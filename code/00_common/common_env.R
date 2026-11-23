# C3 正式分析：公共函数与建模口径
# 运行顺序：c3_01_prep.R -> c3_02_main.R -> c3_03_external.R -> c3_04_sensitivity.R -> c3_05_figures.R
# 所有路径均相对项目根目录 E:/课题9-征稿

suppressPackageStartupMessages({
  library(duckdb)
  library(mice)
  library(geepack)
  library(lme4)
  library(rms)
  library(pROC)
  library(EValue)
  library(sandwich)
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

SEED <- 20261113
DX_LEVELS <- c("COPD/哮喘", "心源性肺水肿", "肺炎/其他呼吸", "脓毒症", "神经系统", "术后/其他")

# 一级协变量（建模用，已按方案缩放）
# 方案 2.7.2：全部一级变量 + 缺失指示变量（APACHE、GCS）
MAIN_TERMS <- c("age10", "female", "dx_group", "apache10", "gcs3", "sbp10", "rr5", "spo2_5",
                "m_apache", "m_gcs")

# 7 个一级协变量（外部验证方向对照用；诊断组因分类不可对齐，不进入方向对照）
DIRECTION_TERMS <- c("age10", "female", "apache10", "gcs3", "sbp10", "rr5", "spo2_5")

# 连续变量缩放与派生（eICU）
prep_eicu <- function(d) {
  d$age10   <- d$age_numeric / 10
  d$female  <- as.integer(d$gender == "Female")
  d$dx_group <- factor(d$dx_group, levels = DX_LEVELS)
  d$apache10 <- d$apachescore / 10
  d$gcs3     <- d$gcs_total / 3
  d$sbp10    <- d$sbp_min / 10
  d$rr5      <- d$rr_max / 5
  d$spo2_5   <- d$spo2_min / 5
  d$m_apache <- as.integer(is.na(d$apachescore))
  d$m_gcs    <- as.integer(is.na(d$gcs_total))
  d$m_ph     <- as.integer(is.na(d$ph_min))
  d$m_paco2  <- as.integer(is.na(d$paco2_max))
  d$m_pf     <- as.integer(is.na(d$pf_min))
  d$m_lactate<- as.integer(is.na(d$lactate_max))
  d$hosp_id  <- factor(d$hospitalid)
  d
}

# 连续变量缩放与派生（MIMIC）
prep_mimic <- function(d) {
  d$age10   <- d$age / 10
  d$female  <- as.integer(d$gender == "F")
  d$dx_group <- factor(d$dx_group, levels = DX_LEVELS)
  d$apache10 <- d$sofa / 2          # SOFA 每 2 分，与 eICU 的「严重度每单位」对齐命名
  d$gcs3     <- d$gcs_min / 3
  d$sbp10    <- d$sbp_min / 10
  d$rr5      <- d$rr_max / 5
  d$spo2_5   <- d$spo2_min / 5
  d$m_apache <- as.integer(is.na(d$sofa))
  d$m_gcs    <- as.integer(is.na(d$gcs_min))
  d$m_ph     <- as.integer(is.na(d$ph_min))
  d$m_paco2  <- as.integer(is.na(d$paco2_max))
  d$m_pf     <- as.integer(is.na(d$pf_min))
  d$m_lactate<- as.integer(is.na(d$lactate_max))
  d
}

# 对一组插补数据集按行子集
subset_list <- function(dl, idx) lapply(dl, function(d) d[idx, , drop = FALSE])

# 主模型公式
main_formula <- function(y = "to_invasive_48h", terms = MAIN_TERMS) {
  as.formula(paste(y, "~", paste(terms, collapse = " + ")))
}

# 聚类稳健 logistic（医院 / 患者层面）
# 点估计与 logistic GEE 相同；稳健方差取 sandwich::vcovCL（HC0），
# 在 independence working correlation 下与 geeglm 的稳健标准误一致。
# 主模型另用 geeglm(exchangeable) 复核，结果见 results/c3/tables/_archive_*。
fit_gee <- function(d, formula, id = "hosp_id") {
  d <- d[complete.cases(d[, all.vars(formula)]), ]
  fit <- glm(formula, data = d, family = binomial)
  V <- sandwich::vcovCL(fit, cluster = d[[id]], type = "HC0")
  structure(list(coefficients = coef(fit), vcov = V), class = "c3fit")
}
coef.c3fit <- function(object, ...) object$coefficients
vcov.c3fit <- function(object, ...) object$vcov

# 普通 logistic（MIMIC 单中心；稳健标准误由 sandwich 提供）
fit_glm <- function(d, formula) {
  d <- d[complete.cases(d[, all.vars(formula)]), ]
  glm(formula, data = d, family = binomial)
}

# Rubin 合并（用于 GEE / glmer 多重插补）
# fits: 每个插补数据集上的模型对象（有 coef 与 vcov）
get_fixef <- function(f) {
  if (inherits(f, "merMod")) lme4::fixef(f) else coef(f)
}

align_fits <- function(fits) {
  bl <- lapply(fits, get_fixef)
  alln <- unique(unlist(lapply(bl, names)))
  b <- sapply(bl, function(x) { v <- rep(NA_real_, length(alln)); names(v) <- alln; v[names(x)] <- x; v })
  if (is.null(dim(b))) b <- matrix(b, ncol = length(fits), dimnames = list(alln, NULL))
  rownames(b) <- alln
  Vl <- lapply(fits, function(f) {
    V <- as.matrix(vcov(f)); out <- matrix(NA_real_, length(alln), length(alln),
                                           dimnames = list(alln, alln))
    out[rownames(V), colnames(V)] <- V; out
  })
  keep <- alln[rowSums(is.na(b)) == 0]
  list(b = b[keep, , drop = FALSE], V = lapply(Vl, function(V) V[keep, keep, drop = FALSE]), names = keep)
}

pool_rubin <- function(fits) {
  m <- length(fits)
  al <- align_fits(fits)
  b <- al$b; keep <- al$names
  V <- al$V
  Ubar <- Reduce(`+`, V) / m
  B <- cov(t(b))
  if (is.null(dim(B))) B <- matrix(B, nrow = nrow(b), ncol = nrow(b))
  Tmat <- Ubar + (1 + 1/m) * B
  se <- sqrt(diag(Tmat))
  lambda <- (1 + 1/m) * diag(B) / diag(Tmat)
  lambda[!is.finite(lambda)] <- 0
  df <- ifelse(lambda > 0, (m - 1) / lambda^2, Inf)
  est <- rowMeans(b)
  tval <- est / se
  p <- 2 * pt(-abs(tval), df = df)
  lo <- est - qt(0.975, df) * se
  hi <- est + qt(0.975, df) * se
  data.frame(term = rownames(b), estimate = est, std.error = se,
             conf.low = lo, conf.high = hi, p.value = p, df = df,
             row.names = NULL)
}

# Rubin 合并，返回完整协方差矩阵（用于联合 Wald 检验）
pool_rubin_full <- function(fits) {
  m <- length(fits)
  al <- align_fits(fits)
  b <- al$b
  Ubar <- Reduce(`+`, al$V) / m
  B <- cov(t(b))
  if (is.null(dim(B))) B <- matrix(B, nrow = nrow(b), ncol = nrow(b))
  Tmat <- Ubar + (1 + 1/m) * B
  list(est = rowMeans(b), T = Tmat, names = rownames(b))
}

# 联合 Wald 检验（H0：指定系数全为 0）
joint_wald <- function(pooled_full, terms) {
  idx <- match(terms, pooled_full$names)
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0) return(list(W = NA, df = 0, p = NA))
  b <- pooled_full$est[idx]
  V <- pooled_full$T[idx, idx, drop = FALSE]
  W <- as.numeric(t(b) %*% solve(V) %*% b)
  df <- length(idx)
  list(W = W, df = df, p = pchisq(W, df, lower.tail = FALSE))
}

# 单个模型提取为统一格式
tidy_fit <- function(fit, model_label = "") {
  b <- get_fixef(fit); V <- as.matrix(vcov(fit))
  ok <- !is.na(b) & rownames(V) %in% names(b)
  b <- b[rownames(V)[ok]]; V <- V[rownames(V)[ok], rownames(V)[ok], drop = FALSE]
  se <- sqrt(diag(V))
  data.frame(model = model_label, term = names(b), estimate = b,
             std.error = se, conf.low = b - 1.96 * se, conf.high = b + 1.96 * se,
             p.value = 2 * pnorm(-abs(b / se)), row.names = NULL)
}

# OR 表（合并后）
or_table <- function(pooled, model_label = "") {
  out <- pooled
  out$model <- model_label
  out$OR <- exp(out$estimate)
  out$OR_low <- exp(out$conf.low)
  out$OR_high <- exp(out$conf.high)
  out[, c("model", "term", "estimate", "std.error", "OR", "OR_low", "OR_high", "p.value")]
}

# 标准化均数差
smd_cont <- function(x, g) {
  a <- x[g == 0]; b <- x[g == 1]
  a <- a[!is.na(a)]; b <- b[!is.na(b)]
  if (length(a) < 2 || length(b) < 2) return(NA_real_)
  (mean(b) - mean(a)) / sqrt((var(a) + var(b)) / 2)
}
smd_bin <- function(p, g) {
  # p: 0/1 指标
  p0 <- mean(p[g == 0], na.rm = TRUE); p1 <- mean(p[g == 1], na.rm = TRUE)
  (p1 - p0) / sqrt((p0 * (1 - p0) + p1 * (1 - p1)) / 2)
}
smd_cat_max <- function(x, g) {
  lev <- levels(factor(x))
  max(sapply(lev, function(l) smd_bin(as.integer(factor(x) == l), g)), na.rm = TRUE)
}

# 基线表：按主结局分组
baseline_table <- function(d, vars_num, vars_cat, group = "to_invasive_48h", labels = NULL) {
  g <- d[[group]]
  ok <- !is.na(g)
  d <- d[ok, ]; g <- g[ok]
  rows <- list()
  rows[[length(rows)+1]] <- data.frame(variable = "N", level = "",
    overall = as.character(nrow(d)), group0 = as.character(sum(g == 0)),
    group1 = as.character(sum(g == 1)), smd = "", missing = "", stringsAsFactors = FALSE)
  for (v in vars_num) {
    x <- d[[v]]
    fmt <- function(z) sprintf("%.1f (%.1f–%.1f)", median(z, na.rm = TRUE),
                               quantile(z, .25, na.rm = TRUE), quantile(z, .75, na.rm = TRUE))
    rows[[length(rows)+1]] <- data.frame(
      variable = v, level = "median (IQR)",
      overall = fmt(x), group0 = fmt(x[g == 0]), group1 = fmt(x[g == 1]),
      smd = sprintf("%.2f", smd_cont(x, g)),
      missing = sprintf("%d (%.1f%%)", sum(is.na(x)), 100*mean(is.na(x))),
      stringsAsFactors = FALSE)
  }
  for (v in vars_cat) {
    x <- factor(d[[v]])
    for (l in levels(x)) {
      ind <- as.integer(x == l)
      rows[[length(rows)+1]] <- data.frame(
        variable = v, level = l,
        overall = sprintf("%d (%.1f%%)", sum(ind == 1, na.rm = TRUE), 100*mean(ind == 1, na.rm = TRUE)),
        group0 = sprintf("%d (%.1f%%)", sum(ind[g==0] == 1, na.rm = TRUE), 100*mean(ind[g==0] == 1, na.rm = TRUE)),
        group1 = sprintf("%d (%.1f%%)", sum(ind[g==1] == 1, na.rm = TRUE), 100*mean(ind[g==1] == 1, na.rm = TRUE)),
        smd = sprintf("%.2f", smd_bin(ind, g)),
        missing = sprintf("%d (%.1f%%)", sum(is.na(d[[v]])), 100*mean(is.na(d[[v]]))),
        stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  if (!is.null(labels)) out$variable <- ifelse(out$variable %in% names(labels), labels[out$variable], out$variable)
  out
}

# Wilson 置信区间
wilson_ci <- function(x, n) {
  if (n == 0) return(c(NA, NA))
  p <- x / n; z <- 1.96
  den <- 1 + z^2/n
  centre <- (p + z^2/(2*n)) / den
  half <- z * sqrt(p*(1-p)/n + z^2/(4*n^2)) / den
  c(centre - half, centre + half)
}

# 从合并系数计算绝对风险（其余协变量取参照/中位数）
risk_at <- function(beta, data_ref, vary, values) {
  out <- sapply(values, function(v) {
    nd <- data_ref; nd[[vary]] <- v
    X <- model.matrix(~ age10 + female + dx_group + apache10 + gcs3 + sbp10 + rr5 + spo2_5,
                      data = nd)
    # 只取截距与已估系数对应的列
    common <- intersect(colnames(X), names(beta))
    lp <- as.numeric(X[, common, drop = FALSE] %*% beta[common])
    mean(plogis(lp))
  })
  data.frame(value = values, risk = out)
}

# 医院层面漏斗图输入：观测率与模型期望率
hospital_observed <- function(d, y = "to_invasive_48h", hosp = "hospitalid") {
  agg <- aggregate(d[[y]], by = list(hospitalid = d[[hosp]]), FUN = function(z)
    c(n = length(z), events = sum(z, na.rm = TRUE)))
  out <- data.frame(hospitalid = agg$hospitalid, n = agg$x[, "n"], events = agg$x[, "events"])
  out$rate <- out$events / out$n
  out
}
