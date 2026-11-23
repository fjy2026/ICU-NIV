# c3_ext2_common.R —— 补强分析（ext2）公共函数
# 只读：data/processed/c3_expanded、results/c3/_data、results/c3_ext/_data、MIMIC 源库（只读）
# 所有新产物写入 results/c3_ext2/
# 用途：校准曲线、可移植性、HACOR NRI/IDI、时间零点敏感性、变量对齐
suppressPackageStartupMessages({
  library(duckdb)
  library(mice)
  library(geepack)
  library(sandwich)
  library(pROC)
  library(ggplot2)
  library(patchwork)
})

SEED <- 20261117
DX_LEVELS <- c("COPD/哮喘", "心源性肺水肿", "肺炎/其他呼吸", "脓毒症", "神经系统", "术后/其他")
MAIN_TERMS <- c("age10","female","dx_group","apache10","gcs3","sbp10","rr5","spo2_5",
                "m_apache","m_gcs")
DIRECTION_TERMS <- c("age10","female","apache10","gcs3","sbp10","rr5","spo2_5")
DIR_LABELS <- c(age10="Age (per 10 years)", female="Female (vs male)",
                apache10="Severity score (per unit)", gcs3="GCS (per 3 points)",
                sbp10="Systolic BP (per 10 mmHg)", rr5="Respiratory rate (per 5/min)",
                spo2_5="SpO2 (per 5%)")
OUT <- "results/c3_ext2"

# ---------------- eICU：准备数据 + 时间零点 ----------------
load_eicu_ext2 <- function() {
  e <- readRDS("results/c3/_data/eicu_prepared_v2.rds")
  con <- dbConnect(duckdb())
  x <- dbGetQuery(con, "SELECT patientunitstayid, niv_start_offset, niv_pre_icu,
    inv_start_strict, inv_start_any, death_offset, hospitaladmitoffset,
    followup_48h_complete, icu_los_days, first_icu_stay
    FROM read_parquet('data/processed/c3_expanded/eicu_niv_expanded.parquet')")
  dbDisconnect(con, shutdown = TRUE)
  stopifnot(!anyDuplicated(x$patientunitstayid), all(e$patientunitstayid %in% x$patientunitstayid))
  x <- x[match(e$patientunitstayid, x$patientunitstayid), ]
  stopifnot(identical(as.numeric(e$patientunitstayid), as.numeric(x$patientunitstayid)))
  # 主结局核对：严格有创开始相对无创开始 0–2880 分钟
  chk <- as.integer(!is.na(x$inv_start_strict - x$niv_start_offset) &
                      (x$inv_start_strict - x$niv_start_offset) >= 0 &
                      (x$inv_start_strict - x$niv_start_offset) <= 48*60 &
                      (is.na(x$death_offset) | x$inv_start_strict < x$death_offset))
  ina <- e$in_analysis == 1
  cat("eICU main outcome reproduced (inv-niv in [0,2880], main analysis set): mismatches =",
      sum(chk[ina] != e$to_invasive_48h[ina]), "\n")
  e$niv_start_offset <- x$niv_start_offset
  e$niv_pre_icu <- x$niv_pre_icu
  e$inv_start_strict <- x$inv_start_strict
  e$death_offset <- x$death_offset
  e$followup_48h_complete <- x$followup_48h_complete
  e$icu_los_days <- x$icu_los_days
  # T0 = ICU 入科：严格有创在入 ICU 后 0–48 h 内开始，且早于死亡（与主结局口径一致）
  e$y_t0_icu <- as.integer(!is.na(e$inv_start_strict) & e$inv_start_strict > 0 &
                             e$inv_start_strict <= 48*60 &
                             (is.na(e$death_offset) | e$inv_start_strict < e$death_offset))
  # 限制在 NIV 于入 ICU 48 h 内开始者（避免结构性 0）
  e$y_t0_icu_restricted <- ifelse(!is.na(e$niv_start_offset) & e$niv_start_offset <= 48*60,
                                  e$y_t0_icu, NA_integer_)
  e
}

# 20 套插补完成数据集（复用主分析 v2 插补，保证与主模型一致）
completed_eicu_ext2 <- function(ei) {
  imp <- readRDS("results/c3/_data/imp_eicu_v2.rds")
  stopifnot(nrow(imp$data) == nrow(ei), imp$m == 20)
  d1 <- complete(imp, 1)
  stopifnot(all(as.character(d1$dx_group) == as.character(ei$dx_group)),
            all(d1$to_invasive_48h == ei$to_invasive_48h))
  lapply(1:imp$m, function(i) {
    d <- complete(imp, i)
    d$to_invasive_48h <- as.integer(as.character(d$to_invasive_48h))
    d$dx_group <- factor(ei$dx_group, levels = DX_LEVELS)
    d$age10 <- d$age_numeric / 10
    d$female <- as.integer(d$gender == "Female")
    d$apache10 <- d$apachescore / 10
    d$gcs3 <- d$gcs_total / 3
    d$sbp10 <- d$sbp_min / 10
    d$rr5 <- d$rr_max / 5
    d$spo2_5 <- d$spo2_min / 5
    d$m_apache <- as.integer(is.na(ei$apachescore))
    d$m_gcs <- as.integer(is.na(ei$gcs_total))
    d$hosp_id <- factor(ei$hospitalid)
    d$hospitalid <- ei$hospitalid
    d$patientunitstayid <- ei$patientunitstayid
    d$y_t0_icu <- ei$y_t0_icu
    d$y_t0_icu_restricted <- ei$y_t0_icu_restricted
    d$niv_start_offset <- ei$niv_start_offset
    d
  })
}

# ---------------- MIMIC：准备数据 + 源库时间零点 ----------------
load_mimic_ext2 <- function() {
  m <- readRDS("results/c3_ext/_data/mimic_ext.rds")
  mi <- subset(m, in_analysis == 1)
  mi$dx_group <- factor(as.character(mi$dx_group_v2), levels = DX_LEVELS)
  mi$female <- as.integer(mi$gender == "F")
  mi$year_index <- as.integer(format(as.POSIXct(mi$niv_start), "%Y")) -
    min(as.integer(format(as.POSIXct(mi$niv_start), "%Y")))
  cu <- as.character(mi$first_careunit)
  cu2 <- ifelse(grepl("Medical/Surgical", cu), "MICU/SICU",
         ifelse(grepl("Medical Intensive", cu), "MICU",
         ifelse(grepl("Surgical|Trauma", cu), "SICU",
         ifelse(grepl("Cardiac Vascular", cu), "CVICU",
         ifelse(grepl("Coronary", cu), "CCU", "Other")))))
  mi$careunit <- factor(cu2)
  mi$pat_id <- factor(mi$subject_id)

  cache <- file.path(OUT, "_data", "mimic_t0_timing.rds")
  if (file.exists(cache)) {
    t0 <- readRDS(cache)
  } else {
    con <- dbConnect(duckdb(), dbdir = "F:/CoredataCenter/MIMIC-IV/mimic.duckdb", read_only = TRUE)
    q <- paste(sprintf("%d", mi$stay_id), collapse = ",")
    ic <- dbGetQuery(con, sprintf("SELECT stay_id, intime FROM mimiciv_icu.icustays WHERE stay_id IN (%s)", q))
    pe <- dbGetQuery(con, sprintf("SELECT stay_id, starttime FROM mimiciv_icu.procedureevents
      WHERE itemid = 225792 AND stay_id IN (%s)", q))
    dbDisconnect(con, shutdown = TRUE)
    pe$starttime <- as.POSIXct(pe$starttime, tz = "UTC")
    agg <- aggregate(starttime ~ stay_id, data = pe, FUN = min)
    names(agg)[2] <- "inv_start"
    t0 <- merge(ic, agg, by = "stay_id", all.x = TRUE)
    saveRDS(t0, cache)
  }
  t0 <- t0[match(mi$stay_id, t0$stay_id), ]
  stopifnot(all(t0$stay_id == mi$stay_id))
  mi$intime <- t0$intime
  mi$inv_start <- t0$inv_start
  mi$niv_after_intime_h <- as.numeric(difftime(mi$niv_start, mi$intime, units = "hours"))
  mi$inv_after_intime_h <- as.numeric(difftime(mi$inv_start, mi$intime, units = "hours"))
  mi$inv_after_niv_h <- as.numeric(difftime(mi$inv_start, mi$niv_start, units = "hours"))
  # 核对：主结局 = 严格有创相对无创 0–48 h
  chk <- as.integer(!is.na(mi$inv_after_niv_h) & mi$inv_after_niv_h >= 0 & mi$inv_after_niv_h <= 48 &
                      (is.na(mi$death_time) | mi$inv_start < mi$death_time))
  cat("MIMIC main outcome reproduced from source: mismatches =", sum(chk != mi$to_invasive_48h), "\n")
  mi$y_t0_icu <- as.integer(!is.na(mi$inv_after_intime_h) & mi$inv_after_intime_h > 0 &
                              mi$inv_after_intime_h <= 48 &
                              (is.na(mi$death_time) | mi$inv_start < mi$death_time))
  mi$y_t0_icu_restricted <- ifelse(!is.na(mi$niv_after_intime_h) & mi$niv_after_intime_h <= 48,
                                   mi$y_t0_icu, NA_integer_)
  mi
}

completed_mimic_ext2 <- function(mi) {
  imp_vars <- c("to_invasive_48h","age","gender","dx_group","sofa","gcs_min","sbp_min","rr_max",
                "spo2_min","ph_min","paco2_max","pf_min","lactate_max","bicarbonate_min",
                "year_index","careunit")
  dimp <- mi[, imp_vars]
  dimp$to_invasive_48h <- factor(dimp$to_invasive_48h, levels = c(0,1))
  dimp$gender <- factor(dimp$gender)
  meth <- make.method(dimp)
  if ("gender" %in% names(meth)) meth["gender"] <- "logreg"
  cache <- file.path(OUT, "_data", "imp_mimic_ext2.rds")
  if (file.exists(cache)) {
    imp <- readRDS(cache)
  } else {
    imp <- mice(dimp, m = 20, maxit = 10, seed = SEED, method = meth, printFlag = FALSE)
    saveRDS(imp, cache)
  }
  lapply(1:imp$m, function(i) {
    d <- complete(imp, i)
    d$to_invasive_48h <- as.integer(as.character(d$to_invasive_48h))
    # dx_group 用插补后的值（MIMIC 有 1 例主诊断缺失，主分析由插补处理）
    d$dx_group <- factor(as.character(d$dx_group), levels = DX_LEVELS)
    d$age10 <- d$age / 10
    d$female <- as.integer(d$gender == "F")
    d$apache10 <- d$sofa / 2
    d$gcs3 <- d$gcs_min / 3
    d$sbp10 <- d$sbp_min / 10
    d$rr5 <- d$rr_max / 5
    d$spo2_5 <- d$spo2_min / 5
    d$m_apache <- as.integer(is.na(mi$sofa))
    d$m_gcs <- as.integer(is.na(mi$gcs_min))
    d$pat_id <- factor(mi$subject_id)
    d$stay_id <- mi$stay_id
    d$y_t0_icu <- mi$y_t0_icu
    d$y_t0_icu_restricted <- mi$y_t0_icu_restricted
    d$niv_after_intime_h <- mi$niv_after_intime_h
    d
  })
}

# ---------------- 建模与预测 ----------------
main_formula <- function(y = "to_invasive_48h", terms = MAIN_TERMS) {
  as.formula(paste(y, "~", paste(terms, collapse = " + ")))
}

fit_exch <- function(d, formula, id = "hosp_id") {
  dd <- d[complete.cases(d[, all.vars(formula)]), ]
  dd$.gee_id <- dd[[id]]
  geeglm(formula, data = dd, family = binomial(), id = .gee_id,
         corstr = "exchangeable", std.err = "san.se")
}

fit_cluster_glm <- function(d, formula, id = "hosp_id") {
  dd <- d[complete.cases(d[, all.vars(formula)]), ]
  fit <- glm(formula, data = dd, family = binomial)
  structure(list(coefficients = coef(fit),
                 vcov = sandwich::vcovCL(fit, cluster = dd[[id]], type = "HC0")),
            class = "c3cl")
}
coef.c3cl <- function(object, ...) object$coefficients
vcov.c3cl <- function(object, ...) object$vcov

# 用固定公式在设计矩阵上算线性预测值，避免 predict.geeglm 的 newdata 兼容问题
lp_matrix <- function(d, formula) {
  X <- model.matrix(formula, data = d)
  X
}
lp_from_coef <- function(X, b) {
  common <- intersect(colnames(X), names(b))
  as.numeric(X[, common, drop = FALSE] %*% b[common])
}

rubin_pool <- function(fits) {
  # fits: list of objects with coef() and vcov()
  b <- sapply(fits, coef)
  V <- lapply(fits, vcov)
  alln <- rownames(b)
  m <- length(fits)
  Ubar <- Reduce(`+`, V) / m
  B <- cov(t(b))
  Tmat <- Ubar + (1 + 1/m) * B
  se <- sqrt(diag(Tmat))
  est <- rowMeans(b)
  df <- (m - 1) / pmax((1 + 1/m) * diag(B) / diag(Tmat), 1e-12)^2
  p <- 2 * pt(-abs(est/se), df = df)
  data.frame(term = alln, estimate = est, std.error = se,
             conf.low = est - qt(.975, df) * se, conf.high = est + qt(.975, df) * se,
             p.value = p, row.names = NULL)
}

# ---------------- 预测性能 ----------------
auc_safe <- function(y, p) {
  ok <- !is.na(y) & !is.na(p)
  if (length(unique(y[ok])) < 2) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(y[ok], p[ok], quiet = TRUE, direction = "<")))
}

brier_safe <- function(y, p) mean((p - y)^2, na.rm = TRUE)

calib_slope <- function(y, p) {
  lp <- qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  fit <- glm(y ~ lp, family = binomial)
  unname(coef(fit)[2])
}
calib_intercept <- function(y, p) {
  # calibration-in-the-large（CITL）：offset 为预测 logit
  lp <- qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  fit <- glm(y ~ 1, offset = lp, family = binomial)
  unname(coef(fit)[1])
}
calib_intercept_lp <- function(y, p) {
  # 斜率模型 y ~ a + b*logit(p) 的截距 a（与既往补充表口径一致，非 CITL）
  lp <- qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  unname(coef(glm(y ~ lp, family = binomial))[1])
}

net_benefit <- function(y, p, pt) {
  n <- length(y)
  tp <- sum(p >= pt & y == 1); fp <- sum(p >= pt & y == 0)
  tp/n - fp/n * (pt/(1 - pt))
}
nb_treat_all <- function(y, pt) {
  prev <- mean(y); prev - (1 - prev) * pt/(1 - pt)
}

wilson_ci <- function(x, n) {
  if (n == 0) return(c(NA, NA))
  p <- x/n; z <- 1.96; den <- 1 + z^2/n
  ctr <- (p + z^2/(2*n))/den
  hw <- z * sqrt(p*(1-p)/n + z^2/(4*n^2))/den
  c(ctr - hw, ctr + hw)
}

# 按风险十分位聚合（用于校准点）
decile_points <- function(y, p, k = 10) {
  qs <- unique(quantile(p, probs = seq(0, 1, length.out = k + 1), na.rm = TRUE))
  grp <- cut(p, breaks = qs, include.lowest = TRUE, labels = FALSE)
  out <- data.frame()
  for (g in sort(unique(grp))) {
    idx <- which(grp == g)
    if (length(idx) < 10) next
    ci <- wilson_ci(sum(y[idx]), length(idx))
    out <- rbind(out, data.frame(mean_pred = mean(p[idx]), obs = mean(y[idx]),
                                 n = length(idx), lo = ci[1], hi = ci[2]))
  }
  out
}

# 聚类 bootstrap（固定预测，只重抽样行）
boot_metrics <- function(y, p, cluster, B = 500, seed = SEED) {
  set.seed(seed)
  cl <- unique(cluster); out <- matrix(NA_real_, B, 4,
    dimnames = list(NULL, c("auc","brier","slope","intercept")))
  for (b in 1:B) {
    samp <- sample(cl, length(cl), replace = TRUE)
    idx <- unlist(lapply(samp, function(z) which(cluster == z)))
    if (length(unique(y[idx])) < 2) next
    out[b, ] <- c(auc_safe(y[idx], p[idx]), brier_safe(y[idx], p[idx]),
                  tryCatch(calib_slope(y[idx], p[idx]), error = function(e) NA),
                  tryCatch(calib_intercept(y[idx], p[idx]), error = function(e) NA))
  }
  as.data.frame(out)
}

# ---------------- NRI / IDI ----------------
continuous_nri <- function(y, p_new, p_ref) {
  up_new <- mean(p_new[y == 1] > p_ref[y == 1])
  down_new <- mean(p_new[y == 1] < p_ref[y == 1])
  down_ref <- mean(p_new[y == 0] < p_ref[y == 0])
  up_ref <- mean(p_new[y == 0] > p_ref[y == 0])
  (up_new - down_new) + (down_ref - up_ref)
}

category_nri <- function(y, p_new, p_ref, cuts = c(0.10, 0.20)) {
  g_new <- as.integer(cut(p_new, c(-Inf, cuts, Inf)))
  g_ref <- as.integer(cut(p_ref, c(-Inf, cuts, Inf)))
  up <- g_new > g_ref; down <- g_new < g_ref
  (mean(up[y == 1]) - mean(down[y == 1])) + (mean(down[y == 0]) - mean(up[y == 0]))
}

idi_calc <- function(y, p_new, p_ref) {
  (mean(p_new[y == 1]) - mean(p_ref[y == 1])) +
    (mean(p_ref[y == 0]) - mean(p_new[y == 0]))
}

# 合并统计量打印
fmt_or <- function(est, lo, hi) sprintf("%.2f (%.2f–%.2f)", exp(est), exp(lo), exp(hi))
