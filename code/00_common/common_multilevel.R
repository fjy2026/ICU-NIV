# c3_ext3_common.R —— ext3 补强分析公共函数
# 只读：data/processed/c3_expanded、results/c3/_data、results/c3_ext/_data、results/c3_ext2/_data
# 全部新产物写入 results/c3_ext3/
# 性质：补充/探索性分析；不改变主分析、主结论与手稿。
suppressPackageStartupMessages({
  library(duckdb)
  library(mice)
  library(lme4)
  library(pROC)
  library(sandwich)
  library(metafor)
  library(ggplot2)
  library(patchwork)
})

SEED <- 20261117
OUT <- "results/c3_ext3"
dir.create(file.path(OUT, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "_data"), recursive = TRUE, showWarnings = FALSE)

DX_LEVELS <- c("COPD/哮喘", "心源性肺水肿", "肺炎/其他呼吸", "脓毒症", "神经系统", "术后/其他")
MAIN_TERMS <- c("age10", "female", "dx_group", "apache10", "gcs3", "sbp10", "rr5", "spo2_5",
                "m_apache", "m_gcs")
DEATH_TERMS <- c("age10", "female", "dx_group", "apache10", "gcs3", "sbp10", "dnr")
SER_LEVELS <- c("Low", "Mid", "High")
BEDS_LEVELS <- c("<100", "100 - 249", "250 - 499", ">= 500", "Unknown")

TILDE <- intToUtf8(126)
fml <- function(lhs, rhs) as.formula(paste(lhs, TILDE, rhs))

# ---------------- eICU：数据与插补完成数据集 ----------------
load_eicu_ext3 <- function() {
  e <- readRDS("results/c3/_data/eicu_prepared_v2.rds")
  ei <- subset(e, in_analysis == 1)
  imp <- readRDS("results/c3/_data/imp_eicu_v2.rds")
  stopifnot(nrow(imp$data) == nrow(ei), imp$m == 20)

  con <- dbConnect(duckdb(shared_home = FALSE))
  hx <- dbGetQuery(con, "SELECT patientunitstayid, hospitalid, hosp_beds_category, hosp_teaching,
    hosp_nivfirst_n, hosp_dnr_rate, death_hosp, dnr_before_niv, dni_before_niv,
    hacor2h_total, hacor24h_total, rox_2h, rox_24h, to_invasive_24h, to_invasive_48h,
    followup_48h_complete, niv_epap_mean, niv_set_rr_mean, unittype, icu_medical, icu_surgical,
    flag_or_pacu_source
    FROM read_parquet('data/processed/c3_expanded/eicu_niv_expanded.parquet')")
  dbDisconnect(con, shutdown = TRUE)
  stopifnot(!anyDuplicated(hx$patientunitstayid), all(ei$patientunitstayid %in% hx$patientunitstayid))
  hx <- hx[match(ei$patientunitstayid, hx$patientunitstayid), ]
  stopifnot(identical(as.numeric(ei$patientunitstayid), as.numeric(hx$patientunitstayid)))

  beds <- ifelse(is.na(hx$hosp_beds_category), "Unknown", as.character(hx$hosp_beds_category))

  prep <- function(d) {
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
    d$hospitalid <- ei$hospitalid
    d$hosp <- factor(ei$hospitalid)
    d$death_hosp <- hx$death_hosp
    d$dnr <- hx$dnr_before_niv
    d$dnr_known <- as.integer(!is.na(d$dnr))
    d$dnr0 <- ifelse(is.na(d$dnr), 0L, as.integer(d$dnr == 1))
    d$dni <- hx$dni_before_niv
    d$teaching <- as.integer(hx$hosp_teaching)
    d$beds <- factor(beds, levels = BEDS_LEVELS)
    d$vol <- hx$hosp_nivfirst_n
    d$vol50 <- d$vol / 50
    d$dnr_rate <- hx$hosp_dnr_rate
    d$hacor2 <- hx$hacor2h_total
    d$hacor24 <- hx$hacor24h_total
    d$rox2 <- hx$rox_2h
    d$rox24 <- hx$rox_24h
    d$y24 <- ei$to_invasive_24h
    d$fu <- ei$followup_48h_complete
    d$icu_medical <- hx$icu_medical
    d$icu_surgical <- hx$icu_surgical
    d$or_pacu <- hx$flag_or_pacu_source
    d$niv_epap <- hx$niv_epap_mean
    d$niv_setrr <- hx$niv_set_rr_mean
    d
  }

  cl <- lapply(1:imp$m, function(i) {
    d <- complete(imp, i)
    stopifnot(all(as.character(d$dx_group) == as.character(ei$dx_group)))
    prep(d)
  })

  hosp <- data.frame(hospitalid = hx$hospitalid, teaching = as.integer(hx$hosp_teaching),
                     beds = factor(beds, levels = BEDS_LEVELS),
                     vol = hx$hosp_nivfirst_n, dnr_rate = hx$hosp_dnr_rate)
  hosp <- hosp[!duplicated(hosp$hospitalid), ]
  rownames(hosp) <- NULL
  list(ei = ei, imp = imp, cl = cl, hosp = hosp, hx = hx)
}

# ---------------- 通用：固定效应与方差 ----------------
get_fixef <- function(f) if (inherits(f, "merMod")) lme4::fixef(f) else stats::coef(f)
get_vcov <- function(f) as.matrix(stats::vcov(f))

# Rubin 合并单参数
pool_scalar <- function(fits, term) {
  b <- vapply(fits, function(f) { x <- get_fixef(f); if (term %in% names(x)) x[[term]] else NA_real_ }, 0)
  v <- vapply(fits, function(f) { V <- get_vcov(f); if (term %in% rownames(V)) V[term, term] else NA_real_ }, 0)
  ok <- is.finite(b) & is.finite(v) & v > 0
  if (sum(ok) < 3) return(data.frame(term = term, est = NA, se = NA, df = NA, p = NA, lo = NA, hi = NA, m = sum(ok)))
  b <- b[ok]; v <- v[ok]; m <- length(b)
  qbar <- mean(b); ubar <- mean(v); bvar <- stats::var(b)
  tv <- ubar + (1 + 1 / m) * bvar
  df <- if (bvar > 0) (m - 1) * (1 + ubar / ((1 + 1 / m) * bvar))^2 else Inf
  se <- sqrt(tv)
  crit <- if (is.finite(df)) stats::qt(0.975, df) else 1.959964
  p <- if (is.finite(df)) 2 * stats::pt(-abs(qbar / se), df) else 2 * stats::pnorm(-abs(qbar / se))
  data.frame(term = term, est = qbar, se = se, df = df, p = p,
             lo = qbar - crit * se, hi = qbar + crit * se, m = m)
}

# Rubin 合并多参数：联合 Wald 检验
pool_joint <- function(fits, terms) {
  terms <- terms[vapply(terms, function(t)
    all(vapply(fits, function(f) t %in% names(get_fixef(f)), TRUE)), TRUE)]
  if (length(terms) == 0) return(list(Qbar = NA, Tv = NA, W = NA, df = 0, p = NA, terms = terms))
  Q <- sapply(fits, function(f) { x <- get_fixef(f); x[terms] })
  if (is.null(dim(Q))) Q <- matrix(Q, nrow = length(terms))
  U <- lapply(fits, function(f) { V <- get_vcov(f); V[terms, terms, drop = FALSE] })
  Qbar <- rowMeans(Q)
  Ubar <- Reduce(`+`, U) / length(U)
  B <- stats::cov(t(Q))
  Tv <- Ubar + (1 + 1 / length(fits)) * B
  W <- as.numeric(t(Qbar) %*% solve(Tv) %*% Qbar)
  df <- nrow(Tv)
  p <- stats::pchisq(W, df = df, lower.tail = FALSE)
  list(Qbar = Qbar, Tv = Tv, W = W, df = df, p = p, terms = terms)
}

# Rubin 合并线性组合 L %*% beta
pool_lincom <- function(fits, L) {
  terms <- names(L)
  if (!all(vapply(terms, function(t) all(vapply(fits, function(f) t %in% names(get_fixef(f)), TRUE)), TRUE))) {
    return(data.frame(est = NA, se = NA, df = NA, p = NA, lo = NA, hi = NA, m = length(fits)))
  }
  Q <- sapply(fits, function(f) { x <- get_fixef(f); x[terms] })
  if (is.null(dim(Q))) Q <- matrix(Q, nrow = length(terms))
  U <- lapply(fits, function(f) { V <- get_vcov(f); V[terms, terms, drop = FALSE] })
  Qbar <- rowMeans(Q)
  Ubar <- Reduce(`+`, U) / length(U)
  B <- stats::cov(t(Q))
  Tv <- Ubar + (1 + 1 / length(fits)) * B
  est <- as.numeric(sum(L * Qbar))
  se <- sqrt(as.numeric(t(L) %*% Tv %*% L))
  bvar <- as.numeric(t(L) %*% B %*% L)
  ubar <- as.numeric(t(L) %*% Ubar %*% L)
  m <- length(fits)
  df <- if (bvar > 0) (m - 1) * (1 + ubar / ((1 + 1 / m) * bvar))^2 else Inf
  crit <- if (is.finite(df)) stats::qt(0.975, df) else 1.959964
  p <- if (is.finite(df)) 2 * stats::pt(-abs(est / se), df) else 2 * stats::pnorm(-abs(est / se))
  data.frame(est = est, se = se, df = df, p = p, lo = est - crit * se, hi = est + crit * se, m = m)
}

# ---------------- 预测、AUC、校准 ----------------
fit_glm_pred <- function(cl, form, subset = NULL, newdata = NULL) {
  ps <- matrix(NA_real_, nrow(if (is.null(newdata)) cl[[1]] else newdata), length(cl))
  for (i in seq_along(cl)) {
    d <- cl[[i]]
    if (!is.null(subset)) d <- d[subset, , drop = FALSE]
    if ("dx_group" %in% names(d)) d$dx_group <- droplevels(d$dx_group)
    f <- glm(form, data = d, family = binomial)
    nd <- if (is.null(newdata)) d else newdata
    ps[, i] <- predict(f, newdata = nd, type = "response")
  }
  rowMeans(ps, na.rm = TRUE)
}

auc_ci_cluster <- function(y, p, cluster, B = 500, seed = SEED) {
  set.seed(seed)
  ok <- is.finite(y) & is.finite(p)
  y <- y[ok]; p <- p[ok]; cluster <- cluster[ok]
  r <- as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE)))
  cl <- unique(cluster)
  vals <- numeric(B)
  for (b in 1:B) {
    s <- sample(cl, length(cl), replace = TRUE)
    idx <- unlist(lapply(s, function(k) which(cluster == k)), use.names = FALSE)
    if (length(unique(y[idx])) < 2) { vals[b] <- NA; next }
    vals[b] <- as.numeric(pROC::auc(pROC::roc(y[idx], p[idx], quiet = TRUE)))
  }
  ci <- stats::quantile(vals, c(0.025, 0.975), na.rm = TRUE)
  c(auc = r, lo = ci[[1]], hi = ci[[2]])
}

auc_diff_cluster <- function(y, p0, p1, cluster, B = 500, seed = SEED) {
  set.seed(seed)
  ok <- is.finite(y) & is.finite(p0) & is.finite(p1)
  y <- y[ok]; p0 <- p0[ok]; p1 <- p1[ok]; cluster <- cluster[ok]
  r <- as.numeric(pROC::auc(pROC::roc(y, p1, quiet = TRUE))) -
       as.numeric(pROC::auc(pROC::roc(y, p0, quiet = TRUE)))
  cl <- unique(cluster)
  vals <- numeric(B)
  for (b in 1:B) {
    s <- sample(cl, length(cl), replace = TRUE)
    idx <- unlist(lapply(s, function(k) which(cluster == k)), use.names = FALSE)
    if (length(unique(y[idx])) < 2) { vals[b] <- NA; next }
    vals[b] <- as.numeric(pROC::auc(pROC::roc(y[idx], p1[idx], quiet = TRUE))) -
               as.numeric(pROC::auc(pROC::roc(y[idx], p0[idx], quiet = TRUE)))
  }
  ci <- stats::quantile(vals, c(0.025, 0.975), na.rm = TRUE)
  c(diff = r, lo = ci[[1]], hi = ci[[2]])
}

calib_metrics <- function(y, p) {
  ok <- is.finite(y) & is.finite(p)
  y <- y[ok]; p <- p[ok]
  lp <- stats::qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  slope <- unname(stats::coef(stats::glm(y ~ lp, family = binomial()))[2])
  citl <- unname(stats::coef(stats::glm(y ~ 1, offset = lp, family = binomial()))[1])
  c(auc = as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE))),
    brier = mean((p - y)^2), slope = slope, citl = citl)
}

nri_continuous <- function(y, p0, p1) {
  ok <- is.finite(y) & is.finite(p0) & is.finite(p1)
  y <- y[ok]; p0 <- p0[ok]; p1 <- p1[ok]
  up <- p1 > p0; dn <- p1 < p0
  e <- sum((up & y == 1) - (dn & y == 1)) / sum(y == 1)
  ne <- sum((dn & y == 0) - (up & y == 0)) / sum(y == 0)
  c(nri = e + ne, events = e, non_events = ne)
}

idi_fun <- function(y, p0, p1) {
  ok <- is.finite(y) & is.finite(p0) & is.finite(p1)
  y <- y[ok]; p0 <- p0[ok]; p1 <- p1[ok]
  (mean(p1[y == 1]) - mean(p0[y == 1])) - (mean(p1[y == 0]) - mean(p0[y == 0]))
}

# ---------------- 决策曲线 ----------------
net_benefit <- function(y, p, pt) {
  n <- length(y)
  tp <- sum(p >= pt & y == 1); fp <- sum(p >= pt & y == 0)
  tp / n - fp / n * (pt / (1 - pt))
}

dca_grid <- function(y, p, pts = seq(0.05, 0.50, by = 0.01)) {
  data.frame(
    threshold = pts,
    model = vapply(pts, function(t) net_benefit(y, p, t), 0),
    all = vapply(pts, function(t) net_benefit(y, rep(1, length(y)), t), 0),
    none = 0
  )
}

# 固定容量策略：取预测风险最高的 top 比例
fixed_capacity <- function(y, p, top = 0.25) {
  n <- length(y); k <- ceiling(n * top)
  sel <- order(p, decreasing = TRUE)[1:k]
  tp <- sum(y[sel] == 1); fp <- sum(y[sel] == 0)
  data.frame(top = top, n_selected = k, selected_pct = 100 * k / n,
             events_captured = tp, sensitivity = tp / sum(y == 1),
             ppv = tp / k, npv = sum(y[-sel] == 0) / (n - k),
             n_monitored_per_case = ifelse(tp > 0, k / tp, NA),
             risk_cutoff = min(p[sel]))
}

# 泊松精确 CI
pois_ci <- function(o, e) {
  lo <- ifelse(o > 0, stats::qgamma(0.025, o) / e, 0)
  hi <- stats::qgamma(0.975, o + 1) / e
  cbind(lo = lo, hi = hi)
}

save_csv <- function(x, name) {
  utils::write.csv(x, file.path(OUT, "tables", name), row.names = FALSE, fileEncoding = "UTF-8")
}
