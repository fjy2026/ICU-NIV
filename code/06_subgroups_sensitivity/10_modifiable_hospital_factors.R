# c3_ext3_02_modifiable.R —— ② 预设可改变候选的增量价值
# 候选一：ΔHACOR（24 h − 2 h），landmark 24 h，仅纳入 24 h 未插管者，结局 24–48 h 转有创
# 候选二：ΔROX（24 h − 2 h），同法
# 候选三：医院 NIV 量 / 教学属性加入主模型后的增量
# 候选四：NIV 参数（压力支持、潮气量、漏气）覆盖检查
source((function(f) { for (p in c(file.path("code/00_common", f), file.path("../00_common", f), file.path("00_common", f), f)) if (file.exists(p)) return(p); f })("common_multilevel.R"))
set.seed(SEED)

L <- load_eicu_ext3()
cl <- L$cl
n <- nrow(cl[[1]])
for (i in seq_along(cl)) {
  cl[[i]]$dhacor <- cl[[i]]$hacor24 - cl[[i]]$hacor2
  cl[[i]]$drox <- cl[[i]]$rox24 - cl[[i]]$rox2
}

# ------------------------------------------------------------------
# 1. ΔHACOR / ΔROX 的 landmark 分析
# ------------------------------------------------------------------
landmark <- cl[[1]]$y24 == 0
cat("Landmark (not intubated by 24 h):", sum(landmark), "admissions;",
    sum(cl[[1]]$to_invasive_48h[landmark]), "events 24-48 h\n")

f_base <- fml("to_invasive_48h", paste(MAIN_TERMS, collapse = " + "))
f_hacor <- fml("to_invasive_48h", paste(c(MAIN_TERMS, "hacor2", "dhacor"), collapse = " + "))
f_rox <- fml("to_invasive_48h", paste(c(MAIN_TERMS, "rox2", "drox"), collapse = " + "))

fit_landmark <- function(var, full_form, subset_extra, label) {
  ok <- landmark & subset_extra
  y <- cl[[1]]$to_invasive_48h[ok]
  cl_sub <- lapply(cl, function(d) d[ok, , drop = FALSE])
  p0 <- fit_glm_pred(cl_sub, f_base)
  p1 <- fit_glm_pred(cl_sub, full_form)
  fits <- lapply(cl_sub, function(d) {
    d$dx_group <- droplevels(d$dx_group)
    tryCatch(glm(full_form, data = d, family = binomial), error = function(e) NULL)
  })
  fits <- fits[!vapply(fits, is.null, TRUE)]
  pool_var <- pool_scalar(fits, var)
  pool_base <- pool_scalar(fits, if (var == "dhacor") "hacor2" else "rox2")
  a0 <- auc_ci_cluster(y, p0, cl[[1]]$hospitalid[ok], B = 500)
  a1 <- auc_ci_cluster(y, p1, cl[[1]]$hospitalid[ok], B = 500)
  d <- auc_diff_cluster(y, p0, p1, cl[[1]]$hospitalid[ok], B = 500)
  list(label = label, n = sum(ok), events = sum(y), coverage = 100 * sum(ok) / n,
       OR = exp(pool_var$est), OR_lo = exp(pool_var$lo), OR_hi = exp(pool_var$hi),
       p = pool_var$p, OR_base = exp(pool_base$est), OR_base_lo = exp(pool_base$lo),
       OR_base_hi = exp(pool_base$hi), p_base = pool_base$p,
       AUC_base = a0[["auc"]], AUC_base_lo = a0[["lo"]], AUC_base_hi = a0[["hi"]],
       AUC_full = a1[["auc"]], AUC_full_lo = a1[["lo"]], AUC_full_hi = a1[["hi"]],
       dAUC = d[["diff"]], dAUC_lo = d[["lo"]], dAUC_hi = d[["hi"]],
       p0 = p0, p1 = p1, y = y, ok = ok)
}

hacor <- fit_landmark("dhacor", f_hacor, !is.na(cl[[1]]$hacor2) & !is.na(cl[[1]]$hacor24),
                      "Delta HACOR (24 h - 2 h), per 1 point")
rox <- fit_landmark("drox", f_rox, !is.na(cl[[1]]$rox2) & !is.na(cl[[1]]$rox24),
                    "Delta ROX (24 h - 2 h), per 1 unit")

# 随访完整子集（敏感性）
hacor_fu <- fit_landmark("dhacor", f_hacor,
  !is.na(cl[[1]]$hacor2) & !is.na(cl[[1]]$hacor24) & cl[[1]]$fu == 1,
  "Delta HACOR, follow-up complete")
rox_fu <- fit_landmark("drox", f_rox,
  !is.na(cl[[1]]$rox2) & !is.na(cl[[1]]$rox24) & cl[[1]]$fu == 1,
  "Delta ROX, follow-up complete")

tab <- data.frame(
  candidate = c(hacor$label, rox$label),
  n = c(hacor$n, rox$n), events = c(hacor$events, rox$events),
  coverage_pct = c(hacor$coverage, rox$coverage),
  OR_per_unit = c(hacor$OR, rox$OR), OR_lo = c(hacor$OR_lo, rox$OR_lo), OR_hi = c(hacor$OR_hi, rox$OR_hi),
  p = c(hacor$p, rox$p),
  OR_baseline = c(hacor$OR_base, rox$OR_base), baseline_lo = c(hacor$OR_base_lo, rox$OR_base_lo),
  baseline_hi = c(hacor$OR_base_hi, rox$OR_base_hi), baseline_p = c(hacor$p_base, rox$p_base),
  AUC_base = c(hacor$AUC_base, rox$AUC_base), AUC_base_lo = c(hacor$AUC_base_lo, rox$AUC_base_lo),
  AUC_base_hi = c(hacor$AUC_base_hi, rox$AUC_base_hi),
  AUC_full = c(hacor$AUC_full, rox$AUC_full), AUC_full_lo = c(hacor$AUC_full_lo, rox$AUC_full_lo),
  AUC_full_hi = c(hacor$AUC_full_hi, rox$AUC_full_hi),
  dAUC = c(hacor$dAUC, rox$dAUC), dAUC_lo = c(hacor$dAUC_lo, rox$dAUC_lo),
  dAUC_hi = c(hacor$dAUC_hi, rox$dAUC_hi)
)
save_csv(tab, "modifiable_hacor_rox.csv")

tab_fu <- data.frame(
  candidate = c(hacor_fu$label, rox_fu$label),
  n = c(hacor_fu$n, rox_fu$n), events = c(hacor_fu$events, rox_fu$events),
  OR_per_unit = c(hacor_fu$OR, rox_fu$OR), OR_lo = c(hacor_fu$OR_lo, rox_fu$OR_lo),
  OR_hi = c(hacor_fu$OR_hi, rox_fu$OR_hi), p = c(hacor_fu$p, rox_fu$p),
  dAUC = c(hacor_fu$dAUC, rox_fu$dAUC), dAUC_lo = c(hacor_fu$dAUC_lo, rox_fu$dAUC_lo),
  dAUC_hi = c(hacor_fu$dAUC_hi, rox_fu$dAUC_hi)
)
save_csv(tab_fu, "modifiable_hacor_rox_followup.csv")

cat("\nDelta HACOR / ROX landmark:\n"); print(tab)
cat("\nFollow-up complete sensitivity:\n"); print(tab_fu)

# ------------------------------------------------------------------
# 2. 医院 NIV 量 / 教学属性的增量
# ------------------------------------------------------------------
f_aug <- fml("to_invasive_48h", paste(c(MAIN_TERMS, "vol50", "teaching"), collapse = " + "))
y <- cl[[1]]$to_invasive_48h
hosp <- cl[[1]]$hospitalid
p_base <- fit_glm_pred(cl, f_base)
p_aug <- fit_glm_pred(cl, f_aug)
a_base <- auc_ci_cluster(y, p_base, hosp, B = 500)
a_aug <- auc_ci_cluster(y, p_aug, hosp, B = 500)
d_aug <- auc_diff_cluster(y, p_base, p_aug, hosp, B = 500)

set.seed(SEED)
nri_vals <- matrix(NA_real_, 500, 3)
cls <- unique(hosp)
for (b in 1:500) {
  s <- sample(cls, length(cls), replace = TRUE)
  idx <- unlist(lapply(s, function(k) which(hosp == k)), use.names = FALSE)
  if (length(unique(y[idx])) < 2) next
  nri_vals[b, ] <- nri_continuous(y[idx], p_base[idx], p_aug[idx])
}
nri_ci <- apply(nri_vals, 2, quantile, c(0.025, 0.975), na.rm = TRUE)
nri_est <- nri_continuous(y, p_base, p_aug)

# 跨医院 5 折交叉验证（沿用 ext2 的均衡分折），比较增量
pred <- readRDS("results/c3_ext2/_data/predictions_calibration.rds")
fold <- pred$fold
stopifnot(length(fold) == n)

cv_pred <- function(form) {
  pm <- matrix(NA_real_, n, length(cl))
  for (i in seq_along(cl)) {
    d <- cl[[i]]
    X <- model.matrix(form, data = d)
    yy <- d$to_invasive_48h
    for (k in 1:5) {
      tr <- fold != k; te <- fold == k
      fit <- glm.fit(X[tr, , drop = FALSE], yy[tr], family = binomial())
      b <- fit$coefficients; b[is.na(b)] <- 0
      pm[te, i] <- plogis(as.numeric(X[te, , drop = FALSE] %*% b))
    }
  }
  rowMeans(pm, na.rm = TRUE)
}
p_cv_base <- cv_pred(f_base)
p_cv_aug <- cv_pred(f_aug)
cv_base <- as.numeric(pROC::auc(pROC::roc(y, p_cv_base, quiet = TRUE)))
cv_aug <- as.numeric(pROC::auc(pROC::roc(y, p_cv_aug, quiet = TRUE)))

vt <- data.frame(
  metric = c("Apparent AUC, base model", "Apparent AUC, + hospital volume and teaching",
             "Apparent delta AUC", "Cross-validated AUC, base model",
             "Cross-validated AUC, + hospital volume and teaching",
             "Cross-validated delta AUC",
             "Continuous NRI (apparent)", "NRI, events", "NRI, non-events"),
  estimate = c(a_base[["auc"]], a_aug[["auc"]], d_aug[["diff"]],
               cv_base, cv_aug, cv_aug - cv_base,
               nri_est[["nri"]], nri_est[["events"]], nri_est[["non_events"]]),
  lo = c(a_base[["lo"]], a_aug[["lo"]], d_aug[["lo"]], NA, NA, NA,
         nri_ci[1, 1], nri_ci[1, 2], nri_ci[1, 3]),
  hi = c(a_base[["hi"]], a_aug[["hi"]], d_aug[["hi"]], NA, NA, NA,
         nri_ci[2, 1], nri_ci[2, 2], nri_ci[2, 3])
)
save_csv(vt, "hospital_volume_teaching_increment.csv")
cat("\nHospital volume / teaching increment:\n"); print(vt)

# ------------------------------------------------------------------
# 3. NIV 参数覆盖检查
# ------------------------------------------------------------------
par_tab <- data.frame(
  variable = c("NIV expiratory positive airway pressure (mean)", "NIV set respiratory rate (mean)",
               "Pressure support", "Tidal volume", "Leak"),
  nonmissing = c(sum(!is.na(cl[[1]]$niv_epap)), sum(!is.na(cl[[1]]$niv_setrr)), NA, NA, NA),
  coverage_pct = c(100 * mean(!is.na(cl[[1]]$niv_epap)), 100 * mean(!is.na(cl[[1]]$niv_setrr)), NA, NA, NA),
  available_in_delivered_data = c("Yes", "Yes", "No", "No", "No")
)
save_csv(par_tab, "niv_parameters_coverage.csv")
cat("\nNIV parameter coverage:\n"); print(par_tab)

saveRDS(list(hacor = hacor, rox = rox, hacor_fu = hacor_fu, rox_fu = rox_fu,
             tab = tab, tab_fu = tab_fu, vt = vt, par_tab = par_tab,
             p_base = p_base, p_aug = p_aug, p_cv_base = p_cv_base, p_cv_aug = p_cv_aug),
        file.path(OUT, "_data", "modifiable.rds"))
cat("\n02 done.\n")
