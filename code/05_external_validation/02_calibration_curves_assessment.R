# c3_ext2_01_calibration.R —— 跨医院 5 折交叉验证校准 + 外部校准
# 输入：results/c3/_data/{eicu_prepared_v2,imp_eicu_v2}.rds、results/c3_ext/_data/mimic_ext.rds
# 输出：results/c3_ext2/tables/calibration_metrics.csv
#       results/c3_ext2/tables/calibration_deciles.csv
#       results/c3_ext2/figures/FigureS4_calibration.png/.pdf
# 说明：补充/敏感性分析；不改变主结论。
source((function(f) { for (p in c(file.path("code/00_common", f), file.path("../00_common", f), file.path("00_common", f), f)) if (file.exists(p)) return(p); f })("common_transportability.R"))
set.seed(SEED)
dir.create(file.path(OUT, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "figures"), recursive = TRUE, showWarnings = FALSE)

e <- load_eicu_ext2()
ei <- subset(e, in_analysis == 1)
cl <- completed_eicu_ext2(ei)
cat("eICU n =", nrow(ei), " events =", sum(ei$to_invasive_48h),
    " hospitals =", length(unique(ei$hospitalid)), "\n")

f_main <- main_formula()

# ---------------- 医院层面 5 折 ----------------
hn <- data.frame(hospitalid = unique(ei$hospitalid))
tab <- table(ei$hospitalid)
hn$n <- as.integer(tab[as.character(hn$hospitalid)])
hn$events <- as.integer(tapply(ei$to_invasive_48h, ei$hospitalid, sum)[as.character(hn$hospitalid)])
set.seed(SEED)
hn <- hn[sample(nrow(hn)), ]
hn$fold <- 0L; tot <- rep(0, 5)
for (i in seq_len(nrow(hn))) {
  f <- which.min(tot); hn$fold[i] <- f; tot[f] <- tot[f] + hn$n[i]
}
fold_map <- setNames(hn$fold, as.character(hn$hospitalid))
ei$fold <- as.integer(fold_map[as.character(ei$hospitalid)])
cat("\n5-fold (hospital-level) sizes and events:\n")
print(data.frame(fold = 1:5,
                 n = as.integer(table(ei$fold)),
                 events = as.integer(tapply(ei$to_invasive_48h, ei$fold, sum)),
                 hospitals = as.integer(table(hn$fold))))

# ---------------- 表观与交叉验证预测 ----------------
# 预测模型与主模型协变量相同，用 logistic 回归（与既往跨中心分析一致；
# GEE 只改变标准误、几乎不改变点估计，且 geeglm 无稳定的 newdata 预测接口）。
fit_pred <- function(d, formula) glm(formula, data = d, family = binomial)
p_app <- matrix(NA_real_, nrow(ei), length(cl))
p_cv  <- matrix(NA_real_, nrow(ei), length(cl))
for (j in seq_along(cl)) {
  d <- cl[[j]]; d$fold <- ei$fold
  X <- lp_matrix(d, f_main)
  fit_all <- fit_pred(d, f_main)
  p_app[, j] <- plogis(lp_from_coef(X, coef(fit_all)))
  for (k in 1:5) {
    tr <- which(d$fold != k); te <- which(d$fold == k)
    fit_k <- try(fit_pred(d[tr, ], f_main), silent = TRUE)
    if (inherits(fit_k, "try-error")) { cat("fold", k, "imp", j, "failed\n"); next }
    b <- coef(fit_k)
    if (any(is.na(b))) cat("fold", k, "imp", j, "has NA coefficients\n")
    p_cv[te, j] <- plogis(lp_from_coef(X[te, , drop = FALSE], b))
  }
}
p_app_e <- rowMeans(p_app, na.rm = TRUE)
p_cv_e  <- rowMeans(p_cv, na.rm = TRUE)
cat("\nCV predictions missing:", sum(is.na(p_cv_e)), "\n")

y_e <- ei$to_invasive_48h

# ---------------- MIMIC：迁移与重拟合 ----------------
mi <- load_mimic_ext2()
cm <- completed_mimic_ext2(mi)
cat("MIMIC n =", nrow(mi), " events =", sum(mi$to_invasive_48h), "\n")

f_mimic <- main_formula(terms = setdiff(MAIN_TERMS, c("m_apache","m_gcs")))
p_tr <- matrix(NA_real_, nrow(mi), length(cl))   # eICU 模型迁移
p_re <- matrix(NA_real_, nrow(mi), length(cm))   # MIMIC 重拟合
# 迁移：eICU 每套插补模型系数 -> MIMIC 协变量
for (j in seq_along(cl)) {
  d <- cl[[j]]
  b <- coef(fit_pred(d, f_main))
  dm <- cm[[1]]  # 协变量取第一套 MIMIC 插补（MIMIC 缺失极少）；结局未被使用
  Xm <- lp_matrix(dm, f_main)
  p_tr[, j] <- plogis(lp_from_coef(Xm, b))
}
# 重拟合：MIMIC 自身模型
for (j in seq_along(cm)) {
  d <- cm[[j]]
  X <- lp_matrix(d, f_mimic)
  p_re[, j] <- plogis(lp_from_coef(X, coef(fit_pred(d, f_mimic))))
}
p_tr_m <- rowMeans(p_tr, na.rm = TRUE)
p_re_m <- rowMeans(p_re, na.rm = TRUE)
y_m <- mi$to_invasive_48h

saveRDS(list(y_e = y_e, p_app = p_app_e, p_cv = p_cv_e,
             fold = ei$fold, hospital = ei$hospitalid,
             y_m = y_m, p_transport = p_tr_m, p_refit = p_re_m,
             subject = mi$subject_id),
        file.path(OUT, "_data", "predictions_calibration.rds"))

# ---------------- 指标 ----------------
metric_row <- function(label, cohort, y, p, cluster) {
  bs <- boot_metrics(y, p, cluster, B = 500)
  q <- function(x) quantile(x, c(.025, .975), na.rm = TRUE)
  data.frame(
    model = label, cohort = cohort, n = length(y), events = sum(y),
    event_rate = round(100*mean(y), 1),
    auc = round(auc_safe(y,p), 3), auc_low = round(q(bs$auc)[1],3), auc_high = round(q(bs$auc)[2],3),
    brier = round(brier_safe(y,p), 4), brier_low = round(q(bs$brier)[1],4), brier_high = round(q(bs$brier)[2],4),
    calib_slope = round(calib_slope(y,p), 3),
    slope_low = round(q(bs$slope)[1],3), slope_high = round(q(bs$slope)[2],3),
    calib_intercept = round(calib_intercept(y,p), 3),
    intercept_low = round(q(bs$intercept)[1],3), intercept_high = round(q(bs$intercept)[2],3),
    slope_model_intercept = round(calib_intercept_lp(y,p), 3),
    stringsAsFactors = FALSE)
}

met <- rbind(
  metric_row("eICU apparent", "eICU", y_e, p_app_e, ei$hospitalid),
  metric_row("eICU 5-fold CV", "eICU", y_e, p_cv_e, ei$hospitalid),
  metric_row("eICU model transported to MIMIC", "MIMIC", y_m, p_tr_m, mi$subject_id),
  metric_row("MIMIC refit", "MIMIC", y_m, p_re_m, mi$subject_id)
)
met$note <- c("In-sample predictions, 20 imputations averaged",
              "Hospital-level 5-fold; each stay predicted once; 20 imputations averaged",
              "eICU coefficients applied to MIMIC covariates (severity coded as SOFA/2)",
              "Same covariates refit in MIMIC; patient-clustered")
write.csv(met, file.path(OUT, "tables", "calibration_metrics.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("\nCalibration metrics:\n"); print(met[, c("model","n","events","auc","brier","calib_slope","calib_intercept")])

# 逐折指标（检查中心间校准漂移）
fold_met <- do.call(rbind, lapply(1:5, function(k) {
  idx <- ei$fold == k
  data.frame(fold = k, n = sum(idx), events = sum(y_e[idx]),
             hospitals = length(unique(ei$hospitalid[idx])),
             auc = round(auc_safe(y_e[idx], p_cv_e[idx]), 3),
             brier = round(brier_safe(y_e[idx], p_cv_e[idx]), 4),
             calib_slope = round(calib_slope(y_e[idx], p_cv_e[idx]), 3),
             calib_intercept = round(calib_intercept(y_e[idx], p_cv_e[idx]), 3))
}))
write.csv(fold_met, file.path(OUT, "tables", "calibration_folds.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("\nPer-fold calibration:\n"); print(fold_met)

# 十分位校准点
dp <- rbind(
  cbind(model = "eICU 5-fold CV", decile_points(y_e, p_cv_e)),
  cbind(model = "eICU model transported to MIMIC", decile_points(y_m, p_tr_m)),
  cbind(model = "MIMIC refit", decile_points(y_m, p_re_m))
)
write.csv(dp, file.path(OUT, "tables", "calibration_deciles.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# ---------------- 图（十分位校准曲线；不画 loess，避免稀疏高预测端的过冲） ----------------
dp$model <- factor(dp$model, levels = c("eICU 5-fold CV", "eICU model transported to MIMIC",
                                        "MIMIC refit"))
lim <- c(0, 0.65)
theme_set(theme_classic(base_size = 11))

pA <- ggplot(subset(dp, model == "eICU 5-fold CV")) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey45") +
  geom_line(aes(mean_pred, obs), color = "#2C6E9B", linewidth = 0.7) +
  geom_pointrange(aes(mean_pred, obs, ymin = lo, ymax = hi), color = "#2C6E9B", size = 0.35) +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = "Predicted probability of invasive ventilation within 48 h",
       y = "Observed proportion", title = "A  eICU, hospital-level 5-fold cross-validation") +
  annotate("text", x = 0.03, y = 0.60, hjust = 0,
           label = sprintf("Calibration slope %.2f\nCalibration intercept %.2f\nBrier score %.3f",
                           met$calib_slope[2], met$calib_intercept[2], met$brier[2]), size = 3.2)

pB <- ggplot(subset(dp, model != "eICU 5-fold CV")) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey45") +
  geom_line(aes(mean_pred, obs, color = model), linewidth = 0.7) +
  geom_pointrange(aes(mean_pred, obs, ymin = lo, ymax = hi, color = model), size = 0.32) +
  scale_color_manual(values = c("eICU model transported to MIMIC" = "#C0504D",
                                "MIMIC refit" = "#4F7942"), name = NULL,
                     labels = c(sprintf("Transported (intercept %.2f)", abs(round(met$calib_intercept[3], 2))),
                                sprintf("Refit (intercept %.2f)", abs(round(met$calib_intercept[4], 2))))) +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = "Predicted probability of invasive ventilation within 48 h",
       y = "Observed proportion", title = "B  MIMIC: transported model vs refit") +
  theme(legend.position = c(0.46, 0.10), legend.background = element_blank(),
        legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8))

fig <- pA | pB
ggsave(file.path(OUT, "figures", "FigureS4_calibration.png"), fig,
       width = 11, height = 5.0, dpi = 300)
ggsave(file.path(OUT, "figures", "FigureS4_calibration.pdf"), fig,
       width = 11, height = 5.0)

writeLines(capture.output(sessionInfo()), file.path(OUT, "_data", "session_info_calibration.txt"))
cat("\ncalibration done\n")
