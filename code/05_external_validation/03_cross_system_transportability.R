# c3_ext2_02_transportability.R —— 可移植性对照表（含 DCA 净获益）
# 输入：results/c3_ext2/_data/predictions_calibration.rds（校准脚本产出的固定预测）
#       以及 eICU/MIMIC 完成数据集（用于补充迁移模型）
# 输出：results/c3_ext2/tables/transportability_summary.csv
# 说明：补充/敏感性分析；不改变主结论。
source((function(f) { for (p in c(file.path("code/00_common", f), file.path("../00_common", f), file.path("00_common", f), f)) if (file.exists(p)) return(p); f })("common_transportability.R"))
set.seed(SEED)
dir.create(file.path(OUT, "tables"), recursive = TRUE, showWarnings = FALSE)

pred <- readRDS(file.path(OUT, "_data", "predictions_calibration.rds"))
e <- load_eicu_ext2(); ei <- subset(e, in_analysis == 1)
mi <- load_mimic_ext2()
cl <- completed_eicu_ext2(ei)
cm <- completed_mimic_ext2(mi)

# ---------------- 补充迁移模型 ----------------
# (1) 去掉严重度评分的 6 变量模型（与既往跨中心分析口径一致）
f6 <- to_invasive_48h ~ age10 + female + gcs3 + sbp10 + rr5 + spo2_5
p6 <- matrix(NA_real_, nrow(mi), length(cl))
for (j in seq_along(cl)) {
  b <- coef(glm(f6, data = cl[[j]], family = binomial))
  Xm <- model.matrix(f6, data = cm[[1]])
  p6[, j] <- plogis(lp_from_coef(Xm, b))
}
p6_m <- rowMeans(p6, na.rm = TRUE)

# (2) 连续变量在各自队列内 z 标准化后的迁移（消除量表/单位差异）
zcols <- c("age10","apache10","gcs3","sbp10","rr5","spo2_5")
zify <- function(d, ref) {
  for (v in zcols) d[[paste0(v,"_z")]] <- (d[[v]] - mean(ref[[v]], na.rm = TRUE)) /
    sd(ref[[v]], na.rm = TRUE)
  d
}
clz <- lapply(cl, zify, ref = ei)
cmz <- lapply(cm, zify, ref = mi)
fz <- to_invasive_48h ~ age10_z + female + dx_group + apache10_z + gcs3_z + sbp10_z + rr5_z + spo2_5_z
pz <- matrix(NA_real_, nrow(mi), length(clz))
for (j in seq_along(clz)) {
  b <- coef(glm(fz, data = clz[[j]], family = binomial))
  Xm <- model.matrix(fz, data = cmz[[1]])
  pz[, j] <- plogis(lp_from_coef(Xm, b))
}
pz_m <- rowMeans(pz, na.rm = TRUE)

# ---------------- 指标与 DCA ----------------
dca3 <- function(y, p, th = c(0.10, 0.20, 0.30)) sapply(th, function(t) net_benefit(y, p, t))
rowf <- function(label, cohort, y, p, note) {
  nb <- dca3(y, p)
  data.frame(model = label, cohort = cohort, n = length(y), events = sum(y),
             event_rate = round(100*mean(y),1),
             auc = round(auc_safe(y,p),3), brier = round(brier_safe(y,p),4),
             calib_slope = round(calib_slope(y,p),3),
             calib_intercept = round(calib_intercept(y,p),3),
             slope_model_intercept = round(calib_intercept_lp(y,p),3),
             nb_0.10 = round(nb[1],4), nb_0.20 = round(nb[2],4), nb_0.30 = round(nb[3],4),
             note = note, stringsAsFactors = FALSE)
}
y_e <- pred$y_e; y_m <- pred$y_m
# MIMIC 可比队列：内科 ICU、非 OR/PACU，且 9 个协变量完整（与稿件 Table S14 同口径）
cov_m <- c("to_invasive_48h","age10","female","dx_group","apache10","gcs3","sbp10","rr5","spo2_5")
cmp <- mi$icu_medical == 1 & mi$flag_or_pacu_source == 0 &
  complete.cases(mi[, cov_m])
cat("MIMIC comparable cohort n =", sum(cmp), " events =", sum(y_m[cmp]), "\n")
core <- rbind(
  rowf("eICU apparent", "eICU", y_e, pred$p_app, "In-sample; 20 imputations averaged"),
  rowf("eICU 5-fold CV", "eICU", y_e, pred$p_cv, "Hospital-level 5-fold; each stay predicted once"),
  rowf("eICU model transported to MIMIC", "MIMIC", y_m, pred$p_transport,
       "Full main model (severity coded APACHE/10 vs SOFA/2)"),
  rowf("MIMIC refit", "MIMIC", y_m, pred$p_refit, "Same covariates refit in MIMIC"),
  rowf("MIMIC comparable cohort: eICU model transported", "MIMIC comparable",
       y_m[cmp], pred$p_transport[cmp], "Medical ICU, non-OR/PACU, complete covariates"),
  rowf("MIMIC comparable cohort: refit", "MIMIC comparable",
       y_m[cmp], pred$p_refit[cmp], "Medical ICU, non-OR/PACU, complete covariates")
)
supp <- rbind(
  rowf("S1: eICU 6-variable model transported to MIMIC", "MIMIC", y_m, p6_m,
       "Prior cross-centre model: severity and diagnostic group removed"),
  rowf("S2: eICU model transported to MIMIC, predictors z-scored within cohort", "MIMIC", y_m, pz_m,
       "Continuous predictors standardised within each cohort")
)
# 参照策略
ref_rows <- rbind(
  data.frame(model = "Treat all", cohort = "eICU", n = length(y_e), events = sum(y_e),
             event_rate = round(100*mean(y_e),1), auc = NA, brier = NA, calib_slope = NA,
             calib_intercept = NA, slope_model_intercept = NA,
             nb_0.10 = round(nb_treat_all(y_e,.10),4), nb_0.20 = round(nb_treat_all(y_e,.20),4),
             nb_0.30 = round(nb_treat_all(y_e,.30),4), note = "Reference strategy", stringsAsFactors = FALSE),
  data.frame(model = "Treat all", cohort = "MIMIC", n = length(y_m), events = sum(y_m),
             event_rate = round(100*mean(y_m),1), auc = NA, brier = NA, calib_slope = NA,
             calib_intercept = NA, slope_model_intercept = NA,
             nb_0.10 = round(nb_treat_all(y_m,.10),4), nb_0.20 = round(nb_treat_all(y_m,.20),4),
             nb_0.30 = round(nb_treat_all(y_m,.30),4), note = "Reference strategy", stringsAsFactors = FALSE),
  data.frame(model = "Treat all", cohort = "MIMIC comparable", n = sum(cmp), events = sum(y_m[cmp]),
             event_rate = round(100*mean(y_m[cmp]),1), auc = NA, brier = NA, calib_slope = NA,
             calib_intercept = NA, slope_model_intercept = NA,
             nb_0.10 = round(nb_treat_all(y_m[cmp],.10),4), nb_0.20 = round(nb_treat_all(y_m[cmp],.20),4),
             nb_0.30 = round(nb_treat_all(y_m[cmp],.30),4), note = "Reference strategy", stringsAsFactors = FALSE)
)
out <- rbind(core, supp, ref_rows)
write.csv(out, file.path(OUT, "tables", "transportability_summary.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("\nTransportability summary:\n"); print(out[, c("model","n","events","auc","brier","calib_slope","calib_intercept","nb_0.10","nb_0.20","nb_0.30")])

# ---------------- DCA 曲线 ----------------
dca_curve <- function(y, p, th) {
  data.frame(threshold = th,
             nb_model = sapply(th, function(t) net_benefit(y, p, t)),
             nb_all = sapply(th, function(t) nb_treat_all(y, t)))
}
th <- seq(0.02, 0.50, by = 0.01)
dA <- rbind(cbind(model = "eICU apparent", dca_curve(y_e, pred$p_app, th)),
            cbind(model = "eICU 5-fold CV", dca_curve(y_e, pred$p_cv, th)),
            cbind(model = "Treat all", data.frame(threshold = th, nb_model =
              sapply(th, function(t) nb_treat_all(y_e, t)), nb_all = NA)))
dB <- rbind(cbind(model = "eICU model transported to MIMIC", dca_curve(y_m, pred$p_transport, th)),
            cbind(model = "MIMIC refit", dca_curve(y_m, pred$p_refit, th)),
            cbind(model = "Treat all", data.frame(threshold = th, nb_model =
              sapply(th, function(t) nb_treat_all(y_m, t)), nb_all = NA)))
dA$model <- factor(dA$model, levels = c("eICU apparent","eICU 5-fold CV","Treat all"))
dB$model <- factor(dB$model, levels = c("eICU model transported to MIMIC","MIMIC refit","Treat all"))
cols <- c("eICU apparent" = "#7FA8C4", "eICU 5-fold CV" = "#2C6E9B",
          "eICU model transported to MIMIC" = "#C0504D", "MIMIC refit" = "#4F7942",
          "Treat all" = "grey25")
theme_set(theme_classic(base_size = 11))
pA <- ggplot(dA[, c("threshold","model","nb_model")], aes(threshold, nb_model, color = model, linetype = model)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = cols, name = NULL) +
  scale_linetype_manual(values = c("solid","solid","dashed"), name = NULL) +
  labs(x = "Threshold probability", y = "Net benefit",
       title = "A  eICU decision curves") +
  coord_cartesian(ylim = c(-0.05, 0.10)) +
  theme(legend.position = c(0.62, 0.80), legend.background = element_blank(),
        legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8))
pB <- ggplot(dB[, c("threshold","model","nb_model")], aes(threshold, nb_model, color = model, linetype = model)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = cols, name = NULL) +
  scale_linetype_manual(values = c("solid","solid","dashed"), name = NULL) +
  labs(x = "Threshold probability", y = "Net benefit",
       title = "B  MIMIC decision curves") +
  coord_cartesian(ylim = c(-0.10, 0.20)) +
  theme(legend.position = c(0.62, 0.80), legend.background = element_blank(),
        legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 7.5))
fig <- pA | pB
ggsave(file.path(OUT, "figures", "FigureS5_dca.png"), fig, width = 10, height = 4.2, dpi = 300)
ggsave(file.path(OUT, "figures", "FigureS5_dca.pdf"), fig, width = 10, height = 4.2)

writeLines(capture.output(sessionInfo()), file.path(OUT, "_data", "session_info_transport.txt"))
cat("\ntransportability done\n")
