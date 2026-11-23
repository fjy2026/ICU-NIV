# c3_ext3_03_decision.R —— ③ 决策分析（DCA / 净获益）
# 策略：常规监测（不强化）vs 对高风险层强化监测
# 这是决策分析，不是影响评估；不能写成“实施后改善结局”。
source((function(f) { for (p in c(file.path("code/00_common", f), file.path("../00_common", f), file.path("00_common", f), f)) if (file.exists(p)) return(p); f })("common_multilevel.R"))
set.seed(SEED)

pred <- readRDS("results/c3_ext2/_data/predictions_calibration.rds")
stopifnot(length(pred$y_e) == length(pred$p_cv), length(pred$y_m) == length(pred$p_refit))

# 预测集：eICU 用 5 折交叉验证预测；MIMIC 用本地重拟合（主）与迁移（敏感性）
sets <- list(
  list(name = "eICU (5-fold CV)", y = pred$y_e, p = pred$p_cv),
  list(name = "MIMIC-IV (local refit)", y = pred$y_m, p = pred$p_refit),
  list(name = "MIMIC-IV (transported eICU model)", y = pred$y_m, p = pred$p_transport)
)

# ------------------------------------------------------------------
# 1. 校准与总体表现
# ------------------------------------------------------------------
cal_tab <- do.call(rbind, lapply(sets, function(s) {
  cm <- calib_metrics(s$y, s$p)
  data.frame(cohort = s$name, n = length(s$y), events = sum(s$y),
             prevalence_pct = 100 * mean(s$y),
             AUC = cm[["auc"]], Brier = cm[["brier"]],
             calibration_slope = cm[["slope"]], CITL = cm[["citl"]])
}))
save_csv(cal_tab, "decision_calibration.csv")
cat("Calibration:\n"); print(cal_tab)

# ------------------------------------------------------------------
# 2. DCA：模型引导 vs 全部强化 vs 不强化
# ------------------------------------------------------------------
pts <- seq(0.05, 0.50, by = 0.01)
dca_list <- lapply(sets, function(s) {
  g <- dca_grid(s$y, s$p, pts)
  g$cohort <- s$name
  g
})
dca <- do.call(rbind, dca_list)
save_csv(dca, "decision_dca_curve.csv")

key <- c(0.10, 0.20, 0.30)
keytab <- do.call(rbind, lapply(sets, function(s) {
  data.frame(cohort = s$name,
             threshold = key,
             model = sapply(key, function(t) net_benefit(s$y, s$p, t)),
             monitor_all = sapply(key, function(t) net_benefit(s$y, rep(1, length(s$y)), t)),
             monitor_none = 0)
}))
keytab$incremental_vs_all <- keytab$model - keytab$monitor_all
save_csv(keytab, "decision_net_benefit_key.csv")
cat("\nNet benefit at key thresholds:\n"); print(keytab)

# ------------------------------------------------------------------
# 3. 固定容量风险分层策略（top 1/3、top 1/4）
# ------------------------------------------------------------------
fc <- do.call(rbind, lapply(sets, function(s) {
  do.call(rbind, lapply(c(1/3, 1/4), function(top) {
    f <- fixed_capacity(s$y, s$p, top)
    f$cohort <- s$name
    f$prevalence_pct <- 100 * mean(s$y)
    f$NB_at_cutoff <- net_benefit(s$y, s$p, f$risk_cutoff)
    f$NB_at_0.20 <- net_benefit(s$y, s$p, 0.20)
    f$NB_monitor_all_at_cutoff <- net_benefit(s$y, rep(1, length(s$y)), f$risk_cutoff)
    f
  }))
}))
fc <- fc[, c("cohort", "top", "n_selected", "selected_pct", "events_captured", "sensitivity",
             "ppv", "npv", "n_monitored_per_case", "risk_cutoff",
             "NB_at_cutoff", "NB_at_0.20", "NB_monitor_all_at_cutoff", "prevalence_pct")]
save_csv(fc, "decision_fixed_capacity.csv")
cat("\nFixed-capacity strategies:\n"); print(fc)

# ------------------------------------------------------------------
# 4. 图：DCA
# ------------------------------------------------------------------
plot_dca <- function(cohort_name) {
  g <- dca[dca$cohort == cohort_name, ]
  f <- fc[fc$cohort == cohort_name & fc$top == 1/3, ]
  ytop <- max(g$model, na.rm = TRUE)
  ggplot(g, aes(threshold)) +
    geom_line(aes(y = model, colour = "Model-guided monitoring"), linewidth = 0.9) +
    geom_line(aes(y = all, colour = "Enhanced monitoring for all"), linewidth = 0.7, linetype = 2) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.5) +
    geom_point(data = f, aes(x = risk_cutoff, y = NB_at_cutoff, shape = "Top tertile"),
               size = 2.6, colour = "#B85042") +
    scale_colour_manual(values = c("Model-guided monitoring" = "#2C6E8F",
                                   "Enhanced monitoring for all" = "grey35"),
                        name = NULL) +
    scale_shape_manual(values = c("Top tertile" = 17), name = NULL) +
    labs(x = "Risk threshold", y = "Net benefit") +
    annotate("text", x = 0.06, y = min(g$all, na.rm = TRUE) * 0.72,
             hjust = 0, vjust = 0, size = 3.0,
             label = sprintf("Top tertile: %.0f%% events\n1 per %.1f monitored",
                             f$sensitivity * 100, f$n_monitored_per_case)) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom",
          plot.title = element_blank())
}

p1 <- plot_dca("eICU (5-fold CV)")
p2 <- plot_dca("MIMIC-IV (local refit)")
p3 <- plot_dca("MIMIC-IV (transported eICU model)")
fig <- (p1 | p2 | p3) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(file.path(OUT, "figures", "fig_decision_dca.png"), fig,
       width = 11, height = 3.8, dpi = 300)
ggsave(file.path(OUT, "figures", "fig_decision_dca.pdf"), fig,
       width = 11, height = 3.8)

saveRDS(list(cal_tab = cal_tab, dca = dca, keytab = keytab, fc = fc, sets = sets),
        file.path(OUT, "_data", "decision.rds"))
cat("\n03 done.\n")
