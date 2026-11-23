# c3temporal_02_landmark24.R —— 敏感性 B：24 小时 landmark 分析
# 输入：scripts/c3temporal_common.R 载入的主分析集与时间顺序补取表
# 输出：results/c3_temporal/tableB_landmark24.csv
#       results/c3_temporal/tableB_landmark24_riskset.csv
#
# 设计（先定死，后看结果）：
#   - 风险集：主分析集；排除 24 h 内已转有创者；排除 24 h 内死亡者；
#     主分析另要求 ICU 停留 >= 24 h（24 h 时仍在 ICU，处于风险中）。
#   - 结局：24–48 h 内转有创（NIV 后 1440–2880 分钟），且转有创早于死亡。
#   - 协变量：主模型的 0–24 h 值；在 24 h landmark 时点均已先于结局发生。
#   - 缺失处理：在 landmark 风险集内重新做 20 套链式方程插补、Rubin 合并；
#     因风险集与主分析不同，不复用主分析的插补。
#   - 估计：exchangeable GEE、医院聚类，与主模型一致。
#   - 敏感性风险集：不要求 ICU 停留 >= 24 h，仅排除 24 h 内转有创/死亡者，
#     出院者在 24–48 h 按无事件计。

source((function(f) { for (p in c(file.path("code/00_common", f), file.path("../00_common", f), file.path("00_common", f), f)) if (file.exists(p)) return(p); f })("common_temporal.R"))
set.seed(SEED)

e  <- load_eicu_temporal()
ei <- subset(e, in_analysis == 1)

# 结局与风险集
ei$y24_48 <- as.integer(ei$to_invasive_48h == 1 & ei$to_invasive_24h == 0)
alive24 <- is.na(ei$tmp_death_rel_min) | ei$tmp_death_rel_min > 1440
rs_base <- ei$to_invasive_24h == 0 & alive24
rs_main <- rs_base & ei$icu_los_days >= 1
rs_alt  <- rs_base

# 与 t_inv_rel_min 交叉核对
chk_ev <- ei$y24_48 == as.integer(!is.na(ei$tmp_t_inv_rel_min) &
                                  ei$tmp_t_inv_rel_min > 1440 &
                                  ei$tmp_t_inv_rel_min <= 2880 &
                                  (is.na(ei$tmp_death_rel_min) |
                                   ei$tmp_t_inv_rel_min < ei$tmp_death_rel_min))
cat("outcome cross-check mismatches:", sum(!chk_ev, na.rm = TRUE), "\n")

imp_vars <- c("y24_48","age_numeric","gender","dx_group","apachescore","gcs_total",
              "sbp_min","rr_max","spo2_min","ph_min","paco2_max","pf_min","lactate_max",
              "bicarbonate_first","hosp_nivfirst_n","hosp_ph_avail_pct")

run_landmark <- function(rs, tag, seed_off = 0) {
  d0 <- ei[rs, ]
  dimp <- d0[, imp_vars]
  dimp$y24_48 <- factor(dimp$y24_48, levels = c(0, 1))
  dimp$gender <- factor(dimp$gender)
  meth <- make.method(dimp)
  if ("gender" %in% names(meth)) meth["gender"] <- "logreg"
  set.seed(SEED + seed_off)
  imp <- mice(dimp, m = 20, maxit = 10, method = meth, printFlag = FALSE)
  dl <- lapply(1:imp$m, function(i) {
    d <- complete(imp, i)
    d$y <- as.integer(as.character(d$y24_48))
    d$dx_group <- factor(d$dx_group, levels = DX_LEVELS)
    d$age10 <- d$age_numeric / 10
    d$female <- as.integer(d$gender == "Female")
    d$apache10 <- d$apachescore / 10
    d$gcs3 <- d$gcs_total / 3
    d$sbp10 <- d$sbp_min / 10
    d$rr5 <- d$rr_max / 5
    d$spo2_5 <- d$spo2_min / 5
    d$m_apache <- as.integer(is.na(d0$apachescore))
    d$m_gcs <- as.integer(is.na(d0$gcs_total))
    d$hosp_id <- factor(d0$hospitalid)
    d
  })
  f <- as.formula(paste("y ~", paste(MAIN_TERMS, collapse = " + ")))
  p <- fit_pool_subset(dl, f)$pool
  n <- nrow(d0); ev <- sum(d0$y24_48)
  key <- c("age10","female","apache10","gcs3","sbp10","rr5","spo2_5")
  lab <- c(age10="Age, per 10 years", female="Female vs male",
           apache10="Severity score, per unit", gcs3="GCS, per 3 points",
           sbp10="SBP, per 10 mmHg", rr5="RR, per 5 breaths/min",
           spo2_5="SpO2, per 5%")
  out <- do.call(rbind, lapply(key, function(tt) {
    r <- or_row(p, tt, lab[tt], tag, n, ev)
    if (is.null(r)) return(NULL)
    r$hospitals <- length(unique(d0$hospitalid))
    r
  }))
  out
}

tabB <- run_landmark(rs_main, "Landmark 24 h (ICU stay >= 24 h)", 1)
write.csv(tabB, file.path(OUT, "tableB_landmark24.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

tabB2 <- run_landmark(rs_alt, "Landmark 24 h (no ICU-LOS restriction)", 2)
write.csv(tabB2, file.path(OUT, "tableB_landmark24_riskset.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

cat("\nRisk set flow:\n")
cat("main cohort:", nrow(ei), "events48:", sum(ei$to_invasive_48h), "\n")
cat("excluded intubated <=24h:", sum(ei$to_invasive_24h == 1), "\n")
cat("excluded death <=24h:", sum(ei$to_invasive_24h == 0 & !alive24), "\n")
cat("excluded ICU LOS <24h:", sum(rs_base & ei$icu_los_days < 1), "\n")
cat("landmark n:", sum(rs_main), "events 24-48h:", sum(ei$y24_48[rs_main]), "\n")
cat("alt risk set n:", sum(rs_alt), "events 24-48h:", sum(ei$y24_48[rs_alt]), "\n")
cat("\nB done.\n")
print(tabB)
