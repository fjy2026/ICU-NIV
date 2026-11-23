# c3_01_prep.R —— 读取正式数据、构建建模变量、输出队列流程与基线表
# 输入：data/processed/c3/eicu_niv_analysis.parquet, mimic_niv_analysis.parquet
# 输出：results/c3/_data/{eicu,mimic}_prepared.rds
#       results/c3/tables/table1_eicu_baseline.csv
#       results/c3/tables/table1_mimic_baseline.csv
#       results/c3/tables/cohort_flow.csv
#       results/c3/tables/event_rates.csv
#       results/c3/tables/missingness.csv
#       results/c3/tables/hospital_summary.csv
# 运行：Rscript scripts/c3_01_prep.R

source((function(f) { for (p in c(file.path("code/00_common", f), file.path("../00_common", f), file.path("00_common", f), f)) if (file.exists(p)) return(p); f })("common_env.R"))
set.seed(SEED)
dir.create("results/c3/_data", recursive = TRUE, showWarnings = FALSE)
dir.create("results/c3/tables", recursive = TRUE, showWarnings = FALSE)

con <- dbConnect(duckdb(shared_home = TRUE))
e_raw <- dbGetQuery(con, "SELECT * FROM read_parquet('data/processed/c3/eicu_niv_analysis.parquet')")
m_raw <- dbGetQuery(con, "SELECT * FROM read_parquet('data/processed/c3/mimic_niv_analysis.parquet')")
DBI::dbDisconnect(con, shutdown = TRUE)

e <- prep_eicu(e_raw)
m <- prep_mimic(m_raw)

saveRDS(e, "results/c3/_data/eicu_prepared.rds")
saveRDS(m, "results/c3/_data/mimic_prepared.rds")

ei <- subset(e, in_analysis == 1)
mi <- subset(m, in_analysis == 1)

# ---- 队列流程（全部由数据计算） ----
e_n1 <- nrow(e)
e_n2 <- sum(e$inv_before_niv == 0, na.rm = TRUE)
e_n3 <- sum(e$inv_before_niv == 0 & e$inv_same_time == 0, na.rm = TRUE)
e_n4 <- sum(e$in_analysis == 1, na.rm = TRUE)
m_n1 <- nrow(m)
m_n2 <- sum(m$inv_before_niv == 0, na.rm = TRUE)
m_n3 <- sum(m$inv_before_niv == 0 & m$inv_same_time == 0, na.rm = TRUE)
m_n4 <- sum(m$in_analysis == 1, na.rm = TRUE)

flow <- data.frame(
  database = c(rep("eICU", 4), rep("MIMIC", 4)),
  step = c("treatment 文本含 non-invasive ventilation",
           "排除严格有创记录早于无创（序贯撤机）",
           "排除有创与无创开始同刻",
           "排除年龄 < 18 岁或年龄缺失（主分析集）",
           "procedureevents itemid = 225794",
           "排除有创（225792）早于无创",
           "排除有创与无创开始同刻",
           "年龄 < 18 岁（外部分析集）"),
  remaining = c(e_n1, e_n2, e_n3, e_n4, m_n1, m_n2, m_n3, m_n4),
  excluded = c(NA, e_n1 - e_n2, e_n2 - e_n3, e_n3 - e_n4,
               NA, m_n1 - m_n2, m_n2 - m_n3, m_n3 - m_n4)
)
write.csv(flow, "results/c3/tables/cohort_flow.csv", row.names = FALSE, fileEncoding = "UTF-8")

# ---- 事件率（Wilson 95% CI） ----
rate_rows <- list()
add_rate <- function(db, cohort, label, x, n) {
  ci <- wilson_ci(x, n)
  data.frame(database = db, cohort = cohort, outcome = label,
             events = x, n = n, rate_pct = round(100*x/n, 1),
             ci_low = round(100*ci[1], 1), ci_high = round(100*ci[2], 1))
}
for (spec in list(
  list(db="eICU", d=ei, cohort="主分析集"),
  list(db="eICU", d=e, cohort="全无创队列"),
  list(db="MIMIC", d=mi, cohort="外部分析集"),
  list(db="MIMIC", d=m, cohort="全无创队列"))) {
  d <- spec$d; n <- nrow(d)
  rate_rows[[length(rate_rows)+1]] <- add_rate(spec$db, spec$cohort, "48h 内转有创（严格）", sum(d$to_invasive_48h, na.rm=TRUE), n)
  rate_rows[[length(rate_rows)+1]] <- add_rate(spec$db, spec$cohort, "24h 内转有创", sum(d$to_invasive_24h, na.rm=TRUE), n)
  rate_rows[[length(rate_rows)+1]] <- add_rate(spec$db, spec$cohort, "任意时间转有创", sum(d$to_invasive_any, na.rm=TRUE), n)
  if ("to_invasive_48h_anydef" %in% names(d))
    rate_rows[[length(rate_rows)+1]] <- add_rate(spec$db, spec$cohort, "48h 内转有创（宽口径）", sum(d$to_invasive_48h_anydef, na.rm=TRUE), n)
  rate_rows[[length(rate_rows)+1]] <- add_rate(spec$db, spec$cohort, "48h 内死亡", sum(d$death_48h, na.rm=TRUE), n)
  rate_rows[[length(rate_rows)+1]] <- add_rate(spec$db, spec$cohort, "院内死亡", sum(d$death_hosp %in% 1), n)
  rate_rows[[length(rate_rows)+1]] <- add_rate(spec$db, spec$cohort, "复合：48h 转有创或院内死亡", sum(d$composite_48h %in% 1), n)
}
rates <- do.call(rbind, rate_rows)
write.csv(rates, "results/c3/tables/event_rates.csv", row.names = FALSE, fileEncoding = "UTF-8")

# ---- 基线表 ----
num_e <- c("age_numeric","apachescore","gcs_total","sbp_min","rr_max","spo2_min",
           "ph_min","paco2_max","pf_min","lactate_max","bicarbonate_first","hacor_score")
cat_e <- c("gender","dx_group")
labels_e <- c(age_numeric="年龄（岁）", apachescore="APACHE IVa 评分", gcs_total="GCS 总分",
  sbp_min="收缩压最低（mmHg）", rr_max="呼吸频率最高（次/分）", spo2_min="SpO2 最低（%）",
  ph_min="pH 最低", paco2_max="PaCO2 最高（mmHg）", pf_min="PaO2/FiO2 最低",
  lactate_max="乳酸最高（mmol/L）", bicarbonate_first="碳酸氢根首次（mmol/L）", hacor_score="HACOR 评分",
  gender="性别", dx_group="入院诊断组")
t1e <- baseline_table(ei, num_e, cat_e, labels = labels_e)
t1e$database <- "eICU"
write.csv(t1e, "results/c3/tables/table1_eicu_baseline.csv", row.names = FALSE, fileEncoding = "UTF-8")

num_m <- c("age","sofa","gcs_min","sbp_min","rr_max","spo2_min",
           "ph_min","paco2_max","pf_min","lactate_max","bicarbonate_min","hacor_score")
cat_m <- c("gender","dx_group")
labels_m <- c(age="年龄（岁）", sofa="首日 SOFA", gcs_min="GCS 最低", sbp_min="收缩压最低（mmHg）",
  rr_max="呼吸频率最高（次/分）", spo2_min="SpO2 最低（%）", ph_min="pH 最低",
  paco2_max="PaCO2 最高（mmHg）", pf_min="PaO2/FiO2 最低", lactate_max="乳酸最高（mmol/L）",
  bicarbonate_min="碳酸氢根最低（mmol/L）", hacor_score="HACOR 评分", gender="性别", dx_group="入院诊断组")
t1m <- baseline_table(mi, num_m, cat_m, labels = labels_m)
t1m$database <- "MIMIC"
write.csv(t1m, "results/c3/tables/table1_mimic_baseline.csv", row.names = FALSE, fileEncoding = "UTF-8")

# ---- 缺失表 ----
miss_rows <- list()
for (spec in list(list(db="eICU", d=ei), list(db="MIMIC", d=mi))) {
  d <- spec$d
  vars <- setdiff(names(d), c("patientunitstayid","stay_id","hadm_id","subject_id","uniquepid",
                              "patienthealthsystemstayid","hospitalid","hosp_id","niv_start",
                              "niv_start_offset","inv_start_strict","inv_start_any","death_offset",
                              "death_time","icu_los_days","hosp_los_days","apacheadmissiondx"))
  for (v in vars) {
    if (!is.numeric(d[[v]]) && !is.character(d[[v]]) && !is.factor(d[[v]])) next
    miss_rows[[length(miss_rows)+1]] <- data.frame(
      database = spec$db, variable = v, n_missing = sum(is.na(d[[v]])),
      pct_missing = round(100*mean(is.na(d[[v]])), 1))
  }
}
write.csv(do.call(rbind, miss_rows), "results/c3/tables/missingness.csv", row.names = FALSE, fileEncoding = "UTF-8")

# ---- 医院层面摘要 ----
hosp <- hospital_observed(ei)
hosp$ph_avail <- sapply(hosp$hospitalid, function(h) mean(!is.na(ei$ph_min[ei$hospitalid == h])))
hosp_summary <- data.frame(
  metric = c("医院数", "每院 NIV-first 例数（中位数）", "每院例数 IQR 下限", "每院例数 IQR 上限",
             "每院例数范围", "医院 48h 转有创率中位数（%）", "转有创率 IQR（%）",
             "转有创率范围（%）", "pH 可得比例中位数（%）", "pH 可得比例范围（%）",
             "例数 < 50 的医院数", "例数 < 25 的医院数"),
  value = c(nrow(hosp),
            sprintf("%.1f", median(hosp$n)),
            sprintf("%.0f", quantile(hosp$n, .25)),
            sprintf("%.0f", quantile(hosp$n, .75)),
            sprintf("%d–%d", min(hosp$n), max(hosp$n)),
            sprintf("%.1f", 100*median(hosp$rate)),
            sprintf("%.1f–%.1f", 100*quantile(hosp$rate, .25), 100*quantile(hosp$rate, .75)),
            sprintf("%.1f–%.1f", 100*min(hosp$rate), 100*max(hosp$rate)),
            sprintf("%.1f", 100*median(hosp$ph_avail)),
            sprintf("%.0f–%.0f", 100*min(hosp$ph_avail), 100*max(hosp$ph_avail)),
            sum(hosp$n < 50), sum(hosp$n < 25))
)
write.csv(hosp_summary, "results/c3/tables/hospital_summary.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(hosp, "results/c3/tables/hospital_observed_rates.csv", row.names = FALSE, fileEncoding = "UTF-8")

cat("prep done\n")
cat("eICU main n =", nrow(ei), " events =", sum(ei$to_invasive_48h), "\n")
cat("MIMIC main n =", nrow(mi), " events =", sum(mi$to_invasive_48h), "\n")
print(hosp_summary)
