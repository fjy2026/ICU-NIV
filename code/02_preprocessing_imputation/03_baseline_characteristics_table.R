# c3_09_baseline_v2.R —— 诊断映射更正后重建基线表
# 输入：results/c3/_data/{eicu,mimic}_prepared_v2.rds
# 输出：results/c3/tables/table1_eicu_baseline.csv
#       results/c3/tables/table1_mimic_baseline.csv
# 说明：变量、分组、SMD 与缺失口径与 c3_01_prep.R 完全一致，只把 dx_group 换成 v2。
# 运行：Rscript scripts/c3_09_baseline_v2.R

source((function(f) { for (p in c(file.path("code/00_common", f), file.path("../00_common", f), file.path("00_common", f), f)) if (file.exists(p)) return(p); f })("common_env.R"))

e <- readRDS("results/c3/_data/eicu_prepared_v2.rds")
m <- readRDS("results/c3/_data/mimic_prepared_v2.rds")
ei <- subset(e, in_analysis == 1)
mi <- subset(m, in_analysis == 1)

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

cat("baseline v2 done\n")
print(t1e[t1e$variable == "入院诊断组", ])
print(t1m[t1m$variable == "入院诊断组", ])
