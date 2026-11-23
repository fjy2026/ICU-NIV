# c3_12_hospital_v2.R —— v2 诊断映射下的医院间异质性（ICC/MOR）
# 输入：results/c3/_data/imp_eicu_v2.rds、eicu_prepared_v2.rds
# 输出：results/c3/tables/hospital_heterogeneity.csv（v2）
# 运行：Rscript scripts/c3_12_hospital_v2.R
# 说明：模型与 c3_02_main.R 的随机截距模型完全一致，只换 v2 诊断映射。

source((function(f) { for (p in c(file.path("code/00_common", f), file.path("../00_common", f), file.path("00_common", f), f)) if (file.exists(p)) return(p); f })("common_env.R"))
set.seed(SEED)

e <- readRDS("results/c3/_data/eicu_prepared_v2.rds")
ei <- subset(e, in_analysis == 1)
imp <- readRDS("results/c3/_data/imp_eicu_v2.rds")

complete_list <- lapply(1:imp$m, function(i) {
  d <- complete(imp, i)
  d$to_invasive_48h <- as.integer(as.character(d$to_invasive_48h))
  d$dx_group <- factor(d$dx_group, levels = DX_LEVELS)
  d$age10 <- d$age_numeric/10; d$female <- as.integer(d$gender == "Female")
  d$apache10 <- d$apachescore/10; d$gcs3 <- d$gcs_total/3
  d$sbp10 <- d$sbp_min/10; d$rr5 <- d$rr_max/5; d$spo2_5 <- d$spo2_min/5
  d$m_apache <- as.integer(is.na(ei$apachescore)); d$m_gcs <- as.integer(is.na(ei$gcs_total))
  d$hosp_id <- factor(ei$hospitalid)
  d
})

f_ri <- as.formula(paste("to_invasive_48h ~", paste(MAIN_TERMS, collapse = " + "), "+ (1 | hosp_id)"))
fits_ri <- lapply(complete_list, function(d) glmer(f_ri, data = d, family = binomial, nAGQ = 0))
vars <- sapply(fits_ri, function(f) as.numeric(VarCorr(f)$hosp_id))
v <- mean(vars)
icc <- v / (v + pi^2/3)
mor <- exp(sqrt(2*v) * qnorm(0.75))
# IOR-85 已删除：原 exp(sqrt(2*v)*qnorm(0.85)) 并非标准 IOR（标准 IOR 为覆盖中间 80% 的区间），
# 纯院间异质性以 ICC 与 MOR 汇总。见 results/c3/ior85_reconciliation.md（2026-11-14）。

out <- data.frame(
  metric = c("随机截距方差", "ICC", "MOR（中位比值比）",
             "主模型参数数", "分析样本量", "主结局事件数", "医院数"),
  value = c(sprintf("%.3f", v), sprintf("%.3f", icc), sprintf("%.2f", mor),
            ncol(model.matrix(~ age10 + female + dx_group + apache10 + gcs3 + sbp10 + rr5 + spo2_5 +
                                m_apache + m_gcs, data = complete_list[[1]])),
            nrow(complete_list[[1]]),
            sum(complete_list[[1]]$to_invasive_48h), length(unique(complete_list[[1]]$hosp_id))))
write.csv(out, "results/c3/tables/hospital_heterogeneity.csv", row.names = FALSE, fileEncoding = "UTF-8")
cat("ICC =", round(icc, 3), " MOR =", round(mor, 2), "\n")
print(out)
