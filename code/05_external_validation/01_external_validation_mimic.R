# c3_11_external_v2.R —— 诊断映射更正后的 MIMIC 外部验证（v2）
# 输入：results/c3/_data/mimic_prepared_v2.rds、eicu_prepared_v2.rds、imp_eicu_v2.rds
# 输出：results/c3/_data/imp_mimic_v2.rds
#       results/c3/tables/mimic_main_model.csv（v2）
#       results/c3/tables/table4_direction.csv（v2；eICU 侧改用最终 GEE）
# 运行：Rscript scripts/c3_11_external_v2.R
# 说明：只把 dx_group 换成 v2 映射、把 eICU 侧从聚类稳健 logistic 换成最终
#       exchangeable GEE；其余变量、模型形式、判定门槛与 c3_03_external.R 完全一致。

source((function(f) { for (p in c(file.path("code/00_common", f), file.path("../00_common", f), file.path("00_common", f), f)) if (file.exists(p)) return(p); f })("common_env.R"))
set.seed(SEED)

m <- readRDS("results/c3/_data/mimic_prepared_v2.rds")
mi <- subset(m, in_analysis == 1)
mi$year_index <- as.integer(format(as.POSIXct(mi$niv_start), "%Y")) - min(as.integer(format(as.POSIXct(mi$niv_start), "%Y")))
cu <- as.character(mi$first_careunit)
cu2 <- ifelse(grepl("Medical/Surgical", cu), "MICU/SICU",
       ifelse(grepl("Medical Intensive", cu), "MICU",
       ifelse(grepl("Surgical|Trauma", cu), "SICU",
       ifelse(grepl("Cardiac Vascular", cu), "CVICU",
       ifelse(grepl("Coronary", cu), "CCU", "Other")))))
mi$careunit <- factor(cu2)

imp_vars <- c("to_invasive_48h","age","gender","dx_group","sofa","gcs_min","sbp_min","rr_max",
              "spo2_min","ph_min","paco2_max","pf_min","lactate_max","bicarbonate_min",
              "year_index","careunit")
dimp <- mi[, imp_vars]
dimp$to_invasive_48h <- factor(dimp$to_invasive_48h, levels = c(0, 1))
dimp$gender <- factor(dimp$gender)
meth <- make.method(dimp)
if ("gender" %in% names(meth)) meth["gender"] <- "logreg"

if (file.exists("results/c3/_data/imp_mimic_v2.rds")) {
  imp <- readRDS("results/c3/_data/imp_mimic_v2.rds")
} else {
  imp <- mice(dimp, m = 20, maxit = 10, seed = SEED, method = meth, printFlag = FALSE)
  saveRDS(imp, "results/c3/_data/imp_mimic_v2.rds")
}

complete_list <- lapply(1:imp$m, function(i) {
  d <- complete(imp, i)
  d$to_invasive_48h <- as.integer(as.character(d$to_invasive_48h))
  d$dx_group <- factor(d$dx_group, levels = DX_LEVELS)
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
  d
})

f_main <- main_formula(terms = setdiff(MAIN_TERMS, c("m_apache","m_gcs")))
fits <- lapply(complete_list, function(d) fit_gee(d, f_main, id = "pat_id"))
pool_m <- pool_rubin(fits)
tab_m <- or_table(pool_m, model_label = "MIMIC adjusted")
tab_m$OR <- round(tab_m$OR, 2); tab_m$OR_low <- round(tab_m$OR_low, 2)
tab_m$OR_high <- round(tab_m$OR_high, 2); tab_m$p.value <- signif(tab_m$p.value, 3)
write.csv(tab_m, "results/c3/tables/mimic_main_model.csv", row.names = FALSE, fileEncoding = "UTF-8")

# eICU 侧：最终 exchangeable GEE（v2 插补）
imp_e <- readRDS("results/c3/_data/imp_eicu_v2.rds")
e <- readRDS("results/c3/_data/eicu_prepared_v2.rds")
ei <- subset(e, in_analysis == 1)
ce <- lapply(1:imp_e$m, function(i) {
  d <- complete(imp_e, i)
  d$to_invasive_48h <- as.integer(as.character(d$to_invasive_48h))
  d$dx_group <- factor(d$dx_group, levels = DX_LEVELS)
  d$age10 <- d$age_numeric / 10; d$female <- as.integer(d$gender == "Female")
  d$apache10 <- d$apachescore / 10; d$gcs3 <- d$gcs_total / 3
  d$sbp10 <- d$sbp_min / 10; d$rr5 <- d$rr_max / 5; d$spo2_5 <- d$spo2_min / 5
  d$m_apache <- as.integer(is.na(ei$apachescore)); d$m_gcs <- as.integer(is.na(ei$gcs_total))
  d$hosp_id <- factor(ei$hospitalid); d
})
fit_exch <- function(d, formula) {
  dd <- d[complete.cases(d[, all.vars(formula)]), ]
  dd$.gee_id <- dd$hosp_id
  geeglm(formula, data = dd, family = binomial(), id = .gee_id,
         corstr = "exchangeable", std.err = "san.se")
}
pe <- pool_rubin(lapply(ce, function(d) fit_exch(d, main_formula())))

dir_terms <- DIRECTION_TERMS
dir_labels <- c(age10="年龄（每 10 岁）", female="女性（比男性）", apache10="严重度评分（每单位）",
                gcs3="GCS（每 3 分）", sbp10="收缩压（每 10 mmHg）", rr5="呼吸频率（每 5 次/分）",
                spo2_5="SpO2（每 5%）")
rows <- list()
for (tt in dir_terms) {
  a <- pe[pe$term == tt, ]; b <- pool_m[pool_m$term == tt, ]
  same_dir <- sign(a$estimate) == sign(b$estimate)
  overlap <- !(a$conf.high < b$conf.low || b$conf.high < a$conf.low)
  verdict <- if (!same_dir) "不一致" else if (overlap) "一致" else "方向一致但量级不同"
  rows[[length(rows)+1]] <- data.frame(
    variable = dir_labels[tt],
    eICU_OR = sprintf("%.2f (%.2f–%.2f)", exp(a$estimate), exp(a$conf.low), exp(a$conf.high)),
    MIMIC_OR = sprintf("%.2f (%.2f–%.2f)", exp(b$estimate), exp(b$conf.low), exp(b$conf.high)),
    direction_same = ifelse(same_dir, "是", "否"),
    ci_overlap = ifelse(overlap, "是", "否"),
    verdict = verdict, stringsAsFactors = FALSE)
}
dir_tab <- do.call(rbind, rows)
n_consistent <- sum(dir_tab$verdict == "一致")
overall <- if (n_consistent >= 5) "支持" else if (n_consistent == 4) "部分支持" else "不支持"
write.csv(dir_tab, "results/c3/tables/table4_direction.csv", row.names = FALSE, fileEncoding = "UTF-8")
cat("MIMIC v2 direction table:\n"); print(dir_tab)
cat("\n一致变量数 =", n_consistent, "/7 ->", overall, "\n")
