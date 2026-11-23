# c3_21_secondary_bloodgas_v2.R —— 含血气与乳酸的次要模型（v2、exchangeable GEE）
# 输入：results/c3/_data/eicu_prepared_v2.rds、imp_eicu_v2.rds
# 输出：results/c3/tables/table2b_secondary_model.csv
# 运行：Rscript scripts/c3_21_secondary_bloodgas_v2.R
# 说明：协变量为主模型全部一级变量＋pH、PaCO2、P/F、乳酸、碳酸氢根及其缺失指示，
#       与 c3_02_main.R 的次要模型一致；估计方法改为 exchangeable GEE、m=20、Rubin。
#       血气缺失由多重插补处理，分析样本仍为 13,986、事件 1,345；另报告血气完整病例数。

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
  d$m_ph <- as.integer(is.na(ei$ph_min)); d$m_paco2 <- as.integer(is.na(ei$paco2_max))
  d$m_pf <- as.integer(is.na(ei$pf_min)); d$m_lactate <- as.integer(is.na(ei$lactate_max))
  d$hosp_id <- factor(ei$hospitalid)
  d$ph10 <- d$ph_min / 0.1 - 73.4
  d$paco2_10 <- d$paco2_max / 10
  d$pf50 <- d$pf_min / 50
  d$lactate1 <- d$lactate_max
  d$bicarb5 <- d$bicarbonate_first / 5
  d
})
fit_exch <- function(d, formula) {
  dd <- d[complete.cases(d[, all.vars(formula)]), ]
  if (is.factor(dd$dx_group)) dd$dx_group <- droplevels(dd$dx_group)
  dd$.gee_id <- dd$hosp_id
  geeglm(formula, data = dd, family = binomial(), id = .gee_id,
         corstr = "exchangeable", std.err = "san.se")
}
sec_terms <- c(MAIN_TERMS, "ph10", "paco2_10", "pf50", "lactate1", "bicarb5",
               "m_ph", "m_paco2", "m_pf", "m_lactate")
f_sec <- as.formula(paste("to_invasive_48h ~", paste(sec_terms, collapse = " + ")))
fits <- lapply(complete_list, function(d) fit_exch(d, f_sec))
pool <- pool_rubin(fits)
tab <- or_table(pool, model_label = "secondary (GEE, exchangeable, +blood gas)")
tab$OR <- round(tab$OR, 2); tab$OR_low <- round(tab$OR_low, 2)
tab$OR_high <- round(tab$OR_high, 2); tab$p.value <- signif(tab$p.value, 3)
tab$n <- nrow(complete_list[[1]])
tab$events <- sum(complete_list[[1]]$to_invasive_48h)
write.csv(tab, "results/c3/tables/table2b_secondary_model.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

cc_bg <- sum(complete.cases(ei[, c("ph_min","paco2_max","pf_min","lactate_max")]))
cat("secondary blood-gas model done; n =", nrow(complete_list[[1]]),
    " events =", sum(complete_list[[1]]$to_invasive_48h),
    " blood-gas complete cases =", cc_bg, "\n")
print(tab[tab$term %in% c("ph10","paco2_10","pf50","lactate1","bicarb5"), ])
