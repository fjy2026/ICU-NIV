# c3_14_secondary_v2.R —— 最终模型（GEE、v2 映射）下的次要结局
# 输入：results/c3/_data/{eicu_prepared_v2,imp_eicu_v2}.rds
# 输出：results/c3/tables/table3_secondary_outcomes.csv（v2）
# 运行：Rscript scripts/c3_14_secondary_v2.R
# 说明：结局集合、协变量与 c3_02_main.R 完全一致，只把估计方法改为 exchangeable GEE、
#       诊断映射改为 v2。仍为探索性分析，未校正多重比较。

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
  for (oc in c("to_invasive_24h","to_invasive_any","death_48h","death_hosp","composite_48h"))
    d[[oc]] <- ei[[oc]]
  d
})
fit_exch <- function(d, formula) {
  dd <- d[complete.cases(d[, all.vars(formula)]), ]
  dd$.gee_id <- dd$hosp_id
  geeglm(formula, data = dd, family = binomial(), id = .gee_id,
         corstr = "exchangeable", std.err = "san.se")
}

outcomes <- c("to_invasive_any","death_48h","death_hosp","composite_48h","to_invasive_24h")
sec <- list()
for (oo in outcomes) {
  fs <- lapply(complete_list, function(d) {
    dd <- d; dd$y <- dd[[oo]]; dd <- dd[!is.na(dd$y), ]
    fit_exch(dd, as.formula(paste("y ~", paste(MAIN_TERMS, collapse = " + "))))
  })
  p <- pool_rubin(fs)
  p$outcome <- oo
  p$events <- sum(complete_list[[1]][[oo]], na.rm = TRUE)
  p$n <- sum(!is.na(complete_list[[1]][[oo]]))
  sec[[oo]] <- p
}
df <- do.call(rbind, sec)
df$OR <- round(exp(df$estimate), 2)
df$OR_low <- round(exp(df$conf.low), 2)
df$OR_high <- round(exp(df$conf.high), 2)
df$p.value <- signif(df$p.value, 3)
write.csv(df[, c("outcome","events","n","term","OR","OR_low","OR_high","p.value")],
          "results/c3/tables/table3_secondary_outcomes.csv", row.names = FALSE, fileEncoding = "UTF-8")
print(df[df$term %in% c("age10","apache10","gcs3","sbp10","rr5","spo2_5"),
         c("outcome","events","term","OR","OR_low","OR_high","p.value")])
