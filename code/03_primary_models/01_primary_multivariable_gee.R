# c3_07_gee_main.R —— 最终主模型：exchangeable GEE × 20 套插补 × Rubin 合并
# 输入：results/c3/_data/eicu_prepared_v2.rds（v2 诊断映射）
# 输出：results/c3/_data/imp_eicu_v2.rds
#       results/c3/tables/table2_main_model.csv（最终主模型）
#       results/c3/tables/table2_gee_vs_cluster.csv（两法对照）
# 运行：Rscript scripts/c3_07_gee_main.R
#
# 方法（在运行前定死）：
#   - 主模型：logistic GEE，working correlation = exchangeable，聚类变量 hospitalid，
#     稳健（sandwich）标准误；m = 20 套链式方程插补，Rubin 规则合并。
#   - 对照：同一批 20 套完成数据集上拟合「logistic + 医院聚类稳健方差（HC0）」，
#     即 independence working correlation 的 GEE；仅用于方法学对照，不作主结果。
#   - 协变量、缩放、缺失指示与方案 2.7.2 完全一致，未增删。

source((function(f) { for (p in c(file.path("code/00_common", f), file.path("../00_common", f), file.path("00_common", f), f)) if (file.exists(p)) return(p); f })("common_env.R"))
set.seed(SEED)

e <- readRDS("results/c3/_data/eicu_prepared_v2.rds")
ei <- subset(e, in_analysis == 1)

imp_vars <- c("to_invasive_48h","age_numeric","gender","dx_group","apachescore","gcs_total",
              "sbp_min","rr_max","spo2_min","ph_min","paco2_max","pf_min","lactate_max",
              "bicarbonate_first","hosp_nivfirst_n","hosp_ph_avail_pct")
dimp <- ei[, imp_vars]
dimp$to_invasive_48h <- factor(dimp$to_invasive_48h, levels = c(0, 1))
dimp$gender <- factor(dimp$gender)
meth <- make.method(dimp)
if ("gender" %in% names(meth)) meth["gender"] <- "logreg"

if (file.exists("results/c3/_data/imp_eicu_v2.rds")) {
  imp <- readRDS("results/c3/_data/imp_eicu_v2.rds")
  cat("reuse existing v2 imputation, m =", imp$m, "\n")
} else {
  imp <- mice(dimp, m = 20, maxit = 10, seed = SEED, method = meth, printFlag = FALSE)
  saveRDS(imp, "results/c3/_data/imp_eicu_v2.rds")
  cat("v2 imputation done, m =", imp$m, "\n")
}

complete_list <- lapply(1:imp$m, function(i) {
  d <- complete(imp, i)
  d$to_invasive_48h <- as.integer(as.character(d$to_invasive_48h))
  d$dx_group <- factor(d$dx_group, levels = DX_LEVELS)
  d$age10 <- d$age_numeric / 10
  d$female <- as.integer(d$gender == "Female")
  d$apache10 <- d$apachescore / 10
  d$gcs3 <- d$gcs_total / 3
  d$sbp10 <- d$sbp_min / 10
  d$rr5 <- d$rr_max / 5
  d$spo2_5 <- d$spo2_min / 5
  d$m_apache <- as.integer(is.na(ei$apachescore))
  d$m_gcs <- as.integer(is.na(ei$gcs_total))
  d$m_ph <- as.integer(is.na(ei$ph_min))
  d$m_paco2 <- as.integer(is.na(ei$paco2_max))
  d$m_pf <- as.integer(is.na(ei$pf_min))
  d$m_lactate <- as.integer(is.na(ei$lactate_max))
  d$hosp_id <- factor(ei$hospitalid)
  for (oc in c("to_invasive_24h","to_invasive_any","to_invasive_48h_anydef",
               "death_48h","death_hosp","composite_48h","composite_any",
               "followup_48h_complete","first_icu_stay","patient_unit_seq",
               "dnr_before_niv"))
    d[[oc]] <- ei[[oc]]
  d
})

# ---------------- 主模型：exchangeable GEE ----------------
f_main <- main_formula()
fit_exch <- function(d, formula, id = "hosp_id") {
  dd <- d[complete.cases(d[, all.vars(formula)]), ]
  dd$.gee_id <- dd[[id]]
  geeglm(formula, data = dd, family = binomial(), id = .gee_id,
         corstr = "exchangeable", std.err = "san.se")
}
fits_gee <- lapply(complete_list, function(d) fit_exch(d, f_main))
pool_gee <- pool_rubin(fits_gee)

# 未校正：单因素 exchangeable GEE（同一批插补）
uni_terms <- c("age10","female","dx_group","apache10","gcs3","sbp10","rr5","spo2_5")
uni <- list()
for (tt in uni_terms) {
  fs <- lapply(complete_list, function(d)
    fit_exch(d, as.formula(paste("to_invasive_48h ~", tt))))
  uni[[tt]] <- or_table(pool_rubin(fs), model_label = "unadjusted")
}
uni_df <- do.call(rbind, uni)
adj_df <- or_table(pool_gee, model_label = "adjusted (GEE, exchangeable)")

tab2 <- rbind(uni_df[, c("model","term","OR","OR_low","OR_high","p.value")],
              adj_df[, c("model","term","OR","OR_low","OR_high","p.value")])
tab2$OR <- round(tab2$OR, 2); tab2$OR_low <- round(tab2$OR_low, 2)
tab2$OR_high <- round(tab2$OR_high, 2); tab2$p.value <- signif(tab2$p.value, 3)
write.csv(tab2, "results/c3/tables/table2_main_model.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

# ---------------- 对照：聚类稳健 logistic（independence GEE） ----------------
fits_cl <- lapply(complete_list, function(d) fit_gee(d, f_main))
pool_cl <- pool_rubin(fits_cl)

a <- pool_gee; b <- pool_cl
cmp <- merge(
  data.frame(term = a$term, OR_gee = exp(a$estimate),
             OR_gee_low = exp(a$conf.low), OR_gee_high = exp(a$conf.high),
             p_gee = a$p.value),
  data.frame(term = b$term, OR_cl = exp(b$estimate),
             OR_cl_low = exp(b$conf.low), OR_cl_high = exp(b$conf.high),
             p_cl = b$p.value),
  by = "term", all = TRUE)
cmp$OR_ratio_gee_over_cl <- cmp$OR_gee / cmp$OR_cl
cmp$sig_diff <- (cmp$p_gee < 0.05) != (cmp$p_cl < 0.05)
num <- sapply(cmp, is.numeric)
cmp[num] <- lapply(cmp[num], function(x) signif(x, 3))
cmp <- cmp[cmp$term != "(Intercept)", ]
write.csv(cmp, "results/c3/tables/table2_gee_vs_cluster.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

# ---------------- 关键数字打印 ----------------
cat("\n=== 最终主模型（exchangeable GEE, m=20）===\n")
print(tab2[tab2$model == "adjusted (GEE, exchangeable)", ])
cat("\n=== 两法对照 ===\n")
print(cmp)

cat("\nn =", nrow(complete_list[[1]]), " events =", sum(complete_list[[1]]$to_invasive_48h),
    " hospitals =", length(unique(complete_list[[1]]$hosp_id)), "\n")
cat("GEE alpha (imputation 1):", fits_gee[[1]]$geese$alpha, "\n")
