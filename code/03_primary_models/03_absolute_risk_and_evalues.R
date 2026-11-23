# c3_13_absrisk_evalue_v2.R —— 最终模型下的绝对风险与 E 值
# 输入：results/c3/_data/{eicu_prepared_v2,imp_eicu_v2}.rds
# 输出：results/c3/tables/absolute_risk.csv、evalues.csv（v2）
# 运行：Rscript scripts/c3_13_absrisk_evalue_v2.R

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
  d$hosp_id <- factor(ei$hospitalid); d
})
fit_exch <- function(d, formula) {
  dd <- d[complete.cases(d[, all.vars(formula)]), ]
  dd$.gee_id <- dd$hosp_id
  geeglm(formula, data = dd, family = binomial(), id = .gee_id,
         corstr = "exchangeable", std.err = "san.se")
}
fits <- lapply(complete_list, function(d) fit_exch(d, main_formula()))
pool <- pool_rubin(fits)
beta <- setNames(pool$estimate, pool$term)

# 绝对风险
d1 <- complete_list[[1]]
ref <- data.frame(
  age10 = median(d1$age10), female = 0,
  dx_group = factor("COPD/哮喘", levels = DX_LEVELS),
  apache10 = median(d1$apache10), gcs3 = 15/3, sbp10 = median(d1$sbp10),
  rr5 = median(d1$rr5), spo2_5 = median(d1$spo2_5),
  m_apache = 0, m_gcs = 0)
risk_sets <- list(
  list(v = "gcs3", label = "GCS 15/12/9", values = c(5, 4, 3), shown = c(15, 12, 9)),
  list(v = "sbp10", label = "收缩压 100/80/60 mmHg", values = c(10, 8, 6), shown = c(100, 80, 60)),
  list(v = "apache10", label = "APACHE IVa 40/60/80", values = c(4, 6, 8), shown = c(40, 60, 80)),
  list(v = "age10", label = "年龄 40/60/80 岁", values = c(4, 6, 8), shown = c(40, 60, 80)))
risk_rows <- list()
for (rs in risk_sets) {
  r <- risk_at(beta, ref, rs$v, rs$values)
  risk_rows[[length(risk_rows)+1]] <- data.frame(
    variable = rs$label, value = rs$shown,
    predicted_risk_pct = sprintf("%.1f", 100*r$risk))
}
write.csv(do.call(rbind, risk_rows), "results/c3/tables/absolute_risk.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

# E 值（四个主关联）
key_labels <- c(age10="年龄/10 岁", apache10="严重度/单位", sbp10="收缩压/10 mmHg", gcs3="GCS/3 分")
ev_rows <- list()
for (tt in names(key_labels)) {
  or <- exp(pool$estimate[pool$term == tt])
  lo <- exp(pool$conf.low[pool$term == tt]); hi <- exp(pool$conf.high[pool$term == tt])
  ev <- EValue::evalues.OR(est = or, lo = lo, hi = hi, rare = FALSE)
  evp <- ev["E-values", "point"]
  evc <- if (or < 1) ev["E-values", "upper"] else ev["E-values", "lower"]
  ev_rows[[length(ev_rows)+1]] <- data.frame(
    variable = key_labels[tt], OR = sprintf("%.2f (%.2f–%.2f)", or, lo, hi),
    E_value_point = sprintf("%.2f", evp), E_value_CI = sprintf("%.2f", evc))
}
write.csv(do.call(rbind, ev_rows), "results/c3/tables/evalues.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
print(read.csv("results/c3/tables/absolute_risk.csv", fileEncoding = "UTF-8"))
print(read.csv("results/c3/tables/evalues.csv", fileEncoding = "UTF-8"))
