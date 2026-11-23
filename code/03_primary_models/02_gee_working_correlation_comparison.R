# c3_10_dxcorr_compare.R —— 分离「诊断映射更正」对主模型的影响
# 输入：results/c3/_data/imp_eicu.rds（v1 映射）、imp_eicu_v2.rds（v2 映射）
# 输出：results/c3/tables/table2_dx_v1_vs_v2.csv
# 说明：同一方法（exchangeable GEE，m=20，Rubin）、同一协变量集，只换 dx_group 映射，
#       因此差值可归因于诊断映射更正。运行：Rscript scripts/c3_10_dxcorr_compare.R

source((function(f) { for (p in c(file.path("code/00_common", f), file.path("../00_common", f), file.path("00_common", f), f)) if (file.exists(p)) return(p); f })("common_env.R"))
set.seed(SEED)

e1 <- readRDS("results/c3/_data/eicu_prepared.rds")       # v1
e2 <- readRDS("results/c3/_data/eicu_prepared_v2.rds")    # v2
ei <- subset(e2, in_analysis == 1)

build <- function(imp, ref) {
  ref <- subset(ref, in_analysis == 1)
  lapply(1:imp$m, function(i) {
    d <- complete(imp, i)
    d$to_invasive_48h <- as.integer(as.character(d$to_invasive_48h))
    d$dx_group <- factor(d$dx_group, levels = DX_LEVELS)
    d$age10 <- d$age_numeric/10; d$female <- as.integer(d$gender == "Female")
    d$apache10 <- d$apachescore/10; d$gcs3 <- d$gcs_total/3
    d$sbp10 <- d$sbp_min/10; d$rr5 <- d$rr_max/5; d$spo2_5 <- d$spo2_min/5
    d$m_apache <- as.integer(is.na(ref$apachescore)); d$m_gcs <- as.integer(is.na(ref$gcs_total))
    d$hosp_id <- factor(ref$hospitalid)
    d
  })
}
fit_exch <- function(d, formula) {
  dd <- d[complete.cases(d[, all.vars(formula)]), ]
  dd$.gee_id <- dd$hosp_id
  geeglm(formula, data = dd, family = binomial(), id = .gee_id,
         corstr = "exchangeable", std.err = "san.se")
}
pool_adj <- function(dl) {
  fs <- lapply(dl, function(d) fit_exch(d, main_formula()))
  list(pool = pool_rubin(fs), alpha = mean(sapply(fs, function(f) as.numeric(f$geese$alpha))))
}

imp1 <- readRDS("results/c3/_data/imp_eicu.rds"); cl1 <- build(imp1, e1)
imp2 <- readRDS("results/c3/_data/imp_eicu_v2.rds"); cl2 <- build(imp2, e2)
r1 <- pool_adj(cl1); r2 <- pool_adj(cl2)

mk <- function(p, tag) data.frame(term = p$term, OR = exp(p$estimate),
                                  low = exp(p$conf.low), high = exp(p$conf.high),
                                  p = p$p.value)
cmp <- merge(mk(r1$pool, "v1"), mk(r2$pool, "v2"), by = "term", all = TRUE, suffixes = c("_v1","_v2"))
cmp <- cmp[cmp$term != "(Intercept)", ]
cmp$OR_change_pct <- 100*(cmp$OR_v2/cmp$OR_v1 - 1)
cmp$sig_diff <- (cmp$p_v1 < 0.05) != (cmp$p_v2 < 0.05)
num <- sapply(cmp, is.numeric)
cmp[num] <- lapply(cmp[num], function(x) signif(x, 4))
write.csv(cmp, "results/c3/tables/table2_dx_v1_vs_v2.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

cat("mean exchangeable alpha: v1 =", round(r1$alpha, 4), " v2 =", round(r2$alpha, 4), "\n")
print(cmp)
