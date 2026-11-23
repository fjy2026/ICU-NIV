# c3_08_dnr_sensitivity.R —— 基线 DNR/DNI 敏感性模型（S10）＋ 诊断映射更正后的 S5b
# 输入：results/c3/_data/eicu_prepared_v2.rds、imp_eicu_v2.rds
# 输出：results/c3/tables/sensitivity_dnr.csv
#       results/c3/tables/sensitivity_summary.csv（追加 S10 与 S5b-v2，保留 S1–S9）
# 运行：Rscript scripts/c3_08_dnr_sensitivity.R
#
# 方法（运行前定死）：
#   - S10：在最终主模型（exchangeable GEE，m=20，Rubin）基础上只增加
#     dnr_before_niv 一个 0/1 协变量；DNR 缺失 14 例在本敏感性分析中排除，
#     因此同时拟合「同一样本、不含 DNR」的对照模型，两者之差即为 DNR 的影响。
#   - S5b-v2：用 v2 诊断映射排除「术后/其他」组，exchangeable GEE 重跑，
#     用于替换 v1 关键词映射下的 S5b 结果（原结果保留在 sensitivity_summary_v1.csv）。

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
  d$dnr_before_niv <- ei$dnr_before_niv
  d
})

fit_exch <- function(d, formula, id = "hosp_id") {
  dd <- d[complete.cases(d[, all.vars(formula)]), ]
  if (is.factor(dd$dx_group)) dd$dx_group <- droplevels(dd$dx_group)
  dd$.gee_id <- dd[[id]]
  geeglm(formula, data = dd, family = binomial(), id = .gee_id,
         corstr = "exchangeable", std.err = "san.se")
}
pool_terms <- function(dl, terms, subset_idx = NULL) {
  if (!is.null(subset_idx)) dl <- lapply(dl, function(d) d[subset_idx, , drop = FALSE])
  fs <- lapply(dl, function(d) fit_exch(d, as.formula(paste("to_invasive_48h ~",
                                                           paste(terms, collapse = " + ")))))
  p <- pool_rubin(fs)
  list(pool = p, fits = fs)
}

# 最终主模型（全样本，用于判定 S5b/S10 是否改变主结论）
main_pool <- pool_terms(complete_list, MAIN_TERMS)$pool
main_or <- setNames(exp(main_pool$estimate), main_pool$term)
main_sig <- setNames(main_pool$p.value < 0.05, main_pool$term)
key4 <- c("age10","apache10","gcs3","sbp10")
key_labels <- c(age10="年龄/10 岁", apache10="严重度/单位", gcs3="GCS/3 分", sbp10="收缩压/10 mmHg")

changed_vs <- function(p, ref_or = main_or, ref_sig = main_sig) {
  ch <- character(0)
  for (tt in key4) {
    or <- exp(p$estimate[p$term == tt]); lo <- exp(p$conf.low[p$term == tt]); hi <- exp(p$conf.high[p$term == tt])
    if (length(or) == 0) next
    if (isTRUE(sign(log(or)) != sign(log(ref_or[tt]))) ||
        (ref_sig[tt] && lo <= 1 && hi >= 1)) ch <- c(ch, key_labels[tt])
  }
  if (length(ch) == 0) "无" else paste(ch, collapse = "、")
}

# ---------------- S10：DNR/DNI ----------------
ok_dnr <- which(!is.na(ei$dnr_before_niv))
mod_no <- pool_terms(complete_list, MAIN_TERMS, ok_dnr)
mod_dnr <- pool_terms(complete_list, c(MAIN_TERMS, "dnr_before_niv"), ok_dnr)

p_no <- mod_no$pool; p_dnr <- mod_dnr$pool
rows <- list()
for (tt in key4) {
  or0 <- exp(p_no$estimate[p_no$term == tt]); lo0 <- exp(p_no$conf.low[p_no$term == tt]); hi0 <- exp(p_no$conf.high[p_no$term == tt])
  or1 <- exp(p_dnr$estimate[p_dnr$term == tt]); lo1 <- exp(p_dnr$conf.low[p_dnr$term == tt]); hi1 <- exp(p_dnr$conf.high[p_dnr$term == tt])
  rows[[length(rows)+1]] <- data.frame(
    variable = key_labels[tt],
    OR_without_DNR = sprintf("%.2f (%.2f–%.2f)", or0, lo0, hi0),
    p_without_DNR = signif(p_no$p.value[p_no$term == tt], 3),
    OR_with_DNR = sprintf("%.2f (%.2f–%.2f)", or1, lo1, hi1),
    p_with_DNR = signif(p_dnr$p.value[p_dnr$term == tt], 3),
    OR_change_pct = sprintf("%+.1f", 100*(or1/or0 - 1)))
}
dnr_eff <- data.frame(
  variable = "DNR（基线，1=有；eICU 仅记录不复苏医嘱）",
  OR_without_DNR = "—", p_without_DNR = NA,
  OR_with_DNR = sprintf("%.2f (%.2f–%.2f)",
                        exp(p_dnr$estimate[p_dnr$term == "dnr_before_niv"]),
                        exp(p_dnr$conf.low[p_dnr$term == "dnr_before_niv"]),
                        exp(p_dnr$conf.high[p_dnr$term == "dnr_before_niv"])),
  p_with_DNR = signif(p_dnr$p.value[p_dnr$term == "dnr_before_niv"], 3),
  OR_change_pct = "—")
out <- rbind(do.call(rbind, rows), dnr_eff)
write.csv(out, "results/c3/tables/sensitivity_dnr.csv", row.names = FALSE, fileEncoding = "UTF-8")
print(out)

# ---------------- S5b-v2：按 v2 映射排除「术后/其他」 ----------------
ok_s5b <- which(ei$dx_group != "术后/其他")
mod_s5b <- pool_terms(complete_list, MAIN_TERMS, ok_s5b)

# ---------------- 更新 sensitivity_summary.csv ----------------
old <- read.csv("results/c3/tables/sensitivity_summary_v1.csv", stringsAsFactors = FALSE,
                fileEncoding = "UTF-8")
old$version <- "v1（聚类稳健 logistic；v1 诊断映射）"
new_rows <- data.frame(
  id = c("S5b-v2", "S10"),
  description = c("排除术后/其他诊断组（v2 映射，GEE 重跑）",
                  "主模型加入基线 DNR/DNI 协变量"),
  n = c(length(ok_s5b), length(ok_dnr)),
  events = c(sum(ei$to_invasive_48h[ok_s5b], na.rm = TRUE),
             sum(ei$to_invasive_48h[ok_dnr], na.rm = TRUE)),
  changed = c(changed_vs(mod_s5b$pool), changed_vs(p_dnr, main_or, main_sig)),
  note = c("替换 v1 映射下的 S5b",
           sprintf("DNR 组 OR %s；排除 DNR 缺失 14 例",
                   sprintf("%.2f", exp(p_dnr$estimate[p_dnr$term == "dnr_before_niv"])))),
  version = "v2（exchangeable GEE；v2 诊断映射）",
  stringsAsFactors = FALSE)
out_sum <- rbind(old, new_rows)
write.csv(out_sum, "results/c3/tables/sensitivity_summary.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
cat("\n=== S5b-v2 / S10 是否改变主结论 ===\n")
print(new_rows[, c("id","n","events","changed")])
