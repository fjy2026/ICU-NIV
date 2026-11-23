# c3_15_subgroups_hacor_v2.R —— v2 诊断映射下的预设亚组、HACOR 增量价值、连续变量 RCS
# 输入：results/c3/_data/eicu_prepared_v2.rds、results/c3/_data/imp_eicu_v2.rds
# 输出：results/c3/tables/table5_subgroups.csv
#       results/c3/tables/subgroup_interaction.csv
#       results/c3/tables/hacor_results.csv
#       results/c3/tables/rcs_tests.csv
#       results/c3/tables/rcs_curves.csv
# 运行：Rscript scripts/c3_15_subgroups_hacor_v2.R

source((function(f) { for (p in c(file.path("code/00_common", f), file.path("../00_common", f), file.path("00_common", f), f)) if (file.exists(p)) return(p); f })("common_env.R"))
set.seed(SEED)

e <- readRDS("results/c3/_data/eicu_prepared_v2.rds")
ei <- subset(e, in_analysis == 1)
imp <- readRDS("results/c3/_data/imp_eicu_v2.rds")

complete_list <- lapply(1:20, function(i) {
  d <- complete(imp, i)
  d$to_invasive_48h <- as.integer(as.character(d$to_invasive_48h))
  d$dx_group <- factor(d$dx_group, levels = DX_LEVELS)
  d$age10 <- d$age_numeric / 10; d$female <- as.integer(d$gender == "Female")
  d$apache10 <- d$apachescore / 10; d$gcs3 <- d$gcs_total / 3
  d$sbp10 <- d$sbp_min / 10; d$rr5 <- d$rr_max / 5; d$spo2_5 <- d$spo2_min / 5
  d$m_apache <- as.integer(is.na(ei$apachescore)); d$m_gcs <- as.integer(is.na(ei$gcs_total))
  d$hosp_id <- factor(ei$hospitalid)
  d
})

# 医院病例量三分位（在医院层面确定，再回填到入住）
hn <- unique(ei[, c("hospitalid","hosp_nivfirst_n")])
q <- quantile(hn$hosp_nivfirst_n, c(1/3, 2/3))
hn$vol_tertile <- cut(hn$hosp_nivfirst_n, breaks = c(-Inf, q, Inf),
                      labels = c("低病例量", "中病例量", "高病例量"), include.lowest = TRUE)
vol_map <- setNames(hn$vol_tertile, hn$hospitalid)

# 亚组定义（在分析集上生成）
ei$age_cat  <- factor(ifelse(ei$age_numeric >= 65, "≥65 岁", "<65 岁"), levels = c("<65 岁", "≥65 岁"))
ei$gcs_cat  <- factor(ifelse(ei$gcs_total >= 13, "GCS ≥13", "GCS <13"), levels = c("GCS ≥13", "GCS <13"))
ei$niv_site <- factor(ifelse(ei$niv_pre_icu == 1, "入 ICU 前开始", "入 ICU 后开始"),
                      levels = c("入 ICU 后开始", "入 ICU 前开始"))
ei$vol_cat  <- factor(vol_map[as.character(ei$hospitalid)], levels = c("低病例量","中病例量","高病例量"))

for (i in 1:20) {
  d <- complete_list[[i]]
  d$age_cat <- ei$age_cat; d$gcs_cat <- ei$gcs_cat
  d$niv_site <- ei$niv_site; d$vol_cat <- ei$vol_cat
  complete_list[[i]] <- d
}

# 稀疏规则：统计各层事件数
strata_defs <- list(
  list(name = "无创适应证", var = "dx_group", drop = "dx_group",
       levels = DX_LEVELS),
  list(name = "年龄", var = "age_cat", drop = "age10", levels = levels(ei$age_cat)),
  list(name = "意识状态", var = "gcs_cat", drop = "gcs3", levels = levels(ei$gcs_cat)),
  list(name = "无创开始时机", var = "niv_site", drop = NULL, levels = levels(ei$niv_site)),
  list(name = "医院病例量", var = "vol_cat", drop = NULL, levels = levels(ei$vol_cat))
)

sub_rows <- list(); inter_rows <- list()
for (sd in strata_defs) {
  for (lv in sd$levels) {
    idx <- !is.na(ei[[sd$var]]) & ei[[sd$var]] == lv
    ev <- sum(ei$to_invasive_48h[idx], na.rm = TRUE); n <- sum(idx)
    terms <- setdiff(MAIN_TERMS, sd$drop)
    if (!is.null(sd$drop) && sd$drop == "dx_group") terms <- setdiff(terms, "dx_group")
    status <- if (ev < 50) "事件数 <50，不报告 OR" else if (ev < 100) "事件数 <100，仅描述" else "校正模型"
    if (ev >= 100) {
      fs <- lapply(complete_list, function(d) {
        dd <- d[idx, ]; fit_gee(dd, as.formula(paste("to_invasive_48h ~", paste(terms, collapse=" + "))))
      })
      p <- pool_rubin(fs)
      keep <- p$term %in% DIRECTION_TERMS
      sub_rows[[length(sub_rows)+1]] <- data.frame(
        subgroup = sd$name, stratum = lv, n = n, events = ev, status = status,
        term = p$term[keep],
        OR = sprintf("%.2f (%.2f–%.2f)", exp(p$estimate[keep]), exp(p$conf.low[keep]), exp(p$conf.high[keep])),
        p = signif(p$p.value[keep], 3), stringsAsFactors = FALSE)
    } else {
      sub_rows[[length(sub_rows)+1]] <- data.frame(
        subgroup = sd$name, stratum = lv, n = n, events = ev, status = status,
        term = DIRECTION_TERMS, OR = NA_character_, p = NA_real_, stringsAsFactors = FALSE)
    }
  }
  # 交互检验：亚组指示变量 × 全部协变量（联合 Wald）
  terms <- setdiff(MAIN_TERMS, sd$drop)
  if (!is.null(sd$drop) && sd$drop == "dx_group") terms <- setdiff(terms, "dx_group")
  form_int <- as.formula(paste("to_invasive_48h ~ (", paste(terms, collapse=" + "), ") * ", sd$var))
  fits_int <- lapply(complete_list, function(d) fit_gee(d, form_int))
  pf <- pool_rubin_full(fits_int)
  inter_terms <- grep(paste0(":", sd$var), pf$names, value = TRUE)
  jw <- joint_wald(pf, inter_terms)
  inter_rows[[length(inter_rows)+1]] <- data.frame(
    subgroup = sd$name, n_terms = jw$df, Wald = round(jw$W, 2),
    p_raw = jw$p, stringsAsFactors = FALSE)
}
sub_tab <- do.call(rbind, sub_rows)
write.csv(sub_tab, "results/c3/tables/table5_subgroups.csv", row.names = FALSE, fileEncoding = "UTF-8")

inter_tab <- do.call(rbind, inter_rows)
inter_tab$p_holm <- p.adjust(inter_tab$p_raw, method = "holm")
write.csv(inter_tab, "results/c3/tables/subgroup_interaction.csv", row.names = FALSE, fileEncoding = "UTF-8")

# ---------------- HACOR ----------------
hc <- ei[!is.na(ei$hacor_score), ]
hc_cc <- hc[complete.cases(hc[, c("to_invasive_48h","age_numeric","gender","dx_group","apachescore",
                                  "gcs_total","sbp_min","rr_max","spo2_min")]), ]
hc_cc$age10 <- hc_cc$age_numeric/10; hc_cc$female <- as.integer(hc_cc$gender=="Female")
hc_cc$apache10 <- hc_cc$apachescore/10; hc_cc$gcs3 <- hc_cc$gcs_total/3
hc_cc$sbp10 <- hc_cc$sbp_min/10; hc_cc$rr5 <- hc_cc$rr_max/5; hc_cc$spo2_5 <- hc_cc$spo2_min/5
hc_cc$m_apache <- as.integer(is.na(hc_cc$apachescore)); hc_cc$m_gcs <- as.integer(is.na(hc_cc$gcs_total))
hc_cc$dx_group <- factor(hc_cc$dx_group, levels = DX_LEVELS)
hc_cc$hacor_score <- as.numeric(hc_cc$hacor_score)

m1 <- glm(to_invasive_48h ~ hacor_score, data = hc_cc, family = binomial)
m2 <- glm(to_invasive_48h ~ hacor_score + age10 + female + dx_group + apache10 + gcs3 +
            sbp10 + rr5 + spo2_5 + m_apache + m_gcs, data = hc_cc, family = binomial)
m3 <- glm(to_invasive_48h ~ age10 + female + dx_group + apache10 + gcs3 +
            sbp10 + rr5 + spo2_5 + m_apache + m_gcs, data = hc_cc, family = binomial)

auc_of <- function(fit, nd) {
  pr <- predict(fit, newdata = nd, type = "response")
  as.numeric(pROC::auc(pROC::roc(nd$to_invasive_48h, pr, quiet = TRUE)))
}
a1 <- auc_of(m1, hc_cc); a2 <- auc_of(m2, hc_cc); a3 <- auc_of(m3, hc_cc)

# bootstrap 乐观度校正与 AUC 差值区间
B <- 200
opt <- matrix(NA, B, 3); diffs <- numeric(B)
for (b in 1:B) {
  idx <- sample(nrow(hc_cc), replace = TRUE)
  db <- hc_cc[idx, ]
  f1 <- try(glm(to_invasive_48h ~ hacor_score, data = db, family = binomial), silent = TRUE)
  f2 <- try(glm(to_invasive_48h ~ hacor_score + age10 + female + dx_group + apache10 + gcs3 +
                  sbp10 + rr5 + spo2_5 + m_apache + m_gcs, data = db, family = binomial), silent = TRUE)
  f3 <- try(glm(to_invasive_48h ~ age10 + female + dx_group + apache10 + gcs3 +
                  sbp10 + rr5 + spo2_5 + m_apache + m_gcs, data = db, family = binomial), silent = TRUE)
  if (any(sapply(list(f1,f2,f3), inherits, "try-error"))) next
  ab <- c(auc_of(f1, db), auc_of(f2, db), auc_of(f3, db))
  ao <- c(auc_of(f1, hc_cc), auc_of(f2, hc_cc), auc_of(f3, hc_cc))
  opt[b, ] <- ab - ao
  diffs[b] <- ab[2] - ab[3]
}
corrected <- c(a1, a2, a3) - colMeans(opt, na.rm = TRUE)
hacor_res <- data.frame(
  model = c("仅 HACOR", "HACOR + 一级协变量", "仅一级协变量"),
  n = nrow(hc_cc), events = sum(hc_cc$to_invasive_48h),
  AUC_apparent = round(c(a1, a2, a3), 3),
  AUC_corrected = round(corrected, 3)
)
hacor_extra <- data.frame(
  model = c("AUC 差（HACOR+协变量 减 仅协变量）", "差值 95% CI 下限", "差值 95% CI 上限",
            "HACOR 可计算（任一窗口）", "其中 2h 窗口", "其中 24h 窗口", "P/F 可计算"),
  value = c(sprintf("%.3f", a2 - a3),
            sprintf("%.3f", quantile(diffs, .025, na.rm = TRUE)),
            sprintf("%.3f", quantile(diffs, .975, na.rm = TRUE)),
            sum(!is.na(ei$hacor_score)), sum(ei$hacor_window == "2h", na.rm = TRUE),
            sum(ei$hacor_window == "24h", na.rm = TRUE), sum(!is.na(ei$pf_min))))
write.csv(hacor_res, "results/c3/tables/hacor_models.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(hacor_extra, "results/c3/tables/hacor_supplement.csv", row.names = FALSE, fileEncoding = "UTF-8")

# ---------------- RCS（完整病例，主模型协变量） ----------------
cc <- ei[complete.cases(ei[, c("to_invasive_48h","age_numeric","gender","dx_group","apachescore",
                               "gcs_total","sbp_min","rr_max","spo2_min")]), ]
cc$age10 <- cc$age_numeric/10; cc$female <- as.integer(cc$gender=="Female")
cc$apache10 <- cc$apachescore/10; cc$gcs3 <- cc$gcs_total/3
cc$sbp10 <- cc$sbp_min/10; cc$rr5 <- cc$rr_max/5; cc$spo2_5 <- cc$spo2_min/5
cc$m_apache <- 0; cc$m_gcs <- 0
cc$dx_group <- factor(cc$dx_group, levels = DX_LEVELS)
cc$y <- cc$to_invasive_48h

rcs_specs <- list(
  list(v = "age_numeric", label = "年龄（岁）", other = c("female","dx_group","apache10","gcs3","sbp10","rr5","spo2_5")),
  list(v = "apachescore", label = "APACHE IVa 评分", other = c("age10","female","dx_group","gcs3","sbp10","rr5","spo2_5")),
  list(v = "sbp_min", label = "收缩压（mmHg）", other = c("age10","female","dx_group","apache10","gcs3","rr5","spo2_5")),
  list(v = "gcs_total", label = "GCS 总分", other = c("age10","female","dx_group","apache10","sbp10","rr5","spo2_5"))
)
rcs_tests <- list(); curves <- list()
for (sp in rcs_specs) {
  f_rcs <- as.formula(paste("y ~ rcs(", sp$v, ", 4) +", paste(sp$other, collapse = " + ")))
  fit <- lrm(f_rcs, data = cc, x = TRUE, y = TRUE)
  an <- anova(fit)
  rn <- trimws(rownames(an))
  p_overall <- an[grep(paste0("^", sp$v), rn)[1], "P"]
  p_nonlin <- an[grep("Nonlinear", rn)[1], "P"]
  rcs_tests[[length(rcs_tests)+1]] <- data.frame(variable = sp$label,
    p_overall = signif(p_overall, 3), p_nonlinear = signif(p_nonlin, 3), n = nrow(cc))
  grid <- seq(quantile(cc[[sp$v]], .01), quantile(cc[[sp$v]], .99), length.out = 100)
  ref <- data.frame(
    age_numeric = median(cc$age_numeric), apachescore = median(cc$apachescore),
    sbp_min = median(cc$sbp_min), gcs_total = median(cc$gcs_total),
    age10 = median(cc$age10), female = 0,
    dx_group = factor("COPD/哮喘", levels = DX_LEVELS),
    apache10 = median(cc$apache10), gcs3 = median(cc$gcs3),
    sbp10 = median(cc$sbp10), rr5 = median(cc$rr5), spo2_5 = median(cc$spo2_5))
  nd <- ref[rep(1, length(grid)), ]
  nd[[sp$v]] <- grid
  pr <- predict(fit, newdata = nd, type = "lp", se.fit = TRUE)
  curves[[length(curves)+1]] <- data.frame(
    variable = sp$label, x = grid,
    risk = plogis(pr$linear.predictors),
    low = plogis(pr$linear.predictors - 1.96*pr$se.fit),
    high = plogis(pr$linear.predictors + 1.96*pr$se.fit),
    ref = median(cc[[sp$v]]),
    p_nonlinear = signif(p_nonlin, 3))
}
write.csv(do.call(rbind, rcs_tests), "results/c3/tables/rcs_tests.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(do.call(rbind, curves), "results/c3/tables/rcs_curves.csv", row.names = FALSE, fileEncoding = "UTF-8")

cat("subgroups / HACOR / RCS done\n")
print(inter_tab)
print(hacor_res)
print(do.call(rbind, rcs_tests))
