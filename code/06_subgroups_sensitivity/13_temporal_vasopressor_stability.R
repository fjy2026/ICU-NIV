# c3temporal_03_vasopressor.R —— 敏感性 C：基线血管活性药加入主模型
# 输入：scripts/c3temporal_common.R 载入的主分析集与 20 套插补；
#       results/c3_ext/_data/mimic_ext.rds（MIMIC 外部验证准备集）
# 输出：results/c3_temporal/tableC_vasopressor.csv
#       results/c3_temporal/tableC_vasopressor_mimic.csv
#
# 设计（先定死，后看结果）：
#   - eICU 两套口径：
#       treatment 口径（全院治疗记录，覆盖 100%）：tx_vaso_pre6h_any / tx_vaso_at_niv_any；
#       infusionDrug 口径（仅 41% 有输液记录）：vaso_pre6h_any；
#         因来源不可得者的 0 不能解释为“未用药”，故
#         (i) 限 vaso_has_infusion_record = 1 的亚组分析；
#         (ii) 全样本 + 来源不可得指示变量的缺失层分析。
#   - 每个模型只在主模型基础上加一个二分类变量，其余协变量与插补不变；
#     eICU 用 exchangeable GEE、医院聚类、Rubin 合并（复用主分析 20 套插补，
#     因血管活性药变量在分析样本中完整可观测）。
#   - MIMIC 用 inputevents 口径 vaso_pre6h_any，沿用外部验证模型
#     （聚类稳健 logistic，患者聚类；m=20 插补）。
#   - 报告四个主关联（年龄、严重度、GCS、SBP）变化、血管活性药本身 OR、覆盖与缺失。

source((function(f) { for (p in c(file.path("code/00_common", f), file.path("../00_common", f), file.path("00_common", f), f)) if (file.exists(p)) return(p); f })("common_temporal.R"))
set.seed(SEED)

e   <- load_eicu_temporal()
imp <- readRDS("results/c3/_data/imp_eicu_v2.rds")
ei  <- subset(e, in_analysis == 1)
tmp_cols <- grep("^tmp_", names(ei), value = TRUE)
dl <- build_complete_list(e, imp, extra_from_e = tmp_cols)

lab <- c(age10="Age, per 10 years", female="Female vs male",
         apache10="Severity score, per unit", gcs3="GCS, per 3 points",
         sbp10="SBP, per 10 mmHg", rr5="RR, per 5 breaths/min",
         spo2_5="SpO2, per 5%")

extract_model <- function(pool, model, n, ev, extra_terms = character(0)) {
  terms <- c("age10","female","apache10","gcs3","sbp10","rr5","spo2_5", extra_terms)
  out <- do.call(rbind, lapply(terms, function(tt) {
    l <- if (tt %in% names(lab)) lab[tt] else tt
    r <- or_row(pool, tt, l, model, n, ev)
    if (is.null(r)) return(NULL)
    r
  }))
  out
}

rows <- list()

# ---- eICU treatment 口径 ----
for (spec in list(
  list(var = "tmp_tx_vaso_pre6h_any", label = "C1 eICU treatment, pre-NIV 6 h"),
  list(var = "tmp_tx_vaso_at_niv_any", label = "C2 eICU treatment, at NIV"))) {
  v <- spec$var
  bad <- is.na(ei[[v]])
  cat(spec$label, "missing:", sum(bad), "\n")
  idx <- which(!bad)
  f <- as.formula(paste("to_invasive_48h ~",
                        paste(c(MAIN_TERMS, v), collapse = " + ")))
  p <- fit_pool_subset(dl, f, idx)$pool
  rows[[length(rows)+1]] <- extract_model(p, spec$label, length(idx),
                                          sum(ei$to_invasive_48h[idx]),
                                          extra_terms = v)
}

# ---- eICU infusionDrug 口径 ----
# (i) 限有输液记录者（先拟合同子集的参照模型，再拟合加药模型，避免把样本选择
#     误当成血管活性药的效应）
idx_i <- which(ei$tmp_vaso_has_infusion_record == 1)
f_ref_i <- as.formula(paste("to_invasive_48h ~", paste(MAIN_TERMS, collapse = " + ")))
p_ref_i <- fit_pool_subset(dl, f_ref_i, idx_i)$pool
rows[[length(rows)+1]] <- extract_model(p_ref_i,
  "C3a eICU infusionDrug subgroup, reference (no vasopressor)",
  length(idx_i), sum(ei$to_invasive_48h[idx_i]))
f_i <- as.formula(paste("to_invasive_48h ~",
                        paste(c(MAIN_TERMS, "tmp_vaso_pre6h_any"), collapse = " + ")))
p_i <- fit_pool_subset(dl, f_i, idx_i)$pool
rows[[length(rows)+1]] <- extract_model(p_i, "C3b eICU infusionDrug, pre-NIV 6 h (recorded subgroup)",
                                        length(idx_i), sum(ei$to_invasive_48h[idx_i]),
                                        extra_terms = "tmp_vaso_pre6h_any")

# (ii) 全样本 + 来源不可得缺失层
dl_ml <- lapply(dl, function(d) {
  d$vaso_ml <- ifelse(d$tmp_vaso_has_infusion_record == 1, d$tmp_vaso_pre6h_any, 0)
  d$vaso_norec <- 1 - d$tmp_vaso_has_infusion_record
  d
})
f_ml <- as.formula(paste("to_invasive_48h ~",
                         paste(c(MAIN_TERMS, "vaso_ml", "vaso_norec"), collapse = " + ")))
p_ml <- fit_pool_subset(dl_ml, f_ml)$pool
rows[[length(rows)+1]] <- extract_model(p_ml, "C4 eICU infusionDrug, missing-record layer",
                                        nrow(ei), sum(ei$to_invasive_48h),
                                        extra_terms = c("vaso_ml","vaso_norec"))

tabC <- do.call(rbind, rows)
write.csv(tabC, file.path(OUT, "tableC_vasopressor.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# ---- MIMIC inputevents 口径 ----
m <- readRDS("results/c3_ext/_data/mimic_ext.rds")
mi <- subset(m, in_analysis == 1)
con <- dbConnect(duckdb(shared_home = TRUE))
mt <- dbGetQuery(con, "SELECT * FROM read_parquet('data/processed/c3_temporal/mimic_niv_temporal.parquet')")
dbDisconnect(con, shutdown = TRUE)
stopifnot(all(mi$stay_id %in% mt$stay_id))
mi$vaso_pre6h_any <- mt$vaso_pre6h_any[match(mi$stay_id, mt$stay_id)]
mi$year_index <- as.integer(format(as.POSIXct(mi$niv_start), "%Y")) -
  min(as.integer(format(as.POSIXct(mi$niv_start), "%Y")))
cu <- as.character(mi$first_careunit)
mi$careunit <- factor(ifelse(grepl("Medical/Surgical", cu), "MICU/SICU",
                      ifelse(grepl("Medical Intensive", cu), "MICU",
                      ifelse(grepl("Surgical|Trauma", cu), "SICU",
                      ifelse(grepl("Cardiac Vascular", cu), "CVICU",
                      ifelse(grepl("Coronary", cu), "CCU", "Other"))))))

imp_vars <- c("to_invasive_48h","age","gender","dx_group","sofa","gcs_min","sbp_min","rr_max",
              "spo2_min","ph_min","paco2_max","pf_min","lactate_max","bicarbonate_min",
              "year_index","careunit")
dimp <- mi[, imp_vars]
dimp$to_invasive_48h <- factor(dimp$to_invasive_48h, levels = c(0, 1))
dimp$gender <- factor(dimp$gender)
meth <- make.method(dimp)
if ("gender" %in% names(meth)) meth["gender"] <- "logreg"
set.seed(SEED + 5)
imp_m <- mice(dimp, m = 20, maxit = 10, method = meth, printFlag = FALSE)
dlm <- lapply(1:imp_m$m, function(i) {
  d <- complete(imp_m, i)
  d$to_invasive_48h <- as.integer(as.character(d$to_invasive_48h))
  d$dx_group <- factor(d$dx_group, levels = DX_LEVELS)
  d$age10 <- d$age / 10
  d$female <- as.integer(d$gender == "F")
  d$apache10 <- d$sofa / 2
  d$gcs3 <- d$gcs_min / 3
  d$sbp10 <- d$sbp_min / 10
  d$rr5 <- d$rr_max / 5
  d$spo2_5 <- d$spo2_min / 5
  d$pat_id <- factor(mi$subject_id)
  d$vaso <- mi$vaso_pre6h_any
  d
})
terms_m <- setdiff(MAIN_TERMS, c("m_apache","m_gcs"))
f_m_ref <- as.formula(paste("to_invasive_48h ~", paste(terms_m, collapse = " + ")))
f_m <- as.formula(paste("to_invasive_48h ~", paste(c(terms_m, "vaso"), collapse = " + ")))
fits_ref <- lapply(dlm, function(d) {
  dd <- d[complete.cases(d[, all.vars(f_m_ref)]), ]
  fit_gee(dd, f_m_ref, id = "pat_id")
})
p_ref <- pool_rubin(fits_ref)
tabC_m <- rbind(
  extract_model(p_ref, "C5a MIMIC reference (no vasopressor)",
                nrow(mi), sum(mi$to_invasive_48h)),
  {
    fits_m <- lapply(dlm, function(d) {
      dd <- d[complete.cases(d[, all.vars(f_m)]), ]
      fit_gee(dd, f_m, id = "pat_id")
    })
    p_m <- pool_rubin(fits_m)
    extract_model(p_m, "C5b MIMIC inputevents, pre-NIV 6 h",
                  nrow(mi), sum(mi$to_invasive_48h), extra_terms = "vaso")
  })
write.csv(tabC_m, file.path(OUT, "tableC_vasopressor_mimic.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# ---- 覆盖与事件率 ----
cov_rows <- list()
for (nm in c("tx_vaso_pre6h_any","tx_vaso_at_niv_any","vaso_pre6h_any","vaso_at_niv_any")) {
  x <- ei[[paste0("tmp_", nm)]]
  cov_rows[[length(cov_rows)+1]] <- data.frame(
    source = "eICU", variable = nm, n_total = nrow(ei), n_observed = sum(!is.na(x)),
    n_positive = sum(x == 1, na.rm = TRUE),
    event_rate_positive = mean(ei$to_invasive_48h[x == 1], na.rm = TRUE),
    event_rate_negative = mean(ei$to_invasive_48h[x == 0], na.rm = TRUE))
}
cov_rows[[length(cov_rows)+1]] <- data.frame(
  source = "eICU", variable = "vaso_has_infusion_record",
  n_total = nrow(ei), n_observed = nrow(ei),
  n_positive = sum(ei$tmp_vaso_has_infusion_record == 1),
  event_rate_positive = mean(ei$to_invasive_48h[ei$tmp_vaso_has_infusion_record == 1]),
  event_rate_negative = mean(ei$to_invasive_48h[ei$tmp_vaso_has_infusion_record == 0]))
cov_rows[[length(cov_rows)+1]] <- data.frame(
  source = "MIMIC", variable = "vaso_pre6h_any", n_total = nrow(mi),
  n_observed = sum(!is.na(mi$vaso_pre6h_any)), n_positive = sum(mi$vaso_pre6h_any == 1, na.rm = TRUE),
  event_rate_positive = mean(mi$to_invasive_48h[mi$vaso_pre6h_any == 1], na.rm = TRUE),
  event_rate_negative = mean(mi$to_invasive_48h[mi$vaso_pre6h_any == 0], na.rm = TRUE))
write.csv(do.call(rbind, cov_rows), file.path(OUT, "tableC_vasopressor_coverage.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

cat("\nC done.\n")
print(tabC)
print(tabC_m)
print(do.call(rbind, cov_rows))
