# c3temporal_01_early_covariates.R —— 敏感性 A：早期（NIV 前 / 0–2 h）协变量替换
# 输入：scripts/c3temporal_common.R 载入的主分析集与 20 套插补
# 输出：results/c3_temporal/tableA_early_covariate.csv
#       results/c3_temporal/tableA_early_covariate_restricted.csv
#
# 设计（先定死，后看结果）：
#   - 每个变体在“早期变量可观测”的同一批患者上做配对比较：
#       参照模型 = 主模型（0–24 h 协变量）在同一子集上重拟合；
#       替换模型 = 仅把对应的 GCS/SBP/RR/SpO2 换成早期值，其余协变量不变。
#     这样 OR 的变化归因于“测量时点”，而不是子集样本构成。
#   - 缺失处理：早期变量采用完整病例（不插补该变量）；其余协变量沿用主分析
#     的 20 套插补（Rubin 合并）。理由：早期变量缺失率高，且若用 24 h 值插补
#     早期值会把主模型的信息带入，反而掩盖时点差异。
#   - 估计：exchangeable GEE、医院聚类，与主模型一致；20 套插补并行拟合。
#   - 变体：A1 = NIV 前 6 h 最近一次；A2 = NIV 后 0–2 h 最差值；
#           A3 = 0–2 h 且转有创前（排除转有创后的测量）。
#   - 附加：限定 t0_niv_before_icu = 0 的联合替换。

source((function(f) { for (p in c(file.path("code/00_common", f), file.path("../00_common", f), file.path("00_common", f), f)) if (file.exists(p)) return(p); f })("common_temporal.R"))
set.seed(SEED)

e   <- load_eicu_temporal()
imp <- readRDS("results/c3/_data/imp_eicu_v2.rds")
ei  <- subset(e, in_analysis == 1)

dl <- build_complete_list(e, imp, extra_from_e = grep("^tmp_", names(ei), value = TRUE))

terms_ref <- MAIN_TERMS
key_terms <- c("age10","apache10","gcs3","sbp10","rr5","spo2_5")
lab <- c(age10="Age, per 10 years", apache10="Severity score, per unit",
         gcs3="GCS, per 3 points", sbp10="SBP, per 10 mmHg",
         rr5="RR, per 5 breaths/min", spo2_5="SpO2, per 5%")
vlab <- c(A1="A1: pre-NIV 6 h, last value", A2="A2: NIV 0-2 h, worst value",
          A3="A3: 0-2 h before intubation, worst value")

# ---- 精简数据（只保留建模列，便于并行导出） ----
base_cols <- c("to_invasive_48h","age10","female","dx_group","apache10",
               "gcs3","sbp10","rr5","spo2_5","m_apache","m_gcs","hosp_id")
early_def <- list(
  A1 = c(gcs3="tmp_pre_gcs_last6h", sbp10="tmp_pre_sbp_last6h",
         rr5="tmp_pre_rr_last6h", spo2_5="tmp_pre_spo2_last6h"),
  A2 = c(gcs3="tmp_post_gcs_worst_0_2h", sbp10="tmp_post_sbp_worst_0_2h",
         rr5="tmp_post_rr_worst_0_2h", spo2_5="tmp_post_spo2_worst_0_2h"),
  A3 = c(gcs3="tmp_post_gcs_worst_0_2h_preintub", sbp10="tmp_post_sbp_worst_0_2h_preintub",
         rr5="tmp_post_rr_worst_0_2h_preintub", spo2_5="tmp_post_spo2_worst_0_2h_preintub")
)
scale_div <- c(gcs3=3, sbp10=10, rr5=5, spo2_5=5)

dl_all <- lapply(dl, function(d) {
  s <- d[, base_cols, drop = FALSE]
  s$before <- d$tmp_t0_niv_before_icu
  for (v in names(early_def)) for (nm in names(early_def[[v]])) {
    s[[paste0(v, "_", nm)]] <- d[[early_def[[v]][nm]]] / scale_div[nm]
  }
  s
})

# ---- 模型任务清单 ----
jobs <- list()
mk_job <- function(v, analysis, replaced, idx, swap_terms, term_map) {
  list(v = v, analysis = analysis, replaced = replaced, idx = idx,
       swap_terms = swap_terms, term_map = term_map)
}
add_var_job <- function(v, analysis, replaced, idx, early_names) {
  swap <- terms_ref
  for (nm in names(early_names)) swap[swap == nm] <- paste0(v, "_", nm)
  mk_job(v, analysis, replaced, idx, swap, setNames(paste0(v, "_", names(early_names)), names(early_names)))
}
for (v in names(early_def)) {
  cols <- unname(early_def[[v]])
  idx_j <- which(complete.cases(ei[, cols, drop = FALSE]))
  # 联合替换
  jobs[[length(jobs)+1]] <- add_var_job(v, "Joint: all four replaced",
                                        "GCS+SBP+RR+SpO2", idx_j, early_def[[v]])
  # 逐一替换
  for (vv in names(early_def[[v]])) {
    idx_s <- which(!is.na(ei[[early_def[[v]][vv]]]))
    jobs[[length(jobs)+1]] <- add_var_job(v, "Single: one replaced", vv, idx_s,
                                          setNames(early_def[[v]][vv], vv))
  }
}
# 限定 NIV 不早于 ICU 入科（仅联合替换）
for (v in names(early_def)) {
  cols <- unname(early_def[[v]])
  idx_j <- which(complete.cases(ei[, cols, drop = FALSE]) & ei$tmp_t0_niv_before_icu == 0)
  j <- add_var_job(v, "Joint, restricted NIV>=ICU", "GCS+SBP+RR+SpO2", idx_j, early_def[[v]])
  j$restricted <- TRUE
  jobs[[length(jobs)+1]] <- j
}

# ---- 并行执行 ----
do_job <- function(job) {
  idx <- job$idx
  n <- length(idx); ev <- sum(dl_all[[1]]$to_invasive_48h[idx], na.rm = TRUE)
  f_ref <- as.formula(paste("to_invasive_48h ~", paste(terms_ref, collapse = " + ")))
  f_sw  <- as.formula(paste("to_invasive_48h ~", paste(job$swap_terms, collapse = " + ")))
  p_ref <- pool_rubin(lapply(dl_all, function(d) fit_exch(d[idx, , drop = FALSE], f_ref)))
  p_sw  <- pool_rubin(lapply(dl_all, function(d) fit_exch(d[idx, , drop = FALSE], f_sw)))
  out <- list()
  for (tt in key_terms) {
    sw_name <- if (tt %in% names(job$term_map)) job$term_map[[tt]] else tt
    rr_ref <- or_row(p_ref, tt, lab[tt], "Reference (main 0-24 h)", n, ev)
    rr_sw  <- or_row(p_sw,  sw_name, lab[tt], "Early value", n, ev)
    if (is.null(rr_ref) || is.null(rr_sw)) next
    out[[length(out)+1]] <- data.frame(
      variant = job$v, variant_label = vlab[job$v], analysis = job$analysis,
      replaced = job$replaced, term = tt, term_label = lab[tt],
      ref_OR = rr_ref$OR, ref_low = rr_ref$OR_low, ref_high = rr_ref$OR_high,
      ref_p = rr_ref$p.value,
      early_OR = rr_sw$OR, early_low = rr_sw$OR_low, early_high = rr_sw$OR_high,
      early_p = rr_sw$p.value,
      pct_change_OR = 100 * (rr_sw$OR / rr_ref$OR - 1),
      n = n, events = ev, stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

n_workers <- 8
cl <- parallel::makeCluster(n_workers)
parallel::clusterEvalQ(cl, { suppressPackageStartupMessages({library(geepack); library(sandwich)}) })
parallel::clusterExport(cl,
  c("dl_all","fit_exch","pool_rubin","align_fits","get_fixef","or_row",
    "do_job","key_terms","lab","terms_ref","vlab"),
  envir = environment())
res <- parallel::parLapply(cl, jobs, do_job)
parallel::stopCluster(cl)

tab <- do.call(rbind, res)
tab_main <- tab[is.na(tab$analysis) | tab$analysis != "Joint, restricted NIV>=ICU", ]
tab_restr <- tab[!is.na(tab$analysis) & tab$analysis == "Joint, restricted NIV>=ICU", ]
write.csv(tab_main, file.path(OUT, "tableA_early_covariate.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(tab_restr, file.path(OUT, "tableA_early_covariate_restricted.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

cat("A done.\nJoint sample sizes:\n")
print(unique(tab_main[tab_main$analysis == "Joint: all four replaced", c("variant","n","events")]))
cat("\nSingle sample sizes:\n")
print(unique(tab_main[tab_main$analysis == "Single: one replaced",
                      c("variant","replaced","n","events")]))
cat("\nRestricted sample sizes:\n")
print(unique(tab_restr[, c("variant","n","events")]))
