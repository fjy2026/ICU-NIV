# c3_ext2_04_t0.R —— 时间零点敏感性：ICU 入科 vs 无创开始
# 输入：扩取数据的时间字段（只读）+ 主分析 20 套插补
# 输出：results/c3_ext2/tables/t0_sensitivity.csv
# 说明：补充/敏感性分析。eICU 用 inv_start_strict 相对 ICU 入科的分钟偏移；
#      MIMIC 用 procedureevents 225792 首次开始时间与 icustays.intime 计算。
#      T0=ICU 入科的“全队列”含 NIV 在 48 h 后才开始者（结构性非事件）；
#      另给“NIV 于入 ICU 48 h 内开始”的限制队列。
source((function(f) { for (p in c(file.path("code/00_common", f), file.path("../00_common", f), file.path("00_common", f), f)) if (file.exists(p)) return(p); f })("common_transportability.R"))
set.seed(SEED)
dir.create(file.path(OUT, "tables"), recursive = TRUE, showWarnings = FALSE)

# ---------------- eICU ----------------
e <- load_eicu_ext2(); ei <- subset(e, in_analysis == 1)
cl <- completed_eicu_ext2(ei)
f_main <- main_formula()

fit_pool <- function(dlist, formula, id = "hosp_id", subset_expr = NULL) {
  fs <- lapply(dlist, function(d) {
    if (!is.null(subset_expr)) d <- d[subset_expr(d), ]
    fit_exch(d, formula, id = id)
  })
  rubin_pool(fs)
}
term_lab <- c("(Intercept)" = "(Intercept)",
  age10 = "Age (per 10 years)", female = "Female (vs male)",
  apache10 = "Severity score (per unit)", gcs3 = "GCS (per 3 points)",
  sbp10 = "Systolic BP (per 10 mmHg)", rr5 = "Respiratory rate (per 5/min)",
  spo2_5 = "SpO2 (per 5%)", m_apache = "APACHE missing indicator",
  m_gcs = "GCS missing indicator",
  "dx_group心源性肺水肿" = "Diagnostic group: cardiogenic pulmonary oedema",
  "dx_group肺炎/其他呼吸" = "Diagnostic group: pneumonia/other respiratory",
  "dx_group脓毒症" = "Diagnostic group: sepsis",
  "dx_group神经系统" = "Diagnostic group: neurological",
  "dx_group术后/其他" = "Diagnostic group: postoperative/other")
label_term <- function(t) ifelse(t %in% names(term_lab), term_lab[t], t)

or_tab <- function(pool, label, n, events) {
  data.frame(database = label[1], definition = label[2], n = n, events = events,
             term = pool$term, term_label = label_term(pool$term), OR = exp(pool$estimate),
             OR_low = exp(pool$conf.low), OR_high = exp(pool$conf.high),
             p = pool$p.value, stringsAsFactors = FALSE)
}

pool_main <- fit_pool(cl, f_main)
pool_t0   <- fit_pool(cl, main_formula("y_t0_icu"))
pool_t0r  <- fit_pool(cl, main_formula("y_t0_icu"), subset_expr = function(d) !is.na(d$y_t0_icu_restricted))

# 核对主模型与已发表主结果一致
pub <- read.csv("results/c3/tables/table2_main_model.csv", fileEncoding = "UTF-8")
pub <- pub[pub$model == "adjusted (GEE, exchangeable)", ]
chk <- merge(pub[, c("term","OR","OR_low","OR_high")],
             data.frame(term = pool_main$term, OR2 = exp(pool_main$estimate)), by = "term")
cat("Main model reproduction max |OR diff| =", max(abs(chk$OR - chk$OR2), na.rm = TRUE), "\n")

n_full <- nrow(ei); ev_t0 <- sum(ei$y_t0_icu)
restr <- !is.na(ei$y_t0_icu_restricted)
eicu_rows <- rbind(
  or_tab(pool_main, c("eICU","Main: T0 = NIV start"), n_full, sum(ei$to_invasive_48h)),
  or_tab(pool_t0,   c("eICU","T0 = ICU admission, full cohort"), n_full, ev_t0),
  or_tab(pool_t0r,  c("eICU","T0 = ICU admission, NIV within 48 h"), sum(restr), sum(ei$y_t0_icu_restricted[restr]))
)

cat("\neICU event counts: main", sum(ei$to_invasive_48h), " T0-ICU full", ev_t0,
    " T0-ICU restricted", sum(ei$y_t0_icu_restricted[restr]), " (n=", sum(restr), ")\n")
cat("NIV started >48 h after ICU admission:", sum(ei$niv_start_offset > 2880, na.rm = TRUE), "\n")
cat("Cross-tab main vs T0-ICU:\n"); print(table(main = ei$to_invasive_48h, t0icu = ei$y_t0_icu))

# ---------------- MIMIC ----------------
mi <- load_mimic_ext2(); cm <- completed_mimic_ext2(mi)
f_mimic <- main_formula(terms = setdiff(MAIN_TERMS, c("m_apache","m_gcs")))
fit_pool_m <- function(dlist, formula, subset_expr = NULL) {
  fs <- lapply(dlist, function(d) {
    if (!is.null(subset_expr)) d <- d[subset_expr(d), ]
    fit_cluster_glm(d, formula, id = "pat_id")
  })
  rubin_pool(fs)
}
pool_m_main <- fit_pool_m(cm, f_mimic)
pool_m_t0   <- fit_pool_m(cm, main_formula("y_t0_icu", setdiff(MAIN_TERMS, c("m_apache","m_gcs"))))
pool_m_t0r  <- fit_pool_m(cm, main_formula("y_t0_icu", setdiff(MAIN_TERMS, c("m_apache","m_gcs"))),
                          subset_expr = function(d) !is.na(d$y_t0_icu_restricted))
restr_m <- !is.na(mi$y_t0_icu_restricted)
mimic_rows <- rbind(
  or_tab(pool_m_main, c("MIMIC","Main: T0 = NIV start"), nrow(mi), sum(mi$to_invasive_48h)),
  or_tab(pool_m_t0,   c("MIMIC","T0 = ICU admission, full cohort"), nrow(mi), sum(mi$y_t0_icu)),
  or_tab(pool_m_t0r,  c("MIMIC","T0 = ICU admission, NIV within 48 h"), sum(restr_m), sum(mi$y_t0_icu_restricted[restr_m]))
)
cat("\nMIMIC event counts: main", sum(mi$to_invasive_48h), " T0-ICU full", sum(mi$y_t0_icu),
    " T0-ICU restricted", sum(mi$y_t0_icu_restricted[restr_m]), " (n=", sum(restr_m), ")\n")
cat("NIV started >48 h after ICU admission:", sum(mi$niv_after_intime_h > 48, na.rm = TRUE), "\n")
cat("Cross-tab main vs T0-ICU (MIMIC):\n"); print(table(main = mi$to_invasive_48h, t0icu = mi$y_t0_icu))

out <- rbind(eicu_rows, mimic_rows)
out[, c("OR","OR_low","OR_high")] <- lapply(out[, c("OR","OR_low","OR_high")], function(x) round(x, 3))
out$p <- signif(out$p, 3)
write.csv(out, file.path(OUT, "tables", "t0_sensitivity.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# 四个主关联的并列表（便于阅读）
four <- c("age10","apache10","gcs3","sbp10")
cat("\nFour main associations (eICU):\n")
print(out[out$database == "eICU" & out$term %in% four, c("definition","term","OR","OR_low","OR_high","p")])
cat("\nFour main associations (MIMIC):\n")
print(out[out$database == "MIMIC" & out$term %in% four, c("definition","term","OR","OR_low","OR_high","p")])

writeLines(capture.output(sessionInfo()), file.path(OUT, "_data", "session_info_t0.txt"))
cat("\nt0 sensitivity done\n")
