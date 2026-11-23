# c3_ext2_05_alignment.R —— 两库变量对齐表
# 输入：eICU/MIMIC 扩取数据与码本（只读）
# 输出：results/c3_ext2/tables/variable_alignment.csv
# 说明：逐项列来源、定义、单位、主分析集缺失率、是否可直接比较；缺失率由数据计算。
source((function(f) { for (p in c(file.path("code/00_common", f), file.path("../00_common", f), file.path("00_common", f), f)) if (file.exists(p)) return(p); f })("common_transportability.R"))
dir.create(file.path(OUT, "tables"), recursive = TRUE, showWarnings = FALSE)

e <- load_eicu_ext2(); ei <- subset(e, in_analysis == 1)
mi <- load_mimic_ext2()

miss <- function(x) sprintf("%d (%.1f%%)", sum(is.na(x)), 100*mean(is.na(x)))
n_e <- nrow(ei); n_m <- nrow(mi)

rows <- list()
add <- function(item, src_e, def_e, unit_e, x_e, src_m, def_m, unit_m, x_m, comparable, note) {
  rows[[length(rows)+1]] <<- data.frame(
    item = item,
    eICU_source = src_e, eICU_definition = def_e, eICU_unit = unit_e,
    eICU_missing_main = miss(x_e),
    MIMIC_source = src_m, MIMIC_definition = def_m, MIMIC_unit = unit_m,
    MIMIC_missing_main = miss(x_m),
    directly_comparable = comparable, note = note, stringsAsFactors = FALSE)
}

add("Age", "patient.age",
    "Admission age; values '> 89' recoded to 90", "years", ei$age_numeric,
    "patients.anchor_age + year(admittime) - anchor_year",
    "Derived admission age", "years", mi$age,
    "Yes (eICU top-coded at 90)", "eICU de-identification top-codes very old patients; MIMIC age is derived and not top-coded in the same way.")

add("Severity score", "apachePatientResult.apachescore (apacheversion = IVa, score > 0)",
    "APACHE IVa admission score", "points", ei$apachescore,
    "mimiciv_derived.first_day_sofa.sofa",
    "First-day SOFA score", "points", mi$sofa,
    "No (different instruments)",
    "The model uses APACHE IVa / 10 in eICU and SOFA / 2 in MIMIC; these are not equivalent scales. SOFA itself contains a GCS component.")

add("Glasgow Coma Scale", "apacheApsVar.eyes + motor + verbal",
    "APACHE admission GCS (sum of three components)", "points (3-15)", ei$gcs_total,
    "mimiciv_derived.first_day_gcs.gcs_min",
    "Minimum GCS during first ICU day", "points (3-15)", mi$gcs_min,
    "Partial", "eICU uses the APACHE admission value; MIMIC uses the first-day minimum. Both scored 3-15 but timing and derivation differ.")

add("Systolic blood pressure", "vitalAperiodic.noninvasivesystolic",
    "Minimum non-invasive systolic BP, first 24 h (20-300 mmHg)", "mmHg", ei$sbp_min,
    "mimiciv_derived.first_day_vitalsign.sbp_min",
    "Minimum systolic BP, first ICU day (all sources)", "mmHg", mi$sbp_min,
    "Partial", "eICU uses non-invasive cuff readings only; MIMIC includes invasive and non-invasive readings.")

add("Respiratory rate", "vitalPeriodic.respiration",
    "Maximum respiratory rate, first 24 h (0-80/min)", "breaths/min", ei$rr_max,
    "mimiciv_derived.first_day_vitalsign.resp_rate_max",
    "Maximum respiratory rate, first ICU day", "breaths/min", mi$rr_max,
    "Partial", "Both are first-day maxima, but derived from different charting systems and device sources.")

add("SpO2", "vitalPeriodic.sao2",
    "Minimum SpO2, first 24 h (0-100%)", "%", ei$spo2_min,
    "mimiciv_derived.first_day_vitalsign.spo2_min",
    "Minimum SpO2, first ICU day", "%", mi$spo2_min,
    "Partial", "Both are first-day minima; eICU uses periodic vital-sign records, MIMIC uses charted vital signs.")

add("NIV definition", "treatment.treatmentoffset",
    "Free-text treatment record containing non-invasive ventilation; earliest offset used as NIV start", "timestamp/offset", ei$niv_start_offset,
    "procedureevents.starttime (itemid = 225794)",
    "Procedure record of non-invasive ventilation; earliest starttime used", "timestamp", mi$niv_start,
    "Partial", "eICU NIV is identified from free-text treatment strings; MIMIC NIV from a procedure item. This is the main source of the different cohort sizes and event rates (eICU 9.6% vs MIMIC 23.1%).")

add("Invasive ventilation outcome", "treatment.treatmentoffset",
    "First strict invasive treatment after NIV start, 0 < delta <= 48 h; NIV text excluded", "binary", ei$to_invasive_48h,
    "procedureevents.starttime (itemid = 225792)",
    "First invasive procedure after NIV start, 0 < delta <= 48 h", "binary", mi$to_invasive_48h,
    "Partial", "Different source systems (free-text treatment vs procedure table). Both use a strict invasive definition. Event rates are not directly comparable.")

add("Time zero", "NIV start offset from ICU admission (treatmentoffset)",
    "First NIV treatment record", "minutes from ICU admission", ei$niv_start_offset,
    "procedureevents.starttime (itemid = 225794)",
    "First NIV procedure start", "absolute timestamp", mi$niv_start,
    "Yes (conceptually)", "Both anchor at NIV start. Sensitivity analysis re-anchors at ICU admission using eICU treatment offsets and MIMIC icustays.intime.")

add("In-hospital death", "patient.hospitaldischargestatus",
    "Expired vs Alive at hospital discharge; discharge offset for timing", "binary / timestamp", ei$death_hosp,
    "admissions.hospital_expire_flag / admissions.deathtime",
    "In-hospital death flag and death time", "binary / timestamp", mi$death_hosp,
    "Yes", "Both define in-hospital death; 48-h death is computed from the same time zero (NIV start).")

out <- do.call(rbind, rows)
out$eICU_n_main <- n_e; out$MIMIC_n_main <- n_m
write.csv(out, file.path(OUT, "tables", "variable_alignment.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("Variable alignment rows:", nrow(out), " (eICU n =", n_e, "; MIMIC n =", n_m, ")\n")
print(out[, c("item","eICU_missing_main","MIMIC_missing_main","directly_comparable")])

writeLines(capture.output(sessionInfo()), file.path(OUT, "_data", "session_info_alignment.txt"))
cat("\nalignment done\n")
