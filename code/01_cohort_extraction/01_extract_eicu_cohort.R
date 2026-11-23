#!/usr/bin/env Rscript
# C3 主队列：eICU-CRD 2.0 无创通气（NIV）ICU 入住正式取数
# 一行 = 一次 ICU 入住（patientunitstayid）。保留全部有无创记录的入住，
# 用 in_analysis 标记主分析集（NIV-first、非同刻、年龄 >= 18 且年龄不缺失），
# 使 S2/S2b/S3 等敏感性分析可在同一张表上完成。
# 只读源数据；不插补、不建模、不做统计检验。
suppressPackageStartupMessages({ library(DBI); library(duckdb) })
source((function(f) { for (p in c(file.path("code/00_common", f), file.path("../00_common", f), file.path("00_common", f), f)) if (file.exists(p)) return(p); f })("trace_logger.R"))

P   <- "F:/CoredataCenter/eICU/parquet"
OUT <- "data/processed/c3"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
con <- dbConnect(duckdb::duckdb())
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

# ---------------------------------------------------------------- 1. 队列与治疗先后
cohort_sql <- sprintf("
WITH niv AS (
  SELECT patientunitstayid, MIN(treatmentoffset) AS t0
  FROM read_parquet('%s/treatment.parquet')
  WHERE lower(treatmentstring) LIKE '%%non-invasive ventilation%%'
  GROUP BY 1
),
tx AS (
  SELECT t.patientunitstayid,
    MIN(CASE WHEN lower(t.treatmentstring) LIKE '%%mechanical ventilation%%'
              AND lower(t.treatmentstring) NOT LIKE '%%non-invasive ventilation%%'
             THEN t.treatmentoffset END) AS inv_start_strict,
    MIN(CASE WHEN lower(t.treatmentstring) LIKE '%%mechanical ventilation%%'
             THEN t.treatmentoffset END) AS inv_start_any,
    MIN(CASE WHEN lower(t.treatmentstring) LIKE '%%mechanical ventilation%%'
              AND lower(t.treatmentstring) NOT LIKE '%%non-invasive ventilation%%'
              AND t.treatmentoffset > n.t0 THEN t.treatmentoffset END) AS t1_strict_any,
    MIN(CASE WHEN lower(t.treatmentstring) LIKE '%%mechanical ventilation%%'
              AND lower(t.treatmentstring) NOT LIKE '%%non-invasive ventilation%%'
              AND t.treatmentoffset > n.t0 AND t.treatmentoffset <= n.t0 + 2880
             THEN t.treatmentoffset END) AS t1_strict_48h,
    MIN(CASE WHEN lower(t.treatmentstring) LIKE '%%mechanical ventilation%%'
              AND t.treatmentoffset > n.t0 THEN t.treatmentoffset END) AS t1_any_any,
    MIN(CASE WHEN lower(t.treatmentstring) LIKE '%%mechanical ventilation%%'
              AND t.treatmentoffset > n.t0 AND t.treatmentoffset <= n.t0 + 2880
             THEN t.treatmentoffset END) AS t1_any_48h,
    MIN(CASE WHEN lower(t.treatmentstring) LIKE '%%mechanical ventilation%%'
              AND lower(t.treatmentstring) NOT LIKE '%%non-invasive ventilation%%'
              AND t.treatmentoffset = n.t0 THEN 1 ELSE 0 END) AS strict_same
  FROM read_parquet('%s/treatment.parquet') t
  JOIN niv n USING (patientunitstayid)
  GROUP BY 1
)
SELECT p.patientunitstayid, p.patienthealthsystemstayid, p.uniquepid, p.hospitalid,
       p.unitvisitnumber, p.hospitaladmitoffset, p.hospitaldischargeoffset,
       p.unitdischargeoffset, p.unitdischargestatus, p.hospitaldischargestatus,
       p.age, p.gender, p.ethnicity, p.apacheadmissiondx,
       n.t0 AS niv_start_offset,
       tx.inv_start_strict, tx.inv_start_any,
       tx.t1_strict_any, tx.t1_strict_48h, tx.t1_any_any, tx.t1_any_48h,
       tx.strict_same
FROM read_parquet('%s/patient.parquet') p
JOIN niv n USING (patientunitstayid)
LEFT JOIN tx USING (patientunitstayid)
", P, P, P)
cohort <- dbGetQuery(con, cohort_sql)
cohort <- mc_checkpoint(cohort, "纳入：有无创通气（NIV）治疗的 ICU 入住", id = "patientunitstayid")

# 年龄（'> 89' 记 90；无法解析记缺失）
cohort$age_numeric <- suppressWarnings(as.numeric(cohort$age))
cohort$age_numeric[grepl("^>\\s*89", cohort$age)] <- 90
cohort$age_missing <- as.integer(is.na(cohort$age_numeric))
cohort$age_lt18    <- as.integer(!is.na(cohort$age_numeric) & cohort$age_numeric < 18)

# 先后与同刻
cohort$inv_before_niv <- as.integer(!is.na(cohort$inv_start_strict) & cohort$inv_start_strict < cohort$niv_start_offset)
cohort$inv_same_time  <- as.integer(!is.na(cohort$inv_start_strict) & cohort$inv_start_strict == cohort$niv_start_offset)
cohort$niv_first      <- 1L - cohort$inv_before_niv
cohort$niv_pre_icu    <- as.integer(cohort$niv_start_offset < 0)
cohort$niv_record     <- 1L
cohort$first_icu_stay <- as.integer(cohort$unitvisitnumber == 1)
cohort$patient_unit_seq <- ave(seq_len(nrow(cohort)), cohort$uniquepid, FUN = function(i) rank(cohort$patientunitstayid[i], ties.method = "first"))

# 三步队列标记：整表保留，只打标记（S2/S2b/S3 需要全样本）
cohort <- mc_checkpoint(cohort, "排除有创记录早于无创（序贯撤机）", id = "patientunitstayid",
                        keep = cohort$inv_before_niv == 0)
cohort <- mc_checkpoint(cohort, "排除有创与无创开始同刻记录", id = "patientunitstayid",
                        keep = cohort$inv_before_niv == 0 & cohort$inv_same_time == 0)
cohort$in_analysis <- as.integer(cohort$inv_before_niv == 0 & cohort$inv_same_time == 0 &
                                   !is.na(cohort$age_numeric) & cohort$age_numeric >= 18)
cohort <- mc_checkpoint(cohort, "纳入：NIV-first、非同刻、年龄 >= 18 岁且年龄不缺失", id = "patientunitstayid",
                        keep = cohort$in_analysis == 1)

# ---------------------------------------------------------------- 2. APACHE / GCS / 生理评分
apache <- dbGetQuery(con, sprintf("
SELECT patientunitstayid,
       MAX(CASE WHEN apacheversion = 'IVa' AND apachescore > 0 THEN apachescore END) AS apachescore,
       MAX(CASE WHEN apacheversion = 'IVa' AND predictedhospitalmortality >= 0
                THEN predictedhospitalmortality END) AS predictedhospitalmortality
FROM read_parquet('%s/apachePatientResult.parquet')
WHERE apacheversion = 'IVa' GROUP BY 1", P))
aps <- dbGetQuery(con, sprintf("
SELECT patientunitstayid, eyes, motor, verbal, intubated
FROM read_parquet('%s/apacheApsVar.parquet')", P))
aps$gcs_total <- ifelse(aps$eyes >= 0 & aps$motor >= 0 & aps$verbal >= 0,
                        aps$eyes + aps$motor + aps$verbal, NA_real_)
aps$intubated_apache <- ifelse(aps$intubated < 0, NA_real_, aps$intubated)
aps <- aps[, c("patientunitstayid", "gcs_total", "intubated_apache")]

# ---------------------------------------------------------------- 3. 0–24h 生命体征（ICU 入科为零点）
vit <- dbGetQuery(con, sprintf("
SELECT patientunitstayid,
  MIN(CASE WHEN heartrate BETWEEN 0 AND 300 THEN heartrate END) AS hr_min,
  MAX(CASE WHEN heartrate BETWEEN 0 AND 300 THEN heartrate END) AS hr_max,
  MAX(CASE WHEN respiration BETWEEN 0 AND 80 THEN respiration END) AS rr_max,
  MIN(CASE WHEN sao2 BETWEEN 0 AND 100 THEN sao2 END) AS spo2_min
FROM read_parquet('%s/vitalPeriodic.parquet')
WHERE observationoffset BETWEEN 0 AND 1440
GROUP BY 1", P))
nibp <- dbGetQuery(con, sprintf("
SELECT patientunitstayid,
  MIN(CASE WHEN noninvasivesystolic BETWEEN 20 AND 300 THEN noninvasivesystolic END) AS sbp_min
FROM read_parquet('%s/vitalAperiodic.parquet')
WHERE observationoffset BETWEEN 0 AND 1440
GROUP BY 1", P))

# ---------------------------------------------------------------- 4. 0–24h 血气与化验（ICU 入科为零点）
labs <- dbGetQuery(con, sprintf("
SELECT patientunitstayid,
  MIN(CASE WHEN labname='pH' AND labresult BETWEEN 6.5 AND 8.0 THEN labresult END) AS ph_min,
  MAX(CASE WHEN labname='paCO2' AND labresult BETWEEN 5 AND 200 THEN labresult END) AS paco2_max,
  MAX(CASE WHEN labname='lactate' AND labresult BETWEEN 0 AND 30 THEN labresult END) AS lactate_max,
  arg_min(CASE WHEN labname='bicarbonate' AND labresult BETWEEN 1 AND 60 THEN labresult END, labresultoffset) AS bicarbonate_first,
  arg_min(CASE WHEN labname='FiO2' AND labresult BETWEEN 21 AND 100 THEN labresult END, labresultoffset) AS fio2_first
FROM read_parquet('%s/lab.parquet')
WHERE labresultoffset BETWEEN 0 AND 1440
GROUP BY 1", P))

# 同一偏移配对计算 PaO2/FiO2：0–24h（ICU 零点）与 0–2h / 0–24h（NIV 零点）
labpiv <- dbGetQuery(con, sprintf("
WITH niv AS (
  SELECT patientunitstayid FROM read_parquet('%s/treatment.parquet')
  WHERE lower(treatmentstring) LIKE '%%non-invasive ventilation%%' GROUP BY 1
),
l AS (
  SELECT lab.patientunitstayid, lab.labresultoffset,
         MAX(CASE WHEN labname='paO2' AND labresult BETWEEN 10 AND 700 THEN labresult END) AS po2,
         MAX(CASE WHEN labname='FiO2' AND labresult BETWEEN 21 AND 100 THEN labresult END) AS fio2,
         MIN(CASE WHEN labname='pH' AND labresult BETWEEN 6.5 AND 8.0 THEN labresult END) AS ph
  FROM read_parquet('%s/lab.parquet') lab JOIN niv USING (patientunitstayid)
  GROUP BY 1,2
)
SELECT l.patientunitstayid, l.labresultoffset, l.po2, l.fio2, l.ph,
       CASE WHEN l.po2 IS NOT NULL AND l.fio2 IS NOT NULL THEN l.po2 / (l.fio2 / 100.0) END AS pf
FROM l
WHERE l.po2 IS NOT NULL OR l.fio2 IS NOT NULL OR l.ph IS NOT NULL
", P, P))

pf_24_icu <- aggregate(pf ~ patientunitstayid, data = subset(labpiv, !is.na(pf) &
                        labresultoffset >= 0 & labresultoffset <= 1440), FUN = min, na.rm = TRUE)
names(pf_24_icu)[2] <- "pf_min"
labpiv$t0 <- cohort$niv_start_offset[match(labpiv$patientunitstayid, cohort$patientunitstayid)]
hacor_lab_2h  <- aggregate(pf ~ patientunitstayid, data = subset(labpiv, !is.na(pf) &
                        labresultoffset >= t0 & labresultoffset <= t0 + 120), FUN = min, na.rm = TRUE)
hacor_lab_24h <- aggregate(pf ~ patientunitstayid, data = subset(labpiv, !is.na(pf) &
                        labresultoffset >= t0 & labresultoffset <= t0 + 1440), FUN = min, na.rm = TRUE)
ph_2h <- aggregate(ph ~ patientunitstayid, data = subset(labpiv, !is.na(ph) &
                        labresultoffset >= t0 & labresultoffset <= t0 + 120), FUN = min, na.rm = TRUE)
ph_24h <- aggregate(ph ~ patientunitstayid, data = subset(labpiv, !is.na(ph) &
                        labresultoffset >= t0 & labresultoffset <= t0 + 1440), FUN = min, na.rm = TRUE)
names(hacor_lab_2h)[2] <- "pf_hacor_2h"; names(hacor_lab_24h)[2] <- "pf_hacor_24h"
names(ph_2h)[2] <- "ph_hacor_2h"; names(ph_24h)[2] <- "ph_hacor_24h"

# ---------------------------------------------------------------- 5. HACOR 的 HR / RR / GCS（NIV 零点）
hacor_vit_exact <- dbGetQuery(con, sprintf("
WITH niv AS (SELECT patientunitstayid, MIN(treatmentoffset) AS t0
  FROM read_parquet('%s/treatment.parquet')
  WHERE lower(treatmentstring) LIKE '%%non-invasive ventilation%%' GROUP BY 1)
SELECT n.patientunitstayid,
  MAX(CASE WHEN v.heartrate BETWEEN 0 AND 300 AND v.observationoffset BETWEEN n.t0 AND n.t0+120 THEN v.heartrate END) AS hr_2h,
  MAX(CASE WHEN v.heartrate BETWEEN 0 AND 300 AND v.observationoffset BETWEEN n.t0 AND n.t0+1440 THEN v.heartrate END) AS hr_24h,
  MAX(CASE WHEN v.respiration BETWEEN 0 AND 80 AND v.observationoffset BETWEEN n.t0 AND n.t0+120 THEN v.respiration END) AS rr_2h,
  MAX(CASE WHEN v.respiration BETWEEN 0 AND 80 AND v.observationoffset BETWEEN n.t0 AND n.t0+1440 THEN v.respiration END) AS rr_24h
FROM niv n JOIN read_parquet('%s/vitalPeriodic.parquet') v USING (patientunitstayid)
GROUP BY 1", P, P))
gcs <- dbGetQuery(con, sprintf("
WITH niv AS (SELECT patientunitstayid, MIN(treatmentoffset) AS t0
  FROM read_parquet('%s/treatment.parquet')
  WHERE lower(treatmentstring) LIKE '%%non-invasive ventilation%%' GROUP BY 1)
SELECT n.patientunitstayid,
  MIN(CASE WHEN c.nursingchartcelltypevallabel='Glasgow coma score'
            AND c.nursingchartcelltypevalname='GCS Total'
            AND TRY_CAST(c.nursingchartvalue AS DOUBLE) BETWEEN 3 AND 15
            AND c.nursingchartoffset BETWEEN n.t0 AND n.t0+120
           THEN TRY_CAST(c.nursingchartvalue AS DOUBLE) END) AS gcs_hacor_2h,
  MIN(CASE WHEN c.nursingchartcelltypevallabel='Glasgow coma score'
            AND c.nursingchartcelltypevalname='GCS Total'
            AND TRY_CAST(c.nursingchartvalue AS DOUBLE) BETWEEN 3 AND 15
            AND c.nursingchartoffset BETWEEN n.t0 AND n.t0+1440
           THEN TRY_CAST(c.nursingchartvalue AS DOUBLE) END) AS gcs_hacor_24h
FROM niv n JOIN read_parquet('%s/nurseCharting.parquet') c USING (patientunitstayid)
GROUP BY 1", P, P))

# ---------------------------------------------------------------- 6. 治疗限制与 NIV 参数
dnr <- dbGetQuery(con, sprintf("
WITH niv AS (SELECT patientunitstayid, MIN(treatmentoffset) AS t0
  FROM read_parquet('%s/treatment.parquet')
  WHERE lower(treatmentstring) LIKE '%%non-invasive ventilation%%' GROUP BY 1)
SELECT n.patientunitstayid,
  MAX(CASE WHEN lower(c.cplgroup) = 'care limitation'
            AND lower(c.cplitemvalue) LIKE 'do not resuscitate%%'
            AND c.cplitemoffset <= n.t0 THEN 1 ELSE 0 END) AS dnr_before_niv
FROM niv n JOIN read_parquet('%s/carePlanGeneral.parquet') c USING (patientunitstayid)
GROUP BY 1", P, P))
nivset <- dbGetQuery(con, sprintf("
SELECT patientunitstayid,
  AVG(CASE WHEN respchartvaluelabel='NIV Setting EPAP'
            AND TRY_CAST(respchartvalue AS DOUBLE) BETWEEN 1 AND 40
           THEN TRY_CAST(respchartvalue AS DOUBLE) END) AS niv_epap_mean,
  AVG(CASE WHEN respchartvaluelabel='NIV Setting Set_RR'
            AND TRY_CAST(respchartvalue AS DOUBLE) BETWEEN 1 AND 60
           THEN TRY_CAST(respchartvalue AS DOUBLE) END) AS niv_set_rr_mean
FROM read_parquet('%s/respiratoryCharting.parquet')
WHERE respchartoffset BETWEEN 0 AND 1440
  AND respchartvaluelabel LIKE 'NIV Setting%%'
GROUP BY 1", P))

# ---------------------------------------------------------------- 7. 合并
apache <- mc_checkpoint(apache, "并入：APACHE IVa 评分", id = "patientunitstayid", side = TRUE)
aps    <- mc_checkpoint(aps, "并入：APACHE 生理变量与 GCS", id = "patientunitstayid", side = TRUE)
vit    <- mc_checkpoint(vit, "并入：0–24h 生命体征", id = "patientunitstayid", side = TRUE)
nibp   <- mc_checkpoint(nibp, "并入：0–24h 无创袖带血压", id = "patientunitstayid", side = TRUE)
labs   <- mc_checkpoint(labs, "并入：0–24h 血气与乳酸", id = "patientunitstayid", side = TRUE)
dnr    <- mc_checkpoint(dnr, "并入：无创前 DNR 记录", id = "patientunitstayid", side = TRUE)
nivset <- mc_checkpoint(nivset, "并入：NIV 设定参数", id = "patientunitstayid", side = TRUE)
hacor_vit_exact <- mc_checkpoint(hacor_vit_exact, "并入：HACOR 心率与呼吸频率", id = "patientunitstayid", side = TRUE)
gcs    <- mc_checkpoint(gcs, "并入：HACOR GCS", id = "patientunitstayid", side = TRUE)

out <- cohort
for (blk in list(apache, aps, vit, nibp, labs, dnr, nivset)) {
  out <- merge(out, blk, by = "patientunitstayid", all.x = TRUE, sort = FALSE)
}
for (blk in list(pf_24_icu, hacor_vit_exact, gcs, ph_2h, ph_24h, hacor_lab_2h, hacor_lab_24h)) {
  out <- merge(out, blk, by = "patientunitstayid", all.x = TRUE, sort = FALSE)
}
stopifnot(!anyDuplicated(out$patientunitstayid))

# ---------------------------------------------------------------- 8. 结局与派生
out$death_hosp <- ifelse(is.na(out$hospitaldischargestatus), NA_integer_,
                         ifelse(out$hospitaldischargestatus == "Expired", 1L, 0L))
out$death_time <- ifelse(!is.na(out$death_hosp) & out$death_hosp == 1,
                         out$hospitaldischargeoffset, NA_real_)
out$death_offset <- out$death_time

# 转有创：t1 严格晚于 t0；死亡早于 t1 时记 0
inv_event <- function(t1, td) {
  as.integer(!is.na(t1) & (is.na(td) | t1 < td))
}
out$to_invasive_48h      <- inv_event(out$t1_strict_48h, out$death_time)
out$to_invasive_24h      <- as.integer(!is.na(out$t1_strict_any) &
                              out$t1_strict_any <= out$niv_start_offset + 1440 &
                              (is.na(out$death_time) | out$t1_strict_any < out$death_time))
out$to_invasive_any      <- inv_event(out$t1_strict_any, out$death_time)
out$to_invasive_48h_anydef <- inv_event(out$t1_any_48h, out$death_time)

out$death_48h <- as.integer(!is.na(out$death_time) & !is.na(out$niv_start_offset) &
                              out$death_time - out$niv_start_offset > 0 &
                              out$death_time - out$niv_start_offset <= 2880)
out$composite_48h            <- as.integer(out$to_invasive_48h == 1 | out$death_hosp == 1)
out$composite_48h_strictdeath <- as.integer(out$to_invasive_48h == 1 | out$death_48h == 1)
out$composite_any            <- as.integer(out$to_invasive_any == 1 | out$death_hosp == 1)

out$followup_48h_complete <- ifelse(is.na(out$unitdischargeoffset) & is.na(out$death_time), NA_integer_,
  as.integer((!is.na(out$unitdischargeoffset) & out$unitdischargeoffset - out$niv_start_offset >= 2880) |
             out$to_invasive_48h == 1 | out$death_48h == 1))

out$icu_los_days  <- round(out$unitdischargeoffset / 1440, 3)
out$hosp_los_days <- round((out$hospitaldischargeoffset - out$hospitaladmitoffset) / 1440, 3)

# 缺失标记
out$m_apache   <- as.integer(is.na(out$apachescore))
out$m_gcs      <- as.integer(is.na(out$gcs_total))
out$m_ph       <- as.integer(is.na(out$ph_min))
out$m_paco2    <- as.integer(is.na(out$paco2_max))
out$m_pf       <- as.integer(is.na(out$pf_min))
out$m_lactate  <- as.integer(is.na(out$lactate_max))
out$m_sbp      <- as.integer(is.na(out$sbp_min))
out$m_rr       <- as.integer(is.na(out$rr_max))
out$m_spo2     <- as.integer(is.na(out$spo2_min))

# 入院诊断 6 组（按 APACHE 入院诊断文本，有序关键词映射）
map_dx <- function(x) {
  x <- ifelse(is.na(x), "", tolower(x))
  grp <- rep("术后/其他", length(x))
  grp[grepl("sepsis", x)] <- "脓毒症"
  grp[grepl("copd|emphysema|asthma|bronchitis", x) & grp == "术后/其他"] <- "COPD/哮喘"
  grp[grepl("chf|congestive heart failure|cardiogenic", x) & grp == "术后/其他"] <- "心源性肺水肿"
  grp[grepl("pneumonia|ards|respiratory|pleural|atelectasis|pneumothorax|hemothorax|pulmonary embol|aspiration|airway obstruction|apnea|restrictive lung", x) & grp == "术后/其他"] <- "肺炎/其他呼吸"
  grp[grepl("coma|cva|cerebrovascular|seizure|intracranial|encephalopath|meningitis|subarachnoid|guillain|myasthenia|neurologic", x) & grp == "术后/其他"] <- "神经系统"
  grp
}
out$dx_group <- factor(map_dx(out$apacheadmissiondx),
                       levels = c("COPD/哮喘", "心源性肺水肿", "肺炎/其他呼吸", "脓毒症", "神经系统", "术后/其他"))

# HACOR：五项按原始量表计分；优先 0–2h，五项不全时用 0–24h
score_hr  <- function(x) ifelse(is.na(x), NA_integer_, ifelse(x < 90, 0L, ifelse(x < 120, 1L, 2L)))
score_ph  <- function(x) ifelse(is.na(x), NA_integer_, ifelse(x >= 7.35, 0L, ifelse(x >= 7.30, 1L, ifelse(x >= 7.25, 2L, ifelse(x >= 7.20, 3L, 4L)))))
score_gcs <- function(x) ifelse(is.na(x), NA_integer_, ifelse(x >= 15, 0L, ifelse(x >= 13, 1L, ifelse(x >= 11, 2L, ifelse(x >= 9, 3L, 4L)))))
score_pf  <- function(x) ifelse(is.na(x), NA_integer_, ifelse(x >= 201, 0L, ifelse(x >= 176, 1L, ifelse(x >= 151, 2L, ifelse(x >= 126, 3L, ifelse(x >= 101, 4L, 5L))))))
score_rr  <- function(x) ifelse(is.na(x), NA_integer_, ifelse(x <= 30, 0L, ifelse(x <= 35, 1L, ifelse(x <= 40, 2L, ifelse(x <= 45, 3L, 4L)))))

mk_hacor <- function(hr, ph, gcs, pf, rr) {
  s <- cbind(score_hr(hr), score_ph(ph), score_gcs(gcs), score_pf(pf), score_rr(rr))
  ifelse(rowSums(is.na(s)) == 0, rowSums(s), NA_integer_)
}
out$hacor_hr_2h  <- score_hr(out$hr_2h);   out$hacor_ph_2h  <- score_ph(out$ph_hacor_2h)
out$hacor_gcs_2h <- score_gcs(out$gcs_hacor_2h); out$hacor_pf_2h <- score_pf(out$pf_hacor_2h)
out$hacor_rr_2h  <- score_rr(out$rr_2h)
out$hacor_2h <- mk_hacor(out$hr_2h, out$ph_hacor_2h, out$gcs_hacor_2h, out$pf_hacor_2h, out$rr_2h)

out$hacor_hr_24h  <- score_hr(out$hr_24h);   out$hacor_ph_24h  <- score_ph(out$ph_hacor_24h)
out$hacor_gcs_24h <- score_gcs(out$gcs_hacor_24h); out$hacor_pf_24h <- score_pf(out$pf_hacor_24h)
out$hacor_rr_24h  <- score_rr(out$rr_24h)
out$hacor_24h <- mk_hacor(out$hr_24h, out$ph_hacor_24h, out$gcs_hacor_24h, out$pf_hacor_24h, out$rr_24h)

out$hacor_score  <- ifelse(!is.na(out$hacor_2h), out$hacor_2h, out$hacor_24h)
out$hacor_window <- ifelse(!is.na(out$hacor_2h), "2h", ifelse(!is.na(out$hacor_24h), "24h", NA_character_))
out$hacor_hr  <- ifelse(!is.na(out$hacor_2h), out$hacor_hr_2h, out$hacor_hr_24h)
out$hacor_ph  <- ifelse(!is.na(out$hacor_2h), out$hacor_ph_2h, out$hacor_ph_24h)
out$hacor_gcs <- ifelse(!is.na(out$hacor_2h), out$hacor_gcs_2h, out$hacor_gcs_24h)
out$hacor_pf  <- ifelse(!is.na(out$hacor_2h), out$hacor_pf_2h, out$hacor_pf_24h)
out$hacor_rr  <- ifelse(!is.na(out$hacor_2h), out$hacor_rr_2h, out$hacor_rr_24h)

# 医院层面：NIV-first 例数与血气可得比例（只在最终分析集内计算）
hosp <- aggregate(cbind(n_in = out$in_analysis,
                        ph_avail = as.integer(out$in_analysis == 1 & !is.na(out$ph_min))) ~ hospitalid,
                  data = out, FUN = sum)
hosp$hosp_ph_avail_pct <- round(100 * hosp$ph_avail / hosp$n_in, 1)
out <- merge(out, hosp[, c("hospitalid", "n_in", "hosp_ph_avail_pct")],
             by = "hospitalid", all.x = TRUE, sort = FALSE)
names(out)[names(out) == "n_in"] <- "hosp_nivfirst_n"

# ---------------------------------------------------------------- 9. 选列、检查点、输出
keep <- c(
  "patientunitstayid","patienthealthsystemstayid","uniquepid","hospitalid",
  "unitvisitnumber","first_icu_stay","patient_unit_seq","hospitaladmitoffset",
  "niv_record","niv_start_offset","niv_pre_icu","inv_record_strict","inv_record_any",
  "inv_start_strict","inv_start_any","inv_before_niv","inv_same_time","niv_first",
  "age","age_numeric","age_lt18","age_missing","in_analysis",
  "to_invasive_48h","to_invasive_24h","to_invasive_any","to_invasive_48h_anydef",
  "death_48h","death_hosp","death_offset","composite_48h","composite_48h_strictdeath","composite_any",
  "followup_48h_complete",
  "gender","ethnicity","apacheadmissiondx","dx_group","apachescore","predictedhospitalmortality",
  "gcs_total","sbp_min","rr_max","spo2_min","ph_min","paco2_max","pf_min","lactate_max",
  "bicarbonate_first","hr_min","hr_max",
  "m_apache","m_gcs","m_ph","m_paco2","m_pf","m_lactate","m_sbp","m_rr","m_spo2",
  "icu_los_days","hosp_los_days",
  "hacor_score","hacor_window","hacor_hr","hacor_ph","hacor_gcs","hacor_pf","hacor_rr",
  "hacor_2h","hacor_24h",
  "dnr_before_niv","niv_epap_mean","niv_set_rr_mean",
  "hosp_nivfirst_n","hosp_ph_avail_pct")
out$inv_record_strict <- as.integer(!is.na(out$inv_start_strict))
out$inv_record_any    <- as.integer(!is.na(out$inv_start_any))
out <- out[, keep]

out <- mc_checkpoint(out, "纳入：C3 eICU 主分析集（NIV-first、非同刻、成人）", id = "patientunitstayid",
                     keep = out$in_analysis == 1)

dbWriteTable(con, "c3_eicu_analysis", out, overwrite = TRUE)
pq <- file.path(OUT, "eicu_niv_analysis.parquet")
dbExecute(con, sprintf("COPY c3_eicu_analysis TO '%s' (FORMAT PARQUET)", pq))

cat("eICU rows:", nrow(out), " unique stay:", length(unique(out$patientunitstayid)), "\n")
cat("in_analysis:", sum(out$in_analysis), " NIV-first:", sum(out$niv_first),
    " inv_before:", sum(out$inv_before_niv), " same:", sum(out$inv_same_time),
    " age<18:", sum(out$age_lt18), " age_missing:", sum(out$age_missing), "\n")
sub <- out[out$in_analysis == 1, ]
cat("MAIN: n=", nrow(sub), " to_invasive_48h=", sum(sub$to_invasive_48h),
    " to_invasive_24h=", sum(sub$to_invasive_24h), " any=", sum(sub$to_invasive_any),
    " death_48h=", sum(sub$death_48h), " death_hosp=", sum(sub$death_hosp), "\n")
cat("MAIN composite_48h=", sum(sub$composite_48h), " strictdeath=", sum(sub$composite_48h_strictdeath),
    " any=", sum(sub$composite_any), " followup_incomplete=", sum(sub$followup_48h_complete == 0, na.rm = TRUE), "\n")
cat("MAIN APACHE missing=", sum(sub$m_apache), " GCS missing=", sum(sub$m_gcs),
    " pH missing=", sum(sub$m_ph), " PF missing=", sum(sub$m_pf), " lactate missing=", sum(sub$m_lactate), "\n")
cat("MAIN HACOR computable=", sum(!is.na(sub$hacor_score)),
    " (2h=", sum(sub$hacor_window == "2h", na.rm = TRUE),
    ", 24h=", sum(sub$hacor_window == "24h", na.rm = TRUE), ")\n")
cat("MAIN DNR before NIV=", sum(sub$dnr_before_niv, na.rm = TRUE), "\n")
cat("Wrote:", pq, "\n")
