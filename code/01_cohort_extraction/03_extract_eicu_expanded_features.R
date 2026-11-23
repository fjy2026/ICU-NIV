#!/usr/bin/env Rscript
# C3 扩取：eICU-CRD 2.0 无创通气（NIV）队列新增变量
# 一行 = 一次 ICU 入住（patientunitstayid）。以既有 c3/eicu_niv_analysis.parquet 为基表，
# 从原始 parquet 只读补取动态变量、ROX/HACOR 组分、医院层面、来源与病房类型、
# 计划/围术期插管线索、DNR/DNI、6/12 h 结局与队列标志。
# 只读源数据；不改 c3/ 旧数据；不插补、不清洗、不建模。
suppressPackageStartupMessages({ library(DBI); library(duckdb) })
source((function(f) { for (p in c(file.path("code/00_common", f), file.path("../00_common", f), file.path("00_common", f), f)) if (file.exists(p)) return(p); f })("trace_logger.R"))

P     <- "F:/CoredataCenter/eICU/parquet"
BASE  <- "data/processed/c3/eicu_niv_analysis.parquet"
OUT   <- "data/processed/c3_expanded"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
con <- dbConnect(duckdb::duckdb())
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

base <- dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s')", BASE))
base <- mc_checkpoint(base, "读入：C3 eICU 既有分析表（15,970 次有无创记录的 ICU 入住）",
                      id = "patientunitstayid")
dbWriteTable(con, "cohort", base, overwrite = TRUE)
cat("base rows:", nrow(base), " in_analysis:", sum(base$in_analysis), "\n")

# ================================================================ 1. 动态变量（ICU 入科零点）
dyn_lab <- dbGetQuery(con, sprintf("
SELECT patientunitstayid,
  arg_min(CASE WHEN labname='lactate' AND labresult BETWEEN 0 AND 30 THEN labresult END, labresultoffset) AS lact_base_first,
  arg_max(CASE WHEN labname='lactate' AND labresult BETWEEN 0 AND 30 THEN labresult END, labresultoffset) AS lact_last_0_24,
  arg_min(CASE WHEN labname='lactate' AND labresult BETWEEN 0 AND 30 AND labresultoffset BETWEEN 1380 AND 1500 THEN labresult END, abs(labresultoffset-1440)) AS lact_at24h,
  arg_min(CASE WHEN labname='pH' AND labresult BETWEEN 6.5 AND 8.0 THEN labresult END, labresultoffset) AS ph_base_first,
  arg_max(CASE WHEN labname='pH' AND labresult BETWEEN 6.5 AND 8.0 THEN labresult END, labresultoffset) AS ph_last_0_24,
  arg_min(CASE WHEN labname='pH' AND labresult BETWEEN 6.5 AND 8.0 AND labresultoffset BETWEEN 1380 AND 1500 THEN labresult END, abs(labresultoffset-1440)) AS ph_at24h
FROM read_parquet('%s/lab.parquet') JOIN cohort USING(patientunitstayid)
WHERE labresultoffset BETWEEN 0 AND 1440
GROUP BY 1", P))

dyn_gcs <- dbGetQuery(con, sprintf("
SELECT patientunitstayid,
  arg_min(CASE WHEN nursingchartcelltypevallabel='Glasgow coma score'
                AND nursingchartcelltypevalname='GCS Total'
                AND TRY_CAST(nursingchartvalue AS DOUBLE) BETWEEN 3 AND 15
               THEN TRY_CAST(nursingchartvalue AS DOUBLE) END, nursingchartoffset) AS gcs_base_first,
  arg_max(CASE WHEN nursingchartcelltypevallabel='Glasgow coma score'
                AND nursingchartcelltypevalname='GCS Total'
                AND TRY_CAST(nursingchartvalue AS DOUBLE) BETWEEN 3 AND 15
               THEN TRY_CAST(nursingchartvalue AS DOUBLE) END, nursingchartoffset) AS gcs_last_0_24,
  arg_min(CASE WHEN nursingchartcelltypevallabel='Glasgow coma score'
                AND nursingchartcelltypevalname='GCS Total'
                AND TRY_CAST(nursingchartvalue AS DOUBLE) BETWEEN 3 AND 15
                AND nursingchartoffset BETWEEN 1380 AND 1500
               THEN TRY_CAST(nursingchartvalue AS DOUBLE) END, abs(nursingchartoffset-1440)) AS gcs_at24h
FROM read_parquet('%s/nurseCharting.parquet') JOIN cohort USING(patientunitstayid)
WHERE nursingchartoffset BETWEEN 0 AND 1440
GROUP BY 1", P))

dyn_sbp_art <- dbGetQuery(con, sprintf("
SELECT patientunitstayid,
  arg_min(CASE WHEN systemicsystolic BETWEEN 20 AND 300 THEN systemicsystolic END, observationoffset) AS sbp_art_base_first,
  arg_min(CASE WHEN systemicsystolic BETWEEN 20 AND 300 AND observationoffset BETWEEN 1380 AND 1500 THEN systemicsystolic END, abs(observationoffset-1440)) AS sbp_art_at24h
FROM read_parquet('%s/vitalPeriodic.parquet') JOIN cohort USING(patientunitstayid)
WHERE observationoffset BETWEEN 0 AND 1440
GROUP BY 1", P))

dyn_sbp_ni <- dbGetQuery(con, sprintf("
SELECT patientunitstayid,
  arg_min(CASE WHEN noninvasivesystolic BETWEEN 20 AND 300 THEN noninvasivesystolic END, observationoffset) AS sbp_ni_base_first,
  arg_min(CASE WHEN noninvasivesystolic BETWEEN 20 AND 300 AND observationoffset BETWEEN 1380 AND 1500 THEN noninvasivesystolic END, abs(observationoffset-1440)) AS sbp_ni_at24h
FROM read_parquet('%s/vitalAperiodic.parquet') JOIN cohort USING(patientunitstayid)
WHERE observationoffset BETWEEN 0 AND 1440
GROUP BY 1", P))

# 血气可得性（供医院层面比例）
bg_any <- dbGetQuery(con, sprintf("
SELECT patientunitstayid,
  MAX(CASE WHEN (labname='pH' AND labresult BETWEEN 6.5 AND 8.0)
             OR (labname='paCO2' AND labresult BETWEEN 5 AND 200)
             OR (labname='paO2' AND labresult BETWEEN 10 AND 700)
           THEN 1 ELSE 0 END) AS bg_any_0_24
FROM read_parquet('%s/lab.parquet') JOIN cohort USING(patientunitstayid)
WHERE labresultoffset BETWEEN 0 AND 1440
GROUP BY 1", P))

# ================================================================ 2. ROX 组分（NIV 零点，±1 h）
rox_spo2_rr <- dbGetQuery(con, sprintf("
WITH c AS (SELECT patientunitstayid, niv_start_offset FROM cohort)
SELECT c.patientunitstayid,
  arg_min(CASE WHEN v.sao2 BETWEEN 0 AND 100 AND v.observationoffset BETWEEN c.niv_start_offset+60  AND c.niv_start_offset+180  THEN v.sao2 END, abs(v.observationoffset-(c.niv_start_offset+120)))  AS spo2_2h,
  arg_min(CASE WHEN v.sao2 BETWEEN 0 AND 100 AND v.observationoffset BETWEEN c.niv_start_offset+660 AND c.niv_start_offset+780 THEN v.sao2 END, abs(v.observationoffset-(c.niv_start_offset+720)))  AS spo2_12h,
  arg_min(CASE WHEN v.sao2 BETWEEN 0 AND 100 AND v.observationoffset BETWEEN c.niv_start_offset+1380 AND c.niv_start_offset+1500 THEN v.sao2 END, abs(v.observationoffset-(c.niv_start_offset+1440))) AS spo2_24h,
  arg_min(CASE WHEN v.respiration BETWEEN 0 AND 80 AND v.observationoffset BETWEEN c.niv_start_offset+60  AND c.niv_start_offset+180  THEN v.respiration END, abs(v.observationoffset-(c.niv_start_offset+120)))  AS rr_2h,
  arg_min(CASE WHEN v.respiration BETWEEN 0 AND 80 AND v.observationoffset BETWEEN c.niv_start_offset+660 AND c.niv_start_offset+780 THEN v.respiration END, abs(v.observationoffset-(c.niv_start_offset+720)))  AS rr_12h,
  arg_min(CASE WHEN v.respiration BETWEEN 0 AND 80 AND v.observationoffset BETWEEN c.niv_start_offset+1380 AND c.niv_start_offset+1500 THEN v.respiration END, abs(v.observationoffset-(c.niv_start_offset+1440))) AS rr_24h
FROM read_parquet('%s/vitalPeriodic.parquet') v JOIN c USING(patientunitstayid)
WHERE v.observationoffset BETWEEN c.niv_start_offset+60 AND c.niv_start_offset+1500
GROUP BY 1", P))

rox_fio2_resp <- dbGetQuery(con, sprintf("
WITH c AS (SELECT patientunitstayid, niv_start_offset FROM cohort),
f AS (
  SELECT r.patientunitstayid, r.respchartoffset AS off,
         TRY_CAST(r.respchartvalue AS DOUBLE) AS fio2,
         CASE when r.respchartvaluelabel='FiO2' then 1
              when r.respchartvaluelabel='FIO2 (%%)' then 2
              when r.respchartvaluelabel='Set Fraction of Inspired Oxygen (FIO2)' then 3 else 4 end AS pri
  FROM read_parquet('%s/respiratoryCharting.parquet') r JOIN c USING(patientunitstayid)
  WHERE r.respchartvaluelabel IN ('FiO2','FIO2 (%%)','O2 Percentage','Set Fraction of Inspired Oxygen (FIO2)')
    AND TRY_CAST(r.respchartvalue AS DOUBLE) BETWEEN 21 AND 100
)
SELECT c.patientunitstayid,
  arg_min(CASE WHEN f.fio2 IS NOT NULL AND f.off BETWEEN c.niv_start_offset+60  AND c.niv_start_offset+180  THEN f.fio2 END, abs(f.off-(c.niv_start_offset+120))*10+pri)  AS fio2r_2h,
  arg_min(CASE WHEN f.fio2 IS NOT NULL AND f.off BETWEEN c.niv_start_offset+660 AND c.niv_start_offset+780 THEN f.fio2 END, abs(f.off-(c.niv_start_offset+720))*10+pri)  AS fio2r_12h,
  arg_min(CASE WHEN f.fio2 IS NOT NULL AND f.off BETWEEN c.niv_start_offset+1380 AND c.niv_start_offset+1500 THEN f.fio2 END, abs(f.off-(c.niv_start_offset+1440))*10+pri) AS fio2r_24h
FROM c LEFT JOIN f USING(patientunitstayid)
GROUP BY 1", P))

rox_fio2_lab <- dbGetQuery(con, sprintf("
WITH c AS (SELECT patientunitstayid, niv_start_offset FROM cohort)
SELECT c.patientunitstayid,
  arg_min(CASE WHEN l.labresult BETWEEN 21 AND 100 AND l.labresultoffset BETWEEN c.niv_start_offset+60  AND c.niv_start_offset+180  THEN l.labresult END, abs(l.labresultoffset-(c.niv_start_offset+120)))  AS fio2l_2h,
  arg_min(CASE WHEN l.labresult BETWEEN 21 AND 100 AND l.labresultoffset BETWEEN c.niv_start_offset+660 AND c.niv_start_offset+780 THEN l.labresult END, abs(l.labresultoffset-(c.niv_start_offset+720)))  AS fio2l_12h,
  arg_min(CASE WHEN l.labresult BETWEEN 21 AND 100 AND l.labresultoffset BETWEEN c.niv_start_offset+1380 AND c.niv_start_offset+1500 THEN l.labresult END, abs(l.labresultoffset-(c.niv_start_offset+1440))) AS fio2l_24h
FROM read_parquet('%s/lab.parquet') l JOIN c USING(patientunitstayid)
WHERE l.labname='FiO2'
GROUP BY 1", P))

# ================================================================ 3. HACOR 组分原始值（NIV 零点）
hacor_lab <- dbGetQuery(con, sprintf("
WITH c AS (SELECT patientunitstayid, niv_start_offset FROM cohort),
l AS (
  SELECT l.patientunitstayid, c.niv_start_offset, l.labresultoffset,
    MAX(CASE WHEN l.labname='paO2' AND l.labresult BETWEEN 10 AND 700 THEN l.labresult END) AS po2,
    MAX(CASE WHEN l.labname='FiO2' AND l.labresult BETWEEN 21 AND 100 THEN l.labresult END) AS fio2,
    MIN(CASE WHEN l.labname='pH' AND l.labresult BETWEEN 6.5 AND 8.0 THEN l.labresult END) AS ph
  FROM read_parquet('%s/lab.parquet') l JOIN c USING(patientunitstayid)
  WHERE l.labresultoffset BETWEEN c.niv_start_offset AND c.niv_start_offset+1440
  GROUP BY 1,2,3
)
SELECT patientunitstayid,
  MIN(CASE WHEN labresultoffset <= niv_start_offset+120 THEN ph END) AS ph_2h,
  MIN(ph) AS ph_24h,
  MIN(CASE WHEN labresultoffset <= niv_start_offset+120 THEN po2 END) AS pao2_2h,
  MIN(po2) AS pao2_24h,
  arg_min(CASE WHEN labresultoffset <= niv_start_offset+120 THEN fio2 END, abs(labresultoffset-(niv_start_offset+60))) AS hacor_fio2_2h,
  arg_min(CASE WHEN fio2 IS NOT NULL THEN fio2 END, abs(labresultoffset-(niv_start_offset+720))) AS hacor_fio2_24h,
  MIN(CASE WHEN labresultoffset <= niv_start_offset+120 AND po2 IS NOT NULL AND fio2 IS NOT NULL THEN po2/(fio2/100.0) END) AS pf_2h,
  MIN(CASE WHEN po2 IS NOT NULL AND fio2 IS NOT NULL THEN po2/(fio2/100.0) END) AS pf_24h
FROM l GROUP BY 1", P))

hacor_vit <- dbGetQuery(con, sprintf("
WITH c AS (SELECT patientunitstayid, niv_start_offset FROM cohort)
SELECT c.patientunitstayid,
  MAX(CASE WHEN v.heartrate BETWEEN 20 AND 300 AND v.observationoffset BETWEEN c.niv_start_offset AND c.niv_start_offset+120  THEN v.heartrate END) AS hr_hacor_2h,
  MAX(CASE WHEN v.heartrate BETWEEN 20 AND 300 AND v.observationoffset BETWEEN c.niv_start_offset AND c.niv_start_offset+1440 THEN v.heartrate END) AS hr_hacor_24h,
  MAX(CASE WHEN v.respiration BETWEEN 0 AND 80 AND v.observationoffset BETWEEN c.niv_start_offset AND c.niv_start_offset+120  THEN v.respiration END) AS rr_hacor_2h,
  MAX(CASE WHEN v.respiration BETWEEN 0 AND 80 AND v.observationoffset BETWEEN c.niv_start_offset AND c.niv_start_offset+1440 THEN v.respiration END) AS rr_hacor_24h
FROM read_parquet('%s/vitalPeriodic.parquet') v JOIN c USING(patientunitstayid)
WHERE v.observationoffset BETWEEN c.niv_start_offset AND c.niv_start_offset+1440
GROUP BY 1", P))

hacor_gcs <- dbGetQuery(con, sprintf("
WITH c AS (SELECT patientunitstayid, niv_start_offset FROM cohort)
SELECT c.patientunitstayid,
  MIN(CASE WHEN n.nursingchartoffset BETWEEN c.niv_start_offset AND c.niv_start_offset+120
            AND n.nursingchartcelltypevallabel='Glasgow coma score'
            AND n.nursingchartcelltypevalname='GCS Total'
            AND TRY_CAST(n.nursingchartvalue AS DOUBLE) BETWEEN 3 AND 15
           THEN TRY_CAST(n.nursingchartvalue AS DOUBLE) END) AS gcs_2h,
  MIN(CASE WHEN n.nursingchartoffset BETWEEN c.niv_start_offset AND c.niv_start_offset+1440
            AND n.nursingchartcelltypevallabel='Glasgow coma score'
            AND n.nursingchartcelltypevalname='GCS Total'
            AND TRY_CAST(n.nursingchartvalue AS DOUBLE) BETWEEN 3 AND 15
           THEN TRY_CAST(n.nursingchartvalue AS DOUBLE) END) AS gcs_24h
FROM read_parquet('%s/nurseCharting.parquet') n JOIN c USING(patientunitstayid)
GROUP BY 1", P))

# ================================================================ 4. 医院层面、来源与病房类型
hosp <- dbGetQuery(con, sprintf("SELECT hospitalid, numbedscategory, teachingstatus, region
FROM read_parquet('%s/hospital.parquet')", P))
names(hosp) <- c("hospitalid", "hosp_beds_category", "hosp_teaching", "hosp_region")

src <- dbGetQuery(con, sprintf("SELECT patientunitstayid, unitadmitsource, unittype, hospitaladmitsource
FROM read_parquet('%s/patient.parquet')", P))

# ================================================================ 5. 计划/围术期插管线索
plan <- dbGetQuery(con, sprintf("
WITH c AS (SELECT patientunitstayid, niv_start_offset FROM cohort),
t AS (
  SELECT t.patientunitstayid, t.treatmentoffset, lower(t.treatmentstring) AS ts
  FROM read_parquet('%s/treatment.parquet') t JOIN c USING(patientunitstayid)
)
SELECT c.patientunitstayid,
  MAX(CASE WHEN t.ts LIKE 'surgery|%%' THEN 1 ELSE 0 END) AS tx_surgery_any,
  MAX(CASE WHEN t.ts LIKE 'surgery|%%' AND t.ts LIKE '%%mechanical ventilation%%'
            AND t.ts NOT LIKE '%%non-invasive ventilation%%' AND t.treatmentoffset > c.niv_start_offset
           THEN 1 ELSE 0 END) AS inv_surgery_path,
  MAX(CASE WHEN t.ts LIKE 'surgery|%%' AND t.ts LIKE '%%mechanical ventilation%%'
            AND t.ts NOT LIKE '%%non-invasive ventilation%%' AND t.treatmentoffset > c.niv_start_offset
            AND t.treatmentoffset <= c.niv_start_offset+2880
           THEN 1 ELSE 0 END) AS inv_surgery_path_48h
FROM c LEFT JOIN t USING(patientunitstayid)
GROUP BY 1", P))

planned_cp <- dbGetQuery(con, sprintf("
SELECT patientunitstayid,
  MAX(CASE WHEN lower(cplgroup)='planned procedures' THEN 1 ELSE 0 END) AS careplan_planned_procedure
FROM read_parquet('%s/carePlanGeneral.parquet')
GROUP BY 1", P))

# ================================================================ 6. DNR / DNI
lim <- dbGetQuery(con, sprintf("
WITH c AS (SELECT patientunitstayid, niv_start_offset FROM cohort)
SELECT c.patientunitstayid,
  MAX(CASE WHEN lower(g.cplgroup)='care limitation' THEN 1 ELSE 0 END) AS has_careplan_record,
  MAX(CASE WHEN lower(g.cplgroup)='care limitation' AND g.cplitemoffset <= c.niv_start_offset THEN 1 ELSE 0 END) AS carelim_any_before_niv,
  MAX(CASE WHEN lower(g.cplitemvalue)='no intubation' AND g.cplitemoffset <= c.niv_start_offset THEN 1 ELSE 0 END) AS dni_before_niv,
  MAX(CASE WHEN lower(g.cplitemvalue)='comfort measures only' AND g.cplitemoffset <= c.niv_start_offset THEN 1 ELSE 0 END) AS comfort_before_niv,
  MAX(CASE WHEN lower(g.cplitemvalue)='do not resuscitate' THEN 1 ELSE 0 END) AS dnr_any,
  MAX(CASE WHEN lower(g.cplitemvalue)='no intubation' THEN 1 ELSE 0 END) AS dni_any,
  MAX(CASE WHEN lower(g.cplitemvalue)='comfort measures only' THEN 1 ELSE 0 END) AS comfort_any
FROM c LEFT JOIN read_parquet('%s/carePlanGeneral.parquet') g USING(patientunitstayid)
GROUP BY 1", P))

# ================================================================ 7. 6 / 12 h 结局（严格有创，t1 > t0 最早）
inv_win <- dbGetQuery(con, sprintf("
WITH c AS (SELECT patientunitstayid, niv_start_offset FROM cohort),
t AS (
  SELECT t.patientunitstayid, MIN(t.treatmentoffset) AS t1
  FROM read_parquet('%s/treatment.parquet') t JOIN c USING(patientunitstayid)
  WHERE lower(t.treatmentstring) LIKE '%%mechanical ventilation%%'
    AND lower(t.treatmentstring) NOT LIKE '%%non-invasive ventilation%%'
    AND t.treatmentoffset > c.niv_start_offset
  GROUP BY 1
)
SELECT patientunitstayid, t1 FROM t", P))

# ================================================================ 8. 合并
for (nm in c("dyn_lab","dyn_gcs","dyn_sbp_art","dyn_sbp_ni","bg_any","rox_spo2_rr",
             "rox_fio2_resp","rox_fio2_lab","hacor_lab","hacor_vit","hacor_gcs",
             "src","plan","planned_cp","lim","inv_win")) {
  assign(nm, mc_checkpoint(get(nm), paste0("并入：", nm), id = "patientunitstayid", side = TRUE))
}
hosp <- mc_checkpoint(hosp, "并入：医院属性表", id = "hospitalid", side = TRUE)

out <- base
for (blk in list(dyn_lab, dyn_gcs, dyn_sbp_art, dyn_sbp_ni, bg_any, rox_spo2_rr,
                 rox_fio2_resp, rox_fio2_lab, hacor_lab, hacor_vit, hacor_gcs,
                 src, plan, planned_cp, lim, inv_win)) {
  out <- merge(out, blk, by = "patientunitstayid", all.x = TRUE, sort = FALSE)
}
out <- merge(out, hosp, by = "hospitalid", all.x = TRUE, sort = FALSE)
stopifnot(!anyDuplicated(out$patientunitstayid))

# --- 动态变量派生 ---
out$lact_delta_at24h    <- out$lact_at24h - out$lact_base_first
out$lact_delta_lastfirst<- out$lact_last_0_24 - out$lact_base_first
out$ph_delta_at24h      <- out$ph_at24h - out$ph_base_first
out$ph_delta_lastfirst  <- out$ph_last_0_24 - out$ph_base_first
out$gcs_delta_at24h     <- out$gcs_at24h - out$gcs_base_first
out$gcs_delta_lastfirst <- out$gcs_last_0_24 - out$gcs_base_first
out$sbp_art_delta_at24h <- out$sbp_art_at24h - out$sbp_art_base_first
out$sbp_ni_delta_at24h  <- out$sbp_ni_at24h - out$sbp_ni_base_first
out$sbp_base_first <- ifelse(!is.na(out$sbp_art_base_first), out$sbp_art_base_first, out$sbp_ni_base_first)
out$sbp_at24h      <- ifelse(!is.na(out$sbp_art_at24h), out$sbp_art_at24h, out$sbp_ni_at24h)
out$sbp_delta_at24h<- out$sbp_at24h - out$sbp_base_first
for (v in c("lact","ph","gcs","sbp")) {
  out[[paste0("m_", v, "_base_first")]] <- as.integer(is.na(out[[paste0(v, "_base_first")]]))
  out[[paste0("m_", v, "_at24h")]]      <- as.integer(is.na(out[[paste0(v, "_at24h")]]))
}

# --- ROX ---
cat("ROX cols:", paste(grep("^(spo2|fio2|rr_|rox)", names(out), value = TRUE), collapse = ", "), "\n")
cat("lens:", length(out[["spo2_2h"]]), length(out[["fio2r_2h"]]), length(out[["fio2l_2h"]]),
    length(out[["rr_2h"]]), "nrow", nrow(out), "\n")
for (t in c("2h","12h","24h")) {
  r <- out[[paste0("fio2r_", t)]]; l <- out[[paste0("fio2l_", t)]]
  out[[paste0("fio2_", t)]]    <- ifelse(!is.na(r), r, l)
  out[[paste0("fio2src_", t)]] <- ifelse(!is.na(r), "resp", ifelse(!is.na(l), "lab", NA_character_))
  out[[paste0("rox_", t)]] <- ifelse(!is.na(out[[paste0("spo2_", t)]]) &
                                       !is.na(out[[paste0("fio2_", t)]]) &
                                       !is.na(out[[paste0("rr_", t)]]) &
                                       out[[paste0("rr_", t)]] > 0,
                                     (out[[paste0("spo2_", t)]] / (out[[paste0("fio2_", t)]]/100)) / out[[paste0("rr_", t)]],
                                     NA_real_)
  out[[paste0("rox_computable_", t)]] <- as.integer(!is.na(out[[paste0("rox_", t)]]))
  for (c0 in c("spo2","fio2","rr","rox")) {
    out[[paste0("m_", c0, "_", t)]] <- as.integer(is.na(out[[paste0(c0, "_", t)]]))
  }
}

# --- HACOR 组分 ---
out$hacor2h_hr <- out$hr_hacor_2h;  out$hacor2h_ph <- out$ph_2h
out$hacor2h_pao2 <- out$pao2_2h; out$hacor2h_fio2 <- out$hacor_fio2_2h
out$hacor2h_pf <- out$pf_2h;  out$hacor2h_gcs <- out$gcs_2h; out$hacor2h_rr <- out$rr_hacor_2h
out$hacor24h_hr <- out$hr_hacor_24h; out$hacor24h_ph <- out$ph_24h
out$hacor24h_pao2 <- out$pao2_24h; out$hacor24h_fio2 <- out$hacor_fio2_24h
out$hacor24h_pf <- out$pf_24h; out$hacor24h_gcs <- out$gcs_24h; out$hacor24h_rr <- out$rr_hacor_24h

score_hr  <- function(x) ifelse(is.na(x), NA_integer_, ifelse(x < 90, 0L, ifelse(x < 120, 1L, 2L)))
score_ph  <- function(x) ifelse(is.na(x), NA_integer_, ifelse(x >= 7.35, 0L, ifelse(x >= 7.30, 1L, ifelse(x >= 7.25, 2L, ifelse(x >= 7.20, 3L, 4L)))))
score_gcs <- function(x) ifelse(is.na(x), NA_integer_, ifelse(x >= 15, 0L, ifelse(x >= 13, 1L, ifelse(x >= 11, 2L, ifelse(x >= 9, 3L, 4L)))))
score_pf  <- function(x) ifelse(is.na(x), NA_integer_, ifelse(x >= 201, 0L, ifelse(x >= 176, 1L, ifelse(x >= 151, 2L, ifelse(x >= 126, 3L, ifelse(x >= 101, 4L, 5L))))))
score_rr  <- function(x) ifelse(is.na(x), NA_integer_, ifelse(x <= 30, 0L, ifelse(x <= 35, 1L, ifelse(x <= 40, 2L, ifelse(x <= 45, 3L, 4L)))))
mk_total <- function(hr, ph, gcs, pf, rr) {
  s <- cbind(score_hr(hr), score_ph(ph), score_gcs(gcs), score_pf(pf), score_rr(rr))
  ifelse(rowSums(is.na(s)) == 0, rowSums(s), NA_integer_)
}
out$hacor2h_total  <- mk_total(out$hacor2h_hr, out$hacor2h_ph, out$hacor2h_gcs, out$hacor2h_pf, out$hacor2h_rr)
out$hacor24h_total <- mk_total(out$hacor24h_hr, out$hacor24h_ph, out$hacor24h_gcs, out$hacor24h_pf, out$hacor24h_rr)
out$hacor2h_computable  <- as.integer(!is.na(out$hacor2h_total))
out$hacor24h_computable <- as.integer(!is.na(out$hacor24h_total))

# --- 医院层面派生 ---
main <- out$in_analysis == 1
hosp_agg <- aggregate(cbind(n_in = as.integer(main),
                            dnr = as.integer(main & !is.na(out$dnr_before_niv) & out$dnr_before_niv == 1),
                            bg  = as.integer(main & !is.na(out$bg_any_0_24) & out$bg_any_0_24 == 1)) ~ hospitalid,
                      data = out, FUN = sum)
hosp_agg$hosp_dnr_rate <- ifelse(hosp_agg$n_in > 0, hosp_agg$dnr / hosp_agg$n_in, NA_real_)
hosp_agg$hosp_bloodgas_avail_pct <- ifelse(hosp_agg$n_in > 0, round(100 * hosp_agg$bg / hosp_agg$n_in, 1), NA_real_)
out <- merge(out, hosp_agg[, c("hospitalid","hosp_dnr_rate","hosp_bloodgas_avail_pct")],
             by = "hospitalid", all.x = TRUE, sort = FALSE)
stopifnot(!anyDuplicated(out$patientunitstayid))

# --- 来源与病房标志 ---
out$src_ed   <- as.integer(!is.na(out$unitadmitsource) & out$unitadmitsource == "Emergency Department")
out$src_floor<- as.integer(!is.na(out$unitadmitsource) & grepl("Floor", out$unitadmitsource))
out$src_or   <- as.integer(!is.na(out$unitadmitsource) & out$unitadmitsource == "Operating Room")
out$src_pacu_recovery <- as.integer(!is.na(out$unitadmitsource) & out$unitadmitsource %in% c("PACU","Recovery Room"))
out$src_other_icu <- as.integer(!is.na(out$unitadmitsource) & out$unitadmitsource %in% c("Other ICU","ICU","ICU to SDU"))
out$src_other_hospital <- as.integer(!is.na(out$unitadmitsource) & out$unitadmitsource == "Other Hospital")
out$src_direct <- as.integer(!is.na(out$unitadmitsource) & out$unitadmitsource == "Direct Admit")
out$hosp_src_or_pacu <- as.integer(!is.na(out$hospitaladmitsource) & out$hospitaladmitsource %in% c("Operating Room","PACU","Recovery Room"))
out$icu_medical  <- as.integer(!is.na(out$unittype) & out$unittype %in% c("MICU","CCU-CTICU","Cardiac ICU","Neuro ICU"))
out$icu_surgical <- as.integer(!is.na(out$unittype) & out$unittype %in% c("SICU","CSICU","CTICU"))
out$icu_medsurg  <- as.integer(!is.na(out$unittype) & out$unittype == "Med-Surg ICU")
out$flag_or_pacu_source <- as.integer(out$src_or == 1 | out$src_pacu_recovery == 1)

# --- 计划/围术期插管线索 ---
out$tx_surgery_any <- ifelse(is.na(out$tx_surgery_any), 0L, out$tx_surgery_any)
out$inv_surgery_path <- ifelse(is.na(out$inv_surgery_path), 0L, out$inv_surgery_path)
out$inv_surgery_path_48h <- ifelse(is.na(out$inv_surgery_path_48h), 0L, out$inv_surgery_path_48h)
out$careplan_planned_procedure <- ifelse(is.na(out$careplan_planned_procedure), 0L, out$careplan_planned_procedure)
out$inv_with_or_source <- as.integer(out$flag_or_pacu_source == 1 & out$to_invasive_48h == 1)
out$inv_planned_context <- as.integer(out$to_invasive_48h == 1 &
  (out$inv_surgery_path_48h == 1 | out$careplan_planned_procedure == 1 | out$inv_with_or_source == 1))

# --- DNR / DNI ---
for (v in c("has_careplan_record","carelim_any_before_niv","dni_before_niv","comfort_before_niv",
            "dnr_any","dni_any","comfort_any")) {
  out[[v]] <- ifelse(is.na(out[[v]]), 0L, out[[v]])
}
out$flag_dnr_before_niv <- as.integer(!is.na(out$dnr_before_niv) & out$dnr_before_niv == 1)

# --- 6 / 12 h 结局 ---
out$t1_strict_after_niv <- out$t1
out$to_invasive_6h <- as.integer(!is.na(out$t1) & out$t1 - out$niv_start_offset > 0 &
                                   out$t1 - out$niv_start_offset <= 360 &
                                   (is.na(out$death_offset) | out$t1 < out$death_offset))
out$to_invasive_12h <- as.integer(!is.na(out$t1) & out$t1 - out$niv_start_offset > 0 &
                                    out$t1 - out$niv_start_offset <= 720 &
                                    (is.na(out$death_offset) | out$t1 < out$death_offset))
out$composite_6h  <- as.integer(out$to_invasive_6h == 1 | out$death_hosp == 1)
out$composite_12h <- as.integer(out$to_invasive_12h == 1 | out$death_hosp == 1)

# --- 队列标志 ---
out$cohort_first_icu      <- as.integer(out$first_icu_stay == 1)
out$cohort_patient_unique <- as.integer(out$patient_unit_seq == 1)
out$cohort_medical_icu    <- out$icu_medical
out$cohort_surgical_icu   <- out$icu_surgical
out$strict_invasive_def   <- 1L
out$cohort_excl_or_pacu <- as.integer(out$in_analysis == 1 & out$flag_or_pacu_source == 0)
out$cohort_excl_dnr     <- as.integer(out$in_analysis == 1 & out$flag_dnr_before_niv == 0)
out$cohort_strict       <- as.integer(out$in_analysis == 1 & out$flag_or_pacu_source == 0 &
                                        out$flag_dnr_before_niv == 0 & out$first_icu_stay == 1)

# ================================================================ 9. 选列、检查点、输出
keep_old <- names(base)
new_cols <- c(
  "lact_base_first","lact_last_0_24","lact_at24h","lact_delta_at24h","lact_delta_lastfirst",
  "ph_base_first","ph_last_0_24","ph_at24h","ph_delta_at24h","ph_delta_lastfirst",
  "gcs_base_first","gcs_last_0_24","gcs_at24h","gcs_delta_at24h","gcs_delta_lastfirst",
  "sbp_art_base_first","sbp_art_at24h","sbp_art_delta_at24h",
  "sbp_ni_base_first","sbp_ni_at24h","sbp_ni_delta_at24h",
  "sbp_base_first","sbp_at24h","sbp_delta_at24h",
  "m_lact_base_first","m_lact_at24h","m_ph_base_first","m_ph_at24h",
  "m_gcs_base_first","m_gcs_at24h","m_sbp_base_first","m_sbp_at24h",
  paste0(rep(c("spo2_","fio2_","fio2src_","rr_","rox_","m_spo2_","m_fio2_","m_rr_","m_rox_","rox_computable_"), each=3), c("2h","12h","24h")),
  "hacor2h_hr","hacor2h_ph","hacor2h_pao2","hacor2h_fio2","hacor2h_pf","hacor2h_gcs","hacor2h_rr","hacor2h_total","hacor2h_computable",
  "hacor24h_hr","hacor24h_ph","hacor24h_pao2","hacor24h_fio2","hacor24h_pf","hacor24h_gcs","hacor24h_rr","hacor24h_total","hacor24h_computable",
  "hosp_beds_category","hosp_teaching","hosp_region","hosp_nivfirst_n","hosp_dnr_rate",
  "hosp_ph_avail_pct","hosp_bloodgas_avail_pct",
  "unitadmitsource","unittype","hospitaladmitsource",
  "src_ed","src_floor","src_or","src_pacu_recovery","src_other_icu","src_other_hospital","src_direct",
  "hosp_src_or_pacu","icu_medical","icu_surgical","icu_medsurg",
  "inv_surgery_path","inv_surgery_path_48h","tx_surgery_any","careplan_planned_procedure",
  "inv_with_or_source","inv_planned_context",
  "dni_before_niv","comfort_before_niv","carelim_any_before_niv","dnr_any","dni_any","comfort_any","has_careplan_record",
  "to_invasive_6h","to_invasive_12h","composite_6h","composite_12h",
  "cohort_first_icu","cohort_patient_unique","cohort_medical_icu","cohort_surgical_icu",
  "flag_or_pacu_source","flag_dnr_before_niv","strict_invasive_def",
  "cohort_excl_or_pacu","cohort_excl_dnr","cohort_strict")
new_cols <- unique(new_cols[new_cols %in% names(out)])
final <- out[, c(keep_old, new_cols)]

final <- mc_checkpoint(final, paste0("输出：C3 eICU 扩取表全表（保留 ", nrow(final), " 行；in_analysis = 1 的主分析集为 ", sum(final$in_analysis), " 行）"),
                       id = "patientunitstayid")

dbWriteTable(con, "c3_eicu_expanded", final, overwrite = TRUE)
pq <- file.path(OUT, "eicu_niv_expanded.parquet")
dbExecute(con, sprintf("COPY c3_eicu_expanded TO '%s' (FORMAT PARQUET)", pq))

# ================================================================ 10. 汇报数字（全部由本脚本算出）
sub <- final[final$in_analysis == 1, ]
cat("\n=== eICU expanded ===\n")
cat("rows:", nrow(final), " unique:", length(unique(final$patientunitstayid)),
    " main:", nrow(sub), "\n")
cat("dynamic 0-24h:", "lact", sum(!is.na(sub$lact_base_first)), "ph", sum(!is.na(sub$ph_base_first)),
    "gcs", sum(!is.na(sub$gcs_base_first)), "sbp", sum(!is.na(sub$sbp_base_first)), "\n")
cat("dynamic 24h:", "lact", sum(!is.na(sub$lact_at24h)), "ph", sum(!is.na(sub$ph_at24h)),
    "gcs", sum(!is.na(sub$gcs_at24h)), "sbp", sum(!is.na(sub$sbp_at24h)), "\n")
for (t in c("2h","12h","24h")) {
  cat("ROX", t, ":", sum(sub[[paste0("rox_computable_", t)]]), " spo2", sum(!is.na(sub[[paste0("spo2_",t)]])),
      " fio2", sum(!is.na(sub[[paste0("fio2_",t)]])), " rr", sum(!is.na(sub[[paste0("rr_",t)]])), "\n")
}
cat("HACOR:", "2h", sum(sub$hacor2h_computable), "24h", sum(sub$hacor24h_computable),
    "either", sum(sub$hacor2h_computable==1 | sub$hacor24h_computable==1), "\n")
cat("hospital fields:", "beds NA", sum(is.na(sub$hosp_beds_category)), "region NA", sum(is.na(sub$hosp_region)),
    "teaching NA", sum(is.na(sub$hosp_teaching)), "hospitals", length(unique(sub$hospitalid)), "\n")
cat("source: OR", sum(sub$src_or), "PACU/recovery", sum(sub$src_pacu_recovery),
    "OR/PACU flag", sum(sub$flag_or_pacu_source), "medical ICU", sum(sub$icu_medical),
    "surgical ICU", sum(sub$icu_surgical), "med-surg", sum(sub$icu_medsurg), "\n")
cat("planned: inv_surgery_path", sum(sub$inv_surgery_path), "48h", sum(sub$inv_surgery_path_48h),
    "tx_surgery_any", sum(sub$tx_surgery_any), "careplan_planned", sum(sub$careplan_planned_procedure),
    "inv_with_or", sum(sub$inv_with_or_source), "inv_planned_context", sum(sub$inv_planned_context), "\n")
cat("DNR/DNI: dnr_before", sum(sub$dnr_before_niv, na.rm=TRUE), "dni_before", sum(sub$dni_before_niv),
    "comfort_before", sum(sub$comfort_before_niv), "dnr_any", sum(sub$dnr_any),
    "dni_any", sum(sub$dni_any), "careplan record 0", sum(sub$has_careplan_record==0), "\n")
cat("outcomes: 6h", sum(sub$to_invasive_6h, na.rm=TRUE), "12h", sum(sub$to_invasive_12h, na.rm=TRUE),
    "24h", sum(sub$to_invasive_24h, na.rm=TRUE), "48h", sum(sub$to_invasive_48h, na.rm=TRUE),
    "composite6h", sum(sub$composite_6h, na.rm=TRUE), "composite12h", sum(sub$composite_12h, na.rm=TRUE), "\n")
cat("cohort flags: first_icu", sum(sub$cohort_first_icu), "patient_unique", sum(sub$cohort_patient_unique),
    "excl_or_pacu", sum(sub$cohort_excl_or_pacu), "excl_dnr", sum(sub$cohort_excl_dnr),
    "strict", sum(sub$cohort_strict), "\n")
cat("Wrote:", pq, "\n")
