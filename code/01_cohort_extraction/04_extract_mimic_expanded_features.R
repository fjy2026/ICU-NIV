#!/usr/bin/env Rscript
# C3 扩取：MIMIC-IV v3.1 无创通气（NIV）外部验证队列新增变量
# 一行 = 一次 ICU 入住（stay_id）。以既有 c3/mimic_niv_analysis.parquet 为基表，
# 从 mimic.duckdb 只读补取动态变量、ROX/HACOR 组分、来源与病房类型、
# 计划/围术期插管线索、DNR/DNI、6/12 h 结局与队列标志。
# 只读源数据；不改 c3/ 旧数据；不插补、不清洗、不建模。
suppressPackageStartupMessages({ library(DBI); library(duckdb) })
source((function(f) { for (p in c(file.path("code/00_common", f), file.path("../00_common", f), file.path("00_common", f), f)) if (file.exists(p)) return(p); f })("trace_logger.R"))

DB  <- "F:/CoredataCenter/MIMIC-IV/mimic.duckdb"
BASE<- "data/processed/c3/mimic_niv_analysis.parquet"
OUT <- "data/processed/c3_expanded"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
con <- dbConnect(duckdb::duckdb(), dbdir = DB, read_only = TRUE)
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

base <- dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s')", BASE))
base <- mc_checkpoint(base, "读入：C3 MIMIC 既有分析表（2,696 次有无创操作的 ICU 入住）", id = "stay_id")
cat("base rows:", nrow(base), " in_analysis:", sum(base$in_analysis), "\n")

T <- "(SELECT o.stay_id, o.hadm_id, o.subject_id, i.intime FROM read_parquet('data/processed/c3/mimic_niv_analysis.parquet') o JOIN main.icustays i USING(stay_id))"

# ================================================================ 1. 动态变量（ICU intime 零点）
dyn_bg <- dbGetQuery(con, sprintf("
WITH t AS %s
SELECT t.stay_id,
  arg_min(CASE WHEN b.lactate BETWEEN 0 AND 30 THEN b.lactate END, b.charttime) AS lact_base_first,
  arg_max(CASE WHEN b.lactate BETWEEN 0 AND 30 THEN b.lactate END, b.charttime) AS lact_last_0_24,
  arg_min(CASE WHEN b.lactate BETWEEN 0 AND 30 AND b.charttime BETWEEN t.intime+INTERVAL 23 HOUR AND t.intime+INTERVAL 25 HOUR THEN b.lactate END, abs(date_diff('minute', b.charttime, t.intime+INTERVAL 24 HOUR))) AS lact_at24h,
  arg_min(CASE WHEN b.ph BETWEEN 6.5 AND 8.0 THEN b.ph END, b.charttime) AS ph_base_first,
  arg_max(CASE WHEN b.ph BETWEEN 6.5 AND 8.0 THEN b.ph END, b.charttime) AS ph_last_0_24,
  arg_min(CASE WHEN b.ph BETWEEN 6.5 AND 8.0 AND b.charttime BETWEEN t.intime+INTERVAL 23 HOUR AND t.intime+INTERVAL 25 HOUR THEN b.ph END, abs(date_diff('minute', b.charttime, t.intime+INTERVAL 24 HOUR))) AS ph_at24h
FROM mimiciv_derived.bg b JOIN t ON b.subject_id=t.subject_id AND b.hadm_id=t.hadm_id
WHERE b.charttime BETWEEN t.intime AND t.intime + INTERVAL 24 HOUR
GROUP BY 1", T))

dyn_gcs <- dbGetQuery(con, sprintf("
WITH t AS %s
SELECT t.stay_id,
  arg_min(CASE WHEN g.gcs BETWEEN 3 AND 15 THEN g.gcs END, g.charttime) AS gcs_base_first,
  arg_max(CASE WHEN g.gcs BETWEEN 3 AND 15 THEN g.gcs END, g.charttime) AS gcs_last_0_24,
  arg_min(CASE WHEN g.gcs BETWEEN 3 AND 15 AND g.charttime BETWEEN t.intime+INTERVAL 23 HOUR AND t.intime+INTERVAL 25 HOUR THEN g.gcs END, abs(date_diff('minute', g.charttime, t.intime+INTERVAL 24 HOUR))) AS gcs_at24h
FROM mimiciv_derived.gcs g JOIN t USING(stay_id)
WHERE g.charttime BETWEEN t.intime AND t.intime + INTERVAL 24 HOUR
GROUP BY 1", T))

dyn_sbp <- dbGetQuery(con, sprintf("
WITH t AS %s
SELECT t.stay_id,
  arg_min(CASE WHEN v.sbp BETWEEN 20 AND 300 THEN v.sbp END, v.charttime) AS sbp_base_first,
  arg_max(CASE WHEN v.sbp BETWEEN 20 AND 300 THEN v.sbp END, v.charttime) AS sbp_last_0_24,
  arg_min(CASE WHEN v.sbp BETWEEN 20 AND 300 AND v.charttime BETWEEN t.intime+INTERVAL 23 HOUR AND t.intime+INTERVAL 25 HOUR THEN v.sbp END, abs(date_diff('minute', v.charttime, t.intime+INTERVAL 24 HOUR))) AS sbp_at24h
FROM mimiciv_derived.vitalsign v JOIN t USING(stay_id)
WHERE v.charttime BETWEEN t.intime AND t.intime + INTERVAL 24 HOUR
GROUP BY 1", T))

# ================================================================ 2. ROX 组分（NIV 零点 ±1 h）
rox_spo2_rr <- dbGetQuery(con, sprintf("
WITH t AS (SELECT stay_id, niv_start FROM read_parquet('data/processed/c3/mimic_niv_analysis.parquet'))
SELECT t.stay_id,
  arg_min(CASE WHEN v.spo2 BETWEEN 0 AND 100 AND v.charttime BETWEEN t.niv_start+INTERVAL 1 HOUR AND t.niv_start+INTERVAL 3 HOUR THEN v.spo2 END, abs(date_diff('minute', v.charttime, t.niv_start+INTERVAL 2 HOUR))) AS spo2_2h,
  arg_min(CASE WHEN v.spo2 BETWEEN 0 AND 100 AND v.charttime BETWEEN t.niv_start+INTERVAL 11 HOUR AND t.niv_start+INTERVAL 13 HOUR THEN v.spo2 END, abs(date_diff('minute', v.charttime, t.niv_start+INTERVAL 12 HOUR))) AS spo2_12h,
  arg_min(CASE WHEN v.spo2 BETWEEN 0 AND 100 AND v.charttime BETWEEN t.niv_start+INTERVAL 23 HOUR AND t.niv_start+INTERVAL 25 HOUR THEN v.spo2 END, abs(date_diff('minute', v.charttime, t.niv_start+INTERVAL 24 HOUR))) AS spo2_24h,
  arg_min(CASE WHEN v.resp_rate BETWEEN 0 AND 80 AND v.charttime BETWEEN t.niv_start+INTERVAL 1 HOUR AND t.niv_start+INTERVAL 3 HOUR THEN v.resp_rate END, abs(date_diff('minute', v.charttime, t.niv_start+INTERVAL 2 HOUR))) AS rr_2h,
  arg_min(CASE WHEN v.resp_rate BETWEEN 0 AND 80 AND v.charttime BETWEEN t.niv_start+INTERVAL 11 HOUR AND t.niv_start+INTERVAL 13 HOUR THEN v.resp_rate END, abs(date_diff('minute', v.charttime, t.niv_start+INTERVAL 12 HOUR))) AS rr_12h,
  arg_min(CASE WHEN v.resp_rate BETWEEN 0 AND 80 AND v.charttime BETWEEN t.niv_start+INTERVAL 23 HOUR AND t.niv_start+INTERVAL 25 HOUR THEN v.resp_rate END, abs(date_diff('minute', v.charttime, t.niv_start+INTERVAL 24 HOUR))) AS rr_24h
FROM mimiciv_derived.vitalsign v JOIN t USING(stay_id)
WHERE v.charttime BETWEEN t.niv_start+INTERVAL 1 HOUR AND t.niv_start+INTERVAL 25 HOUR
GROUP BY 1"))

rox_fio2_chart <- dbGetQuery(con, sprintf("
WITH t AS (SELECT stay_id, niv_start FROM read_parquet('data/processed/c3/mimic_niv_analysis.parquet'))
SELECT t.stay_id,
  arg_min(CASE WHEN TRY_CAST(c.value AS DOUBLE) BETWEEN 21 AND 100 AND c.charttime BETWEEN t.niv_start+INTERVAL 1 HOUR AND t.niv_start+INTERVAL 3 HOUR THEN TRY_CAST(c.value AS DOUBLE) END, abs(date_diff('minute', c.charttime, t.niv_start+INTERVAL 2 HOUR))) AS fio2c_2h,
  arg_min(CASE WHEN TRY_CAST(c.value AS DOUBLE) BETWEEN 21 AND 100 AND c.charttime BETWEEN t.niv_start+INTERVAL 11 HOUR AND t.niv_start+INTERVAL 13 HOUR THEN TRY_CAST(c.value AS DOUBLE) END, abs(date_diff('minute', c.charttime, t.niv_start+INTERVAL 12 HOUR))) AS fio2c_12h,
  arg_min(CASE WHEN TRY_CAST(c.value AS DOUBLE) BETWEEN 21 AND 100 AND c.charttime BETWEEN t.niv_start+INTERVAL 23 HOUR AND t.niv_start+INTERVAL 25 HOUR THEN TRY_CAST(c.value AS DOUBLE) END, abs(date_diff('minute', c.charttime, t.niv_start+INTERVAL 24 HOUR))) AS fio2c_24h
FROM main.chartevents c JOIN t USING(stay_id)
WHERE c.itemid = 223835
  AND c.charttime BETWEEN t.niv_start+INTERVAL 1 HOUR AND t.niv_start+INTERVAL 25 HOUR
GROUP BY 1"))

rox_fio2_bg <- dbGetQuery(con, sprintf("
WITH t AS (SELECT stay_id, hadm_id, subject_id, niv_start FROM read_parquet('data/processed/c3/mimic_niv_analysis.parquet'))
SELECT t.stay_id,
  arg_min(CASE WHEN b.fio2 BETWEEN 21 AND 100 AND b.charttime BETWEEN t.niv_start+INTERVAL 1 HOUR AND t.niv_start+INTERVAL 3 HOUR THEN b.fio2 END, abs(date_diff('minute', b.charttime, t.niv_start+INTERVAL 2 HOUR))) AS fio2b_2h,
  arg_min(CASE WHEN b.fio2 BETWEEN 21 AND 100 AND b.charttime BETWEEN t.niv_start+INTERVAL 11 HOUR AND t.niv_start+INTERVAL 13 HOUR THEN b.fio2 END, abs(date_diff('minute', b.charttime, t.niv_start+INTERVAL 12 HOUR))) AS fio2b_12h,
  arg_min(CASE WHEN b.fio2 BETWEEN 21 AND 100 AND b.charttime BETWEEN t.niv_start+INTERVAL 23 HOUR AND t.niv_start+INTERVAL 25 HOUR THEN b.fio2 END, abs(date_diff('minute', b.charttime, t.niv_start+INTERVAL 24 HOUR))) AS fio2b_24h
FROM mimiciv_derived.bg b JOIN t ON b.subject_id=t.subject_id AND b.hadm_id=t.hadm_id
WHERE b.charttime BETWEEN t.niv_start+INTERVAL 1 HOUR AND t.niv_start+INTERVAL 25 HOUR
GROUP BY 1"))

# ================================================================ 3. HACOR 组分（NIV 零点）
hacor_lab <- dbGetQuery(con, sprintf("
WITH t AS (SELECT stay_id, hadm_id, subject_id, niv_start FROM read_parquet('data/processed/c3/mimic_niv_analysis.parquet'))
SELECT t.stay_id,
  MIN(CASE WHEN b.ph BETWEEN 6.5 AND 8.0 AND b.charttime <= t.niv_start+INTERVAL 2 HOUR THEN b.ph END) AS ph_2h,
  MIN(CASE WHEN b.ph BETWEEN 6.5 AND 8.0 THEN b.ph END) AS ph_24h,
  MIN(CASE WHEN b.po2 BETWEEN 10 AND 700 AND b.charttime <= t.niv_start+INTERVAL 2 HOUR THEN b.po2 END) AS pao2_2h,
  MIN(CASE WHEN b.po2 BETWEEN 10 AND 700 THEN b.po2 END) AS pao2_24h,
  arg_min(CASE WHEN b.fio2 BETWEEN 21 AND 100 AND b.charttime <= t.niv_start+INTERVAL 2 HOUR THEN b.fio2 END, abs(date_diff('minute', b.charttime, t.niv_start+INTERVAL 1 HOUR))) AS hacor_fio2_2h,
  arg_min(CASE WHEN b.fio2 BETWEEN 21 AND 100 THEN b.fio2 END, abs(date_diff('minute', b.charttime, t.niv_start+INTERVAL 12 HOUR))) AS hacor_fio2_24h,
  MIN(CASE WHEN b.pao2fio2ratio BETWEEN 10 AND 1000 AND b.charttime <= t.niv_start+INTERVAL 2 HOUR THEN b.pao2fio2ratio END) AS pf_2h,
  MIN(CASE WHEN b.pao2fio2ratio BETWEEN 10 AND 1000 THEN b.pao2fio2ratio END) AS pf_24h
FROM mimiciv_derived.bg b JOIN t ON b.subject_id=t.subject_id AND b.hadm_id=t.hadm_id
WHERE b.charttime BETWEEN t.niv_start AND t.niv_start+INTERVAL 24 HOUR
GROUP BY 1"))

hacor_vit <- dbGetQuery(con, sprintf("
WITH t AS (SELECT stay_id, niv_start FROM read_parquet('data/processed/c3/mimic_niv_analysis.parquet'))
SELECT t.stay_id,
  MAX(CASE WHEN v.heart_rate BETWEEN 20 AND 300 AND v.charttime <= t.niv_start+INTERVAL 2 HOUR THEN v.heart_rate END) AS hr_hacor_2h,
  MAX(CASE WHEN v.heart_rate BETWEEN 20 AND 300 THEN v.heart_rate END) AS hr_hacor_24h,
  MAX(CASE WHEN v.resp_rate BETWEEN 0 AND 80 AND v.charttime <= t.niv_start+INTERVAL 2 HOUR THEN v.resp_rate END) AS rr_hacor_2h,
  MAX(CASE WHEN v.resp_rate BETWEEN 0 AND 80 THEN v.resp_rate END) AS rr_hacor_24h
FROM mimiciv_derived.vitalsign v JOIN t USING(stay_id)
WHERE v.charttime BETWEEN t.niv_start AND t.niv_start+INTERVAL 24 HOUR
GROUP BY 1"))

hacor_gcs <- dbGetQuery(con, sprintf("
WITH t AS (SELECT stay_id, niv_start FROM read_parquet('data/processed/c3/mimic_niv_analysis.parquet'))
SELECT t.stay_id,
  MIN(CASE WHEN g.gcs BETWEEN 3 AND 15 AND g.charttime <= t.niv_start+INTERVAL 2 HOUR THEN g.gcs END) AS gcs_2h,
  MIN(CASE WHEN g.gcs BETWEEN 3 AND 15 THEN g.gcs END) AS gcs_24h
FROM mimiciv_derived.gcs g JOIN t USING(stay_id)
WHERE g.charttime BETWEEN t.niv_start AND t.niv_start+INTERVAL 24 HOUR
GROUP BY 1"))

# ================================================================ 4. 来源与服务科室
serv_last <- dbGetQuery(con, "
WITH s AS (SELECT hadm_id, curr_service, row_number() OVER (PARTITION BY hadm_id ORDER BY transfertime DESC) rn FROM main.services)
SELECT hadm_id, curr_service AS curr_service_last FROM s WHERE rn = 1")

serv_admit <- dbGetQuery(con, sprintf("
WITH t AS (SELECT o.stay_id, o.hadm_id, i.intime FROM read_parquet('data/processed/c3/mimic_niv_analysis.parquet') o JOIN main.icustays i USING(stay_id)),
s AS (SELECT t.stay_id, s.curr_service,
        row_number() OVER (PARTITION BY t.stay_id ORDER BY s.transfertime DESC) rn
      FROM t JOIN main.services s ON s.hadm_id=t.hadm_id AND s.transfertime <= t.intime)
SELECT stay_id, curr_service AS service_at_icu_admit FROM s WHERE rn = 1", T))

# ================================================================ 5. 计划/围术期插管线索
proc_icd <- dbGetQuery(con, sprintf("
WITH t AS (SELECT o.stay_id, o.hadm_id, CAST(o.niv_start AS DATE) AS niv_date FROM read_parquet('data/processed/c3/mimic_niv_analysis.parquet') o)
SELECT t.stay_id,
  MAX(1) AS proc_icd_any,
  MAX(CASE WHEN abs(date_diff('day', CAST(p.chartdate AS DATE), t.niv_date)) <= 1 THEN 1 ELSE 0 END) AS proc_icd_near_niv
FROM t JOIN main.procedures_icd p USING(hadm_id)
GROUP BY 1", T))

inv_proc <- dbGetQuery(con, sprintf("
WITH t AS (SELECT stay_id, niv_start FROM read_parquet('data/processed/c3/mimic_niv_analysis.parquet'))
SELECT t.stay_id,
  MAX(CASE WHEN p.starttime > t.niv_start AND p.starttime <= t.niv_start+INTERVAL 48 HOUR THEN 1 ELSE 0 END) AS inv_intubation_proc
FROM main.procedureevents p JOIN t USING(stay_id)
WHERE p.itemid = 224385
GROUP BY 1", T))

# ================================================================ 6. DNR / DNI
codest <- dbGetQuery(con, sprintf("
WITH t AS (SELECT stay_id, niv_start FROM read_parquet('data/processed/c3/mimic_niv_analysis.parquet'))
SELECT t.stay_id,
  MAX(CASE WHEN lower(c.value) LIKE '%%dnr%%' OR lower(c.value) LIKE '%%dnar%%' THEN 1 ELSE 0 END) AS dnr_any,
  MAX(CASE WHEN lower(c.value) LIKE '%%dni%%' THEN 1 ELSE 0 END) AS dni_any,
  MAX(CASE WHEN lower(c.value) LIKE '%%comfort%%' THEN 1 ELSE 0 END) AS comfort_any,
  MAX(CASE WHEN c.charttime <= t.niv_start AND (lower(c.value) LIKE '%%dnr%%' OR lower(c.value) LIKE '%%dnar%%') THEN 1 ELSE 0 END) AS dnr_before_niv_new,
  MAX(CASE WHEN c.charttime <= t.niv_start AND lower(c.value) LIKE '%%dni%%' THEN 1 ELSE 0 END) AS dni_before_niv,
  MAX(CASE WHEN c.charttime <= t.niv_start AND lower(c.value) LIKE '%%comfort%%' THEN 1 ELSE 0 END) AS comfort_before_niv,
  MAX(1) AS has_codestatus_record
FROM main.chartevents c JOIN t USING(stay_id)
WHERE c.itemid IN (223758, 228687)
GROUP BY 1", T))

# ================================================================ 7. 6 / 12 h 结局
inv_win <- dbGetQuery(con, sprintf("
WITH t AS (SELECT stay_id, niv_start FROM read_parquet('data/processed/c3/mimic_niv_analysis.parquet'))
SELECT t.stay_id, MIN(p.starttime) AS t1
FROM main.procedureevents p JOIN t USING(stay_id)
WHERE p.itemid = 225792 AND p.starttime > t.niv_start
GROUP BY 1"))

# ================================================================ 8. 队列首次入住
firststay <- dbGetQuery(con, sprintf("
WITH t AS (SELECT o.stay_id, o.subject_id FROM read_parquet('data/processed/c3/mimic_niv_analysis.parquet') o),
a AS (SELECT subject_id, stay_id,
        MIN(intime) OVER (PARTITION BY subject_id) AS pt_first_intime,
        row_number() OVER (PARTITION BY subject_id ORDER BY intime, stay_id) AS patient_unit_seq
      FROM main.icustays)
SELECT t.stay_id, a.pt_first_intime, a.patient_unit_seq,
       CASE WHEN a.patient_unit_seq = 1 THEN 1 ELSE 0 END AS first_icu_stay_new
FROM t JOIN a ON a.stay_id=t.stay_id"))

# ================================================================ 9. 合并
for (nm in c("dyn_bg","dyn_gcs","dyn_sbp","rox_spo2_rr","rox_fio2_chart","rox_fio2_bg",
             "hacor_lab","hacor_vit","hacor_gcs","serv_admit","proc_icd",
             "inv_proc","codest","inv_win","firststay")) {
  assign(nm, mc_checkpoint(get(nm), paste0("并入：", nm), id = "stay_id", side = TRUE))
}
serv_last <- mc_checkpoint(serv_last, "并入：住院末次服务科室", id = "hadm_id", side = TRUE)

out <- base
for (blk in list(dyn_bg, dyn_gcs, dyn_sbp, rox_spo2_rr, rox_fio2_chart, rox_fio2_bg,
                 hacor_lab, hacor_vit, hacor_gcs, serv_admit, proc_icd, inv_proc,
                 codest, inv_win, firststay)) {
  out <- merge(out, blk, by = "stay_id", all.x = TRUE, sort = FALSE)
}
out <- merge(out, serv_last, by = "hadm_id", all.x = TRUE, sort = FALSE)
stopifnot(!anyDuplicated(out$stay_id))

# --- 动态变量派生 ---
out$lact_delta_at24h    <- out$lact_at24h - out$lact_base_first
out$lact_delta_lastfirst<- out$lact_last_0_24 - out$lact_base_first
out$ph_delta_at24h      <- out$ph_at24h - out$ph_base_first
out$ph_delta_lastfirst  <- out$ph_last_0_24 - out$ph_base_first
out$gcs_delta_at24h     <- out$gcs_at24h - out$gcs_base_first
out$gcs_delta_lastfirst <- out$gcs_last_0_24 - out$gcs_base_first
out$sbp_delta_at24h     <- out$sbp_at24h - out$sbp_base_first
out$sbp_delta_lastfirst <- out$sbp_last_0_24 - out$sbp_base_first
for (v in c("lact","ph","gcs","sbp")) {
  out[[paste0("m_", v, "_base_first")]] <- as.integer(is.na(out[[paste0(v, "_base_first")]]))
  out[[paste0("m_", v, "_at24h")]]      <- as.integer(is.na(out[[paste0(v, "_at24h")]]))
}

# --- ROX ---
for (t in c("2h","12h","24h")) {
  c0 <- out[[paste0("fio2c_", t)]]; b0 <- out[[paste0("fio2b_", t)]]
  out[[paste0("fio2_", t)]]    <- ifelse(!is.na(c0), c0, b0)
  out[[paste0("fio2src_", t)]] <- ifelse(!is.na(c0), "chart", ifelse(!is.na(b0), "lab", NA_character_))
  out[[paste0("rox_", t)]] <- ifelse(!is.na(out[[paste0("spo2_", t)]]) &
                                       !is.na(out[[paste0("fio2_", t)]]) &
                                       !is.na(out[[paste0("rr_", t)]]) &
                                       out[[paste0("rr_", t)]] > 0,
                                     (out[[paste0("spo2_", t)]] / (out[[paste0("fio2_", t)]]/100)) / out[[paste0("rr_", t)]],
                                     NA_real_)
  out[[paste0("rox_computable_", t)]] <- as.integer(!is.na(out[[paste0("rox_", t)]]))
  for (cc in c("spo2","fio2","rr","rox")) {
    out[[paste0("m_", cc, "_", t)]] <- as.integer(is.na(out[[paste0(cc, "_", t)]]))
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

# --- 来源与病房标志 ---
out$icu_medical  <- as.integer(!is.na(out$first_careunit) & out$first_careunit %in% c("Medical Intensive Care Unit (MICU)","Coronary Care Unit (CCU)","Cardiac Vascular Intensive Care Unit (CVICU)"))
out$icu_surgical <- as.integer(!is.na(out$first_careunit) & out$first_careunit %in% c("Surgical Intensive Care Unit (SICU)","Trauma SICU (TSICU)","Neuro Surgical Intensive Care Unit (Neuro SICU)"))
out$icu_medsurg  <- as.integer(!is.na(out$first_careunit) & out$first_careunit == "Medical/Surgical Intensive Care Unit (MICU/SICU)")
out$icu_neuro    <- as.integer(!is.na(out$first_careunit) & out$first_careunit %in% c("Neuro Intermediate","Neuro Stepdown"))
out$icu_other    <- as.integer(!is.na(out$first_careunit) & out$first_careunit %in% c("PACU","Intensive Care Unit (ICU)","Medicine","Surgery/Trauma","Neurology","Medicine/Cardiology Intermediate","Med/Surg","Surgery/Vascular/Intermediate"))
out$flag_or_pacu_source <- as.integer((!is.na(out$admission_type) & out$admission_type == "SURGICAL SAME DAY ADMISSION") |
                                        (!is.na(out$first_careunit) & out$first_careunit == "PACU"))
out$src_ed <- as.integer(!is.na(out$admission_type) & out$admission_type %in% c("EW EMER.","DIRECT EMER."))
out$planned_surg_admtype <- as.integer(!is.na(out$admission_type) & out$admission_type %in% c("SURGICAL SAME DAY ADMISSION","ELECTIVE"))
out$surgical_service_at_admit <- as.integer(!is.na(out$service_at_icu_admit) & out$service_at_icu_admit %in% c("SURG","CSURG","NSURG","VSURG","TSURG","PSURG","ORTHO","GYN","GU","TRAUM","ENT"))

# --- 计划/围术期插管线索 ---
for (v in c("proc_icd_any","proc_icd_near_niv","inv_intubation_proc")) {
  out[[v]] <- ifelse(is.na(out[[v]]), 0L, out[[v]])
}
out$inv_planned_context <- as.integer(out$to_invasive_48h == 1 &
  (out$planned_surg_admtype == 1 | out$surgical_service_at_admit == 1 | out$proc_icd_near_niv == 1))

# --- DNR / DNI ---
out$dnr_before_niv_strict <- ifelse(is.na(out$dnr_before_niv_new), 0L, out$dnr_before_niv_new)
for (v in c("dnr_any","dni_any","comfort_any","dni_before_niv","comfort_before_niv","has_codestatus_record")) {
  out[[v]] <- ifelse(is.na(out[[v]]), 0L, out[[v]])
}
out$flag_dnr_before_niv <- as.integer(out$dnr_before_niv_strict == 1)

# --- 6 / 12 h 结局 ---
out$to_invasive_6h <- as.integer(!is.na(out$t1) & as.numeric(difftime(out$t1, out$niv_start, units="secs")) > 0 &
                                   as.numeric(difftime(out$t1, out$niv_start, units="secs")) <= 6*3600 &
                                   (is.na(out$death_time) | out$t1 < out$death_time))
out$to_invasive_12h <- as.integer(!is.na(out$t1) & as.numeric(difftime(out$t1, out$niv_start, units="secs")) > 0 &
                                    as.numeric(difftime(out$t1, out$niv_start, units="secs")) <= 12*3600 &
                                    (is.na(out$death_time) | out$t1 < out$death_time))
out$composite_6h  <- as.integer(out$to_invasive_6h == 1 | out$death_hosp == 1)
out$composite_12h <- as.integer(out$to_invasive_12h == 1 | out$death_hosp == 1)

# --- 队列标志 ---
out$cohort_first_icu <- out$first_icu_stay_new
out$cohort_patient_unique <- as.integer(out$patient_unit_seq == 1)
out$cohort_medical_icu <- out$icu_medical
out$cohort_surgical_icu <- out$icu_surgical
out$strict_invasive_def <- 1L
out$cohort_excl_or_pacu <- as.integer(out$in_analysis == 1 & out$flag_or_pacu_source == 0)
out$cohort_excl_dnr     <- as.integer(out$in_analysis == 1 & out$flag_dnr_before_niv == 0)
out$cohort_strict       <- as.integer(out$in_analysis == 1 & out$flag_or_pacu_source == 0 &
                                        out$flag_dnr_before_niv == 0 & out$cohort_first_icu == 1)

# ================================================================ 10. 选列、检查点、输出
keep_old <- names(base)
new_cols <- c(
  "lact_base_first","lact_last_0_24","lact_at24h","lact_delta_at24h","lact_delta_lastfirst",
  "ph_base_first","ph_last_0_24","ph_at24h","ph_delta_at24h","ph_delta_lastfirst",
  "gcs_base_first","gcs_last_0_24","gcs_at24h","gcs_delta_at24h","gcs_delta_lastfirst",
  "sbp_base_first","sbp_last_0_24","sbp_at24h","sbp_delta_at24h","sbp_delta_lastfirst",
  "m_lact_base_first","m_lact_at24h","m_ph_base_first","m_ph_at24h",
  "m_gcs_base_first","m_gcs_at24h","m_sbp_base_first","m_sbp_at24h",
  paste0(rep(c("spo2_","fio2_","fio2src_","rr_","rox_","m_spo2_","m_fio2_","m_rr_","m_rox_","rox_computable_"), each=3), c("2h","12h","24h")),
  "hacor2h_hr","hacor2h_ph","hacor2h_pao2","hacor2h_fio2","hacor2h_pf","hacor2h_gcs","hacor2h_rr","hacor2h_total","hacor2h_computable",
  "hacor24h_hr","hacor24h_ph","hacor24h_pao2","hacor24h_fio2","hacor24h_pf","hacor24h_gcs","hacor24h_rr","hacor24h_total","hacor24h_computable",
  "curr_service_last","service_at_icu_admit",
  "icu_medical","icu_surgical","icu_medsurg","icu_neuro","icu_other",
  "flag_or_pacu_source","src_ed","planned_surg_admtype","surgical_service_at_admit",
  "proc_icd_any","proc_icd_near_niv","inv_intubation_proc","inv_planned_context",
  "dnr_before_niv_strict","dni_before_niv","comfort_before_niv","dnr_any","dni_any","comfort_any","has_codestatus_record",
  "to_invasive_6h","to_invasive_12h","composite_6h","composite_12h",
  "cohort_first_icu","cohort_patient_unique","cohort_medical_icu","cohort_surgical_icu",
  "flag_dnr_before_niv","strict_invasive_def",
  "cohort_excl_or_pacu","cohort_excl_dnr","cohort_strict")
new_cols <- unique(new_cols[new_cols %in% names(out)])
final <- out[, c(keep_old, new_cols)]
final <- mc_checkpoint(final, paste0("输出：C3 MIMIC 扩取表全表（保留 ", nrow(final), " 行；in_analysis = 1 的外部分析集为 ", sum(final$in_analysis), " 行）"),
                       id = "stay_id")

con_out <- dbConnect(duckdb::duckdb())
dbWriteTable(con_out, "c3_mimic_expanded", final, overwrite = TRUE)
pq <- file.path(OUT, "mimic_niv_expanded.parquet")
dbExecute(con_out, sprintf("COPY c3_mimic_expanded TO '%s' (FORMAT PARQUET)", pq))
dbDisconnect(con_out, shutdown = TRUE)

# ================================================================ 11. 汇报数字
sub <- final[final$in_analysis == 1, ]
cat("\n=== MIMIC expanded ===\n")
cat("rows:", nrow(final), " unique:", length(unique(final$stay_id)), " main:", nrow(sub), "\n")
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
cat("source: OR/PACU", sum(sub$flag_or_pacu_source), "ED", sum(sub$src_ed),
    "medical", sum(sub$icu_medical), "surgical", sum(sub$icu_surgical),
    "med-surg", sum(sub$icu_medsurg), "neuro", sum(sub$icu_neuro), "\n")
cat("service: surgical at admit", sum(sub$surgical_service_at_admit), "curr last missing", sum(is.na(sub$curr_service_last)),
    "admit service missing", sum(is.na(sub$service_at_icu_admit)), "\n")
cat("planned: admtype", sum(sub$planned_surg_admtype), "surg service", sum(sub$surgical_service_at_admit),
    "proc any", sum(sub$proc_icd_any), "proc near", sum(sub$proc_icd_near_niv),
    "intubation proc", sum(sub$inv_intubation_proc), "planned context", sum(sub$inv_planned_context), "\n")
cat("DNR/DNI: dnr_before(merged old)", sum(sub$dnr_before_niv, na.rm=TRUE), "dnr_before(strict)", sum(sub$dnr_before_niv_strict), "dni_before", sum(sub$dni_before_niv),
    "comfort_before", sum(sub$comfort_before_niv), "dnr_any", sum(sub$dnr_any),
    "dni_any", sum(sub$dni_any), "no code status", sum(sub$has_codestatus_record==0), "\n")
cat("outcomes: 6h", sum(sub$to_invasive_6h, na.rm=TRUE), "12h", sum(sub$to_invasive_12h, na.rm=TRUE),
    "24h", sum(sub$to_invasive_24h, na.rm=TRUE), "48h", sum(sub$to_invasive_48h, na.rm=TRUE),
    "composite6h", sum(sub$composite_6h, na.rm=TRUE), "composite12h", sum(sub$composite_12h, na.rm=TRUE), "\n")
cat("cohort flags: first_icu", sum(sub$cohort_first_icu), "patient_unique", sum(sub$cohort_patient_unique),
    "excl_or_pacu", sum(sub$cohort_excl_or_pacu), "excl_dnr", sum(sub$cohort_excl_dnr),
    "strict", sum(sub$cohort_strict), "\n")
cat("Wrote:", pq, "\n")
