#!/usr/bin/env Rscript
# C3 外部验证队列：MIMIC-IV v3.1 无创通气（NIV）ICU 入住正式取数
# 一行 = 一次 ICU 入住（stay_id）。保留全部有无创操作的入住，
# 用 in_analysis 标记外部分析集（NIV-first、非同刻、年龄 >= 18）。
# 只读源数据；不插补、不建模、不做统计检验。
suppressPackageStartupMessages({ library(DBI); library(duckdb) })
source((function(f) { for (p in c(file.path("code/00_common", f), file.path("../00_common", f), file.path("00_common", f), f)) if (file.exists(p)) return(p); f })("trace_logger.R"))

DB  <- "F:/CoredataCenter/MIMIC-IV/mimic.duckdb"
OUT <- "data/processed/c3"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
con <- dbConnect(duckdb::duckdb(), dbdir = DB, read_only = TRUE)
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

# ---------------------------------------------------------------- 1. 队列与有创先后
cohort <- dbGetQuery(con, "
WITH niv AS (
  SELECT stay_id, MIN(starttime) AS t0 FROM procedureevents
  WHERE itemid = 225794 GROUP BY 1
),
inv AS (
  SELECT p.stay_id,
    MIN(CASE WHEN p.starttime <  n.t0 THEN p.starttime END) AS inv_before_min,
    MIN(CASE WHEN p.starttime =  n.t0 THEN p.starttime END) AS inv_same_min,
    MIN(CASE WHEN p.starttime >  n.t0 THEN p.starttime END) AS t1_strict_any,
    MIN(CASE WHEN p.starttime >  n.t0 AND p.starttime <= n.t0 + INTERVAL 48 HOUR
             THEN p.starttime END) AS t1_strict_48h,
    MIN(CASE WHEN p.starttime >  n.t0 AND p.starttime <= n.t0 + INTERVAL 24 HOUR
             THEN p.starttime END) AS t1_strict_24h,
    SUM(CASE WHEN p.starttime = n.t0 THEN 1 ELSE 0 END) AS n_same
  FROM procedureevents p JOIN niv n USING (stay_id)
  WHERE p.itemid = 225792
  GROUP BY 1
)
SELECT i.stay_id, i.hadm_id, i.subject_id, i.first_careunit, i.intime, i.outtime,
       i.los AS icu_los_days,
       a.admittime, a.dischtime, a.deathtime, a.hospital_expire_flag,
       a.admission_type, a.race,
       p.gender, p.anchor_age, p.anchor_year,
       n.t0 AS niv_start,
       i2.inv_before_min, i2.inv_same_min, i2.t1_strict_any, i2.t1_strict_48h,
       i2.t1_strict_24h, i2.n_same
FROM icustays i
JOIN niv n ON n.stay_id = i.stay_id
JOIN admissions a ON a.hadm_id = i.hadm_id
JOIN patients p ON p.subject_id = i.subject_id
LEFT JOIN inv i2 ON i2.stay_id = i.stay_id
")
cohort <- mc_checkpoint(cohort, "纳入：有无创通气（NIV）操作的 ICU 入住", id = "stay_id")

cohort$age <- cohort$anchor_age + as.integer(format(as.Date(cohort$admittime), "%Y")) - cohort$anchor_year
cohort$age_lt18 <- as.integer(!is.na(cohort$age) & cohort$age < 18)
cohort$inv_before_niv <- as.integer(!is.na(cohort$inv_before_min))
cohort$inv_same_time  <- as.integer(!is.na(cohort$inv_same_min) | (!is.na(cohort$n_same) & cohort$n_same > 0))
cohort$niv_first      <- 1L - cohort$inv_before_niv
cohort$niv_record     <- 1L

cohort <- mc_checkpoint(cohort, "排除有创早于无创（序贯撤机）", id = "stay_id",
                        keep = cohort$inv_before_niv == 0)
cohort <- mc_checkpoint(cohort, "排除有创与无创开始同刻", id = "stay_id",
                        keep = cohort$inv_before_niv == 0 & cohort$inv_same_time == 0)
cohort$in_analysis <- as.integer(cohort$inv_before_niv == 0 & cohort$inv_same_time == 0 &
                                   !is.na(cohort$age) & cohort$age >= 18)
cohort <- mc_checkpoint(cohort, "纳入：NIV-first、非同刻、年龄 >= 18 岁", id = "stay_id",
                        keep = cohort$in_analysis == 1)

# ---------------------------------------------------------------- 2. 首日严重度、生命体征、化验
sofa <- dbGetQuery(con, "SELECT stay_id, sofa FROM mimiciv_derived.first_day_sofa")
gcs  <- dbGetQuery(con, "SELECT stay_id, gcs_min FROM mimiciv_derived.first_day_gcs")
vit  <- dbGetQuery(con, "SELECT stay_id, heart_rate_min AS hr_min, heart_rate_max AS hr_max,
  sbp_min, resp_rate_max AS rr_max, spo2_min FROM mimiciv_derived.first_day_vitalsign")
lab  <- dbGetQuery(con, "SELECT stay_id, bicarbonate_min FROM mimiciv_derived.first_day_lab")

# 血气：-6 至 24 小时（ICU 入科为零点）；同一住院内关联
bg24 <- dbGetQuery(con, "
SELECT i.stay_id,
  MIN(CASE WHEN b.ph BETWEEN 6.5 AND 8.0 THEN b.ph END) AS ph_min,
  MAX(CASE WHEN b.pco2 BETWEEN 5 AND 200 THEN b.pco2 END) AS paco2_max,
  MAX(CASE WHEN b.lactate BETWEEN 0 AND 30 THEN b.lactate END) AS lactate_max,
  MIN(CASE WHEN b.pao2fio2ratio BETWEEN 10 AND 1000 THEN b.pao2fio2ratio END) AS pf_min
FROM icustays i
LEFT JOIN mimiciv_derived.bg b
  ON b.subject_id = i.subject_id AND b.hadm_id = i.hadm_id
 AND b.charttime BETWEEN i.intime - INTERVAL 6 HOUR AND i.intime + INTERVAL 24 HOUR
GROUP BY 1")

# ---------------------------------------------------------------- 3. HACOR 五项（NIV 零点）
hacor_vit <- dbGetQuery(con, "
WITH niv AS (SELECT stay_id, MIN(starttime) AS t0 FROM procedureevents
  WHERE itemid = 225794 GROUP BY 1)
SELECT n.stay_id,
  MAX(CASE WHEN vs.heart_rate BETWEEN 20 AND 300
            AND vs.charttime BETWEEN n.t0 AND n.t0 + INTERVAL 2 HOUR THEN vs.heart_rate END) AS hr_2h,
  MAX(CASE WHEN vs.heart_rate BETWEEN 20 AND 300
            AND vs.charttime BETWEEN n.t0 AND n.t0 + INTERVAL 24 HOUR THEN vs.heart_rate END) AS hr_24h,
  MAX(CASE WHEN vs.resp_rate BETWEEN 0 AND 80
            AND vs.charttime BETWEEN n.t0 AND n.t0 + INTERVAL 2 HOUR THEN vs.resp_rate END) AS rr_2h,
  MAX(CASE WHEN vs.resp_rate BETWEEN 0 AND 80
            AND vs.charttime BETWEEN n.t0 AND n.t0 + INTERVAL 24 HOUR THEN vs.resp_rate END) AS rr_24h
FROM niv n LEFT JOIN mimiciv_derived.vitalsign vs USING (stay_id)
GROUP BY 1")

hacor_gcs <- dbGetQuery(con, "
WITH niv AS (SELECT stay_id, MIN(starttime) AS t0 FROM procedureevents
  WHERE itemid = 225794 GROUP BY 1)
SELECT n.stay_id,
  MIN(CASE WHEN g.charttime BETWEEN n.t0 AND n.t0 + INTERVAL 2 HOUR THEN g.gcs END) AS gcs_hacor_2h,
  MIN(CASE WHEN g.charttime BETWEEN n.t0 AND n.t0 + INTERVAL 24 HOUR THEN g.gcs END) AS gcs_hacor_24h
FROM niv n LEFT JOIN mimiciv_derived.gcs g USING (stay_id)
GROUP BY 1")

hacor_bg <- dbGetQuery(con, "
WITH niv AS (SELECT n.stay_id, n.hadm_id, n.subject_id, MIN(p.starttime) AS t0
  FROM procedureevents p JOIN icustays n USING (stay_id)
  WHERE p.itemid = 225794 GROUP BY 1,2,3)
SELECT n.stay_id,
  MIN(CASE WHEN b.ph BETWEEN 6.5 AND 8.0
            AND b.charttime BETWEEN n.t0 AND n.t0 + INTERVAL 2 HOUR THEN b.ph END) AS ph_hacor_2h,
  MIN(CASE WHEN b.ph BETWEEN 6.5 AND 8.0
            AND b.charttime BETWEEN n.t0 AND n.t0 + INTERVAL 24 HOUR THEN b.ph END) AS ph_hacor_24h,
  MIN(CASE WHEN b.pao2fio2ratio BETWEEN 10 AND 1000
            AND b.charttime BETWEEN n.t0 AND n.t0 + INTERVAL 2 HOUR THEN b.pao2fio2ratio END) AS pf_hacor_2h,
  MIN(CASE WHEN b.pao2fio2ratio BETWEEN 10 AND 1000
            AND b.charttime BETWEEN n.t0 AND n.t0 + INTERVAL 24 HOUR THEN b.pao2fio2ratio END) AS pf_hacor_24h
FROM niv n LEFT JOIN mimiciv_derived.bg b
  ON b.subject_id = n.subject_id AND b.hadm_id = n.hadm_id
GROUP BY 1")

# ---------------------------------------------------------------- 4. 主诊断与基线 code status
dx <- dbGetQuery(con, "
WITH niv AS (SELECT stay_id, hadm_id FROM procedureevents WHERE itemid = 225794 GROUP BY 1,2)
SELECT DISTINCT n.hadm_id, d.icd_code, d.icd_version
FROM niv n JOIN diagnoses_icd d USING (hadm_id)
WHERE d.seq_num = 1
")
# 同一住院多条 seq_num=1 时只保留一条（按编码排序，规则固定）
dx <- dx[order(dx$hadm_id, dx$icd_version, dx$icd_code), ]
dx <- dx[!duplicated(dx$hadm_id), ]

cs <- dbGetQuery(con, "
WITH niv AS (SELECT stay_id, MIN(starttime) AS t0 FROM procedureevents
  WHERE itemid = 225794 GROUP BY 1)
SELECT n.stay_id,
  MAX(CASE WHEN c.charttime <= n.t0
            AND (lower(c.value) LIKE '%dnr%' OR lower(c.value) LIKE '%dni%'
                 OR lower(c.value) LIKE '%comfort%') THEN 1 ELSE 0 END) AS dnr_before_niv
FROM niv n LEFT JOIN chartevents c
  ON c.stay_id = n.stay_id AND c.itemid IN (223758, 228687, 229784)
GROUP BY 1")

# ---------------------------------------------------------------- 5. 合并
sofa <- mc_checkpoint(sofa, "并入：首日 SOFA", id = "stay_id", side = TRUE)
gcs  <- mc_checkpoint(gcs, "并入：首日 GCS", id = "stay_id", side = TRUE)
vit  <- mc_checkpoint(vit, "并入：首日生命体征", id = "stay_id", side = TRUE)
lab  <- mc_checkpoint(lab, "并入：首日化验", id = "stay_id", side = TRUE)
bg24 <- mc_checkpoint(bg24, "并入：血气 -6–24h", id = "stay_id", side = TRUE)
hacor_vit <- mc_checkpoint(hacor_vit, "并入：HACOR 心率与呼吸频率", id = "stay_id", side = TRUE)
hacor_gcs <- mc_checkpoint(hacor_gcs, "并入：HACOR GCS", id = "stay_id", side = TRUE)
hacor_bg  <- mc_checkpoint(hacor_bg, "并入：HACOR pH 与 P/F", id = "stay_id", side = TRUE)
cs <- mc_checkpoint(cs, "并入：无创前 code status", id = "stay_id", side = TRUE)
dx <- mc_checkpoint(dx, "并入：本次住院主诊断", id = "hadm_id", side = TRUE)

out <- cohort
for (blk in list(sofa, gcs, vit, lab, bg24, hacor_vit, hacor_gcs, hacor_bg, cs)) {
  out <- merge(out, blk, by = "stay_id", all.x = TRUE, sort = FALSE)
}
out <- merge(out, dx, by = "hadm_id", all.x = TRUE, sort = FALSE)
stopifnot(!anyDuplicated(out$stay_id))

# ---------------------------------------------------------------- 6. 结局与派生
out$death_hosp <- as.integer(out$hospital_expire_flag == 1)
out$death_time <- out$deathtime

inv_event <- function(t1, td) as.integer(!is.na(t1) & (is.na(td) | t1 < td))
out$to_invasive_48h <- inv_event(out$t1_strict_48h, out$death_time)
out$to_invasive_24h <- inv_event(out$t1_strict_24h, out$death_time)
out$to_invasive_any <- inv_event(out$t1_strict_any, out$death_time)
death_secs <- as.numeric(difftime(out$death_time, out$niv_start, units = "secs"))
out$death_48h <- as.integer(out$death_hosp == 1 & !is.na(death_secs) &
                              death_secs > 0 & death_secs <= 48 * 3600)
out$composite_48h <- as.integer(out$to_invasive_48h == 1 | out$death_hosp == 1)
out$composite_48h_strictdeath <- as.integer(out$to_invasive_48h == 1 | out$death_48h == 1)
out$composite_any <- as.integer(out$to_invasive_any == 1 | out$death_hosp == 1)
out$followup_48h_complete <- ifelse(is.na(out$outtime),
  NA_integer_,
  as.integer(as.numeric(difftime(out$outtime, out$niv_start, units = "hours")) >= 48 |
             out$to_invasive_48h == 1 | out$death_48h == 1))
out$hosp_los_days <- as.numeric(difftime(out$dischtime, out$admittime, units = "days"))

out$m_sofa <- as.integer(is.na(out$sofa))
out$m_gcs  <- as.integer(is.na(out$gcs_min))
out$m_ph   <- as.integer(is.na(out$ph_min))
out$m_paco2 <- as.integer(is.na(out$paco2_max))
out$m_pf   <- as.integer(is.na(out$pf_min))
out$m_lactate <- as.integer(is.na(out$lactate_max))
out$m_sbp  <- as.integer(is.na(out$sbp_min))
out$m_rr   <- as.integer(is.na(out$rr_max))
out$m_spo2 <- as.integer(is.na(out$spo2_min))

# 入院诊断 6 组：按本次住院主诊断 ICD-9/10 有序映射
map_icd <- function(code, ver) {
  code <- ifelse(is.na(code), "", toupper(gsub("\\.", "", code)))
  grp <- rep("术后/其他", length(code))
  is9 <- ver == 9
  # 脓毒症
  grp[is9 & grepl("^038", code)] <- "脓毒症"
  grp[is9 & code %in% c("99591","99592")] <- "脓毒症"
  grp[!is9 & grepl("^A4[01]", code)] <- "脓毒症"
  grp[!is9 & grepl("^R652", code)] <- "脓毒症"
  # COPD / 哮喘
  grp[grp == "术后/其他" & is9 & grepl("^49[1-3]|^496", code)] <- "COPD/哮喘"
  grp[grp == "术后/其他" & !is9 & grepl("^J4[1-5]", code)] <- "COPD/哮喘"
  # 心源性肺水肿 / 心衰
  grp[grp == "术后/其他" & is9 & grepl("^428|^4020[13]|^4021[13]|^4029[13]|^404", code)] <- "心源性肺水肿"
  grp[grp == "术后/其他" & !is9 & grepl("^I50|^I110|^I13[02]", code)] <- "心源性肺水肿"
  # 肺炎 / 其他呼吸
  grp[grp == "术后/其他" & is9 & grepl("^48[0-6]|^507|^51[12]|^518|^4151", code)] <- "肺炎/其他呼吸"
  grp[grp == "术后/其他" & !is9 & grepl("^J1[2-8]|^J69|^J8[0-4]|^J9[0-6]|^I26", code)] <- "肺炎/其他呼吸"
  # 神经系统
  grp[grp == "术后/其他" & is9 & grepl("^43[0-8]|^345|^78001|^7803", code)] <- "神经系统"
  grp[grp == "术后/其他" & !is9 & grepl("^I6[0-9]|^G40|^R402|^R56", code)] <- "神经系统"
  grp
}
out$dx_group <- factor(map_icd(out$icd_code, out$icd_version),
                       levels = c("COPD/哮喘", "心源性肺水肿", "肺炎/其他呼吸", "脓毒症", "神经系统", "术后/其他"))

# HACOR 计分
score_hr  <- function(x) ifelse(is.na(x), NA_integer_, ifelse(x < 90, 0L, ifelse(x < 120, 1L, 2L)))
score_ph  <- function(x) ifelse(is.na(x), NA_integer_, ifelse(x >= 7.35, 0L, ifelse(x >= 7.30, 1L, ifelse(x >= 7.25, 2L, ifelse(x >= 7.20, 3L, 4L)))))
score_gcs <- function(x) ifelse(is.na(x), NA_integer_, ifelse(x >= 15, 0L, ifelse(x >= 13, 1L, ifelse(x >= 11, 2L, ifelse(x >= 9, 3L, 4L)))))
score_pf  <- function(x) ifelse(is.na(x), NA_integer_, ifelse(x >= 201, 0L, ifelse(x >= 176, 1L, ifelse(x >= 151, 2L, ifelse(x >= 126, 3L, ifelse(x >= 101, 4L, 5L))))))
score_rr  <- function(x) ifelse(is.na(x), NA_integer_, ifelse(x <= 30, 0L, ifelse(x <= 35, 1L, ifelse(x <= 40, 2L, ifelse(x <= 45, 3L, 4L)))))
mk_hacor <- function(hr, ph, gcs, pf, rr) {
  s <- cbind(score_hr(hr), score_ph(ph), score_gcs(gcs), score_pf(pf), score_rr(rr))
  ifelse(rowSums(is.na(s)) == 0, rowSums(s), NA_integer_)
}
out$hacor_hr_2h <- score_hr(out$hr_2h); out$hacor_ph_2h <- score_ph(out$ph_hacor_2h)
out$hacor_gcs_2h <- score_gcs(out$gcs_hacor_2h); out$hacor_pf_2h <- score_pf(out$pf_hacor_2h)
out$hacor_rr_2h <- score_rr(out$rr_2h)
out$hacor_2h <- mk_hacor(out$hr_2h, out$ph_hacor_2h, out$gcs_hacor_2h, out$pf_hacor_2h, out$rr_2h)
out$hacor_hr_24h <- score_hr(out$hr_24h); out$hacor_ph_24h <- score_ph(out$ph_hacor_24h)
out$hacor_gcs_24h <- score_gcs(out$gcs_hacor_24h); out$hacor_pf_24h <- score_pf(out$pf_hacor_24h)
out$hacor_rr_24h <- score_rr(out$rr_24h)
out$hacor_24h <- mk_hacor(out$hr_24h, out$ph_hacor_24h, out$gcs_hacor_24h, out$pf_hacor_24h, out$rr_24h)
out$hacor_score <- ifelse(!is.na(out$hacor_2h), out$hacor_2h, out$hacor_24h)
out$hacor_window <- ifelse(!is.na(out$hacor_2h), "2h", ifelse(!is.na(out$hacor_24h), "24h", NA_character_))
out$hacor_hr  <- ifelse(!is.na(out$hacor_2h), out$hacor_hr_2h, out$hacor_hr_24h)
out$hacor_ph  <- ifelse(!is.na(out$hacor_2h), out$hacor_ph_2h, out$hacor_ph_24h)
out$hacor_gcs <- ifelse(!is.na(out$hacor_2h), out$hacor_gcs_2h, out$hacor_gcs_24h)
out$hacor_pf  <- ifelse(!is.na(out$hacor_2h), out$hacor_pf_2h, out$hacor_pf_24h)
out$hacor_rr  <- ifelse(!is.na(out$hacor_2h), out$hacor_rr_2h, out$hacor_rr_24h)

# ---------------------------------------------------------------- 7. 选列、检查点、输出
keep <- c(
  "stay_id","hadm_id","subject_id","first_careunit","admission_type",
  "niv_record","niv_start","inv_before_niv","inv_same_time","niv_first",
  "age","age_lt18","in_analysis",
  "to_invasive_48h","to_invasive_24h","to_invasive_any",
  "death_48h","death_hosp","death_time","composite_48h","composite_48h_strictdeath","composite_any",
  "followup_48h_complete",
  "gender","race","dx_group","sofa","gcs_min","sbp_min","rr_max","spo2_min",
  "ph_min","paco2_max","pf_min","lactate_max","bicarbonate_min","hr_min","hr_max",
  "m_sofa","m_gcs","m_ph","m_paco2","m_pf","m_lactate","m_sbp","m_rr","m_spo2",
  "icu_los_days","hosp_los_days",
  "hacor_score","hacor_window","hacor_hr","hacor_ph","hacor_gcs","hacor_pf","hacor_rr",
  "hacor_2h","hacor_24h","dnr_before_niv")

out <- out[, keep]
out <- mc_checkpoint(out, "纳入：C3 MIMIC 外部分析集（NIV-first、非同刻、成人）", id = "stay_id",
                     keep = out$in_analysis == 1)

con_out <- dbConnect(duckdb::duckdb())
dbWriteTable(con_out, "c3_mimic_analysis", out, overwrite = TRUE)
pq <- file.path(OUT, "mimic_niv_analysis.parquet")
dbExecute(con_out, sprintf("COPY c3_mimic_analysis TO '%s' (FORMAT PARQUET)", pq))
dbDisconnect(con_out, shutdown = TRUE)

cat("MIMIC rows:", nrow(out), " unique stay:", length(unique(out$stay_id)), "\n")
cat("in_analysis:", sum(out$in_analysis), " NIV-first:", sum(out$niv_first),
    " inv_before:", sum(out$inv_before_niv), " same:", sum(out$inv_same_time),
    " age<18:", sum(out$age_lt18), "\n")
sub <- out[out$in_analysis == 1, ]
cat("MAIN: n=", nrow(sub), " to_invasive_48h=", sum(sub$to_invasive_48h),
    " 24h=", sum(sub$to_invasive_24h), " any=", sum(sub$to_invasive_any),
    " death_48h=", sum(sub$death_48h), " death_hosp=", sum(sub$death_hosp), "\n")
cat("MAIN composite_48h=", sum(sub$composite_48h), " strictdeath=", sum(sub$composite_48h_strictdeath),
    " any=", sum(sub$composite_any), " followup_incomplete=", sum(sub$followup_48h_complete == 0, na.rm = TRUE), "\n")
cat("MAIN SOFA missing=", sum(sub$m_sofa), " GCS missing=", sum(sub$m_gcs),
    " pH missing=", sum(sub$m_ph), " PF missing=", sum(sub$m_pf), " lactate missing=", sum(sub$m_lactate), "\n")
cat("MAIN HACOR computable=", sum(!is.na(sub$hacor_score)),
    " (2h=", sum(sub$hacor_window == "2h", na.rm = TRUE),
    ", 24h=", sum(sub$hacor_window == "24h", na.rm = TRUE), ")\n")
cat("MAIN DNR/DNI before NIV=", sum(sub$dnr_before_niv, na.rm = TRUE), "\n")
cat("Wrote:", pq, "\n")
