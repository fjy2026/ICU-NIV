#!/usr/bin/env Rscript
# 生成 C3 两份码本：变量定义、单位、源表字段、构造规则、缺失率与取值分布。
# 定义与规则优先取 plan/extraction/c3.csv；码本里的缺失与分布由本脚本在真实 parquet 上计算。
suppressPackageStartupMessages({ library(DBI); library(duckdb) })
con <- dbConnect(duckdb::duckdb())
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

OUT <- "data/processed/c3"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

sheet <- read.csv("plan/extraction/c3.csv", check.names = FALSE, fileEncoding = "UTF-8",
                  stringsAsFactors = FALSE)
names(sheet)[1] <- sub("^\ufeff", "", names(sheet)[1])

units <- c(
  age = "岁", age_numeric = "岁", sbp_min = "mmHg", rr_max = "次/分", spo2_min = "%",
  ph_min = "无", paco2_max = "mmHg", pf_min = "mmHg", lactate_max = "mmol/L",
  bicarbonate_first = "mmol/L", bicarbonate_min = "mmol/L", hr_min = "次/分", hr_max = "次/分",
  heart_rate_min = "次/分", heart_rate_max = "次/分", resp_rate_max = "次/分", rr_max = "次/分",
  apachescore = "分", gcs_total = "分", gcs_min = "分", sofa = "分",
  icu_los_days = "天", hosp_los_days = "天", niv_start_offset = "分钟",
  hospitaladmitoffset = "分钟", death_offset = "分钟", predictedhospitalmortality = "概率",
  hosp_ph_avail_pct = "%", niv_epap_mean = "cmH2O", niv_set_rr_mean = "次/分",
  hacor_score = "分", hacor_2h = "分", hacor_24h = "分", hacor_hr = "分", hacor_ph = "分",
  hacor_gcs = "分", hacor_pf = "分", hacor_rr = "分"
)
extra <- list(
  eICU = data.frame(
    变量名 = c("age", "ethnicity", "hacor_2h", "hacor_24h"),
    中文名 = c("年龄原文", "族裔", "0–2 小时 HACOR", "0–24 小时 HACOR"),
    类型 = c("直接取", "直接取", "复合", "复合"),
    源变量 = c("patient.age", "patient.ethnicity", "vitalPeriodic/lab/nurseCharting", "vitalPeriodic/lab/nurseCharting"),
    计算规则 = c("原值，可为 '> 89'", "原值", "NIV 后 0–2 小时五项计分求和", "NIV 后 0–24 小时五项计分求和"),
    stringsAsFactors = FALSE),
  `MIMIC-IV` = data.frame(
    变量名 = c("hacor_2h", "hacor_24h"),
    中文名 = c("0–2 小时 HACOR", "0–24 小时 HACOR"),
    类型 = c("复合", "复合"),
    源变量 = c("vitalsign/gcs/bg", "vitalsign/gcs/bg"),
    计算规则 = c("NIV 后 0–2 小时五项计分求和", "NIV 后 0–24 小时五项计分求和"),
    stringsAsFactors = FALSE)
)

make_one <- function(df, ds, path) {
  s <- sheet[sheet$数据集 == ds, ]
  x <- extra[[ds]]
  meta <- rbind(
    data.frame(变量名 = s$变量名, 中文名 = s$`中文名`, 类型 = s$类型,
               源变量 = s$源变量, 计算规则 = s$计算规则, stringsAsFactors = FALSE),
    x
  )
  meta <- meta[!duplicated(meta$变量名), ]
  vars <- names(df)
  meta <- meta[match(vars, meta$变量名), ]
  meta$变量名 <- vars
  meta$中文名[is.na(meta$中文名)] <- vars[is.na(meta$中文名)]
  meta$类型[is.na(meta$类型)] <- "直接取"
  meta$源变量[is.na(meta$源变量)] <- ""
  meta$计算规则[is.na(meta$计算规则)] <- "原值"

  res <- data.frame(
    变量名 = vars,
    中文名 = meta$中文名,
    类型 = meta$类型,
    单位 = ifelse(vars %in% names(units), units[vars], ""),
    来源表字段 = meta$源变量,
    构造规则 = meta$计算规则,
    非缺失数 = vapply(df, function(z) sum(!is.na(z)), numeric(1)),
    缺失数 = vapply(df, function(z) sum(is.na(z)), numeric(1)),
    `缺失率(%)` = round(100 * vapply(df, function(z) mean(is.na(z)), numeric(1)), 2),
    取值数 = vapply(df, function(z) length(unique(z[!is.na(z)])), integer(1)),
    最小值 = NA_real_, 最大值 = NA_real_, 取值分布前5 = "",
    check.names = FALSE, stringsAsFactors = FALSE
  )
  for (v in vars) {
    z <- df[[v]]
    if (is.numeric(z) || inherits(z, "difftime")) {
      zz <- suppressWarnings(as.numeric(z)); zz <- zz[is.finite(zz)]
      if (length(zz)) { res$最小值[res$变量名 == v] <- min(zz); res$最大值[res$变量名 == v] <- max(zz) }
    } else {
      tab <- sort(table(as.character(z[!is.na(z)])), decreasing = TRUE)
      res$取值分布前5[res$变量名 == v] <- paste0(names(head(tab, 5)), ":", as.integer(head(tab, 5)), collapse = "; ")
    }
  }
  write.csv(res, path, row.names = FALSE, fileEncoding = "UTF-8")
  cat("wrote", path, ":", nrow(res), "variables\n")
  res
}

e <- dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s/eicu_niv_analysis.parquet')", OUT))
m <- dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s/mimic_niv_analysis.parquet')", OUT))
ce <- make_one(e, "eICU", file.path(OUT, "eicu_niv_analysis_codebook.csv"))
cm <- make_one(m, "MIMIC-IV", file.path(OUT, "mimic_niv_analysis_codebook.csv"))

cat("eICU rows:", nrow(e), " cols:", ncol(e), " MIMIC rows:", nrow(m), " cols:", ncol(m), "\n")
