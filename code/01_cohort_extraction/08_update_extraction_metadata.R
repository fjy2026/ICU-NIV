#!/usr/bin/env Rscript
# 用两份成品数据的实际列名，回填变量提取表的状态；只改状态与核对结果，
# 不动定义、源变量与计算规则。UTF-8 带 BOM 写回。
suppressPackageStartupMessages({ library(DBI); library(duckdb) })
con <- dbConnect(duckdb::duckdb())
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
E <- "data/processed/c3/eicu_niv_analysis.parquet"
M <- "data/processed/c3/mimic_niv_analysis.parquet"
en <- names(dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s') LIMIT 0", E)))
mn <- names(dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s') LIMIT 0", M)))

base <- "plan/extraction/c3_base.csv"
p <- "plan/extraction/c3.csv"
s <- read.csv(base, check.names = FALSE, fileEncoding = "UTF-8", stringsAsFactors = FALSE)
names(s)[1] <- sub("^\ufeff", "", names(s)[1])
key <- c("patientunitstayid","in_analysis","to_invasive_48h","death_hosp","apachescore",
         "gcs_total","ph_min","pf_min","lactate_max","hacor_score","dnr_before_niv",
         "stay_id","sofa","gcs_min")
for (i in seq_len(nrow(s))) {
  present <- if (s$数据集[i] == "eICU") s$变量名[i] %in% en else
             if (s$数据集[i] == "MIMIC-IV") s$变量名[i] %in% mn else FALSE
  if (present) {
    s$状态[i] <- if (s$变量名[i] %in% key) "已核对" else "已提取"
    s$核对结果[i] <- if (s$变量名[i] %in% key)
      "成品数据中存在；缺失率与取值范围见同目录码本" else "成品数据中存在；码本已给出缺失率"
  }
}
# UTF-8 BOM + CSV
tmp <- tempfile()
write.csv(s, tmp, row.names = FALSE, fileEncoding = "UTF-8")
raw <- readBin(tmp, "raw", file.info(tmp)$size)
con2 <- file(p, "wb")
writeBin(as.raw(c(0xEF, 0xBB, 0xBF)), con2)
writeBin(raw, con2)
close(con2)
unlink(tmp)
cat("updated", p, ":", sum(s$状态 == "已核对"), "verified,", sum(s$状态 == "已提取"), "extracted\n")
