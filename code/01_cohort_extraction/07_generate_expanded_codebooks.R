#!/usr/bin/env Rscript
# 生成 C3 扩取两套数据的码本（逐列中文名、来源、定义、缺失与取值）
# 中文名与定义来自 plan/extraction/c3_expanded.csv；旧列沿用 c3 旧码本。
suppressPackageStartupMessages({ library(DBI); library(duckdb) })
con <- dbConnect(duckdb::duckdb())
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

es <- read.csv("plan/extraction/c3_expanded.csv", fileEncoding = "UTF-8-BOM",
               stringsAsFactors = FALSE, check.names = FALSE)
old_e <- read.csv("data/processed/c3/eicu_niv_analysis_codebook.csv", fileEncoding = "UTF-8-BOM",
                  stringsAsFactors = FALSE, check.names = FALSE)
old_m <- read.csv("data/processed/c3/mimic_niv_analysis_codebook.csv", fileEncoding = "UTF-8-BOM",
                  stringsAsFactors = FALSE, check.names = FALSE)
lookup <- function(ds, v, old) {
  hit <- old[old$变量名 == v, ]
  if (nrow(hit)) return(list(cn = hit$中文名[1], type = hit$类型[1], unit = hit$单位[1],
                              src = hit$来源表字段[1], rule = hit$构造规则[1]))
  hit2 <- es[es$数据集 == ds & es$变量名 == v, ]
  if (nrow(hit2)) return(list(cn = hit2$中文名[1], type = hit2$类型[1], unit = "",
                              src = hit2$源变量[1], rule = paste(hit2$定义[1], hit2$计算规则[1])))
  list(cn = "", type = "", unit = "", src = "", rule = "")
}

make_cb <- function(path, ds, old, out) {
  d <- dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s')", path))
  rows <- lapply(names(d), function(v) {
    x <- d[[v]]; lk <- lookup(ds, v, old)
    nmiss <- sum(is.na(x)); n <- length(x) - nmiss
    typ <- if (is.numeric(x)) "数值" else if (inherits(x, "Date") || inherits(x, "POSIXct")) "时间" else "字符"
    if (!nzchar(lk$type)) lk$type <- typ
    if (is.numeric(x) && n > 0) {
      mn <- suppressWarnings(min(x, na.rm = TRUE)); mx <- suppressWarnings(max(x, na.rm = TRUE))
      mn <- if (is.finite(mn)) mn else NA; mx <- if (is.finite(mx)) mx else NA
    } else { mn <- NA; mx <- NA }
    vals <- if (n > 0) {
      tb <- sort(table(x, useNA = "no"), decreasing = TRUE)
      paste(paste0(names(tb)[seq_len(min(5, length(tb)))], "=", as.integer(tb)[seq_len(min(5, length(tb)))]), collapse = "; ")
    } else ""
    data.frame(变量名 = v, 中文名 = lk$cn, 类型 = lk$type, 单位 = lk$unit,
               来源表字段 = lk$src, 构造规则 = lk$rule,
               非缺失数 = n, 缺失数 = nmiss,
               缺失率 = if (length(x) > 0) round(100 * nmiss / length(x), 1) else NA,
               取值数 = length(unique(x[!is.na(x)])), 最小值 = mn, 最大值 = mx,
               取值分布前5 = vals, check.names = FALSE, stringsAsFactors = FALSE)
  })
  cb <- do.call(rbind, rows)
  write.csv(cb, out, row.names = FALSE, fileEncoding = "UTF-8")
  cat("wrote", out, "rows", nrow(cb), "\n")
}

make_cb("data/processed/c3_expanded/eicu_niv_expanded.parquet", "eICU", old_e,
        "data/processed/c3_expanded/eicu_niv_expanded_codebook.csv")
make_cb("data/processed/c3_expanded/mimic_niv_expanded.parquet", "MIMIC", old_m,
        "data/processed/c3_expanded/mimic_niv_expanded_codebook.csv")
