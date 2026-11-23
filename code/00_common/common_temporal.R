# c3temporal_common.R —— C3 时间顺序敏感性分析的公共函数
# 输入：results/c3/_data/eicu_prepared_v2.rds、imp_eicu_v2.rds
#       data/processed/c3_temporal/*.parquet
# 说明：主模型口径完全沿用 scripts/c3_common.R 与 c3_07_gee_main.R：
#       exchangeable GEE、医院聚类、m=20 套插补、Rubin 合并、v2 诊断映射。
# 本文件只做数据拼接与建模封装，不改动任何已有数据或结果。

suppressPackageStartupMessages({
  library(duckdb)
  library(mice)
  library(geepack)
  library(sandwich)
  library(data.table)
})
source((function(f) { for (p in c(file.path("code/00_common", f), file.path("../00_common", f), file.path("00_common", f), f)) if (file.exists(p)) return(p); f })("common_env.R"))

OUT <- "results/c3_temporal"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- 读取主分析准备集与时间顺序补取表，按主键对齐 ----
load_eicu_temporal <- function() {
  e  <- readRDS("results/c3/_data/eicu_prepared_v2.rds")
  con <- dbConnect(duckdb(shared_home = TRUE))
  t  <- dbGetQuery(con, "SELECT * FROM read_parquet('data/processed/c3_temporal/eicu_niv_temporal.parquet')")
  dbDisconnect(con, shutdown = TRUE)
  stopifnot(!anyDuplicated(t$patientunitstayid))
  idx <- match(e$patientunitstayid, t$patientunitstayid)
  stopifnot(all(!is.na(idx)))
  t <- t[idx, , drop = FALSE]
  stopifnot(identical(t$patientunitstayid, e$patientunitstayid))
  # 时间顺序列加前缀，避免与主表同名冲突（主表无同名列，仍显式加前缀便于追溯）
  add <- setdiff(names(t), c("patientunitstayid", "hospitalid", "in_analysis"))
  for (v in add) e[[paste0("tmp_", v)]] <- t[[v]]
  e
}

# ---- 与 c3_07_gee_main.R 完全一致地展开 20 套插补数据 ----
build_complete_list <- function(e, imp, extra_from_e = character(0)) {
  ei <- subset(e, in_analysis == 1)
  stopifnot(nrow(complete(imp, 1)) == nrow(ei))
  lapply(1:imp$m, function(i) {
    d <- complete(imp, i)
    d$to_invasive_48h <- as.integer(as.character(d$to_invasive_48h))
    d$dx_group <- factor(d$dx_group, levels = DX_LEVELS)
    d$age10 <- d$age_numeric / 10
    d$female <- as.integer(d$gender == "Female")
    d$apache10 <- d$apachescore / 10
    d$gcs3 <- d$gcs_total / 3
    d$sbp10 <- d$sbp_min / 10
    d$rr5 <- d$rr_max / 5
    d$spo2_5 <- d$spo2_min / 5
    d$m_apache <- as.integer(is.na(ei$apachescore))
    d$m_gcs <- as.integer(is.na(ei$gcs_total))
    d$hosp_id <- factor(ei$hospitalid)
    oc <- c("to_invasive_24h","to_invasive_any","to_invasive_48h_anydef",
            "death_48h","death_hosp","composite_48h","composite_any",
            "followup_48h_complete","first_icu_stay","patient_unit_seq",
            "dnr_before_niv","icu_los_days")
    for (v in intersect(oc, names(ei))) d[[v]] <- ei[[v]]
    # 时间顺序补取列（按位置对齐；ei 与 complete 行序一致）
    for (v in extra_from_e) d[[v]] <- ei[[v]]
    d
  })
}

# ---- exchangeable GEE（与主模型一致） ----
fit_exch <- function(d, formula, id = "hosp_id") {
  dd <- d[complete.cases(d[, all.vars(formula)]), , drop = FALSE]
  if (is.factor(dd$dx_group)) dd$dx_group <- droplevels(dd$dx_group)
  dd$.gee_id <- dd[[id]]
  geeglm(formula, data = dd, family = binomial(), id = .gee_id,
         corstr = "exchangeable", std.err = "san.se")
}

# ---- 在给定行子集上拟合并 Rubin 合并 ----
fit_pool_subset <- function(dl, formula, idx = NULL, id = "hosp_id") {
  if (!is.null(idx)) dl <- lapply(dl, function(d) d[idx, , drop = FALSE])
  fs <- lapply(dl, function(d) fit_exch(d, formula, id = id))
  list(fits = fs, pool = pool_rubin(fs))
}

# ---- OR 提取 ----
or_row <- function(pool, term, label, model, n, events, extra = NULL) {
  i <- which(pool$term == term)
  if (length(i) == 0) return(NULL)
  out <- data.frame(
    model = model, term = label,
    OR = exp(pool$estimate[i]), OR_low = exp(pool$conf.low[i]),
    OR_high = exp(pool$conf.high[i]), p.value = pool$p.value[i],
    n = n, events = events, stringsAsFactors = FALSE)
  if (!is.null(extra)) out <- cbind(out, extra)
  out
}

fmt_or <- function(or, lo, hi) sprintf("%.2f (%.2f-%.2f)", or, lo, hi)
