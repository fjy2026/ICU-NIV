# c3temporal_04_summary_table.R —— 汇总三组敏感性分析的关键估计
# 输入：results/c3_temporal/tableA_early_covariate.csv
#       results/c3_temporal/tableB_landmark24.csv
#       results/c3_temporal/tableC_vasopressor.csv
#       results/c3_temporal/tableC_vasopressor_mimic.csv
# 输出：results/c3_temporal/summary_all.csv
# 说明：只读取上表已有数值并重排，不重新建模、不写死任何结果数字。

a  <- read.csv("results/c3_temporal/tableA_early_covariate.csv",
               stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")
b  <- read.csv("results/c3_temporal/tableB_landmark24.csv",
               stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")
c  <- read.csv("results/c3_temporal/tableC_vasopressor.csv",
               stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")
cm <- read.csv("results/c3_temporal/tableC_vasopressor_mimic.csv",
               stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")

mk <- function(section, model, term, or, lo, hi, p, n, events, contrast) {
  data.frame(section = section, model = model, term = term,
             OR = round(or, 3), OR_low = round(lo, 3), OR_high = round(hi, 3),
             p = formatC(p, format = "g", digits = 3),
             n = n, events = events, contrast = contrast,
             stringsAsFactors = FALSE)
}

out <- list()

# A: 联合替换（参照与早期值成对）
aj <- a[grepl("^Joint: all four replaced", a$analysis), ]
for (i in seq_len(nrow(aj))) {
  r <- aj[i, ]
  out[[length(out) + 1]] <- mk("A_joint", r$variant, r$term_label,
                               r$ref_OR, r$ref_low, r$ref_high, r$ref_p,
                               r$n, r$events, "reference (0-24 h)")
  out[[length(out) + 1]] <- mk("A_joint", r$variant, r$term_label,
                               r$early_OR, r$early_low, r$early_high, r$early_p,
                               r$n, r$events, "early value")
}

# A: 逐一替换（只保留被替换项）
as <- a[grepl("^Single: one replaced", a$analysis) & a$term == a$replaced, ]
for (i in seq_len(nrow(as))) {
  r <- as[i, ]
  out[[length(out) + 1]] <- mk("A_single", r$variant, r$term_label,
                               r$early_OR, r$early_low, r$early_high, r$early_p,
                               r$n, r$events, "early value")
}

# B
for (i in seq_len(nrow(b))) {
  r <- b[i, ]
  out[[length(out) + 1]] <- mk("B_landmark24", "landmark", r$term,
                               r$OR, r$OR_low, r$OR_high, r$p.value,
                               r$n, r$events, "24-48 h outcome")
}

# C eICU
for (i in seq_len(nrow(c))) {
  r <- c[i, ]
  out[[length(out) + 1]] <- mk("C_eICU", r$model, r$term,
                               r$OR, r$OR_low, r$OR_high, r$p.value,
                               r$n, r$events, "with vasopressor")
}

# C MIMIC
for (i in seq_len(nrow(cm))) {
  r <- cm[i, ]
  out[[length(out) + 1]] <- mk("C_MIMIC", r$model, r$term,
                               r$OR, r$OR_low, r$OR_high, r$p.value,
                               r$n, r$events, "with vasopressor")
}

res <- do.call(rbind, out)
write.csv(res, "results/c3_temporal/summary_all.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
cat("summary rows:", nrow(res), "\n")
