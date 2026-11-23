# c3_ext2_03_hacor_nri.R —— HACOR 在常规协变量之外的增量价值：NRI / IDI
# 输入：results/c3/_data/{eicu_prepared_v2,imp_eicu_v2}.rds（仅取协变量与 HACOR 子集）
# 输出：results/c3_ext2/tables/hacor_nri_idi.csv
# 说明：补充/敏感性分析；HACOR 可计算子集，不能外推全队列。
# 方法：参考模型 = 主协变量模型；新模型 = HACOR + 主协变量（logistic）。
#      连续 NRI（无数值切点）；分类 NRI（预测风险 <10%、10–20%、≥20% 三档，分析前定）；
#      IDI；bootstrap 500 次（按医院整群重抽样）给 95% 百分位区间。
source((function(f) { for (p in c(file.path("code/00_common", f), file.path("../00_common", f), file.path("00_common", f), f)) if (file.exists(p)) return(p); f })("common_transportability.R"))
set.seed(SEED)
dir.create(file.path(OUT, "tables"), recursive = TRUE, showWarnings = FALSE)

ei <- subset(load_eicu_ext2(), in_analysis == 1)
covs <- c("to_invasive_48h","age_numeric","gender","dx_group","apachescore",
          "gcs_total","sbp_min","rr_max","spo2_min")
hc <- ei[!is.na(ei$hacor_score) & complete.cases(ei[, covs]), ]
hc$age10 <- hc$age_numeric/10; hc$female <- as.integer(hc$gender=="Female")
hc$apache10 <- hc$apachescore/10; hc$gcs3 <- hc$gcs_total/3
hc$sbp10 <- hc$sbp_min/10; hc$rr5 <- hc$rr_max/5; hc$spo2_5 <- hc$spo2_min/5
hc$dx_group <- factor(hc$dx_group, levels = DX_LEVELS)
hc$hacor <- as.numeric(hc$hacor_score)
hc$m_apache <- as.integer(is.na(hc$apachescore)); hc$m_gcs <- as.integer(is.na(hc$gcs_total))
hc$hosp_id <- factor(hc$hospitalid)
y <- hc$to_invasive_48h

cat("HACOR non-missing:", sum(!is.na(ei$hacor_score)), "of", nrow(ei),
    sprintf("(%.1f%%)", 100*mean(!is.na(ei$hacor_score))), "\n")
cat("Complete-covariate subset n =", nrow(hc), " events =", sum(y),
    " hospitals =", length(unique(hc$hospitalid)), "\n")

f_ref <- to_invasive_48h ~ age10 + female + dx_group + apache10 + gcs3 +
  sbp10 + rr5 + spo2_5 + m_apache + m_gcs
f_new <- update(f_ref, . ~ . + hacor)

m_ref <- glm(f_ref, data = hc, family = binomial)
m_new <- glm(f_new, data = hc, family = binomial)
p_ref <- predict(m_ref, type = "response")
p_new <- predict(m_new, type = "response")

metrics <- function(y, p_ref, p_new, cuts = c(0.10, 0.20)) {
  up1 <- mean(p_new[y==1] > p_ref[y==1]); dn1 <- mean(p_new[y==1] < p_ref[y==1])
  dn0 <- mean(p_new[y==0] < p_ref[y==0]); up0 <- mean(p_new[y==0] > p_ref[y==0])
  g_ref <- as.integer(cut(p_ref, c(-Inf, cuts, Inf)))
  g_new <- as.integer(cut(p_new, c(-Inf, cuts, Inf)))
  cat_up1 <- mean((g_new > g_ref)[y==1]); cat_dn1 <- mean((g_new < g_ref)[y==1])
  cat_dn0 <- mean((g_new < g_ref)[y==0]); cat_up0 <- mean((g_new > g_ref)[y==0])
  c(nri_cont = (up1 - dn1) + (dn0 - up0),
    nri_events = up1 - dn1,
    nri_nonevents = dn0 - up0,
    nri_cat = (cat_up1 - cat_dn1) + (cat_dn0 - cat_up0),
    idi = (mean(p_new[y==1]) - mean(p_ref[y==1])) + (mean(p_ref[y==0]) - mean(p_new[y==0])),
    auc_ref = auc_safe(y, p_ref), auc_new = auc_safe(y, p_new))
}

pt <- metrics(y, p_ref, p_new)

# ---------------- 聚类 bootstrap ----------------
B <- 500
boot <- matrix(NA_real_, B, length(pt), dimnames = list(NULL, names(pt)))
hosp <- unique(hc$hospitalid); ok <- 0
set.seed(SEED)
for (b in 1:B) {
  samp <- sample(hosp, length(hosp), replace = TRUE)
  idx <- unlist(lapply(samp, function(h) which(hc$hospitalid == h)))
  db <- hc[idx, ]
  fb1 <- try(glm(f_ref, data = db, family = binomial), silent = TRUE)
  fb2 <- try(glm(f_new, data = db, family = binomial), silent = TRUE)
  if (any(sapply(list(fb1, fb2), inherits, "try-error"))) next
  pb1 <- try(predict(fb1, type = "response"), silent = TRUE)
  pb2 <- try(predict(fb2, type = "response"), silent = TRUE)
  if (any(sapply(list(pb1, pb2), inherits, "try-error"))) next
  v <- try(metrics(db$to_invasive_48h, as.numeric(pb1), as.numeric(pb2)), silent = TRUE)
  if (inherits(v, "try-error") || anyNA(v[1:5])) next
  boot[b, ] <- v; ok <- ok + 1
}
cat("Bootstrap replicates used:", ok, "/", B, "\n")
qs <- function(x) quantile(x, c(.025, .975), na.rm = TRUE)
boot_diff <- boot[, "auc_new"] - boot[, "auc_ref"]

out <- data.frame(
  metric = c("Continuous NRI", "  Events", "  Non-events", "Categorical NRI (<10%, 10-20%, >=20%)",
             "IDI", "AUC (reference model)", "AUC (HACOR + covariates)", "AUC difference"),
  estimate = c(pt["nri_cont"], pt["nri_events"], pt["nri_nonevents"], pt["nri_cat"],
               pt["idi"], pt["auc_ref"], pt["auc_new"], pt["auc_new"] - pt["auc_ref"]),
  ci_low = c(qs(boot[,"nri_cont"])[1], qs(boot[,"nri_events"])[1], qs(boot[,"nri_nonevents"])[1],
             qs(boot[,"nri_cat"])[1], qs(boot[,"idi"])[1], qs(boot[,"auc_ref"])[1],
             qs(boot[,"auc_new"])[1], qs(boot_diff)[1]),
  ci_high = c(qs(boot[,"nri_cont"])[2], qs(boot[,"nri_events"])[2], qs(boot[,"nri_nonevents"])[2],
              qs(boot[,"nri_cat"])[2], qs(boot[,"idi"])[2], qs(boot[,"auc_ref"])[2],
              qs(boot[,"auc_new"])[2], qs(boot_diff)[2])
)

out$n <- nrow(hc); out$events <- sum(y); out$bootstrap_ok <- ok
out$estimate <- round(out$estimate, 4); out$ci_low <- round(out$ci_low, 4); out$ci_high <- round(out$ci_high, 4)
write.csv(out, file.path(OUT, "tables", "hacor_nri_idi.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("\nHACOR NRI/IDI:\n"); print(out)
writeLines(capture.output(sessionInfo()), file.path(OUT, "_data", "session_info_hacor.txt"))
cat("\nhacor nri done\n")
