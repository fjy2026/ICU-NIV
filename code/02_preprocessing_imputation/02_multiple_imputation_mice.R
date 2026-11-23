# c3_02_impute.R —— 单独重建 eICU 多重插补对象（供追踪步骤留证）
# 与 c3_02_main.R 的插补块完全一致（m=20, maxit=10, seed=20261113）
source((function(f) { for (p in c(file.path("code/00_common", f), file.path("../00_common", f), file.path("00_common", f), f)) if (file.exists(p)) return(p); f })("common_env.R"))
set.seed(SEED)
e <- readRDS("results/c3/_data/eicu_prepared.rds")
ei <- subset(e, in_analysis == 1)
imp_vars <- c("to_invasive_48h","age_numeric","gender","dx_group","apachescore","gcs_total",
              "sbp_min","rr_max","spo2_min","ph_min","paco2_max","pf_min","lactate_max",
              "bicarbonate_first","hosp_nivfirst_n","hosp_ph_avail_pct")
dimp <- ei[, imp_vars]
dimp$to_invasive_48h <- factor(dimp$to_invasive_48h, levels = c(0, 1))
dimp$gender <- factor(dimp$gender)
meth <- make.method(dimp)
if ("gender" %in% names(meth)) meth["gender"] <- "logreg"
imp <- mice(dimp, m = 20, maxit = 10, seed = SEED, method = meth, printFlag = FALSE)
saveRDS(imp, "results/c3/_data/imp_eicu.rds")
cat("imputation rebuilt: m =", imp$m, "\n")
