# c3_18_sensitivity_v2.R —— S1–S9 全部改用 v2 诊断映射重跑
# 输入：results/c3/_data/eicu_prepared_v2.rds、imp_eicu_v2.rds
# 输出：results/c3/tables/sensitivity_summary.csv、sensitivity_detail.csv、evalues.csv
#       results/c3/_data/imp_s2a_v2.rds、imp_s2b_v2.rds、imp_s4c_v2.rds
# 运行：Rscript scripts/c3_18_sensitivity_v2.R
#
# 估计方法：
#   - S1、S2、S3、S5、S6、S9：exchangeable GEE，医院聚类，m=20，Rubin 合并。
#   - S4a：exchangeable GEE，完整病例（按定义不插补），单次拟合。
#   - S4c：exchangeable GEE，m=50、maxit=5（沿用原设定），Rubin 合并。
#   - S7：随机截距 logistic（glmer, nAGQ=0），m=20，Rubin 合并（S7 本身就是检验
#     随机截距模型，故不使用 GEE）。
#   - S2b-i 与主分析完全相同，直接复用主分析合并结果。
# 每条结果均记录 method、n、events，不写死数字。

source((function(f) { for (p in c(file.path("code/00_common", f), file.path("../00_common", f), file.path("00_common", f), f)) if (file.exists(p)) return(p); f })("common_env.R"))
set.seed(SEED)

e <- readRDS("results/c3/_data/eicu_prepared_v2.rds")
ei <- subset(e, in_analysis == 1)
imp <- readRDS("results/c3/_data/imp_eicu_v2.rds")

fit_exch <- function(d, formula, id = "hosp_id") {
  dd <- d[complete.cases(d[, all.vars(formula)]), ]
  if (is.factor(dd$dx_group)) dd$dx_group <- droplevels(dd$dx_group)
  dd$.gee_id <- dd[[id]]
  geeglm(formula, data = dd, family = binomial(), id = .gee_id,
         corstr = "exchangeable", std.err = "san.se")
}

make_complete <- function(imp_obj, ref) {
  lapply(1:imp_obj$m, function(i) {
    d <- complete(imp_obj, i)
    d$to_invasive_48h <- as.integer(as.character(d$to_invasive_48h))
    d$dx_group <- factor(d$dx_group, levels = DX_LEVELS)
    d$age10 <- d$age_numeric/10; d$female <- as.integer(d$gender == "Female")
    d$apache10 <- d$apachescore/10; d$gcs3 <- d$gcs_total/3
    d$sbp10 <- d$sbp_min/10; d$rr5 <- d$rr_max/5; d$spo2_5 <- d$spo2_min/5
    d$m_apache <- as.integer(is.na(ref$apachescore)); d$m_gcs <- as.integer(is.na(ref$gcs_total))
    d$m_ph <- as.integer(is.na(ref$ph_min)); d$m_paco2 <- as.integer(is.na(ref$paco2_max))
    d$m_pf <- as.integer(is.na(ref$pf_min)); d$m_lactate <- as.integer(is.na(ref$lactate_max))
    d$hosp_id <- factor(ref$hospitalid)
    for (oc in c("to_invasive_24h","to_invasive_any","to_invasive_48h_anydef",
                 "death_48h","death_hosp","composite_48h","followup_48h_complete",
                 "first_icu_stay","patient_unit_seq"))
      d[[oc]] <- ref[[oc]]
    d
  })
}
complete_list <- make_complete(imp, ei)

fit_pool <- function(dl, outcome = "to_invasive_48h", terms = MAIN_TERMS, idx = NULL) {
  if (!is.null(idx)) dl <- lapply(dl, function(d) d[idx, , drop = FALSE])
  fs <- lapply(dl, function(d) {
    dd <- d; dd$y <- dd[[outcome]]; dd <- dd[!is.na(dd$y), ]
    fit_exch(dd, as.formula(paste("y ~", paste(terms, collapse = " + "))))
  })
  pool_rubin(fs)
}

key <- DIRECTION_TERMS
key_labels <- c(age10="年龄/10 岁", female="女性", apache10="严重度/单位", gcs3="GCS/3 分",
                sbp10="收缩压/10 mmHg", rr5="呼吸频率/5", spo2_5="SpO2/5%")

main_pool <- fit_pool(complete_list)
main_or <- setNames(exp(main_pool$estimate), main_pool$term)
main_sig <- setNames(main_pool$p.value < 0.05, main_pool$term)

changed_terms <- function(p) {
  eval_terms <- c("age10","apache10","gcs3","sbp10")
  ch <- character(0)
  for (tt in eval_terms) {
    or <- exp(p$estimate[p$term == tt]); lo <- exp(p$conf.low[p$term == tt]); hi <- exp(p$conf.high[p$term == tt])
    if (length(or) == 0) next
    if (isTRUE(sign(log(or)) != sign(log(main_or[tt]))) ||
        (main_sig[tt] && lo <= 1 && hi >= 1)) ch <- c(ch, key_labels[tt])
  }
  if (length(ch) == 0) "无" else paste(ch, collapse = "、")
}

detail_rows <- list(); summary_rows <- list()
add_sens <- function(id, desc, p, n, events, method, note = "") {
  for (tt in key) {
    detail_rows[[length(detail_rows)+1]] <<- data.frame(
      id = id, description = desc, term = key_labels[tt],
      OR = sprintf("%.2f (%.2f–%.2f)", exp(p$estimate[p$term==tt]),
                   exp(p$conf.low[p$term==tt]), exp(p$conf.high[p$term==tt])),
      p = signif(p$p.value[p$term==tt], 3), n = n, events = events,
      method = method, note = note, stringsAsFactors = FALSE)
  }
  summary_rows[[length(summary_rows)+1]] <<- data.frame(
    id = id, description = desc, n = n, events = events,
    changed = changed_terms(p), method = method, note = note, stringsAsFactors = FALSE)
}
n_main <- nrow(ei); ev_main <- sum(ei$to_invasive_48h)

GEE20 <- "exchangeable GEE，m=20，Rubin"
GEE50 <- "exchangeable GEE，m=50，Rubin"
GEE_CC <- "exchangeable GEE，完整病例（不插补）"
RI20 <- "随机截距 logistic（glmer, nAGQ=0），m=20，Rubin"

# S1 结局口径
add_sens("S1a", "主结局改为任意时间转有创", fit_pool(complete_list, "to_invasive_any"),
         n_main, sum(complete_list[[1]]$to_invasive_any, na.rm=TRUE), GEE20)
add_sens("S1b", "主结局改为复合（48h 转有创或院内死亡）", fit_pool(complete_list, "composite_48h"),
         sum(!is.na(complete_list[[1]]$composite_48h)), sum(complete_list[[1]]$composite_48h, na.rm=TRUE), GEE20)
add_sens("S1c", "主结局改为 24 小时内转有创", fit_pool(complete_list, "to_invasive_24h"),
         n_main, sum(complete_list[[1]]$to_invasive_24h, na.rm=TRUE), GEE20)
add_sens("S1d", "有创定义改宽口径（含复合路径）", fit_pool(complete_list, "to_invasive_48h_anydef"),
         n_main, sum(complete_list[[1]]$to_invasive_48h_anydef, na.rm=TRUE), GEE20)

# S2a 全无创队列（v2 重新插补）
e_all <- subset(e, !is.na(age_numeric))
dall <- e_all[, c("to_invasive_48h","age_numeric","gender","dx_group","apachescore","gcs_total",
                  "sbp_min","rr_max","spo2_min","ph_min","paco2_max","pf_min","lactate_max",
                  "bicarbonate_first","hosp_nivfirst_n","hosp_ph_avail_pct")]
dall$to_invasive_48h <- factor(dall$to_invasive_48h, levels = c(0,1)); dall$gender <- factor(dall$gender)
meth <- make.method(dall); if ("gender" %in% names(meth)) meth["gender"] <- "logreg"
if (file.exists("results/c3/_data/imp_s2a_v2.rds")) {
  imp_all <- readRDS("results/c3/_data/imp_s2a_v2.rds")
} else {
  imp_all <- mice(dall, m = 20, maxit = 10, seed = SEED, method = meth, printFlag = FALSE)
  saveRDS(imp_all, "results/c3/_data/imp_s2a_v2.rds")
}
cl_all <- make_complete(imp_all, e_all)
add_sens("S2a", "队列改为全无创（含序贯撤机）", fit_pool(cl_all),
         nrow(e_all), sum(e_all$to_invasive_48h, na.rm=TRUE), GEE20)

# S2b 同刻记录
add_sens("S2b-i", "同刻记录归为有创在前（等于主分析）", main_pool, n_main, ev_main,
         "exchangeable GEE，m=20，Rubin（同主分析）")
rows2b <- which(e$in_analysis == 1 | e$inv_same_time == 1)
e_s2b <- e[rows2b, ]
d2b <- e_s2b[, c("to_invasive_48h","age_numeric","gender","dx_group","apachescore","gcs_total",
                 "sbp_min","rr_max","spo2_min","ph_min","paco2_max","pf_min","lactate_max",
                 "bicarbonate_first","hosp_nivfirst_n","hosp_ph_avail_pct")]
d2b$to_invasive_48h <- factor(d2b$to_invasive_48h, levels = c(0,1)); d2b$gender <- factor(d2b$gender)
meth <- make.method(d2b); if ("gender" %in% names(meth)) meth["gender"] <- "logreg"
if (file.exists("results/c3/_data/imp_s2b_v2.rds")) {
  imp_2b <- readRDS("results/c3/_data/imp_s2b_v2.rds")
} else {
  imp_2b <- mice(d2b, m = 20, maxit = 10, seed = SEED, method = meth, printFlag = FALSE)
  saveRDS(imp_2b, "results/c3/_data/imp_s2b_v2.rds")
}
cl_2b <- make_complete(imp_2b, e_s2b)
add_sens("S2b-ii", "同刻记录归为 NIV-first", fit_pool(cl_2b),
         nrow(e_s2b), sum(e_s2b$to_invasive_48h, na.rm=TRUE), GEE20)

# S3 患者层面去重
idx3a <- which(ei$first_icu_stay == 1)
add_sens("S3a", "仅本次住院首次 ICU 入住", fit_pool(complete_list, idx = idx3a),
         length(idx3a), sum(ei$to_invasive_48h[idx3a], na.rm=TRUE), GEE20)
idx3b <- which(ei$patient_unit_seq == 1)
add_sens("S3b", "每位患者仅一条（patientunitstayid 最小）", fit_pool(complete_list, idx = idx3b),
         length(idx3b), sum(ei$to_invasive_48h[idx3b], na.rm=TRUE), GEE20)

# S4 缺失处理
cc_ei <- ei[complete.cases(ei[, c("to_invasive_48h","age_numeric","gender","dx_group","apachescore",
                                  "gcs_total","sbp_min","rr_max","spo2_min")]), ]
cc_ei$age10 <- cc_ei$age_numeric/10; cc_ei$female <- as.integer(cc_ei$gender=="Female")
cc_ei$apache10 <- cc_ei$apachescore/10; cc_ei$gcs3 <- cc_ei$gcs_total/3
cc_ei$sbp10 <- cc_ei$sbp_min/10; cc_ei$rr5 <- cc_ei$rr_max/5; cc_ei$spo2_5 <- cc_ei$spo2_min/5
cc_ei$m_apache <- 0; cc_ei$m_gcs <- 0; cc_ei$dx_group <- factor(cc_ei$dx_group, levels = DX_LEVELS)
cc_ei$hosp_id <- factor(cc_ei$hospitalid)
fit_cc <- fit_exch(cc_ei, main_formula(terms = setdiff(MAIN_TERMS, c("m_apache","m_gcs"))))
add_sens("S4a", "完整病例分析（不插补）", tidy_fit(fit_cc),
         nrow(cc_ei), sum(cc_ei$to_invasive_48h), GEE_CC)
add_sens("S4b", "主模型去掉缺失指示变量",
         fit_pool(complete_list, terms = setdiff(MAIN_TERMS, c("m_apache","m_gcs"))),
         n_main, ev_main, GEE20)

# S4c 50 套插补（v2 重新插补，maxit=5 沿用原设定）
dmain50 <- ei[, c("to_invasive_48h","age_numeric","gender","dx_group","apachescore","gcs_total",
                  "sbp_min","rr_max","spo2_min","ph_min","paco2_max","pf_min","lactate_max",
                  "bicarbonate_first","hosp_nivfirst_n","hosp_ph_avail_pct")]
dmain50$to_invasive_48h <- factor(dmain50$to_invasive_48h, levels = c(0,1))
dmain50$gender <- factor(dmain50$gender)
meth50 <- make.method(dmain50); if ("gender" %in% names(meth50)) meth50["gender"] <- "logreg"
if (file.exists("results/c3/_data/imp_s4c_v2.rds")) {
  imp50 <- readRDS("results/c3/_data/imp_s4c_v2.rds")
} else {
  imp50 <- mice(dmain50, m = 50, maxit = 5, seed = SEED, method = meth50, printFlag = FALSE)
  saveRDS(imp50, "results/c3/_data/imp_s4c_v2.rds")
}
cl50 <- make_complete(imp50, ei)
for (i in seq_along(cl50)) cl50[[i]]$ph_avail10 <- cl50[[i]]$hosp_ph_avail_pct * 10
p50 <- fit_pool(cl50, terms = c(MAIN_TERMS, "ph_avail10"))
add_sens("S4c", "50 套插补并加入医院血气可得比例", p50,
         nrow(ei), sum(ei$to_invasive_48h, na.rm=TRUE), GEE50,
         note = "maxit=5（沿用原设定）")

# S5 计划插管
ei$delta_min <- ei$inv_start_strict - ei$niv_start_offset
idx5a <- which(is.na(ei$delta_min) | !(ei$delta_min > 0 & ei$delta_min <= 60))
add_sens("S5a", "排除无创后 60 分钟内即转有创者", fit_pool(complete_list, idx = idx5a),
         length(idx5a), sum(ei$to_invasive_48h[idx5a], na.rm=TRUE), GEE20)
idx5b <- which(ei$dx_group != "术后/其他")
add_sens("S5b", "排除术后/其他诊断组（v2 映射）", fit_pool(complete_list, idx = idx5b),
         length(idx5b), sum(ei$to_invasive_48h[idx5b], na.rm=TRUE), GEE20)

# S6 医院
keep6a <- which(ei$hosp_nivfirst_n >= 50)
add_sens("S6a", "排除 NIV-first 例数 < 50 的医院", fit_pool(complete_list, idx = keep6a),
         length(keep6a), sum(ei$to_invasive_48h[keep6a], na.rm=TRUE), GEE20)

hosp_rate <- hospital_observed(ei)
cand <- hosp_rate[hosp_rate$n >= 20, ]
top3 <- cand$hospitalid[order(-cand$rate)][1:3]
bot3 <- cand$hospitalid[order(cand$rate)][1:3]
for (k in 1:3) {
  idx <- which(!(ei$hospitalid %in% top3[k]))
  add_sens(paste0("S6b-top", k), paste0("剔除转有创率最高的第 ", k, " 家医院（例数≥20）"),
           fit_pool(complete_list, idx = idx), length(idx),
           sum(ei$to_invasive_48h[idx], na.rm=TRUE), GEE20)
  idx <- which(!(ei$hospitalid %in% bot3[k]))
  add_sens(paste0("S6b-bot", k), paste0("剔除转有创率最低的第 ", k, " 家医院（例数≥20）"),
           fit_pool(complete_list, idx = idx), length(idx),
           sum(ei$to_invasive_48h[idx], na.rm=TRUE), GEE20)
}

# S7 随机截距模型
f_ri <- as.formula(paste("to_invasive_48h ~", paste(MAIN_TERMS, collapse=" + "), "+ (1 | hosp_id)"))
fits_ri <- lapply(complete_list, function(d) glmer(f_ri, data = d, family = binomial, nAGQ = 0))
p_ri <- pool_rubin(fits_ri)
add_sens("S7", "随机截距模型替代 GEE", p_ri, n_main, ev_main, RI20,
         note = sprintf("随机截距方差 %.3f",
                        mean(sapply(fits_ri, function(f) as.numeric(VarCorr(f)$hosp_id)))))

# S8 E 值（由最终主模型计算，不单独拟合）
ev_rows <- list()
for (tt in c("age10","apache10","sbp10","gcs3")) {
  or <- main_or[tt]; lo <- exp(main_pool$conf.low[main_pool$term == tt])
  hi <- exp(main_pool$conf.high[main_pool$term == tt])
  ev <- EValue::evalues.OR(est = or, lo = lo, hi = hi, rare = FALSE)
  evp <- ev["E-values", "point"]
  evc <- if (or < 1) ev["E-values", "upper"] else ev["E-values", "lower"]
  ev_rows[[length(ev_rows)+1]] <- data.frame(
    variable = key_labels[tt], OR = sprintf("%.2f (%.2f–%.2f)", or, lo, hi),
    E_value_point = sprintf("%.2f", evp), E_value_CI = sprintf("%.2f", evc))
}
write.csv(do.call(rbind, ev_rows), "results/c3/tables/evalues.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

# S9 排除随访不足 48 小时
idx9 <- which(ei$followup_48h_complete == 1)
add_sens("S9", "排除观察不足 48 小时者", fit_pool(complete_list, idx = idx9),
         length(idx9), sum(ei$to_invasive_48h[idx9], na.rm=TRUE), GEE20)

summary_tab <- do.call(rbind, summary_rows)
detail_tab <- do.call(rbind, detail_rows)
write.csv(summary_tab, "results/c3/tables/sensitivity_summary.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(detail_tab, "results/c3/tables/sensitivity_detail.csv", row.names = FALSE, fileEncoding = "UTF-8")

cat("sensitivity v2 done\n")
print(summary_tab[, c("id","n","events","changed","method")])
