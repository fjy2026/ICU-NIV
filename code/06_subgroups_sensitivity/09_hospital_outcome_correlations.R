# c3_ext3_01_hospital_outcomes.R —— ① 医院升级率（SER）与风险调整院内死亡（SMR）
# 输入：results/c3/_data（只读）、data/processed/c3_expanded（只读）
# 输出：results/c3_ext3/tables/*、figures/fig_hospital_ser_smr.png、fig_patient_interaction.png
# 说明：医院层面关联只能解释为“医院情境”，不能写成“提高升级率能改善生存”。
source((function(f) { for (p in c(file.path("code/00_common", f), file.path("../00_common", f), file.path("00_common", f), f)) if (file.exists(p)) return(p); f })("common_multilevel.R"))
set.seed(SEED)

L <- load_eicu_ext3()
cl <- L$cl
ei <- L$ei
n <- nrow(cl[[1]])
cat("eICU analysis set:", n, "admissions;", length(unique(cl[[1]]$hospitalid)), "hospitals\n")

# ------------------------------------------------------------------
# 1. 主模型预测的升级概率（20 套插补平均）
# ------------------------------------------------------------------
f_esc <- fml("to_invasive_48h", paste(MAIN_TERMS, collapse = " + "))
p_esc <- fit_glm_pred(cl, f_esc)
y_esc <- cl[[1]]$to_invasive_48h
cat("Escalation AUC (apparent, averaged predictions):",
    round(as.numeric(pROC::auc(pROC::roc(y_esc, p_esc, quiet = TRUE))), 3), "\n")

# 医院随机截距（与 SER 并列的第二种口径）
f_ri <- fml("to_invasive_48h", paste(c(MAIN_TERMS, "(1|hosp)"), collapse = " + "))
ri_fits <- lapply(cl, function(d) {
  tryCatch(glmer(f_ri, data = d, family = binomial, nAGQ = 0,
                 control = glmerControl(optimizer = "bobyqa",
                                        optCtrl = list(maxfun = 2e5),
                                        calc.derivs = FALSE)),
           error = function(e) NULL)
})
ri_ok <- !vapply(ri_fits, is.null, TRUE)
cat("Random-intercept models converged:", sum(ri_ok), "/", length(ri_fits), "\n")
ri_fits <- ri_fits[ri_ok]
eb_list <- lapply(ri_fits, function(f) {
  re <- ranef(f)$hosp
  data.frame(hospitalid = as.numeric(rownames(re)), eb = re[[1]])
})
eb <- Reduce(function(a, b) merge(a, b, by = "hospitalid", suffixes = c("", ".x")), eb_list)
eb$eb_mean <- rowMeans(eb[, -1, drop = FALSE], na.rm = TRUE)
v_vec <- vapply(ri_fits, function(f) as.numeric(VarCorr(f)$hosp[1, 1]), 0)
v_ri <- mean(v_vec, na.rm = TRUE)
icc_ri <- v_ri / (v_ri + pi^2 / 3)
mor_ri <- exp(sqrt(2 * v_ri) * qnorm(0.75))
cat(sprintf("Random-intercept variance %.3f; ICC %.3f; MOR %.2f\n", v_ri, icc_ri, mor_ri))

# ------------------------------------------------------------------
# 2. 院内死亡模型与 SMR（年龄、性别、严重度、GCS、SBP、诊断组、DNR）
# ------------------------------------------------------------------
f_death <- fml("death_hosp", paste(DEATH_TERMS, collapse = " + "))
ok_death <- !is.na(cl[[1]]$death_hosp) & cl[[1]]$dnr_known == 1
cat("Death model analysable:", sum(ok_death), "admissions;",
    sum(cl[[1]]$death_hosp[ok_death]), "deaths\n")
p_death <- rep(NA_real_, n)
pmat <- matrix(NA_real_, sum(ok_death), length(cl))
for (i in seq_along(cl)) {
  d <- cl[[i]][ok_death, , drop = FALSE]
  d$death_hosp <- as.integer(d$death_hosp)
  d$dnr0 <- as.integer(d$dnr0)
  fit <- glm(f_death, data = d, family = binomial)
  pmat[, i] <- predict(fit, newdata = d, type = "response")
}
p_death[ok_death] <- rowMeans(pmat, na.rm = TRUE)
cat("Deaths expected (sum p):", round(sum(p_death[ok_death]), 1),
    "; observed:", sum(cl[[1]]$death_hosp[ok_death]), "\n")

# ------------------------------------------------------------------
# 3. 医院层面汇总
# ------------------------------------------------------------------
pat <- data.frame(hospitalid = cl[[1]]$hospitalid,
                  y_esc = y_esc, p_esc = p_esc,
                  y_death = cl[[1]]$death_hosp, p_death = p_death,
                  death_ok = as.integer(ok_death))
agg_e <- aggregate(pat[, c("y_esc", "p_esc")], by = list(hospitalid = pat$hospitalid), FUN = sum)
names(agg_e) <- c("hospitalid", "O_esc", "E_esc")
agg_e$n_esc <- as.numeric(table(pat$hospitalid)[as.character(agg_e$hospitalid)])
agg_e$SER <- agg_e$O_esc / agg_e$E_esc
ci_ser <- pois_ci(agg_e$O_esc, agg_e$E_esc)
agg_e$SER_lo <- ci_ser[, 1]; agg_e$SER_hi <- ci_ser[, 2]

pd <- pat[pat$death_ok == 1, ]
agg_d <- aggregate(pd[, c("y_death", "p_death")], by = list(hospitalid = pd$hospitalid), FUN = sum)
names(agg_d) <- c("hospitalid", "O_death", "E_death")
agg_d$n_death <- as.numeric(table(pd$hospitalid)[as.character(agg_d$hospitalid)])
agg_d$SMR <- agg_d$O_death / agg_d$E_death
ci_smr <- pois_ci(agg_d$O_death, agg_d$E_death)
agg_d$SMR_lo <- ci_smr[, 1]; agg_d$SMR_hi <- ci_smr[, 2]

agg <- merge(agg_e, agg_d, by = "hospitalid", all.x = TRUE)
agg <- merge(agg, L$hosp, by = "hospitalid", all.x = TRUE)
agg <- merge(agg, eb[, c("hospitalid", "eb_mean")], by = "hospitalid", all.x = TRUE)
agg$logSER <- log(agg$SER)
agg$ser_z <- as.numeric(scale(agg$SER))

# SER 与随机截距的一致性
a_all <- agg[is.finite(agg$SER) & is.finite(agg$eb_mean), ]
a_pos <- a_all[a_all$SER > 0, ]
rho_all <- cor(log(a_pos$SER), a_pos$eb_mean, method = "spearman")
a_big <- a_pos[a_pos$n_esc >= 50, ]
rho_big <- cor(log(a_big$SER), a_big$eb_mean, method = "spearman")
ser3_h <- cut(a_all$SER, c(-Inf, quantile(a_all$SER, c(1/3, 2/3)), Inf), labels = SER_LEVELS)
eb3_h <- cut(a_all$eb_mean, c(-Inf, quantile(a_all$eb_mean, c(1/3, 2/3)), Inf), labels = SER_LEVELS)
tab_agree <- table(ser3_h, eb3_h)
agree_pct <- 100 * sum(diag(tab_agree)) / sum(tab_agree)
cat("Hospitals with SER = 0 (no observed escalation):", sum(a_all$SER == 0), "\n")

# 分布汇总
summ_dist <- function(x, lab) {
  x <- x[is.finite(x)]
  data.frame(metric = lab,
             n_hospitals = length(x),
             median = median(x), q1 = quantile(x, 0.25), q3 = quantile(x, 0.75),
             min = min(x), max = max(x))
}
dist_ser <- summ_dist(agg$SER, "SER (observed/expected escalation)")
dist_smr <- summ_dist(agg$SMR, "SMR (observed/expected in-hospital death)")
dist_tab <- rbind(dist_ser, dist_smr)
dist_tab$above1 <- c(sum(agg$SER > 1, na.rm = TRUE), sum(agg$SMR > 1, na.rm = TRUE))
dist_tab$ci_above1 <- c(sum(agg$SER_lo > 1, na.rm = TRUE), sum(agg$SMR_lo > 1, na.rm = TRUE))
dist_tab$ci_below1 <- c(sum(agg$SER_hi < 1, na.rm = TRUE), sum(agg$SMR_hi < 1, na.rm = TRUE))
save_csv(dist_tab, "hospital_distribution.csv")

ser_eb_tab <- data.frame(
  metric = c("Spearman rho, log(SER) vs random intercept, all hospitals",
             "Spearman rho, log(SER) vs random intercept, hospitals with >=50 admissions",
             "Tertile agreement, %", "Random-intercept variance", "ICC", "MOR"),
  value = c(round(rho_all, 3), round(rho_big, 3), round(agree_pct, 1),
            round(v_ri, 3), round(icc_ri, 3), round(mor_ri, 2))
)
save_csv(ser_eb_tab, "ser_vs_random_intercept.csv")

agg_out <- agg[, c("hospitalid", "n_esc", "O_esc", "E_esc", "SER", "SER_lo", "SER_hi",
                   "n_death", "O_death", "E_death", "SMR", "SMR_lo", "SMR_hi",
                   "eb_mean", "teaching", "beds", "vol", "dnr_rate")]
save_csv(agg_out, "hospital_level_ser_smr.csv")

cat(sprintf("SER median %.2f (IQR %.2f-%.2f); SMR median %.2f (IQR %.2f-%.2f)\n",
            median(agg$SER, na.rm = TRUE), quantile(agg$SER, .25, na.rm = TRUE),
            quantile(agg$SER, .75, na.rm = TRUE), median(agg$SMR, na.rm = TRUE),
            quantile(agg$SMR, .25, na.rm = TRUE), quantile(agg$SMR, .75, na.rm = TRUE)))
cat(sprintf("SER vs EB: rho all %.3f, n>=50 %.3f, tertile agreement %.1f%%\n",
            rho_all, rho_big, agree_pct))

# ------------------------------------------------------------------
# 4. 医院层面关联：SMR ~ SER（加权线性回归 / 荟萃回归）
# ------------------------------------------------------------------
d <- agg[is.finite(agg$SMR) & agg$E_death > 0, ]
d$beds <- factor(d$beds, levels = BEDS_LEVELS)
d$vol50 <- d$vol / 50
d$teaching <- as.integer(d$teaching)
sd_ser <- sd(d$SER, na.rm = TRUE)

rob_row <- function(m, term, model, r2 = NA) {
  V <- vcovHC(m, type = "HC1")
  b <- coef(m)[term]; se <- sqrt(V[term, term])
  df <- m$df.residual
  p <- 2 * pt(-abs(b / se), df = df)
  crit <- qt(0.975, df)
  data.frame(model = model, term = term, beta = b, se = se,
             lo = b - crit * se, hi = b + crit * se, p = p,
             beta_per_unit_SER = b / sd_ser, r2 = r2, n_hosp = nobs(m))
}

f0 <- fml("SMR", "ser_z")
f1 <- fml("SMR", "ser_z + vol50 + teaching + beds + dnr_rate")
m0 <- lm(f0, data = d, weights = E_death)
m1 <- lm(f1, data = d, weights = E_death)
assoc <- rbind(
  rob_row(m0, "ser_z", "WLS unadjusted (weights = expected deaths)", summary(m0)$r.squared),
  rob_row(m1, "ser_z", "WLS adjusted (weights = expected deaths)", summary(m1)$r.squared)
)

# 敏感性：权重改为观测死亡数；限制医院规模
m0n <- lm(f0, data = d, weights = O_death)
m1n <- lm(f1, data = d, weights = O_death)
assoc <- rbind(assoc,
  rob_row(m0n, "ser_z", "WLS unadjusted (weights = observed deaths)", summary(m0n)$r.squared),
  rob_row(m1n, "ser_z", "WLS adjusted (weights = observed deaths)", summary(m1n)$r.squared))

for (k in c(10, 25)) {
  dk <- d[d$n_death >= k, ]
  mk0 <- lm(f0, data = dk, weights = E_death)
  mk1 <- lm(f1, data = dk, weights = E_death)
  assoc <- rbind(assoc,
    rob_row(mk0, "ser_z", sprintf("WLS unadjusted, hospitals with >=%d deaths", k), summary(mk0)$r.squared),
    rob_row(mk1, "ser_z", sprintf("WLS adjusted, hospitals with >=%d deaths", k), summary(mk1)$r.squared))
}

# 荟萃回归：log(SMR)，Poisson 方差
d$yi <- log((d$O_death + 0.5) / d$E_death)
d$vi <- 1 / (d$O_death + 0.5)
mr0 <- rma(yi, vi, mods = as.formula(paste(TILDE, "ser_z")), data = d, method = "REML")
mr1 <- rma(yi, vi, mods = as.formula(paste(TILDE, "ser_z + vol50 + teaching + beds + dnr_rate")),
           data = d, method = "REML")
meta_row <- function(mr, model) {
  idx <- match("ser_z", rownames(mr$beta))
  b <- mr$beta[idx, 1]; se <- mr$se[idx]
  data.frame(model = model, term = "ser_z", beta = b, se = se,
             lo = mr$ci.lb[idx], hi = mr$ci.ub[idx],
             p = 2 * pnorm(-abs(b / se)), beta_per_unit_SER = b / sd_ser,
             r2 = if (is.null(mr$R2)) NA_real_ else mr$R2, n_hosp = mr$k)
}
assoc <- rbind(assoc,
  meta_row(mr0, "Random-effects meta-regression, unadjusted"),
  meta_row(mr1, "Random-effects meta-regression, adjusted"))
assoc$model <- as.character(assoc$model)
save_csv(assoc, "hospital_association.csv")

# 散点图数据
save_csv(d[, c("hospitalid", "SER", "SMR", "n_esc", "n_death", "E_death", "ser_z")],
         "hospital_scatter_data.csv")

cat("\nHospital-level association (SMR ~ SER):\n")
print(assoc[, c("model", "beta", "lo", "hi", "p", "r2", "n_hosp")])

# ------------------------------------------------------------------
# 5. 病人层面效应修饰：多水平 logistic（院内死亡）
# ------------------------------------------------------------------
qs_ser <- quantile(agg_e$SER, c(1/3, 2/3))
agg_e$ser3 <- cut(agg_e$SER, c(-Inf, qs_ser, Inf), labels = SER_LEVELS)
qs_ser2 <- median(agg_e$SER)
agg_e$ser2 <- ifelse(agg_e$SER > qs_ser2, "High", "Low")
qs_ser4 <- quantile(agg_e$SER, c(1/4, 2/4, 3/4))
agg_e$ser4 <- cut(agg_e$SER, c(-Inf, qs_ser4, Inf), labels = c("Q1", "Q2", "Q3", "Q4"))

qa <- quantile(ei$apachescore, c(1/3, 2/3), na.rm = TRUE)

prep_patient <- function(d) {
  d$ser3 <- agg_e$ser3[match(d$hospitalid, agg_e$hospitalid)]
  d$ser2 <- agg_e$ser2[match(d$hospitalid, agg_e$hospitalid)]
  d$ser4 <- agg_e$ser4[match(d$hospitalid, agg_e$hospitalid)]
  d$sev3 <- cut(d$apachescore, c(-Inf, qa, Inf), labels = c("S1", "S2", "S3"), include.lowest = TRUE)
  d$gcs_cat <- cut(d$gcs_total, c(-Inf, 8, 12, 15), labels = c("G1", "G2", "G3"), include.lowest = TRUE)
  d$ser3 <- factor(d$ser3, levels = SER_LEVELS)
  d$sev3 <- factor(d$sev3, levels = c("S1", "S2", "S3"))
  d$gcs_cat <- factor(d$gcs_cat, levels = c("G1", "G2", "G3"))
  d$ser_z <- as.numeric(scale(agg_e$SER))[match(d$hospitalid, agg_e$hospitalid)]
  d$ok <- !is.na(d$death_hosp) & d$dnr_known == 1
  d
}
clp <- lapply(cl, prep_patient)

fit_pool_glmer <- function(form, dat_list, terms = NULL, joint = NULL) {
  fits <- lapply(dat_list, function(d) {
    dd <- d[d$ok, , drop = FALSE]
    for (v in c("dx_group", "ser3", "sev3", "gcs_cat", "cutv")) {
      if (!is.null(dd[[v]])) dd[[v]] <- droplevels(dd[[v]])
    }
    tryCatch(glmer(form, data = dd, family = binomial, nAGQ = 0,
                   control = glmerControl(optimizer = "bobyqa",
                                          optCtrl = list(maxfun = 2e5),
                                          calc.derivs = FALSE)),
             error = function(e) NULL)
  })
  okf <- !vapply(fits, is.null, TRUE)
  fits <- fits[okf]
  out <- list(fits = fits, converged = sum(okf), total = length(okf))
  if (!is.null(terms)) out$terms <- do.call(rbind, lapply(terms, function(t) pool_scalar(fits, t)))
  if (!is.null(joint)) out$joint <- pool_joint(fits, joint)
  out
}

# 5.1 SER 三分位 × 严重度三分位
f_sev <- fml("death_hosp", "ser3 * sev3 + age10 + female + apache10 + gcs3 + sbp10 + dx_group + dnr0 + (1|hosp)")
res_sev <- fit_pool_glmer(f_sev, clp,
  terms = c("ser3Mid", "ser3High"),
  joint = c("ser3Mid:sev3S2", "ser3High:sev3S2", "ser3Mid:sev3S3", "ser3High:sev3S3"))
cat("Severity interaction model converged:", res_sev$converged, "/", res_sev$total, "\n")

strat_sev <- do.call(rbind, lapply(c("S1", "S2", "S3"), function(s) {
  L <- c(ser3High = 1)
  if (s == "S2") L <- c(L, "ser3High:sev3S2" = 1)
  if (s == "S3") L <- c(L, "ser3High:sev3S3" = 1)
  z <- pool_lincom(res_sev$fits, L)
  data.frame(stratum = s, est = z$est, se = z$se, p = z$p, lo = z$lo, hi = z$hi,
             OR = exp(z$est), OR_lo = exp(z$lo), OR_hi = exp(z$hi))
}))
strat_sev$stratum_label <- c("APACHE tertile 1 (low)", "APACHE tertile 2", "APACHE tertile 3 (high)")[match(strat_sev$stratum, c("S1", "S2", "S3"))]

# 5.2 SER 三分位 × GCS 分类
f_gcs <- fml("death_hosp", "ser3 * gcs_cat + age10 + female + apache10 + gcs3 + sbp10 + dx_group + dnr0 + (1|hosp)")
res_gcs <- fit_pool_glmer(f_gcs, clp,
  terms = c("ser3Mid", "ser3High"),
  joint = c("ser3Mid:gcs_catG2", "ser3High:gcs_catG2", "ser3Mid:gcs_catG3", "ser3High:gcs_catG3"))
cat("GCS interaction model converged:", res_gcs$converged, "/", res_gcs$total, "\n")

strat_gcs <- do.call(rbind, lapply(c("G1", "G2", "G3"), function(s) {
  L <- c(ser3High = 1)
  if (s == "G2") L <- c(L, "ser3High:gcs_catG2" = 1)
  if (s == "G3") L <- c(L, "ser3High:gcs_catG3" = 1)
  z <- pool_lincom(res_gcs$fits, L)
  data.frame(stratum = s, est = z$est, se = z$se, p = z$p, lo = z$lo, hi = z$hi,
             OR = exp(z$est), OR_lo = exp(z$lo), OR_hi = exp(z$hi))
}))
strat_gcs$stratum_label <- c("GCS <=8", "GCS 9-12", "GCS 13-15")[match(strat_gcs$stratum, c("G1", "G2", "G3"))]

# 5.3 连续交互（SER 每 1 SD × 严重度每 10 分 / GCS 每 3 分）
f_cont <- fml("death_hosp", "ser_z * apache10 + ser_z * gcs3 + age10 + female + apache10 + gcs3 + sbp10 + dx_group + dnr0 + (1|hosp)")
res_cont <- fit_pool_glmer(f_cont, clp, terms = c("ser_z", "ser_z:apache10", "ser_z:gcs3"))
cat("Continuous interaction model converged:", res_cont$converged, "/", res_cont$total, "\n")

# 汇总表
eff_tab <- rbind(
  data.frame(model = "SER tertile x APACHE tertile", term = "Interaction (4 df)",
             est = NA, se = NA, p = res_sev$joint$p, OR = NA, OR_lo = NA, OR_hi = NA),
  data.frame(model = "SER tertile x GCS category", term = "Interaction (4 df)",
             est = NA, se = NA, p = res_gcs$joint$p, OR = NA, OR_lo = NA, OR_hi = NA),
  data.frame(model = "Continuous SER x APACHE", term = "SER per SD x APACHE per 10",
             est = res_cont$terms$est[2], se = res_cont$terms$se[2], p = res_cont$terms$p[2],
             OR = exp(res_cont$terms$est[2]), OR_lo = exp(res_cont$terms$lo[2]), OR_hi = exp(res_cont$terms$hi[2])),
  data.frame(model = "Continuous SER x GCS", term = "SER per SD x GCS per 3",
             est = res_cont$terms$est[3], se = res_cont$terms$se[3], p = res_cont$terms$p[3],
             OR = exp(res_cont$terms$est[3]), OR_lo = exp(res_cont$terms$lo[3]), OR_hi = exp(res_cont$terms$hi[3]))
)
save_csv(eff_tab, "patient_effect_modification.csv")
strat_all <- rbind(
  data.frame(model = "APACHE tertile", strat_sev[, c("stratum_label", "OR", "OR_lo", "OR_hi", "p")]),
  data.frame(model = "GCS category", strat_gcs[, c("stratum_label", "OR", "OR_lo", "OR_hi", "p")])
)
names(strat_all)[2] <- "stratum"
save_csv(strat_all, "patient_stratified_or.csv")

cat("\nInteraction tests:\n"); print(eff_tab)
cat("\nStratified ORs (high vs low SER):\n"); print(strat_all)

# ------------------------------------------------------------------
# 6. 敏感性：排除 DNR、仅内科 ICU、排除 OR/PACU、不同 SER 切点
# ------------------------------------------------------------------
run_subset <- function(label, keep) {
  keep <- as.logical(keep)
  # 升级预测与 SER（子集内重新拟合）
  p_esc2 <- rep(NA_real_, n)
  pm <- matrix(NA_real_, sum(keep), length(cl))
  for (i in seq_along(cl)) {
    d <- cl[[i]][keep, , drop = FALSE]
    d$dx_group <- droplevels(d$dx_group)
    fit <- glm(f_esc, data = d, family = binomial)
    pm[, i] <- predict(fit, newdata = d, type = "response")
  }
  p_esc2[keep] <- rowMeans(pm, na.rm = TRUE)
  pe <- data.frame(hospitalid = cl[[1]]$hospitalid[keep], y = y_esc[keep], p = p_esc2[keep])
  ae <- aggregate(pe[, c("y", "p")], by = list(hospitalid = pe$hospitalid), FUN = sum)
  names(ae) <- c("hospitalid", "O", "E"); ae$SER <- ae$O / ae$E
  ae$ser3 <- cut(ae$SER, c(-Inf, quantile(ae$SER, c(1/3, 2/3)), Inf), labels = SER_LEVELS)

  # 死亡预测与 SMR（子集内重新拟合）
  ok2 <- keep & !is.na(cl[[1]]$death_hosp) & cl[[1]]$dnr_known == 1
  p_death2 <- rep(NA_real_, n)
  pmd <- matrix(NA_real_, sum(ok2), length(cl))
  for (i in seq_along(cl)) {
    d <- cl[[i]][ok2, , drop = FALSE]
    d$dx_group <- droplevels(d$dx_group)
    fit <- glm(f_death, data = d, family = binomial)
    pmd[, i] <- predict(fit, newdata = d, type = "response")
  }
  p_death2[ok2] <- rowMeans(pmd, na.rm = TRUE)
  pd2 <- data.frame(hospitalid = cl[[1]]$hospitalid[ok2], y = cl[[1]]$death_hosp[ok2], p = p_death2[ok2])
  ad <- aggregate(pd2[, c("y", "p")], by = list(hospitalid = pd2$hospitalid), FUN = sum)
  names(ad) <- c("hospitalid", "O", "E"); ad$SMR <- ad$O / ad$E

  ag <- merge(ae, ad, by = "hospitalid", all.x = TRUE)
  ag <- merge(ag, L$hosp, by = "hospitalid", all.x = TRUE)
  ag$vol50 <- ag$vol / 50; ag$teaching <- as.integer(ag$teaching)
  ag$beds <- factor(ag$beds, levels = BEDS_LEVELS)
  ag$ser_z <- as.numeric(scale(ag$SER))
  dg <- ag[is.finite(ag$SMR) & ag$E.y > 0, ]
  if (nrow(dg) < 20) return(NULL)
  m <- lm(fml("SMR", "ser_z + vol50 + teaching + beds + dnr_rate"), data = dg, weights = E.y)
  b <- coef(m)["ser_z"]; se <- sqrt(vcovHC(m, type = "HC1")["ser_z", "ser_z"])
  p <- 2 * pt(-abs(b / se), df = m$df.residual)

  # 病人层面交互（子集内）
  clp2 <- lapply(cl, function(d) {
    d$ser3 <- ae$ser3[match(d$hospitalid, ae$hospitalid)]
    d$sev3 <- cut(d$apachescore, c(-Inf, qa, Inf), labels = c("S1", "S2", "S3"), include.lowest = TRUE)
    d$ser3 <- factor(d$ser3, levels = SER_LEVELS); d$sev3 <- factor(d$sev3, levels = c("S1", "S2", "S3"))
    d$ok <- keep & !is.na(d$death_hosp) & d$dnr_known == 1
    d
  })
  rs <- fit_pool_glmer(f_sev, clp2, joint = c("ser3Mid:sev3S2", "ser3High:sev3S2",
                                               "ser3Mid:sev3S3", "ser3High:sev3S3"))
  data.frame(subset = label, n_hosp = nrow(dg), beta_SMR_SER = b, se = se,
             lo = b - qt(0.975, m$df.residual) * se, hi = b + qt(0.975, m$df.residual) * se,
             p_hospital = p, p_interaction = rs$joint$p, n_glmer = rs$converged)
}

sens <- rbind(
  run_subset("Excluding pre-NIV DNR", cl[[1]]$dnr_known == 1 & cl[[1]]$dnr0 == 0),
  run_subset("Medical ICU only", cl[[1]]$icu_medical == 1),
  run_subset("Excluding OR/PACU source", cl[[1]]$or_pacu == 0)
)
save_csv(sens, "sensitivity_hospital_interaction.csv")
cat("\nSensitivity analyses:\n"); print(sens)

# 不同 SER 切点（全队列）：中位数、四分位
run_cut <- function(cutvar, label) {
  clp2 <- lapply(cl, function(d) {
    d$cutv <- factor(agg_e[[cutvar]][match(d$hospitalid, agg_e$hospitalid)])
    d$sev3 <- cut(d$apachescore, c(-Inf, qa, Inf), labels = c("S1", "S2", "S3"), include.lowest = TRUE)
    d$sev3 <- factor(d$sev3, levels = c("S1", "S2", "S3"))
    d$ok <- !is.na(d$death_hosp) & d$dnr_known == 1
    d
  })
  form <- fml("death_hosp", "cutv * sev3 + age10 + female + apache10 + gcs3 + sbp10 + dx_group + dnr0 + (1|hosp)")
  fits <- lapply(clp2, function(d) {
    dd <- d[d$ok, , drop = FALSE]
    dd$dx_group <- droplevels(dd$dx_group); dd$cutv <- droplevels(dd$cutv); dd$sev3 <- droplevels(dd$sev3)
    tryCatch(glmer(form, data = dd, family = binomial, nAGQ = 0,
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5),
                                          calc.derivs = FALSE)), error = function(e) NULL)
  })
  fits <- fits[!vapply(fits, is.null, TRUE)]
  trm <- grep(":sev3", names(get_fixef(fits[[1]])), value = TRUE)
  j <- pool_joint(fits, trm)
  data.frame(definition = label, n_levels = nlevels(droplevels(clp2[[1]]$cutv)),
             interaction_terms = length(trm), p_interaction = j$p)
}
cut_tab <- rbind(
  run_cut("ser2", "SER above vs below median"),
  run_cut("ser3", "SER tertiles"),
  run_cut("ser4", "SER quartiles")
)
save_csv(cut_tab, "sensitivity_ser_cutpoints.csv")
cat("\nSER cutpoint sensitivity:\n"); print(cut_tab)

# ------------------------------------------------------------------
# 7. 图
# ------------------------------------------------------------------
theme_mc <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        plot.title = element_blank())

# 图 A/B：SER-随机截距一致性；SMR-SER 关联
pa <- ggplot(a_pos, aes(log(SER), eb_mean)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_point(aes(size = n_esc), alpha = 0.65, colour = "#2C6E8F") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#B85042", linewidth = 0.8) +
  scale_size_continuous(range = c(0.6, 4), name = "Admissions") +
  labs(x = "log observed/expected escalation ratio (SER)",
       y = "Random-intercept estimate (log-odds)") +
  annotate("text", x = min(log(a_pos$SER)), y = max(a_pos$eb_mean),
           hjust = 0, vjust = 1, size = 3.2,
           label = sprintf("Spearman rho = %.2f", rho_all)) +
  theme_mc

m0_line <- lm(SMR ~ SER, data = d, weights = E_death)
xr <- range(d$SER[is.finite(d$SER)])
xl <- seq(xr[1], quantile(d$SER, 0.99, na.rm = TRUE), length.out = 100)
yl <- predict(m0_line, newdata = data.frame(SER = xl))
line_df <- data.frame(SER = xl, SMR = yl)
pb <- ggplot(d, aes(SER, SMR)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey60") +
  geom_point(aes(size = E_death), alpha = 0.65, colour = "#2C6E8F") +
  geom_line(data = line_df, colour = "#B85042", linewidth = 0.8) +
  scale_size_continuous(range = c(0.6, 4), name = "Expected deaths") +
  coord_cartesian(xlim = c(0, quantile(d$SER, 0.99, na.rm = TRUE)),
                  ylim = c(0, quantile(d$SMR, 0.99, na.rm = TRUE))) +
  labs(x = "Observed/expected escalation ratio (SER)",
       y = "Observed/expected in-hospital death (SMR)") +
  theme_mc

fig1 <- (pa | pb) + plot_annotation(tag_levels = "A")
ggsave(file.path(OUT, "figures", "fig_hospital_ser_smr.png"), fig1,
       width = 9, height = 3.8, dpi = 300)
ggsave(file.path(OUT, "figures", "fig_hospital_ser_smr.pdf"), fig1,
       width = 9, height = 3.8)

# 图 C/D：分层 OR 森林图
strat_plot <- function(df, p_int, xlab) {
  df$stratum <- factor(df$stratum_label, levels = rev(df$stratum_label))
  ggplot(df, aes(OR, stratum)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
    geom_pointrange(aes(xmin = OR_lo, xmax = OR_hi), colour = "#2C6E8F", size = 0.5) +
    scale_x_log10() +
    labs(x = xlab, y = NULL) +
    annotate("text", x = min(df$OR_lo) * 0.95, y = 0.5, hjust = 0, size = 3.1,
             label = sprintf("P for interaction = %.3f", p_int)) +
    theme_mc
}
pc <- strat_plot(strat_sev, res_sev$joint$p,
                 "OR for death, high vs low SER (log scale)")
pd_plot <- strat_plot(strat_gcs, res_gcs$joint$p,
                      "OR for death, high vs low SER (log scale)")
fig2 <- (pc | pd_plot) + plot_annotation(tag_levels = "A")
ggsave(file.path(OUT, "figures", "fig_patient_interaction.png"), fig2,
       width = 9.5, height = 3.2, dpi = 300)
ggsave(file.path(OUT, "figures", "fig_patient_interaction.pdf"), fig2,
       width = 9.5, height = 3.2)

# 保存关键对象
strip_fits <- function(x) { x$fits <- NULL; x }
saveRDS(list(agg = agg, assoc = assoc,
             res_sev = strip_fits(res_sev), res_gcs = strip_fits(res_gcs),
             res_cont = strip_fits(res_cont), strat_sev = strat_sev, strat_gcs = strat_gcs,
             sens = sens, cut_tab = cut_tab, ser_eb_tab = ser_eb_tab,
             rho_all = rho_all, rho_big = rho_big, agree_pct = agree_pct,
             v_ri = v_ri, icc_ri = icc_ri, mor_ri = mor_ri,
             p_esc = p_esc, p_death = p_death, ok_death = ok_death),
        file.path(OUT, "_data", "hospital_outcomes.rds"))
cat("\n01 done.\n")
