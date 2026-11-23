# Reproducible Analytical Pipeline: Noninvasive Ventilation (NIV) Escalation in Critical Care

**Official Code Repository:** [https://github.com/fjy2026/ICU-NIV](https://github.com/fjy2026/ICU-NIV)  
**Target Journal:** *JMIR Medical Informatics* (JMIR Publications Inc.)  
**Manuscript Title:** Escalation From Noninvasive to Invasive Ventilation in Intensive Care: A Multicentre Cohort Study of Patient Factors, Between-Hospital Variation, and Transportability  
**Authors:** Hangyu He, MD; Wei Li, MD; Tian Qiu, MD; Shuang Li, MD, PhD; Jiayu Fang, MD  

---

## 1. Scientific Overview & Study Aims

Noninvasive ventilation (NIV) is a cornerstone of acute respiratory support in intensive care units (ICUs). However, NIV failure requiring escalation to invasive mechanical ventilation (IMV) occurs in 10%–50% of patients and is associated with substantial peri-intubation morbidity and excess mortality. While clinical decision support systems (CDSS) and machine-learning risk prediction architectures have been proposed, existing algorithms frequently suffer from high risk of bias, lack of multicentre external validation, and catastrophic performance degradation when transported across diverse hospital environments.

This repository provides the core data science, electronic phenotyping, and statistical modeling codebase supporting our multicentre investigation across two major clinical EHR databases:
1. **Derivation Cohort (eICU-CRD v2.0):** 13,986 adult ICU admissions treated with first-line NIV across 186 diverse United States hospitals (1,345 48-hour escalations to IMV, 9.6%);
2. **External Validation Cohort (MIMIC-IV v3.1):** 1,983 adult ICU admissions from a major academic medical center (Beth Israel Deaconess Medical Center, Boston, MA; 459 48-hour escalations to IMV, 23.1%).

### Key Methodological Findings
- **Multicentre Association Structure:** Exchangeable correlation Generalized Estimating Equations (GEE) with robust sandwich variance across 20 multiple-imputed datasets (Rubin's rules) to identify robust physiological and clinical predictors of 48-hour escalation;
- **Quantification of Institutional Practice Variation:** Multilevel mixed-effects logistic regression models revealing marked institutional variation (intraclass correlation coefficient [ICC] 0.169; median odds ratio [MOR] 2.18), demonstrating that institutional escalation propensity acts as a profound algorithmic confounder;
- **Empirical Transportability Evaluation:** Formal testing of cross-system transportability demonstrating severe calibration drift (calibration-in-the-large 2.511) and feature discordance when deploying algorithms across disparate health systems;
- **Temporal Horizon & Intubation Trigger Bias:** Sensitivity analyses examining pre-NIV physiological windows (0–6 h vs 6–24 h) demonstrating that extreme vital signs capture proximate dynamic intubation triggers rather than static baseline risk;
- **Feature Collinearity & Decomposition:** Resolution of negative suppressor confounding between composite scores (SOFA) and individual physiological features (GCS) through non-neurological feature decomposition;
- **Information Sparsity in Passive Surveillance:** Audit of clinical scores (HACOR / ROX) revealing that mandatory arterial blood gas components were absent in 77.3% of multicentre admissions, precluding autonomous passive electronic surveillance.

---

## 2. Standardized Code Repository Architecture

The `code/` directory is strictly focused on pure data science, cohort phenotyping, and statistical modeling, structured into six sequential pipeline stages alongside shared environment utilities:

```
code/
├── 00_common/                           # Shared statistical environment, scaling & modeling formulas
│   ├── common_env.R                     # Core packages, scaling rules, random seeds (20261113), diagnostic tiers
│   ├── common_transportability.R        # External validation, calibration curves & transportability helpers
│   ├── common_multilevel.R              # Multilevel mixed-effects models, ICC, and MOR calculation utilities
│   ├── common_temporal.R                # Pre-NIV physiological window extraction & temporal alignment
│   ├── trace_logger.R                   # R execution tracing and analytical step audit logger
│   └── trace_logger.py                  # Python execution tracing and data extraction logger
│
├── 01_cohort_extraction/                # 1. 人群纳排与特征提取 (Cohort Extraction & Phenotyping)
│   ├── 01_extract_eicu_cohort.R         # eICU-CRD NIV cohort extraction, ventilation timing & inclusion/exclusion
│   ├── 02_extract_mimic_cohort.R        # MIMIC-IV NIV cohort extraction, procedure events & mechanical ventilation
│   ├── 03_extract_eicu_expanded_features.R # Expanded eICU feature extraction (arterial blood gases, GCS subscores)
│   ├── 04_extract_mimic_expanded_features.R # Expanded MIMIC feature extraction (SOFA components, vasopressors)
│   ├── 05_build_extraction_metadata.py  # Automated feature extraction metadata dictionary compiler
│   ├── 06_generate_concept_codebooks.R  # EHR concept dictionary and variable definitions generator
│   ├── 07_generate_expanded_codebooks.R # Codebook compiler for supplemental clinical covariates
│   └── 08_update_extraction_metadata.R  # Cohort flow counts and extraction audit synchronizer
│
├── 02_preprocessing_imputation/         # 2. 数据清洗与多重插补 (Preprocessing & Imputation)
│   ├── 01_data_preprocessing_scaling.R  # Scaling (age/10, SBP/10, RR/5, SpO2/5), missingness flags & baseline prep
│   ├── 02_multiple_imputation_mice.R    # 20-fold multiple imputation via chained equations (MICE, maxit=10)
│   └── 03_baseline_characteristics_table.R # Standardized baseline comparisons and SMD balance diagnostics
│
├── 03_primary_models/                   # 3. GEE 主模型与绝对风险估计 (Primary GEE Models)
│   ├── 01_primary_multivariable_gee.R   # Primary exchangeable GEE logistic model with Rubin's rules pooling
│   ├── 02_gee_working_correlation_comparison.R # Working correlation comparison (Exchangeable vs Independence HC0)
│   └── 03_absolute_risk_and_evalues.R   # Baseline absolute risk differences and E-value sensitivity
│
├── 04_hospital_variation/               # 4. 多水平随机截距模型与中心变异 (Hospital Variation)
│   └── 01_multilevel_random_intercept_hospital.R # Multilevel mixed-effects model for ICC (0.169) and MOR (2.18)
│
├── 05_external_validation/              # 5. 跨系统外部验证与可迁移性 (External Validation)
│   ├── 01_external_validation_mimic.R   # External validation cohort evaluation in MIMIC-IV
│   ├── 02_calibration_curves_assessment.R # Logistic and Loess calibration curves (intercept & calibration slope)
│   ├── 03_cross_system_transportability.R # Transportability diagnostics across 7 alignable clinical features
│   └── 04_decision_curve_analysis.R    # Decision Curve Analysis (DCA) for net benefit across threshold probabilities
│
├── 06_subgroups_sensitivity/            # 6. 亚组分析、动态评分与敏感性分析 (Subgroups & Sensitivity)
│   ├── 01_sensitivity_dnr_exclusion.R   # Sensitivity analysis excluding Do-Not-Resuscitate / palliative admissions
│   ├── 02_secondary_clinical_outcomes.R # Secondary clinical outcomes (ICU mortality, hospital LOS, 28-day mortality)
│   ├── 03_subgroup_analyses_and_hacor.R # Subgroup stratification (diagnostic groups, age tiers, severity strata)
│   ├── 04_multivariable_sensitivity_suite.R # Complete-case analysis, extreme outlier winsorization, alternate cutoffs
│   ├── 05_arterial_bloodgas_subgroup.R  # Blood-gas complete-case subgroup model (PaO2/FiO2, pH, PaCO2)
│   ├── 06_hacor_reclassification_nri_idi.R # Net Reclassification Improvement (NRI) and Integrated Discrimination (IDI)
│   ├── 07_time_zero_pre_icu_niv_sensitivity.R # Time-zero sensitivity: pre-ICU NIV vs ICU-initiated NIV stratification
│   ├── 08_cross_system_feature_alignment.R # Cross-system variable alignment and semantic harmonization checks
│   ├── 09_hospital_outcome_correlations.R # Hospital-level outcome correlation with institutional escalation rate
│   ├── 10_modifiable_hospital_factors.R # Modifiable practice factors (nurse-to-patient ratio, unit volume)
│   ├── 11_temporal_pre_niv_vitals.R     # Pre-NIV physiological window (0-6h prior to initiation) analysis
│   ├── 12_temporal_landmark_24h.R       # 24-hour landmark analysis excluding early (<24h) post-initiation intubations
│   ├── 13_temporal_vasopressor_stability.R # Concomitant vasoactive infusions and cardiovascular stability checks
│   └── 14_temporal_sensitivity_summary.R # Consolidated temporal sensitivity comparison matrix
│
└── README.md                            # Master reproducibility and pipeline execution documentation
```

---

## 3. Step-by-Step Analytical Execution Pipeline

To replicate the statistical findings from raw EHR extracts, run the pipeline sequentially from the project root:

### Stage 1: Cohort Extraction and Phenotyping
Ensure raw PhysioNet eICU-CRD and MIMIC-IV database tables are accessible. Run:
```bash
Rscript code/01_cohort_extraction/01_extract_eicu_cohort.R
Rscript code/01_cohort_extraction/02_extract_mimic_cohort.R
Rscript code/01_cohort_extraction/03_extract_eicu_expanded_features.R
Rscript code/01_cohort_extraction/04_extract_mimic_expanded_features.R
python code/01_cohort_extraction/05_build_extraction_metadata.py
```
*Key Deliverables:* Preprocessed analytical cohort Parquet files and variable extraction dictionaries.

### Stage 2: Data Preprocessing and Multiple Imputation
Standardize continuous variables, flag missing values, and execute 20-fold chained equations imputation:
```bash
Rscript code/02_preprocessing_imputation/01_data_preprocessing_scaling.R
Rscript code/02_preprocessing_imputation/02_multiple_imputation_mice.R
Rscript code/02_preprocessing_imputation/03_baseline_characteristics_table.R
```
*Key Deliverables:* Analysis-ready prepared datasets and 20 completed imputed datasets (`m = 20`, seed = 20261113).

### Stage 3: Primary Multivariable GEE Modeling
Fit the primary hospital-clustered GEE models and evaluate unmeasured confounding:
```bash
Rscript code/03_primary_models/01_primary_multivariable_gee.R
Rscript code/03_primary_models/02_gee_working_correlation_comparison.R
Rscript code/03_primary_models/03_absolute_risk_and_evalues.R
```
*Key Deliverables:* Pooled multivariable odds ratios (Rubin's rules), sandwich robust standard errors, and E-values.

### Stage 4: Multilevel Modeling & Hospital Practice Variation
Fit the mixed-effects random-intercept models to quantify between-hospital variation:
```bash
Rscript code/04_hospital_variation/01_multilevel_random_intercept_hospital.R
```
*Key Deliverables:* Intraclass correlation coefficient (ICC = 0.169), median odds ratio (MOR = 2.18), and empirical Bayes hospital residuals.

### Stage 5: External Validation & Cross-System Transportability
Evaluate algorithm performance, calibration, and transportability in the independent MIMIC-IV cohort:
```bash
Rscript code/05_external_validation/01_external_validation_mimic.R
Rscript code/05_external_validation/02_calibration_curves_assessment.R
Rscript code/05_external_validation/03_cross_system_transportability.R
Rscript code/05_external_validation/04_decision_curve_analysis.R
```
*Key Deliverables:* Calibration-in-the-large (2.511), calibration slope, 7-feature transportability metrics, and decision curves.

### Stage 6: Subgroup Stratifications & Sensitivity Suite
Execute comprehensive sensitivity analyses, temporal horizon evaluations, and clinical score audits:
```bash
Rscript code/06_subgroups_sensitivity/01_sensitivity_dnr_exclusion.R
Rscript code/06_subgroups_sensitivity/02_secondary_clinical_outcomes.R
Rscript code/06_subgroups_sensitivity/03_subgroup_analyses_and_hacor.R
Rscript code/06_subgroups_sensitivity/04_multivariable_sensitivity_suite.R
Rscript code/06_subgroups_sensitivity/05_arterial_bloodgas_subgroup.R
Rscript code/06_subgroups_sensitivity/06_hacor_reclassification_nri_idi.R
Rscript code/06_subgroups_sensitivity/07_time_zero_pre_icu_niv_sensitivity.R
Rscript code/06_subgroups_sensitivity/11_temporal_pre_niv_vitals.R
Rscript code/06_subgroups_sensitivity/12_temporal_landmark_24h.R
```
*Key Deliverables:* Complete sensitivity matrices confirming the robustness of primary clinical conclusions across operational definitions.

---

## 4. Software Dependencies & Runtime Specifications

### R Environment
- **R Version:** >= 4.3.0 (Tested on R 4.4.1 x86_64)
- **Core Packages:**
  - `duckdb` (>= 0.10.0) — High-performance embedded database engine for Parquet reads
  - `mice` (>= 3.16.0) — Multivariate Imputation by Chained Equations
  - `geepack` (>= 1.3.9) — Generalized Estimating Equations with user-defined working correlations
  - `sandwich` (>= 3.1.0) — Model-robust covariance matrix estimators (HC0)
  - `lme4` (>= 1.1-35) — Mixed-effects logistic regression for multilevel hospital modeling
  - `rms` (>= 6.8-0) — Regression modeling strategies & restricted cubic splines
  - `pROC` (>= 1.18.5) — ROC analysis and DeLong non-parametric comparisons
  - `EValue` (>= 4.1.3) — Sensitivity analysis for unmeasured confounding
  - `metafor` (>= 4.6-0) — Random-effects meta-analytic pooling
  - `data.table` (>= 1.15.0), `dplyr` (>= 1.1.4), `tidyr` (>= 1.3.1), `readr` (>= 2.1.5)
- **Installation:**
  ```r
  install.packages(c(
    "duckdb", "mice", "geepack", "sandwich", "lme4", "rms", "pROC",
    "EValue", "metafor", "data.table", "dplyr", "tidyr", "readr"
  ))
  ```

### Python Environment
- **Python Version:** >= 3.10 (Tested on Python 3.12.3)
- **Required Libraries:**
  - `duckdb` (>= 0.10.0), `pandas` (>= 2.2.0), `pyarrow` (>= 15.0.0)
- **Installation:**
  ```bash
  pip install duckdb pandas pyarrow
  ```

---

## 5. Data Access & Governance Protocols

The raw clinical datasets analyzed in this study are third-party electronic health record databases hosted on PhysioNet:
1. **eICU Collaborative Research Database (eICU-CRD v2.0):**  
   Pollard TJ, Johnson AEW, Raffa JD, Celi LA, Mark RG, Badawi O. The eICU Collaborative Research Database, a freely available multi-center database for critical care research. *Sci Data* 2018;5:180178.  
   DOI: [10.1038/sdata.2018.178](https://doi.org/10.1038/sdata.2018.178) | PhysioNet: [https://physionet.org/content/eicudb/2.0/](https://physionet.org/content/eicudb/2.0/)
2. **Medical Information Mart for Intensive Care IV (MIMIC-IV v3.1):**  
   Johnson AEW, Bulgarelli L, Shen L, Gayles A, Shammout A, Horng S, Pollard TJ, Hao S, Moody B, Gow B, Lehman LH, Celi LA, Mark RG. MIMIC-IV, a freely accessible electronic health record dataset. *Sci Data* 2023;10(1):1.  
   DOI: [10.1038/s41597-022-01899-x](https://doi.org/10.1038/s41597-022-01899-x) | PhysioNet: [https://physionet.org/content/mimiciv/3.1/](https://physionet.org/content/mimiciv/3.1/)

### Data Governance Requirements
Under the PhysioNet Credentialed Data Use Agreement, individual-level health records cannot be redistributed in public repositories. To execute the cohort extraction from raw records:
1. Complete CITI Program training in *Human Subjects Research* and *Health Privacy (HIPAA)*;
2. Execute the formal PhysioNet Data Use Agreement for eICU-CRD and MIMIC-IV;
3. Place raw database tables in the designated project directory;
4. Execute `01_cohort_extraction/` scripts to generate analysis-ready cohorts.

---

## 6. Reproducibility & Seed Protocol

To ensure exact numerical reproducibility across stochastic modeling steps:
- **Random Seeds:**
  - Primary analysis and multiple imputation: `SEED <- 20261113`
  - Supplemental external and temporal analyses: `SEED <- 20261117`
  - Bootstrap resampling (1,000 BCa iterations): Pre-set seeds embedded in respective scripts
- **Imputation Stability:**
  - Multiple imputation uses chained equations (`mice`) with `m = 20` completed datasets and `maxit = 10` iterations;
  - All multivariable estimates and odds ratios are pooled across the 20 imputed datasets using Rubin's rules.

---

## 7. Ethics Approvals & Research Compliance

- **eICU-CRD Approval:** Approved by the Institutional Review Board (IRB) of the Massachusetts Institute of Technology (protocol 0403000206).
- **MIMIC-IV Approval:** Approved by the IRBs of Beth Israel Deaconess Medical Center and MIT (protocol 2001P001699).
- **Informed Consent Waiver:** Both databases comprise certified, fully de-identified records adhering to the HIPAA Expert Determination method; individual informed consent was formally waived.
- **Reporting Fidelity:** Conducted in strict accordance with STROBE (Strengthening the Reporting of Observational Studies in Epidemiology) and TRIPOD (Transparent Reporting of a multivariable prediction model for Individual Prognosis Or Diagnosis) guidelines.

---

## 8. License & Inquiries

- **Code License:** Released under the [MIT License](https://opensource.org/licenses/MIT). You are free to inspect, adapt, and reuse this analytical code for academic research with appropriate attribution.
- **Academic Inquiries:**  
  - **Jiayu Fang, MD** (Co-corresponding author), Department of Otorhinolaryngology–Head & Neck Surgery, Zhongnan Hospital of Wuhan University. Email: `366619601@qq.com`  
  - **Shuang Li, MD, PhD** (Co-corresponding author), Department of Otorhinolaryngology–Head & Neck Surgery, Zhongnan Hospital of Wuhan University. Email: `entlishuang@whu.edu.cn`
