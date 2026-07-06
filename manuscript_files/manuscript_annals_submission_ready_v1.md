# Hemodynamic Recovery Capacity After Early Circulatory Support Across Three Critical Care Databases

Target journal: Annals of Intensive Care

Article type: Research Article

Submission status: Condensed submission-ready draft v1. Author-specific declarations, affiliations, ethics wording, database versions, and code repository details require confirmation before journal upload.

## Abstract

### Background

After early circulatory support, restoration of mean arterial pressure may not indicate recovery of organ perfusion. We developed Hemodynamic Recovery Capacity (HRC), a model-based measure of physiologic recovery, to quantify observed recovery relative to expected recovery after early support.

### Methods

We conducted a retrospective multicohort study across MIMIC-IV, eICU, and SICdb. The index time was initiation of early vasoactive circulatory support. Baseline physiology was summarized from 0 to 6 h, response from 6 to 24 h, and outcomes after a 24 h landmark. HRC was constructed from harmonized MAP, vasoactive-support, and urine-output response domains after within-database expected-recovery modeling, cross-fitting, residualization, and standardization. The primary outcome was post-landmark hospital mortality. The key secondary outcome was post-landmark organ nonrecovery. Database-specific estimates were combined using random-effects meta-analysis. MIMIC-IV mechanistic modules evaluated CVP, MPP, cardiogenic shock, and MCS.

### Results

The primary HRC cohorts included 16,757 patients in MIMIC-IV, 9,275 in eICU, and 8,274 in SICdb. Each 1-SD higher HRC was associated with lower post-landmark hospital mortality in MIMIC-IV (OR 0.76, 95% CI 0.73 to 0.80), eICU (OR 0.75, 95% CI 0.71 to 0.80), and SICdb (OR 0.67, 95% CI 0.61 to 0.72), with a pooled OR of 0.73 (95% CI 0.68 to 0.78). HRC was also associated with organ nonrecovery (pooled OR 0.61, 95% CI 0.57 to 0.66). Among patients with restored MAP, absent organ recovery identified a higher-risk subgroup across all databases. In MIMIC-IV, low HRC combined with venous congestion or low MPP had the highest event rates. Findings were directionally consistent across landmark, time-window, missingness, model-specification, and negative-control analyses.

### Conclusions

HRC identified a reproducible dynamic recovery phenotype after early circulatory support across three critical care databases. Lower HRC was associated with mortality and organ nonrecovery and distinguished pressure restoration from broader physiologic recovery. Prospective studies are needed before HRC is used to guide treatment decisions.

Keywords: hemodynamics; circulatory support; vasopressors; organ perfusion; critical care

## Introduction

After vasopressors, fluid resuscitation, or mechanical circulatory support are initiated, clinicians must rapidly decide whether a critically ill patient is recovering or only meeting short-term hemodynamic targets [1-4]. Early circulatory support is intended to restore perfusion, but improvement in macrocirculatory pressure does not necessarily indicate recovery of renal, metabolic, or broader organ function [5-13]. This distinction matters because a patient may reach a conventional MAP target while remaining dependent on high vasoactive support, developing oliguria, accumulating fluid, or failing to clear metabolic stress [1-3,7-9,12-14].

Current bedside assessment of hemodynamic response remains fragmented. MAP, vasopressor dose, lactate, urine output, creatinine, fluid balance, and severity scores are all clinically meaningful, but each captures only part of the response to support [1-3,8,9,12-16]. Static thresholds identify hypotension or organ dysfunction, and raw changes describe whether a variable improved or worsened, but neither approach fully accounts for the patient's initial state or the support intensity required to produce the observed response.

We developed Hemodynamic Recovery Capacity (HRC) as a model-based measure of early physiologic recovery after circulatory support. HRC was designed to quantify the deviation between observed and expected recovery during the 6 to 24 h response window, conditional on baseline physiology, illness severity, and early support intensity. It is not intended to replace severity scores or single bedside markers. Instead, it formalizes a familiar clinical question: after accounting for how sick the patient was and how much support was required, did the patient recover better or worse than expected?

Prior hemodynamic studies have linked arterial pressure, vasopressor exposure, tissue perfusion markers, serial organ dysfunction, and venous congestion to outcomes in critically ill patients [1-16,22-24]. Annals of Intensive Care has also published database-based hemodynamic studies using MIMIC and eICU that evaluated MAP, blood pressure components, and CVP as outcome-associated exposures [30-32]. Whether early post-support recovery can be quantified as a reproducible physiologic phenotype across independent ICU databases remains uncertain.

We derived and externally validated HRC in critically ill adults receiving early circulatory support across MIMIC-IV, eICU, and SICdb. We tested whether lower HRC was associated with post-landmark hospital mortality and organ nonrecovery, evaluated discordance between MAP restoration and organ recovery, and assessed robustness across alternative windows, landmark handling, missingness, model specifications, and falsification analyses. In MIMIC-IV, we further examined CVP, MPP, cardiogenic shock, and MCS modules as mechanistic validation analyses.

## Methods

### Study design and data sources

We conducted a retrospective multicohort study using MIMIC-IV, eICU, and SICdb [17-19]. MIMIC-IV was used as the discovery and mechanistic cohort because it provides granular time-stamped physiologic, treatment, and mechanical circulatory support data. eICU was used as a multicenter external validation cohort. SICdb was used as an independent cross-system validation cohort. The study was designed and reported according to STROBE [20].

### Study population and index time

The primary population consisted of adult ICU patients who received early vasoactive circulatory support within the first 24 h of ICU admission [1-3,21]. For patients with more than one qualifying ICU stay, the first qualifying patient-level ICU stay was used. The index time was the first qualifying vasoactive support initiation within the first ICU day. Patients were required to have baseline measurements, response-window measurements, and post-landmark outcome ascertainment. Patients who died or left the ICU before the 24 h landmark were excluded from the primary HRC cohort because a complete 0 to 24 h recovery exposure could not be measured; they were evaluated in landmark eligibility, IPCW, early-event, and composite sensitivity analyses.

### Analytic windows and harmonization

The analysis used a time-locked landmark design. Baseline physiology was summarized during 0 to 6 h after index time, physiologic response during 6 to 24 h, and outcomes after the 24 h landmark. Variables were harmonized according to physiologic directionality rather than identical raw measurement systems. Core HRC domains were MAP response, vasoactive-support response, and urine-output response. MAP was expressed in mmHg. Urine output was normalized to body weight and time where possible. Vasoactive support was represented as norepinephrine-equivalent dose in MIMIC-IV where feasible and as a harmonized vasoactive burden measure in eICU and SICdb when stable dose equivalence could not be ensured. Lactate, creatinine, CVP, and MPP were used for secondary, construct-validation, or mechanistic analyses according to availability.

### HRC construction

HRC was defined as a model-based residualized measure of physiologic recovery. For each response domain, observed recovery was calculated from the baseline and response windows, with higher values aligned to more favorable recovery. Expected recovery was estimated within each database using multivariable regression models conditioned on baseline physiologic state, illness severity, and early support intensity. Candidate predictors included age, sex where available, database-native illness severity, baseline MAP, baseline vasoactive support burden, and baseline urine output. Post-landmark variables and outcome variables were excluded from expected-recovery models.

Expected-recovery models were trained separately within each database using 5-fold cross-fitting. For each patient, expected recovery was predicted from a model trained without that patient's fold. Domain-specific residual recovery was calculated as observed minus expected recovery and standardized within database. The primary HRC was the equal-weighted mean of standardized residual components, followed by within-database standardization. The framework used the same physiologic structure across databases while allowing database-specific coefficient recalibration. Continuous HRC per 1-SD higher value was the primary exposure; low-HRC was a secondary clinical phenotype.

### Outcomes and statistical analysis

The primary outcome was post-landmark hospital mortality. The key secondary outcome was post-landmark organ nonrecovery, defined from available post-24 h renal, metabolic, and support-related variables informed by organ dysfunction, lactate, and urine-output criteria [8,9,12,13,15,16,21]. Additional analyses evaluated pressure-organ discordance among patients with restored MAP, MIMIC-IV CVP/MPP joint phenotypes, and MIMIC-IV cardiogenic shock/MCS subgroups.

Analyses were performed separately within each database. Logistic regression estimated associations between continuous HRC and mortality or organ nonrecovery, adjusted for baseline covariates measured before or during the 0 to 6 h baseline window. Adjusted dose-response relationships were estimated using natural-spline logistic models. Database-specific estimates were combined using random-effects meta-analysis; fixed-effect estimates were descriptive. Mechanistic analyses were restricted to MIMIC-IV and evaluated CVP/MPP phenotypes and early IABP, Impella, or ECMO support [22-29].

Sensitivity analyses addressed major threats to validity: landmark eligibility, IPCW for landmark eligibility, early death or discharge profiling, early death assigned as worst recovery, alternative 6 to 12 h and 12 to 24 h windows, missingness, measurement frequency, index-time adjustment, fluid balance adjustment, simplified raw-domain HRC, PCA-weighted HRC, no-leakage HRC excluding first-day severity scores, generalized additive expected-recovery models, and low-HRC cutoff sensitivity. A negative-control permutation analysis randomly permuted HRC labels within each database 200 times. Analyses used R version 4.5.1 and Python version 3.9.6.

## Results

### Cohort construction and baseline characteristics

Early vasoactive support was identified in 24,176 adult ICU stays/cases in MIMIC-IV, 22,427 in eICU, and 11,328 in SICdb. After 24 h landmark eligibility and restriction to first qualifying patient-level stay or case, formal modeling cohorts included 17,966 patients in MIMIC-IV, 16,677 in eICU, and 8,700 in SICdb. After core HRC data requirements, the primary HRC cohorts included 16,757 patients in MIMIC-IV, 9,275 in eICU, and 8,274 in SICdb (Figure 1; Table 1).

Baseline characteristics differed across databases. Median age was 67 years in MIMIC-IV, 67 years in eICU, and 70 years in SICdb. Male patients accounted for 60.7%, 56.2%, and 65.0%, respectively. Post-landmark hospital mortality occurred in 2,557 MIMIC-IV patients (15.3%), 1,912 eICU patients (20.6%), and 926 SICdb patients (11.2%). Organ nonrecovery occurred in 5,415 MIMIC-IV patients (37.0%), 3,196 eICU patients (37.4%), and 2,852 SICdb patients (37.2%) among patients with available post-landmark organ-recovery assessment.

### HRC distribution and construct behavior

HRC was standardized within each database after residualization. Median HRC was -0.025 (IQR -0.584 to 0.577) in MIMIC-IV, -0.034 (IQR -0.597 to 0.578) in eICU, and -0.060 (IQR -0.689 to 0.627) in SICdb. Residualized response-domain correlations were modest, indicating that HRC was not driven by a single component. MAP and urine-output residual correlations were 0.121 in MIMIC-IV, 0.141 in eICU, and 0.166 in SICdb. Construct-validation analyses showed expected directional correlations between HRC and raw physiologic recovery signals, including delta MAP, vasoactive-support reduction, and urine-output recovery.

### HRC and post-landmark hospital mortality

Higher HRC was associated with lower post-landmark hospital mortality in all three databases (Figure 2). Each 1-SD higher HRC was associated with lower adjusted odds of mortality in MIMIC-IV (OR 0.76, 95% CI 0.73 to 0.80), eICU (OR 0.75, 95% CI 0.71 to 0.80), and SICdb (OR 0.67, 95% CI 0.61 to 0.72). The random-effects pooled estimate was OR 0.73 (95% CI 0.68 to 0.78). Between-database heterogeneity was present (I2 76.4%), but directionality was consistent.

### HRC and organ nonrecovery

HRC was also associated with organ nonrecovery after the 24 h landmark (Figure 3). Each 1-SD higher HRC was associated with lower odds of organ nonrecovery in MIMIC-IV (OR 0.64, 95% CI 0.62 to 0.67), eICU (OR 0.57, 95% CI 0.54 to 0.60), and SICdb (OR 0.62, 95% CI 0.59 to 0.65). The random-effects pooled estimate was OR 0.61 (95% CI 0.57 to 0.66). Component endpoint analyses were directionally consistent for oliguria, creatinine worsening, lactate nonclearance, persistent hyperlactatemia, and persistent vasopressor use where available.

### Pressure restoration and organ recovery were discordant

Among patients with restored MAP, absence of organ recovery identified a higher-risk subgroup in all three databases (Figure 3). In MIMIC-IV, mortality was 8.3% among MAP-restored patients with organ recovery and 26.4% among MAP-restored patients without organ recovery. Corresponding mortality rates were 13.5% versus 29.0% in eICU and 5.3% versus 16.5% in SICdb.

### MIMIC-IV mechanistic modules

In the MIMIC-IV CVP/MPP module, the reference group without low HRC and without congestion or low MPP had mortality of 6.3% and organ nonrecovery of 29.8%. Patients with congestion or low MPP alone had mortality of 12.0% and organ nonrecovery of 42.6%; patients with low HRC alone had mortality of 9.0% and organ nonrecovery of 39.7%. Patients with both low HRC and congestion or low MPP had the highest event rates, with mortality of 18.4% and organ nonrecovery of 55.4% (Figure 4).

The MCS module included 367 MIMIC-IV patients who received IABP, Impella, or ECMO during the 0 to 24 h window. Higher HRC was directionally associated with lower mortality (OR 0.85, 95% CI 0.69 to 1.06) and lower organ nonrecovery (OR 0.87, 95% CI 0.71 to 1.06), although confidence intervals crossed 1.0. For the combined support or organ nonrecovery endpoint, higher HRC was associated with lower odds of the endpoint (OR 0.70, 95% CI 0.52 to 0.95).

### Sensitivity and falsification analyses

Landmark sensitivity analyses supported the primary findings. Patients excluded before the 24 h landmark accounted for 13.8% of MIMIC-IV patients, 18.6% of eICU patients, and 16.0% of SICdb patients. In IPCW models, HRC remained associated with mortality in MIMIC-IV (OR 0.76, 95% CI 0.73 to 0.80), eICU (OR 0.75, 95% CI 0.71 to 0.79), and SICdb (OR 0.67, 95% CI 0.62 to 0.72). Assigning early deaths the worst recovery value strengthened the association in all three databases.

The association did not depend on the exact response window. Cross-database random-effects estimates remained directionally consistent for the primary 6 to 24 h window (OR 0.74, 95% CI 0.68 to 0.81), early 6 to 12 h window (OR 0.84, 95% CI 0.76 to 0.93), and late 12 to 24 h window (OR 0.74, 95% CI 0.71 to 0.77). Additional analyses addressing measurement frequency, complete-case weighting, multiple imputation, index time, fluid adjustment, no-leakage construction, generalized additive expected-recovery models, PCA weighting, simplified raw-domain HRC, and cutoff sensitivity preserved the direction of association.

In negative-control permutation analyses, permuted HRC estimates centered near OR 1.0 for mortality and organ nonrecovery, whereas observed associations remained separated from the null distributions in all three databases (Figure 4). This indicated that the findings were not reproduced by random assignment of database-specific HRC values.

## Discussion

In this multicohort study of ICU patients receiving early circulatory support, HRC identified a reproducible post-support recovery phenotype across MIMIC-IV, eICU, and SICdb. Higher HRC was consistently associated with lower post-landmark mortality and lower odds of organ nonrecovery. The association was not limited to mortality: HRC also distinguished patients in whom MAP was restored but organ recovery remained absent. In MIMIC-IV, low HRC aligned with higher-risk physiology characterized by venous congestion or low MPP, and the MCS module showed directionally consistent findings despite limited precision.

The main contribution is the framing of early hemodynamic response as a dynamic recovery construct rather than a static severity marker. Conventional ICU risk models summarize how ill a patient is at a given time. HRC instead quantifies how much recovery occurred after accounting for baseline severity, physiologic state, and support intensity. The modest incremental AUC gain after adding HRC to baseline models should therefore be interpreted carefully. HRC is not presented as a stand-alone mortality prediction score, but as a reproducible physiologic response phenotype.

Our findings extend prior hemodynamic literature emphasizing arterial pressure, vasopressor exposure, lactate, urine output, tissue perfusion, and venous congestion in shock [1-14,22-24]. They also extend Annals database-based studies that evaluated MAP, blood pressure components, and CVP as outcome-associated exposures [30-32]. Compared with these studies, the present analysis shifts the focus from single hemodynamic exposures to residualized recovery after support initiation, tests the construct across three databases, and explicitly separates pressure restoration from downstream organ recovery.

The pressure-organ discordance observed here is biologically plausible. Arterial pressure is only one determinant of organ perfusion; renal and systemic recovery may also depend on venous pressure, perfusion pressure gradient, microcirculatory flow, cardiac output, endothelial function, and vasoactive burden [10,11,22-24]. This may explain why HRC, which integrates pressure response, vasoactive-support response, and urine-output response, was more informative than MAP restoration alone.

The MIMIC-IV CVP/MPP analyses support this interpretation but should not be read as proving a single mechanism. CVP and MPP were available only in a subset and may reflect both disease severity and monitoring decisions. Similarly, the MCS module should not be interpreted as device-specific efficacy or failure. Randomized evidence for temporary MCS differs by device, patient selection, and shock context [25-29]. The safer interpretation is that persistent low recovery capacity despite mechanical support may mark refractory shock physiology.

Clinically, HRC may provide a framework for post-resuscitation reassessment. Early circulatory support is often judged by pressure targets, vasopressor trends, lactate, urine output, and clinician gestalt [1-3,7-9,12-14]. HRC formalizes part of that reasoning by asking whether the patient recovered more or less than expected given the starting point and early support intensity. A simplified raw-domain HRC correlated strongly with residualized HRC, suggesting that a bedside version may be possible. However, these findings do not justify using HRC alone to escalate therapy, select mechanical support, or guide withdrawal of care.

This study has limitations. First, it was retrospective and observational; HRC is associated with outcomes but cannot establish that improving HRC would improve survival or organ recovery. Residual confounding by shock etiology, treatment decisions, cardiac function, source control, and unmeasured trajectory may remain. Second, HRC requires survival and data availability through the 24 h landmark. Although landmark audits, IPCW analyses, and early-event composite sensitivities were performed, HRC cannot be directly measured in patients who die or leave the ICU before the response window is complete. Third, database harmonization required physiologic directionality rather than identical raw measurement. Fourth, organ nonrecovery was operationally defined from available variables and was not prospectively adjudicated. Fifth, lactate, CVP, and MPP were not consistently available for the core cross-database HRC definition. Sixth, the MCS module was restricted to MIMIC-IV and underpowered for precise device-specific inference.

In conclusion, HRC identified a reproducible dynamic recovery phenotype after early circulatory support across three ICU databases. Lower recovery capacity was associated with higher mortality and organ nonrecovery, distinguished pressure restoration from organ recovery, and was supported by MIMIC-IV hemodynamic and MCS mechanistic analyses. Future work should test whether prospective HRC monitoring can improve recognition of post-support recovery failure.

## Abbreviations

CI, confidence interval; CVP, central venous pressure; ECMO, extracorporeal membrane oxygenation; HRC, Hemodynamic Recovery Capacity; IABP, intra-aortic balloon pump; ICU, intensive care unit; IPCW, inverse probability weighting for 24 h landmark eligibility; MAP, mean arterial pressure; MCS, mechanical circulatory support; MPP, mean perfusion pressure; OR, odds ratio; PCA, principal component analysis; STROBE, Strengthening the Reporting of Observational Studies in Epidemiology.

## Declarations

### Ethics approval and consent to participate

This study used only de-identified data from MIMIC-IV, eICU-CRD, and SICdb. Access to the databases was obtained through the required database-specific credentialing, training, and data-use agreements where applicable. The original database projects obtained the relevant institutional approvals and consent waivers or exemptions as described in their source documentation and database publications. No direct patient contact occurred and no identifiable patient-level data were accessed by the study investigators. Because the present analysis used only de-identified, publicly or credential-accessible retrospective databases, additional local institutional review board approval was not required.

### Consent for publication

Not applicable. The study used de-identified retrospective database records and did not include identifiable individual patient information.

### Availability of data and materials

The raw patient-level data are not redistributed by the authors because access is governed by the original database licenses and data-use agreements. MIMIC-IV v3.1 is available to credentialed users through PhysioNet at https://physionet.org/content/mimiciv/3.1/; eICU-CRD v2.0 is available through PhysioNet at https://physionet.org/content/eicu-crd/2.0/; and SICdb v1.0.8 is available through PhysioNet at https://physionet.org/content/sicdb/1.0.8/. Derived aggregate results supporting the main findings are provided in the main text, figures, tables, and supplementary materials. Reproduction of patient-level analytic datasets requires independent database access and execution of the study code under the applicable database-use agreements.

### Code availability

Analysis code for cohort construction, variable harmonization, HRC construction, statistical modeling, sensitivity analyses, and figure generation will be made available in a public GitHub repository, with an archived release DOI generated through Zenodo before submission or upon acceptance. The repository will exclude database credentials, local paths, and patient-level derived datasets. Repository URL and Zenodo DOI: to be added after repository release.

### Competing interests

The authors declare that they have no competing interests.

### Funding

This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.

### Authors' contributions

Concept and design: first author and corresponding author. Data access and curation: first author. Statistical analysis: first author. Figure and table generation: first author. Drafting of the manuscript: first author. Critical revision of the manuscript for important intellectual content: all authors. Supervision: corresponding author. All authors read and approved the final manuscript. Author initials will be inserted after the author list is finalized.

### Acknowledgements

Not applicable.

## Figure legends

### Figure 1. Study design, HRC construction, and implementation cohorts

Panel A shows the time-locked analytic design. The index time was initiation of early circulatory support. Baseline physiology was summarized during 0 to 6 h, response during 6 to 24 h, and outcomes after a 24 h landmark. Panel B shows HRC construction from harmonized MAP, vasoactive support burden, and urine output response domains. Expected recovery was estimated within each database conditional on baseline physiologic state, illness severity, and support intensity. Panel C shows the three-database cohort flow.

### Figure 2. HRC and post-landmark hospital mortality

Panel A shows unadjusted mortality across database-specific HRC quartiles. Panel B shows adjusted dose-response relationships between continuous HRC and hospital mortality. Panel C shows database-specific adjusted ORs for mortality per 1-SD higher HRC and the random-effects pooled estimate. ORs below 1 indicate lower odds of mortality with higher recovery capacity.

### Figure 3. HRC and post-support organ recovery

Panel A shows database-specific adjusted ORs for organ nonrecovery and the random-effects pooled estimate. Panel B evaluates pressure-organ discordance among patients with restored MAP. Patients without organ recovery had substantially higher mortality than those with organ recovery, supporting the distinction between pressure restoration and physiologic recovery.

### Figure 4. Mechanistic support and robustness analyses

Panel A shows the joint MIMIC-IV phenotype combining low HRC with venous congestion or low MPP. Panel B shows the MIMIC-IV MCS module. Panel C summarizes key sensitivity analyses, including landmark handling, alternative time windows, fluid adjustment, simplified HRC, missing-data approaches, and no-leakage HRC construction. Panel D shows the negative-control permutation analysis.

## Table

Table 1. Baseline characteristics of the primary HRC cohorts.

## References

1. Monnet X, Messina A, Greco M, Bakker J, Aissaoui N, Cecconi M, et al. ESICM guidelines on circulatory shock and hemodynamic monitoring 2025. Intensive Care Med. 2025;51:1971-2012. doi:10.1007/s00134-025-08137-z.
2. Cecconi M, De Backer D, Antonelli M, Beale R, Bakker J, Hofer C, et al. Consensus on circulatory shock and hemodynamic monitoring: task force of the European Society of Intensive Care Medicine. Intensive Care Med. 2014;40:1795-1815. doi:10.1007/s00134-014-3525-z.
3. Evans L, Rhodes A, Alhazzani W, Antonelli M, Coopersmith CM, French C, et al. Surviving Sepsis Campaign: international guidelines for management of sepsis and septic shock 2021. Intensive Care Med. 2021;47:1181-1247. doi:10.1007/s00134-021-06506-y.
4. Vincent JL, De Backer D. Circulatory shock. N Engl J Med. 2013;369:1726-1734. doi:10.1056/NEJMra1208943.
5. Asfar P, Meziani F, Hamel JF, Grelon F, Megarbane B, Anguel N, et al. High versus low blood-pressure target in patients with septic shock. N Engl J Med. 2014;370:1583-1593. doi:10.1056/NEJMoa1312173.
6. Lamontagne F, Richards-Belle A, Thomas K, Harrison DA, Sadique MZ, Grieve RD, et al. Effect of reduced exposure to vasopressors on 90-day mortality in older critically ill patients with vasodilatory hypotension: a randomized clinical trial. JAMA. 2020;323:938-949. doi:10.1001/jama.2020.0930.
7. Hernández G, Ospina-Tascón GA, Damiani LP, Estenssoro E, Dubin A, Hurtado J, et al. Effect of a resuscitation strategy targeting peripheral perfusion status vs serum lactate levels on 28-day mortality among patients with septic shock: the ANDROMEDA-SHOCK randomized clinical trial. JAMA. 2019;321:654-664. doi:10.1001/jama.2019.0071.
8. Nguyen HB, Rivers EP, Knoblich BP, Jacobsen G, Muzzin A, Ressler JA, et al. Early lactate clearance is associated with improved outcome in severe sepsis and septic shock. Crit Care Med. 2004;32:1637-1642. doi:10.1097/01.CCM.0000132904.35713.A7.
9. Jones AE, Shapiro NI, Trzeciak S, Arnold RC, Claremont HA, Kline JA. Lactate clearance vs central venous oxygen saturation as goals of early sepsis therapy: a randomized clinical trial. JAMA. 2010;303:739-746. doi:10.1001/jama.2010.158.
10. De Backer D, Creteur J, Preiser JC, Dubois MJ, Vincent JL. Microvascular blood flow is altered in patients with sepsis. Am J Respir Crit Care Med. 2002;166:98-104. doi:10.1164/rccm.200109-016OC.
11. Dubin A, Pozo MO, Casabella CA, Pálizas F Jr, Murias G, Moseinco MC, et al. Increasing arterial blood pressure with norepinephrine does not improve microcirculatory blood flow: a prospective study. Crit Care. 2009;13:R92. doi:10.1186/cc7922.
12. Kidney Disease: Improving Global Outcomes (KDIGO) Acute Kidney Injury Work Group. KDIGO clinical practice guideline for acute kidney injury. Kidney Int Suppl. 2012;2:1-138. doi:10.1038/kisup.2012.1.
13. Macedo E, Malhotra R, Bouchard J, Wynn SK, Mehta RL. Oliguria is an early predictor of higher mortality in critically ill patients. Kidney Int. 2011;80:760-767. doi:10.1038/ki.2011.150.
14. Boyd JH, Forbes J, Nakada TA, Walley KR, Russell JA. Fluid resuscitation in septic shock: a positive fluid balance and elevated central venous pressure are associated with increased mortality. Crit Care Med. 2011;39:259-265. doi:10.1097/CCM.0b013e3181feeb15.
15. Vincent JL, Moreno R, Takala J, Willatts S, De Mendonça A, Bruining H, et al. The SOFA (Sepsis-related Organ Failure Assessment) score to describe organ dysfunction/failure. Intensive Care Med. 1996;22:707-710. doi:10.1007/BF01709751.
16. Ferreira FL, Bota DP, Bross A, Mélot C, Vincent JL. Serial evaluation of the SOFA score to predict outcome in critically ill patients. JAMA. 2001;286:1754-1758. doi:10.1001/jama.286.14.1754.
17. Johnson AEW, Bulgarelli L, Shen L, Gayles A, Shammout A, Horng S, et al. MIMIC-IV, a freely accessible electronic health record dataset. Sci Data. 2023;10:1. doi:10.1038/s41597-022-01899-x.
18. Pollard TJ, Johnson AEW, Raffa JD, Celi LA, Mark RG, Badawi O. The eICU Collaborative Research Database, a freely available multi-center database for critical care research. Sci Data. 2018;5:180178. doi:10.1038/sdata.2018.178.
19. Rodemund N, Wernly B, Jung C, Cozowicz C, Koköfer A. The Salzburg Intensive Care database (SICdb): an openly available critical care dataset. Intensive Care Med. 2023;49:700-702. doi:10.1007/s00134-023-07046-3.
20. von Elm E, Altman DG, Egger M, Pocock SJ, Gøtzsche PC, Vandenbroucke JP, et al. The Strengthening the Reporting of Observational Studies in Epidemiology (STROBE) Statement: guidelines for reporting observational studies. PLoS Med. 2007;4:e296. doi:10.1371/journal.pmed.0040296.
21. Singer M, Deutschman CS, Seymour CW, Shankar-Hari M, Annane D, Bauer M, et al. The Third International Consensus Definitions for Sepsis and Septic Shock (Sepsis-3). JAMA. 2016;315:801-810. doi:10.1001/jama.2016.0287.
22. Mullens W, Abrahams Z, Francis GS, Sokos G, Taylor DO, Starling RC, et al. Importance of venous congestion for worsening of renal function in advanced decompensated heart failure. J Am Coll Cardiol. 2009;53:589-596. doi:10.1016/j.jacc.2008.05.068.
23. Damman K, van Deursen VM, Navis G, Voors AA, van Veldhuisen DJ, Hillege HL. Increased central venous pressure is associated with impaired renal function and mortality in a broad spectrum of patients with cardiovascular disease. J Am Coll Cardiol. 2009;53:582-588. doi:10.1016/j.jacc.2008.08.080.
24. Legrand M, Dupuis C, Simon C, Gayat E, Mateo J, Lukaszewicz AC, et al. Association between systemic hemodynamics and septic acute kidney injury in critically ill patients: a retrospective observational study. Crit Care. 2013;17:R278. doi:10.1186/cc13133.
25. van Diepen S, Katz JN, Albert NM, Henry TD, Jacobs AK, Kapur NK, et al. Contemporary management of cardiogenic shock: a scientific statement from the American Heart Association. Circulation. 2017;136:e232-e268. doi:10.1161/CIR.0000000000000525.
26. Baran DA, Grines CL, Bailey S, Burkhoff D, Hall SA, Henry TD, et al. SCAI clinical expert consensus statement on the classification of cardiogenic shock. Catheter Cardiovasc Interv. 2019;94:29-37. doi:10.1002/ccd.28329.
27. Thiele H, Zeymer U, Neumann FJ, Ferenc M, Olbrich HG, Hausleiter J, et al. Intraaortic balloon support for myocardial infarction with cardiogenic shock. N Engl J Med. 2012;367:1287-1296. doi:10.1056/NEJMoa1208410.
28. Thiele H, Zeymer U, Akin I, Behnes M, Rassaf T, Mahabadi AA, et al. Extracorporeal life support in infarct-related cardiogenic shock. N Engl J Med. 2023;389:1286-1297. doi:10.1056/NEJMoa2307227.
29. Møller JE, Engstrøm T, Jensen LO, Eiskjær H, Mangner N, Polzin A, et al. Microaxial flow pump or standard care in infarct-related cardiogenic shock. N Engl J Med. 2024;390:1382-1393. doi:10.1056/NEJMoa2312572.
30. Vincent JL, Nielsen ND, Shapiro NI, Gerbasi ME, Grossman A, Doroff R, et al. Mean arterial pressure and mortality in patients with distributive shock: a retrospective analysis of the MIMIC-III database. Ann Intensive Care. 2018;8:107. doi:10.1186/s13613-018-0448-9.
31. Khanna AK, Kinoshita T, Natarajan A, Schwager E, Linn DD, et al. Association of systolic, diastolic, mean, and pulse pressure with morbidity and mortality in septic ICU patients: a nationwide observational study. Ann Intensive Care. 2023;13:9. doi:10.1186/s13613-023-01101-4.
32. Li DK, Wang XT, Liu DW. Association between elevated central venous pressure and outcomes in critically ill patients. Ann Intensive Care. 2017;7:83. doi:10.1186/s13613-017-0306-1.
