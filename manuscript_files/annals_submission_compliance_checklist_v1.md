# Annals of Intensive Care Submission Compliance Checklist v1

Target journal: Annals of Intensive Care

Article type: Research Article

Prepared: 2026-07-06

## Official Requirements Checked

| Requirement | Target journal requirement | Current manuscript status | Action |
|---|---:|---|---|
| Article type | Original Research / Research Article | Matched | Keep as Research Article |
| Main text length | <= 4,000 words | Approximately 2,810 words | Pass |
| Abstract | Structured, <= 350 words | Approximately 340 words | Pass |
| Keywords | 3 to 5 | 5 keywords selected | Keep concise; avoid too many database-specific keywords |
| References | <= 40 | 32 references | Pass |
| Main display items | <= 5 total tables and/or figures | 4 figures + 1 table | Pass |
| Reporting guideline | STROBE for observational studies | STROBE cited and supplementary map available | Keep STROBE checklist in supplement |
| Data availability | Required in Declarations | Drafted with database links and locked versions: MIMIC-IV v3.1, eICU-CRD v2.0, SICdb v1.0.8 | Pass |
| Ethics | Required | Drafted with de-identified database wording and no additional local IRB requirement | Confirm final author approval of wording |
| Cover letter | Should explain fit, policies, conflicts, author approval, originality | Draft prepared | Add corresponding author details |
| Editable manuscript | DOCX/RTF/LaTeX accepted | Markdown source available; Word conversion still pending | Convert after author fields are filled |

Official sources checked:

- ScienceDirect Guide for Authors: https://www.sciencedirect.com/journal/annals-of-intensive-care/publish/guide-for-authors
- Springer submission guidelines: https://link.springer.com/journal/13613/submission-guidelines

## Published Annals Article Pattern Used

The manuscript package follows the pattern observed in Annals database/hemodynamic studies:

| Benchmark article | Why it matters for this manuscript |
|---|---|
| Vincent et al., 2018, MIMIC-III MAP and mortality | Shows Annals accepts database-based hemodynamic exposure studies with clinically interpretable pressure endpoints |
| Khanna et al., 2023, eICU blood pressure components | Shows Annals accepts large eICU observational hemodynamic studies with mortality and organ-dysfunction outcomes |
| Li et al., 2017, MIMIC CVP and outcomes | Supports the CVP/venous-congestion mechanistic framing |

These articles are cited in the manuscript as refs 30-32. Other Annals database articles were used as format benchmarks but were not forced into the reference list unless they directly supported the scientific argument.

## Locked Main Display Plan

| Display item | Content | File |
|---|---|---|
| Figure 1 | Study design, HRC construction, and cohort flow | `outputs/annals_main_v5/figures/figure1_annals_hrc_concept_flow_v5.tiff` |
| Table 1 | Baseline characteristics across three HRC cohorts | `outputs/annals_main_v5/tables/table1_baseline_characteristics.md` |
| Figure 2 | HRC distribution, mortality spline, mortality forest plot | `outputs/annals_main_v5/figures/figure2_annals_hrc_mortality_v5.tiff` |
| Figure 3 | Organ nonrecovery and pressure-organ discordance | `outputs/annals_main_v5/figures/figure3_annals_organ_nonrecovery_v5.tiff` |
| Figure 4 | CVP/MPP, MCS module, robustness, negative control | `outputs/annals_main_v5/figures/figure4_annals_mechanism_robustness_v5.tiff` |

Total main display items: 5.

## Files Prepared in This Submission Package

| File | Role |
|---|---|
| `manuscript_annals_submission_ready_v1.md` | Condensed Annals-ready manuscript text |
| `title_page_annals_v1.md` | Title page with author placeholders |
| `cover_letter_annals_v1.md` | Annals-targeted cover letter draft |
| `annals_submission_compliance_checklist_v1.md` | This checklist |

## Author-Confirmation Items Still Required

- Author names, affiliations, corresponding author address, email, and ORCID IDs.
- Author contributions by initials.
- Public code repository URL and release identifier.

## Submission Risk Notes

- The HRC construct is novel; claims must stay framed as a reproducible physiologic recovery phenotype, not as a validated treatment-guidance tool.
- The MCS module remains mechanistic and MIMIC-IV-only; avoid device-efficacy language.
- The 24 h landmark remains a methodologic vulnerability, but the manuscript now explicitly supports it with landmark audit, IPCW, early-event composite, and sensitivity analyses.
- Main text should remain compact. Detailed sensitivity results belong in the supplement.
