# Pathway-Level Risk Scores

Code repository for the paper:

> **Rank-Portable Pathway-Level Prognostic Scores from Routine Targeted-Panel
> Sequencing in Five Solid Tumors.**
> *npj Precision Oncology*, under review.

This repository provides the **fixed compact-score formula** and the
analysis code to apply it to new patient cohorts, perform external
validation without refitting, and reproduce the Kaplan–Meier and decision-
curve analyses reported in the paper.

The compact scores were derived from MSK-CHORD 2024 (BRCA, LUAD, PAAD,
PRAD, CRC) and validated on independent GENIE BPC cohorts **without
refitting coefficients or rescaling features**.

---

## Repository structure

```
preprocessing/
    build_pathway_features.py       # mutation + CNA + SV → pathway-level features
    build_msk_cohort.py             # patient filtering + train/test split
    tcga_10_pathway.xlsx            # TCGA 10-pathway gene definitions

compact_score/
    apply_compact_score.py          # apply fixed Cox formula to new cohort
    compact_score_formula.yaml      # per-cancer fixed β coefficients,
                                    # MSK-derived z-score parameters reference,
                                    # MSK-derived partial-hazard cutoffs

ibcga_run_records/
    <cancer>_30run_selected_features.csv    # per-run selected non-fixed features
    cv_cindex_by_run.csv                    # per-run cross-validated C-index
    feature_selection_frequency_summary.csv # aggregate selection statistics
                                            # (mirrors Supplementary Table S1)
    final_feature_selection_decision.csv    # final compact features per cancer
    README.md                               # schemas, naming conventions, caveats

reference_algorithms/
    oax_binary_reference.py         # fixed-K binary OAX reference implementation
                                    # (NOT the full IBCGA engine)
    README_OAX.md                   # algorithm description and citations

external_validation/
    genie_bpc_validation.py         # GENIE BPC no-refit external Cox
    zscore_params/                  # MSK-derived μ/σ per cancer
    README.md                       # reproduction map for Table 4, Figure 4,
                                    # Supplementary Table S4, Supp Figure S12

KaplanMeier_plot/
    5feature_km.py                  # risk-stratified Kaplan–Meier curves

decision_curve_analysis/
    clinical_utility_5features.py   # decision-curve + calibration analysis
    README.md                       # reproduction map for Table 3, Figure 7,
                                    # Supp Figure S6
```

---

## Reproduction map

| Manuscript artifact                                | Folder / script(s)                                                                  |
|----------------------------------------------------|--------------------------------------------------------------------------------------|
| Table 1 (cohort construction)                      | `preprocessing/`                                                                     |
| Table 2 (internal Cox performance)                 | `compact_score/`                                                                     |
| Table 3 (rank-based clinical decision simulations) | `decision_curve_analysis/` + `analyses/clinical_simulation/run_clinical_simulation.py` |
| Table 4 (no-refit GENIE BPC external Cox)          | `external_validation/genie_bpc_validation.py`                                        |
| Figure 1 (overview)                                | (BioRender; no code)                                                                 |
| Figure 2 (fixed β coefficients per cancer)         | `compact_score/` (β from `compact_score_formula.yaml`)                               |
| Figure 3 (internal validation)                     | `KaplanMeier_plot/`                                                                  |
| Figure 4 (external validation)                     | `external_validation/`                                                               |
| Figure 7 (decision-curve + clinical utility)       | `decision_curve_analysis/clinical_utility_5features.py`                              |
| Supp Table S1 (top-15 IBCGA selection statistics)  | `ibcga_run_records/feature_selection_frequency_summary.csv`                          |
| Supp Table S3 (cohort construction)                | `preprocessing/`                                                                     |
| Supp Table S4 (GENIE BPC sensitivity analyses)     | `external_validation/scripts/genie_bpc_sensitivity_cox.py`                           |
| Supp Table S5 (IBCGA per-cancer parameters)        | (parameters documented; see `ibcga_run_records/README.md`)                           |
| Supp Figure S1 (CV C-index stability)              | `ibcga_run_records/cv_cindex_by_run.csv`                                             |
| Supp Figure S2 (post-IBCGA m=1..10 sensitivity)    | `ibcga_run_records/feature_selection_frequency_summary.csv` + Cox refit              |
| Supp Figure S6 (DCA per cancer)                    | `decision_curve_analysis/clinical_utility_5features.py`                              |
| Supp Figure S12 (external quintile calibration)    | `external_validation/scripts/external_combined_figures.py`                           |

---

## Data access

The clinicogenomic datasets used in this study are publicly available
under their own access agreements and are **not redistributed by this
repository**:

| Dataset | Source |
|---|---|
| MSK-CHORD 2024 | cBioPortal study `msk_chord_2024` |
| AACR Project GENIE BPC | Synapse `syn27056172` |

---

## Compact score formula

The fixed Cox formula is documented in
`compact_score/compact_score_formula.yaml`:

```
linear_predictor(s) = Σ  β_i · feature_i(s)
risk_score(s)       = exp(linear_predictor(s))
```

- β coefficients are MSK-derived and **not refit** in external cohorts.
- Feature z-scoring uses MSK-derived μ/σ stored under
  `external_validation/zscore_params/`.
