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
    <cancer>_30run_selected_features.csv    # per-run selected features
    cv_cindex_by_run.csv                    # per-run cross-validated C-index
    feature_selection_frequency_summary.csv # aggregate selection statistics
                                            # (mirrors Supplementary Table S2)
    final_feature_selection_decision.csv    # final compact features per cancer
    README.md

reference_algorithms/
    oax_binary_reference.py         # OAX reference implementation
                                    # (NOT the full IBCGA engine)
    README_OAX.md

external_validation/
    genie_bpc_validation.py         # GENIE BPC no-refit external Cox (Table 3)
    zscore_params/                  # MSK-derived μ/σ per cancer
    scripts/
        external_combined_figures.py    # Figure 6 panels, Supp Figure S11
        genie_bpc_sensitivity_cox.py    # Supp Table S5 (OS-anchor sensitivity)
        figure_style.py                 # shared matplotlib styling
    README.md

KaplanMeier_plot/
    5feature_km.py                  # risk-stratified Kaplan–Meier curves

decision_curve_analysis/
    clinical_utility_5features.py   # Figure 7, Supp Figure S3 (DCA + calibration)
    run_clinical_simulation.py      # Table 4 (trial enrichment, MDT referral,
                                    # fixed-cutoff drift)
    README.md
```

---


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
