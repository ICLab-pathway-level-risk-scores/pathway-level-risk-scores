# Pathway-Level Risk Scores

Code repository for the manuscript:

> **Multi-variant pathway aggregation of sparse targeted-panel alterations for
> prognostic modeling across five solid tumors.**

Figure and table numbers below follow this manuscript.

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
    genie_bpc_validation.py         # GENIE BPC no-refit external Cox (Table 3, Figure 5)
    zscore_params/                  # MSK-derived μ/σ per cancer
    scripts/
        external_combined_figures.py    # external calibration (Supp Figure S11), HR forest
        genie_bpc_sensitivity_cox.py    # Supp Table S4 (OS-anchor sensitivity)
        figure_style.py                 # shared matplotlib styling
    README.md

KaplanMeier_plot/
    5feature_km.py                  # risk-stratified Kaplan–Meier curves

decision_curve_analysis/
    clinical_utility_5features.py   # Supp Figure S3 (decision curves)
    run_clinical_simulation.py      # Supp Table S7, Figure S13 (trial enrichment, MDT referral,
                                    # fixed-cutoff drift)
    README.md
```


## Data access

The clinicogenomic datasets used in this study are publicly available
under their own access agreements and are **not redistributed by this
repository**:

| Dataset | Source |
|---|---|
| MSK-CHORD 2024 | cBioPortal study `msk_chord_2024` |
| AACR Project GENIE BPC | Synapse `syn27056172` |

### Local data layout

Scripts resolve paths relative to the repository root and expect the access-controlled
data (not included here) in the following directories:

| Directory | Contents |
|---|---|
| `raw_data/all_timeline_csv/` | MSK-CHORD 2024 clinical and timeline tables |
| `build_3flag_inputs/` | per-cancer MSK-CHORD alteration tables (`<Cancer>_integrated_genomics.csv`: one row per sample with mutation, CNA and SV calls) and `tcga_10_pathway.xlsx` (copy from `preprocessing/`) |
| `model_input/filter1_data/`, `model_input/sample_id_mapping/` | MSK-CHORD per-cancer cohort tables and train/test sample-ID maps |
| `model_input/filter1_age_fixed/` | written by `preprocessing/build_msk_cohort.py` |
| `external_validation/raw_data/GENIE_BPC/<cohort>/` | GENIE BPC `clinical_sample_nonMSK.csv`, `survival_nonMSK.csv`, `data_mutations.tsv`, `data_cna.tsv`, and `data_sv.tsv` where available (non-MSK samples) |

Run order for the external and benchmark analyses:

```
python3 external_validation/genie_bpc_validation.py               # Table 3, Figure 5 inputs
python3 gene_level_analysis/scripts/benchmark_common_coverage.py  # Table 2
python3 gene_level_analysis/scripts/burden_samemethod_benchmark.py  # Supp Table S11
```

---

## Compact score formula

The fixed Cox formula is documented in
`compact_score/compact_score_formula.yaml`:

```
score(s) = exp( Σ_k β_k · x_k(s) )
```

This is the **partial-hazard score** from the Cox model
(paper terminology; see Supplementary Table S3).

- β coefficients are MSK-derived and **not refit** in external cohorts.
- Feature z-scoring uses MSK-derived μ/σ stored under
  `external_validation/zscore_params/`.
