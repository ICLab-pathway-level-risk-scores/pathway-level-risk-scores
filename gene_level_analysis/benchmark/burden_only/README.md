# Burden-only benchmark (Supplementary Table S10)

Script: `gene_level_analysis/scripts/burden_only_benchmark.py` (seed 42, B = 1000).
It writes to `figures_tables/burden_only_benchmark/`; this folder archives those outputs.

## Question
Does the 4-5-feature compact score carry information beyond the overall pathway-alteration burden?

## Inputs
* fixed MSK-CHORD 80:20 split (`random_state=42`, stratified on overall-survival event) from
  `model_input/filter1_age_fixed/<cancer>/<cancer>_{train,test}.csv`;
* GENIE BPC pathway features written by `external_validation/genie_bpc_validation.py`
  (`external_validation/outputs/GENIE_BPC_<cohort>/`).

Run `external_validation/genie_bpc_validation.py` first.

## Arms (both genomic-only)
* `compact_fixed`: the published compact score (`compact_score/compact_score_formula.yaml`)
* `burden_only`: univariable Cox on `PW_HIT_COUNT` (All-pathway alteration count), fitted on MSK train

Both models are fixed on MSK train and applied unchanged to the MSK held-out split and to GENIE BPC
(no refit, no re-standardization).

## Statistics
* Point estimates: observed full-cohort C-index and ΔC = C(compact) − C(burden-only) on the same patients
  (`*_bootmean` columns give bootstrap means; Table 2 reports bootstrap means, so values can differ by ≤0.001).
* 95% CIs: percentiles of 1,000 paired patient-level bootstrap resamples within each evaluation cohort
  (both arms scored on the same resample). Not adjusted for multiple comparisons; exploratory.
* The bootstrap random-number stream is advanced exactly as in the original run (which also evaluated
  comparators that are not reported), so these files are reproduced exactly.

## Files
`arm_metrics.csv` (C-index and CI per cancer × arm × cohort), `paired_deltaC.csv` (compact − burden-only).
All files are aggregate; no patient-level data.
