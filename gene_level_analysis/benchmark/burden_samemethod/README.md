# Burden-only and same-method benchmarks (Supplementary Table S11)

Script: `gene_level_analysis/scripts/burden_samemethod_benchmark.py` (seed 42, B = 1000).
It writes to `figures_tables/burden_samemethod_benchmark/`; this folder archives those outputs.

## Inputs
Same inputs as the Table 2 benchmark (`benchmark_common_coverage.py`):
* fixed MSK-CHORD 80:20 split (`random_state=42`, stratified on overall-survival event) from
  `model_input/filter1_age_fixed/<cancer>/<cancer>_{train,test}.csv`;
* MSK gene-level events from `build_3flag_inputs/<cancer>_integrated_genomics.csv`;
* GENIE BPC pathway features written by `external_validation/genie_bpc_validation.py`
  (`external_validation/outputs/GENIE_BPC_<cohort>/`), and GENIE BPC gene events rebuilt from
  `external_validation/raw_data/GENIE_BPC/<cohort>/data_{mutations,cna}.tsv`.

Run `external_validation/genie_bpc_validation.py` first.

## Arms (all genomic-only)
* `compact_fixed`: the published compact score (`compact_score/compact_score_formula.yaml`)
* `burden_only`: univariable Cox on `PW_HIT_COUNT` (All-pathway alteration count), fitted on MSK train
* `pathway_EN_capK` / `pathway_EN_free`: Elastic-Net Cox on the pathway candidate pool, at most K_c terms / uncapped
* `gene_EN_capK` / `gene_EN_free`: the same Elastic-Net Cox routine on gene × alteration indicators

Every model (scaler, coefficients, candidate columns, penalty) is fixed on MSK train and applied unchanged to the
MSK held-out split and to GENIE BPC (no refit, no re-standardization).

## Settings
* Elastic-Net: `gene_vs_pathway_headtohead.coxnet_select` (same l1 ratio, penalty path and five-fold CV for both
  representations). With a term cap, the eligible model with the best CV C-index is kept (`cap_mode="cv_best"`),
  so realized term counts can be below K_c (LUAD pathway arm: 2 terms). Without a cap, the CV-chosen penalty.
* Gene candidates: retained-pathway genes with ≥1% MSK-train prevalence, restricted to genes present in the
  corresponding GENIE BPC cohort's mutation or copy-number file (assay-content coverage proxy; no external
  outcomes used). SV indicators are not used externally (SV set to 0, as in the compact-score external procedure).
* `PW_DRIVER_BURDEN_z` is excluded from the pathway pool (not reconstructed for GENIE BPC; no compact score uses
  it). `driver_burden_sensitivity.csv`: adding it back changes internal held-out C by ≤0.003.
* External coverage handling: missing copy-number calls and absent mutation records are coded as unaltered; no
  sample-level panel gene lists are used.

## Statistics
* Point estimates: observed full-cohort C-index and ΔC = C(arm A) − C(arm B) on the same patients
  (`*_bootmean` columns give bootstrap means; Table 2 reports bootstrap means, so values can differ by ≤0.001).
* 95% CIs: percentiles of 1,000 paired patient-level bootstrap resamples within each evaluation cohort
  (both arms scored on the same resample). Not adjusted for multiple comparisons; exploratory.

## Files
`arm_metrics.csv` (C-index and CI per cancer × arm × cohort), `paired_deltaC.csv` (planned contrasts),
`driver_burden_sensitivity.csv`. All files are aggregate; no patient-level data.
