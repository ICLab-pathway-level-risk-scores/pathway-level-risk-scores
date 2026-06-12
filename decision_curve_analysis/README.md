# Decision-curve and clinical-utility analyses

This folder reproduces decision-curve, calibration, trial-enrichment,
multidisciplinary-team (MDT) referral, and fixed-cutoff drift analyses for
the deployed compact pathway scores. All analyses use the fixed compact-score
outputs (`../compact_score/compact_score_formula.yaml`) and cohort-specific
rank thresholds; no coefficients or μ/σ are refit here.

## Files

- `clinical_utility_5features.py` — adjusted Cox forest plot for the compact
  pathway score, clinical-only vs score-only vs combined C-index table,
  per-cancer decision-curve analysis at the better of two candidate
  horizons, and combined-model calibration plots.
- `best_result_per_cancer.csv` — headline summary across cancers.

The trial-enrichment, MDT-referral, surveillance, and fixed-cutoff drift
scenarios are produced by `../analyses/clinical_simulation/run_clinical_simulation.py`,
with output CSVs `scenario1_trial_enrichment.csv`,
`scenario2_surveillance.csv`, `scenario3_mdt_referral.csv`, and
`summary_table.csv`.

## Reproduction map (manuscript figures/tables)

| Manuscript artifact                                     | Script                                                                     | Output                                              |
|----------------------------------------------------------|-----------------------------------------------------------------------------|-----------------------------------------------------|
| Main Table 3 (Q4 trial-enrichment NNS, enrichment ratio) | `../analyses/clinical_simulation/run_clinical_simulation.py`                | `../analyses/clinical_simulation/scenario1_trial_enrichment.csv`, `table4A_scenario1_pivot.csv` |
| Main Figure 7 (decision curve + calibration)             | `clinical_utility_5features.py`                                             | `../figures_tables/clinical_utility_5features/`     |
| MDT-referral PPV/NPV table                               | `../analyses/clinical_simulation/run_clinical_simulation.py`                | `scenario3_mdt_referral.csv`, `table4B_scenario3_pivot.csv` |
| Supplementary Figure S6 (DCA per cancer)                 | `clinical_utility_5features.py`                                             | `../figures_tables/clinical_utility_5features/dca_*.png` |
| Fixed MSK-derived cutoff drift sensitivity               | `../analyses/clinical_simulation/run_clinical_simulation.py`                | `summary_table.csv` (fixed-cutoff rows)             |
| Panel-coverage / SV-availability sensitivity             | `../analyses/clinical_simulation/panel_coverage_sensitivity.py`, `../analyses/clinical_simulation/sv_sensitivity.py` | `panel_coverage_sensitivity.csv`, `sv_sensitivity.csv` |

## Caveats

These analyses are **prognostic-enrichment simulations** at cohort-specific
rank thresholds and do not imply treatment-selection benefit. The
fixed-cutoff drift analysis demonstrates why MSK-derived absolute cutoffs
require local recalibration before any absolute-risk use in an external
cohort, motivating the rank-portable deployment posture adopted in the
paper.
