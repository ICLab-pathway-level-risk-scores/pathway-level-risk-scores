# External validation — GENIE BPC

## Reproduction map

| Manuscript artifact | Script |
|---|---|
| Table 3, Figure 5 | `genie_bpc_validation.py` |
| Supp Figure S11 (calibration) | `scripts/external_combined_figures.py` |
| Supp Table S4 | `scripts/genie_bpc_sensitivity_cox.py` |

All scripts load fixed β from `../compact_score/compact_score_formula.yaml`
and MSK-derived μ/σ from `zscore_params/`. No refit, no rescaling.
