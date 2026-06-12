# External validation — GENIE BPC

## Reproduction map

| Manuscript artifact | Script |
|---|---|
| Table 3 | `genie_bpc_validation.py` |
| Figure 6 | `genie_bpc_validation.py`, `scripts/external_combined_figures.py` |
| Supp Table S4 | `scripts/genie_bpc_sensitivity_cox.py` |
| Supp Figure S12 | `scripts/external_combined_figures.py` |

All scripts load fixed β from `../compact_score/compact_score_formula.yaml`
and MSK-derived μ/σ from `zscore_params/`. No refit, no rescaling.
