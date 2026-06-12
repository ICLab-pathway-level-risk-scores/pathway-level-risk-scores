# External validation — GENIE BPC

## Reproduction map

| Manuscript artifact | Script (`scripts/`) |
|---|---|
| Table 4 | `genie_bpc_validation.py` |
| Figure 4 | `genie_bpc_validation.py`, `external_combined_figures.py` |
| Supp Table S4 | `genie_bpc_sensitivity_cox.py` |
| Supp Figure S12 | `external_combined_figures.py` |

All scripts load fixed β from `../compact_score/compact_score_formula.yaml`
and MSK-derived μ/σ from `zscore_params/`. No refit, no rescaling.
