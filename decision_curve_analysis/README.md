# Decision-curve and clinical-utility analyses

## Reproduction map

| Manuscript artifact | Script |
|---|---|
| Table 3 (trial-enrichment NNS) | `run_clinical_simulation.py` |
| MDT-referral PPV/NPV | `run_clinical_simulation.py` |
| Figure 7 (decision curve + calibration) | `clinical_utility_5features.py` |
| Supp Figure S6 (DCA per cancer) | `clinical_utility_5features.py` |
| Fixed-cutoff drift sensitivity | `run_clinical_simulation.py` |

Both scripts load the deployed compact pathway scores from
`../compact_score/compact_score_formula.yaml`.
