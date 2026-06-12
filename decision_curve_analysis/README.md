# Decision-curve and clinical-utility analyses

Reproduces decision-curve, calibration, trial-enrichment, MDT-referral,
and fixed-cutoff drift analyses using the deployed compact pathway
scores (`../compact_score/compact_score_formula.yaml`).

## Files

- `clinical_utility_5features.py` — adjusted Cox forest plot, C-index
  comparison, decision-curve analysis, and calibration plots.

Trial-enrichment, MDT-referral, and fixed-cutoff drift simulations are
produced by `../analyses/clinical_simulation/run_clinical_simulation.py`.

## Reproduction map

| Manuscript artifact | Script |
|---|---|
| Table 3 (trial-enrichment NNS) | `../analyses/clinical_simulation/run_clinical_simulation.py` |
| MDT-referral PPV/NPV | `../analyses/clinical_simulation/run_clinical_simulation.py` |
| Figure 7 (decision curve + calibration) | `clinical_utility_5features.py` |
| Supplementary Figure S6 (DCA per cancer) | `clinical_utility_5features.py` |
| Fixed-cutoff drift sensitivity | `../analyses/clinical_simulation/run_clinical_simulation.py` |
