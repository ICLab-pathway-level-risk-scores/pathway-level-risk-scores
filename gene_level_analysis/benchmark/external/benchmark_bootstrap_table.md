# Main Table — Pathway-level vs gene-level benchmark (observed C-index, bootstrap 95% CI)

**Conclusion:** at a matched feature budget, pathway-level compact models showed higher point estimates of C-index in all five cancers internally and in four of five cancers externally, with IDC as the external exception; larger gene-level models can recover modest additional discrimination at the cost of substantially more features. We do not claim gene-level models are statistically outperformed (no formal non-inferiority test).

| Cancer | Model | Selected gene-variant features, n | MSK held-out C | GENIE BPC C | Transfer dC |
|---|---|---|---|---|---|
| IDC | Pathway compact | 5 | 0.646 (0.611-0.681) | 0.625 (0.589-0.661) | -0.022 |
| IDC | Gene (matched-size ≤Kc) | 5 | 0.637 (0.603-0.671) | 0.681 (0.652-0.712) | 0.045 |
| IDC | Gene (pathway genes) | 33 | 0.675 (0.639-0.713) | 0.669 (0.633-0.701) | -0.006 |
| LUAD | Pathway compact | 5 | 0.623 (0.599-0.646) | 0.611 (0.582-0.644) | -0.012 |
| LUAD | Gene (matched-size ≤Kc) | 5 | 0.601 (0.578-0.626) | 0.611 (0.580-0.643) | 0.011 |
| LUAD | Gene (pathway genes) | 53 | 0.646 (0.622-0.671) | 0.623 (0.592-0.654) | -0.023 |
| PAAD | Pathway compact | 5 | 0.582 (0.550-0.615) | 0.574 (0.538-0.609) | -0.008 |
| PAAD | Gene (matched-size ≤Kc) | 5 | 0.555 (0.526-0.582) | 0.564 (0.532-0.594) | 0.009 |
| PAAD | Gene (pathway genes) | 8 | 0.574 (0.540-0.607) | 0.577 (0.542-0.612) | 0.003 |
| PRAD | Pathway compact | 4 | 0.641 (0.601-0.684) | 0.629 (0.580-0.683) | -0.012 |
| PRAD | Gene (matched-size ≤Kc) | 4 | 0.595 (0.554-0.634) | 0.614 (0.562-0.667) | 0.019 |
| PRAD | Gene (pathway genes) | 19 | 0.641 (0.599-0.679) | 0.640 (0.586-0.697) | -0.002 |
| CRC | Pathway compact | 5 | 0.559 (0.528-0.588) | 0.505 (0.471-0.540) | -0.054 |
| CRC | Gene (matched-size ≤Kc) | 2 | 0.528 (0.509-0.547) | 0.500 (0.476-0.524) | -0.028 |
| CRC | Gene (pathway genes) | 45 | 0.569 (0.538-0.599) | 0.606 (0.572-0.637) | 0.037 |

*All models derived on the MSK training split; coefficients and preprocessing parameters fixed before application to the MSK held-out test set and GENIE BPC (no external refitting). C-index = Harrell's concordance; 95% CIs are patient-level bootstrap (B=1000) on each cohort. 'Gene (matched-size ≤Kc)' = unaggregated gene-level model constrained to at most Kc selected gene–variant features; 'Gene (pathway genes)' = unaggregated gene–variant indicators from the same retained-pathway gene set. Uno C-IPCW and time-dependent AUC are in Supplementary. Transfer ΔC = external − internal. IDC = invasive ductal carcinoma.*
