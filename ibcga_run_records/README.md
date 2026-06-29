# IBCGA run records

Run-level outputs of the 30 independent IBCGA feature-selection runs
per cancer, plus aggregate selection-frequency statistics and the final
compact feature-selection decision used to assemble the deployed
compact pathway score (`../compact_score/compact_score_formula.yaml`).

Numerical values follow **Supplementary Table S2**.

## Files

| File | Purpose |
|---|---|
| `<cancer>_30run_selected_features.csv` | Per-run selected features (one per cancer; 5 cancers) |
| `cv_cindex_by_run.csv` | Per-(cancer, run) cross-validated Cox C-index |
| `feature_selection_frequency_summary.csv` | Aggregate selection statistics; mirrors Supplementary Table S2 |
| `final_feature_selection_decision.csv` | Final compact features per cancer; matches the deployed YAML |

## Naming convention

The reader-friendly form `"All-pathway alteration count"` (paper / Table S2
/ CSV `feature_name_display`) corresponds to the code-level identifier
`PW_HIT_COUNT` (YAML, scripts). The `PW_` prefix on this global feature
is a historical code convention; it does not mean per-pathway.

Per-pathway features follow `<Pathway> {amplification hit / deletion hit
/ SV hit / mutation rate (z) / any-alteration rate (z) / alt.-class
composite}` in paper form and `PW_<Pathway>_{amp_hit / del_hit / sv_hit
/ mut_rate_z / any_rate_z / zsum}` in code form.

`feature_selection_frequency_summary.csv` columns: `cancer`, `rank`,
`feature_name_paper`, `feature_name_display`, `selected_runs_of_30`,
`selection_frequency`, `mean_main_effect`. The four bookkeeping columns
that earlier mirrored Supplementary Table S2's "Mean |main effect|" and
"Sign stability", the per-run `pr_list`, and the redundant
`feature_name_internal` (still resolvable from `feature_name_display`
via the `<Pathway>_{amp_hit/del_hit/sv_hit/mut_rate_z/any_rate_z/zsum}`
convention plus the `PW_` prefix) have been removed to match the
published table layout. `abs(mean_main_effect)` reproduces the "Mean
main effect" magnitude printed in Supplementary Table S2; the sign of
the orthogonal-array main-effect estimate is preserved.

## Per-cancer IBCGA search parameters (Supplementary Table S6)

| Cancer | r_start | r_end | First G_max | Subsequent G_max |
|---|---|---|---|---|
| IDC  | 20 | 5  | 300 | 100 |
| LUAD | 30 | 15 | 300 | 100 |
| PAAD | 30 | 15 | 300 | 100 |
| PRAD | 20 | 5  | 300 | 100 |
| CRC  | 30 | 15 | 300 | 100 |

`r_start` and `r_end` bound the feature-subset size during the bi-objective
search; `First G_max` is the generation budget for the first IBCGA epoch
and `Subsequent G_max` applies to refinement epochs.

## Caveats

- These records audit the **feature-selection decision**, not the exact
  stochastic trajectory of any single IBCGA run. With
  `random_seed = -1`, individual trajectories are not seed-replayable
  by design.

- External GENIE BPC data did not contribute to candidate-feature
  construction, IBCGA feature selection, IBCGA parameter setting, Cox
  coefficient estimation, standardization parameter estimation, or
  cutoff determination.
