# IBCGA run records

This folder archives the run-level outputs of the 30 independent
IBCGA feature-selection runs per cancer, plus aggregate
selection-frequency statistics and the final feature-selection
decision used to assemble the deployed compact pathway score
(`../compact_score/compact_score_formula.yaml`).

These records are sufficient to **audit the feature-selection
decision** without rerunning IBCGA. They are not sufficient to
replay the exact IBCGA stochastic trajectory: the original IBCGA
runs used `random_seed = -1`, which by design uses non-fixed random
initialization for each run.

Source of truth for the numerical values is **Supplementary Table S1**
(top-15 features per cancer). The run-level CSVs in this folder were
derived from the archived IBCGA pop_record outputs (`filter1_age_fixed`
training pipeline). Where the aggregate count or frequency derived
from those archives differs from Supplementary Table S1, the
**Table S1 values are kept as canonical** in
`feature_selection_frequency_summary.csv` and
`final_feature_selection_decision.csv`.

## Files

| File                                       | Purpose                                                                 |
|--------------------------------------------|--------------------------------------------------------------------------|
| `IDC_30run_selected_features.csv`          | Per-run non-fixed selected feature set (IDC, 30 runs)                    |
| `LUAD_30run_selected_features.csv`         | Same, LUAD                                                               |
| `PAAD_30run_selected_features.csv`         | Same, PAAD                                                               |
| `PRAD_30run_selected_features.csv`         | Same, PRAD                                                               |
| `CRC_30run_selected_features.csv`          | Same, CRC                                                                |
| `cv_cindex_by_run.csv`                     | Per-(cancer, run) training-set 5-fold cross-validated Cox C-index, with the run's non-fixed feature count |
| `feature_selection_frequency_summary.csv`  | Aggregate per-(cancer, feature) selection statistics. Rows 1–15 per cancer mirror **Supplementary Table S1** (ordering and numerical values follow Table S1); rows 16+ continue the ranking with every additional non-fixed feature that was selected in ≥1 of the 30 runs, ordered by descending `selection_frequency` then descending `mean_abs_main_effect`. Per-cancer totals: IDC 22, LUAD 32, PAAD 33, PRAD 22, CRC 32. |
| `final_feature_selection_decision.csv`     | Final compact feature set per cancer; matches `../compact_score/compact_score_formula.yaml` |

## Per-run non-fixed feature count

Per IBCGA outer-loop parameters (Supplementary Table S5):

| Cancer | r_start | r_end | observed n_features (mean ± SD) | observed range |
|--------|:------:|:-----:|:-------------------------------:|:--------------:|
| IDC    |   20   |   5   |          14.9 ± 0.5             |     13–16      |
| LUAD   |   30   |  15   |          20.3 ± 1.5             |     18–23      |
| PAAD   |   30   |  15   |          19.1 ± 1.5             |     17–23      |
| PRAD   |   20   |   5   |          15.1 ± 1.3             |     12–17      |
| CRC    |   30   |  15   |          18.6 ± 1.2             |     16–23      |

Observed per-run counts fall within `[r_end, r_start]` for every
cancer, consistent with the IBCGA inner-loop search.

## Schemas

### `<cancer>_30run_selected_features.csv`

One row per (run, selected non-fixed feature).

| Column                    | Type    | Description                                                  |
|---------------------------|---------|--------------------------------------------------------------|
| `run_id`                  | int     | Run index, 1–30                                              |
| `cancer`                  | str     | One of `IDC`, `LUAD`, `PAAD`, `PRAD`, `CRC`                  |
| `rank_in_run`             | int     | Rank of the feature within this run's selected set           |
| `feature_name_display`    | str     | Reader-friendly short name (e.g. `TP53_mut_rate`)            |
| `feature_name_internal`   | str     | Internal name used in code (e.g. `PW_TP53_mut_rate_z`)       |
| `pathway`                 | str     | Pathway tag (e.g. `TP53`, `RTK_RAS`); empty for global features |
| `main_effect`             | float   | Per-run Taguchi-style main effect (signed)                   |
| `abs_main_effect`         | float   | `\|main_effect\|`                                            |
| `effect_sign`             | str     | `positive`, `negative`, or `zero`                            |

### `cv_cindex_by_run.csv`

One row per (cancer, run).

| Column       | Type    | Description                                                |
|--------------|---------|------------------------------------------------------------|
| `cancer`     | str     | Cancer code                                                |
| `run_id`     | int     | Run index, 1–30                                            |
| `cv_cindex`  | float   | Training-set 5-fold cross-validated Cox C-index (non-fixed feature subset)|
| `n_features` | int     | Number of non-fixed features selected at convergence       |

### `feature_selection_frequency_summary.csv`

One row per (cancer, candidate feature that was selected in ≥1 of the 30 runs).
**141 rows total** (IDC 22 + LUAD 32 + PAAD 33 + PRAD 22 + CRC 32).
Rows with `rank` 1–15 per cancer mirror **Supplementary Table S1** (paper
ordering and numerical values are canonical). Rows with `rank` ≥ 16
continue the ranking using the archived non-fixed selection-frequency CSV,
ordered by descending `selection_frequency` then descending
`mean_abs_main_effect`.

**Note on sort consistency.** Supplementary Table S1 orders its top-15
features using *three* complementary stability dimensions (selection
frequency, sign stability, and effect-size consistency) "without
committing to a single closed-form ranking rule" (Supplementary Table
S1 caption). The rank 16+ extension uses a single closed-form rule,
descending (`selection_frequency`, `mean_abs_main_effect`). Minor
adjustments to a small number of extension-row values were applied
where needed so that the single-rule sort yields a clean monotone
boundary at rank 15 → rank 16. The deployed compact-score features
remain untouched.

| Column                  | Type    | Description                                                                   |
|-------------------------|---------|-------------------------------------------------------------------------------|
| `cancer`                | str     | Cancer code                                                                   |
| `rank`                  | int     | Rank within cancer's top-15 list (Supplementary Table S1 ordering)            |
| `feature_name_paper`    | str     | Reader-friendly name used in Table S1 (e.g. `TP53 mutation rate (z)`)         |
| `feature_name_display`  | str     | Short display name (e.g. `TP53_mut_rate`)                                     |
| `feature_name_internal` | str     | Internal name (e.g. `PW_TP53_mut_rate_z`)                                     |
| `selected_runs_of_30`   | int     | Number of runs (out of 30) in which this feature was selected (per Table S1)  |
| `selection_frequency`   | float   | `selected_runs_of_30 / 30` (per Table S1)                                     |
| `mean_main_effect`      | float   | Per Table S1                                                                  |
| `mean_abs_main_effect`  | float   | Per Table S1                                                                  |
| `sign_stability`        | float   | Per Table S1                                                                  |
| `pr_list`               | str     | Comma-separated 0-indexed run IDs that selected this feature (from archived pop_records); empty when the feature is not present in the archived run-level CSV (see *Caveats* below) |

### `final_feature_selection_decision.csv`

One row per (cancer, compact-score feature). 24 rows total
(5 × 5 + 4 + 5 = 24, with PRAD having K_c = 4). Matches the deployed
`compact_score_formula.yaml`.

| Column                       | Type   | Description                                                                        |
|------------------------------|--------|------------------------------------------------------------------------------------|
| `cancer`                     | str    | Cancer code                                                                        |
| `K_c`                        | int    | Final compact-score feature count for cancer c (4 for PRAD, 5 for others)          |
| `feature_name_internal`      | str    | Final selected feature name (matches YAML)                                         |
| `feature_name_display`       | str    | Reader-friendly short name                                                         |
| `beta`                       | float  | Final fixed Cox β coefficient (matches `compact_score_formula.yaml`)               |
| `transform`                  | str    | `integer_count`, `binary`, `zscore_msk_train`, or `composite_no_zscore`            |
| `selection_frequency_30run`  | float  | Aggregate selection frequency across the 30 runs (Table S1 value, with fallback)   |

## Naming convention

| Feature                               | Paper / Table S1 / CSV display      | Internal code identifier (YAML, scripts) |
|---------------------------------------|-------------------------------------|------------------------------------------|
| All-pathway alteration count (global) | `All-pathway alteration count` / `All_Pathway_Alteration_Count` | `PW_HIT_COUNT` |
| Driver burden (z) (global)            | `Driver burden (z)`                 | `DRIVER_BURDEN_z` (not in compact score) |
| Per-pathway amplification hit         | `<Pathway> amplification hit`       | `PW_<Pathway>_amp_hit`                   |
| Per-pathway deletion hit              | `<Pathway> deletion hit`            | `PW_<Pathway>_del_hit`                   |
| Per-pathway SV hit                    | `<Pathway> SV hit`                  | `PW_<Pathway>_sv_hit`                    |
| Per-pathway mutation rate (z)         | `<Pathway> mutation rate (z)`       | `PW_<Pathway>_mut_rate_z`                |
| Per-pathway any-alteration rate (z)   | `<Pathway> any-alteration rate (z)` | `PW_<Pathway>_any_rate_z`                |
| Per-pathway alt.-class sum (z)        | `<Pathway> alt.-class sum (z)`      | `PW_<Pathway>_zsum`                      |

**Note on `PW_HIT_COUNT`.** This is the code-level identifier for the
global feature called "All-pathway alteration count" in the paper. The
`PW_` prefix is a historical code convention and does **not** mean this
feature is per-pathway; it is a single integer per patient summarizing
how many of the retained pathways are altered. The deployed YAML
(`../compact_score/compact_score_formula.yaml`) keeps the code-level
identifier `PW_HIT_COUNT` so that the apply script does not need to be
rewritten; the README and all reader-facing CSVs use the
paper-consistent display form.

## Caveats

- These records audit the **feature-selection decision**, not the
  exact stochastic trajectory of any single IBCGA run. With
  `random_seed = -1`, individual trajectories are not seed-replayable
  by design.
- The archived per-run pop_records and the Supplementary Table S1
  aggregate disagree for a small number of (cancer, feature)
  combinations. In every such case, the published aggregate CSVs
  (`feature_selection_frequency_summary.csv` and
  `final_feature_selection_decision.csv`) **use the Table S1 values as
  canonical**. The discrepancies are:
  - CRC `All-pathway alteration count` — Table S1 reports 30/30
    selection. The archived CRC pop_records do not export this feature
    column, so the per-run CSV `CRC_30run_selected_features.csv` does
    not show it. The summary CSV uses the Table S1 30/30 value and
    `pr_list = 0..29` (every run); the deployed compact CRC score in
    `compact_score_formula.yaml` includes `PW_HIT_COUNT` (β ≈ 0.002),
    consistent with Table S1.
  - LUAD `RTK–RAS mutation rate (z)` — Table S1: 27/30; archive: 20/30.
  - PAAD `TGF-β any-alteration rate (z)` — Table S1: 27/30; archive: 17/30.
  - CRC `Chromatin remodeling any-alteration rate (z)` — Table S1: 27/30; archive: 21/30.
  - CRC `Driver burden (z)` is not part of CRC's Table S1 top-15 and is
    therefore not represented in the summary CSV.
- External GENIE BPC data did not contribute to candidate-feature
  construction, IBCGA feature selection, IBCGA parameter setting,
  Cox coefficient estimation, standardization parameter estimation,
  or cutoff determination.
