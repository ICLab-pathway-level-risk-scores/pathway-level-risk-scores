# gene_level_analysis

Code and outputs for the gene-level analyses in the manuscript
(Methods: "Gene-level benchmark and score architecture").

## Layout

```
gene_level_analysis/
├── benchmark/                 # Figure 2, Table 2 — gene-level vs pathway-level benchmark
│   ├── internal/              #   MSK-CHORD held-out test set (head-to-head)
│   └── external/              #   GENIE BPC cross-cohort transfer
├── decomposition/             # Figure 4 — constituent-gene decomposition of the
│                              #   compact features (feature -> gene traceback,
│                              #   prevalence x |log HR| contribution)
└── scripts/                   # generating scripts
```

## Figure mapping

| Manuscript | Folder |
|---|---|
| Figure 2a internal C-index, parsimony–performance | `benchmark/internal/` |
| Figure 2b, Table 2 internal-to-external transfer; Supp Table S6 | `benchmark/external/` |
| Supp Table S11 burden-only and same-method benchmarks | `benchmark/burden_samemethod/` |
| Figure 4 constituent-gene decomposition | `decomposition/` |

All arms are evaluated on the same fixed 8:2 split (`random_state=42`, stratified
on `Event_OS`) under fixed MSK-derived coefficients; the compact pathway score is
genomic-only, so the gene-level arms add no clinical covariates. Uncertainty uses
test-set bootstrap (B = 1000). See each subfolder's own `README.md` for details.
