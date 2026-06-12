# External Validation: §4.2 3-Feature Risk Score (5 Cancers, Best Cohort Each)

## Methodology — which logic was used?

**This validation uses `km_3feature_vs_hitcount` logic (cohort-internal rank quantiles).**

| Logic | Cut-point source | Used here? |
|-------|------------------|:---:|
| `km_3feature_vs_hitcount` | external cohort's own rank quantile, `pd.qcut(values.rank, q=q)` | ✅ |
| `km_3feature_vs_hitcount_train_thresholds` | MSK-CHORD train-set quantile cuts applied to test | ❌ |

**MSK partial hazard cutoffs** (BRCA 0.980 / LUAD 1.001 / PAAD 0.948 / PRAD 0.918 / CRC 0.953) are reported as a **secondary** binary high/low row labelled `MSK_cutoff` in each summary CSV, but the primary Q2/Q3/Q4 KM use cohort-internal quantiles.

z-score normalization for `mut_rate` / `any_rate` features uses **MSK train μ/σ** from `zscore_params/{cancer}_zscore_params.csv` (not cohort-internal z).

## Best result per cancer (5-cancer headline)

| Cancer | Cohort | Endpoint | n | events | C-index | 2Q p | 3Q p | 4Q p | Source |
|--------|--------|:---:|:-:|:-:|:---:|:---:|:---:|:---:|---|
| BRCA | METABRIC IDC | OS | 1537 | 891 | 0.571 | **1.7×10⁻⁵** | **8.2×10⁻⁵** | **1.4×10⁻⁴** | brca_metabric (Curtis *Nature* 2012 + Pereira *Nat Commun* 2016) |
| LUAD | OncoSG (Asian) | OS | 299 | 92 | 0.549 | 0.776 | 0.268 | **4.8×10⁻³** | luad_oncosg_2020 (Chen *Genome Med* 2020) |
| PAAD | TCGA GDC | OS | 181 | 98 | 0.577 | **6.2×10⁻³** | **0.011** | **0.047** | paad_tcga_gdc (TCGA reprocessed) |
| PRAD | TCGA PanCancer Atlas | OS | 489 | 9 | **0.774** | 0.32 | 0.41 | 0.47 | prad_tcga_pan_can_atlas_2018 (TCGA *Cell* 2015) |
| CRC | SYSUCC | DFS | 760 | 150 | 0.581 | **7.4×10⁻⁴** | **0.015** | **5.2×10⁻³** | crc_sysucc_2022 (Sun Yat-sen UCC 2022) |

## Layout (best cohorts only)

```
external_validation/
├── README.md
├── raw_data/                              # cBioPortal cohort dumps (input)
│   ├── BRCA_metabric/
│   ├── LUAD_oncosg_2020/
│   ├── PAAD_TCGA_GDC/
│   ├── PRAD_TCGA_pan_can_atlas_2018/
│   └── CRC_sysucc_2022/
├── scripts/                               # one validation script per cohort
│   ├── metabric_external_validation.py    → BRCA
│   ├── luad_oncosg_validation.py          → LUAD
│   ├── paad_tcga_validation.py            → PAAD (locked to GDC only)
│   ├── tcga_external_validation.py        → PRAD (5-cancer aggregate, but only PRAD output retained)
│   └── crc_sysucc_validation.py           → CRC
├── outputs/                               # per-cohort {figures, tables, features}
│   ├── BRCA_metabric/
│   ├── LUAD_oncosg_2020/
│   ├── PAAD_tcga_gdc/
│   ├── PRAD_tcga/
│   └── CRC_sysucc/
├── summary/
│   └── best_result_per_cancer.csv         ← headline 5-cancer summary
└── zscore_params/                         # MSK train μ/σ for mut_rate / any_rate (5 cancers)
```

## Best-result PNG paths

- BRCA: `outputs/BRCA_metabric/figures/metabric_brca_3feature_quantile_km.png`
- LUAD: `outputs/LUAD_oncosg_2020/figures/luad_oncosg_2020_3feature_quantile_km.png`
- PAAD: `outputs/PAAD_tcga_gdc/figures/paad_tcga_gdc_3feature_quantile_km.png`
- PRAD: `outputs/PRAD_tcga/figures/PRAD_tcga_3feature_quantile_km.png`
- CRC: `outputs/CRC_sysucc/figures/crc_sysucc_3feature_quantile_km.png`

## §4.2 features and β coefficients (applied)

| Cancer | Feature 1 (β) | Feature 2 (β) | Feature 3 (β) | Cutoff |
|--------|---|---|---|:---:|
| BRCA | PW_HIT_COUNT (0.102) | PW_TP53_mut_rate_z (0.178) | PW_PI3K_any_rate_z (0.059) | 0.980 |
| LUAD | PW_HIT_COUNT (0.055) | PW_NRF2_zsum (0.106) | PW_TP53_zsum (0.102) | 1.001 |
| PAAD | PW_HIT_COUNT (0.043) | PW_Cell_Cycle_any_rate_z (0.158) | PW_TP53_mut_rate_z (0.065) | 0.948 |
| PRAD | PW_HIT_COUNT (0.121) | PW_TP53_zsum (0.199) | PW_WNT_amp_hit (1.110) | 0.918 |
| CRC | PW_TGF_Beta_zsum (0.050) | PW_RTK_RAS_amp_hit (0.166) | PW_TP53_any_rate_z (0.028) | 0.953 |

## Key methodology notes

- **Sample type**: solid-tumor only (Primary, Metastasis); cfDNA/ctDNA excluded
- **CANCER_TYPE_DETAILED whitelist**: BRCA = IDC; LUAD = LUAD; PAAD = PDAC; PRAD = PRAD; CRC = Colon + Rectal Adenocarcinoma
- **MSK leakage filter**: any cohort with overlapping MSK SAMPLE_IDs is excluded (verified per-cohort; all 5 best cohorts are 0% MSK overlap)
- **METABRIC / SYSUCC**: no SV data → `PW_*_sv_hit = 0` (SV is ~3–6% of HIT_COUNT in MSK, impact minimal)
- **SYSUCC**: no OS_MONTHS → uses DFS endpoint (recurrence/progression)
- **PAAD TCGA GDC**: CNA file uses Entrez_Gene_Id only; mapped to Hugo via mutation file in script

## Re-run

Each script reads from `raw_data/<COHORT>/` (paths in script constants). z-score params loaded from `zscore_params/`.

---

## Reproduction map (manuscript figures/tables)

Primary external validation = **no-refit GENIE BPC** (MSK-contributed
cases excluded).

| Manuscript artifact | Script (`scripts/`) |
|---|---|
| Table 4 (external Cox, 5 cancers) | `genie_bpc_validation.py` |
| Figure 4 | `genie_bpc_validation.py`, `external_combined_figures.py` |
| Supp Table S4 (OS-anchor sensitivity) | `genie_bpc_sensitivity_cox.py` |
| Supp Figure S12 (quintile calibration) | `external_combined_figures.py` |

All scripts load fixed β from `../compact_score/compact_score_formula.yaml`
and MSK-derived μ/σ from `zscore_params/`. No refit, no rescaling.
