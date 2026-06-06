"""genie_bpc_validation.py
===========================
GENIE BPC public release external validation for the 4 cohorts overlapping
our MSK-trained compact pathway score:

    NSCLC (LUAD)  - 1038 non-MSK samples, 976 OS-evaluable, 585 events
    CRC           -  817                ,  789               , 396
    PANC (PAAD)   -  580                ,  573               , 489
    Prostate (PRAD) - 557               ,  555               , 210

For each cohort:
  1. Filter to the relevant ONCOTREE codes
  2. Rebuild PW_* features using MSK pathway gene sets + MSK z-score params
  3. Compute compact risk_score with published MSK β (compact_score_formula.yaml)
  4. Output:
       outputs/GENIE_BPC_<cohort>/<cohort>_pathway_features.csv
       outputs/GENIE_BPC_<cohort>/figures/  (KM + calibration)
       outputs/GENIE_BPC_<cohort>/tables/   (C-index, multivar Cox, calibration)

SV is unavailable for non-MSK panels in GENIE BPC — all PW_*_sv_hit = 0,
mirroring METABRIC / TCGA handling.
"""
from __future__ import annotations

import os
import sys
import warnings
from pathlib import Path


CANCER_LABEL = {"BRCA": "IDC", "LUAD": "LUAD", "PAAD": "PAAD",
                "PRAD": "PRAD", "CRC": "CRC"}

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib-cache")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from figure_style import apply_arial_style
apply_arial_style()
import numpy as np
import pandas as pd
import yaml
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test
from lifelines.utils import concordance_index

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from build_3flag_enhanced import (
    ALIAS_MAP, ALL_PW_KEYS, CHROMATIN_ALL, CHROMATIN_OG, CHROMATIN_TSG,
    DDR_TSG, NONSILENT, PATHWAY_XLSX, parse_10pathway,
)

RAW_BASE = ROOT / "external_validation" / "raw_data" / "GENIE_BPC"
OUT_BASE = ROOT / "external_validation" / "outputs"
FORMULA = ROOT / "external_validation" / "zscore_params" / "compact_score_formula.yaml"
MSK_ZSCORE_DIR = ROOT / "build_3flag_outputs" / "pathonly_3flag_enhanced"

# GENIE BPC cohort → (cancer code, ONCOTREE codes, HIT_COUNT pathways)
COHORTS = {
    "BrCa": {
        "cancer": "BRCA",
        "oncotree": ["IDC"],
        "hit_count_pw": ["TP53", "PI3K", "RTK_RAS", "Cell_Cycle", "Chromatin", "DDR", "MYC"],
        "has_sv": True,
    },
    "NSCLC": {
        "cancer": "LUAD",
        "oncotree": ["LUAD"],
        "hit_count_pw": ["TP53", "PI3K", "RTK_RAS", "Cell_Cycle", "Chromatin", "DDR", "MYC", "NOTCH", "NRF2"],
    },
    "CRC": {
        "cancer": "CRC",
        "oncotree": ["COAD", "READ", "COADREAD"],
        "hit_count_pw": ["TP53", "RTK_RAS", "WNT", "TGF_Beta", "PI3K", "DDR", "MYC", "Chromatin"],
    },
    "PANC": {
        "cancer": "PAAD",
        "oncotree": ["PAAD"],
        "hit_count_pw": ["TP53", "RTK_RAS", "Cell_Cycle", "Chromatin", "TGF_Beta", "MYC", "NOTCH", "HIPPO"],
    },
    "Prostate": {
        "cancer": "PRAD",
        "oncotree": ["PRAD"],
        "hit_count_pw": ["TP53", "PI3K", "Cell_Cycle", "Chromatin", "DDR", "WNT"],
    },
}


def clean_gene(g):
    if pd.isna(g):
        return None
    s = str(g).strip()
    if not s or s.lower() == "nan":
        return None
    return ALIAS_MAP.get(s, s)


def is_nonsilent(v):
    if pd.isna(v):
        return False
    return str(v).strip() in NONSILENT


def parse_event(s):
    if pd.isna(s):
        return np.nan
    s = str(s).upper()
    if s.startswith("1") or "DECEASED" in s or "DEAD" in s:
        return 1.0
    if s.startswith("0") or "LIVING" in s or "ALIVE" in s:
        return 0.0
    return np.nan


def load_msk_zscore(cancer):
    df = pd.read_csv(MSK_ZSCORE_DIR / cancer / f"{cancer}_zscore_params.csv")
    return {row["col"]: (float(row["mean"]), float(row["std"])) for _, row in df.iterrows()}


def zscore_msk(s, params):
    mu, sd = params
    sd = sd if sd >= 1e-9 else 1.0
    return (s - mu) / sd


def pathway_defs():
    pw = parse_10pathway(PATHWAY_XLSX)
    pw["DDR"] = {"OG": set(), "TSG": DDR_TSG, "all": DDR_TSG}
    pw["Chromatin"] = {"OG": CHROMATIN_OG, "TSG": CHROMATIN_TSG, "all": CHROMATIN_ALL}
    return pw


def build_features(cohort_key, cfg, pw):
    raw = RAW_BASE / cohort_key
    samp = pd.read_csv(raw / "clinical_sample_nonMSK.csv", low_memory=False)
    surv = pd.read_csv(raw / "survival_nonMSK.csv", low_memory=False)
    cna = pd.read_csv(raw / "data_cna.tsv", sep="\t", low_memory=False)
    mut = pd.read_csv(raw / "data_mutations.tsv", sep="\t", low_memory=False)

    surv["OS_MONTHS"] = pd.to_numeric(surv["OS_MONTHS"], errors="coerce")
    surv["Event_OS"] = surv["OS_STATUS"].apply(parse_event)
    surv = surv[(surv["OS_MONTHS"] > 0) & surv["Event_OS"].notna()].copy()

    samp = samp[samp["ONCOTREE_CODE"].isin(cfg["oncotree"])].copy()
    samp = samp.merge(surv[["PATIENT_ID", "OS_MONTHS", "Event_OS"]], on="PATIENT_ID", how="inner")
    print(f"  After ONCOTREE filter + OS-evaluable: n={len(samp)}, "
          f"events={int(samp['Event_OS'].sum())}")

    # CNA matrix
    gene_col = cna.columns[0]
    cna[gene_col] = cna[gene_col].map(clean_gene)
    cna = cna.dropna(subset=[gene_col]).drop_duplicates(gene_col).set_index(gene_col)
    cna = cna.apply(pd.to_numeric, errors="coerce").fillna(0).astype(float)
    cna_samples = set(cna.columns)
    samp = samp[samp["SAMPLE_ID"].isin(cna_samples)].copy()
    print(f"  After CNA coverage: n={len(samp)}, events={int(samp['Event_OS'].sum())}")

    sample_ids = samp["SAMPLE_ID"].tolist()
    cna_sub = cna.reindex(columns=sample_ids, fill_value=0)
    covered_genes = set(cna_sub.index)

    # Mutations
    mut["Hugo_Symbol"] = mut["Hugo_Symbol"].map(clean_gene)
    mut = mut.dropna(subset=["Hugo_Symbol", "Tumor_Sample_Barcode"])
    mut = mut.rename(columns={"Tumor_Sample_Barcode": "SAMPLE_ID"})
    mut = mut[mut["Variant_Classification"].apply(is_nonsilent)]
    mut_pairs = set(zip(mut["SAMPLE_ID"], mut["Hugo_Symbol"]))

    # SV (only some cohorts)
    sv_pairs = set()
    sv_path = raw / "data_sv.tsv"
    if cfg.get("has_sv") and sv_path.exists():
        sv = pd.read_csv(sv_path, sep="\t", low_memory=False)
        sid_col = "Sample_Id" if "Sample_Id" in sv.columns else (
            "Sample_ID" if "Sample_ID" in sv.columns else "Tumor_Sample_Barcode")
        for col in ("Site1_Hugo_Symbol", "Site2_Hugo_Symbol"):
            if col in sv.columns:
                pairs = sv[[sid_col, col]].dropna()
                pairs[col] = pairs[col].map(clean_gene)
                pairs = pairs.dropna()
                sv_pairs |= set(zip(pairs[sid_col].astype(str), pairs[col].astype(str)))
        print(f"  SV gene-pair set: {len(sv_pairs)} (samples-with-SV: "
              f"{len({s for s, _ in sv_pairs})})")

    msk_z = load_msk_zscore(cfg["cancer"])

    out = pd.DataFrame(index=sample_ids)
    pathway_binary = {}
    for pathway in ALL_PW_KEYS:
        genes = sorted(set(pw[pathway]["all"]) & covered_genes)
        denom = max(len(genes), 1)
        if genes:
            vals = cna_sub.loc[genes]
            amp = vals.ge(2).any(axis=0)
            dele = vals.le(-2).any(axis=0)
        else:
            amp = pd.Series(False, index=sample_ids)
            dele = pd.Series(False, index=sample_ids)
        mut_count, any_count = [], []
        for sid in sample_ids:
            mgenes = {g for g in genes if (sid, g) in mut_pairs}
            if genes:
                amp_genes = set(cna_sub.index[cna_sub[sid].ge(2)]) & set(genes)
                del_genes = set(cna_sub.index[cna_sub[sid].le(-2)]) & set(genes)
            else:
                amp_genes, del_genes = set(), set()
            mut_count.append(len(mgenes) / denom)
            any_count.append(len(mgenes | amp_genes | del_genes) / denom)
        pk = f"PW_{pathway}"
        out[f"{pk}_amp_hit"] = amp.astype(int).values
        out[f"{pk}_del_hit"] = dele.astype(int).values
        # SV hit: any of this pathway's genes appears as Site1 or Site2 in SV for this sample
        if sv_pairs and genes:
            gene_set = set(genes)
            out[f"{pk}_sv_hit"] = [
                int(any((sid, g) in sv_pairs for g in gene_set)) for sid in sample_ids
            ]
        else:
            out[f"{pk}_sv_hit"] = 0
        out[f"{pk}_mut_rate_raw"] = mut_count
        out[f"{pk}_any_rate_raw"] = any_count
        pathway_binary[pathway] = (
            out[f"{pk}_amp_hit"].astype(bool)
            | out[f"{pk}_del_hit"].astype(bool)
            | (out[f"{pk}_mut_rate_raw"] > 0)
        ).astype(int)

    for pathway in ALL_PW_KEYS:
        pk = f"PW_{pathway}"
        out[f"{pk}_mut_rate_z"] = zscore_msk(out[f"{pk}_mut_rate_raw"],
                                             msk_z.get(f"{pk}_mut_rate", (0.0, 1.0)))
        out[f"{pk}_any_rate_z"] = zscore_msk(out[f"{pk}_any_rate_raw"],
                                             msk_z.get(f"{pk}_any_rate", (0.0, 1.0)))
        out[f"{pk}_zsum"] = (out[f"{pk}_amp_hit"] + out[f"{pk}_del_hit"]
                            + out[f"{pk}_sv_hit"] + out[f"{pk}_mut_rate_z"])

    out["PW_HIT_COUNT"] = pd.DataFrame(
        {p: pathway_binary[p] for p in cfg["hit_count_pw"]}
    ).sum(axis=1)
    out = out.reset_index(names="SAMPLE_ID")
    out = out.merge(samp[["PATIENT_ID", "SAMPLE_ID", "OS_MONTHS", "Event_OS",
                          "ONCOTREE_CODE", "SAMPLE_TYPE_DETAILED",
                          "AGE_AT_SEQUENCING", "INSTITUTION", "SEQ_ASSAY_ID"]],
                    on="SAMPLE_ID")
    out["Event_OS"] = out["Event_OS"].astype(int)
    return out


def compute_score(df, cancer, formula):
    spec = formula["cancers"][cancer]
    feats = [f["name"] for f in spec["features"]]
    betas = np.array([float(f["beta"]) for f in spec["features"]])
    return (df[feats].values * betas).sum(axis=1)


def km_pvalue(df, group):
    use = df[["OS_MONTHS", "Event_OS"]].copy()
    use["group"] = group.values
    use = use.dropna()
    if use["group"].nunique() < 2:
        return float("nan")
    if use["group"].nunique() == 2:
        vals = list(use["group"].unique())
        a = use[use["group"] == vals[0]]; b = use[use["group"] == vals[1]]
        return float(logrank_test(a["OS_MONTHS"], b["OS_MONTHS"], a["Event_OS"], b["Event_OS"]).p_value)
    return float(multivariate_logrank_test(use["OS_MONTHS"], use["group"], use["Event_OS"]).p_value)


def plot_km(ax, df, group, title):
    kmf = KaplanMeierFitter()
    for name in sorted(group.dropna().unique()):
        mask = group == name
        if mask.sum() < 3:
            continue
        kmf.fit(df.loc[mask, "OS_MONTHS"], df.loc[mask, "Event_OS"],
                label=f"{name} (n={mask.sum()})")
        kmf.plot_survival_function(ax=ax, ci_show=False, linewidth=1.6)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Months"); ax.set_ylabel("OS")
    ax.set_xlim(0, min(120, max(36, df["OS_MONTHS"].quantile(0.95))))
    ax.grid(alpha=0.2); ax.legend(fontsize=7)


def calibration_at(df, score, t):
    """Stratify by score quintile, compute predicted survival vs observed at t."""
    from lifelines import KaplanMeierFitter
    df = df.copy(); df["score"] = score
    df = df.dropna(subset=["OS_MONTHS", "Event_OS", "score"])
    if df["score"].nunique() < 5:
        return pd.DataFrame()
    df["q"] = pd.qcut(df["score"].rank(method="first"), 5,
                      labels=["Q1", "Q2", "Q3", "Q4", "Q5"])
    rows = []
    kmf = KaplanMeierFitter()
    # We use a Cox model on score alone for predicted survival proxy
    cph = CoxPHFitter(penalizer=0.1)
    cph.fit(df[["OS_MONTHS", "Event_OS", "score"]], duration_col="OS_MONTHS", event_col="Event_OS")
    for q in ["Q1", "Q2", "Q3", "Q4", "Q5"]:
        sub = df[df["q"] == q]
        if len(sub) < 5:
            continue
        kmf.fit(sub["OS_MONTHS"], sub["Event_OS"])
        observed = 1 - float(kmf.predict(t))
        # predicted event prob at t from cox using mean score in quintile
        s_mean = sub["score"].mean()
        sf = cph.predict_survival_function(pd.DataFrame({"score": [s_mean]}))
        # find S(t)
        if t in sf.index:
            S_t = float(sf.loc[t].iloc[0])
        else:
            # interpolate
            S_t = float(np.interp(t, sf.index.values, sf.iloc[:, 0].values,
                                   left=1.0, right=float(sf.iloc[-1, 0])))
        predicted = 1 - S_t
        rows.append({"q": q, "n": int(len(sub)),
                     "mean_predicted_event": predicted,
                     "observed_event": observed})
    return pd.DataFrame(rows)


def run_cohort(cohort_key, formula, pw):
    cfg = COHORTS[cohort_key]
    cancer = cfg["cancer"]
    print(f"\n===== GENIE BPC {cohort_key} → {cancer} =====")

    out_dir = OUT_BASE / f"GENIE_BPC_{cohort_key}"
    fig_dir = out_dir / "figures"; fig_dir.mkdir(parents=True, exist_ok=True)
    tab_dir = out_dir / "tables";  tab_dir.mkdir(parents=True, exist_ok=True)

    feat = build_features(cohort_key, cfg, pw)
    feat["risk_score"] = compute_score(feat, cancer, formula)
    feat.to_csv(out_dir / f"GENIE_BPC_{cohort_key}_pathway_features.csv", index=False)

    cidx = float(concordance_index(feat["OS_MONTHS"], -feat["risk_score"], feat["Event_OS"]))
    print(f"  C-index = {cidx:.4f}")

    # quantile KM
    summary_rows = []
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), constrained_layout=True)
    for i, q in enumerate([2, 3, 4]):
        labels = [f"Q{j}" for j in range(1, q + 1)]
        group = pd.qcut(feat["risk_score"].rank(method="first"), q=q, labels=labels)
        p = km_pvalue(feat, group)
        summary_rows.append({
            "cohort": f"GENIE_BPC_{cohort_key}", "cancer": cancer,
            "n": len(feat), "events": int(feat["Event_OS"].sum()),
            "quantiles": q, "logrank_p": p, "external_cindex": cidx,
        })
        plot_km(axes[i], feat, group,
                f"{q}Q  n={len(feat)} ev={int(feat['Event_OS'].sum())}  p={p:.2e}")
    fig.suptitle(f"GENIE BPC {cohort_key} (non-MSK, n={len(feat)}) — {CANCER_LABEL.get(cancer, cancer)} compact score quantile KM",
                 fontsize=11)
    fig.savefig(fig_dir / f"genie_bpc_{cohort_key}_quantile_km.png", dpi=220)
    plt.close(fig)

    # fixed MSK cutoff
    cutoff = float(formula["cancers"][cancer]["msk_partial_hazard_cutoff"])
    feat["risk_score_exp"] = np.exp(feat["risk_score"])
    feat["risk_high_msk"] = (feat["risk_score_exp"] > cutoff).astype(int)
    p_msk = km_pvalue(feat, feat["risk_high_msk"].map({0: "Low", 1: "High"}))
    p_med = km_pvalue(feat, (feat["risk_score"] > feat["risk_score"].median()).astype(int).map({0: "Low", 1: "High"}))
    n_high = int(feat["risk_high_msk"].sum())
    print(f"  MSK cutoff KM p={p_msk:.3e} (n_high={n_high}/{len(feat)})  median p={p_med:.3e}")
    summary_rows.append({"cohort": f"GENIE_BPC_{cohort_key}", "cancer": cancer,
                         "n": len(feat), "events": int(feat["Event_OS"].sum()),
                         "quantiles": "MSK_cutoff", "logrank_p": p_msk, "external_cindex": cidx,
                         "n_high_msk": n_high})
    summary_rows.append({"cohort": f"GENIE_BPC_{cohort_key}", "cancer": cancer,
                         "n": len(feat), "events": int(feat["Event_OS"].sum()),
                         "quantiles": "median", "logrank_p": p_med, "external_cindex": cidx})

    # high/low KM figure
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    plot_km(axes[0], feat, feat["risk_high_msk"].map({0: "Low", 1: "High"}),
            f"MSK β-cutoff (>{cutoff:.3f})  p={p_msk:.2e}")
    plot_km(axes[1], feat,
            (feat["risk_score"] > feat["risk_score"].median()).astype(int).map({0: "Low", 1: "High"}),
            f"GENIE median cutoff  p={p_med:.2e}")
    fig.suptitle(f"GENIE BPC {cohort_key} ({CANCER_LABEL.get(cancer, cancer)}): compact score high/low risk KM", fontsize=11)
    fig.savefig(fig_dir / f"genie_bpc_{cohort_key}_high_low_km.png", dpi=220)
    plt.close(fig)

    # Calibration at 24m / 36m
    cal_t = 24 if cancer in ("LUAD", "PAAD") else 36
    cal = calibration_at(feat, feat["risk_score"], cal_t)
    if not cal.empty:
        cal["t_month"] = cal_t
        cal["cancer"] = cancer
        cal["cohort"] = f"GENIE_BPC_{cohort_key}"
        cal.to_csv(tab_dir / f"calibration_{cal_t}m.csv", index=False)
        fig, ax = plt.subplots(figsize=(4.5, 4.5), constrained_layout=True)
        ax.plot([0, 1], [0, 1], "k--", linewidth=0.8)
        ax.scatter(cal["mean_predicted_event"], cal["observed_event"], s=50, color="#234")
        for _, row in cal.iterrows():
            ax.annotate(row["q"], (row["mean_predicted_event"], row["observed_event"]),
                        fontsize=8, xytext=(4, 4), textcoords="offset points")
        ax.set_xlabel(f"Mean predicted event probability at {cal_t}m")
        ax.set_ylabel(f"Observed event probability at {cal_t}m")
        lim = max(cal["mean_predicted_event"].max(), cal["observed_event"].max()) + 0.05
        ax.set_xlim(0, lim); ax.set_ylim(0, lim)
        ax.set_title(f"{cohort_key} ({CANCER_LABEL.get(cancer, cancer)}) calibration @ {cal_t}m\n(5 quintiles of compact score)",
                     fontsize=10)
        fig.savefig(fig_dir / f"calibration_{cal_t}m.png", dpi=220)
        plt.close(fig)

    pd.DataFrame(summary_rows).to_csv(tab_dir / "compact_score_summary.csv", index=False)
    return summary_rows, cal


def main():
    with FORMULA.open() as f:
        formula = yaml.safe_load(f)
    pw = pathway_defs()
    all_summary = []
    all_cal = []
    for ck in COHORTS:
        rows, cal = run_cohort(ck, formula, pw)
        all_summary.extend(rows)
        if not cal.empty:
            all_cal.append(cal)
    all_summary_df = pd.DataFrame(all_summary)
    out_dir = OUT_BASE / "_GENIE_BPC_combined"
    out_dir.mkdir(parents=True, exist_ok=True)
    all_summary_df.to_csv(out_dir / "compact_score_summary.csv", index=False)
    if all_cal:
        pd.concat(all_cal, ignore_index=True).to_csv(out_dir / "calibration_combined.csv", index=False)
    print(f"\n=== Combined summary ===")
    print(all_summary_df[all_summary_df["quantiles"].isin([2, "MSK_cutoff"])].to_string(index=False))


if __name__ == "__main__":
    main()
