#!/usr/bin/env python3
"""
gene_vs_pathway_external.py
===========================
No-refit EXTERNAL head-to-head on GENIE BPC (institutionally independent).

The internal benchmark (gene_vs_pathway_headtohead.py) showed gene-level only
edges the pathway compact score by spending 7-25x more features. This script
tests the decisive question: does that gene-level edge SURVIVE transfer to an
independent cohort, or is it MSK-specific overfitting?

Protocol (mirrors the published compact-score external validation exactly):
  * Fit each gene-level Coxnet on the MSK training split ONLY, then FREEZE
    (coefficients + StandardScaler mean/sd + selected feature columns).
  * On each GENIE BPC cohort, rebuild the SAME per-(gene, alteration) events
    (mut from data_mutations, amp/del from data_cna; SV set to 0 for all
    non-MSK panels, exactly as the compact score does), align to the frozen
    feature columns (genes absent on the external panel -> 0), apply the frozen
    scaler + coefficients. No refitting, no re-scaling.
  * Compact pathway score external risk is taken from the already-validated
    external_validation/outputs/GENIE_BPC_*/*_pathway_features.csv (risk_score).
  * Report external Harrell C-index and Uno C-IPCW, plus paired bootstrap
    Delta-C (compact - gene) on the external cohort.

Run:
    python3 scripts/gene_vs_pathway_external.py --boot 1000
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sksurv.metrics import concordance_index_ipcw
from sksurv.util import Surv
from lifelines.utils import concordance_index

import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
import gene_vs_pathway_headtohead as HH

RAW_BASE = ROOT / "external_validation" / "raw_data" / "GENIE_BPC"
OUT_BASE = ROOT / "external_validation" / "outputs"
OUT_DIR = ROOT / "figures_tables" / "gene_vs_pathway_external"

# cancer -> (GENIE raw dir, GENIE output cohort dir, ONCOTREE codes)
COHORTS = {
    "BRCA": ("BrCa", "GENIE_BPC_BrCa", ["IDC"]),
    "LUAD": ("NSCLC", "GENIE_BPC_NSCLC", ["LUAD"]),
    "PAAD": ("PANC", "GENIE_BPC_PANC", ["PAAD"]),
    "PRAD": ("Prostate", "GENIE_BPC_Prostate", ["PRAD"]),
    "CRC":  ("CRC", "GENIE_BPC_CRC", ["COAD", "READ", "COADREAD"]),
}
NONSILENT = {
    "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
    "In_Frame_Ins", "In_Frame_Del", "Splice_Site", "Translation_Start_Site",
    "Nonstop_Mutation", "Splice_Region",
}
TAU = 24.0


# --------------------------------------------------------------- MSK: fit + freeze
def fit_frozen_gene_model(cancer, gene_set, cap_k=None):
    """Fit a gene-level Coxnet on MSK TRAIN only; return frozen (model, scaler, cols)."""
    train, _test, tr_sids, te_sids = HH.load_split(cancer)
    ev = HH.build_gene_events(cancer, tr_sids + te_sids, gene_set)
    prev = ev.loc[tr_sids].mean(0)
    cols = prev[prev >= HH.PREVALENCE].index.tolist()
    Xtr = ev.loc[tr_sids, cols].reset_index(drop=True)
    scaler = StandardScaler().fit(Xtr.values)
    y_tr = HH.make_surv(train)
    model, n_used = HH.coxnet_select(scaler.transform(Xtr.values), y_tr, cap_k=cap_k)
    return {"model": model, "scaler": scaler, "cols": cols, "n_used": n_used}


# --------------------------------------------------------------- GENIE: gene events
def build_genie_gene_events(raw_dir, sample_ids):
    """Per-(gene, alteration) binary events for GENIE samples. SV set to 0 (non-MSK
    panels lack callable SV, matching the compact-score external protocol)."""
    raw = RAW_BASE / raw_dir
    mut = pd.read_csv(raw / "data_mutations.tsv", sep="\t", low_memory=False)
    mut = mut.rename(columns={"Tumor_Sample_Barcode": "SAMPLE_ID"})
    mut = mut[mut["Variant_Classification"].isin(NONSILENT)]
    mut = mut[mut["SAMPLE_ID"].isin(sample_ids)]
    idx = pd.Index(sample_ids, name="SAMPLE_ID")
    feat = pd.DataFrame(index=idx)
    # mutations
    for g, sub in mut.groupby("Hugo_Symbol"):
        feat[f"g_{g}_mut"] = idx.isin(sub["SAMPLE_ID"].unique()).astype(int)
    # CNA amp/del
    cna = pd.read_csv(raw / "data_cna.tsv", sep="\t", low_memory=False)
    gcol = cna.columns[0]
    cna = cna.dropna(subset=[gcol]).drop_duplicates(gcol).set_index(gcol)
    keep = [s for s in sample_ids if s in cna.columns]
    cna = cna[keep].apply(pd.to_numeric, errors="coerce")
    amp = (cna >= 2)
    dele = (cna <= -2)
    for g in cna.index:
        a = amp.loc[g]
        if a.any():
            feat[f"g_{g}_amp"] = idx.isin(a.index[a]).astype(int)
        d = dele.loc[g]
        if d.any():
            feat[f"g_{g}_del"] = idx.isin(d.index[d]).astype(int)
    return feat.fillna(0)


def apply_frozen(frozen, genie_events):
    X = genie_events.reindex(columns=frozen["cols"], fill_value=0)
    Xs = frozen["scaler"].transform(X.values)
    return frozen["model"].predict(Xs)     # risk (higher = worse)


# --------------------------------------------------------------- metrics
def ext_metrics(dur, evt, risk):
    dur = np.asarray(dur, float); evt = np.asarray(evt, bool); risk = np.asarray(risk, float)
    c = float(concordance_index(dur, -risk, evt))
    y = Surv.from_arrays(event=evt, time=dur)
    tau = min(TAU, float(dur[evt].max()) if evt.any() else TAU)
    try:
        cip = float(concordance_index_ipcw(y, y, risk, tau=tau)[0])
    except Exception:
        cip = float("nan")
    return c, cip


def boot_delta(dur, evt, r_ref, r_oth, B, rng):
    dur = np.asarray(dur, float); evt = np.asarray(evt, bool)
    n = len(dur); d = []
    for _ in range(B):
        idx = rng.integers(0, n, n)
        if evt[idx].sum() < 2:
            continue
        cr = concordance_index(dur[idx], -r_ref[idx], evt[idx])
        co = concordance_index(dur[idx], -r_oth[idx], evt[idx])
        d.append(cr - co)
    d = np.array(d)
    p = 2 * min((d < 0).mean(), (d > 0).mean())
    return float(d.mean()), float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5)), float(min(p, 1))


# --------------------------------------------------------------- driver
def run_cancer(cancer, B, rng):
    raw_dir, out_dir, _onc = COHORTS[cancer]
    pw = pd.read_csv(OUT_BASE / out_dir / f"{out_dir}_pathway_features.csv")
    pw = pw[(pw["OS_MONTHS"] > 0) & pw["Event_OS"].notna()].copy()
    sids = pw["SAMPLE_ID"].tolist()
    dur = pw["OS_MONTHS"].to_numpy(float)
    evt = pw["Event_OS"].astype(bool).to_numpy()
    print(f"\n=== {cancer} (GENIE {out_dir}) n={len(sids)} events={int(evt.sum())} ===", flush=True)

    genie_ev = build_genie_gene_events(raw_dir, sids)

    gs = HH.retained_gene_set(cancer)
    K = int((pd.read_csv(HH.DECISION)
             .pipe(lambda d: d[d["cancer"] == HH.CANCER_ALIAS.get(cancer, cancer)])).shape[0])
    frozen = {
        "gene_restricted": fit_frozen_gene_model(cancer, gs),
        "gene_unrestricted": fit_frozen_gene_model(cancer, None),
        "gene_capped_K": fit_frozen_gene_model(cancer, gs, cap_k=K),
    }

    risks = {"compact_fixed": pw["risk_score"].to_numpy(float)}
    for name, fr in frozen.items():
        risks[name] = apply_frozen(fr, genie_ev)

    rows = []
    for name, r in risks.items():
        c, cip = ext_metrics(dur, evt, r)
        nf = K if name == "compact_fixed" else frozen[name]["n_used"]
        rows.append({"cancer": cancer, "cohort": out_dir, "arm": name,
                     "n_features_msk": nf, "n_ext": len(sids), "events_ext": int(evt.sum()),
                     "ext_c_index": round(c, 4), "ext_c_ipcw": round(cip, 4)})
        print(f"  {name:18s} nf(MSK)={nf:3d}  extC={c:.4f}  extC-IPCW={cip:.4f}", flush=True)

    boots = []
    ref = risks["compact_fixed"]
    for name, r in risks.items():
        if name == "compact_fixed":
            continue
        dc, lo, hi, p = boot_delta(dur, evt, ref, r, B, rng)
        boots.append({"cancer": cancer, "comparison": f"compact_fixed_minus_{name}",
                      "ext_delta_c_mean": round(dc, 4), "ci_lo": round(lo, 4),
                      "ci_hi": round(hi, 4), "boot_p": round(p, 4)})
        sig = "*" if p < 0.05 else " "
        print(f"    ext ΔC vs {name:18s} = {dc:+.4f} [{lo:+.4f},{hi:+.4f}] p={p:.3f}{sig}", flush=True)
    return rows, boots


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--boot", type=int, default=1000)
    ap.add_argument("--cancers", nargs="+", default=list(COHORTS))
    args = ap.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(42)
    all_m, all_b = [], []
    for c in args.cancers:
        m, b = run_cancer(c, args.boot, rng)
        all_m += m; all_b += b
    mdf = pd.DataFrame(all_m); bdf = pd.DataFrame(all_b)
    mdf.to_csv(OUT_DIR / "external_metrics.csv", index=False)
    bdf.to_csv(OUT_DIR / "external_bootstrap_deltaC.csv", index=False)
    json.dump({"metrics": all_m, "bootstrap": all_b},
              open(OUT_DIR / "external.json", "w"), indent=2)
    print("\n" + "=" * 70)
    print(mdf.to_string(index=False))
    print(f"\nwrote {OUT_DIR}")


if __name__ == "__main__":
    main()
