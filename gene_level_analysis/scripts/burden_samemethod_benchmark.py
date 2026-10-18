#!/usr/bin/env python3
"""
burden_samemethod_benchmark.py
==============================
Two reviewer-style controls for the representation claim (BiB revision, 2026-10).

1. Burden-only comparator
   Does the 4-5-feature compact score carry information beyond the overall
   pathway-alteration burden?  Arm "burden_only" = Cox on PW_HIT_COUNT (All-pathway
   alteration count) alone, fitted on MSK train; compared with the fixed compact score.

2. Same-selection-method comparator
   Table 2 compares the IBCGA-selected compact score with Elastic-Net gene models, so a
   difference can come from the selection method as well as the representation.  Here the
   pathway pool is put through the SAME Coxnet routine (gene_vs_pathway_headtohead.
   coxnet_select, same l1 ratio, lambda path, 5-fold CV, prevalence filter for genes) as the
   gene indicators, at the matched budget (<= Kc terms) and without a budget cap.

Every model is fitted on the MSK training split only, frozen (scaler + coefficients +
columns) and applied unchanged to the MSK held-out test split and to the GENIE BPC external
cohorts (genes absent from an external panel -> 0, SV -> 0, exactly as in Table 2).
PW_DRIVER_BURDEN_z is excluded from the pathway pool because it is not reconstructed for
GENIE BPC; no compact score uses it.

Outputs (figures_tables/burden_samemethod_benchmark/):
  arm_metrics.csv      C-index (bootstrap 95% CI) per cancer x arm x cohort
  paired_deltaC.csv    paired bootstrap ΔC for the planned contrasts
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
import gene_vs_pathway_external as EX  # noqa: E402
import gene_vs_pathway_headtohead as HH  # noqa: E402
import benchmark_common_coverage as BC  # noqa: E402  (Table 2 gene-candidate restriction)
sys.path.insert(0, str(ROOT / "external_validation" / "scripts"))
COVERAGE = "observed"   # "observed" = Table 2 rule; "panel" = assay gene lists (set in main)
V4_DIR = ROOT / "external_validation" / "outputs" / "panel_coverage_V4"


def panel_common_genes(raw_dir, sample_ids, min_frac=0.01):
    """Genes on EVERY assay used by >= min_frac of the evaluated GENIE samples."""
    import genie_panel_coverage_sensitivity as PC
    panels = PC.load_panels()
    samp = pd.read_csv(EX.RAW_BASE / raw_dir / "clinical_sample_nonMSK.csv", low_memory=False)
    a = samp.set_index("SAMPLE_ID")["SEQ_ASSAY_ID"].reindex(sample_ids).dropna()
    freq = a.value_counts(normalize=True)
    used = [x for x in freq.index if freq[x] >= min_frac]
    common = set.intersection(*(panels[x] for x in used))
    return common, used
from compact_score.apply_compact_score import apply_score, load_formula  # noqa: E402

OUT = ROOT / "figures_tables" / "burden_samemethod_benchmark"
DISP = {"BRCA": "IDC", "LUAD": "LUAD", "PAAD": "PAAD", "PRAD": "PRAD", "CRC": "CRC"}
B = 1000
EXCLUDE = {"PW_DRIVER_BURDEN_z"}
CONTRASTS = [("compact_fixed", "burden_only"),
             ("pathway_EN_capK", "gene_EN_capK"),
             ("pathway_EN_free", "gene_EN_free"),
             ("compact_fixed", "pathway_EN_capK")]


def boot_ci(t, e, r, rng):
    n = len(t); v = []
    for _ in range(B):
        i = rng.integers(0, n, n)
        if e[i].sum() < 2:
            continue
        v.append(concordance_index(t[i], -r[i], e[i]))
    v = np.array(v)
    obs = float(concordance_index(t, -r, e))
    return obs, float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5)), float(v.mean())


def paired(t, e, ra, rb, rng):
    n = len(t); d = []
    for _ in range(B):
        i = rng.integers(0, n, n)
        if e[i].sum() < 2:
            continue
        d.append(concordance_index(t[i], -ra[i], e[i]) - concordance_index(t[i], -rb[i], e[i]))
    d = np.array(d)
    p = min(1.0, 2 * min((d <= 0).mean(), (d >= 0).mean()))
    obs = float(concordance_index(t, -ra, e) - concordance_index(t, -rb, e))
    return obs, float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5)), float(p), float(d.mean())


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); prev = 1.0
    for rank, i in reversed(list(enumerate(o, 1))):
        prev = min(prev, p[i] * n / rank); q[i] = prev
    return q


def frozen_coxnet(Xtr, Xte, Xext, y_tr, cap):
    """Same rule as Table 2: under a cap, the CV-best eligible model ('cv_best');
    without a cap, the CV-chosen lambda."""
    sc = StandardScaler().fit(Xtr.values)
    m, nu = HH.coxnet_select(sc.transform(Xtr.values), y_tr, cap_k=cap,
                             cap_mode="cv_best" if cap is not None else "densest")
    pr = lambda X: m.predict(sc.transform(X.values))
    return pr(Xte), pr(Xext), nu


def run_cancer(cancer, rng):
    raw_dir, out_dir, _ = EX.COHORTS[cancer]
    key = HH.CANCER_ALIAS.get(cancer, cancer)
    formula = load_formula()
    K = len(formula["cancers"][formula.get("cancer_key_aliases", {}).get(key, key)]["features"])
    train, test, tr, te = HH.load_split(cancer)
    y_tr = HH.make_surv(train)
    t_i = test["OS_MONTHS"].to_numpy(float); e_i = test["Event_OS"].astype(bool).to_numpy()
    ext = pd.read_csv(EX.OUT_BASE / out_dir / f"{out_dir}_pathway_features.csv")
    ext = ext[(ext["OS_MONTHS"] > 0) & ext["Event_OS"].notna()].copy()
    if COVERAGE == "panel":   # per-sample panel-coverage features, same samples and order
        v4 = pd.read_csv(V4_DIR / f"{out_dir}_pathway_features_V4.csv").set_index("SAMPLE_ID")
        keep = ["SAMPLE_ID", "OS_MONTHS", "Event_OS"]
        ext = ext[keep].join(v4.drop(columns=[c for c in keep if c in v4.columns and c != "SAMPLE_ID"],
                                     errors="ignore"), on="SAMPLE_ID")
    t_e = ext["OS_MONTHS"].to_numpy(float); e_e = ext["Event_OS"].astype(bool).to_numpy()

    risks, nterms = {}, {}
    # compact (fixed published formula); sanity: external risk_score == formula applied to ext
    ci, _ = apply_score(test, key, formula)     # log score (monotone in exp risk; same C-index)
    ce, _ = apply_score(ext, key, formula)
    assert np.allclose(ce.to_numpy(float), ext["risk_score"].to_numpy(float), atol=1e-6), "external z/feature mismatch"
    risks["compact_fixed"] = (ci.to_numpy(float), ext["risk_score"].to_numpy(float)); nterms["compact_fixed"] = K
    # burden only: univariable Cox on All-pathway alteration count (MSK train)
    cph = CoxPHFitter().fit(train[["PW_HIT_COUNT", "OS_MONTHS", "Event_OS"]], "OS_MONTHS", "Event_OS")
    b = float(cph.params_["PW_HIT_COUNT"])
    risks["burden_only"] = (b * test["PW_HIT_COUNT"].to_numpy(float), b * ext["PW_HIT_COUNT"].to_numpy(float))
    nterms["burden_only"] = 1
    # pathway pool through the same Coxnet routine as the gene models
    pool = [c for c in train.columns if c.startswith("PW_") and c not in EXCLUDE]
    Xtr, Xte, Xex = train[pool], test[pool], ext.reindex(columns=pool)
    assert not Xex.isna().any().any(), "missing pathway columns in GENIE"
    for name, cap in (("pathway_EN_capK", K), ("pathway_EN_free", None)):
        ri, re, nu = frozen_coxnet(Xtr, Xte, Xex, y_tr, cap)
        risks[name] = (ri, re); nterms[name] = nu
    # gene indicators (identical construction to Table 2)
    gs = HH.retained_gene_set(cancer)
    ev = HH.build_gene_events(cancer, tr + te, gs)
    prev = ev.loc[tr].mean(0); cand = prev[prev >= HH.PREVALENCE].index.tolist()
    if COVERAGE == "panel":
        common, used = panel_common_genes(raw_dir, ext["SAMPLE_ID"].tolist())
        mut_genes = cna_genes = common
        print(f"  panel-common genes: {len(common)} (assays: {', '.join(used)})", flush=True)
    else:
        mut_genes, cna_genes = BC.genie_coverage(raw_dir)      # common coverage, as in Table 2
    cols = [c for c in cand if BC.evaluable(c, mut_genes, cna_genes)]
    print(f"  gene candidates: {len(cand)} -> {len(cols)} after coverage restriction", flush=True)
    gev = EX.build_genie_gene_events(raw_dir, ext["SAMPLE_ID"].tolist()).reindex(columns=cols, fill_value=0)
    for name, cap in (("gene_EN_capK", K), ("gene_EN_free", None)):
        ri, re, nu = frozen_coxnet(ev.loc[tr, cols], ev.loc[te, cols], gev, y_tr, cap)
        risks[name] = (ri, re); nterms[name] = nu

    rows, deltas = [], []
    for arm, (ri, re) in risks.items():
        im, il, ih, ibm = boot_ci(t_i, e_i, np.asarray(ri, float), rng)
        em, el, eh, ebm = boot_ci(t_e, e_e, np.asarray(re, float), rng)
        rows.append({"cancer": DISP[cancer], "arm": arm, "n_terms": nterms[arm],
                     "msk_C": round(im, 4), "msk_lo": round(il, 4), "msk_hi": round(ih, 4),
                     "genie_C": round(em, 4), "genie_lo": round(el, 4), "genie_hi": round(eh, 4),
                     "msk_C_bootmean": round(ibm, 4), "genie_C_bootmean": round(ebm, 4)})
        print(f"  {DISP[cancer]:5} {arm:16} terms={nterms[arm]:3}  MSK={im:.3f} ({il:.3f}-{ih:.3f})  "
              f"GENIE={em:.3f} ({el:.3f}-{eh:.3f})", flush=True)
    for a, bname in CONTRASTS:
        for coh, t, e, idx in (("MSK held-out", t_i, e_i, 0), ("GENIE BPC", t_e, e_e, 1)):
            d, lo, hi, pv, dbm = paired(t, e, np.asarray(risks[a][idx], float), np.asarray(risks[bname][idx], float), rng)
            deltas.append({"cancer": DISP[cancer], "contrast": f"{a} - {bname}", "cohort": coh,
                           "deltaC": round(d, 4), "lo": round(lo, 4), "hi": round(hi, 4),
                           "excludes_0": bool(lo > 0 or hi < 0), "deltaC_bootmean": round(dbm, 4),
                           "boot_p": round(pv, 4)})
    # sensitivity (internal only; separate RNG so the main results above are unchanged):
    # put PW_DRIVER_BURDEN_z back into the pathway pool
    srng = np.random.default_rng(7)
    sens = []
    pool_db = [c for c in train.columns if c.startswith("PW_")]
    for name, cap in (("pathway_EN_capK", K), ("pathway_EN_free", None)):
        sc = StandardScaler().fit(train[pool_db].values)
        m, nu = HH.coxnet_select(sc.transform(train[pool_db].values), y_tr, cap_k=cap,
                                 cap_mode="cv_best" if cap is not None else "densest")
        coef = dict(zip(pool_db, m.coef_[:, 0]))
        r_db = m.predict(sc.transform(test[pool_db].values))
        c_db = boot_ci(t_i, e_i, r_db, srng)[0]
        c_ex = boot_ci(t_i, e_i, np.asarray(risks[name][0], float), srng)[0]
        sens.append({"cancer": DISP[cancer], "arm": name, "terms_with_DB": nu,
                     "DB_selected": bool(abs(coef["PW_DRIVER_BURDEN_z"]) > 1e-8),
                     "msk_C_with_DB": round(c_db, 4), "msk_C_without_DB": round(c_ex, 4),
                     "difference": round(c_db - c_ex, 4)})
        print(f"  [sens] {DISP[cancer]:5} {name:16} with DB terms={nu:3} C={c_db:.3f} | without C={c_ex:.3f} "
              f"| DB selected={sens[-1]['DB_selected']}", flush=True)
    return rows, deltas, sens


def main():
    global COVERAGE, OUT
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--coverage", choices=["observed", "panel"], default="observed")
    COVERAGE = ap.parse_args().coverage
    if COVERAGE == "panel":
        OUT = OUT / "panel_coverage"
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(42)
    R, D, SE = [], [], []
    for c in HH.CANCERS:
        print(f"=== {c}", flush=True)
        r, d, se = run_cancer(c, rng); R += r; D += d; SE += se
    pd.DataFrame(R).to_csv(OUT / "arm_metrics.csv", index=False)
    D = pd.DataFrame(D)
    D["BH_q_all40"] = np.round(bh(D["boot_p"].values), 4)
    D["BH_q_within_contrast"] = np.nan
    for k, g in D.groupby("contrast"):
        D.loc[g.index, "BH_q_within_contrast"] = np.round(bh(g["boot_p"].values), 4)
    D.to_csv(OUT / "paired_deltaC.csv", index=False)
    pd.DataFrame(SE).to_csv(OUT / "driver_burden_sensitivity.csv", index=False)
    print(pd.DataFrame(D).to_string(index=False))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
