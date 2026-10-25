#!/usr/bin/env python3
"""
burden_only_benchmark.py
========================
Burden-only comparator (Supplementary Table S10).

Does the 4-5-feature compact score carry information beyond the overall pathway-alteration
burden?  Arm "burden_only" = Cox on PW_HIT_COUNT (All-pathway alteration count) alone, fitted
on the MSK training split; it is compared with the fixed compact score.

Both arms are fixed on MSK training data and applied unchanged to the MSK held-out test split
and to the GENIE BPC external cohorts.  Point estimates are observed full-cohort C-indices and
paired differences; 95% CIs are percentiles of 1,000 paired patient-level bootstrap resamples.

The bootstrap random-number stream (seed 42) is advanced exactly as in the original run, which
also evaluated comparators that are not reported, so the archived Supplementary Table S10
values are reproduced exactly.

Outputs (figures_tables/burden_only_benchmark/):
  arm_metrics.csv      C-index (bootstrap 95% CI) per cancer x arm
  paired_deltaC.csv    paired bootstrap ΔC, compact score minus burden-only model
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import gene_vs_pathway_external as EX  # noqa: E402
import gene_vs_pathway_headtohead as HH  # noqa: E402
from compact_score.apply_compact_score import apply_score, load_formula  # noqa: E402

OUT = ROOT / "figures_tables" / "burden_only_benchmark"
DISP = {"BRCA": "IDC", "LUAD": "LUAD", "PAAD": "PAAD", "PRAD": "PRAD", "CRC": "CRC"}
B = 1000
N_UNREPORTED_ARMS = 4        # comparators evaluated in the original run but not reported
N_UNREPORTED_CONTRASTS = 3


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
    obs = float(concordance_index(t, -ra, e) - concordance_index(t, -rb, e))
    return obs, float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5)), float(d.mean())


def advance(rng, n):
    """Consume one bootstrap's worth of draws (B resamples of size n)."""
    for _ in range(B):
        rng.integers(0, n, n)


def run_cancer(cancer, rng):
    _, out_dir, _ = EX.COHORTS[cancer]
    key = HH.CANCER_ALIAS.get(cancer, cancer)
    formula = load_formula()
    K = len(formula["cancers"][formula.get("cancer_key_aliases", {}).get(key, key)]["features"])
    train, test, _, _ = HH.load_split(cancer)
    t_i = test["OS_MONTHS"].to_numpy(float); e_i = test["Event_OS"].astype(bool).to_numpy()
    ext = pd.read_csv(EX.OUT_BASE / out_dir / f"{out_dir}_pathway_features.csv")
    ext = ext[(ext["OS_MONTHS"] > 0) & ext["Event_OS"].notna()].copy()
    t_e = ext["OS_MONTHS"].to_numpy(float); e_e = ext["Event_OS"].astype(bool).to_numpy()

    # compact (fixed published formula); sanity: external risk_score == formula applied to ext
    ci, _ = apply_score(test, key, formula)     # log score (monotone in exp risk; same C-index)
    ce, _ = apply_score(ext, key, formula)
    assert np.allclose(ce.to_numpy(float), ext["risk_score"].to_numpy(float), atol=1e-6), "external z/feature mismatch"
    risks = {"compact_fixed": (ci.to_numpy(float), ext["risk_score"].to_numpy(float))}
    nterms = {"compact_fixed": K, "burden_only": 1}
    # burden only: univariable Cox on All-pathway alteration count (MSK train)
    cph = CoxPHFitter().fit(train[["PW_HIT_COUNT", "OS_MONTHS", "Event_OS"]], "OS_MONTHS", "Event_OS")
    b = float(cph.params_["PW_HIT_COUNT"])
    risks["burden_only"] = (b * test["PW_HIT_COUNT"].to_numpy(float), b * ext["PW_HIT_COUNT"].to_numpy(float))

    rows = []
    for arm, (ri, re) in risks.items():
        im, il, ih, ibm = boot_ci(t_i, e_i, np.asarray(ri, float), rng)
        em, el, eh, ebm = boot_ci(t_e, e_e, np.asarray(re, float), rng)
        rows.append({"cancer": DISP[cancer], "arm": arm, "n_terms": nterms[arm],
                     "msk_C": round(im, 4), "msk_lo": round(il, 4), "msk_hi": round(ih, 4),
                     "genie_C": round(em, 4), "genie_lo": round(el, 4), "genie_hi": round(eh, 4),
                     "msk_C_bootmean": round(ibm, 4), "genie_C_bootmean": round(ebm, 4)})
        print(f"  {DISP[cancer]:5} {arm:14} terms={nterms[arm]:2}  MSK={im:.3f} ({il:.3f}-{ih:.3f})  "
              f"GENIE={em:.3f} ({el:.3f}-{eh:.3f})", flush=True)
    for _ in range(N_UNREPORTED_ARMS):
        advance(rng, len(t_i)); advance(rng, len(t_e))

    deltas = []
    for coh, t, e, idx in (("MSK held-out", t_i, e_i, 0), ("GENIE BPC", t_e, e_e, 1)):
        d, lo, hi, dbm = paired(t, e, np.asarray(risks["compact_fixed"][idx], float),
                                np.asarray(risks["burden_only"][idx], float), rng)
        deltas.append({"cancer": DISP[cancer], "contrast": "compact_fixed - burden_only", "cohort": coh,
                       "deltaC": round(d, 4), "lo": round(lo, 4), "hi": round(hi, 4),
                       "excludes_0": bool(lo > 0 or hi < 0), "deltaC_bootmean": round(dbm, 4)})
    for _ in range(N_UNREPORTED_CONTRASTS):
        advance(rng, len(t_i)); advance(rng, len(t_e))
    return rows, deltas


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(42)
    R, D = [], []
    for c in HH.CANCERS:
        print(f"=== {c}", flush=True)
        r, d = run_cancer(c, rng); R += r; D += d
    pd.DataFrame(R).to_csv(OUT / "arm_metrics.csv", index=False)
    pd.DataFrame(D).to_csv(OUT / "paired_deltaC.csv", index=False)
    print(pd.DataFrame(D).to_string(index=False))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
