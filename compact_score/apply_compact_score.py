#!/usr/bin/env python3
"""Apply the published compact pathway-level risk score to a new cohort.

Two entry points:
  1. ``apply_score(df, cancer)`` — given a DataFrame already containing the
     PW_* features used by ``cancer``, return ``(linear_predictor, risk_score)``.
     Useful when the caller built features with their own pipeline.
  2. CLI: ``python apply_compact_score.py --cancer BRCA --features X.csv``
     reads a feature CSV, writes a copy with added ``linear_predictor``,
     ``risk_score`` and ``risk_group_msk_cutoff`` columns, and prints
     C-index + Kaplan-Meier log-rank p (if OS_MONTHS / Event_OS present).

The formula and β are loaded from
``external_validation/zscore_params/compact_score_formula.yaml``.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
FORMULA = ROOT / "external_validation" / "zscore_params" / "compact_score_formula.yaml"


def load_formula(path: Path = FORMULA) -> dict:
    with path.open() as f:
        return yaml.safe_load(f)


def apply_score(df: pd.DataFrame, cancer: str, formula: dict | None = None
                ) -> tuple[pd.Series, pd.Series]:
    formula = formula or load_formula()
    if cancer not in formula["cancers"]:
        raise ValueError(f"unknown cancer {cancer!r}; valid: {list(formula['cancers'])}")
    spec = formula["cancers"][cancer]
    missing = [f["name"] for f in spec["features"] if f["name"] not in df.columns]
    if missing:
        raise KeyError(f"{cancer}: missing features in input df: {missing}")
    linear = sum(float(f["beta"]) * df[f["name"]] for f in spec["features"])
    risk = np.exp(linear)
    return linear, risk


def _cindex(t, score, e):
    from lifelines.utils import concordance_index
    return float(concordance_index(t, -score, e))


def _logrank(df, group):
    from lifelines.statistics import logrank_test
    g = df.assign(_g=group.values).dropna(subset=["OS_MONTHS", "Event_OS", "_g"])
    if g["_g"].nunique() != 2:
        return float("nan")
    vals = list(g["_g"].unique())
    a = g[g["_g"] == vals[0]]; b = g[g["_g"] == vals[1]]
    return float(logrank_test(a["OS_MONTHS"], b["OS_MONTHS"], a["Event_OS"], b["Event_OS"]).p_value)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cancer", required=True, choices=["BRCA", "LUAD", "PAAD", "PRAD", "CRC"])
    ap.add_argument("--features", required=True, type=Path,
                    help="CSV with PW_* feature columns (and optionally OS_MONTHS, Event_OS)")
    ap.add_argument("--out", type=Path, default=None,
                    help="output CSV (default: <features>.with_score.csv)")
    args = ap.parse_args()

    df = pd.read_csv(args.features)
    formula = load_formula()
    spec = formula["cancers"][args.cancer]
    linear, risk = apply_score(df, args.cancer, formula)
    df["linear_predictor"] = linear
    df["risk_score"] = risk
    df["risk_group_msk_cutoff"] = (risk > float(spec["msk_partial_hazard_cutoff"])).astype(int)

    out = args.out or args.features.with_suffix(".with_score.csv")
    df.to_csv(out, index=False)
    print(f"wrote {out}")
    print(f"features used ({args.cancer}):")
    for f in spec["features"]:
        print(f"  {f['name']:32s}  β = {float(f['beta']):+.6f}")
    print(f"msk_partial_hazard_cutoff = {spec['msk_partial_hazard_cutoff']}")

    if {"OS_MONTHS", "Event_OS"}.issubset(df.columns):
        valid = df[["OS_MONTHS", "Event_OS", "risk_score"]].dropna()
        c = _cindex(valid["OS_MONTHS"], valid["risk_score"], valid["Event_OS"])
        p = _logrank(df, df["risk_group_msk_cutoff"].map({0: "Low", 1: "High"}))
        n_high = int(df["risk_group_msk_cutoff"].sum())
        print(f"\nC-index             = {c:.4f}")
        print(f"KM log-rank p (MSK cutoff): {p:.3e}  (n_high={n_high}/{len(df)})")


if __name__ == "__main__":
    main()
