#!/usr/bin/env python3
"""Phase A1 — Clinical Decision Pathway Simulation.

Computes three clinical-decision scenarios per cancer, converting DCA
net-benefit values into deployable decision metrics:

  Scenario 1 — Clinical-trial enrichment (top-quartile screening)
  Scenario 2 — Risk-stratified surveillance scheduling
  Scenario 3 — Multidisciplinary referral at MSK-derived median cutoff

All metrics use the held-out MSK internal test set per cancer.
Compact score formula (fixed β + MSK cutoff) is loaded from YAML.
Event probabilities at each cancer-specific time horizon are estimated
via Kaplan-Meier, properly handling right-censoring.

Outputs:
  - analyses/clinical_simulation/scenario1_trial_enrichment.csv
  - analyses/clinical_simulation/scenario2_surveillance.csv
  - analyses/clinical_simulation/scenario3_mdt_referral.csv
  - analyses/clinical_simulation/summary_table.csv  (Table 4 candidate)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from apply_compact_score import apply_score, load_formula

OUT = ROOT / "analyses" / "clinical_simulation"
OUT.mkdir(parents=True, exist_ok=True)

# Cancer-specific clinical-decision time horizons (matches DCA/calibration)
HORIZON_MONTHS = {
    "BRCA": 36,
    "PRAD": 36,
    "CRC": 36,
    "LUAD": 24,
    "PAAD": 24,
}

# Phase II trial target enrollment for Scenario 1
TRIAL_TARGET = 200

# Risk-stratified surveillance cadences (visits/year) for Scenario 2
SURVEILLANCE_CADENCE = {
    "Q4": 4,  # top quartile, every 3 months
    "Q3": 3,  # every 4 months
    "Q2": 3,  # every 4 months
    "Q1": 2,  # bottom quartile, every 6 months
}
UNIFORM_CADENCE = 3  # every 4 months for all (baseline)


def km_event_probability(times, events, t):
    """KM-based cumulative event probability at time t."""
    if len(times) == 0:
        return np.nan
    kmf = KaplanMeierFitter()
    kmf.fit(times, event_observed=events)
    surv = kmf.survival_function_at_times(t).iloc[0]
    return float(1 - surv)


def metrics_at_threshold(df, score_col, threshold, horizon, horizon_eps=1.0):
    """Compute sensitivity / specificity / PPV / NPV / referral rate
    at a given score threshold and time horizon.

    Patients censored before horizon are excluded from the event-definition
    denominator for sens/spec/PPV/NPV (standard handling). Referral rate
    uses all patients.
    """
    n_total = len(df)
    n_above = int((df[score_col] > threshold).sum())
    referral_rate = n_above / n_total if n_total else np.nan

    # Define event status at horizon: only patients with sufficient follow-up
    df = df.copy()
    df["event_at_horizon"] = np.where(
        (df["Event_OS"] == 1) & (df["OS_MONTHS"] <= horizon), 1,
        np.where(df["OS_MONTHS"] >= horizon - horizon_eps, 0, np.nan)
    )
    cls = df.dropna(subset=["event_at_horizon"]).copy()
    cls["event_at_horizon"] = cls["event_at_horizon"].astype(int)
    cls["above"] = (cls[score_col] > threshold).astype(int)

    tp = int(((cls["above"] == 1) & (cls["event_at_horizon"] == 1)).sum())
    fp = int(((cls["above"] == 1) & (cls["event_at_horizon"] == 0)).sum())
    fn = int(((cls["above"] == 0) & (cls["event_at_horizon"] == 1)).sum())
    tn = int(((cls["above"] == 0) & (cls["event_at_horizon"] == 0)).sum())

    sens = tp / (tp + fn) if (tp + fn) else np.nan
    spec = tn / (tn + fp) if (tn + fp) else np.nan
    ppv = tp / (tp + fp) if (tp + fp) else np.nan
    npv = tn / (tn + fn) if (tn + fn) else np.nan

    return {
        "n_total": n_total,
        "n_above_threshold": n_above,
        "referral_rate": referral_rate,
        "n_classifiable": len(cls),
        "n_events_at_horizon": int(cls["event_at_horizon"].sum()),
        "sens": sens, "spec": spec, "ppv": ppv, "npv": npv,
        "tp": tp, "fp": fp, "fn": fn, "tn": tn,
    }


def scenario1_trial_enrichment(df, score_col, horizon, target=TRIAL_TARGET):
    """Top-quartile screening for trial enrichment."""
    q3 = df[score_col].quantile(0.75)
    overall_event_prob = km_event_probability(df["OS_MONTHS"], df["Event_OS"], horizon)

    high = df[df[score_col] > q3]
    high_event_prob = km_event_probability(high["OS_MONTHS"], high["Event_OS"], horizon)

    rest = df[df[score_col] <= q3]
    rest_event_prob = km_event_probability(rest["OS_MONTHS"], rest["Event_OS"], horizon)

    enrichment = high_event_prob / overall_event_prob if overall_event_prob else np.nan

    # NNS to identify one high-risk patient in Q4 (= 4, since 25%)
    # NNS to identify one expected event in Q4 = 1 / high_event_prob
    nns_one_event_random = 1 / overall_event_prob if overall_event_prob else np.nan
    nns_one_event_q4 = (1 / 0.25) * (1 / high_event_prob) if high_event_prob else np.nan
    # i.e., to screen 1 expected event in Q4 = screen 4 patients to get 1 Q4 patient,
    # then need 1/p events per Q4 patient
    # Equivalent: 4 / high_event_prob

    # Phase II trial: target N high-risk enrolled
    screening_random = target / overall_event_prob if overall_event_prob else np.nan
    screening_q4 = target / (high_event_prob * 0.25) if high_event_prob else np.nan
    # i.e., need target events; if Q4 gives event prob X, and 25% are Q4,
    # screen 1 person → 0.25 * X probability of getting a Q4 event

    # Simpler: to enroll target Q4 patients (event-conditioned), need target / (0.25 * high_event_prob) screened
    # But more useful: to enroll target high-risk patients in Q4 (regardless of confirmed event),
    # need target / 0.25 = 4 * target screened
    enrollment_q4 = target / 0.25  # = 4 × target

    return {
        "n_test": len(df),
        "horizon_months": horizon,
        "overall_event_prob": overall_event_prob,
        "q3_threshold": float(q3),
        "n_q4": len(high),
        "q4_event_prob": high_event_prob,
        "rest_event_prob": rest_event_prob,
        "enrichment_ratio": enrichment,
        "nns_one_event_random": nns_one_event_random,
        "nns_one_event_q4": nns_one_event_q4,
        "n_screened_to_enroll_q4_target": enrollment_q4,
        "expected_events_among_q4_target": target * high_event_prob,
        "expected_events_random_same_size": target * overall_event_prob,
    }


def scenario2_surveillance(df, score_col, horizon):
    """Risk-stratified surveillance scheduling."""
    q1, q2, q3 = df[score_col].quantile([0.25, 0.50, 0.75])
    df = df.copy()
    df["quartile"] = pd.cut(
        df[score_col],
        bins=[-np.inf, q1, q2, q3, np.inf],
        labels=["Q1", "Q2", "Q3", "Q4"],
    )

    out = {}
    for q in ["Q1", "Q2", "Q3", "Q4"]:
        sub = df[df["quartile"] == q]
        event_prob = km_event_probability(sub["OS_MONTHS"], sub["Event_OS"], horizon)
        out[f"{q}_n"] = len(sub)
        out[f"{q}_event_prob"] = event_prob
        out[f"{q}_cadence_per_yr"] = SURVEILLANCE_CADENCE[q]

    # Imaging volume (visits per year per patient)
    n = len(df)
    risk_strat_volume = sum(
        SURVEILLANCE_CADENCE[q] * out[f"{q}_n"] / n for q in ["Q1", "Q2", "Q3", "Q4"]
    )
    uniform_volume = UNIFORM_CADENCE

    # % of horizon-events occurring in each quartile (sensitivity-style)
    df["event_at_horizon"] = np.where(
        (df["Event_OS"] == 1) & (df["OS_MONTHS"] <= horizon), 1,
        np.where(df["OS_MONTHS"] >= horizon, 0, np.nan)
    )
    cls = df.dropna(subset=["event_at_horizon"])
    total_events = int(cls["event_at_horizon"].sum())
    for q in ["Q1", "Q2", "Q3", "Q4"]:
        sub = cls[cls["quartile"] == q]
        events_in_q = int(sub["event_at_horizon"].sum())
        out[f"{q}_pct_of_events"] = events_in_q / total_events if total_events else np.nan

    out["risk_strat_avg_visits_per_yr"] = risk_strat_volume
    out["uniform_avg_visits_per_yr"] = uniform_volume
    out["q4_imaging_share_pct"] = (
        SURVEILLANCE_CADENCE["Q4"] * out["Q4_n"] / sum(
            SURVEILLANCE_CADENCE[q] * out[f"{q}_n"] for q in ["Q1", "Q2", "Q3", "Q4"]
        )
    )
    out["q4_event_share_pct"] = out["Q4_pct_of_events"]
    out["concentration_ratio"] = (
        out["q4_event_share_pct"] / out["q4_imaging_share_pct"]
        if out["q4_imaging_share_pct"] else np.nan
    )
    out["q1_volume_reduction_vs_uniform"] = (
        (UNIFORM_CADENCE - SURVEILLANCE_CADENCE["Q1"]) / UNIFORM_CADENCE
    )
    return out


def scenario3_mdt_referral(df, score_col, threshold, horizon, label=""):
    """MDT referral triggered when compact score > threshold."""
    m = metrics_at_threshold(df, score_col, threshold, horizon)
    m["threshold_label"] = label
    m["threshold_value"] = float(threshold)
    return m


GENIE_MAPPING = {
    "BRCA": "GENIE_BPC_BrCa",
    "LUAD": "GENIE_BPC_NSCLC",
    "PAAD": "GENIE_BPC_PANC",
    "PRAD": "GENIE_BPC_Prostate",
    "CRC": "GENIE_BPC_CRC",
}


def load_external(cancer):
    name = GENIE_MAPPING[cancer]
    p = ROOT / "external_validation" / "outputs" / name / f"{name}_pathway_features.csv"
    if not p.exists():
        return None
    df = pd.read_csv(p)
    if "OS_MONTHS" not in df.columns or "Event_OS" not in df.columns:
        return None
    return df


def run_all():
    formula = load_formula()
    cancers = ["BRCA", "LUAD", "PAAD", "PRAD", "CRC"]

    s1_int, s2_int, s3_int = [], [], []
    s1_ext, s3_ext = [], []
    for cancer in cancers:
        # Internal MSK test set
        test_path = ROOT / "model_input" / "filter1_age_fixed" / f"{cancer}_test.csv"
        df = pd.read_csv(test_path)
        _, score = apply_score(df, cancer, formula)
        df["compact_score"] = score
        horizon = HORIZON_MONTHS[cancer]
        cutoff = formula["cancers"][cancer]["msk_partial_hazard_cutoff"]
        q3 = df["compact_score"].quantile(0.75)
        test_median = df["compact_score"].quantile(0.50)

        r1 = scenario1_trial_enrichment(df, "compact_score", horizon)
        r1["cancer"] = cancer
        r1["cohort"] = "MSK_test"
        s1_int.append(r1)

        r2 = scenario2_surveillance(df, "compact_score", horizon)
        r2["cancer"] = cancer
        r2["cohort"] = "MSK_test"
        s2_int.append(r2)

        for thr, label in [(cutoff, "MSK_median_cutoff"),
                           (test_median, "Test_median"),
                           (q3, "Top_quartile_Q4")]:
            r = scenario3_mdt_referral(df, "compact_score", thr, horizon, label)
            r["cancer"] = cancer
            r["cohort"] = "MSK_test"
            r["horizon_months"] = horizon
            s3_int.append(r)

        # External GENIE BPC
        ext_df = load_external(cancer)
        if ext_df is None:
            continue
        if "risk_score" in ext_df.columns:
            ext_df["compact_score"] = ext_df["risk_score"]
        else:
            try:
                _, score = apply_score(ext_df, cancer, formula)
                ext_df["compact_score"] = score
            except KeyError:
                continue
        ext_df = ext_df.dropna(subset=["OS_MONTHS", "Event_OS", "compact_score"]).copy()
        if len(ext_df) == 0:
            continue
        ext_q3 = ext_df["compact_score"].quantile(0.75)
        ext_median = ext_df["compact_score"].quantile(0.50)

        r1e = scenario1_trial_enrichment(ext_df, "compact_score", horizon)
        r1e["cancer"] = cancer
        r1e["cohort"] = "GENIE_BPC"
        s1_ext.append(r1e)

        for thr, label in [(cutoff, "MSK_median_cutoff"),
                           (ext_median, "Ext_median"),
                           (ext_q3, "Top_quartile_Q4")]:
            r = scenario3_mdt_referral(ext_df, "compact_score", thr, horizon, label)
            r["cancer"] = cancer
            r["cohort"] = "GENIE_BPC"
            r["horizon_months"] = horizon
            s3_ext.append(r)

    s1 = s1_int + s1_ext
    s2 = s2_int
    s3 = s3_int + s3_ext

    df1 = pd.DataFrame(s1)
    df2 = pd.DataFrame(s2)
    df3 = pd.DataFrame(s3)

    df1.to_csv(OUT / "scenario1_trial_enrichment.csv", index=False)
    df2.to_csv(OUT / "scenario2_surveillance.csv", index=False)
    df3.to_csv(OUT / "scenario3_mdt_referral.csv", index=False)

    # ---- Summary table candidates ----

    # Table 4A: Scenario 1 (trial enrichment) internal vs external
    s1_pivot = df1.pivot_table(
        index="cancer",
        columns="cohort",
        values=["overall_event_prob", "q4_event_prob", "enrichment_ratio",
                "nns_one_event_random", "nns_one_event_q4",
                "expected_events_among_q4_target", "expected_events_random_same_size"],
    ).round(3)
    s1_pivot.to_csv(OUT / "table4A_scenario1_pivot.csv")

    # Table 4B: Scenario 3 (MDT) — all thresholds, internal + external
    s3_pivot = df3.set_index(["cancer", "cohort", "threshold_label"])[
        ["referral_rate", "sens", "spec", "ppv", "npv", "n_events_at_horizon"]
    ].round(3)
    s3_pivot.to_csv(OUT / "table4B_scenario3_pivot.csv")

    # Table 4 main: per-cancer headline numbers (internal Q4 + external Q4)
    int_q4 = df3[(df3["cohort"] == "MSK_test") & (df3["threshold_label"] == "Top_quartile_Q4")].set_index("cancer")
    int_msk = df3[(df3["cohort"] == "MSK_test") & (df3["threshold_label"] == "MSK_median_cutoff")].set_index("cancer")
    ext_q4 = df3[(df3["cohort"] == "GENIE_BPC") & (df3["threshold_label"] == "Top_quartile_Q4")].set_index("cancer")
    s1i = df1[df1["cohort"] == "MSK_test"].set_index("cancer")
    s1e = df1[df1["cohort"] == "GENIE_BPC"].set_index("cancer")

    summary_rows = []
    for c in cancers:
        row = {
            "cancer": c,
            "horizon_months": HORIZON_MONTHS[c],
            # Internal
            "int_n": int(s1i.at[c, "n_test"]),
            "int_overall_event_pct": round(s1i.at[c, "overall_event_prob"] * 100, 1),
            "int_q4_event_pct": round(s1i.at[c, "q4_event_prob"] * 100, 1),
            "int_enrichment_ratio": round(s1i.at[c, "enrichment_ratio"], 2),
            "int_nns_q4": round(s1i.at[c, "nns_one_event_q4"], 1),
            "int_mdt_q4_referral_pct": round(int_q4.at[c, "referral_rate"] * 100, 1),
            "int_mdt_q4_sens_pct": round(int_q4.at[c, "sens"] * 100, 1),
            "int_mdt_q4_spec_pct": round(int_q4.at[c, "spec"] * 100, 1),
            "int_mdt_q4_ppv_pct": round(int_q4.at[c, "ppv"] * 100, 1),
            "int_mdt_msk_referral_pct": round(int_msk.at[c, "referral_rate"] * 100, 1),
            "int_mdt_msk_sens_pct": round(int_msk.at[c, "sens"] * 100, 1),
            "int_mdt_msk_spec_pct": round(int_msk.at[c, "spec"] * 100, 1),
        }
        if c in s1e.index:
            row.update({
                "ext_n": int(s1e.at[c, "n_test"]),
                "ext_overall_event_pct": round(s1e.at[c, "overall_event_prob"] * 100, 1),
                "ext_q4_event_pct": round(s1e.at[c, "q4_event_prob"] * 100, 1),
                "ext_enrichment_ratio": round(s1e.at[c, "enrichment_ratio"], 2),
                "ext_nns_q4": round(s1e.at[c, "nns_one_event_q4"], 1),
            })
        if c in ext_q4.index:
            row.update({
                "ext_mdt_q4_referral_pct": round(ext_q4.at[c, "referral_rate"] * 100, 1),
                "ext_mdt_q4_sens_pct": round(ext_q4.at[c, "sens"] * 100, 1),
                "ext_mdt_q4_spec_pct": round(ext_q4.at[c, "spec"] * 100, 1),
                "ext_mdt_q4_ppv_pct": round(ext_q4.at[c, "ppv"] * 100, 1),
            })
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT / "summary_table.csv", index=False)

    print("\n=== Scenario 1 — Trial enrichment (Internal + External) ===")
    print(df1.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print("\n=== Scenario 2 — Risk-stratified surveillance (Internal) ===")
    print(df2.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print("\n=== Scenario 3 — MDT referral (multi-threshold, Internal + External) ===")
    print(df3.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print("\n=== Summary (Table 4 candidate) ===")
    print(summary.to_string(index=False))
    return summary


if __name__ == "__main__":
    run_all()
