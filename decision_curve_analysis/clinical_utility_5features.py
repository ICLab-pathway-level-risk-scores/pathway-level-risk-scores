#!/usr/bin/env python3
"""Clinical utility analyses for the current feature-score model.

Outputs:
- adjusted Cox forest plot for the pathway feature-score
- clinical-only vs score-only vs combined C-index table
- decision curve analysis for two candidate time points per cancer; only the
  better-performing time point is plotted per cancer
- cancer-specific calibration plots for the combined model

Decision curves and calibration use a binary endpoint at the selected time point.
Events by t are cases; patients surviving ≥t are controls;
patients censored before t are excluded from those binary summaries.
"""
from __future__ import annotations

import os
from pathlib import Path
import warnings

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib-cache")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from figure_style import apply_arial_style
apply_arial_style()
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from scipy import stats
from scipy.ndimage import uniform_filter1d

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "model_input" / "filter1_age_fixed"
OUT = ROOT / "figures_tables" / "clinical_utility_5features"
OUT.mkdir(parents=True, exist_ok=True)

CANCERS = ["BRCA", "LUAD", "PAAD", "PRAD", "CRC"]
CANCER_LABEL = {"BRCA": "IDC", "LUAD": "LUAD", "PAAD": "PAAD",
                "PRAD": "PRAD", "CRC": "CRC"}

DURATION = "OS_MONTHS"
EVENT = "Event_OS"

# Candidate DCA/calibration time points (months). The better one is chosen per
# cancer by average net benefit gain of Clinical + feature-score vs Clinical only.
CANDIDATE_TIMEPOINTS = {
    "BRCA": [36, 60],
    "LUAD": [24, 36],
    "PAAD": [12, 24],
    "PRAD": [36, 60],
    "CRC":  [24, 36],
}

FEATURES = {
    "BRCA": ["PW_HIT_COUNT", "PW_TP53_mut_rate_z", "PW_Cell_Cycle_sv_hit", "PW_Chromatin_zsum", "PW_RTK_RAS_amp_hit"],
    "LUAD": ["PW_HIT_COUNT", "PW_TP53_any_rate_z", "PW_NOTCH_any_rate_z", "PW_NRF2_zsum", "PW_RTK_RAS_mut_rate_z"],
    "PAAD": ["PW_HIT_COUNT", "PW_RTK_RAS_amp_hit", "PW_TP53_mut_rate_z", "PW_TGF_Beta_any_rate_z", "PW_MYC_any_rate_z"],
    "PRAD": ["PW_HIT_COUNT", "PW_Chromatin_sv_hit", "PW_Cell_Cycle_zsum", "PW_WNT_zsum"],
    "CRC":  ["PW_HIT_COUNT", "PW_DDR_zsum", "PW_RTK_RAS_any_rate_z", "PW_TGF_Beta_amp_hit", "PW_Chromatin_any_rate_z"],
}

CLINICAL_COMMON = ["age_z", "ecog_z", "SAMPLE_TYPE_Metastasis"]
CLINICAL_EXTRA  = {"BRCA": ["HR_positive"], "CRC": ["MSI_high"]}

CLIN_DIR = Path("/Users/kao/Downloads/all timeline csv")
ID_DIR   = ROOT / "model_input" / "sample_id_mapping"


def load_extra_clinical() -> dict:
    """Return per-sample dict with HR_positive (BRCA) and MSI_high (CRC)."""
    sample  = pd.read_csv(CLIN_DIR / "data_clinical_sample.csv",
                          usecols=["SAMPLE_ID", "PATIENT_ID", "MSI_TYPE"])
    patient = pd.read_csv(CLIN_DIR / "data_clinical_patient.csv",
                          usecols=["PATIENT_ID", "HR"])
    merged = sample.merge(patient, on="PATIENT_ID", how="left")
    merged["HR_positive"] = (merged["HR"] == "Yes").astype(float)
    merged["MSI_high"]    = (merged["MSI_TYPE"] == "Instable").astype(float)
    return merged.set_index("SAMPLE_ID")[["HR_positive", "MSI_high"]].to_dict(orient="index")


def attach_extra(df: pd.DataFrame, cancer: str, ids: np.ndarray, clin_map: dict) -> pd.DataFrame:
    df = df.copy()
    df.insert(0, "SAMPLE_ID", ids)
    for col in CLINICAL_EXTRA.get(cancer, []):
        df[col] = df["SAMPLE_ID"].map(lambda s, c=col: clin_map.get(s, {}).get(c, float("nan")))
    return df


def clinical_cols(cancer: str, df: pd.DataFrame) -> list[str]:
    cols = CLINICAL_COMMON + CLINICAL_EXTRA.get(cancer, [])
    return [c for c in cols if c in df.columns]


def fit_score(train: pd.DataFrame, test: pd.DataFrame, features: list[str]) -> tuple[np.ndarray, np.ndarray, CoxPHFitter]:
    cph = CoxPHFitter(penalizer=0.1)
    cph.fit(train[features + [DURATION, EVENT]], duration_col=DURATION, event_col=EVENT)
    train_score = cph.predict_log_partial_hazard(train[features]).to_numpy()
    test_score = cph.predict_log_partial_hazard(test[features]).to_numpy()
    return train_score, test_score, cph


def fit_model(train: pd.DataFrame, test: pd.DataFrame, cols: list[str], penalizer: float = 0.01) -> tuple[CoxPHFitter, float]:
    cph = CoxPHFitter(penalizer=penalizer)
    cph.fit(train[cols + [DURATION, EVENT]], duration_col=DURATION, event_col=EVENT)
    cindex = cph.score(test[cols + [DURATION, EVENT]], scoring_method="concordance_index")
    return cph, float(cindex)


def lrt_p(full: CoxPHFitter, reduced: CoxPHFitter, df: int = 1) -> float:
    stat = 2 * (full.log_likelihood_ - reduced.log_likelihood_)
    return float(stats.chi2.sf(max(stat, 0), df=df))


def binary_at_t(df: pd.DataFrame, t: int) -> pd.Series:
    """Return 1/0 for known event status at t months; NaN if censored before t."""
    event_by_t    = (df[EVENT].astype(int) == 1) & (df[DURATION] <= t)
    known_control = df[DURATION] >= t
    y = pd.Series(np.nan, index=df.index, dtype=float)
    y[event_by_t]    = 1.0
    y[known_control] = 0.0
    return y


def predict_risk_at(cph: CoxPHFitter, df: pd.DataFrame, cols: list[str], t: int) -> np.ndarray:
    surv = cph.predict_survival_function(df[cols], times=[t]).T.iloc[:, 0].to_numpy()
    return np.clip(1 - surv, 0, 1)


def net_benefit(y_true: np.ndarray, risk: np.ndarray, thresholds: np.ndarray) -> list[float]:
    out = []
    n = len(y_true)
    for pt in thresholds:
        pred = risk >= pt
        tp = np.sum(pred & (y_true == 1))
        fp = np.sum(pred & (y_true == 0))
        out.append(tp / n - fp / n * (pt / (1 - pt)))
    return out


def calibration_bins(y_true: np.ndarray, risk: np.ndarray, n_bins: int = 5) -> pd.DataFrame:
    df = pd.DataFrame({"y": y_true, "risk": risk})
    df["bin"] = pd.qcut(df["risk"].rank(method="first"), q=min(n_bins, len(df)), labels=False)
    rows = []
    for _, g in df.groupby("bin"):
        rows.append({
            "mean_predicted": g["risk"].mean(),
            "observed": g["y"].mean(),
            "n": len(g),
        })
    return pd.DataFrame(rows)


def plot_forest(summary: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(9.8, 4.8))
    plot_df = summary.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(plot_df))
    colors = np.where(plot_df["score_p"] < 0.05, "#c44e52", "#777777")
    for i, row in plot_df.iterrows():
        ax.plot([row["score_hr_low"], row["score_hr_high"]], [i, i], color=colors[i], lw=2)
        ax.scatter(row["score_hr"], i, color=colors[i], s=75, zorder=3)
    ax.axvline(1, color="black", ls="--", lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["cancer"])
    ax.set_xlabel("Adjusted HR for feature-score")
    ax.set_title("Clinical + Feature-Score Cox Model")
    x_right = max(1.1, plot_df["score_hr_high"].max() * 1.55)
    ax.set_xlim(0, x_right)
    for i, row in plot_df.iterrows():
        sig = "***" if row["score_p"] < 0.001 else "**" if row["score_p"] < 0.01 else "*" if row["score_p"] < 0.05 else "ns"
        txt = f"HR={row['score_hr']:.2f} [{row['score_hr_low']:.2f}-{row['score_hr_high']:.2f}] {sig}; ΔC={row['delta_cindex']:+.3f}; LRT p={row['lrt_p']:.2e}"
        ax.text(plot_df["score_hr_high"].max() * 1.05, i, txt, va="center", fontsize=8)
    ax.grid(axis="x", alpha=0.25)
    plt.tight_layout()
    fig.savefig(OUT / "forest_adjusted_feature_score_hr.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT / "forest_adjusted_feature_score_hr.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_dca(dca_rows: pd.DataFrame) -> None:
    for cancer, sub in dca_rows.groupby("cancer"):
        t = int(sub["timepoint_m"].iloc[0])
        fig, ax = plt.subplots(figsize=(6.2, 4.4))
        for model, color, lw, smooth in [
            ("Clinical only",           "#4c72b0", 2.0, True),
            ("Feature-score only",      "#55a868", 1.8, True),
            ("Clinical + feature-score","#c44e52", 2.0, True),
            ("Treat all",               "#777777", 1.4, False),
            ("Treat none",              "#222222", 1.2, False),
        ]:
            m = sub[sub["model"] == model].sort_values("threshold")
            if not len(m):
                continue
            y_vals = m["net_benefit"].to_numpy()
            if smooth and len(y_vals) >= 9:
                y_vals = uniform_filter1d(y_vals, size=9)
            ax.plot(m["threshold"], y_vals, label=model, lw=lw, color=color)
        ax.axhline(0, color="black", lw=0.8)
        ax.set_title(f"{CANCER_LABEL.get(cancer, cancer)}: {t}-month Decision Curve")
        ax.set_xlabel("Threshold probability")
        ax.set_ylabel("Net benefit")
        ax.set_xlim(0.05, 0.60)
        y_min = min(-0.05, sub["net_benefit"].min() - 0.01)
        y_max = max(0.10, sub["net_benefit"].max() + 0.02)
        ax.set_ylim(y_min, y_max)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
        plt.tight_layout()
        fig.savefig(OUT / f"{cancer}_decision_curve_{t}m.png", dpi=200, bbox_inches="tight")
        fig.savefig(OUT / f"{cancer}_decision_curve_{t}m.pdf", bbox_inches="tight")
        plt.close(fig)


def plot_calibration(cal_rows: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, len(CANCERS), figsize=(18, 3.8), squeeze=False)
    axes = axes[0]
    for ax, cancer in zip(axes, CANCERS):
        sub = cal_rows[cal_rows["cancer"] == cancer]
        t = int(sub["timepoint_m"].iloc[0]) if len(sub) else CANDIDATE_TIMEPOINTS[cancer][0]
        ax.plot([0, 1], [0, 1], color="black", ls="--", lw=1)
        if len(sub):
            ax.scatter(sub["mean_predicted"], sub["observed"], s=sub["n"] * 2.2,
                       color="#c44e52", alpha=0.85)
            ax.plot(sub["mean_predicted"], sub["observed"], color="#c44e52", lw=1.5)
        ax.set_title(f"{CANCER_LABEL.get(cancer, cancer)} ({t}m)")
        ax.set_xlim(0, 0.9)
        ax.set_ylim(0, 0.9)
        ax.grid(alpha=0.25)
        ax.set_xlabel("Predicted risk")
    axes[0].set_ylabel("Observed event rate")
    fig.suptitle("Calibration: Clinical + Feature-Score Model (cancer-specific time points)",
                 y=1.03, fontweight="bold")
    plt.tight_layout()
    fig.savefig(OUT / "calibration_combined.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT / "calibration_combined.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    summary_rows = []
    selected_dca_rows = []
    selected_cal_rows = []
    candidate_rows = []

    thresholds = np.arange(0.05, 0.605, 0.005)
    clin_map = load_extra_clinical()

    for cancer in CANCERS:
        print(f"\n{cancer}")
        train_raw = pd.read_csv(DATA / cancer / f"{cancer}_train.csv")
        test_raw  = pd.read_csv(DATA / cancer / f"{cancer}_test.csv")
        ids_tr = pd.read_csv(ID_DIR / f"{cancer}_train_ids.csv")["SAMPLE_ID"].values[:len(train_raw)]
        ids_te = pd.read_csv(ID_DIR / f"{cancer}_test_ids.csv")["SAMPLE_ID"].values[:len(test_raw)]

        train = attach_extra(train_raw, cancer, ids_tr, clin_map)
        test  = attach_extra(test_raw,  cancer, ids_te, clin_map)

        features = FEATURES[cancer]
        clin = clinical_cols(cancer, train)

        tr_score, te_score, score_cph = fit_score(train, test, features)
        train["feature_score"] = tr_score
        test["feature_score"]  = te_score

        cph_clin,     c_clin     = fit_model(train, test, clin)
        cph_score,    c_score    = fit_model(train, test, ["feature_score"])
        cph_combined, c_combined = fit_model(train, test, clin + ["feature_score"])
        lrt = lrt_p(cph_combined, cph_clin)

        s = cph_combined.summary.loc["feature_score"]
        row = {
            "cancer": cancer,
            "selected_timepoint_m": None,
            "n_train": len(train),
            "n_test": len(test),
            "test_events": int(test[EVENT].sum()),
            "clinical_cols": ";".join(clin),
            "features": ";".join(features),
            "cindex_clinical": c_clin,
            "cindex_score": c_score,
            "cindex_combined": c_combined,
            "delta_cindex": c_combined - c_clin,
            "lrt_p": lrt,
            "score_hr": float(s["exp(coef)"]),
            "score_hr_low": float(s["exp(coef) lower 95%"]),
            "score_hr_high": float(s["exp(coef) upper 95%"]),
            "score_p": float(s["p"]),
        }
        summary_rows.append(row)
        print(f"  C-index clinical={c_clin:.3f}, score={c_score:.3f}, combined={c_combined:.3f}, "
              f"ΔC={c_combined-c_clin:+.3f}, LRT p={lrt:.2e}")

        cancer_candidates = []
        for t in CANDIDATE_TIMEPOINTS[cancer]:
            y = binary_at_t(test, t)
            eval_mask = y.notna()
            y_eval = y[eval_mask].to_numpy()
            dca_n = len(y_eval)
            event_rate = float(np.mean(y_eval)) if dca_n else np.nan
            print(f"  {t}m DCA evaluable n={dca_n}, event rate={event_rate:.3f}")
            if dca_n < 50 or len(np.unique(y_eval)) < 2:
                candidate_rows.append({
                    "cancer": cancer,
                    "timepoint_m": t,
                    "n_evaluable": dca_n,
                    "event_rate": event_rate,
                    "mean_nb_gain_vs_clinical": np.nan,
                    "selected": False,
                })
                continue

            risks = {
                "Clinical only":             predict_risk_at(cph_clin,     test.loc[eval_mask], clin, t),
                "Feature-score only":        predict_risk_at(cph_score,    test.loc[eval_mask], ["feature_score"], t),
                "Clinical + feature-score":  predict_risk_at(cph_combined, test.loc[eval_mask], clin + ["feature_score"], t),
            }
            local_rows = []
            nb_by_model = {}
            for model, risk in risks.items():
                nb = net_benefit(y_eval, risk, thresholds)
                nb_by_model[model] = np.array(nb)
                for pt, val in zip(thresholds, nb):
                    local_rows.append({"cancer": cancer, "timepoint_m": t, "model": model,
                                       "threshold": pt, "net_benefit": val,
                                       "n_evaluable": dca_n, "event_rate": event_rate})
            nb_all = [event_rate - (1 - event_rate) * (pt / (1 - pt)) for pt in thresholds]
            for pt, val in zip(thresholds, nb_all):
                local_rows.append({"cancer": cancer, "timepoint_m": t, "model": "Treat all",
                                   "threshold": pt, "net_benefit": val,
                                   "n_evaluable": dca_n, "event_rate": event_rate})
                local_rows.append({"cancer": cancer, "timepoint_m": t, "model": "Treat none",
                                   "threshold": pt, "net_benefit": 0.0,
                                   "n_evaluable": dca_n, "event_rate": event_rate})

            mean_gain = float(np.mean(nb_by_model["Clinical + feature-score"] - nb_by_model["Clinical only"]))
            cal = calibration_bins(y_eval, risks["Clinical + feature-score"], n_bins=5)
            cal["cancer"] = cancer
            cal["timepoint_m"] = t
            cal["n_evaluable"] = dca_n
            cancer_candidates.append({
                "timepoint_m": t,
                "rows": local_rows,
                "calibration": cal.to_dict("records"),
                "mean_gain": mean_gain,
                "n_evaluable": dca_n,
                "event_rate": event_rate,
            })
            candidate_rows.append({
                "cancer": cancer,
                "timepoint_m": t,
                "n_evaluable": dca_n,
                "event_rate": event_rate,
                "mean_nb_gain_vs_clinical": mean_gain,
                "selected": False,
            })

        if cancer_candidates:
            best = max(cancer_candidates, key=lambda x: (x["mean_gain"], x["n_evaluable"]))
            print(f"  selected DCA timepoint={best['timepoint_m']}m "
                  f"(mean NB gain vs clinical={best['mean_gain']:+.4f})")
            selected_dca_rows.extend(best["rows"])
            selected_cal_rows.extend(best["calibration"])
            summary_rows[-1]["selected_timepoint_m"] = best["timepoint_m"]
            summary_rows[-1]["selected_dca_mean_nb_gain_vs_clinical"] = best["mean_gain"]
            for r in candidate_rows:
                if r["cancer"] == cancer and r["timepoint_m"] == best["timepoint_m"]:
                    r["selected"] = True

    summary = pd.DataFrame(summary_rows)
    dca     = pd.DataFrame(selected_dca_rows)
    cal     = pd.DataFrame(selected_cal_rows)
    candidates = pd.DataFrame(candidate_rows)
    summary.to_csv(OUT / "clinical_utility_summary.csv", index=False)
    candidates.to_csv(OUT / "decision_curve_candidate_timepoints.csv", index=False)
    dca.to_csv(OUT / "decision_curve_selected.csv", index=False)
    cal.to_csv(OUT / "calibration_combined.csv", index=False)

    plot_forest(summary)
    if len(dca):
        plot_dca(dca)
    if len(cal):
        plot_calibration(cal)

    print(f"\nOutputs: {OUT}")


if __name__ == "__main__":
    main()
