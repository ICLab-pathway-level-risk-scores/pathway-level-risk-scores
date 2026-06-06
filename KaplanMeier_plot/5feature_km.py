#!/usr/bin/env python3
"""KM plots for the Cox risk-score feature sets without HIT_COUNT comparison.

Uses the same feature sets and test-set rank-quantile grouping logic as
compare_3feature_vs_hitcount_km.py, but only plots the Cox risk score.
The x-axis is chosen per cancer from at-risk counts to avoid sparse tail drops.
"""
from pathlib import Path
import os
import warnings

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib-cache")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from figure_style import apply_arial_style
apply_arial_style()
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "model_input" / "filter1_age_fixed"
OUT = ROOT / "figures_tables" / "5feature_km"
CANCERS = ["BRCA", "LUAD", "PAAD", "PRAD", "CRC"]

CANCER_LABEL = {"BRCA": "IDC", "LUAD": "LUAD", "PAAD": "PAAD",
                "PRAD": "PRAD", "CRC": "CRC"}

FEATURES = {
    "BRCA": ["PW_HIT_COUNT", "PW_TP53_mut_rate_z", "PW_Cell_Cycle_sv_hit", "PW_Chromatin_zsum", "PW_RTK_RAS_amp_hit"],
    "LUAD": ["PW_HIT_COUNT", "PW_TP53_any_rate_z", "PW_NOTCH_any_rate_z", "PW_NRF2_zsum", "PW_RTK_RAS_mut_rate_z"],
    "PAAD": ["PW_HIT_COUNT", "PW_RTK_RAS_amp_hit", "PW_TP53_mut_rate_z", "PW_TGF_Beta_any_rate_z", "PW_MYC_any_rate_z"],
    "PRAD": ["PW_HIT_COUNT", "PW_Chromatin_sv_hit", "PW_Cell_Cycle_zsum", "PW_WNT_zsum"],
    "CRC": ["PW_HIT_COUNT", "PW_DDR_zsum", "PW_RTK_RAS_any_rate_z", "PW_TGF_Beta_amp_hit", "PW_Chromatin_any_rate_z"],
}

COLORS = {
    2: ["#2166ac", "#d6604d"],
    3: ["#2166ac", "#878787", "#d6604d"],
    4: ["#2166ac", "#92c5de", "#f4a582", "#d6604d"],
}


def display_feature(name: str) -> str:
    if name == "PW_HIT_COUNT":
        return "All_Pathway_Alternation_Count"
    if name.startswith("PW_"):
        name = name[3:]
    if name.endswith("_z"):
        name = name[:-2]
    return name


def display_feature_list(features: list[str]) -> str:
    return ", ".join(display_feature(f) for f in features)


def rank_quantile_groups(values, q: int) -> np.ndarray:
    ranked = pd.Series(values).astype(float).rank(method="first")
    return pd.qcut(ranked, q=q, labels=False).astype(int).to_numpy()


def fit_score(train: pd.DataFrame, test: pd.DataFrame, features: list[str]):
    cols = features + ["Event_OS", "OS_MONTHS"]
    cph = CoxPHFitter(penalizer=0.1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cph.fit(train[cols], duration_col="OS_MONTHS", event_col="Event_OS")
    score = cph.predict_partial_hazard(test[features]).to_numpy()
    return cph, score


def logrank_p(duration, event, groups) -> float:
    result = multivariate_logrank_test(duration, groups, event)
    return float(result.p_value)


def km_survival_at(duration, event, groups, q: int, month: int = 24) -> list[float]:
    out = []
    kmf = KaplanMeierFitter()
    for g in range(q):
        mask = groups == g
        kmf.fit(duration[mask], event_observed=event[mask])
        out.append(float(kmf.predict(month)))
    return out


def monotonic_decreasing(vals, tolerance: float = 1e-9) -> bool:
    return all(vals[i] >= vals[i + 1] - tolerance for i in range(len(vals) - 1))


def auto_xlim(duration, score) -> int:
    """Pick a cancer-level x-axis cap based on the sparsest rank-quantile group.

    The cap is the latest 5-month grid point where every group across Q2/Q3/Q4
    still has at least max(8, 4% of group size) samples at risk. It is rounded
    down to a clean 10-month boundary for consistent panel layout.
    """
    duration = np.asarray(duration, dtype=float)
    upper = min(120.0, float(np.nanmax(duration)))
    grid = np.arange(20, np.floor(upper / 5) * 5 + 0.1, 5)
    best = grid[0] if len(grid) else 40
    for t in grid:
        ok = True
        for q in [2, 3, 4]:
            groups = rank_quantile_groups(score, q)
            for g in range(q):
                mask = groups == g
                min_at_risk = max(8, int(np.ceil(mask.sum() * 0.04)))
                if int((duration[mask] >= t).sum()) < min_at_risk:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            best = t
    return int(max(40, np.floor(best / 10) * 10))


def plot_panel(ax, duration, event, score, q: int, x_max: int):
    groups = rank_quantile_groups(score, q)
    p = logrank_p(duration, event, groups)
    kmf = KaplanMeierFitter()
    for g in range(q):
        mask = groups == g
        label = f"Q{g + 1} (n={int(mask.sum())})"
        kmf.fit(duration[mask], event_observed=event[mask], label=label)
        kmf.plot_survival_function(
            ax=ax,
            ci_show=(q == 2),
            color=COLORS[q][g],
            linewidth=2.0,
            alpha=0.92,
        )
    ax.set_title(f"{q}-quantile", fontsize=10)
    ax.set_xlabel("OS months", fontsize=9)
    ax.set_ylabel("Survival probability", fontsize=9)
    ax.set_xlim(0, x_max)
    ax.set_ylim(0, 1.02)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7, loc="best", framealpha=0.72, edgecolor="gray")
    ax.text(
        0.03,
        0.05,
        f"log-rank p={p:.2e}",
        transform=ax.transAxes,
        fontsize=8,
        bbox=dict(boxstyle="round,pad=0.22", facecolor="white", alpha=0.82, edgecolor="gray"),
    )
    return groups, p


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    hr_rows = []

    for cancer in CANCERS:
        train = pd.read_csv(DATA / cancer / f"{cancer}_train.csv")
        test = pd.read_csv(DATA / cancer / f"{cancer}_test.csv")
        features = FEATURES[cancer]
        duration = test["OS_MONTHS"].to_numpy()
        event = test["Event_OS"].to_numpy()

        cph, score = fit_score(train, test, features)
        x_max = auto_xlim(duration, score)

        for feature in features:
            s = cph.summary.loc[feature]
            hr_rows.append(
                {
                    "cancer": cancer,
                    "feature": feature,
                    "HR": float(s["exp(coef)"]),
                    "coef": float(s["coef"]),
                    "p": float(s["p"]),
                }
            )

        fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.6))
        fig.suptitle(
            f"{CANCER_LABEL.get(cancer, cancer)}: Cox risk score Q2/Q3/Q4 KM\n"
            f"Features: {display_feature_list(features)}",
            fontsize=11,
            fontweight="bold",
            y=1.02,
        )

        for ax, q in zip(axes, [2, 3, 4]):
            groups, p = plot_panel(ax, duration, event, score, q, x_max)
            surv24 = km_survival_at(duration, event, groups, q, 24)
            rows.append(
                {
                    "cancer": cancer,
                    "quantile": q,
                    "logrank_p": p,
                    "surv24_by_group": ";".join(f"{v:.3f}" for v in surv24),
                    "surv24_monotonic_decreasing": monotonic_decreasing(surv24),
                    "n_by_group": ";".join(str(int((groups == g).sum())) for g in range(q)),
                    "x_axis_months": x_max,
                    "features": ";".join(features),
                }
            )

        plt.tight_layout()
        fig.savefig(OUT / f"{cancer}_risk_score_km.png", dpi=220, bbox_inches="tight")
        fig.savefig(OUT / f"{cancer}_risk_score_km.pdf", bbox_inches="tight")
        plt.close(fig)

    pd.DataFrame(rows).to_csv(OUT / "risk_score_km_summary.csv", index=False)
    pd.DataFrame(hr_rows).to_csv(OUT / "risk_score_hr_summary.csv", index=False)
    print(OUT)


if __name__ == "__main__":
    main()
