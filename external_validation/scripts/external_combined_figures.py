"""external_combined_figures.py
================================
Build manuscript-ready combined figures pulling together all external cohorts:

  - Combined calibration plot (1 panel per cohort, 5 quintiles each)
  - Combined forest plot (adjusted HR per cohort, ranked by cancer)

Reads from:
  external_validation/outputs/_GENIE_BPC_combined/calibration_combined.csv
  figures_tables/external_multivar_cox/summary.csv
"""
from __future__ import annotations

import os
from pathlib import Path

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

ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "figures_tables" / "external_combined"
OUT_DIR.mkdir(parents=True, exist_ok=True)

CANCER_ORDER = ["BRCA", "LUAD", "PAAD", "PRAD", "CRC"]
CANCER_LABEL = {"BRCA": "IDC", "LUAD": "LUAD", "PAAD": "PAAD",
                "PRAD": "PRAD", "CRC": "CRC"}


def plot_calibration():
    p = ROOT / "external_validation" / "outputs" / "_GENIE_BPC_combined" / "calibration_combined.csv"
    if not p.exists():
        print(f"missing {p}")
        return
    df = pd.read_csv(p)
    cohorts = df["cohort"].drop_duplicates().tolist()
    ncols = min(4, len(cohorts))
    nrows = int(np.ceil(len(cohorts) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.5 * ncols, 3.5 * nrows),
                             constrained_layout=True)
    axes = np.array(axes).reshape(-1)
    for i, cohort in enumerate(cohorts):
        sub = df[df["cohort"] == cohort]
        t = sub["t_month"].iloc[0]
        cancer = sub["cancer"].iloc[0]
        ax = axes[i]
        ax.plot([0, 1], [0, 1], "k--", linewidth=0.8)
        ax.scatter(sub["mean_predicted_event"], sub["observed_event"], s=55, color="#1d3d7b")
        for _, r in sub.iterrows():
            ax.annotate(r["q"], (r["mean_predicted_event"], r["observed_event"]),
                        fontsize=8, xytext=(4, 4), textcoords="offset points")
        lim = max(sub["mean_predicted_event"].max(), sub["observed_event"].max()) + 0.07
        ax.set_xlim(0, lim); ax.set_ylim(0, lim)
        ax.set_xlabel(f"Predicted event @ {t}m")
        ax.set_ylabel(f"Observed event @ {t}m")
        ax.set_title(CANCER_LABEL.get(cancer, cancer), fontsize=10)
        ax.grid(alpha=0.2)
    for j in range(len(cohorts), len(axes)):
        axes[j].axis("off")
    fig.suptitle("GENIE BPC external calibration — compact score quintiles",
                 fontsize=12)
    fig.savefig(OUT_DIR / "calibration_combined.png", dpi=220)
    fig.savefig(OUT_DIR / "calibration_combined.pdf")
    plt.close(fig)
    print(f"wrote {OUT_DIR / 'calibration_combined.png'}")


def plot_forest():
    p = ROOT / "figures_tables" / "external_multivar_cox" / "summary.csv"
    df = pd.read_csv(p)
    df = df[df["model"] == "multivariable"].copy()

    # Map cohort label to cancer
    def cancer_of(label):
        for c in CANCER_ORDER:
            if label.startswith(c):
                return c
        return "?"
    df["cancer"] = df["cohort"].map(cancer_of)
    df["sort_key"] = df["cancer"].map({c: i for i, c in enumerate(CANCER_ORDER)})
    df = df.sort_values(["sort_key", "cohort"]).reset_index(drop=True)

    # Reverse for top-down plotting
    fdf = df.iloc[::-1].reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(9.5, 0.55 * len(fdf) + 1.8), constrained_layout=True)
    y = np.arange(len(fdf))
    # colour-code by cancer
    palette = {"BRCA": "#d6604d", "LUAD": "#4393c3", "PAAD": "#7fbf7b",
               "PRAD": "#9970ab", "CRC": "#f1a340"}
    colors = [palette.get(c, "#234") for c in fdf["cancer"]]
    ax.errorbar(fdf["HR"], y,
                xerr=[fdf["HR"] - fdf["lower95"], fdf["upper95"] - fdf["HR"]],
                fmt="o", ecolor="#789", capsize=3, markersize=7, linestyle="none",
                markerfacecolor="white", markeredgecolor="#234")
    for i, c in enumerate(colors):
        ax.plot(fdf["HR"].iloc[i], y[i], "o", color=c, markersize=7)
    ax.axvline(1.0, color="grey", linestyle="--", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.cohort}\n[{r.cancer}]  n={r.n}, ev={r.events}"
                        for r in fdf.itertuples()])
    ax.set_xlabel("Adjusted HR for compact risk_score (95% CI)")
    ax.set_xscale("log")
    ax.set_title("External multivariable Cox — compact pathway score per cohort\n"
                 "(adjusted for available age/sex/stage/treatment/sample type)",
                 fontsize=10)
    xmax = ax.get_xlim()[1]
    for i, r in enumerate(fdf.itertuples()):
        sig = "**" if r.p < 0.01 else ("*" if r.p < 0.05 else "")
        ax.text(xmax, i,
                f"  HR={r.HR:.2f} ({r.lower95:.2f}–{r.upper95:.2f}), p={r.p:.2g} {sig}",
                va="center", fontsize=8)
    fig.savefig(OUT_DIR / "forest_combined.png", dpi=220)
    fig.savefig(OUT_DIR / "forest_combined.pdf")
    plt.close(fig)
    print(f"wrote {OUT_DIR / 'forest_combined.png'}")

    # also write summary table sorted by cancer
    df[["cohort", "cancer", "n", "events", "HR", "lower95", "upper95", "p",
        "cindex", "covariates"]].to_csv(OUT_DIR / "summary_sorted.csv", index=False)


def main():
    plot_calibration()
    plot_forest()


if __name__ == "__main__":
    main()
