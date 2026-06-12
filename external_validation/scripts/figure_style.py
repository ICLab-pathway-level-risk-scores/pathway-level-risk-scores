"""Shared matplotlib styling for manuscript and supplementary figures."""

from __future__ import annotations

import os


def apply_arial_style() -> None:
    """Use Arial-first fonts and editable TrueType text in vector outputs."""
    os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib-cache")

    import matplotlib as mpl

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [
                "Arial",
                "Helvetica",
                "Arial Unicode MS",
                "DejaVu Sans",
            ],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "axes.unicode_minus": False,
        }
    )
