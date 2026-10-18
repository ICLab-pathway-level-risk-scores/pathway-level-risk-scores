"""Import shim for the public repository layout.

The analysis scripts import the pathway-feature builder under its working-copy module
name, build_3flag_enhanced. In this repository the same code is
preprocessing/build_pathway_features.py; this shim re-exports it unchanged.
"""
import sys
from pathlib import Path

_here = Path(__file__).resolve()
_root = next(p for p in _here.parents if (p / "preprocessing" / "build_pathway_features.py").exists())
sys.path.insert(0, str(_root))

from preprocessing.build_pathway_features import *  # noqa: E402,F401,F403
