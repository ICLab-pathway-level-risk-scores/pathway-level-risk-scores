# Reference algorithm — OAX (Orthogonal Array Crossover)

This folder provides a **reference implementation** of the
orthogonal-array crossover (OAX) operator used by the Inheritable
Bi-objective Combinatorial Genetic Algorithm (IBCGA) for fixed-K binary
feature-subset recombination.

## What is OAX?

OAX is a recombination operator originally developed in the
Intelligent Evolutionary Algorithm (IEA) framework. It uses a Taguchi
orthogonal array to systematically sample a small but
information-rich subset of all possible parental-chunk
recombinations. Compared with single-point or uniform crossover, OAX:

- evaluates only `α` offspring (e.g. 8 with `L8(2^7)`) instead of the
  combinatorial `2^B` parental-chunk assignments,
- estimates a per-factor **main effect** by Taguchi factor analysis,
  which IBCGA uses to bias subsequent search toward informative
  chunks.

In this study, OAX is the binary recombination operator inside
IBCGA. The outer IBCGA loop progressively reduces the selected
feature count from `r_start` to `r_end` (per-cancer values are
reported in Supplementary Table S5). Within each outer cycle, OAX
operates at a **fixed K**, i.e. every individual in the population has
exactly K ones in its binary feature mask.

## What this folder is — and is not

This folder contains a **reference implementation** for transparency:

- `oax_binary_reference.py` — a self-contained Python implementation
  of OAX for fixed-K binary feature-subset recombination, using the
  `L8(2^7)` orthogonal array, with a small worked example.

This folder is **not** the full IBCGA engine used to derive the
reported compact pathway scores. The original IBCGA engine is not
redistributed in this repository.

## Why is it sufficient for reproducibility?

The reported compact pathway scores and all reported validation
analyses are reproduced from:

- `compact_score/compact_score_formula.yaml` — fixed Cox β
  coefficients, MSK training-set z-score parameters, and MSK-derived
  cutoffs;
- `external_validation/` — no-refit GENIE BPC Cox, calibration, and
  sequencing-origin sensitivity analyses;
- `decision_curve_analysis/` — DCA, trial-enrichment, MDT-referral,
  and fixed-cutoff drift analyses;
- `ibcga_run_records/` — archived 30-run IBCGA selected-feature sets
  and aggregate selection-frequency statistics, sufficient to audit
  the feature-selection decision.

Rerunning IBCGA from scratch is not required to reproduce any
reported number.

## References

1. **Ho, S.-Y. & Chen, Y.-C.** "An efficient evolutionary algorithm
   for accurate polygonal approximation." *Pattern Recognition* **34**,
   2305–2317 (2001). DOI: [10.1016/S0031-3203(00)00159-X](https://doi.org/10.1016/S0031-3203(00)00159-X).
   *(Origin paper for OAX; Sections 3.2–3.3 give the algorithm and a
   worked example.)*

2. **Ho, S.-Y., Shu, L.-S. & Chen, J.-H.** "Intelligent evolutionary
   algorithms for large parameter optimization problems."
   *IEEE Trans. Evol. Comput.* **8**, 522–541 (2004).
   DOI: [10.1109/TEVC.2004.835176](https://doi.org/10.1109/TEVC.2004.835176).
   *(IEA framework; OAX as the IEA recombination operator.)*

3. **Ho, S.-Y., Chen, J.-H. & Huang, M.-H.** "Inheritable genetic
   algorithm for biobjective 0/1 combinatorial optimization
   problems and its applications." *IEEE Trans. Syst. Man Cybern. B*
   **34**, 609–620 (2004).
   DOI: [10.1109/TSMCB.2003.817090](https://doi.org/10.1109/TSMCB.2003.817090).
   *(IBCGA: the bi-objective inheritable GA that uses OAX as its
   binary feature-subset recombination operator.)*

## Run

```bash
python reference_algorithms/oax_binary_reference.py
```

The script prints two random fixed-K parents, the 8 OAX-generated
offspring, their toy fitness values, the per-factor Taguchi main
effects, and the indices of the best two offspring (the standard
IBCGA selection-for-survival rule).
