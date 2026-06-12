"""
OAX (Orthogonal Array Crossover) — fixed-K binary feature-subset reference
=========================================================================

This is a *reference implementation* of the orthogonal array crossover
(OAX) operator used inside IBCGA for binary feature-subset recombination
under a fixed selected-feature count K.

This file is provided for transparency. It is NOT the full IBCGA engine
used to derive the reported compact pathway scores. The reported
results are reproduced from the fixed YAML score formula
(`compact_score/compact_score_formula.yaml`) and the validation scripts
under `external_validation/` and `decision_curve_analysis/`; rerunning
IBCGA is not required.

References
----------
1. Ho, S.-Y. & Chen, Y.-C. "An efficient evolutionary algorithm for
   accurate polygonal approximation." Pattern Recognition 34,
   2305-2317 (2001). https://doi.org/10.1016/S0031-3203(00)00159-X
   (OAX algorithm and worked example: Sections 3.2-3.3.)

2. Ho, S.-Y. et al. "Intelligent evolutionary algorithms for large
   parameter optimization problems." IEEE Trans. Evol. Comput. 8,
   522-541 (2004). https://doi.org/10.1109/TEVC.2004.835176
   (IEA / OAX framework.)

3. Ho, S.-Y. et al. "Inheritable genetic algorithm for biobjective 0/1
   combinatorial optimization problems and its applications."
   IEEE Trans. Syst. Man Cybern. B 34, 609-620 (2004).
   https://doi.org/10.1109/TSMCB.2003.817090
   (IBCGA: the bi-objective inheritable GA that uses OAX as its
   binary recombination operator.)

Operator semantics in this study
--------------------------------
For each cancer and outer IBCGA cycle, the population consists of
binary indicator vectors of length n_features. Each individual has
exactly K ones (and n_features - K zeros), representing a selected
feature subset of size K. OAX recombines two such parents and
produces multiple offspring, then selects the best two offspring by
fitness for the next generation. The orthogonal-array structure
samples a small but information-rich subset of all 2^B possible
chunk-assignment combinations.

The outer IBCGA loop decreases K stepwise from r_start to r_end
(per-cancer values are reported in Supplementary Table S5). This
file implements OAX for one fixed K.

Notes on this reference implementation
--------------------------------------
* This implementation uses the L8(2^7) orthogonal array, which is
  appropriate for B = 7 chunks. For more chunks, an L16(2^15) array
  could be substituted.
* This implementation includes a simple "fixed-K repair" step that
  restores exactly K ones in each offspring by random bit flips on
  positions inherited from the parents. The exact repair / feasibility
  policy in the production IBCGA engine may differ.
* This file is self-contained and only depends on numpy.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


# ---------------------------------------------------------------------
# L8(2^7) orthogonal array (Taguchi)
# rows = 8 experimental conditions, columns = 7 two-level factors.
# Each row tells one offspring which parent (0 or 1) to inherit each
# chunk from.
# ---------------------------------------------------------------------
L8 = np.array(
    [
        [0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 1, 1, 1],
        [0, 1, 1, 0, 0, 1, 1],
        [0, 1, 1, 1, 1, 0, 0],
        [1, 0, 1, 0, 1, 0, 1],
        [1, 0, 1, 1, 0, 1, 0],
        [1, 1, 0, 0, 1, 1, 0],
        [1, 1, 0, 1, 0, 0, 1],
    ],
    dtype=np.uint8,
)


# ---------------------------------------------------------------------
# Helper structures
# ---------------------------------------------------------------------
@dataclass
class OAXResult:
    offspring: np.ndarray          # (8, n_features) uint8, all with exactly K ones
    fitness: np.ndarray            # (8,) float, evaluated by user-provided fn
    main_effects: np.ndarray       # (B,) per-factor estimated main effect
    best_two_idx: np.ndarray       # (2,) indices of best two offspring


# ---------------------------------------------------------------------
# Core OAX function
# ---------------------------------------------------------------------
def oax_binary_recombination(
    parent_a: np.ndarray,
    parent_b: np.ndarray,
    K: int,
    fitness_fn,
    n_chunks: int = 7,
    rng: np.random.Generator | None = None,
) -> OAXResult:
    """OAX recombination of two binary parents under fixed-K feasibility.

    Parameters
    ----------
    parent_a, parent_b
        1D binary arrays of length n_features. Each must have exactly K ones.
    K
        Required number of ones in every offspring (fixed-K constraint).
    fitness_fn
        Callable that takes a binary vector (or batch (B, n_features)) and
        returns a scalar fitness (higher is better, e.g. cross-validated
        Cox C-index). Will be called on each offspring.
    n_chunks
        Number of factors B used in the OA. With L8 we use B = 7.
    rng
        Optional numpy random Generator for reproducibility.

    Returns
    -------
    OAXResult containing 8 offspring (one per OA row), their fitness,
    per-factor main effects, and the indices of the best two offspring.

    Notes
    -----
    The differing positions between parent_a and parent_b are partitioned
    into n_chunks contiguous-ish chunks. Each OA row assigns each chunk
    to parent_a (0) or parent_b (1). The resulting offspring is then
    repaired to have exactly K ones.
    """
    if rng is None:
        rng = np.random.default_rng()
    parent_a = np.asarray(parent_a, dtype=np.uint8)
    parent_b = np.asarray(parent_b, dtype=np.uint8)
    assert parent_a.shape == parent_b.shape
    assert parent_a.sum() == K and parent_b.sum() == K, (
        "Both parents must have exactly K ones (fixed-K feasibility)."
    )
    n_features = parent_a.size
    if n_chunks > 7:
        raise ValueError("This reference implementation uses L8, max 7 chunks.")

    # Step 1. identify differing positions between parents
    diff_mask = parent_a != parent_b
    diff_idx = np.flatnonzero(diff_mask)
    same = np.where(diff_mask, 0, parent_a).astype(np.uint8)  # copied to all offspring

    # Step 2. partition differing positions into n_chunks
    # If fewer differing positions than chunks, use degenerate chunks of size 1.
    if diff_idx.size == 0:
        # parents identical: return 8 copies of parent_a as offspring
        offspring = np.tile(parent_a, (8, 1))
    else:
        chunk_assignment = np.array_split(diff_idx, min(n_chunks, diff_idx.size))
        offspring = np.empty((8, n_features), dtype=np.uint8)
        for row_i, oa_row in enumerate(L8):
            child = same.copy()
            for chunk_i, chunk_positions in enumerate(chunk_assignment):
                source = parent_b if oa_row[chunk_i] else parent_a
                child[chunk_positions] = source[chunk_positions]
            offspring[row_i] = _repair_to_K(child, K, rng)

    # Step 3. evaluate fitness
    fitness = np.array([float(fitness_fn(child)) for child in offspring], dtype=float)

    # Step 4. estimate per-factor main effect (Taguchi MED-style)
    main_effects = np.zeros(min(n_chunks, max(1, diff_idx.size)), dtype=float)
    used_factors = main_effects.size
    for f in range(used_factors):
        col = L8[:, f]
        mean_lo = fitness[col == 0].mean() if (col == 0).any() else 0.0
        mean_hi = fitness[col == 1].mean() if (col == 1).any() else 0.0
        main_effects[f] = mean_hi - mean_lo

    # Step 5. select best two offspring by fitness
    best_two_idx = np.argsort(-fitness)[:2]

    return OAXResult(
        offspring=offspring,
        fitness=fitness,
        main_effects=main_effects,
        best_two_idx=best_two_idx,
    )


def _repair_to_K(child: np.ndarray, K: int, rng: np.random.Generator) -> np.ndarray:
    """Restore exactly K ones in `child` by random flips of inherited bits."""
    s = int(child.sum())
    if s == K:
        return child
    child = child.copy()
    if s > K:
        # flip s - K random ones to zeros
        on_pos = np.flatnonzero(child == 1)
        to_flip = rng.choice(on_pos, size=s - K, replace=False)
        child[to_flip] = 0
    else:
        # flip K - s random zeros to ones
        off_pos = np.flatnonzero(child == 0)
        to_flip = rng.choice(off_pos, size=K - s, replace=False)
        child[to_flip] = 1
    return child


# ---------------------------------------------------------------------
# Worked example (mirrors Section 3.3 of Ho & Chen, Pattern Recognition 2001
# in spirit: small problem, easy to inspect)
# ---------------------------------------------------------------------
def _toy_fitness(child: np.ndarray) -> float:
    """Toy fitness: prefers ones in the first half of the feature vector."""
    n = child.size
    weights = np.linspace(1.0, 0.0, n)
    return float(np.dot(child, weights))


def _demo() -> None:
    rng = np.random.default_rng(seed=2024)
    n_features = 20
    K = 5
    # two random fixed-K parents
    pa = np.zeros(n_features, dtype=np.uint8)
    pa[rng.choice(n_features, K, replace=False)] = 1
    pb = np.zeros(n_features, dtype=np.uint8)
    pb[rng.choice(n_features, K, replace=False)] = 1
    print("parent_a:", pa)
    print("parent_b:", pb)
    print(f"|diff positions|: {int((pa != pb).sum())}")

    result = oax_binary_recombination(pa, pb, K, _toy_fitness, n_chunks=7, rng=rng)
    print("\nOffspring (8 rows from L8):")
    for i, child in enumerate(result.offspring):
        marker = " <-- best" if i in result.best_two_idx else ""
        print(f"  row {i}: {child}  fitness={result.fitness[i]:.3f}  ones={int(child.sum())}{marker}")
    print("\nMain effects per OA factor:", np.round(result.main_effects, 3))
    print("Best-two indices:", result.best_two_idx.tolist())


if __name__ == "__main__":
    _demo()
