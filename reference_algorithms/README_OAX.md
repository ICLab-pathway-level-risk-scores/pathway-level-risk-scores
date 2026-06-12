# Reference algorithm — OAX (Orthogonal Array Crossover)

Reference implementation of the orthogonal-array crossover (OAX)
operator used by IBCGA for binary feature-subset recombination.

`oax_binary_reference.py` — self-contained Python implementation using
the `L8(2^7)` orthogonal array, with a worked example.

This is provided for transparency. **It is not the full IBCGA engine**
used to derive the reported compact pathway scores. The deployed
compact scores and reported validation analyses are reproduced from
`../compact_score/compact_score_formula.yaml`, the run records under
`../ibcga_run_records/`, and the validation scripts under
`../external_validation/` and `../decision_curve_analysis/`.

## References

1. **Ho, S.-Y. & Chen, Y.-C.** "An efficient evolutionary algorithm for
   accurate polygonal approximation." *Pattern Recognition* **34**,
   2305–2317 (2001). DOI: [10.1016/S0031-3203(00)00159-X](https://doi.org/10.1016/S0031-3203(00)00159-X). *(OAX algorithm and worked example: Sections 3.2–3.3.)*

2. **Ho, S.-Y., Shu, L.-S. & Chen, J.-H.** "Intelligent evolutionary
   algorithms for large parameter optimization problems."
   *IEEE Trans. Evol. Comput.* **8**, 522–541 (2004).
   DOI: [10.1109/TEVC.2004.835176](https://doi.org/10.1109/TEVC.2004.835176). *(IEA framework.)*

3. **Ho, S.-Y., Chen, J.-H. & Huang, M.-H.** "Inheritable genetic
   algorithm for biobjective 0/1 combinatorial optimization
   problems and its applications." *IEEE Trans. Syst. Man Cybern. B*
   **34**, 609–620 (2004).
   DOI: [10.1109/TSMCB.2003.817090](https://doi.org/10.1109/TSMCB.2003.817090). *(IBCGA.)*

## Run

```bash
python reference_algorithms/oax_binary_reference.py
```
