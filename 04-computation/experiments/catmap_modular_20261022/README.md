# The cat map in modular data, the logarithmic side, and the Rogers–Ramanujan family (opus-2026-10-22-S55)

These scripts support:
* [THM-4762](../../../01-canon/theorems/THM-4762-the-arnold-cat-map-is-the-commutator-of-t-and-s-pisano-bounds-in-modular-data.md);
* [THM-4763](../../../01-canon/theorems/THM-4763-the-cat-map-on-the-logarithmic-side-lucas-traces-and-spectral-radius-phi-squared.md);
* [HYP-9372](../../../05-knowledge/hypotheses/HYP-9372-pisano-saturation-of-the-cat-map-and-the-ising-spin-structure-exception.md);
* [HYP-9373](../../../05-knowledge/hypotheses/HYP-9373-lorentzian-growth-equals-logarithmic-degree-for-c2-cofinite-vertex-algebras.md);
* the [synthesis note](../../../05-knowledge/results/signature_axis_cat_map_rr_mrb_20261022.md).

They need Python 3 with numpy, sympy and mpmath. No supplied PDF or extracted text is stored here.

| script | output | what |
|---|---|---|
| `catmap_mtc.py` | `catmap_mtc.out`, `catmap_census.json` | Pisano facts. Modular data for Fibonacci, Virasoro minimal models, `SU(2)_k`, `SU(3)_k`, lattice discriminant forms (definite and indefinite), `D(Z_N)`, `D(S_3)` and a Deligne product. Census of `C = T S T⁻¹ S⁻¹` |
| `center_check.py` | `center_check.out` | faithfulness of the cat map modulo the centre of `SL2(Z/N)` |
| `log_catmap.py` | `log_catmap.out` | the 3p−1 W(p) torus amplitudes; S, T computed at 40 digits; eigenvalue moduli and the Lucas trace law, `p = 2..6` |
| `rr_family.py` | `rr_family.out` | exact verification to `q^300` of 18 Rogers–Ramanujan-type identities and of the fifth-order mock theta Hecke-type identity |
| `rr_asymptotics.py` | `rr_asymptotics.out` | effective central charges and q → 1 ratios |
| `mrb_checks.py` | `mrb_checks.out` | checks of the supplied MRB documents (FIT bounds, prime defect, Ising-ring identity) and of the Collatz abstract's trailing-ones claim |
