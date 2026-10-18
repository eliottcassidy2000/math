---
id: HYP-9316
title: "First homotopy of the derived arc space of a quasi-homogeneous plane curve: the first weight at which the arc-space Hilbert series leaves the derived Euler product prod(1 - t^(d+n)) / prod((1 - t^(a+j))(1 - t^(b+j))) equals the minimal weight of a nonzero Cramer (determinantal) syzygy among jets f_0..f_K lying in a K-generated monomial ideal of jet variables. This explains kind-pasteur's cusp observation (agreement through weight 33, failure at 34): on the lct contact locus x_0 = x_1 = y_0 = y_1 = y_2 = 0 six jets lie in a five-generator ideal, and 6+7+8+9+10+11 - (2+3+3+4+5) = 34"
status: >
  OPEN, strongly supported.
  - Two pre-registered predictions PASSED. They were committed in 06b7a0fe91 before the runs:
    - E6 y^3 - x^4 (wt 3,4): first deviation at weight 29 (2166 vs 2165);
    - A4 y^2 - x^5 (wt 2,5): first deviation at weight 31 (4717 vs 4716).
  - Matches observed before registration: D4 x^3 + y^3 and x^2 y + y^3 at 10; A3 y^2 - x^4 (wt 1,2) at 16;
    x^4 + y^4 at 11, via the non-contact ideal (x_0^2, y_0^2).
  - kind-pasteur's cusp (wt 2,3): first failure at 34 = the contact-locus Cramer weight.
  - Control: the node xy (log canonical) never deviates (checked to weight 16).
  - Each deviation is by exactly 1 at the first weight, consistent with a single H_1 class.
  - The Cramer syzygy exists by construction once the coefficient matrix has full rank K (checked by random
    evaluation). That it is the FIRST class and is not a Koszul boundary is the conjectural part.
session: opus-2026-10-18-S44
source: 05-knowledge/results/li_spectral_sequence_quantum_transfer_20261018.md
scripts:
  - 04-computation/quantum_transfer_s44/derived_arc_curves_first_syzygy_20261018.py (+ curves_E6_W30.out, curves_A4_W32.out, derived_arc_curves_first_syzygy_*.json)
related: [HYP-9298, kind-pasteur-2026-10-18-S1 note sections 1.3 and 5a, Mustata 2001 (jet schemes of lci singularities), Ein-Mustata-Yasuda]
---

# HYP-9316 — the first homotopy of a curve's derived arc space is a Cramer syzygy

## Setting

Let `f(x, y)` be quasi-homogeneous of weight `d` (`wt x = a`, `wt y = b`).
- Its derived arc space has functions the Koszul complex of the jets `f_n = d^n f` (weight `d + n`) on
  `C[x_j, y_j]` (weights `a + j`, `b + j`).
- The Euler characteristic is `E(t) = Π_(n>=0)(1 - t^(d+n)) / Π_(j>=0)((1 - t^(a+j))(1 - t^(b+j)))`.
- `H_0` is the arc algebra, and `HS(H_0) = E` exactly below the first weight `w_1` with `H_1 ≠ 0`.
- For the cusp (`a, b, d = 2, 3, 6`), `E = Π_(n>=2)(1 - t^n)^(-1) / ((1-t^3)(1-t^4)(1-t^5))`. kind-pasteur found it
  equal to `HS(H_0)` through weight 33 and different at 34 (13436 vs 13435).

## Mechanism

Suppose the `K + 1` jets `f_0, ..., f_K` lie in an ideal `Q = (h_1, ..., h_K)` with `K` homogeneous generators.
- Write `f_l = Σ_j A_(lj) h_j`.
- Cramer's rule gives `Σ_l (-1)^l det(A_l) f_l = 0`, where `A_l` is `A` with row `l` deleted.
- This syzygy has weight `w(Q) = Σ_(l<=K) (d + l) - Σ_j wt(h_j)`, and it is nonzero when `A` has rank `K`.

**Cusp.** On the contact locus `x_0 = x_1 = y_0 = y_1 = y_2 = 0` (arcs with `ord x >= 2`, `ord y >= 3`), `f` vanishes
to order 6. So the six jets `f_0..f_5` lie in a five-generator ideal, and
`w = (6+7+...+11) - (2+3+3+4+5) = 51 - 17 = 34`. This is the contact locus computing
`lct = codim/order = 5/6 < 1`.

## Statement

`w_1 = min w(Q)`, the minimum taken over monomial ideals `Q` in the jet variables that contain `K + 1` consecutive jets,
`K` being the number of generators, and whose coefficient matrix has rank `K`.

## Evidence

| curve | weights | `d` | predicted `w(Q)` (ideal) | observed `w_1` |
|---|---|---|---|---|
| node `xy` | (1,1) | 2 | none | none to weight 16 |
| D4 `x^3 + y^3`; `x^2 y + y^3` | (1,1) | 3 | 10 `(x_0, y_0)` | 10 |
| A3 `y^2 - x^4` | (1,2) | 4 | 16 `(x_0, y_0, y_1)` | 16 |
| `x^4 + y^4` | (1,1) | 4 | 11 `(x_0^2, y_0^2)` | 11 |
| E6 `y^3 - x^4` | (3,4) | 12 | 29 `(x_0^2, y_0)` | **29 (pre-registered)** |
| A4 `y^2 - x^5` | (2,5) | 10 | 31 `(x_0^2, y_0, y_1)` | **31 (pre-registered)** |
| cusp `y^2 - x^3` | (2,3) | 6 | 34 `(x_0, x_1, y_0, y_1, y_2)` | 34 (kind-pasteur) |

The pure contact-locus predictor, which allows only ideals `(x_0..x_(i-1), y_0..y_(k-1))`, is right for D4, A3 and the
cusp. It is too large for `x^4 + y^4` (13 vs 11), E6 (32 vs 29) and A4 (33 vs 31), where powers `x_0^2` are needed.

## Context

- **Log canonical curves.** By Mustață (2001) and Ein-Mustață-Yasuda, jet schemes of a locally complete intersection
  have the expected dimension for all `m` exactly when the singularity is log canonical. For a plane curve this means
  the jets form a regular sequence, `H_(>=1) = 0`, and `HS(H_0) = E` forever (the node). A non-log-canonical curve has
  a contact locus with `codim < order`, which carries surplus jets in a smaller ideal, which gives a Cramer syzygy.
- **The weighting.** `w_1` is the weighted position of the first homotopy class. It depends on both the singularity
  and the grading.
- **Fat points.** For the fat point `x^N` (kind-pasteur S1, section 1.3) the same mechanism, with the principal ideal
  `(x_0^(N-1))` containing `f_0, f_1`, gives the Euler syzygy at weight `2N + 3`.

## Falsifiers

- A curve where the first deviation is below the predicted minimum.
- A curve where the predicted Cramer syzygy is a Koszul boundary, so the first deviation lies above the prediction.
