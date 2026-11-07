---
id: HYP-9490
title: "Gate towers forget their core, except through the cyclic-triangle count: for tournaments T, T' of the same order n, H(Gamma^r T)/H(Gamma^r T') = 1 + 2(c3(T) - c3(T')) rho_(n,r) (1 + o(1)) as r -> infinity, where rho_(n,r) = H(Gamma^r TT_(n-2))/H(Gamma^r TT_n) -> 0 is universal; equivalently the universal tower weights U_(d,r) = (L^r q^d)(1) of the gate's cover-polynomial operator satisfy U_(d-2,r)/U_(d,r) -> 0, with rho depending asymptotically only on the total order N = n + 2r, of size about (log N)^2/N^2"
status: >
  OPEN (CONJECTURE) with an exact reduction.
  PROVED (THM-4900(6)):
  - H(Gamma^r T) = sum_j sw_j(T) U_(n-2j,r), with sw_1 = 2 c3;
  - hence H(Gamma^r T) - H(Gamma^r T') = 2(c3(T) - c3(T')) U_(n-2,r) + sum_(j>=2) (sw_j(T) - sw_j(T')) U_(n-2j,r) exactly.
  So the conjecture is equivalent to an analytic statement about one explicit linear operator L
  (L f(q) = q[f(q+1) + q f(q) - f(q-1)], Session210), with no tournaments left in it.
  FINITE-EXACT evidence:
  - every isomorphism class of orders 4, 5, 6, 7 (4 + 12 + 56 + 456 cores) at total order N = 401/402: the relative
    deviation from the reference tower is kappa(T) eps_N with kappa = c3(T) - c3(reference), residual <= 1.7e-6 against
    eps_N = 3.95e-4; at order 4 the law is exact (P_T = q^4 + 2 c3 q^2);
  - six cores (TT3, C3, TT5, R5, QR7, TT7) to N = 1201.
session: mac-mini-2026-11-06-microcosm
source: 05-knowledge/results/microcosm_gate_doubling_resolvent_20261106.md
scripts:
  - 04-computation/experiments/microcosm_gate_resolvent_20261106/c3_first_order_law.py (+ .out)
  - 04-computation/experiments/microcosm_gate_resolvent_20261106/cover_tower_values.py (+ .out)
related: [THM-4900, Session210 graded-spectrum-gates (operator, Theorem 4/5), THM-078 (odd-cycle form), klein THM-4889/HYP-9481 (levels by cyclic triangles), Session210 structure atlas (q = (H-1)/2)]
---

# HYP-9490 — gate towers forget their core, except through the cyclic triangles

## Statement

Let L act on polynomials by L f(q) = q[f(q+1) + q f(q) − f(q−1)]. This is the gate's action on the ordered
path-cover polynomial (Session210). Put

    U_(d,r) = (L^r q^d)(1) = H(Γ^r TT_d),    ρ_(d,r) = U_(d−2,r)/U_(d,r).

**Conjecture.**
- **(a)** For each fixed d ≥ 2, ρ_(d,r) → 0 as r → ∞.
- **(b)** ρ depends asymptotically only on N = d + 2r: ρ_(d,r)/ρ_(d+2,r−1) → 1, and U_(d,r)/U_(d+2,r−1) → 1.
- **(c)** ρ_(d,r) = Θ((log N)²/N²). Exact values give 2N²ρ = 63.7, 77.8, 86.8 at N = 401, 801, 1201.

**Consequence**, PROVED from (a) and THM-4900(6). Let T, T′ have the same order n. Then

    H(Γ^r T)/H(Γ^r T′) = 1 + 2(c3(T) − c3(T′))·ρ_(n,r)·(1 + o(1)).

- So deep gate towers forget their core. To first order the only trace left is the number of cyclic triangles, i.e.
  the metagraph level of klein's THM-4889.
- With (b), the same holds against the one-vertex tower of the same total order:

      H(Γ^r T)/H(Γ^((N−1)/2) K1) = 1 + 2(c3(T) − ((n−1)/2)²)·ρ_N·(1 + o(1))    (n odd),

  because Γ^((n−1)/2)(K1) has c3 = ((n−1)/2)². For even n, use the tower over TT2.

## Evidence (exact integer arithmetic)

Tower values come from the cover-polynomial value recursion P_Γ(j) = j[P(j+1) + jP(j) − P(j−1)] at integer points,
which reproduces brute-force Hamiltonian-path counts to order 19.

**All cores of a given order, at N ≈ 401** (c3_first_order_law.out):

| Order n | Classes | Reference | Spread of dev/κ | Max \|dev − κ·ε\| |
|---:|---:|---|---|---:|
| 4 | 4 | TT2 tower, c3 = 2 | exact | 1.1e-19 |
| 5 | 12 | K1 tower, c3 = 4 | [3.9595e-4, 3.9619e-4] | 1.65e-7 |
| 6 | 56 | TT2 tower, c3 = 6 | [3.9397e-4, 3.9499e-4] | 5.57e-7 |
| 7 | 456 | K1 tower, c3 = 9 | [3.9429e-4, 3.9754e-4] | 1.74e-6 |

**Large N against the one-vertex tower** (cover_tower_values.out). Relative deviations, with N²·dev in brackets:

| Core | κ | N = 401 | N = 801 | N = 1201 |
|---|---:|---|---|---|
| TT3 | −1 | −3.96e-4 (−63.7) | −1.21e-4 (−77.8) | −6.02e-5 (−86.8) |
| R5 | +1 | +3.96e-4 (+63.7) | +1.21e-4 (+77.8) | +6.02e-5 (+86.8) |
| TT5 | −4 | −1.58e-3 (−254.8) | −4.85e-4 (−311.3) | −2.41e-4 (−347.3) |
| QR7 | +5 | +1.98e-3 (+318.3) | +6.06e-4 (+389.1) | +3.01e-4 (+434.1) |
| TT7 | −9 | −3.56e-3 (−572.9) | −1.09e-3 (−700.3) | −5.42e-4 (−781.5) |

- The ratios to the κ = 1 column are 4.00, 5.00 and 9.00 at every N.
- C3 = Γ(K1) gives exactly zero.

## Heuristic mechanism (HEURISTIC)

H(Γ^r T) is a weighted sum over Motzkin-type walks of the value recursion. Upward steps carry weight j, flat steps
j², downward steps −j, and the walk ends at a point m where P_T is evaluated.
- The typical endpoint m grows with N, so lower coefficients of P_T are suppressed relative to the leading one by
  about m⁻².
- The log factors in (c) would come from the spread of the endpoint distribution.
- A saddle-point analysis of the walk weights should prove (a)–(c).

## Tests

- Extend to N = 4001 with a C implementation of the value recursion, to fix the rate in (c).
- Second order: the residual should be governed by sw_2 = 2c5 + 4·(number of pairs of disjoint 3-cycles), against
  the factor U_(n−4,r)/U_(n,r).
- Prove (a) directly from the walk representation.
