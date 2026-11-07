---
id: HYP-9490
title: "Gate towers forget their core, except through the cyclic-triangle count: for tournaments T, T' of the same order n, H(Gamma^r T)/H(Gamma^r T') = 1 + 2(c3(T) - c3(T')) rho_(n,r) (1 + o(1)) as r -> infinity, where rho_(n,r) = H(Gamma^r TT_(n-2))/H(Gamma^r TT_n) -> 0 is universal; equivalently the universal tower weights U_(d,r) = (L^r q^d)(1) of the gate's cover-polynomial operator satisfy U_(d-2,r)/U_(d,r) -> 0, with rho depending asymptotically only on the total order N = n + 2r, of size about (log N)^2/N^2"
status: >
  OPEN (CONJECTURE (a) and the rate (c)) with an exact reduction; everything else follows from (a).
  PROVED (THM-4900(6)):
  - H(Gamma^r T) = sum_j sw_j(T) U_(n-2j,r), with sw_1 = 2 c3;
  - hence H(Gamma^r T) - H(Gamma^r T') = 2(c3(T) - c3(T')) U_(n-2,r) + sum_(j>=2) (sw_j(T) - sw_j(T')) U_(n-2j,r) exactly.
  So the conjecture is equivalent to an analytic statement about one explicit linear operator L
  (L f(q) = q[f(q+1) + q f(q) - f(q-1)], Session210), with no tournaments left in it.
  PROVED here (Lemmas A-C below):
  - the diffusion form H(Gamma^r T) = [P_T(d/dy) M^r e^y](0), with M = d^2/dy^2 + 2 sinh(y) d/dy the adjoint of L;
  - monotonicity of the weights;
  - (a) implies (b) and the Consequence.
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
  - 04-computation/experiments/microcosm_gate_resolvent_20261106/tower_weights_diffusion.py (+ .out): Lemmas A-C exactly
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

## Proved lemmas (this session)

**Lemma A (the tower is a diffusion; PROVED).** Let

    M = ∂² + 2 sinh(y)·∂ = e^(−2 cosh y)·∂·e^(2 cosh y)·∂    (∂ = d/dy),

and pair polynomials with entire functions by ⟨P, Φ⟩ = [P(∂)Φ](0). Then ⟨LP, Φ⟩ = ⟨P, MΦ⟩ and ⟨P, e^y⟩ = P(1). Hence

    H(Γ^r T) = ⟨L^r P_T, e^y⟩ = [P_T(∂)·M^r e^y](0),    U_(d,r) = (d/dy)^d [M^r e^y](0).

- The macro part of the tower is the r-th power of one Sturm–Liouville operator applied to e^y. The core enters only
  as the constant-coefficient operator P_T(∂), read at the origin.
- On exponentials, M e^(ky) = k e^((k+1)y) + k² e^(ky) − k e^((k−1)y). This is the transpose of Session210's value
  recursion, so M^r e^y = Σ_k W_r(k) e^(ky) with the signed walk weights of the heuristic below.
- r = 1 gives M e^y = e^y + e^(2y) − 1, hence H(Γ T) = H(T) + P_T(2) − P_T(0), which is Session210 (5).

*Proof.* Take P = q^e. Since (q+1)^e − (q−1)^e = 2Σ_(i odd) C(e,i) q^(e−i),

    L q^e = q^(e+2) + Σ_(k≥0) 2C(e, 2k+1) q^(e−2k),

so ⟨L q^e, Φ⟩ = Φ^(e+2)(0) + Σ_k 2C(e,2k+1) Φ^(e−2k)(0). By Leibniz,
∂^e(Φ″ + 2 sinh(y)Φ′)(0) = Φ^(e+2)(0) + 2Σ_i C(e,i) sinh^(i)(0) Φ^(e−i+1)(0), and sinh^(i)(0) = [i odd]. The two
agree. Finally ∂^e e^y = e^y, which is 1 at 0. ∎

**Lemma B (monotonicity; PROVED).** Every U_(e,r) is a positive integer, and:
1. U_(e,r) ≤ U_(e+2,r): at fixed depth, a bigger transitive core carries more paths.
2. U_(e,r) ≥ U_(e+2,r−1): at fixed total order N = e + 2r, a deeper tower carries more paths.
3. 1 ≤ U_(d,r)/U_(d+2,r−1) ≤ 1 + 2^d·ρ_(d+2,r−1).

*Proof.*
- The displayed L q^e has nonnegative coefficients. Splitting off the first step gives
  U_(e,r) = U_(e+2,r−1) + Σ_k 2C(e,2k+1) U_(e−2k,r−1), with nonnegative terms. This gives 2.
- For 1, induct on r. If row r−1 is nondecreasing in e, each term of the recursion at e+2 dominates the matching term
  at e, because C(e+2,·) ≥ C(e,·). The extra terms at e+2 are nonnegative.
- For 3, Σ_k 2C(d,2k+1) = 2^d, and by 1 every U_(d−2k,r−1) ≤ U_(d,r−1) = ρ_(d+2,r−1)·U_(d+2,r−1). ∎

**Corollary C ((a) ⇒ (b); PROVED).** If (a) holds, then for each fixed d both U_(d,r)/U_(d+2,r−1) → 1 and
ρ_(d,r)/ρ_(d+2,r−1) → 1.
- The first limit is Lemma B.3 together with (a) at d + 2.
- For the second, write ρ_(d,r)/ρ_(d+2,r−1) = [U_(d−2,r)/U_(d,r−1)]·[U_(d+2,r−1)/U_(d,r)]. Both factors tend to 1 by
  the first limit, at d − 2 and at d. For d = 2 the first factor is exactly 1, since L·1 = q².

So HYP-9490 is now the single analytic statement (a), plus the rate (c).

**Consequence**, PROVED from (a) and THM-4900(6). Let T, T′ have the same order n. Then

    H(Γ^r T)/H(Γ^r T′) = 1 + 2(c3(T) − c3(T′))·ρ_(n,r)·(1 + o(1)).

- So deep gate towers forget their core. To first order the only trace left is the number of cyclic triangles, i.e.
  the metagraph level of klein's THM-4889.
- With (b), which follows from (a) by Corollary C, the same holds against the one-vertex tower of the same total
  order:

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

**The rate (c) is not settled by the data** (`tower_weights_diffusion.out`, `cover_tower_values.out`). The normalised
value N²ρ/(log N)² is still drifting downward:

| N | 23 | 43 | 83 | 123 | 163 | 203 | 243 | 401 | 801 | 1201 |
|---|---|---|---|---|---|---|---|---|---|---|
| N²ρ/(log N)² | 1.149 | 1.030 | 0.962 | 0.936 | 0.921 | 0.911 | 0.903 | 0.887 | 0.870 | 0.863 |

- For N ≤ 243 the value is ρ_(3,r) = U_(1,r)/U_(3,r) itself; d = 4 agrees to within 1%.
- For N ≥ 401 it is |dev|/2 from the TT3-against-K1-tower column above. If (a) holds, this is the same quantity
  asymptotically (Corollary C).
- A Θ-law with a log-log correction is consistent with these values, and so is a slow decay to 0.

## Heuristic mechanism (HEURISTIC)

H(Γ^r T) is a weighted sum over Motzkin-type walks of the value recursion. Upward steps carry weight j, flat steps
j², downward steps −j, and the walk ends at a point m where P_T is evaluated.
- The typical endpoint m grows with N, so lower coefficients of P_T are suppressed relative to the leading one by
  about m⁻².
- By Lemma A this is literally U_(d,r) = Σ_k W_r(k)·k^d with M^r e^y = Σ_k W_r(k)e^(ky). If the signed weights W_r
  concentrate near k ≈ K_r, then ρ ≈ 1/K_r².
- Weighing flat steps (k²) against up steps (k) over r ≈ N/2 steps puts the height at K_r ≈ N/log N, up to log-log
  corrections. That would give (c).
- The difficulty is the sign of the down steps: W_r is not a positive measure, so concentration needs a real argument
  (a saddle point for the diffusion M, or a positive reformulation).

## Tests

- Extend to N = 4001 with a C implementation of the value recursion, to fix the rate in (c).
- Second order: the residual should be governed by sw_2 = 2c5 + 4·(number of pairs of disjoint 3-cycles), against
  the factor U_(n−4,r)/U_(n,r).
- Prove (a) directly from the walk representation, or from the diffusion M = e^(−2 cosh y)∂e^(2 cosh y)∂ of Lemma A.
  (a) says the Taylor coefficients of M^r e^y at 0 satisfy c_(d−2)/c_d → 0, i.e. M^r e^y oscillates or grows ever
  faster near the origin.
