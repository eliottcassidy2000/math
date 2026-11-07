---
id: HYP-9490
title: "The forgetting rate of gate towers: rho_(d,r) = H(Gamma^r TT_(d-2))/H(Gamma^r TT_d) = Theta((log N)^2/N^2) as r -> infinity, N = d + 2r. Deep gate towers forget their core except through the cyclic triangles, H(Gamma^r T)/H(Gamma^r T') = 1 + 2(c3(T) - c3(T')) rho_(n,r) + O(rho_(n,r)/r), and this is PROVED (THM-4900(7)) with rho <= 1/(r+1); the conjecture is the sharp rate"
status: >
  OPEN (CONJECTURE: the rate only). The first version (2026-11-06) also conjectured (a) rho -> 0, (b) universality in
  N, and the first-order law. These are now PROVED in THM-4900(7):
  - rho_(d,r) <= 1/(r+1) uniformly in d;
  - the weights and rho depend asymptotically only on N;
  - H(Gamma^r T)/H(Gamma^r T') = 1 + 2 Delta c3 rho_(n,r) + E, with |E| <= C(T, T') rho_(n,r)/(r + 1 - sigma(T')).
  Audit H1 corrected the first version's error form 2 Delta c3 rho (1 + o(1)) to 2 Delta c3 rho + o(rho): the old
  form fails when Delta c3 = 0 but a higher sw_j differs (MISTAKE-729).
  FINITE-EXACT evidence for the rate: rho_(3,r) exactly to N = 243, and its exact transform to N = 1201. The
  normalised value N^2 rho/(log N)^2 is still drifting (1.149 at N = 23 down to 0.864 at N = 1201).
session: mac-mini-2026-11-06-microcosm
source: 05-knowledge/results/microcosm_gate_doubling_resolvent_20261106.md
scripts:
  - 04-computation/experiments/microcosm_gate_resolvent_20261106/c3_first_order_law.py (+ .out)
  - 04-computation/experiments/microcosm_gate_resolvent_20261106/cover_tower_values.py (+ .out)
  - 04-computation/experiments/microcosm_gate_resolvent_20261106/tower_weights_diffusion.py (+ .out)
  - 04-computation/experiments/microcosm_gate_resolvent_20261106/tower_forgetting_bound.py (+ .out)
related: [THM-4900 (statements 6-7), Session210 graded-spectrum-gates (operator, Theorem 4/5), Session159 cycle_profiles (1) and profile_operator (1)-(2) (operator, all-order odd-cycle form of P_T, sw_2 = 2c5 + 4 pairs), Session159 lagrange (5), (7) (transitive-core towers as exponential sums), Session150 transfer, THM-078 (odd-cycle form, odd n), klein THM-4889/HYP-9481 (levels by cyclic triangles), Session210 structure atlas (q = (H-1)/2)]
---

# HYP-9490 — the rate at which gate towers forget their core

## Setting

- L f(q) = q[f(q+1) + q f(q) − f(q−1)] is the gate's action on the ordered path-cover polynomial (Session210;
  first in Session150/159).
- U_(d,r) = (L^r q^d)(1) = H(Γ^r TT_d), and ρ_(d,r) = U_(d−2,r)/U_(d,r) for d ≥ 2.
- N = d + 2r is the total order of the tower Γ^r TT_d.

## Conjecture

For each fixed d ≥ 2, as r → ∞,

    ρ_(d,r) = Θ((log N)²/N²).

Weaker forms, each still open:
- ρ_(d,r) = o(1/N). This would improve the proved bound ρ ≤ 1/(r+1) ≈ 2/N.
- ρ_(d,r) = N^(−2+o(1)).
- A lower bound ρ_(d,r) ≥ c/N². None is proved, because the exponential coefficients c_(r,q) below are signed.

## What is proved (THM-4900(7))

- **Operator form.** H(Γ^r T) = [P_T(∂)·M^r e^y](0), with M = ∂² + 2 sinh(y)∂ = e^(−2 cosh y)∂e^(2 cosh y)∂, the
  adjoint of L under ⟨P, Φ⟩ = [P(∂)Φ](0).
  - Equivalently, U_(d,r) = Σ_q c_(r,q) q^d, with Session159's coefficients
    c_(r+1,q) = (q−1)c_(r,q−1) + q²c_(r,q) − (q+1)c_(r,q+1) (lagrange/NOTE.md (7)).
- **Monotonicity and the bound.** U_(e,r) ≥ U_(e+2,r−1) and U_(e,r) ≥ (r+1)·U_(e−2,r). The proof is an induction on r
  with Pascal's rule; the slack term dominates U_(e−2,r) term by term.
- **Forgetting with an explicit error.** For T, T′ of the same order n:

      H(Γ^r T)/H(Γ^r T′) = 1 + 2(c3(T) − c3(T′))·ρ_(n,r) + E,
      |E| ≤ ρ_(n,r)·[σ(T) + σ(T′) + 2|Δc3|·σ(T′)]/(r + 1 − σ(T′)),    σ(X) = Σ_(j≥1)|sw_j(X)|.

- **Universality.** 1 ≤ U_(d,r)/U_(d+2,r−1) ≤ 1 + 2^d/r, and ρ_(d,r)/ρ_(d+2,r−1) → 1. So ρ depends asymptotically
  only on N.

So the conjecture asks only for the true size of ρ. The proved bound 1/(r+1) is linear in 1/N. The data say ρ is
smaller than that by a factor of about N/(log N)².

## Evidence (exact integer arithmetic)

Tower values come from the cover-polynomial value recursion P_Γ(j) = j[P(j+1) + jP(j) − P(j−1)] at integer points.
That recursion reproduces brute-force Hamiltonian-path counts to order 19. The weights come from the first-step
recursion (R) of THM-4900(7).

**The rate** (`tower_weights_diffusion.out`, `cover_tower_values.out`). The normalised value N²ρ/(log N)², at d = 3:

| N | 23 | 43 | 83 | 123 | 163 | 203 | 243 | 401 | 801 | 1201 |
|---|---|---|---|---|---|---|---|---|---|---|
| N²ρ/(log N)² | 1.149 | 1.030 | 0.962 | 0.936 | 0.921 | 0.911 | 0.903 | 0.887 | 0.871 | 0.864 |

- For N ≤ 243 the value is ρ_(3,r) = U_(1,r)/U_(3,r) itself. d = 4 agrees to within 3.0% at N = 24, 1.2% at N = 44,
  and 0.5% from N = 84 on.
- For N ≥ 401 it is |dev|/2 from the TT3-against-K1-tower comparison below. That equals ρ_(3,r)/(1 + 2ρ_(3,r))
  exactly (P_(C3) = q³ + 2q), so it is ρ_(3,r) to relative accuracy 2ρ < 4·10⁻⁴.
- The value is still drifting. A Θ-law with a log-log correction fits, and so does a slow decay to 0.

**The first-order law, as an illustration of THM-4900(7)(c), all cores of a given order, at N ≈ 401**
(`c3_first_order_law.out`):

| Order n | Classes | Reference | Spread of dev/κ | Max \|dev − κ·ε\| |
|---:|---:|---|---|---:|
| 4 | 4 | TT2 tower, c3 = 2 | exact | 1.1e-19 |
| 5 | 12 | K1 tower, c3 = 4 | [3.9595e-4, 3.9619e-4] | 1.65e-7 |
| 6 | 56 | TT2 tower, c3 = 6 | [3.9397e-4, 3.9499e-4] | 5.57e-7 |
| 7 | 456 | K1 tower, c3 = 9 | [3.9429e-4, 3.9754e-4] | 1.74e-6 |

Here κ = c3(T) − c3(reference) and ε = 2ρ. At order 4 the law is exact, since P_T = q⁴ + 2c3·q².

**Large N against the one-vertex tower** (`cover_tower_values.out`). Relative deviations, with N²·dev in brackets:

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

- U_(d,r) = Σ_q c_(r,q)·q^d, and M^r e^y = Σ_q c_(r,q)e^(qy). If the signed weights c_(r,·) concentrate near
  q ≈ K_r, then ρ ≈ 1/K_r².
- In the value-recursion walk, up steps from height q carry weight q, flat steps q², and down steps −q. Weighing flat
  against up steps over r ≈ N/2 steps puts the height at K_r ≈ N/log N, up to log-log corrections. That gives the
  conjectured rate.
- The difficulty is the sign of the down steps. The c_(r,q) are not a positive measure, so concentration needs a real
  argument.
- Two candidates:
  - a saddle point for M. On y = iθ, M becomes −∂_θ² + 2 sin θ·∂_θ, minus the generator of the von Mises diffusion,
    which is self-adjoint and nonnegative on L²(e^(2 cos θ)dθ);
  - a positive reformulation. The degree walk of THM-4900(7)(R) and Session159's level walk both have positive
    weights.
- The proved bound 1/(r+1) compares the slack S with U_(e−2,r) only at ratio ≥ 1 per term. The k = 0 term has ratio 4,
  and the true per-step gain grows with r. A comparison that keeps track of where U_(e−2,r)'s mass sits would be the
  natural next step.

## Tests

- Extend to N = 4001 with a C implementation of the value recursion, to fix the drift.
- Second order: the residual E in THM-4900(7)(c) should be governed by Δsw_2 = 2Δc5 + 4Δ(pairs of disjoint 3-cycles),
  against the factor U_(n−4,r)/U_(n,r).
- Prove ρ = o(1/N) from a two-step induction that also uses the slack of the down jumps.
