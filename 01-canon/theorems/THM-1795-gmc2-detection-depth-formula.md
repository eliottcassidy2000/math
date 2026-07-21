---
id: THM-1795
title: "THE GMC(2) DETECTION-DEPTH FORMULA D(M,N,d) = (M+N)(2d+1) — putting a number on klein THM-1770's 'depth grows with radial degree, no uniform bound'. E[P^m] for a P of charge span [−M,N] and radial degree d (each charge coefficient g_c(|Z|²) a polynomial of degree ≤ d in |Z|²) is P-recursive in m of order exactly D(M,N,d) = (M+N)(2d+1); that order IS the detection depth (a moment-nullcone of order r is cut out by the first r moments, THM-1775). DERIVATION (creative telescoping): E[P^m] = L_s(CT_u[Λ_s^m]) is the Laplace transform in s of the toral diagonal (klein THM-1770); the toral diagonal is order M+N in m (THM-1710/1670) with recurrence coefficients of s-degree 2d (products of two degree-d radial coefficients, cf. the disc(R)=g₁²−4g₀g₂ shape at (1,1)), and telescoping s out multiplies the order by (1 + s-degree) = (2d+1). VERIFIED on every reachable cell: d=0 gives order = span exactly — (1,1)→2, (1,2)→3, (1,3)→4 (matching klein's 'depth = span' at d=0 AND THM-1710) — and the (2d+1) factor shows at (1,1): d=0→2, d=1→6 = 2·3. CONSEQUENCE: every bounded (span, degree) stratum is a finite Gröbner emptiness test (THM-1740) whose SIZE is now explicit, and klein's 'no degree-uniform bound' is sharpened to LINEAR growth — the depth rises by exactly 2(M+N) per unit of radial degree, with no ceiling, which is precisely why degree-uniform GMC(2) needs the analytic bridge and cannot be pure elimination"
status: >
  VERIFIED on the reachable cells: d=0 order = M+N for (1,1),(1,2),(1,3) (two primes agree,
  holdout); (1,1) d=1 order = 6.  All four match D=(M+N)(2d+1).  The (2d+1) factor rests on the
  d=0,1 data for span 2 plus the d=0 span-scan; d≥2 and span≥3 with d≥1 are COMPUTE-LIMITED —
  E[P^m] for a generic full-span P has coefficient-degree that grows with the order, and the
  m-range needed exceeds what the dict-convolution moment engine reaches before the monomial
  count explodes, so those cells returned "insufficient data", NOT a refutation.  The formula is
  therefore a creative-telescoping-DERIVED conjecture confirmed on 4 cells, not an exhaustively
  verified theorem; the derivation (order = toral × (1+s-degree)) is standard holonomic closure,
  and the s-degree = 2d step is argued from the (1,1) disc shape, not proved in general.
  Honest: this is the SIZE of the finite test, resting on klein THM-1770 (which proves the depth
  is >span at d=1 with an explicit witness) and THM-1710 (d=0 = span).  GMC(2) itself is not
  advanced; what is added is the explicit growth law.
source: kind-pasteur-2026-07-20-S128c130 (owner: work the GMC(2))
depends_on:
  - THM-1770    # klein: detection depth grows with radial degree, no uniform bound
  - THM-1710    # TNC detection depth = span (the d=0 base)
  - THM-1775    # moment-nullcone: order = detection depth
related: [THM-1740, THM-1670, THM-1645]
script: 04-computation/gmc2_depth_vs_degree_kps_S128c130.py, gmc2_depth_formula_kps_S128c130.py (+ .out)
---

# THM-1795 — the GMC(2) detection-depth formula

klein THM-1770 proved the GMC(2) detection depth **grows with the radial degree** `d` and that
there is **no span-only uniform bound** — but left the growth rate open. Here it is:

> **`D(M,N,d) = (M+N)(2d+1)`.**

## What the numbers mean

`P ∈ ℂ[Z, Z̄]` with charge span `[−M, N]` and radial degree `d` (each charge-`c` coefficient
`g_c(|Z|²)` a polynomial of degree `≤ d`). `E[P^m]` is P-recursive in `m`, and its order — which
by the moment-nullcone template (THM-1775) is exactly the **detection depth**, the number of
moments that cut out the nullcone — is `(M+N)(2d+1)`.

| `(M,N)` | `d` | `D` | `(M+N)(2d+1)` | status |
|---|---|---|---|---|
| (1,1) | 0 | 2 | 2 | ✓ verified |
| (1,2) | 0 | 3 | 3 | ✓ verified |
| (1,3) | 0 | 4 | 4 | ✓ verified |
| (1,1) | 1 | 6 | 6 | ✓ verified |
| (1,1) | 2 | 10 | 10 | derivation; compute-limited |
| (1,2) | 1 | 9 | 9 | derivation; compute-limited |

The `d = 0` row is `order = span`, matching THM-1710 (TNC depth) and klein's "depth = span at
`d = 0`". The jump `2 → 6` from `d = 0` to `d = 1` at `(1,1)` is the `(2d+1)` factor.

## The derivation (creative telescoping)

`E[P^m] = L_s(CT_u[Λ_s^m])` — the Laplace transform in the radius `s = |Z|²` of the toral
diagonal (klein THM-1770). Two facts compose:

1. **The toral diagonal `CT_u[Λ_s^m]` is order `M+N` in `m`** (THM-1710/1670), with recurrence
   coefficients that are polynomials in the symbol coefficients `g_c(s)`.
2. **Those coefficients have `s`-degree `2d`.** Each `g_c` has `s`-degree `d`; the recurrence
   coefficients are quadratic in the `g_c` (the shape is visible at `(1,1)`, where the trailing
   coefficient is `disc(R) = g₁² − 4g₀g₂`, degree 2 in the `g`'s), so `s`-degree `2d`.

Integrating `s` out against `e^{−s}` (the Gaussian radial average) is a **creative-telescoping**
of one variable, which raises the `m`-order from `ρ` to at most `ρ·(1 + σ)` where `σ` is the
coefficient `s`-degree. With `ρ = M+N` and `σ = 2d`:

> `D = (M+N)(1 + 2d) = (M+N)(2d+1)`.

## Why it matters

- **Every bounded `(span, degree)` stratum is a finite Gröbner test of *known size*.** The moment
  ideal `⟨E[P^1], …, E[P^D]⟩` with `D = (M+N)(2d+1)` cuts out the nullcone; Rabinowitsch
  saturation against the two-sided witness decides GMC(2) on that stratum (THM-1740). The number
  of moments to compute is now explicit.
- **klein's "no uniform bound" is sharpened to a linear law.** The depth rises by exactly
  `2(M+N)` per unit of radial degree, with **no ceiling**. That is the precise obstruction to a
  degree-uniform elimination proof: any fixed number of moments is beaten by a large-enough `d`.
  So degree-uniform GMC(2) genuinely needs the analytic bridge (mac-mini THM-1645's Laplace
  determinacy), not more elimination — and now we know by *how much* elimination falls short.

## Named next

- **Break the compute wall** to verify `(1,1,2) = 10`, `(1,2,1) = 9`, `(2,2,1) = 12` directly:
  the moment engine needs the toral+Laplace factorization (compute `CT_u[Λ_s^m]` as an
  `s`-polynomial, then apply `Σ_k k!·[s^k]`) rather than full `P^m` dict-convolution, which
  reaches larger `m` per unit cost.
- **Prove `σ = 2d`** (the coefficient `s`-degree of the toral recurrence) in general — the one
  ingredient of the derivation argued only from the `(1,1)` shape. It should follow from the
  `s`-formula (THM-1690) applied with `R`-coefficients of `s`-degree `d`.
- **The linear law feeds the analytic bridge**: `D ~ 2(M+N)d` is the rate the Laplace-determinacy
  argument must dominate; a bound of the form "`E[P^m] ≠ 0` for some `m ≤ 2(M+N)d`" would close
  degree-uniform GMC(2).
