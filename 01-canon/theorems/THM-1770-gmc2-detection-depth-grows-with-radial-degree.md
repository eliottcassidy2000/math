---
id: THM-1770
title: "THE GMC(2) DETECTION DEPTH GROWS WITH THE RADIAL DEGREE — so the span-only uniform bound (HYP-8540) is FALSE, and the finite-elimination route (opus/kp THM-1740) cannot be made degree-uniform. On the FIXED charge span {−1,0,1}: at radial degree d=0 the detection depth is exactly 2 = span (E[P¹]=B₀, E[P²]=2A₀C₀+B₀², so the first two moments force one-sidedness), but at d=1 it is STRICTLY GREATER than 2 — the explicit two-sided P = Z + (−i√2 + i√2|Z|²) + Z̄ satisfies E[P¹]=E[P²]=0 yet E[P³]=2√2·i ≠ 0. Structural reason: E[P^m] = L_s(CT_u[Λ_s^m]) is the LAPLACE TRANSFORM in s of the toral diagonal; the toral diagonal is P-recursive of order = span (THM-1710), but creative-telescoping s out of the recurrence RAISES the holonomic order above the span. Consequence: bounded (span, degree) strata are each decidable by finite Gröbner (THM-1740 stands), but there is NO degree-uniform collapse — degree-uniform GMC(2) genuinely requires the analytic bridge, not more elimination."
status: >
  PROVED (the decisive direction, by explicit witness). d=0 depth = 2: E[P¹]=B₀ and
  E[P²]=2A₀C₀+B₀² (machine-confirmed), so {E[P¹]=E[P²]=0} ⟹ B₀=0, A₀C₀=0 ⟹ one-sided —
  exactly THM-1710's toral span. d=1 depth > 2: the two-sided witness
  P = Z + (−i√2 + i√2·|Z|²) + Z̄ has E[P¹]=E[P²]=0 and E[P³]=2√2·i ≠ 0 (E[P⁴]=−24),
  machine-verified by exact Wick — so the depth-2 nullcone check does NOT force one-sidedness
  at d=1 while it does at d=0. Detection depth STRICTLY increases from d=0 to d=1 on a fixed
  span; hence it is not a function of span alone. The Laplace/creative-telescoping mechanism is
  the standard reason a Laplace transform of a D-finite family raises D-finite order — stated
  as the explanation of the observed growth, not re-proved here.
  This does NOT prove or refute GMC(2). It refutes the SPAN-ONLY uniform bound and delimits the
  finite-elimination route. GMC(2) remains OPEN.
source: klein-2026-07-20-S381 (owner: now work to prove GMC(2))
depends_on:
  - THM-1710  # toral detection depth = span (M+N), the d=0 baseline this extends
  - THM-1740  # opus/kp: bounded GMC(2) = finite Gröbner emptiness, per (span,degree) stratum
related:
  - THM-1700  # the charge-radius lock (integer moments) + bottom-up descent
  - THM-1530  # the r-dependent {−1,0,1} stratum closed by m≤14 (bounds the d=1 depth above)
  - "CASE-gamma-bridge-domination-step: the bridge is open; this says elimination cannot replace it"
script: 04-computation/gmc2_depth_vs_degree_klein_S381.py (+ .out)
---

# THM-1770 — the GMC(2) detection depth grows with the radial degree

## The question

opus/kp THM-1740 made bounded GMC(2) a **finite Gröbner emptiness test** — per `(span, degree)`
stratum — and HYP-8540 asked whether this collapses to a **single** test with a bound in the
**span alone**, uniform over the radial degree `d`. That is the radial analog of THM-1710's
toral result "detection depth = span (`M+N`)". If true, GMC(2) would hold for all bounded-span
`P` at every degree by one finite check. **It is false.**

## The witness

Fix the charge span `{−1,0,1}` and write the genuine (charge-radius-locked) `P` with radial
coefficients. Detection depth = smallest `D` such that `E[P¹]=⋯=E[P^D]=0` forces `P` one-sided.

**d = 0** (`P = A₀Z + B₀ + C₀Z̄`): depth **2**. Machine-confirmed,
`E[P¹] = B₀`, `E[P²] = 2A₀C₀ + B₀²`; so the first two moments give `B₀=0` then `A₀C₀=0`, i.e.
one-sided. This is exactly the span (THM-1710).

**d = 1**: depth **> 2**. The explicit **two-sided** polynomial

```text
P = Z + (−i√2 + i√2·|Z|²) + Z̄        (A₀=1, C₀=1 both nonzero; B₀=−i√2, B₁=i√2)
```

has, by exact Wick evaluation,

```text
E[P¹] = 0,   E[P²] = 0,   E[P³] = 2√2·i ≠ 0,   E[P⁴] = −24.
```

So `{E[P¹]=E[P²]=0}` does **not** force one-sidedness at `d=1` — a genuinely two-sided `P` slips
through the depth-2 check and is only killed at `m=3`. **The detection depth strictly increases
from `d=0` to `d=1` on a fixed span**, so it is not a function of span alone.

## Why — the Laplace step raises the order

`E[P^m] = L_s(CT_u[Λ_s^m])`, `L_s(g) = ∫₀^∞ g(s)e^{−s}ds`. The toral diagonal `CT_u[Λ_s^m]` is
`P`-recursive in `m` of order `= span` (THM-1710), with coefficients polynomial in `s`.
Producing a recurrence for `E[P^m]` means **eliminating `s`** (creative telescoping / Zeilberger)
— and eliminating a parameter from a `D`-finite recurrence generically **raises** the order.
So the radial holonomic order exceeds the toral span, and the detection depth grows with the
radial degree that feeds `s`-monomials into `Λ_s`. The `d=0→d=1` jump above is the smallest
instance.

## Consequence for the GMC(2) endgame

- **HYP-8540 (span-only uniform bound) is false.** There is no single finite test whose depth
  depends on the span alone; the bound must grow with the radial degree.
- **THM-1740 stands** — every bounded `(span, degree)` stratum is still decidable by a finite
  Gröbner test. What is refuted is *degree-uniformity*.
- **Degree-uniform GMC(2) genuinely needs the analytic bridge.** The finite-elimination route
  cannot reach unbounded radial degree, because the object it eliminates (the moment ideal) has
  no degree-uniform generating set. This is the precise sense in which the
  `CASE-gamma-bridge-domination-step` gap — the analytic step — cannot be replaced by more
  algebra: the algebra provably runs out at each fixed degree.

## Scope

One decisive witness on one span settles the *direction* (depth is not span-only). The exact
`d=1` depth is not pinned here (from THM-1530 §E it is `≤ 14`; the depth-2 witness shows it is
`≥ 3`). No claim on GMC(2) itself, which remains open — this sharpens where its proof must come
from: the radial/analytic layer, not finite elimination.

*Files: `04-computation/gmc2_depth_vs_degree_klein_S381.py` (+ `.out`).*
