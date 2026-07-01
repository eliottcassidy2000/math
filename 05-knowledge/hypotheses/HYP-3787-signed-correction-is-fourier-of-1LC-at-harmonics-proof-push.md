---
id: HYP-3787
title: THE SIGNED CORRECTION = FOURIER OF 1_{L_C} AT HARMONICS OF THE TEST SPEED -- a quantitative bound that makes the single-huge-speed and <=6-huge-speed cases of the covering-min residual RIGOROUS. Upgrading S74's equidistribution heuristic: the measure of the core's fixed lonely set L_C(r) covered by a test speed w is covered(w) = integral 1_{L_C} 1_{||wt||<r} = 2r|L_C| (DIAGONAL, the equidistribution main term, j=0) + sum_{j!=0} hat(1_{L_C})(jw) ghat(j) (the SIGNED CORRECTION, off-diagonal = the Fourier coefficients of 1_{L_C} sampled at the harmonics jw of the test speed, ghat(j)=sin(2pi j r)/(pi j)). Since L_C is a union of N arcs, |hat(1_{L_C})(m)| <= N/(pi|m|), hence |signed correction| <= sum_{j!=0}(N/(pi|jw|))(1/(pi|j|)) = (N/(pi^2 w)) 2 zeta(2) = N/(3w). So covered(w) <= 2r|L_C| + N/(3w) < |L_C| once w > N/(3(1-2r)|L_C|) => M(C u {w}) >= r RIGOROUSLY (single huge speed, explicit threshold). Multi-patch: union bound covered(H) <= |H| 2r|L_C| + sum N/(3 w_i) < |L_C| when |H| 2r < 1 (|H|<=6) and each w_i > threshold => M(C u H) >= r RIGOROUS for <=6 huge speeds. Bounded speeds <= threshold: lazy-cut ILP (HYP-3782). RESIDUAL narrows to >=7 huge speeds (needs the cross-harmonic terms hat(1_{L_C})(jw_i - j'w_j), inclusion-exclusion). VERIFIED n=14: the Fourier identity matches the direct covered measure to 5-6 digits; the signed correction decays like 1/w (0.067,0.013,0.0013,-0.00002 for w=7,50,182,10000); |hat|<=N/(pi m) confirmed; thresholds 11,59,259 for cores {1..6},{1..9},{1..11}.
status: FOURIER IDENTITY + DECAY BOUND rigorous/verified (Parseval + the finite-arcs 1/m decay); single-huge-speed and <=6-huge-speed cases RIGOROUS (explicit threshold w>N/(3(1-2r)|L_C|)) modulo the elementary |hat|<=N/(pi m). The >=7-huge-speed case (cross-harmonic) is the remaining residual. Combined with the lazy-cut (bounded), this closes most of the covering-min lower bound and reduces the rest to a named Fourier statement. NOT a complete proof; a substantial proof-push improving S74 (heuristic -> quantitative) and the repo's signed-cut thread.
source: mac-mini-2026-06-30-S75
related:
  - HYP-3786   # S74 the equidistribution reduction (this quantifies its correction term)
  - HYP-3784   # S73 the single-patch three-gap scaling (this is the Fourier version, any huge w not just 182k)
  - HYP-3782   # S72 lazy-cut (the bounded regime; this handles speeds > threshold)
  - HYP-2948   # Beurling-Selberg minorants (ghat is the danger-arc Fourier / Fejer kernel)
  - HYP-3571   # the floor 1/(2 zeta(2)); the zeta(2)=sum 1/j^2 here is the same constant
results:
  - 04-computation/signed_correction_fourier_LC_macmini_20260630.py
  - 05-knowledge/results/signed_correction_fourier_LC_macmini_20260630.out
---

# HYP-3787 -- the signed correction is the Fourier of 1_{L_C} at harmonics of the test speed

The owner's seed -- "signed correction is the off-diagonal version of this: the Fourier coefficients of
`1_{L_C}` sampled at harmonics of the test speed" -- turns S74's equidistribution *heuristic* into a
*quantitative* bound, and that bound closes the single-huge-speed and `<=6`-huge-speed cases of the covering-min
residual rigorously.

## The Fourier identity (diagonal + signed correction)
For a test speed `w`, the measure of the core's fixed lonely set `L_C(r)` that `w`'s danger arcs cover is, by
Parseval (`1_{||wt||<r}(t) = g(wt)`, `hat(g)` supported on multiples of `w`):

    covered(w) = 2r|L_C|                                  <- DIAGONAL (j=0): the equidistribution main term
               + sum_{j != 0} hat(1_{L_C})(jw) * ghat(j)  <- SIGNED CORRECTION (off-diagonal)

with `ghat(j) = sin(2 pi j r)/(pi j)` (the Fourier coefficients of a single danger arc, a Fejer/Dirichlet
kernel). The correction is exactly the **Fourier coefficients of `1_{L_C}` sampled at the harmonics `jw` of the
test speed**. Verified numerically: `covered_direct = 2r|L_C| + corr` to 5-6 digits for all `w`.

## The decay bound (why the correction is small)
`L_C` is a union of `N` arcs, so `|hat(1_{L_C})(m)| <= N/(pi|m|)` (verified: `sup_m pi m |hat(m)| <= N`). Hence

    |signed correction| <= sum_{j!=0} (N/(pi|jw|)) (1/(pi|j|)) = (N/(pi^2 w)) * 2 zeta(2) = N/(3w).

(The `zeta(2) = sum 1/j^2` here is the same Basel constant as the floor `1/(2 zeta(2))`, HYP-3571.) So

    covered(w) <= 2r|L_C| + N/(3w) < |L_C|   as soon as   w > N/(3(1-2r)|L_C|).

## The proof push
- **Single huge speed** `w > N/(3(1-2r)|L_C|)` (explicit threshold; `11, 59, 259` for cores `{1..6},{1..9},
  {1..11}`): `covered(w) < |L_C|`, so a lonely time survives => `M(C u {w}) >= r`. RIGOROUS (Fourier decay of
  `1_{L_C}`). This is the S73 three-gap result for *any* huge `w`, not just multiples of `n(n-1)`.
- **Multi-patch, `<= 6` huge speeds** (each `> threshold`): union bound `covered(H) <= |H| 2r|L_C| + sum N/(3
  w_i) < |L_C|` when `|H| 2r < 1` (`|H| <= 6`). RIGOROUS => `M(C u H) >= r`.
- **Bounded speeds** (`<= threshold`): the lazy-cut ILP (HYP-3782).
- **Residual**: `|H| >= 7` huge speeds -- the union bound fails, and one needs the **cross-harmonic** terms
  `hat(1_{L_C})(jw_i - j'w_j)` (inclusion-exclusion / joint equidistribution). This is the last piece.

## What it improves
- **S74** (equidistribution heuristic `covered ~ 2r|L_C|`) becomes a **quantitative identity with an explicit
  error bound `N/(3w)`** and an explicit survival threshold.
- The repo's **"signed correction / signed cut"** thread gets its explicit form: the off-diagonal correction is
  `sum_{j!=0} hat(1_{L_C})(jw) ghat(j)`, decaying `O(N/w)` by the finite-arcs Fourier decay.
- The residual (S74's "effective Erdos-Turan") is now concrete: **the Erdos-Turan bound IS this Fourier decay**;
  for `>= 7` huge speeds it is the cross-harmonic sum.

## Honest scope
The Fourier identity and the decay bound `|signed correction| <= N/(3w)` are rigorous (Parseval + `|hat| <=
N/(pi m)` from `N` arcs), verified numerically. They make the **single-huge-speed** and **`<=6`-huge-speed**
cases rigorous (explicit thresholds), which with the lazy-cut (bounded speeds) closes most of the covering-min
lower bound. The **`>=7`-huge-speed** case (cross-harmonic inclusion-exclusion) remains the residual. This is a
substantial proof-push, not a complete proof; and it needs, per core `C`, the finite quantities `N(L_C)` and
`|L_C|` (finitely many cores, all computable).
