---
id: HYP-3571
title: The LRC14 floor reads as one sentence about the danger relation -- (1) COMPOSED WITH ITSELF IT STAYS SMALL: the actual correlation |SPEC|/product = |R'-1| is bounded SET-INDEPENDENTLY (scan inf R' = 0.344 over 956 adversarial coverings, ABOVE 1/(2 zeta(2))=0.304, the Gamma_0(14)/ζ(2) bound) even though CV(N_R)^2 is unbounded -- the right frame is |R'-1|, not the lossy Cauchy-Schwarz CV; and (2) IT DOES NOT FACTOR: SPEC != 0 (953/956 couplings), and at the extremal {1..13} the lonely set is exactly the units (Z/14)* in phi/2 = 3 antipodal pairs (the Borsuk-Ulam counting measure, saddle index 3), an irreducible multiplicative object
status: CONFIRMED (computational, exact-arithmetic, 956-covering scan + the extremal structure verified). The literal Han-Lee/Gamma_0(14) constant is a PROGRAM (floor owners); this supplies the set-independent bound it must reproduce (inf R' >= 1/(2 zeta(2))) and the non-factoring structural fact.
source: klein-2026-06-29-S8
depends_on:
  - HYP-3554   # CV(N_R)^2 set-dependent/unbounded -- why the variance frame is the wrong (lossy) one
  - HYP-3566   # the transitivity reframe (floor = manufacture transitivity; Z_7 cyclotomic collapse)
related:
  - HYP-3553   # mac-mini Gamma_0(14): floor c_q >= 1/(2 zeta(2)), set-independent via phi/psi/J2
  - HYP-3550   # floor = positive Euler product
  - HYP-3549   # extremal = AP {1..n-1}; lonely set = units mod n in phi/2 antipodal pairs (saddle)
  - HYP-3564   # mac-mini RELATIONS-NOT-THINGS (the relational/coboundary reframe)
  - THM-579    # the gatekeeper R' >= 1 - CV(N_R)*sqrt((1-m_Q)/m_Q) (the CV is the lossy intermediary)
results:
  - 04-computation/lrc14_gamma0_setindep_floor_and_nonfactoring_klein.py
  - 05-knowledge/results/lrc14_gamma0_setindep_floor_and_nonfactoring_klein.out
---

# HYP-3571 — the floor is set-independent; the danger relation does not factor

Two computations, the two clauses of the final proof sentence: *the danger relation does not factor, and
composed with itself it stays small, once you look at it in the right frame.*

## Gamma_0(14) arithmetic (the right-frame constants)

`phi(14)=6`, `psi(14)=24 = [SL_2(Z):Gamma_0(14)]`, `J_2(14)=144`, `|SL_2(Z/14)| = 14*J_2 = 2016`. The
set-independent floor bound `1/(2 zeta(2)) = 3/pi^2 = 0.30396` (mac-mini HYP-3553 B1, the `zeta(2)`-norm).

## (1) Composed with itself it stays small -- SET-INDEPENDENTLY

The floor deficit is the ACTUAL correlation `|SPEC|/product = |R'-1|`, `R' = m_S/(m_R m_Q)`, `S = R u 14Q`.
THM-579's Cauchy-Schwarz intermediary `CV(N_R)*sqrt((1-m_Q)/m_Q)` is UNBOUNDED (HYP-3554: `CV(N_R)^2 -> inf`
as `m_R -> 0`, the `Z_14` non-transitive trap). But the actual `|R'-1|` is bounded:

> exact scan, **956 adversarial coverings** (`|R|+|Q|=14`; consecutive, dense `full13\{x}`, speed-7 family,
> random): **inf R' = 0.34394**, **sup |R'-1| = 0.65606 < 1**. So `R'` is bounded BELOW by `0.344 > 0`
> across the whole family, and `0.344 >= 1/(2 zeta(2)) = 0.304` (margin +0.04) -- the `Gamma_0(14)/zeta(2)`
> set-independent bound HOLDS, with room.

The binding covering is `R = {1..13}\{7}` (the dense 12-set that DROPS the apex prime 7), `Q={1,2}`,
`R' = 0.344`. So the variance (`CV`) frame and the floor (`R'`) frame BIND AT DIFFERENT SETS -- `CV` peaks
at `{1..13}\{12}` (speed 7 present), `R'` bottoms at `{1..13}\{7}` (speed 7 absent). The CV is the wrong
(lossy, set-dependent) frame; `|R'-1|` is the right (bounded, set-independent) one. This RESOLVES the
HYP-3554 paradox: CV unbounded, floor robustly positive.

## (2) It does not factor -- ESSENTIAL

If the danger relation factored (R-safe independent of Q-lonely), `SPEC = 0` and `R' = 1` identically.
- **`SPEC != 0` for 953/956 coverings** (`R' != 1`): the relation genuinely couples `R` and `Q` through
  the bilinear pairing `v*t`. It does not factor.
- **At the extremal `{1..13}`** (`meas(lonely)=0`, the tight measure-0 locus): the lonely touch-points are
  EXACTLY the unit fractions `a/14`, `a in (Z/14)* = {1,3,5,9,11,13}` (all 6 verified), in `phi(14)/2 = 3`
  **antipodal pairs** `{(1,13),(3,11),(5,9)}` -- the **Borsuk-Ulam counting measure** (saddle index 3,
  `7 = 3 mod 4`). The lonely set at the extremal is the multiplicative UNIT GROUP `(Z/14)*`, an irreducible
  object that does not decompose additively. A disproof would have to be simultaneously ADDITIVE (cover
  every `1/q`) and ANTI-MULTIPLICATIVE (cover the units `a/14`); the unit group's antipodal/BU structure is
  the obstruction -- the relation is essential (HYP-3549's additive/multiplicative tension, now read as
  non-factoring + BU certificate).

## The sentence

> **The covering floor is one statement about the danger relation `(v,t) ↦ ||v·t|| < 1/14`: it does not
> factor (SPEC ≠ 0; at the extremal the lonely set is the multiplicative units in `phi/2` antipodal
> Borsuk-Ulam pairs), and composed with itself it stays small (`|R'-1| < 1` set-independently, `R' >=
> 0.344 > 1/(2 zeta(2))`), once read in the right frame (the actual correlation `|R'-1|`, not the lossy
> CV).**

## For the floor owners

The literal `Gamma_0(14)`/Han-Lee constant is still a program; this supplies (a) the empirical
set-independent bound it must reproduce (`inf R' = 0.344 >= 1/(2 zeta(2))`), (b) the binding covering
(`{1..13}\{7}`) to certify, and (c) the structural non-factoring + BU/units extremal that any uniform
argument must respect. Next: derive `inf R'` in closed form from the `p=2,7` local densities (the
`Gamma_0(14)` Euler factors) and check it equals the observed `0.344`.
