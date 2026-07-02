---
id: HYP-3846
title: The sec-7.3 arc x radius LP, n=6 pilot -- with single-anchor data the LP is SEPARABLE (closed form = the arc-localized ladder); the constraint set (monotonicity + coverage + slope transport) is EXACT on the n=6 final window for both tight sets, at Q=6 and Q=31
status: CONFIRMED (pilot; exact rationals, no solver needed at this scope)
source: klein-2026-07-01-S89
script: 04-computation/arc_radius_lp_n6_pilot_klein.py (+ .out)
related:
  - opus S33 doc sec 7.3 (the LP formulation piloted here)
  - HYP-3841 (the ladder), HYP-3842 (the smearing negative -- why coverage-side fractional variables fail)
  - HYP-3834/3835 (c(AP_5) = c({1,3,4,5,9}) = 2/5; final window [1/7, 1/6] exactly linear)
---

# HYP-3846: the arc x radius LP pilot at n=6

## Finding 1 (structure): separability

With data at a single anchor layer r0 (per-arc lonely masses x_{a,0} and per-arc
gap-rate budgets b_a) and constraints (i) monotonicity, (ii) coverage, (iii) slope
transport, the LP min sum_a x_{a,L} decouples per arc:
    LP = sum_a max(0, x_{a,0} - b_a (r_L - r0))
-- the ARC-LOCALIZED LADDER in closed form. A solver only earns its keep with multi-layer
coupling and the universal (covering-axiom-parameterized) version.

## Finding 2 (exactness on the pilot)

n=6, anchor 1/7, target 1/6, grids Q = 6 and Q = 31 (= Phi_6(6), the Gamma_0-aligned
modulus): the LP value equals the true Lambda(1/6) = 0 EXACTLY for both tight sets
(AP_5 {1..5}, anchor 2/35 = the collapse law (2/5)(1-6/7); sporadic {1,3,4,5,9}, anchor
41/630). The n=6 window is defect-free by the band lemma (d=1 band (6d,7d) at d=1 is
(6,7), empty; diam <= 12), so transport budgets at the anchor are the true worst case and
the constraint set (i)-(iii) is exact -- the collapse law is LP-certified per-set at n=6.

## Honest scope + next

Per-set only. The universal version (min over covering configurations with D_{a,l}
constrained by covering axioms alone) remains open -- and HYP-3842's smearing negative
says its COVERAGE side must carry integrality/concentration, not fractional masses; the
lonely-mass-side variables of sec 7.3 with transport constraints are the promising
difference. Next: two-layer coupling (does transport + a second anchor beat the
single-anchor closed form on a set with an interior convex kink?), then the universal
parameterization at n=6 where the census is tiny.
