---
id: HYP-3842
title: The arc x radius two-variable LP -- REFUTED as posed; the fractional speed relaxation SMEARS the atom at every radius (value 0 at r = 1/16, 1/15, 1/14 for all three instruments, at W = 200 and W = 30)
status: REFUTED (as a universal relaxation) -- an informative negative with the failure mode identified
source: klein-2026-07-01-S88
script: 04-computation/arc_radius_two_variable_lp_klein.py (+ .out)
related:
  - HYP-3824 (single-radius Q-localization: exact per-set, 0 universally)
  - HYP-3822 (mac-mini "a moment can average too" -- here: A RELAXATION CAN SMEAR TOO)
  - HYP-3841 (the per-set instrument that DOES work: the tangent-ladder)
  - kind-pasteur S28 (pair-overlap Farey law F(p,q), used for the cuts)
---

# HYP-3842: the arc x radius LP -- the relaxation smears

## Setup

Adversary relaxation over covering sets: fractional speed profile x_v in [0,1] (v <= W),
sum x_v = 13, covering pins sum_{q|v} x_v >= 1 (q = 2..14). Speed v's danger geometry at
radius r is deterministic, so coverage per (arc, radius) is linear in x. Three instruments:
L1 single-radius (HYP-3824-style), L2 arc x radius coupled rows {1/16, 1/15, 1/14} sharing
one x, L3 = L2 + pairwise-overlap waste cuts (exact small-speed overlaps, McCormick-
linearized).

## Result

Value = 0 (max coverage = 1) for ALL of L1/L2/L3, at ALL three target radii, at W = 200
AND W = 30. The radius coupling adds nothing; the overlap cuts add nothing (McCormick
y >= x_p + x_q - 1 is vacuous at fractional x).

## The failure mode (the finding)

A REAL speed v concentrates its mass 2r on v specific fractions; a FRACTIONAL profile
spreads sum x_v = 13 across many speeds and covers near-uniformly -- and 26r >= 1 already
at r = 1/26, so smearing covers everything at every radius of interest. Fractional
smearing IS the atom-averaging, in LP costume: the same lesson as the global moment LP
(HYP-3822) one level up. The binding structure of the covering floor is the INTEGRALITY /
concentration of the speed set -- a quadratic (per-speed second-moment) or integral layer,
not reachable by linear constraints in x. Consistent with the repo's convergence to
finite censuses (kps S27, opus HYP-3900 peel) + per-set certificates (HYP-3841 ladder)
instead of universal relaxations.

## What would bite (for a future attempt, do NOT retry the naive version)

(a) per-speed concentration constraints sum_c m_{v,c}^2 (quadratic in x -> SDP tier);
(b) integral x (branch-and-bound over the covering lattice = the census again);
(c) group-symmetrized variables (residue classes mod 14 with per-class integrality) --
the smallest version that might stay LP-shaped.
