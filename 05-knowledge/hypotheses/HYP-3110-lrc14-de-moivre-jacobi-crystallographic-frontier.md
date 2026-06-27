---
id: HYP-3110
title: LRC14 De Moivre, Jacobi theta, and crystallographic proof-carrier frontier
status: RESERVED / synthesis and formalization lane; not a proof
source: codex-2026-06-27-S263
tangent: T1186
technique: LTI-247
tournament_technique: LTT-145
related:
  - HYP-3109
  - HYP-3108
  - HYP-3107
  - HYP-3106
  - HYP-3105
  - HYP-3104
  - HYP-3103
  - HYP-3073
  - HYP-3063
  - HYP-2614
  - HYP-2613
  - HYP-2309
  - OPEN-Q-108
---

# HYP-3110: LRC14 De Moivre, Jacobi Theta, And Crystallographic Proof-Carrier Frontier

## Reservation Claim

This lane is reserved to test the user's De Moivre quintic, Jacobi theta,
17-wallpaper-group, and 230-space-group prompt against the live LRC14 proof
frontier.

It is not a proof.  The intended contribution is a controlled-forgetting
interface:

```text
de_moivre_quintic_fold
  -> exact finite-depth quintic cancellation detector;
jacobi_theta_tail
  -> signed residue-cusp / support-six theta tail carrier;
wallpaper_17_orbifold
  -> finite 2D crystallographic quotient audit;
space_group_230_orbifold
  -> finite 3D crystallographic quotient audit;
lee_yang_root_curve
  -> HYP-3109 zero-locus / zero-real ear sidecar;
finite_address_packet
  -> HYP-3107 terminal proof interface.
```

The preserved LRC predicate is: a residual row must retain enough coverage,
observer, endpoint-owner, residue-cusp, and finite-address data to imply
`LRC14Statement` through the HYP-3107 frontier.  The destroyed coordinates are
raw runner labels, raw time, and any scalarized root/moment/lattice count.

## Immediate Work

1. Build a dependency-free scout that verifies the De Moivre quintic identity
   as a finite Laurent-polynomial cancellation, records Jacobi-theta and
   crystallographic sidecar fields, and runs Tournament Analysis on the proof
   carriers.
2. Add a Lean-facing ledger with exact counts for the 17 wallpaper groups and
   230 three-dimensional space groups, plus the De Moivre quintic normal-form
   identity as a checked algebraic lemma if the local Lean algebra automation
   accepts it cleanly.
3. Integrate the result with HYP-3108/HYP-3109 rather than duplicating their
   Lee-Yang root-curve work.

## Assumption Challenge

Do not assume the useful vertices are runners, arcs, roots, or group names.
The candidate vertices are proof-carrier sidecars: quintic-fold normal forms,
theta tails, wallpaper orbifold quotients, space-group orbifold quotients,
root-curve packets, observer-gluing certificates, and finite-address exits.

