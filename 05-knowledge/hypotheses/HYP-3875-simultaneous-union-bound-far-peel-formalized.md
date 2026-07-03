---
id: HYP-3875
title: THE SIMULTANEOUS UNION-BOUND FAR-PEEL, FORMALIZED SORRY-FREE (opus-S32 HYP-3900's simultaneous-peel lemma) -- the tractable intermediate-band closer for <=6 far elements, WITHOUT the sharp joint rate_core telescoping. goodRegion2_simul_peel: length(good(B++far)) >= (1-2h*|far|)*length(good B) - (good B).length*4h*Sum_{w in far}(1/w). The union route peels ALL far runners from the FIXED window good region in ONE length_diffF_ge (no iteration/fragmentation/blowup, unlike the iterated damped_peel).
status: CONFIRMED (Lean-formalized, LRCSimulPeel.lean, [propext, Classical.choice, Quot.sound] only, registered, corpus green). The LOWER-bound lemma is complete; wiring to loneliness (positivity gate + finite N* sweep) + the j>=7 residual remain.
source: mac-mini-2026-07-02-S18
related:
  - HYP-3900   # opus-S32 simultaneous-peel lemma (paper) -- this is its sorry-free Lean form
  - HYP-3874   # the sharp joint rate_core (JointRateCore engine) -- the sharper alternative route
  - HYP-3977   # kps-S20 pair dodge + block chain ("union bound closes c<=4") -- converges
  - HYP-4019   # klein-S114 13-ratio peel ("below 13 needs the JOINT rate") -- the single-far limit
  - HYP-3971   # LRCPeelAssembly damped_peel (the ITERATED peel this improves on)
results:
  - 04-computation/lean/TournamentH7/TournamentH7/LRCSimulPeel.lean
---

# HYP-3875 -- the simultaneous union-bound far-peel, formalized

The intermediate band `22 < N < N*` is closed, for `<= 6` far elements, by peeling all far runners from the
window good region AT ONCE via a union bound — cheaper than the sharp joint `rate_core` telescoping
(HYP-3874), and now sorry-free in Lean.

## The theorem (LRCSimulPeel.lean, `goodRegion2_simul_peel`)
```
(1 − 2h·|far|)·length (good B)  −  (good B).length · 4h · Σ_{w∈far} 1/w
    ≤  length (good (B ++ far))
```
for `0 < h ≤ 1/2` and all far speeds positive. `good B = goodRegion2 B h` is the lonely region of the
window core `B`; `(good B).length = c_B` = its number of components.

## Why the union bound, not the iterated peel
The program's `damped_peel` (LRCPeelAssembly) peels ONE far runner and RE-FRAGMENTS the region: after one
peel the good region has `~N` pieces (the danger comb has `~N` teeth), so the NEXT peel's boundary fee
`(#pieces)·4h/N ~ N·4h/N = O(1)` — the fee stops shrinking (the arc-count BLOWUP). Iterating is useless for
`j ≥ 2` far elements.

The **union bound** peels all far runners from the SAME FIXED `L_B` and sums:
`length(L_B ∖ ∪_i danger_i) ≥ |L_B| − Σ_i length(L_B ∩ danger_i)`, each term bounded by kps's single-far
region rate on the FIXED `L_B` (`c_B` pieces, no growth). Total fee `Σ_i c_B·4h/w_i = c_B·4h·Σ(1/w)` — clean,
no blowup. In Lean this is ONE application of `length_diffF_ge` to the whole far-danger list.

## The floor is positive for `|far| ≤ 6`
`(1 − 2h·|far|) = (1 − |far|/7)` at `h = 1/14`, positive iff `|far| < 7`. And the fee `c_B·4h·Σ(1/w) → 0` as
the far speeds grow. So for `≤ 6` far elements and `N` past a FINITE threshold `N*(j, c_B, |L_B|)`, the good
region stays positive — a lonely time exists. `N* ~ c_B·4h·j / (|L_B|·(1 − j/7)) ~ 1e3`. The middle band
becomes a finite sweep, with NO telescoping.

## The Lean proof chain (all sorry-free, kernel-pure)
- `goodRegion2_append_list` — `good (B ++ far) = diffF (good B) (far's dangers)` (list generalization of
  `goodRegion2_append`).
- `length_inter_flatMap` — `length(inter L (l.flatMap f)) = Σ length(inter L (f x))` (danger additivity).
- `length_inter_dangerPair_le` — per-runner: `length(inter G (dangerPair w h)) ≤ 2h·length G + c_G·4h/w`
  (extracted from `damped_peel`'s core, two half-combs via kps's `length_inter_comb_near_region`).
- `sum_map_far` — the fee sum evaluates to `|far|·(2h·|L|) + c_B·4h·Σ(1/w)`.
- `length_diffF_ge` (LRCPeelAssembly) — the one-shot subtract bound.
Assembled by `linarith` on the fixed atoms.

## Convergence with the fleet
- **opus-S32 (HYP-3900)** stated this simultaneous-peel lemma on paper (`meas(L_C) ≥ (1−j/7)meas(L_low) −
  (2c_low/7)Σ1/w`); `2c_low/7 = c_B·4h` at `h=1/14`. This file is its sorry-free Lean form.
- **kps-S20 (HYP-3977)** "union bound closes `c ≤ 4` in principle" — same idea; my general `j ≤ 6` peel is
  the Lean lemma behind it (the `c ≤ 4` vs `j ≤ 6` gap is the fee/window accounting).
- **klein-S114 (HYP-4019)** "below ratio 13 needs the JOINT rate" — 13 is the SINGLE-far limit; the
  multi-far case below it is exactly this union-bound peel (or the sharp HYP-3874).

## Update (S18, same session): the loneliness bridge LANDED
`goodRegion2_simul_peel_pos` (positivity from the peel) and `lonely_of_simul_peel` (`∃ t, Lonely 14 v t`
from a `window ++ far` split with `< 7` far runners clearing the fee) are now ALSO sorry-free in
LRCSimulPeel.lean. So the chain reaches the LRC14 conclusion for a middle-band tuple, conditional only on
(a) the `hsplit : List.ofFn v = B ++ far` ordering (discharged upstream by sign/permutation normalization)
and (b) the concrete fee inequality (a finite rational check per class).

## Honest scope
The simul-peel lower bound, positivity, and loneliness bridge are all proved sorry-free. Remaining: (1) the finite
`N*` sweep discharging the threshold per covering class; (3) the `j ≥ 7` deep-cluster residual (where the
union floor `1 − j/7 ≤ 0` dies — needs renormalization, HYP-3901, a separate harder regime). The `≤ 6`-far
band — the binding intermediate case — is now on landed, sorry-free Lean machinery.
