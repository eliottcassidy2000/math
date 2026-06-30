---
id: HYP-3551
title: Where 7/89 comes from (binding-pair / packing-radius anatomy of M for covering sets) AND a TIGHTER covering set: {1..12,182} has M=14/183 (+7.1%), superseding the recorded tightest 7/89 (+10.1%); the tightest covering set = the densest coverable core {1..12} (M=1/13) + the minimal killer 182=lcm(13,14) that perturbs it least; through the Littlewood lens the LRC is ANTI-Littlewood -- a POSITIVE floor on the simultaneous-approximation product that Littlewood drives to 0
status: VERIFIED (exact pair-scan + 4M-point grid): M({1..11,13,84})=7/89, M({1..12,182})=14/183<7/89; structured search supports 14/183 as the tightest. The Littlewood/anti-Littlewood framing is SYNTHESIS. Does NOT threaten LRC (14/183 > 1/14).
source: mac-mini-2026-06-29-S18
related:
  - HYP-3548   # the two razor-thin lines; refines its GAP-line margin from +10.1% to +7.1%
  - HYP-3550   # the resonance Euler product (the floor); this is the GAP-side companion
  - THM-523    # covering reduction; covering = kill all b<=14
  - HYP-2566   # 'inf M over covering sets > 1/14' (uniform looseness) -- this supports it
external: Littlewood conjecture (simultaneous Diophantine approximation); three-gap; continued fractions
results:
  - 04-computation/seven_over_89_littlewood_macmini_20260629.py
  - 05-knowledge/results/seven_over_89_littlewood_macmini_20260629.out
---

# HYP-3551 -- where 7/89 comes from, and the tighter 14/183

## Where 7/89 comes from (the anatomy)
`M({1..11,13,84}) = 7/89` at `t* = 37/89`, **binding pair (5, 84)**: since `84 == -5 (mod 89)`, the two
runners 5 and 84 sit at the SAME distance `7/89` from the wall at `t*`, so `89 = 84 + 5` is the
**binding modulus**. The numerator `7` is the **packing (covering) radius**: at the optimal rotation
`a=37`, the 13 speeds `*37 mod 89` all avoid the `7`-neighbourhood of `0` (sorted distances
`7,7,8,14,15,22,23,29,30,36,37,38,44`), and `7` is the largest such radius. General law (S15):
`M(S) = j/D` with `D` a **binding-pair sum/difference** and `j` the **covering radius of the speeds in
Z/D`. The set is the covering completion of the **skip-12 core** `{1..11,13}` (`M=1/12`) by the minimal
killer `84 = lcm(12,14)` (kills resonances 12 and 14). The single-large family `{1..11,13,84k}` gives
`M = 7k/(84k+5) -> 1/12`, minimized at `k=1` -> `7/89`. (`89 = F_11` is a coincidence: the `k=2`
sibling is `14/173`, not Fibonacci.)

## A TIGHTER covering set: 14/183 supersedes 7/89
`{1..12, 182}` is primitive covering (`182 = 13*14` kills resonances 13 and 14) with
> `M({1..12,182}) = 14/183 ≈ 0.07650`  (**+7.1%** above `1/14`),
GRID-VERIFIED (4M points + golden refine) -- **tighter** than the previously-recorded tightest
`7/89 ≈ 0.07865` (+10.1%). Binding pair `(1, 182)` (`182 == -1 mod 183`), modulus `183 = 182+1`, radius
`14`. So the recorded gap-line margin of HYP-3548 sharpens from `+10.1%` to `+7.1%`.

## Why the densest core wins (the structural understanding)
The tightest covering set is the **densest coverable core + the minimal killer that perturbs it least**.
Among single-skip cores, `{1..12}` (skip 13) has `M=1/13`, **denser** (closer to `1/14`) than `{1..11,13}`
(skip 12, `M=1/12`); and its required killer `182=lcm(13,14)` is the LARGEST minimal killer (13,14
coprime), so it equidistributes and barely punctures the `1/13` structure -> `14/183` just below `1/13`.
Verified across cores: skip-13 -> 14/183 (tightest), skip-12 -> 7/89, skip-11 -> 14/157, skip-9 ->
14/131, skip-10/8 -> 1/11, 1/9. The covering infimum is bounded below by the densest-coverable-core
value (`~1/13`), well above `1/14` -- **supporting HYP-2566's uniform looseness** (covering sets stay a
fixed margin above `1/14`).

## The Littlewood lens: the LRC is anti-Littlewood
At the bind, `||v_a t*|| = ||v_b t*|| = M`, a **simultaneous approximation** of one wall by both runners;
the Littlewood product `||v_a t*||*||v_b t*|| = M^2`, and `D*M^2 = j^2/D` (e.g. `89*(7/89)^2 = 49/89`).
The **Littlewood conjecture** says `inf_q q||q alpha||||q beta|| = 0` -- the product CAN be driven to 0.
The **LRC floor `M >= 1/14`** is exactly the opposite: a POSITIVE lower bound `M^2 >= 1/14^2` on that
product at the binding -- **the LRC is the anti-Littlewood**, the obstruction to simultaneous
wall-approximation collapsing. The tightest covering set is where this obstruction is WEAKEST (the LRC's
closest approach to a Littlewood vanishing), pinned at `14/183`. NOTE (honest): the "badly-approximable
golden-ratio" intuition is a RED HERRING here -- the binding ratios `L/partner` have LARGE partial
quotients (`84/5=[16,1,4]`, not bounded), so the tightness comes from **core density + minimal
perturbation**, not from bad approximability. The Littlewood connection is conceptual (anti-Littlewood
positive floor), pairing with the S17 result that the floor itself is a positive Euler product.

## What it buys
Corrects the recorded extremal (14/183 < 7/89), sharpens the gap-line margin to +7.1%, and gives the
structural law for the tightest covering set (densest coverable core + minimal far killer), which both
explains the margin and supports the uniform-looseness conjecture (HYP-2566). The anti-Littlewood
framing places the LRC opposite Littlewood: where Littlewood's product vanishes, the LRC's stays
positive -- the same multiplicative-positivity theme as the Euler-product floor (HYP-3550).
