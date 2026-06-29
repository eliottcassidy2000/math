---
id: HYP-3549
title: Across ALL analogous forms of the LRC (Diophantine, covering-combs, view-obstruction torus-geodesic, billiard, lattice/Littlewood-sibling, Fourier, tournament), a disproof candidate is universally a RATIONAL/closed-geodesic object (equidistribution kills irrationals: central cube vol (1-2/n)^{n-1} -> e^-2 > 0), and the boundary's SHAPE is governed by the multiplicative group (Z/n)^*: the extremal lonely set = the phi(n) UNITS mod n in phi(n)/2 antipodal pairs (= the saddle index), touching at multiplicative inverses; covering = ADDITIVE pinning to walls at Farey points 1/q (q<=n); so a disproof must be simultaneously additive (covering) and anti-multiplicative (cover the units) -- the additive/multiplicative tension. Hyperoctahedral B_{n-1} symmetry with the global R-antipodal x->1-x = the complement (klein THM-584). LRC(2p) inherits saddle index phi(2p)/2=(p-1)/2 from apex p (geometric 2-adic fold).
status: VERIFIED patterns n=4..14 (lonely=units, saddle=phi/2, non-covering witness 683/683, vol->e^-2, 2p->p inheritance). The (Z/n)^* meta-pattern and additive/multiplicative framing are SYNTHESIS.
source: mac-mini-2026-06-29-S16
related:
  - HYP-3548   # the two razor-thin lines / disproof structure (this gives its geometric shape)
  - HYP-3547   # apex-7 Mersenne/Heegner (the multiplicative richness of (Z/14)^*=(Z/7)^*)
  - THM-580    # 2-adic descent = the geometric fold T^{2p-1}->T^{p-1}; saddle inheritance
  - THM-581    # saddle index = phi(n)/2; Borsuk-Ulam witness on the antipodal unit pairs
  - THM-584    # complement = antipodal map (= central-cube x->1-x); hyperoctahedral theme
external: Cusick view-obstruction; Littlewood conjecture (simultaneous-approx sibling); Weyl equidistribution
results:
  - 04-computation/lrc_disproof_geometry_across_forms_macmini_20260629.py
  - 05-knowledge/results/lrc_disproof_geometry_across_forms_macmini_20260629.out
---

# HYP-3549 -- the disproof candidate across the LRC's forms; the multiplicative boundary

## The disproof candidate in each analogous form
| form | space | a disproof candidate is... |
|---|---|---|
| Diophantine | speeds `v in Z^{n-1}` | a primitive covering set with `M(S) < 1/n` |
| covering combs | `t in [0,1)` | danger combs (each meas `1/n`... width `2/n`) cover `[0,1)`, no gap |
| **view-obstruction** | torus `T^{n-1}` | the closed geodesic `{t v mod 1}` AVOIDS the central cube `C_n=[1/n,1-1/n]^{n-1}` |
| **billiard** | box `[0,1/2]^{n-1}` | an unfolded trajectory trapped in the boundary slab (within `1/n` of a wall) |
| **lattice / Littlewood-sibling** | `P^{n-2}(Q)` | a rational direction staying within `1/n` of the coordinate-hyperplane arrangement `forall t` |
| Fourier | `SPEC` | `SPEC = -1` (danger-comb correlations saturate) |
| tournament | covering structure | a covering structure with empty lonely set |

## The cross-form patterns (verified n=4..14)
- **[P1] Equidistribution dichotomy.** `vol(C_n)=(1-2/n)^{n-1}` increases to `e^-2=0.1353 > 0`. So a direction with `Q`-independent ratios EQUIDISTRIBUTES (Weyl) and VISITS `C_n` => not a disproof. **Disproof candidates are RATIONAL/integer directions only** (closed geodesics). The Littlewood conjecture is the simultaneous-approximation sibling: a disproof is a direction "badly distributed toward the walls."
- **[P3] Universal extremal = AP/staircase `{1,...,n-1}`.** Its lonely set is EXACTLY the `phi(n)` **units mod n** (proof: `t=a/n` is lonely iff `gcd(a,n)=1`; if `g=gcd>1` then runner `i=n/g` hits `0`). These form `phi(n)/2` **antipodal pairs** `{a, n-a}` = the **saddle index**; runner `a^{-1} mod n` does the touching. (n=14: 6 units `{1,3,5,9,11,13}`, 3 pairs, saddle 3.)
- **[P4] Covering = ADDITIVE pinning at Farey points.** Non-covering `S` (omits all multiples of some `q<=n`) is killed by the `t=1/q` witness (verified 683/683): the geodesic, at the Farey time `1/q`, has every runner `>= 1/q >= 1/n` from a wall. So a disproof geodesic must pass through a wall at EVERY `1/2,1/3,...,1/n` (covering) AND stay in the wall-tube between (`M<1/n`): **pinned + trapped.**
- **[P2/P6] Hyperoctahedral symmetry + the R-antipodal.** `Sym(C_n)=B_{n-1}=(Z/2)^{n-1} >| S_{n-1}` (sign flips `v_i->-v_i` x coordinate permutations). The GLOBAL reflection `x->1-x` (all coords) `= t->-t =` the **complement/reversal R**; the central-cube antipodal map IS the tournament complement (klein THM-584). Disproof candidates are `B_{n-1}`-orbits of rational directions, R-antipodally paired -- the same signed-permutation theme as the metagraph signed cycle index (THM-587).

## The meta-pattern: (Z/n)^* shapes the boundary; additive vs multiplicative
In every form the boundary's SHAPE is the **multiplicative group `(Z/n)^*`**: lonely set = the units, antipodal pairs = `(Z/n)^*/{+-1}` (saddle `= phi(n)/2`), touch = inversion `a -> a^{-1}`. The covering constraint, dually, is **additive** (pinning at the `q`-divisions `1/q`). **A disproof must be at once additive (covering -- hit every `1/q`) and anti-multiplicative (cover the units `a/n`)** -- and these pull against each other (the S15 "covering forces structure" tension, now revealed as the additive/multiplicative tension of `n`). Where `(Z/n)^*` is multiplicatively rich (n=2p, p Mersenne/Heegner -- HYP-3547), the boundary is most rigid. **`(Z/14)^* = (Z/7)^*`**: LRC(14)'s multiplicative boundary IS the apex-7 cyclotomic group.

## The 2p -> p geometric fold
LRC(2p) inherits its boundary shape from LRC(p): `phi(2p)/2 = phi(p)/2 = (p-1)/2` (verified p=3,5,7,11,13), so LRC(14) and LRC(7) have the SAME 3 antipodal pairs. This is the geometric face of the 2-adic descent (THM-580): folding the even coordinates projects the `T^{2p-1}` geodesic to a `T^{p-1}` geodesic, carrying the units `(Z/2p)^* -> (Z/p)^*`. The descent IS a torus fold.

## What it buys
A unified geometric picture of where a disproof can live (rational directions, pinned at Farey points, trapped in the wall-tube, organized by `(Z/n)^*` and the R-antipodal) -- and a sharpened "why hard": the disproof must reconcile the additive covering pins with the multiplicative unit-avoidance, hardest exactly when `(Z/n)^*` is rich (the apex prime). Ties the geometric LRC to the tournament metagraph (shared `B`, shared R) and to Littlewood (the approximation sibling).
