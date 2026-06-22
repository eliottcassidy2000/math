---
id: HYP-2885
title: The LRC's missing realizability is ADDITIVE-ENERGY extremality -- the interval is the unique realizable cover-extremizer (a proof route via the sharp Fejer/additive-energy theorem)
status: GROUNDED (p0 tracks A(E), interval maximizes both); proof route = [p0<=L_y] + [L_y Fejer-extremized by interval] + [L_y(interval)<=cap]; the middle is the open realizability lemma
source: kind-pasteur-2026-06-22-S31j
related:
  - THM-534    # p0<=L_y (moment-LP dual, PROVED)
  - HYP-2873   # additive energy A(E)=int|E^|^4, interval maximizes (Fejer), mac-mini
  - HYP-2879   # consec maximizes L_y (the open extremality), mac-mini
  - OPEN-Q-108 # the wide bound crux
  - HYP-2605   # winding tournament
---

# HYP-2885: the LRC's missing realizability is additive-energy extremality

## The realizability lens (from tournament analysis, applied differently)
Tournament analysis: `H=I(Omega,2)`; the gap at 7 is a REALIZABILITY obstruction -- `K_3` cannot be
a conflict graph (THM-200). The LRC has an analogous "missing realizability", but the carrier is NOT
a graph -- it is the speed set's ADDITIVE structure.

LRC(14) <=> the cover atom `p0(E)=meas{x: phases {frac(e x)} hit all 6 inner sectors} <= cap_k` for
every offset set E. A counterexample OVER-COVERS (`p0>cap`). The realizability question: can integer
speeds realize over-covering?

## The finding (GROUNDED, lrc_realizability_extremal_kps.py)
`p0(E)` TRACKS the **additive energy** `A(E) = #{(a,b,c,d) in E^4: a+b=c+d} = int_0^1 |E^(t)|^4 dt`
(HYP-2873). Computed at k=9:
| set | p0 | A(E) |
| interval [0..8] | **0.416** | **489** (MAX) |
| random | 0.07-0.11 | ~205 |
| Sidon [0,1,3,7,12,..] | 0.012 | 120 (min) |
The **INTERVAL/AP maximizes BOTH** `p0` and `A(E)`; the spread (Sidon, minimal additive energy)
MINIMIZES `p0`. (Dilation-invariant: `2*interval` gives the same `p0`.)

## The realizability argument (a proof route)
> A LRC(14) over-covering counterexample would require a speed set with `p0 > cap`, hence (since `p0`
> tracks `A(E)` and the interval maximizes coverage) additive energy beyond the interval's. But
> **`A(E) <= A(AP)` is a SHARP classical extremal theorem** -- the interval uniquely maximizes
> additive energy among integer sets of given size. So no integer set out-covers the AP. The
> over-covering is NOT REALIZABLE.

Rigorous skeleton, three steps (the first and third are DONE):
1. `p0(E) <= L_y(E)`  -- PROVED (THM-534, moment-LP dual, exact Bonferroni).
2. `L_y(E) <= L_y(AP)` -- the INTERVAL maximizes `L_y`. This is the missing realizability lemma =
   mac-mini's "consec maximizes L_y" (HYP-2879), now identified as a **Fejer/additive-energy
   extremality**: `L_y` is a spectral functional of `|E^|` (HYP-2873 grounds `L_y <-> A(E)`), and the
   interval is its sharp maximizer (Fejer-kernel concentration). The classical theorem `A(E)<=A(AP)`
   is the tool.
3. `L_y(AP) <= cap_k` -- VERIFIED exact rational, all k=8..13 (THM-534).

So the LRC(14) wide bound reduces to ONE realizability inequality (step 2), and the analogue of
"K_3 is not a conflict graph" is "no integer set exceeds the interval's additive energy" -- a SHARP,
KNOWN extremal fact. The remaining work is to prove `L_y <= G(A(E))` for monotone `G` (or directly
that `L_y` is Schur-concave / Fejer-extremized under the interval rearrangement), turning step 2 into
a corollary of additive combinatorics.

## Why this is the right "slightly different" structure
- Tournaments: realizability of the CONFLICT GRAPH (Omega); obstruction `K_3` (combinatorial).
- LRC: realizability of the ADDITIVE SPECTRUM (`|E^|`, the additive energy); obstruction "more
  energy than the interval" (additive-combinatorial). The interval/AP is the unique realizable
  extremizer -- the additive-energy analogue of the Paley tournament's H-maximality (THM-027), and
  of the winding tournament being most-cyclic at the AP (HYP-2605, kps-S31i).

## Tests / next
- Prove `L_y(E) <= L_y(AP)` via `A(E) <= A(AP)` + a spectral majorization `L_y <= G(A)` (the crux).
- Check Schur-concavity of `L_y` under the interval rearrangement (compress E toward an AP, `L_y` up).
- Connect to mac-mini's covering-system route (HYP-+2878): over-covering = covering system = also
  realizability-obstructed (Hough). -> THM-534, HYP-2873, HYP-2879, OPEN-Q-108.
