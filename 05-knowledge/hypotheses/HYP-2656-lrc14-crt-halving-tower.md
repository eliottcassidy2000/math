---
id: HYP-2656
title: LRC14 CRT/halving — odd base is constant, even speeds = dyadic doublings, minimizer CRT-factors
status: CONFIRMED (exact, single-deletion layer)
source: kind-pasteur-2026-06-19
depends_on:
  - HYP-2651
  - THM-541
  - THM-523
related:
  - HYP-2660
  - HYP-2659
  - HYP-2658
  - HYP-2646
  - HYP-2644
  - HYP-2655
  - OPEN-Q-108
---

# HYP-2656 — LRC14 CRT / Halving Structure

## Claim

The `n=14 = 2·7` factorization, with gap `1/14 = (1/2)·(1/7)`, gives a literal
2-adic/7-adic split of the lonely measure. Three exact facts, all on the
single-deletion core layer `[1,13]\{e}`:

1. **The odd sub-core is a constant base.** For every even `e`, the odd speeds of
   `[1,13]\{e}` are the SAME set `{1,3,5,7,9,11,13}`, with constant lonely measure
   ```
   L_odd = 75454/315315 = 0.2392972107.
   ```
   The odd speeds carry NONE of the discrimination between cores; all variation is
   in the even/doubling content (Glaisher: every even speed `= 2^a · (odd)`).

2. **The dyadic levels HALVE the lonely measure** (the "2" in `14 = 2·7`).
   Adding the 2-adic levels of the even speeds to the odd base, by level
   `level_j = {v even : v_2(v) = j}`:
   ```
   odd base                       L = 0.23929721
   + level 1 {2,6,10}  (= 2·odd)  L = 0.11729223   ratio 0.490 ≈ 1/2
   + level 2 {4,12}    (= 4·odd)  L = 0.04519290   ratio 0.385
   + level 3 {8}       (= 8·odd)  L = 0.00000000
   ```
   Level 1 halves the measure almost exactly. The doubling map `z ↦ 2z` (order 3
   mod 7, splitting sectors into QR/NQR) is the same map that halves band widths.

3. **The drop-6 minimizer CRT-FACTORS exactly.** Using the cell decomposition
   `τ = a/7 + ξ/7`, `a = 0..6` (a 7-adic × archimedean split):
   ```
   L(drop-6) = (#surviving cells / 7) · R = (2/7) · (49/1716) = 7/858.
   ```
   The 7-adic factor `2/7` (only cells `a=1,6` survive, a τ↔1−τ mirror pair) is
   pure residue-mod-7 combinatorics; the archimedean factor `R = 49/1716` is the
   within-cell density, equal on both surviving cells by the mirror.

## Why drop-6 wins, through the halving lens

The 4 critical free intervals of the drop-6 minimizer are bounded by exactly the
three largest speeds `{11, 12, 13}`, and the boundary-owner parity is balanced:
**4 even-owner boundaries (all owned by 12 = 2·6) and 4 odd-owner boundaries
(11, 13).** Deleting `6 = 2·3` is the cheapest single deletion precisely because
its dyadic-tower parent `12 = 2²·3 = 2·6` survives and re-covers `6`'s danger
bands. The halving structure is the mechanism: the deepest member of odd-3's
tower `{6,12}` stays in.

The full single-deletion ranking (exact):
```
drop  6 (2·3 ) L = 7/858     = 0.00815851   <- global min, tower {6,12} keeps 12
drop 12 (4·3 ) L = 426/35035 = 0.01215927
drop 10 (2·5 ) L = 1520/63063= 0.02410288
drop  4 (4·1 ) L = 97/4004   = 0.02422577
drop  2 (2·1 ) L = 11/364    = 0.03021978
drop 13 (odd ) L = ...       = 0.03410122   <- best ODD deletion is far worse
...
drop  7 (odd ) L = ...       = 0.08389666   <- worst (removes the mult-of-7 runner)
```
Codex S37 caveat: the phrase "five even deletions" is the low-even block
`{2,4,6,10,12}`.  The remaining even deletion, drop `8=2^3`, is the high-level
even outlier with `L=950/21021=0.04519290`, after the odd drops `13,5,3` in
the fixed-observer safe-measure order.  Thus the dyadic tower explains the low
collar rows, but tower parity alone is not a total scalar order; the CRT cell
and endpoint-owner layers still matter.

## CRT factorization is special, not generic

Checking the clean-product property `R_a` constant across all surviving cells for
every single-deletion core: it holds ONLY for the low-surviving-cell cores
(drops 1,2,3,4,6, with #surv ∈ {1,2}) and FAILS for the high-L cores (drops
5,7,8,9,10,11,12,13 with #surv ∈ {4,6}). The minimizer drop-6 is the unique core
that minimizes BOTH factors simultaneously: minimal 7-adic factor (2/7) AND
minimal archimedean within-cell density among the clean cores. Minimization
happens on both CRT factors at once.

## Evidence

Scripts (all exact `fractions.Fraction`, target `1/14`):
- `04-computation/lrc14_crt_halving_oddcontrol_kps.py` — odd base constant + halving cuts
- `04-computation/lrc14_odd_base_dyadic_tower_kps.py` — exact `L_odd`, dyadic cascade
- `04-computation/lrc14_drop6_evenowner_kps.py` — boundary-owner parity, tower deletion
- `04-computation/lrc14_crt_factor_test_kps.py` — exact CRT product for drop-6
- `04-computation/lrc14_crt_factor_general_kps.py` — factorization generality scan

Outputs in `05-knowledge/results/` (same basenames, `.out`).

## Relationship to existing proof state

This is a re-reading of the proved single-deletion atlas (HYP-2651 / THM-541)
through the `14 = 2·7` lens, and it extends the prior `lrc14_7adic_archimedean_
split_s6.py` (which found the 2-cell support for stranger-cores) to the EXACT
single-deletion minimizer, showing the clean `(2/7)·R` product holds there too.

It does NOT close OPEN-Q-108 (the uniform fattening lemma over ALL primitive
cores). It clarifies architecture: on the single-deletion layer the lonely
measure is `L_odd`-rigid-base + dyadic-halving-refinement, and minimization is a
joint 7-adic/2-adic optimization. The open multi-element / far-element layer
(HYP-2644, HYP-2655) is where the rigid-base picture must be extended.
