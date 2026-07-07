---
source: boxeph-2026-07-07-S2 (Lean formalization)
status: FORMALIZED, kernel-pure (no sorry, no native_decide, axioms = [propext,
  Classical.choice, Quot.sound]). Three new Lean modules bridging opus's Farey roof
  (THM-637) to the tail-diameter good set, and an unconditional AP30 density-floor
  certificate that extends death-star's AP20 (diameter <=19 -> <=29).
tags:
  - lonely-runner
  - LRC14
  - lean
  - formalization
  - farey-roof
  - density-floor
---

# The Farey roof bridges to the good set, and derives the diameter floor

Owner (remote): formalize the LRC(14) proof, pull often, integrate incoming work,
best possible state, "healthy diet from the other agents." This session took two
fresh fleet Lean modules -- opus-S135's `LRCFareyRoof` (the pointwise roof theorems
`zero_gap_empty`/`lemmaC`, THM-637) and mac-mini-S42/monad-S2's `LRCTailDiameter`
(the good-set measure `muGood`, `good_anti`, `good_translate`, and the diameter
chain reducing `μ_θ(E)` to `μ_θ(AP_{D+1})`) -- and supplied the missing LINK between
them, plus an unconditional certificate it unlocks.

## What was missing

`LRCTailDiameter`'s whole route was GREEN *conditional on* `AP76Certificate`
(`muGood (1/7) {0..75} ≥ 2314528732/40290957525`), and `LRCFareyRoof` computed the
roof pointwise -- but nothing connected `roof(x) > θ` to membership in the
measure-theoretic good set `Good θ E = {x | ∃ a, ∀ e∈E, θ < fract(e·x − a)}`. The
certificate could not be discharged because the roof lived in one file and the
measure in another, unbridged.

## The bridge (`LRCFareyRoofBridge`, kernel-pure)

- **`good_of_roof_gt`**: on the open Farey-`k` cell `(p/q, p'/q')`, if
  `roof := (qx−p)+(p'−q'x) > θ`, then `x ∈ Good θ (AP_k)`. Witness: put the θ-arc's
  left end at `a := (q'x−p') + (roof−θ)/2` -- the midpoint offset that seats the
  closed θ-arc STRICTLY inside the empty roof-interval `(q'x−p', qx−p)`. Any orbit
  point in `(a, a+θ]` would land strictly inside that interval, contradicting
  `zero_gap_empty`. (Needs no sign hypothesis on θ.)
- **`muGood_ge_of_cell_interval`**: an in-cell interval `(c,d) ⊆ [0,1]` with roof
  `> θ` throughout contributes its length to `μ_θ(AP_k)` (via `good_of_roof_gt` +
  `measure_mono` + `Real.volume_Ioo`).
- **`muGood_ge_sum_intervals`**: a finite pairwise-disjoint family of such intervals
  contributes the SUM of its lengths (via `measure_biUnion_finset`).

Together these reduce `AP76Certificate` from an analytic measure statement to a
**decidable sum of Farey-cell interval lengths** -- and only cells adjacent to a
`q ≤ 6` node contribute (roof node values are `1/q`, so `roof>1/7` needs
`min(q,q') ≤ 6`), so the sum is over ~24 cells, not the ~1776 of Farey-76.

## The payoff (`LRCAP30Floor`, unconditional)

death-star-S2 proved `μ_{1/7}(AP₂₀) ≥ m_P` by exhibiting TWO empty arcs BY HAND.
Those arcs are *exactly* the roof-superlevel intervals of the two cells adjacent to
`0/1` and `1/1`: for `AP_k`, near-0 cell `(0/1,1/k)` has roof `1−(k−1)x`, giving
`(0, 6/(7(k−1)))`, and the mirror near 1 gives `((7k−13)/(7(k−1)), 1)`, each of
length `6/(7(k−1))`. This session DERIVES both from `good_of_roof_gt` and pushes `k`:

> **`ap30_certificate` (kernel-pure): `μ_{1/7}(AP₃₀) ≥ 6/203 + 6/203 = 12/203 ≈
> 0.0591 ≥ m_P = 14249/252252`.**  `ap30_certificate_icc0` rephrases it on `{0,…,29}`
> via `good_translate`, ready for `TailDiameter.muGood_ge_APD` at diameter `≤ 29`.

Two endpoint intervals clear `m_P` for every `k ≤ 31` (`12/(7(k−1)) ≥ m_P`), so the
*systematic* roof route reaches **diameter ≤ 29** with the same two intervals
death-star used by hand for `≤ 19` -- the extra reach is free once the intervals
come from the roof. Adding the `q=2,3,…` node intervals pushes the floor further
toward the tight `AP₇₆` end (where manual enumeration is infeasible but the roof is
not).

## Why this is the right division of labor

- opus-S135 = the roof pointwise engine (`zero_gap_empty`).
- mac-mini/monad = the good-set measure + diameter chain (`muGood`, `muGood_ge_APD`).
- death-star = the loose small-`D` end by hand (`AP₂₀`).
- **this session = the roof→measure bridge + the systematic derivation**, which is
  what the TIGHT `AP₇₆` end needs (its ~24 intervals cannot be hand-built one arc at
  a time, but each is a one-line `good_of_roof_gt` on its Farey cell).

## Ledger

- **Kernel-pure, built, axiom-audited:** `good_of_roof_gt`,
  `roof_superlevel_subset_good`, `muGood_ge_of_cell_interval`,
  `muGood_ge_sum_intervals`, `ap30_certificate`, `ap30_certificate_icc0` -- all
  `[propext, Classical.choice, Quot.sound]`.
- **Open (handed to the fleet):** the tight `AP₇₆` certificate = a ~24-interval
  instantiation of `muGood_ge_sum_intervals` with the `q ≤ 6` Farey-76 cell data
  (each membership a `good_of_roof_gt` call; the sum is decidable arithmetic). Also:
  whether the post-far-peel single-scale residual has diameter `≤ 29` (then `AP₃₀`
  already closes the `k=13` `hlarge` leg) or needs the full `≤ 75`.
- **Files:** `04-computation/lean/TournamentH7/TournamentH7/LRCFareyRoofBridge.lean`,
  `LRCAP30Floor.lean` (+ root import).
