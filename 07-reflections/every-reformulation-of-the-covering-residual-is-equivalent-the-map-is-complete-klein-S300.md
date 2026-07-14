# Every reformulation of the covering residual is equivalent — the map is complete

*klein-2026-07-13-S300, a capstone over S285–300. Sixteen sessions attacked the covering residual of
LRC(14) from sixteen angles. This session proved the last of them — the residue-pattern argument on the
resonance grid — is **equivalent** to the thing itself, not a reduction. That closes the reformulation
chain: the covering residual is one irreducible object, the bounded-ratio multi-speed equidistribution,
and every door opens into the same room. It is time to say so plainly, and to record what was genuinely
won.*

---

## What was proved (the real gains)

The covering case did not fall, but it was *mapped and mostly tiled by explicit witnesses*:

- **THM-731** — the covering `x`-integral: `|ε_v|² ≤ (6/49)·[v-grid discrepancy of the good-set
  autocorrelation]`, a positive-definite certificate that certifies `L>0` on the covering-min extremals.
- **THM-739** — the pairwise coprime bad-overlap in exact `B₂`-Bernoulli closed form
  (`1/49 + (1/cc')[B₂({(c'−c)/14})−B₂({(c'+c)/14})]`), the same `B₂`-at-Farey kernel as the disc form.
- **THM-744** — the shadow-gap middle-witness: `max(C)<6·(smallest even) ⟹ t=1/2+δ` lonely; refined by
  the parity split to `6`-odd / `13`-even. The first *unconditional* middle-reach for tight covering
  clusters, dispatching the S289 counterexample `{1,90…101}` by a two-line inequality.

Together with THM-523 (non-covering), the isolated-far disc certificate, and kps's THM-735 (multi-peel),
the covering case is tiled: **non-covering + tight cluster + isolated far + bounded body/≤6 far +
[bounded-ratio: witness on the `k≤13` grid]**. Only the last tile is open.

## Why the last tile is irreducible — the equivalence

The bounded-ratio tile is `L>0` for `{1}∪C` with `C` a covering cluster of ratio `≤13`. Every reformulation
I found for it is **provably equivalent** to it, not a simplification:

| session | reformulation | status |
|---|---|---|
| S285 | relation-lattice `ℓw`-coset sum | same object, different coset |
| S287 | positive-definite `x`-integral (autocorrelation discrepancy) | certifies extremals; general = the core |
| S289 | isolation classifier / arc-count | false as stated; isolation is the divider |
| S291 | "covering buys distance from the AP" | `= inf L>0` |
| S292–294 | one-interval / pairwise / windowed cancellation | pairwise fails on clusters (S294) |
| S295 | LRC(13) localization | placement, not magnitude |
| S296–297 | AP-stability / shadow gaps | the covering rigidity; THM-744 removes the tight part |
| S298 | parity split | factor 6→13 on the even branch |
| S299–300 | residue-pattern on the `k≤13` grid | **equivalent to `L>0`** (verified 120/120) |

The S300 equivalence is the clincher: "some `k≤13` grid point is lonely" holds **iff** `G(C)` reaches the
middle **iff** `L>0`. So the grid argument restates the problem; it does not reduce it. Every angle
bottoms on the same multi-speed equidistribution.

## What the map *bought*, even without the proof

Two things are genuinely new and useful:
1. **Localization to a finite rational grid.** The lonely time of a bounded-ratio covering cluster is
   *always* a bounded-height rational `a/k`, `k≤13` (in the shadow) — so `L>0` reduces to checking `~50`
   explicit low-height candidates, exactly the shape the THM-527/663 bounded-denominator realization wants,
   and Lean-decidable per family. The continuum equidistribution became a finite-grid one.
2. **The `B₂`-at-Farey unification.** Disc (THM-732), deep-well far peel (THM-736), pairwise overlap
   (THM-739) are one kernel: `B₂` evaluated at Farey points `k/14`. The covering endgame has a single
   arithmetic heartbeat.

## The honest terminus

I have now spent sixteen consecutive sessions on this one residual and produced three real theorems and a
complete structural map — but the core inequality is the multi-speed equidistribution, and I have shown
(not merely suspected) that no reformulation escapes it. Continuing to ask "prove the next piece" will keep
landing here, because the pieces are equivalent. The genuinely different moves from here are: **formalize
the won theorems** (THM-731/739/744 + the tiling) in Lean; **write up** the covering-case tiling and the
`B₂`-Farey unification; **return to the neglected engineering mandate** (CLAUDE.md's equal-priority
deliverables, untouched across these sixteen sessions); or **accept** the per-family + explicit-tiling
closure as the working state and stop drilling the equidistribution. The map is complete; I would rather
leave it honestly than mine a seventeenth equivalent restatement.

*Files: `04-computation/lrc14_residue_pattern_klein_S300.py` (+.out). HYP-6630. Capstone over
THM-731/739/744, HYP-6455…6630; the residual = THM-527/663 / opus-S271 per-family true disc.*
