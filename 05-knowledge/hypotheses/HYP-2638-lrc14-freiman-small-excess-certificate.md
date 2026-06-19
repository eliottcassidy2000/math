---
id: HYP-2638
title: LRC(14) Freiman small-excess certificate - finite 3k-4 pocket for the additive-energy route
status: CLAIMED stub; evidence pending
source: codex-2026-06-19-S27
depends_on:
  - HYP-2637
  - HYP-2635
  - HYP-2607
  - HYP-2604
related:
  - THM-531
  - THM-535
  - OPEN-Q-108
---

# HYP-2638 - LRC(14) Freiman Small-Excess Certificate

## Claim Being Tested

The additive-energy lead from HYP-2635/HYP-2637 has a sharp first finite
pocket: ordinary sumset excess

```text
exc(E) = |E+E| - (2|E|-1)
```

should be enough to close the low-dimensional side up to the Freiman `3k-4`
threshold.  If `|E+E| <= 3k-4`, equivalently `exc(E) <= k-3`, Freiman's
`3k-4` theorem places `E` inside an arithmetic progression of length at most
`|E+E|-k+1 = k+exc(E)`.  Since the LRC sector functional is invariant under
translation and dilation of the offset set, this turns the infinite
small-excess pocket into a finite normalized enumeration.

The intended certificate is:

1. Enumerate normalized `k`-subsets of `[0, k+e-1]` containing `0`, grouped by
   exact excess `e <= k-3`.
2. Compute exact `L_y(E)` and the sector distribution for each row.
3. Verify that `e=0` is precisely the AP row and every `e>=1` row has a
   strict margin below the AP/cap frontier for `k=8,9,10` (the dangerous
   LRC14 small-cluster sizes).
4. Use this as the finite theorem-shaped pocket inside the broader
   relation-fiber/GAP split: low excess is finite; high excess must either
   peel an uncovered relation vertex or enter a higher-rank GAP/tail estimate.

## Missing Evidence

This file is an early namespace reservation.  The next commit should add the
script, stored output, and exact maxima by excess.  If the finite table does
not preserve a clean AP margin, this hypothesis should be downgraded to a
negative route note rather than promoted.

## Assumption Challenge

The first challenged assumption is that the useful tournament vertices are raw
runners or ordinary pair sums.  For this pocket the vertices should be proof
obligations:

```text
Freiman_3k4_pocket
> exact_excess_layer
> AP_invariance
> relation_fiber_cover
> higher_rank_GAP_tail
> raw_pair_sum_energy
> raw_runner_vertices
```

This quotient preserves the LRC predicate only after recording the exact sector
functional on each finite normal form; it destroys the individual multiplicand
sign and reciprocal-tail data, which must remain in HYP-2636/HYP-2634.
