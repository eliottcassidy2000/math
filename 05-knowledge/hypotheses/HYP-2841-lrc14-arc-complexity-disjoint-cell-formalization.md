---
id: HYP-2841
title: LRC14 arc-complexity disjoint-cell formalization -- the THM-546/HYP-2840 factor-six saving is now a sorry-free finite Lean lemma
status: FORMALIZED SUPPORT LEMMA; analytic B_j(E') and breakpoint-cell instantiation remain open
source: codex-2026-06-22-S87
related:
  - THM-546
  - HYP-2839
  - HYP-2840
  - HYP-2838
---

# HYP-2841 -- disjoint-cell arc-complexity formalization

## Summary

Incoming work created two distinct HYP-2840 files: mac-mini's L_y/arc-complexity
route and kind-pasteur's Vitali-covering route. I left both intact and reserved
HYP-2841 for the formal support lemma common to the THM-546 sharpening:

> If the exactly-one-miss regions `B_j` are pairwise disjoint and all live in a
> common finite breakpoint-cell partition, then their total occupied cell count
> is bounded by the size of that partition. If the partition has at most
> `7 * sumE` cells, the old `6 * 7 * sumE` ledger collapses to `7 * sumE`.

This is the finite Lean core of the claimed factor-six saving in THM-546: the
six sector ledgers should not be summed independently once exactly-one-miss
cells are known to be disjoint.

Concurrent integration: kind-pasteur S31 pushed `LRCMarginalUniform.lean` while
this note was being prepared. That module proves the marginal one-speed sector
measure atom for HYP-2840; this HYP-2841 module proves the adjacent finite
disjoint-cell counting atom. They are complementary pieces of the same
decorrelation/Vitali/far-peel pipeline.

S87b addendum: `LRCMarginalUniform.lean` is now root-imported and `Verify`-
audited next to `LRCArcComplexity.lean`; its unused-hypothesis warnings were
removed by marking the unused interval-side assumptions intentionally unused.
Combined support-atom transcript:
`05-knowledge/results/lrc14_hyp2840_support_atoms_lean_codex_s87b.out`.
After KPS S31c, the cleanest interpretation is that these atoms support the
far-element/decorrelation half of "consec maximizes `L_y`" in the THM-534 route:
`LRCMarginalUniform` supplies one-speed sector measure, while `LRCArcComplexity`
keeps the exactly-one-miss variation ledger at `7*sumE` rather than `42*sumE`.

## Lean status

New module:

`04-computation/lean/TournamentH7/TournamentH7/LRCArcComplexity.lean`

The module defines a finite occupied-cell count
`occupiedCount I B := sum i in I, (B i).card` and proves:

- `occupiedCount_eq_card_biUnion`: pairwise disjoint finite cell families have
  occupied count equal to the cardinality of their union.
- `occupiedCount_le_cells`: if every occupied cell belongs to a common finite
  partition `C`, then the occupied count is at most `C.card`.
- `occupiedCount_le_seven_mul_sum`: if additionally `C.card <= 7 * sumE`, then
  the occupied count is at most `7 * sumE`.

Root and audit wiring:

- `04-computation/lean/TournamentH7/TournamentH7.lean`
- `04-computation/lean/TournamentH7/TournamentH7/Verify.lean`

Build transcript:

- `05-knowledge/results/lrc14_arc_complexity_lean_codex_s87.out`

Verification command:

```bash
cd 04-computation/lean/TournamentH7
lake build TournamentH7.LRCArcComplexity TournamentH7.Verify \
  > /home/claude/math/05-knowledge/results/lrc14_arc_complexity_lean_codex_s87.out 2>&1
rg -n "warning:|error:|sorryAx|declaration uses .sorry|failed" \
  /home/claude/math/05-knowledge/results/lrc14_arc_complexity_lean_codex_s87.out
```

The `rg` scan returned no matches. The new declarations depend only on the
standard/classical Mathlib axioms shown by `#print axioms`
(`propext`, `Classical.choice`, `Quot.sound`).

## What this proves, and what it does not prove

Proved now:

- the finite disjointness-to-counting argument behind `V <= 7 * sumE`;
- the exact support lemma needed once actual `B_j(E')` regions are represented
  as unions of cells from one breakpoint partition;
- the fact that the factor-six loss is not forced by finite bookkeeping.

Not proved yet:

- the construction of the actual analytic regions `B_j(E')`;
- measurability and cell decomposition for those regions;
- the proof that every connected component or arc of `B_j(E')` is controlled by
  occupied cells of the shared breakpoint partition;
- the breakpoint bound `#cells <= 7 * sumE` for the concrete LRC sector-crossing
  partition;
- any direct improvement of `hp0cap` by itself. This is support infrastructure
  for the far-peel/decorrelation route, not a standalone LRC14 proof.

## Use in the LRC14 proof route

After the analytic side supplies:

1. a common breakpoint partition for the six exactly-one-miss sets,
2. pairwise disjointness of the `B_j`,
3. the component/readout inequality from arcs to occupied cells, and
4. the concrete breakpoint count `#cells <= 7 * sumE`,

the old THM-546 estimate

`V(E') <= 42 * sumE`

can be replaced rigorously by

`V(E') <= 7 * sumE`.

With the signed Abel constant this moves the conservative gapped single-peel
cutoff from the several-thousand range to roughly the recorded `~545` scale.
The empirical `V_actual ~ 4 * span` route would shrink the cutoff toward
`~80`, but that stronger structural claim is not formalized here.

## Assumption challenge / Tournament Analysis

I did not use runners or arcs as the tournament vertices for this formal layer.
The useful vertices are proof-obligation/cell families:

`common breakpoint partition > occupied cells > exactly-one-miss sector cells > component readout > far-peel error budget`.

This quotient preserves the LRC-relevant predicate "each exactly-one-miss sector
cell is counted at most once in the far-peel variation ledger." It destroys
speed ownership, interval adjacency, actual component geometry, phase values,
and the identity of the missed sector beyond disjoint cell membership.

Challenged assumption: the six sector ledgers must be paid independently. The
disjoint-cell quotient shows the factor `6` is a bookkeeping artifact once the
exactly-one-miss regions are subordinated to a shared partition.

Tournament fingerprints to record when the concrete `B_j` are built:

- score histogram of cells by how many sector ledgers attempt to claim them
  (target: only `0` or `1`);
- directed edge flips under refinement/coarsening of breakpoint cells;
- SCCs of proof obligations under implication;
- Hamiltonian path count for the proof-obligation ordering above, as a
  stability check on the preferred formal route.

## Next Lean targets

1. Define the concrete sector-crossing breakpoint cells for a finite integer
   list `E`.
2. Define `B_j(E)` as an exactly-one-miss cell family over that partition.
3. Prove pairwise disjointness of the six `B_j`.
4. Prove the concrete partition cardinality bound `#cells <= 7 * sumE`.
5. Bridge occupied cells to the analytic `#arcs(B_j)` or variation count used
   by THM-546.
6. Only then plug the resulting `V <= 7 * sumE` bound into the hp0cap/far-peel
   pipeline.
