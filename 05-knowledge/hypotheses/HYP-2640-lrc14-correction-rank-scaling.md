---
id: HYP-2640
title: LRC(14) correction rank scaling - raw rank is a switch, signed visible quotient is the ruler
status: OPEN
source: codex-2026-06-19-S28
depends_on:
  - HYP-2639
  - HYP-2638
  - HYP-2637
  - HYP-2636
  - HYP-2635
  - HYP-2633
related:
  - HYP-2607
  - HYP-2604
  - THM-531
  - THM-534
  - OPEN-Q-108
---

# HYP-2640 - LRC(14) Correction Rank Scaling

## Claim

The correction terms

```text
p0(E) - M7(k)
L_y(E) - L_y^inf(k)
```

are not controlled linearly by raw low-height relation-lattice rank. Raw
height-2 rank is a switch:

```text
low rank / uncovered fibre      -> independent-tail peel
saturated low-height rank       -> inverse-combinatorics / signed-coimage pocket
```

Once rank saturates, correction size is governed by a signed visible quotient:
observer-coupled fold multiplicity, summand-shell visibility, and mod-27
coimage phase. The rank supplies capacity; the visible quotient supplies
coherence.

The theorem target should therefore not be

```text
correction <= C * raw_relation_rank.
```

It should be the split lemma:

```text
height-2 rank/coverage not saturated
    -> peel a dissociated coordinate and use the independent baseline;

height-2 rank saturated
    -> apply an inverse theorem to get AP / small-excess / GAP / relation-covered structure;
       then show every non-AP saturated row loses either fold multiplicity/coherence
       or signed mod-27 coimage alignment.
```

## Evidence

The atlas
`04-computation/lrc14_relation_rank_correction_scaling_codex_s28.py` stores
its output at
`05-knowledge/results/lrc14_relation_rank_correction_scaling_codex_s28.out`.
It compares exact `p0`, `L_y`, independent baselines, height-2 exact relation
rank, fold and pair-collision motif rank/count, and mod-27 shell ranks for AP,
near-AP, GAP, third-pocket, and dissociated rows.

For `k=8`, the height-2 exact relation rank on the seven nonzero coordinates
is already saturated (`rE=6`, corank `1`) for AP, near-AP, GAP, and the two
wide third-pocket examples:

| row | excess | `rE` | fold count | exact rels | `p0-M7` | `L_y-L_y^inf` |
|---|---:|---:|---:|---:|---:|---:|
| AP | 0 | 6 | 12 | 1786 | 0.302731 | 0.308965 |
| nearAP top | 1 | 6 | 9 | 1496 | 0.249160 | 0.253795 |
| d2 GAP worst | 6 | 6 | 8 | 920 | 0.171140 | 0.175263 |
| third pocket A | 16 | 6 | 3 | 326 | 0.008019 | 0.013547 |
| third pocket B | 16 | 6 | 1 | 402 | 0.008255 | 0.007063 |
| ternary dissociated | 21 | 0 | 0 | 0 | 0.026667 | 0.021535 |

So the raw rank has the right qualitative tail split: the ternary
height-2-dissociated row has rank `0` and sits near the independent baseline.
But inside the relation-rich pocket the rank is constant while the correction
varies by a factor of more than `40`. AP and third-pocket A have the same raw
rank but `L_y` correction `0.308965` versus `0.013547`.

The visible/multiplicity data explain the gap better. AP has `12` fold motifs
at `k=8`; the worst d2 GAP has `8`; third-pocket A has `3`; third-pocket B has
`1`; the ternary row has none. Exact relation count falls in the same direction
(`1786 -> 920 -> 326/402 -> 0`), while raw rank is blind to the drop.

The mod-27 shell ranks also saturate (`mod27_observer` full in most structured
rows), which is another warning: even the finite coimage rank is not enough.
The proof must retain signed mass and phase, not only image dimension.

## Interpretation

The correction is best viewed as a signed trace/transfer through three maps:

```text
integer relation lattice
    -> finite coimage over Z/27Z
    -> seven-sector signed functional.
```

Rank measures the dimension of the domain or image. The LRC correction is the
signed mass that survives the last map. A large kernel can carry huge absolute
mass but tiny signed correction if the observer-visible phase cancels. This is
the same lesson as HYP-2633/HYP-2636 for the reciprocal tail: expose the channel
first, then scalarize.

## Proof Route

1. Use HYP-2637's weighted relation-fibre coverage to prove the low-rank peel:
   an uncovered or low-rank coordinate contributes only independent-tail sized
   correction.
2. Use HYP-2638 for the finite small-excess pocket.
3. Use the KPS Freiman-dimension/GAP penalty for proper GAP pockets.
4. For saturated but non-AP relation-covered rows, use HYP-2639's labelled
   relation hypergraph. The AP is the only row that keeps high fold
   multiplicity and signed coimage coherence at once.
5. Feed the remaining signed tail into HYP-2636's block-frequency transfer /
   HYP-2633's residue-lift summation-by-parts, before absolute values.

This does not prove LRC(14), but it tightens the last target: prove a
visible-coimage coherence deficit for non-AP saturated rank.

## Tournament Analysis

The pairwise observable ranks proof quotients by correction power, decoy
resistance, sign visibility, finite closure, and theorem readiness.

Hamiltonian path:

```text
coherent_visible_relation_rank
> freiman_corank_dimension_proxy
> mod27_observer_gcd3_tail
> fold_motif_rank
> balanced_pair_rank
> raw_height2_relation_rank
> raw_runner_vertices
```

Fingerprint: score histogram `{0:1,1:1,2:1,3:1,4:1,5:1,6:1}`, no directed
3-cycles, one Hamiltonian path.

Assumption challenge: the vertices are not runners. I considered runners, raw
relation vectors, folds, pair-collision shells, mod-27 residues, gcd strata,
Freiman corank, Fourier correction terms, and proof obligations. This quotient
preserves correction-carrying rank; it destroys exact time geometry and
high-height tails. The challenged assumption is that total short-relation rank
linearly predicts the LRC correction.
