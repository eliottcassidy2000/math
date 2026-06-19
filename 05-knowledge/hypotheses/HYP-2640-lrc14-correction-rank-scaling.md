---
id: HYP-2640
title: LRC(14) correction rank scaling - raw rank is a switch, signed/coset quotient is the ruler
status: PARTIALLY-TRUE; two exact scouts complete; proof obligation refined
source: codex-2026-06-19-S28
depends_on:
  - HYP-2639
  - HYP-2638
  - HYP-2637
  - HYP-2636
  - HYP-2635
  - HYP-2633
  - HYP-2606
related:
  - HYP-2607
  - HYP-2604
  - THM-531
  - THM-534
  - OPEN-Q-108
---

# HYP-2640 - LRC(14) Correction Rank Scaling

## Result

The correction terms

```text
p0(E) - M7(k)
L_y(E) - L_y^inf(k)
```

are not controlled linearly by raw low-height relation-lattice rank.  The
right statement is two-stage:

```text
raw relation rank / coverage is a switch;
signed visible/coset phase is the ruler.
```

Low rank or uncovered weighted fibres route to an independent-tail peel.  Once
low-height rank saturates, the correction size is governed by a signed quotient:
observer-coupled fold multiplicity, integer summand-shell visibility,
weighted-fibre coherence, and finite coimage phase.

So the theorem target is not:

```text
correction <= C * raw_relation_rank.
```

It is the split lemma:

```text
height-2 rank/coverage not saturated
    -> peel a dissociated coordinate and use the independent baseline;

height-2 rank saturated
    -> apply inverse combinatorics to get AP / small-excess / GAP /
       relation-covered structure;
       then show every non-AP saturated row loses either fold coherence,
       signed/coset alignment, or weighted-fibre phase before scalarization.
```

## Evidence A: Height-2 Rank Atlas

The concurrent atlas
`04-computation/lrc14_relation_rank_correction_scaling_codex_s28.py` stores
output at
`05-knowledge/results/lrc14_relation_rank_correction_scaling_codex_s28.out`.
It compares exact `p0`, `L_y`, independent baselines, height-2 exact relation
rank, fold and pair-collision motif rank/count, and mod-27 shell ranks for AP,
near-AP, GAP, third-pocket, and dissociated rows.

For `k=8`, the height-2 exact relation rank on the seven nonzero coordinates is
already saturated (`rE=6`, corank `1`) for AP, near-AP, GAP, and the two wide
third-pocket examples:

| row | excess | `rE` | fold count | exact rels | `p0-M7` | `L_y-L_y^inf` |
|---|---:|---:|---:|---:|---:|---:|
| AP | 0 | 6 | 12 | 1786 | 0.302731 | 0.308965 |
| nearAP top | 1 | 6 | 9 | 1496 | 0.249160 | 0.253795 |
| d2 GAP worst | 6 | 6 | 8 | 920 | 0.171140 | 0.175263 |
| third pocket A | 16 | 6 | 3 | 326 | 0.008019 | 0.013547 |
| third pocket B | 16 | 6 | 1 | 402 | 0.008255 | 0.007063 |
| ternary dissociated | 21 | 0 | 0 | 0 | 0.026667 | 0.021535 |

Thus raw rank has the right qualitative tail split: the ternary height-2
dissociated row has rank `0` and sits near the independent baseline.  But
inside the relation-rich pocket the raw rank is constant while the correction
varies by a factor greater than `40`.  AP and third-pocket A have the same raw
rank but `L_y` correction `0.308965` versus `0.013547`.

The visible/multiplicity data explain the gap better.  AP has `12` fold motifs
at `k=8`; the worst d2 GAP has `8`; third-pocket A has `3`; third-pocket B has
`1`; the ternary row has none.  Exact relation count falls in the same
direction (`1786 -> 920 -> 326/402 -> 0`), while raw rank is blind to the drop.

The mod-27 shell ranks also saturate in most structured rows, so finite coimage
rank is not enough either.  The proof must retain signed mass and phase, not
only image dimension.

## Evidence B: Pair-Sum And Weighted Rank Scout

The complementary scout
`04-computation/lrc14_correction_rank_scaling_codex_s28.py` stores output at
`05-knowledge/results/lrc14_correction_rank_scaling_codex_s28.out`.  It
computes exact `L_y`, exact `p0=meas(S_7)`, independent baselines, pair-sum
Freiman relation rank, visible fold rank, hidden pair-sum shell rank, and
bounded weighted relation rank on named rows.

Named-row evidence:

- AP rows maximize correction inside the scan:
  - `AP8`: `Corr_y=0.308964891`, pair rank `6`, visible `6`, hidden `5`.
  - `AP9`: `Corr_y=0.341582294`, pair rank `7`, visible `7`, hidden `6`.
  - `AP10`: `Corr_y=0.351249852`, pair rank `8`, visible `8`, hidden `7`.
- KPS third-pocket rows have small pair-rank correction but full bounded
  weighted rank:
  - `KPS_third_A`: pair rank `3`, weighted rank `6`, `Corr_y=0.013546678`.
  - `KPS_third_B`: pair rank `5`, weighted rank `6`, `Corr_y=0.007062543`.
- The k=8 default bank gives the decisive guardrail against visible-rank-only:
  `E=(0,1,3,5,7,9,11,13)` has total pair rank `5`, visible rank `0`, hidden
  rank `5`, and still has `Corr_y=0.215709575`.  Hidden integer odd-coset AP
  shells can be positive correction payload, not merely cancellation.

Default exact banks:

- k=8, max coordinate 14, `3431` primitive rows: AP is top, but the
  `(total,visible,hidden)=(5,0,5)` packet is third by max correction.
- k=9, max coordinate 13, `1287` primitive rows: AP is top, and lower visible
  ranks still retain large positive corrections.
- k=10, max coordinate 12, `220` primitive rows: all rows in the small bank
  have maximal total rank, so visible/hidden packet data is the observed
  separator.

## Reconciliation With KPS S12

This does not contradict the KPS S12 correction that the mod-27 antipodal
summand graph is inert for small binding clusters.  There are two shell notions
in play:

- KPS's inertness statement is about `{a,27-a}` residue shells for elements
  below `27/2`.
- HYP-2640's hidden ranks are integer pair-sum and weighted-fibre shells inside
  `E+E` and bounded summand fibres.

The binding correction comes from genuine integer relations, with the mod-27
factorization diagnosing why n=14 is hard rather than carrying the small-cluster
correction by itself.

## Interpretation

The correction is best viewed as a signed trace/transfer through several maps:

```text
integer relation lattice
    -> finite coimage / residue shell data
    -> seven-sector signed functional.
```

Rank measures dimension of the domain or image.  The LRC correction is the
signed mass that survives the last map.  A large kernel can carry huge absolute
mass but tiny signed correction if the observer-visible phase cancels.  This is
the same lesson as HYP-2633/HYP-2636 for the reciprocal tail: expose the channel
first, then scalarize.

## Refined Proof Target

```text
low relation rank                    -> independent-limit / peel;
high rank + AP/odd-coset phase       -> finite AP or near-AP envelope;
high rank + cancelling shell phase   -> dimension penalty / signed tail bound;
small pair rank + high weighted rank -> weighted-fibre GAP pocket.
```

A reasonable analytic form is:

```text
Corr(E) <= sum_{typed relation packets P} phase_coeff(P) * rank_capacity(P)
           + independent/peel error,
```

where `phase_coeff(P)` depends on fold visibility, shell parity, coset offset,
coimage phase, and multiplicand visibility, while `rank_capacity(P)` records
only the available relation capacity.

## Proof Route

1. Use HYP-2637's weighted relation-fibre coverage to prove the low-rank peel:
   an uncovered or low-rank coordinate contributes only independent-tail sized
   correction.
2. Use HYP-2638 for the finite small-excess pocket.
3. Use the KPS Freiman-dimension/GAP penalty for proper GAP pockets.
4. For saturated but non-AP relation-covered rows, use HYP-2639's labelled
   relation hypergraph.  The AP is the only row that keeps high fold
   multiplicity and signed coherence at once.
5. Feed the remaining signed tail into HYP-2636's block-frequency transfer /
   HYP-2633's residue-lift summation-by-parts, before absolute values.

This does not prove LRC(14), but it tightens the last target: prove a
visible/coset coherence deficit for non-AP saturated rank.

## Tournament Analysis

The pairwise observable ranks proof quotients by correction power, decoy
resistance, sign visibility, finite closure, and theorem readiness.

Hamiltonian path:

```text
coherent_visible_relation_rank
> signed_coset_rank
> weighted_relation_rank
> freiman_corank_dimension_proxy
> mod27_observer_gcd3_tail
> fold_motif_rank
> balanced_pair_rank
> raw_height2_relation_rank
> raw_additive_energy
> raw_runner_vertices
```

Assumption challenged: tournament vertices are not runners or arcs here.  They
are proof quotients.  The quotient preserves correction power and destroys raw
runner identity, which is intentional because raw runner identity has already
proved too coarse for this proof obligation.
