---
id: HYP-2637
title: LRC(14) relation-fiber / GAP split - weighted summand energy replaces the failed dissociated-stranger dichotomy
status: OPEN proof program; PARTIALLY-CONFIRMED by local relation-fiber scout
source: codex-2026-06-19-S26
depends_on:
  - HYP-2635
  - HYP-2634
  - HYP-2636
  - HYP-2607
  - HYP-2606
  - HYP-2599
related:
  - THM-446
  - THM-442
  - THM-414
  - THM-515
  - OPEN-Q-108
---

# HYP-2637 - LRC(14) Relation-Fiber / GAP Split

## Claim

The failed wide/dissociated dichotomy in HYP-2635 should be replaced by a
weighted summand-fiber split.

For a cluster offset set `E`, define the bounded weighted summand fibers

```text
c in N^k, 0 < sum c_i <= M, c_i <= H,
c |-> sum_i c_i e_i.
```

Collisions in these fibers are exactly bounded-height integer relations

```text
sum_i n_i e_i = 0.
```

Ordinary pair-sum energy is only the `sum c_i=2`, coefficient-one slice.  It
misses relations such as `2x+y=z`, which are visible to the LRC Fourier/orbit
relation lattice and are exactly the obstruction to peeling a dissociated
stranger.

The proposed proof split is:

1. If some nonzero element is not touched by bounded weighted-summand
   relations, peel it and use the dissociated / independent-limit estimate.
2. If every nonzero element is touched, the set has high bounded additive
   relation energy.  Use a Freiman / Balog-Szemeredi-Gowers style inverse lemma
   to place it inside a low-rank GAP, then bound the resulting GAP families via
   AP-orbit invariance and dimension.
3. In reciprocal lift tails, retain the multiplicative sign coordinate:

```text
term sign = sign(residue coefficient) * (-1)^(number of negative relation coefficients).
```

So addition creates relation fibers; multiplication creates denominator parity,
divisibility, and the positive/negative signed lift.

This does not prove LRC(14).  It sharpens the next proof obligation after
HYP-2635 and connects it to HYP-2634's low-height lift-opposition mechanism.

## Computation

Script:

- `04-computation/lrc14_relation_fiber_gap_codex_s26.py`
- output: `05-knowledge/results/lrc14_relation_fiber_gap_codex_s26.out`

The script compares ordinary pair-sum energy with weighted summand-fiber energy
on the two HYP-2635 no-dissociated-stranger examples.

```text
case              span  pair_E  pair_col  wt_E  wt_colfib  maxfib  h2_cover
kp_example_A        35      40         6   913         73       7      7/7
kp_example_B        31      34         3  1075         68       7      7/7
consecutive_8        7      72        22  4177         23      22      7/7
dissociated_powers15625     28         0   273          0       1      0/7
```

`kp_example_A = (0,3,5,16,28,30,33,35)` has a rank-2 scout cover with steps
`(3,5)` and box `(2,7)`.  Its ordinary pair-sum slice sees only a few
collisions, but the weighted relation fibers touch every nonzero vertex.  Sample
relations include:

```text
5 + 30 = 35
5 + 28 = 33
3 + 30 = 33
28 + 35 = 30 + 33
```

The missing relation for the apparent outlier is also weighted-small:

```text
2*16 + 3 = 35.
```

`kp_example_B = (0,4,12,15,20,21,25,31)` likewise has full height-2 coverage,
with sample relations:

```text
4 + 21 = 25
15 + 31 = 21 + 25
4 + 31 = 15 + 20
4 + 12 + 15 = 31
```

By contrast, the dissociated row `(0,1,5,25,125,625,3125,15625)` has no
weighted collisions at this height and no nonzero relation coverage.

## Additive Versus Multiplicative Sign

The same script revisits HYP-2634's QR pair:

```text
a=2: E=(1,2,8,9,15,22)
a=4: E=(1,4,8,11,15,22)
```

Both rows are relation-dense in weighted summand fibers through mass `7`, but
their dangerous named motifs differ.  The sign table makes the
addition/multiplication split explicit:

```text
a2_negative:
  2*2 + 2*22 = 1 + 8 + 9 + 2*15
  denominator parity positive, residue coefficient negative, term negative

universal_positive:
  1 + 2 + 2*22 = 8 + 9 + 2*15
  denominator parity negative, residue coefficient negative, term positive

a4_h3_negative:
  3*4 + 2*15 = 1 + 8 + 11 + 22
  denominator parity positive, residue coefficient negative, term negative

a4_h4_negative:
  4*4 + 3*11 = 4*1 + 8 + 15 + 22
  denominator parity positive, residue coefficient negative, term negative
```

Thus the summand graph supplies the equality, but positive/negative lift mass
is decided by the multiplicative denominator parity and the finite residue
coefficient.  This is the requested even/odd to positive/negative bridge.

## Tournament Analysis

Candidate vertices included runners, gaps, ordinary pair sums, weighted
summand fibers, relation vectors, Fourier modes, residue characters,
denominator products, and proof obligations.

Chosen Hamiltonian path:

```text
weighted_summand_fibres
> height2_relation_coverage
> freiman_gap_scout
> multiplicand_sign_parity
> finite_residue_character
> ordinary_pair_sums
> raw_speed_vertices
```

Fingerprints:

```text
score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1}
directed_3_cycles = 0
SCC_sizes = [1,1,1,1,1,1,1]
hamiltonian_paths = 1
```

The challenged assumption is that ordinary pair sums are enough.  They are not:
the live LRC predicate is preserved by weighted summand fibers plus
multiplicand sign parity, and is destroyed by reducing to raw runners or
ordinary summand edges.

## Next Proof Obligations

1. Define a quantitative bounded relation-coverage threshold:

```text
Cov_H,M(E) = {nonzero e_i touched by some primitive relation from weighted fibers}.
```

2. Prove a peel lemma when `Cov_H,M(E)` is not full.
3. Prove an inverse theorem when `Cov_H,M(E)` is full: bounded relation coverage
   plus enough fibers implies containment in a low-rank GAP with controlled
   rank/volume.
4. Bound `L_y` or `meas(S7)` on those GAP families using AP-orbit invariance and
   dimension.
5. For HYP-2633/HYP-2634 reciprocal tails, split every low-height relation ledger
   by denominator parity and residue-character sign before applying Abel
   summation to the remaining high-height tail.
