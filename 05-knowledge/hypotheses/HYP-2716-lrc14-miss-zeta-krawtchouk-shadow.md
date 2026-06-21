---
id: HYP-2716
title: LRC14 miss-zeta discrepancy is a Krawtchouk top-character budget
status: OPEN; exact identity and scout evidence
source: codex-2026-06-21
depends_on:
  - HYP-2715
  - HYP-2714
  - HYP-2698
  - THM-351
  - THM-352
related:
  - THM-561
  - HYP-2718
  - HYP-2717
  - HYP-2676
  - HYP-2707
  - HYP-2710
  - THM-414
  - OPEN-Q-108
---

# HYP-2716 - Miss-Zeta Krawtchouk Shadow

## Claim

The remaining HYP-2714 multi-block carrier-product gap should be bounded as
one Boolean-cube character, not as a residual-coordinate `L1` discrepancy.

For a split row, write

```text
d_R = z_prod(R) - z_actual(R),     R subset {1,...,6},
W_j = sum_{|R|=j} d_R.
```

Let the binary Krawtchouk polynomial on the 6-cube be

```text
K_h(j) = sum_t (-1)^t binom(j,t) binom(6-j,h-t),
M_h    = sum_j W_j K_h(j).
```

Then

```text
M_h = sum_{|A|=h} sum_R d_R (-1)^|A cap R|
```

is the weight-`h` Walsh/MacWilliams shadow of the miss-zeta discrepancy, and
the actual cover error is exactly the top character

```text
Product(E)-p0(E) = M_6.
```

This is immediate from `K_6(j)=(-1)^j`.  The important point is proof
currency: `M_6` can be small even when the lower shadows `M_0,...,M_5` and the
coordinate `L1` discrepancy are much larger.  A signed Erdos-Turan/Koksma
argument should therefore target the all-six-sector character after routing
low-height resonances, not attempt to dominate every miss-zeta coordinate.

## Exact Evidence

The HYP-2715 scout now reports `W_j` and `M_h` exactly:

```text
04-computation/lrc14_multiblock_miss_zeta_layers_codex_20260621.py
05-knowledge/results/lrc14_multiblock_miss_zeta_layers_codex_20260621.out
```

Representative rows:

```text
two 4-blocks, moderate gap:
  M_6 = Product-p0 = -3495299/477209040
  |M_6| / sum_j |W_j| = 3495299/147256542 ~= 0.02374
  |M_6| / coordinate_L1 = 10485897/667494196 ~= 0.01571

3+3+2 split:
  M_6 = -51629953/57697542630
  |M_6| / sum_j |W_j| ~= 0.00493
  |M_6| / coordinate_L1 ~= 0.00345

five 2-blocks:
  M_6 = -2447628624709/93106921650624
  |M_6| / sum_j |W_j| ~= 0.07536
  |M_6| / coordinate_L1 ~= 0.06305
```

Across the six tested split rows, the aggregate character-weight tournament
orders the shadows as

```text
h=2 > h=0 > h=4 > h=1 > h=3 > h=5 > h=6.
```

Thus the cover-error character `h=6` is the smallest aggregate shadow in this
bank, even though the raw discrepancy has large lower-character content.

## Proof Route

The proposed proof obligation becomes:

```text
|M_6(E)| <= cap_k - Product(E)
```

for the HYP-2714 moderate-span balanced split branch.

A reasonable split is:

1. Keep the exact miss-zeta product factorization before scalarizing.
2. Express the orbit discrepancy in the Boolean character basis.
3. Prove cancellation for the `h=6` character using signed packet estimates,
   not coordinatewise domination.
4. Send low-gap or low-height resonances to the finite HYP-2714 ledger.

This lines up HYP-2715 with HYP-2676: both say that the dangerous term is a
signed packet/character after quotienting, and both lose too much if absolute
values are taken before the correct quotient.

HYP-2717 gives the next analytic refinement: expand this same `M_6` character
in carrier Fourier modes, then split exact carrier relations `n.M=0` from
nonresonant modes `n.M!=0`.  The exact relation modes are unavoidable for
integer carrier vectors, so the missing proof is a high-height Fourier tail
bound plus a finite low-height resonance ledger, not literal full-torus
equidistribution.

## Assumption Challenge and Tournament Analysis

The useful vertices are not runners, arcs, or raw residual masks.  This session
considered runners, carrier blocks, residual masks, residual-size layers,
Boolean characters, and proof obligations.  The chosen quotient has vertices
`h=0..6`, the Krawtchouk character weights.  It preserves the scalar cover
error at `h=6` and destroys sector labels and coordinatewise residual signs.

Pairwise observable: larger aggregate `sum |M_h|` over tested split rows.
Switch/gauge: pass from residual coordinates to the Boolean character shadow.
Tie Hamiltonian path:

```text
h=2 > h=0 > h=4 > h=1 > h=3 > h=5 > h=6.
```

Fingerprint: transitive tournament, score histogram
`{0:1,1:1,2:1,3:1,4:1,5:1,6:1}`, zero directed 3-cycles.

Challenged assumption: the multi-block discrepancy must be controlled as a
positive residual cone.  HYP-2716 replaces that with a top-character budget.

## Status

This is not yet the LRC(14) proof.  It is a sharper exact formulation of the
remaining analytic lemma.  The next work should produce a signed
Erdos-Turan/Koksma packet estimate directly for `M_6`, with finite resonance
exceptions handed to HYP-2714.
