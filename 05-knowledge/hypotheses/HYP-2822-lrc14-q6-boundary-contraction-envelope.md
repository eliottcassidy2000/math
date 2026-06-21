---
id: HYP-2822
title: LRC14 q6 boundary-contraction envelope for the gK8 route
status: OPEN theorem target; exact finite evidence and boundary correction
source: codex-2026-06-22-S79
depends_on:
  - HYP-2820
  - HYP-2816
  - HYP-2815
  - HYP-2814
  - HYP-2812
  - HYP-2811
  - HYP-2810
  - HYP-2808
related:
  - OPEN-Q-108
  - THM-563
  - THM-564
  - HYP-2798
  - HYP-2788
---

# HYP-2822: q6 Boundary-Contraction Envelope

## Claim

The gK8 concentration route should not use a uniform small-`f` slogan
`q6 -> q6/7`.  The correct boundary theorem appears to be a two-scale envelope:

```text
near boundary f>=15:
  q6 contracts, but the exact worst single-far ratio can be 14/15;

as f -> infinity:
  q6 contracts to the equidistribution factor 1/7 per independent far speed;

gK8 proof:
  use the exact boundary envelope plus positive L_yK8 margin at small f,
  then hand large f to the asymptotic/Krawtchouk/R-tail route.
```

This refines HYP-2816 and complements HYP-2820's endpoint-period q6-ratio
bound.  The small-`f` correction is real and large, but it is not a
counterexample to gK8.  It says the proof must keep the endpoint-period /
finite-window address at the boundary and only use the `1/7` law after the far
coordinate is genuinely decorrelated.

## Exact Evidence

Script:

```text
04-computation/lrc14_q6_boundary_contraction_codex_s79.py
```

Stored output:

```text
05-knowledge/results/lrc14_q6_boundary_contraction_codex_s79.out
```

The script validates its q6-only integer-grid engine against full `miss_dist`
on S78 leader rows, then scans two layers:

1. a frontier packet bank of current risk families;
2. an exhaustive bounded-base single-far q6-ratio scan for k=10,11,12,
   with all `B subset [0,14]`, `0 in B`, `|B|=k-1`, and `15<=f<=60`.

Exhaustive single-far result:

```text
k=10: rows=138138, max ratio=14/15
  B=(0,7,8,9,10,11,12,13,14), f=15,
  q6(B)=1/98, q6(B u {15})=1/105,
  L_yK8=11281937/2522520,
  10cap-L=3964063/2522520.

k=11: rows=92092, max ratio=14/15
  B=(0,6,7,8,9,10,11,12,13,14), f=15,
  q6(B)=1/98, q6(B u {15})=1/105,
  L_yK8=6563731/1261260,
  10cap-L=2583869/1261260.

k=12: rows=46046, max ratio=14/15
  B=(0,5,6,7,8,9,10,11,12,13,14), f=15,
  q6(B)=1/98, q6(B u {15})=1/105,
  L_yK8=379829/63063,
  10cap-L=160711/63063.
```

Frontier double-far result:

```text
adjacent top14/top-cluster packets reach q6 ratio 7/8
by adding {15,16} to a bounded base whose q6 interval is controlled by 14.
```

The S78 k=12 odd-bridge/gK8 leader is a different risk mode:

```text
B=(0,2,4,6,8,9,10,11,12,14), far pair {16,18},
q6(B)=1/98, q6(E)=1/126, ratio=7/9,
L_yK8=30494/4851,
10cap-L=11086/4851.
```

Thus the worst q6-ratio row is not the gK8-margin leader, and the gK8-margin
leader is not the worst q6-ratio row.  The proof needs a two-coordinate risk
ledger, not a scalar contraction ratio alone.

## Proposed Theorem Shape

The finite evidence suggests the following exact boundary lemma.

```text
If a bounded base B has its all-missed atom q6(B) controlled by a largest
bounded speed b<=14, then adding a single far speed f>b gives the boundary
envelope

    q6(B u {f}) / q6(B) <= b/f <= 14/15.

The equality family is the top-cluster or top14 base with b=14 and f=15.
```

For two adjacent far speeds the same dominant-interval heuristic predicts

```text
q6(B u {15,16}) / q6(B) <= 14/16 = 7/8,
```

which is exactly what the current frontier packet bank sees.

The missing rigor is to prove when the all-missed atom is controlled by the
dominant initial interval and to route non-dominant multi-interval q6 packets
through the finite atlas.  This is much sharper and more honest than trying to
force a false small-`f` `1/7` bound.

## Post-Rebase Seam With HYP-2820

Incoming mac-mini S23 work clarifies the method boundary.  HYP-2820 proves the
single-far q6-ratio statement by the THM-563 endpoint-period method: the base
set `Q6(B)` is fixed, so `f*(q6(B u {f})-(1/7)q6(B))` is periodic in `f`.
That is exactly the mechanism HYP-2822 should use for the single-far boundary
envelope.

The adjacent two-far `7/8` row should not be forced through the same period-max
tool.  For `{M,M+g}` the second far coordinate creates `M`-dependent
breakpoints, so the correction is almost-periodic/R-tail data rather than a
fixed-base periodic sequence.  Thus HYP-2822 splits into two proof obligations:

1. single far: endpoint-period q6 envelope, aligned with HYP-2820;
2. two or more fars: finite boundary atlas plus the HYP-2817/HYP-2818/HYP-2819
   R-tail ledger, with the structural dichotomy/covolume route deciding when a
   row belongs to the single-far side.

## Tournament Analysis

The script's tournament vertices are proof packets, not runners:

```text
single_consecutive, single_top14, single_top_cluster,
double_consecutive, double_top14, double_top_cluster,
double_even_ap, double_s78_leader_base.
```

The pairwise observable is the risk tuple `(minimum gK8 margin, maximum q6
ratio)`.  The switch orients toward smaller margin, then larger q6 ratio, with
lexicographic tie path.

The resulting order separates two risks:

- top14/top-cluster single-far packets maximize q6 ratio at `14/15`;
- S78 odd-bridge packets can be gK8-margin leaders while having smaller q6
  ratios such as `7/9`.

This preserves the LRC predicate `L_yK8<=10cap` and destroys the misleading
single scalar predicate "q6 ratio close to 1/7".

## Status

No LRC14 proof is claimed.  HYP-2822 contributes a concrete finite correction
to the gK8 concentration proof order:

1. prove the boundary q6 envelope (`14/15`, `7/8`, and its finite-atlas
   variants);
2. combine that envelope with exact gK8 margin, not q6 alone;
3. use asymptotic `1/7` equidistribution and the `12*zeta(3)` R-tail only after
   the boundary window is discharged.
