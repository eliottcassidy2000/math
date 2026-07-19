---
id: THM-1134
title: Five-residue multiplier Kakeya cone and the exact all-scale r6 step-two family
status: PROVED / FINITE-EXACT — a ten-orbit affine certificate gives a six-unit gap after a nonzero multiplier on Z/13; an all-792-core rectangle atlas plus a 12,771-row exact complement closes the step-two five-comb family at every legal scale
source: codex-2026-07-18-S73 with death-star-S58 analytic handoff
depends_on: [THM-1128, THM-1132]
related: [THM-1102, THM-1121, THM-1126]
verification:
  - 04-computation/lrc14_r6_five_residue_kakeya_exact_codex_S73.py
  - 05-knowledge/results/lrc14_r6_five_residue_kakeya_exact_codex_S73.out
  - 04-computation/lrc14_r6_step2_rectangle_atlas_exact_codex_S73.py
  - 05-knowledge/results/lrc14_r6_step2_rectangle_atlas_exact_codex_S73.out
  - 04-computation/lrc14_r6_step2_finite_complement_exact_codex_S73.cpp
  - 05-knowledge/results/lrc14_r6_step2_finite_complement_exact_codex_S73.out
  - 04-computation/lean/TournamentH7/TournamentH7/LRCSharpCombArithmetic.lean
---

# THM-1134 — multiplier Kakeya cone and r6 step-two family

This theorem records two rigorous advances on the open uniform `r=6` sharp
horn.  The first repairs the apparent five-centre failure of the fixed
`t=1/13` chart by allowing all nonzero thirteenth-grid charts.  The second
uses fixed core-safe rectangles plus a finite exact complement to close the
entire step-two family containing the sampled maximum in THM-1132.

## 1. The five-residue multiplier lemma

For every subset `A` of `Z/13Z` with at most five elements, some nonzero
multiplier `u` makes the largest cyclic gap of `uA` at least six grid units.
It is enough to prove the five-element case: extend a smaller set to five
elements, choose the chart there, and then delete the added points; deleting
points cannot decrease the largest cyclic gap.

The statement has a compact exhaustive certificate.  Translation and
nonzero multiplication preserve cyclic gap multisets up to choosing the
chart.  The `1,287` five-subsets split into exactly ten affine orbits:

```text
representative       orbit size   witness u   gap word
0 1 2 3 4                  78          1      1 1 1 1 9
0 1 2 3 5                 156          1      1 1 1 2 8
0 1 2 3 6                 156          1      1 1 1 3 7
0 1 2 3 7                 156          2      1 1 2 2 7
0 1 2 3 8                  78          2      2 1 1 2 7
0 1 2 4 5                 156          1      1 1 2 1 8
0 1 2 4 7                 156          1      1 1 2 3 6
0 1 2 4 10                156          1      1 1 2 6 3
0 1 2 5 6                 156          1      1 1 3 1 7
0 1 2 6 9                  39          2      2 2 1 7 1
```

The orbit sizes sum to `1,287`, and the exact referee constructs every orbit,
checks disjointness/coverage, and replays each displayed gap word.  The lower
bound six is sharp: `A={0,1,2,4,7}` has maximum largest gap exactly six over
all twelve multipliers.

Equivalently, the eight-point complement of every five-subset of `Z/13Z`
contains a five-term affine interval.  This is the useful additive form of
the finite lemma.

## 2. A uniform five-comb Kakeya cone

Let a seven-speed core `P` lie in `{1,...,12}` and let

```text
k1<k2<k3<k4<k5.
```

Choose an integer `B in [k1,k5]`, put `a_i=k_i-B`, and define

```text
A=max_i |a_i|,             M=max(A,80).
```

Apply the multiplier lemma to the set of distinct residues `{a_i mod 13}`
(repetitions only reduce its size) and choose
`t0=u/13`.  Every core speed is at distance at least `1/13` from an integer at
`t0`.  On

```text
I=[t0-epsilon,t0+epsilon],   epsilon=10/(273M),
```

the core remains `1/14`-safe because `12 epsilon<=1/182`.  The six-unit
centre gap at `t0` has vertical safe width

```text
6/13-1/7=29/91.
```

Every offset centre moves by at most `A epsilon`.  Moreover
`2A epsilon<=20/273<1/13`, so distinct thirteenth-grid centre clusters cannot
cross or enter the selected gap on `I`.  The intersection of that gap over
`I` therefore contains a fixed vertical arc `X` of length

```text
|X|>=29/91-2A epsilon>=67/273.                         (1)
```

If

```text
B>=17M,                                                   (2)
```

then

```text
B|I|>=340/273=1+67/273.
```

Consequently the slope-`B` needle `x=Bt mod 1` has a complete preimage of
`X` inside `I`.  That preimage is contained in a component of the five-comb
safe set, so the component has length at least

```text
67/(273B)>1/(7k5).                                      (3)
```

Thus a sixth killer larger than `k5` cannot cover that component.

With `Delta=k5-k1`, choose the upper integer midpoint
`B=k1+ceil(Delta/2)`.  Then (2) follows from the transparent sufficient
condition

```text
k1>=max(1360,8(Delta+1)).                               (4)
```

This is the first arbitrary-shape five-comb cone.  It does not close uniform
`r=6`; the finite atlas and shapes outside (4) remain separate obligations.

The constants are sharp within this fixed universal-chart/no-crossing
rectangle scheme.  Put `y=2M epsilon` and `r=B/M`.  Core safety gives
`y<=M/1092`, order preservation requires `y<1/13`, the vertical width is at
most `29/91-y`, and a complete crossing requires

```text
(r+1)y>=1+29/91=120/91.
```

At `r=16` this forces `y>1/13`, so no such chart works.  At `r=17` the least
possible value is `y=20/273`, and core safety then forces `M>=80`.  Thus both
`17` and `80` are explained by the two competing walls rather than selected
by a numerical search.

### A complementary separated-ratio gate

The cone controls near-diagonal shapes.  A direct extension of THM-1097's
mass/component argument gives a rigorous gate in the opposite direction.
Let `J` be a nonwrapping core-safe interval of length `ell`, let
`k1<...<k5=K`, and put

```text
Q5=ell(14K-7 sum_i ki)-6K sum_i(1/ki)-47.              (5)
```

Then `Q5>0` forces a surviving component longer than `1/(7K)`.

Indeed, the sharp one-comb discrepancy and the union bound give survivor
mass

```text
A>=2ell/7-(6/49)sum_i(1/ki).                           (6)
```

A tooth of `D_ki` meeting `J` has its centre index in an interval of length
`ki ell+1/7`, so there are at most `ki ell+8/7` such teeth.  Five combs have
at most `ell sum_i ki+40/7` meeting teeth.  Their union has no more components
than teeth, and its complement in connected `J` has at most one additional
component.  Hence

```text
C<=ell sum_i ki+47/7.                                  (7)
```

The exact identity between the right sides is

```text
7K*(right side of (6))-(right side of (7))=Q5/7.
```

Thus `Q5>0` gives `7KA>C`, so one of the `C` survivor components has length
greater than `1/(7K)`.  Endpoint-only teeth are harmless overcounts, and
open/closed conventions do not change mass or positive-length component
count.  This gate does not by itself cover every shape; paired with the
multiplier cone it makes the remaining intermediate-ratio obstruction
explicit rather than calling the whole tail “uniform”.

## 3. The step-two five-comb family is all-scale

### An explicit rectangle for the formerly sampled worst core

Take

```text
P={1,2,4,7,9,11,12},
(k1,...,k5)=(b,b+2,b+4,b+6,b+8),
B=b+4,
offsets=(-4,-2,0,2,4).
```

The fixed interval

```text
J=[15/98,9/56],               |J|=3/392                (8)
```

is core-safe.  On it `2t in [15/49,9/28]`, and the cyclic centre gap from
`0` to `2t` contains the fixed vertical arc

```text
X=[1/14,15/49-1/14],          |X|=8/49.                (9)
```

For every `b>=148`,

```text
B|J|=(b+4)3/392>=57/49=1+|X|,
```

so a complete needle preimage lies in `J`.  Its length satisfies

```text
8/(49(b+4))>1/(7(b+8)).                                (10)
```

The covering-killer ray is legal only for `b>13 max(P)=156`; hence (8)--(10)
prove it for every legal `b>=157`.  This replaces the finite `b<=4000`
telemetry and the unquantified drift paragraph in THM-1132 by an exact
all-scale rectangle proof.

### Uniform rectangle atlas over all cores

The same fixed-polygon idea can be made uniform over all `C(12,7)=792`
seven-speed cores.  In the centred coordinates

```text
B=b+4,               offsets=(-4,-2,0,2,4),
```

split `[0,1]` at every possible core-tooth endpoint and every collision of
two offset centres.  On each live order cell, and for each labelled cyclic
gap and inward side, both bounding walls are affine.  If the endpoint safe
width is `w` and the inward retreat rate is `rho`, a one-sided interval of
length `d` carries a fixed vertical arc of width `w-rho*d`.  Maximising the
exact crossing certificate over those finite candidates gives

```text
B>=168                                                   (11)
```

for every core.  The unique worst certificate is

```text
P={1,2,3,4,7,9,11},
I=[13/168,13/154],       |I|=13/1848,
X=[9/22,13/22],          |X|=2/11,
168|I|=13/11=1+|X|.                                      (12)
```

Thus every `b>=164` has a complete safe preimage of length
`|X|/B>=1/(7B)>1/(7(B+4))=1/(7(b+8))`.

### Independent finite complement

A separate C++ endpoint-subtraction referee checks every legal row with
`b<=164`, without using the rectangle atlas.  Across all 792 cores this is

```text
12,771 rows,             failures of 7(b+8)L>1: 0.       (13)
```

The finite-bank extremal is

```text
P={1,2,4,7,9,11,12}, b=158,
killers=(158,160,162,164,166),
L=67/61992,             7(b+8)L=5561/4428.              (14)
```

The atlas and finite referee overlap at `b=164`; together they cover every
legal base with no boundary convention gap.  Optimized, unoptimized, and
sanitized C++ runs are byte-identical; normal and optimized Fraction-atlas
runs are byte-identical.

## 4. Structural reading

The key challenged assumption is that the chart `t=1/13` must stay fixed.
The proof vertices are residue patterns together with a selectable
multiplicative chart.  The affine-orbit quotient preserves the maximum-over-
charts gap predicate, while a single residue order or a tournament on the
five killers destroys that freedom.  Once a chart is chosen, the proof-
bearing Kakeya carrier is again the labelled cyclic gap with endpoint owners,
offset slopes, and the slope-`B` needle.  Its naked boundary-order tournament
is transitive and loses every metric inequality used above.

The normal and optimized multiplier-certificate runs are byte-identical to
the frozen output.  SHA-256 values for the multiplier certificate, all-core
rectangle atlas, and finite complement (source/output pairs) are

```text
7d6b12d3ac58f9e582c007d27a925e1c9f0f81138e854cf633ca50614fecee7d
4f192dcee76fe7cbebe7a835958e39f74e71ea86258acd64a1f651dd76465880
f9901493e3386fd8f58f9d429b79bbdeb7da9447fa9ba87c209443a44e01fbea
7b526996945c20d302ee28ece52069f10d53c899d2a6ceef75533bc5e275d08b
ff0fc4cd717fa01ae722f117e697bca5c7b5028371f0bc26a3ca1a93dcc09fc6
eaeefdcdbefce78e61d5d299e858d8914ff50da64b5eb9e705fb8ee4dc94d7c4
```
