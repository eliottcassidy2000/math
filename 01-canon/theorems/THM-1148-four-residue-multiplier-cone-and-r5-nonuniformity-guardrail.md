---
id: THM-1148
title: Four-residue multiplier cone, Q4 separated tail, and exact-transfer cone for the open r=5 stratum
status: PROVED three uniform four-comb cones and exact nonuniformity counterexample; FINITE-EXACT affine-orbit/core/ray ledgers. Uniform r=5 remains OPEN
source: codex-2026-07-18-S74
depends_on:
  - THM-1128  # earlier thirteen-grid rectangle gate
  - THM-1133  # exact all-scale span-at-most-30 bank
  - THM-1137  # corrected exact interval transfer Phi
related:
  - THM-1097  # sharp one-comb discrepancy
  - THM-1140  # corrected four-comb reconnaissance
  - THM-1141  # actual-mean proposal refuted here
  - THM-1147  # exact two-comb endpoint law
artifact: 03-artifacts/lrc14-r5-nonuniformity-refutation-and-four-residue-cone-codex-20260718.md
script: 04-computation/lrc14_r5_nonuniformity_refutation_multiplier_cone_codex_20260718.py
output: 05-knowledge/results/lrc14_r5_nonuniformity_refutation_multiplier_cone_codex_20260718.out
---

# THM-1148 -- three exact cones and the surviving four-comb shape

Write

```text
D_k={t in R/Z : ||kt||<1/14}.
```

The uniform `r=5` problem asks whether, after an eight-speed core and four
ordered killers `k1<k2<k3<k4`, some surviving component is longer than
`1/(7k4)`.  Such a component cannot be covered by the fifth danger comb.

This theorem proves three complementary cones, corrects the proposed
actual-mean argument, and gives an explicit infinite shape outside the union.
It does not close uniform `r=5`.

## 1. The ordinary `4/3` nonuniformity claim is false

Take

```text
P={1,2,3,4,5,6,7,8},
J=[1/14,13/112],
(k1,k2,k3,k4)=(108,109,110,111).
```

`J` is a maximal core-safe component and the killers are legal because
`108>13 max(P)`.  Exact endpoint subtraction gives five survivor lengths

```text
319/55944, 305/55944, 97/18648, 277/55944, 13/3024.
```

Hence

```text
mass A       =955/37296,
mean A/5     =191/37296,
longest L    =319/55944,
L/(A/5)      =638/573<4/3.                              (1)
```

This refutes THM-1141's proposal when “mean” means the actual mean survivor
length.  It does not refute the desired metric conclusion; here

```text
7k4L=319/72>1.                                          (2)
```

The missing datum is the tooth-cluster count.  Nearby teeth coalesce, so low
gap variance can coexist with a large mean.  On a cyclic chart of
circumference `1/k1` containing one full tooth from each comb, if the four
teeth form at most three connected clusters, the longest complementary gap
is already greater than `1/(7k4)`.  Indeed the total safe length is at least

```text
1/k1-(1/7)sum_i 1/ki,
```

and for three clusters

```text
7k4/k1-k4 sum_i(1/ki)
  >=4k4/k1-1>3.
```

Boundary fragments and charts containing multiple teeth of one comb require
their own ledger, so this last observation is local rather than a global
four-comb theorem.

## 2. Sharp four-residue multiplier lemma

> **Lemma A.**  For every nonempty subset `A` of `Z/13Z` with at most four elements,
> some nonzero multiplier `u` makes the largest cyclic gap of `uA` at least
> seven grid units.  Seven is sharp.

It is enough to check four-sets.  The 715 labelled four-subsets split into
seven affine orbits:

```text
representative  orbit size  witness u  gap word
0 1 2 3              78          1     1 1 1 10
0 1 2 4             156          1     1 1 2 9
0 1 2 5             156          1     1 1 3 8
0 1 2 6             156          2     2 2 8 1
0 1 3 4              78          1     1 2 1 9
0 1 3 9              52          2     2 3 1 7
0 1 3 11             39          1     1 2 8 2.
```

The orbit sizes sum to 715.  The sixth orbit proves sharpness.  The exact
minimum best gaps at cardinalities `1,2,3,4` are respectively

```text
13,12,9,7.                                              (3)
```

## 3. Multiplier-Kakeya cone

Let `P subset {1,...,12}` be any core and let

```text
k1<k2<k3<k4,
B in [k1,k4],
ai=ki-B,
A=max_i |ai|,
M=max(A,84).
```

> **Theorem B.**  If `B` is an integer and
>
> ```text
> B>=15M,                                                (4)
> ```
>
> then the core-safe set after deleting the four danger combs has a
> component of length at least
>
> ```text
> 11/(73B)>1/(7k4).                                     (5)
> ```

**Proof.**  Apply Lemma A to `{ai mod 13}` and choose `t0=u/13` carrying a
seven-unit cyclic centre gap.  Put

```text
epsilon=14/(365M),       I=[t0-epsilon,t0+epsilon].
```

Every core speed is safe on `I`, since at `M=84`

```text
1/13-12epsilon=339/4745>1/14.                           (6)
```

The initial vertical safe width between the four offset danger combs is
`7/13-1/7=36/91`.  Both bounding walls drift by at most

```text
2Aepsilon<=28/365<1/13,                                 (7)
```

so their grid order cannot cross.  The common gap has width

```text
36/91-28/365=10592/33215>11/73.                         (8)
```

Choose a fixed vertical arc `X` of length `11/73` there.  By (4),

```text
B|I|>=15(28/365)=84/73=1+|X|.                           (9)
```

Therefore the slope-`B` needle has a complete preimage of `X` inside `I`.
That preimage is safe for the core and all four killers and has length
`|X|/B`.  Since `B<=k4` and `77>73`, (5) follows.  ∎

A transparent midpoint corollary, with `Delta=k4-k1`, is

```text
k1>=max(1260,7(Delta+1)).                              (10)
```

The optimal integer choice of `B` within this chart is

```text
Delta<=84:        (B,M)=(k4,84),
85<=Delta<=168:   (B,M)=(k1+84,84),
Delta>=169:       (B,M)=(k1+ceil(Delta/2),ceil(Delta/2)). (11)
```

Thus the exact piecewise gate is respectively

```text
k4>=1260,       k1>=1176,       k1>=14ceil(Delta/2).    (12)
```

The constants `15` and `84` are jointly optimal **within this fixed
universal, complete-preimage, no-grid-crossing scheme**.  Ratio 14 cannot
supply the more-than-`8/7` winding budget before order preservation fails;
at ratio 15, core safety forces the least integer floor to be 84.  This is a
scheme optimality statement, not a global optimality claim about all Kakeya
constructions.

## 4. The exact separated-ratio gate Q4

Let `J` be any nonwrapping core-safe interval of length `ell` and put
`K=k4`.  Define

```text
Q4=ell(21K-7sum_i ki)-6Ksum_i(1/ki)-39.                 (13)
```

> **Theorem C.**  If `Q4>0`, the four-comb complement in `J` has a component
> longer than `1/(7K)`.

**Proof.**  Write `A` for the actual survivor mass and `C` for the actual
number of survivor components.  THM-1097's sharp one-comb discrepancy and a
union bound give the lower estimate

```text
A>=A0:=3ell/7-(6/49)sum_i(1/ki).                        (14)
```

A tooth of `D_ki` meeting `J` has its centre in a scaled interval of length
`ki ell+1/7`, so at most `ki ell+8/7` such teeth meet `J`.  Four combs and
the two ends of the connected window therefore give the real upper bound

```text
C<=C0:=ell sum_i ki+39/7.                               (15)
```

Direct expansion yields

```text
7K A0-C0=Q4/7.                                          (16)
```

When `Q4>0`, `7KA>=7KA0>C0>=C`.  In particular `A>0`, so `C>=1`.
Hence the mean length of the `C` actual components is greater than
`1/(7K)`, so one component is longer than that. ∎

Removing the `K` self-terms gives the useful equivalent form

```text
Q4=ell(14K-7(k1+k2+k3))
      -6K(1/k1+1/k2+1/k3)-45.                          (17)
```

Along a ratio ray `m(a,b,c,d)`, every shape with `2d>a+b+c` enters this cone
at an explicit finite scale.  At equality the gate remains a fixed negative
rational, so the slope wall is genuine for this method.

## 5. Corrected exact-transfer cone

The live correction in THM-1137 is essential.  The old `7/6` recursion used
in early THM-1140 notes is false; the exact normalized transfer is

```text
Phi(x)=min(6/7,(x-1/7)/2),       x>=1.                  (18)
```

For a core `P`, let `ell(P)` be its longest safe-component length.  An exact
atlas of all 495 eight-subsets gives

```text
min_P ell(P)(13max(P)+1)=72/35>13/7,                    (19)
```

attained by `P={1,2,6,7,8,9,10,11}` with `ell=1/70`.
Thus the first legal killer puts (18) in its saturated regime.  Define

```text
c1=6/7,
xi=(ki/k_(i-1))c_(i-1),
ci=Phi(xi),                           i=2,3,4.           (20)
```

If each `xi>=1`, the transfers are legal; if the resulting `c4>1/7`, the
desired four-comb component exists.

> **Corollary D.**  Three adjacent ratios `ki/k_(i-1)>=9/5` suffice.

At equality, the exact inputs and outputs are

```text
inputs:   54/35,       63/50,        3519/3500,
outputs:   7/10,      391/700,       3019/7000>1/7.     (21)
```

The optimal common ratio within recurrence (20) is the positive root

```text
r_*=1.797111878...,
6r_*^3-r_*^2-2r_*-28=0.                                (22)
```

`9/5` is a transparent rational corollary and strengthens the corrected
coarse `7/3` cone of THM-1140.

## 6. Exact combined residual

For a fixed core, the currently proved gates are

```text
(i)   k4-k1<=30                                      [THM-1133],
(ii)  the piecewise multiplier gate (12),
(iii) Q4>0,
(iv)  the exact Phi path (20) succeeds.
```

They are not exhaustive.  The first primitive ratio ray left by their
asymptotic predicates is

```text
(k1,k2,k3,k4)=m(3,4,5,6),       m>=53.                  (23)
```

It is legal for every core, has span `3m>30`, misses the multiplier cone,
and lies exactly on the Q4 slope wall:

```text
Q4=-6*6(1/3+1/4+1/5)-45=-366/5.                        (24)
```

The transfer begins with `c1=6/7`, then

```text
x2=(4/3)(6/7)=8/7,       c2=1/2,
x3=(5/4)(1/2)=5/8<1,                                  (25)
```

so it cannot continue.  This is an infinite **proof-method residual**, not a
counterexample to the four-comb theorem.

In normalized coordinates `(x,y,z)=(k1,k2,k3)/k4`, the asymptotic residual
obeys

```text
x<=7/8,       x+y+z>=2,       and the Phi path fails.   (26)
```

The equality `x=7/8` retains an integer-midpoint parity sidecar.  Among all
3,646,069 primitive rays with `1<=a<b<c<d<=100`, the strict multiplier,
strict Q4, and exact-transfer union covers 3,054,412; 591,657 remain
(`16.2273%`).  This census diagnoses ratio space and does not discharge the
unbounded quantifier.

The ray (23), rather than an undifferentiated “clustered majority”, is the
next exact target.  A successful new invariant must cross both a Q4 slope
equality and a failed owner-bearing transfer.

## 7. Verification and tournament audit

Normal and optimized runs are byte-identical.  Frozen SHA-256 values are

```text
artifact 771407d6c072a4438626183ea6c87eaacebba0d21bf1e1267777ec35ad81aaa9
script   cc375186ed1959b7c9f5a1cc5232edd80c6d47cc46d17cee31cc3789846ab17f
output   5759498be1ea709b41959c70cf85cb77832a3a9b150fff65d9bcaea255422181
```

Ordering exact wall events gives a transitive tournament with score
histogram `0,...,N-1`, no directed cycles, singleton SCCs, and one tie-grouped
Hamiltonian path.  It forgets lengths, endpoint owners, overlap clusters,
and wall slopes.

The useful finite carrier instead starts with the cyclic residue word under
the multiplier switch `a -> ua`, labels its Hamiltonian cycle by integer
gaps, and retains the wall-drift and needle-slope budgets.  In interval
transfer, the vertices are metric proof states rather than runners.  The
faithful packet is

```text
(cyclic residue word, gap labels, overlap-cluster count,
 endpoint owner, normalized interval width, wall slope).
```

This challenged vertex choice explains both the negative and positive
results: the naked pair tournament loses (1), while the labelled seven-gap
carrier proves Theorem B.
