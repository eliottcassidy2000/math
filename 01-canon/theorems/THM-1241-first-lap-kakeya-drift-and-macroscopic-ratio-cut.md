---
id: THM-1241
title: First-lap Kakeya drift forces a macroscopic cut in every six-comb cover
status: PROVED analytic theorem / COMPUTER-EXACT constant and one-comb-envelope audit.  Every pivot d_h>=7c/6 in a six-comb cover pays sum_i|d_i-d_h|>d_h/14.  In particular d_6/c>7/6 and sum_i(d_6-d_i)>d_6/14, hence d_6/d_1>70/69 and some adjacent ratio d_(j+1)/d_j>211/210.  A full six-speed tangent stalk with coalescing normalized ratios cannot occur
source: codex-2026-07-19 tangent-stalk session
depends_on: []
related: [THM-1176, THM-1198, THM-1199, THM-1232, THM-1233]
script: 04-computation/lrc14_first_lap_kakeya_drift_thm1241.py
output: 05-knowledge/results/lrc14_first_lap_kakeya_drift_thm1241.out
lean: 04-computation/lean/TournamentH7/TournamentH7/LRCFirstLapKakeyaDrift.lean
formalization_sha256: 2213a9430262b1820bf40c21f9b8cf31171ecfba2fc62dd47e385ec9f06f997a
---

# THM-1241 -- the six-comb ratio packet cannot coalesce

Let

```text
D_s={t in R: ||st||<1/14}
```

and let `c` be a positive integer and `k` an integer.  Put

```text
G=[(14k+1)/(14c),(14k+13)/(14c)]                    (1)
```

and assume it is a complete closed safe gap of speed `c`.  Suppose that the
six distinct faster positive integer speeds

```text
c<d_1<d_2<d_3<d_4<d_5<d_6                           (2)
```

have strict danger combs covering `G`.  Every pivot at or above one full
normalized lap satisfies

```text
d_h>=7c/6  implies  sum_(i=1)^6 |d_i-d_h|>d_h/14.    (2a)
```

There is always at least one such pivot, and in particular

```text
d_6/c>7/6,                                            (3)
sum_(i=1)^6 (d_6-d_i)>d_6/14.                         (4)
```

Consequently

```text
(d_6-d_1)/d_6>1/70,          d_6/d_1>70/69,          (5)
```

and at least one edge of the ordered speed path obeys

```text
d_(j+1)/d_j>211/210.                                  (6)
```

This is a scale-free separation theorem inside THM-1232/1228's compact
ratio box.  It rules out the proposed full-cluster tangent degeneration:
bounded offsets, or even `d_6-d_1=o(d_6)`, cannot lift to six-comb covers.

## 1. Slow-gap normalization and the exact one-comb load

Write

```text
t=(14k+1)/(14c)+(6/(7c))x,             0<=x<=1,
L_i=6d_i/(7c),
phi_i=d_i(14k+1)/(14c) mod 1.                         (7)
```

The normalized `i`th comb is

```text
U_i={x in [0,1]: ||phi_i+L_i x||<1/14}.               (8)
```

For later reference, the exact phase-free Lebesgue-load envelope is

```text
Q(L)=sup_phi |{x in [0,1]:||phi+Lx||<1/14}|
    =[floor(L)/7+min(frac(L),1/7)]/L.                 (9)
```

Indeed, after `u=phi+Lx`, an interval of length
`L=m+s`, `m=floor(L)`, contains `m` complete periods of the danger arc and
one residual interval of length `s`.  The full periods contribute `m/7`,
and the residual intersects a circular arc of length `1/7` in at most
`min(s,1/7)`.  Translation attains that upper bound, proving (9).

Suppose first that `L_6<=1`.  Since every `d_i>c`, one has

```text
6/7<L_i<=1.                                           (10)
```

An interval of length at most one meets the periodic danger set in total
length at most `1/7`; equivalently, (9) gives

```text
|U_i|<=1/(7L_i)<1/6.                                  (11)
```

Thus `sum_i |U_i|<1`, so their union cannot cover `[0,1]`.  This proves

```text
L_6>1,                                                (12)
```

which is exactly (3).

## 2. Freeze one complete pivot lap

Fix any index `h` with `L_h>=1`.  The initial interval

```text
0<=x<=1/L_h                                           (13)
```

lies inside `[0,1]`.  Put `u=L_hx`, so `u` traverses one complete copy of
the circle.  Define the signed relative drifts

```text
eta_i=L_i/L_h-1=d_i/d_h-1.                            (14)
```

On this first pivot lap the `i`th danger condition is

```text
||phi_i+u+eta_i u||<1/14.                             (15)
```

The circle distance is 1-Lipschitz.  Hence every point satisfying (15)
also satisfies

```text
||phi_i+u||<1/14+|eta_i|.                             (16)
```

In other words, freeze the six arcs at the start of the lap and enlarge the
`i`th radius by exactly its total relative drift during that lap.  If the
moving arcs (15) cover the lap, the six static open arcs (16) cover the
circle.

If some radius in (16) is at least `1/2`, then

```text
|eta_i|>=1/2-1/14=3/7>1/14,                          (17)
```

and (2a) follows immediately.  Otherwise every arc is a proper circle arc,
with exact length

```text
min(1,2(1/14+|eta_i|))=1/7+2|eta_i|.                 (18)
```

The `min(1,2r)` form records the large-radius branch rather than silently
using the short-arc formula outside its range: at radius `1/2` the strict
ball has full measure but misses its antipode, while larger radii give the
whole circle.

A finite family of nonempty proper **open** arcs covering the connected
circle (take representatives `0<=u<1`) has total length strictly greater
than one.  Measure gives at least
one.  If equality held, all pairwise overlaps would have measure zero; any
nonempty overlap of open arcs has positive measure, so the arcs would be
pairwise disjoint.  A disjoint union of at least two nonempty open sets
cannot be the connected circle.  Therefore (18) gives

```text
1<sum_i(1/7+2|eta_i|)
 =6/7+2sum_i |eta_i|,
sum_i |eta_i|>1/14.                                   (19)
```

Substituting (14) into (19) is precisely the all-pivot law (2a).  Equation
(12) permits the fastest pivot `h=6`, for which every `eta_i<=0`; this
specialization is (4).

The proof uses arbitrary real phases.  It therefore retains every arithmetic
phase orbit as a special case, without assuming equidistribution, a stable
residue word, or bounded endpoint-owner complexity.

## 3. The macroscopic path cut

Since `d_i>=d_1`,

```text
sum_i(d_6-d_i)<=5(d_6-d_1).                           (20)
```

Equations (4) and (20) give (5).  There is a sharper way to expose where
the separation sits.  Put

```text
g_j=d_(j+1)-d_j,                       1<=j<=5.        (21)
```

Telescoping with multiplicity gives

```text
sum_(i=1)^6(d_6-d_i)=sum_(j=1)^5 j g_j.              (22)
```

If every `g_j<=d_6/210`, then the right side of (22) is at most

```text
(1+2+3+4+5)d_6/210=d_6/14,                           (23)
```

contrary to (4).  Thus some `g_j>d_6/210`; since `d_j<d_6`,

```text
d_(j+1)/d_j=1+g_j/d_j>1+1/210=211/210,               (24)
```

which proves (6).

There is a useful carrier-scale refinement.  Let

```text
r=#{i:d_i>=7c/6},
M=max_(1<=j<=5)(d_(j+1)-d_j).                         (24a)
```

The qualifying pivots form a fastest suffix.  For pivot `h`, its absolute
drift invoice has gap-weight sum

```text
W_6=1+2+3+4+5=15,
W_5=1+2+3+4+1=11,
W_4=1+2+3+2+1=9.                                    (24b)
```

Since `d_h/14>=c/12`, the all-pivot law gives the strict branch ledger

```text
r=1    => M>c/180,
r=2    => M>c/132,
r>=3   => M>c/108.                                  (24c)
```

Here `r>=1` follows from (3).  These are single-pivot consequences, not a
claim that the joint system of all pivot invoices has been optimized.  They
turn the qualitative five-way cut into a progressively larger beat-scale
edge as more speeds cross one full normalized lap.

## 4. Compactness consequence and adjusted target

Let a sequence of putative covers have carriers tending to infinity and
pass to any convergent ratio subsequence

```text
x_i=lim d_i/c.                                        (25)
```

Equations (3)--(6) survive weakly in the limit.  In particular,

```text
x_6>=7/6,
x_6-x_1>=x_6/70>=1/60,                               (26)
```

and after passing once more to a fixed adjacent index, some

```text
x_(j+1)-x_j>=x_6/210>=1/180.                          (27)
```

Thus the remaining compact object is not one six-fold tangent stalk.  It
has at least one macroscopic cut, leaving at most five-speed tangent
subclusters on its two sides.  First-order offsets remain relevant inside
those smaller clusters, but they cannot carry the whole obstruction.

This adjusts the post-THM-1233 target: stratify by the cut index in (24),
then attach the arithmetic phase/offset stalk only within each side.  The
five possible cut indices are a finite branch set.  The all-pivot law (2a)
adds another dispersion invoice whenever an earlier ratio crosses `7/6`,
without creating a new phase case.

## 5. Kakeya and Tournament Analysis audit

The protected slow gap is the Kakeya needle.  On its first fastest lap the
faithful vertices are not runners but the six moving danger arcs

```text
(static centre phi_i, radius 1/14, drift eta_i).      (28)
```

This carrier preserves the cover predicate on (13).  Freezing the motion
destroys its exact wall chronology but retains a safe outer approximation;
the lost motion is paid exactly by the radius increments `|eta_i|`.
Raw wall counts would be unfaithful because arbitrarily fine teeth can add
walls without increasing the relative first-lap drift.

For the required pairwise tournament, take observable `d_j-d_i`, use strict
speed order as the switch/gauge, and use the ordered speeds as the tie
Hamiltonian path.  The tournament is transitive, has score histogram
`0,1,2,3,4,5`, no directed cycles, six singleton SCCs, and unique path

```text
d_1 -> d_2 -> d_3 -> d_4 -> d_5 -> d_6.             (29)
```

The new information is the edge-weighted path invoice (22).  It forces one
path edge to be macroscopic, whereas the unweighted tournament sees only
the same total order in every instance.  We challenged runners, teeth,
fixed circle sections, section boundaries, wall-crossing events, residues,
cover arcs, Fourier modes, overlap components, and proof obligations as
vertices.  Moving first-lap arcs preserve the local Kakeya predicate;
weighted path edges are the smallest quotient that preserves the resulting
dispersion consequence.

## 6. Exact replay and honest frontier

The companion dependency-free `Fraction` referee independently reconstructs
the one-comb phase arrangement, checks (9) on a rational bank, derives every
constant in (2a)--(27), and prints the diagonal alias windows that the scalar
load envelope alone would leave.  Normal and optimized runs are required to
be byte-identical.  The Lean companion kernel-checks the six-arc arithmetic
invoice, the weighted path identity, the diameter and ratio conversions,
and the five-way additive edge disjunction; it keeps the open-circle cover
lemma as the explicit paper input.  The targeted Lean build is green with no
proof placeholders.  Frozen SHA-256 hashes are

```text
source  795fd008f4c146746171b768ac5b9eae887e69e8a66bf7042e1d23e939064338
output  e983ba1e187f3bb437f5a8f04cdaf635872247fceb8923343c0a16984098dcf4
```

THM-1241 excludes all full-cluster coalescing limits and creates a finite
five-way macroscopic-cut stratification.  It does not exclude separated
ratio packets, classify the phase stalks inside the resulting subclusters,
extend the finite-carrier BV bank past its current frontier, prove universal
six-comb noncoverage, or prove LRC(14).
