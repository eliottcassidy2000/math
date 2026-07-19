---
id: THM-1198
title: Universal five-comb noncoverage from a six-bin dual density
status: PROVED / COMPUTER-EXACT.  A six-step rational probability density has closed-tooth load at most 7/36 for every real slope L>=6/7 and every phase.  Therefore five arbitrary faster danger combs leave dual mass at least 1/36 and normalized Lebesgue measure at least 1/42, hence physical length at least 1/(49c) in every c-slow gap.  Consequently every six-comb cover has d_1<13c/6 and six disjoint private regions of length at least 1/(49c), collapsing the THM-1176 toothpick ladder at its first rung
source: codex-2026-07-18-S76 five-comb universal-density agent
depends_on: []
related: [THM-1094, THM-1176, THM-1178]
script: 04-computation/lrc14_five_comb_universal_six_bin_density_codex_20260718.py
output: 05-knowledge/results/lrc14_five_comb_universal_six_bin_density_codex_20260718.out
lean: 04-computation/lean/TournamentH7/TournamentH7/LRCFiveCombDual.lean
---

# THM-1198 -- universal five-comb six-bin dual density

At radius `1/14`, write

```text
D(phi,L)={x in [0,1]: ||phi+Lx||<=1/14}.             (1)
```

The closed convention in (1) makes the result stronger than the usual open
danger-comb statement and does not change any integral.  Define the symmetric
step function

```text
       3/4,    0<=x<1/6 or 5/6<x<=1,
f(x)= 13/12,   1/6<=x<1/3 or 2/3<x<=5/6,             (2)
       7/6,    1/3<=x<=2/3.
```

Its six bin heights are

```text
(3/4,13/12,7/6,7/6,13/12,3/4),                       (3)
```

and

```text
integral_0^1 f(x)dx=1.                                (4)
```

The central result is a phase-free operator norm.

> **Theorem 1 (six-bin dual density).**  For every real `L>=6/7` and every
> real phase `phi`,
>
> ```text
> A(L,phi):=integral_(D(phi,L)) f(x)dx <=7/36<1/5.    (5)
> ```
>
> The constant is exact:
>
> ```text
> sup_(L>=6/7,phi mod 1) A(L,phi)=7/36.               (6)
> ```
>
> Equality occurs only at `L=6/7`; there it occurs precisely for
> `phi mod 1 in [1/2,9/14]`.

The proof has a small exact arrangement on the bounded-frequency piece and a
one-line bounded-variation tail.  Both parts are given below so the finite
certificate has an explicit analytic boundary.

## 1. Exact compact arrangement

First suppose

```text
6/7<=L<=3,        y=1/L,        z=phi/L,              (7)
```

where the phase is reduced to `0<=phi<=1`.  Thus

```text
1/3<=y<=7/6,                 0<=z<=y.                 (8)
```

The two endpoints of the tooth centred at the integer `n` are

```text
x=(n+-1/14)y-z.                                      (9)
```

On this compact parameter domain, `phi+Lx` lies in `[0,4]`; hence only
`n=0,1,2,3,4` can meet `[0,1]`.  Formula (2) changes only at the seven bin
boundaries `j/6`.  Draw the seventy rational event lines

```text
(n+-1/14)y-z=j/6,
n=0,...,4,  j=0,...,6,                               (10)
```

together with the four boundary lines in (8).

On every cell of this line arrangement, the order of every tooth endpoint
relative to every bin boundary is fixed.  Each tooth/bin overlap length is
therefore the difference of two fixed affine endpoint expressions, or zero.
Consequently `A` is continuous and affine on each cell.  Its maximum on the
compact polygon is attained at an arrangement vertex.

The dependency-free `Fraction` referee intersects all pairs of the `74`
lines.  Among `2,701` raw line pairs it finds exactly `101` distinct feasible
vertices.  At every vertex it evaluates the integral in two independent
forms:

1. direct tooth-by-six-bin intersection;
2. the superlevel identity

   ```text
   f=3/4+(1/3)1_[1/6,5/6]+(1/12)1_[1/3,2/3].         (11)
   ```

The two exact evaluations agree at every vertex.  The largest values are

```text
7/36, 55/288, 3/16, 13/72, 5/28, 8/45, 11/63,...    (12)
```

and the only maximizing vertices are

```text
(y,z)=(7/6,7/12), (7/6,3/4).                         (13)
```

The segment between them lies on the boundary `y=7/6` and is an arrangement
edge on which the affine load is constant.  There is no maximizing vertex
off that boundary.  Thus the compact maximum is `7/36`, with equality only
when `L=6/7`.  Dividing the two endpoint `z` values by `y=7/6` gives the phase
interval in Theorem 1.

There is also a direct equality check.  At `L=6/7`, a complete danger
preimage has length `1/6`.  Its load reaches `7/36` exactly when it lies inside
the height-`7/6` plateau `[1/3,2/3]`.  For the tooth centred at `1`, its left
endpoint is

```text
13/12-(7/6)phi.                                       (14)
```

Requiring (14) to lie in `[1/3,1/2]` is exactly
`phi in [1/2,9/14]`.

### 1.1 The exact phase-free load envelope through `L=3`

The same arrangement yields more than its global maximum.  Define

```text
P(L)=sup_(phi mod 1) A(L,phi).                         (14a)
```

For fixed `y`, a maximum in `z` occurs at `z=0`, `z=y`, or on one of the
event lines (10), because the load is piecewise affine in `z`.  Intersecting
those phase lines with the full arrangement gives `175` affine load
segments.  Their exact upper envelope has twelve pieces:

```text
L interval             P(L)
[6/7,68/63]             1/(6L)
[68/63,8/7]             3/4-9/(14L)
[8/7,6/5]               3/(14L)
[6/5,48/35]             5/18-5/(42L)
[48/35,3/2]             11/(42L)
[3/2,12/7]              2/9-1/(14L)
[12/7,2]                13/(42L)
[2,244/119]             1/24+19/(84L)
[244/119,15/7]          3/4-103/(84L)
[15/7,12/5]             8/(21L)
[12/5,18/7]             5/18-2/(7L)
[18/7,3]                3/(7L).                       (14b)
```

The non-arrangement breakpoint `244/119` is the exact crossing of two affine
phase candidates; the referee inserts every such segment intersection before
taking the upper envelope.  Thus (14b) is an exact vertical-slice theorem,
not a phase grid.  For example it gives

```text
P(1)=1/6, P(8/7)=3/16, P(13/7)=1/6,
P(2)=13/84, P(15/7)=8/45.                             (14c)
```

The last value is a useful guardrail: `P(L)` is not monotone, and in
particular `sup_(L>=13/7)P(L)=1/6` is false already at `L=15/7`.

## 2. Analytic BV tail

Put

```text
h(u)=1_{||u||<=1/14}-1/7.                             (15)
```

On the period `[-1/2,1/2]`, start a primitive `H` at zero.  It has slope
`-1/7` outside `[-1/14,1/14]` and slope `6/7` inside.  Its two extrema are
`-3/49` and `3/49`, so it closes periodically and

```text
||H||_infinity=3/49.                                  (16)
```

The interior variation of (2) is

```text
TV_(0,1)(f)=2(1/3+1/12)=5/6.                         (17)
```

The zero extension also pays the two endpoint jumps `3/4`, hence

```text
TV_R(f 1_[0,1])=3/4+5/6+3/4=7/3.                    (18)
```

Stieltjes integration by parts, using
`d H(phi+Lx)=L h(phi+Lx)dx`, now gives

```text
|A(L,phi)-1/7|
 <= ||H||_infinity TV_R(f 1_[0,1])/L
 =1/(7L).                                             (19)
```

For `L>=3`, therefore,

```text
A(L,phi)<=1/7+1/(7L)<=4/21<7/36.                     (20)
```

Equations (7)--(20) prove Theorem 1 and its exact equality statement.  Notice
that the tail cutoff has a genuine rational margin:

```text
7/36-4/21=1/252.                                      (21)
```

## 3. Five arbitrary phases cannot cover

> **Theorem 2 (universal real five-comb noncoverage).**  Let
> `L_1,...,L_5>=6/7` be arbitrary real slopes and let
> `phi_1,...,phi_5` be arbitrary phases.  Then the five closed sets
>
> ```text
> D(phi_i,L_i), i=1,...,5,                            (22)
> ```
>
> do not cover `[0,1]`.  More precisely, for
>
> ```text
> U=[0,1] minus union_i D(phi_i,L_i),                 (22a)
> ```
>
> one has
>
> ```text
> integral_U f>=1/36,             |U|>=1/42.          (22b)
> ```

**Proof.**  If they covered, then their indicator sum would be at least one
on `[0,1]`.  Integrating against the nonnegative probability density `f` and
using Theorem 1 would give

```text
1<=sum_i A(L_i,phi_i)<=5(7/36)=35/36<1,              (23)
```

a contradiction.  More quantitatively, the union bound gives

```text
integral_U f=1-integral_(union_i D_i)f
 >=1-sum_i A(L_i,phi_i)>=1/36.                        (23a)
```

Since `f<=7/6`,

```text
1/36<=integral_U f<=(7/6)|U|,                        (23b)
```

which gives `|U|>=1/42`.  This argument permits repeated slopes, independent
arbitrary phases, real rather than integer frequencies, and closed rather
than open teeth. ∎

## 4. Exact slow-gap normalization

Let `c>0`, and let

```text
G=[t_0,t_0+6/(7c)]                                    (24)
```

be any complete safe gap of the speed-`c` danger comb.  Write

```text
t=t_0+(6/(7c))x.                                      (25)
```

For a frequency `d>=c`,

```text
||dt||<=1/14
 iff ||phi+(6d/(7c))x||<=1/14,                       (26)
phi=dt_0 mod 1.
```

Thus its normalized slope is

```text
L=6d/(7c)>=6/7.                                       (27)
```

The phases arising from integer slow gaps are constrained rationals, but
Theorem 2 already allows every real phase.  This checks the arbitrary-phase
relaxation in the safe direction: it enlarges the family being ruled out and
loses no original counterexample.

> **Corollary 3 (uniform `r=5` closure).**  No five frequencies
> `d_i>=c`, integer or real, have danger combs covering a complete `c`-slow
> gap.  Their common survivor inside that gap has length at least
>
> ```text
> (6/(7c))(1/42)=1/(49c).                             (27a)
> ```

This closes the uniform `r=5` slow-gap problem left open in
THM-1094/1176.  It is not a finite carrier statement or a consecutive-shape
statement; it is uniform in the scale, phase, and five independent ratios.

## 5. The six-comb first-tooth law

> **Corollary 4.**  Suppose
>
> ```text
> c<d_1<...<d_6                                      (28)
> ```
>
> and the six faster danger combs cover a complete `c`-slow gap.  Then
>
> ```text
> d_1<13c/6.                                          (29)
> ```

**Proof.**  An interval of length at least

```text
1/d_1+6/(7d_1)=13/(7d_1)                             (30)
```

contains a complete `d_1`-slow gap: choose the first slow-gap left endpoint
after the interval's left endpoint, paying less than one `d_1` period, and
then append the gap length.  If `d_1>=13c/6`, then

```text
|G|=6/(7c)>=13/(7d_1).                                (31)
```

Hence `G` contains a complete `d_1`-slow gap.  With the closed-tooth
convention, the `d_1` danger comb is absent on the **interior** of that gap
but may contain its two endpoints.  The other five closed danger combs cover
the interior.  Their finite union is closed, so it also covers the closure,
contrary to Corollary 3.  Equivalently, the two endpoints are null for the
dual-mass proof. ∎

For integer speeds (29) has the exact form

```text
6d_1<=13c-1.                                          (32)
```

THM-1176's former finite toothpick ladder said that a six-comb cover must
satisfy at least one of

```text
d_1/c<13/6,  d_2/d_1<13/6,  d_3/d_2<4/3.             (33)
```

Corollary 4 proves that the first alternative always holds.  In the H-drift
coordinates `x_i=c/d_i`, every six-comb cover therefore begins with

```text
x_1>6/13,                                             (34)
```

in addition to THM-1176's harmonic pressure
`sum_i x_i>1` and THM-1178's rational seam surplus.

> **Corollary 5 (six private needles).**  Under the six-cover hypotheses of
> Corollary 4, define the unique-provider region
>
> ```text
> Q_i=G minus union_(j!=i) D_(d_j).                   (34a)
> ```
>
> Then the six `Q_i` are pairwise disjoint, `Q_i` is contained in `D_(d_i)`,
> and
>
> ```text
> |Q_i|>=1/(49c)                         for every i. (34b)
> ```
>
> In normalized coordinates each also has dual mass at least `1/36`.

**Proof.**  Delete comb `i`.  The other five have normalized slopes at least
`6/7`, so Theorem 2 gives dual survivor mass at least `1/36` and Corollary 3
gives physical survivor length at least `1/(49c)`.  Because all six combs
cover, that survivor lies in comb `i`; hence it is the unique-provider region
(34a), up to the irrelevant endpoint convention.  A point cannot have two
different unique providers, so the six regions are disjoint. ∎

There is a second immediate structural form.  Let

```text
C(x)=#{i:x in D(phi_i,L_i)}                            (34c)
```

for a normalized six-cover.  Theorem 1 gives

```text
integral f(x)(C(x)-1)dx
 =sum_i A(L_i,phi_i)-1<=6(7/36)-1=1/6.               (34d)
```

Therefore the region covered at least twice has `f`-mass at most `1/6`, and
the region with a unique provider has `f`-mass at least `5/6`.  Thus every
putative six-cover is forced to be an almost-partition under the same dual
measure, not merely a harmonic near-equality.  Since `f>=3/4`, the normalized
Lebesgue measure of the multiply-covered region is at most

```text
(1/6)/(3/4)=2/9.                                      (34e)
```

Hence the union of the six private regions has normalized length at least
`7/9`, or physical length at least

```text
(6/(7c))(7/9)=2/(3c).                                 (34f)
```

There is also an exact functional drift.  Define the phase-free majorant

```text
Pbar(L)=P(L),                    6/7<=L<=3,
       =1/7+1/(7L),             L>3,                  (34g)
```

where the compact part is (14b) and the tail is (19).  If six combs cover,
integrating their indicator sum against `f` forces

```text
sum_(i=1)^6 Pbar(6d_i/(7c))>=1.                       (34h)
```

Unlike a scalar harmonic inequality, (34h) retains the oscillatory tooth
geometry at twelve exact ratio intervals.  The tail branch reads
`Pbar(6d/(7c))=1/7+c/(6d)`.  The majorant jumps upward at the chosen analytic
cutoff `L=3`; this is harmless and honest--(14b) is the exact envelope only
through `3`, while (19) is a uniform rather than exact tail.

## 6. Carrier and Tournament Analysis audit

The protected slow gap is a one-dimensional Kakeya needle, but the primary
dual object is not a runner configuration.  It is the probability measure
`f dx` whose mass cannot be captured by five arbitrary periodic tooth
preimages.  The compact proof vertices are the endpoint/bin crossing events
(10).

A runner-order tournament is not meaningful for the proof inequality.  There
is no pairwise observable: Theorem 1 is a universal **one-comb operator
norm**, and Theorem 2 sums five scalar bounds.  Sorting runners gives a
transitive tournament, but it discards the phase, the clipped tooth endpoints,
and their weighted overlap lengths.  Thus no switch/gauge, cycles, SCCs, or
Hamiltonian-path count can sharpen (5).  This records the challenged
assumption explicitly: the faithful vertices in the core proof are
arrangement events and dual-mass obligations, not runners or arcs.

Corollary 5 does create a meaningful derived tournament carrier.  Take the
six private needles `Q_i` as vertices, use the pairwise observable

```text
w_i-w_j,                  w_i=integral_(Q_i)f,         (35)
```

and switch it to an orientation by sign.  Resolve ties by `inf Q_i` in the
chronological coordinate on the gap, then by index (the relatively open
`Q_i` need not have a first point).  The resulting
tournament is transitive: score histogram `0,1,2,3,4,5`, no directed cycles,
six singleton SCCs, and one tie-resolved Hamiltonian path.  This ranking keeps
the six lower bounds `w_i>=1/36` but destroys the shapes and locations of the
private regions, so the faithful carrier is the labelled six-region family
itself.  The first-tooth law (29), harmonic pressure, seam surplus, Fano/gcd
periods, and beat redundancy can now be attached as obligations to those
vertices; none should replace the phase-bearing arrangement carrier used in
the proof of (5).

## 7. Reproducibility and scope

The exact script has no nonstandard dependencies.  Normal and optimized
Python runs are byte-identical to the stored output.  Frozen SHA-256 hashes
are

```text
source  4a498cc295ed7a1f3511ca889beecca07e8565c7bef5492ac46e14dc9ec3c73f
output  38cfd0e710bc9eb8d5dfa2d4d5ef1aefee005d419b59f6a307a93a14925fd772
```

The imported Lean module `LRCFiveCombDual` kernel-checks the six-bin mass,
variation and margin arithmetic, the abstract five-load contradiction, the
`1/42`, `1/(49c)`, `2/9`, `7/9`, and `2/(3c)` conversions, the six-load
overlap-surplus consumer, the integer first-tooth form (32), the `L=15/7`
nonmonotonicity guardrail, and the abstract phase-free functional drift.  The
arrangement maximum and BV integration-by-parts inequality remain explicitly
analytic inputs.  Targeted `lake build TournamentH7.LRCFiveCombDual` succeeds
with no `sorry` or `native_decide`.

This theorem completely closes universal five-comb coverage of a slow gap,
gives the scale-free survivor floor `1/(49c)`, and gives the first-rung plus
private-needle laws for any six-comb cover.  It does **not** by itself exclude
all six-comb covers, prove the remaining crown extraction needed for LRC(14),
or turn the scale-sensitive seam surplus into a uniform contradiction.  Those
are the honest residuals.
