---
id: THM-1168
title: Exact component-aware shape erosion closes every primitive proportional four-comb residual of height at most 64
status: PROVED -- finite exact closed-endpoint census of all 95,336 primitive THM-1148 residual shapes (a,b,c,d), d<=64. A universal eroded-radius inequality closes 95,335 shapes; the unique failed inequality (49,50,51,56) is closed at its sole exceptional scale by an exact 165-core component-address bank. Uniform r=5 beyond this finite primitive-shape atlas remains OPEN
source: codex-2026-07-18-S76
depends_on:
  - THM-1148  # exact core floor and the proportional residual predicate
  - THM-1159  # first exact primitive-shape erosion certificate
related:
  - THM-1137  # exact Phi gate used to define the residual census
  - THM-1147  # exact two-comb endpoint law
  - THM-1167  # gap-local route refuted; this theorem retains the component address
script: 04-computation/lrc14_r5_shape_erosion_height64_referee_codex_S76.py
output: 05-knowledge/results/lrc14_r5_shape_erosion_height64_referee_codex_S76.out
---

# THM-1168 -- exact proportional shape erosion through height 64

Write

```text
D_k={t in R/Z : ||kt||<1/14}.
```

The danger comb is open and its safe complement is closed.  This theorem
extends THM-1159 from the first primitive ray `(3,4,5,6)` to every primitive
proportional shape of height at most 64 which remains after THM-1148's three
scale-free gates.  It also records the exact point at which a core-independent
shape radius first stops being sufficient.

Uniform `r=5` for primitive shapes of height greater than 64 remains open.

## 1. The shape invariant

Let

```text
A=(a,b,c,d),                    1<=a<b<c<d,
gcd(a,b,c,d)=1,
S_A={x : ||ax||,||bx||,||cx||,||dx|| >= 1/14},
delta_A=1/(7d).
```

The positive-length components of `S_A` are closed rational intervals.  If
`[u,v]` is such a component, the starts of closed `delta_A`-intervals which
stay in it form

```text
[u,v-delta_A]       when v-u>=delta_A,                (1)
empty               otherwise.
```

Equality in (1) produces a legitimate singleton start.  A safe wall which is
isolated between two unsafe cells is retained in the wall telemetry but
cannot carry a positive interval.  Let `T_A` be the union of all intervals in
(1), let `G_A` be the largest cyclic complementary gap of `T_A`, and put

```text
H_A=G_A+delta_A.                                      (2)
```

> **Shape-needle lemma.**  Every real interval of length strictly greater
> than `H_A` contains a closed subinterval of length `delta_A` whose image
> modulo one lies in `S_A`.

Indeed, an interval of length `L` has an allowed-start interval of length
`L-delta_A`.  When `L>H_A`, this is longer than every cyclic gap of the closed
set `T_A`, so it meets `T_A`.  This is the THM-1159 argument with the primitive
shape left as a parameter.  The strict inequality avoids a closed interval
placed exactly across a largest complementary gap.

## 2. Uniform core-to-shape comparison

Let `P` be an eight-subset of `{1,...,12}`, put `M=max(P)`, and let `ell(P)`
be the length of a longest component of its closed safe set.  THM-1148 proves

```text
ell(P)>=72/[35(13M+1)].                               (3)
```

For a proportional four-comb prefix

```text
(k1,k2,k3,k4)=m(a,b,c,d),                             (4)
```

legality is `am>13M`.  Hence the least possible scale for a fixed `M` is

```text
m_a(M)=floor(13M/a)+1.
```

Define the exact universal phase floor

```text
C_a=min_{8<=M<=12} 72 m_a(M)/[35(13M+1)].             (5)
```

The integer inequality `am>=13M+1` also gives the simpler, slightly weaker
bound

```text
m ell(P)>=72/(35a).                                   (6)
```

> **Universal radius criterion.**  If `H_A<C_a`, every legal scale of the
> ray (4) leaves a closed interval of length
>
> ```text
> delta_A/m=1/(7dm)=1/(7k4)                           (7)
> ```
>
> safe for the core and all four killers.

Choose a longest core-safe component `J`.  Under the phase map `x=mt`, its
lift has length `m ell(P)>=C_a>H_A`.  The shape-needle lemma supplies a closed
`delta_A`-interval in `S_A`; pulling it back gives (7) inside `J`.

This criterion is component-aware: it uses the whole actual core component
as a phase needle.  It does **not** assert that every individual `k1`-gap is
good, the statement refuted in THM-1167.

## 3. Exact residual census

The finite shape census uses exactly THM-1148's proportional residual
predicate.  A primitive shape is retained when all three of

```text
8a>7d,                         [strict multiplier cone]
2d>a+b+c,                      [strict Q4 slope]
exact three-step Phi transfer  [THM-1137]
```

fail.  The scale-dependent span-at-most-30 gate is not needed in this shape
census.  For every retained shape with `d<=64`, the referee reconstructs all
open danger-wall events, closed safe components, eroded start components,
cyclic gaps, `H_A`, and the exact floor (5).

There are exactly

```text
95,336 primitive residual shapes with d<=64.
```

The universal criterion closes `95,335`.  The cumulative counts and the
first rays are:

```text
height cap d       40       45       50       55       56       64
residual count   13548    22259    34386    51233    55000    95336

shape          delta_A      G_A        H_A       C_a-H_A
(3,4,5,6)       1/42      23/168      9/56        21/40
(3,5,6,7)       1/49       5/49       6/49       138/245
(4,5,6,7)       1/49       9/49      10/49        76/245
(3,6,7,8)       1/56       9/56       5/28        71/140
(4,5,7,8)       1/56      71/784     85/784     1591/3920
(4,6,7,8)       1/56       4/49      39/392      813/1960
(5,6,7,8)       1/56      13/112     15/112      111/400.
```

Thus the two next rays after THM-1159, `(3,5,6,7)` and `(4,5,6,7)`, and all
four height-eight residual rays close with large exact margins.  Across the
entire positive bank the tightest row is

```text
A=(45,46,47,52),
H_A=17/364,
C_45=216/4585,
C_45-H_A=97/238420>0.                                 (8)
```

The extension from height 56 through height 64 contains exactly 40,336
additional residual shapes.  All 40,336 satisfy the universal radius
criterion; none needs a new component-address bank.  Its tightest row is

```text
A=(53,59,60,62),
H_A=55/1736,
C_53=48/1225,
C_53-H_A=2279/303800>0.                               (8a)
```

The exhaustive quantifier matters: these are finite exact ray statements,
not sampled evidence about the unbounded ratio space.

## 4. The unique radius obstruction is component-address closed

Exactly one shape fails `H_A<C_a`:

```text
A_*=(49,50,51,56),
delta_*=1/392,
G_*=2/49,
H_*=17/392,
C_49=3/70,
C_49-H_*=-1/1960.                                    (9)
```

Its shape complex has 412 distinct no-tie walls, 124 positive safe
components, 76 eroded start components, and two complementary gaps attaining
`G_*`.  Equation (9) is a failure of the universal-placement criterion, not
a failure of the ray.

The scale split is exact.

- `m<=2` is illegal even for `M=8`, since `49*2<=13*8`.
- At `m=3`, legality is equivalent to `M<=11`; hence the possible cores are
  exactly the `binom(11,8)=165` eight-subsets of `{1,...,11}`.
- For `m>=4`, (3) gives

  ```text
  m ell(P)>=4*72/[35(13*12+1)]
           =288/5495
           =H_*+2783/307720>H_*,                     (10)
  ```

  so the universal shape-needle argument applies again.

It remains only to keep the actual component address at `m=3`.  Erode every
closed core component by `delta_*/3=1/1176`, obtaining the core-start set
`U_P`.  A start `t` works exactly when

```text
t in U_P intersect union_{j=0}^2 (T_*+j)/3.           (11)
```

The referee constructs all intervals in (11) with exact rational endpoints.
The result is

```text
cores checked                         165
cores with no intersection              0
closed component intersections       8456
intersections per core              36..82
minimum best overlap              47/39984
minimum positive overlap          1/119952.            (12)
```

The weakest core is

```text
P=(1,3,4,5,6,7,8,11),
```

with a witness overlap

```text
[1639/2352,13955/19992]
```

of length `47/39984`.  Every one of the 76 primitive shape-start components
participates in at least one core witness.  Choosing any point in an
intersection (11) makes the closed interval of length `1/1176` core-safe,
and multiplication by three maps it to a closed `1/392`-interval in `S_*`.
Thus the four killers `(147,150,153,168)` are safe on the same interval.
Together with (10), this closes every legal scale of the exceptional ray.

## 5. Appending the fifth killer

For any closed interval produced above, its length is exactly

```text
1/(7k4).
```

If `k5>k4`, every connected component of the open danger set `D_k5` has
length `1/(7k5)<1/(7k4)`.  A connected closed interval of the latter length
cannot lie in one such open tooth, so it contains a point safe for `k5` as
well.  Even equality of lengths would be harmless because an open interval
cannot contain a closed interval of the same length.  Hence every ray in the
stated atlas is `1/14`-lonely after an arbitrary ordered fifth killer.

## 6. Relation to the gap-local refutation

THM-1167 exhibits consecutive high-speed quadruples with individual
`k1`-gaps whose longest survivor is below `1/(7k4)`.  Nothing here revives
that false quantifier.

The objects are different:

```text
refuted object:  minimum over every individual k1-gap;
THM-1168 object: full cyclic primitive-shape start complex T_A,
                 intersected with an actual core-component phase needle.
```

For 95,335 shapes the needle is longer than the largest cyclic gap of `T_A`,
which forces an intersection without declaring any one `k1`-gap good.  At
the unique short-needle row, the proof explicitly retains the core address
and checks (11).  This is precisely the information THM-1167 says a positive
extension must preserve.

## 7. Tournament and alternate-carrier audit

Ordering the four runner vertices by speed always gives the same transitive
tournament:

```text
scores=(0,1,2,3), cycles=0, four singleton SCCs,
Hamiltonian paths=1.
```

It cannot distinguish any closed shape from the open frontier.  For wall
vertices, the pairwise observable is chronological position after cutting
the circle at zero.  The gauge orders a simultaneous tie by `(speed, side)`;
this produces a transitive tournament and its tie-resolved Hamiltonian path.
The exceptional shape has 412 distinct no-tie wall vertices, scores
`0,...,411`, no cycles, singleton SCCs, and one chronological Hamiltonian
path.  For shapes with simultaneous events, the `(speed,side)` gauge supplies
the tie Hamiltonian path; naked chronology is a total preorder rather than a
tournament.

The proof-bearing carrier is instead

```text
(wall-owner event word,
 closed safe-component lengths,
 delta-eroded cyclic start complex and weighted gap word,
 core component and phase scale,
 weighted core-start/shape-lift incidence graph).       (13)
```

Candidate vertices explicitly challenged are runners, combs, wall events,
safe components, eroded starts, core components, lifted shape starts, and
incidence edges.  Runner and wall tournaments destroy metric lengths.
Component order without the cyclic gap weights destroys `H_A`.  Primitive
shape starts without core addresses fail exactly at (9).  Core components
without the phase lift destroy (11).  The bipartite incidence object in (13),
not a naked tournament, is the minimal carrier for the exceptional repair.

## 8. Exact replay and honest frontier

The dependency-free referee uses `Fraction` throughout and no `assert`, so
optimized mode checks the same conditions.  It independently audits the
event sweep by direct midpoint tests on every leading shape and the unique
exception.  Normal and `python -O` outputs are byte-identical to the frozen
output.

```text
04-computation/lrc14_r5_shape_erosion_height64_referee_codex_S76.py
SHA-256 253e5e28248a3d5875018fa3b168a7b74c47b45bf86b85e0d9a17e076cfa7d6e

05-knowledge/results/lrc14_r5_shape_erosion_height64_referee_codex_S76.out
SHA-256 97f7ff416e9e815ac9b8dbbefca19debed1e76b9c836b50b8142ee6dbf6fe2eb
```

This theorem closes all 95,336 primitive proportional THM-1148 residual
rays with `d<=64`.  It does not classify any primitive shape with `d>64`,
does not prove a uniform bound for `H_A`, and does not prove uniform `r=5`.
