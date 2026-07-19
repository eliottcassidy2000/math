---
id: THM-1267
title: Centered-spoke survivor protrusion law
status: PROVED.  In every hypothetical six-comb cover of a complete c-safe gap, the full five-comb survivor on the slowest fast speed's centered safe component is forced into one endpoint tail.  Exact six-bin endpoint density and the THM-1241-separated one-comb envelope give d1/c<563/270, hence 270d1<=563c-1.  Paper topology, optimization-safe exact referee, and sorry-free Lean arithmetic consumer are supplied
source: codex-2026-07-19 H-drift / centered-component synthesis
depends_on: [THM-1198, THM-1240, THM-1241]
related: [THM-1199, THM-1236, THM-1244, THM-1252, THM-1274, THM-1283]
script: 04-computation/lrc14_centered_survivor_protrusion_thm1267_referee.py
output: 05-knowledge/results/lrc14_centered_survivor_protrusion_thm1267_referee.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCCenteredSurvivorProtrusion.lean
script_sha256: 8c9c88e52743301b0b7c7e9c64b343a0b87b704100b46f80af574be3004ac43b
output_sha256: 1e2821930e319c545439a6e3ae0e58bb87928dc5f720c434c9e737a6c13f7c68
formalization_sha256: b53309299dc978e4d064f141fd83b48c57dac09e16a6bcdef60cb30c1fbc2ded
---

# THM-1267 -- centered-spoke survivor protrusion law

## 1. Statement

For a positive integer speed `v`, write

```text
D_v={t: ||vt||<1/14},             Dbar_v={t: ||vt||<=1/14}. (1)
```

Let

```text
G=[(14k+1)/(14c),(14k+13)/(14c)]                     (2)
```

be a complete closed safe gap of the positive integer speed `c`.  Suppose

```text
c<d1<d2<d3<d4<d5<d6                                (3)
```

and the six strict danger combs `D_(d1),...,D_(d6)` cover `G`.  Put `d=d1`,

```text
t0=(k+1/2)/c,                 q=c+d.                  (4)
```

Choose either nearest integer `p` to `qt0`, and define

```text
epsilon=p-qt0,       rho=|epsilon|<=1/2,       t=p/q. (5)
```

Let `S` be the complete closed `d`-safe component containing `t`.  Then its
center and its part outside `G` satisfy the exact laws

```text
center(S)=t0+epsilon/d,                              (6)
|S minus G|=max(0,(rho+3/7)/d-3/(7c)).              (7)
```

Normalize length on `S` and set

```text
ell=|S minus G|/|S|.
```

Then

```text
ell=max(0,1/2+7rho/6-d/(2c)).                        (8)
```

The five closed combs belonging to `d2,...,d6` have a THM-1198 survivor on
`S`, and every such survivor lies in the single endpoint tail `S minus G`.
This first gives

```text
ell>1/27,
rho>[27(d/c)-25]/63,
d/c<113/54.                                          (9)
```

Using THM-1241's macroscopic ratio separation for the `d6` load sharpens the
same chain to

```text
ell>11/270,
rho>[135(d/c)-124]/315,
d/c<563/270.                                         (10)
```

Since `c,d` are positive integers, the final strict inequality is equivalently
the useful integer cut

```text
270d1<=563c-1.                                       (11)
```

## 2. Exact center and endpoint-tail geometry

Nearest-integer selection gives

```text
t-t0=epsilon/q,              |epsilon|<=1/2.         (12)
```

Because `q>2c`,

```text
|t-t0|<=1/(2q)<1/(4c)<3/(7c),                        (13)
```

so `t` lies in the interior of `G`.  Moreover

```text
ct=k+1/2+c epsilon/q,
dt=p-ct=p-k-1/2-c epsilon/q.                         (14)
```

The correction in the second line has absolute value strictly below `1/4`.
Consequently

```text
p-k-3/4<dt<p-k-1/4,
floor(dt)=p-k-1.                                     (15)
```

Writing `h=p-k-1`, the complete component through `t` is therefore

```text
S=[(14h+1)/(14d),(14h+13)/(14d)].                    (16)
```

Its center is

```text
(h+1/2)/d=(p-k-1/2)/d=t0+epsilon/d,                  (17)
```

where the last equality uses

```text
p=(c+d)t0+epsilon=k+1/2+dt0+epsilon.                 (18)
```

This proves (6), including its sign.  It also reproves that `t` is interior
to `S`: equations (14)--(15) put its `d`-depth above `1/4`.

The radii of `G` and `S` are respectively

```text
a=3/(7c),                       b=3/(7d),             (19)
```

with `b<a`.  The two intervals meet at the common interior point `t`.
Therefore the shorter interval `S` can cross at most one endpoint of `G`;
`S minus G` is empty or one left/right endpoint interval.  For two such
intersecting intervals its length is

```text
max(0,|center(S)-center(G)|+b-a),                    (20)
```

which is exactly (7).  The sign of `epsilon` chooses the right or left tail.
Finally `|S|=6/(7d)`, and multiplying (7) by `7d/6` gives (8).

This is an exact full-component calculation.  It does not use the smaller
protected intersection `K=G intersect S` from THM-1244 and it loses no
exterior tail.

## 3. Why the five-comb survivor must protrude

Normalize `S` affinely to `[0,1]`.  The five remaining speeds have normalized
slopes

```text
L_j=6d_j/(7d)>6/7,                 j=2,...,6.         (21)
```

Let `U` be the complement in `[0,1]` of their five **closed** normalized
danger combs.  THM-1198 measures `U` with the probability density whose six
successive heights are

```text
f=(3/4,13/12,7/6,7/6,13/12,3/4).                    (22)
```

There is an exact containment

```text
U subset S minus G.                                  (23)
```

Indeed, throughout the closed safe component `S`, the `d=d1` comb is not
strictly dangerous.  If a point of `S intersect G` also avoided all five
closed combs `d2,...,d6`, then none of the six open danger combs would cover
that point, contrary to the assumed cover of `G`.  Using closed combs for
`U` only strengthens this implication.  Endpoint conventions introduce no
exception; for the later integral argument the finitely many interval and
density boundaries are null in any case.

For every `L>6/7`, THM-1198's equality classification gives

```text
P(L)=sup_phi integral_(D(phi,L)) f <7/36.            (24)
```

Thus the union bound and (21) imply

```text
integral_U f>1-5(7/36)=1/36.                         (25)
```

Let `E` denote the normalized endpoint interval corresponding to `S minus G`,
so `|E|=ell`.  There is one necessary case split which must not be suppressed:

- If `ell>=1/6`, then immediately `ell>1/27`.
- If `ell<1/6`, the whole of `E` stays inside the left or right endpoint
  sixth.  By the symmetry in (22), `f=3/4` almost everywhere on either one,
  and hence

  ```text
  1/36<integral_U f<=integral_E f=(3/4)ell.          (26)
  ```

Both branches give

```text
ell>1/27.                                             (27)
```

Because this is positive, the `max` in (8) is on its affine branch.  With
`r=d/c`, equations (8) and (27) rearrange to

```text
rho>(27r-25)/63.                                     (28)
```

Combining (28) with `rho<=1/2` gives `54r<113`, proving every assertion in
(9).  Notice that the endpoint density improves the earlier pure-Lebesgue
constant `89/42` to `113/54` before any Kakeya separation is spent.

## 4. The separated one-comb envelope

THM-1241 proves

```text
d6/d>70/69.                                          (29)
```

Consequently the last of the five loads in (21) obeys

```text
L_6=6d6/(7d)>(6/7)(70/69)=140/161=20/23.             (30)
```

Only this `d6` load receives (30).  The other four loads are bounded by the
global `7/36`; no unproved separation is assigned to them.

The exact THM-1198 envelope through `L=3`, restricted by (30), has the
following twelve branches.  The last column is the maximum relevant on that
branch after imposing `L>20/23`.

| interval | `P(L)` | relevant maximum |
|---|---:|---:|
| `[6/7,68/63]` | `1/(6L)` | open supremum `23/120` at `L downarrow 20/23` |
| `[68/63,8/7]` | `3/4-9/(14L)` | `3/16` |
| `[8/7,6/5]` | `3/(14L)` | `3/16` |
| `[6/5,48/35]` | `5/18-5/(42L)` | `55/288` |
| `[48/35,3/2]` | `11/(42L)` | `55/288` |
| `[3/2,12/7]` | `2/9-1/(14L)` | `13/72` |
| `[12/7,2]` | `13/(42L)` | `13/72` |
| `[2,244/119]` | `1/24+19/(84L)` | `13/84` |
| `[244/119,15/7]` | `3/4-103/(84L)` | `8/45` |
| `[15/7,12/5]` | `8/(21L)` | `8/45` |
| `[12/5,18/7]` | `5/18-2/(7L)` | `1/6` |
| `[18/7,3]` | `3/(7L)` | `1/6` |

The first formula is decreasing and

```text
1/[6(20/23)]=23/120=161/840.                         (31)
```

The threshold is excluded, so the value is not attained there.  Every later
candidate in the table is strictly smaller; the closest is

```text
23/120-55/288=1/1440.                                (32)
```

For `L>=3`, THM-1198's periodic-primitive bounded-variation estimate gives

```text
P(L)<=1/7+1/(7L)<=4/21<23/120,                       (33)
23/120-4/21=1/840.
```

Approaching `20/23` from above inside the first exact branch proves the
reverse supremum inequality.  Hence the exact all-frequency statement is

```text
sup_(L>20/23) P(L)=23/120,                            (34)
```

and every individual load on the open domain is strictly below this value.

Apply (34) to `d6` and the global THM-1198 bound to `d2,...,d5`.  Their total
weighted load is strictly below

```text
4(7/36)+23/120=349/360.                              (35)
```

Therefore

```text
integral_U f>11/360.                                 (36)
```

Repeat the endpoint split exactly as in Section 3.  If `ell>=1/6`, then
`ell>11/270`.  If `ell<1/6`, equations (23), (22), and (36) give

```text
11/360<integral_U f<=(3/4)ell,
ell>11/270.                                          (37)
```

Substitution into (8) yields

```text
rho>[135(d/c)-124]/315.                              (38)
```

Finally `rho<=1/2` turns (38) into

```text
270(d/c)<563,                                        (39)
```

which proves (10)--(11).

## 5. Functional `H`-drift interpretation

The gain in this theorem is functional, not another pair-overlap invoice.
THM-1198 centers the tooth indicator by `1/7` and takes its periodic primitive
`H`; its bounded amplitude and the variation of `f` control the large-`L`
drift in (33).  The exact finite arrangement supplies the same load as the
piecewise function `P(L)`.  THM-1267 retains that dependence on `L` instead
of replacing it everywhere by the scalar `7/36`.

The complete conversion chain is

```text
THM-1241 ratio drift
 -> L6>20/23
 -> one-comb load defect
 -> survivor f-mass
 -> endpoint-tail length
 -> centered error rho
 -> slowest-fast ratio d1/c.                         (40)
```

Numerically, the separated load cap gains exactly

```text
7/36-23/120=1/360                                   (41)
```

over the universal bound.  The endpoint density `3/4` converts this into
exactly one extra `1/270` of normalized protrusion:

```text
1/27=10/270            becomes            11/270.    (42)
```

This is the positional face of the `H`-drift: a frequency-side load defect
must appear as a literal endpoint displacement of the centered safe
component.  It is complementary to THM-1252's seam/toothpick detuning.  No
gcd quantum, Hunter forest, or pair edge is charged here, so these credits
remain available downstream.

The conversion is also recursively suggestive.  A later theorem may repeat
the same operation on a protruding component with a newly separated top load.
THM-1267 proves only the single centered rung; it does not assert that such an
iteration closes or is self-similar without further phase control.

## 6. Alternate carriers and tournament audit

We explicitly challenged runners, gaps, fixed sixths, section boundaries,
tooth wall-crossing events, residues, cover arcs, Fourier modes, pair circuits,
and proof obligations as vertices.

On the five noncarrier runners `d2,...,d6`, use normalized-slope comparison as
the pairwise observable and increasing speed as the switch/gauge.  The result
is the transitive tournament with score histogram

```text
(0,1,2,3,4),
```

zero directed cycles, five singleton SCCs, and the unique tie Hamiltonian path

```text
d2,d3,d4,d5,d6.                                      (43)
```

This quotient preserves rank and identifies `d6` as the one load to which
THM-1241 applies.  It destroys the five phases, actual values of `P(L_j)`, the
orientation of the endpoint tail, and the containment `U subset S minus G`.
It is therefore a useful loss audit but not the proof object.

Using the six density bins as vertices preserves the endpoint weight `3/4`
but forgets which comb pays each load.  Using wall-crossing events preserves
the finite arrangement behind `P` but forgets survivor position.  The faithful
carrier is instead the hybrid collection

```text
(one oriented endpoint-tail obligation;
 five phase-bearing weighted-load obligations).      (44)
```

Its preserved predicate is precisely: the common five-comb survivor has
enough `f`-mass and must lie in that tail.  The challenged assumption is that
a pairwise tournament should carry every LRC reduction.  Here the decisive
object is a one-body functional budget coupled to interval topology.

## 7. Paper topology provider, formalization, and referee

The paper proof provides the topological/measure layer:

1. nearest-integer selection and the floor identification (15);
2. identification of `S` as the complete closed safe component;
3. the one-endpoint form and length of `S minus G`;
4. the exact containment (23);
5. transport of THM-1198's density to `S`, the union bound, and the two
   `ell>=1/6` / `ell<1/6` integrations.

The Lean arithmetic consumer is sorry-free.  It proves the center identity,
the half-integer strip, physical-to-normalized protrusion including the
`max 0` branch, both endpoint quantiles and their explicit case splits, the
basic and separated five-load budgets, all relevant envelope-candidate gaps,
the strict first-piece comparison, both exact `rho` implications, both ratio
caps, rational cross-multiplication, and the integer rounding in (11).  Its
`#print axioms` audit contains only standard quotient/classical foundations;
there are no theorem-specific axioms and no `native_decide`.

The standard-library exact referee uses `Fraction` and an explicit `require`
rather than Python assertions.  It checks:

- `371,752` signed-gap centered geometry rows, including tie choices and both
  contained and protruding branches;
- the exact mirror positive control `(c,d,k,p)=(1,2,0,1 or 2)`, where the two
  tails have opposite orientation and `ell=1/12`;
- normalization, symmetry, endpoint masses, and both inverse quantiles of the
  six-bin density;
- the partition, continuity, monotonic endpoint maximum, and exact maximum of
  every one of the twelve `P(L)` pieces;
- the BV tail, the open supremum (34), every displayed constant, and the
  integer strict-rounding consequence;
- the alternate-carrier tournament loss audit.

Ordinary and optimized Python outputs are byte-identical.  Frozen hashes are

```text
referee  8c9c88e52743301b0b7c7e9c64b343a0b87b704100b46f80af574be3004ac43b
output   1e2821930e319c545439a6e3ae0e58bb87928dc5f720c434c9e737a6c13f7c68
Lean     b53309299dc978e4d064f141fd83b48c57dac09e16a6bcdef60cb30c1fbc2ded
```

THM-1267 does **not** prove six-comb noncoverage, the sporadic branch empty,
or LRC(14).  It proves a strict new structural condition on every hypothetical
six-cover: its slowest fast speed must satisfy (11), while its centered
five-comb survivor occupies a quantitatively nontrivial exterior endpoint
tail.  Any closing argument must combine this positioned tail with an
independent contradiction, invoice, or recursive carrier obstruction.

THM-1274 makes the first such composition.  In a slowest-rooted blocker
two-cycle, the protrusion sign fixes the reverse descent digit, and the
blocker wall facing that endpoint continues through only the other five owner
labels.  It therefore yields a closest return with factor greater than `3/2`
or reaches the carrier-gap endpoint within five tooth occurrences.  The
exterior survivor itself remains outside the local six-cover word; converting
that terminal landing into a neighboring-carrier descent or endpoint tax is
performed in THM-1283.  The terminal endpoint owner crosses the adjacent
carrier tooth in a proper exterior seam, and subtracting that tooth segment
from the compulsory survivor tail gives the strict residue/gcd refinement
`270d_1x+45d_1Q<=563cx-1`.  Following the resulting mixed-owner child remains
open.
