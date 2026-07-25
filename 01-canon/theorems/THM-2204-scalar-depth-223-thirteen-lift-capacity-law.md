---
id: THM-2204
title: "Scalar depth-(2,2,3) exclusion and thirteen-lift capacity law"
status: >
  PROVED + VERIFIED-EXACT. In the scalar five-unit/three-blocker branch, the
  actual blocker valuation profile (2,2,3) is empty. On the primitive
  13^4 guard-safe annulus, an exhaustive exact audit of all 3,081 unordered
  depth-two blocker pairs and all 13,182 unit sign classes gives a unique
  minimum conditional top-five capacity margin of 1,830. More generally,
  the thirteen coefficient lifts of one unit class obey an exact all-depth
  family-sum law: their capacities sum to 2W-W_a. This scalar recurrence
  does not determine the lifted top-five order statistic or retain sharp
  recursive control. In the
  audited universe one fixed-sum lift family has capacities ranging from
  zero to 2,460; the missing recursive state is its labelled guard-hole
  correlation vector, equivalently its nonconstant Fourier data. The
  THM-2205/2207 subsequently close (1,1,3)/(1,2,3); profiles of deepest
  valuation at least four remain open.
source: codex-2026-07-24-scalar-depth223-lift-capacity
depends_on:
  - THM-2192-scalar-five-plus-three-root-sheet-chord-invoice
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
related:
  - THM-2138-all-depth-unit-annulus-extremal-law
  - THM-2197-scalar-chord-coverage-has-a-boolean-deficiency-quotient
  - THM-2200-convex-semigroup-and-finite-place-support-hole-trichotomy
  - THM-2201-cyclic-root-fibre-hasse-jet-transition-carrier
  - THM-2207-scalar-depth-123-labelled-guard-hole-exclusion
script: 04-computation/lrc14_scalar_depth223_thirteen_lift_capacity_thm2204.py
output: 05-knowledge/results/lrc14_scalar_depth223_thirteen_lift_capacity_thm2204.out
script_sha256: 9c16e1e7a69834f9304c877f1627232f374688aa7f8d2372c260dfdbfa056ac8
output_sha256: 0b89f80f6f82c63bbc6700239f133ff1e4c9ed3b2fe6389bc604f6b8de229bb8
hash_basis: working-tree bytes (LF)
---

# THM-2204 -- scalar depth-(2,2,3) exclusion and thirteen-lift capacity law

Use

```text
D_a={t in R/Z:||at||<1/14},
C_H={t in R/Z:||Ht||>1/7}.                           (1)
```

In the scalar `5+3` branch of THM-2192 and THM-2198, suppose

```text
C_H subset union_(i=1)^5 D_(q_i)
             union union_(j=1)^3 D_(c_j)             (2)
```

almost everywhere.  The coefficients `H,q_1,...,q_5` are positive
thirteen-units, the three actual blockers are positive multiples of
thirteen, and after relabelling their valuations satisfy

```text
1<=lambda_1<=lambda_2<lambda_3,
lambda_j=nu_13(c_j).                                 (3)
```

This theorem has two parts.  First, it proves the exact capacity identity
obeyed by every thirteen-lift family.  Second, it uses the full labelled
capacities, rather than only their family sum, to eliminate

```text
(lambda_1,lambda_2,lambda_3)=(2,2,3).                (4)
```

## 1. The exact thirteen-lift family

Let

```text
N=13Q,                  Q=13^m,
m>=1,
```

and normalize the guard to one.  A primitive image phase is a residue
`r mod Q` with `13` not dividing `r`.  Its thirteen primitive lifts are

```text
z=r+kQ mod N,                    k in F_13.           (5)
```

Label a sheet by

```text
s=-k mod 13.
```

Let `B(r)` be the set of labels for which the lift in (5) belongs to
`C_1`.  Thus `|B(r)|` is nine or ten.

Fix a unit residue `a mod Q`.  Its thirteen coefficient lifts are

```text
q_l=a+lQ mod N,                  l in F_13.           (6)
```

They all have the same inverse

```text
delta=a^(-1) mod 13.
```

For a primitive phase put

```text
I_a(r)=Z intersection
       (ar/Q-13/14,ar/Q+13/14).                       (7)
```

The rooted integer-window law of THM-2198 gives the exact active-sheet set

```text
A_(q_l)(r)
 ={delta(d+lr) mod 13:d in I_a(r)}.                  (8)
```

Indeed

```text
q_l r/Q=ar/Q+lr,
q_l=a mod 13,
```

so lifting the coefficient translates every endpoint by the
phase-dependent sheet displacement `delta*l*r`.

Let `P` be any finite set of primitive image phases; in applications it is
the set on which all fibre-constant blocker bits are off.  Define

```text
C_l(a;P)=sum_(r in P)
          |A_(q_l)(r) intersection B(r)|,            (9)

W(P)=sum_(r in P)|B(r)|,

W_a(P)=sum_(r in P intersection D_a)|B(r)|.          (10)
```

Here `D_a` in (10) is evaluated at `r/Q`.

### Thirteen-lift capacity-sum law

For every `a` and `P`,

```text
sum_(l in F_13) C_l(a;P)=2W(P)-W_a(P).               (11)
```

For a fixed phase `r`, both `delta` and `r` are nonzero modulo thirteen.
As `l` runs through `F_13`, the displacement `delta*l*r` therefore runs
through every sheet exactly once.  Each endpoint in (8) visits every
sheet, and hence visits `B(r)` exactly `|B(r)|` times.  Consequently

```text
sum_l |A_(q_l)(r) intersection B(r)|
       =|I_a(r)| |B(r)|.                             (12)
```

Because `Q` is a power of thirteen, there are no `1/14` torsion
endpoints, and

```text
|I_a(r)|=
 1 if r/Q in D_a,
 2 otherwise
 =2-1_(D_a)(r/Q).                                    (13)
```

Summing (12) and using (13) proves (11).

Equation (11) is an exact recursive equality, not a heuristic average.
It is valid at every depth and for every residual phase set `P`.

## 2. The missing coordinate in the scalar recurrence

Let

```text
F(r)=F_13\B(r),
E_a(P)=sum_(r in P)|I_a(r)|,                         (14)

b_(a,l)(P)=sum_(r in P)
 |{delta(d+lr):d in I_a(r)} intersection F(r)|.      (15)
```

Then

```text
C_l(a;P)=E_a(P)-b_(a,l)(P).                          (16)
```

The family-sum law fixes the sum of the thirteen coordinates in (16), but
the largest lifted capacity is

```text
max_l C_l=E_a-min_l b_(a,l).                         (17)
```

Thus a recursive upper bound for the lifted order statistic requires a
lower bound on the smallest labelled guard-hole correlation.  Neither
`W`, `W_a`, nor their combination in (11) determines that information.
Equivalently, sharp recursive control requires the full thirteen-vector

```text
(b_(a,l))_(l in F_13),
```

or enough of its nonzero Fourier coefficients to control its minimum.

The finite audit below shows that this is a real concentration phenomenon,
not merely a logical possibility.  On the residual associated with the
diagonal shallow pair `(4,4)`, the thirteen sign-class lifts over base
class `a=1098` have `(sign representative,capacity)` rows

```text
((1098, 2460), (1099, 2454), (3295, 2456),
 (3296, 2452), (5492, 2454), (5493, 2452),
 (7689, 2450), (7690, 2456), (9886, 2456),
 (9887, 2450), (12083,2452), (12084,2452),
 (14280,   0)).                                      (18)
```

Their fixed family sum is

```text
29444,
```

but their spread is

```text
2460-0=2460.                                         (19)
```

This is the largest within-family spread in the complete depth-four
profile audit.  In particular, replacing the lift vector by its scalar
sum or average loses the lifted top-five order statistic and any sharp
recursive bound based on it.  The trivial bound by the whole family sum
of course remains valid.

The geometric reason is precise.  For each fixed phase and endpoint,
varying `l` sweeps one affine thirteen-sheet needle.  Equation (11) counts
the total incidences of those needles with the guard-safe sheets.  The
capacity maximum instead detects whether one common lift parameter avoids
the guard holes across many different phases.  The destroyed information
is phase alignment.

## 3. Exact elimination of `(2,2,3)`

Assume (4).  Take the primitive layer

```text
N=13^4=28561,                    Q=13^3=2197.
```

Multiplication of primitive numerators by the guard unit `H` normalizes the
guard to one and replaces every coefficient by its product with `H^(-1)`
modulo `N`.  This bijection preserves all blocker valuations.  Define

```text
U={z mod N:13 does not divide z and 7||z||_N>N},
|U|=18830,                                           (20)
```

where

```text
||z||_N=min(z mod N,-z mod N).
```

The unique deepest blocker is safe everywhere on `U`.  Indeed, writing it
as `13^3w` with `w` a thirteen-unit gives

```text
13^3w*z/13^4=wz/13,
```

a nonzero thirteenth root whose norm is at least
`1/13>1/14`.

The two depth-two blockers have the form `13^2u`.  Their masks depend only
on the unit part `u mod 13^2`, up to sign.  There are

```text
phi(13^2)/2=78                                       (21)
```

such sign classes.  Distinct actual blockers can reduce to the same mask,
so the exhaustive pair universe includes all diagonal pairs:

```text
C(78+1,2)=3081.                                      (22)
```

The unit masks are indexed up to sign modulo `N`; there are

```text
phi(13^4)/2=13182                                    (23)
```

classes.

For a unit label `q mod N` and a depth-two unit part `u mod 13^2`, put

```text
S_q={z in U:14||qz||_N<N},
T_u={z in U:14||13^2uz||_N<N}.                       (24)
```

For each shallow pair define

```text
R_(u,v)=U\(T_u union T_v).                           (25)
```

Order the `13182` conditional capacities

```text
|S_q intersection R_(u,v)|
```

decreasingly as

```text
c_1(u,v)>=...>=c_13182(u,v).                         (26)
```

The exact exhaustive result is

```text
|R_(u,v)|-sum_(i=1)^5 c_i(u,v)>=1830                (27)
```

for all `3081` shallow pairs.  The minimum is unique in the fixed
sign-representative convention.  Its row is

```text
(u,v)=(11,79),
|R_(u,v)|=14484,                                     (28)

((c_i,unit label))_(i=1)^5
 =((2614,1858),(2614,1860),(2504,929),
   (2504,930),(2418,6)).
```

Thus the five-capacity sum is `12654`, leaving the margin

```text
14484-12654=1830.                                    (29)
```

Reduce the five actual unit coefficients modulo `N` and sign.  Repeated
mask classes do not enlarge their union, so discard repetitions.  The
union of the remaining at most five masks inside `R_(u,v)` has size at
most the sum of their individual capacities, and hence at most the five
largest entries in (26).  Equation (27) therefore leaves at least one
residue in `R_(u,v)` safe from all five unit masks.  It is also safe from
the two shallow blockers by (25), from the deepest blocker by the
thirteenth-root argument, and lies in the guard-danger set by (20).

There are no torsion endpoints because a power of thirteen is coprime to
seven and fourteen.  The uncovered residue satisfies all inequalities
strictly and thickens to an uncovered interval.  This contradicts the
almost-everywhere cover (2), proving that (4) is impossible.

## 4. Exact audit

The companion checker uses the root-image disintegration at `N=13^4`.
There are `2028` primitive image phases modulo `13^3`; their guard-safe
root fibres have the exact histogram

```text
1450 fibres of size 9,
 578 fibres of size 10,
1450*9+578*10=18830.                                 (30)
```

A depth-two blocker is constant on each root fibre.  For a unit label `q`,
the checker stores two exact phase bitsets,

```text
{r:h_q(r)>=1},             {r:h_q(r)>=2},
```

where `h_q(r)` is the number of guard-safe, `q`-dangerous roots.  Since a
unit root window has one or two sheets, its capacity on any residual phase
set is the sum of two bit counts.

The checker verifies:

```text
all 78 direct depth-two masks equal their divided image masks;
all 3081 shallow pairs, including all 78 diagonal pairs;
all 13182 unit sign classes on every pair;
the unique minimum row (28);
the full pair-table SHA256 digest
  b82b4f22fe2242197fc8ee587a99acd026120d52001d9093e8b465557cc57024;
an independent direct-torsion reconstruction of the worst residual
  and its five reported capacities in (28);
all 3081*1014=3124134 instances of the family-sum law (11);
the maximum-spread witness (18)--(19).
```

Normal and optimized Python transcripts are byte-identical.  Reproduce
with

```text
python3 04-computation/lrc14_scalar_depth223_thirteen_lift_capacity_thm2204.py
python3 -O 04-computation/lrc14_scalar_depth223_thirteen_lift_capacity_thm2204.py
```

## 5. Connection and scope ledger

The recursive carrier is

```text
source:
  one depth-m unit class and a residual set of primitive image phases;
map:
  lift the coefficient through a+l*13^m for l in F_13;
target:
  a thirteen-vector of labelled sheet-time capacities;
preserved by the family sum:
  total endpoint visits and weighted base danger;
destroyed by scalar averaging:
  phase-aligned guard-hole avoidance and the lifted order statistic;
needed sidecar:
  the labelled hole-correlation vector b_(a,l), or nonzero Fourier data;
cheapest hostile test:
  the spread-2460 family (18).                        (31)
```

The affine needle picture in Section 2 is useful only as an incidence
model.  It is not a tournament: different lift parameters can tie, and the
target is a capacity order statistic rather than a pairwise orientation.
Nor is its scalar family sum a knot invariant: it forgets precisely the
phase alignment needed to reconstruct the lifted movie.

THM-2201 gives a finite faithful coordinate system for the missing data:
phase by phase, the two endpoint sheets and the guard-hole set have
invertible Hasse jets, from which each cyclic correlation summand in
`b_(a,l)` can be reconstructed.  Exact integer accumulation of those
phasewise correlations recovers the thirteen-vector in Section 2.
Augmentation or the scalar family sum keeps only its constant coordinate.
The jet is therefore a faithful carrier for the sidecar, not an inequality
controlling its minimum.  Likewise, THM-2200's independent-product support
pumps classify static matching holes but do not preserve this dynamic
phase alignment.

This theorem closes only the actual valuation profile `(2,2,3)` and does
not establish a uniform all-depth top-five recurrence. THM-2205 and
THM-2207 separately close `(1,1,3)` and `(1,2,3)` using their larger typed
sign-class universes. Profiles with deepest valuation at least four remain
open. QED.
