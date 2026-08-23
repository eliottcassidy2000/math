---
id: THM-3786
title: "Quadratic etale surface irregular two-by-three Euler-support no-go"
status: >
  PROVED + VERIFIED-EXACT + PENDING INDEPENDENT HOSTILE AUDIT.  On the
  THM-3783 quadratic etale surface, every Darboux pair with two Euler
  weights in one output and three in the other must have common-step
  arithmetic-progression supports.  All irregular two-by-three cells are
  empty: the contribution-bucket census leaves a unique commuting hub, and
  its sign synchronization, constant-profile, and weight-zero boundaries
  reduce to a forbidden homogeneous bracket, the proved two-by-two no-go,
  or an identically zero scalar bucket.  The common-step AP lane remains
  open except for THM-3783's already-proved aligned orientation.
source: jc_sparse_direct_search / contribution-bucket hub peel, 2026-08-23
audit: >
  SELF-AUDITED PROOF CANDIDATE.  The pair-sum census is exact for all
  positive gaps.  The companion independently checks 4,096 bounded gap
  triples, a necessary constant/nonconstant commutation atlas over weights
  [-64,64] and gaps [1,18], every surviving zero-weight family, the
  negative-arm exit, the upper synchronized identity, and the central
  reduction to equal-gap two-by-two support.  Normal, optimized, and frozen
  outputs agree byte for byte.  Independent hostile audit remains due,
  especially for the all-degree hub-lemma case split.
depends_on:
  - THM-3783-quadratic-tower-etale-surface-maximal-polynomial-observable
related:
  - THM-3781-common-step-three-by-three-danielewski-darboux-nonentry
script: 04-computation/jc2_quadratic_etale_surface_irregular_two_by_three_thm3786.py
output: 05-knowledge/results/jc2_quadratic_etale_surface_irregular_two_by_three_thm3786.out
script_sha256: 4fcc6e552fc79a252018555e9ae8b341cdec117006adcefda97795be4c2d0be2
output_sha256: 5461057220838202cd7a92f5a947efbb7cbcb65f0fde8dc0de502564cc886764
semantic_sha256: 917995306a060c8d6c27869d549ceb6cb6fa8538c72b5374ada5aa6636a2493c
hash_basis: raw LF bytes
---

# THM-3786 -- irregular two-by-three supports have a commuting hub

**PROVED + VERIFIED-EXACT + PENDING INDEPENDENT HOSTILE AUDIT.**  The
complete two-by-two no-go in THM-3783 leaves two-by-three Euler support as
the first sparse construction lane on the smooth surface

```text
S=k[r,z,g]/(r^3g-z^2+r/4).                            (1)
```

This theorem closes every irregular placement at once.  Only supports on
one common-step arithmetic progression survive the contribution-bucket
gate.  This is a structural reduction, not a bounded search.

Let `F,G in S` have exact active Euler-weight supports

```text
supp(F)={a,a+d},
supp(G)={b,b+e,b+e+f},             d,e,f>0.           (2)
```

If

```text
{F,G}=1,                                               (3)
```

then necessarily

```text
d=e=f.                                                 (4)
```

Consequently all irregular `2 by 3` cells, and by output swap all irregular
`3 by 2` cells, are empty.  No assertion is made here about every
common-step cell.

## 1. Inherited homogeneous calculus

Retain THM-3783's source chart

```text
s=1+2x^4y,
r=x^2,                 z=xs/2,       g=(s^2-1)/(4x^4).               (5)
```

A nonzero homogeneous element of weight `u` is uniquely `x^u p(s)`, where

```text
p(-s)=(-1)^u p(s),
(s^2-1)^ceil(-u/4) divides p(s) when u<0.              (6)
```

Thus a negative-weight profile is nonconstant and vanishes at `s=+/-1`.
A constant profile is possible only in an even nonnegative weight; at
weight zero it is a global scalar.

For profiles `p,q`, put

```text
W_(u,v)(p,q)=u p q'-v p' q.                            (7)
```

Then

```text
{x^u p(s),x^v q(s)}=2x^(u+v+3)W_(u,v)(p,q).           (8)
```

THM-3783 proves that no single homogeneous bracket is a nonzero scalar and
that no Darboux pair has at most two active weights in each output.

## 2. Exact collision census

Subtract `a+b` from the six pair sums in `(2)`.  The addresses are

```text
0, e, e+f, d, d+e, d+e+f.                            (9)
```

The bottom and top are unique.  Positivity leaves exactly three possible
equalities among the other addresses:

```text
d=e,                 d=f,                 d=e+f.      (10)
```

If the scalar occupied a singleton address, that address alone would be a
forbidden homogeneous scalar bracket.  Hence a collision is necessary.  If
two relations in `(10)` hold, positivity forces

```text
d=e=f.                                                (11)
```

Thus a non-common-step cell either has no collision and is dead, or has
exactly one type

```text
L: d=e, f!=d;          U: d=f, e!=d;          C: d=e+f.                (12)
```

In each type one `G` component occurs only in singleton buckets with both
components of `F`.  It must commute with both and will be called the hub.

## 3. Hub lemma, including constants and weight zero

Suppose components of weights `u,c` have profiles `p,h` and commute.  Then

```text
u p h'-c p'h=0.                                       (13)
```

If `h` is nonconstant and `c!=0`, either `u=0` and `p` is constant, or `p`
is nonconstant and `uc>0`.  Indeed, at every root,

```text
u ord_alpha(h)=c ord_alpha(p),                        (14)
```

and positive orders force equal signs.  If `c=0` and `h` is nonconstant,
`(13)` forces `u=0`; two distinct spoke weights cannot both satisfy this.

If `h` is constant, regularity gives `c>=0` even.  For `c>0`, equation
`(13)` forces `p` constant.  For `c=0`, the hub is a global scalar and may
be deleted by target translation.  Therefore a nonconstant hub synchronizes
all nonzero spoke signs, with a zero spoke forced constant; a positive
constant hub forces both spoke profiles constant.  These alternatives are
exhaustive and use no division by a profile.

## 4. Lower collision

Assume `d=e`, `f!=d`, put `A=a+d`, and normalize the scalar bucket.  Then

```text
b=-3-A.                                               (15)
```

Its two addresses are `(F_a,G_(b+d))` and `(F_A,G_b)`.  The top `G`
component is the hub, of weight

```text
c=b+d+f=-3-a+f.                                      (16)
```

The four singleton equations are

```text
{F_a,G_b}=0,              {F_A,G_(b+d)}=0,
{F_a,G_c}=0,              {F_A,G_c}=0.                (17)
```

For a nonconstant hub, two positive spokes contradict the fixed negative
weight `b=-3-A`.  Two negative spokes have `A<=-1`, hence

```text
c>=d+f-2>=0,                                          (18)
```

with equality only at the excluded common-step case `d=f=1`.  Straddling
signs are impossible.  The boundaries `a=0` and `A=0` respectively meet
the fixed negative weight `b+d=-3` or force `d+f<=2`, again only the
excluded case.

For a constant hub, positive `c` forces both `F` profiles constant; the
negative singleton partners in `(17)` then force both spoke weights zero,
impossible.  Hence `c=0`, so `f=a+3` and `a>=-2`.  The first two equations
in `(17)` force `a<=0`, `A<=0`, with a zero profile constant.  The complete
integer list is

```text
(a,d,f)=(-2,1,1), (-2,2,1), (-1,1,2).                (19)
```

The first is common-step.  In either remaining cell `A=0`; its scalar
summand vanishes, while the other scalar summand has two negative profiles
and vanishes at `s=1`.  Type `L` is empty.

## 5. Upper collision

Assume `d=f`, `e!=d`, put `A=a+d`, and normalize

```text
b=-3-A-e.                                             (20)
```

The bottom component `G_b` is the hub.  The scalar addresses are

```text
(F_a,G_(-3-a))              and              (F_A,G_(-3-A)),          (21)
```

and the singleton equations are

```text
{F_a,G_b}=0,              {F_A,G_b}=0,
{F_a,G_(-3-A)}=0,         {F_A,G_(-3-a)}=0.           (22)
```

A positive nonconstant spoke cannot share the negative hub.  For two
negative spokes, the last two equations in `(22)` leave only

```text
a=-3, d=1 or 2;                   or                  a=-2, d=1.       (23)
```

At `a=-3`, the top `G` weight is zero and its singleton equation makes it a
global constant.  Its scalar summand vanishes, leaving a forbidden
homogeneous scalar bracket.

At `a=-2,d=1`, the spoke weights are `-2,-1`.  If their profiles are `p,q`
and the hub profile is `h`, the two hub equations give

```text
p'/p=2q'/q,                         p=lambda q^2.      (24)
```

The other singleton equations make the weight `-2` and `-1` profiles of
`G` proportional to `p` and `q`.  Both summands in `(21)` therefore commute
separately, and the scalar bucket is zero.

The remaining zero-spoke boundary is `A=0`.  That spoke is a global
constant, so one scalar summand vanishes and the other is a forbidden
homogeneous bracket.  A weight-zero constant hub would force
`A=-3-e<0` while a singleton partner has weight `e>0`; a positive constant
hub contradicts `(20)`.  Type `U` is empty.

## 6. Central collision

Assume `d=e+f`, put `A=a+d`, and normalize

```text
b=-3-A.                                               (25)
```

The middle `G` component is the hub, of weight

```text
c=b+e=-3-a-f.                                        (26)
```

The scalar addresses are `(F_a,G_(-3-a))` and `(F_A,G_(-3-A))`; the
singleton equations are

```text
{F_a,G_(-3-A)}=0,         {F_a,G_c}=0,
{F_A,G_c}=0,              {F_A,G_(-3-a)}=0.           (27)
```

For a nonconstant hub, the sign lemma and endpoint equations leave only:

1. `A=0`, when one scalar summand vanishes and the other is a forbidden
   homogeneous bracket; or
2. `a=-3,d=2` (so `e=f=1`), when the top `G` component is a global constant
   and the same reduction applies.

There is no remaining nonzero-sign regime: two negative spokes would need
`-3-a<0`, but `a>-3` and `d=e+f>=2` make `A>=0`.

For a constant hub, positive weight contradicts the spoke signs.  At weight
zero, `(26)` gives

```text
a=-f-3,                         A=e-3.                (28)
```

If `e<3`, the last singleton in `(27)` has opposite nonzero weights.  For
`e>=3`, delete the global constant hub.  The remaining supports are

```text
F: {-f-3,e-3},                    G: {-e,f}.           (29)
```

Both gaps equal `e+f`, and both middle sums equal `-3`.  This is literally
the complete two-by-two cell excluded by THM-3783.  Type `C` is empty.

## 7. Consequence and controls

The no-collision cells and the three one-collision types exhaust every
irregular support, proving `(4)`.  Output swap changes only the bracket sign,
so the same conclusion holds for `3 by 2` support.

The exact companion verifies all buckets for `1<=d,e,f<=16`; keeps the
common-step double collision as a positive open-boundary control; runs a
necessary constant/nonconstant sign atlas over weights `[-64,64]` and gaps
`[1,18]`; and checks the lower arm exit, upper synchronization, and central
equal-gap reduction exactly.  The finite atlas is a hostile control, not the
universality proof; the all-degree proof is Sections 2--6.

Common-step AP cells remain open here, apart from THM-3783's named aligned
orientation.  **QED, conditional only on independent hostile audit.**
