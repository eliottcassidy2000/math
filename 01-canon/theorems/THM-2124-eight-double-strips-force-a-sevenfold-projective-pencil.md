---
id: THM-2124
title: "Eight double strips over F_13 force a sevenfold projective pencil"
status: >
  PROVED. For every odd prime power q>=7, the minimal (q+1)/2 double strips
  can cover F_q^2 only in one projective pencil, while (q+3)/2 double strips
  have normal-direction partition (k) or (k-1,1). Thus eight double strips
  over F_13 force seven normals to be projectively equal. In the rank-two LRC
  torus model, a guard divisible by 13 and eight terminal characters nonzero
  modulo 13 can cover almost everywhere only if at least seven terminal
  residues lie in one pencil. In the (7,1) branch, the aligned seven already
  cover; under THM-2097's LRC specialization hypotheses they are all
  rationally guard-proportional, and THM-2128 excludes that branch by an exact
  mod-169 seven-comb descent. THM-2131's digit-section lift excludes the
  all-aligned (8) branch under the same rank-two LRC hypotheses. Thus the
  guard-blocker/no-terminal-blocker lane is empty, but LRC(14) remains open.
source: codex-2026-07-22-LRC-F13-double-strip-pencil
depends_on:
  - THM-2114
related:
  - THM-2120
  - THM-2122
  - THM-2123
  - THM-2125
  - THM-2128
  - THM-2131
---

# THM-2124 -- eight double strips force a sevenfold pencil

## 1. Finite-plane near-pencil theorem

Let `q>=7` be an odd prime power and put

```text
k=(q+3)/2.                                             (1)
```

For `i=1,...,k`, let `ell_i` be a nonzero linear form on `F_q^2`, let
`A_i` be a subset of `F_q` of size at most two, and set

```text
S_i={v in F_q^2:ell_i(v) in A_i}.                      (2)
```

If

```text
F_q^2=union_i S_i,                                    (3)
```

then the multiplicity partition of the projective normals `[ell_i]` is

```text
(k) or (k-1,1).                                       (4)
```

In particular, eight unions of at most two parallel affine lines cover
`F_13^2` only if seven of their eight normal directions are equal.

### The vanishing polynomial

Padding a set `A_i` can only enlarge `S_i`, so choose distinct
`a_i,b_i in F_q` with `A_i subset {a_i,b_i}`. The nonzero polynomial

```text
P(x,y)=product_i (ell_i(x,y)-a_i)(ell_i(x,y)-b_i)      (5)
```

vanishes at every point of `F_q^2` and has total degree

```text
2k=q+3.                                                (6)
```

We use the following elementary bounded-degree form of the finite-field
vanishing ideal. If a polynomial `F` of total degree at most `q+3` vanishes
on `F_q^2`, then

```text
F=U(x,y)(x^q-x)+V(x,y)(y^q-y),   deg U,deg V<=3.       (7)
```

Indeed, divide first by `x^q-x` and then by `y^q-y`. The final remainder has
degree less than `q` in each variable and still vanishes on `F_q^2`. Fixing
one coordinate and using the one-variable root bound twice makes that
remainder zero. The total-degree division also gives the bounds in (7).

Let

```text
R(x,y)=product_i ell_i(x,y).                           (8)
```

Taking the homogeneous part of degree `q+3` in (7) gives cubic homogeneous
forms `A,B` such that

```text
R(x,y)^2=x^q A(x,y)+y^q B(x,y).                       (9)
```

The shape on the right is invariant under `GL_2(F_q)`: after a linear
coordinate change, Frobenius gives

```text
(a x+b y)^q=a x^q+b y^q.                              (10)
```

Thus (9) may be used in coordinates adapted to any projective root of `R`.

### A repeated root is almost the whole pencil

Let a projective root of `R` have multiplicity `m`, move it to `[0:1]`, and
write

```text
r(s)=R(s,1),  a(s)=A(s,1),  b(s)=B(s,1).              (11)
```

Then

```text
r(s)^2=s^q a(s)+b(s),              deg a,deg b<=3.     (12)
```

If `m>=2`, the left side is divisible by `s^4`. Reducing (12) modulo `s^q`
shows that the degree-at-most-three polynomial `b` is zero. Hence `r^2` is
divisible by `s^q`, so

```text
2m>=q,
m>=(q+1)/2=k-1.                                       (13)
```

Thus every projective root of `R` is either simple or has multiplicity at
least `k-1`.

### The all-simple case is impossible

Suppose all `k` roots are simple. Since `P^1(F_q)` has `q+1>k` points, a
linear coordinate change can place a nonroot at infinity. All `k` roots of
`r(s)=R(s,1)` then belong to `F_q`, and `r` is a nonzero polynomial of degree
`k`.

For `s in F_q`, equation (12) and `s^q=s` give

```text
r(s)^2=s a(s)+b(s).                                   (14)
```

The right side is represented by a polynomial of degree at most four. It
vanishes at the `k>=5` distinct roots of `r`, so it is the zero polynomial.
Equation (14) would then make `r` vanish at all `q` field elements. But
`deg r=k<q`, a contradiction.

There is therefore a repeated root. By (13) it has multiplicity at least
`k-1`; since `deg R=k`, the only possible partitions are exactly (4). QED.

### The sharp minimal layer

The same argument also explains why seven is the critical number at `q=13`.
Put

```text
r=(q+1)/2.                                             (15)
```

First, fewer than `r` double strips cannot cover by point count: each has at
most `2q` points, while

```text
(r-1)2q=q(q-1)<q^2.                                   (16)
```

If `r` double strips do cover, their padded vanishing polynomial has degree
`q+1`. Equation (7) now has sidecars of degree at most one, and its top part
is

```text
R(x,y)^2=x^q A_1(x,y)+y^q B_1(x,y).                   (17)
```

At any projective root of multiplicity `m>=1`, adapted coordinates give

```text
r_0(s)^2=s^q a_1(s)+b_1(s),       deg a_1,deg b_1<=1. (18)
```

The left side is divisible by `s^2`, so `b_1=0`. Hence `2m>=q`, and oddness
of `q` gives `m>=(q+1)/2=r`. This one root consumes the full degree of `R`.
Thus a minimal cover has direction partition `(r)`. In particular,

```text
seven double strips cover F_13^2 only if all seven are parallel. (19)
```

## 2. Guard-blocker corollary for LRC(14)

Let `Gamma` be a rank-two character lattice and

```text
K=Hom(Gamma,R/Z).
```

Let the guard and eight terminal characters be

```text
g=13u,                 c_1,...,c_8 in Gamma,           (20)
```

where `u` is nonzero and every `c_i mod 13` is nonzero. Suppose that, outside
a null set, the terminal danger bands cover the guard-safe region:

```text
{X:||g.X||>1/7}
 subset union_i {X:||c_i.X||<1/14}.                   (21)
```

Write `T=K[13]`, so `T` is a two-dimensional vector space over `F_13`. Let
`E` be the null exceptional set in (21), enlarged by any boundary sets, and
put

```text
E^#=union_(z in T)(E-z).                              (22)
```

This is still null. The guard-safe set has positive measure, so choose
`X` in it but outside `E^#`. For every `z in T`,

```text
g.(X+z)=13u.X=g.X,                                    (23)
```

and `X+z` is outside `E`. Hence all 169 points of the root fiber `X+T` must
be covered by the eight terminal bands.

For fixed `i`, the map

```text
z |-> c_i.z,       T -> (1/13 Z)/Z                    (24)
```

is a nonzero `F_13`-linear form because `c_i mod 13` is nonzero. The thirteen
possible values of `c_i.(X+z)` are equally spaced on the circle. An arc of
length `1/7<2/13` contains at most two of them. Therefore

```text
{z in T:||c_i.(X+z)||<1/14}                           (25)
```

is a union of at most two parallel affine lines in `T`, with projective
normal `[c_i mod 13]`.

The eight sets (25) cover `T`. Applying the `q=13` case of Section 1 proves

```text
there is a one-dimensional subspace L <= Gamma/13 Gamma
such that #{i:c_i mod 13 belongs to L}>=7.             (26)
```

Thus THM-2123's surviving guard-`13`-blocker branch is not arbitrary: when no
terminal is itself a `13`-blocker, it is forced into an eightfold or
seven-plus-one terminal pencil modulo thirteen.

### The seven-plus-one branch is high-vertical

Suppose the direction partition is `(7,1)`, and call the seven aligned
terminals `c_1,...,c_7`. On any clean guard-safe fiber, their seven strips
already cover all of `T`. Indeed, if they missed one affine line in their
common direction, the exceptional transverse double strip would meet that
thirteen-point line in at most two points and could not finish the cover.
Thus, outside a null set,

```text
{X:||g.X||>1/7}
 subset union_(i=1)^7 {X:||c_i.X||<1/14}.             (27)
```

Now assume the inherited LRC specialization hypotheses of THM-2097: the
characters admit a primitive specialization on which the guard value is
positive odd and the seven terminal values are positive and pairwise
distinct. If `g,c_1,...,c_7` spanned rational rank two, THM-2097 would give a
strict torus point outside all seven dangers but inside the guard-safe set,
contradicting (27). Therefore

```text
c_i in Q g,                    i=1,...,7.              (28)
```

If the guard and all eight terminals span the ambient rational rank two, the
exceptional terminal is consequently the sole rationally transverse
terminal. THM-2128 subsequently proves that this seven-high-vertical-plus-one-
transverse handoff is impossible.

## 3. Exact scope and the next obstruction

The theorem does **not** say that either pencil pattern covers the torus.
It records a necessary condition on one finite root fiber. In particular:

1. affine line offsets were used to obtain the cover but disappear from the
   final projective conclusion;
2. congruence modulo thirteen forgets integral height, content above one
   factor of thirteen, and the drift of danger residues as `13X` moves;
3. if a terminal also vanishes modulo thirteen, its fiber danger set is not a
   double strip and this corollary does not apply; those branches are routed
   through THM-2120, THM-2122, THM-2123, and the complementary THM-2125 gate;
4. ranks nine through twelve do not have the critical identity
   `2k=q+3`, so the same binary-form argument does not give (4) unchanged.

The theorem itself leaves a dynamic question: can eight mod-`13`-parallel
integral characters keep choosing covering line offsets on almost every
guard-safe root fiber? THM-2131 answers no by retaining the affine offsets on
a periodic next-digit section and lifting the pencil through every power of
thirteen. Together with THM-2128, this empties both patterns in the guard-
blocker/no-terminal-blocker lane; simultaneous terminal blockers are outside
the present finite-plane model.

## 4. Kakeya duality, assumption challenge, and Tournament Analysis

The useful vertex set here is not the runners. It can be taken to be the
eight doubled affine needles, their projective normal directions, the 169
root-fiber points, or the eight linear factors of `R`. Passing from needles
to factors preserves the union-cover predicate exactly through the vanishing
of (5). Passing onward to projective directions destroys the affine offsets,
but the Frobenius gap in (9) recovers the near-pencil multiplicity.

The intrinsic pairwise observable is only

```text
det(ell_i,ell_j)=0 or not.                             (29)
```

It is a symmetric parallel/transverse relation, not a canonical tournament.
One can normalize representatives and orient nonparallel pairs by the
quadratic character of their determinant, breaking parallel ties along the
label Hamiltonian path `1,...,8`. But rescaling one representative by a
nonsquare flips all incident edges. That switch preserves neither the cover
predicate nor the affine offsets, so its score histogram, cycles, SCCs, edge
flips, and Hamiltonian-path count are gauge artifacts and are rejected here.
Equation (4) says only that the intrinsic tie graph contains `K_8` or `K_7`.

This challenges the assumption that a pair-overlap graph or tournament can
detect the critical finite-plane obstruction. The obstruction is genuinely
global: it is the degree-`q+3` product of all doubled needles inside the
vanishing ideal of the entire plane. The tournament becomes useful again only
after it is enriched by the moving affine offsets, which are precisely the
data needed to attack the surviving sevenfold-pencil drift.
