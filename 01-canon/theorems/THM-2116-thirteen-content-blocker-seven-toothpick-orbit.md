---
id: THM-2116
title: "One 13-content blocker forces an exact seven-toothpick orbit cover"
status: >
  PROVED as a structural reduction. In a rank-eight two-torus cover, assume
  the guard is nonzero modulo 13, exactly one terminal is 13 times a character
  independent of the guard, and each of the other seven terminals is
  transverse to the guard modulo 13. Every guard-kernel orbit of order 13 on
  which the guard and blocker are safe must, under the stated pointwise cover,
  be covered by seven terminal subsets of sizes one or two. Hence either six
  two-point toothpicks plus one singleton partition the orbit, or seven
  toothpicks cover it with exactly one doubled point. For the almost-everywhere
  cover supplied by THM-2098, the same dichotomy holds for almost every base
  orbit. Two singleton phases on a pointwise-covered safe orbit give an
  immediate escape. This is the first exact consumer of THM-2114's mandatory
  13-content blocker. THM-2120 subsequently excludes both extremal patterns
  throughout this branch; neither theorem closes every rank-eight branch or
  proves LRC(14).
source: codex-2026-07-22-LRC-thirteen-blocker-toothpick-orbit
depends_on:
  - THM-2098
  - THM-2114
related:
  - THM-2095
  - THM-2105
  - THM-2120
---

# THM-2116 -- the first 13-blocker becomes seven finite toothpicks

Let `Gamma` be the saturated rank-two character lattice, put

```text
Gamma^vee=Hom(Gamma,Z),              K=Gamma^vee tensor (R/Z),
```

and regard `K` as the corresponding connected two-torus.  Let

```text
g,c_*,c_1,...,c_7 in Gamma                              (1)
```

be a guard and eight terminal characters.  Assume

```text
c_*=13u in Gamma,
g,u are Q-independent,
g mod 13 is nonzero,
c_i mod 13 is not proportional to g mod 13
                                  for i=1,...,7.         (2)
```

The last condition says that the seven residual terminals are transverse to
the guard over `F_13`.  Suppose their terminal-danger bands cover the guard
complement:

```text
{X:||g.X||>1/7}
 subset union_(i in {*,1,...,7}) {X:||c_i.X||<1/14}.    (3)
```

## 1. The guard-kernel needle

Choose a nonzero vector

```text
v in Gamma^vee/13 Gamma^vee
```

which spans the kernel of the pairing map

```text
Gamma^vee/13 Gamma^vee -> F_13,       v |-> <g,v>,      (3a)
```

and choose an integral lift, still called `v`.  For any `X_0 in K`, form the
order-thirteen orbit

```text
O(X_0)={X_0+jv/13:j in F_13}.                           (4)
```

The points are distinct because `v mod 13` is nonzero.  Moreover

```text
g.(X_0+jv/13)=g.X_0,
c_*.(X_0+jv/13)=c_*.X_0+u.v*j=c_*.X_0                  (5)
```

in `R/Z`.  Thus the guard and the 13-content blocker are constant on this
finite Kakeya needle.

There are safe base points.  Independence of `g,u` makes the integer torus
homomorphism

```text
X |->(g.X,u.X)
```

full-rank and hence surjective, even when its determinant is not `+/-1`.  For
example, choose

```text
g.X_0=1/2,                     u.X_0=1/52.              (6)
```

Then `||g.X_0||=1/2` and `||c_*.X_0||=1/4`.  More
generally, the following argument applies to every base point satisfying

```text
||g.X_0||>1/7,                 ||c_*.X_0||>1/14.        (7)
```

## 2. Every residual terminal is a singleton or a toothpick

Fix `i in {1,...,7}` and put

```text
r_i=c_i.v mod 13.                                      (8)
```

Transversality in (2) says `r_i!=0`.  Consequently the thirteen values of
`c_i` on (4) are a permutation of a translated uniform grid:

```text
{c_i.X_0+j/13:j in F_13}.                              (9)
```

An open circle interval of radius `1/14` has length `1/7`, strictly between
one and two spacings of this grid:

```text
1/13 < 1/7 < 2/13.                                    (10)
```

It therefore contains either one or two of the points in (9).  Define

```text
D_i(X_0)={j in F_13:
          ||c_i.(X_0+jv/13)||<1/14}.                   (11)
```

Then

```text
|D_i(X_0)| in {1,2}.                                   (12)
```

The two-point case is a colored toothpick: in the `c_i`-coordinate the two
points are adjacent, while in the common orbit coordinate their difference
is `+/-r_i^(-1) mod 13`.  The step label must be retained; the seven
toothpicks need not be edges of one common 13-cycle.

There is also an exact phase test.  If
`s_i=13(c_i.X_0) mod 1`, the roots in (9) solve `13z=s_i`.  Direct endpoint
counting gives

```text
|D_i(X_0)|=1  iff  ||s_i||<=1/14,
|D_i(X_0)|=2  iff  ||s_i||> 1/14.                       (13)
```

At equality the second grid point lies on the strict terminal boundary and
is excluded, which is why the first line is closed.

## 3. The two extremal orbit patterns

Take a base point satisfying (7).  Equations (5) make both `g` and `c_*`
safe on all of `O(X_0)`.  If (3) holds, the remaining seven danger sets must
cover all thirteen orbit labels:

```text
F_13=union_(i=1)^7 D_i(X_0).                            (14)
```

Let `s` be the number of singleton sets in (12).  Their total incidence is

```text
sum_i |D_i|=s+2(7-s)=14-s.                              (15)
```

Covering thirteen points forces `s<=1`.  This yields the exact dichotomy.

1. If `s=1`, (15) equals thirteen.  Every orbit point has multiplicity one:
   six toothpicks and the singleton form a disjoint partition of `F_13`.
2. If `s=0`, (15) equals fourteen.  Exactly one orbit point has multiplicity
   two and the other twelve have multiplicity one.  In particular there is no
   triple point.

Thus

```text
two singleton phases on one base point satisfying (7)
  imply that (3) is false.                              (16)
```

Equation (13) turns (16) into a concrete simultaneous-phase target:

```text
find X_0 with guard and blocker safe and
||13c_i.X_0||<=1/14, ||13c_j.X_0||<=1/14
for some i!=j.                                         (17)
```

This is a two-window problem on the original torus, not an unbounded speed
search.

There is one measure-theoretic distinction to retain. THM-2098 obtains its
strict-band cover only up to a null exceptional set `E`. In that setting,
replace the pointwise hypothesis (3) by its almost-everywhere version and
discard

```text
union_(j in F_13) (E-jv/13).                           (17a)
```

This is still null. Every remaining base point has all thirteen translates
covered, so (14)--(15) and the same dichotomy hold for almost every base
point satisfying (7). Accordingly, one isolated two-singleton base point
refutes the pointwise hypothesis (3), whereas a positive-measure set of such
base points refutes the almost-everywhere cover inherited from THM-2098.

## 4. Frontier effect and Tournament Analysis

THM-2114 says a rank-eight through rank-eleven cover must first pay a
13-content blocker.  This theorem consumes only the rank-eight branch in which
that payment is a unique terminal blocker independent of the guard and the
other seven columns satisfy (2).  A guard-content blocker is separate; ranks
nine through eleven leave eight through ten residual bands and different
incidence ledgers. Under the stated rank-eight assumptions, the arithmetic
payment becomes one of only two exact patterns on thirteen points, pointwise
under (3) and almost everywhere under THM-2098's weaker cover. The next
lawful attacks are:

```text
force two singleton windows by (17),
exclude the six-toothpick-plus-singleton partition using step labels r_i,
or exclude the one-double seven-toothpick pattern by higher intersections. (18)
```

The challenged assumption is that all size-two residual dangers are edges of
one fixed cycle.  They are adjacent only in their own terminal coordinate;
the common orbit sees seven possibly different steps `r_i^(-1)`.  Forgetting
those colors creates a false parity argument.

Candidate tournament vertices are the thirteen orbit points, the seven
colored toothpicks, singleton phases, collision points, and cover obligations.
Orienting a toothpick by the numerical order of its endpoints gives a tie
Hamiltonian path but has no invariant meaning under translation or step
change.  Score histograms, directed cycles, SCCs, and edge flips forget the
singleton-versus-double switch in (13).  The faithful carrier is the colored
chord hypergraph

```text
(F_13; D_1,...,D_7; r_1,...,r_7)
```

together with the continuous base phase `X_0`. This theorem proves the
reduction and the two incidence ledgers. THM-2120 assembles the full base-
phase torus and uses its finite character kernels to eliminate both. QED.
