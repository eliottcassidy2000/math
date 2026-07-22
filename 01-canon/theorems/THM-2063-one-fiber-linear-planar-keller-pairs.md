---
id: THM-2063
title: "One-fiber-linear planar Keller pairs are tame"
status: >
  PROVED. If one nonzero member of the output pencil of a complex planar
  Keller map is affine on the fibers of a linear source projection,
  coefficient descent removes the other component by a polynomial in that
  member and leaves one of two explicit triangular normal forms. The inverse
  is displayed. Hence every hypothetical JC(2) counterexample has fiber
  degree at least two in every source direction for every nonzero output-pencil
  member. This completely classifies the stated class; it does not prove or
  disprove JC(2).
source: codex-2026-07-21-JC2-one-fiber-linear
depends_on: []
related:
  - THM-1330
  - THM-1345
  - THM-1340
  - THM-2045
  - THM-2046
  - MISTAKE-205
script: 04-computation/jc2_one_fiber_linear_referee_codex_20260721.py
output: 05-knowledge/results/jc2_one_fiber_linear_referee_codex_20260721.out
script_sha256: d26bc835fe5effd65ce327f188634610445b67b6cc5dfad3f2be4cfc40884943
output_sha256: b164d0266bf3ec84fa897c08947a0c3463b20e7a3f2006531a199ac7ba5ac46b
hash_basis: normalized repository blobs (LF)
---

# THM-2063 -- one-fiber-linear planar Keller pairs are tame

Let

```text
F=(P,Q): C^2 -> C^2,       J(P,Q)=P_x Q_y-P_y Q_x=kappa in C*.
```

Suppose that, after invertible affine changes in the source and target, the
first component is affine in one source variable:

```text
P=A(x)y+B(x).                                             (1)
```

Then `F` is a tame polynomial automorphism. More precisely, exactly one of the
following two normal forms applies.

```text
A=a in C*:  P=ay+B(x),
             Q=H(P)-(kappa/a)x+b;                        (2)

A=0:         P=alpha x+beta,
             Q=(kappa/alpha)y+H(P),                      (3)
```

where `H` and `B` are arbitrary one-variable polynomials and the displayed
constants required to be nonzero are nonzero. Conversely every map in (2) or
(3) has Jacobian `kappa`.

## 1. Coefficient descent

First assume `A!=0` and expand

```text
Q=sum_(j=0)^n q_j(x)y^j,       q_n!=0.                  (4)
```

If `n>=1`, the coefficient of `y^n` in the Jacobian is

```text
n A'(x)q_n(x)-A(x)q_n'(x).                              (5)
```

Since the Jacobian is constant, (5) vanishes. In `C(x)`,

```text
(q_n/A^n)'=0,
```

so characteristic zero gives `q_n=c_n A^n` for a constant `c_n`. The
polynomial `c_n P^n` has the same leading `y`-coefficient, while
`J(P,c_nP^n)=0`. Replace `Q` by `Q-c_nP^n`. Its `y`-degree falls without
changing the Jacobian. Finite descent gives

```text
Q=H(P)+R(x).                                             (6)
```

Now

```text
kappa=J(P,R)=-A(x)R'(x).                                (7)
```

Both factors in (7) are polynomials and their product is a unit. Therefore
`A=a in C*` and `R=-(kappa/a)x+b`, proving (2).

If `A=0`, then `P=B(x)` and

```text
kappa=B'(x)Q_y(x,y).                                    (8)
```

The factors in (8) are units in `C[x,y]`. Hence `B'=alpha in C*` and
`Q_y=kappa/alpha`; its remaining `x`-polynomial can be written as `H(P)`
because `P=alpha x+beta`. This is (3). The impossible case `A=B'=0` would
give zero Jacobian.

## 2. Explicit inverses and tameness

For (2), with target coordinates `(u,v)=(P,Q)`,

```text
x=(a/kappa)(H(u)+b-v),
y=(u-B(x))/a.                                           (9)
```

For (3),

```text
x=(u-beta)/alpha,
y=(alpha/kappa)(v-H(u)).                               (10)
```

These are polynomial inverses. The proof itself factors the map into affine
maps and elementary shears: repeatedly subtract a polynomial in `P` from the
second target coordinate, then use (7) or (8). Thus `F` is tame, not merely
invertible. Direct substitution also proves the converse statements.

## 3. The all-directions curvature gate

Return to an arbitrary planar Keller pair `(P,Q)`. Let `R=sP+tQ` be any
nonzero member of its output pencil. Complete `(s,t)` to a matrix in
`GL_2(C)`, making `R` the first target coordinate. Let `(ell,m)` be any affine
linear coordinate system on the source. If

```text
deg_m R<=1,                                             (11)
```

then the transformed pair satisfies (1), so the original pair is an
automorphism. Consequently every hypothetical counterexample to `JC(2)` must
satisfy the simultaneous gate

```text
deg_m(sP+tQ)>=2                                        (12)
```

for every nonzero output direction `(s,t)` and every linear source fiber
direction `m` completed to affine coordinates.

This is stronger than saying that both displayed components have degree at
least two in the current `y`: target cancellation and every linear source
foliation are included. It is also affine-invariant. It does not cover a
general Laurent-monomial direction, because a `GL_2(Z)` change with negative
entries is not a polynomial source automorphism.

## 4. What this does and does not classify

The theorem completely classifies the one-fiber-linear stratum, including all
degrees of the unrestricted shear polynomials `B,H`. In the current planar
counterexample atlas it supplies an empty class:

```text
hypothetical counterexamples
  subset {Keller pairs satisfying (12)}.
```

Together with THM-1345, any survivor must also lose every nontrivial
`C*`-equivariance after affine conjugacy. Together with the valid universal
parts of THM-1330, it must be nonproper, have function-field degree at least
three, and ramify only at infinity. These are necessary gates, not an
enumeration of the surviving maps. Neither the inverse Jelonek problem nor
the planar Jacobian conjecture is closed here.

## 5. Assumption challenge and tournament analysis

The useful vertices are not just the two output polynomials. They are pairs

```text
(output-pencil direction, source-fiber direction),
```

labelled by the exact fiber degree in (11). Quotienting to total degrees
destroys the gate: a high-degree polynomial can become affine on a special
linear foliation, and two high-degree outputs can have a low-degree pencil
combination. The challenged assumption is therefore that a planar search may
fix its displayed coordinates or its leading total degree.

A tournament is not the proof carrier. The incidence relation `deg_m R<=1`
is bipartite and invariant under independent projective changes on the source
and target direction lines; forcing a direction between two incidences loses
that symmetry. One may orient candidates by smaller minimum fiber degree and
break ties by Newton support size, but this only schedules a search. The
coefficient descent (5)--(7) is the certificate.

## 6. Exact computational referee

The companion audit checks (5) symbolically for six successive `y`-degrees,
then generates exact rational instances of both normal forms, verifies their
Jacobians, and composes the displayed inverses in both directions. It also
checks that a nonconstant trial `A` can pass the top-coefficient equation but
fails the final unit equation, guarding against stopping the descent one layer
too early. The stored output ends in `PASS`. QED.
