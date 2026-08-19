---
id: THM-3554
title: "Punctured Kummer normal form of the fixed Keller collision surface"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED WITH
  TYPING/TOOLING REPAIRS.  The curved collision component C=0 of the fixed
  THM-1300 three-dimensional Keller map is G_m x A^1.  Its restriction to
  the target plane F3=0 is, after explicit Laurent source and polynomial
  target automorphisms, exactly (s,b)->(b,4s^2).  Thus it is a degree-two
  finite etale Kummer cover of the complement of a parabola and carries the
  known orbit collision.  Its natural affine-plane completion ramifies, and
  no everywhere-etale A^2 filling can preserve this quadratic function-field
  extension.  This is a punctured planar near-counterexample, not JC(2).
source: kps-s188
depends_on:
  - THM-1300-jacobian-counterexample-dixmier-A3-explicit
related:
  - THM-3543-torus-quotient-ramification-square-no-go
  - THM-3545-catalan-self-intersection-keller-thickening-boundary
  - THM-3546-invariant-graph-keller-descent-criterion
companion: 04-computation/jacobian_punctured_kummer_collision_surface_kps_s188.py
output: 05-knowledge/results/jacobian_punctured_kummer_collision_surface_kps_s188.out
script_sha256: a4301507f47fd18d45097d8b8bf079652de5640b1c36be439fd190ea8d1683df
output_sha256: 789ed0276d8e2ad5e1e8d358ef56e2a6e7815bae52c98cc0045664150a81b5e1
hash_basis: LF-normalized bytes
---

# THM-3554 -- punctured Kummer normal form of the fixed Keller collision surface

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED WITH
TYPING/TOOLING REPAIRS.**  The two nonfixed points in THM-1300's triple
collision already lie on a two-dimensional etale noninjective map.  The
missing ingredient for a planar counterexample is completely visible: its
source is a punctured plane, and filling that puncture restores the quadratic
branch divisor.

The field has characteristic different from two.

## 1. The collision surface is a Laurent plane

For THM-1300's map `F=(F_1,F_2,F_3)`, put

```text
C=2-3xy-x^2z.                                           (1)
```

Then `F_3=xC`.  On the curved component `S=V(C)`,

```text
x(3y+xz)=2.                                             (2)
```

Hence `x` is a unit in `k[S]`.  Set

```text
s=x^(-1),                 v=xy.                         (3)
```

Equation `(1)` and its inverse substitutions give

```text
x=s^(-1),        y=vs,        z=s^2(2-3v),              (4)
```

and therefore an explicit ring isomorphism

```text
k[S] ~= k[s,s^(-1),v].                                  (5)
```

In particular, `S` is smooth and isomorphic to `G_m x A^1`, not to `A^2`.
The nonconstant unit `x` is already an algebraic obstruction to `S` being a
coordinate hypersurface: a polynomial plane has only constant units.  Thus
THM-3546's coordinate-hypersurface descent theorem does not silently turn
this surface into a planar Keller counterexample.

## 2. Exact restricted map

Substituting `(4)` into the first two components of `F` gives

```text
alpha=F_1|S=s^2(v+1)(v+2),
beta =F_2|S=2s(2v+3),
F_3|S=0.                                                (6)
```

The tangential Jacobian in these Laurent coordinates is

```text
det d(alpha,beta)/d(s,v)=-2s^2,                         (7)
```

which is a unit of `k[s,s^(-1),v]`.  Equivalently, in the chart `(x,y)` on
`S` it is `2/x^3`, again a unit on `S`.  Thus the restriction is etale even
though its ordinary determinant is not a constant when written in a
non-volume-preserving Laurent chart.

Now make the source Laurent automorphism

```text
(s,v) -> (s,b),             b=beta=4sv+6s,              (8)
```

whose Jacobian is `4s` and whose inverse is

```text
v=(b-6s)/(4s).                                          (9)
```

On the target plane make the polynomial automorphism

```text
(alpha,beta) -> (b,delta),
b=beta,                    delta=beta^2-16alpha.        (10)
```

Its Jacobian is `16`, and its inverse is
`alpha=(b^2-delta)/16`.  From `(6)`,

```text
delta=4s^2.                                             (11)
```

Consequently the restricted map has the exact normal form

```text
G_m x A^1 -> A^1 x G_m,
(s,b)       -> (b,4s^2).                                (12)
```

The companion verifies `(1)`--`(12)`, including the ambient determinant
`det JF=-2`, by exact symbolic arithmetic under normal and optimized Python
execution.

## 3. Kummer cover and the collision

On coordinate rings, `(12)` is

```text
k[b,delta,delta^(-1)] -> k[b,s,s^(-1)],
delta                  -> 4s^2.                        (13)
```

The target ring extension is free of rank two with basis `{1,s}`; the
derivative `8s` is a unit.  Hence `(12)` is a finite etale cover of degree
two.  Over an algebraic closure its deck involution is `s->-s`.  In the
coordinates `(s,v)` from `(3)`, keeping `b` fixed gives

```text
iota(s,v)=(-s,-v-3).                                   (14)
```

The two THM-1300 orbit points are

```text
(s,v)=(1,-3/2),          (-1,-3/2),                    (15)
```

and both map to

```text
(alpha,beta)=(-1/4,0).                                  (16)
```

Thus the collision is exactly the nontrivial Kummer deck orbit.  Put

```text
D: beta^2-16alpha=0,             U=A^2 minus D.         (17)
```

The morphism `S->U` is finite etale and surjective.  Accordingly, its
topological image and its image on algebraically closed points are exactly
`U`.  Two qualifications are load-bearing:

- as a morphism to the ambient target `A^2`, its scheme-theoretic image is
  the closure of `U`, hence all of `A^2`; and
- its image on `k`-rational points can be smaller than `U(k)`.  For example,
  over `Q` the point `(alpha,beta)=(-1/8,0)` has `delta=2`, but a preimage
  would require `s^2=1/2`.

Over `D` the two geometric square-root sheets would have to meet.

## 4. The puncture is load-bearing

The same formulas extend to the natural affine closure

```text
A^2_(s,b) -> A^2_(b,delta),       (s,b)->(b,4s^2),      (18)
```

but the Jacobian is `-8s`; the extension ramifies on `s=0`.  This is not an
artifact of the chosen closure.  Over any field of characteristic different
from two, there is no dominant etale morphism

```text
H:A^2 -> A^2_(b,delta)                                  (19)
```

whose induced extension of function fields is, over `k(b,delta)`, the
quadratic extension

```text
k(b,delta)(s),                 s^2=delta/4.             (20)
```

Indeed, dominance makes `H^*delta` a nonconstant polynomial.  It is not a
unit, so its zero divisor is nonempty.  Because `H` is etale, the pullback of
the smooth reduced divisor `delta=0` is reduced; along each of its irreducible
components `P`,

```text
ord_P(H^*delta)=1.                                      (21)
```

But `(20)` holds in `k(A^2)`, so

```text
ord_P(H^*delta)=2 ord_P(s),                             (22)
```

which is even.  Equations `(21)` and `(22)` contradict each other.  Therefore
the exact Kummer field extension cannot be filled by an everywhere-etale
polynomial plane map.

The same argument applies verbatim to a smooth integral affine source `X`
with only constant global units.  If an etale dominant `X->A^2` induced
`(20)`, then `delta` could not pull back to a unit, so it would have a reduced
height-one zero divisor; `(22)` would again make every multiplicity even.
The source `G_m x A^1` escapes this obstruction for exactly the advertised
reason: `s` is a nonconstant unit and the divisor `s=0` has been removed.

## 5. What this changes in the planar search

THM-3543 found the opposite two-dimensional shadow of the same geometry.  The
categorical torus quotient forgets `s`, contracts the whole line `C=0`, and
pays the Jacobian square `C^2`.  The present restriction keeps `s`; it keeps
the map etale and the collision finite, but pays with the nonconstant unit
`s` and a missing divisor.  The two failures are dual:

```text
forget the character s  -> polynomial plane, but ramified contraction;
retain the character s  -> etale collision, but punctured source/target. (23)
```

Together with THM-3545's nonterminating Catalan thickening, this isolates a
three-way counterexample pincer.  A genuine planar Keller counterexample must
simultaneously retain a finite collision, avoid divisorial ramification, and
replace the punctured Laurent source by `A^2` without returning to the exact
quadratic extension `(20)`.

Concrete live ways to escape `(22)` are therefore sharply limited:

1. introduce mixed transverse corrections that change the Kummer function
   field before the boundary is filled;
2. use a nonproper affine modification whose missing divisor maps to infinity
   rather than to the finite branch parabola; or
3. pass to higher sheet number with no finite branch divisor, so the
   noninjectivity is carried entirely by the asymptotic set.

These are search architectures, not constructions.  The theorem excludes
only the exact quadratic filling `(20)`.  It does not exclude higher-degree
mixed thickenings, other surfaces in the three-dimensional map, rational
modifications with a different function field, or planar Keller
counterexamples in general.
