---
id: THM-3879
title: "Rational torus sextic C3 packet and one-place tradeoff"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE AWAITING INDEPENDENT HOSTILE AUDIT.  An
  explicit rational sextic, obtained as the dual of a trinodal quartic, has
  exactly six A2 cusps and four A1 nodes and admits a torus equation
  4Q2^3-27Q3^2=0.  Its irreducible Cardano cubic supplies a connected
  codimension-one-unramified C3 layer over the quadratic resolvent, exactly
  the global packet absent from THM-3874.  No projective line leaves one
  normalization place; the best line leaves two.  This is a sharp
  counterexample laboratory, not a planar Jacobian counterexample.
source: root / post-THM-3874 global-cusp-gluing reframe, 2026-08-23
audit: >
  PROVISIONAL EXACT PROOF CANDIDATE.  The assertion-free companion
  reconstructs the trinodal quartic and dual map, torus pullbacks, bidual
  inverse, complete singular Groebner packet, inner transversality, outer
  Hessian gates, Cardano discriminant and irreducibility gate, Kummer norm,
  and the saturated sixth-power line systems.  Normal and optimized runs
  must byte-match the frozen transcript.  Independent audit must recheck
  geometric irreducibility, the complete projective singular locus, the
  height-one Kummer divisor, and the exact projective-line scope.
depends_on:
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
related:
  - THM-3851-tricuspidal-quartic-rank-two-two-place-tradeoff
  - THM-3874-three-cusp-quadratic-k3-affine-class-group
script: 04-computation/jc2_rational_torus_sextic_c3_one_place_tradeoff_thm3879.py
output: 05-knowledge/results/jc2_rational_torus_sextic_c3_one_place_tradeoff_thm3879.out
script_sha256: e3ee7cd5cba7e9653602463ee4220071b5e5e93fbe28404468ccf1eaea4994bc
output_sha256: d868045804a4d9b2ed08a49f499729660dd64503abb7945c774c82f9cd1048f8
semantic_sha256: 892b5f28dbe930f764abf25c72b663b36bd23b84519728380ec5894cbe3adb01
hash_basis: raw LF bytes
---

# THM-3879 -- the missing C3 packet exists, but it costs the one-place chart

**PROVED + VERIFIED-EXACT CANDIDATE AWAITING INDEPENDENT HOSTILE AUDIT.**
Work over an algebraically closed field `k` of characteristic zero.  Put

```text
X=ST(S^2-2T^2),
Y=ST(S^2-T^2),
Z=(S^2-T^2)(S^2-2T^2).                                    (1)
```

The degree-four map `[S:T] |-> [X:Y:Z]` is a rational trinodal quartic.
Its dual curve has normalization

```text
A=-(S^2-T^2)^2(S^2+2T^2),
B= (S^2-2T^2)^2(S^2+T^2),
C=2S^3T^3.                                                  (2)
```

Define

```text
Q_2=2A^2+3AB+B^2+11C^2,
Q_3=C(4A^2+5AB+2B^2+14C^2),
F=4Q_2^3-27Q_3^2.                                          (3)
```

Then `F` is the irreducible equation of the image of `(2)`.  The sextic
`V(F)` is rational and has exactly

```text
6 A2 cusps + 4 A1 nodes.                                   (4)
```

The six cusps are the transverse intersections `Q_2=Q_3=0`; in particular
they lie on the contact conic `Q_2=0`.  The cubic

```text
u^3-Q_2u-Q_3                                                (5)
```

is irreducible over `k(P2)` and has discriminant `F`.  Its quadratic
resolvent therefore carries a connected cyclic cubic layer which is
unramified at every height-one valuation.  This is the precise global
three-divisible packet that the affine class-group calculation of THM-3874
proved absent for the three-cusp quintic.

The gain has an exact boundary.  For every nonzero line `ell` in the plane
of `(2)`, the pullback `nu^*ell` has at least two support points on `P1`.
The line `C=0` attains two:

```text
nu^*C=2S^3T^3.                                               (6)
```

Thus no projective affine chart of this embedded sextic has affine-line
normalization.  Its best chart has normalization `G_m`.  The construction
supplies the desired `C3` gluing and rationality simultaneously, but not the
one-place infinity forced on a planar Keller branch.

## 1. The trinodal quartic and its dual

Let

```text
q_x=S^2-T^2,             q_y=S^2-2T^2,             q_z=ST.   (7)
```

The three quadratics have pairwise disjoint reduced zero sets, and `(1)` is

```text
[X:Y:Z]=[q_yq_z:q_xq_z:q_xq_y].                             (8)
```

The two roots of `q_x` map to `[1:0:0]`, the roots of `q_y` to `[0:1:0]`,
and the roots of `q_z` to `[0:0:1]`.  Each pair has distinct tangent
directions, giving three ordinary nodes.  The map is birational: on its
generic image,

```text
T/S=XY/((2Y-X)Z).                                           (9)
```

Since a generic line pulls back with degree four, the image is a quartic.
Its arithmetic genus is three, so the three nodes exhaust its singularities.

Taking the cross product of the two homogeneous derivative rows gives

```text
(partial_S(X,Y,Z)) cross (partial_T(X,Y,Z))=4(A,B,C),       (10)
```

with `(A,B,C)` exactly as in `(2)`.  There is no common zero in `(2)`, so it
defines a degree-six morphism to the dual plane.

## 2. Torus equation and birationality

Put

```text
H=ST(S^4+S^2T^2+2T^4).                                    (11)
```

Direct substitution gives the stronger pair of identities

```text
Q_2(A(S,T),B(S,T),C(S,T))=3H^2,
Q_3(A(S,T),B(S,T),C(S,T))=2H^3.                            (12)
```

Hence `F` vanishes on the image.  There is no hidden covering degree.  Along
the parametrization, the gradient of `F` is a common nonzero rational factor
times

```text
[X:Y:Z].                                                     (13)
```

Equation `(9)` therefore recovers `[S:T]` rationally from a generic point of
the dual image.  The map `(2)` is birational.  A generic target line pulls
back to a binary sextic, so the image has degree six.  Its irreducible
equation divides the degree-six polynomial `F`; consequently it equals `F`
up to a scalar.  This proves geometric irreducibility as well as the claimed
normalization.

## 3. The complete singular packet

In the chart `C=1`, write `x=A,y=B`.  Exact elimination of

```text
F=F_x=F_y=0                                                  (14)
```

has one row linear in `x` and eliminant

```text
(y^4-13y^2+128)^2 (y^4+272y^2+64).                         (15)
```

The two factors are squarefree and coprime.  The first gives the four finite
inner points

```text
y^4-13y^2+128=0,            x=(y^3-13y)/16.                 (16)
```

They are exactly `Q_2=Q_3=0`.  The Jacobian of `(Q_2,Q_3)` is a unit modulo
their reduced coordinate ring, so `Q_2,Q_3` are analytic parameters.  In
those parameters the local equation is

```text
4Q_2^3-27Q_3^2=0,                                          (17)
```

an `A2` cusp.

The second factor gives four outer singularities

```text
y^4+272y^2+64=0,            x=y(y^2+136)/192.               (18)
```

The determinant of the affine Hessian of `F` is a unit at every point in
`(18)`, so all four are ordinary `A1` nodes.

On the line `C=0`, in the chart `B=1`,

```text
F=4(A+1)^3(2A+1)^3.                                        (19)
```

The only projective singularities there are `A=-1,-1/2`.  At both, the
Jacobian of `(Q_2,Q_3)` is nonzero, so they are the remaining two `A2`
cusps.  Equations `(15)-(19)` prove that no projective singular point has
been omitted.  Finally

```text
p_a(6)=10=6 delta(A2)+4 delta(A1),                           (20)
```

consistent with the explicit rational normalization.

## 4. The Cardano and Kummer packets

The discriminant of `(5)` is exactly `F`.  It is also irreducible.  Indeed a
rational root of the monic cubic would, by Gauss' lemma and the homogeneous
degree comparison, homogenize to a linear form

```text
L=aA+bB+cC.                                                  (21)
```

At `C=0`, the root identity reads

```text
L(L^2-(A+B)(2A+B))=0.                                      (22)
```

The displayed product of two distinct lines is not a square, so `(22)`
forces `a=b=0`.  For `L=cC`, comparison of the `A^2C` and `ABC` coefficients
would require simultaneously

```text
2c+4=0,                         3c+5=0,                     (23)
```

which is impossible.  Thus `(5)` is irreducible.  Since its irreducible
discriminant is nonsquare, its Galois group is `S3`.

Let `K=k(P2)`, choose `rho in k` with `rho^2=27`, and put

```text
E=K(w),                 w^2=27Q_3^2-4Q_2^3=-F,
gamma=rho Q_3+w.                                              (24)
```

Then

```text
Norm_(E/K)(gamma)=4Q_2^3.                                  (25)
```

Over the generic point of `Q_2=0`, the quadratic cover splits into the two
divisors `w=+/-rho Q_3`.  On the divisor where `gamma=0`, the conjugate
factor is a unit, so `(25)` gives

```text
div(gamma)=3D                                                (26)
```

at every affine height-one valuation.  Cardano identifies `E(cuberoot
gamma)/E` with the `C3` subextension of the `S3` normal closure of `(5)`.
Irreducibility makes it connected.  Thus `(26)` is a genuine nontrivial
codimension-one-unramified cyclic layer, not merely a square/norm identity.

The cubic `(5)` is globally monogenic and already normal.  In the chart
`C=1`, its singular ideal has the exact reduced support

```text
48u+y^2+88=0,
192x-y^3-136y=0,
y^4+272y^2+64=0,                                           (27)
```

so its singular locus consists of four closed points.  A hypersurface is
`S2`, and regularity in codimension one therefore proves normality.
THM-3801 prevents this particular completion from containing a constant-unit
degree-three etale plane atlas.  The point of the example is different:
unlike THM-3874, its
quadratic resolvent really does possess the required cyclic layer, so future
searches must now solve the atlas and infinity problems rather than recreate
local cusp torsion.

## 5. Why every line leaves at least two places

The three coordinate forms in `(2)` expand as

```text
A=-S^6+3S^2T^4-2T^6,
B= S^6-3S^4T^2+4T^6,
C=2S^3T^3.                                                  (28)
```

Suppose a nonzero line pulled back to a divisor supported at one point.
Its degree-six binary form would be

```text
aA+bB+cC=kappa(alpha S+beta T)^6,             kappa!=0.     (29)
```

The left side has zero `S^5T` and `ST^5` coefficients.  In characteristic
zero, `(29)` forces `alpha beta=0`.  If the right side is a multiple of
`S^6`, the `S^4T^2,S^3T^3,S^2T^4` rows force respectively
`b=c=a=0`; the `S^6` row then contradicts `kappa!=0`.  The `T^6` case is
identical.  Hence the one-place line system is empty.

By `(6)`, the line `C=0` has support exactly `{S=0,T=0}`.  The minimum number
of normalization places at infinity is therefore two.  This proves the
claimed sharp projective-line tradeoff.

## 6. Scope and reproduction

This theorem constructs no Keller atlas and no Jacobian counterexample.  It
does not exclude nonlinear target-coordinate changes, different rational
torus sextics, higher-degree branch curves, or nonmonogenic cubic orders.
It identifies a precise next design target: preserve a nontrivial splitting-
conic/C3 packet while collapsing the normalization boundary from two places
to one without making the finite completion monogenic.

Run

```bash
python3 04-computation/jc2_rational_torus_sextic_c3_one_place_tradeoff_thm3879.py
python3 -O 04-computation/jc2_rational_torus_sextic_c3_one_place_tradeoff_thm3879.py
```

and compare both streams byte-for-byte with
`05-knowledge/results/jc2_rational_torus_sextic_c3_one_place_tradeoff_thm3879.out`.
The companion uses exact rational polynomial arithmetic and has no inactive
`assert` gates.
