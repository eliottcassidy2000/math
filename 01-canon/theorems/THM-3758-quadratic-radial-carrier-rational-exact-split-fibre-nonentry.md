---
id: THM-3758
title: "Quadratic radial-carrier rational-exact split-fibre nonentry"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  A four-parameter
  family Q=XA(XT)+X^2B(XT), with A linear, B'=beta A, and
  delta=2BA'-AB' a nonzero constant, consists of smooth reducible
  noncoordinates with explicit rational constant-Jacobian mates.  Every
  rational mate belongs to one primitive torsor, whose principal parts on
  the two components of Q=0 are opposite.  No target-only correction can
  regularize both, so no member has a polynomial mate.  This is a rationally
  exact near-counterexample family, not a planar Keller pair.
source: root + jc_sparse_direct_search / 2026-08-23
audit: >
  PASS.  An independent hostile audit rederived the smoothness Bezout
  identity, quadratic generic-fibre equation and discriminant, exact rational
  primitive, geometrically integral constant field, residual-divisor
  irreducibility after localization, and the opposite signed principal parts.
  Gamma is unrestricted in both theorem and companion.  Normal, optimized,
  and frozen output agree; script/output/semantic hashes and CHECKS=87 match.
depends_on: []
related:
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
  - THM-3551-one-ray-planar-jacobian-mate-no-go
  - THM-3598-danielewski-rational-exact-polar-graph-family-and-classification
  - THM-3754-affine-variable-euclidean-descent-classification
  - THM-3755-composite-monomial-generic-fibre-residue-obstruction
script: 04-computation/jc2_quadratic_radial_carrier_split_fibre_thm3758.py
output: 05-knowledge/results/jc2_quadratic_radial_carrier_split_fibre_thm3758.out
script_sha256: b9d0a9bd1de2237750ce760e3db47a4c6a03c1a4c584bcbb4a6c1ba3b6c19ec5
output_sha256: 696b216444d3b373401f31d5865bf867de161b6cc191fdbd91817f9b28ac0a01
semantic_sha256: ca31fa26678b1998450d1006a90fbcfe319096ad6b6fbc17f2006acc64fc7ecb
hash_basis: raw LF bytes
---

# THM-3758 -- rational exactness reaches a split-fibre wall

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  The residue
obstructions in THM-3551 and THM-3755 stop a rational primitive on the
generic fibre.  The present family passes that gate: its Hamiltonian time
form is rationally exact.  It still fails globally, because one primitive
has opposite poles on two vertical components and target shears cannot
distinguish them.

Let `k` be an algebraically closed field of characteristic zero.  Choose

```text
a0,a1,beta in k*,                  gamma in k,
delta=2a1 gamma-beta a0^2 !=0,                         (1)
```

and define

```text
z=XT,
A(z)=a0+a1z,
B(z)=gamma+beta a0z+(beta a1/2)z^2,
Q(X,T)=X A(z)+X^2 B(z).                                (2)
```

Then:

1. `Q` has no critical point in `k^2`;
2. `Q=X[A(z)+XB(z)]` has a reducible zero fibre and is not a coordinate;
3. for every `c in k*`, `Q` has the rational Jacobian mate

```text
Y=A(z)+2XB(z),
P0=-c/(2Q) [z+Y/(a1+2beta Q)];                         (3)
```

4. every rational mate is `P0+H(Q)`, with `H in k(Q)`;
5. no rational mate is regular on `k^2`, so no polynomial mate exists.

A normalized integral member is

```text
(a0,a1,beta,gamma)=(1,-2,-1,0),
Q=X-2X^2T-X^3T+X^4T^2.                                (4)
```

It is a smooth degree-six noncoordinate with a rational mate of Jacobian
one, but it is not a Jacobian-conjecture counterexample.

## 1. One Bezout constant proves smoothness and coprimality

The parameterization in `(2)` was chosen so that

```text
B'=beta A,
2BA'-AB'=2a1B-beta A^2=delta.                         (5)
```

On the axis `X=0`, direct differentiation gives `Q_X=a0`, so the gradient is
nonzero.  On `X!=0`, `(X,T)->(X,z)` is an etale chart because

```text
J_(X,T)(X,z)=X.                                       (6)
```

In that chart a critical point would satisfy

```text
Q_X|z=A+2XB=0,                 Q_z|X=X(A'+XB')=0.      (7)
```

But twice `B` times the second bracket minus `B'` times the first bracket is
the nonzero constant `delta` in `(5)`.  Thus `(7)` is impossible and `Q` is
smooth everywhere.

The same constant records the zero-fibre geometry:

```text
Res_z(A,B)=a1 delta/2 !=0.                             (8)
```

Hence `A,B` are coprime.  The residual factor

```text
R=A(z)+XB(z)                                           (9)
```

is primitive and linear in `X` in `k[z][X]`, so it is irreducible there.
Localizing `k[X,T]` at `X` identifies it with `k[X,X^-1,z]`; irreducibility
survives that localization, and `R` is not divisible by `X`.  Therefore

```text
Q=X R                                                  (10)
```

has exactly the two distinct reduced components needed below.  In
particular `Q` is reducible and cannot be a coordinate polynomial.

The hypotheses have sharp geometric roles.  If `a0=0`, the origin is
critical.  If `delta=0`, the two equations in `(7)` are dependent: for any
`z` with `A(z)!=0`, taking `X=-a1/[beta A(z)]` solves both.  Thus the
nonzero conditions in `(1)` are not genericity decoration.

## 2. The quadratic generic fibre has an exact time form

Put `L=Q` and retain the chart `(X,z)`.  Define

```text
Y=A+2XB.                                               (11)
```

Eliminating `X` from the quadratic equation `L=XA+X^2B` gives

```text
Y^2=A^2+4BL=:D(z,L).                                  (12)
```

The derivative constraint in `(5)` makes this double cover integrable:

```text
D_z=2AA'+4LB'=2A(a1+2beta L),
dY=(a1+2beta L) A dz/Y.                               (13)
```

Its discriminant as a quadratic in `z` is

```text
disc_z(D)=-8 delta L(a1+2beta L).                     (14)
```

It is nonzero over `k(L)`, so `D` is not a rational square, the generic
fibre is geometrically integral, and the constant field of its tangent
derivation is exactly `k(L)`.

For any rational function `P`, differentiation at fixed `L` gives

```text
J_(X,T)(P,Q)=-XY dP/dz |_L.                           (15)
```

The identities

```text
Y+A=2(A+XB)=2L/X,
1/(XY)=1/(2L)(1+A/Y)                                  (16)
```

therefore turn `J(P,Q)=c` into

```text
dP/dz |_L=-c/(2L)(1+A/Y).                             (17)
```

Equation `(13)` integrates `(17)` exactly and gives `(3)`.  Conversely two
rational solutions differ by an element of the constant field `k(L)`.  This
proves both rational existence and the complete primitive torsor
`P0+k(Q)`.

This is a materially stronger near miss than a deep finite coefficient jet:
the generic Hamiltonian response has actually been integrated.

## 3. Opposite vertical principal parts forbid a regular mate

It remains to put back the divisor information lost on the generic fibre.
Along the component

```text
D_X: X=0                                               (18)
```

one has `z=0`, `Y=a0`, and `Q` has order one because `R=a0` there.  The
quantity in brackets in `(3)` specializes at `Q=0` to

```text
z+Y/a1=a0/a1.                                         (19)
```

Along the other component

```text
D_R: R=A+XB=0,                                        (20)
```

the coordinate `X` is generically nonzero, `Y=A+2XB=-A`, and again `Q`
has order one.  Now the same bracket is

```text
z-A/a1=-a0/a1.                                        (21)
```

Consequently `P0` has simple `1/Q` principal coefficients

```text
D_X: -c a0/(2a1),                 D_R: +c a0/(2a1).    (22)
```

They are nonzero and opposite.  A correction `H(Q)` with a pole of order at
least two creates an uncancellable higher pole on both components.  A
regular `H` cancels neither.  If `H` has a simple pole `h_-1/Q`, it adds the
same coefficient `h_-1` on both components, so it cannot cancel both values
in `(22)`.  Therefore no element of the complete torsor `P0+k(Q)` is regular
on `k^2`.  This proves the polynomial no-go.

The obstruction is not a generic-fibre residue: that debt was paid in
`(13)`.  It is the sidecar lost when the two components of one target value
are collapsed to the single parameter `L`.

## 4. Exact controls and construction lesson

Reproduce with

```bash
python3 -B 04-computation/jc2_quadratic_radial_carrier_split_fibre_thm3758.py
python3 -B -O 04-computation/jc2_quadratic_radial_carrier_split_fibre_thm3758.py
```

The assertion-free companion verifies `(5)`, `(8)`, `(12)--(17)`, the
signed principal parts `(19)--(22)`, five admissible parameter controls,
the `delta=0` torus-critical and `a0=0` axis-critical boundaries, the
normalized support `(4)`, and 52 full constant-Jacobian linear systems
through total mate degree twelve.  Every bounded system has augmented-rank
gap one.  These are exact hostile controls; Sections 1--3 prove the whole
parameter family.  **QED.**

The design lesson is precise.  Moving from one composite invariant to two
carrier powers can make the time form rationally exact.  The remaining debt
is vertical: a future construction must either make the primitive's
principal parts agree on every component of a special fibre, or provide a
second response channel whose correction is not constant on that fibre.
