# Positive moment tensors force coefficient-phase dispersion

**Status: PROVED GENERAL LEMMA / PROMOTED AS
[THM-3327](../01-canon/theorems/THM-3327-positive-moment-tensor-covering-arc-and-factorial-phase-dispersion.md).**  This note
extracts the phase mechanism used in THM-3310 from its cyclic-quartic chart
and strengthens it for the standard factorial monomial basis.  It proves no
factorial conjecture.

## Fixed-order positive-tensor sector lemma

Let `Lambda` be a complex-linear functional on a commutative algebra and fix
nonzero elements `u_1,...,u_r`.  For one integer `m>=2`, assume

```text
Lambda(u_1^k1 ... u_r^kr)>0
for every k_i>=0 with sum k_i=m.                            (1)
```

Put `P=sum_i c_i u_i`, omitting zero coefficients.  If

```text
Lambda(P^m)=0,                                              (2)
```

then the shortest closed covering-arc width of the coefficient phases
`{arg c_i}` is strictly greater than `pi/m`.  This is not the usual pairwise
circular distance, which is always at most `pi`.

Indeed, if all phases lie in a closed arc of width `pi/m`, rotate `P` so that
`0<=arg c_i<=pi/m`.  The multinomial expansion is

```text
Lambda(P^m)=
 sum_(k_1+...+k_r=m) multinomial(m;k) c_1^k1...c_r^kr
                         Lambda(u_1^k1...u_r^kr).           (3)
```

Every summand in `(3)` lies in the closed upper half-plane.  If their sum is
zero, every summand lies on its boundary.  The pure terms force every occupied
coefficient phase to be `0` or `pi/m`.  If both endpoints occur, the mixed
term `c_i^(m-1)c_j` has phase strictly between `0` and `pi`; if only one
endpoint occurs, all terms lie on one nonzero real ray.  Both alternatives
contradict `(2)`.

The strict statement fails at `m=1`: two opposite coefficient rays can cancel.
The correct order-one conclusion is only that a null vector cannot lie in an
open semicircle.

The constant `pi/m` is sharp in the class of positive tensors.  With two
basis elements, choose the mixed tensor entries positive but tending to zero;
the moment polynomial tends to `1+z^m`, whose first root has argument `pi/m`.
Thus no uniform larger sector is available without quantitative tensor data.

## Factorial and simplex specialization

For the factorial functional

```text
L(x^alpha)=alpha!,                                         (4)
```

condition `(1)` holds for every finite standard-monomial support and every
`m`: a product of monomials is a monomial and its factorial is positive.
Equivalently, `(4)` is integration against the product exponential measure on
the positive orthant.  On a fixed homogeneous degree, THM-3018 gives the same
positivity through the simplex integral.

Consequently, for every `m>=2`,

```text
L(P^m)=0  =>  coefficient-phase covering width(P)>pi/m.    (5)
```

The first two factorial moments give a much stronger joint statement.  If the
standard-monomial coefficient phases of `P` lie in a closed semicircle and
`L(P)=0`, positivity of the weights `L(x^alpha)` forces every occupied phase
to one of the two endpoints.  After one global rotation, `P` therefore has
real coefficients.  But the product-exponential representation gives

```text
L(P^2)>0                                                    (6)
```

for every nonzero real polynomial `P`.  Hence

```text
L(P)=L(P^2)=0
  => coefficient-phase covering width(P)>pi.               (7)
```

The homogeneous simplex version of `(6)` is simply the strict positivity of
the integral of a nonzero real polynomial square.  Thus `(7)` applies to both
full `FC(n)` and homogeneous `HFC(n)` in their standard monomial coordinates.
Moreover `L(P)=0` writes the origin as a strict convex combination of the
occupied unit phase vectors, with weights `L(x^alpha)|c_alpha|`.  Equation
`(6)` excludes the collinear opposite-ray case.  Hence the origin lies in the
two-dimensional interior of their phase polygon: at least three noncollinear
coefficient rays are necessary.  This immediately re-excludes support one and
two, but does not exclude a three-ray polygon surrounding the origin.

## Relation to the cyclic quartic

THM-3310 works in the Fourier basis

```text
zbar, z^2, z zbar^2, z^3 zbar, zbar^4.
```

That basis is not the standard positive monomial basis, so `(7)` must not be
transported into it.  What THM-3310 verifies is exactly the fixed-order
hypothesis `(1)` at `m=3`: all `35` cubic tensor entries are positive.  The
strict `pi/3` coamoeba barrier is therefore an instance of `(1)--(3)`, while
its exact lopsided radii use additional quantitative tensor information.

## Symmetry-vacuity boundary

There is a useful stopping rule before transporting the standard-basis gate.
Suppose an order-`q` automorphism permutes a real monomial basis, preserves
`L`, and `P` is an eigenvector for a primitive character `chi`.  Then

```text
L(P^j)=chi^j L(P^j),
```

every occupied monomial orbit has full length `q`, and its coefficients carry
all `q` equally spaced phases.  Their shortest covering width is
`2*pi*(q-1)/q`.  For `q>=3`, `L(P)=L(P^2)=0` and the standard-basis `>pi` gate
are therefore automatic.  At `q=2`, total width `pi` forces all antipodal
orbit axes to coincide, making `P` global-phase real and hence
`L(P^2)!=0`; otherwise the width is already greater than `pi`.

The cyclic quartic has `q=3`.  Its standard barycentric phase gate is thus
vacuous, not a new support-five exclusion.  Passing to the Fourier basis
destroys the automatic positivity guarantee, so the order-three tensor must
be reverified there—as THM-3310 does.  The live leverage is its quantitative
Fourier tensor on support five and other unresolved cells, not a second use of
the standard-basis semicircle test.

## Typed boundary and next use

```text
source:     labelled coefficient vector and a positive order-m tensor;
map:        multinomial expansion into a half-plane;
preserved:  vanishing/nonvanishing of the selected moment;
lost:       coefficient magnitudes and later-moment compatibility;
sidecar:    quantitative tensor weights or the next moment ideal;
test:       phase chamber first, saturation/elimination second.
```

The lemma is a pre-elimination gate.  It cannot prove that every phase chamber
outside the forbidden sector contains a zero, nor can it identify Fourier
coordinates with standard monomial coordinates.  Its immediate use is to
prune unresolved projective phase charts before exact saturation and to
require every full factorial candidate to wrap strictly around the coefficient
origin—after first checking whether symmetry has made the gate automatic.
