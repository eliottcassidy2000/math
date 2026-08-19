---
id: THM-3560
title: "Jelonek Euler gate and monomial target-shear no-go"
status: >
  PROVED + VERIFIED-EXACT.  For every target coordinate surface T in the
  fixed THM-1300 map, constructible Euler integration through the exact
  3/1/0 fibre stratification gives chi(F^-1(T))=3-2chi(T intersect L)-
  chi(T intersect E).  Hence a coordinate pullback must satisfy
  2chi(T intersect L)+chi(T intersect E)=2.  For every m>=1 and nonzero
  lambda, the target shear c+lambda*b^m has Jelonek section normalized by
  s^2=4+3lambda*b^(m+1), with one finite point removed.  Its section Euler
  characteristic is -m, its omitted-curve intersection has m+1 points, and
  its full pullback has Euler characteristic m+2, never that of A2.  The
  specialized core cubic has no rational root by the a-infinity valuation,
  so the pullback is irreducible.  Thus no monomial b-shear pullback has any
  source-coordinate factor, in any degree.
source: kps-s188
depends_on:
  - THM-1300-jacobian-counterexample-dixmier-A3-explicit
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
related:
  - THM-1335-trisection-modulus-master-identities-trace-polynomiality
  - THM-3546-invariant-graph-keller-descent-criterion
  - THM-3559-affine-target-coordinate-pullback-no-go
companion: 04-computation/jacobian_jelonek_euler_monomial_shear_kps_s188.py
output: 05-knowledge/results/jacobian_jelonek_euler_monomial_shear_kps_s188.out
script_sha256: a1f757e1bb83db9475d8e6d771cb8c7b88eefd4a55c766bafee5db6b246b3c60
output_sha256: 8aec39457572d62841f9faa25ed836e0ec9b8d758c7712adee804ff8e4f92675
hash_basis: LF-normalized bytes
---

# THM-3560 -- Jelonek Euler gate and monomial target-shear no-go

**PROVED + VERIFIED-EXACT.**  The fixed three-dimensional counterexample has
enough global fibre information to test nonlinear target hypersurfaces before
attempting coordinate recognition.  The resulting Euler gate closes an
all-degree target-shear family and unexpectedly passes through the
equianharmonic and lemniscatic curves at its first two nonlinear levels.

All varieties are over `C`.  Write `chi` for compactly supported topological
Euler characteristic.

## 1. The constructible Euler gate

Use target coordinates `(a,b,c)`.  THM-2473 proves that the fixed Keller map
`F:A^3->A^3` has the exact fibre-count stratification

```text
3 points  on A^3 minus V(L),
1 point   on V(L) minus E,
0 points  on E,                                             (1)
```

where

```text
L=27a^2c^2-18abc+b^3c+16a-b^2,                         (2)

E={ (4/(27t^2),4/(3t),t) : t in C* }.                  (3)
```

Let `T` be any closed target subvariety.  Put

```text
D=(T intersect V(L))_red,          e=(T intersect E)_red.  (4)
```

Euler integration of the constructible fibre-count function `(1)` gives

```text
chi(F^(-1)(T))
 =3 chi(T minus D)+chi(D minus e)
 =3chi(T)-2chi(D)-chi(e).                               (5)
```

No local triviality over the lower strata is being assumed: `(5)` is the
standard additivity formula for a quasi-finite constructible map, whose
fibre Euler characteristic is its finite cardinality.

If `T` is a target coordinate hypersurface, then `T~=A^2` and `chi(T)=1`.
If its complete pullback is also a source coordinate hypersurface, then
`F^(-1)(T)~=A^2`.  Equation `(5)` therefore imposes the necessary condition

```text
2chi(D)+chi(e)=2.                                       (6)
```

In particular,

```text
chi(e) is even.                                         (7)
```

When `T` meets `E` in finitely many reduced points, `(7)` is a literal even
intersection-support rule.  Multiplicities do not enter: Euler integration
sees the reduced constructible strata.

## 2. Monomial target shears

For `m>=1` and `lambda!=0`, consider the triangular target coordinate

```text
g_m(a,b,c)=c+lambda b^m                                (8)
```

and its coordinate surface

```text
T_m=V(g_m)~=A^2.                                       (9)
```

Every `T_m` contains the collision value `(-1/4,0,0)`.  Substituting
`c=-lambda b^m` into `(2)` gives the affine Jelonek section

```text
D_m: Q_m(a,b)=0,

Q_m=27lambda^2 a^2 b^(2m)
    +(16+18lambda b^(m+1))a
    -b^2(1+lambda b^(m+1)).                            (10)
```

As a quadratic in `a`, its discriminant is the exact cube

```text
disc_a(Q_m)=4(4+3lambda b^(m+1))^3.                    (11)
```

This identity is the curve-level shadow of THM-1335's global
cube-plus-square trisection identity.

## 3. Normalization and Euler characteristic

Put

```text
n=m+1,                P(b)=4+3lambda b^n.              (12)
```

The polynomial `P` is squarefree.  The normalization of `D_m` is

```text
C_m: s^2=P(b),                                         (13)
```

with the single affine point `(b,s)=(0,-2)` removed.  One exact normalization
coordinate is

```text
s=(54lambda^2 b^(2m)a+16+18lambda b^n)/(2P).           (14)
```

Indeed, on `Q_m=0`, direct expansion gives

```text
(54lambda^2 b^(2m)a+16+18lambda b^n)^2=4P^3.           (15)
```

Conversely, the quadratic formula recovers `a` from `(b,s)`.  At `b=0`, the
branch `s=2` extends to the smooth point `(a,b)=(0,0)`; rationalizing the
quadratic formula makes this extension explicit.  The branch `s=-2` has a
pole and is the one removed point.  At a root of `P`, both quadratic branches
join in one unibranch cusp; normalization is bijective there, so no Euler
correction is introduced.

The smooth projective model of `(13)` has

```text
genus floor((n-1)/2),
one point at infinity if n is odd,
two points at infinity if n is even.                   (16)
```

After deleting those infinity points and `(0,-2)`, `(16)` gives in both
parities

```text
chi(D_m)=1-n=-m.                                        (17)
```

This also identifies two familiar boundary cases:

- `m=2`: `s^2=4+3lambda b^3`, an equianharmonic elliptic curve (`j=0`);
- `m=3`: `s^2=4+3lambda b^4`, a lemniscatic elliptic curve (`j=1728`).

Thus the Bernoulli-lemniscate and CM-signature objects appear here not as
decorative analogies but as the normalizations of the first nonlinear
Jelonek sections.

## 4. The omitted curve and the all-degree obstruction

On the omitted curve `(3)`, equation `(8)` becomes

```text
t+lambda(4/(3t))^m=0,

t^(m+1)=-lambda(4/3)^m.                                (18)
```

The right side is nonzero, so `(18)` has exactly `m+1` distinct points:

```text
chi(e_m)=m+1.                                           (19)
```

For even `m`, this is already odd and violates the parity gate `(7)`.  For
all `m`, substitute `(17)` and `(19)` into `(5)`:

```text
chi(F^(-1)(T_m))
 =3-2(-m)-(m+1)
 =m+2.                                                  (20)
```

Since `m+2!=1`, the pullback hypersurface

```text
V(F3+lambda F2^m)                                      (21)
```

is never isomorphic to `A^2`.  In fact it is irreducible.  THM-2473's core
fibre cubic over the function field `C(a,b)` of `T_m` specializes to

```text
L_m x^3+(4+3lambda b^(m+1))x+2lambda b^m.             (22)
```

At the discrete valuation `a=infinity`, `v(a)=-1`, its three terms at an
alleged rational root of valuation `k in Z` have valuations

```text
-2+3k,                         k,                      0.  (23)
```

For `k<=0`, the first is the unique minimum.  For `k>=1`, the constant term
is the unique minimum.  A vanishing valued sum must attain its minimum at
least twice, so `(22)` has no root in `C(a,b)` and is irreducible.  Since `F`
is quasi-finite, every irreducible component of the hypersurface `(21)`
would dominate `T_m`; `(22)` therefore proves that `(21)` has only one
component.

Thus `F3+lambda F2^m` is irreducible and not a coordinate.  It has no
source-coordinate factor.  This closes the complete monomial `b`-shear
family in every degree.

The theorem does not rule out a mixed polynomial shear `c+phi(a,b)` or a
nonlinear target coordinate not triangular in `c`.  Those are the surviving
descent cells.

## 5. Exact verification

Reproduce with

```bash
python3 04-computation/jacobian_jelonek_euler_monomial_shear_kps_s188.py
python3 -O 04-computation/jacobian_jelonek_euler_monomial_shear_kps_s188.py
```

The ordinary and optimized transcripts agree.  The companion verifies the
universal discriminant identity by replacing `b^(m+1)` with one symbol, and
then independently checks `(11)`, `(15)`, `(17)--(20)` for `1<=m<=12`,
covering both parity classes and the two CM elliptic cases.  It also checks
the valuation minimum across `-24<=k<=24`; the two-case proof after `(23)`
covers every integer.  The all-degree claims are displayed algebra, not an
extrapolation from finite controls.

**QED.**
