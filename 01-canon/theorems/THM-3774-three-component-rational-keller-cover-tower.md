---
id: THM-3774
title: "Three-component rational Keller cover tower and cuspidal degree-three near miss"
status: >
  PROVED + VERIFIED-EXACT + PENDING INDEPENDENT HOSTILE AUDIT.  For every
  m>=1 there is a smooth noncoordinate polynomial U with three pairwise
  disjoint zero-fibre components and a rational mate P satisfying J(U,P)=1.
  The generated rational target field has exact degree m+1, with minimal
  trinomial t^(m+1)-mPt+mU and an explicit irreducible binomial
  discriminant.  The log-canonical blowdown W=UP is polynomial and has
  component spectrum (1,0,0), so the complete rational-mate torsor contains
  no polynomial mate.  At m=2 the cover is an irreducible non-Galois cubic
  with S3 closure and A2 cuspidal discriminant.  This is a rational Keller
  near miss, not a polynomial Keller map or a planar Jacobian counterexample.
source: root + jc_zero_debt_lift / nonradial equalizer and cover-tower session, 2026-08-23
audit: >
  SELF-AUDITED PROOF CANDIDATE.  The all-m chart identity and off-exponent
  error are checked symbolically.  Direct source identities,
  irreducibility specializations, component spectra, minimal relations,
  discriminants, nearest hostile exponents, and degeneration controls are
  checked for m=1,...,7; gradient resultants are checked for m=1,...,5.
  Normal and optimized runs byte-match.  Independent hostile audit remains
  due, especially for the Laurent-binomial irreducibility and generic
  transposition-inertia arguments.
depends_on:
  - THM-3770-vertical-principal-part-equalizer-and-log-canonical-dressing-gate
related:
  - THM-3758-quadratic-radial-carrier-rational-exact-split-fibre-nonentry
  - THM-3772-variable-flank-three-charge-rational-exact-classification
  - THM-3773-quadratic-rational-keller-cover-degree-and-target-word-nonpolynomialization
script: 04-computation/jc2_three_component_rational_keller_cover_tower_thm3774.py
output: 05-knowledge/results/jc2_three_component_rational_keller_cover_tower_thm3774.out
script_sha256: b1a60faee125b9d3a82a220026c26623ce4a546f0de6d4feef9b075c988e96de
output_sha256: dd8ce3d9d2a5253652693d61e09350ae8c7b05a7f9741872da31b04c9fa909ff
semantic_sha256: 7afa16470c6cb9f0cdc68c83c258dae8dd6223f80b42f5f6717cef38ab489b96
hash_basis: raw LF bytes
---

# THM-3774 -- a rational Keller tower reaches the cubic cusp boundary

**PROVED + VERIFIED-EXACT + PENDING INDEPENDENT HOSTILE AUDIT.**  The
vertical-equalizer route does not stop at the quadratic near miss of
THM-3758.  It contains rational Keller seeds of every cover degree.  The
first degree not excluded by the Galois Keller theorem already has the
degree-three, non-Galois, and cuspidal-discriminant anatomy sought in the
planar counterexample search.  Its entire unpaid invoice is one component
value on one reducible fibre.

Let `k` be an algebraically closed field of characteristic zero and let
`m>=1` be an integer.  In `k[x,y]` put

```text
A=1+xy,
B=1+x^(2m+1) A^m,
U=xAB,
P=((m+1)B-1)/(mx),
W=UP=AB((m+1)B-1)/m.                                 (1)
```

Then

```text
J(U,P)=1,                         J(U,W)=U.            (2)
```

The zero fibre of `U` is the disjoint union of the three smooth irreducible
curves

```text
x=0,                              A=0,                 B=0,              (3)
```

and the zero-component spectrum of the log-canonical pair is

```text
Spec_0(U,W)=(1,0,0)               on (x,A,B).          (4)
```

All rational response-one mates of `U` are exactly

```text
P+H(U),                           H in k(T),            (5)
```

and none is polynomial.  Thus `U` is a smooth noncoordinate polynomial
with a rational, but no polynomial, constant-Jacobian mate.

The rational map `(U,P)` has exact function-field degree

```text
[k(x,y):k(U,P)]=m+1.                                      (6)
```

Writing target coordinates as `(L,p)` and using an indeterminate `T`, its
minimal polynomial is

```text
f_m(T)=T^(m+1)-mpT+mL.                                  (7)
```

Its discriminant is

```text
disc_T(f_m)
 =(-1)^(m(m+1)/2) m^m
   [(m+1)^(m+1)L^m-m^(m+1)p^(m+1)].                    (8)
```

For `m>=2` the extension in `(6)` is not Galois.  For `m=2`, specifically,

```text
f_2(T)=T^3-2pT+2L,
disc_T(f_2)=-4(27L^2-8p^3),                            (9)
```

and the Galois closure has group `S_3`.  The branch curve
`27L^2=8p^3` is an `A_2` cusp.  This cubic row is the centerpiece of the
tower: all generic cover invariants have crossed the degree-three boundary,
while the polynomial obstruction remains exactly the single axis entry in
`(4)`.

## 1. The all-degree Jacobian cancellation

The useful coordinates are `(x,A)`.  Since `A=1+xy`,

```text
J_(x,y)(F,G)=x J_(x,A)(F,G).                            (10)
```

Set `eta=B-1=x^(2m+1)A^m`.  At fixed `A`, the four derivatives in `(1)`
are

```text
U_x=A[1+2(m+1)eta],          U_A=x[1+(m+1)eta],
P_x=[-1+2(m+1)eta]/x^2,      P_A=(m+1)eta/(xA).        (11)
```

The two quadratic terms and the two linear terms cancel separately:

```text
x(U_xP_A-U_AP_x)
 =(m+1)eta[1+2(m+1)eta]
  -[1+(m+1)eta][-1+2(m+1)eta]
 =1.                                                      (12)
```

This proves the first identity in `(2)`.  The second follows without another
expansion:

```text
J(U,UP)=U J(U,P)=U.                                    (13)
```

The exponent `2m+1` is forced inside this ansatz, rather than guessed from
the first members.  For an arbitrary positive integer `q`, define

```text
B_q=1+x^qA^m,
U_q=xAB_q,
P_q=((m+1)B_q-1)/(mx).                                (14)
```

The same chart calculation gives the exact hostile identity

```text
J(U_q,P_q)-1
 =-[(m+1)(q-2m-1)/m] x^q A^m B_q
 =-[(m+1)(q-2m-1)/m] x^(q-1)A^(m-1)U_q.              (15)
```

Thus the nonconstant response vanishes if and only if `q=2m+1`.  The tower
exponent is the unique balance between the `m`-fold `A` charge and the two
Jacobian weights carried by `x`.

## 2. The zero fibre is smooth, irreducible by pieces, and disjoint

The factors in `(3)` do not meet.  On `x=0`, both `A` and `B` equal one;
on `A=0`, `B` equals one.  It remains to justify that the displayed factors
are genuinely the irreducible components.

The only nonobvious factor is `B`.  Localizing at `x` identifies

```text
k[x,x^(-1),y]=k[x,x^(-1),A].                          (16)
```

After also localizing at `A`, the polynomial `B` is the Laurent binomial

```text
1+x^(2m+1)A^m.                                        (17)
```

The exponent vector `(2m+1,m)` is primitive because
`gcd(2m+1,m)=1`.  A unimodular monomial change therefore sends `(17)` to
`1+z`, so it is prime in the two-variable Laurent ring.  If a factor of
`B` in the ring in `(16)` became a Laurent unit, it would be a monomial
`c x^a A^b`.  The nonzero constant term of `B` as a polynomial in `A`
forces `b=0`; hence it was already a unit in `(16)`.  Thus `B` is
irreducible there.  Finally, a factor in `k[x,y]` that becomes a unit after
localizing at `x` is a power of `x` up to a scalar, but `B|_(x=0)=1`.
Therefore `B` is irreducible in `k[x,y]`.

There are no critical points of `U`.  On `x!=0`, a critical point of `U`
would make the left side of `J(U,P)=1` vanish.  On the remaining line,

```text
U_x|_(x=0)=1.                                          (18)
```

Consequently `U` is smooth, and its three disjoint irreducible components
are smooth as well.  More concretely, the first component is `A^1`, the
second is `G_m`, and `(17)` realizes the third as a translate of a primitive
one-dimensional subtorus, hence again as `G_m`.  In particular, every
component has constant field `k`.  Since a coordinate polynomial has an
irreducible zero fibre, the three-component fibre also proves that `U` is
not a coordinate.

Reduction of the polynomial `W` in `(1)` on the three components gives

```text
W|_(x=0)=1,                 W|_(A=0)=0,                W|_(B=0)=0.       (19)
```

This proves `(4)`.  Equivalently, because `P=W/U` and every factor of `U`
has multiplicity one, `(1,0,0)` is the vector of first principal
coefficients of `P` relative to the common vertical parameter `U`.

## 3. A birational source chart exposes every cover degree

Put

```text
s=xA,                             t=xs=x^2A.           (20)
```

Since

```text
x^(2m+1)A^m=t^(m+1)/s,                                     (21)
```

the two rational target coordinates become

```text
U=s+t^(m+1),
P=s/t+[(m+1)/m]t^m=U/t+t^m/m.                         (22)
```

Conversely, on a dense open set,

```text
s=U-t^(m+1),       x=t/s,       A=s^2/t,       y=(A-1)/x.                (23)
```

Hence

```text
k(x,y)=k(U,t).                                         (24)
```

Equation `(22)` yields the relation `(7)`.  It is irreducible over
`k(p,L)`.  Indeed, a factorization of the monic polynomial `f_m` over that
fraction field descends, by Gauss's lemma, to a factorization into positive-
degree monic polynomials in `k[p,L][T]`.  Specializing `p=0` would factor

```text
T^(m+1)+mL,                                             (25)
```

but `(25)` is Eisenstein at the prime `L`.  This contradiction proves
irreducibility.  Equations `(23)--(25)` prove the degree formula `(6)`.

The same chart identifies the complete Hamiltonian constant field.  A
direct calculation gives

```text
J(t,U)=xt !=0.                                          (26)
```

For `R in k(U,t)`, therefore,

```text
J(R,U)=xt (partial R/partial t).                        (27)
```

Characteristic zero now gives

```text
ker_(k(x,y)) J(-,U)=k(U).                              (28)
```

Because `J(U,P)=1`, subtracting two rational response-one mates and using
`(28)` proves the torsor statement `(5)`.  It also shows that the generic
`U`-fibre is geometrically integral and rational, with parameter `t`; the
degree in `(6)` is target-cover degree, not hidden genus.

## 4. The discriminant tower and the first non-Galois row

The standard discriminant formula for the monic trinomial `T^n+aT+b`, with
`n=m+1`, `a=-mp`, and `b=mL`, gives exactly `(8)`.  Its nonconstant branch
factor

```text
D_m=(m+1)^(m+1)L^m-m^(m+1)p^(m+1)                     (29)
```

is irreducible: the exponent vector `(m,m+1)` is primitive, so the same
Laurent-binomial argument as in Section 2 applies.

The normalization of `(29)` is visible without elimination.  At its generic
point set

```text
p=(m+1)r^m/m,                  L=r^(m+1).              (30)
```

Then `r` is a double root of `f_m`, since

```text
f_m(r)=f_m'(r)=0,             f_m''(r)=m(m+1)r^(m-1) !=0.               (31)
```

It is the unique repeated root generically.  Indeed, another one would have
the same `m`th and `(m+1)`st powers as `r`, and coprimality forces equality;
equivalently `r=(m+1)L/(mp)` in the function field of `(29)`.  Thus inertia
in the Galois closure contains one transposition.

Irreducibility of `f_m` makes the monodromy action transitive.  If the
degree-`m+1` extension itself were Galois, this action would be regular, so
every nonidentity element would act without fixed points.  For `m>=2`, the
transposition in `(31)` fixes `m-1` of the roots, a contradiction.  This
proves the asserted non-Galoisness.  The quadratic row `m=1` is, as expected,
Galois.

For `m=2`, an irreducible cubic has transitive Galois group `A_3` or `S_3`.
The discriminant in `(9)` is not a square because its irreducible branch
factor occurs to the first power.  Its Galois closure therefore has group
`S_3`.  After scaling the target coordinates, `27L^2-8p^3=0` is the standard
`v^2=u^3` cusp, proving the `A_2` assertion.  More generally, `(29)` gives
the coprime `(m,m+1)` branch cusp for every `m>=2`; no full symmetric-group
claim is made beyond the cubic row.

## 5. The one-entry vertical obstruction is complete

The rational function in `(1)` has the transparent form

```text
P=1/x+[(m+1)/m]x^(2m)A^m.                             (32)
```

Thus it has no horizontal pole and its only affine pole is on `x=0`, inside
the fibre `U=0`.  Suppose `P+H(U)` were polynomial.  A finite pole of `H`
at a value `lambda!=0` would create a pole on `U=lambda`, where `(32)` is
regular, so `H` can have no such pole.  At `lambda=0`, a pole of order at
least two would create the same higher-order pole on all three reduced
components, while `P` has order at most one.  A simple correction `c/U`
changes the principal coefficient vector to

```text
(1+c,c,c),                                             (33)
```

which is never zero.  Hence no `H` regularizes `P`.  By `(5)`, this excludes
every rational mate that could be polynomial, not merely a chosen ansatz.
It is also the direct instance of THM-3770's necessary-and-sufficient
vertical equalizer.

The multiplication

```text
P |-> W=UP                                             (34)
```

is an affine blowdown: it clears the axis pole and changes constant response
to the polynomial log-canonical law `J(U,W)=U`.  It does not erase the debt;
the lost regularity survives exactly as `(4)`.  Two of three critical
components already collapse to `(0,0)`, while the axis maps to `(0,1)`.
Equalizing that last address without making `U` critical would, by THM-3770,
turn `(W-c)/U` into an actual polynomial Keller mate.

## 6. Sharp boundaries and scope

Three controls show what is load-bearing.

1. Replacing `2m+1` by any `q!=2m+1` produces the explicit nonconstant
   response in `(15)`.
2. Replacing `A=1+xy` by `A=xy` can equalize the apparent axis value only by
   making the origin critical: both first derivatives of `U` vanish there.
3. Deleting the nonlinear term, so `B=1`, recovers the two-component seed
   `U=x(1+xy), P=1/x`; it preserves rational exactness but loses the new
   degree tower.

At `m=1`, `(7)` is the conic `T^2-pT+L=0` with discriminant `p^2-4L`, so the
tower contains the earlier quadratic geometry.  At `m=2`, it crosses to the
cubic `S_3` cusp in `(9)`.  The exact companion verifies the symbolic all-m
identities, seven direct members, five smoothness resultants, fourteen
nearest off-exponent hostiles, and both degeneration mechanisms.

Nothing here asserts that `(U,P)` is polynomial: equation `(32)` proves the
opposite.  Nothing here is a planar Jacobian counterexample.  The theorem
instead isolates a much narrower remaining construction problem: preserve
the cubic rational cover and its smooth noncoordinate first coordinate while
changing the single exceptional value `(0,1)` to `(0,0)` without introducing
a critical point, a horizontal pole, or a new component mismatch.  **QED.**
