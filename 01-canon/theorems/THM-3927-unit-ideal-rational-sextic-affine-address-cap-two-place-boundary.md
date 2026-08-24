---
id: THM-3927
title: "Unit-ideal rational sextic passes the affine address cap but has two infinity places"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. Over an
  algebraically closed field of characteristic zero, an explicit binary
  cubic over k[A,C] has unit coefficient ideal, irreducible S3 generic
  fibre, smooth normal Delone--Faddeev surface, scalar units, class group
  Z^2, and no scalar-unit index-form value. Its irreducible rational sextic
  discriminant has exactly two affine collision fibres, two addresses in
  each, both A5 contacts, so the target-branch affine address cap is two.
  Nevertheless its affine normalization is P1 minus two points, so
  THM-3841's Jelonek gate rejects it. Independently, the ramification curve
  identifies both affine address pairs inside the completion and is not
  unibranch, violating THM-3920. Its class and one vertical prime form a
  class-group basis, but their natural etale deletion has compactly
  supported Euler characteristic 5. Compressing the target by B=AC produces
  a rational one-place octic;
  after the forced index-A saturation, its normal maximal cubic order is
  explicitly monogenic. This is a sharp address-cap-positive near miss, not
  a Keller atlas or a counterexample to JC(2).
source: jc_zero_debt_lift / post-THM-3920 binary-cubic address-cap search, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (jc-cohn3709, 2026-08-23). The audit
  reconstructed generic irreducibility and nonmonogenicity, the projective
  normalization, both infinity valuations, the complete collision/genus
  ledger, the Nagata class and unit calculation, different valuations,
  ramification identifications, Euler ledger, and compressed maximal
  overorder. It repaired a substantive mechanism/scope error: two
  projective infinity places invoke THM-3841, whereas THM-3920 applies
  independently because the ramification curve identifies the two affine
  branch pairs and is non-unibranch. The assertion-free companion
  byte-matches in normal and optimized mode; frozen output, raw hashes, and
  documentation checks pass.
depends_on:
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
  - THM-3841-deleted-ramification-three-puncture-jelonek-nonentry
  - THM-3920-affine-plane-boundary-unibranch-depressed-cubic-chart-obstruction
  - THM-3922-affine-plane-open-boundary-basis-class-group-obstruction
related:
  - THM-3906-degree-six-common-zero-normal-cubic-two-place-boundary
  - THM-3907-unit-ideal-nonmonogenic-cubic-six-place-boundary
  - THM-3926-unit-ideal-cubic-primitive-class-genus-two-boundary
script: 04-computation/jc2_unit_ideal_rational_sextic_address_cap_thm3927.py
output: 05-knowledge/results/jc2_unit_ideal_rational_sextic_address_cap_thm3927.out
script_sha256: 90cc89a19a7593bded7d07216accbe12fc21cbedf5eda6ef9f12491abc73a533
output_sha256: b35a13f3b8ff83e226a9d6c205fafcc0feb77d183399ce90bac6cb91989811d6
semantic_sha256: 67524a2bbcb3b40b473373aaffba08ca57c23978935ef9937e070e34f11bb034
hash_basis: raw LF bytes
---

# THM-3927 -- the affine cap passes, but infinity still splits

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of characteristic zero. Compactly supported
Euler characteristics are asserted after specialization to `k=C`.

Let `R=k[A,C]` and consider the binary cubic

```text
Psi=A(1+27A)U^3+C U^2V-(1+24A)UV^2+8A V^3.                (1)
```

Let `S` be its Delone--Faddeev cubic algebra and `X=Spec S`. This theorem
shows that `(1)` clears three invoices that had previously failed in
different examples: it is normal and genuinely nonmonogenic, its branch is
rational, and no affine target-branch fibre has more than two normalization
addresses. It then fails three separate global invoices: polynomial
uniruledness at infinity, boundary unibranchness inside the completion, and
the Euler characteristic of the natural deletion.

## 1. Unit coefficient ideal, nonmonogenicity, and the S3 order

Write the coefficient row of `(1)` as

```text
(a,b,c,d)=(A(1+27A), C, -(1+24A), 8A).                    (2)
```

It is unimodular because

```text
c+3d=-1.                                                   (3)
```

This does not make the cubic order monogenic. If the index form represented
a scalar unit, there would be `U,V in k[A,C]` and `lambda in k*` with
`Psi(U,V)=lambda`. Specializing at `A=0` gives

```text
U V(CU-V)=lambda                                           (4)
```

in `k[C]`. All three factors on the left would be units. The first two then
make `U,V` nonzero constants, while `CU-V` remains nonconstant, a
contradiction. By the cubic index-form criterion of THM-3801, `S` is not
monogenic over `R`.

Dehomogenize at `V=1`:

```text
f(t)=A(1+27A)t^3+C t^2-(1+24A)t+8A.                       (5)
```

Solving for the target-linear coefficient gives

```text
C=-[A(1+27A)t^3-(1+24A)t+8A]/t^2.                         (6)
```

The numerator and denominator in `(6)` are coprime over `k(A)` and the
rational map in `t` has degree three. Hence `(5)` is irreducible over
`k(A,C)`. Its discriminant is

```text
Delta=
 -1259712A^6+1399680A^5-93312A^4C+240192A^4
 -7344A^3C+14688A^3+576A^2C^2-144A^2C+396A^2
 -32AC^3+48AC^2+4A+C^2.                                  (7)
```

Moreover

```text
Delta(A,0)=-4A(27A+1)
  (11664A^4-13392A^3-1728A^2-72A-1),                     (8)
```

so `Delta` has an odd prime valuation and is not a square in `k(A,C)`.
Thus the generic Galois group is `S3`.

In the basis `(1,omega,theta)`, the Delone--Faddeev multiplication is

```text
omega^2=C omega-A(1+27A)theta+A(1+27A)(1+24A),
omega theta=-8A^2(1+27A),
theta^2=-8AC+8A omega+(1+24A)theta.                        (9)
```

The ideal of these three relations, together with all `2 x 2` minors of
their `3 x 4` Jacobian matrix, has reduced Groebner basis `[1]`. The
Delone--Faddeev algebra is finite flat of rank three and its generic fibre
is a field, so `(9)` defines an integral surface. The Jacobian calculation
therefore proves

```text
X is smooth, hence S is the normal maximal R-order.         (10)
```

## 2. Rational sextic normalization and two places at infinity

Eliminating `C` from `f=f_t=0` gives the repeated-root incidence

```text
27A^2t^3+A(t^3+24t-16)+t=0.                               (11)
```

Its discriminant as a quadratic in `A` is

```text
(t^2-8t+4)(t^2+4t-8)^2.                                  (12)
```

Parametrize the nonsquare conic factor by

```text
t=-4(s+2)/((s-1)(s+1)),
y=-2(s^2+4s+1)/((s-1)(s+1)),
y^2=t^2-8t+4,
s=(y-2)/t,
y=[54At^3+t^3+24t-16]/(t^2+4t-8).                       (13)
```

The last two formulas recover `(y,s)` rationally from a generic repeated-root
incidence point. Equations `(5),(11)` then give the birational normalization

```text
A(s)=-(s-1)^2(s+2)/108,
C(s)=(s-1)(s+1)(s^4+6s^3-22s-21)/(72(s+2)).               (14)
```

The recovery formula in `(13)` proves generic injectivity. In homogeneous
coordinates `[S:T]` on `P1`, `(14)` extends to

```text
[S:T] |-> [
 -2(S-T)^2(S+2T)^2T^2 :
  3(S-T)(S+T)(S^4+6S^3T-22ST^3-21T^4) :
  216(S+2T)T^5 ].                                         (15)
```

The three sextics in `(15)` have no common zero. Their image lies on the
homogenization of `(7)`, and the recovery formula makes its degree six.
Consequently `(7)` is the irreducible equation of a rational plane sextic.

At `Z=0`, that equation is

```text
Delta_h(A,C,0)=-1259712A^6,                               (16)
```

so the projective curve has only the support `[0:1:0]` at infinity.
However `(15)` has two preimages of that point:

```text
s=-2,                         s=infinity.                  (17)
```

In the target chart `x=A/C, z=Z/C`, with `epsilon=s+2` at the first place
and `r=1/s` at the second, the leading terms are

```text
s=-2:       x=(2/9)epsilon^2+...,     z=-(8/3)epsilon+...,
s=infinity: x=-(2/3)r^2+...,          z=72r^5+....          (18)
```

Thus the first branch is smooth, the second is a `(2,5)` cusp, and their
tangent lines are distinct. There is one projective infinity point but
exactly two normalization places. Equivalently, the affine normalization is
`P1` minus two points rather than `A1`. THM-3841's deleted-divisor/Jelonek
argument therefore forbids this curve as a planar nonproperness component.
This is a projective-puncture obstruction, not an application of
THM-3920's affine finite-fibre address cap.

## 3. Exhaustive affine collisions: the exact cap is two

Let `r != s` be affine normalization parameters with the same `A`-value,
and put

```text
u=r+s,                    v=rs.                            (19)
```

Dividing `A(s)-A(r)` by `s-r` gives

```text
s^2+sr+r^2=3,             v=u^2-3.                        (20)
```

Set

```text
N(s)=(s-1)(s+1)(s^4+6s^3-22s-21).                        (21)
```

On the same-`A` correspondence `(20)`, the numerator of the `C`-difference
has the exact reduction

```text
[N(s)(r+2)-N(r)(s+2)]/(s-r)
       =-(u^2+2u-2)^3.                                    (22)
```

The omitted denominator is harmless: `(s+2)(r+2)` reduces to `3` at the
roots of `u^2+2u-2`. Therefore all off-diagonal affine collisions, and no
others, are the two roots

```text
u^2+2u-2=0.                                                (23)
```

For either root, the two addresses are the roots of

```text
z^2-u z+u^2-3=0,                                          (24)
```

whose discriminant reduces to `6(u+1) != 0`. The two target points are
distinct because

```text
A=(u-2)/36.                                                (25)
```

Equivalently, the four collision addresses are the four simple roots of

```text
s^4+2s^3-10s-11.                                          (26)
```

This also excludes a hidden third address in either fibre. Hence

```text
maximum number of affine normalization addresses = 2.      (27)
```

The affine parametrization is immersive: `A'=-(s-1)(s+1)/36` is coprime
to the numerator of `C'`. Since `(23)` is squarefree and `(22)` is its
cube, the two smooth branches in each collision fibre have intersection
multiplicity three. Both singularities are therefore of type `A5`, with
delta invariant three.

This list is exhaustive even as a singularity list. A plane sextic has
arithmetic genus ten. The two affine `A5` points contribute `3+3`; at
infinity the `(2,5)` cusp contributes two and its transverse intersection
with the smooth branch contributes two. Thus

```text
3+3+2+2=10=p_a,                                            (28)
```

leaving no unrecorded singularity or collision.

## 4. The class invoice passes, but the natural open has Euler five

The element `theta` obeys

```text
P(T)=T^3-(1+24A)T^2+8ACT+64A^3(1+27A).                   (29)
```

On `A theta !=0`, equations `(9),(29)` solve both missing coordinates:

```text
C=-[theta^3-(1+24A)theta^2+64A^3(1+27A)]/(8A theta),
omega=-8A^2(1+27A)/theta.                                 (30)
```

Hence

```text
S_(A theta)=k[A,A^(-1),theta,theta^(-1)],                 (31)
```

a factorial two-torus. Its four boundary primes are

```text
P00=(A,omega,theta),
P01=(A,omega,theta-1),
PC =(A,theta,omega-C),
Q  =(A+1/27,theta,omega-C).                               (32)
```

Their complete divisor relations are

```text
div(A)    =P00+P01+PC,
div(theta)=P00+2PC+Q.                                     (33)
```

For the multiplicity-two entry, `omega theta=-8A^2(1+27A)` and `omega=C`
is generically a unit at `PC`. At `P00`, writing `theta=Au` gives the two
residual roots `u=8C` and `u=0`, separating `P00` from `PC`; at `Q` the
factor `1+27A` is simple. These observations also show that `(32)` lists
every component of the chart boundary.

Nagata localization presents `Cl(S)` by the four prime classes modulo the
two rows

```text
          P00 P01 PC Q
A          1   1  1 0
theta      1   0  2 1.                                    (34)
```

The Smith invariants are `(1,1)`. With `g=[PC]` and `q=[Q]`, one obtains

```text
Cl(S)=Z g direct_sum Z q,
[P00]=-2g-q,       [P01]=g+q.                             (35)
```

The same boundary rows determine the global units. Every unit restricts on
the torus to `lambda A^m theta^n`; its coefficients at `P01` and `Q` force
`m=n=0`. Thus

```text
S*=k*.                                                     (36)
```

On the chart `(31)`, the different is generated by

```text
D_theta=P'(theta)=3theta^2-2(1+24A)theta+8AC.             (37)
```

Let `E` be the irreducible ramification divisor. Put `theta=Au`. Then

```text
P(Au)/A^2=-u^2+8Cu+A(u^3-24u^2+64+1728A),                (38)
```

and `D_theta/A` restricts to `-8C` on the `u=8C` branch `P00` and to `8C`
on the `u=0` branch `PC`. The different is a unit at `P01` and is
generically a unit at `Q`. Therefore

```text
div(D_theta)=E+P00+PC,
[E]=g+q=[P01].                                             (39)
```

In particular `([E],[Q])` is a basis of `Cl(S)`. Deleting `E union Q`
produces a smooth quasi-finite etale generically degree-three open

```text
U=X minus (E union Q)                                      (40)
```

with `Cl(U)=0` and `Gamma(U,O_U)^*=k*`. Thus this natural deletion passes
the full THM-3922 class-and-unit invoice.

It is still not an affine plane. The normalization of `E` is
`P1 minus {-2,infinity}`, with coordinates

```text
theta(s)=(s-1)^3(s+1)/54,
omega(s)=(s-2)(s-1)(s+1)(s+2)^2/108.                     (41)
```

At precisely the two collision pairs `(23)`, the extra coordinates also
coincide: on `(20)` their difference quotients reduce to

```text
[theta(s)-theta(r)]/(s-r)
  =-(u-2)(u^2+2u-2)/54,
[omega(s)-omega(r)]/(s-r)
  =-(u-1)(u+1)(u^2+2u-2)/108.                             (42)
```

Thus `E` is not smooth: it identifies two distinct pairs of its
normalization points. Its two singular points, indexed by the roots of
`u^2+2u-2`, are

```text
(A,C,omega,theta)=
((u-2)/36, (5u-4)/12, (5u-4)/36, (2u-1)/9).               (43)
```

Each identification gives two normalization branches through one point of
the irreducible boundary curve `E`. Thus `E` is not unibranch, and
THM-3920 independently forbids an affine-plane open in this completion.
This source-boundary failure is stronger than merely observing that the
target branch obeys the numerical address cap `(27)`.

Also `Q meet E` is the single point `(-1/27,0,0,0)`, corresponding to
`s=-1`.

The torus `(31)` has Euler characteristic zero. The three affine lines over
`A=0` have only the intersection `P00 meet PC` at `C=0`, so their union has
Euler characteristic two; the disjoint line `Q` contributes one. Hence

```text
chi_c(X)=3.                                                (44)
```

The normalization of `E` has Euler characteristic zero and each of its two
pair-identifications lowers Euler characteristic by one, giving
`chi_c(E)=-2`. Inclusion--exclusion with the unique point `E meet Q` now
gives

```text
chi_c(U)=3-[(-2)+1-1]=5 !=1=chi_c(A2).                    (45)
```

The natural class-basis deletion is therefore Euler-hostile even before
the two-place infinity obstruction is invoked.

## 5. One-place compression pays for infinity with monogenicity

The most instructive attempted repair is the target compression

```text
B=AC.                                                      (46)
```

It preserves the function field away from `A=0`. Clearing the denominator
in `(1)` gives the binary row over `k[A,B]`

```text
(a0,b0,c0,d0)=
(A^2(1+27A), B, -A(1+24A), 8A^2).                         (47)
```

The reduced branch equation is

```text
Delta_*=A^2 Delta(A,B/A)
=-1259712A^8+1399680A^7+240192A^6-93312A^5B+14688A^5
 -7344A^4B+396A^4-144A^3B+4A^3+576A^2B^2
 +48AB^2-32B^3+B^2.                                      (48)
```

Its polynomial normalization is

```text
A(s)=-(s-1)^2(s+2)/108,
B(s)=-(s-1)^3(s+1)(s^4+6s^3-22s-21)/7776.                (49)
```

The coordinate degrees are `(3,8)`. Homogenization has no base point and
`Z=T^8`, so the octic `(48)` has exactly one normalization place at
infinity, namely `s=infinity`. The compression really repairs the
two-place defect.

But the order from `(47)` has discriminant

```text
disc(S0)=A^2 Delta_*.                                      (50)
```

It is not maximal. If `(1,omega,theta0)` is its standard basis, put
`e=theta0/A`. The index-`A` overorder has multiplication

```text
omega^2=A^3(1+27A)(1+24A)+B omega-A^3(1+27A)e,
omega e=-8A^3(1+27A),
e^2=-8B+8omega+(1+24A)e.                                  (51)
```

The last relation solves

```text
omega=[e^2-(1+24A)e+8B]/8,                                (52)
```

and `e` obeys the monic cubic

```text
e^3-(1+24A)e^2+8Be+64A^3(1+27A)=0.                       (53)
```

Its discriminant is `64 Delta_*`. The hypersurface `(53)` is irreducible
by the same degree-three rational-map argument as `(6)`. Its singular locus
is only `(A,B,e)=(0,0,0)`: `partial P/partial B=8e` first forces `e=0`,
`partial P/partial e` then forces `B=0`, and

```text
gcd(A^3(1+27A), d/dA[A^3(1+27A)])=A^2                   (54)
```

forces `A=0`. A hypersurface is Cohen--Macaulay, and this singular point
has codimension two, so Serre's criterion makes `(53)` normal. It is
therefore the maximal order in the compressed generic field.

Equivalently its index form is

```text
I(x,y)=-A^3(1+27A)x^3-Bx^2y+(1+24A)xy^2-8y^3,            (55)
```

and

```text
I(0,1)=-8 in k*.                                           (56)
```

Thus the very saturation that removes the artificial index divisor makes
the one-place maximal order monogenic. The compression does not preserve
the desired nonmonogenic object; it exchanges the infinity-place debt for a
generator.

## 6. Exact scope

This example proves that rationality plus the affine address cap is not the
same as the one-place condition. It also proves that a primitive
ramification class completing to a boundary basis is not sufficient for an
affine-plane chart. The missing data live at two different boundaries:
normalization places at projective infinity, and the singular/equivariant
geometry of the deleted ramification curve.

The construction is not a Keller map. No pair of polynomial coordinates on
an affine source, constant Jacobian identity, or affine-plane atlas is
produced. THM-3841 rejects the two-puncture target branch; THM-3920
independently rejects the non-unibranch ramification boundary; and the
one-place compression has a monogenic maximal order in a different design
class. Consequently this theorem neither supplies a counterexample to
`JC(2)` nor proves `JC(2)`. **QED.**
