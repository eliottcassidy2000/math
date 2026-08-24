---
id: THM-3915
title: "Rational decic cube resolvent, normal nonmonogenic cubic, and etale Euler tariff"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  An explicit rational
  one-place decic has finite packet one ordinary eightfold point, two A2
  cusps, and two A5 contacts.  Its quadratic double plane has class group
  Z^3 direct-sum Z/3: the three-class is a mixed boundary--ADE class, while
  the pure split-boundary residue is primitive.  An explicit irreducible S3
  cubic power order has discriminant 432 A^10 Delta; a global index-A^5
  basis gives its normal finite-flat, rational, globally nonmonogenic maximal
  order with discriminant 6912 Delta.  Thus A=0 is pure order debt and the
  sole reduced field-discriminant divisor is the decic.  The maximal etale
  locus has compact-support Euler characteristic 14 and is not A2; this
  Euler argument alone does not exclude a proper plane open.  Subsequent
  THM-3916 proves that the contracted common-zero divisor has genus two and
  excludes every same-field plane open on which the target functions form a
  Keller pair.  This candidate is closed; JC(2) remains open.
source: decic_resolvent_scout / rational cube-discriminant lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS.  The audit rederived the complete
  branch/collision packet, conic-bundle Picard quotient, global integral
  basis and multiplication table, index form, normality bridge, rational
  chart, singular-fibre census, and Euler ledger.  In particular, it checked
  that discriminant exponent one gives maximality only after the displayed
  finite-free/S2 and irreducibility gates, and that the Euler calculation
  excludes only the maximal etale locus itself, not a proper plane open.
  Normal and optimized replays byte-match the frozen output, raw hashes and
  documentation pass, and no repair was required.  The earlier multi-path
  cross-check also derived the branch/collision and local DVR packets
  were derived independently from the parametrization and from the conic
  cube identity.  A separate conic-bundle lattice audit recovered Picard
  rank 20, the removed-lattice determinant/Smith packet, the exact endpoint
  quotient Z^3 plus Z/3, and the mixed class.  A separate global-order audit
  produced the same integral basis, multiplication table, binary index form,
  A^5 debt, rational chart, and singular-fibre Euler packet.  The companion
  freezes all of these in 91 active gates; normal and optimized replays
  byte-match the frozen output and documentation passes.  The genus-two
  contracted-divisor idea is intentionally left for a successor audit and
  is not used to claim a no-atlas theorem here.  THM-3916 subsequently pays
  that invoice; it is a successor, not a dependency of this proof.
related:
  - THM-3912-even-one-place-split-boundary-a2-three-torsion-design-sieve
  - THM-3913-moving-triple-root-one-place-decic-normal-nonmonogenic-cubic
  - THM-3914-decic-boundary-three-class-degree-one-isotropic-divisor
  - THM-3916-positive-genus-collapsed-valuation-keller-obstruction
  - THM-3917-quintic-parameter-rational-collapsed-cubic
script: 04-computation/jc2_rational_decic_cube_resolvent_index_debt_thm3915.py
output: 05-knowledge/results/jc2_rational_decic_cube_resolvent_index_debt_thm3915.out
script_sha256: 8408e35292a674be3731fd7c95ecfdc7d07d5324c81d8803d1e23922e8884341
output_sha256: 5e1d4cd168597e8975580c00cb1e6cfd06c6d1d8a9657f5b57a056983866ecf1
semantic_sha256: 238246047160a73f2de01d4ff62109d7cac03434ad31fdc58ae0378054030c6b
hash_basis: raw LF bytes
---

# THM-3915 -- a rational decic pays the three-class and normal-order invoices

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an
algebraically closed field `k` of characteristic zero.  Put

```text
G=t^8+6t^6+6t^4-19t^2-24,
A=tG,                         C=(t^2+1)G.                    (1)
```

The image is the irreducible rational decic `Gamma=V(Delta)`, where

```text
Delta=-A^10+30A^8C-16A^8+5A^6C^3-309A^6C^2
      -9A^4C^5-246A^4C^4-3A^2C^7-29A^2C^6+C^9+24C^8.      (2)
```

Its projective closure has one smooth place at infinity and contact ten
with the infinity line.  Its complete finite singularity packet is

```text
one ordinary 8-fold point + two A2 + two A5.                 (3)
```

Let

```text
Q_2=Spec k[A,C,W]/(W^2-Delta).                              (4)
```

Then

```text
Q_2^*=k^*,                    Cl(Q_2)=Z^3 direct-sum Z/3.    (5)
```

The unique three-torsion line is not the pure degree-ten boundary class of
THM-3914.  It mixes that boundary residue with all four finite ADE radical
directions.

There is also an explicit cubic.  Define

```text
P=(C^2-4A^2)(C^2+A^2)^2,

Q=2A^10-C^9+3A^2C^7+9A^4C^5-5A^6C^3-30A^8C,               (6)

f(Z)=Z^3-3PZ-2Q.                                            (7)
```

The polynomial `f` is irreducible over `k(A,C)` and has Galois group `S3`.
Its power-basis order has discriminant

```text
disc_Z(f)=432 A^10 Delta.                                   (8)
```

The factor `A^10` is exactly square-index debt.  The integral closure is an
explicit finite-free normal cubic domain `B` with

```text
disc(B/k[A,C])=6912 Delta,                                  (9)
```

and the inclusion `k[A,C,Z]/(f) subset B` has index ideal `(A^5)`.  The
maximal order `B` is globally nonmonogenic.  Its fraction field is rational.

Finally, if `X=Spec B` and `U_et subset X` is the maximal open on which
`X -> A^2_(A,C)` is etale, then

```text
e_c(U_et)=14.                                               (10)
```

Thus `U_et` itself is not an affine plane.  Inside this theorem, this is an
exact boundary tariff rather than a no-atlas theorem: a hypothetical plane
atlas could have been a proper open subset with complement Euler
characteristic `13`.  Subsequent THM-3916 proves that the contracted divisor
has genus two and excludes a same-field affine-plane Keller model for these
target functions.  No Keller map is constructed and `JC(2)` remains **OPEN**.

## 1. Normalization, the forbidden third cusp, and the genus ledger

The exact resultant is

```text
Res_t(A-tG,C-(t^2+1)G)=Delta.                               (11)
```

The homogeneous normalization map is

```text
[T:S] |-> [TS G_h:(T^2+S^2)G_h:S^10].                     (12)
```

It has no base point.  Away from `A=0`, its ratio is

```text
C/A=t+t^-1.                                                 (13)
```

Thus two generic parameters with the same image are either equal or
reciprocal.  The reciprocal equality is not an identity, as the collision
factor below shows, so `(12)` is generically one-to-one.  Its image has
degree ten and is irreducible; consequently `(12)` is the normalization and
`Delta` is absolutely irreducible.

There is no hidden third critical parameter.  Indeed

```text
A'=G+tG',                    C'=2tG+(t^2+1)G',              (14)
```

and the determinant of the coefficient matrix in `(G,G')` is `1-t^2`.
Away from `t=+/-1`, simultaneous vanishing would force `G=G'=0`, contrary
to

```text
disc(G)=-605291616000000 !=0.                               (15)
```

The imposed relations

```text
G(1)+G'(1)=0,                  G(-1)-G'(-1)=0               (16)
```

give the two critical targets

```text
t=1  |-> (-30,-60),            t=-1 |-> (30,-60).          (17)
```

At each address the determinant of the second and third derivative vectors
is `180000`, so both are `A2` cusps.

The reciprocal mechanism is completely visible.  Put `s=t+t^-1` and
`a(t)=tG(t)`.  Then

```text
a(t)-a(t^-1)=(t-t^-1)^3 (s^2+1)^3,                         (18)

t^10G(t)-t^8G(t^-1)
 =(t-1)^3(t+1)^3(t^4+3t^2+1)^3.                            (19)
```

The quartic in `(19)` is squarefree, avoids `G` and `t=+/-1`, and its roots
form the two reciprocal pairs with `s=+i,-i`.  In the coordinates `(s,A)`,
the two immersed graph branches differ to exact order three.  Hence their
targets

```text
(10i,-10),                         (-10i,-10)                (20)
```

are `A5` singularities, not additional cusps.

The eight simple roots of `G` map to the origin.  Their tangent directions
are distinct because two slopes `r+r^-1` agree exactly for equal or
reciprocal roots, while

```text
gcd(G(t),t^8G(t^-1))=1.                                    (21)
```

Thus the origin is ordinary eightfold.  An independent singular-ideal
elimination has final support polynomial

```text
C^11(C+10)^5(C+60)^2,                                      (22)
```

and the three fibre gcds are

```text
C=0: A^7,             C=-10: A^2+100,
C=-60: A^2-900.                                             (23)
```

This already lists the five targets.  Completeness also follows from the
genus ledger

```text
p_a(Gamma)=36
 =28 (ordinary 8-fold)+2 (two A2)+6 (two A5).               (24)
```

At infinity the degree-ten form of `(2)` is `-A^10`, while the derivative
of the homogenization in `Z` is `1` at `[0:1:0]`.  Therefore this is the
unique infinity point, it is smooth, and its line contact is ten.

## 2. The conic cube and the mixed three-class

On a pencil line `C=sA`, equation `(2)` becomes

```text
Delta(A,sA)=-A^8(A^2-u(s)A+v(s)),                           (25)

u=s(s^8-3s^6-9s^4+5s^2+30),
v=-24s^8+29s^6+246s^4+309s^2+16.                           (26)
```

The discriminant of the residual quadratic is the perfect cube

```text
u^2-4v=(s^2-4)^3(s^2+1)^6
      =((s^2-4)(s^2+1)^2)^3.                               (27)
```

The surface `Q_2` is normal.  Indeed it is a hypersurface and hence `S2`;
because `Delta` is reduced and its singular locus is the finite set found
in `(22),(23)`, the singular locus of the double plane is finite.  Thus it
is also `R1`, and Serre's criterion applies.

This explains both the local packet and the global three-class.  The pencil
gives a conic bundle on a smooth projective resolution `S` of `(4)`.  Its
four degenerate fibres have exponents

```text
n_1,n_2,n_3,n_4=3,3,6,6.                                  (28)
```

The first two contain the `A2` chains and the last two the `A5` chains.
The two infinity curves `B_+,B_-` are sections.  The conic bundle is
rational, and its Picard group is the unimodular rank-twenty lattice with
the following integral basis:

```text
B_+, F,
E_(i,1),...,E_(i,n_i-1), L_(i,n_i)       for i=1,...,4.    (29)
```

Here `F` is a fibre, and the `i`th resolved degenerate fibre is the chain

```text
L_(i,0)-E_(i,1)-...-E_(i,n_i-1)-L_(i,n_i),                 (30)
```

with `B_+` meeting `L_(i,0)` and `B_-` meeting
`L_(i,n_i)`.  In `(29)`, the first two classes have Gram matrix
`[[-4,1],[1,0]]`; each displayed vertical block has diagonal
`-2,...,-2,-1` and adjacent intersections one.  The full determinant is
`-1`.

Put

```text
w=B_+-B_-.                                                  (31)
```

Pairing with the basis `(29)` gives the exact class identity

```text
w=-9F+sum_i (sum_(j=1)^(n_i-1) j E_(i,j)+n_i L_(i,n_i)).  (32)
```

For each chain define its mod-three radical representative

```text
rho_i=sum_(j=1)^(n_i-1) (j mod 3) E_(i,j),                 (33)
```

where the residues are chosen in `{0,1,2}`.  Then

```text
rho_i^2=-2n_i,

tau=(w-rho_1-rho_2-rho_3-rho_4)/3                         (34)

   =-3F+sum_i (sum_(j=1)^(n_i-1) floor(j/3)E_(i,j)
                                  +(n_i/3)L_(i,n_i))
```

is an honest integral Picard class, with

```text
tau^2=(-18-6-6-12-12)/9=-6.                                (35)
```

This is the promised globalization.  It is necessarily mixed.  Indeed the
pure boundary vector pairs as

```text
w L_(i,n_i)=-1,                                             (36)
```

so `w/3` is not integral.  Thus the order-three permission of the boundary
Gram determinant does not globalize by itself.

For completeness, blow up the ordinary eightfold point in the base.  If
`D` is the inverse image of its exceptional line, then `D` is a genus-three
double cover of `P1`,

```text
D^2=-2.                                                     (37)
```

The boundary block is

```text
Gram(B_+,B_-)=[-4  5],                                      (38)
                    [ 5 -4]
```

and the complete removed lattice has orthogonal blocks

```text
R=[-2] orthogonal_sum K_10 orthogonal_sum A2(-)^2
      orthogonal_sum A5(-)^2.                               (39)
```

It has rank `17`, determinant `5832=2^3 3^6`, and nonunit Smith factors

```text
3,3,6,6,18.                                                 (40)
```

The abstract discriminant form has `49` isotropic order-three lines; that
number alone would grossly overstate globalization.  The actual quotient is
much smaller and exact.  In `(29)`, killing `B_+,B_-`, every middle ADE
curve, and

```text
D=B_++B_--F                                                 (41)
```

kills `B_+,F,E_(i,j)` and leaves the four far endpoints with the sole
relation

```text
3L_(1,3)+3L_(2,3)+6L_(3,6)+6L_(4,6)=0.                    (42)
```

Consequently

```text
Cl(Q_2)=Z^4/<(3,3,6,6)>=Z^3 direct-sum Z/3.                (43)
```

The image of `tau` is `(1,1,2,2)` and generates the unique torsion line.
Nondegeneracy of `(39)` also makes the unit kernel zero, proving `(5)`.

## 3. The cubic field and removal of the A5 index debt

The Cardano identity behind `(6)` is

```text
P^3-Q^2=4A^10 Delta.                                       (44)
```

It gives `(8)`.  The cubic `f` is irreducible by a cheap hostile
specialization.  At `C=0`,

```text
f=Z^3+12A^6Z-4A^10.                                       (45)
```

If this cubic had a root in `k(A)`, monicity and integral closedness of
`k[A]` would make the root a polynomial `r(A)`.  The equation also makes
`r` divide `4A^10`, so Gauss's lemma gives `r=cA^j`.  The three exponents in

```text
c^3A^(3j)+12cA^(j+6)-4A^10=0
```

are `3j,j+6,10`.  No nonnegative integer `j` makes their largest value occur
twice: the only pairwise equalities give `j=3`, `j=4`, or a noninteger, and
in the first two cases the remaining exponent is larger.  Thus `(45)` has
no root and is irreducible.  Gauss's lemma over the UFD `k[A,C]` also says
that any factorization of the monic polynomial `(7)` over `k(A,C)` has monic
factors in `k[A,C][Z]`; their degrees survive `C=0`.  Hence `(7)` is
irreducible.
Its discriminant squareclass is the nonsquare irreducible `Delta`, so its
generic Galois group is `S3`.

The normal order is explicit.  In `k(A,C)[Z]/(f)` put

```text
L=C^3-CA^2,                        N=CL-4A^4,

e_1=(Z^2+LZ-2L^2)/A^4,

e_2=(2CL^2-16A^4L-NZ-CZ^2)/A^5.                           (46)
```

Then

```text
B=k[A,C] 1 + k[A,C] e_1 + k[A,C] e_2                      (47)
```

is an algebra.  Its exact multiplication table is

```text
e_1^2=
 -24C(A^2-C^2)(A^2-12C)
 +(A^2C-12A^2-36C^2)e_1
 +A(A^2-12C)e_2,

e_1e_2=
 8A(2A^4+3A^2C^2-6A^2C-3C^4-18C^3)
 -AC(C-12)e_1
 -(A^2C+12A^2+12C^2)e_2,

e_2^2=
 -8C(4A^4+3A^2C^2+24A^2C-3C^4-72C^3)
 +(16A^2+C^3+24C^2)e_1
 +AC(C+60)e_2.                                             (48)
```

The determinant of `(1,e_1,e_2)` relative to `(1,Z,Z^2)` is

```text
-4/A^5.                                                     (49)
```

Conversely `Z,Z^2` have polynomial coordinates in `(47)`, so this is a
genuine overorder with index ideal `(A^5)`.  Equations `(8),(49)` give
`(9)`.

This also supplies the requested local DVR audit.  At the generic point of
`A=0`, where `C` is a unit, put

```text
m=C^3-CA^2-4A^4/C,
beta=(Z^2+mZ-2m^2)/A^5.                                    (50)
```

The basis `(1,Z,beta)` is closed under multiplication.  Equivalently, for
the two roots reducing to `C^3`, set `Z=m+A^5x`; after division by `A^10`
the residual equation modulo `A` is

```text
3C^3x^2-4-96/C=0,

or 3C^4x^2-4(C+24)=0.                                     (51)
```

It is separable and irreducible over `k(C)`.  Hensel separately lifts the
simple root `-2C^3`, so the maximal local algebra is one unramified linear
factor plus one unramified quadratic factor.  The field discriminant
exponent at `A` is zero and the power-order index exponent is five.  Thus
`A=0` is not a branch divisor.

Away from `Delta`, `(9)` is a unit and `B` is etale.  At the generic point
of the irreducible divisor `Delta`, the discriminant exponent is one.  An
enlargement of this cubic order over that DVR would lower the exponent by
twice a positive index length, which is impossible.  Thus `B` is the maximal
order at every height-one prime: off `Delta` it is etale over a DVR, and over
`Delta` its height-one localizations are the DVR factors of that maximal
order.  Hence `B` is `R1`.  It is finite free over the regular surface
`k[A,C]`, hence Cohen--Macaulay and `S2`.  Serre's criterion makes `B`
normal.  This proves that `(47)` is the global integral closure and that the
sole reduced field-discriminant divisor is `V(Delta)`.

## 4. Why the normal order is not monogenic

For `alpha=u e_1+v e_2`, the determinant of
`(1,alpha,alpha^2)` in the basis `(47)` is the binary index form

```text
I(u,v)=
 A(A^2-12C)u^3
 -3(A^2C+4A^2-4C^2)u^2v
 +3AC(C+12)uv^2
 -(16A^2+C^3+24C^2)v^3.                                  (52)
```

Its discriminant is `6912 Delta`, and its four coefficient polynomials have
gcd one.  Nevertheless every coefficient vanishes at `(A,C)=(0,0)`.
Therefore `I(u,v)` lies in the maximal ideal `(A,C)` for every
`u,v in k[A,C]`; it can never be a unit.  Adding a scalar to `alpha` does
not change the generated order or the index form.  Hence no element
generates `(47)` as a `k[A,C]`-algebra: the maximal order is globally
nonmonogenic.

At the origin all products in `(48)` vanish after reduction, so the fibre
is

```text
k direct-sum m,                         m^2=0,
dim_k m=2.                                                   (53)
```

This is the local embedding-dimension shadow of the same obstruction.

## 5. Rationality and the non-Keller natural chart

The cubic field pays the rationality invoice as well.  On `A!=0` put

```text
s=C/A,                         z=Z/A^3,
p(s)=(s^2-4)(s^2+1)^2,
h(s)=s^9-3s^7-9s^5+5s^3+30s.                               (54)
```

Equation `(7)` becomes

```text
z^3-3p(s)z+2h(s)-4A=0.                                    (55)
```

Thus, with

```text
F(s,z)=z^3-3p(s)z+2h(s),                                  (56)
```

one has

```text
A=F/4,                         C=sF/4,                     (57)
```

and conversely `s=C/A,z=Z/A^3`.  Therefore

```text
Frac(B)=k(s,z).                                             (58)
```

The polynomial map `(s,z) |-> (A,C)` in `(57)` has Jacobian

```text
Jac_(s,z)(A,C)=-F F_z/16.                                  (59)
```

It is not Keller.  In particular the exhibited rational chart does not
solve the atlas problem; it exposes the exact Jacobian divisor that a
different chart would have to absorb.

## 6. Exact Euler tariff and remaining atlas debt

The affine normalization of `Gamma` is `A1`, so it has compact-support Euler
characteristic one.  The ordinary eightfold point identifies eight
normalization addresses and each `A5` identifies two; the two cusps are
unibranch.  Hence

```text
e_c(Gamma)=1-(8-1)-2(2-1)=-8.                              (60)
```

There are five singular target points, so

```text
e_c(Gamma_sm)=-13,
e_c(A2 minus Gamma)=1-(-8)=9.                              (61)
```

Over the complement, `U_et -> A2 minus Gamma` is etale of degree three and
contributes `27`.  At a smooth branch point the cubic fibre has type
`k direct-sum k[epsilon]/(epsilon^2)`: its simple factor is the unique
unramified companion point.  Thus over `Gamma_sm` that companion remains in
`U_et` and contributes `-13`.

No singular target contributes an etale point.  At the origin, `(53)` is
one non-etale point.  At the other four targets `A` is a unit, so the maximal
and power orders agree; direct substitution gives `P=Q=0`, hence the fibre
is the one-point algebra `k[Z]/(Z^3)`.  Therefore

```text
e_c(U_et)=27-13=14.                                        (62)
```

This proves `U_et` is not `A2`.  The Euler computation itself does not exclude
a proper open plane: such an open would need complement Euler characteristic
`13`.  THM-3916 subsequently tests the intrinsic common-zero valuation,
proves its residue curve has genus two, and rules out the relevant Keller
plane model.  That successor result is not silently imported into the proof
of the present Euler statement.

## 7. Scope and replay

This theorem constructs a rational branch, genuine quadratic-resolvent
three-class, irreducible rational `S3` field, and normal nonmonogenic
finite-flat cubic order with sole reduced discriminant `Delta`.  It does not
construct a polynomial-plane etale atlas, a Keller pair, or a Jacobian
counterexample; THM-3916 subsequently proves that no such Keller model exists
for this candidate.  The exact three-class is mixed, so THM-3914's pure-boundary
hyperbolic plane does not describe it.  The abstract `49` isotropic lines of
the removed discriminant form do not give an index-`27` saturation; the
actual quotient `(43)` has exactly one torsion line.

Reproduce all exact algebra, order, lattice, and Euler gates with

```bash
python3 04-computation/jc2_rational_decic_cube_resolvent_index_debt_thm3915.py
python3 -O 04-computation/jc2_rational_decic_cube_resolvent_index_debt_thm3915.py
```

Both streams must byte-match the output named in the metadata.  **QED.**
