---
id: THM-3869
title: "Three-cusp square-residual Cardano peel and residual line ramification"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The
  THM-3854 three-cusp quintic has an exact Cardano square residual.  Its
  irreducible S3 cubic admits an explicit finite-free normal maximal order
  which is globally nonmonogenic, but its field discriminant still contains
  the line 9x+4 with multiplicity two.  Thus the construction passes the
  cubic-domain and nonmonogenic completion gates but remains ramified along
  an extra affine-line component.  Deleting that unique totally ramified
  divisor makes 9x+4 a nonconstant unit, so the etale locus has no dominant
  polynomial-plane atlas.  Trace shifts and every scalar change of
  the two square/cube representatives fail to remove that component.
  Nonconstant coefficient changes, a Keller atlas, and JC(2) remain OPEN.
source: jc_quartic_c3_construct / three-cusp inverse-discriminant lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root, 2026-08-23).  The audit reconstructed
  the primitive Newton segment and its factor-degree consequence; checked
  the order at every height-one prime and the finite-free S2 bridge; derived
  the corrected index determinant directly; typed both split-line Kummer
  valuations; and separately checked that the scalar-shift specialization
  and coefficient cases exhaust all scalar pairs.  It also used the local
  fundamental equality to verify the nonconstant-unit atlas obstruction.
  The companion verifies the Cardano and
  normalization identities, all 81 basis associators, the trace matrix and
  corrected binary index form, the local Eisenstein equation, the split
  quadratic-resolvent carrier, the complementary bicubic identity, trace
  invariance, and the complete scalar-shift square obstruction.  Normal and
  optimized runs must byte-match the frozen transcript.
depends_on:
  - THM-3854-integrated-three-cusp-quintic-s5-natural-completion-obstruction
  - THM-3864-integrated-three-cusp-conductor-seminormal-three-direction-gate
related:
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
  - THM-3862-russell-finite-completion-nonmonogenic-branch-contract
  - THM-3865-one-place-inverse-discriminant-resolvent-class-group
script: 04-computation/jc2_three_cusp_cardano_line_ramification_thm3869.py
output: 05-knowledge/results/jc2_three_cusp_cardano_line_ramification_thm3869.out
script_sha256: 0c7b7c75030c1f84c705e14785b3f2ef6c422fff10a12bacf8212f967a9f8e4d
output_sha256: 8f4b36e6dbf8ff9fd3ad77e18f6ed74f8ee6dfdd20279a4c9ad8cc14da4642a2
semantic_sha256: 3ba97aae53c83071a76d456a738b5eee13a93b2ea14de112aa4e959db8f8af25
hash_basis: raw LF bytes
---

# THM-3869 -- an exact nonmonogenic cubic, with one line still unpaid

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
Work over an algebraically closed field `k` of characteristic zero.  Put

```text
B=k[x,y],                         K=Frac(B),

Delta=81x^5+90x^4+25x^3+30x^2y^2+30xy^2-y^4+8y^2,           (1)
a=x+1,                  ell=9x+4,
q=y^2-(15x^2+15x+4).                                      (2)
```

Then

```text
Delta=a^3 ell^2-q^2.                                       (3)
```

Consequently the polynomial representatives

```text
P=a ell^2,                         Q=q ell^2                 (4)
```

satisfy the exact square-residual Cardano identity

```text
P^3-Q^2=Delta ell^4.                                       (5)
```

The depressed cubic

```text
f(T)=T^3-3PT+2Q                                             (6)
```

is irreducible over `K` and has generic Galois group `S3`.  More strongly,
if `tau` is a root and

```text
z=tau^2/ell,                                                (7)
```

then

```text
O=B direct_sum B tau direct_sum B z                         (8)
```

is the integral closure of `B` in `K(tau)`.  It is finite free, normal, and
globally nonmonogenic.  The three discriminants which must not be conflated
are

```text
disc_T(f)=disc_B(1,tau,tau^2)=108 Delta ell^4,               (9)
disc_B(1,tau,z)=108 Delta ell^2,                             (10)
Disc_div(K(tau)/K)=V(Delta)+2V(ell).                         (11)
```

Thus one factor `ell^2` in `(9)` is the square of the power-order index,
but the other `ell^2` in `(10)-(11)` is genuine field ramification.  At the
generic point of `Delta` the cubic has one ramified and one unramified sheet,
of type `(2,1)`.  At `ell=0` it is tamely and totally ramified of degree
three.  This is an explicit normal, nonmonogenic, degree-three completion
model with two rational branch components, not a Keller map or a planar
Jacobian counterexample.

## 1. The identity and its normalization meaning

Equation `(3)` is direct expansion.  On the THM-3854 normalization

```text
x=t^4-2t^2,                         y=3t^5-5t^3,              (12)
```

one has

```text
a=(t^2-1)^2,
ell=9t^4-18t^2+4,
q=(t^2-1)^3(9t^4-18t^2+4).                                 (13)
```

Hence, with

```text
h=(t^2-1)(9t^4-18t^2+4),                                  (14)
```

the representatives in `(4)` pull back exactly to

```text
P(t)=h(t)^2,                         Q(t)=h(t)^3.             (15)
```

This is not a fourth seminormal direction omitted by THM-3864.  For its
first basis element

```text
h_1=t^2(t^2-1)(9t^2-14),                                   (16)
```

one has the literal equality

```text
h=h_1-4(x(t)+1).                                            (17)
```

Thus `[h]=[h_1]` in the three-dimensional defect `S/R` of THM-3864.  Their
common cusp-derivative coordinate is

```text
(h'(0),h'(1),h'(-1))=(0,-10,10).                            (18)
```

The gain in `(5)` comes from changing the polynomial representatives of an
existing direction, not from discovering a new direction.

There is also a complementary exact polarization.  Put

```text
F=15a^2-15a+4,                  M=9a-4=9x+5,
b=1-a=-x.                                                   (19)
```

Then

```text
F^2-a^3 ell^2=b^3M^2,                                      (20)
Delta=y^2(2F-y^2)-b^3M^2.                                  (21)
```

This explains why the neighboring factor `9x+5` recurs in the odd
seminormal control of THM-3854, but it is not a second polynomial Cardano
peel.  On `(12)`,

```text
b=t^2(2-t^2),                                               (22)
```

and `2-t^2` has two simple roots.  Thus `b` is not a square in `k[t]`.
Identity `(20)` trades the `ell` payment for `M` plus a nonsquare carrier;
it does not remove the payment.

## 2. Irreducibility and the `S3` field

Localize `B` at the height-one prime `(ell)`.  This is a DVR with residue
field `k(y)`.  Both

```text
a mod ell=5/9,                    q mod ell=y^2-8/27         (23)
```

are nonzero.  The nonzero coefficient points of `(6)` have valuations

```text
(T^0,T^1,T^3):                  (2,2,0).                     (24)
```

The lower Newton polygon is the single edge from `(0,2)` to `(3,0)`;
the `T` coefficient lies strictly above it.  Its horizontal and vertical
lengths are coprime.  The Newton-polygon factor theorem therefore forces
every nonconstant factor to have degree divisible by three.  The cubic is
irreducible over the completed local field, hence over `K`.

By `(5)`,

```text
disc_T(f)=108 Delta ell^4.                                  (25)
```

THM-3854 proves `Delta` irreducible.  Its valuation in `(25)` is odd, while
`ell^4` is a square.  Therefore the cubic discriminant is not a square in
`K`.  An irreducible cubic with nonsquare discriminant has Galois group
`S3`.

## 3. The corrected maximal order

Inside the cubic field define `z` by `(7)`.  The following relations are
immediate from `(6)`:

```text
tau^2=ell z,
tau z=3a ell tau-2q ell,
z^2=3a ell z-2q tau.                                       (26)
```

They prove that the free `B`-module `(8)` is closed under multiplication.
Linear independence follows from that of `1,tau,tau^2` over `K`, so `(8)`
is a finite-free cubic domain.  The trace matrix in its displayed basis is

```text
      [ 3          0             6a ell       ]
      [ 0          6a ell^2     -6q ell       ]
      [ 6a ell    -6q ell        18a^2 ell^2  ].             (27)
```

Its determinant is exactly `(10)`.

We now prove that this is the maximal order, rather than merely a better
order.  Away from `Delta ell=0`, determinant `(10)` is a unit and the order
is etale.  At the generic point of `Delta`, `ell` is a unit and `(8)` equals
the power order.  Its discriminant valuation is one.  In the DVR
discriminant-index formula

```text
v(disc O)=v(disc O_bar)+2 length(O_bar/O),                   (28)
```

the left side is one, so the nonnegative even index contribution is zero.

At the generic point of `ell`, `q` is a unit.  Eliminating `tau` from
`(26)` gives

```text
z(z-3a ell)^2-4q^2 ell=0,                                   (29)
```

or

```text
z^3-6a ell z^2+9a^2 ell^2z-4q^2 ell=0.                      (30)
```

This is Eisenstein at `ell`.  Hence the local ring generated by `z` is a
DVR, with `z` as a uniformizer.  The last relation in `(26)` recovers

```text
tau=(3a ell z-z^2)/(2q),                                    (31)
```

so it is exactly the localization of `(8)`.

Thus `(8)` is regular in codimension one.  It is Cohen--Macaulay, hence
`S2`, because it is finite free over the regular ring `B`.  Serre's
`R1+S2` criterion proves normality.  Therefore `(8)` is the integral
closure of `B` in the cubic field and `(10)` is the field-discriminant
ideal, up to the scalar unit `108`.

The basis change from `(1,tau,tau^2)` to `(1,tau,z)` has determinant
`1/ell`.  Equivalently, the power order has divisorial index `(ell)` in the
normal order.  This accounts for exactly two of the four copies of `ell`
in `(9)` and proves the distinction `(9)-(11)`.

## 4. Exact ramification and the cyclic layer

At `Delta`, the field-discriminant valuation is one, so tame cubic inertia
is a transposition: one sheet has `e=2` and the companion sheet has `e=1`.
At `ell`, Eisenstein equation `(30)` gives one sheet with `e=3`; its tame
different exponent is `3-1=2`, exactly the exponent in `(10)`.

The same fact is visible in the `C3` layer of the `S3` closure.  Let

```text
w^2=Delta,                         i^2=-1.                    (32)
```

Over `ell=0`, equation `(3)` becomes `w^2=-q^2`, so the line splits on the
quadratic surface into

```text
E_+: (ell,w-iq),                    E_-: (ell,w+iq).          (33)
```

A Cardano Kummer carrier is

```text
alpha=-Q+i ell^2w=ell^2(-q+iw),                              (34)
```

and

```text
Norm(alpha)=P^3.                                             (35)
```

On one component in `(33)`, the parenthesis in `(34)` is a unit; on the
other it has order two because its product with the conjugate is

```text
q^2+w^2=q^2+Delta=a^3 ell^2.                                (36)
```

Consequently the two valuations of `alpha` are

```text
(2,4),                         equivalently (2,1) mod 3.     (37)
```

Both are nonzero modulo three.  The cyclic cubic layer is therefore
ramified of degree three over both split primes above `ell`.  The line in
`(11)` is intrinsic field ramification, not an order-index illusion.

The two branch components nevertheless have the correct elementary
normalizations for the THM-3862 completion contract: `Delta` has affine
normalization `A1` by THM-3854, and `ell=0` is a line.  Their projective
closures share the unique point `[0:1:0]` at infinity.  The precise failure
is the presence of the second component, not its curve type.

There is a stronger atlas obstruction.  The fundamental equality at the
generic point of `ell` is exhausted by the one prime `E` of ramification
index three.  Hence

```text
div_O(ell)=3E.                                               (37a)
```

Let `U` be the maximal etale locus of `Spec O -> Spec B`.  It omits `E`, so
`ell` restricts to a nonconstant unit on `U`.  A dominant morphism
`A2 -> U` would pull that unit back to a scalar; its image would then lie in
one level curve of `ell`, contradicting dominance.  Thus this completion
cannot be the finite target completion of a polynomial-plane Keller atlas.

## 5. Global nonmonogenicity

For an arbitrary element modulo translation by `B`, write

```text
theta=u tau+v z,                         u,v in B.            (38)
```

The determinant of `(1,theta,theta^2)` relative to `(1,tau,z)` is the
binary index form

```text
I(u,v)=ell u^3-3a ell uv^2+2qv^3.                            (39)
```

There is no `u^2v` term: the scalar term `-4q ell uv` in `theta^2` does
not enter this determinant.  The coefficient ideal of `(39)` is exactly

```text
(ell,q),                                                     (40)
```

which is proper; its zero set consists of the two points

```text
x=-4/9,                         y^2=8/27.                    (41)
```

Every value of `(39)` therefore lies in a proper ideal and cannot be a
unit of `B`.  The index-form criterion for a finite-free cubic algebra
proves that no element generates `(8)` over `B`.  The normal completion is
globally nonmonogenic.

This is the positive part of the near miss: the algebra passes exactly the
nonmonogenicity gate that defeated the simpler two-place cubic of THM-3847.
It still fails the desired discriminant passport by `(11)`.

## 6. Trace shifts and all scalar representative shifts

A source-root translation `T |-> T+s`, with `s in B`, does not change the
polynomial discriminant, the cubic field, or its field-discriminant divisor.
It cannot remove `V(ell)` from `(11)`.

There is also a complete exact obstruction to the cheapest changes of the
square/cube representatives.  For constants `c,d in k`, put

```text
P_c=P+c Delta,                       Q_d=Q+d Delta,
R_(c,d)=(P_c^3-Q_d^2)/Delta.                              (42)
```

Then `R_(c,d)` is a square in `B` if and only if

```text
c=d=0,                         R_(0,0)=ell^4.                (43)
```

Thus every scalar change either destroys the square residual or retains
the same paid line.

Here is a coefficient proof.  On `ell=0`, put

```text
q_0=y^2-8/27.
```

Exact specialization gives

```text
R_(c,d)|_(ell=0)=q_0^2(d^2+c^3q_0^2).                       (44)
```

If `c!=0` and the right side is a square, comparison of the three
coefficients of a quadratic square in `y^2` forces `d=0`.  If `c=0` and
`d!=0`, write the full residual as a quadratic in `Y=y^2`; its discriminant
is

```text
4d^4a^3ell^2,                                               (45)
```

which is nonzero, so it is not a square.

It remains to exclude `c!=0,d=0`.  The residual is a quartic in `Y` with
nonzero `Y^0` row, so any square root is even in `y` and has the form

```text
r_2Y^2+r_1Y+r_0.
```

Choose `rho in k*` with `rho^2=c^3`.  The `Y^4,Y^3,Y^2` rows force
`r_2,r_1,r_0`; the `Y` row then agrees automatically, while the final
constant row misses by exactly

```text
ell^4(3ca^2+4)/4.                                           (46)
```

This is a nonzero polynomial for every scalar `c`.  Equations `(44)-(46)`
exhaust all scalar pairs and prove `(43)`.

Nonconstant changes

```text
P+A(x,y)Delta,                         Q+B(x,y)Delta          (47)
```

are not classified here.  They are the first honest coefficient lane which
can change the field rather than merely its root coordinate.

## 7. Scope and reproduction

THM-3869 constructs a genuine normal nonmonogenic `S3` cubic over the target
plane.  It does **not** construct an etale cubic completion with sole branch
`Delta`, an affine-plane source atlas, a Keller map, or a Jacobian
counterexample.  The exact remaining design problem is to replace `(4)` by
nonconstant representatives or a different binary cubic so that the
maximal-order field discriminant loses `ell` while retaining `Delta`.

Run

```bash
python3 04-computation/jc2_three_cusp_cardano_line_ramification_thm3869.py
python3 -O 04-computation/jc2_three_cusp_cardano_line_ramification_thm3869.py
```

and compare both streams byte-for-byte with
`05-knowledge/results/jc2_three_cusp_cardano_line_ramification_thm3869.out`.
The companion uses exact rational polynomial arithmetic and has no inactive
`assert` gates.
