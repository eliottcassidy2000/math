---
id: THM-2692
title: "Degree-twenty-two full-support terminal toric order-twelve closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In the inherited
  polynomial exact-square-prefix, genuine nonsplit degree-twenty-two branch,
  the remaining full-support BCDEW chart is empty.  The order-eleven Hensel
  equations cut the faithful sign-quotient torus to a reduced degree-seven
  carrier over each factor field.  Exact order twelve kills the carrier over
  the degree-five root field; a finite-free good-reduction/Nakayama
  certificate at p=103 kills it
  over the irreducible degree-ten unordered-pair field.  Hence the universal
  full-support quintic eliminant is absolutely irreducible, and the inherited
  fixed-place Kummer and y=0 arguments close the chart.  Together with the
  axis, plane, support-three, and support-four theorems, this empties every
  coefficient-support stratum in this inherited degree-twenty-two branch.
  Split/even short edges, integral raising, branches outside the inherited
  reduction, JC(2), and DC(2) remain open.
source: root/jc2655-referee/jc2692-logic-audit-2026-07-28
depends_on:
  - THM-2411-degree-twenty-two-first-flux-pole-divisor-square-class-reduction
  - THM-2636-degree-twenty-two-BCD-triple-spectral-square-Kummer-closure
  - THM-2671-degree-twenty-two-complete-support-three-weighted-hensel-closure
  - THM-2683-degree-twenty-two-complete-support-four-terminal-hensel-toric-closure
related:
  - MISTAKE-308
script: 04-computation/jc2_degree22_full_support_order12_toric_hensel.py
output: 05-knowledge/results/jc2_degree22_full_support_order12_toric_hensel.out
script_sha256: 1ef693fc76b162a8607904c9ca9b5d431fe7b61ccc556081784e9590f80c36bb
output_sha256: 0bbb9cada5bc0e99e8d154c2d1c06489adde69ab6f7ae3575139a84dbd339d2d
hash_basis: LF-normalized bytes
---

# THM-2692 -- the full-support degree-twenty-two chart is empty

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2683 leaves only the coefficient torus (subscripts distinguish the fixed
coefficients from normalized scale variables below)

```text
B_0 C_0 D_0 E_0 W_0 != 0.                                      (1)
```

inside THM-2411's inherited genuine nonsplit degree-twenty-two branch.  The
linear terminal-rank method cannot close this last chart: its order-eleven
matrix has rank four on a four-dimensional toric image, so the expected
intersection is finite rather than empty.  The decisive object is a finite
degree-seven carrier containing that scheme over each factor field.  The next
Hensel order removes every point of both carriers.

The pair-field calculation is performed by good reduction, but not by the
invalid inference that an arbitrary empty special fibre has empty generic
fibre.  Monicity makes the carrier finite free; every reconstruction
denominator is a unit; Nakayama's lemma then lifts the unit ideal from the
special fibre to characteristic zero.

## 1. Fixed section and exhaustive factor types

Use the normalized variables of THM-2636/2671.  Choose the constant
`rho in C^*` with `rho^2=B_0` and put

```text
c=C_0/rho^3, d=D_0/rho^4, e=E_0/rho^5, w=W_0/rho^6,
t=rho/y.                                                      (2)
```

Thus the normalized weighted coefficients `b_j` are

```text
(b_2,b_3,b_4,b_5,b_6)=(t^2,c t^3,d t^4,e t^5,w t^6),
c d e w != 0.                                                (2a)
```

The pinned sixty-term universal eliminant, divided by its constant leading
coefficient

```text
-88239118492602,
```

is a polynomial

```text
P(t,v) in Q[c,d,e,w,t,v],
deg_v P=5,       universal deg_t P=10,
                  specialized deg_t P<=10.                  (3)
```

monic in `v`.  Its fixed section is the same irreducible squarefree quintic
used in THM-2636, THM-2671, and THM-2683:

```text
L_5(v)=
 155624547606v^5+3215383215v^4-1700698560v^3
 +58124770v^2-855470v+2583.                                 (4)
```

Suppose a full-support specialization of `(3)` factors over `C(t)`.  Monicity
and Gauss integrality over the integrally closed ring `C[t]` give monic
factors in `C[t,v]`.  Their positive `v`-degrees survive at `t=0` and add to
five.  After exchanging the factors, the smaller fixed factor therefore
selects either one root or one unordered pair of roots of `(4)`.

The exact degree-five root field `K_1` and the separately verified
irreducible degree-ten unordered-pair field `K_2` exhaust these choices under
their complex embeddings.  Irreducibility of `(4)` alone would not imply
pair transitivity; the irreducibility of the degree-ten pair modulus is a
required input.

For the terminal calculation below, rename the four scale parameters

```text
(c,d,e,w) as (C,D,E,W).                                   (4a)
```

These uppercase letters are the constant scales in `(2)`, not the weighted
coefficients `c t^3,d t^4,e t^5,w t^6` themselves.

## 2. Orders eleven and twelve are necessary

Over either `K_i`, write the unique formal Hensel factors as

```text
Q=sum_(n>=0) q_n(v)t^n,       S=sum_(n>=0) s_n(v)t^n,
P=Q S.                                                         (5)
```

The fixed factors are coprime because `(4)` is squarefree.  The companion
implements the exact recursion

```text
f_n=r_n-sum_(i=1)^(n-1)q_i s_(n-i),
q_n=rem_(q_0)(f_n (s_0 mod q_0)^(-1)),
s_n=(f_n-q_n s_0)/q_0,                                      (6)
```

checks every division and reconstructed product coefficient, and lifts
through order twelve.

For an actual polynomial factorization, degree additivity and `(3)` imply

```text
q_11=s_11=q_12=s_12=0.                                     (7)
```

This remains necessary if a parameter specialization lowers the `t`-degree.
The complete order-eleven support is

```text
(E,EW,DE,C,CW,CD,CD^2,C^2E,C^3),                           (8)
```

and the five-by-nine coefficient matrix has rank four in both `K_1` and
`K_2`.  The order-twelve system has five rows and fourteen monomials; its
matrix also has rank four.  Rank alone is not the conclusion at either
order.

## 3. The faithful sign quotient has degree seven

On `(1)`, divide `(8)` by `C` and put

```text
a=E/C,                         h=C^2.                       (9)
```

The terminal monomial vector becomes

```text
m=(a,aW,aD,1,W,D,D^2,ah,h).                               (10)
```

This change is lossless for the equations.  Intrinsically,

```text
C[C^+-1,D^+-1,E^+-1,W^+-1]^((C,E)->(-C,-E))
  =C[a^+-1,D^+-1,W^+-1,h^+-1].                            (11)
```

Thus `(9)` is the quotient by the free simultaneous-sign action, not an
ad hoc reparameterization.  Every order-twelve monomial has even total
`C,E` exponent and descends through `(11)`.  Conversely, every torus point
of the quotient lifts over `C`: choose a square root `C^2=h` and set `E=aC`.

Eliminating `(a,D,W,h)` from the graph of `(10)` gives eight quadratic
binomials.  A regular lower triangulation of the faithful exponent
configuration has maximal simplices

```text
01257, 01458, 01578, 03458, 12567, 14568, 15678.            (12)
```

Every determinant in `(12)` is one.  Hence the projective toric image has
degree seven.  In the pre-quotient scale coordinates `(C,D,E,W)`, the same
determinants are two and the affine exponent-difference lattice has index two.  Those
ambient determinants must be divided by the lattice index; reporting degree
fourteen would count the discarded sign twice.

## 4. The exact order-eleven carrier

Choose four independent terminal rows.  Each row is linear in `W,h`, with
coefficients polynomial in `(a,D)`.  Any pair of rows gives a two-by-two pivot
`Delta_ij(a)`.  The six pivots satisfy

```text
gcd(Delta_01,Delta_02,Delta_03,Delta_12,Delta_13,Delta_23)=1. (13)
```

Thus their charts cover every terminal point.  On a fixed pivot chart, solve
the chosen rows for `W,h`.  The other two rows give compatibility polynomials

```text
F(a,D)=G(a,D)=0,              deg_D F=deg_D G=2.            (14)
```

In all six charts and over both factor fields, exact elimination gives, up to
a nonzero scalar in the factor field,

```text
Res_D(F,G) .= Delta_ij(a)^2 H_7(a),
deg Res_D=11,              deg Delta_ij=2,
deg H_7=7.                                                     (15)
```

Equivalently, `(15)` is a literal equality after monic normalization of the
resultant, squared pivot, and `H_7`.

The removed degree-four factor is exactly the squared pivot-collapse
artifact.  Within each factor field, the monic saturated carrier `H_7` is
identical in all six charts, is squarefree, and is coprime to every
corresponding pivot.  The root-field and pair-field polynomials need not be
the same.  Their common degree agrees with the intrinsic toric degree in
Section 3.

Over `K_1`, cancelling the `D^2` terms in `(14)` in every chart gives a
linear subresultant

```text
U(a)D+V(a).                                                  (16)
```

Here `U` is a unit modulo `H_7`, so `(16)` reconstructs `D`; the pivot then
reconstructs `W,h`.  Substitution satisfies all four selected terminal rows,
and the elements

```text
a, D, W, h                                                   (17)
```

are units in the root-field carrier algebra.

For `K_2`, the exact pre-reduction bridge uses less.  A terminal point lies
in some pivot chart by `(13)`; `(15)` then forces `H_7(a)=0`.  Since
`gcd(H_7,Delta_01)=1`, that point lies in pivot chart `01`.  Section 6
certifies the chart-`01` subresultant, reconstruction, and torus units by
finite-free good reduction.  Thus every possible order-eleven trajectory,
including every pivot boundary, is contained in one of the two finite
degree-seven torus carriers.

## 5. Exact root-field closure at order twelve

Over `K_1`, the companion performs every quotient operation exactly in

```text
A_1=K_1[a]/(H_7).                                           (18)
```

It reconstructs `(D,W,h)` independently on all six pivot charts and reduces
the five order-twelve rows from `(7)`.  Every residue has degree six, and in
every chart the successive gcd degrees with `H_7` are

```text
0,0,0,0,0.                                                  (19)
```

In fact the first order-twelve row already generates the unit ideal with
`H_7`.  Hence the root-field terminal scheme is empty.

## 6. Pair-field closure by finite-free good reduction

The order-eleven carrier over `K_2` is first constructed exactly in
characteristic zero.  All six pivot-dependent saturated carriers agree; all
pivots and compatibility polynomials are constructed exactly, and `(13)`
holds.  The chart-`01` reconstruction-unit and final quotient-algebra decisions
are delegated to good reduction.

Write `K_2=Q(alpha)`.  At the prime `p=103`, the irreducible degree-ten
modulus has the simple residue root

```text
alpha -> 61 in F_103.                                       (20)
```

Every coefficient denominator occurring in the exact carrier, the chart-`01`
data `F,G,Delta_01,U_01`, the terminal rows, and the reconstruction is a
`103`-unit.  Let

```text
O=(Z_(103)[alpha]/(P_10))_(103,alpha-61),
A=O[a]/(H_7).                                               (21)
```

The simple-root condition makes `O` unramified at this degree-one prime.
Because `H_7` is exact and monic, `A` is finite free of rank seven over `O`.
Reduction preserves degree seven and squarefreeness.  The reductions of

```text
a,D,W,h,Delta_01,U_01                                       (22)
```

are units.  Since `103A` lies in the Jacobson radical of the finite
`O`-algebra `A`, the same elements are units in `A`; no generic point can
disappear through a `103`-adic reconstruction pole.  Exact coprimality of
`H_7` and `Delta_01` means this one pivot chart contains the entire carrier;
the other five charts provide exact hostile agreement checks.

Let `I` be the ideal in `A` generated by the five reconstructed order-twelve
residues.  The exact finite-field calculation gives five degree-six residues
and aggregate gcd degrees

```text
0,0,0,0,0.                                                  (23)
```

Therefore

```text
I+103A=A,
(A/I)/103(A/I)=0.                                           (24)
```

The `O`-module `A/I` is finite, so `(24)` and Nakayama's lemma give `A/I=0`.
Tensoring with `K_2` proves the unit-ideal identity in the abstract pair field
and hence that the characteristic-zero terminal scheme is empty under every
complex embedding.  This is why one degree-one prime suffices.

Monicity and the unit gates in `(21)--(22)` are load-bearing.  Without them,
an empty special fibre need not imply an empty generic fibre; `103x-1=0` is
the elementary hostile pattern that the finite-free argument excludes.

## 7. Absolute irreducibility and the Kummer closure

Any proper factor of `(3)` selects a root or unordered root pair at `t=0`.
Its Hensel lift would therefore give an order-eleven/order-twelve terminal
point over `K_1` or `K_2`.  Sections 5--6 exclude both possibilities.  Hence

```text
P_(c,d,e,w)(t,v) is absolutely irreducible
for every c d e w != 0.                                    (25)
```

Now apply THM-2671's fixed-place mechanism.  Let `C_(c,d,e,w)` be the smooth
projective normalization of `P_(c,d,e,w)(t,v)=0`, proved irreducible in
`(25)`.  At its five points above
`t=0`, squarefreeness of `(4)` makes `t` a uniformizer and the first flux
makes `zeta` a unit.  The retained Kummer function is the element

```text
mathscr H=rho^3 zeta/t^3 in C(C_(c,d,e,w))^*.               (26)
```

It has valuation `-3` at all five points.  Thus it is nonsquare, its
connected quadratic cover ramifies at those five places and, by parity, at
least six places in total.  Riemann--Hurwitz gives genus at least two.

If `t` were nonconstant along a Keller trajectory, absolute irreducibility
would give an embedding

```text
phi:C(C_(c,d,e,w)) -> C(x).
```

Only after this pullback do the physical equations say
`phi(mathscr H)=Z=T^2`.  The physical element `T` therefore lifts `phi` to the
connected cover `X^2=mathscr H`, producing a nonconstant map from `P^1` to a
curve of genus at least two, which is impossible.  If `t` were constant,
then constant nonzero `rho` makes `y` constant; `v`, hence `u=v y^2`, then
`zeta,Z,T,q` are successively constant, contradicting the genuine nonsplit
deck.

The boundary `y=0` is handled before choosing `t`.  Full support has `B_0!=0`,
so THM-2671's uniform `B!=0` boundary quintic applies; its `u^5` coefficient
is the nonzero constant `-2264031`.  It makes `u`, then `Z,T,q`, constant and
again contradicts the deck.  Therefore the full-support chart `(1)` is empty.

## 8. Consequence, reproduction, and boundary

The degree-twenty-two coefficient strata have now all been exhausted:

```text
support 0--2: THM-2411 and the axis/plane closure chain;
support 3:    THM-2671;
support 4:    THM-2683;
support 5:    this theorem.                                (27)
```

Thus no trajectory remains in the inherited polynomial exact-square-prefix,
genuine nonsplit degree-twenty-two branch.

Run

```bash
python3 04-computation/jc2_degree22_full_support_order12_toric_hensel.py
python3 -O 04-computation/jc2_degree22_full_support_order12_toric_hensel.py
```

Both modes byte-match the declared output.  The companion pins the audited
THM-2636 engine by SHA-256, reconstructs the universal eliminant, both factor
fields, every Hensel coefficient through order twelve, the eight toric
quadrics and seven-simplex triangulation, all six pivot charts, both saturated
carriers, every unit, and the exact `p=103` residue certificate.  Independent
audits reconstructed the factor taxonomy, sign quotient, normalized lattice
degree, finite-carrier/Nakayama implication, Kummer step, and `y=0` boundary.

This theorem does not address split/even short-edge descent, integral
`2`-adic order raising, a branch outside the inherited reduction, another
terminal degree, `JC(2)`, or `DC(2)`.

QED.
