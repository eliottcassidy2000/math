---
id: THM-2747
title: "Highest-odd reduced boundary divisor and one-ended factorization atlas"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Every complete-
  intersection scheme in the full split degree-22 Faber family with a
  nonzero odd seed is reduced.  Its normalized infinity divisor consists of
  five fixed transverse residual points plus 3*gcd(22-j,6) G2 points, and the
  third response has a pole at every one.  Any hypothetical physical image
  component is rational, one-ended, and of weighted degree 1, 2, or 6; its
  response composition is pure-power, and the two lowest odd rows cannot use
  a G2 end.  THM-2745 separately closes all such physical components;
  THM-2755 subsequently closes the complementary all-even zero-flux edge.
  Other reduced degrees, the upstream reduction, JC(2), and DC(2) remain
  open.
source: odd-faber-component-boundary-scout-2026-07-28
audit: thm2694-full-lift-fibre-scout-2026-07-28 (independent scheme/divisor/one-ended audit, exact normal/-O replay, and hash verification: ACCEPT)
depends_on:
  - THM-2719-full-split-odd-faber-generic-normalization-genus-four-hundred-nineteen
  - THM-2723-split-exact-square-prefix-rational-primitive-pole-capacity
  - THM-2741-highest-odd-faber-response-pole-capacity-closure
related:
  - THM-2745-highest-odd-faber-componentwise-exact-prefix-closure
  - THM-2748-one-response-pole-monomial-composition-boundary
  - THM-2755-all-even-zero-flux-componentwise-global-regularity-closure
script: 04-computation/jc2_degree22_highest_odd_boundary_divisor_thm2747.py
output: 05-knowledge/results/jc2_degree22_highest_odd_boundary_divisor_thm2747.out
script_sha256: 482d14f3bac9ac46586c2b950fc777f458d0a2ea7d5e9b386a26217b4af21aa3
output_sha256: 0bb604be2f6b9c57756ee0b58171f000fbfc5911d1a73951dec23964278d3411
hash_basis: LF-normalized bytes
---

# THM-2747 -- highest-odd reduced boundary divisor and one-ended factorization atlas

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2745 excludes every physical highest-odd response component.  The present
theorem records the stronger scheme-theoretic atlas behind that closure.  The
apparent nonreduced escape is not real: the fixed top intersection and the
odd Newton germs make every component generically reduced.  Even before the
exact-prefix contradiction is used, a hypothetical physical trajectory can
land only on a rational component with one response pole, hence on one
weighted component of degree `1`, `2`, or `6`.

## 1. Statement

Work over `C` and use the full chosen-sheet split family

```text
F_23=Phi_22+sum_k a_k h^(22-k)Phi_k-lambda h^23,
G_24=Psi_22+sum_k a_k h^(22-k)Psi_k-W h^24             (1)
```

in `P(1,2,3,4)_[h,d,q,s]`.  Suppose at least one odd coefficient is
nonzero, let `j` be the highest such index, and put

```text
r=22-j,                    g=gcd(r,6).                 (2)
```

Let `C` be the complete-intersection scheme `(F_23,G_24)`.  Then:

1. `C` is a reduced pure curve.  It may be reducible.
2. On the normalization of `C`, the boundary `h=0` has five fixed simple
   residual points and `3g` points over `P_infty=[0:1:0:0]`.
3. Writing these divisors as `D_res` and `D_G2`, respectively,

```text
div(h)=D_res+(6/g)D_G2,

div_infty(R_aff)=25D_res+((150-6r)/g)D_G2,             (3)

deg div(h)=23,
deg div_infty(R_aff)=575-18r.                          (4)
```

Here every point occurs once in `D_res` or `D_G2`; the displayed
coefficients are its local orders.  In particular the third response has a
pole at every normalization point above infinity.

The degrees in `(3)`--`(5)` are weighted-stack `O(1)` degrees; equivalently,
they are one twelfth of the ordinary Cartier `O(12)` degree.  On the
normalization the central `mu_2` swaps two distinct index-cover branches, so
one coarse G2 branch retains order `6/g`; that order is not divided by two.

If one irreducible component `X` carries a physical polynomial Keller
trajectory, then `X` is rational and contains exactly one of those boundary
points.  Its weighted `O(1)` degree is therefore

```text
1                           if it contains a residual point,
6/g in {6,2}                if it contains a G2 point. (5)
```

For `r=19` and `r=21`, the second alternative in `(5)` is impossible.  Thus
the two lowest odd rows can escape THM-2741 only through a degree-one
residual-boundary component.

## 2. The five residual infinity points

At `h=0`, the two equations are the fixed top forms `Phi_22=Psi_22=0` in
`P(2,3,4)_[d,q,s]`.  If `q=0`, the second equation is `224s^6=0`, so
`s=0` and the only weighted-projective point is the common point `P_infty`.
Thus every other point lies in `q!=0`.

On that chart use the weight-zero invariants

```text
t=s^3/q^4,                         u=ds/q^2.           (6)
```

After setting `q=1` and clearing the nonzero power of `s`, the exact two
equations are

```text
fbar=-84u^2+(3+280t)u-21t-84t^2,

gbar=-224u^3+3360u^2t-336ut-3360ut^2
     +3t+560t^2+224t^3.                                (7)
```

Their `u`-resultant is `224t f(t)`, where

```text
f(t)=20141047808t^5-14386462720t^4+1089822720t^3
     -21288960t^2-35910t+81.                           (8)
```

The factor `t` is only the clearing artifact: `f(0)=81`, and the original
top equations have no point there.  Reduction modulo `47` makes `(8)`
irreducible.  Its exact discriminant is

```text
-2^86 3^16 5^12 7^15 13^2 19^2 29,                   (9)
```

so it has five distinct roots.  The Groebner basis of `(7)` makes `u`
unique over each root.  Consequently `(8)` parametrizes exactly five simple
coarse weighted-projective points.  The residual `mu_3` action on the `q=1`
cover is free because `s!=0`; its fifteen cover representatives form five
three-element orbits, not fifteen coarse points.

The invariant response numerator is, up to the nonzero scalar `33/1024`,

```text
s^2 R_22=
84t^2u-70t^2-280tu^2+105tu-3t+84u^3-9u^2.            (10)
```

The exact Groebner basis of `(fbar,gbar,(10))` is `(u,t)`.  Since every
genuine residual point has `t!=0`, `R_22` vanishes at none of them.  The
top intersection is therefore transverse there, `h` has order one, and
`R_aff=R_25/h^25` has pole order `25`.

## 3. The length eighteen G2 point

On the `d=1` index cover at `P_infty`, the two degree-six tangent forms are
the interlaced `G2=I_2(6)` quadratures

```text
Phi_(22,6)=-(231/128)q s(q^2-3s^2)(3q^2-s^2),

Psi_(22,6)=-(231/256)(q^2-s^2)
                         (q^4-14q^2s^2+s^4).          (11)
```

Their stripped resultant is `-2^24 s^36`, so the index-cover intersection
length is `36`.  The central `mu_2` quotient gives weighted length `18`.
Together with the five residual points this is exactly weighted Bezout:

```text
23*24/(2*3*4)=18+5=23.                                (12)
```

For the highest odd gap `r`, the complete weighted initial system is

```text
Psi_(22,6)=0,
Phi_(22,6)+c h^r=0,                  c!=0.             (13)
```

Along each of the six distinct `Psi` lines, `(13)` becomes

```text
h^r=c_alpha s^6,                    c_alpha!=0.        (14)
```

It has `g` distinct reduced branches, and the central involution freely
pairs the six line directions.  Hence the coarse normalization has `3g`
points.  On every one,

```text
ord(h)=6/g,
ord(R_25)=6r/g,
ord_infty(R_aff)=(150-6r)/g.                           (15)
```

Equations `(10)` and `(15)` prove the divisor formulas `(3)`--`(4)`.

## 4. Why no nonreduced member survives

The top forms are coprime, so `(F_23,G_24)` is a regular sequence for every
specialization.  Thus `C` is Cohen--Macaulay and pure of dimension one; in
particular it has no embedded zero-dimensional components.

The five residual boundary germs are smooth by Section 2.  At `P_infty`,
each binomial branch in `(14)` occurs with exponent one.  At its punctured
generic point, the derivative of `Psi_(22,6)` transverse to its simple line
and the derivative of `h^r-c_alpha s^6` along the binomial give rank two.
These are branch multiplicities, not the larger intersection order
`ord(h)=6/g`: every local minimal branch is generically reduced.

No one-dimensional component is contained in `h=0`, because the top
intersection there is zero-dimensional.  Every projective component meets
`h=0`.  Hence every minimal component of `C` is generically reduced.
A one-dimensional Cohen--Macaulay ring is `S_1`; `S_1` plus generic
reducedness implies reducedness.  This proves statement 1 and removes the
nonreduced boundary left open by THM-2741.

## 5. A physical component is rational and one-ended

Let `X` be the normalization of an irreducible component carrying a physical
trajectory.  THM-2723 gives

```text
U R_source'=kappa!=0,                                  (16)
```

so the source map is nonconstant.  It extends to a finite surjective map
`gamma:P1->X`.  Therefore `X` is rational.  More importantly, THM-2723's
rational-primitive classification says that `R_source` has at most one pole
point on `P1`.

Every boundary point on `X` is a response pole by `(3)`.  If `X` contained
two, surjectivity would give two disjoint nonempty pole fibres upstairs,
contradicting the source capacity.  Every projective component has at least
one boundary point, so `X` has exactly one.  Its weighted degree is the
coefficient of that point in `div(h)`, proving `(5)`.

This also gives the exact component-count boundary.  The whole reduced curve
has at most `5+3g` irreducible components.  A formal singleton partition of
the boundary has degree multiset

```text
1,1,1,1,1, (6/g repeated 3g times),                   (17)
```

whose sum is `23`.  Thus pole counting and weighted degree alone do not rule
out the remaining one-ended factorization.  Formula `(17)` is a sharp
combinatorial hostile, not a claim that such a Faber member exists.

## 6. Monomial ramification and the last two odd rows

Let the unique response pole of `X` have order `P`, so

```text
P=25                              on a residual component,
P=(150-6r)/g                     on a G2 component.    (18)
```

Choose a coordinate `T` on `X=P1` with the boundary point at infinity.  Then
`R_X(T)` is a polynomial of degree `P`.  If `U` were constant, THM-2723
would make `R_source` affine-linear, impossible because
`deg(R_source)=P deg(gamma)>=8`.

Thus, after translating the source and putting `y=1/(x-a)`, THM-2723 gives

```text
U=u_0(x-a)^m,
R_source=s_0+s_1 y^(m-1).                              (19)
```

The boundary has one totally ramified source fibre, so `T=T(y)` is a
polynomial of degree `e=deg(gamma)`.  Since

```text
R_X(T(y))-s_0=s_1 y^(m-1),                             (20)
```

every root fibre on the left must be the single point `y=0`.  Surjectivity
of a complex polynomial then forces

```text
R_X(T)-s_0=c(T-T_0)^P,
T(y)-T_0=c' y^e,
m-1=eP.                                                (21)
```

This is the exact monomial-composition boundary; arbitrary response
polynomials do not survive.

Finally, on a G2 branch

```text
ord(q_aff)=(r-18)/g.                                  (22)
```

For `r=19` and `r=21`, this is `+1`.  If such a branch were the only
boundary point of `X`, then `q_aff` would be regular on the projective curve
and vanish at that point, hence would be identically zero.  But its local
tangent law is `q=alpha s+...` with `alpha!=0`, a contradiction.  This
proves the final assertion of Section 1.

## 7. Scope and relation to the physical closure

The divisor partition `(17)` and pure-power composition `(21)` are a sharp
abstract stopping boundary after the polynomial exact-prefix lift equations
are forgotten: divisor and pole counting alone permit degree-`1`, degree-`2`,
and degree-`6` singleton components.  THM-2745 retains those lift equations
and excludes every such physical component.  Thus this theorem is a reusable
scheme/divisor atlas, not a competing physical closure route.

It does not justify the split degree-22 gauge for an arbitrary Keller pair.
THM-2755 subsequently closes the all-even zero-first-flux edge and therefore
finishes this chosen-sheet reduced-degree-22 family, but neither the upstream
quartic/split/degree reduction, other reduced degrees, `JC(2)`, nor `DC(2)`
follows.

## 8. Exact companion

Run

```text
python 04-computation/jc2_degree22_highest_odd_boundary_divisor_thm2747.py
python -O 04-computation/jc2_degree22_highest_odd_boundary_divisor_thm2747.py
```

and compare both transcripts with

```text
05-knowledge/results/jc2_degree22_highest_odd_boundary_divisor_thm2747.out.
```

The companion reconstructs the Faber rows, the invariant residual quotient,
the irreducible squarefree quintic, the response nonvanishing ideal, the
weighted length, every branch/pole invoice, and the sharp singleton-degree
hostile.  It uses explicit exceptions rather than Python `assert`.
