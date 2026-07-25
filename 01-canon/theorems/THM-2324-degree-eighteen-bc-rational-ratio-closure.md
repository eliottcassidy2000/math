---
id: THM-2324
title: "Degree-eighteen B--C ratio-bank closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In the genuine
  nonsplit polynomial exact-square-prefix degree-eighteen branch of
  THM-2262/2297, all nine B--C ratio points in THM-2311's two-sparse bank
  are empty. The two rational points have absolutely irreducible
  genus-three trigonal spectra with ten simple branches and one ordinary
  nonvertical node. The two-point quadratic orbit has genus four: ten
  simple branches and one smooth totally ramified cubic fibre. The
  five-point quintic orbit has genus three: ten simple branches and one
  ordinary nonvertical node. Infinity is unramified throughout. Every
  rational Keller trajectory is therefore constant and gives the
  inherited nonsplit-deck contradiction. Together with THM-2314, THM-2316,
  and THM-2320, this leaves only the eight B--W ratios in the exactly
  two-sparse bank. This does not prove JC(2).
source: codex-2026-07-25-degree18-bc-bank-closure
depends_on:
  - THM-2262-degree-eighteen-trigonal-spectral-discriminant-reduction
  - THM-2297-degree-eighteen-target-translation-normal-form
  - THM-2311-degree-eighteen-two-sparse-weighted-ratio-bank
related:
  - THM-2314-degree-eighteen-bd-linear-ratio-closure
  - THM-2316-degree-eighteen-cd-linear-ratio-closure
  - THM-2320-degree-eighteen-dw-ratio-bank-closure
script: 04-computation/jc2_degree18_bc_ratio_bank_closure_thm2324.py
output: 05-knowledge/results/jc2_degree18_bc_ratio_bank_closure_thm2324.out
script_sha256: 52b1ebe8bf9ff1a5e48f113c3defa37ae9242b54b15008ae050e425121407697
output_sha256: 68c354ef57052cdf6f41ef34a9dc2924e92a144cbca46b05c61bb8e512ca40f5
hash_basis: working-tree bytes (LF)
---

# THM-2324 -- the full B--C ratio bank has positive genus

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2311 reduces the exactly two-sparse degree-eighteen branch to a finite
weighted-projective bank. On the `B`--`C` line, with

```text
t=C^2/B^3,
```

the four ratio factors are

```text
P_1(t)=2000+15309t,
P_2(t)=125+1134t,

P_3(t)=3321125-161754894t-2812385772t^2,

P_4(t)
 =410644531250000
  -18114791748046875t
  -545436228093750000t^2
  -4951165276923468750t^3
  -18946967714644599000t^4
  -26529827304546537363t^5.                         (1)
```

Their degrees `1+1+2+5=9` are the full off-axis `B`--`C` bank. This
theorem treats all four Galois orbits.

## 1. Weighted representatives

Use THM-2297's target-translation normal form

```text
(B,C,D,W),                         weights (2,3,4,5),

G_0(u,y;B,C,D,W)=0,               wt(u,y)=(2,1).    (2)
```

At a rational ratio `t`, put

```text
B=C=1/t.                                             (3)
```

Then `C^2/B^3=t`. For either algebraic orbit, let `alpha` be a root of
the relevant factor and put

```text
(B,C)=(alpha,alpha^2).                               (4)
```

Again `C^2/B^3=alpha`. These choices avoid denominators in the number-field
curves and represent the required weighted orbits.

## 2. Absolute irreducibility

Every specialized spectral curve is cubic in `u`, with the nonzero
constant leading coefficient `-26040609`. If such a curve were reducible
over an algebraic closure of its constant field, it would have a root in
that field's rational function field in `y`. After dividing by the leading
coefficient, the root is integral over the polynomial ring in `y`; that
ring is integrally closed, so the root is a polynomial.

If its degree were `d>2`, the cubic term would have `y`-degree `3d`,
strictly larger than the possible degrees

```text
2d+2,             d+4,             6                (5)
```

of the remaining terms. Hence every possible root is

```text
u=ay^2+by+c.                                        (6)
```

There is in fact a uniform coefficient proof that avoids specialization.
Let

```text
L_infinity(a)
 =1127-138915a+1607445a^2-26040609a^3.
```

The `y^6` coefficient after (6) is `L_infinity(a)`, and the `y^5`
coefficient is `b L'_infinity(a)`.  The infinity cubic is separable, so
these two equations force `b=0`.  With `b=0`, the `y^3` coefficient is
exactly

```text
-435456 C.
```

Every point in the off-axis `B`--`C` bank has `C!=0`, a contradiction.
Thus absolute irreducibility is already forced by the interaction between
the separable infinity fibre and the odd `C`-term.

As an independent exact control, for each rational point substituting (6)
and equating the seven coefficients of `1,y,...,y^6` gives the reduced
Groebner basis

```text
{1}.                                                 (7)
```

For an algebraic orbit, perform the same substitution in the universal
curve

```text
G_0(u,y;t,t^2,0,0)
```

and adjoin `P_3(t)` or `P_4(t)`. Exact Buchberger reduction over
`Q[a,b,c,t]` again gives

```text
{1}.                                                 (8)
```

Thus both the coefficient argument and the Gröbner computation show that
every curve below is absolutely irreducible, and its projective
normalization is a connected degree-three cover of the `y`-line.

## 3. The two rational points have genus three

The rational ratios and representatives are

```text
t_1=-2000/15309,              B=C=-15309/2000,
t_2=-125/1134,                B=C=-1134/125.         (9)
```

After nonzero constant rescaling and primitive integral normalization, the
two curves are

```text
H_1
 =3321506250u^3-205031250u^2y^2+48427561125u^2
  +17718750uy^4-2790065250uy^2+156905298045u
  -143750y^6+76271625y^4-425152800y^3
  -5811307335y^2+41841412812y,                      (10)

H_2
 =1660753125u^3-102515625u^2y^2+28697814000u^2
  +8859375uy^4-1653372000uy^2+110199605760u
  -71875y^6+45198000y^4-251942400y^3
  -4081466880y^2+29386561536y.                      (11)
```

For the unscaled curves, exact `u`-discriminant factorization gives

```text
Delta_1
 =-1628150074335205281/3906250000
   *(10y-81)^2(200y^2-19683)h_8(y),                 (12)

h_8
 =18400000000y^8+298080000000y^7+1653372000000y^6
  -11479125600000y^5-441659357460000y^4
  -3389154437772000y^3+2287679245496100y^2
  +111181211331110460y+450283905890997363,

Delta_2
 =-6668902704477000830976/244140625
   *(5y-54)^2(5y+54)h_9(y),                         (13)

h_9
 =44921875y^9+485156250y^8-455625000y^7
  -59049000000y^6-1142598150000y^5
  -2008846980000y^4+77484097800000y^3
  +234311911747200y^2-1446039226782720y
  -15617223649253376.
```

The degree-ten residuals in (12)--(13) are squarefree and coprime to their
displayed exceptional linear factors. Thus each curve has ten distinct
simple discriminant roots away from the exceptional value.

At those exceptional values the complete fibres are

```text
G_1(u,81/10)
 =-26040609/1000000
   *(100u+81)^2(100u+891),                           (14)

G_2(u,54/5)
 =-26040609/3125
   *(5u+36)(25u+36)^2.                              (15)
```

The repeated points are

```text
(u,y)=(-81/100,81/10),       (-36/25,54/5).
```

At both, `G=G_u=G_y=0`. If the tangent cone is

```text
A U^2+B UV+C V^2,
```

the exact pairs `(A,B^2-4AC)` are

```text
(-2109289329/10,
 -106778435362398485784/15625),

(-3749847696/25,
  1199902527419435384832/15625).                    (16)
```

Both points are ordinary nodes. Since `A!=0`, neither tangent is the
vertical fibre; both normalization branches are unramified over `y`.
The simple third point is unramified as well. Consequently only the ten
simple residual roots ramify, each with index two.

## 4. The quadratic orbit has genus four

Modulo nineteen, the monic reduction of `P_3` is

```text
t^2-7t+5,
```

which has no root in `F_19`. Hence `P_3` is irreducible over `Q`. Let

```text
K_2=Q(alpha),                         P_3(alpha)=0.
```

In `K_2[y]`,

```text
Delta_alpha(y)
 =-153384762202971019112448
   *(y-Y)^2h_(10,alpha)(y),                          (17)

Y=(69984alpha-4075)/8919,                           (18)
```

where `h_(10,alpha)` is monic of degree ten and

```text
gcd(h_(10,alpha),h'_(10,alpha))=1,
gcd(h_(10,alpha),y-Y)=1.                            (19)
```

Put

```text
R=5(237718152alpha+3321125)/2867360391.             (20)
```

The exceptional fibre is

```text
G_alpha(u,Y)=-26040609(u-R)^3,                      (21)
```

while

```text
partial_y G_alpha(R,Y)
 =128(16490073897927alpha-235328275250)
   /60214568211
 !=0.                                               (22)
```

The nonzero statement follows because the numerator is a nonzero
degree-at-most-one element of the quadratic field. Thus the curve is
smooth at `(R,Y)`, with one ramification point of index three. The ten
roots in (19) contribute ten more, for total ramification twelve and

```text
2g-2=3*(-2)+12=6,

g=4.                                                 (23)
```

Every identity is in `K_2`, so both complex conjugates have genus four.

## 5. The quintic orbit has genus three

The monic reduction of `P_4` modulo thirteen is

```text
t^5-3t^4+2t^3-3t^2-t-4.                            (24)
```

Rabin's exact finite-field test gives

```text
gcd(P_4,t^13-t)=1,
P_4 divides t^(13^5)-t                              (25)
```

after reduction modulo thirteen. Since five is prime, (25) proves that
`P_4` is irreducible. Let

```text
K_5=Q(alpha),                         P_4(alpha)=0.
```

Exact Euclidean reduction in `K_5[y]` gives

```text
Delta_alpha(y)
 =-153384762202971019112448
   *(y-Y)^2h_(10,alpha)(y),                         (26)
```

with

```text
deg h_(10,alpha)=10,
gcd(h_(10,alpha),h'_(10,alpha))=1,
gcd(h_(10,alpha),y-Y)=1.                            (27)
```

The exact common root of the exceptional fibre and its `u`-derivative is
`R in K_5`. Put

```text
S=(40/21)alpha+(5/81)Y^2-2R.                       (28)
```

Then

```text
G_alpha(u,Y)=-26040609(u-R)^2(u-S),                (29)

G_alpha(R,Y)
 =partial_u G_alpha(R,Y)
 =partial_y G_alpha(R,Y)=0.                         (30)
```

The exact companion reduces `Y`, `R`, and the tangent cone to the
five-element power basis of `K_5`. Its `U^2` coefficient and discriminant
are both nonzero degree-at-most-four elements. Therefore `(R,Y)` is an
ordinary node, and the nonzero `U^2` coefficient makes both tangent
branches nonvertical. They and the simple point `S` are unramified over
`y`.

Only the ten simple roots in (27) ramify. Thus

```text
2g-2=3*(-2)+10=4,

g=3.                                                 (31)
```

All five embeddings of `K_5` inherit the same exact factorization and
local type, so all five curves have genus three.

## 6. Infinity is unramified in every orbit

Use the weighted infinity chart

```text
z=1/y,                         v=u/y^2.              (32)
```

The `B` and `C` terms have lower weighted order, so at `z=0` every curve
has the same cubic

```text
L_infinity(v)
 =1127-138915v+1607445v^2-26040609v^3.              (33)
```

Its discriminant is

```text
Disc_v L_infinity
 =-153384762202971019112448 !=0.                    (34)
```

There are three distinct smooth points above infinity, with `z` a local
parameter at each. Infinity is unramified. Equations (12)--(34) therefore
list every ramification contribution used in the genus counts.

## 7. Every B--C Keller trajectory is impossible

A putative survivor in the retained branch supplies

```text
(u,y) in C(x)^2,                  G_0(u,y)=0.         (35)
```

If this pair is nonconstant, it extends to a nonconstant morphism from
`P^1` to the relevant projective normalization. Riemann--Hurwitz forbids
a nonconstant map from genus zero to any of the genus-three or genus-four
curves above. Thus `u,y` are constant.

Undo the constant weighted curve isomorphism used in (3) or (4). The wall
`y=0` is already empty by THM-2262. In the original normalization the
inherited first-flux identity

```text
Z=T^2=-2N_2/(5103y)                                 (36)
```

then makes `Z`, nonzero `T`, and `q` constant. The genuine nonsplit deck
fixes the algebraically closed constant field but sends `q` to `-q`,
contradicting `q!=0`. All nine ratios in (1) are empty.

## 8. Exact gain and stopping boundary

The complete `B`--`C` bank of THM-2311 is empty:

```text
9 -> 0.                                             (37)
```

Together with the audited `B`--`D`, `C`--`D`, and `D`--`W` closures, the
exactly two-sparse bank shrinks from `31` to

```text
BW:8,

total:8.                                            (38)
```

The connection ledger is

```text
source:
  all nine ratios in THM-2311's B--C bank;

map:
  use B=C=1/t at the rational points and
  (B,C)=(alpha,alpha^2) over the quadratic and quintic fields;

preserved:
  labelled weighted orbit, Galois orbit, spectral isomorphism class,
  normalization, branch indices, and genus;

temporarily forgotten:
  the scaled third flux, Keller one-form, and whole-Faber sidecar;

why restoration is unnecessary:
  genera three and four already make every rational spectral trajectory
  constant;

hostile control:
  B=C=1 has a squarefree degree-twelve branch discriminant, detecting an
  accidentally generic or identically repeated calculation.           (39)
```

This theorem closes only the `B`--`C` bank inside the genuine nonsplit
polynomial exact-square-prefix degree-eighteen branch. The eight `B`--`W`
ratios, every three-/four-sparse singular stratum, split/even descent,
other Newton edges, `JC(2)`, and `DC(2)` remain open.

## 9. Exact reproduction

Run

```bash
python3 04-computation/jc2_degree18_bc_ratio_bank_closure_thm2324.py
python3 -O 04-computation/jc2_degree18_bc_ratio_bank_closure_thm2324.py
```

Both executions are byte-identical to the stored output. The companion
verifies both rational representatives and primitive curves, the complete
rational factorizations (12)--(15), every residual squarefreeness and
coprimality statement, both rational tangent cones, irreducibility of the
quadratic and quintic ratio polynomials, the field factorizations
(17), (26), the smooth cubic fibre (20)--(22), the quintic node
(28)--(30), four unit Groebner bases, separable infinity, every
ramification count, and a squarefree hostile control. The
Riemann--Hurwitz and deck steps are the mathematical proof above, not
computer assumptions. QED.

The independent hostile audit also reconstructed the two rational
discriminants and tangent cones using pure fraction arithmetic, checked
both quadratic conjugates over `F_53` and all five quintic conjugates over
`F_353`, reran the ordinary companion byte-for-byte, and verified the
ratio-field irreducibility, infinity, ramification, deck, and scope
ledgers.  It found the coefficient proof in Section 2 independently of
the four Gröbner certificates.
