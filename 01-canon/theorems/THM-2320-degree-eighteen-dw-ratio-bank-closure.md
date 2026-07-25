---
id: THM-2320
title: "Degree-eighteen D--W ratio-bank closure"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT. In the genuine
  nonsplit polynomial exact-square-prefix degree-eighteen branch of
  THM-2262/2297, all four D--W ratio points in THM-2311's two-sparse bank
  are empty. The rational point has an absolutely irreducible genus-four
  trigonal spectrum: ten simple branches plus one smooth totally ramified
  cubic fibre. The other three points form one irreducible cubic orbit.
  Uniformly over its ratio field, their spectra are absolutely irreducible
  of genus three: ten simple branches plus one ordinary node whose two
  normalization branches are unramified. Infinity is unramified
  throughout. Every rational Keller trajectory is therefore constant and
  gives the inherited nonsplit-deck contradiction. The D--W bank is empty;
  this does not prove JC(2).
source: codex-2026-07-25-degree18-dw-bank-closure
depends_on:
  - THM-2262-degree-eighteen-trigonal-spectral-discriminant-reduction
  - THM-2297-degree-eighteen-target-translation-normal-form
  - THM-2311-degree-eighteen-two-sparse-weighted-ratio-bank
related:
  - THM-2314-degree-eighteen-bd-linear-ratio-closure
  - THM-2316-degree-eighteen-cd-linear-ratio-closure
script: 04-computation/jc2_degree18_dw_ratio_bank_closure_thm2320.py
output: 05-knowledge/results/jc2_degree18_dw_ratio_bank_closure_thm2320.out
script_sha256: 42b35245b8965ffeb37f351ccbc39a7f5be808c9e3056e1c7b3a969a72acf48a
output_sha256: f20e8f98b4d4a3022a31662b1bef611b01d0d30b8734b3d6588ff13d21ab7760
hash_basis: working-tree bytes (LF)
---

# THM-2320 -- the full D--W ratio bank has positive genus

**PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT.**

THM-2311 reduces the exactly two-sparse degree-eighteen branch to a finite
weighted-projective bank. On the `D`--`W` line, its factors are

```text
935886848+430565625t,

36028797018963968+17932072576352256t
 -1448500838400000t^2+56162900390625t^3,

t=W^4/D^5.                                           (1)
```

The linear root has a rational weighted representative and an absolutely
irreducible genus-four spectral normalization. The three roots of the
cubic form one Galois orbit; a single computation over their cubic ratio
field proves that every conjugate normalization has genus three.

## 1. The rational weighted representative

Use THM-2297's target-translation normal form

```text
(B,C,D,W),                         weights (2,3,4,5),

G_0(u,y;B,C,D,W)=0,               wt(u,y)=(2,1).    (2)
```

Put

```text
t=-935886848/430565625,
r=1/t=-430565625/935886848,

(B,C,D,W)=(0,0,r,r).                                (3)
```

Then

```text
W^4/D^5=1/r=t.                                      (4)
```

THM-2311 proves that this ratio classifies the off-axis orbit on the
weighted `D`--`W` line. Substitution in `G_0`, followed by multiplication
by the nonzero scalar `7311616/49`, gives

```text
H(u,y)
 =-3885692518656u^3
  +239857562880u^2y^2
  -20728431360uy^4
  +3632067084375u
  +168167168y^6
  -134521003125y^2
  +403563009375y
 =0.                                                 (5)
```

This constant rescaling changes neither the affine curve nor its
normalization.

## 2. Absolute irreducibility at the rational point

The leading `u` coefficient in (5) is a nonzero constant. If `H` were
reducible over `C[y,u]`, its cubic `u`-degree would give a root in `C(y)`.
After division by the leading coefficient that root is integral over
`C[y]`; since `C[y]` is integrally closed, the root belongs to `C[y]`.

If its `y`-degree were `d>2`, substitution in (5) would leave the cubic
term with degree `3d`, strictly larger than the possible degrees

```text
2d+2,             d+4,             6                (6)
```

of the other terms. Hence every possible polynomial root has the form

```text
u=ay^2+by+c.                                        (7)
```

Substitution of (7) into (5), followed by equating the seven coefficients
of `1,y,...,y^6`, gives an exact coefficient ideal over `Q`. Its reduced
Groebner basis is

```text
{1}.                                                 (8)
```

The ideal remains the unit ideal after extending constants to `C`.
Therefore no root (7) exists, `H` is absolutely irreducible, and its
projective normalization is a connected degree-three cover of the
`y`-line.

## 3. Ten simple branches and one smooth cubic fibre

For the unscaled normalized curve `G_0` at (3), the exact
`u`-discriminant is

```text
Delta(y)
 =-1628150074335205281/97719251621562548224
   *(104y-405)^2 h_10(y),                            (9)
```

where

```text
h_10(y)
 =851140463828047757312y^10
  +6629074766353064263680y^9
  +38722720389995062886400y^8
  +201060278948051288064000y^7
  +609926274714380083200000y^6
  +999177077047752345600000y^5
  -1467500022145297296000000y^4
  -26582121110232749880000000y^3
  -59732448939525919050000000y^2
  -282296044395397652343750000y
  -549662971058346390380859375.                     (10)
```

Exact Euclidean reduction gives

```text
gcd(h_10,h_10')=1,
gcd(h_10,104y-405)=1.                               (11)
```

Thus the ten roots of (10) are distinct simple discriminant roots. Each
gives one normalization point of ramification index two, for total
ramification contribution ten.

At the remaining finite value

```text
y_0=405/104
```

the complete fibre is

```text
G_0(u,y_0)
 =-26040609/1265319018496*(10816u-3375)^3.          (12)
```

At

```text
u_0=3375/10816
```

direct differentiation gives

```text
partial_y G_0(u_0,y_0)
 =-692110561078125/95051008 !=0.                    (13)
```

The curve is smooth there. The first two `u`-derivatives vanish and the
third does not, so locally

```text
y-y_0=unit*(u-u_0)^3+higher terms.                  (14)
```

This is one point of ramification index three and contribution two. It
explains the squared factor in (9).

## 4. Infinity is unramified

Use the weighted chart

```text
z=1/y,                         v=u/y^2.              (15)
```

At `z=0`, independently of `D` and `W`,

```text
z^6G_0(v/z^2,1/z)
 =L_infinity(v)+O(z),

L_infinity(v)
 =1127-138915v+1607445v^2-26040609v^3.              (16)
```

Its discriminant is

```text
Disc_v L_infinity
 =-153384762202971019112448 !=0.                    (17)
```

The normalization therefore has three distinct smooth points over
infinity and `z` is a local parameter at each. Infinity contributes no
ramification. Equations (9)--(14) give total ramification

```text
10+2=12.
```

Riemann--Hurwitz for the connected degree-three cover yields

```text
2g-2=3*(-2)+12=6,

g=4.                                                 (18)
```

## 5. The cubic ratio orbit

The remaining ratios are the roots of

```text
P(t)
 =56162900390625t^3
  -1448500838400000t^2
  +17932072576352256t
  +36028797018963968.                               (19)
```

Its monic reduction modulo thirteen is

```text
t^3-4t^2+5t-1,
```

which has no root in `F_13`. Hence (19) is irreducible over `Q`. Its
discriminant is the nonzero integer

```text
-1239334538118889697076867375203545322038089020866560000000000000.
                                                               (20)
```

Let

```text
K=Q(alpha),                         P(alpha)=0.
```

The denominator-free representative

```text
(B,C,D,W)=(0,0,alpha^3,alpha^4)                    (21)
```

satisfies

```text
W^4/D^5=alpha.
```

The degree argument in Section 2 remains valid over an algebraic closure
of `K`. A possible factor would again give a root `ay^2+by+c`. To verify
all three conjugates at once, substitute this root into the universal
curve

```text
G_0(u,y;0,0,t^3,t^4)
```

and adjoin `P(t)` to its seven coefficient equations. Exact Buchberger
reduction over `Q[a,b,c,t]` gives

```text
{1}.                                                 (22)
```

Thus the curve is absolutely irreducible over every embedding of `K`.

## 6. The doubled factor is an ordinary node

In `K[y]`, the exact discriminant has the form

```text
Delta_alpha(y)
 =-153384762202971019112448
   *(y-Y)^2 h_(10,alpha)(y),                         (23)
```

where `h_(10,alpha)` is monic of degree ten and

```text
Y
 =(3790040625alpha^2-25880297472alpha-343597383680)
   /42471522304.                                    (24)
```

Exact Euclidean reduction gives

```text
gcd(h_(10,alpha),h'_(10,alpha))=1,
gcd(h_(10,alpha),y-Y)=1.                            (25)
```

Thus there are ten distinct simple branch values away from `Y`.

Put

```text
R
 =(-7058559375alpha^2+104983822336alpha+206158430208)
   /25517520000,

S=(5/81)Y^2-2R.                                     (26)
```

At the exceptional fibre,

```text
G_alpha(u,Y)=-26040609(u-R)^2(u-S),                 (27)
```

and

```text
G_alpha(R,Y)
 =partial_u G_alpha(R,Y)
 =partial_y G_alpha(R,Y)=0.                         (28)
```

For local coordinates `U=u-R`, `V=y-Y`, write the tangent cone as

```text
A U^2+B UV+C V^2.
```

The following exact power-basis elements of `K` are respectively `A` and
`B^2-4AC`:

```text
A
 =6561/202520000
  *(611312821875alpha^2
    -11216021225472alpha
    -19310172962816),

B^2-4AC
 =-945539748965690376192/57572832660675048828125
  *(136966554945917353125alpha^2
    -4257922067489564393472alpha
    -7854173444688698146816).                       (29)
```

Both are nonzero: their displayed numerators have degree at most two and
are nonzero, while `alpha` has degree three. Therefore `(R,Y)` is an
ordinary node. Since `A!=0`, neither tangent is the vertical fibre
`V=0`; both branches on the normalization are unramified over the
`y`-line. The third point `u=S` is also unramified.

Consequently only the ten roots in (25) ramify, each with index two.
Infinity is still governed by (16)--(17) and is unramified. Hence

```text
sum_P(e_P-1)=10,

2g-2=3*(-2)+10=4,

g=3.                                                 (30)
```

Every identity in this section holds in `K`, so all three complex
embeddings give the same genus-three conclusion.

## 7. Every D--W Keller trajectory is impossible

A putative survivor in the retained branch supplies

```text
(u,y) in C(x)^2,                  G_0(u,y)=0.         (31)
```

If this pair is nonconstant, it extends to a nonconstant morphism from
`P^1` to the relevant projective normalization. Riemann--Hurwitz forbids
a nonconstant map from genus zero to genus four or genus three. Thus
`u,y` are constant.

Undo the constant weighted curve isomorphism used to choose (3) or (21).
The wall `y=0` is already empty by THM-2262. In the original normalization
the inherited first-flux identity

```text
Z=T^2=-2N_2/(5103y)                                 (32)
```

then makes `Z`, nonzero `T`, and `q` constant. The genuine nonsplit deck
fixes the algebraically closed constant field but sends `q` to `-q`,
contradicting `q!=0`. All four ratios in (1) are therefore empty.

## 8. Exact gain and stopping boundary

The complete `D`--`W` bank of THM-2311 is empty:

```text
4 -> 0.                                             (33)
```

Together with THM-2314's full `B`--`D` closure and THM-2316's full
`C`--`D` closure, the exactly two-sparse bank shrinks from `31` to

```text
BC:9,                         BW:8,

total:17.                                           (34)
```

The connection ledger is

```text
source:
  all four ratios in THM-2311's D--W bank;

map:
  use D=W=1/t at the rational point and
  (D,W)=(alpha^3,alpha^4) over the cubic orbit;

preserved:
  labelled weighted orbit, spectral isomorphism class, normalization,
  branch indices, and genus;

temporarily forgotten:
  the scaled third flux, Keller one-form, and whole-Faber sidecar;

why restoration is unnecessary:
  genera four and three already make every rational spectral trajectory
  constant;

hostile control:
  D=W=1 has squarefree branch discriminant, detecting an accidentally
  generic or identically repeated calculation.                        (35)
```

This theorem closes only the `D`--`W` bank inside the genuine nonsplit
polynomial exact-square-prefix degree-eighteen branch. The `BC` and `BW`
banks, every three-/four-sparse singular stratum, split/even descent,
other Newton edges, `JC(2)`, and `DC(2)` remain open.

## 9. Exact reproduction

Run

```bash
python3 04-computation/jc2_degree18_dw_ratio_bank_closure_thm2320.py
python3 -O 04-computation/jc2_degree18_dw_ratio_bank_closure_thm2320.py
```

Both executions are byte-identical to the stored output. The companion
verifies the rational weighted representative, integral curve (5),
complete factorization (9)--(10), squarefreeness and coprimality (11),
smooth cubic fibre (12)--(13), cubic ratio-field irreducibility, field
factorization (23), node coordinates (24), (26), tangent cone (29), the
rational and global unit Groebner bases, separable infinity, ramification
arithmetic, and a squarefree hostile control. The Riemann--Hurwitz and deck
steps are the mathematical proof above, not computer assumptions. QED.
