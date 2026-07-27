---
source: codex-2026-07-25-degree18-bc-algebraic-floor
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Within the genuine nonsplit polynomial exact-square-prefix degree-eighteen
  branch of THM-2262/2297, finite places 13 and 19 certify that every
  root of all four THM-2311 B--C ratio factors defines a
  geometrically irreducible trigonal curve with ten simple finite branch
  values and normalization genus at least three.  The genus/deck
  argument empties the full nine-point B--C bank.  The separate rational
  nodal packet refines two of these genus floors to exact nodal
  genus-three models.  The B--W/C--D/D--W banks, higher-support strata,
  split/even descent, JC(2), and DC(2) remain open.
depends_on:
  - THM-2262-degree-eighteen-trigonal-spectral-discriminant-reduction
  - THM-2297-degree-eighteen-target-translation-normal-form
  - THM-2311-degree-eighteen-two-sparse-weighted-ratio-bank
related:
  - jc2-degree18-bc-rational-nodal-genus3-closure-opus-20260725
  - jc2-degree18-bd-quartic-nodal-atlas-opus-20260725
  - THM-2345-degree-eighteen-common-root-wall-saturation
script: 04-computation/jc2_degree18_bc_algebraic_modular_floor.py
output: 05-knowledge/results/jc2_degree18_bc_algebraic_modular_floor.out
script_sha256: 262d99129916bb6904507322b824830e3c6ab140c6c996f76a0d0d31f12b2d4b
output_sha256: d9970040d8d818334ee43caf71f61e91c2e6ca137aadf742d70375470d5c26c9
hash_basis: LF-normalized working-tree bytes
---

# Finite places force genus at least three on the full B--C bank

## 1. The four parameter packets

THM-2311's `B`--`C` bank has the two rational linear ratios

```text
t=-2000/15309,                  t=-125/1134,        (1a)
```

and the roots of

```text
p_2(t)
 =3321125-161754894t-2812385772t^2,                 (1)

p_5(t)
 =410644531250000
  -18114791748046875t
  -545436228093750000t^2
  -4951165276923468750t^3
  -18946967714644599000t^4
  -26529827304546537363t^5.                         (2)
```

The rationalizing coordinate from the companion note applies uniformly.
Over `K=Q(t)` for any root of (1a), (1), or (2), normalize `B=1`,
choose `C^2=t` over the algebraic closure, put `x=Cy`, and use

```text
H_t(u,x)
 =-26040609t^3u^3
  +(49601160t^3+1607445t^2x^2)u^2
  +(-20995200t^3-2857680t^2x^2-138915tx^4)u
  -5598720t^3x+777600t^2x^2-435456t^2x^3
  +78120tx^4+1127x^6.                              (3)
```

The issue is to prove simultaneously for every conjugate parameter that
(3) is geometrically irreducible and has enough simple branch values.
A good finite place supplies both facts without factoring a degree-ten
polynomial over a number field.

## 2. The four finite-place certificates

The exact finite-field script finds:

```text
p_linear,1:
  prime 13,
  monic reduction 3+t,
  residue field F_13,
  irreducible fibre x=2;

p_linear,2:
  prime 13,
  monic reduction 7+t,
  residue field F_13,
  irreducible fibre x=1;

p_2:
  prime 19,
  monic reduction 5+12t+t^2,
  residue field F_(19^2),
  irreducible fibre x=3;

p_5:
  prime 13,
  monic reduction
    9+12t+10t^2+2t^3+10t^4+t^5,
  residue field F_(13^5),
  irreducible fibre x=3.                            (4)
```

The displayed reductions are irreducible: the script tests every monic
factor through half the parameter degree.  Thus (1) and (2) are
irreducible over `Q`, and each line in (4) is one good finite place of
its entire conjugacy packet.

At all four places:

```text
the leading u^3 coefficient is a unit;

deg_x Disc_u(H_t)=12;

deg_x gcd(Disc_u(H_t),partial_x Disc_u(H_t))=1;

the leading cubic at x=infinity is separable;

the displayed specialization has no residue-field root.               (5)
```

The last test is exact: for the residue-field size `Q`, the script
computes the gcd of the cubic with `u^Q-u`.  A cubic with no root is
irreducible.

## 3. Lifting the discriminant floor

Write

```text
Delta_t(x)=Disc_u H_t(u,x).                         (6)
```

THM-2311's exact repeated-branch resultant contains all four packets.
The changes `x=Cy` and the multiplication by `t^3` rescale the branch
polynomial and its base coordinate only by nonzero constants, so they
preserve the property of having a common factor with its derivative.
Therefore, over every characteristic-zero parameter field,

```text
deg gcd(Delta_t,Delta_t')>=1.                       (7)
```

At a good finite place, gcd degree can increase under reduction but
cannot decrease: divide by the unit leading coefficient, take the monic
common divisor over the local fraction field, and use integrality over
the DVR before reducing.  The degree-one residue gcd in (5) therefore
gives the reverse inequality.  Hence

```text
deg gcd(Delta_t,Delta_t')=1                         (8)
```

in characteristic zero for every root in all four packets.

Because `deg Delta_t=12`, equation (8) says exactly that there is one
double discriminant root and ten simple roots.  Each simple root is a
smooth index-two branch point of the normalization: the cubic's leading
coefficient is constant and nonzero, so there is no lost-degree
exception.  Whatever happens over the remaining double value can only
add nonnegative ramification.

## 4. Geometric irreducibility

The irreducible residue cubic at the displayed fibre proves irreducibility over
`K`: after dividing by the unit leading coefficient it is monic over
the local DVR, so any factorization over `K` would reduce without losing
`u`-degree.  The same constant-leading-coefficient specialization
argument then makes (3) irreducible in `K(x)[u]`.

As in the rational packet, this is automatically geometric.  The
`K`-rational point

```text
(u,x)=(0,0)
```

is smooth because

```text
H_x(0,0)=-5598720t^3!=0.                            (9)
```

If the `K`-irreducible curve split into conjugate geometric components,
the rational smooth point would lie on all of them and would be
singular.  Thus the normalization is one connected degree-three cover
of the `x`-line.

The infinity certificate in (5) also gives three distinct unramified
points there after putting `v=u/x^2` and `r=1/x`.  This is useful
boundary control, although the genus floor below only needs the ten
simple finite values.

## 5. Genus floor and Keller contradiction

Let `R` be the ramification degree of the normalized trigonal cover.
The ten simple finite values give

```text
R>=10.
```

Riemann--Hurwitz yields

```text
2g-2=-6+R>=4,

g>=3.                                               (10)
```

A rational Keller trajectory into this normalization is therefore
constant.  Constancy of `u,x` makes `y=x/C` constant over the algebraic
closure.  If `y!=0`, THM-2262's first flux makes `T^2,T,q` constant and
contradicts the genuine nonsplit deck.  If `y=0`, THM-2262 Section 4
closes the exceptional center.

Thus both linear roots, both roots of (1), and all five roots of (2) are
empty in the scoped branch.  This empties every one of THM-2311's nine
`B`--`C` ratios.  The separate rational nodal packet records sharper
local geometry at the two linear points but is not needed for bank
closure.

## 6. Information ledger

```text
source:
  all four B--C parameter packets;

map:
  one inert good finite place per packet, carrying both the branch gcd
  and the small-integral-fibre irreducibility test;

preserved:
  parameter conjugacy, discriminant degree, common-factor degree,
  cover degree, and the rational smooth point;

destroyed:
  exact characteristic-zero locations of the double and simple branch
  values;

sidecar:
  THM-2311 supplies the characteristic-zero lower bound gcd-degree>=1;
  good reduction supplies gcd-degree<=1;

hostile boundary:
  reduction alone cannot prove a repeated root exists in characteristic
  zero, and the resultant alone cannot cap the gcd degree; both halves
  are load-bearing;

next target:
  combine with the audited B--W/C--D/D--W packet and move to the first
  three-sparse face.                                   (11)
```

Run

```text
python3 04-computation/jc2_degree18_bc_algebraic_modular_floor.py
python3 -O 04-computation/jc2_degree18_bc_algebraic_modular_floor.py
```

Both executions must match

```text
05-knowledge/results/jc2_degree18_bc_algebraic_modular_floor.out
```

after LF normalization.
