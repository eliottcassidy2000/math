---
id: THM-2470
title: "Degree-twenty-two E-W plane coprime-weight ramification closure"
status: >
  PROVED + VERIFIED-EXACT. In the open first-flux chart of the genuine
  nonsplit polynomial exact-square-prefix degree-twenty-two branch, the
  complete coefficient plane B=C=D=0 is empty. For E,W nonzero, the
  root-free invariants lambda=W^5/E^6 and r=W/(Ey) give a 13-term
  plane curve of bidegree (10,5). Its Newton triangle reduces every
  factorization to one linear or quadratic v-factor, and both exact
  coefficient ideals are unit, so every physical fibre is absolutely
  irreducible. The v-discriminant is
  constant*lambda^8*W_5^2*K_30, where W_5 is exactly the excluded
  first-flux wall. A complete exceptional-parameter factorization,
  together with a Sylvester-valuation bound and ramification parity,
  leaves total ramification at least twenty in every fibre.
  Riemann--Hurwitz gives normalization genus at least six. No rational
  trajectory survives. Together with the axes, this closes the sixth
  of ten support-two planes. It does not close the other four planes,
  higher mixed strata, JC(2), or DC(2).
source: codex-2026-07-27-degree-twenty-two-EW-ramification
depends_on:
  - THM-2411-degree-twenty-two-first-flux-pole-divisor-square-class-reduction
  - THM-2423-degree-twenty-two-W-axis-genus-two-and-origin-cusp-closure
  - THM-2425-degree-twenty-two-CDE-axis-hyperelliptic-closure
related:
  - THM-2429-degree-twenty-two-CW-plane-hyperelliptic-family-closure
  - THM-2437-degree-twenty-two-DW-plane-quartic-ramification-closure
  - THM-2463-degree-twenty-two-BD-plane-square-lift-closure
  - THM-2468-degree-twenty-two-BW-plane-square-lift-closure
  - THM-2469-degree-twenty-two-CD-plane-coprime-weight-ramification-closure
script: 04-computation/jc2_degree22_ew_plane_ramification_thm2470.py
output: 05-knowledge/results/jc2_degree22_ew_plane_ramification_thm2470.out
script_sha256: 87e71280db14fefb0057e2c2cc88920bb091a72d9d662689a6e480a1dd30299e
output_sha256: f431869559a28acedcd5a513f7917a98a1d0f704bc44e9ff89939c19559f9899
hash_basis: working-tree bytes (LF)
---

# THM-2470 -- the degree-twenty-two E-W plane is empty

**PROVED + VERIFIED-EXACT.**

The `C,D` closure of THM-2469 showed that a coprime pair of coefficient
weights can be attacked without a retained power cover: use a weight-zero
ratio and a root-free weight-one coordinate, then read genus directly from
the base curve.  The same operation is especially sparse for weights five
and six.  The exact conclusion is

```text
genuine degree-22 trajectory,
mathcal A!=0,
B=C=D=0
    => contradiction.                                           (1)
```

Thus the complete `E,W` support-two coefficient plane is empty.

## 1. Coprime-weight invariant coordinates

Use the target-translated coordinates of THM-2411:

```text
y=11s,                   u=dT,                   Z=T^2,

wt(y,u,Z,B,C,D,E,W)=(1,2,3,2,3,4,5,6).              (2)
```

The cases `E=0` and `W=0` are the closed `W` and `E` axes of THM-2423
and THM-2425.  Suppose from now on that

```text
E,W!=0.                                                (3)
```

First take `y!=0` in `C(x)` and define the constant weighted ratio and the
moving root-free coordinate

```text
lambda=W^5/E^6 in C*,             r=W/(Ey).           (4)
```

Bezout arithmetic for weights five and six gives the lossless identities

```text
E/y^5=r^5/lambda,                  W/y^6=r^6/lambda.  (5)
```

As usual put

```text
v=u/y^2,                           zeta=Z/y^3.         (6)
```

Multiplying the first two normalized fluxes of THM-2411 by `lambda` gives

```text
F_1
 =1331lambda(63-1089v)zeta
   +4[-3748096r^5
      +lambda(922383v^2-25410v+63)]
 =0,                                                        (7)

F_2
 =lambda[15944049zeta^2+(-162339408v+2236080)zeta
          -1190488992v^3+147581280v^2-1219680v+672]
   -239878144r^5-1319329792r^6
 =0.                                                        (8)
```

The open first-flux chart is

```text
63-1089v!=0,                 equivalently v!=7/121.   (9)
```

Thus (7) reconstructs `zeta` uniquely.

## 2. The 13-term eliminant

Exact elimination gives

```text
Res_zeta(F_1,F_2)=255104784lambda P_lambda(r,v),       (10)
```

where

```text
P_lambda
 =14048223625216r^10

  +lambda(-10865422960128v^2+1257156375552v
           -36364027392)r^6

  +lambda(4938828618240v^2-571434716160v
           +3935500800)r^5

  -63lambda^2 L_5(v),                                  (11)
```

and

```text
L_5(v)
 =155624547606v^5+3215383215v^4-1700698560v^3
   +58124770v^2-855470v+2583.                         (12)
```

In particular the fixed section survives:

```text
P_lambda(0,v)=-63lambda^2L_5(v).                       (13)
```

Since the active weights are coprime, there is no nontrivial power cover
above this section.  The required genus will instead come from (11).

## 3. Uniform absolute irreducibility

For every `lambda!=0`, the polynomial `P_lambda` is absolutely irreducible
in `C[r,v]`.

Its Newton polygon is

```text
N=conv{(0,0),(10,0),(0,5)}=5 Delta,

Delta=conv{(0,0),(2,0),(0,1)}.                        (14)
```

Every lattice Minkowski summand of this triangle is, up to translation,
`k Delta` for an integer `0<=k<=5`: a proper summand has the same three
edge directions, and their primitive lattice lengths force the common
integer scale.  The coefficients of `v^5` and of the constant monomial are
nonzero for `lambda!=0`.  Gauss normalization therefore makes every factor
monic in `v`, removes translations, and excludes point summands.

If (11) factored, one factor would have `v`-degree one or two.  The complete
degree-one ansatz is

```text
v+ar^2+br+c.                                           (15)
```

Substitution in (11), division in `v`, and collection in `r` gives the unit
coefficient ideal

```text
I_1=(1) in Q[a,b,c,lambda].                            (16)
```

The complete degree-two ansatz is

```text
v^2+(ar^2+br+c)v+dr^4+er^3+fr^2+gr+h.                (17)
```

Its exact remainder coefficients likewise give

```text
I_2=(1) in Q[a,b,c,d,e,f,g,h,lambda].                 (18)
```

The exact companion reconstructs both reduced Groebner bases.  Hence neither
factor can occur at any physical ratio, proving absolute irreducibility.
Let `C_lambda` be the smooth projective normalization.  Projection to the
`r`-line has degree five.

## 4. Exact ramification divisor

Define the quintic wall polynomial

```text
W_5(r,lambda)=234256r^5-105lambda.                    (19)
```

The squared factor below is not usable branch mass.  Direct substitution
shows

```text
P_lambda(r,7/121)=256W_5(r,lambda)^2,                 (20)
```

so its roots are exactly intersections with the excluded first-flux wall.

The exact quintic discriminant factors as

```text
Disc_v(P_lambda)
 =c lambda^8 W_5(r,lambda)^2 K_30(r,lambda),          (21)
```

for `c in Q*`.  Normalize `K_30` to be primitive with

```text
K_30(0,lambda)=3469890498046875lambda^6.              (22)
```

It has degree thirty in `r`, degree six in `lambda`, and 22 nonzero terms.
Its leading coefficient is

```text
[r^30]K_30
 =5164096645133820805624281694208
    (24057lambda-1225000).                            (23)
```

Put

```text
A_1=56133lambda+327680000.                            (24)
```

Exact factorization over `Q[lambda]` gives

```text
Disc_r(K_30)
 =c' lambda^174 A_1^3 P_4^2 Q_4^3 R_6,              (25)
```

where `P_4,Q_4,R_6` are the unique primitive positive-leading irreducible
factors of degrees `4,4,6` occurring with exponents `2,3,1`.  Equations
(21)--(25) and the exact companion specify their coefficients exactly; no
numerical root approximation enters the proof.

The collision with the excluded wall is independently exact:

```text
Res_r(K_30,W_5)
 =c'' lambda^30 A_1^3 S_2(lambda),                    (26)
```

where

```text
S_2
 =1605730621565435904lambda^2
   -4390842515926741875lambda
   -6871947673600000000.                              (27)
```

The six polynomials

```text
24057lambda-1225000,
A_1,P_4,Q_4,R_6,S_2                                  (28)
```

are squarefree and pairwise coprime.  Thus degree drop, root collision, and
wall collision never overlap, except for the deliberately shared `A_1`
between (25) and (26).

At the sole degree-drop ratio

```text
lambda=1225000/24057,                                 (29)
```

`K_30` becomes a squarefree degree-29 polynomial and remains coprime to
`W_5`.

## 5. Uniform ramification floor

We use the elementary Sylvester bound from THM-2469, repeated here to keep
the dependency local.

**Lemma.**  Let `k_s(r)` be a degree-`n` polynomial over a
characteristic-zero DVR whose leading coefficient is a unit.  If its
discriminant has valuation `e`, then

```text
deg gcd(k_0,partial_r k_0)<=e.                        (30)
```

Indeed, the resultant `Res(k_s,partial_r k_s)` is the determinant of the
Sylvester matrix.  Its reduction has nullity equal to the gcd degree.  In
Smith normal form, each lost rank contributes at least one uniformizer to
the determinant.

Every nonzero exceptional factor in (25) is squarefree and has exponent at
most three.  Coprimality with (23) makes the leading coefficient a unit at
all such parameters.  Therefore every degree-30 exceptional fibre has

```text
d=deg gcd(K_30,partial_rK_30)<=3.                     (31)
```

If its nonsimple roots have multiplicities `m_i>=2`, then

```text
d=sum_i(m_i-1),

sum_i m_i=d+#{i}<=2d.                                 (32)
```

Consequently at least

```text
30-2d>=24                                             (33)
```

roots of `K_30` are simple.  At an `A_1` root, at most the five roots of
`W_5` must also be discarded, leaving at least nineteen simple roots away
from the wall.  At a `P_4,Q_4`, or `R_6` root there is no wall collision.
At an `S_2` root, `K_30` is squarefree and at most five roots lie on the
wall.  At (29), all 29 finite roots are simple and none lies on the wall.
The generic case is squarefree and wall-disjoint.  Equations (23)--(28)
exhaust every `lambda in C*`, so uniformly

```text
at least 19 simple finite roots of K_30 lie off W_5.   (34)
```

Each such root contributes one unit of tame ramification to the normalized
degree-five cover.  Indeed, (21) has valuation one there; normalization can
change an order discriminant only by twice the index valuation, so the order
is already maximal and the local ramification is simple.

The total ramification `R` of a degree-five map satisfies

```text
2g(C_lambda)-2=-10+R.                                 (35)
```

Both terms on the left and `-10` are even, hence `R` is even.  The nineteen
visible simple contributions in (34) therefore force

```text
R>=20,                    g(C_lambda)>=6.             (36)
```

This parity step is useful precisely because the wall has odd degree.

## 6. Genus and trajectory closure

The rational pair `(r,v)` from (4),(6) would give a nonconstant rational map
from `P^1` to the positive-genus normalization `C_lambda`, impossible.
Hence `r,v` are constant.  Since `r=W/(Ey)` is nonzero, `y` and then `u`
are constant.  Equation (7) reconstructs constant `zeta`, so `Z,T,q` are
constants.  This contradicts the genuine deck, which fixes the constant
field and sends the nonzero `q` to `-q`.

It remains to treat `y=0`.  The original first flux reduces to

```text
-1449459uZ-14992384E=0.                               (37)
```

The open chart gives `u!=0`, and `E!=0`, so

```text
uZ=k:=-14992384E/1449459 in C*.                       (38)
```

The second flux is

```text
15944049Z^2-1319329792W-1190488992u^3=0.             (39)
```

Substitute `Z=k/u` and multiply by `u^2`.  This is the constant-field
quintic

```text
1190488992u^5+1319329792Wu^2-15944049k^2=0.          (40)
```

Thus `u`, and then `Z,T,q`, are constants, giving the same deck
contradiction.  Together with the two axes, this proves (1).

## 7. Scope and structural lesson

This theorem closes the sixth of ten support-two planes.  The open pairs are

```text
(B,C), (B,E), (C,E), (D,E).                           (41)
```

All have coprime weights.  The `C,D` and `E,W` proofs now show two distinct
Newton shapes for the replacement operation: a quadrilateral may force
small factors through edge arithmetic, while a dilated primitive triangle
forces them through homothetic summands.  In both cases the fixed `L_5`
section survives, but genus is certified by the moving discriminant after
the excluded wall is removed.

Higher mixed strata, branches outside the inherited reduction, split/even
short edges, and integral order raising remain open.  Nothing here proves
`JC(2)` or `DC(2)`.

## 8. Exact verification

Run

```bash
python3 04-computation/jc2_degree22_ew_plane_ramification_thm2470.py
python3 -O 04-computation/jc2_degree22_ew_plane_ramification_thm2470.py
```

The companion reconstructs (7)--(13), the Newton triangle and both
exhaustive factor ideals, (19)--(29), every squarefreeness and coprimality
control, the degree-drop fibre, the branch/genus floors, and the `y=0`
equations.  Normal and optimized transcripts byte-match the stored output.
All truth-bearing checks use explicit exceptions and remain active under
optimized Python.  **QED.**
