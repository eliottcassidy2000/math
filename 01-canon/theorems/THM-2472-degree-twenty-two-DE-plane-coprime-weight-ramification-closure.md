---
id: THM-2472
title: "Degree-twenty-two D-E plane coprime-weight ramification closure"
status: >
  PROVED + VERIFIED-EXACT. In the open first-flux chart of the genuine
  nonsplit polynomial exact-square-prefix degree-twenty-two branch, the
  complete coefficient plane B=C=W=0 is empty. For D,E nonzero, the
  root-free invariants lambda=E^4/D^5 and r=E/(Dy) give a 16-term
  plane curve of bidegree (10,5). Its Newton triangle reduces every
  factorization to one linear or quadratic v-factor, and both exact
  coefficient ideals are unit, so every physical fibre is absolutely
  irreducible. The v-discriminant is
  constant*lambda^7*W_5^2*K_30, where W_5 is exactly the excluded
  first-flux wall. A complete exceptional-parameter factorization,
  together with a Sylvester-valuation bound and ramification parity,
  leaves total ramification at least twenty in every fibre.
  Riemann--Hurwitz gives normalization genus at least six. No rational
  trajectory survives. Together with the axes, this closes the seventh
  of ten support-two planes. It does not close the other three planes,
  higher mixed strata, JC(2), or DC(2).
source: codex-2026-07-27-degree-twenty-two-DE-ramification
depends_on:
  - THM-2411-degree-twenty-two-first-flux-pole-divisor-square-class-reduction
  - THM-2425-degree-twenty-two-CDE-axis-hyperelliptic-closure
related:
  - THM-2429-degree-twenty-two-CW-plane-hyperelliptic-family-closure
  - THM-2437-degree-twenty-two-DW-plane-quartic-ramification-closure
  - THM-2463-degree-twenty-two-BD-plane-square-lift-closure
  - THM-2468-degree-twenty-two-BW-plane-square-lift-closure
  - THM-2469-degree-twenty-two-CD-plane-coprime-weight-ramification-closure
  - THM-2470-degree-twenty-two-EW-plane-coprime-weight-ramification-closure
script: 04-computation/jc2_degree22_de_plane_ramification_thm2472.py
output: 05-knowledge/results/jc2_degree22_de_plane_ramification_thm2472.out
script_sha256: f94615ebfc8c00d994b60e935eb8c8786aa60bb89aee73a1880cfb0b119d627d
output_sha256: 84b2fa686c3dd02ac8e667145b34d61566201f5cfc98dcea1d48fd7a63ffbbe7
hash_basis: working-tree bytes (LF)
---

# THM-2472 -- the degree-twenty-two D-E plane is empty

**PROVED + VERIFIED-EXACT.**

THM-2469 and THM-2470 replaced the exhausted retained-power operation by a
root-free coordinate and exact ramification on the base curve.  Adjacent
weights four and five provide the same lossless coordinates as weights five
and six, but both coefficient contributions remain visible in the excluded
wall.  The exact conclusion is

```text
genuine degree-22 trajectory,
mathcal A!=0,
B=C=W=0
    => contradiction.                                           (1)
```

Thus the complete `D,E` support-two coefficient plane is empty.

## 1. Coprime-weight invariant coordinates

Use the target-translated coordinates of THM-2411:

```text
y=11s,                   u=dT,                   Z=T^2,

wt(y,u,Z,B,C,D,E,W)=(1,2,3,2,3,4,5,6).              (2)
```

The cases `D=0` and `E=0` are the closed `E` and `D` axes of THM-2425.
Suppose from now on that

```text
D,E!=0.                                                (3)
```

First take `y!=0` in `C(x)` and define

```text
lambda=E^4/D^5 in C*,             r=E/(Dy).           (4)
```

Bezout arithmetic for weights four and five gives the lossless identities

```text
D/y^4=r^4/lambda,                  E/y^5=r^5/lambda.  (5)
```

Put

```text
v=u/y^2,                           zeta=Z/y^3.         (6)
```

Multiplying the first two normalized fluxes of THM-2411 by `lambda` gives

```text
F_1
 =1331lambda(63-1089v)zeta
   +4[511104r^4-3748096r^5
      +lambda(922383v^2-25410v+63)]
 =0,                                                        (7)

F_2
 =lambda[15944049zeta^2+(-162339408v+2236080)zeta
          -1190488992v^3+147581280v^2-1219680v+672]
   +(-1978994688v+16355328)r^4-239878144r^5
 =0.                                                        (8)
```

The open first-flux chart is

```text
63-1089v!=0,                 equivalently v!=7/121.   (9)
```

Thus (7) reconstructs `zeta` uniquely.

## 2. The 16-term eliminant

Exact elimination gives

```text
Res_zeta(F_1,F_2)=255104784lambda P_lambda(r,v),       (10)
```

where

```text
P_lambda
 =14048223625216r^10-3831333715968r^9+261227298816r^8

  +lambda(4938828618240v^2-571434716160v
           +3935500800)r^5

  +lambda(-16298134440192v^3+1077562607616v^2
           +38961457920v-987452928)r^4

  -63lambda^2 L_5(v),                                  (11)
```

and

```text
L_5(v)
 =155624547606v^5+3215383215v^4-1700698560v^3
   +58124770v^2-855470v+2583.                         (12)
```

Thus the fixed section remains

```text
P_lambda(0,v)=-63lambda^2L_5(v).                       (13)
```

The coprime active weights leave no retained power cover.  Genus must be
read from (11) itself.

## 3. Uniform absolute irreducibility

For every `lambda!=0`, the polynomial `P_lambda` is absolutely irreducible
in `C[r,v]`.

Its Newton polygon is the same dilated primitive triangle as in THM-2470:

```text
N=conv{(0,0),(10,0),(0,5)}=5 Delta,

Delta=conv{(0,0),(2,0),(0,1)}.                        (14)
```

Every lattice Minkowski summand is, up to translation, `k Delta` for an
integer `0<=k<=5`.  For `lambda!=0`, the `v^5` and constant coefficients are
nonzero.  Gauss normalization therefore makes factors monic in `v`, removes
translations, and excludes point summands.  Any factorization has a factor
of `v`-degree one or two.

The complete degree-one ansatz is

```text
v+ar^2+br+c.                                           (15)
```

Substitution in (11), division in `v`, and collection in `r` gives

```text
I_1=(1) in Q[a,b,c,lambda].                            (16)
```

The complete degree-two ansatz is

```text
v^2+(ar^2+br+c)v+dr^4+er^3+fr^2+gr+h.                (17)
```

Its exact remainder coefficient ideal is likewise

```text
I_2=(1) in Q[a,b,c,d,e,f,g,h,lambda].                 (18)
```

The exact companion reconstructs both reduced Groebner bases.  Hence (11)
is absolutely irreducible at every physical ratio.  Let `C_lambda` be its
smooth projective normalization.  Projection to the `r`-line has degree
five.

## 4. Exact ramification divisor

Define

```text
W_5(r,lambda)=234256r^5-31944r^4-105lambda.           (19)
```

Direct substitution shows

```text
P_lambda(r,7/121)=256W_5(r,lambda)^2,                 (20)
```

so this square records only intersections with the excluded first-flux
wall, not usable branch mass.

The exact quintic discriminant factors as

```text
Disc_v(P_lambda)
 =c lambda^7 W_5(r,lambda)^2 K_30(r,lambda),          (21)
```

for `c in Q*`.  Normalize `K_30` to be primitive with

```text
K_30(0,lambda)=-24289233486328125lambda^7.            (22)
```

It has degree thirty in `r`, degree seven in `lambda`, and 31 nonzero terms.
Its leading coefficient is

```text
[r^30]K_30
 =3755706651006415131363113959424
    (11790625lambda+2519424).                         (23)
```

Put

```text
A_1=1203125lambda+1944.                               (24)
```

Exact factorization over `Q[lambda]` gives

```text
Disc_r(K_30)
 =c' lambda^190 A_1^3 P_6^2 Q_6^3 R_9,              (25)
```

where `P_6,Q_6,R_9` are the unique primitive positive-leading irreducible
factors of degrees `6,6,9` occurring with exponents `2,3,1`.  Equations
(21)--(25) and the exact companion specify their coefficients exactly.

The collision with the excluded wall is

```text
Res_r(K_30,W_5)
 =c'' lambda^29 A_1^3 S_3(lambda),                    (26)
```

where

```text
S_3
 =679985152000000000000lambda^3
   -128456167991962809796875lambda^2
   -48607662933281282414472lambda
   -5052710781009410264064.                           (27)
```

The six polynomials

```text
11790625lambda+2519424,
A_1,P_6,Q_6,R_9,S_3                                  (28)
```

are squarefree and pairwise coprime.  The only overlap between root
collision and wall collision is the displayed `A_1`.

At the sole degree-drop ratio

```text
lambda=-2519424/11790625,                              (29)
```

`K_30` becomes a squarefree degree-29 polynomial and remains coprime to
`W_5`.

## 5. Uniform ramification floor

We repeat the local Sylvester bound used in THM-2469/2470.

**Lemma.**  Let `k_s(r)` be a degree-`n` polynomial over a
characteristic-zero DVR whose leading coefficient is a unit.  If its
discriminant has valuation `e`, then

```text
deg gcd(k_0,partial_r k_0)<=e.                        (30)
```

The resultant `Res(k_s,partial_r k_s)` is the determinant of the Sylvester
matrix.  Its reduction has nullity equal to the gcd degree; Smith normal
form makes every lost rank contribute at least one uniformizer.

Every nonzero exceptional factor in (25) is squarefree and appears with
exponent at most three.  Coprimality with (23) makes the leading coefficient
a unit there.  Hence each degree-30 exceptional fibre has

```text
d=deg gcd(K_30,partial_rK_30)<=3.                     (31)
```

If the nonsimple roots have multiplicities `m_i>=2`, then

```text
d=sum_i(m_i-1),

sum_i m_i=d+#{i}<=2d.                                 (32)
```

At least

```text
30-2d>=24                                             (33)
```

roots of `K_30` are therefore simple.  At an `A_1` root, at most the five
roots of `W_5` must be discarded, leaving at least nineteen simple roots
off the wall.  At a `P_6,Q_6`, or `R_9` root there is no wall collision.  At
an `S_3` root, `K_30` is squarefree and at most five roots lie on the wall.
At (29), all 29 finite roots are simple and wall-disjoint.  The generic case
is squarefree and wall-disjoint.  Equations (23)--(28) exhaust every
`lambda in C*`, so uniformly

```text
at least 19 simple finite roots of K_30 lie off W_5.   (34)
```

At each such root, (21) has discriminant valuation one.  Normalization
changes an order discriminant only by twice the index valuation, so the
order is already maximal and tame local ramification is simple.

For the degree-five projection,

```text
2g(C_lambda)-2=-10+R.                                 (35)
```

The total ramification `R` is even.  Thus the nineteen visible contributions
force

```text
R>=20,                    g(C_lambda)>=6.             (36)
```

## 6. Genus and trajectory closure

The rational pair `(r,v)` from (4),(6) would give a nonconstant rational map
from `P^1` to the positive-genus normalization `C_lambda`, impossible.
Hence `r,v` are constant.  Since `r=E/(Dy)` is nonzero, `y` and then `u`
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
15944049Z^2-1978994688Du-1190488992u^3=0.            (39)
```

Substituting `Z=k/u` and multiplying by `u^2` gives the constant-field
quintic

```text
1190488992u^5+1978994688Du^3-15944049k^2=0.          (40)
```

Thus `u`, and then `Z,T,q`, are constants, giving the same deck
contradiction.  Together with the two axes, this proves (1).

## 7. Scope and structural lesson

This theorem closes the seventh of ten support-two planes.  The open pairs
are

```text
(B,C), (B,E), (C,E).                                  (41)
```

The three coprime-weight closures now expose a stable mechanism.  The
Bezout coordinate preserves the fixed `L_5` section; a dilated primitive
Newton triangle gives an exhaustive small-factor test; the moving
discriminant separates as excluded wall squared times a degree-30 divisor;
and a coarse collision-valuation bound already leaves enough branch mass.
The coefficients and exceptional factors change, but the proof invoice does
not.

Higher mixed strata, branches outside the inherited reduction, split/even
short edges, and integral order raising remain open.  Nothing here proves
`JC(2)` or `DC(2)`.

## 8. Exact verification

Run

```bash
python3 04-computation/jc2_degree22_de_plane_ramification_thm2472.py
python3 -O 04-computation/jc2_degree22_de_plane_ramification_thm2472.py
```

The companion reconstructs (7)--(13), the Newton triangle and both
exhaustive factor ideals, (19)--(29), every squarefreeness and coprimality
control, the degree-drop fibre, the branch/genus floors, and the `y=0`
quintic.  Normal and optimized transcripts byte-match the stored output.
All truth-bearing checks use explicit exceptions and remain active under
optimized Python.  **QED.**
