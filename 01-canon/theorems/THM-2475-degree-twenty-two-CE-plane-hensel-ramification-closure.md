---
id: THM-2475
title: "Degree-twenty-two C-E plane Hensel-ramification closure"
status: >
  PROVED + VERIFIED-EXACT. In the open first-flux chart of the genuine
  nonsplit polynomial exact-square-prefix degree-twenty-two branch, the
  complete coefficient plane B=D=W=0 is empty. For C,E nonzero, the
  root-free invariants lambda=E^3/C^5 and r=E^2/(C^3y) give an 18-term
  plane curve of bidegree (10,5). Its Newton triangle reduces every
  factorization to one linear or quadratic v-factor. The linear ideal
  is unit. For a quadratic factor, the squarefree fixed L_5 section and
  the r,r^2 coefficient gap force its first two Hensel deformations to
  vanish; the two r^3 equations are incompatible on all ten quadratic
  divisors of L_5. Thus every physical fibre is absolutely irreducible.
  The v-discriminant is constant*lambda^22*W_5^2*K_30, where W_5 is
  exactly the excluded first-flux wall. A complete exceptional-parameter
  factorization, a Sylvester-valuation bound, and ramification parity
  force normalization genus at least six in every fibre. No rational
  trajectory survives. Together with the axes, this closes the eighth
  of ten support-two planes. It does not close the B-C or B-E plane,
  higher mixed strata, JC(2), or DC(2).
source: codex-2026-07-27-degree-twenty-two-CE-hensel-ramification
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
  - THM-2472-degree-twenty-two-DE-plane-coprime-weight-ramification-closure
script: 04-computation/jc2_degree22_ce_plane_hensel_ramification_thm2475.py
output: 05-knowledge/results/jc2_degree22_ce_plane_hensel_ramification_thm2475.out
script_sha256: 8b2a8abf48064ac9ebe347bf10f1b83053433a910c5b0cf6a91fe7fb7dc67d12
output_sha256: 5e5e3828ad4a132f954200221a3bb3ea9b7e0a39e06843ea09025849f3d0ce19
hash_basis: working-tree bytes (LF)
---

# THM-2475 -- the degree-twenty-two C-E plane is empty

**PROVED + VERIFIED-EXACT.**

The root-free ramification method of THM-2469/2470/2472 extends to the
nonadjacent coprime weights three and five, but its efficient
irreducibility proof changes.  Instead of a large unrestricted quadratic
factor elimination, the fixed squarefree section gives a finite etale
Hensel problem whose first nonzero deformation is already obstructed.  The
exact conclusion is

```text
genuine degree-22 trajectory,
mathcal A!=0,
B=D=W=0
    => contradiction.                                           (1)
```

Thus the complete `C,E` support-two coefficient plane is empty.

## 1. Coprime-weight invariant coordinates

Use the target-translated coordinates of THM-2411:

```text
y=11s,                   u=dT,                   Z=T^2,

wt(y,u,Z,B,C,D,E,W)=(1,2,3,2,3,4,5,6).              (2)
```

The cases `C=0` and `E=0` are the closed `E` and `C` axes of THM-2425.
Suppose

```text
C,E!=0.                                                (3)
```

First take `y!=0` in `C(x)` and define

```text
lambda=E^3/C^5 in C*,             r=E^2/(C^3y).      (4)
```

The Bezout identity `2*5-3*3=1` gives the lossless formulas

```text
C/y^3=r^3/lambda^2,                E/y^5=r^5/lambda^3. (5)
```

Put

```text
v=u/y^2,                           zeta=Z/y^3.         (6)
```

Multiplying the first two normalized fluxes of THM-2411 by `lambda^3`
gives

```text
F_1
 =1331lambda^3(63-1089v)zeta
   +4[lambda(2342560v-58080)r^3-3748096r^5
      +lambda^3(922383v^2-25410v+63)]
 =0,                                                        (7)

F_2
 =lambda^3[15944049zeta^2+(-162339408v+2236080)zeta
             -1190488992v^3+147581280v^2-1219680v+672]
   +lambda(-206145280zeta+449771520v-1239040)r^3
   -239878144r^5
 =0.                                                        (8)
```

The open first-flux chart is

```text
63-1089v!=0,                 equivalently v!=7/121.   (9)
```

Thus (7) reconstructs `zeta` uniquely.

## 2. The 18-term eliminant

Exact elimination gives

```text
Res_zeta(F_1,F_2)=255104784lambda^3 P_lambda(r,v),    (10)
```

where

```text
P_lambda
 =14048223625216r^10-580505108480lambda r^8

  +lambda^2(-5487587353600v^2+634927462400v
             -12368716800)r^6

  +lambda^3(4938828618240v^2-571434716160v
             +3935500800)r^5

  +lambda^4(-4938828618240v^3+537420744960v^2
             -12593602560v+146361600)r^3

  -63lambda^6 L_5(v),                                  (11)
```

and

```text
L_5(v)
 =155624547606v^5+3215383215v^4-1700698560v^3
   +58124770v^2-855470v+2583.                         (12)
```

The fixed section and the new coefficient gap are

```text
P_lambda(0,v)=-63lambda^6L_5(v),

[r]P_lambda=[r^2]P_lambda=0.                          (13)
```

The quintic `L_5` is squarefree.

## 3. Uniform absolute irreducibility

For every `lambda!=0`, the polynomial `P_lambda` is absolutely irreducible
in `C[r,v]`.

Its Newton polygon is

```text
N=conv{(0,0),(10,0),(0,5)}=5 Delta,

Delta=conv{(0,0),(2,0),(0,1)}.                        (14)
```

As in THM-2470/2472, every nontrivial lattice Minkowski summand is `k Delta`
up to translation.  The nonzero constant and `v^5` coefficients normalize
factors to be monic in `v` and remove translations.  Any factorization has
a factor of `v`-degree one or two.

For a linear factor the complete ansatz is

```text
v+ar^2+br+c.                                           (15)
```

Exact division gives the unit coefficient ideal

```text
I_1=(1) in Q[a,b,c,lambda].                            (16)
```

For a quadratic factor, write

```text
q(r,v)=q_0(v)+rq_1(v)+r^2q_2(v)+r^3q_3(v)+... .      (17)
```

At `r=0`, the monic quadratic

```text
q_0=v^2+cv+h                                           (18)
```

divides `L_5`.  Because `L_5` is squarefree, `q_0` is coprime to its cubic
cofactor.  At orders one and two, (13) gives

```text
q_0s_i+q_i s_0=0.                                     (19)
```

Thus `q_0` divides `q_i`; since `deg_v q_i<2`, one has

```text
q_1=q_2=0.                                            (20)
```

The Newton support `2 Delta` allows no `v r^3` monomial, so

```text
q_3=e in C.                                           (21)
```

This reduces the whole quadratic-factor question to the first nonzero
Hensel equation.  Let

```text
rem_v(L_5,v^2+cv+h)=J_1(c,h)v+J_0(c,h),               (22)
```

where

```text
J_1
 =155624547606c^4-3215383215c^3-466873642818c^2h
   -1700698560c^2+6430766430ch-58124770c
   +155624547606h^2+1700698560h-855470,

J_0
 =155624547606c^3h-3215383215c^2h-311249095212ch^2
   -1700698560ch+3215383215h^2-58124770h+2583.        (23)
```

Put `x=lambda^2 e` and define

```text
A=-669650058c^2+9223830c+446433372h+2439360,
N=112442880c^2+12235520c-112442880h+286720,

B=27009219006c^3-558041715c^2-108036876024ch
  -295162560c+1116083430h-10087770,
M=13605588480ch+1480497920h-403200.                  (24)
```

The `v r^3` and `r^3` remainder coefficients, after removing nonzero scalar
factors, are exactly

```text
Ax+N,                         Bx+M.                  (25)
```

Exact Groebner reduction gives

```text
(J_1,J_0,Ax+N,Bx+M)=(1) in Q[x,c,h].                 (26)
```

Thus none of the ten quadratic divisors of the fixed squarefree quintic can
even begin the required polynomial Hensel lift.  Linear and quadratic
factors are both impossible, proving absolute irreducibility.  Let
`C_lambda` be the smooth projective normalization.  Projection to the
`r`-line has degree five.

## 4. Exact ramification divisor

Define

```text
W_5(r,lambda)=234256r^5-4840lambda r^3-105lambda^3.   (27)
```

Direct substitution shows

```text
P_lambda(r,7/121)=256W_5(r,lambda)^2,                 (28)
```

so this square is exactly the excluded first-flux wall.

The quintic discriminant factors as

```text
Disc_v(P_lambda)
 =c lambda^22 W_5(r,lambda)^2 K_30(r,lambda),         (29)
```

for `c in Q*`.  Normalize `K_30` to be primitive with

```text
K_30(0,lambda)=-209858977321875lambda^20.             (30)
```

It has degree thirty in `r`, degree twenty in `lambda`, and 39 nonzero
terms.  Its leading coefficient is

```text
[r^30]K_30
 =19752025963219313237822537728
    (19370043lambda^2-12500).                         (31)
```

Put

```text
A_2=881552179200lambda^2-6987070464lambda+15353125.   (32)
```

Exact factorization gives

```text
Disc_r(K_30)
 =c' lambda^580 A_2^3 P_8^2 Q_8^3 R_12,             (33)
```

where `P_8,Q_8,R_12` are the unique primitive positive-leading irreducible
factors of degrees `8,8,12` occurring with exponents `2,3,1`.

The wall collision divisor is

```text
Res_r(K_30,W_5)
 =c'' lambda^90 A_2^3 S_4(lambda),                    (34)
```

where

```text
S_4
 =8665219524371609991249920000lambda^4
  -16387179390132006347512532172800lambda^3
  -775731974166149050267099107687lambda^2
  -12728746868025855794767200000lambda
  -86790242512199520000000000.                        (35)
```

The six polynomials

```text
19370043lambda^2-12500,
A_2,P_8,Q_8,R_12,S_4                                 (36)
```

are squarefree and pairwise coprime.  The only overlap between a root
collision and a wall collision is the displayed `A_2`.

There are two degree-drop ratios, the roots of

```text
19370043lambda^2-12500=0.                             (37)
```

At both, the `r^29` coefficient is nonzero.  The degree-30 discriminant and
the wall resultant are also nonzero there.  The standard leading-degree
specialization identity

```text
Disc_30(K)|_[r^30]=0=[r^29](K)^2 Disc_29(K)           (38)
```

therefore shows that both degree-29 fibres are squarefree and wall-disjoint.

## 5. Uniform ramification floor

Use the Sylvester bound from THM-2469/2470/2472: over a
characteristic-zero DVR, discriminant valuation `e` bounds

```text
deg gcd(k_0,partial_rk_0)<=e.                         (39)
```

Every nonzero exceptional factor in (33) is squarefree and has exponent at
most three; (31),(36) make the leading coefficient a unit there.  Thus each
degree-30 exceptional fibre has gcd degree at most three.  Its nonsimple
roots consume total multiplicity at most six, leaving at least

```text
30-6=24                                               (40)
```

simple roots of `K_30`.

At an `A_2` root, at most the five roots of `W_5` must be discarded, leaving
nineteen simple roots off the wall.  At a `P_8,Q_8`, or `R_12` root there is
no wall collision.  At an `S_4` root, `K_30` is squarefree and at most five
roots lie on the wall.  The two degree-29 fibres are squarefree and
wall-disjoint by (38).  The generic case is also squarefree and
wall-disjoint.  Hence uniformly

```text
at least 19 simple finite roots of K_30 lie off W_5.   (41)
```

At each, (29) has discriminant valuation one.  The order is already maximal
because normalization changes discriminant valuation by twice the index,
so tame local ramification is simple.

For the degree-five projection,

```text
2g(C_lambda)-2=-10+R.                                 (42)
```

The total ramification `R` is even.  Nineteen visible contributions force

```text
R>=20,                    g(C_lambda)>=6.             (43)
```

## 6. Genus and trajectory closure

The rational pair `(r,v)` from (4),(6) would give a nonconstant rational map
from `P^1` to the positive-genus normalization, impossible.  Hence `r,v`
are constant.  Since `r=E^2/(C^3y)` is nonzero, `y` and then `u` are
constant.  Equation (7) reconstructs constant `zeta`, so `Z,T,q` are
constants.  This contradicts the genuine deck.

It remains to treat `y=0`.  The original first flux is

```text
-1449459uZ+9370240Cu-14992384E=0.                    (44)
```

The open chart gives `u!=0`; thus

```text
Z=alpha+beta/u,

alpha=9370240C/1449459,
beta=-14992384E/1449459.                              (45)
```

The second flux is

```text
15944049Z^2-206145280CZ-1190488992u^3=0.             (46)
```

After (45) and multiplication by `u^2`, this is the nonzero constant-field
quintic

```text
1190488992u^5
 +206145280C(alpha u^2+beta u)
 -15944049(alpha u+beta)^2=0.                         (47)
```

Thus `u`, and then `Z,T,q`, are constants, again contradicting the deck.
Together with the two axes, this proves (1).

## 7. Scope and structural lesson

This theorem closes the eighth of ten support-two planes.  Only

```text
(B,C), (B,E)                                           (48)
```

remain open.

The fixed section is more than a source of possible branch points.  When a
sparse moving eliminant has an initial coefficient gap, squarefreeness makes
each factorization of the fixed section lift uniquely.  A proposed small
Newton factor can then fail at its first allowed deformation order.  Here
this replaces a stalled nine-variable quadratic-factor elimination by the
three-variable unit ideal (26).  The operation is likely reusable on the two
remaining planes and on higher-support strata with sparse boundary sections.

Higher mixed strata, branches outside the inherited reduction, split/even
short edges, and integral order raising remain open.  Nothing here proves
`JC(2)` or `DC(2)`.

## 8. Exact verification

Run

```bash
python3 04-computation/jc2_degree22_ce_plane_hensel_ramification_thm2475.py
python3 -O 04-computation/jc2_degree22_ce_plane_hensel_ramification_thm2475.py
```

The companion reconstructs (7)--(13), the Newton triangle, the line ideal,
the fixed-section Hensel obstruction (22)--(26), the complete discriminant
and wall divisors (27)--(38), all coprimality controls, the branch/genus
floors, and the `y=0` quintic.  Normal and optimized transcripts byte-match
the stored output.  All truth-bearing checks use explicit exceptions and
remain active under optimized Python.  **QED.**
