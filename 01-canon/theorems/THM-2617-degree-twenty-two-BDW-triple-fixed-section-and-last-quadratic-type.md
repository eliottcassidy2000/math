---
id: THM-2617
title: "Degree-twenty-two B-D-W triple fixed section and last quadratic type"
status: >
  PROVED CANDIDATE REDUCTION + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT
  PENDING. In the inherited genuine nonsplit degree-twenty-two branch on
  the open first-flux chart, the first support-three stratum B,D,W nonzero
  and C=E=0 has a 34-term two-parameter eliminant of bidegree (5,5). Its
  Newton polygon is the full degree-five triangle and its fixed p=0 section
  is the same squarefree quintic L_5 as on the closed B-D and B-W planes.
  Every line factor and two of the three possible quadratic top types are
  excluded by exact coefficient ideals. Reducibility is reduced to one
  explicit quadratic top type: the product of two moving-cubic roots. If
  this last type is absent, the physical square lift is connected with at
  least six branch places and genus at least two, closing the entire B-D-W
  stratum. A 324-parameter rational hostile sweep finds no reducible fibre,
  but uniform irreducibility and the stratum closure remain OPEN. No other
  mixed stratum, split/even descent, 2-adic raising, JC(2), or DC(2) follows.
source: klein-2026-07-28-degree22-bdw-triple
depends_on:
  - THM-2411-degree-twenty-two-first-flux-pole-divisor-square-class-reduction
  - THM-2463-degree-twenty-two-BD-plane-square-lift-closure
  - THM-2468-degree-twenty-two-BW-plane-square-lift-closure
related:
  - THM-2480-degree-twenty-two-BC-plane-hensel-ramification-closure
script: 04-computation/jc2_degree22_bdw_triple_last_quadratic_thm2617.py
output: 05-knowledge/results/jc2_degree22_bdw_triple_last_quadratic_thm2617.out
script_sha256: 9bb6239e12e443aea02ac7e04ffe08dca612ad35e9f6078d61c4f39fcf7ade27
output_sha256: 1434b9726b8c9a5fb33397ed22ce5ec28278bf938d87943f2f2696a5caea5468
hash_basis: working-tree bytes (LF)
---

# THM-2617 -- the B-D-W triple reduces to one quadratic at infinity

**PROVED CANDIDATE REDUCTION + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT
PENDING.**

The ten support-two coefficient planes are closed, but a support-three torus
has two independent weight-zero moduli.  The first low-complexity triple is

```text
B,D,W != 0,                    C=E=0,                         (1)
```

because its weights `(2,4,6)` admit the direct invariants `D/B^2,W/B^3`
and retain the physical square lift used on the `B,D` and `B,W` planes.

The calculation gives a sharp reduction, not a closure.  Four of the five
possible small top-factor types are impossible.  The remaining type has an
explicit two-root formula and coefficient ideal; that ideal is now the exact
highest-leverage algebraic target on this stratum.

## 1. Inherited coordinates and the two moduli

Use THM-2411's target-translated degree-twenty-two coordinates

```text
y=11s,                    u=dT,                    Z=T^2,

wt(y,u,Z,B,C,D,E,W)=(1,2,3,2,3,4,5,6).             (2)
```

Work in the genuine nonsplit branch with

```text
mathcal A=616B-1089u+63y^2 != 0.                    (3)
```

On (1), the identically-zero boundary `y=0` makes the first flux

```text
N_1=1331(616B-1089u)Z.                              (4)
```

The open chart then forces `Z=0`, contrary to `Z=T^2` and `T=q^2!=0`.
Hence `y!=0` in `C(x)`.

Define

```text
lambda=D/B^2 in C*,              mu=W/B^3 in C*,

p=B/y^2,                         v=u/y^2,
zeta=Z/y^3.                                               (5)
```

Then

```text
D/y^4=lambda p^2,                W/y^6=mu p^3.             (6)
```

Dividing the first two normalized fluxes by `y^5,y^6` gives

```text
F_1
 =2044416 lambda p^2-2981440pv+819896p zeta+24640p
  +3689532v^2-1449459v zeta-101640v+83853zeta+252=0,

F_2
 =-1319329792mu p^3-1978994688lambda p^2v
  +16355328lambda p^2+1443016960pv^2-71554560pv
  +65591680p zeta+98560p-1190488992v^3+147581280v^2
  -162339408v zeta-1219680v+15944049zeta^2
  +2236080zeta+672=0.                                      (7)
```

The open chart is

```text
616p-1089v+63 != 0,                                     (8)
```

so `F_1` reconstructs `zeta` uniquely.

## 2. The 34-term eliminant and fixed section

Exact elimination gives

```text
Res_zeta(F_1,F_2)=2^4 11^6 R_(lambda,mu)(p,v),            (9)
```

where `R` has exactly 34 terms and bidegree `(5,5)`.  Its `(p,v)` support is
every lattice point of

```text
Delta_5=conv{(0,0),(5,0),(0,5)}.                         (10)
```

The fixed section and first moving coefficient are independent of both
moduli:

```text
R(0,v)=-567L_5(v),                                       (11)

[p]R
 =11088(18649222647v^4+563356398v^3-136161300v^2
          +3108490v-17451),                              (12)
```

with

```text
L_5(v)
 =155624547606v^5+3215383215v^4-1700698560v^3
   +58124770v^2-855470v+2583.                             (13)
```

The quintic `L_5` is squarefree, and the polynomials in (11)--(12) are
coprime.  Thus `p=0` meets every fibre in five distinct smooth points, each
with valuation exactly one.

The top homogeneous form is

```text
-38974342(56p-99v)^2
 (256mu p^3+384lambda p^2v-280pv^2+231v^3).              (14)
```

This combines, without losing either modulus, the top factors from the
closed `B,D` and `B,W` planes.

## 3. Exhaustive small-factor inventory

If a total-degree-five polynomial is reducible over an algebraic closure,
one factor has degree one or two.  Since the Newton polygon is the full
triangle, the top form of that small factor is a product of one or two of
the five lines in (14).  There are exactly two line types and three
quadratic types:

```text
L1: one wall line,
L2: one root of the moving cubic;

Q1: wall times wall,
Q2: wall times one moving root,
Q3: two moving roots.                                    (15)
```

This inventory is exhaustive over the algebraic closure; no rationality of
the individual cubic roots is assumed.

## 4. Exact exclusion of four types

For `L1`, substitute

```text
p=(99/56)v+b.                                            (16)
```

All five coefficients of `R` in `v` generate the unit ideal in
`Q[b,lambda,mu]`.

For `L2`, choose one cubic root `h`.  Since `mu,h!=0`, its equation solves

```text
mu=(-384lambda h^2+280h-231)/(256h^3).                   (17)
```

Substitution `p=hv+b` gives a coefficient ideal whose Groebner basis
contains

```text
b^3,              b^2h,              bh^2,              h^3. (18)
```

Hence `h=0`, impossible because the cubic in (14) has constant term `231`.

For a monic quadratic factor use the complete ansatz

```text
p^2+(av+b)p+cv^2+dv+e.                                  (19)
```

For `Q1`, its top coefficients are

```text
a=-99/28,                    c=9801/3136.                (20)
```

The eleven remainder coefficients generate the unit ideal in
`Q[b,d,e,lambda,mu]`.

For `Q2`, use

```text
a=-(99/56+h),                 c=(99/56)h,                (21)
```

together with (17).  The exact remainder ideal in
`Q[b,d,e,lambda,h]` contains `h^3`, again impossible.  Thus `L1,L2,Q1,Q2`
are absent for every `lambda,mu in C*`.

## 5. The sole remaining quadratic type

Choose one root `h` of the moving cubic and put

```text
D_0=-384lambda h^2+280h-231=256mu h^3.                  (22)
```

The monic quadratic whose roots are the other two cubic roots has top form

```text
Q_3^top
 =p^2+[h(280h-231)/D_0]pv-[231h^2/D_0]v^2.             (23)
```

Both `h` and `D_0` are nonzero on the physical torus.  Add the lower terms
`bp+dv+e` from (19), substitute (17), and divide `R` by this quadratic.
The resulting nine explicit coefficient equations in

```text
b,d,e,lambda,h                                           (24)
```

are the complete remaining irreducibility obligation.

The two top-next equations are linear in `b,d`.  Their determinant factors
exactly as

```text
966306 h
 (384lambda h^2-560h+693)
 (384lambda h^2-280h+231)^5
 (114048lambda h^2-25088h^2-44352h+68607)^2.             (25)
```

The first and third displayed factors are already nonzero by (22).  Away
from the other two factors, (25) solves `b,d` uniquely and leaves a smaller
ideal in `e,lambda,h`; on either exceptional factor it gives two separate
boundary ideals.  None of those three ideals is declared empty here.

This is the exact stopping boundary:

```text
R_(lambda,mu) reducible
  iff a Q3 factor satisfying the nine equations exists.               (26)
```

## 6. Why the last factor would close the stratum

Suppose the `Q3` ideal is empty.  Sections 3--5 then make every physical
fibre `R_(lambda,mu)=0` absolutely irreducible.  Let `C_(lambda,mu)` be its
smooth projective normalization.

Restore the physical square class

```text
Y=y/sqrt(B),                    Y^2=1/p.                (27)
```

At the five points (11), `p` has valuation one.  Hence `1/p` is not a
square, the double cover is connected, and it ramifies at all five points.
The branch divisor of a quadratic extension has even degree, so it has at
least six places.  Riemann--Hurwitz gives

```text
g(double cover)>=-1+6/2=2.                              (28)
```

A rational Keller trajectory would give a map from `P^1` to this positive-
genus normalization, hence a constant map.  Then `p,v,Y,y,u` are constant;
the open first flux reconstructs `zeta`, and then `Z,T,q` are constant,
contradicting the genuine deck involution.  Thus emptiness of the one ideal
in (24) closes all of (1).

## 7. Finite hostile probe

As a hostile control only, the companion takes the 18 distinct nonzero
rationals

```text
{n/d : d in {1,2,3}, -4<=n<=4, n!=0}                  (29)
```

for each of `lambda,mu`, factors all `324` specialized bivariate
eliminants exactly over `Q`, and finds

```text
reducible fibres: 0.                                    (30)
```

This is **FINITE-EXACT** evidence.  It is not used in (26), does not exclude
algebraic parameter values, and is not a proof of uniform irreducibility.

## 8. Reproduction and scope

Run

```bash
python3 04-computation/jc2_degree22_bdw_triple_last_quadratic_thm2617.py
python3 -O 04-computation/jc2_degree22_bdw_triple_last_quadratic_thm2617.py
```

The companion reconstructs (7)--(14), the exhaustive top inventory, all
four exact exclusions, the last quadratic chart and determinant, the square-
lift genus invoice, the `y=0` boundary, and the frozen hostile universe.
All truth-bearing checks use `require` and remain active under optimized
Python.

The exact next target is the `Q3` remainder ideal in (24), split into the
generic determinant chart and the two genuine exceptional factors of (25).
Uniform irreducibility and the `B,D,W` closure remain **OPEN**.

Nothing here closes another support-three stratum, a branch outside the
inherited nonsplit exact-square-prefix reduction, split/even short edges, or
integral `2`-adic order raising.  It proves neither `JC(2)` nor `DC(2)`.

**QED for the reduction, pending independent hostile audit.**
