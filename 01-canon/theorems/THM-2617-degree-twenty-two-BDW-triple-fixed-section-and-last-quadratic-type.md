---
id: THM-2617
title: "Degree-twenty-two B-D-W triple pair-field Hensel closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In the inherited
  genuine nonsplit degree-twenty-two branch on
  the open first-flux chart, the support-three stratum B,D,W nonzero and
  C=E=0 has a 34-term two-parameter eliminant of bidegree (5,5). Its fixed
  p=0 section is the squarefree quintic L_5. Every line factor is excluded
  by exact coefficient ideals. Any quadratic factor chooses two roots of
  L_5; the ten unordered pairs form one reduced irreducible degree-ten
  pair field. Exact Hensel lifting in p uniquely solves the factor through
  order p^2, but its two p^3 remainder equations differ by a degree-nine
  element coprime to the pair-field modulus. Thus every physical fibre is
  absolutely irreducible. The physical square lift is connected with at
  least six branch places and genus at least two, so the entire B-D-W
  stratum is closed. No other mixed stratum, split/even descent, 2-adic
  raising, JC(2), or DC(2) follows.
source: klein-2026-07-28-degree22-bdw-triple
depends_on:
  - THM-2411-degree-twenty-two-first-flux-pole-divisor-square-class-reduction
  - THM-2463-degree-twenty-two-BD-plane-square-lift-closure
  - THM-2468-degree-twenty-two-BW-plane-square-lift-closure
related:
  - THM-2480-degree-twenty-two-BC-plane-hensel-ramification-closure
script: 04-computation/jc2_degree22_bdw_pair_hensel_closure_thm2617.py
output: 05-knowledge/results/jc2_degree22_bdw_pair_hensel_closure_thm2617.out
script_sha256: f4b92bfb3076ca85ea2da17f8821a0c787cca18cdea0ee5f7994486d94e35a0e
output_sha256: e57cba24899521acadb150360e67a49f4f36c2bf03d15d9b6623c03ff1278ea9
independent_script: 04-computation/jc2_degree22_bdw_triple_last_quadratic_thm2617.py
independent_output: 05-knowledge/results/jc2_degree22_bdw_triple_last_quadratic_thm2617.out
independent_script_sha256: c49330245c970481d3cb36a7a5354af0eaa3196b4c0ee499eb665fb264248bd9
independent_output_sha256: c5733177f0b883f4fd6c534394c22dba0608f6b8d89e33c6562de11d7dd41120
hash_basis: working-tree bytes (LF)
---

# THM-2617 -- the B-D-W triple closes over its fixed-section pair field

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The ten support-two coefficient planes are closed, but a support-three torus
has two independent weight-zero moduli.  The first low-complexity triple is

```text
B,D,W != 0,                    C=E=0,                         (1)
```

because its weights `(2,4,6)` admit the direct invariants `D/B^2,W/B^3`
and retain the physical square lift used on the `B,D` and `B,W` planes.

The calculation now closes this triple.  The first companion reduces the
only difficult top type to one compatibility divisor.  A second, structurally
different calculation starts from the fixed `p=0` section instead: a quadratic
factor must choose an unordered pair of its five roots.  Those ten choices form
one exact degree-ten field, and the attempted Hensel lift fails uniformly at
order `p^3`.

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

## 5. The sole remaining quadratic type and its determinant walls

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

The first and third displayed factors are already nonzero by (22).  Name the
other two factors

```text
A=384lambda h^2-560h+693,
C=114048lambda h^2-25088h^2-44352h+68607.              (26)
```

Both exceptional walls are empty.

On `A=0`, substitute

```text
lambda=(560h-693)/(384h^2).                             (27)
```

The Groebner basis of all nine specialized remainder equations contains
`h^2`.  This contradicts `h!=0`.

On `C=0`, substitute

```text
lambda=(25088h^2+44352h-68607)/(114048h^2).             (28)
```

The two top-next equations become parallel.  Up to a nonzero rational unit,
the gcd of their two augmented consistency minors is

```text
h^6(8h-99)^3(32h-99)(56h-99)^2(64h-99)^4.              (29)
```

Since `h!=0`, a solution can occur only at

```text
h in {99/8,99/32,99/56,99/64}.                         (30)
```

At each of these four exact rational values, the complete specialized ideal
in `Q[b,d,e]` is the unit ideal.  Hence `C=0` is empty too.

It remains to work on `AC!=0`, where (25) solves `b,d` uniquely.  The next
two equations are linear in `e`.  Their exact resultant factors as

```text
h^3 A^2(384lambda h^2-280h+231)^3 C^3 K(lambda,h),      (31)
```

where the residual primitive polynomial `K` has bidegree

```text
(deg_lambda K,deg_h K)=(7,15)
```

and exactly 69 terms.  All displayed factors except `K` are inverted on
this chart.  Therefore every surviving `Q3` factor must lie over `K=0`.
This is a necessary reduction, not a claim that every point of `K=0` lifts
to a factor.

For continuing the residual calculation, the localized root-scaled chart

```text
Lambda=lambda h^2,       beta=b/h,
delta=d/h^2,             epsilon=e/h^2                 (32)
```

puts all nine saturated coefficient equations in total degree at most eight.

This was the determinant-chart stopping boundary:

```text
R_(lambda,mu) reducible
  iff a Q3 factor satisfying the nine equations exists;
every such factor satisfies A!=0, C!=0, K(lambda,h)=0.  (33)
```

The pair-field argument below bypasses the residual `K`-curve and proves
that the nine-equation scheme is empty.

## 6. Fixed-section pair field and Hensel obstruction

Write

```text
R(p,v)=r_0(v)+r_1(v)p+...+r_5(v)p^5,
r_0=-567L_5.                                             (H1)
```

Suppose `R` has a quadratic factor over `C`.  Because `mu!=0`, the top form
(14) has no `v` factor.  The quadratic side can therefore be scaled uniquely
to be monic in `p`:

```text
Q=p^2+q_1(v)p+q_0(v),       deg(q_1)<=1, deg(q_0)<=2.    (H2)
```

Its `v^2` coefficient at `p=0` is nonzero.  Hence

```text
q_0=c qhat,       c!=0,       qhat=v^2+alpha v+beta.     (H3)
```

Equation `r_0=q_0s_0` says exactly that `qhat` divides `L_5`.  Dividing
`L_5` by `qhat` and taking the two remainder coefficients defines the
unordered-pair scheme.  Its exact lexicographic basis has the triangular
form

```text
constant*beta+B_9(alpha),             P_10(alpha),       (H4)
```

where

```text
P_10(alpha)
=290627997810865923974832 alpha^10
 -24018842794286439997920 alpha^9
 -8783750360720454296760 alpha^8
 +471786003250475288760 alpha^7
 +105739913176342048395 alpha^6
 -2189869887890543364 alpha^5
 -556279469540931180 alpha^4
 -5665718237324640 alpha^3
 +836795172789840 alpha^2
 +29407080274880 alpha
 +321069284032.
```

The coefficient of `beta` in (H4) is the nonzero integer
`2466637614990092544000`.  Exact factorization proves that `P_10` is
irreducible over `Q`; it is also squarefree.  Consequently

```text
A=Q[alpha]/(P_10)                                       (H5)
```

is a reduced degree-ten field, `beta` is a fixed element of `A`, and its ten
geometric embeddings are precisely the ten unordered pairs of the five
distinct roots of `L_5`.  There are no missing components: a monic quadratic
divisor of a squarefree quintic is exactly the product attached to one such
pair.  The complementary quotient `L_5/qhat` has degree three, so exchanging
the names “factor” and “cofactor” in a `2+3` factorization does not escape
this atlas; one simply selects its quadratic side.

Now retain the scale rather than normalizing it away.  Put

```text
x=1/c,       q_1=cg,       g=g_a v+g_b,
shat=L_5/qhat.                                           (H6)
```

Coefficient comparison in `p` gives the following successive divisibility
conditions in `A[v]`:

```text
r_1+567g shat                         =qhat t_1,
r_2-g t_1+567x shat                   =qhat t_2,
r_3-g t_2-x t_1                       =qhat t_3.          (H7)
```

The two remainder coordinates of the first line form a nonsingular `2x2`
system for `g_a,g_b`.  Since `r_2` is affine in `lambda`, those of the second
line form a nonsingular `2x2` system for `x,lambda`.  The reduced numerators
of both determinants and of the forced `x` have degree nine, ten terms, and
are coprime to `P_10`; they are units at every root pair, and `x!=0`.

This scaling is checked directly, not inferred from the remainders.  With

```text
q_0=c qhat, q_1=cg,
s_0=-567x shat, s_1=xt_1, s_2=xt_2,                    (H8)
```

the companion multiplies `Q` by its provisional cofactor and recovers
`r_0,r_1,r_2` exactly.

The coefficient `r_3` is affine in `mu`.  One remainder coordinate in the
third line of (H7) determines `mu`; its pivot is again a degree-nine unit.
After that substitution the other coordinate is, up to a nonzero rational
scalar, the following primitive polynomial:

```text
O_9(alpha)
=1130885334253871586766989948754019527574252106308468112 alpha^9
 -175171904604809806522915642021425303092320752846661248 alpha^8
 -23987998932486549476841913451129298773242820286819048 alpha^7
 +4155941615827072747203528754630926010782491728169112 alpha^6
 +180146909639856808193786189847345565305728206396137 alpha^5
 -34078508926274194253641821924511345688894165032362 alpha^4
 -563188439587546864344955993805021039226394561752 alpha^3
 +99164223051444387326517957329206207387270288048 alpha^2
 +305846347722090566604280525587712300295540688 alpha
 -70610888837110659252936402830409827008993312.          (H9)
```

Exact Euclidean division gives

```text
gcd(P_10,O_9)=1.                                        (H10)
```

Thus `O_9` is nonzero at every one of the ten root pairs.  The third Hensel
condition can never hold, so `R` has no quadratic factor for any physical
`lambda,mu`.  Sections 3--4 already exclude every line factor.  Since a
reducible total-degree-five polynomial has a factor of degree at most two,
every physical `R_(lambda,mu)` is absolutely irreducible.

## 7. The square lift closes the stratum

By Section 6 every physical fibre `R_(lambda,mu)=0` is absolutely
irreducible.  Let `C_(lambda,mu)` be its smooth projective normalization.

Restore the physical square class

```text
Y=y/sqrt(B),                    Y^2=1/p.                (34)
```

At the five points (11), `p` has valuation one.  Hence `1/p` is not a
square, the double cover is connected, and it ramifies at all five points.
The branch divisor of a quadratic extension has even degree, so it has at
least six places.  Riemann--Hurwitz gives

```text
g(double cover)>=-1+6/2=2.                              (35)
```

A rational Keller trajectory would give a map from `P^1` to this positive-
genus normalization, hence a constant map.  Then `p,v,Y,y,u` are constant;
the open first flux reconstructs `zeta`, and then `Z,T,q` are constant,
contradicting the genuine deck involution.  Thus emptiness of the one ideal
in (24), now proved by the pair-field obstruction, closes all of (1).

## 8. Finite hostile probe

As a hostile control only, the companion takes the 18 distinct nonzero
rationals

```text
{n/d : d in {1,2,3}, -4<=n<=4, n!=0}                  (36)
```

for each of `lambda,mu`, factors all `324` specialized bivariate
eliminants exactly over `Q`, and finds

```text
reducible fibres: 0.                                    (37)
```

This is **FINITE-EXACT** evidence.  It is not used in the exact exclusions
or in (33); by itself it does not exclude algebraic parameter values and is
not a proof of uniform irreducibility.  Section 6 supplies that proof.

## 9. Reproduction and scope

Run

```bash
python3 04-computation/jc2_degree22_bdw_pair_hensel_closure_thm2617.py
python3 -O 04-computation/jc2_degree22_bdw_pair_hensel_closure_thm2617.py
python3 04-computation/jc2_degree22_bdw_triple_last_quadratic_thm2617.py
python3 -O 04-computation/jc2_degree22_bdw_triple_last_quadratic_thm2617.py
```

The closure companion independently reconstructs the eliminant and fixed
section, proves the ten-point pair algebra is one reduced degree-ten field,
performs the three Hensel steps in that field, checks every division unit,
reconstructs the product through `p^2`, and proves (H10).  The reduction
companion reconstructs (7)--(14), the exhaustive top inventory, all four
top-type exclusions, the last quadratic chart and determinant, both
exceptional-wall eliminations, the generic compatibility divisor, the
localized degree-eight chart, the square-lift genus invoice, the `y=0`
boundary, and the frozen hostile universe.  All truth-bearing checks use
`require` and remain active under optimized Python.

An independent hostile audit rederived the coefficient recurrences (H7)
directly from a monic quadratic times cubic product, checked that `mu!=0`
makes every top line nonvertical and hence makes the line/quadratic inventory
exhaustive, and verified that exchanging factor and cofactor still selects
the unique quadratic side.  It separately checked the reduced ten-point
pair scheme, every field division, the direct product identities through
`p^2`, and the coprimality of `O_9` with `P_10`.  Normal and optimized runs
of both exact companions byte-match their stored transcripts and hashes.
Finally, the audit checked that the five simple `p=0` points are distinct
odd branch places, so connectedness and the genus lower bound do not rely on
a singular model or on counting multiplicity as separate places.

Uniform irreducibility and the `B,D,W` closure are **PROVED**.  The residual
`K(lambda,h)=0` chart is now a superseded intermediate description of an
empty factor scheme, not an open dependency.

Nothing here closes another support-three stratum, a branch outside the
inherited nonsplit exact-square-prefix reduction, split/even short edges, or
integral `2`-adic order raising.  It proves neither `JC(2)` nor `DC(2)`.

**QED.**
