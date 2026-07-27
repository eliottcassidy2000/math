---
id: THM-2518
title: "Perron inverse-branch owner-word cospan recovery"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Perron powers
  disintegrate every entrywise moment exactly over
  ordered inverse-branch tuples of T^L, or equivalently over the deck-arrow
  components of the fibre product T x_(T^L) T.  For rational BV response
  densities the owner-gated m-th Perron moment converges uniformly to rho A^m
  at an explicit O(13^-L) rate.  Applied to the THM-2513 first/square vector,
  one genuine nonidentity cospan carries all 72 mixed colours and all 5,184
  vector cut coefficients.  More strongly, for every prescribed ordered pair
  of distinct nonzero first-collision residues there is a same-owner
  recurrence-supported triangle: a Boolean old owner at x, two distinct
  response predecessors, and one common future owner--word event at T^L x.
  A distinguished response arm retains the old deep diagonal-zero lift, and
  expansion before the deep sum is literally Boolean.  The selected higher
  addresses may differ by residue, are not covariant without the role-wise
  ancestry carries, and need not preserve old source/deep phase labels.
  Square self-cospans C_d are antipodally equal, and recurrence support does
  not imply nonzero owner-loop drift, a present-time owner/arrival
  identification, a row exclusion, or LRC(14).
source: codex-2026-07-27-perron-inverse-branch-cospan
depends_on:
  - THM-2449-coprime-owner-anova-and-delta-replica-boundary
  - THM-2478-delayed-owner-handoff-graft-and-deep-sheet-rebase-boundary
  - THM-2508-affine-cut-bundle-covariance-and-carry-permutation
  - THM-2513-anchored-first-or-second-moment-spectrum-and-pair-space-boundary
related:
  - THM-2299-rooted-current-service-energy-and-base-phase-no-go
  - THM-2471-owner-first-collision-weighted-root-service-and-temporal-atom-boundary
  - THM-2515-haar-self-correlation-disintegration-and-rational-shift-recovery
  - THM-2516-anchored-moment-simplex-disintegration-and-owner-star-recovery
  - THM-2517-cubic-anchored-spectrum-flt3-and-three-time-boolean-lift
script: 04-computation/lrc14_perron_collision_cospan_thm2518_2519_referee.py
output: 05-knowledge/results/lrc14_perron_collision_cospan_thm2518_2519_referee.out
script_sha256: b09f9e8940914729583057cf5c92833a0f8dd85673b47bf1b2d5d71afe895e39
output_sha256: 2818c04f4ef56f352e0365504d7b2f8bf194bf8bdb0b04e4268a5340acba7898
hash_basis: working-tree bytes (LF)
---

# THM-2518 -- Perron inverse branches turn the moment simplex into a lawful cospan

**PROVED.**

THM-2513 proves that the first/square response vector has no primitive mixed
zero.  THM-2516 puts that vector on an owner-rooted edge--triangle, but its
free translation arms need not be lawful ancestry edges.  The correct
replacement is already intrinsic to the expanding map

```text
T(x)=13x mod 1.
```

Two points are related at depth `L` when they have the same image under
`T^L`.  Thus the natural binary object is not a cosmetic tournament on
runners.  It is the finite equivalence-relation groupoid

```text
R_L={(x,x'):T^Lx=T^Lx'}
   =T x_(T^L) T.                                                (1)
```

Its arrow components are the deck needles

```text
x' = x+d/13^L,                  d in Z/13^L Z.                  (2)
```

The identity arrow is `d=0`; a unit `d` joins two histories which first
coalesce exactly at depth `L`.  A positive owner--word event evaluated at
their common future point is automatically the same Boolean factor on both
ends.  This is the missing structure behind the free shifts.

## 1. Exact finite fibre-product identity

Let `Y` be a nonempty finite set and let

```text
X=Y x C_N,                         pi(y,r)=y.                    (3)
```

For a finite label set `J`, functions `F_j:X->Q`, and `G:Y->Q`, put

```text
(P_N F_j)(y)=1/N sum_(r in C_N)F_j(y,r).                        (4)
```

For every integer `m>=1`,

```text
1/|Y| sum_y G(y)(P_NF_j(y))^m

 =1/N^m sum_(r_1,...,r_m)
    1/|Y| sum_y G(y) product_(q=1)^m F_j(y,r_q).                (5)
```

This is just expansion of the product, but its coordinate meaning matters:
all `m` arms in every summand lie in one fibre of `pi` and retain the same
future factor `G(y)`.

There is an exact cyclic-grid model.  Take

```text
Y=C_D,                     X=C_(DN),

pi(k mod DN)=k mod D.                                           (6)
```

The inverse branches over `a in C_D` are `a+Dr`, `r in C_N`.  Hence (5)
is an identity of rational finite tables, with no limiting or genericity
qualification.

For `m=2` and `N>=2`, deleting the diagonal branch pairs gives the pointwise
identity

```text
1/[N(N-1)] sum_(r!=s)F_j(y,r)F_j(y,s)

 =N/(N-1)(P_NF_j(y))^2
  -1/(N-1)P_N(F_j^2)(y).                                      (7)
```

Thus distinct arms have the same independent-pair limit whenever the Perron
averages converge to the global mean.  The diagonal is not needed.

## 2. Circle Perron identities and the sheet-summed needle

Let `N>=1` and define the Perron operator of `x->Nx mod 1` by

```text
(P_N f)(y)=1/N sum_(r=0)^(N-1) f((y+r)/N).                     (8)
```

Let `F_j:T->R` be bounded measurable functions and let `G:T->R` be
integrable.  Write

```text
A_j=integral_T F_j,
rho=integral_T G.                                               (9)
```

The exact circle form of (5) is

```text
M_N^[m](j)
 =integral_T G(y)(P_NF_j(y))^m dy

 =1/N^m sum_(r_1,...,r_m)
    integral_T G(y) product_q F_j((y+r_q)/N)dy.                (10)
```

The terms on the second line are ordered inverse-branch charts.  An
absolute branch label depends on the chosen cut of the covering.  Summing
the simultaneous branch index first retains only deck differences and gives
a global circle object.  For

```text
d=(d_2,...,d_m) in C_N^(m-1),
```

put

```text
L^N_d(j)
 =integral_T G(Nx)F_j(x)
    product_(q=2)^m F_j(x+d_q/N) dx.                            (11)
```

Because `N(x+d_q/N)=Nx` modulo one, every arm in (11) has exactly the same
future point and the same factor `G(Nx)`.  Splitting the `x`-circle into its
`N` inverse-branch cells proves

```text
1/N^(m-1) sum_d L^N_d(j)=M_N^[m](j).                           (12)
```

At degrees one and two, use the notation

```text
S^N(j)=integral_T G(Nx)F_j(x)dx,

C^N_d(j)=integral_T G(Nx)F_j(x)F_j(x+d/N)dx.                  (13)
```

Then

```text
S^N(j)=integral_T G P_NF_j,

1/N sum_d C^N_d(j)=integral_T G(P_NF_j)^2.                    (14)
```

This is also the precise bridge to THM-2515.  Let

```text
I_r=[r/N,(r+1)/N),
x_r(y)=(y+r)/N,

S_r(j)=integral_T G(y)F_j(x_r(y))dy,

T_(r,s)(j)=integral_T G(y)F_j(x_r(y))F_j(x_s(y))dy.            (15)
```

For `d=s-r mod N`, change variables on `I_r` to obtain

```text
S_r(j)=N integral_(I_r)G(Nx)F_j(x)dx,

T_(r,s)(j)
 =N integral_(I_r)G(Nx)F_j(x)F_j(x+d/N)dx.                    (16)
```

Consequently

```text
S^N=1/N sum_r S_r,

C^N_d=1/N sum_r T_(r,r+d).                                   (17)
```

Equations (16)--(17) retain the exact ordered-pair refinement requested by
the ancestry stalk, while `C_d` is the sheet-summed, gauge-safer needle.
Only the difference survives a simultaneous reindexing of the covering.

There is an unavoidable antipodal symmetry:

```text
C^N_(-d)=C^N_d.                                                (18)
```

Indeed, translate `x` by `-d/N` and use periodicity of `G(Nx)`.  Hence this
square cospan is intrinsically unoriented.

## 3. Quantitative BV recovery of independent moments

For a periodic BV function `f`, the shifted `N`-point rectangle rule gives

```text
||P_Nf-integral_T f||_infinity<=Var(f)/N.                      (19)
```

To prove (19), use one sample `(y+r)/N` in each interval
`[r/N,(r+1)/N)`.  The error on that interval is at most its oscillation
divided by `N`; the sum of the interval oscillations is at most `Var(f)`.

Assume now that

```text
0<=G<=1,                     rho>0,
0<=F_j<=M_j,                 V_j=Var(F_j).                       (20)
```

For every `m>=1`, equations (19)--(20) give

```text
|M_N^[m](j)-rho A_j^m|
 <=rho m M_j^(m-1)V_j/N.                                     (21)
```

In particular,

```text
m=1:  error<=rho V_j/N,
m=2:  error<=2rho M_jV_j/N,
m=3:  error<=3rho M_j^2V_j/N.                                 (22)
```

The last estimate is the quantitative cubic interface reserved for
THM-2517.  It uses only
`|z^m-a^m|<=mM^(m-1)|z-a|` on `[0,M]`.

The ordered distinct-pair average in (7) also has an explicit bound.  Put

```text
Q^off_N(j)
 =1/[N(N-1)] sum_(r!=s)
   integral_T G(y)F_j(x_r(y))F_j(x_s(y))dy.                    (23)
```

For `N>=2`, (7), (19), and `0<=P_N(F_j^2),A_j^2<=M_j^2` imply

```text
|Q^off_N(j)-rho A_j^2|
 <=rho(2M_jV_j+M_j^2)/(N-1).                                 (24)
```

Thus one ordered pair with `r!=s` suffices once a limiting square
coefficient is nonzero.  This pair-level conclusion is exact in a fixed
branch chart; it is not yet a covariant choice of an absolute ancestry
sheet.

For a global nonidentity arrow, one can prescribe its first-collision
residue.  From now on take

```text
N=13^L,                       L>=1,
M=N/13,

D_u={d in C_N:d=u mod 13},    u in F_13^*.                    (25)
```

For fixed `x`, the points `x+d/N`, `d in D_u`, form a translated
`M`-point grid.  Therefore

```text
B^N_(u,j)(x)=1/M sum_(d in D_u)F_j(x+d/N),

||B^N_(u,j)-A_j||_infinity<=13V_j/N.                          (26)
```

Average (13) over `D_u`.  Since

```text
integral_T G(Nx)F_j(x)dx=S^N(j),
```

equations (19) and (26) give, uniformly in `u!=0`,

```text
|1/M sum_(d in D_u)C^N_d(j)-rho A_j^2|
 <=14rho M_jV_j/N.                                            (27)
```

The constant `14` consists of the `13`-grid error for the second arm and
the one-grid mixing error for the first.  Hence every prescribed nonzero
collision residue, not merely some distinct branch pair, recovers the square
moment asymptotically.

For completeness, the cubic version of (12) is

```text
L^N_(d,e)(j)
 =integral_T G(Nx)F_j(x)F_j(x+d/N)F_j(x+e/N)dx,

1/N^2 sum_(d,e)L^N_(d,e)=M_N^[3].                             (28)
```

Only `3N-2` of the `N^2` difference pairs have a repeated arm
(`d=0`, `e=0`, or `d=e`).  Each normalized mixed coefficient of one such
table has absolute value at most `rho(max_j M_j)^3`.  Thus the total repeated
contribution to the normalized average is at most

```text
(3N-2)/N^2 rho(max_j M_j)^3.                                  (29)
```

Combined with the cubic bound in (22), a nonzero mixed coefficient of
`A^(circ 3)` survives on one triple of pairwise distinct inverse branches
for every sufficiently large `N`.

If `F_j` and `G` are rational step functions with rational endpoints, then
every table in (10)--(18), (23), and (28) is rational: affine pullback keeps
rational endpoints, finite products keep rational values, and integration
of the common rational refinement is rational.  This makes the Galois
propagation below exact.

## 4. The first/square vector on every collision residue

Now let

```text
J=F_7 x F_13
```

and let the rational nonnegative table `A` satisfy the THM-2513 hypotheses:

```text
A_(ell,0)=a 1_(ell=0),                 a>0,

s -> A_(0,s) is nonconstant.                                  (30)
```

Suppose `A_j=integral F_j` with rational nonnegative BV step densities.
For a mixed character define

```text
Vhat_A(kappa,b)
 =(Ahat(kappa,b),(A^(circ 2))hat(kappa,b)).                    (31)
```

THM-2513 says

```text
Vhat_A(kappa,b)!=(0,0)
               for every kappa in F_7^*, b in F_13^*.         (32)
```

Let `G` be any positive rational Boolean future event.  For example, the
THM-2478 owner--word block is

```text
G(y)=1_(E_j)(y)1_Q(T^K y),                 rho=integral G>0.   (33)
```

Fix one mixed pair `(kappa_0,b_0)` and put

```text
delta=||Vhat_A(kappa_0,b_0)||_2>0,

E_1(N)=rho/(91N) sum_j V_j,

E_2(N)=14rho/(91N) sum_j M_jV_j.                              (34)
```

Choose `L` in either prescribed parity class so large that

```text
sqrt(E_1(N)^2+E_2(N)^2)<rho delta/2,       N=13^L.             (35)
```

For every `u in F_13^*`, equations (22) and (27) show that

```text
1/|D_u| sum_(d in D_u)
  ( Shat^N(kappa_0,b_0), Chat^N_d(kappa_0,b_0) )

```

has Euclidean norm at least `rho delta/2`.  Hence, for every prescribed
`u`, there is some

```text
d_u in D_u                                                    (36)
```

such that

```text
(Shat^N(kappa_0,b_0),Chat^N_(d_u)(kappa_0,b_0))!=(0,0),

||(Shat^N,Chat^N_(d_u))||_2>=rho delta/2.                     (37)
```

Equation (17) gives a still finer fixed-chart statement: the vector in
(37) is the average over `r` of

```text
(Shat_r(kappa_0,b_0),
 That_(r,r+d_u)(kappa_0,b_0)).                                 (37a)
```

Thus one ordered pair `(r,r+d_u)`, with distinct branches, has norm at
least the right side of (37).  Its two tables are rational, so the same
ordered pair Galois-saturates all mixed colours.  The sheet-summed `d_u`
statement remains the canonical one: the absolute `r` in (37a) depends on
the chosen branch chart.

Every selected table is rational.  Galois conjugation over
`Q(zeta_91)` therefore turns (37) into

```text
(Shat^N(kappa,b),Chat^N_(d_u)(kappa,b))!=(0,0)

for all kappa,b!=0                                             (38)
```

on that same cospan.

ANOVA-centre the two tables, transpose their interactions to row-zero
`F_13 x F_7` arrays, and apply THM-2508 coordinatewise.  Its exact formula is

```text
(Psi^S_(tau,a_0)(alpha,beta),
 Psi^C_(tau,a_0)(alpha,beta))

 =91K_(alpha tau,beta)
   (Shat^N(beta a_0,-alpha),
    Chat^N_(d_u)(beta a_0,-alpha)).                            (39)
```

The geometric factor is nonzero whenever

```text
tau,alpha in F_13^*,             a_0,beta in F_7^*.
```

Thus all

```text
12*12*6*6=5,184                                               (40)
```

primitive vector cut coefficients survive on the selected nonidentity
cospan.

Equation (36) is an **all-residue availability atlas**.  The twelve
addresses `d_u` may have different higher digits.  It is not one common
`u`-indexed Boolean stalk, and it does not assert a nonzero DFT across the
collision-residue coordinate.

## 5. Stronger same-owner recurrence triangle

The abstract future cospan can retain the old Boolean owner as well.  At the
fixed THM-2449 clock, take its genuine pre-deep-sum owner event

```text
H(x)=H_own(x)
 =U_(0,0)(x)d_(j,0)(x)Q^epsilon_0(Rx),

h=integral_T H>0.                                             (41)
```

Choose in (33) the future owner--word event for the **same** owner label
`j`.  Let `u,v in F_13^*` be prescribed and distinct.  For

```text
d in D_u,                       e in D_v,
```

define the edge--triangle response vector

```text
mathcal S_d(j_0)
 =integral_T H(x)G(Nx)F_(j_0)(x+d/N)dx,

mathcal T_(d,e)(j_0)
 =integral_T H(x)G(Nx)
    F_(j_0)(x+d/N)F_(j_0)(x+e/N)dx.                            (42)
```

Here `j_0=(ell,s)` is a response-table index, not the fixed owner label in
`E_j`.

Let

```text
w_N=integral_T H(x)G(Nx)dx,

epsilon_N=Var(H)Var(G)/(12N).                                 (43)
```

The two-BV covariance estimate of THM-2478 gives

```text
|w_N-hrho|<=epsilon_N.                                        (44)
```

Average (42) over the higher digits of `d` and `e`.  With the grids from
(26),

```text
1/|D_u| sum_d mathcal S_d(j_0)
 =integral H G(Nx)B^N_(u,j_0),

1/(|D_u||D_v|) sum_(d,e)mathcal T_(d,e)(j_0)
 =integral H G(Nx)B^N_(u,j_0)B^N_(v,j_0).                     (45)
```

Since `H,G<=1`, equations (26) and (44) yield the audited bounds

```text
|average_d mathcal S_d(j_0)-hrho A_(j_0)|
 <=13V_(j_0)/N+M_(j_0)epsilon_N,

|average_(d,e) mathcal T_(d,e)(j_0)-hrho A_(j_0)^2|
 <=26M_(j_0)V_(j_0)/N+M_(j_0)^2epsilon_N.                    (46)
```

For the mixed pair fixed in Section 4, put

```text
E^H_1(N)
 =1/91 sum_(j_0)(13V_(j_0)/N+M_(j_0)epsilon_N),

E^H_2(N)
 =1/91 sum_(j_0)(26M_(j_0)V_(j_0)/N
                  +M_(j_0)^2epsilon_N).                       (47)
```

Take `L` in the desired parity class so large that

```text
sqrt((E^H_1(N))^2+(E^H_2(N))^2)<hrho delta/2.                 (48)
```

The vector average of

```text
(mathcal Shat_d(kappa_0,b_0),
 mathcal That_(d,e)(kappa_0,b_0))
```

over `D_u x D_v` is then nonzero and has norm at least
`hrho delta/2`.  Consequently one ordered pair `(d,e)` satisfies

```text
||(mathcal Shat_d,
   mathcal That_(d,e))(kappa_0,b_0)||_2
 >=hrho delta/2.                                               (49)
```

Rational Galois propagation and (39) give all `72` mixed vector colours and
all `5,184` primitive vector cuts on this same triangle.

The geometry is exact.  Since

```text
d=u mod 13,              e=v mod 13,

u,v,u-v!=0,                                                   (50)
```

the three points

```text
x,                       x+d/N,                  x+e/N          (51)
```

are pairwise distinct and first coalesce under `T` exactly at depth `L`.
Their penultimate relative residues are `0,u,v`.  Moreover,

```text
T^L(x+d/N)=T^L(x+e/N)=T^Lx,                                  (52)
```

so every entry in (42) has a representation containing:

```text
the old Boolean owner event H at x;
two distinct full response densities on sibling histories,
  expanded into lawful packet atoms in Section 6;
the same future owner E_j at T^Lx;
the same future terminal word Q at T^(L+K)x.                   (53)
```

This is a genuine **same-owner-recurrence-supported multiprehistory
signal**.  It replaces THM-2516's free arms by one owner-labelled recurrence
triangle with prescribed distinct first-collision colours.

## 6. Exact old-deep lift and Boolean expansion

Write the lawful THM-2449 response density before its old deep-label sum as

```text
H_ell(q,s,t)=integral_T h_(ell,q,s,t)(x)dx,

F_(ell,s)(x)=sum_q h_(ell,q,s,0)(x),                           (54)
```

where the retained moving safe factor gives the pointwise diagonal zero

```text
h_(ell,t,s,t)(x)=0                  almost everywhere.         (55)
```

For the selected recurrence edge and triangle of Section 5, define

```text
mathcal H^S_ell(q,s,t)
 =integral_T H(x)G(Nx)
    h_(ell,q,s,t)(x+d/N)dx,

mathcal H^T_ell(q,s,t)
 =integral_T H(x)G(Nx)
    h_(ell,q,s,t)(x+d/N)F_(ell,s)(x+e/N)dx.                    (56)
```

Then exactly

```text
sum_q mathcal H^S_ell(q,s,0)=mathcal S_d(ell,s),

sum_q mathcal H^T_ell(q,s,0)=mathcal T_(d,e)(ell,s),

mathcal H^S_ell(t,s,t)=mathcal H^T_ell(t,s,t)=0.               (57)
```

Therefore the THM-2449 deep-frequency proof applies to each surviving
channel without approximation: a nonzero mixed coefficient forces some
nonzero old-deep character on the distinguished `d`-arm.  The deep index
`q in F_13` and the depth-`L` branch addresses `d,e in C_N` are different
typed coordinates and are not identified here.

Finally expand the last factor in (56):

```text
F_(ell,s)(x+e/N)
 =sum_(q')h_(ell,q',s,0)(x+e/N).                              (58)
```

Each resulting summand is a product of literal Boolean packet factors on
the inverse-branch fibre product.  The live summed density `F` may take the
values one and two; Booleanity is asserted only after the expansion (58),
not for `F` itself.

## 7. Exact sheet obstruction and stopping boundary

The cospan is lawful as a product of complete packet events evaluated at
their own physical points.  That does **not** make its selected address
sheet-free.  For the standard LRC danger indicator, let a phase bank be

```text
{d(cx-q/p):q in C_p}
```

and, for `d!=0`, write

```text
delta_x=d/N,               v=nu_13(d),       0<=v<L,

k_col=L-v.                                                        (59)
```

Translation by `delta_x` preserves this bank as a set, by a permutation of
its labels, if and only if

```text
p c d/N in Z.                                                  (60)
```

For the septimal source bank, (60) is equivalent to

```text
k_col<=nu_13(c_j).                                             (61)
```

For the old deep thirteen-bank it is equivalent to

```text
k_col<=nu_13(c_3)+1.                                          (62)
```

A unit needle has `k_col=L`.  Perron recovery requires arbitrarily large
`L`, while the right sides of (61)--(62) are fixed.  Divisibility can restore
the old phase banks only by moving the collision back to a bounded valuation
horizon.  Thus the missing carry/sheet coordinate is an exact incompatibility,
not merely an omitted continuity argument.

The full Perron aggregate is covariant, and the deck difference is invariant
under simultaneous scalar reindexing.  A selected absolute branch pair or
difference is nevertheless not covariant under the full THM-2365 factor-wise
target action: distinct packet roles acquire distinct depth-`L` carries.
The per-role ancestry residues of THM-2478 remain necessary.  In particular:

- the twelve residue witnesses in Section 4 can use different higher
  addresses;
- source parity can be fixed by taking the cofinal even or odd subsequence,
  but parity does not identify the other role-wise carries;
- the old deep label in (56) is retained, not rebased to the future owner;
- source-time and arrival-time atom projections are still different; and
- the common future owner--word factor does not turn two alternative
  histories into one present-time Boolean current.

There is also a sharp charge boundary.  Equation (18) already prevents the
square cospan from orienting `d` against `-d`.  More basically, recurrence
support does not control a signed owner-loop drift.  On two recurrence atoms
with canonical owner contributions `+1,-1`, the aggregate owner drift is
zero, while a Boolean response table supported only on the first atom may be
the single delta entry at `(ell,s)=(0,0)` and hence has nonzero mixed spectrum
and cut bundle.  A coupling identity to the signed drift, not another common
support factor, is required.

Thus THM-2518 closes the following former seam:

```text
first/square spectrum
  -> arbitrary owner-supported translations
```

by the stronger lawful object

```text
first/square spectrum
  -> old-owner / two-prehistory / common-future-owner-word cospan
     with prescribed distinct first-collision residues
     and an exact old-deep lift.                                (63)
```

It does not prove nonzero owner-loop charge, a role-wise sheet rebase, a
same-time owner/arrival equality, a scalar-row exclusion, or LRC(14).
**QED.**

## 8. Exact finite-tower referee

Run

```bash
python3 04-computation/lrc14_perron_collision_cospan_thm2518_2519_referee.py
python3 -O 04-computation/lrc14_perron_collision_cospan_thm2518_2519_referee.py
```

Both runs reproduce the stored transcript byte-for-byte.  The independent
`Fraction` model checks the finite Perron and groupoid identities, ancestry-
sheet/needle disintegration, every fixed last digit, exact first-collision
timing, antipodal symmetry, the anchored recurrence edge and triangle
averages, and both positive and hostile controls.  A separate line-by-line
audit rederived every BV constant, selection inequality, Galois/cut factor,
deep lift, and valuation boundary used above.
