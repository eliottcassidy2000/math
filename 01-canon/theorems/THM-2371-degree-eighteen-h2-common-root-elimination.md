---
id: THM-2371
title: "Degree-eighteen H2 common-root elimination and unconditional cube descent"
status: >
  PROVED + VERIFIED-EXACT. Every point of THM-2357's degree-ten
  H_2 S_4^2 coefficient locus satisfies q_5(1)!=0 and
  Res_y(p_3,q_5)!=0. The singular-order wall and all three exact
  resultant components are eliminated by forced-root quotients and
  univariate subresultant certificates. Consequently THM-2360's
  Laurent-UFD factorization B=c(s-r_0)U_3^3 is unconditional on the
  whole H_2 stratum. This closes the two auxiliary-hypothesis gaps, not
  the H_2 stratum itself: the final coprime linear-times-cube coefficient
  locus remains open. The H_4 stratum, JC(2), and DC(2) also remain open.
source: codex-2026-07-25-degree-eighteen-h2-common-root
depends_on:
  - THM-2332-degree-eighteen-genus-zero-square-class-and-dessin-trap
  - THM-2338-degree-eighteen-deep-common-root-wall-hurwitz-quartet
  - THM-2345-degree-eighteen-common-root-wall-saturation
  - THM-2357-degree-eighteen-h2-moving-root-reduction
  - THM-2360-degree-eighteen-quadratic-ring-cube-descent
related:
  - THM-2347-degree-eighteen-double-zero-wall-saturation
  - THM-2359-degree-eighteen-perfect-quartic-wall-closure
script: 04-computation/jc2_degree18_h2_common_root_elimination_thm2371.py
output: 05-knowledge/results/jc2_degree18_h2_common_root_elimination_thm2371.out
script_sha256: 53e224d562d9ea843ea9f9bae8a8fa95f981c417a3770cef5c9a9dfd656cb6a3
output_sha256: 51b55ac02ffdb539941c09d296e9500cf8a157015cfd9aa428158094ee3e2863
hash_basis: working-tree bytes (LF)
---

# THM-2371 -- eliminate the H2 common-root residual

**PROVED + VERIFIED-EXACT.**

THM-2357 moves the unique three-cycle value to `y=1` and reduces the
mixed degree-eighteen branch to

```text
R_10=4(y-1)p_3^3+49q_5^2=H_2 S_4^2,              (1)

H_2 squarefree,                  deg(H_2,S_4)=(2,4).
```

THM-2360 proves a rigid Laurent cube factorization under

```text
q_5(1)!=0,                       Res_y(p_3,q_5)!=0. (2)
```

The two conditions in (2) were genuine gaps: THM-2360 gives an exact
nonmaximal-order hostile showing that neither follows from an abstract
`H_2 S_4^2` identity. This theorem uses the special degree-eighteen
coefficient pattern to prove both conditions. It therefore makes the
cube descent unconditional on the `H_2` stratum.

It does **not** show that the resulting cube coefficient equations are
empty.

## 1. The moving-root coefficient atlas

On THM-2357's incidence surface, direct division gives

```text
p_3
 =35(y+1)(54B+7y^2+7),                            (3)

q_5
 =7y[
    1620By^2+1620By+2430B+26244C(y+1)
    +77(y^4+y^3+y^2+y)+182
   ].                                              (4)
```

In particular,

```text
lc(p_3)=245,             lc(q_5)=539,

lc(R_10)=73060029.                                 (5)
```

All three leading coefficients are constant and nonzero. Hence every
subresultant used below specializes without a degree-drop exception.

The distinguished evaluations are

```text
p_3(0)=35(54B+7),        q_5(0)=0,

p_3(-1)=0,               q_5(-1)=-14(1215B+91),

p_3(1)=140(27B+7),       q_5(1)=14K,               (6)

K:=245+2835B+26244C.
```

The weighted action before setting the three-cycle value to one sends

```text
p_3(Y) -> lambda^3 p_3(Y/lambda),

q_5(Y) -> lambda^5 q_5(Y/lambda).                  (7)
```

Thus `q_5(r)` scales by `lambda^5`, while the `(3,5)` resultant scales
by `lambda^15`. Their vanishing is invariant under THM-2357's
normalization. It is therefore enough to work in (3)--(4).

We repeatedly use the elementary multiplicity fact

```text
G=H_d T_m^2,                     H_d squarefree

implies T_m divides gcd(G,G').                         (8)
```

If `alpha` is a root of `G` of order at least two, then
`alpha` is a root of `T_m`, even if it is also a root of `H_d`.
Consequently a forced double root can be divided from both sides of
the square-class identity without a coprimality assumption.

## 2. The singular-order wall is empty

The first condition in (2) fails exactly on

```text
K=0,

C=-(245+2835B)/26244.                              (9)
```

On this line,

```text
R_10=343(y-1)T_9,                                  (10)
```

where

```text
T_9
 =213003y^9+639009y^8
  +7938(720B+161)y^7
  +490(28755B+3808)y^6
  +105(466560B^2+209790B+19061)y^5
  +21(3936600B^2+1162350B+83447)y^4
  +(78732000B^3+114434775B^2+25401600B+1539041)y^3
  +3(78732000B^3+40441275B^2+6495930B+333739)y^2
  +1500(54B+7)^3y
  +500(54B+7)^3.                                  (11)
```

The endpoint value is

```text
T_9(1)=32000(27B+7)^3.                             (12)
```

Suppose first that `27B+7!=0`. Then `y=1` is a simple root of
`R_10`. In (1) it must belong to `H_2`, not `S_4`. Dividing that
linear factor gives

```text
T_9=constant * H_1 S_4^2,

deg gcd(T_9,T_9')>=4.                              (13)
```

Let

```text
Sres_3(T_9,T_9')=A_3(B)y^3+A_2(B)y^2+A_1(B)y+A_0(B),
```

and take the primitive integer parts of `A_3,A_2`. They have degrees
`15,15`. Their reductions modulo `11` are

```text
abar_3
 =-B^15+B^14-2B^13-2B^12+2B^11-4B^10+4B^9+2B^8
   -3B^6+B^5-B^4+2B^3+4B^2+3B+3,

abar_2
 =-3B^15+2B^14-2B^13-B^11+B^10-3B^9+4B^8+B^7+B^6
   +4B^5+2B^4-3B^3+4B^2-5B-2.                    (14)
```

Both degrees are preserved modulo `11`, and an exact Sylvester
determinant gives

```text
Res_B(abar_3,abar_2)=6 mod 11.                     (15)
```

Therefore the integer resultant of the characteristic-zero primitive
polynomials is nonzero. They have no common complex root, so the
degree-three subresultant can never vanish identically. This contradicts
(13).

At the omitted value

```text
B=-7/27,
```

the residual factors exactly as

```text
R_10
 =7^6(y-1)^4
   (621y^6+3726y^5+8721y^4+10396y^3
    +7536y^2+3000y+500).                          (16)
```

The sextic in (16) is squarefree: its resultant with its derivative is
`2 mod 7`. Consequently

```text
gcd(R_10,R_10')=(y-1)^3,

deg gcd(R_10,R_10')=3<4.                           (17)
```

This also contradicts (1) and (8). Thus

```text
q_5(1)!=0                                           (18)
```

at every `H_2 S_4^2` point.

## 3. The complete common-root resultant

An exact Sylvester determinant using (3)--(4) gives

```text
Res_y(p_3,q_5)
 =7061881225000
  (54B+7)(1215B+91) Psi(B,C),                     (19)
```

where

```text
Psi
 =767400804B^4-172777374B^3+1750329B^2
  +65086642152BC^2+1562436540BC+6260436B
  +16874314632C^2+117021996C+405769.               (20)
```

The polynomial `Psi` is irreducible over `Q`. Its discriminant as a
quadratic in `C` is

```text
Disc_C(Psi)
 =-4821232752(54B+7)^3(513B-91)^2.                 (21)
```

Equations (6) explain the first two components: they are the common
roots `0` and `-1`. The remaining factor is the norm of `q_5` over the
quadratic factor of `p_3`. We now eliminate all three factors in (19).

### 3.1. The zero-root component

THM-2357 gives

```text
126D-25B^2=-(35/972)(54B+7).                       (22)
```

Hence `54B+7=0` is precisely the common-root wall of THM-2345.
Every low-square-class point on that wall lies on the deeper wall.
THM-2338 proves that the deep wall has exactly three square-class
candidates and that all three are in the `H_4`, not `H_2`, stratum.
Therefore

```text
54B+7!=0                                            (23)
```

at every point of (1). This use is purely at the square-class level;
no Keller flux consequence is needed.

### 3.2. The antipodal component

Put

```text
B=-91/1215.
```

Then

```text
R_10=(2401/729)(y+1)^2 E_8,                        (24)
```

with

```text
E_8
 =22182741y^8-214326y^6+2946308904Cy^5
  -11638431y^4-1696359672Cy^3
  +(502096953744C^2-8207696)y^2-1344364.           (25)
```

Because `p_3(-1)=q_5(-1)=0`, the root `-1` has order at least two in
`R_10`. Equations (1) and (8) make `y+1` divide `S_4`, so

```text
E_8=constant * H_2 T_3^2,

deg gcd(E_8,E_8')>=3.                              (26)
```

Let the degree-two subresultant be

```text
Sres_2(E_8,E_8')=e_2(C)y^2+e_1(C)y+e_0(C).
```

The involution `(y,C)->(-y,-C)` makes `e_2` even and `e_1` odd.
After taking primitive parts, write

```text
e_2(C)=f(C^2),                 e_1(C)=C g(C^2).
```

The value `C=0` is excluded separately: the reduction of `f(0)` modulo
`17` is `-6`, so it is nonzero in characteristic zero. For `C!=0`, put

```text
X=C^2,
```

so the two necessary equations become `f(X)=g(X)=0`, of degrees `5`
and `4`.
Their degree-preserving reductions modulo `17` are

```text
f_5=3X^5-2X^4+3X^3-7X^2+4X-6,

g_4=X^4-3X-2.                                      (27)
```

The exact modular resultant is

```text
Res_X(f_5,g_4)=5 mod 17.                           (28)
```

Thus the characteristic-zero equations are coprime and (26) is
impossible. The complete antipodal component is empty in the `H_2`
stratum.

### 3.3. The moving quadratic-root component

Suppose `Psi(B,C)=0` away from the two factors already treated. Since
the roots of `p_3` are `-1` and the two roots of

```text
7y^2+54B+7,
```

there is a common root `a` in the quadratic pair. The excluded factors
give

```text
a!=0,-1.
```

Solving `p_3(a)=q_5(a)=0` gives the exact rational normalization

```text
B=-7(a^2+1)/54,

C=
 7(19a^4+19a^3+64a^2+19a+19)
 /[26244(a+1)].                                    (29)
```

This parametrization loses no affine point of `Psi=0` in the present
chart. Indeed, a common root exists by (19); it cannot be `0` or `-1`,
and the coefficient of `C` in `q_5(a)` is nonzero for
`a!=0,-1`, so (29) follows. Conversely, substitution makes `Psi`
vanish identically.

The apparent exceptional parameters are already accounted for:

```text
a=0       gives 54B+7=0;

a=-1      would give B=-7/27 but q_5(-1)=3136, so
          it is only a point at infinity of (29), not an affine point;

a=1       gives (B,C)=(-7/27,245/13122), hence K=0. (30)
```

At the common root, (8) makes `y-a` divide `S_4`. Put

```text
(a+1)^2 R_10=7^6(y-a)^2 U_8.                      (31)
```

The quotient is

```text
U_8
 =621(a+1)^2y^8
  +1242(a+1)^3y^7
  -27(a+1)^2(11a^2-92a+11)y^6
  -2(a+1)(709a^4+1006a^3-110a^2+1006a+709)y^5
  -(a+1)^2(139a^4+2836a^3-1415a^2+2836a+139)y^4
  -2a(a+1)(139a^4+1203a^3-1008a^2+1203a+139)y^3
  +a^2(361a^4+1012a^3+5398a^2+1012a+361)y^2
  +1000a^3(a+1)^3y
  +500a^4(a+1)^2.                                 (32)
```

Its leading coefficient is `621(a+1)^2`, nonzero in this chart.
Since `(a+1)^2` is a square,

```text
U_8=constant * H_2 T_3^2,

deg gcd(U_8,U_8')>=3.                              (33)
```

Write

```text
Sres_2(U_8,U_8')=rho_2(a)y^2+rho_1(a)y+rho_0(a).
```

The exact maximal endpoint factors of the first two primitive
coefficients are

```text
rho_2=a^2(a+1)^12 A_36(a),

rho_1=a^3(a+1)^13 B_34(a).                        (34)
```

Both stripped polynomials are reciprocal. Since `a!=0`, put

```text
z=a+a^(-1),

A_36(a)=a^18 Ahat_18(z),

B_34(a)=a^17 Bhat_17(z).                           (35)
```

The degree-preserving reductions modulo `11` are

```text
Ahat_18
 =-z^18+3z^17+5z^16-2z^15+z^14+2z^13-3z^12-z^11
   +3z^10+z^9-2z^7-z^6+2z^5-5z^4-z^3-2z^2-2z+1,

Bhat_17
 =2z^17-5z^16-3z^15+4z^13+3z^12+4z^11+3z^10-5z^9
   +z^8-2z^7-5z^6+z^5+3z^4-4z^3-z^2+4z-4.        (36)
```

A custom Sylvester determinant gives

```text
Res_z(Ahat_18,Bhat_17)=9 mod 11.                   (37)
```

Thus `Ahat_18` and `Bhat_17`, and hence `rho_2,rho_1` away from
`a=0,-1`, have no common complex zero. The necessary subresultant in
(33) cannot vanish. This eliminates the final component `Psi=0`.

Combining (19), (23), (28), and (37) proves

```text
Res_y(p_3,q_5)!=0                                  (38)
```

on the complete `H_2 S_4^2` locus.

## 4. The cube descent is now unconditional

Equations (18) and (38) supply exactly the two hypotheses of THM-2360.
After sending the two roots of `H_2` to `-1,+1`, put

```text
t^2=x^2-1,

s=x+t,                         s^(-1)=x-t.
```

Then every degree-eighteen `H_2` point satisfies

```text
B_10(s)
 :=s^5(7q_5+S_4t)
 =c_0(s-r_0)U_3(s)^3,                             (39)

deg U_3=3,                       c_0 r_0 U_3(0)!=0.
```

Here

```text
r_0+r_0^(-1)=2x_0,

P_6(s)=s^3p_3((s+s^(-1))/2),

B_10 B_10^vee=-4X_0P_6^3,

X_0=(s^2-2x_0s+1)/2.                              (40)
```

Repeated roots of `p_3` and an overlap `p_3(x_0)=0` remain allowed, as
in THM-2360. The common-root theorem removes nonmaximal overlap between
`p_3,q_5`; it does not add squarefreeness of `p_3`.

The smallest current coefficient realization of the remaining problem
is obtained by writing

```text
U_3=s^3+u_2s^2+u_1s+u_0
```

after absorbing its leading coefficient into `c_0`, and imposing

```text
[s^j] B_10
 =c_0 [s^j](s-r_0)(s^3+u_2s^2+u_1s+u_0)^3,

j=0,...,10,                                        (41)
```

together with (40) and the affine pullback of the structured
degree-eighteen coefficients (3)--(4). Equivalently, before normalizing
`H_2`, THM-2357's top-down comparison leaves the six exact equations

```text
E_5=...=E_0=0 in Q[B,C,s_3,s_2],                  (42)
```

with `E_5` linear in `C`. This theorem proves that every
squarefree-`H_2` solution of (42) lies in the localization

```text
K * Res_y(p_3,q_5) !=0.                            (43)
```

The lawful next target is the coprime coefficient system (39)--(43).
Its emptiness is **OPEN**. Thus THM-2371 closes both auxiliary-hypothesis
gaps but does not close the `H_2` stratum.

The `H_4` branch is also untouched except on walls already closed by
THM-2338, THM-2345, THM-2347, and THM-2359. In particular, no
off-wall `H_4 S_4^2` reduction follows from the quadratic Laurent UFD.
Degree eighteen, `JC(2)`, and `DC(2)` remain open.

## 5. Exact reproduction

Run

```bash
python3 04-computation/jc2_degree18_h2_common_root_elimination_thm2371.py
python3 -O 04-computation/jc2_degree18_h2_common_root_elimination_thm2371.py
```

Both transcripts are byte-identical to

```text
05-knowledge/results/jc2_degree18_h2_common_root_elimination_thm2371.out
```

The companion reconstructs (3)--(6), the complete resultant
(19)--(21), the forced-root quotients (10), (24), and (31), every
subresultant degree profile, the exceptional factorization (16), the
exact factor orders (34), the parity and reciprocal quotients, and all
displayed nonzero modular resultants. Each main modular resultant is computed
both by SymPy and by a separate Gaussian determinant of the Sylvester
matrix. Positive square-class controls and raw-wall hostile controls
are included. No executable check uses Python `assert`.

Independent hostile audit is pending. QED.
