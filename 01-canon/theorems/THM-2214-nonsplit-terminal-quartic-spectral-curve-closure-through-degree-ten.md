---
id: THM-2214
title: "Nonsplit terminal quartic spectral-curve closure through reduced degree ten"
status: >
  PROVED + VERIFIED-EXACT + HOSTILE-AUDITED. Let P=H^2+L be a polynomial
  exact square-prefix quartic in one source fibre, with H quadratic, L
  linear, and nonsquare leading coefficient for H over C(x). If a polynomial
  Keller mate has reduced fibre degree 2, 6, or 10, then degree two is tame
  by THM-2071 and degrees six and ten are impossible. In degree ten the two
  constant Faber fluxes cut out a spectral cubic. A smooth cubic admits no
  nonconstant rational P^1-trajectory; singular normalization, the Keller
  one-form, and polynomiality of H eliminate the nodal and cuspidal cases.
  Thus a nonsplit terminal survivor has reduced degree at least 14. Split
  degree at least 6, nonsplit degree at least 14, other short edges, JC(2),
  and DC(2) remain open.
source: codex-2026-07-24-nonsplit-terminal-spectral-curve
depends_on:
  - THM-2071-planar-keller-quadratic-pencil-normal-form
  - THM-2129-quartic-faber-three-coefficient-boundary-classification
  - THM-2180-quartic-first-pole-stratum-is-empty
  - THM-2189-nonsplit-quartic-deck-forces-the-remaining-pole-congruence
  - THM-2202-uniform-all-degree-quartic-pole-closure
related:
  - THM-2194-uniform-degree-six-quartic-pole-closure
  - THM-2206-quadratic-deck-augmentation-and-hamiltonian-characteristic-incompatibility
script: 04-computation/jc2_nonsplit_terminal_degree10_spectral_curve_thm2214.py
output: 05-knowledge/results/jc2_nonsplit_terminal_degree10_spectral_curve_thm2214.out
script_sha256: f17728dae43d5160c6fc38ad7a858cf6f99b5d46b0c04634a8674696e5d3f18f
output_sha256: 879b9f8bc93139de3659eb13d4f3f5cccd936202a1b3e5a73ace9a1b18db78d3
hash_basis: working-tree bytes (LF)
---

# THM-2214 -- nonsplit terminal quartics through degree ten

This theorem concerns only the nonsplit quartic deck with reduced mate degree
`2`, `6`, or `10`. It does **not** handle the split deck, reduced degree at
least `14`, arbitrary quartic source fibres, `JC(2)`, or `DC(2)`.

The main new object is not another pole coefficient. Once the quadratic
approximate root is polynomial, the two constant Faber fluxes cut out a
one-dimensional **spectral curve** in the three translation-invariant
coordinates of the exact square prefix. In reduced degree ten that curve is a
plane cubic. A smooth cubic cannot carry a nonconstant `C(x)` trajectory.
The singular normalization leaves one rational parameter, but the Keller
one-form and polynomiality of the approximate root eliminate every such
trajectory.

## 1. Setup and exact statement

Let

```text
R=C[x],                 K=C(x),
```

where `C` is algebraically closed of characteristic zero. Suppose

```text
P=H^2+L,

H=V z^2+B z+C_0 in R[z],       V!=0,
L=A z+E in R[z],                                      (1)
```

and let `Q in R[z]` satisfy

```text
J_(x,z)(P,Q)=kappa in C*.                             (2)
```

Assume that `V` is not a square in `K`. Reduce `Q` by target shears and
suppose its remaining `z`-degree is

```text
n in {2,6,10}.                                        (3)
```

Then no noninvertible Keller pair survives. More precisely:

* if `n=2`, the reduced mate is `H` up to a target shear, so `(P,H)` is a
  Keller pair with a quadratic member and THM-2071 closes it;
* if `n=6` or `n=10`, assumptions (1)--(3) are contradictory.

Thus a nonsplit terminal quartic survivor, if one exists after THM-2202, has
reduced mate degree at least `14`.

The proof starts with THM-2180/2189's constant-coefficient Faber normal form,
uses THM-2189's nonsplit deck parity, and uses THM-2129's exact Hamiltonian
identity. It does not reuse the finite-pole argument of THM-2194: that theorem
makes `H` polynomial in degree six, while the present argument starts from
that polynomiality and closes the subsequent terminal branch. THM-2202
supplies the required polynomial `H` uniformly in every twice-odd degree.

## 2. What the exact square prefix does to every even Faber seed

For `s>=1`, put `m_s=4s-2` and define the polynomial

```text
R_s(P,H)
 =sum_(k=0)^(s-1) binom(s-1/2,k)
      H^(2s-1-2k)(P-H^2)^k.                          (4)
```

Writing `L=P-H^2`, the binomial expansion at fibre infinity gives

```text
E_(m_s)=R_s(P,H)+S_s,             deg_z(S_s)<=s-2,   (5)
```

where `E_(m_s)=Pol_z(P^(s-1/2))`. The first omitted term is

```text
binom(s-1/2,s) L^s/H,                                (6)
```

whose degree at infinity is `s-2`; every subsequent term loses three more
fibre degrees.

The large expression (4) has a small derivative. At fixed `P`,

```text
partial_H R_s
 =binom(s-1/2,s-1)L^(s-1).                           (7)
```

This follows by differentiating (4): all intermediate terms telescope, or
equivalently by putting `X=H/P^(1/2)` and checking

```text
d/dX sum_(k<s) binom(s-1/2,k)
 X^(2s-1-2k)(1-X^2)^k

 =binom(s-1/2,s-1)(1-X^2)^(s-1).
```

Consequently

```text
J(P,R_s(P,H))
 =-binom(s-1/2,s-1)L^(s-1) J(H,L).                  (8)
```

This is the exact `C[P,H]` compression: the top Faber seed does not expose
its nominal degree after taking a Jacobian. The correction `S_s` is the
sidecar that cannot be discarded.

For degrees two, six, and ten that sidecar is explicit. The coordinate
calculation in the next section gives

```text
E_2 =H,

E_6 =R_2(P,H)+(3/8)T,

E_10=R_3(P,H)+(5/16)T(L-2s_0),                      (9)
```

where `T=A^2/V` and `s_0` is defined in (12). Thus degree six really is
`C[P,H]` plus an `x`-only remainder. Degree ten adds exactly one copy of
the linear remainder, with coefficient `T`. The ratio `A^2/V` is the first
missing module coordinate.

## 3. Intrinsic centered coordinates

First dispose of `A=0`. Then `L=E(x)` and every omitted binomial term in
(4) has negative fibre degree, so every even Faber seed belongs to
`C[P,H]`. Deck parity removes all odd seeds, hence the reduced `Q` belongs
to `C[P,H]`. Equation (2) becomes

```text
kappa=Q_H(P,H) J(P,H)
     =Q_H(P,H) E'(x) H_z.                            (10)
```

The right side cannot be a unit: if `E'=0` it is zero, and otherwise `H_z`
has fibre degree one. Hence `A!=0`.

Express `H` in the linear coordinate `L`:

```text
a=V/A^2,
s_0=AB/(2V)-E,
d=C_0-B^2/(4V).                                     (11)
```

Then the square completion is exact:

```text
H=a(L+s_0)^2+d.                                      (12)
```

Adjoin `u` with `u^2=a`. Because `V` is nonsquare and `A^2` is a square,
this is a genuine quadratic field. Put

```text
w=u(L+s_0),
q=1/u,
T=q^2=1/a=A^2/V.                                    (13)
```

Then

```text
H=w^2+d,
L=q w-s_0,
P=w^4+2d w^2+q w+(d^2-s_0).                         (14)
```

If `U^2=V`, one may take `u=U/A`; then

```text
w=Uz+B/(2U),                 partial w/partial z=U. (15)
```

The deck involution sends `(u,w,q)` to `(-u,-w,-q)` and fixes
`d,s_0,T`. Equation (2) becomes

```text
J_(x,w)(P,Q)=kappa/U.                                (16)
```

Here and below a prime denotes the derivation of the coefficient field with
`w` held fixed.

The exact divisions

```text
Pol_L(L^2/H)=T,
Pol_L(L^3/H)=T(L-2s_0)                               (17)
```

prove (9) directly from the binomial series. All later terms have negative
fibre degree.

## 4. The exact Laurent bank

For a depressed quartic

```text
P=w^4+p w^2+q w+r
```

write

```text
P^(m/4)=E_m+c_(m,1)w^-1+c_(m,2)w^-2+c_(m,3)w^-3+...,

Phi_m=4c_(m,1),
Psi_m=4c_(m,2),
R_m=4c_(m,3)+p c_(m,1).                             (18)
```

THM-2129 gives

```text
J(P,E_m)=(w^2+p/4)Phi_m'+w Psi_m'+R_m'.             (19)
```

Substitute `p=2d`, `q^2=T`, and `r=d^2-s_0`. The fibre derivative recurrence

```text
P partial_w(P^(m/4))
 =(m/4)(partial_w P)P^(m/4)
```

gives

```text
degree 2:
 Phi_2=2q,
 Psi_2=-2s_0,
 R_2=-dq;

degree 6:
 Phi_6=-3q s_0,
 Psi_6=-(3/2)(dT-s_0^2),
 R_6=q(6d s_0-T)/4;

degree 10:
 Phi_10=-(5/4)q(dT-3s_0^2),
 Psi_10=-(5/32)(-24dT s_0+T^2+8s_0^3),
 R_10=(5/8)q(d^2T-3d s_0^2+T s_0).                 (20)
```

These are polynomial identities, not a generic calculation.

For an even reduced Faber combination, (19) and (16) first give

```text
Phi_Q', Psi_Q'=0.                                    (21)
```

The constant field of the algebraic function field remains `C`. The deck
fixes constants, while `Phi_Q` is anti-invariant. Therefore

```text
Phi_Q=0,                  Psi_Q=psi_0 in C,
R_Q'=kappa/U.                                         (22)
```

This is the terminal use of nonsplit parity. On a split deck `Phi_Q` may be
a nonzero constant, so the argument below does not apply.

## 5. Reduced degree two

After normalization and a target shear,

```text
Q=E_2=H.
```

Thus `J(P,H)=kappa`. The pair has a quadratic member, and THM-2071 gives
the tame normal form. This case does not require the remaining sections.

## 6. Reduced degree six

Normalize the coefficient of `E_6` and write

```text
Q=J_0(P)+E_6+alpha E_2,             alpha in C.      (23)
```

The first equation in (22) and (20) gives

```text
s_0=2alpha/3.                                         (24)
```

The second says that `dT` is constant. Substitution into the third
observable makes every `d`-term cancel:

```text
R_Q=-qT/4=-q^3/4.                                    (25)
```

Since `q^2=T` and `q=A/U`, (22) becomes

```text
A T'=-8kappa/3.                                      (26)
```

Formula (9) and the polynomiality of `P,H,Q` show that

```text
(3/8)T=Q-J_0(P)-R_2(P,H)-alpha H
```

belongs to `C[x]`. Hence both `A` and `T` are polynomials. Their product in
(26) is a nonzero constant, so `A` is constant and `T` is nonconstant
affine-linear. But

```text
V=A^2/T                                                (27)
```

cannot then be a polynomial. This contradiction closes reduced degree six.

## 7. Reduced degree ten and the spectral cubic

Normalize the top coefficient and write

```text
Q=J_0(P)+E_10+alpha E_6+beta E_2,
alpha,beta in C.                                     (28)
```

The first flux in (22) gives

```text
dT=3s_0^2-(12alpha/5)s_0+8beta/5.                   (29)
```

Put

```text
y=5s_0-2alpha.                                       (30)
```

Eliminating `dT` from `Psi_Q=psi_0` gives the exact cubic

```text
125T^2
 =64y^3+(640beta-192alpha^2)y
   +128alpha^3-640alpha beta-800psi_0.               (31)
```

The third observable collapses at the same time:

```text
R_Q=qTy/8.                                           (32)
```

Differentiating (32), using `q^2=T` and `q=A/U`, turns the last equation in
(22) into

```text
A(2T y'+3y T')=16kappa.                              (33)
```

Equations (31) and (33) are the complete degree-ten terminal carrier.

### 7.1 The smooth cubic is impossible

If the cubic (31) is nonsingular, its smooth projective completion has genus
one. A rational pair `(y,T) in K^2` defines a rational map from `P^1` to
that completion. Properness extends it across every pole. A nonconstant map
is impossible by Riemann--Hurwitz:

```text
-2=deg(f)(2*1-2)+Ramification>=0.
```

Thus `y` and `T` are constant. This contradicts (33). Every survivor must
therefore lie on the discriminant of (31).

### 7.2 Exact normalization of the singular cubic

Write the right side of (31) as

```text
64y^3+p y+r_0.
```

On the discriminant there is an `e in C` such that

```text
64y^3+p y+r_0=64(y-e)^2(y+2e).                      (34)
```

The trajectory cannot be the singular point `(y,T)=(e,0)`, because `T` is
nonzero. Define

```text
rho=T/(y-e).                                         (35)
```

Then (34) has the rational normalization

```text
y=125rho^2/64-2e,
T=rho(125rho^2/64-3e).                               (36)
```

Substitution into (33) gives

```text
A (F_e(rho))'=kappa,                                 (37)

F_e(rho)
 =40625rho^5/65536
  -1625e rho^3/1024
  +9e^2 rho/8.
```

Indeed

```text
F_e'(rho)
 =(203125rho^4-312000e rho^2+73728e^2)/65536.        (38)
```

The parameter `rho` is nonconstant; otherwise (33) fails.

### 7.3 Rational integrability forces the triple cusp

We use the elementary rational-primitive lemma from THM-2071. If

```text
A in C[x]\{0},       S in C(x),       A S'=kappa,    (39)
```

then either:

* `A` is constant and `S` is affine-linear; or
* after translating `X=x-xi`,

  ```text
  A=a_0 X^m,       m>=2,
  S=gamma_0+gamma_1 X^(1-m),     gamma_1!=0.         (40)
  ```

For completeness, a rational derivative has no simple pole. If `A` had
`h` distinct roots and degree `D`, a rational primitive would have map
degree `D-h`, while its fibre at infinity has multiplicity `D-1`; hence
`D-1<=D-h` and `h=1`. Direct integration gives (40).

Apply this to `S=F_e(rho)`.

If `A` is constant, `F_e(rho)` is affine-linear. A finite pole of `rho`
would give a pole of order multiplied by five, so `rho` is a polynomial.
Then

```text
deg F_e(rho)=5 deg rho
```

cannot equal one.

Suppose `A` is nonconstant and put `t=X^-1`. Equation (40) says

```text
F_e(rho)=gamma_0+gamma_1 t^(m-1).                    (41)
```

Again `rho=R(t)` is a nonconstant polynomial. Differentiating gives

```text
F_e'(R(t))R'(t)
 =gamma_1(m-1)t^(m-2).                              (42)
```

Every zero on the left must therefore lie at `t=0`. If `e!=0`, the even
quartic (38) has at least two distinct roots. A nonconstant polynomial
`R` assumes each of those values at a finite point; (42) would force both
values to equal `R(0)`, impossible. An even quartic with only one distinct
root would be a scalar multiple of `(rho-rho_0)^4`; its cubic coefficient
forces `rho_0=0`, and then its constant coefficient forces `e=0`.

Hence

```text
e=0.                                                 (43)
```

Now `F_0` is a nonzero scalar multiple of `rho^5`. Equation (41) can be a
fifth power only if `gamma_0=0`, and therefore

```text
m-1=5k,              rho=lambda X^-k,
k>=1, lambda!=0.                                    (44)
```

### 7.4 The triple cusp restores an unavoidable pole in `H`

At `e=0`, the singularity conditions and (29), (36) give

```text
beta=3alpha^2/10,

T=125rho^3/64,
s_0=25rho^2/64+2alpha/5,
d=15rho/64,
a=1/T=64/(125rho^3).                                (45)
```

Return to the polynomial linear remainder `L=Az+E`. From (12),

```text
H=a(Az+E+s_0)^2+d.                                   (46)
```

At `X=0`, equation (44) makes `rho` a pole of order `k`. The terms involving
the polynomial `E` are regular:

```text
ord_X(aE^2)>=3k,
ord_X(2aE s_0)>=k.
```

The unique polar part of the constant coefficient in (46) is

```text
a s_0^2+d
 =(5/64)rho+(15/64)rho+regular
 =(5/16)rho+regular.                                 (47)
```

It cannot cancel. This contradicts the starting hypothesis `H in R[z]`.
Thus the singular cubic has no terminal Keller trajectory either.

Sections 7.1--7.4 close reduced degree ten.

## 8. Exact surviving frontier

The theorem proves the following strict improvement of the terminal ledger:

```text
nonsplit quartic deck
 + polynomial exact square prefix
 + reduced mate degree 2, 6, or 10
 -> tame at degree 2, empty at degrees 6 and 10.
```

It leaves:

```text
nonsplit reduced degree >=14,
split reduced degree >=6,
and the unrelated terminal short-edge branches.      (48)
```

The spectral-curve construction extends formally to every reduced degree:
`Phi_Q=0` and `Psi_Q=constant` cut a curve in `(d,s_0,T)`, while `R_Q`
supplies its Keller one-form. What is not yet proved is a uniform
classification of the higher-degree discriminant strata. The degree-ten
cubic shows the likely recursive object: normalize each singular flux curve,
pull back the Keller one-form, and compare its rational primitives with the
polynomiality of `H=a(L+s_0)^2+d`.

## 9. Exact referee

Run

```bash
python3 04-computation/jc2_nonsplit_terminal_degree10_spectral_curve_thm2214.py
python3 -O 04-computation/jc2_nonsplit_terminal_degree10_spectral_curve_thm2214.py
```

The script independently checks:

* the polynomial-part divisions in (9) and (17);
* every row of the Laurent bank (20);
* the same Laurent coefficients through a separately implemented finite
  multinomial expansion, so the Faber recurrence is not its own oracle;
* the degree-six reduction (24)--(26);
* the cubic and third-observable eliminations (29)--(32);
* the singular normalization and one-form (36)--(38); and
* the coefficient `5/16` in (47).

Both modes have byte-identical output stored in
`05-knowledge/results/jc2_nonsplit_terminal_degree10_spectral_curve_thm2214.out`.
The script is an exact symbolic referee only. Riemann--Hurwitz, rational
integrability, and the pole argument supply the mathematical quantifiers.
