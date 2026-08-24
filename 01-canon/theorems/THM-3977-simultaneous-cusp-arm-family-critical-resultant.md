---
id: THM-3977
title: "Simultaneous cusp-arm family critical resultant"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. On the
  height-two exact-volume completion, every member of the one-parameter
  simultaneous cusp/arm jet family has an affine critical point in its
  first coordinate. The exact derivative resultant is x^24 times a degree-12
  residual whose leading and constant coefficients are both nonzero for
  a!=0. An independent u=x^2t normalization gives a degree-12 critical
  address equation and reconstructs an actual affine critical point from
  every root. Thus no companion, including the displayed boundary-jet
  companion, can give a polynomial Darboux pair. The endpoint a=0 first
  coordinate is also critical, although the companion itself is undefined.
  The obstruction survives both the delayed x^5 boundary payment and the
  actual two-term formally lifted first coordinate: their reduced residual
  resultants have degree 43 with nonzero leading and constant terms. A
  separately audited lowest-seam supplement treats
  A_(c,r)=c p+y^2+r x for every c!=0: r!=0 is critical, while r=0 is a
  source-plane submersion whose generic fibre forbids every nonzero
  invariant rational Jacobian value.
source: jc-degree6-one-place / post-THM-3973 simultaneous boundary-jet hostile, 2026-08-24
audit: >
  PASS (root / jc-cohn3709, 2026-08-24). The audit independently rederived
  both completion-chart jet restrictions and the normalized D-chart volume
  ratio, checked the full derivative resultant and the finite-root argument
  using the stable t-leading rows, and reconstructed the intrinsic
  degree-twelve u-cover including all excluded colors and the Bezout x
  address. It separately checked the a=0 endpoint and both degree-43 lifted
  resultants after quotient-ring reduction modulo 6L^2-1. Normal, optimized,
  and frozen outputs byte-match at CHECKS=44; all hashes agree. No formal
  completion obstruction or unrestricted first-coordinate claim was used.
  A separate root/Huygens audit rederived the generalized (c,r) resultant,
  the r=0 submersion endpoint, and the six-residue generic-fibre obstruction;
  its 25-gate normal and optimized runs byte-match their independently frozen
  companion.
depends_on:
  - THM-3973-exact-volume-simple-cubic-determinantal-affine-plane-completion
related:
  - THM-3975-danielewski-one-arm-modification-cubic-control-and-hyperelliptic-no-mate
  - THM-3976-rational-compression-quadratic-pseudoplane-intersection
script: 04-computation/jc2_simultaneous_cusp_arm_critical_resultant_thm3977.py
output: 05-knowledge/results/jc2_simultaneous_cusp_arm_critical_resultant_thm3977.out
script_sha256: 0ba04625c793c316e5ef60969efd64ce32ec4c4fdc281f9da60feac80f95f29b
output_sha256: 533f1aa2576842f4b4254931176c479dd4bacf21a493c54e70596a1d0106535c
semantic_sha256: 5aee4e0e5b0ed5c710ff3f2dcc204f9afeade620bc350ceeafee8f3613c0ce63
hash_basis: raw LF bytes
---

# THM-3977 -- simultaneous boundary jets force an interior critical point

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over
an algebraically closed field `k` of characteristic zero. Fix

```text
6L^2=1,                         a in k^*,                       (1)
```

and on the source plane put

```text
z=1+x^2t,             p=zt,             y=xzt^2.              (2)
```

These generate the height-two completion `B_2` of THM-3973. Consider

```text
A_a=y^2+2Lx+ap,                                                (3)

C_a=y^3+3Lxy+(3a/2)py(1-z)-xz/a
       +x^2(-3aL/2+3a^2y/8).                                  (4)
```

The pair has the intended cusp and arm jets, but it is never Darboux:

```text
for every a!=0, there is (x_0,t_0) in A2 with
(partial_x A_a)(x_0,t_0)=(partial_t A_a)(x_0,t_0)=0.           (5)
```

Consequently `J(A_a,G)` vanishes at `(x_0,t_0)` for every source polynomial
`G`; in particular `J(A_a,C_a)!=1`. This closes the whole displayed family
at its first coordinate, not merely at the forced low-order companion.

## 1. The simultaneous cusp and arm jets are genuine

Let `X_2=Spec(B_2)`. Its boundary and the second component of `div(x)` are

```text
D=V(x,z,p) isomorphic to A1_y,
L_1=V(x,z-1,y) isomorphic to A1_p,                       (6)
```

where `L_1` lies inside the source open `A2_(x,t)`. Direct restriction gives

```text
(A_a,C_a)|D   =(y^2,y^3),
(A_a,C_a)|L_1 =(ap,0).                                  (7)
```

Thus `D` maps with the standard cusp parametrization while `L_1` maps to
the target line `C=0` (the cusp tangent at the origin). The source
calculation gives

```text
J(A_a,C_a)-1 in x k[x,t,a,a^-1,L]/(6L^2-1),             (8)
```

so the Jacobian is one along the interior arm `L_1`. Statement `(8)` alone
does **not** control `D`, because `t` has a pole there.

For the actual boundary calculation use the smooth `(x,y)` chart. The
determinantal relations give

```text
z(z-1)^2=x^3y,                    p=xy/(z-1).            (9)
```

Since `z` has `x`-order three along `D`, equations `(3)--(4)` become

```text
A_a =y^2+x(2L-ay)                                      mod x^3,
C_a =y^3+x(3Ly-3ay^2/2)
          +x^2(-3aL/2+3a^2y/8)                        mod x^3. (10)
```

The coefficient of `x` in their `(x,y)` Jacobian is

```text
(2L-ay)(3L-3ay)
 -4y(-3aL/2+3a^2y/8)
 +a(3Ly-3ay^2/2)
=6L^2=1.                                                (11)
```

THM-3973 gives

```text
dx wedge dt =x/((z-1)(3z-1)) dx wedge dy.               (12)
```

The denominator restricts to one on `D`, so `(10)--(12)` prove that the
Jacobian ratio also extends with value one along `D`. The simultaneous jet
invoice is therefore real; its failure occurs in the affine interior.

## 2. Exact derivative resultant

In source coordinates,

```text
A_a=2Lx+a(t+x^2t^2)+x^2t^4+2x^4t^5+x^6t^6.             (13)
```

Sylvester elimination in `t`, reduced by `6L^2=1`, gives

```text
Res_t(partial_x A_a,partial_t A_a)=(32/3)x^24 H_a(x),   (14)

H_a(x)=3888Lx^12+864ax^9+144La^2x^6-256x^5
       +36La^5x^4+48a^3x^3+9a^6x-36La^4.               (15)
```

For audit, before reducing powers of `L`, the same factorization is

```text
96x^24(15552L^5x^12+3456L^4ax^9-1024L^4x^5
       +96L^3a^2x^6+32L^2a^3x^3+4La^5x^4
       -4La^4+a^6x).                                    (16)
```

For `a!=0`, the residual has

```text
lc_x(H_a)=3888L!=0,              H_a(0)=-36La^4!=0.     (17)
```

It therefore has a nonzero root `x_0` in `k`. At every `x_0!=0`, the
`t`-leading coefficients of the two derivatives are respectively

```text
6x_0^5                         and                     6x_0^6. (18)
```

Their degrees remain six and five, so the resultant root is not a common
root at projective `t=infinity`: it gives a finite `t_0` satisfying `(5)`.
This already proves the theorem.

## 3. The intrinsic degree-twelve critical-address cover

There is a second proof that explains the residual degree. On `x!=0` put

```text
u=x^2t,              f=u^4(u+1)^2,              g=u(u+1). (19)
```

Then

```text
A_a=x^-6 f+a x^-2 g+2Lx.                                (20)
```

Since `partial_t=x^2 partial_u`, at a zero of `partial_u A_a` the fixed-`t`
and fixed-`u` `x` derivatives agree. Thus the critical equations are exactly

```text
2u^3(u+1)(3u+2)+a x^4(2u+1)=0,
-6u^4(u+1)^2-2a x^4u(u+1)+2Lx^7=0.                     (21)
```

Away from the four displayed colors they give

```text
x^4=X_4:=-2u^3(u+1)(3u+2)/(a(2u+1)),
x^7=X_7:=-u^4(u+1)^2/(L(2u+1)).                         (22)
```

Compatibility `X_4^7=X_7^4`, using `L^4=1/36`, is precisely

```text
Phi_a(u):=9a^7(u+1)(2u+1)^3+32u^5(3u+2)^7=0.            (23)
```

This polynomial has degree twelve and constant term `9a^7`. Moreover,

```text
Phi_a(-1)=32,
Phi_a(-1/2)=-1/128,
Phi_a(-2/3)=-a^7/9,                                    (24)
```

so every root avoids `0,-1,-1/2,-2/3`. Choose any root `u_0`. Both rows in
`(22)` are then nonzero and compatible. Since `gcd(4,7)=1`, the explicit
choice

```text
x_0=X_4^2/X_7
   =-4L u_0^2(3u_0+2)^2/(a^2(2u_0+1)),
t_0=u_0/x_0^2                                           (25)
```

satisfies `x_0^4=X_4` and `x_0^7=X_7`, hence `(21)`. This reconstructs an
actual affine critical point rather than only detecting one by elimination.

## 4. The actual formal boundary lift is still globally critical

The boundary calculation can be continued formally, but that does not repair
the global submersion. Distinguish the delayed one-term test and the actual
two-term lifted first coordinate:

```text
A_a^(5)=A_a+(5L/6)x^5(1-z),

A_a^(full)=A_a+x^2(1-z)/4+(5L/6)x^5(1-z).              (26)
```

For `star in {5,full}`, let `Q_a^(star)` be the unique representative of
`x^-24 Res_t(partial_x A_a^(star),partial_t A_a^(star))` having `L`-degree
at most one modulo `6L^2-1`. The complete coefficient rows are frozen by the
exact companion and its semantic hashes. Their decisive endpoints are the
same:

```text
deg_x Q_a^(star)=43,
lc_x Q_a^(star)=-93750,
Q_a^(star)(0)=-384La^4.                                  (27)
```

The correction terms are only linear in `t`, so the two derivative leading
rows remain `6x^5` and `6x^6`. For every `a!=0`, each residual in `(27)` has
a nonzero root, and the stable-degree resultant argument again gives a finite
affine critical point. In particular the actual formally lifted first
coordinate cannot occur in a Darboux pair, regardless of how its companion
is coordinated. This is exact evidence that formal boundary lifting and
global submersion are orthogonal obligations. It is not an obstruction to a
formal Darboux normal form in either completed boundary chart: only the three
displayed global polynomial first coordinates `(3)` and `(26)` are tested
here.

## 5. Sharp endpoint and scope

The restriction `a!=0` is needed only because `(4)` contains `1/a`. The
first coordinate at `a=0` remains critical: take

```text
u=-2/3,             x^7=16/(243L),             t=u/x^2. (28)
```

Then both equations `(21)` with `a=0` vanish. Thus neither deleting nor
specializing the arm coefficient repairs the submersion failure.

What survives is exact and useful: the completion admits simultaneous
cusp/arm boundary jets with Jacobian residue one on both divisors. What
fails is the global submersion condition, detected by the interior
degree-twelve address cover `(23)`. No arbitrary first coordinate on `B_2`,
no further higher-order or formal modification beyond the two explicitly
tested coordinates in `(26)`, no arbitrary Darboux pair, and no consequence
for unrestricted `JC(2)` is claimed.

## 6. Generalized lowest-seam supplement

The first-coordinate obstruction is not tied to the normalization
`6L^2=1`. For arbitrary `c in k^*` and `r in k`, put

```text
A_(c,r)=c p+y^2+r x.                                    (29)
```

It has the same two boundary colors,

```text
A_(c,r)|D=y^2,                  A_(c,r)|L_1=c p.        (30)
```

Direct elimination gives

```text
Res_t(partial_x A_(c,r),partial_t A_(c,r))
  =96x^24 H_(c,r)(x),                                   (31)

H_(c,r)=c^6x+2c^5rx^4-2c^4r+8c^3r^2x^3
        +12c^2r^3x^6+216cr^4x^9+486r^5x^12
        -64r^4x^5.                                      (32)
```

If `r!=0`, then `H_(c,r)(0)=-2c^4r!=0`, while its linear
coefficient is `c^6!=0`. Over the algebraically closed field it therefore
has a nonzero root. At every nonzero `x`, the two derivative polynomials
retain `t`-degrees six and five with leading rows `6x^5` and `6x^6`.
Thus that root reconstructs a finite affine source critical point. The
primary family `(3)` is the specialization `(c,r)=(a,2L)`.

At the omitted seam `r=0`, the exact resultant is instead

```text
Res_t(partial_x A_(c,0),partial_t A_(c,0))=96c^6x^25.   (33)
```

The only possible `x`-coordinate is zero, where
`partial_t A_(c,0)=c`; hence `A_(c,0)` is a source-plane submersion.
It still has no rational mate with nonzero invariant Jacobian. On the
generic fibre `A_(c,0)=P`, with `Y=y`, one has

```text
k(P)(X_P)=k(P)(Y),
E_P(Y)=(P-Y^2)^3-c^3Y^2,                                (34)

Disc_Y(E_P)=64P^3c^12(27P^2+4c^3)^2.                   (35)
```

The sextic is squarefree over `k(P)` and coprime to `P-Y^2`. If
`J(A_(c,0),Q)=h(P)` with `h(P)!=0`, restriction to the generic fibre would
make the differential

```text
-c h(P)(P-Y^2)dY/E_P(Y)                                 (36)
```

exact. It has six simple poles with nonzero residues, which is impossible
for the differential of a rational function. This closes the entire
lowest-seam family `(29)`: nonzero `r` fails by criticality, whereas
`r=0` fails by logarithmic exactness debt.

The independent companion
`04-computation/jc2_generalized_cusp_arm_seam_thm3977_independent.py`
has 25 optimization-safe gates. Its raw hashes are
`8fc3507c40e31a54117ba96c7f56261f49386a0a61c5769800d3ae5ef56b9451`
for the script and
`1534633377e5722a13a933db89da66a22cb480e22f0f042acd753642c0b78c72`
for the output; the semantic hash is
`8aadeddb4635f490ea0acde0694b539d41f6145ac66e93df36af5543c8e6ae47`.

**QED.**

## Reproduction

```bash
python3 04-computation/jc2_simultaneous_cusp_arm_critical_resultant_thm3977.py
python3 -O 04-computation/jc2_simultaneous_cusp_arm_critical_resultant_thm3977.py
python3 agents/check_docs.py
```
