---
id: THM-2429
title: "Degree-twenty-two C-W plane hyperelliptic family closure"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE AUDIT.
  In the open first-flux chart of the genuine nonsplit polynomial
  exact-square-prefix degree-twenty-two branch, the complete coefficient
  plane B=D=E=0 is empty. For C=0 this is THM-2423. For C!=0, the sole
  weighted ratio lambda=W/C^2 turns the first two fluxes into a quadratic
  whose completed square is a degree-five hyperelliptic family. Generic
  fibres have genus two. Its one rational double-root ratio, seven
  algebraic double-root ratios, and one degree-drop ratio all normalize to
  genus one. Hence no rational trajectory survives. This closes the C-W
  plane, not other mixed strata, degree twenty two, JC(2), or DC(2).
source: klein-2026-07-26-degree-twenty-two-cw-plane-triage
depends_on:
  - THM-2129-quartic-faber-three-coefficient-boundary-classification
  - THM-2214-nonsplit-terminal-quartic-spectral-curve-closure-through-degree-ten
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
  - THM-2411-degree-twenty-two-first-flux-pole-divisor-square-class-reduction
  - THM-2423-degree-twenty-two-W-axis-genus-two-and-origin-cusp-closure
related:
  - THM-2425-degree-twenty-two-CDE-axis-hyperelliptic-closure
  - THM-2428-degree-twenty-two-B-axis-trigonal-ramification-closure
script: 04-computation/jc2_degree22_cw_plane_hyperelliptic_family_thm2429.py
output: 05-knowledge/results/jc2_degree22_cw_plane_hyperelliptic_family_thm2429.out
script_sha256: 05fc2b14318f41ed8e11a78a4d619a731d6f1edb785c042e8e6da5f68f481a70
output_sha256: c5f1cce5764d2caf4f04122f6813fb89e8575f7a0448c66f65d1abacf54b453f
hash_basis: working-tree bytes (LF)
---

# THM-2429 -- the degree-twenty-two C-W plane is empty

**PROVED CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE AUDIT.**

The hostile-audited THM-2425 closes the `C`-axis, and candidate THM-2423
closes the `W`-axis. The axes are actually boundary fibres of one uniformly
positive-genus family. The exact conclusion here is

```text
genuine degree-22 trajectory,
mathcal A!=0,
B=D=E=0
    => contradiction.                                           (1)
```

Thus this is the first complete support-at-most-two coefficient plane in
the open degree-twenty-two chart.

## 1. The sole weighted ratio

Use the target-translated coordinates of THM-2411,

```text
y=11s,                   u=dT,                   Z=T^2,

wt(y,u,Z,B,C,D,E,W)=(1,2,3,2,3,4,5,6).                 (2)
```

If `C=0`, (1) is exactly the full `W`-axis of THM-2423. Hence suppose
`C!=0`. The coefficient plane has one weighted projective ratio,

```text
lambda=W/C^2 in C.                                      (3)
```

First take `y!=0` in `C(x)` and put

```text
v=u/y^2,            zeta=Z/y^3,            p=C/y^3.
                                                               (4)
```

Then `W/y^6=lambda p^2`. Dividing `N_1,N_2` by `y^5,y^6`
gives

```text
f_1
 =11979(7-121v)zeta
  +4[922383v^2-25410v+63+(2342560v-58080)p]
 =0,                                                        (5)

f_2
 =15944049zeta^2-162339408zeta v+2236080zeta
  -1190488992v^3+147581280v^2-1219680v+672
  +(-206145280zeta+449771520v-1239040)p
  -1319329792lambda p^2
 =0.                                                        (6)
```

The open chart says

```text
121v-7!=0,                                               (7)
```

so (5) reconstructs `zeta` uniquely. Exact elimination gives

```text
Res_zeta(f_1,f_2)=255104784 R_lambda(v,p),               (8)

R_lambda=A_2(v,lambda)p^2+A_1(v)p+A_0(v),                (9)
```

where

```text
A_2
 =-82458112[
     131769lambda v^2-15246lambda v+441lambda
     +66550v^2-7700v+150
   ],

A_1
 =-464640(10629366v^3-1156639v^2+27104v-315),

A_0
 =-63(155624547606v^5+3215383215v^4-1700698560v^3
      +58124770v^2-855470v+2583).                       (10)
```

The coefficient `A_2` is not identically zero for any `lambda`, and
`gcd(A_1,A_0)=1`. Hence (9) is always genuinely quadratic and has no
vertical identity fibre.

## 2. The completed-square family

The quadratic discriminant is

```text
Disc_p(R_lambda)
 =-809588736(121v-7)^2 H_lambda(v),                      (11)
```

where

```text
H_lambda
 =(35949270496986lambda+18156197220700)v^5
 +(742753522665lambda-1682717215850)v^4
 +(-392861367360lambda-8503492800)v^3
 +(13426821870lambda+370417300)v^2
 +(-197613570lambda+3642100)v
 +(596673lambda-337050).                                (12)
```

Choose `c in C*` with `c^2=-809588736` and set

```text
r=(2A_2p+A_1)/[c(121v-7)].                              (13)
```

On (9), the completed-square identity gives

```text
r^2=H_lambda(v).                                        (14)
```

The only question is whether a special ratio can lower the normalization
to genus zero.

## 3. Complete exceptional-ratio classifier

The exact parameter discriminant is

```text
Disc_v(H_lambda)
 =-2^38 3^7 5^12 7^2 11^40
   (231lambda+100) P_7(lambda),                         (15)
```

where

```text
P_7
 =5385425879119671lambda^7
  +10493769091223880lambda^6
  +8643047959054740lambda^5
  +3838851484116880lambda^4
  +956820534904080lambda^3
  +120517531825152lambda^2
  +3892886965440lambda
  -418060512000.                                        (16)
```

The septic is squarefree, and

```text
P_7(-100/231)=-4185216000,

P_7(-50/99)=-294697041920000/2187.                      (17)
```

There are four exhaustive cases.

### 3.1 Generic ratio

If `lambda` is neither `-50/99` nor a root of (15), then
`H_lambda` is a squarefree quintic. The smooth projective model of (14)
has genus two.

### 3.2 The degree-drop ratio

At

```text
lambda=-50/99,
```

the quintic leading coefficient vanishes, but

```text
H_(-50/99)
 =-1600/3 [
    3858459858v^4-356083761v^3+12020261v^2
    -193963v+1197
   ].                                                   (18)
```

The quartic in brackets is squarefree, so the normalization has genus one.

### 3.3 The rational double-root ratio

At the linear factor in (15),

```text
lambda=-100/231,

H_(-100/231)
 =50(121v-3)^2
   (3543122v^3-2562175v^2+91476v-1323).                 (19)
```

The cubic is squarefree and has value `-576` at `v=3/121`. Dividing `r`
by `121v-3` normalizes (14) to a squarefree cubic double cover, again of
genus one.

### 3.4 The seven algebraic double-root ratios

Let `lambda` be a root of `P_7`. Equations (15)--(17) show that this is a
simple zero of the parameter discriminant, distinct from both ratios above,
and the fibre still has degree five. The singular locus of the quintic
discriminant consists of polynomials with a triple root or at least two
double roots. A simple transverse zero cannot lie there. Consequently

```text
H_lambda=(v-alpha)^2 G_3(v)
```

with `G_3` squarefree and `G_3(alpha)!=0`. Dividing `r` by `v-alpha`
normalizes (14) to `rho^2=G_3(v)`, a genus-one curve.

Thus every value of the weighted ratio has normalization genus at least one.

## 4. Constant map and the y=0 boundary

If `v` were nonconstant, (13)--(14) would give a nonconstant morphism from
`P^1` to the positive-genus normalization classified above, impossible by
Riemann--Hurwitz. Hence `v` is constant. Since (9) has no vertical identity
fibre, `p` is algebraic over the algebraically closed constant field and is
therefore constant. Since `C!=0`, equation (4) makes `y`, then `u,zeta,Z,T`
and `q`, constant. This contradicts the genuine deck, which fixes the
constant field and sends `q` to `-q`.

It remains to treat `y=0`. The open chart gives `u!=0`. The first flux
reconstructs

```text
Z=640C/99,
```

and the second becomes

```text
12800C^2+25344W+22869u^3=0.                             (20)
```

Thus `u,Z,T,q` are constants, giving the same deck contradiction. This
proves (1).

## 5. Scope and structural lesson

This theorem closes the entire `C,W` coefficient plane in
`mathcal A!=0`, including its two axes. It is the first support-two plane,
not a classification of the other nine coordinate planes or general mixed
strata. Branches outside the inherited reduction, split/even short edges,
and integral order raising remain open. Nothing here proves degree twenty
two, `JC(2)`, or `DC(2)`.

The structural object is the weighted projective ratio, not the raw pair of
coefficients. The family is stratified by fibre type:

```text
generic quintic genus 2
  -> eight double-root ratios of genus 1
  -> one degree-drop quartic of genus 1,
```

but never reaches genus zero. This is a precise instance of retaining the
global sidecar after a local quotient: the parameter discriminant classifies
every exceptional fibre, while the degree-drop ratio must be checked
separately because it is invisible to the fixed-degree discriminant.

## 6. Exact verification

Run

```bash
python3 04-computation/jc2_degree22_cw_plane_hyperelliptic_family_thm2429.py
python3 -O 04-computation/jc2_degree22_cw_plane_hyperelliptic_family_thm2429.py
```

The companion checks the weighted reductions (5)--(6), exact resultant
(8)--(10), completed-square family (11)--(12), complete parameter
discriminant (15)--(17), the rational exceptional factorization (19), the
degree-drop quartic (18), all squarefreeness and separation controls, and
the `y=0` equation (20). All truth-bearing checks use explicit exceptions
and remain active under optimized Python.

Normal, optimized, and stored transcripts byte-match after LF
normalization. The declared hashes are over the working-tree bytes.
