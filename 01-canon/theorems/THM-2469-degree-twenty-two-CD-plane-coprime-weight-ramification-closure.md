---
id: THM-2469
title: "Degree-twenty-two C-D plane coprime-weight ramification closure"
status: >
  PROVED + VERIFIED-EXACT. In the open first-flux chart of the genuine
  nonsplit polynomial exact-square-prefix degree-twenty-two branch, the
  complete coefficient plane B=E=W=0 is empty. For C,D nonzero, the
  root-free invariants lambda=D^3/C^4 and r=D/(Cy) give a 19-term
  plane curve of bidegree (8,5). Its Newton polygon reduces every
  factorization to one linear or quadratic v-factor; both exact
  coefficient ideals force lambda=0, so every physical fibre is
  absolutely irreducible. The v-discriminant is
  constant*lambda^4*W_4^2*K_30, where W_4 is exactly the excluded
  first-flux wall. A complete exceptional-parameter factorization,
  together with a Sylvester-valuation bound, leaves at least twenty
  simple finite branch values away from the wall in every fibre.
  Riemann--Hurwitz gives normalization genus at least six. No rational
  trajectory survives. Together with the axes, this closes the fifth
  of ten support-two planes. It does not close the other five planes,
  higher mixed strata, JC(2), or DC(2).
source: codex-2026-07-27-degree-twenty-two-CD-ramification
depends_on:
  - THM-2411-degree-twenty-two-first-flux-pole-divisor-square-class-reduction
  - THM-2425-degree-twenty-two-CDE-axis-hyperelliptic-closure
related:
  - THM-2429-degree-twenty-two-CW-plane-hyperelliptic-family-closure
  - THM-2437-degree-twenty-two-DW-plane-quartic-ramification-closure
  - THM-2463-degree-twenty-two-BD-plane-square-lift-closure
  - THM-2468-degree-twenty-two-BW-plane-square-lift-closure
script: 04-computation/jc2_degree22_cd_plane_ramification_thm2469.py
output: 05-knowledge/results/jc2_degree22_cd_plane_ramification_thm2469.out
script_sha256: 8d38b9f43776076b76cd8cd799294633a97b57d30397f0dbeb2035d8e39c4298
output_sha256: 1bc5eb1e296dc9b2b8047d69a8cbf8c3392809f58d8eb26be13d690cabda2721
hash_basis: working-tree bytes (LF)
---

# THM-2469 -- the degree-twenty-two C-D plane is empty

**PROVED + VERIFIED-EXACT.**

The four previously closed mixed planes are exactly the four pairs among
weights `(2,3,4,5,6)` having nontrivial gcd. On those planes a retained power
class can create a branched cover. The six pairs left after THM-2468 have
coprime weights, so that operation has degree one and supplies no genus.

The `C,D` plane is nevertheless rigid. Bezout arithmetic for weights three
and four gives a root-free quotient, and the base curve itself has uniformly
large ramification. The exact conclusion is

```text
genuine degree-22 trajectory,
mathcal A!=0,
B=E=W=0
    => contradiction.                                           (1)
```

Thus the complete `C,D` support-two coefficient plane is empty.

## 1. Coprime-weight invariant coordinates

Use the target-translated coordinates of THM-2411:

```text
y=11s,                   u=dT,                   Z=T^2,

wt(y,u,Z,B,C,D,E,W)=(1,2,3,2,3,4,5,6).              (2)
```

The cases `C=0` and `D=0` are the closed `D` and `C` axes of THM-2425.
Suppose from now on that

```text
C,D!=0.                                                (3)
```

First take `y!=0` in `C(x)` and define the constant weighted ratio and the
moving root-free coordinate

```text
lambda=D^3/C^4 in C*,             r=D/(Cy).           (4)
```

These satisfy the lossless identities

```text
C/y^3=r^3/lambda,                  D/y^4=r^4/lambda.  (5)
```

As usual put

```text
v=u/y^2,                           zeta=Z/y^3.         (6)
```

Multiply the first two normalized fluxes of THM-2411 by `lambda`. They become

```text
F_1
 =1331lambda(63-1089v)zeta
  +4[(2342560v-58080)r^3+511104r^4
      +lambda(922383v^2-25410v+63)]
 =0,                                                        (7)

F_2
 =lambda[15944049zeta^2+(-162339408v+2236080)zeta
          -1190488992v^3+147581280v^2-1219680v+672]
  -206145280r^3zeta+(449771520v-1239040)r^3
  +(-1978994688v+16355328)r^4
 =0.                                                        (8)
```

The open first-flux chart is

```text
63-1089v!=0,                 equivalently v!=7/121.   (9)
```

Thus (7) reconstructs `zeta` uniquely.

## 2. The 19-term eliminant

Exact elimination gives

```text
Res_zeta(F_1,F_2)=255104784lambda P_lambda(r,v),       (10)
```

where

```text
P_lambda
 =261227298816r^8+79159787520r^7

  +(-5487587353600v^2+634927462400v-12368716800)r^6

  +lambda(-16298134440192v^3+1077562607616v^2
           +38961457920v-987452928)r^4

  +lambda(-4938828618240v^3+537420744960v^2
           -12593602560v+146361600)r^3

  -63lambda^2 L_5(v),                                    (11)
```

and

```text
L_5(v)
 =155624547606v^5+3215383215v^4-1700698560v^3
   +58124770v^2-855470v+2583.                         (12)
```

The old fixed section survives:

```text
P_lambda(0,v)=-63lambda^2L_5(v).                       (13)
```

Unlike THM-2463/2468, no nontrivial power cover lies above it because the
active weights are coprime. The proof must extract genus from (11) itself.

## 3. Uniform absolute irreducibility

For every `lambda!=0`, the polynomial `P_lambda` is absolutely irreducible
in `C[r,v]`.

Its Newton polygon is

```text
N=conv{(0,0),(8,0),(6,2),(0,5)}.                      (14)
```

The four primitive edge directions have lengths `(8,2,3,5)`. A Minkowski
summand has edge lengths

```text
(a,b,c,d)=(b+2c,b,c,b+c),

0<=b<=2,                    0<=c<=3.                  (15)
```

The coefficient of `v^5` in (11) is the nonzero constant
`-9804346499178lambda^2`. Hence any factorization descends by Gauss to monic
factors in `C[r,v]`; no nonconstant factor can have `v`-degree zero. Since
the two positive `v`-degrees add to five, one factor has `v`-degree at most
two.

For degree one, (15) leaves `(b,c)=(1,0)` or `(0,1)`. Both supports lie in
the single maximal ansatz

```text
v+ar^2+br+c.                                           (16)
```

Substitution in (11), followed by collection in `r`, gives the reduced
coefficient ideal

```text
1331c^2+154c+3,
a,
55b-363c-21,
lambda.                                                (17)
```

Thus a linear factor forces `lambda=0`.

For degree two, the only pairs in (15) are `(2,0),(1,1),(0,2)`. Every such
support lies in

```text
v^2+(ar^2+br+c)v+dr^4+er^3+fr^2+gr+h.                (18)
```

The exact remainder ideal for (18) has reduced basis

```text
a, b, 121c+14, d, e, 3025f+144,
6655g+96, 1331h-3, lambda.                            (19)
```

Again divisibility forces `lambda=0`. The factor degrees are exhausted, so
(11) is absolutely irreducible for every physical ratio. Let `C_lambda` be
its smooth projective normalization. Projection to the `r`-line has degree
five.

## 4. Exact ramification divisor

Define the quartic

```text
W_4(r,lambda)=31944r^4+4840r^3+105lambda.             (20)
```

The squared factor below is not branch mass. Direct substitution gives

```text
P_lambda(r,7/121)=256W_4(r,lambda)^2,                 (21)
```

so its roots are exactly intersections with the excluded first-flux wall.

The exact quintic discriminant factors as

```text
Disc_v(P_lambda)
 =c lambda^4 W_4(r,lambda)^2 K_30(r,lambda),          (22)
```

for `c in Q*`. Normalize `K_30` to be primitive with

```text
K_30(0,lambda)=655809304130859375lambda^10.           (23)
```

It has degree thirty in `r`, degree ten in `lambda`, and 46 nonzero terms.
Its leading coefficient is

```text
[r^30]K_30
 =176357374671601011051986944000000
    (9504lambda+4375).                                (24)
```

Put

```text
A_1=7392lambda-625,

A_2=90112lambda-5625,

A_3=3372171264lambda^2-218592000lambda+383828125.     (25)
```

Exact factorization over `Q[lambda]` gives

```text
Disc_r(K_30)
 =c' lambda^290 A_1 A_2^2 A_3^3
      P_10^3 Q_10^2 R_12,                            (26)
```

where `P_10,Q_10,R_12` are the unique primitive positive-leading
irreducible factors of degrees `10,10,12` occurring with exponents `3,2,1`.
This characterization, (22)--(26), and the exact companion uniquely specify
all their coefficients; no numerical root approximation enters the proof.

The collision with the excluded wall is independently exact:

```text
Res_r(K_30,W_4)
 =c'' lambda^30 A_3^3 S_4(lambda),                    (27)
```

where

```text
S_4
 =6125963995157747283975979008lambda^4
  +7723423480030158582978912000lambda^3
  +4548030222783870066786609375lambda^2
  +1419376581390195945000000000lambda
  +17937798138268750000000000.                        (28)
```

The eight polynomials

```text
9504lambda+4375,
A_1,A_2,A_3,P_10,Q_10,R_12,S_4                       (29)
```

are squarefree and pairwise coprime. Thus degree drop, root collision, and
wall collision never overlap, except for the deliberately shared `A_3`
between (26) and (27).

At the sole degree-drop ratio

```text
lambda=-4375/9504,                                    (30)
```

`K_30` becomes a squarefree degree-29 polynomial and remains coprime to
`W_4`.

## 5. The Sylvester valuation bound

We need only a coarse but uniform consequence of the large factors in (26).

**Lemma.** Let `k_s(r)` be a degree-`n` polynomial over a characteristic-zero
DVR whose leading coefficient is a unit. If its discriminant has valuation
`e`, then

```text
deg gcd(k_0,partial_r k_0)<=e.                        (31)
```

Indeed, the resultant `Res(k_s,partial_r k_s)` is the determinant of the
Sylvester matrix. Its reduction has nullity equal to the gcd degree. In Smith
normal form, each lost rank contributes at least one uniformizer to the
determinant, proving (31).

Every nonzero exceptional factor in (26) is squarefree and appears with
exponent at most three. Coprimality with (24) makes the leading coefficient a
unit there. Therefore every degree-30 exceptional fibre has

```text
d=deg gcd(K_30,partial_rK_30)<=3.                     (32)
```

If the nonsimple roots have multiplicities `m_i>=2`, then

```text
d=sum_i(m_i-1),

sum_i m_i=d+#{i}<=2d.                                 (33)
```

Consequently at least

```text
30-2d>=24                                             (34)
```

roots of `K_30` are simple. At an `A_3` root, at most the four roots of
`W_4` must also be discarded. Hence even there at least twenty simple roots
remain away from the wall. At an `S_4` root, `K_30` is squarefree and again
at most four roots are lost. At (30) all 29 finite roots are simple and none
is on the wall. These cases and the generic fibre exhaust every
`lambda in C*` by (26)--(29). Thus uniformly

```text
at least 20 simple finite roots of K_30 lie off W_4.   (35)
```

## 6. Genus and trajectory closure

At a simple root of `K_30` away from `W_4`, equation (22) gives discriminant
valuation one for the degree-five `v`-polynomial. Its leading coefficient is
constant and nonzero. The order-to-normalization index changes discriminant
valuation by an even integer, so the order is already maximal there; tame
characteristic zero then gives one simple ramification point.

Let `g_lambda` be the genus of `C_lambda`. Riemann--Hurwitz and (35) give

```text
2g_lambda-2
 =5(-2)+sum_P(e_P-1)
 >=-10+20=10,

g_lambda>=6.                                           (36)
```

The rational pair `(r,v)` from (4),(6) would give a nonconstant map from
`P^1` to this positive-genus normalization, impossible. Hence `r,v` are
constant. Since `r=D/(Cy)` is nonzero, `y` and then `u` are constant.
Equation (7) reconstructs constant `zeta`, so `Z,T,q` are constants. This
contradicts the genuine deck, which fixes the constant field and sends the
nonzero `q` to `-q`.

It remains to treat `y=0`. The original first flux reduces to

```text
N_1=u(-1449459Z+9370240C).                            (37)
```

The open chart gives `u!=0`, so

```text
Z=9370240C/1449459 in C*.                             (38)
```

The second flux is then the constant-field cubic

```text
15944049Z^2-206145280CZ
 -1978994688Du-1190488992u^3=0.                      (39)
```

Thus `u`, and then `T,q`, are constants, giving the same deck contradiction.
Together with the two axes, this proves (1).

## 7. Scope and structural lesson

This theorem closes the fifth of ten support-two planes. The open pairs are

```text
(B,C), (B,E), (C,E), (D,E), (E,W).                    (40)
```

All have coprime weights. The retained-power-class operation that closed
THM-2463/2468 is therefore exhausted on the degree-twenty-two coefficient
axes. The replacement operation is the root-free Bezout coordinate: form a
weight-zero ratio and a weight-one moving coordinate, retain the fixed `L_5`
section, but extract genus from the base curve's exact ramification divisor.
Higher mixed strata, branches outside the inherited reduction, split/even
short edges, and integral order raising remain open. Nothing here proves
`JC(2)` or `DC(2)`.

## 8. Exact verification

Run

```bash
python3 04-computation/jc2_degree22_cd_plane_ramification_thm2469.py
python3 -O 04-computation/jc2_degree22_cd_plane_ramification_thm2469.py
```

The companion reconstructs (7)--(13), the Newton polygon and both exhaustive
factor ideals, (20)--(30), all squarefreeness and coprimality controls, the
degree-drop fibre, the branch/genus floors, and the `y=0` equations. Normal and
optimized transcripts byte-match the stored output. All truth-bearing checks
use explicit exceptions and remain active under optimized Python. **QED.**
