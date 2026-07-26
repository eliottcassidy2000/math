---
id: THM-2425
title: "Degree-twenty-two C, D, and E axis hyperelliptic closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  In the open first-flux chart of the genuine nonsplit polynomial
  exact-square-prefix degree-twenty-two branch, each of the three
  one-sparse coefficient axes C, D, and E is empty. After the natural
  weighted quotient, first-flux reconstruction makes the second flux
  quadratic in the axis coordinate. Completing the square removes only
  the excluded first-flux wall and gives squarefree hyperelliptic
  residuals of degrees 5, 4, and 3, hence genera 2, 1, and 1. A rational
  trajectory cannot map nontrivially to any of them. This closes the C,
  D, and E axes, not the B axis, mixed strata, degree twenty two, JC(2),
  or DC(2).
source: klein-2026-07-26-degree-twenty-two-cde-axis-triage
depends_on:
  - THM-2129-quartic-faber-three-coefficient-boundary-classification
  - THM-2214-nonsplit-terminal-quartic-spectral-curve-closure-through-degree-ten
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
  - THM-2411-degree-twenty-two-first-flux-pole-divisor-square-class-reduction
related:
  - THM-2247-nonsplit-terminal-quartic-degree-fourteen-closure
  - THM-2406-degree-eighteen-H4-weighted-pole-deep-wall-collapse
  - THM-2423-degree-twenty-two-W-axis-genus-two-and-origin-cusp-closure
  - THM-2428-degree-twenty-two-B-axis-trigonal-ramification-closure
script: 04-computation/jc2_degree22_cde_axis_hyperelliptic_thm2425.py
output: 05-knowledge/results/jc2_degree22_cde_axis_hyperelliptic_thm2425.out
script_sha256: 12d3e8a98ffa8498a92715bf72d087f2be54758e1d83ee85d43155e466f70d63
output_sha256: a29e8d0401f1aed45ae6845c77d3ae1bc0c61e3d5d9d77416be7cc5c41a6c59e
hash_basis: working-tree bytes (LF)
---

# THM-2425 -- the degree-twenty-two C, D, and E axes are empty

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2411 closes the divisor on which the first Faber flux loses its `Z`
coefficient, and THM-2423 closes the complete `W`-axis in the complementary
chart. This theorem closes three more one-sparse axes. Its exact conclusion is

```text
genuine degree-22 trajectory,
mathcal A!=0,

(B,D,E,W)=(0,0,0,0), C!=0
or
(B,C,E,W)=(0,0,0,0), D!=0
or
(B,C,D,W)=(0,0,0,0), E!=0

    => contradiction.                                           (1)
```

The proof uses only the two hostile-audited fluxes of THM-2411 and the
constant field. No third-flux or polynomial-sidecar assertion is needed.

## 1. Weighted first-flux quotient

Retain the normalized coordinates

```text
y=11s,             u=dT,             Z=T^2,

wt(y,u,Z,B,C,D,E,W)=(1,2,3,2,3,4,5,6)                (2)
```

and the equations `N_1=N_2=0` from THM-2411, equations (12)--(16).
On each axis in (1),

```text
mathcal A=9y^2(7-121u/y^2).                            (3)
```

First suppose `y!=0` in `C(x)` and put

```text
v=u/y^2,                 zeta=Z/y^3,

p_C=C/y^3,               p_D=D/y^4,               p_E=E/y^5.
                                                               (4)
```

Write

```text
k_0=922383v^2-25410v+63,

g_0
 =15944049zeta^2-162339408zeta v+2236080zeta
  -1190488992v^3+147581280v^2-1219680v+672.            (5)
```

After division by `y^5` and `y^6`, respectively, the two fluxes on an axis
`X in {C,D,E}` are

```text
f_X=11979(7-121v)zeta+4k_X=0,
g_X=0,                                                     (6)
```

where

```text
k_C=k_0+(2342560v-58080)p_C,
g_C=g_0+(-206145280zeta+449771520v-1239040)p_C;

k_D=k_0+511104p_D,
g_D=g_0+(-1978994688v+16355328)p_D;

k_E=k_0-3748096p_E,
g_E=g_0-239878144p_E.                                    (7)
```

The open chart `mathcal A!=0` says

```text
121v-7!=0.                                                (8)
```

Thus (6) reconstructs `zeta` uniquely. Eliminating it introduces no
component and gives one quadratic

```text
R_X(v,p_X)=A_(2,X)(v)p_X^2+A_(1,X)(v)p_X+A_(0,X)(v)=0.   (9)
```

The exact resultants are nonzero constant multiples of `R_C,R_D,R_E`,
with multipliers

```text
255104784,             2295943056,             255104784. (10)
```

## 2. The three exact quadratics

For the `C`-axis,

```text
A_(2,C)
 =-5487587353600v^2+634927462400v-12368716800,

A_(1,C)
 =-4938828618240v^3+537420744960v^2
   -12593602560v+146361600,

A_(0,C)
 =-9804346499178v^5-202569142545v^4+107144009280v^3
   -3661860510v^2+53894610v-162729.                     (11)
```

For the `D`-axis,

```text
A_(2,D)=29025255424,

A_(1,D)
 =-1810903826688v^3+119729178624v^2
   +4329050880v-109716992,

A_(0,D)
 =-1089371833242v^5-22507682505v^4+11904889920v^3
   -406873390v^2+5988290v-18081.                        (12)
```

For the `E`-axis,

```text
A_(2,E)=14048223625216,

A_(1,E)=4938828618240v^2-571434716160v+3935500800,

A_(0,E)
 =-9804346499178v^5-202569142545v^4+107144009280v^3
   -3661860510v^2+53894610v-162729.                     (13)
```

In every row,

```text
gcd(A_(2,X),A_(1,X),A_(0,X))=1.                         (14)
```

Consequently no constant value of `v` makes (9) an identity in `p_X`.

## 3. Completing the square exposes genera two, one, and one

The three quadratic discriminants factor as

```text
Disc_p(R_C)
 =-40479436800(121v-7)^2 H_C(v),

Disc_p(R_D)
 = 116101021696(121v-7)^2 H_D(v),

Disc_p(R_E)
 = 1011472101015552(121v-7)^2 H_E(v),                   (15)
```

where

```text
H_C
 =363123944414v^5-33654344317v^4-170069856v^3
   +7408346v^2+72842v-6741,

H_D
 =1929229929v^4+42517464v^3-790614v^2-203280v+2485,

H_E
 =37202781v^3+6720219v^2-134673v+497.                   (16)
```

Their exact discriminants are

```text
Disc(H_C)
 = 2^40 3^9 5 7^2 11^40 23 63113,

Disc(H_D)
 =-2^27 3^6 5^2 7^3 11^24 17,

Disc(H_E)
 = 2^12 3^4 5^4 7 11^12 59.                            (17)
```

They are therefore squarefree. Choose `c_X in C*` whose square is the
corresponding constant in (15). On (9), define

```text
r_X=(2A_(2,X)(v)p_X+A_(1,X)(v))
       /[c_X(121v-7)].                                  (18)
```

Equations (8), (9), and the elementary completed-square identity give

```text
r_C^2=H_C(v),             r_D^2=H_D(v),
                           r_E^2=H_E(v).                 (19)
```

The smooth projective models in (19) have genera `2,1,1`, respectively.
If `v` were nonconstant, the rational functions `(v,r_X)` would give a
nonconstant morphism from `P^1` to a positive-genus curve, impossible by
Riemann--Hurwitz.

Hence `v` is constant. By (14), equation (9) is then a nonzero polynomial
of degree at most two in `p_X`; thus `p_X`, algebraic over the algebraically
closed constant field, is constant. Since `X in C*`, equation (4) makes
`y`, then `u,zeta,Z,T`, constant. The genuine deck fixes the constant field
but sends `q` to `-q`, contradicting `q^2=T!=0`.

## 4. The boundary y=0

It remains to justify the division in (4).

On the `C`-axis, `mathcal A!=0` and `y=0` give `u!=0`. The first flux gives

```text
Z=640C/99,
```

and the second becomes

```text
12800C^2+22869u^3=0.                                   (20)
```

Thus `u,Z,T,q` are constants, contradicting the genuine deck.

On the `D`-axis, the first flux at `y=0` is a nonzero constant multiple of
`uZ`. Since `u!=0`, it forces `Z=0`, contradicting `T!=0`.

On the `E`-axis, the first flux gives

```text
Z=-1024E/(99u),
```

and the second becomes

```text
22869u^5-32768E^2=0.                                   (21)
```

Again every coordinate is constant and the deck is contradicted. This
finishes all three axes in (1).

## 5. Scope and the ranked surviving target

This theorem closes exactly the `C`, `D`, and `E` one-sparse axes in the
open first-flux chart. Together with THM-2423, four of its five weighted
coordinate axes are empty. At the time of this proof the `B`-axis remained:
its analogous equation is cubic in `p_B=B/y^2`, with discriminant

```text
constant
  (131769v^2-20570v+189)^2
  K_9(v),                                               (22)
```

where `K_9` is squarefree of degree nine. This is a trigonal, not
hyperelliptic, normalization problem. Candidate THM-2428 subsequently uses
its nine simple finite branch places and Riemann--Hurwitz to close that last
axis without locating every infinity place. Mixed coefficient strata,
branches outside the inherited reduction, split/even short edges, and
integral order raising remain open. Nothing here proves degree twenty two,
`JC(2)`, or `DC(2)`.

The useful cross-thread operation is again quotient plus retained sidecar:
the weighted quotient preserves the full first-flux reconstruction
coefficient `(121v-7)`; completing the square is legal precisely because
that coefficient is retained and excluded by (8). No tournament, knot,
additive-graph, or carry object is identified with the hyperelliptic curves.

## 6. Exact verification

Run

```bash
python3 04-computation/jc2_degree22_cde_axis_hyperelliptic_thm2425.py
python3 -O 04-computation/jc2_degree22_cde_axis_hyperelliptic_thm2425.py
```

The companion checks the weighted reductions (5)--(7), exact resultants
(9)--(13), coefficient gcds (14), all three completed-square factorizations,
the residual discriminants and excluded-wall controls, and the three
`y=0` hostile boundaries. Every truth-bearing check uses explicit
exceptions and remains active under optimized Python.

Normal, optimized, and stored transcripts byte-match after LF
normalization. The declared hashes are over the working-tree bytes.

## 7. Independent hostile audit

An independent audit reran the normal and optimized companions and
byte-compared both with the stored transcript; the declared source and output
hashes match. It independently substituted the three weighted axes into the
hostile-audited THM-2411 fluxes, reconstructed all three resultants, and
rederived

```text
Disc_p(R_X)=constant (121v-7)^2 H_X(v)
```

with residual degrees `5,4,3`. Separate Euclidean gcds confirmed that every
`H_X` is squarefree and that every coefficient triple has gcd one, so no
vertical identity fibre was lost.

The audit also attacked the two proof boundaries. At `y=0`, it recovered
either `Z=0` or the displayed constant-field equations (20)--(21). On
`y!=0`, it checked that completing the square divides only by the
first-flux coefficient already excluded in (8), and that constant `v`
cannot leave a varying `p_X`. No algebraic, geometric, deck, scope, or
reproducibility defect remains. **QED.**
