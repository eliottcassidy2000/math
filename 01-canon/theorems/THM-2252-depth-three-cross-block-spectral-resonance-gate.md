---
id: THM-2252
title: "Depth-three cross-block spectral resonance gate"
status: >
  PROVED + VERIFIED-EXACT + HOSTILE-AUDITED. Every scalar (3,4,5)
  survivor in the five-unit/three-blocker branch has a nonzero bounded
  relation aH+b_1u_1+13b_2u_2+169b_3u_3=0, where |a|<=15 and
  |b_j|<=33320. Each b_j retains a two-digit representation
  b_j=169n_j+m_j with |n_j|,|m_j|<=196. The proof groups each core's
  unavoidable two-step resonance before applying cross-block Fourier
  orthogonality. In the absence of a relation between the four aggregate
  blocks, positive Jackson smoothing and the exact two-step danger law give
  carrier mass below 433/3125 < 961/6930. This places the last depth-three
  profile on a finite union of explicit relation-and-carry planes; it does
  not exclude those planes or prove LRC(14).
source: codex-2026-07-25-depth-three-cross-block-resonance
depends_on:
  - THM-2145-two-block-spectral-crossing-and-6-plus-7-carry
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
  - THM-2222-scalar-transfer-parity-tower-and-four-checkpoint-survivor-reduction
related:
  - THM-2163-radix-relation-carry-descent
  - THM-2223-four-checkpoint-selberg-relation-and-carry-gate
  - THM-2239-unrestricted-multicore-signed-dual-profile-exclusion
  - THM-2243-composite-union-transfer-dual-depth-three-profile-exclusion
  - THM-2250-depth-three-pair-incidence-partition-reduction
script: 04-computation/lrc14_depth_three_cross_block_resonance_gate_thm2252.py
output: 05-knowledge/results/lrc14_depth_three_cross_block_resonance_gate_thm2252.out
script_sha256: 813c5e9c79f634dfcf30f4ee5fa146f1cfe36bfc17b15912f5f60c15afdc7192
output_sha256: d062cfd3b45a4262e79fa30b379825504a3e435e420bdf47b09f3619f1abc068
hash_basis: working-tree bytes (LF)
---

# THM-2252 -- depth-three survival forces cross-block resonance

The arbitrary-coupling Bellman bounds in THM-2239 and THM-2243 deliberately
forget that all current bits use one root digit. This theorem attacks that
loss from the Fourier side. It first retains each normalized core's exact
two-step self-incidence, then asks whether the four resulting blocks can
correlate. If not, the positive carrier is too small.

## 1. The hostile positive carrier

On `T=R/Z`, put

```text
D_b={x:||bx||<1/14},       C_H={x:||Hx||>1/7}.       (1)
```

Work in THM-2198's scalar five-unit/three-blocker branch with profile

```text
(lambda_1,lambda_2,lambda_3)=(3,4,5),

c_j=13^(lambda_j)u_j,       13 does not divide H u_1u_2u_3. (2)
```

For the signed five-unit residual `R`, THM-2198 gives

```text
p:=integral R_+ >=delta_5:=961/6930.                 (3)
```

The even-checkpoint support inclusions in THM-2222 give

```text
support(R_+) subset P:=C_H intersection U_0 intersection U_2, (4)

U_0=D_(13^3u_1) union D_(13^4u_2) union D_(13^5u_3),
U_2=D_(13u_1)   union D_(13^2u_2) union D_(13^3u_3). (5)
```

Define the three lower checkpoint frequencies

```text
v_1=13u_1,       v_2=13^2u_2,       v_3=13^3u_3.    (6)
```

Then

```text
U_2=union_j D_(v_j),       U_0=union_j D_(169v_j).  (7)
```

Thus a scalar survivor must satisfy

```text
measure(P)>=p>=delta_5.                              (8)
```

## 2. Group before taking Fourier support

Use THM-2145's normalized squared-Fejer kernel

```text
J_N=F_N^2/integral F_N^2
```

with `N=99`, and let

```text
K=2N-2=196.                                         (9)
```

For a circular interval `I`, write

```text
q_I=J_N*1_I.
```

The positive-kernel construction gives

```text
0<=q_I<=1,
integral q_I=measure(I),
FourierSupport(q_I) subset [-K,K],
||q_I-1_I||_1<=eta_N.                               (10)
```

THM-2145's exact integral Jackson coefficients and `pi<355/113` give the
strict rational estimate

```text
eta_N<2129/500000.                                  (11)
```

Apply (10) to the guard interval `C={x:||x||>1/7}` and the danger interval
`D={x:||x||<1/14}`. Put

```text
g(x)=q_C(Hx),

A_j(x)=1-q_D(169v_jx),
B_j(x)=1-q_D(v_jx),                                 (12)

Q(x)=g(x)(1-product_j A_j(x))(1-product_j B_j(x)).
```

The map in the last line is `[0,1]`-valued and has coordinate Lipschitz
constant one in each of its seven inputs. Haar invariance and a coordinate
telescope therefore give

```text
|measure(P)-integral Q|<=7eta_N.                    (13)
```

The grouping in (12) is load-bearing. The pair `A_jB_j` has Fourier
frequencies

```text
b_jv_j,       b_j=169n_j+m_j,
|n_j|,|m_j|<=K.                                     (14)
```

Its zero-frequency coefficient includes the genuine internal relation

```text
169*(+1)+(-169)=0.                                  (15)
```

Since `K>=169`, treating the two checkpoints as separate Fourier blocks
would always detect (15) and produce a vacuous "resonance." We instead
absorb every such internal relation into the constant coefficient of the
single core block `A_jB_j`.

## 3. Cross-block factorization

Call a tuple `(a,b_1,b_2,b_3)` an admissible cross-block frequency if

```text
|a|<=K,
b_j=169n_j+m_j,       |n_j|,|m_j|<=K.               (16)
```

Suppose there is no nonzero admissible tuple satisfying

```text
aH+b_1v_1+b_2v_2+b_3v_3=0.                          (17)
```

Expand `Q` as

```text
Q
 =g
  -g product_j A_j
  -g product_j B_j
  +g product_j(A_jB_j).                             (18)
```

Every Fourier mode in a term of (18) has the form on the left of (17).
Under the supposition, a zero total frequency requires zero aggregate
frequency in every block. Hence Fourier orthogonality factors the constant
coefficient block by block, while retaining the full internal constant
coefficient of each `A_jB_j`.

The one-block means are

```text
integral g=5/7,
integral A_j=integral B_j=6/7.                      (19)
```

Integer dilation preserves Haar measure, so the common paired mean is

```text
r_N:=integral
 (1-q_D(169x))(1-q_D(x)).                           (20)
```

Equations (18)--(20) give the exact factorized smooth value

```text
integral Q
 =(5/7)(1-2(6/7)^3+r_N^3).                         (21)
```

This is not an independence assumption about the two appearances of one
core: their complete Jackson-smoothed correlation is the number `r_N`.
Only correlations between distinct aggregate blocks have been excluded by
(17).

## 4. The exact two-step baseline

For the unsmoothed danger process, THM-2222's backward root law is

```text
P(X_current=1 | X_future=0)=2/13,
P(X_current=1 | X_future=1)=1/13.                   (22)
```

Its stationary distribution is `(6/7,1/7)` on `(0,1)`. Squaring the
transition matrix gives

```text
P(X_0=0 | X_2=0)=145/169,

r_0:=P(X_0=0,X_2=0)
    =(6/7)(145/169)
    =870/1183.                                      (23)
```

The product telescope and (10) imply

```text
|r_N-r_0|<=2eta_N.                                  (24)
```

Combining (13), (21), and (24), and using monotonicity of the cubic on
`[0,infinity)`, gives

```text
measure(P)
 <=7eta_N
   +(5/7)(1-2(6/7)^3+(r_0+2eta_N)^3).               (25)
```

At the rational cap in (11), the right side of (25) is

```text
5017879304282732213894543
------------------------------------
36216151278125000000000000

 =0.138553632210874...                              (26)

 <433/3125
 <961/6930,

961/6930-433/3125=487/4331250>0.                    (27)
```

This contradicts (8). Therefore (17) has a nonzero admissible solution.

For perspective, the completely cross-block-factorized unsmoothed value is

```text
(5/7)(1-2(6/7)^3+(870/1183)^3)
 =1144584995/11589168409
 =0.098763341303344....                             (28)
```

The large gap between (28) and `delta_5` is what pays the explicit Jackson
error without discarding the exact internal two-step correlation.

## 5. The triangular relation-and-carry plane

Substitute (6) into the relation just forced:

```text
aH+13b_1u_1+13^2b_2u_2+13^3b_3u_3=0.               (29)
```

Since `H` is a thirteen-unit, reduction modulo thirteen gives `13|a`.
Write `a=13a_0` and divide (29) by thirteen:

```text
a_0H+b_1u_1+13b_2u_2+169b_3u_3=0,                  (30)

|a_0|<=15,
|b_j|<=170K=33320,
b_j=169n_j+m_j,
|n_j|,|m_j|<=196.                                   (31)
```

The aggregate vector `(a_0,b_1,b_2,b_3)` is nonzero. Because all four
weighted speeds in (30) are positive, at least two aggregate coefficients
are nonzero and both signs occur.

Equation (30) is the promised cross-block resonance gate. It is stronger
than a generic coefficient-height relation: the powers `1,13,169` retain
the valuation profile, and (31) retains the two checkpoint digits whose
internal carry was grouped in Section 2. There are only finitely many such
coefficient packets, so `(3,4,5)` is reduced to a finite union of explicit
integer hyperplanes in `(H,u_1,u_2,u_3)`.

## 6. Failure boundary and typed handoff

The exact referee also checks that `N=98`, degree `194`, does not cross the
target with the same equal-degree Jackson/L1 ledger. Thus `N=99` is
certificate-minimal for this proof budget, not optimal among all kernels or
all resonance theorems.

The theorem's typed connection is

```text
source:
  the (3,4,5) positive residual carrier C_H intersection U_0 intersection U_2;

target:
  a finite packet of triangular relation-and-carry planes (30)--(31);

map:
  group each core's two-step orbit, smooth by a positive Jackson kernel,
  and factor Fourier constant terms between aggregate blocks;

preserved:
  the guard mass, exact profile scales, exact two-step danger
  autocorrelation, both checkpoint digits, and coefficient signs;

destroyed:
  Fourier phase inside a nonzero relation plane, the five terminal unit
  dangers, the negative residual carrier, owners, and exact common-root
  incidence;

next sidecar:
  combine the finite relation packet with THM-2250's pair-incidence
  partition or an exact relation-plane Bellman calculation.              (32)
```

In particular, (30) does not say that `a_0` is nonzero; a survivor may have
a bounded relation among the three core blocks alone. It also does not
make the relation planes finite sets of cores. Those are the precise
remaining obligations.

There is a decisive hostile control when all normalized danger combs
coincide. If `u_1=u_2=u_3=u`, then

```text
(a_0,b_1,b_2,b_3)=(0,13,-1,0)                       (33)
```

satisfies (30) identically. Thus THM-2250's reduction to a common danger
comb cannot be combined with this theorem to infer a finite list of ratios
`H/u`: one must force a relation outside the diagonal syzygy module, or
close the common-core branch by a different signed/incidence argument.
This is exactly the relative rank-harvesting boundary, not a bookkeeping
technicality.

## 7. Exact audit

Run

```bash
python3 04-computation/lrc14_depth_three_cross_block_resonance_gate_thm2252.py
python3 -O 04-computation/lrc14_depth_three_cross_block_resonance_gate_thm2252.py
```

Both modes reproduce the stored transcript byte for byte. The companion
checks the exact Jackson coefficient formula and rational-`pi` error cap,
the two-step transition matrix and `870/1183` correlation, (26)--(28), the
degree-194 hostile boundary, the internal-relation control (15), all
coefficient heights, the diagonal hostile syzygy (33), and a digest of the
load-bearing constants.

This theorem reduces the last depth-three scalar profile to bounded
resonance planes. It does not close those planes, the 165 first-depth-one
profiles, the non-scalar branch, or LRC(14). QED.
