---
id: THM-3649
title: "Russell-cylinder Q-star actual lift and fourth-stable closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For the fixed
  degree-seven Q-star fold, a deterministic exact target-ring
  certificate gives J_0=1 and J_1=J_2=0.  Composing this nonvacuity witness
  with PROVED THM-3642's arbitrary-two-form identity forces
  lambda(J_4)=5440/81.  The fixed compiler, quadratic fold, and Q-star
  polynomial are the exact scope; no claim about the whole degree-seven
  zero-debt line, J_3 solvability, unrestricted Keller families, or JC(2)
  is made.
source: root/actual-zero-debt-lifts/2026-08-21
depends_on:
  - THM-3642-russell-cylinder-zero-debt-actual-lift-and-fourth-stable-closure
related:
  - THM-3641-russell-cylinder-principal-noneven-curvature-debt-boundary
script: 04-computation/jc2_russell_cylinder_qstar_actual_lift_fourth_stable_closure_thm3649.py
output: 05-knowledge/results/jc2_russell_cylinder_qstar_actual_lift_fourth_stable_closure_thm3649.out
script_sha256: ce396345bdf4dbaae68791f0de7ecbacbf3554aa29b49b168a4f443a2c2638db
output_sha256: edc93e3a3706b9e8f3ebc4ced7732a9613cb87cf0de0688aa5f2fe8e1ccf1620
hash_basis: raw LF bytes
audit: >
  PASS -- independent hostile reconstruction checked the target monomial,
  block, row, and hash conventions; representative-bound delta data; both
  modular prime controls; the selected 776-by-776 rational solve; the full
  J1=J2 residuals; and composition with THM-3642.  Normal, optimized, and
  stored transcripts are byte-identical at 826 active gates; hashes, AST0,
  docs, diff, LF, and whitespace checks pass.
---

# THM-3649 -- Russell-cylinder Q-star actual lift and fourth-stable closure

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
This theorem fills the fixed-`Q_*` nonvacuity gap left deliberately open in
THM-3642.  It constructs one actual target-ring lift through `J_2`, then uses
THM-3642's already-proved arbitrary-two-form identity to show that every such
equality-wall pair has nonzero fourth-stable debt.

All rings and points are over `C`; the deterministic certificate is rational.
The frozen certificate and its composition passed independent hostile audit
before promotion.

## 0. Fixed compiler, fold, and notation

Use

```text
R=C[b,c,e]/(c^2e-b(b+4)),

D=1+x^2q,
b=(D-1)(D+2)^2,       c=xD(D+2),       e=q(D+3),       (1)
```

and the quadratic fold

```text
(x,t) |-> (x,Q_*(x)+t^2,w=t),                           (2)
```

where

```text
Q_*=-x^7-(27/4)x^6+3x^5+18x^4-3x^3
       -(27/2)x^2+x-3/4.                                (3)
```

At `x=-1,0,1`, respectively,

```text
Q_*=(-3,-3/4,-3),
Q_*'=(-9/2,1,9/2),
Q_*''=(-27/2,-27,-27/2).                                (4)
```

Thus `Q_*` is an ordinary principal-chart control and lies on the zero
second-debt wall

```text
5Q_*''(-1)+13Q_*''(1)+243=0.                            (5)
```

For target coefficients in `R`, put

```text
F^#=c+wF_1+w^2F_2+w^3F_3,
G^#=e+wG_1+w^2G_2+w^3G_3.                              (6)
```

Write the pulled-back source Jacobian as

```text
Jac_(x,t)(F^#,G^#)=sum_(n>=0)t^nJ_n(x).                 (7)
```

The claim is that one deterministic choice in `(6)` satisfies the exact
polynomial identities

```text
J_0=1,                    J_1=0,                    J_2=0. (8)
```

No `J_3=0` claim is included.

## 1. Exact target-ring certificate for `J_0=1`

Let `gamma` denote restriction through `q=Q_*(x)` and let

```text
delta(K)=partial_q K(b(x,q),c(x,q),e(x,q))|_(q=Q_*(x)). (9)
```

The raw target monomials are `b^i c^j e^k`, enumerated in nested increasing
`(i,j,k)` order.  Their restriction-degree weights are

```text
(deg gamma(b),deg gamma(c),deg gamma(e))=(27,19,16).    (10)
```

The cutoffs below are chosen reconstruction cutoffs, not minimality claims.
At cutoff `160` there are `141` monomials.  In ascending `x`-coefficient
rows, the `J_0` system has shape `179 x 282` and columns

```text
[c' gamma(m)]_m,                  [-e' gamma(m)]_m.      (11)
```

The first block is the `G_1` coefficient vector and the second is the `F_1`
vector.  Exact rational RREF has rank `178`; setting every free variable to
zero produces

```text
c'gamma(G_1)-gamma(F_1)e'=1.                            (12)
```

Target-vector hashes bind every nonzero coefficient to its raw target
monomial as semicolon-separated
`i,j,k:numerator/denominator` strings.  Polynomial hashes bind all descending
rational coefficients, including internal zeros.  The selected data are

```text
F_1: target support 89,
     target hash af74cd37692b68ba30226421cc3b89e6efc80c6a43efe999fbc862ef4f4b155c,
     gamma degree 160,
     gamma hash d32dcfd698d46929a03a10fa9adb9118b1abf3a6c013a21ba5bd3088c5189491,
     delta degree 153,
     delta hash 226a82c7753d044927287af3f023b6e29292397047153bffcb7ff7dc3baf7566;

G_1: target support 86,
     target hash 93c4166a2d6472efd799d57e65eda841745149b33dc059528c9925f01bdac657,
     gamma degree 157,
     gamma hash 1c8115e53b579e096ed4297f7cb96efc29612c9126897d230056db4f1e643af7,
     delta degree 150,
     delta hash 3514c6aaff656d5c5e5321247df26288381b1451889fd65c2c52143a89d563ef. (13)
```

Thus `(12)` is an actual target-ring identity, not merely a restriction-ring
ansatz with unspecified representatives: the target coefficient vectors and
their representative-dependent `delta` images are frozen together.

## 2. Exact coupled lift through `J_2`

At restriction cutoff `301` there are `744` raw target monomials.  Order the
four coefficient blocks as

```text
[F_2,G_2,F_3,G_3].                                      (14)
```

Put every ascending `J_1` coefficient row before every ascending `J_2` row.
The row envelopes are `320+461=781`, and the operator has `4*744=2976`
columns.  Modulo the audit prime `1000003`, the operator and augmented ranks
are both `776`.  Augmented RREF selects the operator pivots; transpose RREF of
those columns selects independent rows.  With comma-separated decimal-index
hashes, the frozen skeleton is

```text
pivot columns:
a9697a7a900ca6ddb98a14981d01847456e04cfb48932e8c300d07065c7084d6,

pivot rows:
449e1c433c85fe550e6dda35664a1584b70fe5c4da4846003173ade228fbbfef. (15)
```

The modular calculation in `(15)` is only a pivot selector.  Rebuilding the
selected `776 x 776` square over `Q`, solving it exactly, and substituting
into the **full** rational polynomial system gives `J_1=J_2=0`.  No rank is
claimed for the full rational operator.  An independent prime `1000033`
replay also gives operator and augmented ranks `776/776`.

The chosen higher target vectors and restrictions are

```text
F_2: support 230,
     target hash 6a7e392c760108412806fa8a5edfa97786e6224ba20bd56b8761b5b3d10f2c7c,
     gamma degree 301,
     gamma hash bc865a499cf172e313db29f53af28c8ff817f3eeb0876cd78aacc75eb258aec9;

G_2: support 228,
     target hash f1940ce8dea1d9780caef84075a583e97208c6db8f21f80f39c1be0f8c8cf598,
     gamma degree 298,
     gamma hash 101fc1dbfe686a9eba0c7ce84b6e9d1a07230344526f86ae06019f44c9dc77b4;

F_3: support 230,
     target hash cafa02f3f4692af50c323406c7fdeae9a1afd759480969913b296b23d2c665cd,
     gamma degree 301,
     gamma hash 818ea15d5b7dd9dcf333ba58aaf64ac02fa47607c97aa82115782673e5164f0c;

G_3: support 86,
     target hash 5be70524a0c9dfeeb5ef3f8c26e839df241c490ffa04e8dcbe8d6b56d30ebf3a,
     gamma degree 300,
     gamma hash b895ac5bc31a323d6de69e347425e5393ef19e40bc79f5c467263f6e40f9bdb4. (16)
```

Equations `(12)` and the full residual replay after `(16)` prove `(8)`.

## 3. Composition with the proved fourth-stable identity

THM-3642 proves, for **every** target two-form on this fixed compiler and
`Q_*` fold,

```text
lambda(P)=(5P(-1)-18P(0)+13P(1))/18,

lambda(J_4)=
 (2300/81)J_0(-1)-(1/9)J_0'(-1)-(5/81)J_0''(-1)
 +(3140/81)J_0(1)-(7/9)J_0'(1)-(13/81)J_0''(1)
 +(4/9)J_2(-1)+(5/27)J_2'(-1)
 -(4/9)J_2(1)-(13/27)J_2'(1).                         (17)
```

The identity requires neither `J_1=0` nor `J_3=0`.  Substituting `(8)` gives

```text
lambda(J_4)=2300/81+3140/81=5440/81 !=0.                (18)
```

Therefore the fixed `Q_*` second-stable equality wall is genuinely inhabited
by an actual rational target-ring pair, but no actual target pair on this
fixed fold with `J_0=1,J_2=0` can have identically zero `J_4`.  In particular,
the certificate in Sections 1--2 cannot extend to a constant source
Jacobian.  This conclusion is a direct composition with PROVED THM-3642;
the present theorem does not re-prove the arbitrary-two-form identity.

## 4. Exact scope and stopping boundary

This theorem fixes compiler `(1)`, quadratic fold `(2)`, and polynomial
`Q_*` in `(3)`.  It does not claim

- minimality of the cutoffs `160` or `301`;
- a lift for every member of the degree-seven zero-second-debt line;
- `J_3=0` or a classification of the first nonzero stable coefficient;
- an all-order target pair or a polynomial Keller map on `A^2`; or
- any proof or refutation of `JC(2)`.

The two-dimensional Jacobian conjecture remains **OPEN**.

## 5. Reproduction

The companion uses active `require` gates rather than Python assertions.  It
prints and gates the target, restriction, representative-`delta`, and pivot
hashes; checks the second modular prime; solves the selected square over
`Q`; verifies the full `J_0,J_1,J_2` identities; and requires zero AST
assertion nodes.

```bash
python3 04-computation/jc2_russell_cylinder_qstar_actual_lift_fourth_stable_closure_thm3649.py
python3 -O 04-computation/jc2_russell_cylinder_qstar_actual_lift_fourth_stable_closure_thm3649.py
```

The promoted package has byte-identical normal, optimized, and stored
transcripts; its hashes, document checks, and independent hostile
reconstruction all pass.
