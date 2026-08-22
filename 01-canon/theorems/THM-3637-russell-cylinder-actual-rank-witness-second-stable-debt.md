---
id: THM-3637
title: "Russell-cylinder actual-rank witness second-stable debt"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For the minimal
  rho=0,u=1 Hermite polynomial and the fixed zero-stable restriction classes
  gamma(F_0)=c, gamma(G_0)=e, the THM-3635 actual-ring determinant witness
  extends through the first stable Jacobian identity.  Nevertheless every
  normalized determinant witness in those restriction classes that reaches
  J_1=0 has the same nonzero second-stable retained cokernel -2072/81.  This is
  independent of all unbounded first-jet gauges, all higher restriction lifts,
  and all target representatives modulo ker(gamma).  It excludes only this
  fixed cell and makes no global Keller or JC(2) claim.
source: root / audit_thm3632 second-stable continuation, 2026-08-21
audit: >
  PASS -- independent hostile reconstruction verified the stable expansion,
  all unbounded retained gauges, the order-three graph-ideal representative
  invariance, every vertical first-coefficient contribution, the exact
  first-order lift, and byte-identical normal, optimized, and stored 59-gate
  transcripts.
depends_on:
  - THM-3634-russell-cylinder-quadratic-fold-first-stable-jet-rank-debt
  - THM-3635-russell-cylinder-retained-curve-jet-plane-actual-rank-witness
related:
  - THM-3624-russell-cylinder-noneven-fold-weighted-cokernel-boundary
  - THM-3630-russell-cylinder-noneven-formal-survivor-no-finite-jet-bound
  - THM-3632-russell-cylinder-formal-pair-algebraization-triple-fibre-obstruction
script: 04-computation/jc2_russell_cylinder_actual_rank_witness_second_stable_debt_thm3637.py
output: 05-knowledge/results/jc2_russell_cylinder_actual_rank_witness_second_stable_debt_thm3637.out
script_sha256: cc59af1d4188d96072dd3edd1f5125f78fafd576b26b61d662f105c14c6fdca6
output_sha256: 32b91e0927708fc1c471fdbd2c680f20306a59920bd7e0db73427a59c7d5ed5d
hash_basis: raw LF bytes
---

# THM-3637 -- Russell-cylinder actual-rank witness second-stable debt

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This theorem
continues the actual rank-two/rank-two point of THM-3635 through one more
stable order.  Its positive and negative statements have different scopes:

```text
positive:  one exact actual-ring witness extends through J_1=0;

negative:  every J_0=1 witness with zero-stable restrictions (c,e)
           that reaches J_1=0 has the same nonzero quotient at J_2.       (1)
```

The positive construction uses the exact target representatives `F_0=c` and
`G_0=e`.  The negative statement fixes only their **restriction classes**;
it does not fix their representatives in the target ring.

All rings and closed points are over `C`.  The exact companion works over
`Q` and then base-changes.

## 0. Setup and statement

Use the exponent-two Danielewski surface and compiler

```text
R=C[b,c,e]/(c^2e-b(b+4)),

D=1+x^2q,
b=(D-1)(D+2)^2,        c=xD(D+2),        e=q(D+3).       (2)
```

Fix the minimal THM-3635 Hermite polynomial

```text
Q=x^5+(9/2)x^4-2x^3-(27/4)x^2+x-3/4,                  (3)
```

and define

```text
gamma:R -> C[x],
gamma(K)=K(b(x,Q(x)),c(x,Q(x)),e(x,Q(x))),

S=gamma(R),

delta(K)=partial_q K(b(x,q),c(x,q),e(x,q)) |_(q=Q(x)). (4)
```

As in THM-3635, the same letters `b,c,e` below also denote their restricted
polynomials in `S` when no confusion is possible.

Let

```text
F^#=sum_(i>=0) w^i F_i,       G^#=sum_(i>=0) w^i G_i,
F_i,G_i in R,                                             (5)
```

and pull back by the quadratic fold

```text
(x,t) |-> (x,Q(x)+t^2,w=t).                              (6)
```

Suppose only that

```text
gamma(F_0)=c,                  gamma(G_0)=e.              (7)
```

Then there is no such pair whose source Jacobian is identically `1`.
More precisely, write `A=gamma(F_1), B=gamma(G_1)`.  If the zero-stable
identity

```text
c'B-Ae'=1                                                 (8)
```

holds and the first-stable identity also vanishes, then the second-stable
coefficient has the exact retained quotient

```text
lambda(J_2)=-2072/81 != 0,                                (9)

lambda(P)=(5P(-1)-18P(0)+13P(1))/18.                     (10)
```

This assertion quantifies over

- every `A,B in S` satisfying `(8)`, with no degree bound;
- every target representative of `F_0,G_0,F_1,G_1` in its prescribed
  restriction class; and
- every choice of `F_2,G_2,F_3,G_3 in R`.

On the positive side, the canonical degree-`94/91` witness of THM-3635 has
an exact choice of `F_2,G_2` for which the first-stable coefficient vanishes,
using `F_0=c,G_0=e`.  Thus the obstruction in `(9)` genuinely occurs one
stable order after the determinant witness, not already at `J_1`.

## 1. The first three stable identities

Put

```text
U=gamma(F_0),                 V=gamma(G_0),
A=gamma(F_1),                 B=gamma(G_1),

C=delta(F_0)+gamma(F_2),      E=delta(G_0)+gamma(G_2),
D=delta(F_1)+gamma(F_3),      H=delta(G_1)+gamma(G_3).  (11)
```

The critical quadratic substitution gives

```text
f=U+tA+t^2C+t^3D+O(t^4),
g=V+tB+t^2E+t^3H+O(t^4).                                (12)
```

Consequently

```text
J_0=U'B-AV',

J_1=2U'E+A'B-AB'-2CV',

J_2=3U'H+2A'E+C'B-AE'-2CB'-3DV'.                       (13)
```

All `delta(F_1),delta(G_1)` terms in `D,H` are load-bearing.  They will
contribute nontrivially to `(9)`.

## 2. Retained tangent and cokernel data

At `x=-1,0,1`, the restricted curve has the common target value
`p=(b,c,e)=(0,0,-3)`.  Since the `b`-derivative of
`c^2e-b(b+4)` at `p` is `-4`, the surface is smooth there and `(c,e+3)` are
regular local coordinates.  The curve tangent columns and the vertical
compiler columns are

```text
tau_c=c'= ( 3,3,3),           tau_e=e'=(-9,4,9),

n_c=delta(c)=( 2,0,-2),       n_e=delta(e)=(-2,4,-2),

(delta(c))'=(-9,0,-9).                                  (14)
```

The three pairwise tangent determinants are `(39,15,54)`, so this is the
ordinary retained triple of THM-3635.  Every `P in S` has one common retained
value and

```text
dP=(P'(-1),P'(0),P'(1)) in
T=span{tau_c,tau_e}.                                    (15)
```

The row in `(10)` annihilates both generators of `T`.  Two consequences will
be used repeatedly:

```text
lambda(dP)=0                         for P in S,

lambda(c'Y-e'X)=0                    for X,Y in S.       (16)
```

The second identity follows because `X,Y` have common retained values.

### Representative-invariance lemma

The restriction-class quantifier in the theorem is stronger than fixing
`F_0=c,G_0=e` as target functions.  Let `K in ker(gamma)`.  Its germ at `p`
vanishes on all three retained branches.  Its linear initial form therefore
vanishes on three distinct tangent lines and is zero.  Its quadratic initial
form also vanishes on those lines; a nonzero binary quadratic has at most two
projective zero directions.  Hence

```text
K in m_p^3.                                              (17)
```

The `q`-variation in `(4)` is a tangent direction to the smooth surface at
the retained point.  Differentiating `(17)` once in that direction leaves
branch order at least two:

```text
delta(K)(x_i)=0,             (delta(K))'(x_i)=0,
                             x_i=-1,0,1.                 (18)
```

Thus replacing `F_0` or `G_0` by another representative of the same
restriction class changes neither the values nor first derivatives from
`delta(F_0),delta(G_0)` that occur in `J_2`.  Replacing `F_1,G_1` by other
representatives changes neither retained value of `delta(F_1),delta(G_1)`.
There is no hidden representative gauge in `(9)`.

## 3. Every unbounded zero-stable gauge

Let `A_0,B_0` be the common retained values of `A,B`.  Evaluating `(8)` at
the three retained points gives

```text
3B_0-A_0(-9,4,9)=(1,1,1),

A_0=0,                         B_0=1/3.                 (19)
```

By `(15)`, write uniquely

```text
dA=p tau_c+q tau_e,            dB=r tau_c+s tau_e.      (20)
```

The second derivative of `c` at the retained points is

```text
c''=(157/2,0,-221/2).                                    (21)
```

Differentiate `(8)` and use `(19)`:

```text
c''/3+3dB-tau_e*dA=0,                                    (22)
```

where `*` is coordinatewise multiplication.  The `3 x 4` affine system
`(22)` has rank three.  Its complete solution is

```text
p=sigma-7/6,              q=-58/195,
r=-3658/1755,             s=sigma,       sigma in C.    (23)
```

This is the full unbounded gauge, not a bounded target-monomial search.
Arbitrary high-degree homogeneous Bezout corrections can change `sigma`,
but no other retained first-jet coordinate survives.

The retained values of the vertical target derivatives are fixed by the
same gradients:

```text
delta(F_1)|_triple=p n_c+q n_e,
delta(G_1)|_triple=r n_c+s n_e.                          (24)
```

Equation `(18)` proves that `(24)` is independent of the chosen target
representatives of `A,B`.

## 4. The first-stable identity fixes two values

Put

```text
X=gamma(F_2),                  Y=gamma(G_2),             (25)
```

and let `x_0,y_0` be their common retained values.  By `(14),(18)`, the
retained values of `C,E` are respectively

```text
n_c+x_0(1,1,1),                n_e+y_0(1,1,1).           (26)
```

Evaluate `J_1=0` from `(13)` and use `(19)`:

```text
2 tau_c*(n_e+y_0 1)+(1/3)dA
                   -2(n_c+x_0 1)*tau_e=0.              (27)
```

This rank-two consistent system has the unique solution

```text
x_0=-29/585,                  y_0=-sigma/6-137/36.       (28)
```

No other coefficient of `X,Y` will enter the retained quotient of `J_2`:
their derivatives lie in `T` and are killed by `lambda`.

For completeness, `J_1=0` is globally solvable for the canonical THM-3635
witness with the exact representatives `F_0=c,G_0=e`.  In the complete
`C[z]`-module basis of `S`, the filtered piece `S_(<=174)` has dimension
`134`.  The exact coefficient matrix of

```text
(X,Y) |-> c'Y-e'X                                      (29)
```

on `S_(<=174)^2` has shape `189 x 268`, rank `188`, and augmented rank `188`
for the `J_1` right-hand side.  Setting every free RREF variable to zero gives

```text
deg X=173,                         deg Y=91,

sha256(X)=13403da7b39e92a0ae35b690a07184d673de63361c8bae243c7041e3e45081ef,

sha256(Y)=6432e36c32a516a2204f315142eee0dba83fd0cbc0400f213038ef5a864a5e62.
                                                                  (30)
```

Every module basis element is the restriction of an element of `R`, so
these `X,Y` have target lifts `F_2,G_2`.  Exact substitution gives `J_1=0`.

## 5. The second-stable retained debt

Write `Z=gamma(F_3), W=gamma(G_3)`.  Their retained values are common, so
`lambda` kills the term

```text
3(c'W-Ze').                                               (31)
```

At the retained triple, `A=0`, `B=1/3`.  Equations `(16),(18)` also kill
the `S`-derivative and representative-gauge parts of `C'B`.  What remains
from `J_2` has five contributions.  Using `(14),(23),(24),(28)`, their exact
`lambda`-values are

```text
term                         lambda-value

3c' delta(G_1)              -2(4r+27s)
2A'E                        -12(3p+4q)
C'B                         -3
-2CB'                       (4/3)(4r+27s)
-3 delta(F_1)e'             18(3p+4q).                  (32)
```

In particular, the two vertical `F_1,G_1` rows are not optional.  Summing
`(32)` gives

```text
lambda(J_2)=(54p+72q-8r-54s-9)/3.                       (33)
```

Substitution of the complete affine line `(23)` cancels `sigma` and yields

```text
lambda(J_2)=-2072/81.                                    (34)
```

This contradicts `J_2=0` and proves the theorem.  As an active hostile
control, deleting the two vertical `F_1,G_1` contributions changes the
computed quotient to `3415/81`; their combined contribution is `-1829/27`.

## 6. Quantifier and stopping boundary

The conclusion is exactly a no-extension theorem for the one fixed cell

```text
Q=Q_1,          gamma(F_0)=c,          gamma(G_0)=e,
                 normalized source Jacobian=1.          (35)
```

It fixes restriction classes, not target representatives.  It is also
stronger than excluding only the one displayed THM-3635 pair: every
unbounded `A,B in S` satisfying `(8)` is covered by `(23)`.

The theorem does **not** prove

- a second-stable obstruction for other zero-stable pairs `U,V`;
- an obstruction for every scalar normalization while keeping `(c,e)` fixed;
- termination or algebraization failure for a general formal-in-`w` survivor;
- nonexistence of all quadratic-fold global pairs; or
- the two-dimensional Jacobian conjecture.

A global pair outside `(35)` would require a new argument.  Conversely, the
positive `J_1` lift in `(30)` is only a finite stable jet, not a polynomial
Keller map.

## 7. Exact verification

The companion verifies the three stable identities, compiler and vertical
jets, retained cokernel, graph-ideal order-three lemma, reconstruction of the
THM-3635 witness, exact first-order extension, full affine retained gauge,
all five second-order quotient terms, an active vertical-term mutation, and
an assertion-free AST.

Reproduce with

```bash
python3 04-computation/jc2_russell_cylinder_actual_rank_witness_second_stable_debt_thm3637.py
python3 -O 04-computation/jc2_russell_cylinder_actual_rank_witness_second_stable_debt_thm3637.py
```

Both streams must be byte-identical to

```text
05-knowledge/results/jc2_russell_cylinder_actual_rank_witness_second_stable_debt_thm3637.out
```

The frozen companion reports zero Python assertion statements.
