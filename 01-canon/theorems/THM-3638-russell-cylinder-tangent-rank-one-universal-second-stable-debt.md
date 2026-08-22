---
id: THM-3638
title: "Russell-cylinder tangent-rank-one universal second-stable debt"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For the minimal rho=0,u=1 Hermite polynomial, every normalized
  zero-stable determinant witness whose two horizontal restriction
  derivatives span only one retained tangent direction has the same
  nonzero second-stable cokernel -2072/81.  The calculation is symbolic in
  every zero-stable Hessian and first-stable gradient gauge, is invariant
  under target representatives modulo ker(gamma), and has no target-degree
  bound.  A fixed actual-ring rank-one control exists at t=0 but its minimal
  degree-102 witness already fails J_1 locally.  General retained-tangent
  rank two is not treated here, and JC(2) remains OPEN.
source: root / general_rank_pair_frontier rank-stratum continuation, 2026-08-21
audit: >
  PASS -- an independent hostile audit reconstructed the complete rank-one
  output normalization, symbolic Hessian invoice, representative invariance,
  universal quotient, degree-101/102 actual-ring boundary, local J_1 hostile,
  and byte-identical normal, optimized, and stored 44-gate transcripts.
depends_on:
  - THM-3634-russell-cylinder-quadratic-fold-first-stable-jet-rank-debt
  - THM-3635-russell-cylinder-retained-curve-jet-plane-actual-rank-witness
  - THM-3637-russell-cylinder-actual-rank-witness-second-stable-debt
related:
  - THM-3630-russell-cylinder-noneven-formal-survivor-no-finite-jet-bound
  - THM-3632-russell-cylinder-formal-pair-algebraization-triple-fibre-obstruction
script: 04-computation/jc2_russell_cylinder_tangent_rank_one_universal_second_stable_debt_thm3638.py
output: 05-knowledge/results/jc2_russell_cylinder_tangent_rank_one_universal_second_stable_debt_thm3638.out
script_sha256: 7514f2b6b1c306d3b08aa50753be5b042d20ef2057a81748d565b0b95991b380
output_sha256: d31d23c54e3ec10c67215835993feea12ea2085547a60ba115cb2cd90c6292ff
hash_basis: raw LF bytes
---

# THM-3638 -- Russell-cylinder tangent-rank-one universal second-stable debt

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
This theorem closes the whole retained-tangent-rank-one stratum left after
THM-3634.  Its rank is local and must not be confused with the global ranks
from that theorem:

```text
closed here: rank span{dU,dV}=1 in the retained tangent plane T;

not treated here: rank span{dU,dV}=2 for arbitrary zero-stable pairs U,V. (1)
```

The horizontal restrictions `U,V` and the first-stable restrictions `A,B`
may still each have global rank two.  Indeed Section 6 gives an actual-ring
rank-two/rank-two `t=0` control whose retained horizontal tangent rank is one.

All rings and closed points are over `C`.  The companion works exactly over
`Q` and then base-changes.

## 0. Setup and theorem

Use the exponent-two Danielewski surface and compiler

```text
R=C[b,c,e]/(c^2e-b(b+4)),

D=1+x^2q,
b=(D-1)(D+2)^2,       c=xD(D+2),       e=q(D+3).       (2)
```

Fix the minimal THM-3635 Hermite polynomial

```text
Q=x^5+(9/2)x^4-2x^3-(27/4)x^2+x-3/4,                  (3)
```

and put

```text
gamma(K)=K(b(x,Q(x)),c(x,Q(x)),e(x,Q(x))),
S=gamma(R),

delta(K)=partial_q K(b(x,q),c(x,q),e(x,q)) |_(q=Q(x)). (4)
```

Write `z=e+3`.  Let

```text
F^#=sum_(i>=0) w^i F_i,          G^#=sum_(i>=0) w^i G_i,
F_i,G_i in R,                                                    (5)
```

and pull back through the quadratic fold

```text
(x,t) |-> (x,Q(x)+t^2,w=t).                              (6)
```

Define the zero- and first-stable restrictions

```text
U=gamma(F_0),       V=gamma(G_0),
A=gamma(F_1),       B=gamma(G_1).                        (7)
```

Suppose their zero-stable determinant is a nonzero constant and their
retained derivative columns have rank one:

```text
U'B-AV'=kappa in C*,

dim_C span{dU,dV}=1,                                    (8)

dP=(P'(-1),P'(0),P'(1)).                                (9)
```

Then no choice of the higher `F_i,G_i` gives a constant source Jacobian.
After a constant output normalization taking `kappa` to `1`, the exact
second-stable coefficient always satisfies

```text
lambda(J_2)=-2072/81 != 0,                              (10)

lambda(P)=(5P(-1)-18P(0)+13P(1))/18.                   (11)
```

The assertion has no degree bound.  It quantifies over every target
representative of every restriction class and every higher target
coefficient.

## 1. Stable expansion

Put

```text
C=delta(F_0)+gamma(F_2),       E=delta(G_0)+gamma(G_2),
D=delta(F_1)+gamma(F_3),       H=delta(G_1)+gamma(G_3). (12)
```

The critical substitution `q=Q+t^2` gives

```text
f=U+tA+t^2C+t^3D+O(t^4),
g=V+tB+t^2E+t^3H+O(t^4).                               (13)
```

Thus the first three Jacobian coefficients are

```text
J_0=U'B-AV',

J_1=2U'E+A'B-AB'-2CV',

J_2=3U'H+2A'E+C'B-AE'-2CB'-3DV'.                      (14)
```

The terms `delta(F_1),delta(G_1)` in `D,H` are part of the decomposable
target-pair calculation; they are not optional arbitrary-two-form data.

## 2. Retained geometry and representative invariance

At `x=-1,0,1`, the restricted curve has the common target point
`(b,c,z)=(0,0,0)`.  The local surface coordinates are `(c,z)`, and the
retained tangent and vertical columns are

```text
tau_c=( 3,3,3),                 tau_z=(-9,4,9),

n_c=delta(c)=( 2,0,-2),         n_z=delta(z)=(-2,4,-2),

n_c'=(-9,0,-9),                 n_z'=(3,0,-3),

c''=(157/2,0,-221/2).                                  (15)
```

The three tangent directions in the `(c,z)` plane are

```text
(3,-9),                  (3,4),                  (3,9). (16)
```

Their pairwise determinants are `(39,15,54)`, so this is the ordinary
triple of THM-3635.  Its retained tangent plane is

```text
T=span{tau_c,tau_z}.                                    (17)
```

Every `P in S` has one common retained value and `dP in T`.  The row
`lambda` in `(11)` annihilates both generators of `T`.

The target-representative quantifier is load-bearing.  If `K in ker(gamma)`,
its germ at the common point vanishes on all three branches in `(16)`.  Its
linear initial form is therefore zero.  Its quadratic initial form is also
zero, since a nonzero binary quadratic has at most two projective zero
directions.  Hence

```text
K in m^3.                                               (18)
```

The operator `delta` differentiates once in a tangent direction of the
smooth surface.  Consequently

```text
delta(K)(x_i)=0,             (delta(K))'(x_i)=0,
                             x_i=-1,0,1.                (19)
```

Thus zero-stable representative changes do not affect the value or first
branch derivative of the vertical terms used below.  First-stable
representative changes do not affect their retained vertical values.

## 3. Complete rank-one output normalization

Evaluating `J_0=kappa` at the retained triple shows

```text
B_0 dU-A_0 dV=kappa (1,1,1),                            (20)
```

where `A_0,B_0` are the common retained values.  Under the rank-one
hypothesis, the one-dimensional span in `(20)` is therefore exactly the
all-ones line.

Choose a constant output matrix whose first row sends `(dU,dV)` to
`tau_c=3(1,1,1)` and whose independent second row kills it.  Scale the
second row to normalize the transformed determinant to `1`.  This gives

```text
dU=tau_c,                    dV=0,
J_0=1,                       B_0=1/3.                  (21)
```

The residual determinant-one shear

```text
(F^#,G^#) |-> (F^#-3A_0 G^#,G^#)                       (22)
```

preserves `(21)` and changes the first-stable value from `A_0` to zero.
Hence every rank-one cell has the complete normalization

```text
dU=tau_c,          dV=0,          A_0=0,          B_0=1/3. (23)
```

Because `(tau_c,tau_z)` are independent, the local target gradients of
`F_0,G_0` in `(c,z)` are respectively `(1,0)` and `(0,0)`.

## 4. The symbolic Hessian invoice

Write the arbitrary local Hessian of `F_0` as

```text
H_U=[[u_cc,u_cz],[u_cz,u_zz]],                          (24)
```

and write the target gradient of `G_1` as `(b_c,b_z)`.  Then

```text
dB=b_c tau_c+b_z tau_z.                                 (25)
```

No other zero-stable second derivative or first-stable gradient is fixed.
At branch `i`, with tangent `v_i` from `(16)`, the second derivative of `U`
is

```text
U_i''=c_i''+v_i^T H_U v_i.                              (26)
```

Differentiate `J_0=1`.  Under `(23)` the complete three-row identity is

```text
U''/3+3dB=0.                                            (27)
```

This system has rank three and is equivalent to

```text
u_cc=-3b_c-3658/585,

u_cz=7/4-(3/2)b_z,

u_zz=58/65.                                             (28)
```

Thus `b_c,b_z` remain wholly arbitrary.  Equation `(28)` is a symbolic
identity on the full Hessian/gradient space, not a bounded target-monomial
search.

For comparison, before the residual shear `(22)`, let `A_0=a` and let
`v_zz` be the `zz` Hessian entry of `G_0`.  Applying `lambda` to the
differentiated determinant gives the exact compatibility

```text
195 a v_zz=65 u_zz-58.                                  (29)
```

This includes the apparently degenerate boundary `v_zz=0,u_zz=58/65`; it
does not create an escape.

## 5. Universal second-stable quotient

Apply `lambda` to `J_2` in `(14)` under `(23)`.  Every higher restriction
coefficient has common value and derivative in `T`.  Moreover
`delta(G_0)` has retained value zero because the target gradient of `G_0`
is zero.  Term by term:

- `3U'H` contributes `9 lambda(delta(G_1))`;
- `2A'E` vanishes because `E` has common retained value and `lambda(dA)=0`;
- `C'B` contributes `(1/3)lambda((delta(F_0))')`;
- `AE'` vanishes because `A_0=0`;
- the common part of `C` is killed against `dB`, leaving
  `-2lambda(delta(F_0)*dB)`; and
- `DV'` vanishes because `dV=0`.

Therefore the complete quotient is

```text
lambda(J_2)=9 lambda(delta(G_1))
            +(1/3)lambda((delta(F_0))')
            -2 lambda(delta(F_0)*dB).                  (30)
```

In the local coordinates, its three ingredients are

```text
delta(G_1)_i=(b_c,b_z) dot (n_c,i,n_z,i),

(delta(F_0))'_i=v_i^T H_U (n_c,i,n_z,i)+n_c,i',

delta(F_0)|_triple=n_c.                                (31)
```

Substituting the full affine solution `(28)` into `(30)` cancels both
`b_c` and `b_z` and gives

```text
lambda(J_2)=-2072/81.                                   (32)
```

The unnormalized calculation is an independent covariance check:

```text
lambda(J_2)=40(117a v_zz-39u_zz-17)/81,

lambda(J_2)+2072/81
  =-(8/27)(65u_zz-58-195a v_zz).                       (33)
```

Equation `(29)` again forces `(32)`.  This proves the theorem for every
rank-one cell, including the `v_zz=0` boundary.

As an active hostile control, increasing only the left-branch entry of
`n_c'` by one changes the quotient from `-2072/81` to `-4129/162`.

## 6. An actual-ring `t=0` control and its earlier `J_1` failure

To show that retained tangent rank one is not vacuous inside the global
rank-two/rank-two cell, take

```text
U=c,                         V=z^2=(e+3)^2.             (34)
```

Then `[U],[V]` are globally independent modulo constants, while

```text
dU=tau_c,                         dV=0.                 (35)
```

Use the complete THM-3635 `C[z]`-module basis of `S`.  For the symmetric
source-degree filtration, the exact map

```text
Phi_d:S_(<=d) direct_sum S_(<=d) -> C[x],
Phi_d(A,B)=U'B-AV'                                      (36)
```

has

```text
d=101: dim S_(<=d)=61, matrix 125 x 122,
       ranks (Phi_d,[Phi_d|1])=(122,123);

d=102: dim S_(<=d)=62, matrix 126 x 124,
       ranks (Phi_d,[Phi_d|1])=(124,124).              (37)
```

Thus there is no solution through degree `101`, and there is one unique
solution under the symmetric cutoff `102`.  It has

```text
deg A=93,                         deg B=102,

A(-1)=A(0)=A(1)=-29/195,
B(-1)=B(0)=B(1)=1/3,                              (38)

sha256(A)=a17fb9620b7d2ee308f8ca1f311205ac619cc4e70886e641d7ed4b0863d5c73a,

sha256(B)=d43b93e89fcf1aa61bcfce9b3f497ff2d5c723cadec9c7313e52a1238a2ddf26.
                                                                  (39)
```

The two sidecars are globally independent.  This is a genuine actual-ring
rank-two/rank-two `J_0` witness, but it does **not** supply a positive
`J_1` control.

Indeed write

```text
dA=a_c tau_c+a_z tau_z,           dB=b_c tau_c+b_z tau_z. (40)
```

The exact witness has `b_z=7/6`.  At the retained triple,
`delta(G_0)=delta(z^2)=0` and `V'=0`, so a necessary condition from
`J_1=0` is

```text
a_z=3A_0 b_z=-203/390.                                  (41)
```

The unique witness in `(38)` violates `(41)`; the reduced-rational mismatch
has hash

```text
a48a97b7b9d3d930c9472c610004657e5728cf0dfa0297429a9738482860decc. (42)
```

Hence its `J_1` equation is not even locally solvable at the retained
triple, and therefore is not globally solvable in `S`.  This earlier failure
is only about the one minimal witness `(38)`.  The universal theorem above
says that any other rank-one witness which does reach `J_1=0` still fails at
`J_2` by `(32)`.

## 7. Stopping boundary

The proved candidate payload is

```text
fixed fold Q_1 + retained tangent rank one:
    universal second-stable debt -2072/81;

fixed actual U=c,V=(e+3)^2 minimal witness:
    exact J_0 survivor, but J_1 already locally impossible;

arbitrary retained tangent rank two:
    not part of this theorem;

global quadratic-fold pair / Keller / JC(2):
    OPEN.                                                (43)
```

In particular, this theorem neither promotes a finite local survivor to a
global pair nor turns the formal multi-germ of THM-3630 into one polynomial
fold.  It proves no automorphism theorem and no Jacobian-conjecture claim.

## 8. Exact verification

The companion verifies the three stable identities, compiler and retained
jets, ordinary triple, tangent cokernel, graph-ideal order-three lemma,
complete output normalization equations, the full symbolic Hessian invoice,
the unnormalized shear-covariance factorization, the actual-ring
degree-101/102 boundary and unique witness, the local `J_1` hostile control,
coefficient hashes, and an assertion-free AST.

Reproduce with

```bash
python3 04-computation/jc2_russell_cylinder_tangent_rank_one_universal_second_stable_debt_thm3638.py
python3 -O 04-computation/jc2_russell_cylinder_tangent_rank_one_universal_second_stable_debt_thm3638.py
```

Both streams must be byte-identical to

```text
05-knowledge/results/jc2_russell_cylinder_tangent_rank_one_universal_second_stable_debt_thm3638.out
```

The frozen companion reports zero Python assertion statements.
