---
id: THM-3624
title: "Russell-cylinder non-even-fold weighted-cokernel boundary"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE / PENDING INDEPENDENT AUDIT.
  For non-even quadratic folds, the common-collision first-jet escape is an
  explicit rational surface with one zero-common subcurve removed.  Every
  nonzero escape carries a weighted vertical cokernel, but on the principal
  branch its next invoice couples two independent side jets.  An explicit
  non-even polynomial survives the enlarged arbitrary-target-two-form system
  through total source degree four.  This refutes extension of the THM-3619
  scalar staircase mechanism only; it proves neither a no-Keller theorem nor
  a Keller example for non-even folds.
source: root / noneven_fold_probe parity-hostile continuation, 2026-08-21
audit: PENDING -- provisional theorem and exact companion require hostile audit
depends_on:
  - THM-3612-russell-cylinder-even-fold-nongraph-collision-jet-rigidity
related:
  - THM-3619-russell-cylinder-even-fold-higher-jet-staircase
  - THM-3605-russell-cylinder-graph-slice-puncture-no-filling
script: 04-computation/jc2_russell_cylinder_noneven_fold_weighted_cokernel_thm3624.py
output: 05-knowledge/results/jc2_russell_cylinder_noneven_fold_weighted_cokernel_thm3624.out
script_sha256: 9064f9c7d497e86891c1f6a1ea68457262144b9966636bcbc8da5d3698fca598
output_sha256: 8d0b12b04fa089a00a2fb1bc02534e9be16d440a33d91922b83ec8493bb28ada
hash_basis: raw LF bytes
---

# THM-3624 -- Russell-cylinder non-even-fold weighted-cokernel boundary

**RESERVED / PROVISIONAL PROOF CANDIDATE / PENDING INDEPENDENT AUDIT.**
This theorem identifies exactly where evenness enters the THM-3619 collision
staircase.  The first-jet classification is complete.  The higher calculation
is deliberately finite and negative in scope: it exhibits a non-even
polynomial that the same enlarged arbitrary-two-form obstruction does not
close through total source degree four.

All rings, germs, derivatives, and differential forms are over `C`.

## 0. Setup and precise scope

Retain the THM-3612 compiler and stabilized Russell map

```text
D=1+x^2q,
c=xD(D+2),
b=(D-1)(D+2)^2,
e=q(D+3),                         c^2e=b(b+4),              (1)

B=b+cw,
C=c,
Y=ce+(2b+4)w+cw^2,
S=((b+2)(e+3w^2)+cw(3e+w^2))/8.                            (2)
```

Let `Q in C[x]` be arbitrary, not necessarily even, subject only to

```text
Q(0)=-3/4,                     Q(-1)=Q(1)=-3,              (3)
```

and take the quadratic stable-coordinate fold

```text
E_Q(x,t)=(x,Q(x)+t^2,t).                                  (4)
```

The three source points

```text
p_-=(-1,0),                 p_0=(0,0),
p_+=(1,0)                                                     (5)
```

have the common target

```text
y_0=(B,C,Y,S)=(0,0,0,-3/4).                               (6)
```

Write the three unconstrained first derivatives as

```text
u=Q'(0),                  v_-=Q'(-1),        v_+=Q'(1).    (7)
```

The results below classify the common first-order pullback of a target
two-form, identify its vertical top-symbol cokernel, and give one finite
survivor.  They do **not** prove or refute the existence of a regular target
pair with nonzero constant Jacobian on a non-even fold.

## 1. Exact non-even tangent table

Put

```text
Z=S+3/4.                                                    (8)
```

At `y_0`, the target equation `CY=B(B+4)` has differential `-4dB`,
so `(dC,dY,dZ)` is a cotangent basis.  Direct differentiation of `(1)--(4)`
gives

```text
          dC                         dY                         dZ
p_-       (12+2v_-)dx               (-36-6v_-)dx+4dt          -(v_-+9)dx/2
p_0       3dx                        -9dx+4dt                   u dx
p_+       (12-2v_+)dx               (-36+6v_+)dx+4dt           (9-v_+)dx/2.
                                                                    (9)
```

For a target cotangent two-form

```text
Omega_0=A dC wedge dY+K dC wedge dZ+R dY wedge dZ,        (10)
```

the three coefficients on `dx wedge dt` are

```text
j_-=(48+8v_-)A+2(v_-+9)R,
j_0=12A-4uR,
j_+=(48-8v_+)A+2(v_+-9)R.                                (11)
```

The coefficient `K` is invisible at first order.

## 2. Complete first-jet nonzero locus

There exist `A,R`, not both zero, for which the three values in `(11)` are
one common **nonzero** scalar if and only if

```text
4u(v_-+v_+)+4v_-v_+-27v_-+27v_+-162=0,                  (12)
```

and the following two equations do not hold simultaneously:

```text
(4u-3)v_+-24u+27=0,
(4u+3)v_-+24u+27=0.                                      (13)
```

Thus `(12)` is the tangent-match surface and `(13)` is exactly its
zero-common subcurve.

Indeed, equality of the side values with the middle value is the homogeneous
system

```text
[18-4v_+    v_+-9+2u] [A] = [0],
[18+4v_-    v_-+9+2u] [R]   [0].                         (14)
```

Its determinant is `-2` times the left side of `(12)`.  On `(12)`, the
unique projective null direction gives zero middle value precisely when the
two minors obtained by adjoining `[12,-4u]` vanish.  Those minors are `4`
and `-4` times the two expressions in `(13)`.

Every point of the nonzero locus has the following projective
parametrization.  Put `rho=R/A`.

For finite `rho` with `rho!=+-4`,

```text
v_+ = (18+rho(2u-9))/(4-rho),
v_- =-(18+rho(2u+9))/(4+rho),
u rho !=3.                                                 (15)
```

The exceptional projective values are

```text
rho=0:       v_-=-9/2,       v_+=9/2,       u arbitrary;
rho=4:       u=9/4,          v_-=-9,         v_+ arbitrary;
rho=-4:      u=-9/4,         v_+=9,          v_- arbitrary;
rho=infinity:
             v_+=9-2u,       v_-=-9-2u,      u!=0.        (16)
```

On the finite generic branch the common value is `4A(3-u rho)`; on the two
`rho=+-4` branches it is `-24A`; and on the infinite branch it is `-4uR`.
Hermite interpolation realizes every prescribed value-and-first-derivative
packet at `-1,0,1` by a polynomial.  Hence `(12)--(16)` are genuine
polynomial-fold tangent strata, not merely independent formal germs.

## 3. The weighted vertical top symbol

Write the linear target germs from `(9)` uniformly as

```text
C_i=a_i xi,              Y_i=-3a_i xi+4t,
Z_i=k_i xi,                                                (17)
```

where, in branch order `(-,0,+)`,

```text
a=(12+2v_-, 3, 12-2v_+),
k=(-(v_-+9)/2, u, (9-v_+)/2).                             (18)
```

For homogeneous target coefficient functions of degree `N`, only the pure
`Y^N` coefficient can contribute to the three pure vertical `t^N` rows.
Those three values lie in

```text
4^(N+1) span{a,k}.                                        (19)
```

On the nonzero locus `(12)--(13)`, the vectors `a,k` are independent.  If
they were dependent, the nonzero common vector in `(11)` would force their
span to be the constant vector.  Then `a=(3,3,3)` would give
`(v_-,v_+)=(-9/2,9/2)`, while `k=(-9/4,u,9/4)` could not be constant.

Consequently the unique pure-vertical top-symbol cokernel weight is

```text
                         W=a cross k.                     (20)
```

On the principal `rho=0` branch this becomes, up to a nonzero scalar,

```text
W_u=(9-4u,-18,9+4u).                                      (21)
```

Only at the even central tangent `u=0` does `(21)` reduce to the THM-3619
second difference `(1,-2,1)`.

## 4. The degree-two invoice on the principal branch

Remain on

```text
v_-=-9/2,                         v_+=9/2.                 (22)
```

Use the regular local target parameters `(c,epsilon,w)=(c,e+3,w)` before the
Russell-cylinder isomorphism.  Allow an arbitrary formal target two-form

```text
P dc wedge d epsilon+K dc wedge dw+R d epsilon wedge dw. (23)
```

This strictly enlarges the closed decomposable forms `dF wedge dG`.
Let

```text
r_-=Q''(-1),              r_0=Q''(0),       r_+=Q''(1).   (24)
```

Exact source-degree-two elimination gives the following necessary invoice
for every surviving constant target column:

```text
(9-4u)r_-+(9+4u)r_+=-243.                                (25)
```

The central jet `r_0` cancels completely.  One sparse left quotient of the
`18 x 30` arbitrary-two-form matrix evaluates on the normalized constant
column as

```text
-16/(3(4u+9))
  ((4u-9)r_--(4u+9)r_+-243),                              (26)
```

which proves `(25)` for `u!=-9/4`.  At `u=9/4`, an additional exact left-null
row evaluates to `36`, so no degree-two survivor exists at all.  Source
reflection `x -> -x`, together with the target involution
`(b,c,e,w) -> (b,-c,e,w)`, exchanges the side branches and sends `u -> -u`.
It therefore closes `u=-9/4` as well.  Thus every possible degree-two
survivor has `u!=+-9/4` and obeys `(25)`.

When evenness is imposed, `u=0` and `r_-=r_+`; equation `(25)` collapses to

```text
r_-=r_+=-27/2,                                             (27)
```

the THM-3612/THM-3619 scalar condition.  Without parity, `(25)` is only one
equation on two independent side jets.

## 5. An explicit polynomial finite-jet survivor

Define the non-even polynomial

```text
Q_h(x)=1/5408 (
  44069x^11+112059x^10-154749x^9-406377x^8
 +188107x^7+513081x^6-82835x^5-230931x^4
 +5408x-4056).                                             (28)
```

Its first four jets at the collision abscissae are

```text
x=-1:  (Q_h,Q_h',Q_h'',Q_h''')=(-3,-9/2,0,0),
x= 0:  (Q_h,Q_h',Q_h'',Q_h''')=(-3/4,1,0,0),
x= 1:  (Q_h,Q_h',Q_h'',Q_h''')
        =(-3,9/2,-243/13,10449/169).                       (29)
```

The second jets satisfy `(25)` with `u=1`.  Source degree three adds no
constant-column incompatibility for this packet.  At source degree four the
selected third jets satisfy the exact hostile relation

```text
65Q_h'''(-1)-169Q_h'''(1)+10449=0.                        (30)
```

For a direct full-system check, let

```text
T_N=(C[c,epsilon,w]_(degree<=N))^3,
S_N=direct_sum_(i=-,0,+) C[xi,t]_(degree<=N),              (31)
```

and let `P_N:T_N -> S_N` be pullback of `(23)` along the three branch germs.
Let `tau_N` be the normalized constant column: value `12` in the three
constant rows and zero elsewhere.  Exact arithmetic gives

| `N` | matrix shape | `rank P_N` | `rank[P_N|tau_N]` |
|---:|---:|---:|---:|
| 0 | `3 x 3` | 2 | 2 |
| 1 | `9 x 12` | 7 | 7 |
| 2 | `18 x 30` | 15 | 15 |
| 3 | `30 x 60` | 26 | 26 |
| 4 | `45 x 105` | 40 | 40 |

Thus one arbitrary formal target two-form jet is compatible with the
constant column through total source degree four.  This is a genuine finite
survivor of the enlarged collision-two-form obstruction.

## 6. What survives and what does not

The theorem establishes three exact facts:

1. the complete non-even first-jet nonzero locus `(12)--(16)`;
2. the unique weighted pure-vertical cokernel `(20)`, with principal weight
   `(21)`; and
3. the explicit polynomial finite survivor `(28)--(31)`.

It follows that the **scalar** THM-3619 staircase does not extend verbatim to
non-even folds.  Its unweighted side recurrence depends on parity: after
parity is removed, the first new invoice couples independent side jets, and
the polynomial `Q_h` survives the same enlarged system through degree four.

Nothing here asserts an all-order formal survivor.  In particular:

- no regular pair `F,G` is constructed;
- the compatible arbitrary two-form jet is not shown closed or decomposable;
- no Keller map or Jacobian-conjecture counterexample is produced; and
- no all-regular-target-pair no-Keller theorem for non-even folds is proved or
  refuted.

A different higher-order mechanism may still close every non-even fold.  That
question remains **OPEN**.

## 7. Exact reproduction

The deterministic companion verifies every collision value and tangent sign,
the determinant and zero-common minors, every projective stratum, the weighted
top-symbol cancellation, the sparse degree-two quotient, the exceptional
`u=9/4` debt and reflection gate, all displayed jets of `Q_h`, and every rank
in the table.

Reproduce with

```bash
python3 04-computation/jc2_russell_cylinder_noneven_fold_weighted_cokernel_thm3624.py
python3 -O 04-computation/jc2_russell_cylinder_noneven_fold_weighted_cokernel_thm3624.py
```

Both streams must be byte-identical to
`05-knowledge/results/jc2_russell_cylinder_noneven_fold_weighted_cokernel_thm3624.out`.
