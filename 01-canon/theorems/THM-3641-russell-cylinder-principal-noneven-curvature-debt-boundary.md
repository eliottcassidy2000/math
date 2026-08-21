---
id: THM-3641
title: "Russell-cylinder non-even curvature debt atlas and zero-debt boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On every
  ordinary retained triple in the full projective first-jet atlas, the
  decomposable quadratic-fold second-stable cokernel has an exact all-cell
  derivative identity.  Its constant term is an explicit affine function of
  the three curvatures of Q.  The principal chart recovers exactly the
  THM-3624 arbitrary-two-form invoice, and every ordinary slope cell has
  polynomial zero-debt controls.  Zero debt removes only this retained J_2
  quotient; it supplies no actual J_0 witness, J_1 lift, global J_2 identity,
  Keller pair, or JC(2) consequence.
source: root/general-q-curvature-boundary/2026-08-21
audit: >
  Independent hostile reconstruction rederived the generic finite-rho,
  rho=+/-4, infinite, and principal identities; checked every ordinary-locus
  and representative quantifier; and verified the polynomial controls and
  Hermite boundary.  Normal, optimized, and stored transcripts are
  byte-identical at 104 active gates; hashes, AST0, docs, and diff checks pass.
depends_on:
  - THM-3624-russell-cylinder-noneven-fold-weighted-cokernel-boundary
  - THM-3635-russell-cylinder-retained-curve-jet-plane-actual-rank-witness
  - THM-3639-russell-cylinder-all-retained-cells-universal-second-stable-debt
related:
  - THM-3630-russell-cylinder-noneven-formal-survivor-no-finite-jet-bound
  - THM-3637-russell-cylinder-actual-rank-witness-second-stable-debt
script: 04-computation/jc2_russell_cylinder_principal_noneven_curvature_debt_boundary_thm3641.py
output: 05-knowledge/results/jc2_russell_cylinder_principal_noneven_curvature_debt_boundary_thm3641.out
script_sha256: f82fc5480a63b7cc562086cd10ea92005c4398d3720a922216c4c7b88deda9db
output_sha256: 8db1298888729dee2c1f9db720af1219504f723fb518c9683f7d062ccb612262
hash_basis: raw LF bytes
---

# THM-3641 -- Russell-cylinder non-even curvature debt atlas and zero-debt boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This theorem
replaces the single fixed collision polynomial in THM-3639 by the
full projective first-jet atlas of THM-3624.  On every ordinary retained
triple it gives an exact identity

```text
lambda(J_2)=D_Q+mu_- J_0'(-1)+mu_0 J_0'(0)+mu_+ J_0'(1).       (1)
```

The identity is uniform across zero-stable tangent ranks one and two and is
division-free in that rank parameter.  The scalar `D_Q` is the retained
curvature debt.  Unlike the fixed value `-2072/81` for the minimal polynomial
of THM-3639, it can vanish.  The generic, exceptional, and infinite charts
below give its complete ordinary-triple atlas.

All rings, germs, and closed points are over `C`.  The companion verifies the
displayed identities over `Q` and then base-changes.

## 0. Compiler, fold, and exact scope

Use the exponent-two Danielewski surface and compiler

```text
R=C[b,c,e]/(c^2 e-b(b+4)),

D=1+x^2 q,
b=(D-1)(D+2)^2,        c=xD(D+2),        e=q(D+3).            (2)
```

Let `Q in C[x]` obey

```text
Q(-1)=Q(1)=-3,                         Q(0)=-3/4,              (3)
```

and write

```text
v_-=Q'(-1),          u=Q'(0),          v_+=Q'(1),
r_-=Q''(-1),         r_0=Q''(0),       r_+=Q''(1).             (4)
```

Put `z=e+3`, restrict `(2)` through `q=Q(x)`, and define

```text
gamma(K)=K(b(x,Q(x)),c(x,Q(x)),e(x,Q(x))),
delta(K)=partial_q K(b(x,q),c(x,q),e(x,q))|_(q=Q(x)).          (5)
```

For actual target coefficients `F_i,G_i in R`, set

```text
F^#=sum_(i>=0) w^i F_i,                 G^#=sum_(i>=0) w^i G_i
```

and pull back by the quadratic fold

```text
(x,t) |-> (x,Q(x)+t^2,w=t).                                  (6)
```

With

```text
U=gamma(F_0),        V=gamma(G_0),
A=gamma(F_1),        B=gamma(G_1),

C=delta(F_0)+gamma(F_2),       E=delta(G_0)+gamma(G_2),
L=delta(F_1)+gamma(F_3),       H=delta(G_1)+gamma(G_3),       (7)
```

the first stable coefficients of the source Jacobian are

```text
J_0=U'B-AV',

J_1=2U'E+A'B-AB'-2CV',

J_2=3U'H+2A'E+C'B-AE'-2CB'-3LV'.                           (8)
```

The candidate assumes the nonzero first-jet cell condition of THM-3624 and
normalizes the common retained values of `J_0` to one:

```text
J_0(-1)=J_0(0)=J_0(1)=1.                                   (9)
```

For the zero-debt consequence one strengthens `(9)` to the polynomial
identity `J_0=1`.  No `J_1=0` hypothesis is needed for `(1)`.

## 1. Universal compiler jets

At `x=-1,0,1`, the restricted curve has common target point
`(b,c,z)=(0,0,0)`.  Its tangent columns in the regular local coordinates
`(c,z)` are

```text
tau_c=(2(v_-+6), 3, 2(6-v_+)),
tau_z=(-2(v_-+9), 4u, 2(9-v_+)).                          (10)
```

The vertical compiler columns and their branch derivatives are

```text
n_c=delta(c)=(2,0,-2),             n_z=delta(z)=(-2,4,-2),

m_c=(delta(c))'=(-18-2v_-,0,-18+2v_+),
m_z=(delta(z))'=(12+2v_-,0,-12+2v_+).                    (11)
```

The retained curve curvatures are

```text
c''=(2(r_--v_-^2-18v_--54),
     0,
     2(v_+^2-18v_++54-r_+)),

z''=(2(v_-^2+12v_-+9-r_-),
     4r_0+9/8,
     2(v_+^2-12v_++9-r_+)).                              (12)
```

The mechanism behind the atlas is the branchwise symplectic identity

```text
det((tau_c,i,tau_z,i),(n_c,i,n_z,i))=12,
                                      i=-,0,+.             (13)
```

The first-jet surface is

```text
4u(v_-+v_+)+4v_-v_+-27v_-+27v_+-162=0,                  (14)
```

with the zero-common subcurve removed as in THM-3624.  On this locus the
constant vector belongs to `span{tau_c,tau_z}`.

## 2. The all-cell cancellation mechanism

Choose linear target coordinates `(y,eta)` for which the three retained
tangents are

```text
t_i=(1,h_i),                                               (15)
```

and write the vertical columns and the needed derivative/curvature as

```text
n_i=(p_i,q_i),          m_i=(delta(y))'_i,          s_i=y''_i. (16)
```

Equation `(13)` says that `q_i-h_i p_i` is one common scalar.  Let `lambda`
be the row normalized by

```text
lambda(1)=lambda(h)=0,                    lambda_0=-1.   (17)
```

On an ordinary triple the three `h_i` are distinct.  There is then a unique
row `mu` with moments

```text
mu(1)=lambda(p),
mu(h)=lambda(3q-2hp),
mu(h^2)=lambda(hq).                                      (18)
```

The common value of `q-hp` makes `(18)` simultaneously cancel every Hessian
and first-stable-gradient gauge in `J_2` against the derivatives
`E_i=J_0'(x_i)`.  It leaves exactly

```text
D_Q=lambda(m)-mu(s),                                     (19)
```

which proves the abstract identity `(1)`.

The all-cell quantifier follows exactly as in THM-3639.  A constant
`SL_2(C)` output change puts every rank-two cell into

```text
grad(U)=(1,0),         grad(V)=(0,R),
(A_0,B_0)=(0,1),       R!=0,                             (20)
```

and every rank-one cell into the same form with `R=0`.  The calculation of
`(1)` is polynomial in `R`, so it does not divide by the tangent rank.

Ordinarity is load-bearing for the representative quantifier.  Three
distinct tangent lines imply

```text
ker(gamma) subset m_p^3,                                 (21)
```

so a vertical derivative of a representative change has branch order at
least two.  At a nonordinary tangent collision `(21)` can fail.  The formulas
below therefore state theorem candidates only on their displayed ordinary
loci; no representative-independent result is claimed on the omitted
boundaries.

## 3. Generic finite projective chart

For finite `rho!=+-4`, put

```text
M=3-u rho,                         K=9-u rho,             (22)

v_-=-(18+rho(2u+9))/(4+rho),
v_+=(18+rho(2u-9))/(4-rho).                              (23)
```

The nonzero cell condition is `M!=0`.  Use

```text
y=(4c-rho z)/(4M),                  eta=z.                (24)
```

Then

```text
h=(-4K/(4+rho), 4u, 4K/(4-rho)),

h_0-h_-=4(4u+9)/(rho+4),
h_+-h_0=4(9-4u)/(4-rho),
h_+-h_-=32K/(16-rho^2).                                 (25)
```

Thus this chart is ordinary exactly when

```text
K(4u-9)(4u+9)!=0.                                       (26)
```

The rows in `(1)` are

```text
lambda=((rho+4)(9-4u)/(8K),
        -1,
        (4-rho)(9+4u)/(8K)),

mu=((rho+4)^2(9-4u)/(16MK),
    rho/M,
    -(rho-4)^2(9+4u)/(16MK)).                           (27)
```

The curvature debt is

```text
D_Q=-N/(32 K M^2),                                      (28)

N=(rho+4)^3(9-4u)r_-
  -32 rho^2 K r_0
  +(4-rho)^3(9+4u)r_+
  +27(576-320rho u+225rho^2-9rho^3u).                  (29)
```

Under `J_0=1`, zero retained debt is exactly the affine curvature equation

```text
N=0.                                                     (30)
```

## 4. Exceptional finite charts

The projective values `rho=+-4` are not obtained by substituting into
`(23)`: each has a free side slope and may still be ordinary.

### 4.1 The `rho=4` chart

Here

```text
u=9/4,                 v_-=-9,                 v_+=v,    (31)
```

and the chart is ordinary exactly when `v notin {9,9/2}`.  With
`y=(z-c)/6`, `eta=z`, one has

```text
lambda=((2v-9)/(2(v-9)), -1, -9/(2(v-9))),
mu=(-(2v-9)/(3(v-9)), -2/3, 0),                         (32)

D_Q=[32(v-9)r_0-16(2v-9)r_-+27(117-29v)]/[72(v-9)].    (33)
```

### 4.2 The `rho=-4` chart

Here

```text
u=-9/4,                v_+=9,                  v_-=v,    (34)
```

and the chart is ordinary exactly when `v notin {-9,-9/2}`.  With
`y=-(c+z)/6`, `eta=z`, one has

```text
lambda=(9/(2(v+9)), -1, (2v+9)/(2(v+9))),
mu=(0, 2/3, (2v+9)/(3(v+9))),                           (35)

D_Q=[32(v+9)r_0-16(2v+9)r_+-27(29v+117)]/[72(v+9)].    (36)
```

In each exceptional chart, zero debt is the vanishing of the displayed
numerator.

## 5. Infinite projective chart

For `rho=infinity`,

```text
u!=0,             v_-=-9-2u,             v_+=9-2u.      (37)
```

Use `y=z/(4u)`, `eta=c`.  The triple is ordinary exactly when
`u!=+-9/4`, and

```text
lambda=((4u-9)/(8u), -1, (4u+9)/(8u)),

mu=((9-4u)/(16u^2), -1/u, -(4u+9)/(16u^2)),             (38)

D_Q=[(9-4u)r_-+32u r_0-(9+4u)r_+-243u]/(32u^3).        (39)
```

Again zero debt is exactly the vanishing of the numerator in `(39)`.

## 6. Principal chart and inheritance from THM-3624

At `rho=0`, equations `(23)` give

```text
v_-=-9/2,                         v_+=9/2.               (40)
```

The generic identity reduces to

```text
lambda_u(P)=((9-4u)P(-1)-18P(0)+(9+4u)P(1))/18,

lambda_u(J_2)=D_Q+(9-4u)J_0'(-1)/27
                    -(9+4u)J_0'(1)/27,                 (41)

D_Q=-2((9-4u)r_-+(9+4u)r_++243)/81.                   (42)
```

Therefore

```text
D_Q=0  iff  (9-4u)r_-+(9+4u)r_+=-243.                 (43)
```

Equation `(43)` is **exactly** the degree-two arbitrary-target-two-form
invoice of THM-3624.  That theorem had already located the principal
curvature boundary and supplied a finite enlarged-form survivor.  The new
content of this candidate is the decomposable all-cell identity `(1)`, its
generic and exceptional projective atlas, and the identification of `(43)`
as the zero set of the decomposable retained `J_2` quotient.

For the minimal THM-3635 polynomial

```text
Q_1=x^5+(9/2)x^4-2x^3-(27/4)x^2+x-3/4,                (44)
```

one has `(u,r_-,r_+)=(1,65/2,97/2)`, so `(42)` returns

```text
D_(Q_1)=-2072/81,                                      (45)
```

recovering THM-3639.

## 7. Exact zero-debt polynomial controls

The THM-3624 control

```text
Q_h=1/5408 (
 44069x^11+112059x^10-154749x^9-406377x^8
+188107x^7+513081x^6-82835x^5-230931x^4
+5408x-4056)                                            (46)
```

has

```text
(v_-,u,v_+)=(-9/2,1,9/2),
(r_-,r_0,r_+)=(0,0,-243/13),                            (47)
```

and hence `D_(Q_h)=0`.  THM-3624 proves that this polynomial survives the
larger arbitrary-two-form system through total source degree four.  It does
not prove that the survivor is closed or decomposable.

A smaller control on the same slope packet is

```text
Q_6=Q_1-(259/36)x^2(x^2-1)^2

   =-(259x^6-36x^5-680x^4+72x^3+502x^2-36x+27)/36.    (48)
```

Its curvatures are

```text
(r_-,r_0,r_+)=(-451/18,-251/9,-163/18),                (49)
```

so `D_(Q_6)=0`.  The six value-and-slope conditions have a unique polynomial
of degree at most five, namely `Q_1`, whose debt is nonzero.  Thus degree six
is minimal inside this fixed value-and-slope packet.

A symmetric-endpoint-curvature control is

```text
Q_*= -x^7-(27/4)x^6+3x^5+18x^4-3x^3-(27/2)x^2+x-3/4

    =Q_even+x(1-x^2)^3,

Q_even=-(27/4)x^6+18x^4-(27/2)x^2-3/4.                 (50)
```

It has

```text
(v_-,u,v_+)=(-9/2,1,9/2),
(r_-,r_0,r_+)=(-27/2,-27,-27/2),                        (51)
```

and again `D_(Q_*)=0`.  All three non-even controls have the same ordinary
retained tangent directions as `Q_1`, with pairwise determinants

```text
(39,15,54).                                             (52)
```

For every ordinary slope cell in Sections 3--5, the zero-debt equation is
one nontrivial affine equation on `(r_-,r_0,r_+)`.  Hermite interpolation at
`-1,0,1` realizes every prescribed value, slope, and curvature packet by a
polynomial of degree at most eight.  Hence every displayed ordinary slope
cell has genuine polynomial zero-debt controls, not merely independent
formal branch jets.

## 8. Strict stopping boundary

The theorem's positive conclusion is only

```text
conditional on an actual J_0=1 cell:  lambda(J_2)=D_Q,

D_Q=0:  this one retained three-value quotient vanishes.             (53)
```

In particular, the controls `Q_h,Q_6,Q_*` do **not** by themselves provide

- an actual-ring quadruple `U,V,A,B in gamma_Q(R)` with `J_0=1`;
- a solution of the first-stable identity `J_1=0`;
- an identity `J_2=0` as a polynomial in `x`;
- an all-order decomposable target pair;
- a Keller map on `A^2`; or
- any proof or refutation of the two-dimensional Jacobian conjecture.

At the retained three-value level, `D_Q=0` removes the unique ordinary-triple
cokernel because the remaining higher restriction derivative ranges over the
tangent plane.  That local statement does not discharge global target-ring
algebraization or any source point away from the retained triple.

The first-jet tangent-collision boundaries excluded by `(26)`, `(31)`,
`(34)`, and `(37)` are also **OPEN** for a representative-independent
analogue.  `JC(2)` remains **OPEN**.

## 9. Exact verification

The companion checks the stable expansion, all generic compiler jets,
branchwise determinant `12`, complete projective slope atlas and ordinary
boundaries, exhaustive rank-one/rank-two gauge cancellation, generic and
exceptional formulas `(27)--(39)`, principal recovery `(41)--(45)`, every
displayed polynomial jet, degree-six minimality inside the fixed slope
packet, the degree-eight Hermite isomorphism, and an assertion-free AST.
Normal, optimized, and stored transcripts are required to be byte-identical
before audit promotion.

```bash
python3 04-computation/jc2_russell_cylinder_principal_noneven_curvature_debt_boundary_thm3641.py
python3 -O 04-computation/jc2_russell_cylinder_principal_noneven_curvature_debt_boundary_thm3641.py
```

The stored transcript and both execution modes are byte-identical.  The
independent hostile audit reconstructed every projective chart and confirmed
the stated ordinary-locus and zero-debt scope before promotion.
