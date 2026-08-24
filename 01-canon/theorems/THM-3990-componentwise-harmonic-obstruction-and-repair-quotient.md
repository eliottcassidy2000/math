---
id: THM-3990
title: "Componentwise harmonic obstruction and repair quotient"
status: >
  PROVED. A Laplacian correction can turn a smooth defect into a strictly
  positive third-order coefficient on a compact residual zero locus exactly
  when the defect has positive average on every connected component (after
  one common choice of sign). It can in fact flatten the defect to those
  component averages. Combined with a nonnegative quadratic coefficient and
  a uniform fourth-order remainder, this gives a compact stratified
  positivity certificate. The finite-graph analogue has the same real
  component-average quotient, but its integral cokernel may retain torsion;
  the triangle Laplacian has Smith form (1,3,0). Thus a global average is
  sufficient only on a connected residual locus, and a real Poisson audit is
  not an integral gluing audit.
source: root / Hopf-product--Hopf-S6 cross-frontier session, 2026-08-24
audit: >
  PASS (root, 2026-08-24). The proof checks both directions of the Poisson
  criterion, the two-region compactness estimate, the disconnected mixed-sign
  hostile, and the integral triangle hostile 3(1,-1,0)=L(1,-1,0). The theorem
  is independent of the two 2026 preprints which motivated it; their global
  claims remain under their separate source audits.
depends_on: []
related:
  - THM-3955-node-cotangent-normalization-kernel-and-conductor-torsion
  - THM-3957-triple-normal-crossing-cotangent-conductor-kernel-and-normalization-cokernel
  - THM-3989-cusp-log-laurent-conductor-and-nondividing-depth-reduction
---

# THM-3990 -- componentwise harmonic obstruction and repair quotient

**PROVED.** There are three exact layers: the smooth Poisson quotient, its
compact third-jet consequence, and the finite/integral analogue. The point is
not that every repair problem is linear. The point is that once the allowed
correction has been identified with the image of an operator, the surviving
class lies in its cokernel and must be tested there.

## 1. Smooth componentwise repair

Let

```text
S=S_1 disjoint_union ... disjoint_union S_c                 (1)
```

be a compact smooth Riemannian manifold without boundary, written as its
connected components, and let `w in C^infinity(S)`. Put

```text
mu_alpha=vol(S_alpha)^(-1) integral_(S_alpha) w.            (2)
```

Use the Laplacian convention for which the integral of `Delta chi` on every
closed component is zero. For a fixed sign `epsilon in {+1,-1}`, the following
are equivalent:

```text
(i)  there is chi in C^infinity(S) with
     epsilon*(w-Delta chi)>0 everywhere;

(ii) epsilon*mu_alpha>0 for every alpha.                    (3)
```

When these conditions hold, `chi` may be chosen so that

```text
w-Delta chi=mu_alpha on S_alpha.                            (4)
```

In particular, if `S` is connected, the single condition
`integral_S w !=0` permits one of the two signs. If `S` is disconnected, a
nonzero total integral is not enough: the component averages must all have
the same nonzero sign.

### Proof

If `(i)` holds, integration over `S_alpha` gives

```text
epsilon*vol(S_alpha)*mu_alpha
 =integral_(S_alpha) epsilon*(w-Delta chi)>0,               (5)
```

which proves `(ii)`. Conversely, `w-mu_alpha` has integral zero on each
component. The standard solvability statement for the scalar Laplacian on a
closed connected manifold therefore gives a smooth solution of

```text
Delta chi=w-mu_alpha                         on S_alpha.    (6)
```

The componentwise solutions combine to a smooth function on `S`, and `(4)`
and `(ii)` imply `(i)`.

The mechanism is exactly the Fredholm alternative:

```text
coker(Delta:C^infinity(S)->C^infinity(S))
     = span{1_(S_1),...,1_(S_c)}.                           (7)
```

Thus the averages are not estimates. They are the complete smooth repair
quotient.

## 2. Compact third-jet positivity

Let `X` be compact, let `rho:X->R_(>=0)` be continuous, and put

```text
S={rho=0}.                                                  (8)
```

Assume that `S` is a closed smooth submanifold as in Section 1. Let `q` be a
continuous function satisfying

```text
q|S=0,                    q>=kappa*rho                      (9)
```

for some `kappa>0`. Suppose an allowed correction indexed by
`chi in C^infinity(S)` produces a continuous third coefficient `r_chi` on
`X` whose restriction is

```text
r_chi|S=w-Delta chi.                                      (10)
```

Finally suppose the corrected family has a uniform expansion

```text
F_chi(x,s)=s^2 q(x)+s^3 r_chi(x)+R_chi(x,s),
|R_chi(x,s)|<=C_chi |s|^4                                 (11)
```

for `x in X` and sufficiently small `|s|`.

If all component averages in `(2)` have a common nonzero sign `epsilon`,
then one can choose `chi` and `t_0>0` such that

```text
F_chi(x,epsilon*t)>0       for every x in X and 0<t<=t_0. (12)
```

### Proof

Choose `chi` by `(4)`. Since the finitely many numbers `mu_alpha` have sign
`epsilon`, compactness and continuity give a neighborhood `U` of `S` and a
constant `m>0` such that

```text
epsilon*r_chi>=m                  on U.                   (13)
```

On `U`, equations `(9)--(11)` give

```text
F_chi(x,epsilon*t)>=m t^3-C_chi t^4>0                    (14)
```

for small `t>0`. On the compact set `X\U`, the function `rho` has a positive
minimum `rho_0`, while `r_chi` is bounded by some `B`. Hence

```text
F_chi(x,epsilon*t)
 >=kappa*rho_0 t^2-Bt^3-C_chi t^4>0.                     (15)
```

Shrinking `t_0` proves `(12)`.

This statement is deliberately an order-three certificate. If one component
average is zero, a fourth-order term may still repair that component. The
theorem says only that no strictly signed third-order Laplacian repair exists
there; it does not turn failure of this certificate into failure of the
underlying positivity problem.

## 3. Finite graphs and the integral sidecar

Let `G` be a finite undirected graph with Laplacian `L=D-A`, and let
`V_1,...,V_c` be its connected components. For `v in R^V` put

```text
bar(v)_alpha=|V_alpha|^(-1) sum_(i in V_alpha) v_i.       (16)
```

Then

```text
im_R(L)={z: sum_(i in V_alpha)z_i=0 for every alpha}.     (17)
```

Consequently there is `x in R^V` with

```text
v-Lx=bar(v)_alpha*1                    on V_alpha,        (18)
```

and a common strict sign can be produced exactly when all component averages
have that sign. This is the finite version of `(3)--(4)`.

Over the integers, however,

```text
coker_Z(L)=Z^c direct_sum T                              (19)
```

for a finite torsion group `T` which real linear algebra erases. The triangle
is the smallest hostile. Its Laplacian is

```text
L_C3=[[ 2,-1,-1],
      [-1, 2,-1],
      [-1,-1, 2]],                                      (20)
```

with Smith form `(1,3,0)`. The vector `v=(1,-1,0)` has zero average and is
therefore in `im_R(L_C3)`, but it is not in the integral image. Indeed,

```text
L_C3*(1,-1,0)=3*(1,-1,0),                               (21)
```

while solving `L_C3*x=v` modulo the constant kernel forces a coordinate
difference to be nonintegral: for example, the zero-sum solution is
`x=(1/3,-1/3,0)`. Thus `v` is a nonzero order-three class in the integral
cokernel.

## 4. Typed applications and boundaries

### Brendle--Hung curvature preprint

In arXiv:2608.19068v1, the claimed positive-curvature construction on
`S^2 times S^2` has a nonnegative Cheeger--Mueter background, a quadratic
coefficient which is positive off a residual two-torus `Sigma`, and a cubic
coefficient `V^(3)|Sigma`. The paper's final correction changes that
coefficient by `-Delta_Sigma chi`. Since its `Sigma` is connected, Section 1
explains why the single number

```text
mu=vol(Sigma)^(-1) integral_Sigma V^(3)                  (22)
```

is exactly the final smooth obstruction, not merely a convenient statistic.
Section 2 is the abstract compactness mechanism behind the final sign choice.

This paragraph is an interpretation of the preprint's stated construction,
not a proof of its load-bearing symbolic identities or of its main theorem.
Version 1 and its attached notebook have separate audit obligations; this
theorem does not promote them.

### Hopf/S6 manuscript and planar Jacobian cusp

The `S6` manuscript's clutch maps are integral, so their surviving data must
be computed in an integral cokernel rather than by real ranks or averages.
Likewise, THM-3989's cusp-log target shears are allowed repairs, while its
coefficient law and scalar moment are surviving constraints. These are typed
instances of the rule

```text
defect -> quotient by the actual allowed repairs -> test the surviving class.
                                                                  (23)
```

They are not applications of the graph Laplacian theorem: their operators,
domains, and consumers differ. In particular, equal Smith rank, equal Euler
characteristic, or one global scalar does not replace branch orientation,
integrality, conductor incidence, or the full Laurent coefficient system.

## 5. Failure boundary and cheapest controls

1. **Disconnected hostile.** Take two components and a function with averages
   `+1` and `-1`. No common sign is possible, even if component volumes are
   changed so that the total integral is nonzero.
2. **Zero-average hostile.** A positive fourth coefficient can repair a zero
   cubic average; failure of `(3)` is not a no-positivity theorem.
3. **Integral hostile.** Equation `(21)` is invisible over `R` but survives in
   Smith form. A continuous Poisson solve is not an integral gluing solve.
4. **Nonlinear boundary.** If admissible corrections do not form the image of
   the stated operator, its cokernel is the wrong quotient. The operator must
   be derived from the native perturbation or gluing before this theorem is
   invoked.
