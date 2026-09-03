---
id: THM-4401
title: "Rank-two Poisson Casimir leaves as punctured cubic incidences"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Every nonzero Casimir
  reduction of the three-dimensional core in Long's rank-two Poisson map is
  G_m x A^1 and is exactly the simple-marked-root locus of an explicit cubic.
  Its natural A^2 incidence completion is finite flat of degree three but
  ramifies on the deleted derivative divisor.  At Casimir value zero the
  source splits as an automorphic A^2 component plus the punctured quadratic
  Kummer component of THM-3554.  A connected fixed-(T,S) slice and the
  collision-carrying Lagrangian base change are likewise punctured.  These
  identifications do not prove or disprove JC(2), DC(2), or stable
  cancellation.
source: root + clean-room independent referee / JC2 and arXiv:2608.23777 continuation session, 2026-09-03
depends_on:
  - THM-4397-rank-two-poisson-counterexample-symplectic-gauge-equivalence
  - THM-3554-punctured-kummer-collision-surface-normal-form
  - THM-3555-catalan-thickening-universal-cubic-root-cover
  - THM-3561-rational-keller-danielewski-polynomial-completion
  - THM-2045-the-smooth-factorized-R-family-has-no-planar-jacobian-mate
related:
  - THM-1300-jacobian-counterexample-dixmier-A3-explicit
  - THM-2044-explicit-rank-two-poisson-counterexample-by-symplectic-suspension
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
  - THM-3546-invariant-graph-keller-descent-criterion
primary_source: https://arxiv.org/abs/2608.23777v1
primary_script: 04-computation/poisson_rank2_casimir_leaf_cubic_incidence_thm4401.py
primary_output: 05-knowledge/results/poisson_rank2_casimir_leaf_cubic_incidence_thm4401.out
primary_script_sha256: 9d6525790ee51f2132b3d2e96e1a39643da006a625f4a79c81a520de26985507
primary_output_sha256: 0120c6e4727fb5510a6c2ae53ed4761c8adc83bacdb773b997bfabc967362d05
independent_audit_script: 04-computation/poisson_rank2_casimir_leaf_cubic_incidence_thm4401_independent_audit.py
independent_audit_output: 05-knowledge/results/poisson_rank2_casimir_leaf_cubic_incidence_thm4401_independent_audit.out
independent_audit_script_sha256: 9ed19c0bf3c3925e25457ab0dff52ee9b089038e3600acef9338dcdc2827f3ef
independent_audit_output_sha256: ef47b91091543294e5f3c99590f7b241ef739ea7e879c1475a37182407a2c324
hash_basis: raw LF bytes
audit: >
  PASS.  The primary certificate reconstructs the core directly from the
  preprint and checks 36 exact identities, including the Laurent leaves,
  cubic incidence, discriminant, zero-level split, connected fixed-(T,S)
  slice, and rank-four planar Poisson firewall.  The clean-room audit imports
  no primary code and independently checks 51 identities, including inverse
  ring maps, marked/unmarked-root factorization, the unique triple-root
  target, quadratic-discriminant scaling, Kummer coordinate Jacobians, and
  the exact three-point zero fibre.  Normal, optimized, and fixed-hashseed
  runs are byte-identical and reproduce both frozen outputs.
---

# THM-4401 -- rank-two Poisson Casimir leaves as punctured cubic incidences

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  THIS IS AN EXACT
PUNCTURED-PLANE DESCRIPTION OF NATURAL TWO-DIMENSIONAL REDUCTIONS OF A
FOUR-DIMENSIONAL POISSON COUNTEREXAMPLE.  IT IS NOT A PLANAR KELLER MAP AND
DOES NOT SETTLE `JC(2)` OR `DC(2)`.**

Let `k` be a characteristic-zero field.  In the polynomial source
coordinates `(x,y,beta,D)` of Christopher D. Long's
[*An Explicit Counterexample to the Rank-Two Poisson
Conjecture*](https://arxiv.org/abs/2608.23777v1), put

```text
w     =1+xy,
alpha =2-3xy-x^2 beta,
R     =2x-3x^2y-x^3 beta=x alpha,
S     =y+3xw^2 beta+3xy^2(4+3xy),
T     =-(w^3 beta+y^2w(4+3xy))/2.                      (1)
```

THM-4397 identifies this presentation with the proved rank-two Poisson map
of THM-2044.  In these coordinates its point map is

```text
(x,y,beta,D) |-> (R,T,D,S),                            (2)
```

and

```text
{D,R}=1,                 {S,T}=1,
{R,S}={R,T}={D,S}={D,T}=0.                             (3)
```

The core map `(R,T,S)` has determinant

```text
det partial(R,T,S)/partial(x,y,beta)=1.                (4)
```

Thus `R` is a Casimir of the three-variable Jacobian-Poisson core.  On the
full fourfold, reduction at `R=rho` quotients the translation generated in
the free `D` direction and gives

```text
f_rho:X_rho=Spec A_rho -> A^2_(S,T),
A_rho=k[x,y,beta]/(R-rho).                              (5)
```

Equivalently, every global cross-section `D=d` of that translation induces
the same two-dimensional map `(S,T)`.

## 1. Every nonzero leaf is a Laurent plane

Fix `rho in k*`.  Since `x alpha=rho`, both `x` and `alpha` are units in
`A_rho`.  Define

```text
t=x^(-1),                 u=t(1+xy)=t+y.                (6)
```

Here `u` is a new leaf coordinate, not the temporary notation `xy` used in
some presentations of the core.  The inverse formulas are

```text
x=t^(-1),       y=u-t,       beta=5t^2-3tu-rho t^3.    (7)
```

Substitution in both directions proves the ring isomorphism

```text
A_rho ~= k[t,t^(-1),u].                                (8)
```

On this Laurent plane, `(1)` becomes

```text
S=2t+4u-3rho u^2,
T=(rho u^3-u^2-ut)/2,                                  (9)

det partial(S,T)/partial(t,u)=-t.                      (10)
```

The sign in `(10)` uses target order `(S,T)`; reversing the target order
gives `+t`.  Since `t` is a unit, `f_rho` is everywhere etale.  Its source is
not `A^2`: its coordinate ring has the nonconstant unit `t`.

## 2. Exact cubic incidence and the deleted derivative divisor

Introduce the marked-root polynomial

```text
P_rho(U)=rho U^3-2U^2+S U+4T.                          (11)
```

Equations `(9)` give

```text
P_rho(u)=0,                    P_rho'(u)=2t,
t=(S-4u+3rho u^2)/2.                                  (12)
```

More precisely, the maps defined by `(9)`, `U |-> u`, and the final formula
in `(12)` are mutually inverse and give

```text
k[t,u] ~= k[S,T,U]/(P_rho(U)).                         (13)
```

For `rho!=0`, the right side is free of rank three over `k[S,T]`, with basis
`1,U,U^2`.  Hence `(13)` is the natural finite-flat degree-three completion
of `(8)`.  The completed map has Jacobian `-t`; equivalently, a marked root
is ramified exactly when `P_rho'(u)=2t=0`.  Therefore

```text
X_rho = {t!=0}
      = the simple-marked-root locus in the cubic incidence. (14)
```

This is affinely the universal cubic root cover of THM-3555.  Indeed, put

```text
v=u-2/(3rho),
p=S/rho-4/(3rho^2)=2t/rho-3v^2,
q=4T/rho+2S/(3rho^2)-16/(27rho^3).                     (15)
```

The source and target changes in `(15)` are polynomial automorphisms for the
fixed scalar `rho in k*`, and `(11)` becomes

```text
v^3+pv+q=0,                 q=-v^3-pv.                 (16)
```

The deleted line `t=0` becomes the universal ramification parabola

```text
p+3v^2=0.                                               (17)
```

### 2.1 Point image and finiteness

Over an algebraic closure of `k`, every cubic in `(11)` other than the unique
triple-root cubic has a simple root.  Consequently the set-theoretic
geometric image of the punctured map is

```text
A^2_(S,T) minus {(4/(3rho),-2/(27rho^2))}.              (18)
```

At the missing point,

```text
P_rho(U)=rho(U-2/(3rho))^3.                             (19)
```

Off the discriminant, the map is finite etale of degree three.  Over the
smooth part of the discriminant, the double root is deleted but the third,
simple marked root remains.  The point in `(18)` is the cusp where all three
roots coincide and no simple marked root remains.  In particular, the
punctured map is not finite over the full target plane: over an algebraic
closure its point image is not closed.

Two image qualifications are load-bearing:

- the scheme-theoretic image means its closure, which is all of `A^2` because
  the map is dominant; and
- over a non-algebraically-closed field, the image on `k`-rational points may
  be smaller than `(18)`, because a cubic need not have a simple root in `k`.

### 2.2 Discriminant and the unmarked-root sidecar

The exact cubic discriminant is

```text
Disc_U(P_rho)
 =4(S^2-rho S^3+32T-36rho ST-108rho^2T^2).             (20)
```

Its pullback to `(t,u)` factors as

```text
4t^2(9rho^2u^2-8rho t-12rho u+4).                      (21)
```

After pulling the incidence polynomial back to `k[t,u][U]`, division by its
marked factor gives

```text
P_rho(U)/(U-u)
 =rho U^2+(rho u-2)U+2(t+u-rho u^2).                   (22)
```

The second factor in `(21)` is exactly the discriminant of the quadratic in
`(22)`.  It records collision of the two *unmarked* roots.  It is not a
second ramification divisor for the marked source point: when `t!=0`, the
marked derivative remains `2t`, even if `(22)` has a double root.  On `t=0`
that quadratic discriminant becomes `(3rho u-2)^2`; its simultaneous zero
with the marked ramification line is the triple-root point `(18)`.

## 3. The zero leaf is an automorphic plane plus a Kummer leaf

At `rho=0`, the factorization `R=x alpha` is comaximal because

```text
1=alpha/2+x(3y+x beta)/2.                              (23)
```

Chinese remainders therefore give the scheme-level decomposition

```text
A_0 ~= k[y,beta] x k[t,t^(-1),u],
X_0 ~= A^2 disjoint-union (G_m x A^1).                 (24)
```

On the component `x=0`, the reduced map is

```text
(y,beta) |-> (S,T)=(y,-beta/2-2y^2).                   (25)
```

It is a polynomial automorphism with inverse

```text
y=S,                    beta=-2T-4S^2,                 (26)
```

and its `(S,T)` Jacobian is `-1/2`.

On the component `alpha=0`, the same Laurent coordinates as `(6)` give

```text
S=2t+4u,                T=-(u^2+ut)/2,
det partial(S,T)/partial(t,u)=-t.                      (27)
```

The marked-root equation has dropped to the quadratic

```text
P_0(U)=-2U^2+S U+4T.                                   (28)
```

Its **standard quadratic discriminant** and pullback are

```text
Delta_2=S^2+32T,                 Delta_2|source=4t^2.  (29)
```

There is a harmless but important scaling distinction.  Substituting
`rho=0` into the cubic discriminant formula `(20)` gives

```text
Disc_cubic_formula|_(rho=0)=4 Delta_2=16t^2,           (30)
```

not the standard quadratic discriminant itself; the leading cubic
coefficient has vanished.

Set

```text
b=S,                       delta=S^2+32T.              (31)
```

Using `(t,b)` as Laurent source coordinates turns `(27)` into the exact
THM-3554 Kummer normal form

```text
(t,b) |-> (b,4t^2).                                    (32)
```

The Jacobian of `(32)` is `-8t`.  This agrees with `(27)`: the source change
`(t,u)->(t,b)` has determinant `4`, while the target change
`(S,T)->(b,delta)` has determinant `32`.  Thus the curved component is finite
etale of degree two over `delta!=0`, and adjoining `t=0` restores its
quadratic ramification.

At the distinguished target `(R,T,S)=(0,1/8,0)`, the complete reduced fibre
is

```text
(x,y,beta)=(0,0,-1/4),
             (1,-3/2,13/2),
             (-1,3/2,13/2).                            (33)
```

The first point lies on the automorphic plane component.  The other two are
the Kummer pair `t=+/-1`.  Hence the zero reduction retains the exact `1+2`
collision, but its source is disconnected; its only `A^2` component is
invertible, and its noninjective component is Laurent.

## 4. A connected punctured triple slice

There is a complementary connected two-dimensional slice.  Fix
`(T,S)=(tau,sigma)` with `tau in k*`.  The formula for `T` in `(1)` makes
`w=1+xy` a unit on the core curve.  Put

```text
h=x/w,                   d=1-sigma h-6tau h^2.         (34)
```

Direct elimination and its inverse give

```text
y=sigma+6tau h,                         wd=1,
R=2h-sigma h^2-4tau h^3,                x=h/d,
beta=-2tau d^3-4y^2d^2-3hy^3d.          (35)
```

Therefore the full slice, including the free source coordinate `D`, is

```text
Spec k[h,d^(-1),D] -> A^2_(R,D),
(h,D) |-> (2h-sigma h^2-4tau h^3,D),
det partial(R,D)/partial(h,D)=2d.                      (36)
```

It is etale because `d` is a unit, but it is not an affine plane because
`d` is a nonconstant unit.  At `(tau,sigma)=(1/8,0)`,

```text
d=1-3h^2/4,
R=2h-h^3/2=-h(h-2)(h+2)/2.                             (37)
```

The three retained roots `h=0,+/-2` give the collision, while the missing
points are `h=+/-2/sqrt(3)` over an algebraic closure.  The unique extension
of `(36)` to the chosen affine completion `A^2_(h,D)` has Jacobian `2d` and
ramifies on the two restored vertical lines.  This is a connected punctured
triple cover, not a planar Keller counterexample.

## 5. Descent and completion firewalls

The preceding identities separate four operations that must not be
conflated.

### 5.1 Literal four-output Poisson descent is impossible

In target order `(R,T,D,S)`, the required bracket matrix in `(3)` is

```text
[ 0  0 -1  0 ]
[ 0  0  0 -1 ]
[ 1  0  0  0 ]
[ 0  1  0  0 ],                                        (38)
```

which has determinant one and rank four.  For four functions on a planar
symplectic algebra, the bracket matrix is `J Omega_2 J^T` with `J` of size
`4 x 2`, so its fraction-field rank is at most two.  Thus all four canonical
outputs cannot simply be realized inside one polynomial plane.  The genuine
reductions above first remove or fix a canonical pair.

### 5.2 The chosen incidence completion is unique and ramified

Once the Laurent chart `(t,u)` and the output functions `(9)` are fixed,
their polynomial extension to `k[t,u]` is unique: two such extensions that
agree on the dense open set `t!=0` agree everywhere.  Its Jacobian is `-t`.
Thus direct filling of this chosen puncture restores the marked-root
ramification divisor and cannot give a Keller map.

This statement does **not** exclude a different affine modification, a
different polynomial-plane model of the function field, a deformation that
changes the cubic extension, or a construction that sends the ramified
valuation to infinity.

### 5.3 The collision-carrying Lagrangian base change is punctured

In the target of `(2)`, the plane

```text
L={D=0,S=0},                       coordinates (R,T),  (39)
```

is Lagrangian.  Its source inverse image is exactly `{D=0,S=0}`.  Since the
ambient map is etale, this scheme-theoretic base change is etale and contains
the three-point fibre `(33)`.  Under the identification

```text
(R,T,S)=(F_3,-F_1/2,F_2)                              (40)
```

with the THM-1300 core, it is precisely the middle-coordinate-zero surface
analyzed in THM-3561.  That theorem proves that the source is an affine plane
minus a nonempty curve.  Its rational planar Keller pair retains the triple
collision, while its polynomial filling maps to a non-plane Danielewski
surface.

Thus even the canonical collision-carrying Lagrangian base change is a
punctured planar near-counterexample, not a map `A^2->A^2`.  An arbitrary
Lagrangian restriction does not inherit a planar Jacobian determinant merely
because the pulled-back two-form vanishes; the inverse-image/base-change
condition in `(39)` is the decisive extra structure.

### 5.4 Stable equivalence does not imply planar cancellation

The preprint gives an ordinary polynomial factorization

```text
Phi_pt = target permutation o (G x id_A1) o source automorphism, (41)
```

where `G=(R,T,S)` is the three-dimensional Keller core.  THM-4397 further
proves that this four-dimensional map and THM-2044 are symplectically
right-left equivalent.  Neither statement supplies an equivalence

```text
G ~= F x id_A1                                      (42)
```

for a planar map `F`.  If `(42)` existed, then `(41)` would make the full map
equivalent to `F x id_A2`; constant Jacobian and noninjectivity would descend
to `F`, producing a planar Jacobian counterexample.  No stable-cancellation
theorem and no equivalence `(42)` is proved here.

THM-2045 separately shows that, in the original `(x,q)` variables, the
displayed first coordinate

```text
R=x(2-3xq)                                             (43)
```

has no polynomial planar Jacobian mate.  This excludes the literal
coordinate-preserving destabilization, not an arbitrary left-right
equivalence or an unrelated planar Keller map.

Finally, the reductions `(8)`, `(24)`, `(36)`, and `(39)` do not satisfy the
polynomial-coordinate-hypersurface hypothesis of THM-3546: they are Laurent,
disconnected, or punctured.  Their failures do not exclude another
noncanonical invariant surface.

## 6. Transferable mechanism and scope

The four-dimensional Poisson construction does not cancel the universal
cubic's finite ramification inside an affine plane.  It makes the derivative
a unit by deleting its divisor:

```text
marked-root derivative =2t,
source leaf             ={t!=0},
chosen completion Jac   =-t.                           (44)
```

For planar JC, this isolates a precise construction target: export the
marked-root ramification to nonproper escape at infinity while keeping the
source ring a polynomial plane, hence without introducing nonconstant units.
The present map achieves the first objective only on a Laurent source.

This theorem proves no assertion about arbitrary Casimir choices, arbitrary
Lagrangian surfaces, alternate affine modifications, quantization to `A_2`,
`DC(2)`, or `JC(2)`.  Long's Appendix-B Weyl construction remains in `A_4`,
and THM-4397's exact gauge equivalence does not change that rank.

## 7. Exact audit and replay

The primary certificate reconstructs `(1)` directly from the preprint and
checks 36 exact gates.  It verifies the defining formula identities in
`(6)--(37)`, the discriminant, the fixed-slice inverses, and the determinant
of `(38)`.  The clean-room audit
was frozen before comparison with the primary.  It imports no theorem
candidate and independently verifies 51 gates, including both ring-map
directions in `(8)` and `(13)`, the quotient quadratic `(22)`, the geometric
triple-root target, the CRT identity, all discriminant scalings, Kummer
coordinate-change Jacobians, and the three points `(33)`.

Reproduce the frozen outputs with

```bash
python3 -B 04-computation/poisson_rank2_casimir_leaf_cubic_incidence_thm4401.py
python3 -B -O 04-computation/poisson_rank2_casimir_leaf_cubic_incidence_thm4401.py
PYTHONHASHSEED=314159 python3 -B \
  04-computation/poisson_rank2_casimir_leaf_cubic_incidence_thm4401.py

python3 -B \
  04-computation/poisson_rank2_casimir_leaf_cubic_incidence_thm4401_independent_audit.py
python3 -B -O \
  04-computation/poisson_rank2_casimir_leaf_cubic_incidence_thm4401_independent_audit.py
PYTHONHASHSEED=271828 python3 -B \
  04-computation/poisson_rank2_casimir_leaf_cubic_incidence_thm4401_independent_audit.py
```

The normal, optimized, and fixed-hashseed runs are byte-identical within each
certificate and match the stored outputs.  The raw-LF hashes are

```text
primary script  9d6525790ee51f2132b3d2e96e1a39643da006a625f4a79c81a520de26985507
primary output  0120c6e4727fb5510a6c2ae53ed4761c8adc83bacdb773b997bfabc967362d05
audit script    9ed19c0bf3c3925e25457ab0dff52ee9b089038e3600acef9338dcdc2827f3ef
audit output    ef47b91091543294e5f3c99590f7b241ef739ea7e879c1475a37182407a2c324
```
