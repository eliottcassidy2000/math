---
id: THM-2191
title: "Catalytic localization of the Gordian metric"
status: >
  PROVED (abstract metric-monoid localization, maximality, diagonal-kernel
  formula, catalytic-capacity identity, amplification inequality, and
  homogenization bounds) + CITED APPLICATION (classical prime-decomposition
  cancellation and Owens--Strle's singular-concordance computations for
  knots). The common-translation envelope of any
  commutative nonexpansive metric monoid is its greatest
  translation-invariant pseudometric below the original metric. For knots,
  integrality and connected-sum cancellation make this envelope an attained
  genuine metric d_cat. Its root norm u_cat is the minimum diagonal of
  THM-2176's min-plus kernel, and u(K)-u_cat(K) is exactly the maximum
  directional catalytic saving. Moreover
  max{u_hash(K),c*(K)} <= u_cat(K) <= u(K), with
  sup_phi |phi(K)| <= u_hash(K).
  Hence positive translation catalysis forces a strict connected-sum
  homogenization gap; in particular a one-crossing saving for a knot of
  unknotting number three forces u_hash<=2. Owens--Strle's equality
  c*(9_10)=u(9_10)=3 closes the formerly unresolved reverse 4_1/9_10
  direction as pure bypass. The conditional 10_6 seed and uncalibrated
  pretzel families remain; no positive knot catalyst is produced here.
source: codex-2026-07-24-catalytic-gordian-localization
depends_on:
  - THM-2176-gordian-continuation-profile-and-interaction-cocycle
related:
  - THM-2183-order-join-is-an-exact-tournament-metric-product
  - THM-2188-finite-phase-bank-and-pairwise-overlap-no-go
  - THM-2195-transitive-quotients-exactly-control-universal-substitution-products
external:
  - "Horst Schubert, Die eindeutige Zerlegbarkeit eines Knotens in Primknoten, Sitzungsberichte der Heidelberger Akademie der Wissenschaften 1949/3, 57--104, DOI 10.1007/978-3-642-45813-2."
  - "Brendan Owens and Saso Strle, Immersed disks, slicing numbers and concordance unknotting numbers, Communications in Analysis and Geometry 24 (2016), 1107--1138, DOI 10.4310/CAG.2016.v24.n5.a8, arXiv:1311.6702."
---

# THM-2191 -- catalytic localization of the Gordian metric

THM-2176 separates strict connected-sum subadditivity into two mechanisms:
translation catalysis and geodesic bypass. Its published knot audit finds
many bypasses but no certified catalyst. The right next object is not another
binary relation. It is the metric obtained by allowing the same auxiliary
summand on both sides and then optimizing it away.

The construction is universal, and its comparison with connected-sum
homogenization turns any future catalytic witness into an asymptotic
self-power shortcut.

## 1. Common-translation localization

Let `(M,+,0)` be a commutative monoid with a metric `d` satisfying

```text
d(a+c,b+d)<=d(a,b)+d(c,d).                           (1)
```

In particular every simultaneous translation is nonexpansive:

```text
d(a+c,b+c)<=d(a,b).                                  (2)
```

Define

```text
d_cat(a,b)=inf_(c in M) d(a+c,b+c).                  (3)
```

> **Localization theorem.** The function `d_cat` is a
> translation-invariant pseudometric, satisfies the same joint
> nonexpansivity law as `d`, and is the greatest translation-invariant
> pseudometric bounded above by `d`.

### Pseudometric and triangle inequality

Nonnegativity, symmetry, and `d_cat(a,a)=0` follow immediately from (3).
For arbitrary catalysts `r,s`, (2) and the triangle inequality give

```text
d_cat(a,e)
 <=d(a+r+s,e+r+s)
 <=d(a+r+s,b+r+s)+d(b+r+s,e+r+s)
 <=d(a+r,b+r)+d(b+s,e+s).                            (4)
```

Taking the infimum independently over `r` and `s` proves

```text
d_cat(a,e)<=d_cat(a,b)+d_cat(b,e).                   (5)
```

### Exact translation invariance

For any `t`, translating a candidate in (3) and using (2) gives

```text
d_cat(a+t,b+t)<=d_cat(a,b).                          (6)
```

Conversely, the catalysts of the form `t+c` are a subset of all catalysts,
so

```text
d_cat(a+t,b+t)
 =inf_c d(a+(t+c),b+(t+c))
 >=inf_w d(a+w,b+w)
 =d_cat(a,b).                                        (7)
```

Thus

```text
d_cat(a+t,b+t)=d_cat(a,b).                           (8)
```

Commutativity is load-bearing here. In a noncommutative monoid, a one-sided
version must retain the side of the context, exactly as in MISTAKE-249's
correction to THM-2176.

### Compatibility with addition

For arbitrary `r,s`, (1) gives

```text
d_cat(a+c,b+d)
 <=d(a+c+r+s,b+d+r+s)
 <=d(a+r,b+r)+d(c+s,d+s).                            (9)
```

Taking the two infima yields

```text
d_cat(a+c,b+d)<=d_cat(a,b)+d_cat(c,d).               (10)
```

Consequently `d_cat=0` is a monoid congruence, and the metric quotient is a
metric monoid on which every translation is an isometry.

### Universal maximality

Let `rho` be any translation-invariant pseudometric with `rho<=d`.
For every `c`,

```text
rho(a,b)=rho(a+c,b+c)<=d(a+c,b+c).                   (11)
```

Taking the infimum over `c` proves

```text
rho(a,b)<=d_cat(a,b).                                (12)
```

This proves maximality. The localization does not guess a preferred
catalyst: it is forced by the universal property.

## 2. The catalytic root norm and THM-2176's kernel

Put

```text
ell(x)=d(x,0),
ell_cat(x)=d_cat(x,0).                               (13)
```

For THM-2176's faithful min-plus kernel

```text
P_x(a,b)=d(x+a,b),                                   (14)
```

equation (3) gives the exact diagonal formula

```text
ell_cat(x)
 =inf_y d(x+y,y)
 =inf_y P_x(y,y).                                    (15)
```

Thus `ell_cat` is the infimal diagonal, or min-plus trace, of the full
operation kernel.

Recall THM-2176's directional catalytic term

```text
C_y(x)=ell(x)-d(x+y,y)>=0.                           (16)
```

Equations (15)--(16) imply

```text
ell(x)-ell_cat(x)=sup_y C_y(x).                      (17)
```

Therefore the following are equivalent:

```text
some y has C_y(x)>0;
ell_cat(x)<ell(x);
the identity context is not a diagonal minimizer of P_x.                  (18)
```

No attainment assumption is needed for (18): if an infimum is strictly
below `ell(x)`, some member of the defining set is already strictly below
`ell(x)`.

Inside the localized metric every directional catalytic term vanishes:

```text
d_cat(x+y,y)=d_cat(x,0)=ell_cat(x).                  (19)
```

The interaction defect

```text
ell_cat(x)+ell_cat(y)-ell_cat(x+y)                   (20)
```

can nevertheless remain positive. In THM-2176's split, the localization
removes translation catalysis but does not force geodesic bypass to vanish.

## 3. A catalyst amplifies to every self-power

Let

```text
r=d(x+y,y).                                          (21)
```

For every positive integer `n`, translate the path `x+y -> y` successively
by `0,x,...,(n-1)x`, then concatenate the translated paths in descending
order from `nx+y` to `y`. Each translated step costs at most `r`, so

```text
d(nx+y,y)<=nr.                                       (22)
```

The two end attachments cost at most `ell(y)` each:

```text
d(nx,nx+y)<=ell(y),
d(y,0)=ell(y).                                       (23)
```

Hence

```text
ell(nx)<=nr+2ell(y).                                 (24)
```

Subadditivity of `ell(nx)` gives the connected-sum homogenization

```text
ell_hash(x)=lim_(n->infinity) ell(nx)/n.             (25)
```

Dividing (24) by `n` and then optimizing over `y` proves

```text
ell_hash(x)<=ell_cat(x)<=ell(x).                     (26)
```

This is stronger than the pointwise calibration obstruction in THM-2176:
every finite catalytic saving propagates to a linear asymptotic saving.
Specifically, if

```text
d(x+y,y)<=ell(x)-s,                                  (27)
```

then

```text
ell(nx)<=n(ell(x)-s)+2ell(y),
ell_hash(x)<=ell(x)-s.                               (28)
```

The catalyst may be complicated, but its cost is only the bounded endpoint
term `2ell(y)`.

## 4. Additive invariants and the four-level obstruction chain

Let `phi:M->R` be additive and 1-Lipschitz for `d`. For every catalyst `c`,

```text
|phi(a)-phi(b)|
 =|phi(a+c)-phi(b+c)|
 <=d(a+c,b+c).                                       (29)
```

Therefore `phi` is also 1-Lipschitz for `d_cat`:

```text
|phi(a)-phi(b)|<=d_cat(a,b).                         (30)
```

Applying the same invariant to `nx` and dividing by `n` also gives

```text
|phi(x)|<=ell_hash(x).                               (31)
```

Combining (26) and (31),

```text
sup_(additive 1-Lipschitz phi) |phi(x)|
 <=ell_hash(x)
 <=ell_cat(x)
 <=ell(x).                                           (32)
```

If one additive invariant calibrates `x`, meaning `|phi(x)|=ell(x)`, every
inequality in (32) is equality. This recovers THM-2176's conclusion that a
calibrated first summand cannot be catalyzed.

If `d` is integer-valued and `ell(x)=m`, positive catalysis forces

```text
ell_cat(x)<=m-1,
ell_hash(x)<=m-1.                                    (33)
```

Thus the strict inequality

```text
ell_hash(x)>m-1                                      (34)
```

rules out every positive catalyst, even when no single additive invariant
calibrates `x`.

## 5. Knot specialization

Let `M` be oriented knot types under connected sum and let `d=d_G` be
Gordian distance. Define

```text
d_cat(K,L)=min_J d_G(K#J,L#J),
u_cat(K)=d_cat(K,U).                                 (35)
```

The infimum in (3) is a minimum here: the nonempty set of possible Gordian
distances is a set of nonnegative integers. If `d_cat(K,L)=0`, a minimizing
`J` satisfies

```text
K#J is isotopic to L#J.                              (36)
```

Classical unique prime decomposition of knots cancels `J` from (36), so
`K` and `L` are isotopic. Hence:

> **Catalytic Gordian metric.** `d_cat` is a genuine integer-valued metric
> on knot types, connected sum acts by isometries, and `d_cat` is the
> greatest connected-sum-translation-invariant metric below `d_G`.

No cancellativity hypothesis was needed for the abstract pseudometric
theorem. It enters only here, to prove nondegeneracy.

The minimum-diagonal and capacity formulas become

```text
u_cat(K)=min_J d_G(K#J,J)=min_J P_K(J,J),             (37)

u(K)-u_cat(K)=max_J C_J(K).                           (38)
```

Thus an actual positive knot catalyst exists exactly when

```text
u_cat(K)<u(K).                                       (39)
```

Mirror is an isometric monoid automorphism and `J -> mirror(J)` permutes all
catalysts, so

```text
d_cat(mirror(K),mirror(L))=d_cat(K,L),
u_cat(mirror(K))=u_cat(K).                           (40)
```

### The singular-concordance metric is a second universal floor

Owens--Strle define

```text
d_*(K,L)
 =minimum number of double points in a normally immersed
  concordance from K to L,                           (40a)

c*(K)=d_*(K,U).                                      (40b)
```

Owens--Strle's orientation convention identifies an immersed concordance
from `K_1` to `K_2` with an immersed disk bounded by `-K_1#K_2`. Therefore

```text
d_*(K_1,K_2)=c*(-K_1#K_2),
d_*(K,U)=c*(-K)=c*(K),                               (40b')
```

where the last equality follows because concordance inversion preserves the
unsigned number of double points. This identifies (40b) with their
four-ball crossing number without an orientation ambiguity.

A crossing-change trace is such an immersed concordance, so

```text
d_*(K,L)<=d_G(K,L).                                  (40c)
```

The function `d_*` is a pseudometric on knot types and a metric on smooth
concordance classes. It is exactly connected-sum-translation invariant.
One direction follows by adding the product annulus of a common summand.
For the reverse direction, an immersed concordance between `K#J` and `L#J`
may be summed with the oriented concordance inverse `-J` (the mirror with
reversed orientation); the ribbon, hence slice, knot `J#(-J)` is then
removed by smooth concordances at both ends, without adding double points.
Hence

```text
d_*(K#J,L#J)=d_*(K,L).                               (40d)
```

By the universal maximality (12),

```text
d_*(K,L)<=d_cat(K,L),
c*(K)<=u_cat(K).                                     (40e)
```

This argument is more general than a real-valued calibration: every
translation-invariant pseudometric below Gordian distance supplies a lower
floor for `d_cat`.

For the connected-sum homogenized unknotting number `u_hash` of THM-2176,
(32) and (40e) read

```text
sup_phi |phi(K)|<=u_hash(K),

max{u_hash(K),c*(K)}
 <=u_cat(K)
 <=u(K),                                             (41)
```

where the supremum ranges over additive real-valued Gordian
1-Lipschitz invariants.

All fully additive-calibrated published symbiont families audited in
THM-2176 have

```text
u_hash(K)=u_cat(K)=u(K)                              (42)
```

in each calibrated direction, so their strict connected-sum defects remain
pure bypass.

### The apparent `9_10` catalyst is already impossible

THM-2176 left the reverse direction of the `4_1,9_10` seed open because its
ordinary signature and Rasmussen invariant give only the lower bound two.
Owens--Strle prove the stronger equality

```text
c*(9_10)=u(9_10)=3.                                  (43)
```

Equations (40e) and (43) force

```text
3=c*(9_10)<=u_cat(9_10)<=u(9_10)=3,
u_cat(9_10)=3.                                       (44)
```

Since adding any common summand never costs more than `u(9_10)`, (44)
actually gives the uniform identity. Equivalently, for every `J` the direct
sandwich is

```text
3=d_*(9_10,U)
 =d_*(9_10#J,J)
 <=d_G(9_10#J,J)
 <=d_G(9_10,U)=3.                                  (45)
```

In particular,

```text
d_G(9_10#4_1,4_1)=3,
C_(4_1)(9_10)=0.                                     (46)
```

Thus the `4_1,9_10` symbiosis is pure geodesic bypass in both directions.
This closes the only unconditional reverse candidate in THM-2176's
published seed table.

Owens--Strle's table gives `c*(10_6)=2`; the newer source audited in
THM-2176 retains `u(10_6) in {2,3}`. Therefore the conditional
`3_1,10_6` reverse seed is the sharp remaining boundary:

```text
u(10_6)=2  => c*(10_6)=u_cat(10_6)=u(10_6),

u(10_6)=3 and d_G(10_6#3_1,3_1)=2
  => c*(10_6)=u_cat(10_6)=u_hash(10_6)=2.            (47)
```

The second line uses (28), (40e), and `u(3_1)=1`: the upper bound
`u_hash(10_6)<=2` comes from amplification, while the opposite bound is
supplied by the additive signature calibration `|sigma(10_6)|/2=2`.
More generally, every one-crossing catalyst for a knot with unknotting
number three forces its entire self-power sequence below slope three:

```text
u(nK)<=2n+O(1).                                      (48)
```

The next decisive knot targets are therefore to decide `u(10_6)`, determine
the displayed translated distance if its value is three, and test the
uncalibrated pretzel families against `c*` or another
translation-invariant metric floor.

## 6. Hostile boundary examples

The hypotheses and converses above are sharp.

### Cancellativity does not prevent collapse

Take `M=Z_(>=0)` under addition and

```text
d(m,n)=|2^(-m)-2^(-n)|.                              (47)
```

This is a metric. Writing `x=2^(-m)`, `y=2^(-n)`,
`u=2^(-p)`, and `v=2^(-q)`, all four numbers lie in `[0,1]`, and

```text
|xu-yv|
 <=u|x-y|+y|u-v|
 <=|x-y|+|u-v|
 =|2^(-m)-2^(-n)|+|2^(-p)-2^(-q)|,                  (48)
```

so (1) holds. But

```text
d(m+c,n+c)=2^(-c)d(m,n),
d_cat(m,n)=0                                        (49)
```

for every `m,n`. The monoid is cancellative, yet the localization collapses
and its off-diagonal infimum is not attained.

### A homogenization gap is not sufficient for catalysis

Again take `M=Z_(>=0)`, now with the discrete metric. Every translation is an
isometry, hence

```text
d_cat=d.                                             (50)
```

For `x=1`,

```text
ell(n)=1,
ell_hash(1)=0<1=ell_cat(1)=ell(1).                   (51)
```

Thus a strict homogenization gap is necessary for catalysis by (26), but is
not sufficient. This example also shows that additive real-valued
1-Lipschitz invariants need not recover `d_cat`: additivity forces
`phi(n)=n phi(1)`, while discrete 1-Lipschitzness forces `phi(1)=0`.

### Integrality without cancellation is not enough

Let `M={0,a}` with `a+a=a` and the discrete metric. Addition is
nonexpansive, but translation by the absorbing element sends both points to
`a`, so

```text
d_cat(0,a)=0.                                        (52)
```

In the knot application, integrality supplies attainment and cancellation
then prevents collapse.

### The LRC safe-set monoid has a null absorber

The same construction gives an exact transfer boundary for LRC. Let the
commutative monoid consist of finite speed sets under union, and put

```text
d_mu(A,B)=measure(G_A symmetric_difference G_B).    (LRC.1)
```

This is a pseudometric, and union is nonexpansive because

```text
G_(A union C)=G_A intersection G_C.                 (LRC.2)
```

In fact the full joint law (1) holds: the symmetric difference of
`G_A intersection G_C` and `G_B intersection G_D` is contained in the union
of `G_A symmetric_difference G_B` and `G_C symmetric_difference G_D`.

However, the complete arithmetic-progression packet

```text
P={1,2,...,13}
```

has `measure(G_P)=0`. Therefore, for every `A,B`,

```text
d_mu(A union P,B union P)=0,
inf_C d_mu(A union C,B union C)=0.                  (LRC.3)
```

So unrestricted common-context localization collapses identically on the
positive-Haar safe-measure observable. The obstruction is not
noncancellativity alone but the presence of a null absorbing context for that
a.e. predicate. This does not collapse weak nonemptiness or THM-2047's Euler
data: isolated AP points are invisible to `d_mu`. A useful measure-level LRC
quotient must be graded by remaining runner slots or forbid already terminal
contexts; THM-2174's fixed-core two-state first-period flag is one such
restricted carrier. Catalytic localization itself supplies no LRC pump or
proof.

## 7. Typed frontier

The construction changes the knot search from an unstructured pair hunt to
three nested targets:

```text
source:       Gordian metric plus connected sum;
operation:    add one common catalyst to both endpoints;
envelope:     greatest translation-invariant metric below d_G;
one-body norm:u_cat(K)=min_J P_K(J,J);
dual floor:   additive 1-Lipschitz invariants;
homog. floor: connected-sum homogenization u_hash;
metric floor: four-ball crossing number c*;
catalysis:    strict gap u_cat(K)<u(K);
cheap no-go:  max{u_hash(K),c*(K)}>u(K)-1;
closed seed:  c*(9_10)=u(9_10)=3;
first target: conditional 10_6 boundary, then uncalibrated pretzels.  (53)
```

This theorem does not produce a positive knot catalyst, compute
`u_cat(10_6)` in the possible `u(10_6)=3` branch, or claim that
`u_hash=u_cat` generally. It deliberately avoids the phrase “stable
unknotting number,” which is used for a different twist-family construction
in the literature. “Catalytic Gordian metric” is descriptive repo-local
terminology, not a priority claim; related abstract constructions occur in
the theory of subinvariant metrics, metric monoids, and group completion.
Likewise, “singular-concordance metric” is descriptive here; Owens--Strle
call `d_*` the crossing number distance.

QED.
