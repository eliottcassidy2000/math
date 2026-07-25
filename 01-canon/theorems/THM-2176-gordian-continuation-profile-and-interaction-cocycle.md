---
id: THM-2176
title: "Gordian continuation profile, rooted order, min-plus kernel, and interaction cocycle"
status: >
  PROVED (abstract metric-monoid, universal continuation, cocycle,
  rooted-order, catalysis/bypass, and homogenization statements) +
  CITED APPLICATION (Brittenham--Hermiller crossing-change certificates and
  the classical signature lower bound). The connected-sum continuation
  profile is the largest connected-sum congruence contained in equality of
  unknotting number. The interaction defect is a normalized symmetric
  nonnegative 2-coboundary, monotone on product Gordian cones. The
  T(2,7)-mirror counterexample is pure geodesic bypass, not translation
  contraction, and connected-sum homogenization remains nonadditive. This
  theorem does not compute the still-unknown exact unknotting number of the
  counterexample or give a finite complete knot invariant.
source: codex-2026-07-24-knot-relations
depends_on: []
related:
  - THM-840-hamming-five-continuation-congruence-boundary
  - THM-853-closed-cf-return-semigroup
  - THM-1862-order-join-reduction-principle
  - THM-1975-the-path-cover-polynomial-is-the-refined-compositional-invariant
  - THM-2162-signed-endpoint-cocycle-and-bv-component-split
  - THM-2174-endpoint-phase-scale-obstruction
external:
  - "Mark Brittenham and Susan Hermiller, Unknotting number is not additive under connected sum, arXiv:2506.24088v2."
  - "Mark Brittenham and Susan Hermiller, Unknotting number and connected sums: The knots 4_1 and 5_1, arXiv:2601.18757v1."
---

# THM-2176 -- the Gordian continuation profile and interaction cocycle

The failure of additivity of unknotting number is not repaired by guessing a
second scalar. The exact operation-ready object is a continuation profile,
and its first compression is a positive interaction cocycle.

The proofs below first work in an abstract setting. This makes clear which
claims use knot theory and which follow solely from metric composition.

## 1. Metric monoids and the min-plus regular representation

Let `(M,+,0)` be a commutative monoid with a metric `d` satisfying

```text
d(x+y,x'+y') <= d(x,x')+d(y,y').                    (1)
```

Put

```text
ell(x)=d(x,0).                                       (2)
```

For every `x in M`, translation `T_x(a)=x+a` is 1-Lipschitz, the translations
commute, and

```text
T_(x+y)=T_x after T_y,
ell(x+y)<=ell(x)+ell(y).                             (3)
```

These are immediate from (1).

There is a stronger exact representation. Define the min-plus kernel

```text
P_x(a,b)=d(x+a,b)                                    (4)
```

and tropical convolution

```text
(P_x tensor P_y)(a,c)
  =inf_(b in M) [P_x(a,b)+P_y(b,c)].                 (5)
```

Then

```text
P_x tensor P_y=P_(x+y).                              (6)
```

Indeed, for every intermediate `b`,

```text
d(x+y+a,c)
 <=d(x+y+a,y+b)+d(y+b,c)
 <=d(x+a,b)+d(y+b,c).                                (7)
```

This proves that the left side of (6) is at least the right side. Taking
`b=x+a` makes the two sides equal. The representation is faithful: the row
`P_x(0,-)=d(x,-)` has its unique zero at `x`.

Evaluating (6) at `(0,0)` gives the exact common-intermediate formula

```text
ell(x+y)=inf_(a in M) [d(x,a)+ell(y+a)].             (8)
```

Restricting the infimum to a landmark dictionary gives a certified upper
bound. The ordinary separate-cost estimate is the single landmark `a=0`.
Thus "search for a common intermediate" is tropical matrix multiplication,
not an analogy.

### Knot specialization

Take `M` to be isotopy classes of oriented knots in `S^3`, `+` to be connected
sum, `0=U` the unknot, and `d=d_G` the Gordian distance. Crossing-change
sequences for two pairs may be performed in disjoint summands, so (1) holds:

```text
d_G(K#L,A#B)<=d_G(K,A)+d_G(L,B).                     (9)
```

Consequently connected sum has the faithful min-plus representation

```text
P_K tensor P_L=P_(K#L),
u(K)=P_K(U,U),                                       (10)
```

where `u(K)=d_G(K,U)`.

## 2. The universal continuation profile

The next statement does not require a metric. Let `M` be any monoid and let
`f:M->Y` be any observable. Define

```text
Pi_f(x):M->Y,            Pi_f(x)(z)=f(x+z),           (11)
x equiv_f y              iff Pi_f(x)=Pi_f(y).         (12)
```

> **Universal continuation theorem.** `equiv_f` is the largest
> connected-sum congruence contained in `ker(f)`.

### Proof

Putting `z=0` in (12) shows that `equiv_f` is contained in `ker(f)`. If
`x equiv_f y`, then for every `a,z`,

```text
f((x+a)+z)=Pi_f(x)(a+z)=Pi_f(y)(a+z)=f((y+a)+z),     (13)
```

so `x+a equiv_f y+a`; hence `equiv_f` is a congruence.

Conversely, let `~` be any congruence contained in `ker(f)`. If `x~y`, then
`x+z~y+z` for every continuation `z`, and containment in `ker(f)` gives
`f(x+z)=f(y+z)`. Thus `x equiv_f y`. This proves maximality. QED.

Equivalently, `Pi_f` is the coarsest refinement of `f` on which every future
monoid operation has a deterministic update:

```text
Pi_f(x+a)(z)=Pi_f(x)(a+z).                            (14)
```

This is the infinite Myhill--Nerode form of THM-840's operation-kernel
criterion and THM-853's finite Moore refinement.

For knots, set

```text
Pi_K(J)=u(K#J).                                      (15)
```

The profile is not claimed to be finite or a complete invariant of knot
isotopy. It is exactly the minimal operation-closed refinement of the scalar
`u`.

It is strictly finer than `u`. Let

```text
K=T(2,7),                 Kbar=mirror(K).             (16)
```

Then `u(K)=u(Kbar)=3`. The signature bound and separate unknotting give

```text
Pi_K(K)=u(K#K)=6,                                   (17)
```

whereas Brittenham--Hermiller give

```text
Pi_Kbar(K)=u(Kbar#K)<=5.                             (18)
```

Thus equality of unknotting number is not a connected-sum congruence, and the
continuation `K` is an explicit one-query splitter of the mirror pair.

## 3. The interaction defect is a positive 2-coboundary

Define

```text
sigma(x,y)=ell(x)+ell(y)-ell(x+y).                   (19)
```

Subadditivity makes `sigma>=0`. It is normalized and symmetric:

```text
sigma(x,0)=0,
sigma(x,y)=sigma(y,x).                               (20)
```

It also satisfies the exact cocycle identity

```text
sigma(x,y)+sigma(x+y,z)
 =sigma(y,z)+sigma(x,y+z).                           (21)
```

Both sides are

```text
ell(x)+ell(y)+ell(z)-ell(x+y+z).                     (22)
```

In normalized bar cochains of the commutative monoid with trivial real
coefficients, `sigma=delta ell`. It is therefore a **2-coboundary**, not a
nontrivial cohomology class. Its pointwise nonnegativity is the additional
metric content.

For knots,

```text
sigma(K,L)=u(K)+u(L)-u(K#L).                         (23)
```

Brittenham--Hermiller's "symbiont" relation is precisely the support
`sigma(K,L)>0`, but the magnitude and cocycle law contain more information
than that binary shadow.

For any finite list, put

```text
Sigma(x_1,...,x_r)
 =sum_i ell(x_i)-ell(sum_i x_i).                     (24)
```

Along any binary parenthesization tree, `Sigma` is the sum of `sigma` over
the internal merges. This follows by telescoping, and (21) is the elementary
rotation which proves independence of the tree. In particular the support of
symbiosis is not an arbitrary graph. For every triple,

```text
sigma(x,y)>0
  => [sigma(y,z)>0 or sigma(x,y+z)>0],               (25)
```

because the right side of (21) must be positive.

## 4. The rooted Gordian order and poison-cone monotonicity

Define the rooted geodesic relation

```text
a <=_0 x
  iff ell(x)=d(x,a)+ell(a).                          (26)
```

Thus `a` lies on a geodesic from `x` to the root. This is a partial order.
Reflexivity is immediate. If both `a<=_0 x` and `x<=_0 a`, adding the two
equalities forces `d(a,x)=0`, proving antisymmetry. If
`b<=_0 a<=_0 x`, then

```text
ell(x)-ell(b)=d(x,a)+d(a,b).                         (27)
```

The triangle inequality gives `d(x,b)` at most the right side, while the
reverse triangle inequality gives

```text
ell(x)-ell(b)<=d(x,b).                               (28)
```

Hence equality holds in both and `b<=_0 x`. For an integer graph metric,
`ell` grades the order.

The defect is coordinatewise monotone on this order:

```text
a<=_0 x and b<=_0 y
  => sigma(a,b)<=sigma(x,y).                         (29)
```

Indeed, (1), (26), and the triangle inequality give

```text
ell(x+y)
 <=d(x+y,a+b)+ell(a+b)
 <=ell(x)-ell(a)+ell(y)-ell(b)+ell(a+b),             (30)
```

which rearranges to (29).

Consequently:

```text
{sigma=0} is a down-set,
{sigma>0} is an up-set                                (31)
```

in the product rooted order. This is the abstract reason a single
Brittenham--Hermiller shortcut propagates to every pair having its two knots
as descendants on minimal unknotting paths.

## 5. Translation catalysis versus geodesic bypass

The defect has an exact directional split:

```text
sigma(x,y)=C_y(x)+B_y(x),                            (32)

C_y(x)=ell(x)-d(x+y,y),                              (33)
B_y(x)=d(x+y,y)+ell(y)-ell(x+y).                     (34)
```

Both terms are nonnegative. Translation nonexpansivity proves `C_y(x)>=0`;
the triangle path `x+y -> y -> 0` proves `B_y(x)>=0`.

The two mechanisms are different:

```text
C_y(x)>0: adding y contracts the x-to-root distance;
B_y(x)>0: the separated-summand intermediate y is bypassed by a shorter
          root path.                                 (35)
```

A useful calibration lemma excludes the first mechanism. Suppose
`phi:M->R` is additive and 1-Lipschitz:

```text
|phi(p)-phi(q)|<=d(p,q),                             (36)
```

and `ell(x)=|phi(x)|`. Then for every `y`,

```text
d(x+y,y)>=|phi(x+y)-phi(y)|=|phi(x)|=ell(x),         (37)
```

while translation gives the reverse inequality. Hence

```text
d(x+y,y)=ell(x),             C_y(x)=0.               (38)
```

So every interaction involving a calibrated first summand is pure bypass in
this directional split.

### The first counterexample is a folded Gordian square

For `K,Kbar` in (16), let `X=K#Kbar`. Signature is additive, changes sign
under mirroring, and changes by at most two under one crossing change.
Therefore half the signature difference lower-bounds Gordian distance.
Since the signatures of `K,Kbar,X` are `6,-6,0` up to the global sign
convention,

```text
d_G(X,Kbar)>=3.                                      (39)
```

Unknotting the `K` summand gives the reverse bound, so

```text
d_G(X,Kbar)=3=u(K),             C_Kbar(K)=0.          (40)
```

But `u(X)<=5`, and hence

```text
B_Kbar(K)=6-u(X)>=1.                                 (41)
```

The same holds after swapping `K,Kbar`. More completely,

```text
d_G(U,K)=d_G(U,Kbar)=d_G(K,X)=d_G(Kbar,X)=3,
d_G(K,Kbar)=6,
d_G(U,X)<=5.                                         (42)
```

Thus the known example shortens a diagonal while leaving every separated
side at full calibrated length. It is a **geodesic bypass**, not a catalyst
which makes one summand intrinsically cheaper. The exact value `u(X)` remains
open in the cited paper, so (41) is a lower bound, not an exact defect.

## 6. Connected-sum homogenization does not restore additivity

For every `x`, the sequence

```text
a_n=ell(n x)                                         (43)
```

is subadditive. Fekete's lemma therefore gives the connected-sum
homogenization

```text
ell_hash(x)
 =lim_(n->infinity) ell(n x)/n
 =inf_(n>=1) ell(n x)/n.                             (44)
```

It is positively homogeneous and subadditive:

```text
ell_hash(m x)=m ell_hash(x),
ell_hash(x+y)<=ell_hash(x)+ell_hash(y).               (45)
```

For knots we call (44) **connected-sum homogenized unknotting number**. We
avoid the bare phrase "stable unknotting number", which is already used for
a different twist-family construction in the literature.

Signature calibration gives

```text
ell_hash(K)=ell_hash(Kbar)=3.                         (46)
```

On the other hand, repeating the five-change certificate gives

```text
ell_hash(K#Kbar)<=u(K#Kbar)<=5.                       (47)
```

Therefore

```text
ell_hash(K)+ell_hash(Kbar)-ell_hash(K#Kbar)>=1.       (48)
```

Homogenization removes bounded one-body fluctuations but does **not** repair
connected-sum additivity.

## 7. The native binary relations are not tournaments

Because the Gordian graph is an integer geodesic metric, its threshold
relations

```text
R_m={(x,y):d_G(x,y)<=m}                              (49)
```

satisfy exact Boolean composition

```text
R_m after R_n=R_(m+n).                               (50)
```

The forward inclusion is the triangle inequality. For the reverse inclusion,
choose a vertex at the appropriate distance along a geodesic. Connected sum
respects the budget grading by (9).

The three native relation layers have different types:

```text
R_m:                symmetric graded relation;
<=_0:               rooted partial order;
{sigma>0}:          symmetric interaction relation with a weight sigma.
                                                               (51)
```

None is a tournament. Orienting knots only by `u` produces a transitive rank
shadow with ties; (17)--(18) show that a tied mirror pair can have different
continuation behavior. Forcing a total orientation would delete the exact
coordinate responsible for the counterexample.

The tournament comparison is operation-specific:

- under order-join, THM-1862 proves `H(A▷B)=H(A)H(B)`, so `log H` has zero
  interaction defect and its continuation profile collapses to its scalar;
- under cyclic substitution, scalar `H` does not compose on the verified
  finite universe of THM-1975, while its path-cover profile does there;
- the old Reidemeister-R3 tournament analogy is only a hostile warning:
  3-cycle reversal preserves `H` at `n=4` but already fails at `n=5`.

The general rule is therefore not "replace knots by a tournament." It is:

```text
retain the operation-indexed kernel or continuation profile;
take a tournament shadow only after an intrinsic antisymmetric comparator
has been proved target-preserving.                    (52)
```

## 8. Scope and next finite target

This theorem supplies an exact refinement architecture, not an algorithm for
all unknotting numbers. The full profile (15) is infinite and may identify
non-isotopic knots. A useful finite program must choose a continuation
dictionary `D` and retain

```text
Pi_K restricted to D,
sigma(K,-) restricted to D,
certified crossing-change witnesses and mirror/signature sidecars.         (53)
```

The first controls are forced:

```text
D contains U and T(2,7);
K and mirror(K) must split as in (17)--(18);
product Gordian cones must respect (29);
every proposed cocycle table must satisfy (21).       (54)
```

For LRC, THM-2162 has already realized the legitimate common-intermediate
transfer by retaining a whole safe core before integrating the next comb.
THM-2174 now supplies the matching continuation-congruence warning: endpoint
phase labels without signed magnitude preserve extensive algebraic state but
not the target measure. The universal theorem above explains exactly what a
successful replacement must do, but proves no new LRC(14) case by itself.

QED.
