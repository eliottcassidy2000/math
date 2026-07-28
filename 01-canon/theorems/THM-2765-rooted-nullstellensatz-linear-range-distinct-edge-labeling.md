---
id: THM-2765
title: "Rooted Nullstellensatz linear-range distinct-edge labeling"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  Every forest with n vertices and m>0 edges has an injective integer
  vertex labeling in {0,...,n+2m-3} whose absolute edge differences are
  pairwise distinct.  In particular every n-vertex tree has such a labeling
  in {0,...,3n-5}, a range strictly below three times the graceful range.
  The proof roots every component, orders parents before children, and uses
  the coefficient-one top monomial of the vertex Vandermonde times the
  squared edge-difference discriminant.  Its exponent ladder is exactly one
  injectivity channel plus two mirror-sign channels.  This is a
  range-relaxed theorem, not the Graceful Tree Conjecture.
source: a4-resolvent-next-gate-scout-2026-07-28
depends_on:
  - THM-2270-simultaneous-balanced-cut-relation-and-six-uniform-orientation
  - THM-2761-graph-edge-sum-discriminant-codegree-factorization-and-graceful-sign-gauge
related:
  - THM-2280-centered-polynomial-grid-avoidance-and-bounded-generic-keller-fibre
  - THM-2597-six-vertex-bicycle-modular-abelianization-cycle
  - THM-2606-affine-v4-parity-channels-partial-cubes-and-feuerbach-origin
script: 04-computation/rooted_forest_nullstellensatz_range_thm2765.py
output: 05-knowledge/results/rooted_forest_nullstellensatz_range_thm2765.out
script_sha256: 8fec623eb88546f808de1c96abc72ee92ba56edbaec44a5f8d98d17f1b0192e7
output_sha256: b5ae340356e0b07f226164ce2a5ed8e22d2ae327e0a3bae82a012b85361bdc74
hash_basis: LF-normalized bytes
---

# THM-2765 -- rooted leading monomials give a linear graceful range

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

THM-2761 turns graceful labeling into nonvanishing of an explicit product,
but it does not construct a point outside that hyperplane arrangement.  A
forest has one extra coordinate that a general graph lacks: after choosing
roots, its vertices can be introduced one at a time with exactly one new
edge.  Under that order, the top monomial of the graceful-obstruction product
is triangular.  The coefficient-grid lemma from THM-2270 then evaluates it
nontrivially on a linear-size integer box.

The factor `3` has a literal source.  Each old vertex contributes one
injectivity factor, while each old edge contributes the two sign mirrors
needed to distinguish absolute edge differences.  This is an honest
`1+2` factor count, not a claimed action of `C_3` or `PSL_2(Z)`.

## 1. Rooted forest polynomial

Let `F` be a finite simple forest with

```text
n=|V(F)|,             c=number of components,       m=|E(F)|=n-c.  (1)
```

Choose one root in each component.  Order the vertices

```text
v_1,...,v_c,v_(c+1),...,v_n                              (2)
```

so that every root comes first and every nonroot vertex comes after its
parent.  Write `p(j)<j` for the parent index of `v_j`.  The `k`th nonroot
vertex is `v_(c+k)`, and its oriented edge difference is

```text
Y_k=X_(c+k)-X_(p(c+k)),                 1<=k<=m.       (3)
```

Over `Q[X_1,...,X_n]` form

```text
Phi_F(X)
 =product_(1<=i<j<=n)(X_j-X_i)
  product_(1<=h<k<=m)(Y_k^2-Y_h^2).                    (4)
```

At an integer point `ell`, the first product is nonzero exactly when the
vertex labels are injective.  The second is nonzero exactly when the
absolute edge differences are pairwise distinct.  Injection also makes every
edge difference nonzero.  Thus

```text
Phi_F(ell)!=0                                             (5)
```

is precisely the desired range-relaxed graceful certificate.

## 2. The coefficient-one leading monomial

Use lexicographic order with

```text
X_n > X_(n-1) > ... > X_1.                              (6)
```

The first product in `(4)` has leading monomial

```text
product_(j=1)^n X_j^(j-1).                              (7)
```

For `h<k`, every variable in `Y_h` has index below `c+k`, and the parent of
`v_(c+k)` also has smaller index.  Hence

```text
LM(Y_k^2-Y_h^2)=X_(c+k)^2                               (8)
```

with coefficient one.  Multiplicativity of leading monomials gives

```text
LM(Phi_F)
 =product_(j=1)^c X_j^(j-1)
  product_(k=1)^m X_(c+k)^(c+3k-3),                    (9)
```

again with coefficient exactly one.  Its exponent sum is

```text
binom(n,2)+2binom(m,2)=deg(Phi_F),                      (10)
```

so `(9)` is a full-total-degree monomial of the kind required by the
coefficient form of the Combinatorial Nullstellensatz.

For a tree, `c=1`, `m=n-1`, and `(9)` simplifies to

```text
LM(Phi_F)=product_(j=2)^n X_j^(3j-5).                  (11)
```

The exponent increments are exactly three after the first edge: one new
vertex-injectivity load plus the two mirror factors against each old edge.

## 3. Coefficient-grid extraction

Let `t_j` be the exponent of `X_j` in `(9)`, and take

```text
S_j={0,1,...,t_j}.                                      (12)
```

The successive one-variable Lagrange interpolation identity used in
THM-2270 expresses the coefficient of `product_j X_j^(t_j)` as a weighted
sum of the values of `Phi_F` on `product_j S_j`.  Because that coefficient is
one, `Phi_F` cannot vanish on the entire grid.  Therefore some grid point
`ell` satisfies `(5)`.

If `m>0`, the largest exponent in `(9)` is the one at the last nonroot:

```text
max_j t_j=c+3m-3=n+2m-3.                               (13)
```

Every coordinate of the extracted point is nonnegative and at most `(13)`.
We have proved:

> **Rooted-forest range theorem.** Every `n`-vertex forest with `m>0` edges
> admits an injective labeling
>
> ```text
> ell:V(F)->{0,1,...,n+2m-3}                            (14)
> ```
>
> whose absolute edge differences are pairwise distinct.

For an edgeless forest, `0,1,...,n-1` is immediate and also follows from the
Vandermonde part alone.  For a tree, `(14)` becomes

```text
ell:V(F)->{0,1,...,3n-5}.                              (15)
```

Since a tree has `m=n-1` edges, its graceful target range is `{0,...,n-1}`;
the range in `(15)` is `3m-2<3m`.  Thus every tree has a strict
factor-three range-relaxed graceful labeling.

## 4. Why both channels are necessary

The two products in `(4)` cannot be silently merged or dropped.

First take the three-vertex path rooted at its middle vertex and label the
root and leaves by

```text
(1,0,2).                                                (16)
```

The two oriented differences are `(-1,1)`, so they are distinct, but their
absolute values coincide.  This is the same mirror-sign hostile isolated by
THM-2761: replacing `Y_k^2-Y_h^2` by only `Y_k-Y_h` is false.

Conversely, on the rooted path of length three, labels

```text
(0,1,3,0)                                               (17)
```

have absolute edge differences `(1,2,3)` but repeat the endpoint label zero.
Thus the edge-difference discriminant alone does not imply vertex
injectivity; the Vandermonde is load-bearing.

The forest hypothesis is also structural.  It makes every nonroot introduce
exactly one new edge, so the factors in `(8)` are triangular.  A general graph
can introduce several backward edges at one vertex; formula `(9)` then does
not describe its leading monomial.

## 5. Exact verification

Run

```bash
python 04-computation/rooted_forest_nullstellensatz_range_thm2765.py
python -O 04-computation/rooted_forest_nullstellensatz_range_thm2765.py
```

Both executions byte-match the stored `12`-line transcript.  The companion
uses explicit exceptions and no truth-bearing Python assertions.  It expands
the target coefficient exactly for all `153` recursive-tree parent arrays on
two through six vertices and obtains coefficient one in every case.  It then
finds a valid point in the stated coefficient grids for all `5,913` recursive
parent arrays on two through eight vertices, checking injectivity, distinct
absolute edge differences, and the range.  It separately checks all `35`
forest exponent profiles through eight vertices and both sharp channel
hostiles.

## 6. Boundary ledger

```text
PROVED HERE (candidate):  coefficient-one rooted-forest top monomial;
                          exact 1+2 exponent decomposition;
                          every forest with m>0 has range n+2m-3;
                          every n-vertex tree has range 3n-5;
                          oriented/absolute and edge/injection hostiles.

NOT PROVED:               optimality of n+2m-3 or 3n-5;
                          a labeling in the graceful range 0,...,n-1;
                          the Graceful Tree Conjecture;
                          the same monomial formula for cyclic graphs;
                          a C3, modular-group, partial-cube, Keller, or LRC
                          action or consequence.                         (18)
```

QED (candidate).
