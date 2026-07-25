---
id: THM-2183
title: "Order-join is an exact product for tournament reversal distance"
status: >
  PROVED + VERIFIED-EXACT. For equal corresponding factor orders,
  unlabeled arc-reversal distance satisfies
  d_iso(A join X,B join Y)=d_iso(A,B)+d_iso(X,Y).
  A strict two-edge uncrossing eliminates every block-mixing bijection.
  Thus common order-join factors cancel metrically, fixed block-size joins
  embed Cartesian products isometrically, and no SCC or marked-cut
  hypothesis is needed. The full unlabeled tournament flip graph is still
  not a partial cube: its four-vertex quotient already contains a triangle.
source: codex-2026-07-24-knot-relations
depends_on: []
related:
  - THM-1862
  - THM-1936
  - THM-1975
  - THM-2176
script: 04-computation/tournament_order_join_metric_product_thm2183.py
output: 05-knowledge/results/tournament_order_join_metric_product_thm2183.out
script_sha256: 80424d559ac9a1f473ac9c239f8f70795220cd51c7342617bb336f8e8a933f70
output_sha256: dfdcfcedb27fb531e73ab281d8cf80acd3238d2ca89127bf78514bdefa4175e8
hash_basis: working-tree bytes (LF)
---

# THM-2183 -- order-join is an exact tournament metric product

For tournaments `T,S` of the same order, define unlabeled arc-reversal
distance by

```text
d_iso(T,S)
 =min_(bijections pi:V(T)->V(S))
    #{unordered {u,v}:
       (u->_T v) xor (pi(u)->_S pi(v))}.             (1)
```

For tournaments on disjoint vertex sets, `A join X` denotes the order-join
in which every vertex of `A` beats every vertex of `X`.

## 1. Exact product theorem

> **Theorem.** Let `A,B` have `a` vertices and let `X,Y` have `b` vertices.
> Then
>
> ```text
> d_iso(A join X,B join Y)
>   =d_iso(A,B)+d_iso(X,Y).                          (2)
> ```

**Upper bound.** Take optimal bijections `A->B` and `X->Y`. Their union is
a bijection of the joins. It pays the two internal distances, and every
cross arc agrees. Thus the left side of (2) is at most the right side.

**Block-mixing uncrossing.** Put

```text
P=A join X,                  Q=B join Y.             (3)
```

Fix any bijection `pi:V(P)->V(Q)`. If it mixes the two blocks, equality of
corresponding block sizes implies that there are vertices

```text
a_0 in A,       x_0 in X,
d=pi(a_0) in Y, e=pi(x_0) in B.                     (4)
```

Let `pi'` swap only the two images `d,e`. The pair `{a_0,x_0}` is a mismatch
before the swap: `a_0->_P x_0`, whereas `e->_Q d`. It agrees afterward, so
this pair saves one reversal.

For every third vertex `z`, abbreviate

```text
p=1[a_0->_P z],       r=1[x_0->_P z],
q=1[d->_Q pi(z)],     s=1[e->_Q pi(z)].             (5)
```

The cost of the two incident pairs before and after the swap is

```text
C =|p-q|+|r-s|,
C'=|p-s|+|r-q|.                                     (6)
```

There are four membership cases:

```text
z in A, pi(z) in B:   (r,q)=(0,0),
  C=p+s,                 C'=|p-s|<=C;

z in X, pi(z) in Y:   (p,s)=(1,1),
  C=2-q-r,               C'=|r-q|<=C;

z in A, pi(z) in Y:   (r,s)=(0,1),
  C'-C=q-p-|p-q|<=0;

z in X, pi(z) in B:   (p,q)=(1,0),
  C'-C=r-s-|r-s|<=0.                               (7)
```

Thus `pi'` costs at least one fewer reversal than `pi`. It also removes one
crossed pair of vertices from (4). Repeating the swap eliminates every
block crossing while strictly decreasing the cost at each step.

It follows that no optimal bijection mixes the blocks. Every block-respecting
bijection has cost

```text
cost(pi|_A:A->B)+cost(pi|_X:X->Y)
  >=d_iso(A,B)+d_iso(X,Y).                           (8)
```

Together with the upper bound, this proves (2). QED.

## 2. Cancellation and Cartesian-product patches

Taking `A=B` or `X=Y` in (2) gives both cancellation laws:

```text
d_iso(A join X,A join Y)=d_iso(X,Y),
d_iso(X join A,Y join A)=d_iso(X,Y).                 (9)
```

Consequently, on isomorphism classes, the fixed-size map

```text
([A],[X]) |-> [A join X]                            (10)
```

is an isometric embedding of the Cartesian product of the two quotient flip
graphs. Iterating (2) gives an isometric product patch for every fixed ordered
block-size composition. In particular, if two tournaments have ordered
strong-component decompositions with the same component-size sequence, their
distance is the sum of the distances between corresponding components.

This strengthens THM-1862's invariant algebra:

```text
c_3: additive;
H:   multiplicative;
d_iso: l1-additive.                                  (11)
```

The proof uses no strong connectivity, uniqueness of a displayed cut, root,
or marked-block data. It also covers a wholesale block-swap candidate when
`a=b`: the constant one-way cut makes that mixing strictly uncrossable.

The matching factor orders in the statement are load-bearing. No claim is
made when `|A|!=|B|` or `|X|!=|Y|`. Empty factors are the trivial boundary
case.

## 3. Exact negative control for knot interaction

The binary exchange in (7) is a Monge law: every attempt to exploit the
unlabeled quotient by mixing the two constant-cut blocks can be locally
improved. This is precisely what is absent from connected-sum crossing
changes. If one retains a connected-sum sphere, the decomposition-respecting
Gordian cost is the sum of summand costs; after forgetting the sphere, mixed
crossing changes may create the positive interaction defect of THM-2176.

For tournament order-join, forgetting the displayed cut has **zero** such
defect by (2). For knots, the Brittenham--Hermiller examples prove that the
analogous quotient defect can be positive. The comparison identifies the
missing structure--a constant-cut uncrossing law--without forcing knots into
a tournament relation.

## 4. The ambient quotient graph is not a partial cube

The isometric product patches do not make the full unlabeled tournament flip
graph a partial cube. Encode a labelled `n`-vertex tournament by bits on

```text
(0,1),(0,2),...,(n-2,n-1),                           (12)
```

where zero orients the smaller label toward the larger. On four vertices the
isomorphism classes represented by masks

```text
0, 2, 4                                               (13)
```

are pairwise at `d_iso` distance one. They form a triangle. Since every
partial cube is bipartite, the quotient flip graph is not a partial cube.
Exact orbit enumeration gives four vertices and five edges for this
four-vertex quotient.

This is the useful boundary:

```text
global relabelling quotient:      has bypass triangles;
each fixed order-join patch:      exact isometric product.      (14)
```

The companion audit checks (2) for all labelled factor quadruples with
`1<=a,b<=3` and independently reconstructs the four-vertex quotient and
triangle (13).

QED.
