---
id: THM-2761
title: "Graph edge-sum discriminant codegree factorization and graceful sign gauge"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For any finite
  simple graph H, the discriminant of the polynomial
  whose roots are edge sums x_u+x_v factors exactly into vertex-difference
  powers 2*codeg_H(u,v) and squared additive-energy differences over disjoint
  edge pairs.  For a bipartite graph, the vertex sign gauge x_A=ell and
  x_B=-ell turns the roots into oriented edge differences.  Their
  discriminant detects oriented collisions; a second mirror-plus factor gives
  the discriminant of their squares and detects collisions up to sign.  Hence
  an injection into {0,...,m} is graceful iff that squared-difference
  discriminant is nonzero.  The path labels (0,1,2) are the sharp hostile:
  oriented differences (-1,1) have discriminant 4 but equal absolute value.
  This is a universal graph identity/criterion, not a proof of the Graceful
  Tree Conjecture or a Keller/LRC result.
source: a4-resolvent-next-gate-scout-2026-07-28
audit: >
  root/2026-07-28 (independent edge-pair partition, codegree exponents,
  K4/tree specializations, bipartite gauge, mirror-factor necessity,
  graceful iff scope, and normal/-O/LF-hash replay: ACCEPT)
depends_on: []
related:
  - THM-2595-modular-v4-affine-lift-dichotomy-and-six-vertex-tournament-no-go
  - THM-2606-affine-v4-parity-channels-partial-cubes-and-feuerbach-origin
  - THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration
  - THM-2756-opposite-edge-projectors-parity-cancellation-and-integral-clutch
  - THM-2758-quartic-pair-sum-sextic-resolvent-pullback-and-discriminant-square
script: 04-computation/graph_edge_sum_discriminant_graceful_thm2761.py
output: 05-knowledge/results/graph_edge_sum_discriminant_graceful_thm2761.out
script_sha256: 111af7809545c44f61648833a321d31911c898b615491aa7cd26d2cdcb86e1b4
output_sha256: 1233863e17342e7af3cb743791f948aa4048a6077d3ac987ac8906b0a327b659
hash_basis: LF-normalized bytes
---

# THM-2761 -- edge-sum discriminants split into codegree and matching factors

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2758 factors the discriminant of the six pair sums of four quantities.
The complete graph is not essential.  For every simple graph, two edges are
either adjacent or disjoint.  An adjacent edge-sum difference cancels the
common endpoint and becomes a vertex difference; a disjoint edge pair retains
a genuine four-vertex additive-energy difference.  The number of times one
vertex difference occurs is exactly a graph codegree.

This gives a useful graceful-graph reframe.  A bipartite sign gauge converts
edge sums into oriented label differences.  But oriented distinctness is not
absolute distinctness: a second, sharply necessary mirror factor detects
`d_e=-d_f`.  The result reduces graceful labeling to two explicit nonvanishing
products after the vertex injection is chosen; it does not construct that
injection.

## 1. Universal edge-sum discriminant factorization

Let `H=(V,E)` be a finite simple graph over a field `K`, assign `x_v in K` to
each vertex, and put

```text
s_e=x_u+x_v                         for e={u,v},
P_H(T)=product_(e in E)(T-s_e).                              (1)
```

For distinct vertices define the codegree

```text
c_H(u,v)=|N_H(u) intersect N_H(v)|.                         (2)
```

Also define the disjoint-edge energy factor

```text
D_H(x)=product_( {e,f} subset E, e intersect f=empty )
       (s_e-s_f)^2.                                         (3)
```

Then

```text
disc(P_H)
 =product_( {u,v} subset V)(x_u-x_v)^(2c_H(u,v)) D_H(x).   (4)
```

To prove `(4)`, partition the unordered pairs of distinct graph edges.  If

```text
e={w,u},                    f={w,v},
```

then

```text
s_e-s_f=x_u-x_v.                                           (5)
```

For fixed `{u,v}`, the possible common endpoints `w` are exactly
`N_H(u) intersect N_H(v)`, so `(5)` occurs `c_H(u,v)` times.  Squaring the
Vandermonde contributes exponent `2c_H(u,v)`.  Every remaining edge pair is
disjoint and contributes exactly one factor to `(3)`.  These cases exhaust
the discriminant product.

Consequently `P_H` is separable iff neither of the following occurs:

```text
x_u=x_v for a pair with c_H(u,v)>0;
s_e=s_f for a disjoint edge pair.                           (6)
```

If all vertex values are distinct, only the second, four-vertex additive-
energy obstruction remains.

## 2. Complete graphs, K4, and trees

For `H=K_n`, every vertex pair has codegree `n-2`, so `(4)` becomes

```text
disc(P_Kn)=disc(product_v(T-x_v))^(n-2)
            product_(disjoint e<f)(s_e-s_f)^2.             (7)
```

At `n=4`, there are exactly three disjoint edge pairs.  Their three differences
are the opposite-matching factors `d_1,d_2,d_3` of THM-2758.  Hence `(7)` is

```text
disc(P_K4)=disc(f)^2(d_1d_2d_3)^2,                         (8)
```

recovering that theorem's `disc(f)^2T^2` factorization.

For a tree, `c_H(u,v)=1` exactly when `u` and `v` are at distance two, and is
zero otherwise.  Thus its adjacent contribution is the product of squared
vertex differences over distance-two pairs; the remaining obstruction is
entirely the disjoint-edge energy product.  This is an exact decomposition of
the collision problem, not an existence proof for a labeling avoiding it.

## 3. Bipartite sign gauge

Suppose `H` is bipartite with fixed bipartition `V=A disjoint_union B`.  For an
integer vertex labeling `ell`, set

```text
x_a= ell(a)                 for a in A,
x_b=-ell(b)                 for b in B.                    (9)
```

Every edge has a unique orientation `a -> b`, and its edge sum is the oriented
label difference

```text
s_ab=ell(a)-ell(b)=d_ab.                                  (10)
```

Changing the global bipartition gauge negates every `d_e`, leaving all
discriminants below unchanged.  Thus `(4)` becomes an exact factorization of
the oriented-difference discriminant.  If `ell` is injective, every codegree
factor is nonzero: vertices with a common neighbor lie in the same bipartition
class, where `(9)` preserves label inequality.

The bipartite hypothesis is exact.  A vertex sign assignment that converts
every edge sum into a difference is precisely a two-coloring; an odd cycle has
no global gauge of the form `(9)`.

## 4. The mirror factor and an exact graceful criterion

Write `m=|E|`, let

```text
Q_H(Y)=product_(e in E)(Y-d_e^2),
M_H(ell)=product_(e<f)(d_e+d_f)^2.                         (11)
```

Since

```text
d_e^2-d_f^2=(d_e-d_f)(d_e+d_f),
```

one has the second exact factorization

```text
disc(Q_H)=disc(P_H) M_H(ell).                             (12)
```

The two factors have distinct jobs:

```text
disc(P_H)!=0       iff the oriented differences are distinct;
M_H(ell)!=0        iff no two oriented differences are negatives;
disc(Q_H)!=0       iff the absolute differences are distinct. (13)
```

Now assume

```text
ell:V -> {0,1,...,m}
```

is injective.  Every edge has nonzero absolute difference in `{1,...,m}`.
There are `m` edges, so `m` distinct absolute differences must be exactly
`{1,...,m}`.  Therefore

```text
ell is graceful  iff  disc(Q_H)!=0
                 iff  disc(P_H)M_H(ell)!=0.               (14)
```

This criterion applies in particular to every tree, because every tree is
bipartite.  The Graceful Tree Conjecture becomes the assertion that some
injection makes both explicit factors in `(14)` nonzero; `(14)` does not prove
such an injection exists.

## 5. Sharp absolute-value hostile

The mirror factor cannot be dropped.  Take the path `0-1-2`, with bipartition

```text
A={0,2},                       B={1},
```

and injective labels

```text
(ell(0),ell(1),ell(2))=(0,1,2).                            (15)
```

The oriented differences are

```text
d_01=-1,                         d_21=1.                  (16)
```

Thus `disc(P_H)=(-1-1)^2=4` is nonzero, but `M_H=0` and
`disc(Q_H)=0`: the two absolute differences are both one.  Oriented
distinctness alone is strictly weaker than graceful distinctness.

The nearby positive control `(0,2,1)` has oriented differences `(-2,-1)` and
absolute differences `(2,1)`, so it is a graceful labeling of the same path.

## 6. Exact verification

Run

```bash
python 04-computation/graph_edge_sum_discriminant_graceful_thm2761.py
python -O 04-computation/graph_edge_sum_discriminant_graceful_thm2761.py
```

Both executions byte-match the stored `13`-line transcript
`05-knowledge/results/graph_edge_sum_discriminant_graceful_thm2761.out`.  The
companion uses explicit exceptions and no truth-bearing Python assertions.  It
enumerates all `33,866` labelled simple graphs on `2` through `6` vertices and
classifies all `871,926` edge-pair occurrences as `499,398` adjacent and
`372,528` disjoint, checking the codegree inventory and `(4)` on two exact
integer patterns per graph.  It verifies the `K4` specialization.  Finally it
checks `(12)--(14)` on all `764` injections for the `23` connected labelled
bipartite graphs through four vertices, including `158` graceful labelings and
the two path controls above.

## 7. Boundary ledger

```text
PROVED HERE:              universal graph edge-sum discriminant factorization;
                          codegree multiplicities for adjacent edge pairs;
                          disjoint-edge additive-energy residual;
                          K_n, K4, and tree specializations;
                          exact bipartite edge-sum/difference gauge;
                          oriented/absolute mirror factorization;
                          graceful iff for a fixed bounded injection;
                          sharp P3 absolute-value hostile.

NOT PROVED:               existence of a graceful injection for any new graph;
                          the Graceful Tree Conjecture;
                          a nonbipartite global sign gauge;
                          equivalence with a tournament orientation;
                          a Keller, LRC, partial-cube, or N-body consequence;
                          JC(2), DC(2), or LRC(14).                       (17)
```

QED.
