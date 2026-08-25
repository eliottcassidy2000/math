---
id: THM-4073
title: "Even-graph diameter-layer exact cycle distance and commutative weighted lift"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. For every n>=3,
  tree-diameter layer 2<=D<=n-1, and 3<=m<=n, the simple orbit-graph distance
  from the empty Eulerian graph to a single m-cycle is exactly
  ceil((m-2)/(D-1)). In particular the empty-class degree is D-1 and both its
  eccentricity and the graph diameter are at least ceil((n-2)/(D-1)). The
  orbit projection has exact path lifting but is not generally a covering.
  Cycle-length convolutions lie in the Bose--Mesner algebra of a commuting
  symmetric translation association scheme before quotienting and give
  commuting multiplicity-weighted operators after quotienting; their
  off-diagonal Boolean supports are the relations R_(n,k) and need not
  commute, already at n=4.
source: codex-frontier-synthesis-creative-20260825c / even-graph metric lane
audit: >
  PASS. The primary cyclic-order/anchored-triangle implementation performs
  84,024,922 exact generator gates and checks all 55 triples (n,D,m) through
  n=7. The independent path scans 33,864 edge subsets, recognizes cycles as
  connected 2-regular supports, directly canonicalizes all Eulerian graphs
  through n=6 with 745,164 relabel gates, constructs the simple orbit graphs,
  and checks all 30 triples there using 433,754 transition gates. It also
  recovers the n=6 hostile plateau diam(E_6^(3))=diam(E_6^(4))=3. The primary
  audit computes the n=4 weighted commuting matrices and verifies that their
  Boolean off-diagonal supports do not commute. Normal and optimized outputs
  byte-match; both scripts have zero assert nodes and zero floating literals.
depends_on:
  - THM-4069-even-graph-basis-dependence-and-canonical-cycle-envelope
script: 04-computation/even_graph_diameter_layer_metric_thm4073.py
output: 05-knowledge/results/even_graph_diameter_layer_metric_thm4073.out
script_sha256: d34296b85509a82403b9352a5febab4f0f410488c515e115f24cde7df0d97100
output_sha256: a52befd643d15f22330755a70b5fb6131f81ef855f7c82adc3191d32e4ef316f
independent_audit_script: 04-computation/even_graph_diameter_layer_metric_thm4073_independent_audit.py
independent_audit_output: 05-knowledge/results/even_graph_diameter_layer_metric_thm4073_independent_audit.out
independent_audit_script_sha256: aefa9732c1d07a145aa99f3791233cc5ac544593869ac3ab98159cd5b565c32e
independent_audit_output_sha256: 726ce944a589cee0caa0e84037de9a5f9460d36d9e3db922992ac606dec901ab
hash_basis: raw LF bytes
---

# THM-4073 -- exact radial metric in every even-graph diameter layer

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** THM-4069
classifies every spanning-tree image by one integer, the tree diameter. This
theorem determines a sharp metric profile inside every such layer. It also
identifies the multiplicity sidecar that restores the commutative structure
destroyed by the simple isomorphism-class edge image.

## 1. Labelled lift and simple orbit graph

Retain THM-4069's binary cycle space

```text
Z_n={F subset E(K_n): every vertex has even degree}.              (1)
```

Let `C_(n,k)` be the set of all labelled simple `k`-cycles in `K_n`. For
`2<=D<=n-1`, put

```text
C_(n,<=D+1)=union_(k=3)^(D+1) C_(n,k),
Gamma_(n,D)=Cay(Z_n,C_(n,<=D+1)).                                (2)
```

The symmetric group acts on `(2)` by graph automorphisms. If
`q:Z_n -> Z_n/S_n` is graph isomorphism, then the simple orbit image of
`Gamma_(n,D)` is exactly

```text
E_n^(D)=union_(k=3)^(D+1) R_(n,k)                                (3)
```

from THM-4069. Thus `(3)` is the image of every spanning tree of diameter
`D`, and `D=n-1` is the canonical all-cycle envelope.

The orbit map has the following path-lifting property. Given an orbit path

```text
[F_0]--[F_1]--...--[F_r]                                        (4)
```

and a chosen representative `F'_0` of `[F_0]`, lift inductively. If the next
edge is witnessed by

```text
G -> G xor C,       G in [F_i], C in C_(n,<=D+1),                (5)
```

choose `sigma in S_n` with `sigma(G)=F'_i`. Then

```text
F'_i -> F'_i xor sigma(C)=sigma(G xor C)                         (6)
```

lifts that edge, because the generator set is `S_n`-invariant. Consequently

```text
dist_(E_n^(D))([F],[G])
 =min_{F' in [F], G' in [G]} dist_(Gamma_(n,D))(F',G').          (7)
```

Projection of a labelled path gives one inequality in `(7)`, and the lift
`(6)` of a shortest orbit path gives the other.

This is exact path lifting, not a graph covering: distinct labelled
neighbors can collapse to one orbit neighbor.

## 2. Exact distance to a single cycle

For every `3<=m<=n`, let `[C_m]` denote the isomorphism class of one
`m`-cycle and `n-m` isolated vertices. Then

```text
boxed:
dist_(E_n^(D))([empty],[C_m])=ceil((m-2)/(D-1)).                 (8)
```

### 2.1 Lower bound: a shortest XOR family is edge-connected

By `(7)`, an orbit path of length `r` gives labelled simple cycles
`Q_1,...,Q_r`, each of length at most `D+1`, whose symmetric difference is a
labelled `m`-cycle `C`. Repeatedly discard any nonempty zero-XOR subfamily.
This can only shorten the family, so it is enough to bound an
inclusion-minimal family

```text
Q_1 xor ... xor Q_t=C.                                          (9)
```

Form the intersection graph on the `Q_i`, joining two when they share an
**edge**. This graph is connected. Indeed, suppose it had components
`I_1,...,I_s`, and put

```text
H_j=xor_(i in I_j) Q_i.                                         (10)
```

Different components use disjoint edge sets. Hence the supports of the
`H_j` are disjoint and their union is the edge set of `C`; in particular
each `H_j` is an Eulerian edge-subgraph of `C`. The only Eulerian edge
subgraphs of one simple cycle are `empty` and the whole cycle. Minimality of
`(9)` excludes `H_j=empty`, while two edge-disjoint `H_j` cannot both equal
`C`. Thus `s=1`.

Order the `Q_i` along a rooted spanning tree of their edge-intersection
graph. After `Q_1`, each new `Q_i` shares an edge, hence at least two vertices,
with the preceding union. Therefore

```text
m <= |V(union_i Q_i)|
  <= |Q_1|+sum_(i=2)^t(|Q_i|-2)
   = 2+sum_(i=1)^t(|Q_i|-2)
  <= 2+t(D-1).                                                   (11)
```

This gives `t>=ceil((m-2)/(D-1))`, and hence the lower bound in `(8)`.
The zero-XOR deletion is load-bearing: without it an irrelevant dependent
component could be disjoint from the component producing `C`.

### 2.2 Upper bound: bounded polygon fans

Put

```text
t=ceil((m-2)/(D-1)).                                             (12)
```

Choose positive integers

```text
a_1+...+a_t=m-2,       1<=a_i<=D-1.                             (13)
```

Such a composition exists by `(12)`. Label the cyclically ordered vertices
of `C` as `v_0,...,v_(m-1)`, put `b_i=a_1+...+a_i` and `b_0=0`, and take the
fan cells

```text
Q_i=(v_0,v_(b_(i-1)+1),v_(b_(i-1)+2),...,v_(b_i+1),v_0).        (14)
```

Each `Q_i` has length `a_i+2<=D+1`. Consecutive fan chords occur twice and
cancel, while every boundary edge of `C` occurs once, so

```text
Q_1 xor ... xor Q_t=C.                                          (15)
```

This supplies the matching path and proves `(8)`.

## 3. Degree, eccentricity, and the diameter hostile

From the empty graph, every `k`-cycle generator reaches the single class
`[C_k]`; different lengths have different edge counts. Hence

```text
boxed: deg_(E_n^(D))([empty])=D-1.                              (16)
```

The terminal case `m=n` in `(8)` gives

```text
ecc_(E_n^(D))([empty]) >= ceil((n-2)/(D-1)),
diam(E_n^(D))              >= ceil((n-2)/(D-1)).                 (17)
```

These are lower bounds, not a diameter formula. Strict edge inclusion does
not force strict diameter decrease. The independent exact orbit scan finds

```text
diam(E_6^(3))=diam(E_6^(4))=3,                                  (18)
```

even though `E_6^(3)` is a proper spanning subgraph of `E_6^(4)`. Thus the
cycle-distance profile is a sharp radial statistic while the global farthest
class can lie elsewhere.

## 4. The commuting structure lives in multiplicities

There is an all-`n` association-scheme lift, but it does **not** say that the
simple adjacency matrices of the `R_(n,k)` commute.

Let `O` range over the `S_n`-orbits in the additive group `Z_n`, and define a
relation on labelled Eulerian graphs by

```text
(F,G) in S_O iff F xor G in O.                                  (19)
```

The relations `(19)` form a symmetric translation association scheme. They
partition `Z_n x Z_n`; intersection numbers depend only on the orbit of
`F xor G`; and their adjacency operators commute because they are convolution
by orbit sums in the abelian group `Z_n`. The orbit `O=[C_k]` gives the
labelled `k`-cycle operator

```text
(A_k f)(F)=sum_(C in C_(n,k)) f(F xor C).                        (20)
```

Restrict `(20)` to `S_n`-invariant functions. In the basis of graph-isomorphism
classes this is the multiplicity-weighted matrix

```text
M_k([F],[G])
 =#{C in C_(n,k):[F xor C]=[G]}.                                (21)
```

It follows for all `k,l` that

```text
M_k M_l=M_l M_k.                                                 (22)
```

Moreover

```text
sum_[G] M_k([F],[G])=|C_(n,k)|=n!/(2k(n-k)!),                   (23)
|[F]| M_k([F],[G])=|[G]| M_k([G],[F]).                          (24)
```

Thus the `M_k` are simultaneously diagonalizable after the reversible
orbit-size weighting. Their labelled Fourier eigenvalues are explicit. A
character is indexed by an edge mask `H` modulo the cut space
`Z_n^perp`, and

```text
lambda_k(H)=sum_(C in C_(n,k)) (-1)^|E(H) intersect E(C)|.       (25)
```

The off-diagonal positive support of `M_k` is exactly `R_(n,k)`. But deleting
diagonal orbit loops, multiplicities, and then replacing every positive entry
by one is not an algebra homomorphism. At `n=4`, in the orbit order
`[empty],[C_3],[C_4]`, the exact matrices are

```text
M_3 = [[0,4,0],        M_4 = [[0,0,3],
       [1,0,3],               [0,3,0],
       [0,4,0]],              [1,0,2]].                          (26)
```

They commute. Their simple off-diagonal supports are respectively the path
`[empty]--[C_3]--[C_4]` and the single edge `[empty]--[C_4]`, whose adjacency
matrices do not commute. This is the exact loss ledger:

```text
source:     labelled translation association scheme;
map:        quotient by graph isomorphism, then Boolean edge support;
preserved:  reachability and the orbit-distance formula (7);
destroyed:  local multiplicity, regularity, loops, and commutative products;
sidecar:    the integer matrices M_k.                            (27)
```

## 5. Exact audits and boundary

The primary audit generates simple cycles from cyclic orders, generates
`Z_n` from an anchored-triangle basis, and runs full labelled BFS for every
layer through `n=7`. Its 55 distance profiles are:

| `n` | profiles as `D:(d(C_3),...,d(C_n))` |
|---:|---|
| 3 | `2:(1)` |
| 4 | `2:(1,2); 3:(1,1)` |
| 5 | `2:(1,2,3); 3:(1,1,2); 4:(1,1,1)` |
| 6 | `2:(1,2,3,4); 3:(1,1,2,2); 4:(1,1,1,2); 5:(1,1,1,1)` |
| 7 | `2:(1,2,3,4,5); 3:(1,1,2,2,3); 4:(1,1,1,2,2); 5:(1,1,1,1,2); 6:(1,1,1,1,1)` |

The independent audit does not use a cycle-space basis or cyclic-order cycle
generation. It scans every edge subset through `n=6`, recognizes cycles by
connected `2`-regular support, canonicalizes by direct permutation, builds
the simple orbit graphs, and checks `(8)`, `(16)`, and `(17)` there. It also
prints every edge count, empty eccentricity, and graph diameter, including
the hostile plateau `(18)`.

This theorem does not determine the full diameter, chromatic number, clique
number, or unweighted spectrum of `E_n^(D)`. It gives an exact family of
distances, a sharp marked-vertex degree, and the correct weighted algebra in
which cycle lengths commute.
