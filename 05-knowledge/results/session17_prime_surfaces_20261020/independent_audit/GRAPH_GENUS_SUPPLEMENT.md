# Graph genus: equality over prime fields and a strict extension-field bound

20 October 2026. **PROVED**, as a separate consequence of Theorem FF-SURFACE
in `../FINITE_FIELD_SURFACES.md`. This supplement does not alter the approved
six-file surface packet or its earlier independent audit. No novelty claim,
new tournament interpretation, or knot interpretation is made.

## 1. Exact object and statement

Let q=p^r be an odd prime power. Let G_q be the simple graph on the nonzero
vectors of F_q^2 modulo multiplication by -1, with adjacency determined by
determinant +1 or -1. Write gamma(G_q) for its **minimum orientable graph
genus**: the least genus of a connected closed orientable surface in which
the graph embeds without crossings. This is a graph invariant, distinct
from the topology of the particular all-triangle complex K_q.

The audited surface theorem establishes that G_q is connected and q-regular,
that every edge belongs to exactly two graph triangles, and that

```
V=(q^2-1)/2,    E=q(q^2-1)/4.
```

At every vertex, the graph induced by its neighbors is a disjoint union of
q/p copies of C_p. It is the link graph of the all-triangle complex K_q.
When q=p, that complex is a closed orientable surface of genus

```
g_p=1+(p^2-1)(p-6)/24.
```

### Corollary FF-GRAPH-GENUS

For every odd prime p,

```
gamma(G_p)=1+(p^2-1)(p-6)/24.
```

For every proper finite-field extension q=p^r with r>1,

```
gamma(G_q) >= 2+(q^2-1)(q-6)/24.
```

In particular gamma(G_9)>=12. The second conclusion is only a lower bound;
it does not assert the existence of an embedding attaining that bound.

## 2. Cellular reduction and the Euler bound

It suffices to consider cellular embeddings, meaning embeddings for which
every complementary face is an open disk. Indeed, take a closed regular
neighborhood N of a graph embedded in a closed orientable surface. Since
the graph is connected, N is connected. Cap each boundary component of N
with a disk. The resulting closed orientable surface contains a cellular
embedding of the same graph and has genus at most that of the original
surface. For completeness, if N has genus h and b boundary circles and its
complement has c components with total genus k, the original surface has
genus h+k+b-c>=h; every complementary component has a boundary, so b>=c.
Thus a minimum-genus embedding can be chosen cellular.

Fix an orientation of a cellular embedding. At each vertex v, the incident
neighbors have one cyclic order, represented by a permutation rho_v that
is a single cycle of length deg(v). Facial boundary walks are obtained by
the corresponding successor rule on directed edges. One consistent
convention is that the directed edge u->v is followed by
v->rho_v(u). Every directed edge occurs in exactly one facial walk.

Every facial walk here has length at least three. A length-one walk would
require a loop. In a simple graph, a length-two walk must immediately return
along the same edge; the successor rule would then fix that neighbor at a
vertex, which is impossible for a cyclic order of length q>=3. This argument
allows repeated vertices and repeated undirected edges in longer walks;
it does not assume that all facial boundaries are simple cycles.

If F_emb is the number of embedding faces, counting directed edges gives
3 F_emb<=2E. Euler's formula therefore gives

```
2-2g = V-E+F_emb <= V-E/3,
g >= (E-3V+6)/6 = 1+(q^2-1)(q-6)/24.
```

Equality holds only if every facial walk has length three. Because the graph
is simple, each such walk is a graph triangle with three distinct vertices.
Conversely, a cellular embedding all of whose faces are triangular attains
the displayed Euler bound.

When q=p, the audited K_p supplies just such an orientable triangular
embedding. Its genus is g_p, so the lower and upper bounds agree. This proves
the prime-field assertion, including p=3 and p=5.

## 3. Why no extension-field triangular embedding exists

Suppose all faces of a cellular embedding of G_q were triangular. For each
vertex v and each neighbor u, the corner represented by the transition

```
u -> v -> rho_v(u)
```

would lie in a triangular face. Hence u and rho_v(u) would be adjacent in
the graph induced by the neighbors of v. Since rho_v is a single q-cycle,
its successive neighbor pairs would give a Hamilton cycle in that induced
graph. But this graph is q/p disjoint cycles C_p. When q>p it is disconnected,
and cannot contain a Hamilton cycle. This is a contradiction.

This local proof does not assume in advance that an embedding must use each
graph triangle at most once or must realize K_q. In particular it does not
silently identify an arbitrary cellular embedding with the fixed complex.

### The repeated-triangle loophole, explicitly excluded

Could two different faces use the same graph triangle? They cannot traverse
it in the same direction, since then they would share directed edges, each
of which belongs to just one facial boundary walk. If they traverse it in
opposite directions, then at any one of its vertices v, with the other two
vertices denoted u and w, the two facial transitions force

```
rho_v(u)=w and rho_v(w)=u.
```

Thus {u,w} would be a two-cycle of rho_v, contradicting that rho_v is a single
cycle on q>2 neighbors. This rules out the proposed use of one triangle on
both sides. The degree condition matters: the three-cycle graph by itself
has degree two and its spherical embedding has two faces on the same
triangle. That hostile example does not satisfy the degree hypothesis here.

A longer facial boundary that revisits a triangle is not a length-three
face and already makes the Euler inequality strict. No simplicity claim
about such longer facial walks is needed or asserted.

Alternatively, with this loophole excluded, a triangular embedding would
use all graph triangles exactly once. Indeed, there are exactly 2E/3 graph
triangles because each edge is in two, and exactly 2E/3 embedding faces.
The resulting face incidence would be precisely K_q, whose vertex links
are disconnected when q>p. It therefore could not be a surface embedding.
The Hamilton-cycle argument above proves the same obstruction directly.

## 4. Integer strictness and literal controls

The Euler lower-bound expression

```
g_formal(q)=1+(q^2-1)(q-6)/24
```

is an integer for every odd integer q. The factor q^2-1 is divisible by eight.
For divisibility by three, if 3 does not divide q then 3 divides q^2-1; if
3 divides q then it divides q-6. These facts give divisibility by 24.

For q>p, Section 3 rules out equality in the Euler bound. Since both graph
genus and g_formal(q) are integers, strict inequality implies

```
gamma(G_q)>=g_formal(q)+1.
```

The following are exact substitutions, not graph-genus enumeration:

| field order | V | E | g_formal | proved graph-genus conclusion |
|---:|---:|---:|---:|---|
| 3 | 4 | 6 | 0 | gamma(G_3)=0 |
| 5 | 12 | 30 | 0 | gamma(G_5)=0 |
| 7 | 24 | 84 | 3 | gamma(G_7)=3 |
| 9 | 40 | 180 | 11 | gamma(G_9)>=12 |
| 25 | 312 | 3900 | 495 | gamma(G_25)>=496 |

In particular, the formal Euler value 11 at q=9 is neither the genus of
K_9, which is not a surface, nor the minimum orientable genus of G_9.
The exact extension-field graph genera are not determined by this argument.

## 5. Dependency and scope record

The only family-specific dependency is Theorem FF-SURFACE, particularly
its degree/count calculation, neighbor-link calculation, two-triangles-per-
edge assertion, and prime-field orientable surface construction. The exact
approved dependency is `FINITE_FIELD_SURFACES.md`, SHA256
`ff26f4df81c7fed9b6025a5db84e5c601fba8af064fe57e61f4bdb1405510a8b`.
It is covered by producer manifest SHA256
`bf5919e41985f34c815c7358c59c42fb1337048ee1200813ac84c0cf0c41ba8c`.

The remaining ingredients are the regular-neighborhood reduction, rotation
description of an orientable cellular embedding, and Euler's formula,
with the relevant uses and equality cases proved above. No spectral,
modular-curve, normalization-homology, or exact genus classification theorem
is imported. No extension-field genus upper bound or genus-attainment claim
is made. All previously approved producer and audit files are unchanged.
