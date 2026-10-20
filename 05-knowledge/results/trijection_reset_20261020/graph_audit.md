# Independent audit: the marked-path trijection needs its order

**PROVED:** the explicit coordinate maps, their inverse, symmetry laws, predicate
transfers and quotient failures below. **FINITE-EXACT:** the eight order-four
rows and twenty three-edge graph objects, independently checked by
[`independent_check.py`](independent_check.py), with
[`independent_check.json`](independent_check.json) as its receipt. These are
elementary deductions and a repackaging of the existing tiling model; no
literature novelty claim is made.

## Inheritance and scope

The closest proved mechanism is the path-normalized tiling definition in
[`01-canon/definitions.md`](../../../01-canon/definitions.md), as used in
[THM-4729, setting and L2](../../../01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md).
The quotient fibre formula is already
[LEM-003, Aut acts freely on Hamiltonian paths](../../../01-canon/theorems/LEM-003-aut-acts-freely-on-ham-paths.md).
The relevant corrected near miss is the 2026-10-18 entry in
[`MISTAKES.md`](../../../01-canon/MISTAKES.md): changing vertex labels acts on
tournament reversal bits by `x -> P_g x xor epsilon_g`, not just a permutation.
The hostile below first occurs at four source vertices. The essential sidecar
is the order induced by the marked directed path. Nothing here uses the
provisional pure-blue classification in THM-4730 or THM-4729 L4.

The direct predecessor for the three-edge graph question is
[`defect_boundary_20261020/web.md`, sections 4-6](../defect_boundary_20261020/web.md):
it already proves the `12,4,4` classification, separates complementation
from tetrahedral opposition, and explains why unsigned line graphs identify
the claw and triangle while signed clique products separate them. The
[`tournament_future_boundary_20261020` inheritance ledger](../tournament_future_boundary_20261020/inheritance.md)
also records that equal incidence rank three loses the location/support of
the kernel. The classification below is reproduced as an explicit control,
not claimed as a new result. The retained path order is a different sidecar
from incidence rank or signed line-graph data.

## Exact three-way map

For `n >= 1`, let `M_n` be the isomorphism classes of pairs `(T,P)` with `T` a
tournament and `P=(v_1,...,v_n)` a specified directed Hamiltonian path. Pair
isomorphisms preserve the path. Label `v_i` by its path position `i`, so every
adjacent arc is `i -> i+1`. For each nonadjacent pair `i<j`, put

```text
x_ij = 1 iff j -> i,
G has edge {i,j-1} iff x_ij=1,
U has arc (j-1) -> i iff x_ij=1; otherwise i -> (j-1).
```

Thus `G` is a simple graph and `U` a tournament, each on the *ordered* fixed
set `[n-1]`. For `n=1`, both targets have zero vertices. The free slots
`{(i,j): 1<=i<j-1<=n-1}` correspond bijectively to the unordered pairs of
`[n-1]` by `(i,j) -> {i,j-1}`. Reading every graph edge or tournament arc
therefore restores all the bits, and adjoining the forced forward path
restores `(T,P)`. These are genuine bijections

```text
M_n  <->  simple graphs on fixed ordered [n-1]
     <->  tournaments on fixed ordered [n-1],
|M_n| = 2^binom(n-1,2).
```

The third object is a full tournament on **one fewer vertex**, not the same
unmarked tournament `T`. In particular it need not retain the transitivity
or strong connectivity of `T` as an unmodified predicate.

The graph on the original `[n]` carrying free slots is `K_n minus P_n`;
compressing its *edge indices* does not identify any vertices of `T`. At
`n=4` that carrier has edges `13,14,24` and is itself the path `3-1-4-2`.
It is a particular copy of `P4`, whereas each bit word selects a subset of
its three edges. It is not an arbitrary three-edge graph on four vertices.

## Gauges, symmetries and exact losses

Changing the original names of the `v_i` changes nothing after path
normalization. Changing the chosen path changes the coordinates. On a fixed
unmarked tournament `T`, forgetting `P` has fibre `H(T)/|Aut(T)|`, by
LEM-003: an automorphism fixing the ordered spanning path fixes every vertex.
The graph order cannot instead be forgotten. Three one-edge graphs on `[3]`
are isomorphic as ordinary graphs but lift to three different tournament
types; see the table below. For source orders at most three, the target has
at most two vertices and this failure cannot occur. Thus order four is the
minimal witness.

Let `r(a)=n-a` on `[n-1]`. The operation

```text
(T,P) -> (T^op, reverse(P))
```

induces `G -> r(G)` and `U -> r(U^op)`. Indeed
`x'_ij=x_(n+1-j,n+1-i)`, and compressing gives
`{i,j-1} -> {n-(j-1),n-i}`. A path-reversing anti-automorphism of `T`
therefore exists exactly when `G` is invariant under this specific order
reflection. Existence of some ordinary graph automorphism is a weaker,
different assertion.

Flipping all non-path arcs gives graph complementation and `U -> U^op`.
It retains the forced path in `T`, so it is not full tournament converse
on `T`. Arbitrary graph relabellings act by pure edge permutation; arbitrary
tournament relabellings in reversal coordinates have the additional gauge
mask `epsilon_g`. The two unlabelled quotients cannot be identified by
silently dropping this mask.

## Intransitivity and strong connectivity as graph predicates

Define `x_(i,i+1)=0`. For `i<j<k`, the triangle on those original vertices is
cyclic exactly when

```text
(x_ij, x_jk, x_ik) is (0,0,1) or (1,1,0).
```

This follows by listing the two directed cycle orientations. It is an
order-sensitive predicate, not existence of an ordinary graph triangle.
In the unrestricted fixed-order encoding of a three-vertex tournament,
a single long reversal yields the directed cycle, while three reversals
give a transitive tournament; the latter is not a path-normalized source
with the same path order. In the marked source `T`, transitivity holds
exactly when `G` is empty: a backward chord
closes a directed cycle with the corresponding segment of the fixed path,
whereas no backward chord leaves the transitive order.

There is also an exact **interval-cover** criterion. An edge `{a,b}` of `G`,
with `a<b`, represents the backward chord `(b+1) -> a`, which crosses
exactly the path cuts `k=a,...,b`. For `n>=2`,

```text
T is strong  iff  union_({a,b} in E(G), a<b) [a,b] covers every k in [n-1].
```

Proof: if a cut is uncovered, all crossing arcs point forward, preventing
return across it. Conversely the vertices reachable from `n` form a suffix,
because the forward path lets every reachable vertex reach every later
vertex. If this suffix starts at `r>1`, no backward arc crosses cut `r-1`,
contrary to coverage. Hence `n` reaches every vertex; every vertex reaches
`n` along the fixed path, proving strong connectivity. This is the same
cut mechanism underlying THM-4729 Proposition F, expressed in compressed
graph coordinates. In fact the reachability argument needs only the fixed
directed path and the recorded backward arcs.

For `n>=3`, the marked path closes to a Hamiltonian cycle exactly when
`{1,n-1}` is an edge of `G`. This is sufficient for the interval-cover
criterion but not necessary: at `n=4`, edges `12,23` cover all cuts without
edge `13`. Thus the order data transfers several exact tournament predicates
while graph connectedness by itself is not the transferred condition.

## The complete order-four table

The bit order is `(x13,x14,x24)`; compressed edges are respectively
`12,13,23`. Scores are the sorted outdegrees of `T`. The strong four-vertex
tournament is unique up to isomorphism.

| Bits | Edges of G on [3] | T type | Scores | Cyclic triangles | H/Aut | U type |
|---|---|---|---|---:|---:|---|
| 000 | none | transitive | 0,1,2,3 | 0 | 1/1 | transitive |
| 100 | 12 | C3 -> sink | 0,2,2,2 | 1 | 3/3 | transitive |
| 010 | 13 | strong | 1,1,2,2 | 2 | 5/1 | C3 |
| 110 | 12,13 | strong | 1,1,2,2 | 2 | 5/1 | transitive |
| 001 | 23 | source -> C3 | 1,1,1,3 | 1 | 3/3 | transitive |
| 101 | 12,23 | strong | 1,1,2,2 | 2 | 5/1 | C3 |
| 011 | 13,23 | strong | 1,1,2,2 | 2 | 5/1 | transitive |
| 111 | 12,13,23 | strong | 1,1,2,2 | 2 | 5/1 | transitive |

The four unmarked tournament fibres have sizes `1,1,1,5`, as inherited from
LEM-003 and the existing tiling model. This is separate from deletion-map
fibres `4,2,1` in other path-family constructions.

## The three three-edge graph types

This classification and its complement action are inherited from
[`defect_boundary_20261020/web.md`, section 4](../defect_boundary_20261020/web.md).
A disconnected simple graph on four vertices and three edges must be
`K3 + K1`: a partition into parts of size at most two permits at most two
edges. A connected such graph is a tree, and its possible degree sequences
are `(3,1,1,1)` and `(2,2,1,1)`, giving `K1,3` and `P4`. On fixed labels:

```text
P4:      4!/2 = 12,
K1,3:    4 choices of centre,
K3 + K1: 4 choices of isolated vertex.
```

These exhaust `binom(6,3)=20`. Complementation preserves the `P4` type and
exchanges the other two types; on labelled objects it makes ten disjoint
pairs, six path/path pairs and four star/triangle pairs. This is the exact
complement action on this universe; it alone does not provide the separately
requested `3-4-3` or `4-2-4` decompositions of ten.

## Reproduction and controls

From the repository root:

```text
python 05-knowledge/results/trijection_reset_20261020/independent_check.py --output 05-knowledge/results/trijection_reset_20261020/independent_check.json
```

The implementation enumerates all path-normalized tournaments at `n=1..4`
and all twenty labelled three-edge graphs on `[4]`; no extra filters apply.
Positive controls are inverse-map round trips, transitive and cyclic
triangles, graph complement degree reversal, and `H/Aut` counted from all
vertex permutations. Hostile controls are the three one-edge graph lifts.
Direct arc reachability independently checks interval coverage; direct
triangle enumeration also agrees with the score formula
`C(n,3)-sum_v C(outdegree(v),2)`. The proof establishes the all-order maps
and criteria; the finite check does not replace those proofs.
Every check raises an explicit `RuntimeError` on failure and remains active
under `python -O`; normal and optimized runs produce identical receipts.
