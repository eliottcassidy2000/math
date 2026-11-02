# Session 141: functional capacity separates three unit-distance questions

**Status: elementary theorems proved here; exact finite controls independently
recomputed. No novelty, Collatz convergence, global unit-distance extremality,
or Hadwiger–Nelson improvement is claimed.** This packet recovers prior work at
Git commit `c190e11b7795ed4e2d2f3aeb667e86d08238cdbf`.

## 1. The recovered question and the preservation boundary

The repository's small-value “unit-distance” work means the **Erdős maximum**

\[
u(n)=\max_{P\subset\mathbb R^2,\ |P|=n}
       |\{\{p,q\}\subset P:\|p-q\|=1\}|.
\]

This differs from the **Hadwiger–Nelson** question about coloring every edge of
the full plane unit-distance graph. It also differs from embedding an abstract
functional graph, or routing its edges through a fixed honeycomb lattice.
“Planar unit-distance realization” means points in the plane; edges may cross.
It does not mean a planar graph.

For standard Collatz, use C(n)=n/2 when n is even and C(n)=3n+1 otherwise.
After forgetting orientations, deleting loops, and collapsing repeated edges,
every finite subgraph is a pseudoforest (at most one cycle per component).
Indeed, charge each edge to a distinct directed source: on every vertex subset
there are at most as many edges as vertices. Two cycles joined by a path would
violate that count. The known 1–4–2–1 triangle makes the full graph exactly
three-colorable. This does not assume convergence.

A **topological routing** replaces each edge by a finite path, with disjoint
interiors and no branch vertex in an interior. Its selected route graph is a
subdivision. A mere map with paths sharing interiors is not covered by the
subdivision statements below. Moreover, a selected honeycomb route graph need
not be induced: unused honeycomb edges between selected vertices may create
extra Euclidean unit pairs. Those are extra geometric constraints, not
additional Collatz transitions. The full induced unit-distance completion
still lies in the bipartite degree-three honeycomb graph.

## 2. Exact loss when compressing a graph to one successor system

**Theorem 1 (classical pseudoforest-rank formula, reproved here).** For a finite simple graph G with n
vertices, m edges, and t tree components (isolated vertices count as trees),
the minimum number of edges to delete to obtain a pseudoforest is

\[
\varepsilon(G)=m-n+t.
\]

Equivalently it is the minimum loss when retaining edges that admit an
orientation with outdegree at most one. Arbitrary finite edge subdivisions
preserve this defect exactly.

**Proof.** Within each original tree component on v vertices, at most v−1
edges can be retained. Within each cyclic component at most v can be retained
in a pseudoforest, even if deletions disconnect it. Thus at most n−t edges
can remain. Attain this by retaining the whole of each tree, and a spanning
tree plus one additional edge in each cyclic component. Orient every tree
toward a root; orient a unique cycle cyclically and its attached trees toward
that cycle. This gives the required outdegree bound. Conversely an outdegree
at-most-one orientation supplies at most |W| edges on each W and therefore
cannot contain a component with two cycles.

Subdividing an edge adds one vertex and one edge and changes no component
from tree to cyclic or conversely. Hence m−n and t are unchanged. Iterate. ∎

The ordinary cycle rank is β(G)=m−n+c, with c all components, so
ε(G)=β(G)−(c−t). A functional component is allowed one independent cycle,
not an arbitrary number of cycles merely because a picture uses hexagons.
In particular a subdivided finite Collatz subgraph still has m≤n; replacing
one transition by many unit honeycomb steps cannot generate edge density.
Its density is (m+s)/(n+s), not m/n, where s is the number of inserted vertices.

## 3. A precise meaning of the previous 3n threshold

A **functional layer** is an edge set admitting an outdegree-at-most-one
orientation. A k-layer representation partitions all edges into k such sets;
the orientations may differ between layers. These are arbitrary partial
successor maps, not necessarily arithmetic Collatz maps.

**Theorem 2 (Hakimi capacity criterion, reproved here).** For an integer k≥0,
G has k functional layers if and only if

\[
|E(G[W])|\le k|W|\qquad\text{for every }W\subseteq V(G).
\]

This is the classical bounded-outdegree orientation / pseudoarboricity
criterion, associated with S. L. Hakimi, *On the degrees of the vertices of a
directed graph*, Journal of the Franklin Institute 279(4), 290–308 (1965),
https://doi.org/10.1016/0016-0032(65)90340-6 . The Hall proof below is included
for self-containment; the abstract theorem is not new. The scoped contribution
here is its application to the previously studied unit-distance threshold.

**Proof.** Necessity follows by charging each layer's edges in W to their
sources. For sufficiency form a bipartite graph: a left vertex for every edge
e, and k right-hand slots (v,1),…,(v,k) for every original vertex v. An edge
is adjacent to all slots of its two endpoints. For any collection F of left
vertices, let W be all their endpoints. Its neighbor set has size k|W|, and
|F|≤|E(G[W])|≤k|W|. Hall's matching theorem therefore assigns each edge a
distinct endpoint slot. Orient it from that endpoint and assign it to that
slot's layer. Each layer has at most one outgoing edge per vertex. ∎

**Corollary 3 (extremal threshold equals a dynamical capacity threshold).** For an integer k≥1 put

\[
N_k=\min\{n:u(n)>kn\}.
\]

When this set is nonempty, N_k is exactly the minimum number of vertices in a
plane unit-distance graph that cannot be represented by k functional layers.
For n<N_k every induced subset W has at most u(|W|)≤k|W| edges, so Theorem 2
applies. At n=N_k an extremal configuration violates its full-set inequality.
Thus the old 3n crossover is exactly the first failure of **three** successor
layers. This is a capacity theorem, not a numerical analogy to a kissing
number, and it does not identify those layers with Collatz dynamics.

Using the dated AMP baseline specified below gives N_1=4, N_2=11, and
25≤N_3≤28. In particular, n below the global density crossover cannot hide a
smaller induced density obstruction. Conversely m≤kn alone does **not**
suffice for arbitrary graphs: every subset must be checked.

## 4. Small discrete values and quantitative obstruction

Source: Alexeev–Mixon–Parshall, *The Erdős unit distance problem for small point
sets*, arXiv:2412.11914v2, 12 February 2025, Theorem 1 and Table 1:
https://arxiv.org/html/2412.11914v2 . Their computer-assisted extremality proof
is **cited, not rerun**. These are a dated verified baseline, not a claim that
no later improvement exists.

| n | AMP baseline u(n) | edges necessarily lost to one functional system, at least | cycle rank, at least | maximum degree, at least |
|---:|---:|---:|---:|---:|
| 15 | 37 | 22 | 23 | 5 |
| 20 | 54 | 34 | 35 | 6 |
| 24 | 68…72 | 44 | 45 | 6 |
| 25 | 72…78 | 47 | 48 | 6 |
| 26 | 76…84 | 50 | 51 | 6 |

For each row, with lower bound L, the last three columns follow respectively
from ε≥m−n≥L−n, β=m−n+c≥L−n+1, and Δ≥ceil(2m/n). They apply to any realization
with at least L unit pairs, whether or not it is globally extremal. They do
not require connectedness. For a specified graph the first loss is exactly
m−n+t, which can exceed the tabulated lower bound.

The baseline is exact through n=21; u(22)≤61, u(23)≤66, u(24)≤72, while
u(28)≥85>84. These facts give 25≤N_3≤28. The exact small entries have
u(n)≤n for n≤3 and u(4)=5, and u(n)≤2n for n≤10 and u(11)=23, giving the
first two capacity thresholds.

**Two distinct honeycomb barriers.** First, none of these dense graphs is a
subdivision of a functional graph, by Theorem 1. Second, none can itself have
any subdivision embedded as a honeycomb subgraph: subdivision retains every
original branch degree, but the displayed Δ lower bounds exceed three.
This latter obstruction holds even without asking for functional dynamics.

For comparison, every finite honeycomb vertex set of size n≥4 has at most
floor((3n−6)/2) edges. A connected cyclic planar graph of girth at least six
satisfies this by Euler's formula and face lengths; sum over cyclic components
and tree components. For a forest use m≤n−1≤3n/2−3. Thus the respective upper
bounds are 19,27,33,34,36. These are safe bounds, **not claims of sharp finite
honeycomb edge maxima**. A chosen route graph has the stronger m≤n bound;
its induced honeycomb completion need not.

## 5. Density can grow without increasing chromatic difficulty

**Theorem 4 (faithful Cartesian products).** Given nonempty finite point sets A,B with
their full unit-distance graphs G,H, all but finitely many relative rotations
Rθ produce an injective sum set A+RθB whose full unit-distance graph is G□H.
Consequently χ(G□H)=max(χ(G),χ(H)).

**Proof.** A collision between two labeled sums is u+Rθv=0, with
u=a−a' and v=b−b'. If exactly one difference vanishes, collision is impossible;
if both are nonzero, at most one angle modulo 2π can satisfy it. For unwanted
unit pairs with both coordinates changing, the squared norm is

\[
|u+R_\theta v|^2=|u|^2+|v|^2+2u\cdot R_\theta v.
\]

Because u,v are nonzero, the final term is a nonconstant sinusoid (amplitude
2|u||v|), so equality to one excludes at most two angles modulo 2π. There are
finitely many labeled pairs. Pairs differing in one coordinate have exactly
the intended inherited distances. Avoid all forbidden angles.

For k=max(χ(G),χ(H)), take proper colorings with labels in Z/kZ and color
(g,h) by their sum. Along a Cartesian edge only one summand changes, by a
nonzero residue, proving χ≤k. Each factor occurs as a subgraph, proving ≥. ∎

**Corollary 5.** For every integer d≥1, H_d=K3□…□K3 (d factors) has a faithful plane unit-distance
realization with n=3^d, m=d·3^d, degree 2d, χ=3, and exactly d functional
layers. Orient each coordinate triangle cyclically: its coordinate supplies
one outgoing edge at every vertex, and the d layers partition all edges.
The whole-set inequality m≤kn rules out fewer than d layers.

Thus induced plane unit-distance graphs of chromatic number three have
**unbounded edge density**. This family has m=n log_3(n); it does not give a bound
u(n)≥c n^(1+ε) for any fixed constants c,ε>0. In particular H_3 has 27 points, 81 edges,
and only three colors; H_4 has 81 points, 324 edges, and still only three colors.
The value 81 at 27 is a realized lower bound, not a proof of u(27)=81. Neither
high density nor many functional layers alone produces a new Hadwiger–Nelson
lower bound. A faithful honeycomb route graph has at most two colors; its
suppressed original Collatz graph has three. Subdivision is not a
chromatic-lower-bound transfer operation.

## 6. Exact geometric controls and recovered corrections

`check.py` uses only integer and Fraction arithmetic. It verifies all pairs,
not just claimed edges, for H_d with d=1,2,3,4. The triangle is
(0,0),(1,0),(1/2,√3/2); factor rotations are respectively
(1,0),(3/5,4/5),(5/13,12/13),(8/17,15/17). Coordinates lie in Q(√3)^2.
No tolerance or numerical optimizer is used.

It also independently recounts the repo's corrected lattice 21-point witness. With
w=e^(iπ/3), α=2+w, barα=3−w, T={0,1,w}, and
W={0,±1,±w,±(w−1)}, the 21 points

\[
(\alpha p+\bar\alpha q)/\sqrt7\quad(p\in T,q\in W)
\]

are distinct and have exactly 57 unit pairs. Before scaling, the exact norm
is Q(a,b)=a²+ab+b² and the unit shell is Q=7. Their minimum Q-distance is 1,
so after scaling some pairs are closer than one. This is **not** a penny
configuration. The old THM-431 statement that every triangular-lattice set
on 21 points has at most 47 unit pairs is false; its own October 6 correction
provides this counterexample. The restriction to the nearest-neighbor shell
cannot be silently extended to all scaled lattice layers. The same source's
unqualified “first crossover at 28” is not established by lower bounds; the
valid dated bracket is 25…28. The subsequent THM-461 deletion argument is
compatible with this bracket. We do not import its informal rigidity or
numerical-confluence interpretations.

The checker additionally enumerates every labeled simple graph on 0…5 vertices,
checks the exact one-layer deletion defect against an independent matching
algorithm, checks the capacity criterion for k=1,2,3 against every vertex
subset, and checks single-edge subdivision invariance. Deliberately invalid
collisions, rotations, counts, and capacity claims must fail. Explicit
exceptions rather than Python assertions enforce all acceptance checks.

## 7. Replay, portability, and limits

Run from this directory:

```
PYTHONDONTWRITEBYTECODE=1 python3 check.py > /tmp/unit_distance.normal.json
PYTHONDONTWRITEBYTECODE=1 python3 -O check.py > /tmp/unit_distance.optimized.json
cmp /tmp/unit_distance.normal.json /tmp/unit_distance.optimized.json
cmp /tmp/unit_distance.normal.json replay.normal.json
sha256sum -c SHA256SUMS
```

`PUBLIC_FILES.json` is the exact export allowlist. `SHA256SUMS` hashes all other
allowlisted files, including the selector. `SOURCE_PINS.json` records exact
repository identities and the primary external source. No repository-wide
search of blob contents, forbidden central-history read, or source mutation
is needed for replay.

Finite controls support the proofs; they do not replace them. No full AMP
enumeration, broad search for better small values, or Lean build of any
external Hadwiger–Nelson claim was performed. Later web claims about sharp
small values were not audited or used. No current sharp plane chromatic
interval is needed by any theorem here. The proper induced Euclidean Collatz
placement and unbounded honeycomb route-length results remain Session140
inputs, not newly proved conclusions of this packet.
