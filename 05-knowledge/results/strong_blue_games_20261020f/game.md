# The complementary Hamiltonian-path / spanning-walk game

**PROVED:** with Player 2 choosing the starting vertex, the three variants
that mention a Hamiltonian path admit a winning coloring exactly when
`n >= 4`; the walk/walk variant never does. If both colors must be connected,
the Hamiltonian/Hamiltonian variant has witnesses for every `n >= 7`.
**FINITE-EXACT:** the complete graph atlas excludes such a connected pair
through order 6 and gives exactly one complementary pair at order 7.

## 1. Inheritance and exact rules

The useful inherited mechanism is the path-cover sidecar, rather than a
scalar path count. **THM-1975**,
[the path-cover polynomial is the refined compositional invariant](../../../01-canon/theorems/THM-1975-the-path-cover-polynomial-is-the-refined-compositional-invariant.md),
concerns directed substitution and is marked VERIFIED; no all-order claim
from that file is used below. Section 3 proves the needed undirected join
criterion independently. The canonical hostile is the claw: it is connected,
has a spanning walk, and has no Hamiltonian path. The corrected near miss is
to confuse a vertex-spanning walk with an edge-exhausting Euler trail. The
least-used coordinate is the number of components after deleting a small
vertex set.

The live board is: complement connectivity; path-cover number; separator
components; endpoint choice; shortest walk length; marked torus winding.
The first three control the game; the final coordinate requires extra
geometric data, as documented in [torus.md](torus.md).

Fix `n >= 1` and an undirected simple complete graph on the **same full set
of n vertices**. Player 1 selects edges G; the unselected edges form its
complement Gbar. Isolated vertices remain part of both graphs. For a graph X,

```text
H(X): X has a Hamiltonian path, visiting all vertices exactly once;
W(X): X has a walk visiting all vertices, with repeats allowed.
```

Player 1's prescribed predicate is A and Player 2's is B, each H or W.
The final selection is legal when `not A(G)` and wins when also
`not B(Gbar)`. Player 2 chooses the start. This is a static edge-selection
game, not an additional requirement that Player 1 select edges by walking.
At n=1 the zero-edge traversal exists, so Player 1 has no legal selection.

## 2. Complete answer to the four base variants

| Forbidden in Player 1's color | Required in Player 2's color | Winning selection exists |
| --- | --- | --- |
| H | H | exactly n >= 4 |
| W | H | exactly n >= 4 |
| H | W | exactly n >= 4 |
| W | W | never |

**Proof.** A spanning walk exists exactly when the graph is connected.
Necessity is immediate; sufficiency follows by traversing a spanning tree,
repeating edges when needed. If G is disconnected, two vertices in different
components are adjacent in Gbar, and two vertices in the same component
have a two-edge route through a vertex in another component. Thus Gbar is
connected. This proves the walk/walk impossibility.

For n>=4, take `G=K_(n-1) disjoint-union K_1`. Its complement is the star
`K_(1,n-1)`. The first graph is disconnected; the star has at least three
leaves. A Hamiltonian path can have at most two leaves of the ambient graph,
since each leaf must be a path endpoint. Thus this coloring wins H/H and W/H.
Swap the colors to win H/W.

For n<=3, every nontraceable graph has a traceable complement. At n=2 the
only nontraceable graph is empty. At n=3, nontraceability means at most one
edge, and its complement has at least two edges, hence a spanning path.
The n=1 boundary was stated above. Since H implies W, this also excludes
the two mixed variants below order 4. This proves all four thresholds.

At n=4 the only winning H/H isomorphism classes are the triangle plus an
isolate and the claw, exchanged by complement. No large census is needed
for the all-order statement.

**Start boundary.** Under the different rule that Player 2 must succeed from
every prescribed start, H/H and W/H have threshold 3: one edge plus an
isolated vertex leaves a P3, which cannot be traversed Hamiltonianly starting
at its middle vertex. This does not defeat the actual rule, since Player 2
can start at an endpoint. Allowing repeats makes a connected graph traversable
from every start, so the walk-target cases do not change.

## 3. Exact path-cover criterion for the disconnected-color case

For an undirected graph Y, let `p(Y)` be the minimum number of disjoint paths
covering all its vertices, with single vertices allowed. Let X be the complete
join of nonempty induced blocks `X_1,...,X_r`: every pair in different blocks
is adjacent. Put `n_i=|X_i|`, `p_i=p(X_i)`, and `N=sum n_i`.

**PROVED.**

```text
X has a Hamiltonian path  iff  p_i <= N-n_i+1 for every i.       (1)
```

Necessity: the vertices of block i appear along a Hamiltonian path in at
most `N-n_i+1` runs. These runs are a path cover of that block.

For sufficiency, initially select a minimum path cover in every block.
Any number `k_i` between `p_i` and `n_i` is possible by splitting paths.
If the counts `p_i` satisfy `max p_i <= sum(other p_j)+1`, keep them.
Otherwise there is a unique dominant block i. Split paths in the other
blocks until their total count is `p_i-1`. This is possible by (1), since
their total vertex capacity is `N-n_i >= p_i-1`. The new counts are balanced:
each is at most one plus the sum of the others.

A multiset of labels with largest multiplicity m can be put in a word with
no repeated adjacent label exactly when the other multiplicities sum to
at least m-1. For an explicit construction, place the m copies of a largest
label first. Distribute all other labels cyclically into m gaps after these
copies, putting each label's consecutive allocations into distinct gaps
(each such label occurs at most m times). The first m-1 gaps are nonempty;
within a gap every label is distinct. The resulting word has the required
property. Replace each label occurrence by one path from its block's cover.
All joins between consecutive paths are present, yielding a Hamiltonian path.

Apply (1) to Gbar when G is disconnected: the components of G are its
canonical join blocks after complementation. Consequently the W/H winning
condition has the exact form

```text
G disconnected, and some component C of G satisfies
p(Gbar[C]) > n-|C|+1.                                      (2)
```

The H/W classification is the color-swapped statement. Formula (2) is a
structural iff, although computing a general path-cover number remains a
substantive task. If G is a disjoint union of cliques, each complement block
is independent and `p(Gbar[C])=|C|`; the remaining complete multipartite graph
is traceable iff its largest part is at most `ceil(n/2)`.

The source of this map is the full coloring; its target is the list
`(|C|,p(Gbar[C]))`. It preserves the traceability decision for the complete
join, but loses the number of Hamiltonian paths and their endpoint identities.
Retaining internal path covers is sufficient to construct a traversal.

## 4. Connected in both colors: separator obstructions at order seven

Every traceable graph X obeys

```text
c(X-S) <= |S|+1                                             (3)
```

for every vertex set S, because deleting |S| vertices from a spanning path
leaves at most |S|+1 path segments covering the remaining vertices. This is
necessary, not asserted sufficient for arbitrary graphs.

For every `a>=5`, start with `K_a`, remove one edge uv, and add two new
vertices x,y adjacent only to v. Call the result Q; it has `n=a+2>=7` vertices.
Both Q and Qbar are connected. In Q, deleting v leaves the two new isolates
and the connected old clique minus v: three components after deleting one
vertex, violating (3). In Qbar, the old a vertices induce exactly one edge
uv and the new vertices form an edge. Each old vertex other than v is
adjacent to both new vertices; v connects through u. Deleting x,y leaves
`a-1>=4` components after deleting two vertices, again violating (3).

Thus both colors can fail H even while connected, for every n>=7.
The complete atlas through n=7 has no such graph through n=6 and exactly
two at n=7, a complementary pair (indices 615 and 743). The construction
above is isomorphic to index 743 when a=5. This establishes the sharp
threshold with a finite-exact lower-bound check, rather than claiming a new
analytic classification of all small connected graphs.

When both colors must be connected, W/H and W/W have no legal Player 1
selection; H/W has no winning selection because Player 2's graph is connected.
Only H/H retains a nontrivial obstruction.

## 5. Lengths and the Euler boundary

For a connected n-vertex graph, a spanning tree gives a spanning walk of
length at most `2(n-1)`. Length exactly n-1 is possible precisely when a
Hamiltonian path exists. For a tree T the sharp free-endpoint minimum is

```text
ell(T)=2(n-1)-diam(T);
ell(T; prescribed start s)=2(n-1)-ecc_T(s).                 (4)
```

Indeed, for endpoints s,t each edge on their unique path must be traversed
at least once, and every other tree edge at least twice. Traverse the
off-path branches out and back to attain this bound, then maximize the
s--t distance. For the n>=3 star, the shortest spanning walk has length
`2n-4`, despite the absence of a Hamiltonian path when n>=4.

An Euler trail is another predicate: it must use **every edge exactly once**.
For a connected graph it exists iff the number of odd-degree vertices is
zero or two. If there are two, the start must be one of them. Necessity
comes from entry/exit pairing; sufficiency follows by splicing closed
edge tours, adding and then deleting an edge between the two odd vertices
when necessary. This does not change the four rules above: those ask to
span vertices, not exhaust edges. The explicit torus-grid construction in
[torus.md](torus.md) uses Euler circuits only after declaring the added
orientation and embedding.

## 6. Incoming sharp budgets and the full counting profile

Two companion packets arrived during final integration and were read before
use: [sharp budgets](../kernels_games_detours_20261020f/game.md) and
[rooted budgets / counting](../complement_path_games_20261020/README.md).
Their free-start feasibility agrees with section 2. They strengthen it to
the following minimum numbers of red edges, for every n>=4:

| Variant | Minimum red edges |
|---|---:|
|H/H|n-1|
|H/W|n-1|
|W/H|2n-5, except n=6 gives6 and n=8 gives10|
|W/W|impossible|

The second packet supplies a self-contained degree-obstruction proof,
independently checked here. A nontraceable graph has, for some
`1<=k<=floor(n/2)`, k vertices of degree at most k-1. For a connected
nontraceable graph k>=2. Complementation therefore costs at least

```
g_n(k)=k(n-k)-binom(k,2).
```

This concave quadratic is minimized at an endpoint of the allowed k
interval. The values are n-1 for the unrestricted case, and
`min(2n-5,g_n(floor(n/2)))` for the connected-blue case. The latter simplifies
to the W/H row above. The extremal constructions in the companion proof
have disconnected red, so they satisfy the game's extra legality condition.
This is a proved specialization of classical extremal graph mechanisms,
not a numerical extrapolation or a uniqueness classification.

The same companion extends section 3 from support to counts. Let q_i(c)
count unordered covers of block i by c **oriented** paths; nonsingletons
have two orientations and singletons one. Let w(c_1,...,c_r) count label
words with c_i copies of label i and no equal adjacent labels. Then

```
H_ordered(join_i X_i)
 = sum_c w(c_1,...,c_r) product_i [c_i! q_i(c_i)].          (5)
```

Cut a Hamiltonian path into maximal block runs to obtain the bijection.
The factorial assigns the distinct oriented cover paths to occurrences
of their block label. This is the missing full profile behind the
minimum-path-cover criterion, and it keeps counts separate from existence.

Finally define the edge-colouring polynomial

```
F_G(t)=sum_(vertex permutations pi) t^(red consecutive pairs of pi).
```

Its constant coefficient is the ordered blue HP count; its coefficient of
t^(n-1) is the ordered red HP count. H/H winning means that both vanish.
Complementation gives `F_Gbar(t)=t^(n-1)F_G(1/t)`. Turning a blue edge e red
changes F by `(t-1)D_e(t)`, where D_e counts only permutations using e,
weighted by their other red consecutive pairs. These identities partition
the actual permutation words and need no tournament orientation gauge.
The companion's tests, rather than this packet's game.json, check its rooted
join refinements and minimum budgets on all labelled graphs through six.

### The response polynomial is an invertible path-cover transform

The two incoming constructions combine to give an exact identity. For one
graph G, use the same oriented path-cover counts q_G(c) as above. Then

```
F_G(t) = sum_(c=1)^n c! q_G(c) (t-1)^(n-c),
q_G(c) = [x^(n-c)] F_G(1+x) / c!.                         (6)
```

Proof: expand the product of `1+(t-1)` over the red consecutive pairs in
each permutation. A selected subset of consecutive edges forms a spanning
linear forest, with isolated vertices included. Its c components give
exactly an oriented path cover; contracting them leaves c! possible block
orders. Conversely every such cover and block order produces exactly that
selection in a permutation. A c-component forest has n-c selected edges.
This proves (6) without assuming the selected runs are maximal.

Combining (6) with colour reciprocity gives the explicit complement transform

```
c! q_Gbar(c)
 = sum_(d=c)^n (-1)^(n-d) binom(d-1,c-1) d! q_G(d).       (7)
```

It is an involution, since complementation twice is the identity. In
particular the ordered blue HP count is the alternating evaluation

```
H_ordered(Gbar) = sum_(d=1)^n (-1)^(n-d) d! q_G(d).        (8)
```

Thus the response polynomial and full path-cover profile contain exactly
the same information. The minimum path-cover number is only their support
boundary. The source/target are two coordinates on the same undirected
graph invariant; this does not identify that graph with a tournament.
[cover_response.py](cover_response.py) independently compares forest
enumeration with permutation weights and checks the complement transform.
Its complete universe is all1,099 labelled graphs on1..5 vertices, with
59,809 forest-subset examinations and matching normal/optimized receipts.
The claw and triangle-plus-isolate are a sharp information-loss control:
both have `q=(0,6,6,1)` and `F=12t+12t^2`, yet they are nonisomorphic
complements. A profile fixed by the complement involution does not imply
that the graph is self-complementary.

## 7. Reproduction and scope

Run from the repository root:

```powershell
python3 05-knowledge/results/strong_blue_games_20261020f/game.py
python3 -O 05-knowledge/results/strong_blue_games_20261020f/game.py
```

[game.json](game.json) records the script hash, all class counts, and explicit
edge lists / graph6 encodings for the minimal connected pair. The complete
universe is NetworkX's 1,253 Read--Wilson atlas graphs through seven vertices;
its completeness is documented by the
[official graph_atlas_g reference](https://networkx.org/documentation/stable/reference/generated/networkx.generators.atlas.graph_atlas_g.html).
The checker omits the empty n=0 graph from the game. It independently checks
all 208 graphs through n=6 by literal permutations, verifies the join
criterion on all 256 nonempty atlas graphs whose complement is disconnected,
and checks (4) by shortest-walk state-space BFS on all 25 atlas trees,
including all 142 prescribed starts. The structural hostile families are
also checked through n=30. No sampled-edge claim or tournament canonicalizer
is used.

The game is undirected. In a tournament every vertex pair has one directed
arc, so this two-color complement operation is neither tournament converse
nor the earlier Hamiltonian-path reversal/complement operation. Transferring
a game obstruction to the tournament metagraph requires an additional map;
none is asserted here.
