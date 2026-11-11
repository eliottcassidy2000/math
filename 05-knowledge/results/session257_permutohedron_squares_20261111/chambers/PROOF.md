# Square-sum Hamiltonian paths as singleton braid chambers

**Session257. PROVED in this note**, with finite exact small-graph controls.
The underlying facts about graphical zonotopes and linear extensions are
standard elementary constructions, not literature-priority claims. The new
session contribution is the explicit square-sum/tournament comparison and
its operation/loss ledger, together with the local gluing charts in the
[companion note](../permutohedron/PROOF.md).

## 1. The graph quotient of the braid arrangement

Let G=(V,E) be a finite simple undirected graph, |V|=n>=2. A word
sigma=(v_1,...,v_n) orders V. Orient every G edge from its earlier to its
later vertex, obtaining the acyclic orientation D_sigma. Let L(D) be the
set of words inducing D. Equivalently these are linear extensions of the
reachability partial order of D. Every acyclic orientation has such words.

The permutohedron Pi_V has these n! words as vertices and adjacent swaps as
edges. Delete its edges whose swapped label pair belongs to E, retaining
only NONEDGE swaps. Then:

**Theorem 1.** The connected components of this remaining graph are exactly
L(D), for D ranging over acyclic orientations of G. Moreover the following
are equivalent for a word sigma:

1. sigma is an oriented Hamiltonian path of G;
2. sigma is an isolated vertex of the nonedge-swap graph;
3. L(D_sigma)={sigma};
4. D_sigma has a total reachability order.

**Proof.** An adjacent swap changes the induced orientation exactly when the
swapped pair is an edge of G. Thus each remaining component lies in one fiber.
To connect two linear extensions, move the first letter of the desired word
to the beginning of the other word by adjacent swaps, then induct. Every
crossed letter is incomparable to it: being a predecessor would contradict
the desired first position, and being a successor would contradict the
current linear extension. Incomparable vertices have no G edge, so these
are permitted swaps. This proves fiber connectivity.

A word is isolated precisely when each consecutive pair is an edge of G,
which is precisely Hamiltonicity. A connected fiber contains an isolated
vertex iff it consists of that vertex. Finally a finite partial order has a
unique linear extension iff it is total: an incomparable pair admits
extensions in both relative orders. Alternatively, every consecutive pair
in its unique extension must be directly comparable by a G edge; a longer
D-path between consecutive positions would require an intermediate position.
QED.

For an undirected graph every path has two oriented words, so the number of
singleton fibers is 2h(G), where h counts unoriented Hamiltonian paths.
No existence is obtained merely by constructing the quotient.

## 2. Graphical zonotope and exact cone mass

Set

    Z_G = sum_{ {u,v} in E } [e_u,e_v]  in R^V.

For a generic vector x with distinct coordinates, maximizing <x,z> selects
one endpoint from every edge segment according to the larger x coordinate.
The maximizing vertex records D_x, up to the harmless global choice of
orienting toward larger rather than smaller values. Distinct acyclic
orientations have different vertices: choose an x inducing one of them;
any changed segment endpoint strictly lowers the value.

The relative interior of the normal cone at z_D is the union of the braid
chambers for words in L(D), together with their internal walls. Its
Gaussian cone probability, or spherical probability modulo the all-ones
line, is

    mu(C_D) = |L(D)| / n! .

Indeed the i.i.d. continuous exchangeable coordinates of a standard Gaussian
have every strict order with equal probability. Internal hyperplanes have
measure zero. This also works when G is disconnected: cones have additional
lineality, but the Gaussian chamber partition remains meaningful.

**Corollary 2.** Oriented Hamiltonian paths of G correspond exactly to the
vertices of Z_G whose normal cone consists of a single braid chamber,
equivalently whose cone probability is the smallest possible 1/n!.

Thus Q_n Hamiltonicity can be phrased as the existence of a minimally sized
normal cone in its labelled graphical zonotope. This is a precise reformulation,
not a faster search or an unproved claim that such a cone must exist. For Q25,
the inherited ten unoriented paths give exactly20 singleton cones, each of
mass1/25!. The [new local chart certificate](../permutohedron/certificate.json)
records all ten paths without enumerating25! words.

## 3. Two different native uses of the SAME permutohedron

Session210's `session210/coherent-permutohedron/PROOF.md` associates a
coherently oriented Pi_V to a tournament T. The direction of an adjacent-swap
edge depends only on the comparison of its swapped labels in T. Its sinks
are precisely directed Hamiltonian words of T. This is an inherited theorem.

For the undirected square-sum graph the present construction instead keeps
nonedge swaps UNORIENTED, and its singleton components are the Hamiltonian
words. The two maps are now explicit:

- **Tournament:** comparison on every pair -> coherent direction on every
  permutohedral edge -> sink set.
- **Square-sum:** arithmetic adjacency subset -> deletion of its pair classes
  from the swap graph -> singleton-component set.

They are not interchangeable. For a tournament every permutation pair class
is oriented; for a square graph absent adjacency must be forbidden, not
oriented by an arbitrary tie-break. Turning nonedges into artificial arcs
changes the count and creates paths that never existed in the square graph.

Both constructions admit the common *local legality mask*: at a permutation,
inspect each consecutive pair. For T it must point forward; for Q_n its sum
must be square. But a common Boolean predicate format supplies neither an
H-preserving graph map nor a Chan–Pak realization. The statement keeps the
actual comparison/adjacency data, exactly as Session256 keeps cut order and
named ports rather than only a scalar count.

## 4. The abstract polytope is not enough

Take the path P4 and star K1,3 on four vertices. Both have chromatic polynomial
q(q-1)^3, exactly8 acyclic orientations, and a graphical zonotope combinatorially
a3-cube. Nevertheless their braid-refinement fiber spectra are:

    P4:   two fibers of size1, four of size3, two of size5;
    K1,3: six fibers of size2, two of size6.

Thus h(P4)=1 and h(K1,3)=0. The first spectrum contains minimal cones, the
second does not. Chromatic data, total number of orientations and the abstract
face lattice all forget the crucial alignment of the normal fan with the
labelled braid chambers. These exact spectra are independently produced by
`check.py` through all24 words for each graph.

The correct geometric state is therefore Z_G WITH its labelled braid-fan
refinement, not just its combinatorial polytope type. For the Monster-action
construction in [the filtration note](../filtration/PROOF.md), the analogous
lost data are the orbit labels and singleton endpoint marks. These are two
different, precisely identified losses, not a claim of one universal sidecar.

## 5. The original colors and the local path-cover atlas

On Q25 use the original birth-role edge colors recovered in
[the source atlas](../history/README.md). Once a singleton chamber has been
identified, attach its word's color monomial or ordered color trace. Reversal
pairs the singleton chambers and reverses the color word. This is a well-defined
decoration on actual solutions; it does NOT identify the whole cone partition.
Session40 already proves the full multiaffine birth polynomial and shows
ordered color word plus endpoints can lose identity at36. Nothing here
relabels color information as sufficient incidence data.

For an actual spanning linear forest F, every Hamiltonian word containing F
is a concatenation of its path components, each traversed in one of its two
orientations (singletons have one). This is the exact local chart built in
[the permutohedron note](../permutohedron/PROOF.md): block position sums yield
a weighted permutohedron, while real endpoint seams select the feasible
vertices in its orientation sheets. The chart is a constrained family of
singleton chambers of the global quotient. Refining F by cuts enlarges the
possible interleavings and gives the proved coordinate projection there.

One must not conclude that the global feasible set is connected by adjacent
swaps because its ambient object is a permutohedron. The Q25 experiment has
all10 unoriented feasible chart vertices frozen under both elementary block
swaps and single-block reversals, while two-edge surgery connects them.
The difference between geometric adjacency and a useful path operation is
an explicit result of the chart, not a defect of the representation.

## 6. Finite controls and next questions

Run `python3 -B check.py`. It exhausts all permutations of Q2,...,Q8, P4,
K1,3, C4, K5 and the empty graph on5 vertices. Every fiber's connectivity,
nonedge-swap criterion and singleton/Hamilton equivalence are checked; the
fiber mass is n!. Explicit RuntimeError checks remain active under python -O.
The general proof is independent of these finite controls.

Open targets specific to this formulation:

1. **Construct minimal cones without full enumeration.** Can the reflection
   and path-cover grammars select a singleton cone with controlled named
   endpoints and original birth-color trace? Existence at large n is already
   cited; the target is extra structure, not rediscovering general existence.
2. **Compare surgery to fan moves.** Characterize the smallest refinement of
   a path-cover chart that makes a chosen r-edge exchange into a local move.
   The companion common-component theorem gives a precise starting point.
3. **Operation-compatible comparison.** Find an actual map from a selected
   singleton-cone bank to coherent tournament sinks with a controlled change
   of H. Merely sharing Pi_V or both having Hamilton words does not do this.

Chan–Pak remains OPEN. The Monster shell geometry, square-target triple
compiler, positive two-seam kernel and this chamber quotient are compatible
pieces of a research program, not a proved equivalence of their large
problems.
