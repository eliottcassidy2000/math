# Graceful exchange fibers: a three-component criterion and an unbounded buffer obstruction

**2026-10-19. PROVED elementary statements; FINITE-EXACT controls; global connectivity OPEN.** No external-priority claim. This continues [the exchange packet](graceful_square_exchange_20261019.md), especially its interval-matrix formulation, upper-tail completion, and shortest-edge buffer. It does not refute the packet's global GT1 conjecture, because the obstruction below freezes selected larger colors.

Artifacts: [checker](../../04-computation/experiments/graceful_fiber_obstructions_20261019.py), [output](graceful_fiber_obstructions_20261019.out).

## 1. The state and the inheritance

Put m=n−1. Vertices are the fixed labels 0,...,m, and edge {a,b} has color |a−b|. A graceful tree selects exactly one edge of every color 1,...,m and is connected. Its code is (a_1,...,a_m), where the color-d edge is {a_d,a_d+d}. An exchange changes one coordinate, preserving the tree.

The closest mechanism is the fundamental-cycle exchange criterion in the preceding packet. Its upper-tail completion theorem says that every acyclic selection of the colors k+1,...,m extends using colors k,...,1. That is an existence theorem for completions, **not** connectedness of their exchange graph. The canonical hostile from general rainbow trees is the properly three-colored K4 with four isolated rainbow stars. Here the new hostile lies inside the actual integer-difference coloring.

The active objects are fixed upper-tail components, low-color contact graphs, cuts in numerical order, and the distinction between nonempty fibers and connected fibers. No tournament is imposed on ties: the native carrier is the graph of contacts between components.

## 2. An exact quotient for the shortest-edge buffer

Fix all edges except color one. Removing that edge gives two components, with 0 and m in the same component because edge {0,m} is compulsory. Every adjacent pair of integer labels crossing the cut is an allowed replacement. There are a positive even number of such pairs. Thus each nonempty color-one fiber is an even clique.

Equivalently, contract {0,m}. For m≥3 the entire length-one path becomes a cycle C_m. For m=2 it is a two-edge multicycle; multiplicity must be retained. The allowed positions are the boundaries of a cut in this cyclic order.

This proves that no graceful tree on n≥3 is isolated, but different cliques need not be connected while higher colors are frozen.

## 3. Complete classification when two colors remain free

Assume m≥3. Fix an acyclic selection F with one edge of each color 3,...,m. It has three connected components A,B,C. For j=1,2 define Q_j on {A,B,C}: two components are adjacent when some original labels, one in each, differ by j. Discard loops and, for the support graph only, parallel edges.

**Theorem 1.** The graceful completion fiber with F fixed is disconnected if and only if Q_1 and Q_2 are the same three-vertex path. In this case it has exactly two connected components. Otherwise the fiber is connected.

Equivalently, disconnection occurs precisely when some pair of the three fixed components has no pair of labels at distance one or two.

**Proof.** Q_1 is connected because the full length-one graph is a spanning path. Q_2 is connected because the full length-two graph consists of the even and odd residue chains, and F contains a color-three edge joining these chains. Contracting F therefore leaves a connected support graph in both cases. Each Q_j is consequently either a path on three vertices or a triangle.

A completion chooses one contact edge from Q_1 and one from Q_2. It is a spanning tree exactly when their component-pair types are distinct. If both supports are the same path, the two types must be assigned to the two colors in opposite orders. A one-color move cannot change that assignment without making both edges parallel, so there are two components. Within each assignment, all actual original edges of each type can be varied independently; this is a Cartesian product of cliques, hence connected.

For the converse, name the three possible pair types a,b,c. If both supports are triangles, the ordered unequal pairs form a connected six-cycle under changing one coordinate. If one is a triangle and the other is a path, the four feasible pairs form a connected path. If they are different paths, relabel so the supports are {a,b} and {a,c}; their three feasible pairs (a,c),(b,c),(b,a) form a path. Again replacing a pair type by its actual edge choices only expands a vertex into a connected product of cliques and preserves the displayed adjacencies. This proves both directions and the exact count in the disconnected case. ∎

### A minimal-size finite hostile

At n=7, fix

    (a_3,a_4,a_5,a_6)=(0,0,1,0),
    F={03,04,16,06}.

Its components are H={0,1,3,4,6}, B={2}, C={5}. There are no contacts B–C at distance one or two, while both colors connect each of B,C to H. The full fiber consists of six codes (a_1,a_2), split exactly into

    {(1,3),(2,3)}
    {(4,0),(4,2),(5,0),(5,2)}.

The first component assigns color one to H–B and color two to H–C; the second reverses the assignment. Exhaustive checks of every acyclic tail show no disconnected two-free-color fiber for 4≤n≤6. At n=3 the two completions are directly adjacent; n=2 has only one color. Thus this order-seven witness is minimal. The theorem, rather than extrapolation from this example, supplies the all-size classification.

## 4. No fixed band of short colors can serve as a universal buffer

The preceding obstruction extends to arbitrarily large buffers.

**Theorem 2.** For every k≥2 and every m≥2k+2, there is an acyclic fixed choice of all colors k+1,...,m whose graceful completion fiber has at least k connected components under one-edge, same-color exchanges of the free colors 1,...,k.

**Construction.** Put b=m−1 and choose

    F = { {0,j}: k+1≤j≤m, j≠m−1 } ∪ { {1,m} }.

The edge {1,m} has the omitted color m−1. This forest has one nontrivial component

    H={0,1} ∪ ({k+1,...,m}\{b})

and isolated vertices 2,...,k and b. Thus it has k+1 components and admits graceful completions by the upper-tail completion theorem.

Every neighbor of b at an allowed free distance 1,...,k belongs to H: its lower neighbors are at least m−1−k≥k+1, and the only upper neighbor in range is m, at distance one. Consequently, after contracting F, b has only H as a possible neighbor. Every completion must attach b to H with exactly one edge. Two such edges would create a cycle through the already connected H; zero would disconnect b.

The color d of that unique edge is invariant under any permitted single-color exchange. Changing another color to an edge incident with b creates a cycle; changing b's current edge retains its color. Hence different values of d give different exchange components.

All k values are realized explicitly. For each d∈{1,...,k}, add

    color d:       {b−d,b},
    color j<d:     {1,j+1},
    color j>d:     {0,j}.

The latter edges attach each of 2,...,k separately to H, while the first attaches b. The resulting graph is a tree with every color exactly once. This proves at least k components, without any finite search. ∎

The theorem is about a fixed band of edge colors and single-edge moves, not about the size of a reversed path segment or the number of vertices touched by a relabeling. Allowing several colors to change simultaneously, or changing a larger fixed color, may escape the obstruction. The bound m≥2k+2 is sufficient; no optimal threshold is asserted.

In the smallest family instances, exact completion counts are:

| free colors k | n | completions | components |
|---:|---:|---:|---:|
| 2 | 7 | 6 | 2 |
| 3 | 9 | 26 | 4 |
| 4 | 11 | 150 | 8 |

Only the lower bound k is proved in general. These three component counts are not evidence of an asserted exponential law.

## 5. What this changes in the proof strategy

The global one-color graceful-tree exchange conjecture GT1 remains open here. The examples prevent an induction that treats a nonempty completion fiber as automatically connected, or that assumes a universally adequate fixed set of short buffer colors. A successful global path may have to revisit previously selected larger colors.

The concrete next target is **adaptive unlocking**: given two components of a frozen-tail fiber, characterize a larger-color exchange that connects suitable representatives, and find a well-founded global measure for a sequence of such unlocks. A proposed measure must survive the preceding packet's non-star local minima and this note's attachment-color invariant. Neither no-isolated-vertices nor greedy completion suffices.

The three-component test provides a useful exact interface: source data are the fixed forest and the two numeric contact relations; the map contracts forest components; the preserved predicate is validity and connectivity of two-color completions; lost data are multiplicities and original endpoint positions, restored by clique fibers. This uses the repository's component/cut methods without asserting an equivalence to a tournament or to the Graceful Tree Conjecture itself.

## 6. Exact checks and provenance

Run from repository root:

    python3 04-computation/experiments/graceful_fiber_obstructions_20261019.py
    python3 -O 04-computation/experiments/graceful_fiber_obstructions_20261019.py

The outputs agree. The checker enumerates all 2,879 acyclic upper tails leaving two colors free at n=4,...,10, verifies the support criterion against explicit state-graph components, and finds 950 disconnected fibers. It checks the exact six-state order-seven witness, 627 positive certificates from the unbounded construction, and complete fibers for the three displayed examples. The unbounded theorems follow from the proofs, not the sampled family certificates.

Rainbow-tree reconfiguration is an established subject; see the [Maezawa–Yamaguchi 2024 lecture abstract](https://www.en.sdu.edu.cn/info/1158/8226.htm). The abstract announces a sufficient condition but does not give it, so comparison with that condition is still pending. No global priority claim or major-conjecture closure is made.
