# Two cubes, six squares, and explicit native metagraph tree distributions

Status: FINITE-EXACT structural identification and explicit certificates at tournament order4; PROVED inequalities from those certificates. Independent of the attached preprint's global logarithmic theorem. No optimum or literature-priority claim.

## 1. Exact graph and cyclic core

Use the labelled path-preserving complement metagraph: tournaments on {0,1,2,3} are encoded by6 bits in lexicographic unordered-pair order; a set bit means smaller->larger. For every directed Hamilton path P keep its3 arcs and reverse the other3. The resulting undirected graph G has64 vertices and96 edges. Counts are H=1 at24 vertices, H=3 at16, H=5 at24. The24 H=1 vertices are leaves, each attached to a different H=5 vertex.

Delete these leaves. The core C has40 vertices and72 edges. Partition its edges into A (48 edges between H=3 and H=5 vertices) and B (24 edges between two H=5 vertices). Exact enumeration in check.py gives the following structural certificate:

- The A-subgraph has two20-vertex components. Each is a cube with every edge subdivided once:8 H=3 vertices and12 H=5 midpoints.
- The B-subgraph consists of six disjoint4-cycles on the24 midpoints.
- Contract these six squares. Each of the16 low vertices attaches to three different squares. There are eight such triples, each twice; they are exactly the triangular faces of an octahedron on the six contracted vertices.
- The shortest detour avoiding an A edge has length5, while the shortest detour avoiding a B edge has length3.
- The48 automorphisms consisting of vertex relabellings and optional tournament converse act transitively on A and on B.

Every assertion is checked on all vertices/edges. The cubic bipartite suppressed graph of each A component has bipartition4+4 and degree3, hence is K4,4 minus a perfect matching, i.e. the cube. The output lists both components, all six squares and the three opposite pairs identifying the octahedron. The A/B orbit test checks each of the48 automorphisms against the full64-vertex adjacency relation. This is a finite exact structural proof by exhaustion, not an extrapolation to larger tournament orders. Since A components have no bridges and each B edge belongs to a square, this is indeed the bridgeless core.

This connects the native operation to elementary polyhedral geometry without changing its adjacency to arc flips. Contracting the squares loses ports and distances; the uncontracted incidence data must be retained for metric claims.

## 2. Two unconditional lower bounds

All original edges have length and capacity1. Let S_* be the minimum over distributions on actual spanning trees of maximum expected distance stretch. Let C_* be the corresponding minimum maximum expected fundamental-cut load. The two distributions are allowed to differ.

**Theorem.** For G, and equivalently its bridgeless core,

    S_* >=19/8,       C_* >=79/36.

**Proof.** Every spanning tree includes the24 leaves, and its core restriction has39 edges. Let b be its number of B edges and a=39-b its number of A edges. Since the A-subgraph has two components, b>=1. Each included edge has stretch1. The detour certificate forces stretch>=5 on excluded A edges and >=3 on excluded B edges. Therefore

    sum_(e in A) d_T(e) >=a+5(48-a)=84+4b,
    sum_(e in B) d_T(e) >=b+3(24-b)=72-2b.

For a distribution write x=E b. Its maximum expected edge stretch is at least

    max(7/4+x/12, 3-x/12) >=19/8.

The last minimum is attained algebraically at x=15/2; this is NOT a claim that any tree distribution achieves equality in all preceding estimates.

Equivalently, assign adversarial demand1 to every A edge and demand2 to every B edge (total demand96), independently of the unit lengths. Every spanning tree then has weighted distance cost at least (84+4b)+2(72-2b)=228, giving228/96=19/8. This is an explicit finite weighted-response/minimax dual certificate in precisely the paper's demand-versus-length framework: the coefficient2 cancels the unknown number of square edges. It is stronger than a mere unweighted average-stretch estimate and holds tree by tree.

For capacities, double-count edge/path incidences to obtain

    sum_(f in C) E L_T(f) =sum_(e in C) E d_T(e)
                         >=156+2x>=158.

The maximum over the72 core edges is at least158/72=79/36. Leaf edges have load1 and cannot lower this maximum. These lower bounds use no imported embedding theorem. QED.

The all-pairs stretch objective equals the edge objective in the unit graph: one direction includes edges, and the other sums expected bounds along a graph shortest path.

## 3. Explicit upper bounds from symmetry and BFS

**Theorem.** Actual spanning-tree distributions give

    19/8 <=S_*<=157/48,
    79/36 <=C_*<=13/4.

Build BFS trees of the core starting from masks14,19,51, processing neighbours in increasing numeric order; restore all24 leaf edges. For a fixed tree define s_A,s_B as the sums of tree distances over original edges of types A,B, and l_A,l_B as sums of canonical loads over tree edges of types A,B (missing tree edges contribute zero). Exact certificate values are:

    root14: (s_A,s_B,l_A,l_B)=(160,74,160,74),
    root19: (s_A,s_B,l_A,l_B)=(128,122,196,54),
    root51: (s_A,s_B,l_A,l_B)=(168,66,156,78).

The saved output contains all39 core edges of each tree. Distances and loads are computed by explicit path incidence with integer arithmetic.

For a seed tree choose uniformly one of the48 relabelling/converse automorphisms and take its image. Orbit transitivity makes every A-edge expected stretch equal s_A/48 and every B-edge equal s_B/24; similarly expected loads are l_A/48 and l_B/24. All bridge edges remain exactly1.

For distance, mix the root14 orbit distribution with weight29/32 and root19 with weight3/32. The expected sums are157 on A and157/2 on B, giving equal expected stretch157/48 on both edge orbits. The mixture uses at most96 trees; repeated images can be merged. This explicit support claim is not the paper's stronger general support compression bound.

For capacity, the root51 orbit distribution gives156/48=78/24=13/4 on both orbits, with support at most48 trees. Canonical capacities imply all-cut domination for each tree and expected cut upper bound13/4. The resulting fractional tree-path routing has the same13/4 congestion guarantee by the elementary incidence argument. No claim is made about an individual sampled tree meeting the bound, or about electrical energy.

The distance mixture and capacity distribution differ intentionally. These constants are finite-order upper bounds, not asserted optimal. In particular the O(log N) preprint is unnecessary at this order.

## 4. Inheritance, tests, and next questions

Inherited: the exact K_P definition and degree-H identity; session207's distinction between labelled degrees and class-line degrees; session208 native-tree physics and bridge-factorization (leaves are exact; cut domination does not imply energy domination). This packet adds the two-cube/six-square structure and explicit quantitative native distributions rather than another general bridge reduction.

Anchor: the cyclic core's distortion. Niche: polyhedral incidence. Wildcard: symmetry averaging. The decisive tests are full order4 enumeration, orbit equality and integer tree-path counts; all pass. The octahedral contraction is useful for topology but is not a distance-preserving quotient. No claim about the ordinary isomorphism-class/converse-merged metagraph follows without separately carrying multiplicities and ports.

Open: determine exact S_* and C_*, improve the lower bounds beyond local detours, or find an order-general analogue of this polyhedral structure. The eigenvalue/energy counterexample from the prior packet does not obstruct the valid L1 and cut bounds proved here.

Reproduce: `python3 check.py > RESULT.json` from this directory, using only Python's standard library. The universe is all64 labelled tournaments and all24 Hamilton orders, with every free nonpath orientation. No random or floating computation is used in this certificate.
