# Independent bridge-extension and literal-detour audit

Verdict: PASS. The selected 18-by-18 family has 324 distinct tile edges, and every one is a nonbridge in the full converse-merged tournament metagraph. This does not classify every bridge in the full order-10 metagraph.

## Independent computation

`independent.py` does not import or execute the producer's Python modules. It reads the two certificate JSON files, regenerates the full n=5 spine alphabet, and computes canonical types by exhaustive permutations, graph bridges by edge-ID-aware Tarjan traversal, SCCs by transitive closure, Hamiltonian counts by subset dynamic programming, and literal flips/relabelings with integer adjacency bits. Its checks use explicit exceptions, not assertions. Normal and optimized runs give identical results. The producer's eight-file SHA256SUMS verifies.

All 64 base incidence records match independently regenerated values. The base has 12 ordinary source types and 33 exit-rooted target types. Its ordered-edge symmetric square has 78 source vertices, 561 target vertices, and 4096 edges. Exactly 19 edges are bridges. The four source-leaf tiles are 767, 991, 1021, 1023. The remaining 15 are target-leaf edges in this boundary, not automatically global bridges.

The original selected 18-by-18 boundary has 324 edges and exactly eight bridges, reproducing its 316 already-certified nonbridges without relying on older cycle files. All six local alternate paths have correct endpoints, are simple, and avoid the removed edge. Both remaining literal detours pass every check, giving eight repairs and therefore 324 global nonbridges.

## Why the symmetric square models the boundary

A spine tile e defines a tournament X_e. The actual word (e,h) is the ordered join X_e -> X_h with its displayed Hamiltonian spine. Nonspine complementation gives internal masks e XOR 366 and h XOR 366, all backward cross arcs except the forward seam 4 -> 5. This target is strong because its displayed Hamiltonian path closes by 9 -> 0.

Put f=tau(h), where tau takes converse and reverses all five vertex labels. The source's normalized address is the unordered pair of ordinary types of X_e and X_f. This is exact: the ordered-join cut of size five is canonical from the linear order of strong components, and converse swaps the two normalized factors. The target's normalized address is the unordered pair of exit-rooted types of C(X_e) and C(X_f). The code independently checks that this agrees with complementing/conversing the second actual block and rooting it at its entry vertex.

The target's one-outgoing-arc cut of size five is unique. Indeed, for two such cuts A,B, the cut identity gives

  delta(A)+delta(B) = delta(A intersection B)+delta(A union B)+|A minus B| |B minus A|.

If A,B cross and their union is proper, strongness makes the first two right-hand terms at least one, forcing the crossing product to vanish, a contradiction. For disjoint equal-size-five cuts, A union B is the entire ten-vertex set and their crossing product is 25, also impossible. Thus two different size-five cuts cannot occur. An isomorphism must preserve the cut and its unique seam endpoints; a converse-isomorphism exchanges the two normalized rooted blocks. Conversely, rooted block isomorphisms assemble to a target isomorphism. Therefore address equality is a structural classification here, not merely a numerical invariant.

Ordered products are deliberately retained as distinct multigraph edges: (e,f) and (f,e) have the same endpoints but are not collapsed. The bridge result below depends on that convention.

## General multigraph bridge theorem, including degeneracies

Let B be any finite bipartite multigraph, with distinguishable edge instances and arbitrary isolated vertices. Form a graph with unordered pairs of left vertices and unordered pairs of right vertices, and a distinct edge for every ordered pair of B-edge instances. Then an edge is a bridge if and only if it is a diagonal e squared and at least one endpoint of e has base degree one, where degrees count multiplicity.

For distinct edge instances e,f, (e,f) and (f,e) are parallel, so neither is a bridge. For e=i--j, if i or j has degree one, the corresponding doubled vertex has degree one, proving e squared is a bridge. This includes an isolated base edge, whose square is an isolated single-edge component. Base degree-zero vertices contribute no incident edges and cause no exception.

Otherwise choose f different from e at i and g different from e at j. If either is parallel to e, its diagonal supplies another edge between the endpoints of e squared. If neither is parallel, write f=i--l and g=k--j with l different from j and k different from i. The three-edge path

  {i,i} --(e,f)-- {j,l} --(f,g)-- {i,k} --(e,g)-- {j,j}

avoids e squared and has distinct successive vertices. Hence e squared is not a bridge. This covers all parallel-edge cases, including f=g parallel to e. The independent program additionally tests all 81 bipartite 2-by-2 multigraphs with each of the four endpoint multiplicities in {0,1,2}, including the empty graph.

## Literal 4+6 repairs

For original word (689,689), T=34051591637519 has H(T)=585. The Hamiltonian path (1,2,3,4,5,6,7,8,9,0) has nonpath flip U=23708219209456. Its Hamiltonian path (1,2,3,4,7,5,6,8,9,0) has nonpath flip V=30137802689039. Converse V, then relabel old-to-new by (8,5,7,6,9,1,2,4,3,0), to obtain the canonical mixed target.

For original word (721,661), T=35013664310799 has H(T)=221. The same first path has nonpath flip U=22746146536176. Its Hamiltonian path (1,2,3,4,6,7,5,8,9,0) has nonpath flip V=30651051279887. Relabel old-to-new by (0,3,2,1,4,5,9,7,6,8), without converse, to obtain the same canonical mixed target.

The canonical mixed target has mask335 on vertices0..4, converse mask842 on vertices5..9, all second-to-first cross arcs except4->9. The boundary word (689,661) joins its merged target class to BB. Both U tournaments have strong-component canonical types (4,4),(6,81), in that order. Thus U differs from BB, whose SCC sizes are5,5. The target T and mixed V are strong and have distinct exact rooted addresses, while the two source nodes are nonstrong. Consequently all four vertices of BB--T--U--V--BB are distinct. Every path is explicitly verified Hamiltonian; every nonpath flip and every relabeling is exact, not just an invariant match.

The supplied path enumeration search counts (two alternative classes, return after14/13 neighbors) are not needed for these certificates and were not independently reproduced. No minimality is claimed.

## Remaining scope

The 13 still-unclassified target-boundary-leaf diagonal tiles are 657,659,665,697,729,913,921,977,985,987,989,1009,1017. They are neither assumed globally bridging nor assumed globally nonbridging. None is required to certify the selected324 family.

The four source-leaf bridges are genuinely global: every Hamiltonian path in their ordered-join source traverses the two size-five factors in sequence, and the full64 boundary includes every such spine tile. Thus their source has no omitted incident tile edges. This argument is separate from, and does not apply to, the target-only leaves.
