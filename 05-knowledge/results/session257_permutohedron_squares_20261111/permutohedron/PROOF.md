# A weighted permutohedral atlas for square-sum path surgery

**Status:** PROVED elementary constructions and refinement/surgery laws below; FINITE-EXACT Q25 atlas census. No literature-priority claim. No Chan–Pak, LRC, Monster, or general square-sum existence claim.

## 0. Inheritance and the actual new object

The closest inherited carrier is THM-207, `01-canon/theorems/THM-207-permutohedron-transitive.md`: adjacent order exchanges form a permutohedron. The newer `session210/coherent-permutohedron/PROOF.md` equips that carrier with a tournament-coherent orientation whose sinks are Hamiltonian paths; `session211/score-permutohedron/PROOF.md` identifies information lost by the score projection. Those are tournament constructions. Here the graph is undirected, the blocks are fixed **path components**, and existence is a seam predicate on orientation sheets over a **weighted** permutohedron. We do not import the tournament sink theorem into this different setting.

The closest graph assembly mechanism is `05-knowledge/results/session32_square_monster_20261023/ports/PORT_TRANSFER.md`: keep named ports and their pairing, not just the three proposed supervertices. Its corrected near miss is the branched A-block, which cannot be visited in one contiguous run. The hostile for the current proposal is stronger: even after correctly choosing actual path blocks, every feasible elementary geometric move can disappear. The finite experiment below finds precisely this at Q25.

The underused sidecar is the *orientation of each path piece*. The concept board is: path-cover components / weighted braid geometry / endpoint seam matrices / cut refinement / long signed interval reversal / owner birth colors. The anchor is exact surgery; the niche is functorial zonotope geometry; the wildcard is whether a geometrically adjacent state is a useful graph exchange. The answer to that last question is negative in the main finite example.

The owner's three edge colors are inherited birth-role classes, **not** Q12's three components, and not the square labels4,9,16. Their source is `05-knowledge/results/session40_marked_ports_square_sum_20261023/square-sum/README.md`, section2. These three distinctions are retained throughout.

## 1. Exact path-cover charts (PROVED)

Let G be a simple undirected graph on N>=2 vertices. Let F be a spanning linear forest contained in G, with k nonempty path components B_1,...,B_k. Fix one reference order for every component. A singleton has only one orientation; each larger component has two. An **oriented chart state** consists of an orientation epsilon_i for every component and an ordering pi of the k component labels. Write l_i,r_i for the resulting first and last vertices.

Concatenate the oriented components in pi order. This is an oriented Hamiltonian path of G containing F if and only if every seam

    {r_(pi_j), l_(pi_(j+1))} belongs to E(G), j=1,...,k-1.

Moreover, concatenation is a bijection between feasible states and oriented Hamiltonian paths containing F.

Proof: an existing path component of F must be a contiguous segment in any Hamiltonian superpath: its internal vertices already have degree2, and its two ends are its only possible entry/exit locations. A singleton is immediate. Contract each contiguous segment to obtain its unique order and orientation. Conversely the seam test makes the concatenation a spanning path. These maps are inverse. Global reversal identifies states in pairs for N>=2, giving unoriented paths. In a leaf-rooted chart, requiring a specified leaf first selects exactly one state per unoriented path.

For Q_n the seam test is exactly

    r_i+l_j is an integer square.

The graph source and port labels remain explicit. Merely knowing |B_i| does not determine this test.

### A positive exact polynomial compiler

Give each edge a commuting weight t_e. For fixed orientations put

    M_ij = t_{ {r_i,l_j} } if the seam is in G, and 0 otherwise.

Then the oriented path polynomial over this chart is

    (product_(e in F) t_e) * sum_epsilon sum_pi product_(j=1)^(k-1) M_(pi_j,pi_(j+1)).

Unoriented counting divides by2, or avoids division by fixing an obligatory leaf first. All coefficients before the orientation quotient are nonnegative integers. This is a subset/last-block dynamic program: D(S,j) is the polynomial for orders using S and ending in j; append an unused block i with multiplier M_ji. With variable orientations the state also retains the last block's orientation, and each unused block is appended in every distinct orientation. Initial endpoint restrictions are imposed at the first block.

This is a fixed-forest gluing formula, not a new general birth formula: Session40 already proves the full edge-marked vertex-arrival polynomial. Substituting a common variable per owner color or per square label is a specialization, and destroys individual edge incidences.

## 2. The actual zonotope: block-position sums (PROVED)

Put w_i=|B_i|>0. In an expanded concatenated word let p(v) be the position of vertex v, numbered1,...,N. Define

    y_i = sum_(v in B_i) p(v) - w_i(w_i+1)/2.

If block i follows blocks of total size L_i, then y_i=w_i L_i. In particular these coordinates forget internal orientation but remember the order of distinct blocks.

**Weighted carrier theorem.** The convex hull of these k! vectors is exactly

    Z(w) = sum_(i<j) w_i w_j [e_i,e_j].

Its normal fan is the braid fan, its dimension is k-1, and its graph is the permutohedron: two vertices are adjacent exactly when the corresponding block orders differ by an adjacent transposition. Swapping adjacent blocks a,b changes the vector by

    w_a w_b (e_a-e_b)

when a moves after b.

Proof: each unordered pair contributes w_iw_j to the later block's coordinate. Thus every y is the endpoint selection of the displayed segments prescribed by a total order. For any linear functional c, the zonotope support function is

    h_Z(c)=sum_(i<j) w_iw_j max(c_i,c_j).

Sorting blocks by increasing c realizes this value among the displayed y. Therefore the convex hulls agree. Since every weight is positive, the support-maximizer changes exactly at c_i=c_j; the resulting normal fan is the braid arrangement. Its full chambers correspond to total orders and chamber adjacency to adjacent swaps. The coordinate difference follows directly. k=1 is the one-point degenerate case.

For w=(3,4,5), the six vertices are

    (0,12,35), (0,32,15), (12,0,35),
    (27,0,20), (15,32,0), (27,20,0).

Their hexagon is an honest weighted permutohedron. This is not a Pythagorean consequence: the construction works for every list of positive block sizes. For Q12 its actual three components have these sizes, but there are no intercomponent seams, hence **none** of the six orders or orientation sheets is Hamiltonian. The geometry exists even when the feasible locus is empty.

### Orientation sheets and the exact induced exchange graph

For each fixed orientation vector, the feasible orders form a vertex subset of Z(w). The graph allowing only feasible adjacent block swaps is exactly the **induced** subgraph of its permutohedral graph on that subset. This assertion does not identify the edges of the convex hull of feasible vertices; that smaller polytope can acquire new edges. Nor is the feasible subset necessarily connected, convex, or nonempty.

All orientation sheets lie above the same geometric Z(w). Distinct sheets can have the same coordinate vertex and different feasibility. Example: in Q15 split the familiar path

    8,1 | 15,10,6,3,13,12,4,5,11,14,2,7,9.

Its seam is1+15=16. Reverse the first block without changing its place: the new seam8+15=23 fails. The projected position-sum vector has not changed. Named oriented endpoints are therefore an essential feasibility sidecar, even at order15.

## 3. Cut refinement is an exact zonotope projection (PROVED)

Cut some edges of F to obtain a finer path forest F'. Let fine-block weights be u_a, grouped into sets I_i whose union is coarse block B_i; thus w_i=sum_(a in I_i)u_a. Define the affine map

    (Lz)_i = sum_(a in I_i) z_a - sum_(a<b in I_i) u_a u_b.

**Refinement theorem.** L maps the entire finer zonotope Z(u) onto Z(w). When a fine-block order keeps every coarse group contiguous, its vertex maps to the corresponding coarse-order vertex. Arbitrary interleavings need not map to vertices.

Proof: each fine generator u_a u_b[e_a,e_b] with a,b in different groups maps to the parallel segment u_a u_b[e_i,e_j]. Summing all such generators between groups yields w_iw_j[e_i,e_j]. Each same-group generator maps to the single point u_a u_b e_i, canceled by the displayed translation. This proves equality of images. For a contiguous group, its fine coordinates contain both contributions from all earlier coarse groups and exactly sum_(a<b in I_i)u_au_b from internal relative order. Subtraction leaves w_i times the total earlier coarse weight, the coarse vertex formula.

For an explicit nonvertex image, split a weight2 block A into two singleton pieces and take another singleton B. The fine order A1,B,A2 projects to(1,1), the midpoint of the coarse segment with vertices(0,2),(2,0). The projection faithfully records interleaving rather than pretending it is a coarse order.

Geometric projection alone does not authorize a path gluing. Old feasible paths embed in the refined chart by retaining the removed edges as required seams, in their original incidence order and orientation. Conversely refined states may glue the new pieces differently and escape the old chart. Keeping the removed-edge set identifies precisely which refined paths belong to the old fiber. Thus cutting enlarges the admissible combinatorial search while preserving a canonical geometric quotient.

This is the missing compatibility between scales: a family of real zonotopes and actual graph fibers, rather than the word “permutohedron” attached to unrelated permutations.

## 4. Surgery distance chooses its own chart (PROVED)

For two unoriented Hamiltonian paths P,Q in the same labelled graph, put

    r=|E(P) minus E(Q)|, F=E(P) intersection E(Q).

Then F is a spanning linear forest with exactly r+1 components. Both paths are feasible states of its chart. Every r-edge replacement is therefore exactly a reorder/reorientation of r+1 common path pieces; conversely any two chart superpaths differ only in their seams, at most k-1 removed edges.

Proof: F is a subforest of the path P and has N-1-r edges, hence N-(N-1-r)=r+1 components. Apply section1. The converse follows because every component's internal edges remain common.

A Pósa endpoint rotation is especially simple in these coordinates. If P is split just after pivot v, with following vertex w and far endpoint x, replace edge v-w by v-x, provided v-x is in G. This reverses the suffix and reverses the orientations of **all** its constituent blocks. At a coarse two-piece chart it is one orientation change; in a fine chart it can reverse a long block interval. It is not generally one edge of the fixed fine permutohedron. The exact geometric displacement is

    sum_(a before b within the reversed interval) w_a w_b (e_a-e_b).

All block pairs inside the interval reverse their relative order; other pair contributions are unchanged. The orientation flips are additional sidecar changes invisible in this displacement.

Thus two notions of smallness differ: small *graph-edge surgery* and short *fixed-chart adjacent-swap distance*. Refinement can turn one graph-edge move into a long geometric move, while a small pair-adapted chart exposes it again.

## 5. Q25: the maximal universal forest gives a frozen fine atlas

The following census is FINITE-EXACT, not an all-order statement. An independent vertex DFS from the global leaf18 enumerates all ten unoriented Q25 paths. Their common edge set has14 edges. Its11 path components, in fixed reference orientations, are

    (1), (2,23), (3), (4,21,15,10), (5,20,16,9),
    (6,19,17,8), (7,18), (11,25,24), (12), (13), (14,22).

This forced-edge information is inherited in substance from Session40's ledger; the geometric chart and exchange experiment are the new use. It gives a10-dimensional weighted permutohedron with weights

    (1,2,1,4,4,4,2,3,1,1,2).

There are seven nonsingleton blocks and hence128 orientation sheets before fixing the starting leaf. A separate contracted-block search using only endpoint seam tests reconstructs exactly the same ten paths.

**New finite obstruction.** In the leaf18-first feasible atlas there are **zero** legal adjacent block swaps and **zero** legal individual block reversals. This remains true if global path reversal is allowed as a gauge: the root block has length2, so no one such move can send the first endpoint18 to the last endpoint while preserving a Hamiltonian state. Every elementary chart state is isolated.

Yet the true edge-surgery graphs are not discrete:

- one removed/added edge:8 edges, connected component sizes6 and4;
- at most two removed/added edges:19 edges, one component of size10;
- at most three removed/added edges:31 edges, one component of size10.

The section4 theorem explains how this can happen: the surgery exchanges reverse/reorder *groups* of the11 small blocks. Pair-adapted common-forest charts have only2 blocks for one-edge exchanges or3 blocks for two-edge exchanges. A fixed globally forced forest is canonical but not automatically a good local-move chart.

### Owner color specialization, retained as a control

Color1 is the eight edges born at15,16,17,23,25. Color2 is the fourteen edges born at13,14,19,20,21,22,24. Color3 is the remainder, including7-18. These are exactly Session40's owner birth classes, not a proper edge coloring and not forced/optional/unused status.

The unoriented Q25 path polynomial, after collapsing edge variables to x,y,z by these classes, is

    7 x^7 y^11 z^6 + 2 x^7 y^10 z^7 + x^6 y^12 z^6.

The polynomial is a specialization of the inherited edge-marked polynomial. It does not separate paths: seven paths share the first monomial. The certificate stores the per-path profile along with full paths, component orders, and orientations, so the lost information remains available. This is a concrete positive compiler output, not a new birth-law theorem.

For the square-label coloring, each seam carries q=r_i+l_j, and THM-4510's inherited reflection law R_q(v)=q-v reconstructs vertices along an oriented word. These arithmetic labels coexist with the birth colors; neither is interchangeable with the port pairing. The experiment independently reconstructs every path using its square-label reflection word.

## 6. Source / map / invariant / loss / sidecar

- **Source:** actual labelled graph G, a fixed spanning linear forest F, and Hamiltonian superpaths of F.
- **Target:** orientation sheets over Z(w), with a seam-feasible induced vertex subset in each sheet.
- **Map:** contract path components, retain their orientation/order, then take block position sums minus triangular offsets.
- **Preserved:** all vertices and F edges in the full chart; Hamiltonian feasibility iff the explicit seam predicate; adjacent component order exchanges in each fixed sheet; exact weighted polynomial; refinement projection.
- **Destroyed by geometry alone:** named endpoint labels, component orientations, edge-square and birth-role colors, adjacency restrictions, and internal graph alternatives outside F.
- **Restoring sidecar:** oriented endpoint dictionary, actual seam graph, F itself and removed-edge labels under refinement. For variable local forests use Session32's complete terminal-pairing signatures to choose the path pieces first.
- **Hostiles:** Q12 has a full geometric hexagon and no feasible path; Q15 has identical projection but one feasible and one infeasible orientation; Q25 has a completely frozen elementary atlas despite many one-/two-edge graph surgeries.

Session256's positive-state/observer lesson transfers here as a design rule: preserve the finite response needed by the next operation. The concrete response is an endpoint compatibility matrix and, under changing local covers, a terminal matching. It is not the tournament eight-coordinate response or cyclic Gram kernel. No scalar count-preserving source switch is inferred.

## 7. Open work and reproduction

**OPEN:** choose a small family of adaptive cut charts whose seam-feasible coarse moves provide efficient navigation without enumerating all Hamiltonian paths. Section4 gives a pairwise atlas after knowing both paths, not an algorithm for finding the next useful chart or a uniform connectivity theorem.

**OPEN:** characterize when a fixed forest's feasible induced permutohedral subgraphs are connected, or find useful sufficient conditions compatible with square-sum endpoint arithmetic. Q25 prohibits a blanket assertion even for the forest of all forced edges.

Run `python3 -B reproduce.py` and `python3 -O reproduce.py` in this directory. Both enumerate the same complete finite universes and regenerate `certificate.json`. The code uses explicit checks, not removable assertions. It independently compares vertex DFS with contracted-block seam search, checks all signed adjacent swaps and single-block flips on the ten feasible states, all45 path pairs for edge surgery, the owner-color polynomial, square-reflection decoding, zonotope support controls for3-4-5, and every contiguous fine order among120 refinements of weights(1,2,4,2,3). The unbounded claims rest on the proofs, not the finite support tests.

Read-only source inheritance: THM-207 and THM-4510 by full filenames; `graceful_square_exchange_20261019.md`; `monster_levels_square_sum_7_21_20261022.md`; Session32 ports; Session40 square-sum birth/marking note; Sessions210/211 geometry; Session256 SYNTHESIS. No theorem ID reserved and no Git mutation by this worker. The separate parent controls publication.
