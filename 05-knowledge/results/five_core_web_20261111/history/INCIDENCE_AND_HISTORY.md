# From old boundary incidence to the golden thirteen

**Status: PROVED IN NOTE** for the elementary incidence and book identities;
**FINITE-EXACT** for the displayed small controls; **OPEN** for any all-target
coverage implication. This is a research contribution, not an independently
audited canon promotion. Chan–Pak remains open.

The most productive old source was not another Chan–Pak attempt. The July20
spectral collision in [THM-1780, H-leaves-the-spectral-ladder-at-n6](../../../../01-canon/theorems/THM-1780-H-leaves-the-spectral-ladder-at-n6.md)
is also a collision of the **entire golden triangle-overlap graph**. Its
five-vertex H13 cell is the later physical unit-switch source in
[THM-5117, physical-cut-switches-reversal-cancellation-and-primitive-sparse-cubes](../../../../01-canon/theorems/THM-5117-physical-cut-switches-reversal-cancellation-and-primitive-sparse-cubes.md).
Thus the old obstruction already contains an operation that the new construction
needs. The new contribution is this incidence interpretation and its connection
to the two infinite continuations of the H9 cores, not discovery of H13 or69/71.

## 1. The intrinsic incidence theorem

Let R(T) have a vertex for every cyclic triangle of T; two vertices are adjacent
when the triangles share a directed arc. Let B(T) be the bipartite incidence graph
between those triangles and all tournament arcs belonging to at least one of them.

**Proposition.** Every clique of R(T) with at least two vertices consists of
triangles containing one common directed arc. Every maximal clique with at least
two vertices is exactly the set of cyclic triangles containing one such arc.
Consequently R(T) is diamond-free, and each of its vertices belongs to at most
three maximal cliques of size at least two.

**Proof.** Two adjacent triangles can be written abc and abd, with common arc
a→b. Their cyclicity gives b→c,b→d,c→a,d→a. A third triangle meeting both
in an arc but not containing ab would have vertex set acd or bcd. The first has
a sink a and the second a source b, so neither is cyclic. Thus a third common
neighbour must contain ab. All later vertices of a clique obey the same argument.
The two triangles have only one common arc, so two distinct maximal cliques
cannot share two vertices. Each triangle has only three arcs. QED.

Equivalently R(T) is the intersection graph of the **linear, 3-uniform,
Helly hypergraph** whose ground elements are arcs and whose hyperedges are
cyclic triangles. It need not be an ordinary line graph: one triangle can have
three independent neighbours, one on each of its three arcs.

There is an exact reconstruction of the unlabelled B(T) from R(T). Make one
arc-node for each maximal clique of size at least two and join it to its members.
At a triangle-node v lying in r such cliques, attach 3−r private arc-leaves.
This restores B(T), including the multiplicity of private arcs. What it does
**not** restore is the identification of arc endpoints as tournament vertices.

For S5a, the shared-arc incidence is a path through three triangle-nodes.
For S5b, all three triangle-nodes meet one arc-node, a star. Both B graphs are
trees. The affine A2 label of the latter's overlap K3 is a valid adjacency/Cartan
observation; it does not mean that its common-arc nerve has a topological hole.

This is an exact continuation of the old
[endpoint-collisions-as-incidence-not-adjacency-s95](../../../../07-reflections/endpoint-collisions-as-incidence-not-adjacency-s95.md)
move: three coincident preimages are a common boundary incidence, not necessarily
an independent triangular circuit.

## 2. The same golden P4 carries thirteen and seventeen

Let Q have rows `(18,8,3,5,14)`. Its triangles are

    013, 024, 034, 123.

Their arc-overlap graph is P4, in the order `123–013–034–024`. Hence its
Perron eigenvalue is φ and a Perron vector is `(1,φ,φ,1)`. The six-vertex
chain T6 from Session252 has triangles `012,123,234,345` and the same P4.
But their odd-cycle packing polynomials are

    Z_Q(x)=1+6x,
    Z_T6(x)=1+6x+x².

They both have four cyclic triangles and two directed pentagons. The difference
is one pair of **vertex-disjoint** triangles in T6; every pair in Q intersects.
At x=2 this gives H(Q)=13 and H(T6)=17. The full vertex-intersection matrices
of the triangle families expose the lost coordinate: T6 has a zero off-diagonal
entry, while Q has none. Arc-overlap records only entries equal to2.

Adding a source to Q gives `(18,8,3,5,14,31)` and preserves H. The resulting
order-six tournament and the isomorphic copy of T6 with rows
`(38,8,2,5,15,30)` are **exactly the July20 THM-1780 pair**. They also have equal
ordinary adjacency characteristic polynomials. Thus retaining both adjacency
spectrum and the entire triangle arc-overlap graph still loses H.

This is not merely another negative check. The same Q admits the two actual
cut switches from THM-5117. Append a sink, then an apex beating precisely the
specified original vertices:

| Cut | Resulting rows | H |
|---|---|---:|
| `{0,1}` | `(50,40,99,101,110,64,3)` | 71 |
| `{1,2}` | `(114,40,35,101,110,64,6)` | 69 |

Q itself is in the small-parent construction: deleting original vertex1,2 or4
leaves H5, and deleting0 or3 leaves H3. A concrete next experiment is therefore
to **flatten this two-stage small-parent construction into paired five-core
addresses while retaining the switch's marked response**, rather than adding
Q as an unexplained sixth atom. Equal terminal H alone would not flatten the
needed interface.

The corrected [THM-4803, intransitivity-is-the-relation-modulo-its-modular-tree-and-a000568-from-the-prime-cycle-index-census](../../../../01-canon/theorems/THM-4803-intransitivity-is-the-relation-modulo-its-modular-tree-and-a000568-from-the-prime-cycle-index-census.md)
supplies an inherited positive boundary: the *labelled cyclic-triple hypergraph*
determines the whole path-cover profile. Its published realization theorem and
the corrected in-repo path-cover argument are inherited here, not re-audited.
So there is a precise ladder: arc-overlap → arc incidence → vertex-labelled
triangle incidence. The last carrier retains H; the first two need not.

## 3. The two H9 branches continue into tribonacci and Fermat families

Let T_n be the chain tournament with forward consecutive arcs and all other
arcs backward, as in Session252. Let B_n=C3[1,1,n−2], with the large block
transitive. B5 is S5b up to relabelling; T5 is S5a.

The cyclic triangles of B_n all contain the common arc between its singleton
blocks, so R(B_n)=K_(n−2). Every directed cycle contains those two vertices and
a nonempty subset of the transitive block, traversed in its unique order.
It is odd exactly when that subset has odd cardinality. There are2^(n−3)
such subsets and no pair of disjoint cycles. Therefore

    Z_Bn(x)=1+2^(n−3)x,     H(B_n)=1+2^(n−2), n≥3.

This is the familiar reversed-extreme-arc/Fermat family, now connected directly
to the common-arc continuation of the H9 bank. At n=2^j+2 its H is the Fermat
number F_j, regardless of whether F_j is prime.

The inherited interval-packing recurrence gives T_n's tribonacci H. The branches
agree longer than their overlap geometry does:

| n | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---:|---:|---:|---:|---:|---:|
| H(T_n) | 3 | 5 | 9 | 17 | 31 | 57 |
| H(B_n) | 3 | 5 | 9 | 17 | 33 | 65 |

At six vertices the coincidence is a genuine fugacity-2 cancellation:

    Z_T6(x)−Z_B6(x)=x(x−2).

The chain's one disjoint triangle pair compensates for the book's two extra
odd cycles. At five vertices both full polynomials are `1+4x`. Thus the
two equal-H9 cores are indistinguishable by every fugacity, while their first
extensions become distinguishable by every nonzero fugacity other than2.
This is a direct operational use of the old
[fugacity-axis-five-points](../../../../07-reflections/fugacity-axis-five-points.md)
idea, beyond matching its five evaluation labels.

For all n≥7, H(B_n)>H(T_n). Indeed D_n=H(B_n)−H(T_n) obeys

    D_n=D_(n−1)+D_(n−2)+D_(n−3)+2^(n−5)−2,

with D3=D4=D5=D6=0. The forcing is positive from7 onward. This is a proved
comparison of these two named continuations, not a tournament extremal theorem.

More generally the transitive-block blowup C3[a,b,c] has cyclic triangles
indexed by `(i,j,k)` and

    R(C3[a,b,c]) = K_a □ K_b □ K_c.

Sharing an arc means exactly changing one coordinate. Hence the overlap carrier
includes Hamming graphs, grids and cubes: C3[1,2,2] gives C4 and H11;
C3[2,2,2] gives Q3 and H45. The common-arc star, cyclic square, and cube are
actual tournament carriers, not artificially oriented comparison graphs.

## 4. Distant mechanisms recovered and what they now suggest

The search went beyond the existing Session252 atlas. It read substantial
March, May, June, July and August material, then followed repaired descendants
when a useful object appeared. These are idea provenance routes, not blanket
endorsements of every historical assertion.

| Old source/mechanism | Exact transfer or actionable question |
|---|---|
| March16 [base42-erdos-straus](../../../../07-reflections/base42-erdos-straus.md): factor an equation to expose divisors | On restricted flank words where g is affine, complete the compiler's bilinear rectangle and search divisor constraints. A residue cover needs a positive bounded lift; the old general claims are not used. |
| March16 [the-ladder-counts-itself](../../../../07-reflections/the-ladder-counts-itself.md), [the-fractal-in-the-spectrum](../../../../07-reflections/the-fractal-in-the-spectrum.md): carries as an observable | Derive exact address increments before fitting scale recurrences. The current g recurrence carries the required binary phase explicitly. |
| May05 [lucas-summand-graph-zeckendorf-geometry](../../../../07-reflections/lucas-summand-graph-zeckendorf-geometry.md): local conflict language | Keep its no-adjacent-support automaton as inspiration; do not assume an arbitrary conflict graph is realizable as a tournament's odd-cycle graph. The H13/T6 example identifies the missing vertex incidence. |
| May17 [fugacity-axis-five-points](../../../../07-reflections/fugacity-axis-five-points.md) | The exact x(x−2) collision above turns the old evaluation comparison into a new two-family identity. |
| June06 [polarized-delta-fields-band-gaps-and-numerical-semigroups-s699](../../../../07-reflections/polarized-delta-fields-band-gaps-and-numerical-semigroups-s699.md) | Treat its cofinite-spectrum wording as the desired conclusion, not an input. The useful object is the gradient of H on *physical* moves; the69/71 switch supplies one literal unit gradient. |
| July20 [THM-1780](../../../../01-canon/theorems/THM-1780-H-leaves-the-spectral-ladder-at-n6.md) | Its old cospectral obstruction is the new golden-incidence obstruction and the later unit-switch source. |
| July20 [shear-duality-sums-count-elements-products-count-orientations-opus-S420](../../../../07-reflections/shear-duality-sums-count-elements-products-count-orientations-opus-S420.md) | Track two scales: source-word counts and target values. Their growth functions need not imply the same coverage. |
| [farey-fibonacci-additive-basis-carrier-s169](../../../../07-reflections/farey-fibonacci-additive-basis-carrier-s169.md) | Distinguish four representation economies: smoothing, bounded arity, a normal form, and an address. The five-core claim needs actual acceptance/normalization, not only many representations. |
| [pollock-tetrahedral-defect-pair-descent](../../../../07-reflections/pollock-tetrahedral-defect-pair-descent.md) | A five-summand cover becomes absence of **pairs of four-sum defects at a prescribed displacement**. For five-core charts, examine correlated failures under a physical address switch, rather than proving any one chart covers a tail. |
| August21 [AMM dyadic Fibonacci near miss](../../../../07-reflections/amm-dyadic-fibonacci-near-miss-and-dual-front-two-coordinate-reframe-kps-S188.md) and [Sturmian phase repair](../../../../07-reflections/amm-dyadic-sturmian-phase-cocycle-repair-kps-S189.md) | A two-operator word requires phase plus boundary state. Any decreasing inverse must retain both flank parity and the actual source-interface state. |
| August16 [frozen five ancestry gate](../../../../07-reflections/lrc_endpoint_ufull_frozen_five_common_ancestry_gate_20260816.md) | Five scalar marginals can all be correct while a common realization is missing. For the two H9 cores, transfer the marked joint response before identifying equivalent numbers. |

The Pollock transfer is especially compatible with new address surgery.
Let S_i be the exact image of core/cut chart i. A physical rewrite with target
displacement δ supplies correlated queries `N∈S_i` and `N+δ∈S_j`. Instead of
trying to make S_i cofinite, seek a no-simultaneous-defect theorem for the
joint pair. A useful hypothesis must specify the actual rewrite, its guards,
and its well-founded descent; δ=2 at69/71 is a control, not a uniform rule.

## 5. Connection contract and next decisive tests

| Source → target | Preserved | Lost / required sidecar | Next test |
|---|---|---|---|
| Tournament → triangle arc-overlap | Common-arc pair incidence | Vertex identifications and longer cycles | H13/T6 already separates H |
| Arc-overlap → reconstructed B(T) | Entire abstract triangle–arc incidence | Arc endpoints | Keep triangle vertex-intersection matrix; test whether it suffices in larger examples |
| C3[a,b,c] → Hamming graph | Three block-coordinate changes | Internal ordering/cut contact | Lift a core cut into the three factor coordinates |
| S5a/S5b → chain/book continuations | Explicit birth family and exact H | Equality under other future operations | Compare full seven-coordinate responses at H17 |
| July H13 cell →69/71 switch | Actual two-edge apex exchange | Full target-placement control | Flatten both switched descendants into compatible five-core addresses |
| Pollock paired defects → chart failure correlation | A real displacement under a guarded rewrite | The rewrite's descent/coverage obligation | Test homometric address pairs for simultaneous inverse failures |

Reproduction: `python3 05-knowledge/results/five_core_web_20261111/history/incidence_probe.py`.
It uses NetworkX for finite graph isomorphisms; HP is counted independently by
subset DP and checked against literal odd-cycle packing for the named families.
The clique control universe is all1,099 labelled tournaments of orders1–5.
The collision scout enumerates all those, then order-six codes0–83 and stops at
the first unequal-H collision. It is not an order-six class census. The abstract
proofs above do not depend on the scout.
