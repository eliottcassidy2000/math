# Complete balanced order-ten native-C blue wiring

24 October 2026. **FINITE-EXACT**, relative to the explicitly frozen exhaustive
almost-regular order-ten atlas. This packet does not enumerate the full black
carrier, does not determine the full metagraph, and makes no all-order claim.

## Statement

Let B be the set of self-converse ordinary isomorphism classes of tournaments
on ten vertices with score multiset 4^5 5^5. Use actual fixed-Hamiltonian-path
presentation words, not arbitrary representatives of fibres. A blue word is
fixed by reversal of the path coordinates, with the arc direction reversed.
Join two blue words by native C, which reverses all nonpath arcs. On class nodes
retain every internal C-pair, with multiplicities and loops. Then:

- B has 229 classes, all strong, drawn from the 13,333-class frozen balanced atlas.
- There are 41,699 blue ports: 22,904 internal and 18,795 boundary ports.
- The induced multigraph has 11,452 edges: 11,428 nonloops and 24 loops.
- It has exactly one connected component, no bridges, and no articulation
  vertices. Its cycle rank is 11,452 - 229 + 1 = 11,224.
- The 18,795 boundary edges hit exactly 786 distinct ordinary classes, all
  self-converse with score profile 3,4,4,4,4,5,5,5,5,6 and refinement error zero.
- Exactly 20 source nodes have positive refinement error, all error 40; the
  remaining 209 have error zero. The induced 20-node core is already connected
  and has 90 edges. Every one of the 209 zero-error nodes has a direct neighbor
  in this core. Thus the core is a connected dominating set.
- The internal edge partition is 90 core-core, 1,857 core-zero, and 9,505
  zero-zero. This is not a merger of formerly disconnected balanced core
  components: four isolated nodes in the older 24-node positive-error graph
  are nonbalanced, with scores 3^5 6^5, and are excluded here.

Ordinary and converse-merged node labels agree on these self-converse classes.
Neither assertion counts black edges. Boundary classes are only the targets
of this layer, not all classes with the displayed score profile; no complete
boundary-to-boundary graph is claimed.

## 1. Frozen classification and exact canonicalization

SOURCE_PINS.json identifies every inherited file by SHA256, source packet and
artifact. The atlas freeze is manifest
0462baf160610115c8b2436c2c57cce1964d0439ce49e86b9b43ce00f2bda3f2.
The scalar-control freeze is
7431d29c032d92aa6da1807c436bb839ed9c84909872de31d7f64b09d4963d43.
The checker validates every input tournament and its balanced scores. It
reconstructs all automorphisms and the entire atlas converse involution, and
checks total labelled mass 48,251,508,480 against the frozen receipt.
Completeness of the starting 13,333-class atlas is inherited; this packet
independently verifies its rows and uses a second complete blue-word carrier
below rather than claiming a fresh unrestricted tournament census.

Start color refinement with outdegree, repeatedly retaining the old color and
outneighbor counts into every color. Stable cells encountered here have sizes
1, 3 or 5. The checker verifies internal regularity. On a size-three cell its
unique outgoing neighbor is the successor; on a size-five cell the unique
outneighbor with one common internal outneighbor is the successor. It checks
that this is a single cycle and that the internal arcs are exactly the cyclic
regular tournament. Consequently an isomorphism must preserve the ordered
stable cells and act as a rotation on each intrinsic cycle. Trying every
rotation product (at most 25, verified) and comparing the full adjacency matrix
is an exact canonical form, including all cross-cell correlations. All
minimizers give the full automorphism group. Apply the same construction to
the converse to recover all anti-automorphisms and test the anti-involutions.
Nontrivial relabeling controls cover every output source and boundary class.

The 229 source automorphism orders are 1:220, 3:4, 5:4, 9:1. Actual stable-cell
patterns are 1^10:200, 1^4 3^2:5, and 5^2:24. Refinement error is
sum over automorphism vertex orbits O of (size of the stable cell of O - |O|),
not a blue degree or Hamiltonian-path count. The 786 boundary classes have
orders 1:778, 3:6, 9:2 and all have error zero, checked directly from their
complete automorphism orbits.

## 2. Complete blue carrier, without full black enumeration

Index the 36 nonpath pairs (i,j), j >= i+2, in lexicographic order. A word sets
i -> j when its bit is one; all consecutive arcs i -> i+1 are fixed. For a
Hamiltonian path P of T, encode its full normalized matrix this way. Two paths
of the same T give the same word exactly when an automorphism takes one path
to the other. The automorphism action on complete ordered paths is free.

A blue path induces the anti-involution q(P_i)=P_(9-i). Conversely, for each
fixed-point-free anti-involution choose one member from each of its five
pairs, order the chosen half, and append its reversed q-image. Retain exactly
Hamiltonian paths. This exhausts blue paths, with no duplicate path across
anti-involutions because the path determines q. The 243 anti-involutions need
243 * 2^5 * 5! = 933,120 candidate halfpaths. Every resulting word is checked
against its literal path witness and occurs exactly |Aut(T)| times within
its class. Different ordinary classes have disjoint word sets.

Independently the reflection (i,j) -> (9-j,9-i) partitions the 36 bit positions
into 20 orbits. Enumerating its entire 2^20 fixed-word cube in Gray order and
filtering scores 4/5 produces exactly the same 41,699 words, as sets. This is
an independent completeness control on the entire blue balanced carrier; no
2^36 unrestricted black enumeration is used. An independently frozen scalar
half-path polynomial is recomputed and agrees with all three port totals.

## 3. Native C, endpoint transport, and boundary identification

C(w) = (2^36-1) XOR w. It is a free involution on words and commutes with the
reflection, so it preserves blueness. Literal reversal of all nonpath arcs
is separately checked. If d is the source score and d' the target score, then

    d'(i) = 9 - d(i) + [i=0] - [i=9].

The correction comes from leaving the mandatory path edges unchanged instead
of reversing them. A balanced blue word has complementary endpoint degrees.
Its C-target is balanced precisely for the high-low endpoint case (5,4).
The low-high case (4,5) produces exactly the stated boundary score profile.
Membership of C(w) in the full source word set is tested directly and agrees
with this criterion. Each internal unordered C-pair is counted once by w<C(w),
even when both words belong to the same class. Thus loops have two ports.
Boundary words are canonicalized by the exact procedure above; their full
canonical matrices, canonicalizing orders and automorphism data identify all
786 distinct targets. Boundary C-pairs have exactly one balanced endpoint,
so boundary ports and boundary edges agree in count.

Every node satisfies blue_ports = internal_degree + boundary_edges, with loops
counted twice in internal degree. The global identity is
41,699 = 2 * 11,452 + 18,795. There are 9,243 occupied unordered internal pairs
including loops; multiplicities have histogram 1:7320, 2:1654, 3:252, 4:17.

## 4. Topology and inherited controls

Exact word ownership embeds the older 20 balanced positive-error nodes and all
90 induced edges, with every old literal edge checked. Four nonbalanced old
classes are explicitly excluded. The induced core is connected; multi-source
breadth-first search from that core finds all other nodes at distance one.
The explicit spanning tree consists of 19 core edges and 209 attachments;
each edge carries an actual native-C word witness in topology.json.

Connected components and low-link bridge/articulation tests operate on the
actual recovered adjacency. A multigraph bridge can only occupy a nonloop
pair of multiplicity one. The implementation checks low-link results against
exhaustive vertex deletion for every node and exhaustive eligible edge-pair
deletion. The baseline component count is explicit, so disconnected graphs
cannot incorrectly pass a test that assumes connectedness. Loops and parallel
edges are handled separately and tested in adversarial small fixtures.

## 5. Reproducibility and limits

normal.json is the full finite certificate: 229 source records, 786 boundary
records, 41,699 literal port records, exact adjacency multiplicities, component
partition, histograms, and control counts. The independent record reader
verify_certificate.py checks source paths, blue reflection, C-pair witnesses,
boundary canonicalizing orders, incidence conservation and topology. It rejects
four deliberately corrupted certificates and passes five topology fixtures.
All checks use explicit exceptions, not removable Python assertions. Normal
and optimized (-O) full regeneration must agree byte-for-byte. The independent
auditor's work is separate from this producer packet.

The results establish this finite blue induced graph and its boundary only.
They do not establish black connectivity, an all-size wiring theorem, or a
full classification of order-ten tournaments.

An independent review packet, session52-audit-balanced-blue-wiring, additionally
constructs topology-certificate.json and verify_topology.py: a spanning
2-connected subgraph consisting of a triangle and 226 two-edge ears, using
455 edges. This is an independent compact proof route for the reported
connectivity, absence of articulations and absence of bridges. That review
artifact is not required to run the self-contained producer above.
