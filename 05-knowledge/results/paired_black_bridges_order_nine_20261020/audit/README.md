# Independent audit of the order-nine bridge census

## Verdict

Status: this independent audit is complete. Final root review and publication are separate steps; neither approval nor successful publication by those steps is asserted here.

No mathematical defect or result-changing implementation defect was found in the audited claims. The exact order-nine census, the seven paired-black bridges, their leaf status, and the exhaustive short-cycle coverage all passed independent checks. The new score-based construction also proves the existence of an all-black embedded quotient quadrilateral at every order at least eight.

These are distinct claims: the all-order **existence** construction does not prove that **every** nonleaf edge has a short cycle, nor that every bridge is a leaf at all orders.

## 1. Full classification-map audit

`verify_map.cpp` does not call or copy the producer's incremental path-encoding routine. It constructs each representative's adjacency matrix, enumerates directed Hamiltonian paths, and re-encodes every completed path from its full ordered vertex list. It then explicitly forms the converse-merged fibre and compares it with the stored map.

It checked:

- All 268,435,456 map entries have valid node IDs.
- Every node's map count and least mask equal the recorded mass and representative.
- All 137,563,290 Hamiltonian paths of the 97,144 representative tournaments were re-enumerated.
- Each exact merged fibre has the recorded size, and every member maps to that node. Together with the map counts, this rules out missing members, spurious members, incorrect merges, and incorrect splits.
- The ordinary fibre size, Hamiltonian-path count, automorphism order, converse status, fixed-reflection count, strong connectivity, cyclic-triangle count, and degree formula agree with the producer.

Results: 97,144 merged classes, 191,536 ordinary classes, 2,752 self-converse classes, and 90,397 strongly connected merged classes. See `map_audit.json`.

As usual, passage from all fixed-Hamiltonian-path presentations to all tournament isomorphism classes uses the theorem that every tournament has a directed Hamiltonian path.

The order-seven and order-eight maps match the preceding session's independent Python census entry-for-entry: respectively 32,768 and 2,097,152 entries. All common node fields and graph summaries also agree. See `control_comparison.json`.

## 2. Independent full-graph bridge audit

`verify_graph.cpp` enumerates all masks and accepts exactly the minimum of the complement/reflection orbit, rather than using the producer's half-range filtering loop. It independently reconstructs all 67,092,480 paired-black multigraph edges.

It then runs an iterative Tarjan traversal directly on the full 62,645,527-pair nonloop support graph. This is independent of both the producer's online bridge DSU and its sparse two-forest bridge calculation. A simple-support bridge is retained only when its multigraph multiplicity is one. Loops do not enter the support traversal and can never become bridges.

The result is exactly the same seven bridges. Their unordered endpoint pairs are:

- (58, 96545), mask 64
- (2964, 96655), mask 4096
- (42880, 96554), mask 131072
- (43209, 96524), mask 2097216
- (48851, 97113), mask 16777280
- (50272, 97114), mask 16781312
- (48911, 97138), mask 67108928

Every bridge has a degree-one endpoint of mass two, blue count zero, and non-self-converse type. Both endpoints of every listed bridge are non-self-converse. The support has seven components. See `graph_audit.json`.

### Sparse certificate

`verify_certificates.py` independently checked every supplied certificate edge's canonical mask and map endpoints, verified that its two recorded edge sets are forests, and ran a separate iterative edge-identity Tarjan traversal on the certificate. It found 194,267 edges, comprising forests of sizes 97,137 and 97,130, and exactly the same seven bridges.

The sparse-certificate argument is valid for multigraphs. The first greedy forest spans each component of the original graph; the second spans the graph of edges not selected by the first. If a cut has at least two original edges, then either the first forest already contributes at least two, or the remaining-edge graph has an edge across that cut and its spanning forest contributes another. Thus cuts of size zero, one, or at least two are distinguished correctly. Loops are irrelevant, while parallel copies remain separate edge identities. The code's parent-edge check, rather than a parent-vertex check, handles parallel copies correctly.

## 3. Exhaustive short-cycle coverage

Each row of the 27,524-row quadrilateral file was independently checked for:

1. Four distinct vertices.
2. Presence of all four required support edges.
3. Multiplicity one and nonleaf status of the designated edge.
4. No duplicate designated edge across rows.

The direct graph audit independently classified every support pair. Any single nonleaf edge without a common neighbor had to match a previously verified quadrilateral row. Every row was used exactly once. The final edge partition is:

| Category | Number of multigraph edges |
|---|---:|
| Loops | 4,499 |
| Edges in parallel pairs or larger bundles | 8,461,023 |
| Remaining edges with triangle witnesses | 58,599,427 |
| Remaining edges with quadrilateral witnesses | 27,524 |
| Leaf edges | 7 |
| Uncovered | 0 |
| Total | 67,092,480 |

Consequently every nonleaf paired-black edge at order nine belongs to a cycle of length at most four, with loops and digons counted as lengths one and two. This is an exhaustive finite statement, not a universal-order theorem. Independent checks of the preceding atlas establish the same statement for every order from three through eight; see `lower_order_cycles.json`.

The producer's square search is sound: it searches a-b, a-x, x-y, y-b with x distinct from a and b, and excludes a and b from y. Since the support has no loops, y cannot equal x. Thus an accepted square cannot secretly reuse a vertex.

## 4. Blue graph

An independent enumeration of all reflection-fixed masks and edge-identity Tarjan traversal confirms 32,768 blue edges, 31,832 support pairs, 94,393 components including isolated vertices, and 24 bridges. All 24 bridges have a blue-degree-one endpoint. See `certificates_audit.json`.

This order-nine result does not by itself establish a universal blue leaf-bridge classification.

## 5. Commuting-square criterion and hostile examples

Let P and Q be arc-disjoint directed Hamiltonian paths of a tournament T. Reversing either path leaves the other directed, and the two reversals commute. They therefore give a four-step closed walk after passing to tournament isomorphism and converse classes.

For this to certify a particular paired-black quotient edge e as a nonbridge, it is sufficient that all four moves are black and that e occurs exactly once among the four quotient edge identities. Deleting e then leaves a walk connecting its endpoints; reducing that walk gives a cycle containing e of length at most four. Vertex identifications or repeated other edges do not invalidate that argument.

Arc-disjointness alone is not sufficient, and even arc-disjointness plus all-black status is not sufficient. Independent exhaustive permutation enumeration at order five found five disjoint-path witnesses across four nonloop, nonleaf target edges, and no valid commuting-square certificate for those targets. Two useful counterchecks, with zero-based vertex labels and P=(0,1,2,3,4), are:

- Mask 2, Q=(0,2,4,1,3): the four masks are (2,20,10,52). The last two moves are blue.
- Mask 14, Q=(1,3,0,4,2): all four masks (14,42,25,19) are black, but their quotient edge identities are (14,14,19,19). The target occurs twice, not once.

See `commuting_n5.json`. The producer's `ls.count(t)==1` test is therefore essential, not cosmetic.

## 6. New all-order score certificate

The proposed stronger sufficient criterion is valid.

Write d(v) for outdegree and c(T) for the number of cyclic triangles. Reversing a directed Hamiltonian path P changes only its endpoint scores: the start loses one and the end gains one. From

c(T) = binomial(n,3) - sum_v binomial(d(v),2),

the change is alpha=d(start P)-d(end P)-1.

Suppose P and Q are arc-disjoint and their four endpoints are distinct. Let beta be the analogous increment for Q. Reversal of Q does not change either endpoint score of P, and conversely. The four states consequently have cyclic-triangle counts

c, c+alpha, c+alpha+beta, c+beta.

They are pairwise distinct if alpha and beta are nonzero and alpha is neither beta nor -beta. Cyclic-triangle count is invariant under both tournament isomorphism and converse, so the four states remain distinct quotient vertices.

A blue presentation supplies a path-reversing anti-isomorphism. Hence the endpoint scores sum to n-1. Its reversal increment is therefore 2d(start)-n, of parity n. Requiring alpha and beta to have parity n-1 excludes blue moves in all four states. Thus these assumptions give an embedded all-black quotient quadrilateral.

### Four-vertex extension

Choose one arc A->B of P and one arc C->D of Q whose four distinct vertices avoid the four global endpoints. Introduce x0,x1,x2,x3 and replace the selected arcs by:

- P: A,x0,x1,x2,x3,B.
- Q: C,x2,x0,x3,x1,D.

The two internal paths prescribe all six pairs among the new vertices exactly once. Require both old global starts to dominate the block, and the block to dominate both old global ends. The four attachment requirements and these endpoint requirements are disjoint, so no orientation conflict occurs. Other old-new edges can be oriented arbitrarily.

Both path increments increase by four. With positive unequal initial increments, nonzeroness, inequality, and the required parity are preserved. The new arcs x1->x2 of P and x0->x3 of Q provide four distinct internal attachment vertices for the next extension.

The independently checked seeds are:

| Seed order | alpha | beta | Initial attachment (A,B,C,D) | Four cyclic-triangle counts |
|---:|---:|---:|---|---|
| 8 | 3 | 1 | (4,5,6,3) | (15,18,19,16) |
| 9 | 2 | 4 | (2,3,7,5) | (15,17,21,19) |
| 10 | 1 | 3 | (1,2,8,6) | (32,33,36,35) |
| 11 | 2 | 4 | (2,3,6,4) | (42,44,48,46) |

`verify_score_construction.py` independently decoded the seed masks, checked their scores and paths, and performed 25 extensions per seed, covering every order 8 through 111. All four moves were explicitly checked to be black and all four counts distinct at every step. The symbolic induction proves existence at every order n>=8; the finite tests are additional implementation checks. See `score_construction_audit.json`.

## 7. Reid reference and proof scope

The primary-source indexed abstract is K. B. Reid, *Equivalence of n-tournaments under k-path reversals*, abstract 691-05-29, Notices of the American Mathematical Society, January 1972, p. A-38. Theorem 1 announces that every tournament of order n>=2 can reach the transitive tournament by Hamiltonian-path reversals, with isomorphic copies allowed at each step. Theorem 2 announces the corresponding statement for any permitted path length.

Primary source: https://www.ams.org/journals/notices/197201/197201FullIssue.pdf

The indexed primary-source text was checked. Direct PDF retrieval returned HTTP 403. This audit did not inspect a full published proof of Reid's theorem and makes no such claim.

The announced Hamiltonian-path theorem implies connectivity of the combined metagraph, since a path reversal in a normalized presentation maps to the complement presentation after a converse identification. It does not by itself imply connectivity of the black-only graph or classify its bridges.

### Mixed-anchor consequence

The final manuscript's mixed-anchor theorem is valid, conditional only on the cited combined-connectivity theorem. An active Q component with no self-converse vertex has no incident blue edge, so it would be a whole connected component of the combined graph. The transitive tournament is a pure-blue Q isolate outside that active component, a contradiction. Any self-converse vertex in an active Q component is mixed, since a pure-blue vertex has Q degree zero. Every non-self-converse vertex is active, so the equivalent black-path formulation also follows. Contracting Q components preserves combined connectivity and leaves only blue edges between components. None of these arguments establishes that the active Q component is unique.

## 8. Final manuscript review

The final version of `BRIDGES_THROUGH_NINE_AND_LOCAL_CYCLE_LAWS.md` was reviewed. Its finite numerical claims agree with the independent audits. Its seven leaf representatives were additionally decoded into literal strong-component words; all are nonpalindromic words in singleton and directed-three-cycle blocks. There are exactly seven Q leaves, no mixed Q leaves, and six pure-blue Q isolates. See `leaf_family_audit.json`.

The score-separated square theorem, indefinitely repeatable four-vertex insertion, and mixed-anchor consequence are valid with the stated hypotheses. The manuscript correctly distinguishes:

- finite exhaustive edge coverage through order nine;
- an all-order sufficient local nonbridge certificate;
- all-order existence of an embedded black quadrilateral for every n>=8;
- a cited combined-connectivity theorem and its mixed-anchor consequence;
- the remaining unproved general bridge and short-cycle assertions.

The producer constructor was also run read-only and every one of its 101 reported orders was checked using the independent score/path/blackness verifier. See `producer_construction_audit.json`.

### Imported deficit-two result and exact pendant-bridge series

The later section-seven update imports, rather than reproves, the companion packet's independently approved theorem D(z)=z^4 and its complete leaf series

L(z)=(2z^4-z^5-z^7)/[(1-z-z^3)(1-z^2-z^6)].

Sources: sibling `session15-deficit-centers/EXACT_LEAF_SERIES.md`, with final approval status supplied by `PUBLIC_README.md` and `AUDIT_COMPLETION.md`. The companion `PUBLICATION_MANIFEST.json` hash was verified as `0123f01652018d17f62b6bf25f3acfa71fa743e55e798edb4daf1f7dfaf9a55c`. This audit did not repeat the companion's all-order center proof or residual census.

The new bridge consequence is valid. An NS leaf has a nonpalindromic rigid strong-component word, so it has more than one component. A mixed leaf of order greater than four has a nonempty rigid wrapper around its S4 center, so it too is nonstrong. For any nonstrong tournament of order at least three, a directed Hamiltonian path starts in the first strong component and ends in the last; its first vertex therefore beats its last. Their edge is not a consecutive path edge. Complementing the free arcs reverses that edge and preserves the directed path, producing a directed Hamiltonian cycle. Thus every complement partner is strong.

At orders n>=5 a leaf's partner therefore cannot be another leaf. At n=4 exactly one bridge has two leaf endpoints. The leaf-to-pendant-edge count consequently has exactly one double count, giving the complete pendant-bridge series L(z)-z^4. If E(z) counts bridges with no leaf endpoint, the exact identity is Bridge_Q(z)=L(z)-z^4+E(z), with nonnegative integer coefficients and E(z)=O(z^10). The latter vanishing range is finite census evidence; E(z)=0 is not established.

`verify_pendant_consequence.py` checks the imported rational leaf series against 201 companion coefficients and confirms the pendant correction against every exact bridge count from n=3 through n=9. See `pendant_consequence_audit.json`. This checks the consequence and its algebra, not the imported center theorem again.

A final producer-source cleanup replaced the flattened two-dimensional-array initialization with per-row initialization and made minor formatting and signed-size fixes. These remove portability concerns without changing the mathematical enumeration; the audit's full-map verification remains applicable.

## 9. Reproduction and limitations

Run the audit from the session15 stage directory using `bash audit/run_audit.sh`. The full stage needs `census_n9.json`, the generated `census_n9.map`, both supplied certificate files, the seed files, the sibling session14 atlas for control checks, and the companion session15 deficit-center packet for the imported leaf-series check. If using a compact publication, first decompress the node atlas and sparse certificate and regenerate the map with the producer census program; the 1 GB map is not intended as a public payload. `expected_nodes.tsv` is regenerated from the node JSON by the audit driver. C++ assertions must remain enabled. The checks use the supplied little-endian 32-bit map and independently recompute its mathematical interpretation; they are not a machine-checked proof of the C++ compiler or runtime. Both producer and audit currently assume the standard 32-bit `int`/map layout on this machine. All graph cardinalities fit the declared integer types at n=9.

No producer files, repository checkout, or remote destination were modified. All audit output is confined to this `audit/` directory. The only source-access failure was direct retrieval of the AMS PDF; the indexed primary-source announcement was sufficient to verify the cited statement, but not its proof.
