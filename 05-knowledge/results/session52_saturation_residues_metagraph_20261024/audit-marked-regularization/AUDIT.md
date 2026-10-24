# Independent audit: marked Type-I regularization

24 October 2026. **PASS** for the candidate's MR1–MR9 and elementary count corollary, with the scopes below.
The audited candidate is the complete 13-file manifest-defined copy in
`producer/`, whose manifest SHA-256 is
`30ca3999e8c4abffba1b0641d87bdc25840ff0d5c1fd3fa94cadbb4385e2966e`.
No publication was performed by this auditor.

## Contribution and source boundary

**CITED / KNOWN:** the map from T to the regular tournament D(T) of order
2n+1 is the construction in Lemma 4.3 of Hatem, Mohamed, Müller and Siniora,
[*Cofinality of Regular Tournaments*, arXiv:2605.04846v1](https://arxiv.org/html/2605.04846v1#S4).
The auditor separately opened the primary source and checked its edge recipe.
It has two isomorphic, same-oriented copies, reversed off-diagonal cross
outcomes, matching arcs U_i→V_i, and V_i→z→U_i. The source's b/w name change
for the added vertex does not affect the recipe. The HTML record's displayed
version and manuscript dates differ; no chronology inference is needed.

The packet correctly distinguishes repository THM-4557's skew-Hadamard
doubling: that is a different construction, with a converse second copy,
a different cross rule, and a doubly regular input hypothesis. Neither its
construction nor its automorphism proof is substituted for the current map.
The present structural extractions are proved, but no priority claim for
those statements is established by this audit. The source pins are not a
cryptographic snapshot of the external arXiv website.

## Independent proof review

1. **Regularity.** A vertex U_i has d_i internal out-arcs and n−d_i cross
   out-arcs. V_i has n−1−d_i cross out-arcs, d_i internal out-arcs, and one
   to z. The root has n out-arcs. The displayed matrix is consistent with
   the defining arc conventions and every score equals n.

2. **Rooted maps.** Fixing z fixes its out-copy U and in-copy V setwise.
   The restriction to U is a source isomorphism f. The U-out-neighbor
   profile of V_i is the open in-neighborhood of source i. These profiles
   are pairwise different: if i→j, coordinate i belongs to the latter
   profile, but not the former. Therefore the map on V must also be f.
   Conversely every diagonal extension preserves all arcs. This establishes
   the full rooted isomorphism claim, not merely equality of group orders.
   Swapping the two copies gives the claimed converse intertwiner.

3. **Actual arc transport.** A source pair controls the four cross/same-copy
   edges on its two matched pairs. If i→j, their orientation is the directed
   square U_i→U_j→V_i→V_j→U_i. Distinct source pairs control disjoint edges.
   Any subset of source reversals therefore becomes a simultaneous union
   of edge-disjoint square reversals. The remaining 3n edges are fixed
   root/matching triangles. The labelled Hamming metric is multiplied by
   four; no unmarked edit-metric statement follows.

4. **Paths and cycles.** Read a source path backwards, replacing each vertex
   by its U,V matched pair and prefixing z. Root and matching arcs are
   automatic; each arc between pairs is equivalent to the corresponding
   forward source path arc. This proves an exact equivariant bijection with
   the distinguished cycles, not with all target Hamiltonian cycles.
   An ordered spanning path has trivial automorphism stabilizer.

5. **Native complement boundary.** None of the nonpath source blocks meets
   the lifted cycle. Their union has 2(n−1)(n−2) edges; full target noncycle
   complement reverses (2n+1)(n−1). The former support is contained in the
   latter, hence the two operations differ on exactly 5(n−1) labelled
   target edges. They agree only at n=1. Keeping the target regular is not
   sufficient to identify these different operations or their port data.

6. **Canonical root.** Each U_i,V_i pair has zero common out-neighbors.
   The root has d_i common out-neighbors with U_i and n−1−d_i with V_i.
   Therefore z is the unique all-positive-common-out vertex exactly when
   the source has no source and no sink; otherwise that predicate selects
   no vertex. This proves full automorphism preservation and unmarked
   ordinary/converse class-faithfulness on that subdomain. It does not
   characterize every source whose full automorphism group is preserved.

7. **Recentring.** For T=1+Q and S=Q+1, the candidate's map sends old U_s
   to new z, old V_s to new U_t, old z to new V_t, and fixes both Q copies.
   The old U_s neighborhoods are Q+1 on both sides. The Q blocks, root
   arcs, diagonal arcs and four source/sink-to-Q edge types each match the
   new recipe. Thus the isomorphism is all-size, not an inference from
   matching numerical invariants.

8. **Quotient obstruction.** Q=C3 gives the ordinary collision at source
   order four. Those two sources are converse, so it alone says nothing
   negative about merged faithfulness. Q=1+C3 gives sources with ordered
   SCC-size sequences (1,1,3) and (1,3,1); neither that sequence nor its
   reversal identifies the two. Their doubles of order eleven are
   explicitly isomorphic. Q=C3+S4 supplies the second order-17 witness.
   Generic exhaustive class comparisons below establish that four and
   five really are the respective first collision orders for this map.

9. **Elementary count corollary.** For n≥2, ordinary classes with a source
   and with a sink each number A_(n−1), with intersection A_(n−2), proving
   N_n=A_n−2A_(n−1)+A_(n−2). A self-converse tournament having either extreme
   has both; removing the unique pair leaves a self-converse tournament.
   Conversely adding a source and a sink to such a tournament reverses this
   bijection, so J_n=S_n−S_(n−2), including n=2 with A_0=S_0=1. Converse
   orbits of the nonextreme classes number (N_n+J_n)/2. MR6 injects ordinary,
   self-converse and converse-merged classes into the corresponding regular
   classes at order 2n+1. These are lower bounds requiring source counts as
   inputs, not an A000568 recurrence. Independent enumerations through n=5
   confirm all three quantities and all formula boundaries.

## Independent finite computation

`audit.py` imports no candidate code. It represents tournaments as bit rows
and constructs the lift by XORing copy parity on each source pair. Its exact
isomorphism solver propagates bitset candidate domains, using only degrees
and neighbor-subtournament score multisets for initial pruning. Every
complete map is verified on every ordered pair. On every source through
order five, its source automorphism results are additionally checked against
all permutations. Paths are generated by recursive directed extension,
not by the candidate's permutation filtering or subset dynamic program.

The complete labelled universe is 1,099 tournaments, orders one through
five. Both normal and optimized runs independently check:

- every full target automorphism group and every rooted stabilizer;
- canonical-root conditions and full group preservation on nonextreme inputs;
- all 10,650 single-source-edge block reversals and block disjointness;
- all 7,887 source Hamiltonian paths, distinguished cycles and complements;
- 166 ordinary source-class pairs for rooted and rooted-converse faithfulness;
- all 1,100 recenterings with Q of orders zero through five;
- complete source/target matrices and recenter permutations for all three
  explicit collisions, including the larger source-order-eight witness.

The independently recovered ordinary source/image class counts are
1/1, 1/1, 2/2, 4/3, 12/9. The converse-merged counts are
1/1, 1/1, 2/2, 3/3, 10/9. The candidate's full automorphism histograms and
all corresponding scalar rows agree. These are image counts of this map,
not counts of all regular tournaments.

The candidate's further sixty deterministic samples at source orders six
through ten are rerun as part of the frozen candidate replay, not claimed
as a separately exhaustive independent census.

## Hostile controls and precise limitations

Mathematical negative controls independently reject a converse second copy,
unrestricted full-automorphism equality, unrestricted canonical-root
recognition, equality with all Hamiltonian cycles, equality with native
full-target noncycle complement, and unrestricted ordinary/merged injection.
For transitive order two the target has automorphism group of order five,
versus one in the source, and two rooted Hamiltonian cycles versus one
matched distinguished cycle.

The n=5 collision also confirms a useful non-necessity point: both the source
and target automorphism groups have order three despite source extremes.
In fact full automorphism preservation alone does not imply unmarked
class-faithfulness. Two different source presentations can identify their
roots with different vertices in one target isomorphism class without any
automorphism moving the root of either individual presentation.

`compare.py` rejects nine deliberately corrupted producer receipts: inflated
class counts, erased merged loss, omitted orders, wrong group histograms,
corrupt matrices, a corrupt recenter map, incorrect converse scope, and
all-cycle/distinguished-cycle conflation, and an inflated regular lower bound. `replay.py` also corrupts a real
relocated producer receipt and confirms that manifest verification fails.
No Python assert statements implement these checks, so -O does not erase them.

## Reproducibility and verdict

The complete candidate copy is pinned in this audit's manifest. Normal and
-O direct replays agree byte-for-byte with both the independent receipts and
the producer's original frozen receipts. `replay.py` copies only the actual
manifest-declared files plus MANIFEST.json into a fresh temporary directory,
then runs independent audit, producer checker, and comparison from a
separate working directory with PYTHONPATH removed and bytecode writes
disabled. Both modes must match the frozen byte-level outputs. No sibling
checkout, package, network, or canonical-label atlas is a dependency.

**Verdict:** accept the structural proofs, attribution correction, exact
finite counts, explicit witnesses, and carefully limited complement bridge.
Do not infer a general A000568 counting formula, transport of native full
C-lines or port clocks, or reconstruction of the unrestricted merged
metagraph after its root/block data is discarded.
