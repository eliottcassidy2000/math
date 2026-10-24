# Complete blue edge-neighborhood of the order-ten defect nodes

24 October 2026. **FINITE-EXACT**, with an exact small-cell canonicalization
proof. This packet closes the boundary identification left open in the frozen
induced-wiring packet. It reconstructs every blue edge incident to the24
positive-error self-converse nodes, but not the other edges between their
zero-error neighbors and not any black edges.

## 1. Result and frozen dependence

The source wiring packet is pinned to manifest SHA256
`cd89cad144786c67f076c235409c90e3b4de21e0a5f81926139184068f1d676d`.
It already pins the audited converse-error packet, and supplies the complete
4200 blue ports on the24 defect nodes, paired into90 internal edges and4020
boundary edges. Every boundary target is self-converse with refinement error
zero; the source packet did not identify which boundary targets were isomorphic.

**Exact result.** The4020 boundary edges meet exactly1031 distinct zero-error
self-converse nodes. The blue edge-neighborhood therefore contains1055 vertices
(24 defect nodes and1031 zero-error neighbors) and4110 blue edges incident to
the defect set (90 internal and4020 boundary). Its complete boundary adjacency
and a labelled canonicalization witness for every boundary edge are recorded
in `normal.json`.

The1031 boundary classes have actual stable-cell patterns and full automorphism
orders:

- 1016 classes with ten singleton stable cells, all rigid;
- 11 classes with3+3+1+1+1+1 stable cells: ten have automorphism order3, one
  has order9;
- 4 classes with5+5 stable cells, all with automorphism order5.

Every class has error zero by a direct full automorphism/orbit calculation in
its complete small-cell rotation group, not merely by inheritance from the
boundary argument. Of the4020 distinct target presentation words,3959 have
ten singleton cells,45 have3+3+1+1+1+1, and16 have5+5.

There are3460 occupied positive/zero node pairs. The boundary edge multiplicity
histogram on these pairs is `1:2977, 2:413, 3:63, 4:7`, summing to4020 edges.
The largest number of edges from the defect set incident to one zero-error
neighbor is17. This is its degree into the specified defect set, not its full
blue degree.

The graph consisting only of these4110 incident edges has two components:

- the four crossing-degree-one defect nodes and238 zero-error neighbors;
- the twenty crossing-degree-two defect nodes and793 zero-error neighbors.

These are not asserted to be components of the full blue graph: zero/zero
edges were not reconstructed and can connect the two sets or lead elsewhere.

## 2. Why target isomorphism reduces to at most25 rotations

Every boundary target is blue, hence self-converse. The proved anti-involution
stable-cell lemma in the frozen converse packet says that every actual stable
cell of an even-order self-converse tournament is paired with another cell of
the same odd size. At order ten this leaves exactly

    1^10,    3^2+1^4,    5^2.

This restriction is on the **actual stable partition**, not an arbitrarily
chosen equitable block presentation. The checker computes it independently by
iterative directed refinement. Colors are named by sorting their complete
signatures at each round, so actual cell order is invariant under directed
isomorphism.

Every size-three cell is the directed triangle. Every size-five cell is the
unique regular five-tournament C5. Their cyclic coordinates can be recovered
intrinsically:

- in a triangle, successor is the unique internal outneighbor;
- in C5, successor is the source of the directed edge between the vertex's
  two internal outneighbors.

An automorphism must preserve these successor cycles. Conversely each cyclic
rotation preserves the internal tournament. Thus an isomorphism between two
boundary targets, if it exists, preserves the intrinsically ordered stable
cells and acts by an independent cyclic coordinate shift on each nonsingleton
cell. There are respectively1,9,25 such choices in the three cases.

For each target T, the checker traverses those choices and minimizes the
**entire normalized directed adjacency tuple**, including every crossing
entry between cells. The result is an exact ordinary-isomorphism canonical
key. No parity-only crossing invariant, gain kernel, independent marginal,
or incomplete holonomy replaces the complete matrices.

The smallest labelled vertex is used only to choose a provisional start on
each intrinsic cycle. All rotations are minimized over, so this gauge choice
disappears from the final key. The checker tests idempotence, converse
invariance, and two nontrivial global relabellings of every resulting class.
The normalizer does not enumerate arbitrary10! vertex permutations, and its
completeness follows from the actual-cell and intrinsic-cycle proof above.

## 3. Exact edge aggregation and zero-error verification

For each boundary line in the frozen receipt, decode its target36-bit word
using the fixed base Hamiltonian path. Its canonical vertex order is retained
as an explicit witness, so normalizing the full target matrix by that order
must reproduce its target class key. Grouping these keys gives the1031 classes.
Aggregating by the ordered pair (source defect node, target class key) gives
the entire boundary multigraph, with all4020 distinct edges retained.

Each candidate automorphism on a canonical target is one of the same at most25
cell rotations. The checker tests every such map against every directed arc,
forms its full vertex-orbit partition, and computes

    e(T) = sum_O ( |C(O)| - |O| ).

It obtains zero for every boundary class. This is a complete automorphism
search in this proved small-cell universe, not a partial subgroup test.
There is no additional converse quotient: every boundary target is already
self-converse, as independently witnessed by the positional reversal of its
blue normalized word.

No edge is lost by the target automorphism groups. The distinct presentation
words, not the target Hamiltonian paths, are the line endpoints. Complement
is injective on the distinct source words. Thus all4020 target words are
distinct even when some belong to a target class with nontrivial automorphisms.
Class canonicalization groups endpoints while preserving their actual edge
multiplicity.

The per-source boundary totals agree with the frozen wiring receipt. Summing
all target incidence counts gives4020, and adding the frozen90 internal
positive/positive lines gives4110. Edges joining two zero-error nodes are
outside this enumeration; the word “complete” refers only to all blue edges
touching the24 positive-error nodes.

## 4. Structural reason the two source families share no boundary class

The inherited endpoint-score identity already explains the four degree-one
nodes' isolation from the positive-induced graph. It also shows that degree-one
and degree-two sources cannot have the same boundary neighbor class. The
possible sorted target degree multisets are respectively

    degree-one source: 2,3^4,6^4,7   or   3^4,4,5,6^4;
    degree-two source: 3,4^4,5^4,6   or   4^5,5^5.

They are pairwise distinct. Thus no target class can be adjacent to both source
families. The literal adjacency computation additionally shows that the four
degree-one sources become connected through their238 common-neighborhood
vertices, while the already connected degree-two sources have793 neighbors.
This proves the stated two-component decomposition only for the graph of
edges incident to the defect set, not for the full blue graph.

## 5. Reproduction and stopping boundary

Run `check.py` under normal Python and `-O`; the JSON receipts are byte-identical.
The source files are checked both against their local SHA256 pins and against
the exact immutable wiring manifest. Real-copy relocation reproduces both
modes. Every boundary edge has a canonical full-matrix witness, and every
boundary class has its canonical rows, actual stable cells, automorphism order,
vertex orbits, and defect-set incidence degree in the receipt.

The task stops with this complete edge-neighborhood. There is no further
claim about zero-error/zero-error wiring, black wiring, the full order-ten
metagraph, or a numerical A000568 total.
