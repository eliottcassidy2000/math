# Independent balanced-blue wiring audit

Scope: the complete native-C **blue-only** induced multigraph on all 229 self-converse balanced order-ten classes, with its exact boundary. This is not an enumeration or topology statement about the black carrier or the full metagraph.

## Independent methods

`audit.py` imports no producer code. It reads the frozen 13,333-row atlas and searches generic arc-compatible bijections against the converse of every row. Exactly 229 rows are self-converse. It independently enumerates their full automorphism groups and anti-involutions and checks the reported refinement errors. Thus it does not depend on the producer's special cyclic-cell canonicalizer.

A direct binary-array enumeration of all 2^20 reflection assignments, independent of the producer's Gray traversal and half-path construction, yields exactly the same 41,699 balanced blue words. Every source path, complement word, target owner and target path/canonical-order record is checked literally. Edge counters are rebuilt from the records. Boundary target classes are checked for pairwise nonisomorphism (using invariant buckets then generic bijections), full automorphism order, self-converseness and zero refinement error.

The exact counts are 22,904 internal ports, 11,452 internal C-lines (24 loops and 11,428 nonloops), 18,795 boundary lines and 786 boundary classes. The graph has 229 nodes, one connected component, no bridges, no articulation nodes and cycle rank 11,224. All 20 positive-error nodes form a connected dominating set: all 209 zero-error nodes have a direct core neighbor. Edge multiplicities partition as 90 core-core, 1,857 core-zero and 9,505 zero-zero.

## Small topology certificate

`topology-certificate.json` contains one triangle and 226 two-edge ears, totaling 455 distinct edges and spanning all 229 nodes. Each ear introduces one new vertex joined to two distinct old vertices. Starting from a triangle, this preserves two-vertex-connectivity. The resulting spanning simple subgraph proves connectedness and absence of articulation vertices and bridges in the blue multigraph. `verify_topology.py` checks all certificate edges against the literal adjacency receipt. The complete independent audit separately checks every adjacency entry against literal C-word records and verifies connectedness after all 229 possible single-vertex deletions.

The ear certificate is logically a graph-topology certificate, not by itself an atlas-completeness proof. Its connection to the full 229-class carrier is supplied by the independent atlas, cube and literal-record checks.

## Controls and limits

Normal and optimized Python runs use explicit exceptions, not assertions. Eight hostile mutations test omitted ports, source path corruption, complement corruption, target-order corruption, edge multiplicity corruption, omitted SC classes, false articulation claims and false target errors. All are rejected in both modes. Four additional topology-certificate corruptions (duplicate ear, omitted ear, degenerate triangle and invalid ear endpoint) are rejected in both modes. The full producer checker is rerun from its manifest-only copied closure under normal and optimized Python, yielding byte-identical JSON equal to the frozen receipt.

Completeness of the original almost-regular atlas remains the inherited frozen-atlas dependency. This audit establishes the stated carrier exactly relative to that atlas, checks pinned source bytes, and independently exhausts the blue reflection cube. It does not re-enumerate every labelled almost-regular tournament or every black presentation. External runtime dependency: NumPy for the independent vectorized binary cube; all other audit modules use Python's standard library.

## Distance geometry

Independent breadth-first searches in `distances.py`, in normal and optimized Python, find unweighted simple-support pair distances 1: 9,219; 2: 16,886; 3: 1 (26,106 unordered pairs total). The blue graph has diameter three; the unique distance-three pair has atlas indices 63 and 86. Both have simple neighbor degree 17 and their neighbor sets are disjoint. Full neighbor lists are recorded in `distances-normal.json`. No prime-family interpretation is asserted.
