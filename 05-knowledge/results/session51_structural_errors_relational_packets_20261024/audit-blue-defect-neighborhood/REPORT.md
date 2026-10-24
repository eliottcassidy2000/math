# Independent audit: complete blue neighborhood of the defect slice

24 October 2026. **PASS, finite-exact for every blue edge incident to the 24
positive-error self-converse order-ten nodes.** No full-metagraph claim.

Candidate manifest SHA256:
`141cbcb8e2696238c2db7db40e4bcd9570438d791a4c405e662427671db8f6da`.
The underlying independently reconstructed wiring is frozen at audit manifest
`53b674adb896c5f689e38b82fcbefd3a21ea0dfe9c19e5d96dba8bb43f806c9c`;
its producer wiring manifest is
`cd89cad144786c67f076c235409c90e3b4de21e0a5f81926139184068f1d676d`.
All exact copied source hashes are retained in `SOURCE_PINS.json`.

## Result

The 4020 distinct boundary presentation words belong to exactly **1031 distinct
ordinary target classes**, all self-converse and all of refinement error zero.
Their complete stable partitions and full automorphism orders are:

- 1016 classes with 1^10 stable cells, all rigid.
- 11 classes with 1^4+3^2 stable cells: ten of automorphism order 3 and one
  of order 9.
- Four classes with 5^2 stable cells, all of automorphism order 5.

The graph of all edges incident to the positive slice has 1055 vertices and
4110 edges: 90 internal positive/positive edges and 4020 positive/zero edges.
The boundary adjacency occupies 3460 positive/zero pairs, with multiplicity
histogram 1:2977, 2:413, 3:63, 4:7. The maximum degree *from the positive slice*
of one target class is 17; this is not its full blue degree.

The incident-edge graph has two connected components:

- Four degree-one positive nodes and 238 zero-error neighbors.
- Twenty degree-two positive nodes and 793 zero-error neighbors.

These are components of this **incident-edge graph**, not necessarily of the
full blue graph, nor of the vertex-induced graph on these 1055 vertices.
Zero/zero edges and all black edges have not been enumerated here.

## Independent exact class identification

`independent.py` imports no producer implementation and does not use intrinsic
cycles or a rotation-only search. It starts from all 4020 boundary words in
the independently audited wiring receipt, decodes their full 10-by-10
matrices, and recomputes actual stable refinement from one color.

It canonically orders the actual stable cells by retaining the old cell and
sorting complete outcount signatures at every round. It then minimizes the
**entire directed adjacency matrix over every within-cell permutation**, not
merely rotations. Search sizes are 1, 36, and 14400 on the three actual stable
patterns. They occur respectively on 3959, 45, and 16 boundary presentations.
This is an exhaustive full-arc isomorphism canonicalizer: any isomorphism
must preserve the isomorphism-invariant ordered stable cells, and every
permitted within-cell bijection is tested against the complete matrix.
It does not infer equivalence from degree profiles or stable-cell sizes.

The independent canonical representation and class labels differ from the
producer's. After exact label alignment, all 1031 class identifications,
all 3460 boundary adjacency multiplicities, and all 4020 actual edge witnesses
agree. Every producer canonical vertex-order witness is explicitly applied
to the original target matrix and reproduces its claimed full target matrix.
The independent receipt retains its own full canonical matrices, vertex
orders, complete membership lists, and a bijection from producer class IDs.

For every target class, the audit enumerates **all** within-stable-cell
permutations and keeps exactly those preserving every arc. Thus its full
automorphism group is not assumed to be a subgroup of a guessed rotation
group. Actual orbits and the defect

    e(T) = sum_O ( |C(O)| - |O| )

are recomputed directly, giving zero for every target. The producer's actual
stable cells, vertex orbits, automorphism orders, and incidence counts all
match. Every boundary normalized word has literal positional anti-reflection;
the converse of every independent class matrix also has the same exhaustive
canonical key. There is no further ordinary/converse class quotient.

## Review of the producer's small-cell proof

The producer's faster search is mathematically complete. The inherited
anti-involution stable-cell lemma pairs equal odd stable-cell sizes at even
order. At order ten, the possible stable-size multisets are exactly 1^10,
1^4+3^2, and 5^2. A three-cell is the directed triangle; a regular five-cell
is C5. Their intrinsic successor relation is, respectively, the unique
internal outneighbor or the source among two internal outneighbors. Internal
isomorphisms therefore preserve those oriented cycles and differ only by
rotations. Sorting complete refinement signatures fixes the cell order.
Minimizing the full matrix over 1, 9, or 25 rotation choices gives a complete
class invariant in this restricted universe. The independent exhaustive
search verifies the result without relying on that optimization.

Target automorphisms do not reduce the 4020 edge multiplicities. The line
endpoints are distinct normalized presentation words, not arbitrary target
Hamiltonian paths; complement is injective on all source presentation words.
The audit preserves every source word and every actual complementary target.

## Structural separation and scope

The inherited actual-path endpoint-degree identity makes the two source
sectors' possible target score multisets disjoint:

    degree-one source: 2,3^4,6^4,7 or 3^4,4,5,6^4;
    degree-two source: 3,4^4,5^4,6 or 4^5,5^5.

Thus no target class can meet both source sectors. This proves separation,
while the exact adjacency computation additionally proves connectivity within
each of the two stated incident-edge components. It does not exclude paths
between them through unenumerated zero/zero edges.

The result closes boundary identification for this blue defect slice. It is
not a full order-ten tournament census, does not determine A000568 totals,
does not recover zero/zero or black adjacency, and does not claim that the
1055-vertex graph is induced on its vertex set.

## Reproduction and adversarial controls

Normal and -O independent receipts agree byte-for-byte. The validator makes
literal copied directories under a separate temporary root and runs both
independent and producer checkers in normal and optimized Python from `/tmp`.
It verifies strict candidate manifests and equality with the original frozen
receipts, then confirms that every original producer file is unchanged.

Both normal and optimized modes reject changed pinned source data, a semantic
bug restricting exhaustive canonicalization to the identity permutation, and
a wrong actual canonical vertex-order witness. The last hostile deliberately
repins a disposable corrupted receipt so rejection tests the full-arc semantic
check rather than only its SHA256 guard. Candidate manuscript corruption is
also rejected by its strict manifest. No production source was mutated.
All checks use explicit exceptions and remain active under -O.

The complete portable independent receipt is `normal.json`; the actual-copy
and hostile-test record is `validation.json`. No publication was performed.
