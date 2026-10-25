# Session57: square-sum paths and square-magic incidence

25 October 2026. Second verified checkpoint of an ongoing multiwave session. The question is not dismissed as analogy: we test candidate maps, retain what they preserve, and repair the carrier when an implication fails. No integer magic-square solution or external priority is claimed.

## Wave 1: identify the exact objects

- [Carrier theorem](carriers/THEOREM.md): a square-valued magic grid is exactly an injective embedding of the eight-line grid hypergraph into the fixed-total triple hypergraph H_T. Equivalently, use eight additive three-step paths synchronized on the original cell labels. The line-intersection graph's 22/28 density comes from the grid alone, not arithmetic. [Independent audit](audit-carriers/REPORT.md).
- [Eight source-pinned inherited routes](inheritance/CONNECTION_BOARD.md): common-center opposite-square gaps, Pythagorean square-sum petals, named ports, fixed-endpoint path certificates and their precise boundaries. The exact short-gap map sends d in D_e with d<e²/2 to the positive square-sum AP triangle (2e²-4d,2e²,2e²+4d). At e=1105 there are 13 gaps, five short gaps, but no required additive diamond.
- [Character/graph separation](character-graphs/THEOREM.md): a Parker field can have a square-sum Hamiltonian path on nine squares (F17), while a valid magic grid can have an isolated vertex in its entry-value square-sum graph (F59). In characteristic 3 the graph on a magic entry set collapses to K9 or the empty graph according to the square class of -1. A fixed affine magic plane has exactly 432 magic arrangements among 9! permutations. [Independent audit](audit-character-graphs/AUDIT.md).

The graph uses distinct vertices and admits an edge when their sum is a field square; the zero-sum edge convention is always stated. Hamiltonian paths in the reported graph counts are counted up to reversal. Extension-field arithmetic is not integer arithmetic modulo the field's cardinality.

## The wording distinction

An ordinary integer magic square whose entries are squares has center e² and common total 3e². For e>0 that total is not itself an integer square. Opposite entries sum to 2e², also nonsquare. The separate problem in which different line sums are allowed to be different squares must therefore not be silently identified with the equal-total problem.

Finite-field square status does not make canonical representatives integer squares. The earlier Session 55 no-wrap and complete modular classifications remain intact; discussion of Cain's 2019 conjectures refers to that source's historical status, not a reopening of the inherited repository results.

## Waves 2–3: a positive exact path bridge

- [Decorated K4,4 equivalence](biclique-bridge/THEOREM.md): reflection-matched parts and a parallelogram of offsets make the sixteen edge sums exactly the nine magic values with known multiplicities. Square scaling and a gauge give eight distinct positive integer vertices conditionally on an integer magic square. A bare biclique is insufficient; vertices need not themselves be squares.
- [Exact path certificates](biclique-path-counts/README.md): K4,4 has 576 unoriented Hamiltonian paths and 72 cycles. Two paths suffice to inspect all nine formal labels; exactly 204 optimal pairs have nine union edges and five shared edges. Three paths are needed to cover all sixteen identified edges. [Independent audit](audit-biclique-paths/AUDIT.md).
- [Nine-edge linear frame](two-path-frame/THEOREM.md): the two path cycles provide only two linear constraints; four additional reflection/parallelogram constraints are exactly what remains to recover the six magic equations. The shared-coordinate conditions cannot be inferred from the path skeleton alone.
- [Smaller AP-biclique bridge](ap-biclique-bridge/REPORT.md): two arithmetic-progression parts of K3,3 reduce the carrier to six vertices and two extra diagonal conditions. Two Hamiltonian paths sharing one edge cover all nine identified edges. Their 36 optimal pairs form six partner hexagons. The positive integer normalization makes the selected graph induced, using residue classes1 and3 modulo4. Choudhry's semimagic permutation is explicitly credited. [Independent reduced-frame audit](audit-reduced-frames/AUDIT.md).

The two-path statement is a certificate with its numerical decoration retained, not an assertion that arbitrary intersecting paths imply magic. Minimum coverage means minimum within the declared edge/label-inspection model, not algebraic independence of square predicates. Over odd fields the collision-free gauge count is exactly q-9 when the nine labels are distinct; zero labels require the zero-allowed convention.

## Verification and continuation

All mathematical claims have independent proof/computational audits. Assertion-based producer path/AP scripts are validated normally only; their -O execution is not described as effective assertion validation. The corresponding independent exception-safe audits pass both modes. Path-writing scripts are relocated only in publication copies with receipts; frozen source lineage is retained.

The next waves test bounded integer bicliques, published all-parameter families, the linear/algebraic-geometric locus, and the separate unequal-square-total interpretation. The earlier Session55 modular classifications remain closed in their stated scopes. The actual nondegenerate integer magic-square existence problem remains open, and the session is still active.
