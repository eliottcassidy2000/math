# Independent audit: regular two-complement transport and face incidence

**Verdict: PASS for corrected V3. The initial literal-table audit error is corrected and recorded in `CORRECTION_HISTORY_V3.md`.**

The corrected V3 endpoint/icosahedron packet and the unchanged Gram addendum are valid in their stated scopes. The V1/V2 table contained one invalid displayed path, missed by the initial audit and caught by direct inspection by the final reviewer before publication. The all-order assertions have direct proofs; the order-five classification has a structural proof and an explicit finite adjacency certificate. Independent finite checks are supporting evidence, not substitutes for the all-order arguments. No all-order component classification, ordinary Hamiltonian-path spectrum closure, or black-metagraph connectivity result follows.

## 1. Frozen material and scope

Audited proofs:

- `regular_return_theorems.md` (corrected V3), SHA256 `db2bd384ca2ad4251870723f2e8cea7a872c65ec487ac10fb0a35e84bbef001b`. The original V1/V2 hash was `ea9c2b020b1ed632682b023873806069153fd0eae8655038a23e0955d7925729`.
- `gram_return_addendum.md`, SHA256 `52532f05483cd2b6cf437ea28e6179a2449e3cfa929a2689fbef98be73209046`.

Audited producer programs:

- `check_regular_transport.py`, SHA256 `72f7b7d10ac6464e90f5f3e665f95f3bec67429e215cc38bf861f238cbc67700`.
- `check_gram_return.py`, SHA256 `15852575f9190de42b24a63f3e44b3301a1475a3bab5593d5dfe98931b50ce69`.

The original frozen files and V1/V2 packages remain preserved unchanged. The V3 theorem changes only the mistyped table row; the Gram addendum and both producer programs remain byte-identical. `FROZEN_SHA256SUMS` is a historical record of original staging paths and hashes, not the current V3 package checksum. Use `PUBLICATION_SHA256SUMS_V3_1` for the corrected package.

The objects are labelled tournaments, raw Hamiltonian vertex words, and a graph derived from **two** path-complement steps with an unrestricted intermediate tournament. Path-pair multiplicities, graph degrees, raw H, and automorphism-normalized path counts remain separate throughout.

## 2. Endpoint criterion, parity, and six-edge equality

For any Hamiltonian path from s to t, complementing the nonpath arcs sends the score vector d to `(n-1)1-d+e_s-e_t`. Applying this twice gives `d-e_s+e_t+e_u-e_v`. At positive order greater than one, equality with d is equivalent to `(s,t)=(u,v)`. There is no regularity assumption in this fixed-labelled-score statement. The order-one case is harmless and was checked separately.

The net reversed unordered edges are exactly the symmetric difference of the two path edge sets. A shared edge cannot be traversed in opposite directions, because the first complement preserves the first path. Thus the support size is twice the number of changed successor assignments. Its orientation in the source is balanced for a same-score return.

Closing both paths by the common artificial arc from end to start produces two n-cycles. One changed successor is impossible; exactly two changed successors split the cycle. Hence every nonidentity return has support at least six. At equality, the three cuts divide the first path into four nonempty directed blocks. The fixed start and end force the first and last blocks to remain at the ends. The only nonidentity block order is `A,C,B,D`. Two singleton middle blocks would traverse their common edge backwards and are inadmissible. This proves both necessity and the legal-exchange converse.

The boundary-vertex description correctly distinguishes a directed six-cycle from two directed triangles sharing one vertex. Two vertex-disjoint triangles are excluded. Both permitted types actually occur in regular tournaments: the receipt contains a shared-vertex example at order five and both types at order seven. A six-cycle example in the seven-vertex carousel is

    P=(0,1,3,4,5,6,2), Q=(0,1,5,6,3,4,2).

Its changed arcs form `1->3->6->2->4->5->1`.

The parity obstruction is correctly fixed-label only. Relabelling changes orientation parity by the permutation sign. It does not descend to ordinary isomorphism, and therefore provides no obstruction after ordinary isomorphism and converse merging.

## 3. Exact order-five graph and return weights

The intrinsic map selecting the dominant outneighbour of each vertex is injective. Its lack of fixed points and two-cycles makes it a five-cycle; regularity forces the remaining chords into one of two cyclic orientations, and the definition selects the forward-square orientation. This proves all regular five-tournaments are carousels, with exactly 24 labelled instances and automorphism group of order five.

For two regular five-tournaments, an eight-edge disagreement is impossible: complementing one would produce a balanced nonempty two-edge disagreement, which cannot occur in a tournament. A nontrivial return therefore has six changed arcs. Its support consists of two edge-disjoint directed triangles meeting in one vertex. The displayed legal block exchange and its rotations realize all five possibilities from a carousel.

All twelve codes, all rows of the displayed icosahedral adjacency certificate, both parity components, and the claimed converse-antipode relation were independently verified. Every vertex has five neighbours, 15 identity presentations, ten nonidentity presentations, and total presentation mass 25. Every neighbour has weight two. The corrected V3 table has now been parsed directly from the proof and checked against literal enumeration. The initial audit incorrectly endorsed the V1/V2 table without comparing its displayed words to the enumerated receipt. The invalid word `(0,1,3,4,2)` is replaced by `(0,2,4,1,3)`, whose only matching-endpoint continuation is itself. The counts and all subsequent conclusions remain unchanged.

The conclusion is precisely two labelled icosahedra. The ordinary-isomorphism quotient has one vertex. It is not an identification of the original tournament metagraph with an icosahedron.

## 4. Single-arc and directed-triangle H responses

The single-arc identity correctly splits each affected Hamiltonian word into a prefix ending at one flipped endpoint and a suffix beginning at the other. The endpoint kernels exclude the opposite endpoint, so the relevant induced digraphs are unchanged by the reversal.

For a directed triangle, deleting all three internal arcs is essential. Hamiltonian words use zero, one, or two triangle arcs, never all three. The zero-arc classes cancel. The first cyclic kernel sum counts exactly the one-arc difference. Two triangle arcs necessarily form one consecutive three-vertex block, and the second sum counts exactly the two-arc difference. The two terms were independently checked **separately**, not only by verifying their sum.

The scope is every tournament, with no regularity assumption. The formulas supply endpoint-refined information; score and triangle-count preservation alone do not fix these responses.

## 5. The all-score-fibre Gram theorem

The binary incidence definition is justified. Agreement edges between a source T and its path complement S are exactly the chosen directed Hamiltonian path. Consequently the path witness is unique, including the trivial one-vertex case, and the row mass of B is raw H.

For a fixed labelled score vector d, every common intermediate S of T and U gives a unique two-step presentation. The U-to-S path remains directed in S and is an S-to-U path under the same complement. Its endpoints match the T-to-S endpoints by the score identity. Conversely every two-step same-score return yields that common intermediate. Thus `W=BB^T` counts the intended multiplicities exactly, with `W[T,T]=H(T)`.

Symmetry, positive semidefiniteness, the intersection bound, and the three-by-three principal-minor inequality follow with their stated scopes. The orthonormal orbit compression warning is necessary and correct: its diagonal includes within-orbit return weights, not merely original diagonal entries. At order five the compressed scalar is 25, whereas raw H is 15 and H/|Aut| is 3.

## 6. Actual intermediate tournaments realize the faces

The degree-three intermediate columns are not an abstract factorization manufactured from the known graph. The path-orbit calculation identifies them as actual tournaments of form `t -> s -> C3 -> t`. Their three matching-endpoint paths are the three rotations of the interior directed triangle. There are exactly `5*4*2=40` such labelled tournaments.

The reconstruction proof is sound. Across the three regular neighbours, `s->t` is the unique unanimous arc. Every spoke occurs in exactly one selected path, and every interior-cycle arc in exactly two. The majority orientation on the remaining three vertices recovers the original C3. Thus the unordered triple of regular neighbours recovers the entire labelled intermediate tournament, proving injectivity. The forty distinct graph triangles exhaust the twenty faces in each icosahedron.

The independent checker reconstructed all forty intermediate tournaments in full, rather than checking only the unanimous arc. It also checked the onto statement against **every** graph triangle. There are 280 intermediate columns in total: 240 private columns, ten incident to each regular source, and forty face columns. Hence, on each sheet,

    W = 10I + F F^T = 15I + 2A.

Both stochastic kernels are correct and were independently checked entrywise with rational arithmetic. Uniform path-pair sampling gives `W/25`; uniform first-path sampling followed by a uniform legal continuation gives `(7/9)I+(2/45)A`.

## 7. Exact finite spectrum

The integer identity `A^2=3I+2J-2C` is correct. The even antipodal quotient is K6. On the odd subspace the square is 5I, and zero trace gives equal multiplicities of the two square-root eigenvalues. The resulting spectrum of `15I+2A` is exactly the one in the addendum.

An independent integer characteristic-polynomial calculation checked

    det(xI-A) = (x-5)(x+1)^5(x^2-5)^3,
    det(xI-W) = (x-25)(x-13)^5(x^2-30x+205)^3.

No floating eigensolver or numerical tolerance was used.

## 8. Reproducible independent universes

`independent_audit.py` uses sets of directed arcs and literal vertex words, not the producer's bitmask endpoint DP. Its endpoint kernel inserts a prescribed block into every outside word at every cut. It checks:

- Every labelled tournament at orders 1 through 5: 1, 2, 8, 64, 1024 states.
- Every single arc at those orders: 0, 2, 24, 384, 10240 checks.
- Both triangle-response components for every directed triangle: 0, 0, 2, 64, 2560 checks.
- Every second-word permutation relative to a canonical first path through order 7, and every matching-endpoint second word through order 10. Reverse shared edges are rejected explicitly as incompatible, not silently omitted. These are path-pair geometry universes, not regular-tournament censuses.
- Every regular five-tournament and every return presentation, including graph certificate, weights, automorphisms, parity components, and converse antipodes.
- All fourteen triangle switches in each of the three supplied regular order-seven representatives.

`independent_gram_audit.py` uses those independent arc-set routines and no producer code. It enumerates every labelled score fibre through order 5. At order five this means 291 score fibres, 7680 path incidences, 7736 ordered same-score state pairs, and 9600 same-score return presentations. It additionally checks all forty full face reconstructions, all 576 entries of the regular order-five factorization and both stochastic kernels, and both exact characteristic polynomials.

Both independent programs use explicit exceptions. Normal and optimized Python runs produced byte-identical receipts. Both producer programs were replayed from audit-local copies; normal and optimized outputs matched all supplied receipts exactly. The producer Gram checker imports the producer transport checker, so those two files must remain beside one another. The independent Gram checker similarly requires `independent_audit.py` beside it.

From the publication package directory:

    python3 check_regular_transport.py
    python3 -O check_regular_transport.py
    python3 check_gram_return.py
    python3 -O check_gram_return.py
    python3 independent_audit.py > independent_audit.replay.json
    python3 -O independent_audit.py > independent_audit.optimized.replay.json
    python3 independent_gram_audit.py > independent_gram_audit.replay.json
    python3 -O independent_gram_audit.py > independent_gram_audit.optimized.replay.json

Compare each replay with its corresponding supplied JSON receipt.

The V3 source-aware regression checker `check_literal_examples.py` additionally parses and validates the proof's displayed path table, explicit return pair, triangle lists and pairs, all certificate codes, the Gram addendum's explicit intermediate masks/code/cycle/neighbours, the audit's seven-vertex example, and the SVG's 24 labels and 60 edges. It tests the original mistyped word as a hostile control. See its receipts and `CORRECTION_HISTORY_V3.md`.

## 9. Inheritance and non-consequences

The frozen packet correctly credits `spectrum_basins_20261020e/regular.md` for classical score-fibre connectivity, the cycle-to-triangle compiler, the regular order-seven census, one-step endpoint defects, and the fixed-label parity boundary. It does not promote these to new claims. The inherited H-values 171, 175, and 189 correctly refute an intermediate-value inference from score-fibre connectivity and increment gcd two.

Additional close precursor, recorded here without changing frozen bytes: `05-knowledge/results/merged_metagraph_endpoint_transport_20261020/endpoint-transport-and-complement-leakage.md`, Section 4, especially equation (4.2), already identifies the almost-regular shell's high-to-low endpoint kernel for returns to regular nodes. Its Section 5 separates true half-edge complement from a compressed node-level two-step surrogate. The present labelled path-pair graph and exact common-intermediate incidence refine that boundary; general Gram positivity itself is standard.

This audit makes no exhaustive external priority claim. It neither re-proves nor independently updates the literature status of the global raw-H spectrum conjecture. The inherited source audit states it OPEN; nothing in either audited proof changes that status.
