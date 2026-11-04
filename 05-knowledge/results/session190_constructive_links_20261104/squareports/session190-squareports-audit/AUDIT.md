# Independent audit: square-sum component ports

Session 190, 4 November 2026. **Verdict: PASS within the stated scope.**

The frozen producer packet is identified by `../session190-squareports/FROZEN_MANIFEST.json`, SHA-256 `a3d9b37e0bc6fcdc353b64934bbc27756bcae7f23b9219ff61cdf2f06c83b91f`. All nine file hashes in that manifest were checked. The proof hash is `634fb598634282b8c17cad9a3c8ce117a2498d1c97fbcced04da2d9af1cd33ae`. No producer correction was required by this audit.

## Mathematical review

1. **All-component transfer.** The encoded endpoint pair determines the selected degree of every old vertex and which exposed vertices share a component. This is precisely the information needed when a new vertex selects zero, one, or two incident edges. Selecting two endpoints in one old component creates a cycle, whereas two different components merge. Deleting the new vertex uniquely recovers the old edge subset and selected neighbor set. Aggregating forests with the same state therefore preserves every coefficient. No component ordering or path orientation is counted.

2. **Singleton multiplicity.** A pair `(v,v)` denotes one actual endpoint vertex, not two stubs. It offers one one-edge choice, and joining two singleton components offers one unordered two-edge choice. The implementation's `b != a` guard and distinct-component condition enforce these rules. The singleton old graph has one Hamiltonian path, with endpoint incidence one for the extension formula. The empty old graph is handled by the empty profile and singleton transition rather than by the nonempty endpoint formulas.

3. **Deletion commutation.** For an old vertex `v`, adjoining `z` to `G-v` with neighbors `S-{v}` gives exactly `(G+z)-v`. This proves the card update independently of path counting. Retaining separate deletion profiles is essential: a label absent from an undeleted endpoint state could be internal in many different places, so erasing labels is not justified. The new-vertex deletion card is exactly the old graph. The same reasoning extends to a fixed old deletion set.

4. **Sharp horizon and pruning.** Every arrival decreases component number by at most one; hence more than `r+1` components cannot finish in `r` arrivals. The proposed chain of `r+1` disjoint edges plus exactly `r` joining vertices is a valid strictness example. Every old forest has at least `r+1` components, and the final graph is one path. This proves worst-case necessity of that component level, not minimality among arbitrary encodings or its necessity for every square-sum instance. Retired singleton, excess retired endpoints, and sealed-component pruning are sound because no future edge can change an old internal vertex's degree or touch a retired endpoint. For a retained penultimate state, every component has an endpoint adjacent to the final vertex, so its one- or two-component coefficient can indeed finish without hidden further constraints.

5. **Three-query response.** For four distinct named ports, prescribed old endpoints force the new vertex to have degree two. Exactly one pairing is incompatible with each probe; the other two each contribute once. The response matrix has eigenvalues `2,-1,-1`, determinant `2`, and binary kernel spanned by `(1,1,1)`. The inverse formulas, even-sum condition, and triangle inequalities follow. The stated converse is correctly restricted to a formal nonnegative integer vector rather than arbitrary graph realizability. Two disjoint edges realize each basis vector, and the specified sector of `K4` realizes `(1,1,1)`. This is not the planar cycle-gluing response matrix.

6. **Endpoint-matrix insufficiency.** The two five-vertex paths have the same sole endpoint pair `(0,4)`. Inserting the new vertex at the existing edge `1-2` works in the first graph. In the second, forced leaf edges and the new degree-two vertex saturate vertices 1 and 2, disconnecting vertex 3. The example is valid and is stronger than scalar-count insufficiency.

7. **Human square-sum certificate.** The displayed `Q25-3` path visits exactly its 24 vertices and every adjacent sum is a square. In `Q26-3`, leaves 18 and 22 force the two endpoints. Vertex 9 forces edge `7-9` in addition to `7-18`, so vertex 7 is saturated. Vertices 13 and 26 force `23-13` and `23-26`, saturating vertex 23. Interior vertex 2, whose neighbors are exactly `{7,14,23}`, then has only edge `2-14` available. This contradiction is a human proof, not merely a zero from search. All eight claimed local adjacency lists were independently regenerated.

## Independent finite verification

`independent_check.py` generates square-sum edges by enumerating squares and complements, rather than importing the producer's adjacency code. Its Hamiltonian counter explores all start vertices and fixes reversal by comparing endpoints; it imports neither the producer DFS nor its path pruning. Its all-component control enumerates edge subsets and traverses their components, without using the producer's union-find verifier. Only transfer functions under test are imported.

Both normal and optimized Python runs passed and produced identical `INDEPENDENT.*.json` files, SHA-256 `b81ecb0b578d00382459cd6cff690507464251e7c47aea4a94ab95532c7f7584`:

- 1,100 labeled graphs of orders 0 through 5: complete component profiles and pruned endpoint matrices.
- 1,099 vertex-addition tests with the new label smaller than all old labels.
- 4,306 deletion-commutation tests.
- 3,072 physical four-port queries on all five-vertex hosts, testing a fifth vertex internal to the two paths.
- Explicit singleton and same-component cycle controls.
- Strict horizon witnesses for every `r=1,...,10`.
- Every endpoint matrix for `Q1,...,Q30` and all 134 deletion-card endpoint matrices for `Q24,Q25,Q26,Q29,Q30`.
- Every reported positive witness and the human certificate's local data.

The producer's complete `verify.py` and `verify_structural.py` were also independently rerun in normal and optimized modes. Both full replays equal the frozen `RESULTS.json` byte for byte, SHA-256 `9c9b23281548bfcee9199c5ce1eab26e340f43ae0880911ab8107bd89a5e787c`; both structural replays equal `STRUCTURAL.json`, SHA-256 `067d12b349fd5c5fdc18b9ab2baea9b0a1b9c45b13dfe902e9bc9dc1caa60fe7`. The duplicate replay files are local evidence and excluded from the audit publication allowlist.

## Scope and inheritance

The cited sections of `session190-planar/PROOF.md` and `session32-square-sum/README.md` were inspected. They do contain the deletion-profile direction and the marked two-path-cover/endpoint-aware suppression mechanisms attributed to them. The present packet explicitly treats those as inherited. It gives elementary transfer/response/horizon proofs plus bounded exact endpoint/deletion data, including the structural `Q25-3` to `Q26-3` loss. It does not prove all-order endpoint control, deletion robustness, an efficient general algorithm, or the external infinite square-sum construction. No global literature-priority claim is certified.

## Reproduction

From the workspace parent directory:

    python3 session190-squareports-audit/independent_check.py
    python3 -O session190-squareports-audit/independent_check.py
    python3 session190-squareports/verify.py
    python3 -O session190-squareports/verify.py
    python3 session190-squareports/verify_structural.py
    python3 -O session190-squareports/verify_structural.py

`FROZEN_AUDIT_MANIFEST.json` lists the complete producer and audit publication allowlist, with hashes. `prototype.py`, `explore.py`, bytecode, and duplicate local replay outputs are excluded.
