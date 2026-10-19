# Independent review: twin-signature parity and tree orientation capacity

2026-10-19. **PASS for the mathematical statements and finite scopes described below.**

Reviewed [the proof](../TWIN_SIGNATURE_PARITY.md), its stated dependencies, final verifier, and minimality certificates. No author proof or verifier was modified.

## Proof review

1. The nine-vertex tree has three interchangeable leaf neighbors at vertex0. Every orientation gives two equal direction signatures, hence a free order-two action on its injective embeddings into every directed host. Its induced S(1,2,3) subtree is asymmetric. This refutes the claimed converse in HYP-9173(b), retaining the earlier sufficient implication.
2. Independent-twin signature permutations preserve every source arc and act freely on injective embeddings. The independence assumption and host/injectivity scopes are explicit and necessary.
3. The exact minimum factorial valuation formula is correct. The recurrence `V_k(m)=q+V_k(q)`, `q=max(0,ceil((m-k)/2))`, follows from binary splitting and monotonicity. Every minimizing multiplicity pattern is realizable by orientation, so the claimed sharpness is correctly limited to this subgroup budget.
4. The rooted binomial recurrence correctly requires asymmetric child orientations with distinct edge-direction/type signatures. The two central-edge factors correctly distinguish isomorphic from nonisomorphic halves; orientation prevents an automorphism from exchanging the central edge endpoints.
5. Every nontrivial oriented-tree automorphism group has an involution by branch recursion. This avoids the false inference that every nontrivial finite group has even order.
6. The arbitrary-directed-host converse follows by taking the source itself as host. The asymmetric P5 tournament hostile is correctly proved by deletion–reversal and Rédei parity and independently checked below. These two host classes are never conflated.

## Independent computations

`independent_tree_capacity.py` uses a separately written leaf-extension tree generator, all-root tree canonicalization, explicit degree-class permutations for automorphism groups, and direct group action on every orientation. It does not import the author's implementation.

- All unlabelled trees through order8: counts1,1,1,2,3,6,11,23.
- Every rooted capacity and every unrooted center formula agrees with direct orientation-orbit counts.
- Every nontrivial oriented automorphism group in this universe has even order.
- The nine-vertex counterexample has undirected automorphism order6 and zero asymmetric orientation classes.
- Independent min-plus optimization verifies the factorial budget at every1<=k<=8 and0<=m<=64.
- Receipt: `capacity_receipt.json`;5,464 explicit checks.

`independent_embedding_witnesses.py` separately enumerates all connected induced subtrees and verifies that the eight supplied witness pairs cover precisely all nonhereditarily symmetric trees through order8: one at7, seven at8, none earlier. It counts embeddings by complete permutations and reproduces all listed count pairs:

    (66,81), (572,213), (165,302), (200,323),
    (46,255), (18,427), (4,245), (145,310).

Each pair has different parity. Combined with the inherited hereditary-symmetry sufficient theorem, these finite certificates establish minimality at order9.

For the asymmetric P5 orientation, independent complete enumeration of all1,024 labelled five-tournaments gives the embedding-count spectrum:

    count4:120 hosts; count6:160; count8:600; count10:144.

Its own directed host has one embedding. Receipt: `embedding_receipt.json`.

## Scope limits

The independent census above stops at order8, plus the named order9 witness. The author's complete orientation check through9 and H/Z census through11 were inspected and are separately reproducible, but not independently regenerated above order8 here. The H/Z table is correctly not described as a classification of good trees. No external historical-priority claim is certified.

Reproduce from this directory:

    python3 independent_tree_capacity.py
    python3 independent_embedding_witnesses.py

The embedding verifier reads the supplied `minimality_witnesses.json` beside its own file; that certificate is an exact copy of the author's frozen certificate, with the receipt recording its SHA-256 hash.
