# Independent audit: ordered root amalgamation and the C3 triangle branch

**PASS on the corrected snapshot.** One singleton wording correction was made before final binding; no counting identity or containment theorem changed.

Frozen source bindings:

- ROOT_AMALGAMATION.md: b773da29cfb8298f3df4b99a741dccad7618d1bc84165f69e9066d86b819f072.
- verify_root_amalgamation.py: ccd82862bae7dd1b73be1412f74f6e3f4e1d339aaa3eba711ee1c998401c6bfc.
- root-amalgamation-certificate.json: fdf718daf8267cc399ecab4a56a9982c8eeb6ebe94ec9dfbf6935c86d42d3e15.
- Named canonical-grammar dependency START_THREE.md: 5b0010444d0339741c1ad3ece83e3a30511457d9c3e9f74f037af18510887753, already independently reviewed in this session.

## Mathematical review

The general Hamilton identity has four genuinely disjoint cases after cutting the word at the shared root. A path outside the root cannot return from the second component to the first. The two cases that use a nonempty initial part from the second component and a nonempty terminal part from the first component are therefore counted by exactly the stated reserved-cover responses. Subtracting d removes the reserved singleton-root path, rather than an arbitrary empty-path convention. The other three cases are ordinary endpoint-restricted words. This checks the meaning and multiplicity of every term.

The starting-two-cover identity independently splits according to whether the reserved root-starting path uses any first-component vertex. In the positive case, the second root deletion occupies two labelled slots without another root condition. In the zero case, the entire first deletion is in the other path. The ending identity is the converse. The deletion product and endpoint products are correct with the explicit singleton exceptions. In particular either singleton amalgamation is literally the identity graph, and the formulas reduce accordingly.

For the C3 templates, when both sources are nontrivial, a cycle cut at the shared root traverses the first deletion, the triangle, and the second deletion in that order. Each exterior endpoint path is unique and the triangle supplies three words. With only the second source singleton the count is 3-j; with only the first source singleton it is j. With both singleton the four-vertex core has cycle count zero or one, never three. Thus exactly the two stated one-singleton boundaries are allowed.

The profile table and all four matrices are consistent with the general identities. One convenient independent algebraic reduction is to first compute the profile X of V amalgamated with the core. Then H of the final graph is

    (R_X-d_X)(U_2-d_2)+(H_X-e_X)d_2
        +(D_2-E_2)d_X+e_X E_2.

Here d_X=3d_1, e_X=e_core d_1 and R_X=12U_1+(R_core-12)d_1. Substituting the four core constants gives the displayed matrices, including every negative coefficient. Their negatives do not weaken physicality because the underlying partition is a nonnegative count. The source-swap converse map reverses the triangle labels and interchanges the two exterior blocks, so its j-to-3-j rule and matrix transpose identity have the correct direction.

## Canonical triangle-tail containment and correction

The cited Type-II grammar has deletion blocks U, P, triangle, Q in dominance order. The root has exactly one entry into U and that entry starts one Hamilton path. Hence V consisting of the shared root, U and P is an actual start-one source, regardless of its other root incidences or whether it is weak.

After any ordinary one-entry suffix and the degree-one closure, all new vertices form a bundle G with precisely its first entry g pointing back to the old root. All other cross arcs are forward. The cycle-cut identity forces e_g(G)=1. Thus the shared root, Q and G form an actual ending-one source: a word ending at the shared root must put Q before all of G and must finish G at g. The two assembled sources have the precise uniform cross incidences of the root-amalgamation template, not just equal numerical responses.

The original proof hash b16d52228a2f17790966f377956a1ddf4d96160272a27844a43d76c448d6ff8c called G nontrivial. An empty ordinary suffix leaves G equal to the single closing vertex. The author corrected this to **nonempty**, explicitly including that singleton case. This does not affect e_g(G)=1 or the nontriviality of the second source, which also contains the shared root. The final binding above includes this correction. The independent checker explicitly tests 42 such closing-only bundles.

For a direct degree-one closure with a triangle in the deletion tail, the same grouping yields j=0, with the second source singleton exactly when there is no post-triangle singleton tail. A nonempty post-triangle tail together with the shared root is ending-one because every tail vertex points to the root. Taking converses supplies the opposite allowed boundary.

These arguments prove containment of the specified canonical triangle-tail branch and the direct triangle-tail closure. They do not prove uniqueness of C3 presentations, cover the mixed/triple-entry branch, or establish a complete four-template C3 classification. The source preserves these distinctions.

## Independent controls

The standalone C++ checker imports no author code. It uses forward subset dynamic programming for all path subsets, endpoint paths and rooted cycles, whereas the source enumerates directed words. Its root-normalized adjacency construction also differs from the source's retained-label construction.

All checks pass:

- 285 rooted labelled inputs through order four and all 81,225 ordered general amalgamation pairs, including singleton inputs;
- 72,900 template/source pairs covering every singleton boundary, including the forbidden ones;
- 72,092 allowed physical C3 pairs and their Hamilton, endpoint, deletion and bilinear values;
- 210 literal canonical Type-II reconstructions after actual ordinary suffixes and degree-one closures;
- 42 empty-suffix singleton bundles, plus 18 direct degree-one triangle-tail closures, including the singleton-source boundary.

The source was copied before running its certificate-writing verifier. Normal and optimized source runs reproduce the frozen certificate exactly. Independent optimized and unoptimized assertion-enabled C++ runs agree byte-for-byte.

Reproduce from this audit directory:

    g++ -std=c++20 -O2 check.cpp -o check
    ./check > independent-normal.json
    g++ -std=c++20 -O0 -D_GLIBCXX_ASSERTIONS check.cpp -o check-debug
    ./check-debug > independent-debug.json
    python3 source/verify_root_amalgamation.py > source-normal.json
    python3 -O source/verify_root_amalgamation.py > source-optimized.json

Hashes:

- check.cpp: 2fa4ae7cfd2f1b98929e214bd01c15a4f7ddd072e8bd9a5260d015c6e06facbc.
- Both independent receipts: 91e6f74f9934b542cfd8ce183ab45a185a1e72c25a4d0e25a066f87a91a79f47.
- Both source replay receipts: fdf718daf8267cc399ecab4a56a9982c8eeb6ebe94ec9dfbf6935c86d42d3e15.

No unresolved mathematical issue remains at the stated scope.
