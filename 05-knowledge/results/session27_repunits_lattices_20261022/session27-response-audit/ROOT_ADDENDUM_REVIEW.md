# Final root integration addendum audit

Verdict: PASS. This review binds new root artifacts separately; it does not amend any frozen response proof.

## Direct odd-cycle parity character

Reviewed ODD_CYCLE_PARITY_CHARACTER.md at SHA256 0b125bd8058aaefd8e16ca81284aab4f81fb0cf8122b4c77ae08d8a8ad3769dc.

The direct cancellation proof is valid. Sum over candidate oriented simple cycles with cyclic rotations identified, then over independently varied edge orientations. Any cycle omitting a varied edge cancels in characteristic two. A cycle containing all edges of the alternating path R must traverse that path contiguously. Closing it without u would have even length 2k−2 and is excluded; using u gives the sole possible odd length 2k−1. The fixed last-base→u edge excludes the reverse traversal, and the closing u→v arc exists exactly when Q_uv=1. Hence the XOR equals Q_uv. At k=2, the only odd candidate is the triangle; a purported two-vertex closure cannot be an odd simple cycle, so the boundary is covered.

The identification of odd-cycle parity with the second bit of H is the classical modulo-four theorem, independently checked earlier against [Grinberg–Stanley, Theorem 7.1](https://arxiv.org/html/2307.05569v1). All odd simple cycles, not merely triangles, must be counted. Rotations are identified, not reversed orientations. The argument retains all port addresses and does not imply that a single uncontextualized cycle count reconstructs a tournament.

## Exact first score-to-H failure

Reviewed SCORE_ROOT_COLLISION.md at SHA256 5b4fe17656af6ebe2ad908639b2ae101b2f1e04ce9fc53bf756467f3de84b36c.

The identity c3=binom(n,3)−sum_v binom(d_v,2) follows by counting transitive triples at their unique source. For n≤4, only a single triangle can occur in an odd-cycle collection, so H=1+2c3. Thus the score multiset, and a fortiori the labelled score vector, determines H through order four. The empty case is unique and causes no exception.

For outgoing masks [12,17,2,22,5] and [10,20,1,22,5], independent full-permutation enumeration confirms H=13 and 11. Independently counting all ordered distinct length-l tuples and dividing by l for rotations gives c3=4 in both and c5=2 versus 1. Both exact labelled score vectors are (2,2,1,3,2). The only changed arcs form the cyclic triangle on labels 0,1,2. Hence n=5 is the exact first failure, not merely an upper bound on where failure occurs.

Further exhaustive independent controls cover all 1,099 labelled tournaments of orders 1–5. In each, c3 agrees with the score formula and H=1+2(c3+c5), with nonexistent lengths interpreted as zero. All score-multiset classes through order four have a single H. These controls use independent direct permutation enumeration, not the root subset DP. root_addendum_independent.py and JSON record them; optimized and normal outputs agree.

The eight-vertex same-E8-root hostile is unchanged from the prior successful audit. The now-portable score_root_hostile.py has inline helpers and no private workspace imports. Its output is byte-identical to the previously audited JSON. Its helper's unhandled empty-graph convention is irrelevant to this script's fixed seven/eight-vertex use; this review makes no general API certification claim for that helper.

## Source binding and boundary

ROOT_ADDENDUM_SOURCE_SHA256SUMS binds these two Markdown artifacts, both root scripts and JSON files. ROOT_ADDENDUM_SHA256SUMS binds the independent controls and this review. The old SOURCE_SHA256SUMS intentionally records the earlier root hostile script version; this separate manifest records the portable replacement while preserving the earlier historical audit. None of the frozen response proofs were changed. No repository or Git writes were made.
