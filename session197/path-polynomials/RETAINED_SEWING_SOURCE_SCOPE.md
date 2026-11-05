# Source-scope clarification for retained-root sewing

This note refines attribution without changing any theorem, verifier, bounded census or frozen audit in RETAINED_ROOT_SEWING.md. It records an additional named source read after that packet was frozen:

https://github.com/eliottcassidy2000/math/blob/ee07a5aac9cdfd8d3834fa834df4fb357fcaf769/session195/combined-growth/PROOF.md

That S195 proof already establishes a positive minimal four-coordinate interface for arbitrary moving-root gadgets together with stationary singleton and universal-sink births. Its rank-four continuation matrix is the same determinant16 certificate as the retained-root packet, up to interchanging the last two context rows. Accordingly, the rank certificate in the retained-root packet is an inherited/adapted witness for the narrower stationary grammar, not a priority claim for discovering four coordinates or that numerical minor.

Write K for the S195 stationary singleton birth and D for its stationary universal-sink birth. An old-root attachment of a transitive right block Tr_m marked at position j, counting from0, is exactly the actual graph operation

    D^j, K, D^(m-1-j),

in that order. The j earlier new vertices precede the marked backward-edge vertex in the transitive block, and the remaining vertices follow it. Every new vertex except the K vertex loses to the old graph; the K vertex beats only the retained root. This is equality of labelled graphs when new vertices are labelled in creation order, not just equality of scalar counts.

In particular the two retained transitive-edge letters are KD and DK. From the singleton, the four inputs empty,KD,DK,KDKD have exactly the S195 combined-growth states empty,A,AA,AAD. On these states, the retained continuations empty,K,KK,DK have scalar rows equal to the older empty,A,KA,AA rows. This identifies the determinant witness precisely.

The distinct extension in RETAINED_ROOT_SEWING.md is stationary attachment of an arbitrary nontransitive right block, with the old mark retained, together with its explicit arbitrary-block transfer, target test and complete finite bank calculation. S195 combined-growth's stationary moves only add individual vertices; its arbitrary blocks move the mark. The reviewed retained-only order-at-most5 bank reaches67 and173 without the separate order6 blocks used by the earlier enlarged moving/stationary grammar. This is a difference in the specified construction language and finite bank, not a new claim that these two scalar values were previously unknown as tournament counts.

The all-order scope and finite cutoffs in the published packet are unchanged. In particular neither the positive state, the interval through4999 nor the rank witness establishes the full Hamilton-path spectrum.

Run `python3 verify_retained_sewing_source_scope.py` and its `-O` form. The checker verifies the exact labelled transitive-block decomposition on every rooted tournament through order4 and blocks through order5, and verifies the source state/matrix correspondence. Its only imports are the frozen retained-root checker and its frozen helper.
