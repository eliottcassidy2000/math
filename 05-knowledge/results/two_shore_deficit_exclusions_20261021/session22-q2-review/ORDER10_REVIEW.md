# Independent finite order10 audit

Date: 2026-10-21 UTC. Verdict: PASS within the complete marked q=2,m=5 frame universe. This is separate from the all-m analytic insertion theorem.

The finite universe contains exactly 4*4*8=128 frames: two independent five-shore cyclic orientations and three free cross bits, with the other two bits forced by the marked Hamilton cycle. Every hypothetical order10 deficit-six center occurs in this universe by the inherited q<=8 and odd-quotient reduction. Marked-frame duplicates do not affect an emptiness result and are not counted as ordinary classes. Step-two automorphism and the odd automorphism-group order make the exact index q=2.

Reviewed checker source SHA256: 34fc39b5819d44f346250ee8d84bb2118a20cdaee13b128d1ed526f59f6c27f7. Fresh copied normal/-O runs each pass 1,088 always-active checks and have identical nonvolatile receipt content. Real timing fields differ.

Source inspection confirms the producer independently compares subset-DP and directed-word DFS Hamilton counts, degree-compatible map backtracking and exhaustive shore-permutation products, and direct word reversal tests versus anti-involution half-word enumeration. The shore restriction is complete because the two full degrees differ; automorphisms preserve the shores and anti-automorphisms interchange them. Every blue word has a unique positional reversal, preventing duplicate credit across anti-involutions.

An independent reviewer implementation recomputes all 128 Hamilton counts using a backward endpoint DP. It also recomputes blue counts by a distinct pair-mask DP: for each fixed-point-free involutory anti-map, choose one member of each swapped pair in an oriented half-path, ending with the required middle edge x→a(x). The anti-map property supplies the mirrored suffix arcs. This avoids the producer's explicit half-word permutation enumeration. Reported map permutations are checked against every source edge; map-set completeness is supported by the producer's two inspected exhaustive algorithms, not asserted as a third independent map enumeration.

The recomputation agrees on 64 self-converse marked frames and on the deficit histogram: 1528 (8), 1628 (8), 2908 (12), 2926 (12), 3064 (12), 3114 (12). No deficit is six. The minimum is 1528. The independent script also validates the structural packet's twelve literal construction samples; together it records 1,614 explicit checks, including rejection of one intentionally corrupted receipt.

This finite zero result is independently consistent with the stronger analytic theorem excluding q=2,m>=5. No n9 replay, higher-order finite census, ordinary-class histogram, or metagraph reconstruction is claimed.

Files: independent.py, independent_results.json, replay_results.json, source_pins.json. Copied normal/-O replay directories can be omitted from publication because original producer receipts retain their provenance. Only the review directory was written; no Git or remote writes occurred.
