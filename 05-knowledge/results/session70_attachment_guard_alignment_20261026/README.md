# Session70: legal guards reconstruct tournaments; minimal carriers depend on the operation

2026-10-26. **PROVED** all-order statements below, with explicitly cited classical inputs where needed. **FINITE-EXACT** checks support but do not replace the proofs. No literature-priority assertion, unrestricted prime classification, or solution of the Graceful Tree Conjecture is claimed.

## Main closure: yes, at every finite size

**The complete set of named legal re-marking guards reconstructs the tournament.** Equivalently, a labelled tournament is uniquely determined by its set of directed spanning Hamilton paths.

A useful exact intermediate carrier is the union of all arcs occurring in those paths. It has precisely the same Hamilton-path language. Although some tournament arcs are absent from that union, the omitted orientations can be decoded from vertex-deleted reachability inside strong components. The proof uses the classical Busch–Jacobson–Reid structure theorem for nontraceable arcs and Camion's theorem; the source statements, hypotheses, and proofs were inspected.

For fixed-path presentations X_n, this closes the prior C5 all-order question: every guard-respecting quotient has exactly 2^M valid presentation states, M=(n−1)(n−2)/2. The exact fixed-length cost for these valid states is M bits; explicitly totalizing failed operations adds a separate failure sink. Adding complement or attachment operations cannot make that quotient smaller.

**The operation names matter.** One-step elementary relocation guards already collide at order four. Full left/right relocation traces do reconstruct, using the inherited relocation-connectivity theorem. Anonymous fresh re-marking and weighted aggregate responses are different experiments, not counterexamples.

- [Guard reconstruction proof, decoder, and local-guard boundary](session70-guards/REPORT.md)
- [Independent proof, primary-source verification, and decoder audit](session70-guard-audit/REPORT.md)
- [Cross-domain synthesis and remaining frontier](session70-root/SYNTHESIS.md)

## Why alignment memory grows

The inherited 2^k alignment fibre at order 7k+2 has identical separate rooted endpoint types and colour. This session proves its common strong source is **rigid** for every k, so the loss is not caused by source symmetry.

Every conjunction of hidden block bits is an actual simultaneous named re-marking guard. The resulting Boolean zeta matrix has determinant one and rank 2^k over every field. Keeping all named moments through block order t leaves exactly 2^k−sum_{j≤t}binom(k,j) missing linear coordinates for arbitrary laws. Even/odd parity mixtures agree on all proper-block marginals but differ on a native global guard.

An individual alignment still needs only k bits, read by singleton guards. Exponential linear-law dimension is NOT an exponential deterministic-bit lower bound. This recovers the useful finite-moment mechanism from prior Vitali work without importing infinite nonmeasurability.

- [Native construction, proofs, and historical recovery](session70-vitali/REPORT.md)
- [Independent initial audit](session70-bridge-audit/AUDIT-INITIAL.md)

## Square sums, graceful trees, Petersen, K5, and K3,3

The graph lane proves two other exact, operation-relative minima:

- **Unrestricted forest attachments:** precisely (2r−1)!! pairings of 2r distinct named endpoints. Different pairings are separated by a connector that merges paths in one state and closes a forbidden cycle in the other. Fixed square-sum hosts may allow fewer contexts; unrestricted minimality is not silently transferred.
- **Graceful overwrite guards:** in the inherited factorial family, the minimal validity-only quotient is the k! color-to-component permutations, genuinely forgetting the 2^k left/right signs. This requires no-op requests to be allowed and retains the exact declared operation semantics.

The recovered arbitrary frozen-color width-two graceful connectivity claim was already refuted by a five-state fiber split 3+2; the restricted upper-tail theorem survives. Neither is the Graceful Tree Conjecture.

Petersen contracts along its spokes to K5, but Hamiltonicity is lost unless local ports, repeated visits, and global connectedness are retained. The full local response yields six 2-factors, each two five-cycles. Deleting one vertex and suppressing three degree-two vertices gives K3,3: exactly two of its six Hamilton cycles use all three marked suppressed edges and lift to spanning nine-cycles. K3,3 is gracefully labelable; K5 fails an odd-difference count. These distinguish planarity, Hamiltonicity, and gracefulness rather than conflating them.

- [Graph proofs, exact controls, and recovered square-sum/graceful sources](session70-graph-labeling/REPORT.md)
- [Independent graph and graceful-guard audit](session70-bridge-audit/AUDIT-GRAPH-SUPPLEMENT.md)

## What is inherited and what remains open

[Prior strong-atom C1-C5](../strong_atom_continuation_vitali_20261026/README.md) already supplied the aligned closed carrier, original k-bit obstruction, and operation-kernel criterion. [Session69 P2-P3](../session69_prime_ports_power_aliquot_20261026/session69-tournament/REPORT.md) already closed ordinary, uncoloured all-arity response composition. This session preserves those results and closes the specific full named-guard reconstruction question.

The unrestricted modular-prime atlas, efficient/economical aggregate composition, optimal query complexity, and minimal quotients for weaker or restricted operation menus remain open. Incoming HYP-9416 concerns fresh unnamed successor multisets; its separate finite claims are not promoted to all-order theorems here.

## Verification

The packet contains curated native sources, exact provenance manifests, producer and independent-audit normal/optimized receipts, hostile controls, and a portable packet verifier. Downloaded full third-party PDFs/OCR are not redistributed. Classical-source URLs and SHA-256 identities are recorded in the guard audit.

Run `python3 -B verify_packet.py --replay` from this directory. Its manifest records the precise replay universe; no external network or unpublished checkout is needed for the frozen replay. Publication receipts distinguish exact remote byte verification from CI status; an empty check list is not called a CI pass.
