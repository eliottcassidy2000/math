# Independent audit: target realization and physical interfaces

**PASS.** Reviewed the stable `ports` packet, including its exact family-completeness bound, named adjacency construction, endpoint compiler proof, field codecs, and distinction between count realization and actual named-source continuation.

## Independent computation

`verify.py` does not call the author's six-state compiler. It fixes the last vertex before the root and first vertex after the root, enumerating all possible actual neighbors (or empty segments). Once those two boundary vertices and the root's path slot are fixed, each remaining ordered vertex independently chooses an allowed segment. The product of the numbers of allowed segment choices counts precisely the compatible covers. Endpoint restrictions intersect these per-vertex choices. Summing boundary choices gives every endpoint and directed joint. This is a different computational decomposition from the author's scan-state recurrence.

For all six representations, independently recomputed and matched:

- Every literal row, root cut, nonroot transitive order and strongness condition.
- H=4611841.
- All five physical fields E1,S1,deletion-H,E2,S2 at all24 vertices.
- All552 directed joints per representation.
- All golden vectors, F49/charge and Frobenius coordinates, and F343 orbits/degrees.

The total is144 endpoint records and3312 directed-joint records. Polynomial long division independently recomputes F343 powers using x^3+x^2+1 over F7. No F49 or F343 coordinate is used as an actual graph identity.

The full family census was repeated using the independent low-bit joint recurrence j(2n)=j(n), j(2n+1)=j(n)+n+(n mod2), rather than the author's high-bit scan. The completeness bound is valid: j≤f implies f≤(N−1)/(2^t−1), and no t≥23 can occur for N=4611841. Exactly six chart pairs occur, all order24. Their actual named rows yield exactly four distinct matrices. Therefore order24 is minimal within this declared attachment family only.

Normal and optimized audit outputs and receipts byte-match. The author's hostile controls were rerun: both H462915, cut256 responses1061775/1395279, full graph automorphism orders81/27, no transitive deletion in either graph. This rerun confirms reproducibility of those inherited controls; the new independent boundary-product audit applies to the almost-transitive targets, not to the hostile graphs.

## Clarifications resolved

The author changed the opening from potentially ambiguous 'six strong tournaments' to four distinct labelled graphs with six chart representations. The author also identified fixed-source branching uniqueness as inherited and linked the new all-length classification in `other/PROOF.md`, rather than presenting the earlier finite-scan boundary as the current frontier.

No unresolved mathematical or computational defect was found. The packet explicitly does not prove a legal native guarded transformation from the frozen source, an unrestricted minimum tournament order, or a returning two-digit grammar. Equal numerical H is not treated as a physical named-lineage continuation.

Reproduce with `python verify.py` and `python -O verify.py`. See `receipt.json`.
