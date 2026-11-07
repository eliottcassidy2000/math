# Session 210 permutohedron carrier: independent review

Status: PASS. Reviewed the final frozen proof, source manifest, verifier, and outputs. Both normal and optimized isolated replays exactly match RESULT.json, and supplied SHA-256 checks pass.

The coordinate convention is consistent: positions of labels are coordinates, so an adjacent swap of labels a,b has vector direction ±(e_a−e_b). Tournament preference gives a coherent direction, and all-pairs feedback energy strictly decreases. Sinks are exactly directed Hamiltonian words. Face product and incidence formulas correctly use forward ADJACENCIES, not all forward pairs. The alternating identity and its even-order signed corollary have the correct factors.

The convex-linear/transitive equivalence is correct. Dropping coherence permits precisely every sink count from 1 to n!/2 for n>=2: a perfect matching bounds independent sinks, and distance-to-a-chosen-parity-subset realizes each bound. The n=1 boundary is separate.

For the extreme-reversal family, the increasing-run path classification is complete, energy is m+inv(A,B), and the Gaussian-binomial recurrence has the correct exponent. The quadratic local/global gap is exact; the weighted fixed-order unbounded gap remains a local sink because every pair penalty stays positive.

The field map is a genuine bijection, not a numerical match: reverse coordinate order, use A as pivot columns, and place independent field entries precisely at inversions a>b with a in A and b outside A. These are all and only the free RREF slots, so every subspace appears exactly once. The scalar family formula and classical Schubert-cell machinery are appropriately attributed, and decorations versus bare paths are sharply distinguished.

`independent_permutohedron_checks.py` independently checks every labelled tournament through order five, every ordered-partition face, adjacent-swap gradients, the first nonglobal local sink, extreme-reversal energy distributions through order eight, and actual finite-field row spans for q=2,3 and m<=4. It imports no worker code. `independent_cover_bridge.py` separately verifies the link between this incidence distribution and the gate packet's cover polynomial.

One draft wording issue was corrected: an arbitrary set-theoretic bijection between two countably infinite spectra is not an open problem. The frozen conclusion correctly states the unresolved surjectivity of H↦(H−1)/2 onto the nonnegative integers other than 3 and 10. No completeness or literature-priority claim is established here.

Approved proof SHA-256: 471a51a30e91c57282f8bb7b5748b00c871eab3d556c7b5b63158d1419eb452f.
