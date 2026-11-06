# Independent balanced-atom and oracle audit

PASS, 2026-11-06. Reconstructed the n13 code independently, checked strong connectivity by all-source reachability, and counted paths using a recursive suffix recurrence. All13 start-count values sum to759375=15^5 and agree with the main verifier. Dominance joins give the claimed rate upper13/5, with no global minimum or optimality assertion.

The D_uv identity is correct for both pair orientations: split a path at its consecutive u,v positions. Its prefix/suffix vertex partition is unique, and neither factor contains the reversed arc. Hence the difference is D_vu-D_uv when u->v. The normalization sum over existing arcs equals(n−1)H by counting path/arc incidences. Independently checked10650 reversals across all labelled tournaments of orders2..5; each induced-subset endpoint count was computed by direct permutations, and flipped totals by full-path permutation enumeration. All passed.

Complexity bounds are sound: both endpoint tables cost O(n²2^n) arithmetic operations, all pair convolution sums cost O(n²2^n), and tables store O(n2^n) integers. Their values and convolution totals count at most n! vertex words, giving O(n log n) bit length. The proof properly distinguishes arithmetic from bit complexity and single-arc search moves from native path-complement adjacency. Search is bounded, nonexhaustive, and not a convergence certificate. Pending n14 discovery is not part of this checkpoint. No main preprint theorem is assumed.

Independent scripts/outputs: independent_witness.py and independent_witness.out; independent_delta.py and independent_delta.out. Standalone copies have been renamed without altering mathematics. Main verify.py reviewed: both endpoint recurrences and normalization computation are correct.
