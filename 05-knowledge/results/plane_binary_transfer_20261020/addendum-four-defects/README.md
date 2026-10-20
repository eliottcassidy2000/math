# Addendum: at least four monochromatic norm triangles at q17

**Status:** FINITE-EXACT, independently replayed with optimization-safe verifier gates, plus elementary counting. This is a separate addendum to the immutable reviewed base packet. It strengthens the base's bound of two; it does not alter that base or determine the exact minimum. No further minimum search was run.

## The strengthened statement

Let `G` be the norm-one graph on `F_17[sqrt(3)]` and let `M(s)` be its number of monochromatic unordered graph triangles in a binary vertex coloring. Let `C(s)` be its number of cut edges. Then

    M(s) >= 4,             C(s) <= 1730.                     (A1)

The radial positive-definite kernel in the base has `rho=-1/3` on all 2601 graph edges and hence predicts cut expectation 1734. It violates the valid global Boolean cut inequality by four:

    Σ_(unit edges {i,j}) (1-rho_ij)/2 <= 1730.                (A2)

For every distinct three-point template `T` in this finite norm plane, the averaged eight-point forcing inequality from the base now gives

    D_T >= D_eq(a)/6 >= 4/(1734*6) = 1/2601.                 (A3)

There are 5202 ordered translation-rotation placements, so every such template has at least two monochromatic ordered placements. These remain finite-field statements. They imply neither a new real-plane chromatic lower bound nor a faithful real-plane realization of G.

## Exact root reduction and proof rule

The base independently certifies that `M=0` is impossible. If a coloring had `M<=2`, it would therefore have a monochromatic triangle. The translation/norm-one-rotation group acts transitively on the 1734 unordered triangles. Global color reversal and a group element would make the fixed triangle `[0,17,155]` monochromatic in color zero. The verifier independently enumerates that orbit and checks it equals the entire norm-triangle set.

The new proof tree rules out every coloring of this fixed-root form with at most two monochromatic triangles:

- A branch assigns both possible colors to one unassigned vertex.
- A forced assignment names two distinct already monochromatic triangles. Once that budget is exhausted, a different triangle with two same-colored assigned vertices forces its remaining vertex to the opposite color.
- A conflict leaf names three distinct fully monochromatic triangles.

Thus no positive-defect coloring has `M<=2`. The base's identity `M=1734-C` and degree18 cut parity show that M is even. Therefore `M>=4`, proving A1. Expectations give A2; the independently checked eight-point comparison gives A3.

No optimization algorithm's UNSAT label is trusted in place of the tree: every force, branch, conflict and root symmetry is replayed separately from the generator.

## Certificate and replay

The compressed certificate is `results/q17_no_two_defects_certificate.json.gz`:

- 29,269 proof nodes
- 14,635 conflict leaves
- 565,949 forced assignments
- 11,330,989 uncompressed JSON bytes
- 2,477,596 compressed bytes
- SHA-256: `251cf994140ed97daf9fc963dab8544c3c6018c51cad8f39f095edd3087e23fe`

Run from this addendum's directory:

    python3 scripts/verify_q17_two_defects.py
    python3 -O scripts/verify_q17_two_defects.py
    python3 review/independent_iterative_replay.py

The base zero-defect certificate must also have been replayed:

    python3 ../scripts/verify_q17.py

The new verifier independently constructs the graph by exact modular norm adjacency and checks the triangle orbit. It has explicit exception gates, not removable assertion statements, with strict integer types and bounds on vertex/color/edge references. Negative Python indices cannot alias positive edge IDs. It accepts an optional input certificate path for hostile tests. Independent normal and optimized runs passed; a deliberately corrupted forced color was rejected in both modes. The final range-hardened verifier was also rechecked. The independent iterative replay in `review/` imports neither the generator nor the primary verifier and reconstructs all triangles itself; its retained output agrees on every proof count. See `review/q17_four_defect_addendum_review.md` for the review and negative-index corruption test.

To regenerate the proof tree:

    python3 scripts/certify_q17_two_defects.py

The copied field/template helper is the same small generator used in the base. All runtime dependencies are in the Python standard library. The gzip representation is deterministic (`mtime=0`). This is a finite exact external certificate, not a Lean proof.

## Remaining boundary

The sharp minimum may exceed four. No equality coloring is asserted. The fixed global separating inequality is explicit, but a small-support or short elementary certificate would still be valuable. This addendum ends the requested bounded refinement; it does not start an unrestricted minimum search.
