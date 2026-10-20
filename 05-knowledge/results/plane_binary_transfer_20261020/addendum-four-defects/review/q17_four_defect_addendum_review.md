# Independent review: q17 four-defect addendum

Reviewed 2026-10-20 UTC. This is a bounded verification addendum to the frozen base. No new minimum search was run, no base file was changed, and no checkout/Git writes were made.

## Verdict

PASS. The supplied certificate excludes every coloring with at most two monochromatic norm-one triangles after the valid fixed-monochromatic-root reduction. Together with the separately replayed zero-defect exclusion and even-cut parity, this certifies:

- At least four monochromatic unordered norm-one triangles: M >= 4.
- At most 1730 cut edges in the norm-one graph: Cut <= 1730.
- For every distinct three-point template T, D_T >= 1/2601, hence at least two monochromatic ordered placements among the 5202 placements.

These are lower/upper bounds, not a determination of the exact minimum or maximum cut. No real-plane bound or formal kernel verification is claimed.

## Certificate identity and replay

Certificate: results/q17_no_two_defects_certificate.json.gz

- Compressed SHA-256: 251cf994140ed97daf9fc963dab8544c3c6018c51cad8f39f095edd3087e23fe
- Compressed size: 2,477,596 bytes
- Decompressed JSON size: 11,330,989 bytes
- Proof nodes: 29,269
- Conflict leaves: 14,635
- Forced assignments: 565,949
- Maximum branch depth: 143
- Root orbit: all 1,734 reconstructed norm-one triangles

The final verifier passed both ordinary Python and python -O. An independently implemented iterative replay also reconstructed the triangles directly from norm adjacency and checked every propagation, branch, and conflict, obtaining the same counts. A separate strict type/range audit found every recorded index canonical.

The frozen base zero-defect verifier was separately replayed in ordinary Python and passed with certificate SHA-256 3980c13e87296452da52c9d0f8a6da1a5bfbbf1130e9d8a19f1b00ec0052d9e6. All 28 entries of BASE_MANIFEST.json matched their recorded byte counts and hashes.

## Why the inference rules are valid

A force cites two distinct already fully monochromatic triangles b1 and b2, and a third distinct triangle e. The vertex being forced is unassigned; the other two vertices of e are already assigned the same color. Under the hypothesis that the total monochromatic count is at most two, the remaining vertex of e must receive the opposite color. Otherwise e would be a third monochromatic triangle. Previously fixed assignments never change, so the two budget witnesses remain valid throughout propagation.

A conflict leaf explicitly supplies three distinct fully monochromatic triangles. Each branch splits on an unassigned vertex and includes both possible colors. Consequently the finite tree refutes every extension of its root assignment under the budget-two hypothesis. The argument does not infer UNSAT from a time limit or a solver status; it replays the complete local proof.

## Why the root symmetry reduction is valid

The separate base certificate rules out zero monochromatic triangles. Thus any putative coloring with at most two defects has a monochromatic triangle. The new verifier explicitly enumerates the orbit of its fixed root triangle under all translations and all norm-one multiplications and checks that it equals the entire 1,734-triangle set. These transformations are bijections preserving the hypergraph. A suitable inverse motion therefore puts a monochromatic triangle at the fixed root, and global color complementation makes its three colors zero. The new rooted refutation then applies.

The zero-defect input is essential: the rooted certificate alone would not rule out a coloring with no monochromatic triangle.

## Derivation of the strengthened bounds

Every edge in the 18-regular norm-one graph belongs to exactly two triangles. The binary three-sign identity therefore gives, for every coloring,

    M = 1734 - Cut.

Every cut has even size, because Cut = 18|S| - 2 e(S). Hence M is also even. The new certificate and root reduction give M > 2, so parity improves this to M >= 4 and Cut <= 1730.

For any nonzero side vector a, scaling the coloring transfers the same equilateral lower bound to side a. Each unordered equilateral triangle has three ordered placements, so

    D_eq(a) >= 4/1734.

The frozen, independently checked eight-point forcing inequality gives D_T >= D_eq(a)/6. Therefore

    D_T >= 4/(6*1734) = 1/2601.

Multiplying by 289*18 = 5202 gives at least two monochromatic ordered placements for every distinct template.

## Verifier hardening and rejection tests

The initially inspected version lacked canonical range/type checks. Python negative indices could alias positive edge IDs while still appearing distinct as integers, so this was reported before packaging. The final verifier now enforces exact integer types, nonnegative index bounds, canonical root data, record/list shapes, distinct budget witnesses, and three distinct conflict IDs. Its proof checks use explicit conditional exceptions and remain active under python -O.

Two separate corrupted-certificate copies were tested against the final verifier, in both ordinary and optimized mode:

1. Flip the first recorded forced assignment color: rejected with a certificate-check failure.
2. Set the first branch vertex to -1: rejected as a noncanonical certificate index.

All corruptions were confined to temporary copies. The supplied certificate and frozen base were unchanged.

Final reviewed verifier SHA-256: c30be72d9b8425c06145cd2734bd50e9a1bd54a3e1292b4c39cb9c807741d753
