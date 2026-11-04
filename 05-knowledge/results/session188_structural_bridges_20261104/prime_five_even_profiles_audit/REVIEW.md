# Independent audit: prime-five even-child local responses

PASS. No mathematical correction is needed. This review is bound to the frozen source pins retained in `normal.json`.

## Analytic review

The normalization domain contains every actual even child of order at least four; it is deliberately an upper domain, not a realizability assertion. Expanding around the actual TT4 profile is legitimate because all perturbations are multiples of eight. In each first-order Hamilton term all other children occur once. The first two intermediate return segments for three or five occurrences of the marked child are nonempty, disjoint in label support, and have the same child endpoints. Their exchange is fixed-point-free even when it moves a port, and preserves every join. Thus the q3/q5 derivative parities vanish in each response separately. This argument is Hamilton-specific and is not silently extended to ordered covers.

The full affine formula preserves the common odd product of child Hamilton counts. The displayed Q10 and Q41 gaps cannot vanish modulo sixteen. For Q76 the two gaps are complementary zero/eight values. The complete labeled prime classification has three disjoint orbits of sizes 120, 120, 24, so no five-name prime quotient is omitted. Relabeling expands each chosen bit onto the original physical child set and does not collapse physical names. Profile-preserving replacement is valid because all queries are uniform on each child.

This proves local quotient-converse separation only. It is neither unrestricted reconstruction nor pairwise orbit injectivity; an ambient claim requires the separate all-cover analysis.

## Independent algorithm

The audit imports no producer code. `interpolate.cpp` counts Hamilton paths by forward subset/endpoint dynamic programming on actual transitive-child inflations. Exact integer triangular inversion of the Stirling profile matrix recovers every run multiplicity coefficient needed for the baseline and all fifteen first-order channels. Divisibility and nonnegativity are checked at each inversion. All physical graphs have at most eleven vertices, so unsigned 32-bit path totals are bounded by 11! and cannot overflow; interpolation numerators use signed 64-bit arithmetic.

There are 56 baseline multiplicity vectors of total valuation at most three and five additional one-TT5 controls, giving 61 actual physical size assignments for each query/orientation. The complete run bank is therefore independently recovered from 374,784 physical counts. Every one of the 6,144 full affine response vectors is compared directly, not through hashes.

The independent Python checker additionally:

- Enumerates all 1,024 labeled five-name tournaments, tests every proper nonsingleton subset for the module property, and canonicalizes all 264 prime cases over every permutation.
- Performs a complete ten-variable Boolean transform of every response/channel, checks forbidden higher bidegrees vanish, compares all 256 retained coordinates including cubic/quartic entries and zeros, and reconstructs all 1,024 full queries.
- Reconstructs and retains all 32 complete contrast vectors for each quotient and verifies every explicit battery formula. State digests are supplementary, never the equality test.
- Reconstructs the exact twenty-name old graphs and physical port masks for all eight retained twenty-two-name controls, and independently recovers their residues via the small-physical interpolation. This last check is not claimed to be a second direct twenty-two-name enumeration.

Normal and optimized Python runs complete successfully and produce byte-identical receipts. Both rerun the physical interpolation rather than trusting the earlier exploratory text dump.

## Reproduction and allowlist

From this directory run `python audit.py > normal.json` and `python -O audit.py > optimized.json`, then compare the files. A C++17 compiler is required; no package installation is used. Source and inherited theorem pins are checked before computation.

Publish exactly: `REVIEW.md`, `audit.py`, `interpolate.cpp`, `normal.json`, `optimized.json`, `SHA256SUMS`. The exploratory `interpolation.txt` and `interpolation.stderr` are unnecessary and excluded. Full lossless bank data are already retained in both receipts.
