# Validation and scope

## Exact finite claims

1. The C++ tournament representatives are generated internally by complete permutation orbits, not trusted from a catalogue.
2. An independent Python orbit calculation verifies disjoint coverage of all 32,768 tournament masks at `k=6` and all 2,097,152 at `k=7`, including every orbit-stabilizer identity.
3. At `k=6`, forward subset DP, direct supporting-subcube summation, and the Walsh transform agree at every 1,835,008 reduced frame.
4. At `k=7`, the ordinary-monomial subset-zeta algorithm and the Fourier/Walsh algorithm agree at every 956,301,312 reduced frame. This is a complete pointwise comparison, not a sample or moment-only test.
5. Complete labelled histograms, positive oddness, first moments, and the inherited second-moment formula agree exactly. The order 15 square sum uses unsigned 128-bit arithmetic.
6. Independent Python DP and direct enumeration of potential signed prefixes verify the two maximum witnesses, their full tournaments and anti-involutions.
7. All frame automorphisms and matching completion orbits are explicitly computed. The global maximum counts are 12 and 24 classes; toggling a fixed matching coordinate pairs their parity shores.
8. A separate reviewer recompiled and reran the computations in an isolated directory. The APPROVED report, separately written verifier, and independent results are retained with this packet.

No integer-overflow assumption is hidden in the exhaustive transforms. At `k=7` there are only `7!*2^6=322560` candidate signed prefixes for one tournament. In the scaled Walsh calculation, each prefix contributes at most 64 to the absolute sum of all its Fourier coefficients. The total absolute coefficient sum is therefore at most 20,643,840, below signed 32-bit range; every intermediate butterfly entry is bounded by that same sum. The ordinary-monomial calculation has the same bound because each prefix contributes at most 64 coefficients of absolute value 1. Histogram populations, first moments, and the square sum use explicitly wider types.

The producer Python verifiers raise explicit exceptions and remain active under `python3 -O`. The separately authored independent audit uses assertions and must be run without `-O`. C++ mismatch checks return nonzero status rather than relying on removable assertions.

## Asymptotic claim

`verify_local_design.py` enumerates the entire signed-permutation distribution on the regular four-index frame. It checks every signed path with 1, 2, 3 edges and every disjoint pair of edges, not just representative expectations. It also verifies a concrete affine 16-point K4 design and 12,397 exact rational termwise inequalities supporting the finite error-bound derivation.

The all-order proof is the nonnegative expansion in `SQRT_E_LOWER_BOUND.md`. Finite tests do not substitute for that proof. Existence of the required almost-complete four-block packings is explicitly cited to a primary block-design theorem, and its citation was checked against the publisher's statement.

The resulting statements are lower bounds:

    liminf M_k/k! >= sqrt(e),
    M_k/k! >= sqrt(e)(1-101/(8k)) for k>=4,
    M_k/k! >= sqrt(e)(1-13/(8k)) when k=4 (mod 12).

They do not establish a limiting ratio, a `sqrt(e)` upper bound, or a general formula for `M_k`.

## Search records

The records for `k=8,9,10` are nonexhaustive. Their values are independently counted, but global optimality is not claimed. Searches and restricted extensions are clearly separated from exact censuses. Repeating a large number of local searches without improvement is not used as a proof.

## Provenance

This is an add-only continuation of the approved session14 hollow-frame extrema, moment, and concentration packets. It leaves those packets unchanged. The staged package contains no remote-writing code and makes no repository mutation. Publication, if authorized by the parent, is handled separately.
