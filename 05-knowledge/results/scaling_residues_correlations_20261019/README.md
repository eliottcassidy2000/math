# Scaling, residues, and correlations

Read [the mathematical note](SCALING_RESIDUES_AND_CORRELATIONS.md) for the
proofs and claim boundaries, and [the independent review](REVIEW.md) for the
proof audit and complementary controls.

The package contains five connected results:

- Sharp worst-case lattice precision and an exact continuum-limit hostile.
- Exact prefix-residue transition, total-variation error, and affine wheel factors.
- Relative transfer from consecutive integer windows to actual quadratic cusp windows.
- A rough twin-center Goldbach application, using the existing short-interval sieve.
- A complete Boolean-lattice scaling-gauge criterion, with the four-antichain mantissa hostile.

The actual cusp threshold-width result comes from the
[parabolic escape theorem](../parabolic_escape_precision_20261019/parabolic_escape_precision_20261019.md).
Full-wheel digit mixing is kept distinct from power-depth divisor sieving.
No distribution assertion about the digits of pi is assumed.

## Reproduction

Run from this directory:

    python verify_scaling_correlations.py
    python -O verify_scaling_correlations.py
    python independent_checks.py
    python -O independent_checks.py

The first verifier uses only Python's standard library and passes 148,202
active exact guards. The independent verifier requires SymPy and passes
178,019 active exact guards. Ordinary and optimized runs produce identical
receipts for each verifier. The receipts record the finite universes and
script hashes.

The inverse-square and perturbed-weight tests check transfer algebra, not
the analytic threshold theorem, sieve onset, or an actual-prime assertion.
No floating-point output is used in this packet.

## Review and publication provenance

The reviewed proof and this version differ only in six hyperlink edits.
[publication_provenance.json](publication_provenance.json) records both proof
hashes and the verified link targets. The review’s mathematical content and the independent checker are unchanged;
the review’s presentation was updated to match this package. Its optional suggestion about an additional singular-series
normalization has not been adopted in this checkpoint.

`MANIFEST.sha256` covers the published files. Verify it with
`sha256sum -c MANIFEST.sha256`.
