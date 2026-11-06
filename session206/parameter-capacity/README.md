# Session 206: exact O/Z tail parameter capacity

Status: independently audited all-order theorem. This package is standalone evidence for the parameter-count asymptotic and its finite-bank necessary condition, not a finite-bank covering construction.

For a fixed physical edge-grammar source s=(h,R,e,d), with R>d>0, let A_s(X) count encoded O/Z tail inputs m whose scalar output is at most X, with multiplicity. The theorem is

A_s(X) = X log(R/(R-d))/d + O_s(sqrt(X)).

Thus a fixed finite bank covering all sufficiently large odd scalars must have the sum of these coefficients at least 1/2. A strict sum below 1/2 excludes that particular bank. No conclusion follows at equality, and the root coefficient log(2) exceeds 1/2.

## Files

- PARAMETER_CAPACITY.md: complete mathematical proof, unchanged from the independently reviewed draft.
- verify_capacity.py and RESULT.json: reproducible numerical controls for 55 physically generated source states; 10,000 exact digit-identity checks, plus input counts at three finite thresholds. Counts are inputs, not distinct scalar outputs.
- audit/REVIEW.md: frozen independent mathematical review, including a separate block-decomposition derivation.
- audit/AUTHOR_REPLAY.json: independently replayed author verifier output.
- audit/verify_block_identity.py and audit/BLOCK_RESULT.json: independent exact finite controls for the reviewer's block identity.
- SHA256SUMS: complete payload checksums, excluding this manifest itself.

## Reproduction

From the package directory, with Python 3:

    python verify_capacity.py > /tmp/capacity-replay.json
    cmp RESULT.json /tmp/capacity-replay.json
    python audit/verify_block_identity.py > /tmp/capacity-block-replay.json
    cmp audit/BLOCK_RESULT.json /tmp/capacity-block-replay.json
    sha256sum -c SHA256SUMS

The asymptotic is proved symbolically; no finite computation certifies an asymptotic error bound. The exact script checks are regression controls only. The displayed floating-point logarithms are illustrative finite comparisons, not inputs to the theorem.

## Scope

This does not establish a distinct-output density, eventual coverage, universal finite-bank obstruction, or an asymptotic uniform over banks growing with X. Physical seed attainability is essential to use the arithmetic map as a tournament construction. The prior exact tail formula and positive-width inverse bound are recalled in the proof; the new proof of the capacity coefficient uses only elementary digit identities and residue-class counting.
