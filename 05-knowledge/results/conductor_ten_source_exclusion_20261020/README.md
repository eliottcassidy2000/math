# Genus-five Keller sources and the conductor-twelve bound

20 October 2026. **PROVED; independently audited PASS in the stated scope.**
Prepared for publication; not yet published.

For a complex planar Keller pair with one positive-y-degree coordinate monic
and C(x)(P,Q)=C(x)(y), generic source-algebra codimension five is impossible.
Together with the audited predecessor, a nonautomorphism has codimension at
least six and generic normalization-conductor degree at least twelve.

The exact genus-five plane possibilities are <2,11> and <4,6,7>, with
terminal coordinate degrees (2,11) and (4,6). The symmetric semigroup
<5,6,7,8> is not planar. Explicit immersed realizations prevent conflating
plane admissibility with polynomial realizability. The proof retains all
scalar Jacobian coefficients and treats both branches of its first integrals.

This is not an unrestricted planar Jacobian theorem, a special-fiber-uniform
bound, or a function-field-degree bound for the planar map.

## Read and reproduce

- [Full proof](CONDUCTOR_TEN_SOURCE_EXCLUSION.md).
- [Independent review](INDEPENDENT_REVIEW.md).
- [Exact checker](independent_checks.py), [ordinary results](results.json),
  [optimized results](results-optimized.json): 88 controls, byte-identical outputs.
- [Primary sources](PRIMARY_SOURCES.md) and [dependency scope](DEPENDENCIES.md).
- [Dependency hashes](DEPENDENCY_HASHES.json), [audit summary](REVIEW_SUMMARY.json),
  [copy provenance](PROVENANCE.json), and [file hashes](SHA256SUMS.txt).

    python3 independent_checks.py
    python3 -O independent_checks.py

The [audited predecessor bound ten](../conductor_eight_source_exclusion_20261020/CONDUCTOR_TEN_SYNTHESIS.md)
is the sole earlier exclusion package needed for the final bound-twelve
synthesis. Relative links specify the intended repository layout; this
prepared packet does not assert that its files are already published.
