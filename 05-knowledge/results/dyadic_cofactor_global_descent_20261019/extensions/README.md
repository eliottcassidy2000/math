# Transient heights and inverse basins of the dyadic cofactor map

These two independently audited extensions build on the [global cofactor classification](../COFACTOR_MAP_GLOBAL_CLASSIFICATION.md), preserving its prime-domain restriction.

- [Transient height](TRANSIENT_HEIGHT_SCALING.md): for the maximum M(n) along a complete auxiliary-map orbit, limsup (M(n)-3n)/n^(2/3)=3. The resonant starts n=d²(d-1), for dyadic d>=8, have the exact maximum 3d³-2. The asymptotic upper bound comes from a universal excursion envelope.
- [Inverse basins](INVERSE_BASIN_DIVISOR_LATTICES.md): the factorization F(cp)-1=(c-1)(p+2) gives exact divisor-index descriptions for preimages of fixed points and powers of two, with the native odd-prime restrictions retained.

These are statements about the auxiliary cofactor map. No general aliquot height or inverse-basin conclusion is inferred. The original classification proof and its original review remain unchanged. Historical priority is not established.

## Reproduce

All three programs use only Python 3 and its standard library. From this directory, run each normally and under `python -O`:

- `python -B verify_inverse_basins.py`: compare stdout with `verify_inverse_basins.out` or `verify_inverse_basins.optimized.out` (3,100 author cases).
- `python -B review/audit_transient_height.py`: compare stdout and generated JSON with the transient-height review receipts. Includes two million local excursion starts, 500,011 cached full-orbit nodes, and 254 exact resonant maxima.
- `python -B review/audit_inverse_basins.py`: compare stdout and generated JSON with the inverse-basin review receipts. Includes 50,000 complete inverse targets and the fixed/power divisor lattices at large exponents.

Both [height](review/TRANSIENT_HEIGHT_INDEPENDENT_AUDIT.md) and [inverse-basin](review/INVERSE_BASINS_INDEPENDENT_AUDIT.md) audits pass. The checkers use explicit failures that remain active under optimization. Finite controls supplement the universal proofs. `SOURCE_PROVENANCE.json` records reviewed and published hashes, the status-only proof updates, presentation cleanup, and package-relative checker paths.
