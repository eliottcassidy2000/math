# Minimum-order spectrum framework

**Main result:** every attained h>1 has a tournament witness on at most
`ceil(22 log_{80267} h)` vertices. This is an unrelated-witness compression,
not an induced-core bound or full spectrum completeness theorem.

- [Proof and min-plus first-appearance framework](THEORY.md).
- [Relative gap amplification and a density-zero hostile](COFINALITY.md).
- [Standalone verifier](verify.py), [normal receipt](verification.json),
  [optimized receipt](verification.optimized.json), [hashes](SHA256SUMS).

The proof depends explicitly on inherited finite positive strong-tournament
banks and the exact Moon–Busch strong floor. The script independently checks
new integer bounds, packing, capped min-plus closure, actual small tournaments,
and hostile examples; it does not regenerate the order-eleven bank.

Universal realization of all positive odd values except 7 and 21 remains OPEN.
