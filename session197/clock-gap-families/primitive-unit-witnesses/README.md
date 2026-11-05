# Primitive unit ternary sources

Two explicit strong, rigid, module-indecomposable sources have clone families 260·3^m−1 and 520·3^m+1. The packet includes actual prime counts 2339 and 7019, and 7p counts 1561=7·223 and 189539=7·27077. These are finite witnesses, not an unbounded primality or spectrum-coverage theorem. Module indecomposability applies only to the original sources, not their nontrivial clone expansions.

- `PROOF.md`: construction, exact all-order run identity, structural certificates and scope.
- `verify.py`: standalone standard-library verifier. Run normally and with Python `-O`; redirect output to compare with the included receipts.
- `audit/REPORT.md`: independent PASS review bound to source hashes.
- `audit/check.cpp`: standalone C++17 verifier using literal exterior permutations, walk inclusion-exclusion, module checks, separate automorphism backtracking and a prime sieve. Compile and run with no arguments.
- `audit/INDEPENDENT.json` and `audit/INDEPENDENT_DEBUG.json`: matching independent receipts.

There are no executable file dependencies. The proof's earlier four-seed theorem is corroborating context, already published at [the reviewed clock-gap family checkpoint](https://github.com/eliottcassidy2000/math/blob/ef10dc937a51f36154189a26965e8624e4895eb8/session197/clock-gap-families/README.md). This packet's own explicit run counts and matching bound establish its two all-order identities without importing that verifier.

All files are additive. No navigation index is changed by this checkpoint.
