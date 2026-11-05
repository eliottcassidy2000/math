# Exact converse orbits of the quadratic boundary family

The designated converse involution on a deleted binary mirror fixes exactly h+q-e Hamilton words. Every member of the quadratic C=1 family is rigid, so its class contributes exactly2^(k+1)+1 normalized self-converse presentations, while its ordinary Hamilton-word mass is2(2^k+1)^2-1. The proof keeps ordinary classes, path-pointed fibres and converse orbits distinct.

Read `PROOF.md` and its independent `REVIEW.md`. Compile and run the standalone C++17 `verify.cpp`; its optimized/debug outputs are adjacent. The independent standalone `audit/check.cpp` enumerates full Hamilton words, involution orbits and normalized graph presentations; its receipt is `audit/receipt.txt`.

The physical boundary dependency is the reviewed parent-directory `../PROOF.md`. The exact normalization theorem is copied in `dependencies/CONVERSE_MASS.md`. Frozen proofs and review retain their original scratch-era names for hash stability. Neither C++ checker requires runtime dependencies or network access. `RELEASE_FILES.sha256` pins every released file except itself.
