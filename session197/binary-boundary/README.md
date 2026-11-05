# Binary mirror boundary and quadratic C=1 constructions

The exact physical deletion formula is H(F0)=h^2+q^2-e_start^2. The positive clone formula has no such correction, so the zero-length boundary must be handled separately. The packet also proves the cycle-square law, its two-vertex exception, and an explicit strong self-converse C=1 family with H=2(2^k+1)^2-1, including a literal14-vertex witness for2177=7*311.

Read `PROOF.md` and the independent `REVIEW.md`. Both preserve their exact reviewed bytes. Compile `verify.cpp` with a C++17 compiler and run it; no project, external library, data file or network access is required. The optimized and assertion-enabled unoptimized outputs are the two adjacent receipts. The independent standalone `audit/check.cpp` uses literal word enumeration and walk inclusion-exclusion; its output is `audit/receipt.txt`.

The frozen proof retains its original dependency names. Their exact mathematical texts are copied as `dependencies/BINARY_RIGIDITY.md` and `dependencies/ENDPOINT_RECIPROCITY.md`. They are not runtime dependencies of either standalone verifier. The new zero-length formula is proved directly by a single-edge flip, independently of extrapolating the positive-length identity.

The separately investigated converse-orbit refinement is not included in this checkpoint. `RELEASE_FILES.sha256` pins all released files except itself.
