# Physical clock families and exact gap cancellation

Two independently reviewed all-order packets:

- [Four primitive unit-clock families](two-module-unit-clocks/PROOF.md), with [PASS review](two-module-unit-clocks/REVIEW.md). Primitive amplitudes60,100,156,260 support explicit binary/ternary module updates. Actual prime601/937 and7p721/3241 graphs are included.
- [Arbitrary-gap formula and converse](transitive-crossword-gaps/PROOF.md), with [PASS review](transitive-crossword-gaps/REVIEW.md). In the physical nonzero skeleton1122, the binary coefficient vanishes exactly when both internal gaps have size0 or1. The144-graph tensor certificate proves this for unbounded gap lengths, including deletion boundaries.

Each subdirectory is self-contained. Compile and run `verify.cpp` with C++17; optimized/debug receipts are included. The independent `independent/independent.cpp` uses walk inclusion-exclusion rather than either source counting method, and has its own two receipts. No external source code or graph bank is needed. Frozen proof/verifier/review bytes are preserved.

Original scratch-style paths in the proofs identify the same session companions. The four-seed proof cited by the arbitrary-gap proof is the adjacent subdirectory here. The earlier primitive60 bridge is provenance, not a dependency of either standalone tensor proof.

Scope: these are constructive families and an exact converse inside one specified word skeleton. They do not classify every cross word or unit clock, prove a global mod8 rule, establish full-spectrum coverage, or claim infinitely many prime values.
