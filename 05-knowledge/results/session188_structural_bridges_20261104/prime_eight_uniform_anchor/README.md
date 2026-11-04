# Prime-eight uniform anchor

For the explicit prime quotient `(238,4,240,118,195,146,34,74)`, the fixed block query `(1,0)` has full augmented cover contrast `(4,0,0,...)` modulo eight for arbitrary nontrivial child tournaments. The same query lifts through every fixed ambient. A converse-complement pair detects every original-module family containing this root oddly. This is a single-quotient theorem, not general prime reconstruction.

Run `python3 verify.py > normal.json` and `python3 -O verify.py > optimized.json`, then compare. The complete certificate is two queries times three cover coordinates times 129 source/target run vectors, contracted coefficientwise over all sixteen normalized profile parameters. Literal controls independently count complete augmented tournaments and retain shared cubic and quartic finite differences. See `PROOF.md` for exact cut conventions, normalization, ambient lift and boundaries.

Independent review passed. The separate `prime_eight_uniform_anchor_audit` reconstructs every kernel from 516 literal physical cover profiles and independently checks the larger named controls. Mathematical sources and paired receipts are frozen.

Frozen source allowlist: `PROOF.md`, `README.md`, `verify.py`, `covers.cpp`, `physical.cpp`, `normal.json`, `optimized.json`, `SHA256SUMS`. No exploratory files or executables are included. Source receipt dependencies are sibling proof packets; their hashes are retained in each receipt.
