# Vertex-transitive extreme-child addendum

Proved: every module of a finite vertex-transitive tournament is strong and vertex-transitive; original-module XOR families preserve all complete ordered-cover and one-marked endpoint profiles. This extends the W/C3-hub root-family separator to arbitrary vertex-transitive extreme children and arbitrary nonempty middle children, using at most m+20 fixed queries. It does not claim full orbit injectivity for arbitrary children or an arbitrary strong outer-ambient lift.

Run `python verify.py > normal.json` and `python -O verify.py > optimized.json`, then compare. A standard C++17 compiler builds the literal counter in a temporary directory. Nothing is installed.

Controls include every labelled tournament through order five for the regular-module checks; automorphism witnesses for each small vertex-transitive case; all sixteen module-XOR families of the nine-name nested child C3[C3,C3,C3], with full physical marked profiles retained; eight arbitrary ambient balanced-module controls; and128 full root-containing families. The nested-extreme case uses this nine-name child, and the crossing-middle case uses overlapping TT3 intervals. All actual physical query values are retained in both receipts.

Exact publication allowlist: PROOF.md, README.md, verify.py, count.cpp, normal.json, optimized.json, SHA256SUMS. Dependencies are the four mathematical sources hash-bound in each receipt. The verifier itself imports no sibling code. Independent review passed in the separate vertex_transitive_extremes_audit packet, including complete marked-profile enumeration, all restoration controls and literal count replay.
