# Pair-free cyclic-child quartics and W-prime separation

The packet proves an exact disjoint-quartic kernel for a directed-triangle module by squarefree path-matrix elimination, plus a general same-port zero identity for arbitrary modules. Ordered exterior identities then prove strictly positive marked quartics for every alternating-hub W quotient with a C3 hub and endpoint-balanced extreme children satisfying the stated positive endpoint condition. All other children may be arbitrary nonempty tournaments.

In particular the all-C3 substitutions form an unbounded fully pair-free family: every original-module XOR is visible, and all 2^(2r+2) quotient/child-orientation states have distinct full named responses. The comparison is the whole core or ordered SCC padding. Arbitrary strong outer-ambient lifting and arbitrary prime quotients are not claimed.

Normal and optimized runs agree. They check 2,508 distinct literal augmented tournaments, with complete interpolation values retained: 30 local kernel controls, 12 general-module same-port zero controls, 3 balanced quotient-root controls, 6 ordered-hub physical controls and 40 complete root-containing families at the two smallest W orders. The positive factorized formula is checked through 40 chain children and on mixed profiles. These are controls for the algebraic proof, not an unbounded inference.

The factorized verifier needs child profiles only through h=3: every kernel term has at most three exterior path slots, and each slot meets a transitive-chain child at most once. The generic exterior path-algebra checks do not use this truncation. Physical counts are bounded at seventeen names, so exact uint64 arithmetic is safe. The elementary blow/flip/attachment primitives are imported from ../arbitrary_third_root_family/verify.py, whose hash is pinned.

Reproduce with `python verify.py > normal.json`, `python -O verify.py > optimized.json`, then compare. Both independent reviews are maintained separately: named-response reviews Sections 1–3, and circuit reviews Sections 4–6.

Exact publication allowlist: PROOF.md, verify.py, kernel.py, count.cpp, README.md, normal.json, optimized.json, SHA256SUMS. Exclude explore.py, exploration.json and literal_exploration.json as superseded scratch. No executable is published.
