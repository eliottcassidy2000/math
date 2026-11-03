# Independent audit

**Result: PASS.** No substantive mathematical or implementation correction was required.

## Proof review

- The SCC-ordered ordinal sum gives an induced Cartesian embedding. A between-block reversal merges an interval of nontrivial strong components; a strength-destroying internal reversal splits one. Neither can produce an extra edge within the displayed family.
- The wrapper has order 5n+10 and exactly the stated separated degree ranges. A reversal touching either guard creates a degree outside the wrapper's required guard constants, without entering another valid band. This establishes inducedness, not merely injectivity.
- Uniform module substitution preserves the entire path-cover profile. The segment-cover bijection preserves outside incidences and the abstract path-incidence graph, so it preserves both acyclicity and the number of paths.
- The alternating-group sheet and sign-holonomy statements are correct with the explicitly stated convention of deleting isomorphism-preserving flips.
- The plane realizations are injective and induced unit-distance realizations. Crossing edges are permitted; no abstract-planarity conclusion is drawn.

## External theorem check

The live text of [Tomon, arXiv:2610.03517v2](https://arxiv.org/html/2610.03517v2) was checked on 2026-11-03. Theorem 1.5 / 5.2 applies to every finite graph H and supplies an induced subgraph of a Cartesian power with the requested girth and chromatic number. Lemma 5.1 applies to induced unit-distance graphs in dimension at least two. This packet starts from an induced plane spindle, so both hypotheses apply directly. Arbitrary girth remains an external existence input, not a finite certificate established by these checks.

## Reproduction

- `check.py` and `check_chain.py` passed normally and under `python3 -O`; paired stdout was byte-identical.
- The seed checks verify all 21 pairs over 720 relabellings, all exact geometric distances, two profile algorithms, 1,224 cube flips, 780 wrapper flips, and 122,880 parity tests.
- The G3 certificate checks all 45 pairs over 5,040 relabellings and agrees under two profile algorithms.
- Freshly compiling `census.cpp` reproduced the supplied order-5, order-6, and order-7 JSON files byte-for-byte.
- `search_induced.py` reproduced every positive and negative result. Normal and optimized records agree after omitting elapsed-time fields.
- `audit_chain_independent.py` uses a separate structural search: choose the shared root, enumerate every adjacent neighbor pair and common nonneighbor tip, then exhaustively join compatible diamonds. It independently confirms G2 is absent from M5, G3 from M6, and G4 from M7. Positive controls find G2 in M6 and G3 in M7. Its normal and optimized output is byte-identical. Results are in `audit_chain_independent.json`.

The absence claims remain finite order-specific results; nothing beyond order seven is inferred except the stated induced-subchain implication.
