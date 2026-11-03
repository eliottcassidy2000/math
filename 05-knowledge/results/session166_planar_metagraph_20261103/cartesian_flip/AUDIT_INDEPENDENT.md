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

The absence claims remain finite order-specific results. No nonexistence beyond order seven is inferred except the stated induced-subchain implication.


## Additional order-eight witness

The later G5 witness also passes independent audit. Both normal and optimized executions of `check_chain8.py` pass with identical output. The separate `audit_chain8_independent.py` does not import the author's graph or profile helpers. It computes the permutation action directly on bits, reproduces all 120 minimum quotient distances over 40,320 permutations per source, confirms all sixteen tournaments are strong, and checks every complete profile using endpoint subset DP followed by a first-block recurrence.

The induced graph has exactly 29 edges on 16 classes. Removing O,t1,t3 leaves the five components {t0,a0,b0}, {a1,b1}, {t2,a2,b2}, {a3,b3}, {t4,a4,b4}. This certifies nontraceability, since deleting three vertices from a Hamiltonian path can leave at most four nonempty components. An explicit proper four-coloring passes; the diamond-forced three-color contradiction establishes chromatic number four.

Nontraceability belongs to this metagraph subgraph, not its tournament vertices. Their directed Hamilton-path counts in certificate order are 541,427,497,517,415,497,517,401,509,479,415,497,517,427,497,517. All are positive odd integers. The prior exclusions establish that order eight is minimal for both G4 and G5. No order-eight census or negative search is needed to certify this positive witness.

Independent code and output: `audit_chain8_independent.py` and `audit_chain8_independent.json`. Normal and optimized output is byte-identical.
