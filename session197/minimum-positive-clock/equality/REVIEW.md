# Independent review: unrestricted positive-amplitude 72 equality

**PASS. No correction requested.** This review applies only to the exact frozen files below, under `session197-prime-atoms/positive-equality-global/`.

- `PROOF.md`: `02363cc76958663bae028df5bfc0d15f03c138866cd5881d3cbebb2d9b590934`.
- `verify.py`: `0a7068f130c5db9c3b889086d8fac7c752ae00021df20d491ef8506ede2b0082`.
- `verify.cpp`: `8c2bc8a9cabd6194106f53d79b5063b124d4d53e72254272a3aeead7355fac30`.
- `receipt-normal.json`: `107eaca0607f7d4fa41bfa7396f74d2306fe4278a23943b124cfce741c08d877`.
- `strong-bank.jsonl`: `3d4ddbe02807c7c6bf962147c582e0c6707eb4d221134a715cdd8a84bf730951`.

## Mathematical audit

The run expansion and exact matching-degree theorem reduce an eventual pure positive-unit ternary identity to the degree-three run coefficients. Those dependencies were independently checked in my preceding minimum-amplitude audits. The SCC product outside the root core divides both 217 and 649, whose gcd is 1, so every exterior SCC is a singleton. This permits independent prefix and suffix lengths, not merely symmetric padding.

The strong-order cutoff is correct and uses the necessary sharp boundary: the uniform inequality rules out order 12 and above since 5^11>217^3, while the exact order 11 Busch minimum 225 separately rules out order 11. The uniform inequality alone would not suffice. Thus the complete strong bank through 10 closes the all-order reduction.

The strong-deletion generation is the same exact procedure whose elementary maximal-strong-subset proof I audited in the preceding minimum packets. Every mixed extension of a strong tournament is strong, all eligible predecessor cores satisfy the same induced-path cap, and the degree-class canonicalizer enumerates every permitted permutation exactly. The C++ early cap rejects only after an induced subset already has too many Hamilton paths, so it is safe. Arrays and integer bounds are sufficient for all source and clone orders used.

The two finite clone equations, at matching rank 2, are equivalent to a2=3a3 and a1=3a3+1. The complete finite result leaves one strong rooted core. An isomorphism preserves the ordered SCC list, hence the arbitrary padding lengths. The core is rigid and its clock root is unique. Therefore the n−8 unlabelled and (n−8)n! labelled formulas, and the prescribed-root formula (n−8)(n−1)!, follow exactly.

## Reproduction and independent checks

- Source normal and Python-optimized runs both reproduce receipt hash 107eaca0607f7d4fa41bfa7396f74d2306fe4278a23943b124cfce741c08d877 byte-for-byte.
- My independently written exterior-gap polynomial DP, which has no explicit revisited-root state, checks all 2,456 retained graphs and 2,817 candidate roots. It gives 31 first-equation survivors and exactly one complete hit, order 9 root 4, with full run coefficients (73,217,216,72). Optimized and assertion-enabled unoptimized C++ agree.
- A separate adjacency-compatible backtracking isomorphism check identifies that hit with the sole positive-A72 class in my independently reviewed complete rooted catalog through order 9. The map from this core to that catalog is [4,5,6,0,8,7,1,2,3]. It also independently finds only the identity unmarked automorphism.
- The supplied Python checks independently recount every bank graph, every candidate clone, the full labelled universe through 6, all roots of the equality core, and ten unequal-prefix/suffix padding controls. I reviewed these checks and their assertions.

Independent source and receipts are in `audits/positive-equality/`: `gap_check.cpp`, `gap.json`, `check_catalog.py`, `catalog-check.json`. The frozen source was not altered.

The conclusion classifies only fixed-root transitive-clone sources at amplitude 72. It does not claim any unrestricted scalar-spectrum or arbitrary graph-changing-clock classification.
