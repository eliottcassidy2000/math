# Final bounded incoming-center review

2026-10-21 UTC. Verdict: **PASS**, with exact finite scopes below. This review ran only existing immutable certificate consumers, not any discovery or census producer. It made no Git or remote writes. Exact file bytes are pinned in `source_pins.json`; the parent-reported incoming commit is provenance, not an independently checked Git identity.

## Completeness and structural arguments

1. The orbit-run strengthening is valid. Full automorphism orbits are unions of rotation shores. An anti-automorphism normalizes the full automorphism group and permutes its vertex orbits with cardinalities preserved, without needing to preserve individual rotation shores. A class of one shore gives a four-vertex run; a class of two through four shores gives a five-vertex run by taking the terminal in another shore of that same class. The inherited half-supported four-cycle lemma supplies four distinct black orbits including mates; proper-shore deletion supplies the strong complement. Entry cases can be handled in the opposite tournament and reversed back. Both original and mate runs are absent from all marked-cycle cuts. Adding the cut orbits gives at least eight black orbits. Hence, for the inherited order-25 and order-35 parameter rows, vertex transitivity is necessary.

2. The prime-square group lemma is sound. A Sylow p-subgroup containing a Sylow subgroup of the point stabilizer has index p squared over its stabilizer and is transitive. A central element of order p is semiregular by faithfulness and transitivity. Its p blocks carry a transitive p-group of order p. An element inducing their p-cycle has permutation order p or p squared: no p-power cycle on p squared points can have greater length. In the latter case it is a regular p-squared cycle. In the former, it and the central element generate an independent commuting C_p x C_p, transitive and regular. Thus every vertex-transitive order-25 tournament occurs among the two abelian Cayley group presentations. Twelve inverse pairs per presentation give exactly 8,192 necessary labelled frames. No degree-35 regular-subgroup conclusion is drawn.

3. At order twelve, the exact cut index satisfies q <= 8 and m is an odd divisor of twelve. Therefore m=3 and q=4, and rotation by four is an actual automorphism in marked-cycle coordinates. Its 22 unordered-edge orbits each have size three; four prescribed cycle orbits leave 18 independent bits. All candidates occur in these 262,144 labelled frames; source duplication up to isomorphism is harmless. This is a necessary-frame universe, not an isomorphism census.

4. Order fourteen independently has m=7,q=2. Thirteen edge orbits minus two forced cycle orbits leave eleven bits and 2,048 frames. The consumer's two cyclic shores plus cross-difference construction has exactly this universe. This finite exclusion agrees with the parent's earlier analytic order-fourteen exclusion and does not remove another order from that already-reduced frontier.

## Literal consumer audit

Order25: all 8,192 sequential group/mask addresses are required, and each has seven exact Hamilton permutations. Complete ordered-adjacency keys are equal if and only if the unique word-to-word bijection is a full source automorphism. Distinct keys therefore prove seven distinct orbits. The positional reversal is tested on every vertex pair and must fail the anti-map test for every word. Thus the bank proves Delta >= 7, sufficient to exclude six without requiring parity; self-converseness upgrades the bound to eight. Inversion is independently checked as an anti-automorphism for each source. Discovery-node metadata does not establish any exclusion.

Order12: the bank decoder enforces bounded decompression, exact size/header/hash, no compressed trailing data, and a source-ordered scan with exact final offset. Six bytes encode twelve four-bit vertex labels; the exact permutation and native-edge checks exclude invalid or repeated labels. Every histogram failure is only used in the valid necessary direction: self-converseness implies complementary score multiplicities. Passing that test does not establish self-converseness. Every remaining frame receives eight native words with distinct score strings, each failing the necessary reverse-complement condition for blue words. Score inequality is used only to prove orbit inequality, never conversely. Targets are deterministically reconstructed and the complete expanded stream is committed against the producer stream. This is additional source/port binding, not a substitute for the actual word checks. Public APIs reject aliases, wrong sources, invalid labels, duplicate score words and incorrect targets.

Order14: the complete source-indexed lossless bank supplies eight words and targets per frame. Each is checked by the frozen order12 generic score-port API, with an additional direct reversal test on every word. The imported module does not run the order12 complete replay during the order14 check.

The description “no search” applies to exclusion-bank consumption. The checkers do include tiny exhaustive positive controls (orders at most seven); these are not discovery searches for the order12/14/25 exclusions.

## Independent replay results

All six commands used `--verify`; none used `--write` or a producer. Normal and optimized results agree exactly with saved receipts:

- Order12: 5,357,643 checks; LF receipt SHA256 `fe54716a7b6e6fac30a41e118e9a59e2e80cdea0da6e7602932cfa0697779da2`; 262,144 frames partitioned into 111,104 histogram rejections and 151,040 eight-word certificates (1,208,320 words).
- Order25: 163,904 body checks; LF receipt SHA256 `1cd59dbfdc5ff5c296fc69aa7f5327a41991e237ba2c291d991697e2792ad0c8`; 8,192 frames and 57,344 black words. The 8,176 distinct labelled adjacency matrices are not isomorphism classes.
- Order14: 24,601 checks; LF receipt SHA256 `2b150fd9d93ea708cd201aef0eaa41fa5d5205daed5b1dba41757cf0da3c329d`; 2,048 frames and 16,384 words.

The six logs retain exact stdout and shell timing. The complete run took roughly nineteen seconds wall time in parallel. An initial attempt to use `/usr/bin/time` failed before Python ran because that binary was absent; the successful commands used the shell's `time` keyword instead.

## Promotion and limits

The incoming order12 and order25 exclusions are eligible for promotion. Starting from the parent's nine-order frontier in which fourteen was already analytically excluded, they leave exactly **seven** unresolved orders: **15,18,20,21,24,28,35**. Equivalently, remove 12,14,25 from the frozen ten-order list.

The degree bound remains 35. The inherited minimal eventual recurrence order remains 18 and the uniform safe onset remains 48. No remaining coefficient, exact deficit histogram, ordinary center census, order35 Cayley completeness, or metagraph wiring is established here. The earlier all-order inputs remain inherited at their stated scope; this final bounded review does not rerun their producers or re-credit missing historical computations.
