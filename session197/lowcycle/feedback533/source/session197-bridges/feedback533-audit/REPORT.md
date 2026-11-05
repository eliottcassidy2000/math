# Independent feedback-number separation at H=533

**PASS, no correction.** Exactly six strong H533 tournaments have C<=2, and all six have C=2, order11 and feedback vertex number3. The complete strong-order reduction, fixed-port interpretation and literal changing-port witness pass.

## Frozen source binding

Under `session197-lowcycle/`:

- `FEEDBACK533.md`: `24c37bda3f05d45bb98f685f570046313ae1a464573932b145b7d7a58f83a3cb`.
- `verify_feedback533.py`: `09393e5309b479090bb6cdf0b0d8f37a7259c81f005f39f9d1a5d7b728ba7983`.
- `verify_feedback533.cpp`: `4ee6b8270f13cfa0085826ea08462a0d0160a60d8e698cc7604df8b6d24369db`.
- `feedback533_certificate.json`: `5ba699d571361c0280f176df39711d9e4ed147c8f8f9d1c2de09dd6d85f24c43`.

## All-order and scope audit

The named Busch strong-tournament bound gives 5^(n-1)<=H^3. The exact comparison 533^3=151419437<244140625=5^12 excludes all strong orders n>=13. The proof does not bound the total order of a weak tournament or make an unrequested weak-spectrum conclusion.

The complete C<=3 bank through12 was previously independently regenerated and matched in `session197-bridges/C3_COMPLETE_CENSUS_REVIEW.md` (hash `0b6d87063bd9a655925a31189962fcd00931ea03d72de31849e9b1223b7d73d3`). Its proved strong singleton-extension completeness is inherited here, not inferred from a stable finite list. Restricting that complete bank to C<=2 cannot lose any strong candidate. Orders at most3 are separately harmless because H<=3. The 12,077 count is the order4–12 bank; it intentionally excludes the triangle.

Nonnegative endpoint-DP saturation at534 is safe for comparison with533. The independent recount below uses no saturation, so it also checks this numerical implementation. Pair deletion leaves a directed triangle in every case, while a displayed triple deletion has a completely verified transitive ordering. These certify feedback number exactly3, not merely failure of one chosen feedback pair.

Deleting two permanently retained ports leaves transitive birth order. Conversely a feedback set of size at most two can be enlarged to two initial vertices at these orders; deleting more vertices preserves transitivity, and the remaining vertices can be appended in their transitive order with the required port incidences. Thus the fixed-port interpretation is exact in the stated scope. Changing ports is a different language; the supplied word reconstructs every arc of an actual H533 core. The proof does not rule out H533 at feedback number two with C>=3, and does not assert a least omitted value.

## Independent complete filter and graph checks

The independent filter reads the **independently generated and differently canonicalized** bank `session197-bridges/c3-independent-records.csv`, not the source census CSV. It computes exact uncapped start-vertex counts by removing the first vertex from a subset, rather than forwarding endpoint counts with a cap. It recounts all12,077 C1/C2 graphs and finds exactly six H533 graphs. All are order11,C2. Optimized and unoptimized assertion-enabled runs agree.

A separate Python checker matches these six graphs bijectively to the six source representatives by complete adjacency-compatible degree-constrained isomorphism search. Each has exactly one such match and one matching bijection. It independently obtains H=533 by inclusion-exclusion over unrestricted walks of length10, and C=2 by literal anchored cycle DFS. It tests all2^11 vertex subsets of each graph for transitivity using the exact tournament score criterion, giving feedback number3 in all six cases. Thus12,288 induced-subset tests independently verify the feedback claims.

All330 source pair-deletion triangles and all six upper transitive orders are checked literally. The changing-port word is replayed as constraints on the target graph's original labels, not by copying the source's generated row construction. It checks all55 arcs, the initial forward edge8->9, and final ordered ports[1,0].

## Reproduction artifacts

In this review directory:

- `filter.cpp`: `c46ac7834663d73c78ef04ad5da6f7d5f70ba123d493cd57f785cb167c8707df`.
- Both independent hit CSVs: `b15b1071197f3e31ca94f1a1fe0121c9ff5919a04a8b2ad4db1b6e0c30f42c83`.
- Both independent inventory receipts: `5acb0b097a8236e364c9e8efaad5ef6a317513861026968c694101b7a974da18`.
- `check.py`: `683a621b2a63723edf6bb7c888271a30c15e602494d006cf77e9a53d4736c0fb`.
- Normal/optimized `controls.json`: `81f071dfaba257db3cf65fa19a21fd0756c294e4530bfe088d3fda32ea202a00`.

The controls receipt records the exact independent bank hash, all six independent representatives, their isomorphisms and feedback witnesses. The source verifier and all named dependencies were copied into `source/` before replay because the wrapper writes its certificate. Normal and optimized source runs reproduce the frozen certificate byte-for-byte. No source file was changed.
