# Independent targeted review of the digit-one transport

Date: 2026-11-04. Verdict: **PASS for the stated source-dependent theorem and exact update; no remaining mathematical correction identified.** This is not certification of the entire S10 theorem, its large census, or literary priority.

## Exact scope and source identities

Reviewed `../session189-cutarc-incoming/PROOF.md`, especially sections 2–4; the pinned THM4868 statement; relevant S10 section 2 proofs (insertion, canonical run, asymmetry, admissible parents and run transitions); and the earlier marked-correspondence `REPAIR.md`. Current PROOF SHA256: `d64b370e969b8c2c9608bd40d02e6d9eba60519149e4761828d965067d9b7a87`. Checked source Git blob hashes against SOURCES.json: all five match, namely S74 `54c614fe0ff4ebec0ca49783d43a6bbbeba8ff89`, THM4868 `45289f07626f507ed8b0193d20bfdf9fcabfe06c`, S10 `e7544ce7ef7b41692a824ea9d5b1b1e5987c0a5e`, prior repair `6bcead86383ed61d6e13e199ee0288d2893e2568`, and prior repair review `cac156222db2fc22f44a54085d900ffee0642b45`. The last two entries and PROOF section 1 consistently pin publication to commit `80147be9d0693a3b1b0436d3fe4ac3781d492bd6`.

S10 section 9 is indeed still `(Filled after the independent audit.)`. Its own audit-status claim must not be upgraded on this evidence. The all-order bridge depends on THM-4868(C,E,G), with (H) supporting the stated run-transition tests. No literature source was independently retrieved. No large incoming census or all-tournament census was rerun.

## Mathematical checks

1. **Canonical marking and Phi.** For a strong sole-exit tournament R, q lies in the first strong component of R-x; a spanning path ending at q places it in the last. Positive readiness therefore forces a strong deleted core. Conversely a strong reversed core supports both directions of reachability for its sole-exit attachment. The core reversal is a named involution and preserves exactly the claimed cycle/path identity. The condition on the deleted core is essential and is present. Two degree-one exits cannot both have strong deletions, because deleting the head of their mutual arc strands its tail as a sink. THM-4868(C) supplies existence in UH order at least four, hence a canonical, relabeling-compatible unmarked bijection. The triangle exception is correctly excluded.

2. **Persistent mark and induced heredity.** Each admitted p differs from x: y is a different run endpoint, y- is outside the run, and the optional y-- is critical by the sole-entry lemma. In U^p-x, the new vertex has both an in-neighbor p and another out-neighbor in the old strong core. The old canonical x therefore remains noncritical and degree one. Reversing just the core changes precisely the new edges claimed, yielding N+(w)={x,p}. Canonical UH deletion recovers a ready induced parent at the same physical x. This proves existence of a canonical ready deletion, not that every vertex deletion is ready. Uniqueness of the UH decomposition plus asymmetry distinguishes admissible pairs up to isomorphism, so no duplicate ready child classes occur. The order-four seed and order-five boundary are consistent.

3. **Exact path update.** Paths with w last contribute H-e_x-e_p. For a successor x or p, the ordered split covers contribute F_x+F_p, with exactly the forbidden prefix endpoints removed by G_(p,x)+G_(x,p). Same-endpoint exclusions are impossible because the two parts are disjoint. Empty prefixes are properly included by H(empty)=1. This identity holds for any distinct x,p with the specified attachment; admissibility is needed only to preserve readiness. The added direct derivation of F_x=H+e_q(R-x)=H+1 is correct. No closed recursive state or spectrum coverage follows, and neither is claimed.

4. **Six-vertex obstruction.** Independent permutation enumeration confirms every displayed row, the two unique cycles, H(U)=15 for each, run cardinality 3, four cut arcs, two sole entries/exits, canonical x=2 and q=3, and H(R)=33 versus 29 with digit one. Thus the enlarged unweighted state including H(U) really fails to determine H(R).

## Reproduction and small improvements completed

- Replayed the final incoming checker, SHA256 `a898c73981ec112594c2ff431c056c3cd62cef616f5b640b905b63f05a84868e`, normally and with `-O`. The two receipts in this directory are identical: 88 generated classes, 87 insertions, orders 4–8. Class counts are 1,3,8,21,55. This is a generated-family check, not independent completeness evidence.
- Added `check_witnesses.py` here. It uses direct permutations, not imported incoming checker functions, to check the witnesses and all five distinct-p two-outneighbor updates on each ready witness (including nonadmissible p). Its normal/optimized receipts `witnesses.json` and `witnesses.optimized.json` are identical.
- During review the author made the F_x derivation explicit and added named witness regression assertions to the incoming checker. Both are verified in the final files. The earlier checker already emitted those exact witnesses, so lack of explicit named assertions was a hardening opportunity, not a disproved claim.

No further correction is required for the stated scope. The structural source remains the all-order dependency; these finite tests and this transport review do not certify S10's full census or other unrelated theorems.
