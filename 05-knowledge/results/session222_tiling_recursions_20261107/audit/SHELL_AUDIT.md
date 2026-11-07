# Session222 shell audit: accepted with scope

Reviewed shell REPORT.md and verify.py, independent audit.py/RESULT.json, and producer monotonicity report. Geometry and source/sink proofs are sound. The explicit permutation sends shifted Q_m plus its L onto anchored Q_m plus two diagonals. Both additions have2n−1 cells. This is a cell/Boolean-coordinate bijection, not a tournament or H-preserving map.

Independent endpoint-DP enumerated all1,098 framed tilings of orders3..6: zero L-padding preserves H in every case; anchored zero two-diagonal growth changes H in1,7,63,1023 cases respectively. The decisive minimal witness is a single backward arc: H3=3; order5 zero-leg extension H=3; zero-diagonal extension H=9. Independent geometric map checks cover n3..100.

The producer's additional monotonicity experiment reports all32,768 n7 masks:32,767 strict increases, only zero equal, minimum positive gap6;3,000 seeded nonzero n8 samples also increase. These extra ranges are author-tested, not independently rerun in this audit. They support no all-order monotonicity claim. The single-defect doubling formula is all-order proved separately by the independent assignment of interior vertices to the two increasing runs.

Active checks remain enabled under Python -O. Reproduce independent results with `python3 -O audit/audit.py`. Artifact RESULT.json states the exact universe. The historical and global bit conventions are specified in the shell report, and no recovery of absent historical images is claimed. 
