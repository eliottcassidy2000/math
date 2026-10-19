# Square-sum bipartite cleanup and the first K₃,₃

Read [the proof and cited bridge](SQUARE_SUM_BIPARTITE_CLEANUP.md). The affine maps a↦25a+1 and b↦25b−1 preserve cross-square adjacency and eliminate every within-part edge. The cited positive sumset theorem therefore gives induced K₃,m for every m, with varying anchors and no imported total-height estimate.

The minimum positive-label ceiling for K₃,₃ is exactly 673, both as a subgraph and induced. The displayed witness is checked directly, and two separately written scans exhaust all possibilities below 673. This finite threshold is separate from the arbitrary-m cited theorem.

Run from this directory:

    python3 -B verify_bipartite_bridge.py
    python3 -B -O verify_bipartite_bridge.py
    python3 -B review-root/check_first_k33.py
    python3 -B -O review-root/check_first_k33.py

Each normal/optimized output pair agrees exactly with the corresponding receipts. The main independent suite has 269,721 explicit checks; the root scanner supplies a separate complete lower-bound census and witness search. Neither computation alone proves arbitrary-m existence.
