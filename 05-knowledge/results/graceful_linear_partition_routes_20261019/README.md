# Linear routes for staircase partition systems

Read [the proof](LINEAR_PARTITION_DIAMETER.md) and [independent audit](review-algebra/REVIEW.md). In the stated capacity-one partition subclass, width-two diameter is at most 2r−3 for r≥3 and at most one for r=1,2. The explicit event-stack router uses O(N+r²) work in the given partition representation, retaining actual parallel elements and loops.

For graceful upper-tail fibers this applies when the entire simplified contact graph is a tree. The finite arithmetic census did not find a new nonstar or unequal-support arithmetic example; the strict subclass extension is demonstrated for abstract partition systems. No general-matroid linear bound, optimality, or new arithmetic realization is asserted.

Run each script normally and with optimization from this directory:

    python3 -B verify_partition_routes.py
    python3 -B -O verify_partition_routes.py
    python3 -B arithmetic_partition_controls.py
    python3 -B -O arithmetic_partition_controls.py
    python3 -B deep_relay_control.py
    python3 -B -O deep_relay_control.py
    python3 -B review-algebra/independent_partition_audit.py
    python3 -B -O review-algebra/independent_partition_audit.py

All four normal/optimized receipt pairs match. Main route controls have 2,737,046 checks; the independently written exhaustive audit has 1,880,058 checks. The rank-1,100 event-stack stress emits 2,197 moves beyond the recursion ceiling, while its true width-two distance is one: this is an implementation stress, not a lower bound.
