# Positive insertion-cut geometry

[Proof and scope](positive_cut_geometry.md).

The exact insertion response is a positive directed-cut function. This supplies mincut optimization, compatible lower-band cuts, sharper symmetry bounds, and a divisor-cut compiler for circulant parents. These are statements about one parent and its possible added vertex; cross-parent spectrum completeness remains OPEN.

Reproduce primary and clean-room checks:

    python3 -B audit.py > primary.replay.json
    python3 -B -O audit.py > primary.optimized.replay.json
    cmp primary.replay.json audit.out
    cmp primary.optimized.replay.json audit_optimized.out
    python3 -B symmetry_audit.py > symmetry.replay.json
    python3 -B -O symmetry_audit.py > symmetry.optimized.replay.json
    cmp symmetry.replay.json symmetry_audit.out
    cmp symmetry.optimized.replay.json symmetry_audit_optimized.out

The primary audit imports probe.py for its exact kernel; the symmetry audit reconstructs the response separately and does not import the primary referee.
