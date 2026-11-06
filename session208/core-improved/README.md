# Improved explicit native tree distributions

Independently audited feasible upper bounds at tournament order4:

- Expected stretch at most267/83.
- Expected fundamental-cut congestion at most1063/330.

The two distributions use the same two seed-tree symmetry orbits with different rational weights. Each uses at most96 trees. These improve the bounds in ../core-geometry; the earlier capacity distribution still has the smaller support bound48. No exact optimum is claimed.

Run `python3 /path/to/this/packet/independent_check.py` from any directory. It uses only Python's standard library, reconstructs the native graph, verifies both seed trees, and directly verifies every edge expectation over all96 image trees. Output matches independent.out. WITNESS.json contains the complete finite certificate; no external file is needed for verification. PROOF.md cites the prior packet for shared definitions and structural context.
