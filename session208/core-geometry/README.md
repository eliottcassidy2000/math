# Exact order-four native metagraph geometry

The actual path-complement graph has a bridgeless core made of two subdivided cubes joined through six disjoint squares. Exact symmetry-averaged tree distributions yield:

- Expected stretch between19/8 and157/48 (lower and upper bounds on the optimum).
- Expected fundamental-cut congestion between79/36 and13/4.

These results are independent of the new preprint's global logarithmic theorem. See PROOF.md for precise scope and the distinction between the two distributions. Run `python3 check.py`; only Python's standard library is needed. RESULT.json records all finite certificates. Neither optimum is claimed determined.

Independent exact replay: `python3 /path/to/packet/independent_check.py` works from any working directory and verifies all144 automorphism-image trees and every edge expectation. It uses only the standard library. Main and independent scripts print results and do not alter the packet.
