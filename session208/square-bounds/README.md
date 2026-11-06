# Square-aware native metagraph bounds

For the actual labelled order-four path-complement graph, the certified intervals are:

- Optimal expected stretch:109/40 to211/66.
- Optimal expected fundamental-cut congestion:11/4 to269/84.

Neither exact optimum is determined. New lower bounds preserve each square's boundary subset; new upper bounds are explicit rational mixtures of symmetry images of saved spanning trees. The proofs are independent of the new preprint's global logarithmic theorem.

Run `python3 /path/to/packet/check.py` from any directory. Python's standard library suffices. It checks all finite shortest-path tables, four seed trees,192 image trees and every edge's rational expectation. The saved certificate does not require rerunning the heuristic search or numerical MILP diagnostic. See PROOF.md for scope and open separation-oracle obligations.

Independent replay: `python3 /path/to/packet/independent_check.py` uses tournament-first graph enumeration and Floyd-Warshall shortest paths, then checks all192 tree images. It needs only the standard library and sibling WITNESS.json; output matches independent.out.
