# Independent square-aware bounds audit

PASS, 2026-11-06. Independently reconstructed native adjacency by tournament-first enumeration. Checked each of the180 proper-square-subset cases using Floyd-Warshall all-pairs shortest paths, independently of the primary heap-based routine. The unit table[20,16,12,6] and auxiliary-length table[84,68,52,30] hold for every actual subset, including both two-edge shapes. Independently checked96 A-edge detours across the two length settings, giving5 and22.

The lower proofs are valid. For a fixed spanning tree, deleting a square's unselected edges from the whole core gives a supergraph of the tree, so its shortest distances are legitimate lower bounds. A full square cannot be selected by a tree; thus j<=3 and3s<=b. Summing the disjoint square contributions yields the stated expressions. The7/6 demand certificate proves1308/480=109/40 for expected stretch. Auxiliary lengths4/5 with original capacities1 and exact path/load adjoint give858/312=11/4 for congestion. No equality/attainment follows or is claimed. Leaves cannot create shortcuts and are exact bridge channels.

Independently validated all four saved39-edge trees, their exact path/load sum vectors, and all48 relabelling/converse automorphisms. Recounted every edge's path and load on all192 image trees. Exact rational mixtures produce211/66 expected stretch and269/84 expected load on every core edge. Support at most96 for each mixture is valid. Different distance/load seed pairs and mixture weights are correctly retained.

The two-coordinate finite-library support argument is sound, but does not provide a global separation oracle. The timed-out MILP and stronger proposed universal inequalities are explicitly not used as evidence. Exact optima remain open. No main logarithmic preprint theorem is assumed.

Portable independent_check.py reads sibling WITNESS.json, uses only the standard library, and prints independent.out. It uses an independent graph enumeration and Floyd-Warshall lower-table computation; no solver or floating arithmetic is needed.
