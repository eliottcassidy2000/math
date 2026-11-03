# Square-sum graph: audited core

The simple graph G_n has vertices 1 through n and edges joining distinct vertices with square sum. [The audit](AUDIT.md) describes the inherited components, highlighted pairs, centrality, cuts, and the first triangle; [the proof](DISTANCE_PROOF.md) establishes the sharp distance threshold for 2 and 4. Claims are proved or finite-exact within the stated scope, not statements about an unspecified graph, a triangle interpretation of a three-way split, or arbitrary centrality measures.

Run `python check.py` and `python -O check.py`, or verify the frozen packet and replay with `python verify_frozen_packet.py` and `python -O verify_frozen_packet.py`.
