# Independent audit of the exponent-three graph bridge

2026-10-25. **PROVED elementary proof audit / FINITE-EXACT native controls.**
Reviewer: the session's divisor-ideal agent, independent of the producer
of [graph_bridge.md](graph_bridge.md). No mandatory correction found.

The complement-order obstruction is valid for any Boolean universe: a
monotone C-equivariant map sends x<=Cx to the empty set. Its explicit
forced-empty list on D_k is correct. Therefore injective monotonicity
fails for k>=4, and order reflection fails for every k>=2. The fixed-pair
obstruction rules out every even k without using order. These quantifiers,
together with the positive k=1 and k=3 maps, prove both claimed iff
classifications; the k=3 map deliberately has extra comparisons.

Each k=3 threshold is monotone and complemented by C; prime exchange
permutes its coordinates exactly as actual graph vertex reversal. The
casewise inverse from b,c and the three a-thresholds is exhaustive. The
four-vertex lower bound follows from14>8, independent of a census.

The tournament predicates also pass proof audit. The closing arc is read
by edge03. Without that arc, the leftmost cut needs b=1 and the rightmost
needs c=1; their two length-three backward intervals then cover all cuts.
Marked blue is fixed-path reversal, not arbitrary self-converse existence.
An empty backward-edge set is exactly the transitive tournament in its
unique Hamilton-path order; any backward edge creates a directed cycle
with the forward path segment. This proves the transitive iff p claim.
The final sign-versus-monotonicity obstruction uses the genuine pair
p|Cp and opposite nonzero f_k signs, so does not extrapolate the small
eight-state census.

[The independent checker](audit_graph_bridge.py) imports no producer code.
It reconstructs literal tournament arcs from the six displayed edge
predicates, checks strong connectivity by reachability from every vertex,
checks blue symmetry pairwise, and tests all directed triples for
transitivity. It also checks the closing arc, the g-coordinate difference,
the actual graph reflection, complement, injection, and all196 ordered
source pairs. Its14 sources supply140 tournament triples. Normal and
optimized outputs agree in [the saved receipt](audit_graph_bridge_receipt.json).

    python 05-knowledge/results/reflection_bridge_20261025/audit_graph_bridge.py
    python -O 05-knowledge/results/reflection_bridge_20261025/audit_graph_bridge.py

The small native check supports the specified k=3 map. The all-k
classification rests on the audited proofs above. This audit does not
claim a full independent reread of the inherited eight-state census.
