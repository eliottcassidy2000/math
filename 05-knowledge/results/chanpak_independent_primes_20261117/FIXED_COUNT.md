# Fixed-count SCC rigidity under arc edits

Let H be the Hamiltonian path count. In the labeled tournament metagraph with single-arc reversals, restrict to one value of H. The vertex sets, relative order, and individual H counts of all nonsingleton strong components are constant along every path of this restricted graph. In fact each connected component is a Cartesian product of the fixed-H strong internal flip components and a permutohedron for each maximal run of singleton SCCs.

## Single flip proof, with the actual absolute change

Suppose an old arc a→b is between SCCs, and choose the ordinal cut A→B separating its endpoints. Write E_A(a) and S_B(b) for full paths ending at a and starting at b. Write W_A(a) for ordered covers of A whose second path starts a (with first path optional) and Z_B(b) for analogous covers with the first path ending b (with second optional). One arc reversal b→a changes the count by

    ΔH = W_A(a) Z_B(b) − E_A(a) S_B(b).

The lost paths are exactly the concatenations using a→b. A newly gained path uses b→a and its block segments have shape A₁ B₁ b→a A₂ B₂, with empty outer segments allowed (B₁ ends b, A₂ starts a; the notation for the arc is a boundary, not an extra occurrence). These are counted by the two marked covers independently. No other transition can use the removed arc, since a and b already lie in the inner segments.

Cutting each old A path just before a gives W_A(a)≥H(A); cutting each B path after b gives Z_B(b)≥H(B). The lost factor is at most H(A)H(B). Equality requires every A path to end at a and every B path to start b. These requirements force a to be the sink of A and b the source of B: if either has a contrary arc, extend that directed two-vertex path to a full path by insertions preserving relative order. Conversely, when they are a sink and a source the only change is an adjacent swap of two singleton SCCs, preserving H. This proves the cross-SCC claim. If an internal flip split an old SCC while preserving H, reversing that flip from the split tournament would be a neutral cross-SCC flip, which cannot merge anything. So remaining neutral flips stay inside a still-strong component. By the component product formula, a neutral internal flip holds its individual H fixed.

The factors in the Cartesian product are independent: inside a given nonsingleton component one can only use strong fixed-count flips; in a singleton run, adjacent swaps independently realize all permutations of that run. Singleton vertices cannot cross a nonsingleton SCC at fixed H.

## Even simultaneous exterior flips cannot merge old components at fixed H

Use the positive odd-cycle pack formula for H. If the edits are only between old SCCs, every old odd directed cycle stays intact because it was internal to an old SCC. Similarly every old spanning pack survives. If the new tournament has a directed triangle meeting multiple old SCCs, its single-cycle pack (all other vertices singleton) proves strict count increase.

Conversely, absence of such a triangle forces every vertex outside an old SCC to see it uniformly: mixed adjacency would partition that strongly connected old SCC into two parts with all crossing arcs in a single direction, otherwise a triangle through the outside vertex exists. Applying the same argument with the blocks reversed shows that the orientations between any two old SCCs are uniform by blocks. The quotient has no directed triangle and hence is transitive. Thus old SCCs do not merge. This proves that simultaneous edits of cross-component arcs cannot change the counts of those components without an actual rise in H if a merger occurs.

In particular neither a neutral single-edge walk nor any purely cross-component batch edit can transfer a factor 7 from one existing prime-count component to the other at total 7pq. An attempted fixed-target argument needs internal changes in a genuinely coordinated multi-arc edit, or a path allowed to leave the target value.
