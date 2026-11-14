# What can hide inside k chronological paths?

For a tournament T, let κ(T) be the least number of directed Hamilton paths whose arcs together meet all directed cycles. In a rank order, call arcs going from a later to an earlier vertex backward. Define λ(T) by minimizing, over rank orders, the largest backward indegree or backward outdegree at any vertex.

**Theorem.** T is an induced subtournament of a tournament with κ≤k if and only if λ(T)≤k. For k≥1 and a rank realizing λ(T)≤k with e backward arcs, one can take a supertournament on at most n+1+k(n−1)−e vertices. This theorem concerns induced structure, not preservation of Hamilton path count.

## Proof: turn colored backward chains into full Hamilton paths

Necessity follows by restricting the common rank order of the larger certificate. Every backward arc of T survives in one of the k annotated paths, supplying at most k incoming and k outgoing backward arcs per vertex.

For sufficiency take the backward digraph of the given order. Split each vertex into a tail copy and a head copy to get a bipartite graph of maximum degree ≤k. Color its edges with k matching colors (Kőnig's line-coloring theorem; equivalently, complete to a regular bipartite multigraph and repeatedly remove a perfect matching). On the original vertices, each color is a directed linear forest: indegree and outdegree are at most one, and directed cycles are impossible since every edge decreases in rank. Include isolated vertices as singleton paths.

For each color i, list all paths of its forest in any order. Put a private fresh bridge between each successive pair of paths, forming one list of all old vertices and its private bridges. Add one further shared vertex z. Put all new vertices earlier than all old vertices in the rank order, with z the last of the new ones. Prepare a spanning path for each color i: first list all bridges *of other colors* in the common rank order, then z, then i's own list of forest paths with its private bridges.

Orient the augmented tournament by the old arcs on old pairs, and in increasing rank on all other pairs except where a prepared path requires the opposite arc. This instruction is consistent. Adjacent old vertices on a prepared path are forest edges; bridge–bridge adjacencies are always rank-forward; z leaves for an old vertex in the forward direction. Every adjacency of a private bridge to an old vertex occurs only on the path that owns that bridge. Thus no conflicting requirements arise. Every rank-backward old arc belongs to a forest and every newly backward arc was required by a prepared path. Those k literal Hamilton paths certify the bound. The number of bridges for a forest of n vertices and e_i arcs is n−1−e_i; summing and adding z gives the size bound. For k=0, λ=0 means T itself is transitive.

## Intrinsic cost versus induced cost

The rank parameter has an elementary upper bound λ(T)≤floor((n−1)/2). Take an order with no improving move that relocates any single vertex all the way to the start or the end, where improvement means increasing the number of forward arcs. For a vertex at position j, moving it to the start shows that at most floor((j−1)/2) of its outgoing arcs to the preceding vertices can be backward. Moving it to the end shows at most floor((n−j)/2) of its incoming arcs from later vertices can be backward. Such an order is constructed by repeatedly improving while possible. The bound is sharp for a regular tournament of odd order, since its first vertex in any rank has all its incoming arcs backward.

The two measures can differ by an unbounded additive amount. For the balanced apex B_r of the [rectangle theorem](BALANCED_TWO_RUN.md), λ(B_r)=ceil(r/2). In any rank, if the backward degree bound is k<r, then both sides of the forced rectangle X-after-apex to Y-before-apex have at least r−k vertices. Their degrees cannot exceed k, so k≥ceil(r/2). Conversely use the rectangle rank with k=ceil(r/2); the two apex fans each have size k and the rectangle has each side of size at most k. For r≥2 the intrinsic number is κ(B_r)=r−floor(sqrt(r)), higher by r/2−O(sqrt(r)). Thus κ≤k is not hereditary, even though its hereditary closure has the simple description above.

## Why the construction does not solve the spectrum

The added bridges provide missing transitions between disjoint backward chains, but they may create extra odd cycles and thus extra Hamilton paths. The induced-hull proof does not neutralize them. One positive design problem is now exact: formulate which private bridges can be made count-neutral without losing the prepared paths, or replace old arcs so that their additional packing weights become a controlled correction. The [previous induced-extension tariff](../chanpak_collatz_exact_addresses_20261114/TRIANGLE_TARIFF.md) is a necessary obstruction to any such neutralization that allows extra cycles while keeping a disjoint bank of cyclic components. It does not apply directly to a genuine changing-core synthesis.
