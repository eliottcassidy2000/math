# Slot-shadow collisions and sharp second-order reconstruction

This strengthens Session319's order-22 endpoint/shadow collision. Let an ordered q-slot cover of a tournament be a partition into q labelled directed paths, allowing empty slots. Its count is Q_T(q). Put a variable y_e on each unordered pair e, multiplying y_e whenever a slot path uses the pair consecutively. The resulting polynomial Z_T(q;y) is squarefree in all y_e. Its mixed partials at y=1 count covers simultaneously using prescribed unordered pairs. Write M_T(q;e) and M_T(q;e,f) for first and second moments. This is a path-use observable, not a cycle-packing weight. Ordinary Hamilton paths are the q=1 case.

## Two exact invariances

Let T be an ordinal sum of vertex-transitive blocks U_1→…→U_k, with singletons permitted. Replace any independent subset of blocks by their converses.

1. **Every labelled induced-subset scalar polynomial Q is preserved.** Within any slot the path progresses through the block order, entering each block at most once. For any labelled induced subset S, independently reversing each occupied slot segment in a switched block is a bijection of covers: joins to other blocks remain valid. This also preserves each slot segment's vertex set, size and all its undirected internal edges, giving equality of the entire multivariate polynomial when only internal edge variables are marked (even independently by slot). Vertex transitivity is not needed for this part.
2. **The entire first-moment edge vector on the full tournament is preserved at every q.** Internal marks are covered by the bijection above. For a seam joining u∈U_i to v∈U_j with i<j, choose the slot containing it. Its segments in blocks i and j must be nonempty, and those in every strictly intermediate block must be empty. The last point of block i's segment must be u and the first point of block j's segment must be v. In each block, choices are independent. Transitivity makes these endpoint counts flat on the full block: among covers in which this slot is nonempty, its specified start or end count is [Q_U(q)−Q_U(q−1)]/|U|. These two counts are equal and unchanged by converse. Hence for q≥1 the exact seam moment is
```
q · [Q_i(q)-Q_i(q-1)]/|U_i| · [Q_j(q)-Q_j(q-1)]/|U_j|
  · ∏_(i<l<j) Q_l(q-1) · ∏_(l<i or l>j) Q_l(q).
```
At q=1 it reduces to the adjacent seam law, since nonempty blocks have Q(0)=0. The statement is for the full blocks: on arbitrary induced subsets, endpoint flatness can fail and the first-moment vector need not be invariant. Indeed a three-vertex wedge already detects that failure.

For the Session319 Cayley block U of size 21, adjoining a terminal singleton gives non-converse-isomorphic tournaments U→z and U^op→z. Thus these two match every labelled induced-subset scalar polynomial, all full-tournament first moments, all internal-edge response orders, and the separate Hamilton endpoint vectors (even the joint first/last count matrix, since the last is forced). Using k blocks and a terminal singleton gives 2^k distinct converse-merged classes preserving all these invariants.

## The second-order threshold

For any tournament on n≥3 and any distinct u,v,w, put e={u,v},f={v,w}. A cover using both edges must use them in one slot and consecutively around v, because v occurs once. Thus its count is identically zero if the orientations at v are both incoming or both outgoing. Otherwise the three vertices form one forced directed run through v; putting this run and every other vertex in separate slots yields exactly (q)_(n−2) covers, where the expression is a falling factorial. All other covers use fewer than n−2 nonempty slots. Consequently
```
[q^(n−2)] M_T(q;e,f) = 1 if v is mixed on (e,f), else 0.
```
For a finite one-shot version, q=n−2 and all slots **required nonempty** gives exactly (n−2)! in the mixed case and zero otherwise. For q=n−2 with empties allowed, the support test (zero versus positive) works but the exact positive count need not be (n−2)!.

The wedge statuses determine, at each v, the orientations of all incident edges relative to one chosen incidence. If two tournaments on the same labels share these statuses, for each v all arcs incident to v are either unchanged or all reversed. The status must agree at both ends of every edge, so the complete graph forces a global common choice. Therefore the quadratic observable determines the tournament up to global converse. This is a sharp first/second-order separation for these observables, since the 22-vertex pair matches the entire first-order data.

In fact two slots with empties allowed suffice for the universal positive-support test: put the forced triple into one and any Rédei path on all remaining vertices into the other. This needs no asymptotics or extraction of polynomial coefficients.

The 22-vertex pair is already separated at q=1 by a second mixed moment: take u,v as the last two vertices of any Hamilton path in U, followed by z. The two edges {u,v} and {v,z} occur together in that path, but could not occur together in U^op→z, which has v→u and v→z. For the general one-slot result, see [Wedge reconstruction](WEDGE_RECONSTRUCTION.md).

If a marked start in a specified slot is allowed, one edge mark can instead recover its direction: require that slot start at u and use {u,v}. It must start u→v, since u is used only once. At q≥n−1 all remaining vertices can occupy singleton slots, giving a positive test iff u→v. This contrasts with unpolarized first-order edge moments and separate endpoint marginals, which do not see this correlation.

No new Chan–Pak coverage is asserted. The ordered-slot definition and prior Cayley witness are inherited; the stronger collision and second-order threshold are the present claims.

## Why 21 is minimal for a vertex-transitive non-self-converse block

This is a minimality theorem for the mechanism, not for slot-shadow collisions. A vertex-transitive tournament has odd order because its automorphism group has odd order: any involution swaps a pair, which would reverse the arc on that pair. For prime order p, Cauchy's theorem gives a p-cycle and hence a regular cyclic subgroup. Any tournament with a regular abelian automorphism subgroup is a Cayley tournament for that subgroup; inversion reverses all arcs. Apart from the singleton and prime orders, the only odd orders below 21 are 9 and 15.

At order 9 choose a transitive Sylow 3-subgroup P: take a Sylow subgroup containing a Sylow 3-subgroup of one point stabilizer and use the index 9. Inside the standard Sylow subgroup C3 wr C3 of S9 (after relabelling), the block rotation quotient of P is C3. If P contains a 9-cycle, use it. Otherwise a block rotator g has order 3, and the intersection of P with the base C3^3 is a nonzero invariant linear subspace of F3^3, since |P|≥9. Conjugation by g cyclically permutes coordinates. On any nonzero vector in the invariant subspace, repeatedly apply g−1 until the last nonzero iterate; it lies in the fixed diagonal, since (g−1)^3=0. The resulting nonzero diagonal z commutes with g. Their generated C3×C3 acts regularly, because g moves the blocks and z is regular inside them.

At order 15 the odd automorphism group is solvable by Feit–Thompson. A solvable primitive permutation group has prime-power degree: a minimal normal subgroup is elementary abelian, is transitive by primitivity, and is regular because it is abelian and the action is faithful. Thus there is a nontrivial block system of r blocks of size p with {r,p}={3,5}. The quotient is an odd transitive subgroup of S_r and is cyclic regular: for r=5 its order divides 15, and order 15 would be cyclic, but S5 contains no element of order 15; r=3 is immediate. The kernel K of the block action is consequently the stabilizer of each block and is transitive on it. Its restriction to each block is again the odd cyclic regular group of order p, by the same argument. Faithfulness embeds K as a subspace of F_p^r with each coordinate projection surjective.

If the projection of K to every pair of block coordinates has dimension 2, K is transitive on all cross pairs between those two blocks. Thus the orientation between any two blocks is uniform, so the tournament is a cyclic prime-order quotient inflated by identical prime-order vertex-transitive tournaments. Apply the inversion anti-isomorphism to the quotient and to each internal block (using fixed isomorphisms between them). The full tournament is self-converse.

Otherwise some coordinate pair of K has one-dimensional projection. Conjugating by a lift of the cyclic quotient gives dependence for all cyclic translates of that pair. These pairs form a connected r-cycle because r is prime. A value of one coordinate thus determines every coordinate, and K has dimension 1. Choose a quotient rotator g of order r by Cauchy's theorem (any order-r element is outside K). Its conjugation on K is trivial because r does not divide p−1. The resulting commuting groups of orders p and r act regularly on all 15 vertices, so their abelian Cayley tournament is self-converse. This completes the minimality proof for the vertex-transitive block mechanism.
