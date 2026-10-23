# Static rooted motifs need a linear cutoff on the cycle-phase family

**PROVED, scoped observer lower bound. FINITE-EXACT rooted-canonical controls for g=3,4,5.** This separate addendum does not change the previously frozen five-file marked-cell packet.

## 1. Observer and conclusion

Let rho_k(T,v) be the complete vector of counts of all singly-rooted induced tournament types on at most k vertices, with v the distinguished root. Initialize EVERY vertex with its own vector rho_k(T,v), then perform arbitrary numbers of ordinary directed vertex-refinement rounds, retaining the old color and the multisets of outgoing and incoming neighbor colors.

For every g>=3, the two nonzero-phase tournaments T_(g,+1) and T_(g,-1) from GIRTH_PHASE_HOSTILES.md have identical rho_g vectors at each corresponding actual vertex. Their ordinary unindividualized refinement records remain identical at every depth. Nevertheless they are distinct strong self-converse tournaments of order 3g.

Consequently, a class separator for ALL strong self-converse tournaments of order n whose only observer is the complete static singly-rooted motif profile through cutoff k(n), followed by ordinary vertex refinement, must satisfy

    k(3g) >= g+1       for every g>=3.                       (1)

This is a linear necessary cutoff along these orders, not a claim of sufficiency, a computational lower bound, or a restriction on observers retaining actual marks, multi-vertex overlaps, path flags, or complete adjacency.

## 2. Root-fixed motif bijections

Fix an actual root u=(r,x) and a size t<=g. Compare subsets containing u in the two phase tournaments. The stronger statement with the cell identity of every motif vertex retained is also true.

If a subset uses fewer than all g cells, the support induced by those cells is a forest. Choose the phase-removing cell rotations along each tree, taking rotation zero in the component at the root cell. The gauge equations are solvable on a forest with that prescribed root value. Thus the resulting isomorphism between the unions of occupied cells fixes u exactly. Applied to subsets with any specified cell occupancy, it bijects their rooted induced motif types. It need not fix every other actual vertex, nor is such extra fixing claimed.

If a subset uses all g cells, then t<=g forces t=g and one vertex per cell. On such a transversal, apply

    (i,y) -> (i,-y+2x).                                    (2)

This fixes u=(r,x). It changes the closing phase from +1 to -1, preserves all path matchings and uniform blocks, and preserves every induced arc because the transversal contains no internal C3 arc. Equation (2) is involutive, so it is a bijection of the corresponding rooted transversals.

Summing the two cases proves equality of every rooted count at each fixed actual vertex. The bijections are allowed to depend on the occupied cell set. They do not assemble into a global tournament isomorphism; the nonzero cycle phase is exactly why no such legal global rotation exists.

## 3. Why ordinary refinement cannot recover the loss

The simultaneous cyclic shift (i,x)->(i,x+1) is an automorphism of each tournament. Therefore every rho_g profile is constant on each triple D_i. The pointwise equality above identifies these cellwise profiles between the two tournaments.

The partition into D_i is equitable and its directed quotient is phase-independent: internal outgoing degree is one, support blocks have row sum one or two according to direction, and all other blocks are uniform. Suppose a refinement round has a color constant on each D_i and the same corresponding color in the two tournaments. The outgoing and incoming counts into every color are sums of the same fixed cell-to-cell degrees. The next round consequently has the same two properties. Induction proves indistinguishability at every depth.

For the FULL rooted profile through any cutoff at least two, rooted two-vertex counts include the outdegree. Since the cell degrees are distinct, the initial color partition is already exactly the intrinsic triple partition. It is stable immediately. The argument does not depend on this extra observation: any coarser phase-blind cellwise initial profile also remains unable to distinguish the two graphs.

Combining this with their proved nonisomorphism gives (1). Allowing an unbounded number of these ordinary refinement rounds does not reduce the necessary motif cutoff.

## 4. Static rooting is not actual individualization

A root in rho_k(T,v) is a temporary distinguished position inside each counted induced subgraph. Computing this vector for every v gives a static vertex tag. It does NOT distinguish one actual vertex globally during neighbor refinement.

The latter is a different observer. GIRTH_PHASE_HOSTILES.md proves that an actual distinguished vertex plus the intrinsic cell colors makes these same examples discrete within 1+floor(g/2) rounds. The root-fixed local bijections in Section2 and the global repair by a mark are therefore compatible, not competing claims.

## 5. Reconciliation with the incoming prime-cyclic bound

The incoming bc3c88ec033 packet, `05-knowledge/results/metagraph_five_port_phase_frontier_20261023.md` Section3 and its routed `metagraph_five_census_hostile_20261023/README.md`, already proves fixed-cutoff failure and a necessary asymptotic scale sqrt(2 log2 p) for the observer on PRIME-CYCLIC tournaments.

The new construction is a different subclass: strong self-converse tournaments of composite order 3g with intrinsic triple cells. It gives the explicit linear boundary (1) for a separator on the general strong self-converse class. It neither asserts that the examples are prime-cyclic nor improves the prime-specific bound within that subclass. This is a repository-relative structural corollary, not an external literature-priority claim. The bounded inheritance search checked the incoming rooted-five/phase frontier and the relevant minimal-incidence and marked-cell routes for static rooted motifs, cutoff growth, and prime-cyclic scope.

## 6. Exact controls and replay

The independent standard-library checker `check_rooted_cutoff.py` builds both phase tournaments afresh and imports no code from the frozen producer. For g=3,4,5 it enumerates every subset containing every actual root through size g, canonicalizing with all permutations that fix the root position. Full profile equality is compared exactly, not by hashes. It checks every cell's profile uniformity, obtains the stable ordinary-refinement partition, and verifies (2) entrywise on all rooted full transversals. The receipt retains the full profiles at cell representatives and per-vertex profile digests for inspection.

    python marked-cells/check_rooted_cutoff.py
    python -O marked-cells/check_rooted_cutoff.py

Both modes reproduce `rooted_cutoff_receipt.json` with 3,612 explicit checks, including 81,324,1,215 root-fixed transversals at g=3,4,5 respectively. The all-g conclusion is proved by the root-fixed bijections and equitable induction, not inferred from these finite tests.
