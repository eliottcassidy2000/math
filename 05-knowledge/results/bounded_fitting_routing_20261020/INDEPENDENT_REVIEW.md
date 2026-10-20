# Independent review of bounded Fitting routing

Reviewer: circulant-construction worker. Date: 2026-10-20.
Target: `session10-routing/BOUNDED_FITTING_ROUTING.md`, sections 1–5 as read at 02:10 UTC.
Verdict: **PASS** for the stated finite theorem and binary-word specialization.

1. **Quotient indegree bound.** For any G-orbit [v] mapping to [y], transport v by a group element so that its image is exactly y. This is legitimate because D commutes with G. Since ker D has two elements, there are at most two predecessor G-orbits. The two vector solutions can merge only into one predecessor orbit, never produce an additional predecessor orbit.
2. **Fitting coordinates.** N=ker D^e and W=im D^e are complementary at the stabilization exponent, and G preserves both because it commutes with D. Kernel dimension one implies a single nontrivial nilpotent Jordan block, although this stronger structural observation is unnecessary for routing.
3. **Periodicity after quotienting.** The potentially dangerous implication is correct: D^k v=g v gives D^(k ord(g))v=v by commutation, so its nilpotent component vanishes. Conversely D acts as a permutation on the finite set W, so every orbit represented there is periodic.
4. **Stable versus unstable predecessor.** At any y∈W the two preimages are w∈W and w+k∉W with 0≠k∈ker D⊆N. G-invariance of W prevents their identification in the quotient. Thus every periodic vertex is a branch and has exactly the stated stable and unstable predecessors.
5. **Cut paths and inverse.** Removing one edge into each branch leaves indegree≤1 and outdegree≤1, without creating sources. A source cannot enter a retained cycle; a terminal's backward path cannot enter a cycle. Both endpoint correspondences are inverse. A vertex cannot be selected as predecessor for two branches because it has one outgoing edge.
6. **Bound.** Sources lie outside im D, hence outside W. All paths reach W within e original edges, and their first edge entering W is cut. The source-to-branch and inverse paths therefore have at most e edges. Detecting the final source can require an additional solver call; the revised text correctly bounds predecessor edges rather than an exact number of solver calls.
7. **Cyclic words.** For m=e q with q odd, the stable test is zero XOR on every residue class modulo e. This follows either from the repeated (x+1) factor or directly from D^e=I+R^e and the fact that constants miss that subspace. Solving Dz=y and canonicalizing the two cyclic words justify O(e m²) time with elementary quadratic least rotation. No fast arbitrary-group canonicalization is being claimed.

This review does not review an Artin–Schreier extension not present in the target file at the stated read time. It also does not establish the prime-circulant sign criterion; that is proved and independently checked in the companion circulant packet. No mathematical defect was found.
