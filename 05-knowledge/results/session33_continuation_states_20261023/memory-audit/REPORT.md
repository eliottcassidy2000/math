# Independent audit: fixed-path metagraph involution memory

Verdict: PASS. One implementation correction was requested and applied: producer originally used Python assertions, which disappear under -O; frozen producer now uses explicit checks. No mathematical correction required. Earlier composition-audit artifacts were not changed. No Git mutation or attachment material.

## Independent exact results

The audit constructs every fixed-path tiling for n=3,4,5,6. It canonicalizes using every vertex permutation and converse, without the producer's degree-sort pruning and without importing producer functions. Node names differ but partitions agree, with an explicit label transport for n<=5. It constructs A and fibre sizes independently and compares them after that transport.

- n3: 2 tilings,2 nodes,2 ordered-pair states,2 colored states; defect rank0 and linear dimension2.
- n4: 8 tilings,3 nodes,5 ordered-pair states,5 colored states; defect rank2 and linear dimension5.
- n5: 64 tilings,10 nodes,38 ordered-pair states,38 colored states; defect rank9 and linear dimension19.
- n6: 1024 tilings,34 nodes,356 ordered-pair states,364 colored states. All eight supplied ordered mixed-color witnesses pass. In particular521 and27 have equal current/next merged-node classes but different colors. Canonical labels need not equal the producer's labels; the underlying equivalence is checked.

For n3..5, independently reconstructed exact rational checks include defect Gram identity, rank[J,PJ], support components/bipartiteness, memory powers k=0..6, projected times0..9, the two-time moment-Gram rank, arbitrary signed nonuniform forcing powers0..5, and distinct same-node sources with different next output at n4,n5. Ten-time observation words give precisely the ordered-pair quotient. Color is checked invariant under complement on every tiling at every tested order.

The rational two-time moment Gram [[F,A],[A,F]] has rank r+rank(F-A F^-1 A), providing a second direct certificate of the minimal dilation dimension, independent of the producer's rank[J,PJ] calculation. No irrational square-root approximation or floating threshold is used.

## Mathematical audit

1. Actual operation and deterministic state.

The operation complements free-edge bits while retaining the Hamiltonian path arcs, not full converse. Output merges tournament isomorphism and converse classes. Since tau²=id, all node-observation words alternate current and next nodes. Any adequate deterministic state must recover both; the ordered pair recovers all later outputs. Thus its refinement minimality is exact, not merely a cardinality estimate. Adding invariant color gives the colored pair and needs an extra distinction exactly at mixed supports. Nothing here licenses deletion, extension, or another operation without new refinement.

2. Signed linear state.

The common nullspace of U* and U*P is (V+PV)^perp, so no linear realization for all signed initial sources can have dimension below dim(V+PV). This subspace is P-invariant, establishing sufficiency. Projection of PV off V has Gram I-K²; hence dimension is r+rank(I-K²). The inherited support theorem makes this 2r-c-b. Loops make components nonbipartite. This is not a count of positive coordinate states.

3. Positive stochastic state lower bound.

The new extremality argument is valid and stronger than merely warning about positivity. With nonnegative column-stochastic O,T and probability encodings p_x for every pure tiling, O p_x=e_i forces each column used by p_x to output e_i. Likewise OT p_x=e_j forces every used OT column to output e_j. Distinct ordered pairs cannot share such a hidden coordinate, and each pair needs at least one coordinate. Thus h>=#ordered pairs, attained by the deterministic swap system. A nonlinear assignment x->p_x and omission of full-state intertwining do not weaken this proof. Joint color/node outputs give the analogous colored lower bound. The restriction to every pure tiling is essential: uniform-fibre initial states give mixed next outputs and do not establish this bound. Exact n5 separation38 positive versus19 signed is sound.

4. Minimal selfadjoint involution compression.

K is a selfadjoint contraction. Delta=sqrt(I-K²) commutes with K, and its range is K-invariant. On the observed space plus range(Delta), the displayed block [[K,Delta],[Delta,-K]] is a selfadjoint involution. Every involution with observed block K has BB*=I-K², forcing auxiliary dimension at least rank(Delta). The actual cyclic restriction P|(V+PV) attains this and is identified by polar decomposition. This is not a power dilation of K: even compressed powers are I rather than K². Nor does the construction give a positive coordinate transition, fixed-point-free tiling permutation, integral fibres, color, or actual tournament representatives.

5. Memory and initialization.

P²=I yields BB*=I-K² and BD=-KB. Iteration gives BD^kB*=(-1)^k K^k(I-K²). The order written in the request is equivalent because K commutes with its polynomial I-K². Eliminating hidden coordinates gives the stated finite convolution plus BD^t z0. Its forcing is (-K)^t eta, with eta=Bz0=y1-Ky0. Therefore K determines the kernel sequence but not arbitrary-source initial forcing. The producer correctly states this caveat. Inserting fibre projection at each step changes literal complement into rerandomized dynamics K^t (or similar row-stochastic T^t); equality for every time occurs exactly at zero defect.

## Inheritance and status

Read the actual boundary_support_20261020_metagraph.md Proposition0 and current THM-840/THM-796 frontmatter. The Gram and r-c-b formulas are already PROVED there and explicitly inherited. THM-840 has PROVED operation-kernel and minimal-linear-sidecar statements plus finite-exact replays. THM-796 has PROVED incidence/recursive/color statements and finite-exact coupling, non-lumpability and continuation census. The new claims do not silently upgrade these finite computations to all-order connectivity or classification.

The Halmos construction is classical scaffolding; this audit proves the displayed finite-dimensional identity directly rather than independently reviewing its literature citation. No historical-priority claim is certified.

## Deliverables and reproducibility

Run python3 memory-audit/audit.py and python3 -O memory-audit/audit.py from the workspace parent. All audit checks use explicit exceptions. Files are audit.py, receipt.json, receipt.normal.json, normal.out, optimized.out, SHA256SUMS. Normal and optimized receipts must be byte-identical. n6 matrix ranks are intentionally not claimed or computed; the n6 result is the complete node/pair/color census and mixed witnesses.
