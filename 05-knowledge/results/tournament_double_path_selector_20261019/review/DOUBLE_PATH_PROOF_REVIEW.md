# Independent review: double-path selector and forest primitive floor

2026-10-19. **PASS.** Reviewed `session6-padding-tree-followup/DOUBLE_PATH_SELECTOR.md`, SHA-256 `66840bc2a0cbff9540993e3e4a2ccb45400806cf07238e9c767286863c01cedc`.

The final producer change only updates the status header to acknowledge this completed review; the mathematical text is unchanged.

This proof uses the explicitly listed inherited parity results: all oriented trees through six vertices are parity rigid, every oriented P4 has odd parity, the relevant symmetric seven-trees are good, and the deletion–reversal and component-shuffle formulas. This review verifies the new deductions from these inputs, rather than re-proving the inherited classification.

## 1. Universal selector

The argument is correct with all isolated vertices retained. For the asymmetric seven-spider A plus an isolate, deleting an A edge gives rigid components of orders a, 7-a, and 1. Its shuffle coefficient is 8 binom(7,a), always even; therefore the common defect F is orientation independent.

For D with ordered edges 01,02,03,14,15,26,47, deleting either critical leaf edge 03 or 15 gives F. Deleting the central edge gives P4 plus P4, with shuffle coefficient 70. Deleting 02 or 14 gives orders two and six, with shuffle coefficient 28. Deleting 26 or 47 gives a symmetric seven-tree plus an isolate, with shuffle coefficient eight. Thus all noncritical orientation derivatives vanish pointwise on every host.

The mask-zero anchor is genuinely odd. Its two rooted oriented P4 halves are isomorphic. On any unordered partition into two four-sets, identify the embeddings of either half through this rooted isomorphism. Each pair of embeddings has two assignments to the labelled halves, and exactly one satisfies the oriented central edge. The contribution is the product of two odd P4 counts, hence odd. There are 35 unordered partitions. This proves

    phi_(D,m) = 1 + (bit_2(m) + bit_4(m)) F.

The half-swap is the only nonidentity undirected automorphism. It reverses the central edge, so it fixes no orientation. It preserves the critical-bit sum. Consequently there are exactly 32 universally odd oriented isomorphism classes and 32 nonconstant classes, all with trivial oriented automorphisms.

## 2. Nonconstancy and independent counts

A new verifier, `verify_double_path_independent.py`, imports no producer code. It counts all D orientations by literal permutations, separately checks selected counts through degree-ordered injective backtracking, and enumerates all 8! maps for undirected automorphisms. Results:

- Transitive eight-host: F = 480 and emb(D, mask 5) = 165.
- Host [210,188,169,81,164,201,22,72]: F = 617 and emb(D, mask 5) = 302.
- Every orientation and every first derivative on both hosts agrees with the selector.
- The automorphism group has order two, all oriented orbits have size two, and the split is 32+32.

Thus F is genuinely nonconstant. The finite host checks are witnesses and implementation controls; the universal assertion follows from the proof, not from sampling.

## 3. All-order forest primitive floor

The new all-order statement is valid. Iterating source evaluation on an oriented n-forest leaves a sum, with actual integer multiplicities, of spanning r-forest embedding functions. For r <= 6 every component is parity rigid, so each forest function is constant by the multinomial shuffle formula. The Taylor projector P annihilates positive-degree constants because its value on a degree-r constant is (1-1)^r. Therefore

    P h^(n-r) phi_S = 0,  1 <= r <= min(6,n).

These are exactly the primitive coordinates in the integral reconstruction basis reduced modulo two. No bound on n is used. Evaluation on any r-host with n-r successive universal sources therefore equals the transitive n-host value for r <= 6. This conclusion is limited to source-stack replacements; the text correctly avoids arbitrary-subtournament or minor invariance.

Sharpness at core order seven is independently verified for edges 01,02,03,14,25,46 with mask 52: the counts are 66 on the transitive host and 81 on [118,44,96,85,102,8,34]. The latter host has no universal source. Both direct permutations and injective backtracking agree. Thus the seven threshold cannot be increased. The directed-triangle counterexample correctly shows why the forest hypothesis is needed.

## 4. The degree-seven primitive of D

For the seven-spider A, h phi_A is constant because every eligible source deletion is a six-vertex forest. Evaluating at a transitive six-host identifies that constant with c = phi_A(transitive seven-host). Hence kappa_A = phi_A + c is primitive. Since D_1 of the degree-seven constant is eight times the degree-eight constant, it vanishes modulo two, giving F = D_1 kappa_A.

In the integral source-free reconstruction basis, D_1 is injective on this primitive summand. Thus D has exactly the displayed possible degree-seven defect and no degree-eight primitive component; all lower nonempty primitive components vanish. The next ordinary padding kills this defect because d squared is zero. This deduction correctly explains the inherited nine-host phenomenon without claiming a general tree classification.

The wording now correctly says the orientation functions **take** two affine values. Their linear span can have four vectors, so the previous informal use of “span two values” would have been misleading.

## 5. Credited Hamiltonian-path consumer

The final coordinate translation is also correct. A universal source must be the first vertex of a directed Hamiltonian path, so repeated h leaves the same path-count function in lower degrees. Applying P at a source-free core gives the alternating subset sum of the path counts. Substituting the explicitly inherited Rédei–Berge odd-cycle-collection formula, inclusion-exclusion cancels every collection that leaves an unused vertex. The surviving terms are exactly spanning odd-cycle covers, with weight two to the number of cycles. For a nonempty core this weight is even. The text properly credits this as a coordinate interpretation of an existing theorem and makes no new 7/21 exclusion claim.

## 6. Reproduction

    python3 verify_double_path_independent.py > double_path_independent_receipt.json
    python3 -O verify_double_path_independent.py > double_path_independent_receipt_optimized.json

Both completed successfully and the JSON receipts are byte-identical. The script uses explicit failure checks, not Python assertions; it has no input paths or producer imports. The receipt retains both complete orientation histograms, automorphisms, class split, and sharp seven-core counts.
