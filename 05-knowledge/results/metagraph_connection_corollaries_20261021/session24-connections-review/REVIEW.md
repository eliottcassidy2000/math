# Independent metagraph connection audit

Status: APPROVED mathematically and computationally after the producer corrected the boundary-port singleton scope. Approval binds only the exact producer bytes in APPROVED_SHA256SUMS.

Read-only inspection included repository instructions and current truth-surface guidance, current THM-2477-rigid-strong-classification-through-sixty-two.md, THM-2454-pure-blue-alphabet-and-center-law.md, and THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md. The published converse_exception_classification_20261021 packet has the approved all-order classification; older local provisional prose is not the governing status. No repository mutation, Git operation, or silent producer edit was performed.

## Mathematical gates

1. Signed exceptionality forces ordinary vertex transitivity and therefore strong connectivity. Combining the current strong rigid alphabet {1,C3} and strong pure-blue alphabet {1,C3,R5} with the seed controls gives exactly E intersect R={C3}, E intersect P={C3,R5}. The faithful singleton action is not exceptional. These are corollaries of classifications, not new proofs of their dependencies. The P7 and transitive controls demonstrate incomparability on the correct distinct carriers.

2. The bottom triangle copies in A[C3] are canonical, even when A is composite. Here is a direct hypothesis audit that avoids confusing them with maximal strong modules: they are exactly the cyclic three-vertex modules. A cyclic triangle crossing copies cannot have two vertices in one copy, because any other copy has uniform orientation to those two. If it has one vertex in each of three copies, its cyclic quotient forces every leftover vertex of any one copy to see the vertices in the other two in opposite directions. It consequently fails to be a module. Thus all automorphisms and anti-automorphisms preserve the bottom triangle partition. All independent local rotations and all quotient automorphisms lift, proving Aut(A[C3])=C3^|A| semidirect Aut(A). This reasoning applies to every nonempty A, exactly as in the strengthened packet. Isomorphisms and anti-isomorphisms recover the quotient, so triangle substitution injects converse-merged node classes and preserves self-converse versus non-self-converse status. At a non-self-converse merged pair THM-4731 gives degree=h; the ordinary path-orbit injection therefore proves degree monotonicity there too. At A=K1 and A=T2 the self-converse formula gives zero downstairs, so these small boundaries cause no exception.

3. Block-contiguous path substitution gives a well-defined injection on ordinary orbits. All three triangle Hamilton paths are rotation-equivalent, and any equivalence upstairs descends to an equivalence of projected base paths. Given an anti-automorphism of A, lift it using local triangle anti-isomorphisms: reversing the lifted path then agrees on quotient block order with lifting the reversed base path. Anti-map choice disappears modulo ordinary automorphisms. This proves equivariance. Injectivity forbids a nonfixed base orbit from mapping to a fixed orbit, so two-cycles inject. The paired-degree inequality does not follow merely by subtracting two unrelated monotone quantities, and the producer correctly avoids that fallacy.

4. THM-4731 Statement 3 uses the fundamental paired-black MULTIGRAPH, with loops contributing two and self-converse degree (h-b)/2. All exceptions are self-converse. The independently audited path-cover receipts give C3[C3] degree14 and R5[C3] degree80646; the seed P7 has degree3. Every exceptional tower is either C3/R5 themselves, a positive-height tower over one of them, or a tower over P7. Equivariant injection therefore proves degree0 exactly for C3/R5 and degree at least3 otherwise. This is not distinct-neighbor degree and makes no bridge exclusion or full adjacency claim.

5. The port defect is real: at a boundary between consecutive blocks of sizes r,s, one of rs cross arcs is retained and all others reversed. Uniformity fails iff rs>1. The original phrase “more than one vertex in total” incorrectly included two singleton blocks; the producer corrected this to at least one nonsingleton, equivalently product greater than one. Each nonsingleton block with a neighboring block loses modularity because the external boundary endpoint distinguishes its path endpoint. For the stated nontrivial triangle substitution all blocks have such a neighbor. The quotient-first construction retains uniformity. Different Hamilton counts certify the claimed failure even up to isomorphism; module loss alone would not have certified that stronger assertion.

## Reproduction

Producer check.py was replayed normally and with python3 -O. Both pass all 15 checks and independently count paths by subset DP and all 9! vertex permutations in each of the four cases:

- C3[C3]: 3159.
- Fixed-path complement of the lift: 2983.
- Quotient complemented before substitution: 27.
- Both factors complemented before substitution: 1.

The DP permits updates from zero-valued states but does not create spurious paths; every nonzero endpoint is present in its mask. Diagonal zero handling and preservation of both directed entries of the path pair are correct. Explicit exception guards survive -O.

The audit-owned hostile.py imports no producer code and passes 78 checks in both modes. It exhausts all cyclic three-vertex modules of inflations of C3, R5, P7, composite C3[C3], off-family T3, K1, T2, and the non-self-converse quotient C3→1; they are precisely the bottom blocks. It tests every adjacent block size pair r,s=1,...,4, including the singleton-singleton exception and each one-sided nonsingleton case. Finally it independently enumerates all 9! permutations to find the full 81-element automorphism group and 3159 Hamilton paths of C3[C3], canonizes its 39 path orbits, and finds 11 converse-fixed orbits, 14 paired cycles, and the unique fixed block-contiguous orbit. This finite control supplements rather than replaces the all-order injection proof.

Reproduction commands:

    python3 session24-connections/check.py
    python3 -O session24-connections/check.py
    python3 session24-connections-review/hostile.py
    python3 -O session24-connections-review/hostile.py

Normal and optimized receipts are byte-identical for each script. Publication allowlist: the two producer files in APPROVED_SHA256SUMS, plus this review, hostile.py, four replay receipts, APPROVED_SHA256SUMS and REVIEW_SHA256SUMS. No unrelated scratch data or cache is approved.

The strengthened all-node injection and black-degree monotonicity are approved. The safe sentence “does NOT assert a homomorphism” is important: the explicit noncommutation witness refutes natural lifted-path commutation, and is not by itself a proof that no unmarked graph homomorphism can exist. No such stronger negative is claimed.
