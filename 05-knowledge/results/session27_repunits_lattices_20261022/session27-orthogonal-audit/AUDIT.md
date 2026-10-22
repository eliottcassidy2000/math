# Independent audit: split orthogonal avoidance certificate

Verdict: PASS, with no required mathematical fixes. This verdict is source-bound to the SHA-256 digests in `independent_results.json`, especially ORTHOGONAL_AVOIDANCE.md digest 06d429392539fb063a7f657c11540c6c5b94c650d8ea204957730143e4ca6de7. No source files or repositories were modified.

## Proof and scope checks

- Section 1 explicitly assumes an odd prime power, a split nondegenerate quadratic space in dimension 2r, and uniform sampling of all maximal totally isotropic subspaces, including both families. The quadratic-to-symmetric-bilinear equivalence is valid in odd characteristic. It does not apply the alternating/symplectic formula.
- The local replacement of u has Q(u-Q(u)w/beta(w,u))=0 and retains nonzero beta(w,u), so it produces a nondegenerate hyperbolic plane. Witt splitting gives the stated split quotient W-perp/W and the bijection of maximal spaces above W. The stated use of Witt decomposition is standard and sufficient.
- The singular-point count and incidence recurrence are algebraically correct. The recurrence begins with N_0=1, giving N_1=2 and N_r=2 product(i=1 to r-1)(q^i+1). For d=0 the product is empty; for d=r its j=0 term supplies exactly the factor 2. A non-totally-isotropic space cannot be contained in L.
- Section 2 correctly uses the FULL generated span poset, including nonisotropic elements. The subset-sign coefficients satisfy precisely its Mobius recurrence. They are not silently replaced by those of the ambient subspace lattice. Nonisotropic terms vanish by containment probability zero. Pruning them during coefficient updates is safe because a superspace of a nonisotropic subspace remains nonisotropic. For a retained isotropic W, every subspace below W is isotropic, so no discarded lower term can affect its recurrence.
- C(B) is exactly the cardinality of the avoiding subset of the finite sampling population, so C(B)>0 and A(B)>0 are equivalent to existence. The bound requires STRICT |B|<q^(r-1)+1. Equality alone is inconclusive, as stated.
- The graph counterexample correctly distinguishes pairwise orthogonality from linear rank: all three pairs are orthogonal in both configurations, while their triple spans have dimensions 2 and 3. The coefficient 2 for the dependent triple's two-space is correct.
- Exact extension-field arithmetic is explicitly required for nonprime q. Canonical input representatives, provenance, and completeness obligations are correctly separated from the counting procedure. Complexity is not misstated as polynomial.
- Section 6 is explicitly conditional on a split 48-dimensional form. The independently calculated threshold is 7^23+1=27,368,747,340,080,916,344. It does not claim an actual complete forbidden-set size, prove the split hypothesis for an unspecified construction, certify minimum 10, or prove a rank-96 extremality theorem. The assertion about the supplied manuscripts remains a reference to the separate review, not an independently audited manuscript conclusion here.

## Independent controls

`independent_control.py` uses actual vector-set closure over F3, rather than the source scripts' RREF representation. It enumerates singular projective points, then obtains each isotropic subspace by adjoining orthogonal vectors and deduplicating full vector sets. This gives dimension-level counts:

- r=1: 1, 2.
- r=2: 1, 16, 8.
- r=3: 1, 130, 520, 80.

For every isotropic subspace in these levels it directly checks the claimed containment probability, including d=0 and d=r. For r=2 it calculates rank-sensitive inclusion-exclusion weights independently from vector-set closure and bilinear tests, applies a subset zeta transform, and compares every one of 65,536 answers with direct avoidance among the eight enumerated maximal spaces. All agree. The full avoidance histogram equals the source result exactly: 0:41503, 1:17288, 2:5188, 3:1208, 4:284, 5:48, 6:16, 8:1.

For r=3 it directly counts 60 and 61 avoiders for the dependent and independent K3 configurations among its independently constructed 80 maximal isotropic spaces. Thus the probabilities are 60/80 and 61/80.

Source-script inspection also confirms the enumeration of 40 projective points and 130 planes in F3^4, the general RREF enumeration of 33,880 three-spaces in F3^6, and the 102 actual span-poset Mobius checks. Source scripts were copied into `reproduction/` for execution without altering the audited originals; reproduction outputs and run log reside there.

## Limit of this audit

This verifies the conditional finite-geometry theorem, its executable examples, and the numerical threshold. It does not establish any missing lattice-input completeness, actual gluing lattice, split discriminant form, or rank-96 minimum claim. No novelty assessment beyond the document's own explicit non-novelty disclaimer is claimed.

Both copied source scripts completed successfully, and both reproduced JSON outputs compare byte-for-byte identical with the audited source outputs (`cmp` returned success without differences).
