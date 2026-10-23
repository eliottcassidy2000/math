# Integral augmentation lattice: independent audit

Verdict: PASS. Both publication housekeeping fixes were applied and inspected: the attribution now says “proposed identity,” and the order-seven census script locates its companion relative to __file__ instead of an absolute scratch path. No remaining correction is required.

## All-order identities and signs

For odd n=2m+1>=3 and B=(e_i-e_last), the explicit K_ij=M_ij-M_i,last+M_j,last is correct. Its even-size reduction mod2 is I+J, hence invertible. This agrees exactly with the audited strengthened root theorem: regularity is unnecessary for the specific sum-zero augmentation space. The integral lattice statement adds arithmetic information that rational equivalence forgets.

For the stated ordered B, its maximal minor after deleting row i equals (-1)^(n-1-i)=(-1)^i because n is odd. Pfaffian Cauchy–Binet therefore gives Pf(B^T M B)=sum_i(-1)^i Pf(M deleting i), with exactly the stated signs. The skew cofactor vector w_i=(-1)^i Pf(M deleting i) satisfies adj(M)=w w^T, not its negative. The singular determinant lemma gives det(M+J)=1^T adj(M)1. Since M+J=I+2A, every determinant identity follows with NO missing n factor. Finally adj(-M)=(-1)^(n-1)adj(M)=adj(M), so q_M(0)=detK for q_M(z)=1^T adj(zI-M)1.

An orthonormal real augmentation basis has a different volume normalization: det(B^T B)=n. The note correctly fixes the integral basis rather than importing an orthonormal determinant, which would introduce a factor. Under integral alternating congruence the positive elementary divisors repeat twice in ordinary SNF and their product is |PfK|. They describe the integral commutator pairing; no general equality with Hamiltonian data follows.

For n>=5, at least two of the n-1 old vertices lie on the same side of the last vertex. Their K entry is exactly an M entry, hence +/-1. Thus the gcd of all K entries is1 and d1=1, an all-order proof. It does not extend to cyclic n3, whose d1=3.

## Signed Redei distinction and minimal obstruction

The product |PfK| is determined by det(I+2A), hence adjacency spectrum. It is not the signed sum of Hamiltonian permutation signs. Independently enumerating literal permutations gives:
- cyclic5: PfK=-5, SNF alternating(1,5), H15, signedR=-5;
- Paley7: PfK=-49, SNF alternating(1,7,7), H189, signedR=-147.
These confirm the important49 versus147 hostile without relying on a possibly ambiguous signed-count convention.

For outmasks [18,4,9,19,6] and [30,4,8,18,6], independent SymPy SNF gives ordinary diagonal(1,1,1,1) in both cases. Literal permutations give respectively(H,R)=(11,-3),(5,-1), while direct small-Pfaffian formulas give Pf=-1 in both. Hence the complete alternating type fails to determine H. At n3 the transitive and cyclic cases have distinct profiles(1),(3), and n1 has only one tournament. The n5 collision is therefore minimal among odd orders. No even-order analogy is implicitly included.

## Census scope and independent reproduction

The independent script uses explicit Pf4/Pf6 formulas, not the producer's memoized Pfaffian implementation; normal-form computation, not a principal-minor surrogate, is used for all64 fixed-path n5 presentations and the named examples. It checks the full signed cofactor identity, all maximal B-minor signs, adj(M)=w w^T, and determinant equalities in those cases. Hamiltonian and signed counts are obtained by literal permutation enumeration, not endpoint DP.

For n7, an independently implemented complete32768-presentation census reproduces all profile multiplicities. Profiles are(1,1,d) for d in{1,3,5,7,9,11,13,15,17,19,21,23,25,27,31,33,35}, and(1,7,7), the latter occurring9 times. No two distinct profiles share the same product. Every tournament has a Hamiltonian path, so relabelling along such a path places it in this fixed-path universe; presentation counts are not isomorphism-class counts. This justifies coverage and the negative strict-refinement finding through order7, without extrapolating it to larger orders.

Read THM-174's existing square identity and the cited signed-Redei sources THM-1936/1966. The note appropriately presents an augmentation-lattice interpretation and concrete insufficiency witnesses, not a new square theorem or new finite-group classification. The statement that the integral classification is finer “in principle” is explicitly not claimed as a tournament-family strict-refinement witness.

Independent files: check_augmentation.py and check_augmentation.out. All tests passed. The producer's augmentation.py output was rerun and matches augmentation.out byte-for-byte. No Git operations or repository edits by this auditor.

## Final SHA256 pins

- AUGMENTATION-LATTICE.md: 62f9950ed573ebca98720fd8955ee47c1cce17800f934c9982e18bcda10af182
- augmentation.py: c0828455486925d0a2d8d4236195297ea2facc7f0a291374fdaae2bf41480387
- augmentation.out: d38b73d2f7829c6b4d8bc985aebbbdddf7fb00ec4db05a1e475b77f60e19c83f
- augmentation7.py: 3fb15d1e7ccf1e04a8d837ca8e4556d807e7acefd24e8a92397cacbbe6b719d4
- augmentation7.out: 8a79718cbe11632847e3baf514f344d642ed9ddcf1556e519e1dd199358124a6

All five paths above are in session35-recursions. The corrected augmentation7.py was also rerun to completion and its output matches augmentation7.out byte-for-byte.
