# Separate modular and half-query improvement audit

Verdict: PASS for fixed-extra contrast and reconstruction from every modulus M≥3. This is a same-framework improvement; the earlier ±2 contrast remains valid.

Fix the extra edge from the last base vertex z to omitted core vertex u. Sum only over the d=2k−3 alternating-path edge signs. The unused-extra factor of two disappears from the core-endpoint words: uR contributes Q_uv and reverse(R)u contributes −(1−Q_uv). Of the two base-endpoint words, Ru contributes +1 and u reverse(R) is impossible. Therefore the exact signed contrast is 2Q_uv, either 0 or 2. All other words still cancel by an unused sign.

The number of responses per arc is 2^d=2^(2k−3), each using k−1 base vertices; at most binom(k,2)2^(2k−3) response evaluations suffice in total. The query family is independent of the unknown Q. k=2 works with two one-port queries. k=0,1 have unique cores and require none.

Reduction modulo any integer M≥3 commutes with this finite signed sum. The two possible residues 0 and 2 are distinct for every such M, including composite M and even M such as 4 and 6. Reconstruction uses a two-value lookup, not division by 2 or field linear algebra. In contrast, reducing the earlier symmetric ±2 experiment modulo 4 loses its sign, so that earlier experiment alone would not justify the even-modulus claim.

Modulo 2, the contrast vanishes. The stronger impossibility of using any full attachment responses modulo 2 follows from the separate classical Rédei parity theorem: every tournament has an odd number of directed Hamiltonian vertex words. Thus every such response is 1 modulo 2; for k≥2 there are multiple labelled cores and no information to separate them. The uniqueness cases k=0,1 are exceptions to an unqualified “never reconstructs” phrase. Finite parity checks do not replace this general theorem.

## Independent controls

modular_independent.py imports neither producer code nor subset DP. It enumerates every Hamiltonian vertex permutation for every randomized assignment, checks the per-word survivor classification, and compares the signed total to 2Q_uv.

- 4 exhaustive k=2 contrasts.
- 192 exhaustive k=3 contrasts, including every fixed attachment setting.
- 20 deterministic k=4 samples varying core, ordered endpoints, remaining-vertex order, and fixed attachment arcs.
- Every contrast checked modulo each M=3,...,24.
- All 2,184 individual response counts in these experiments are odd.

Normal and optimized outputs agree. The proof works for all M≥3; the finite modulus list is an implementation control only. No optimal port count, polynomial running time, practical oracle, or universal two-port reconstruction follows.

## Final proof, source, and second-bit extension

Reviewed frozen MODULAR_FAITHFULNESS.md, hash 2ace16955a840be4634770cc57f61ebdeba7432012c16159546128575540eab5, including both fixed-edge directions and the second-bit XOR theorem. Both additional conclusions are valid. Reversing the fixed edge gives 2Q_uv−2; summing the two experiments recovers the old ±2 identity.

For the second-bit result, write H=1+2r as an integer equality. The sum of the character over its nonempty Boolean cube is zero. Thus the exact contrast gives sum χr=Q_uv, and reducing this integer equality modulo 2 gives XOR_x ((H(T(x))−1)/2 mod 2)=Q_uv. This is legitimate integer cancellation before reduction, not division by a nonunit in Z/4Z. With the opposite fixed endpoint, XOR is 1−Q_uv. The independent permutation script additionally checks these second-bit XORs for all its 216 forward-fixed contrasts and accumulates residues from each individual response, not only from the final sum.

The parity reference was independently opened and verified: [Grinberg–Stanley, arXiv:2307.05569v1, Theorem 6.1](https://arxiv.org/html/2307.05569v1) states the oddness of the number of Hamiltonian paths of every tournament, including the empty convention. Its Theorem 7.1 also supplies the stated modulo-four refinement; no priority claim for that classical background is implied.

Producer normal and optimized scripts were independently rerun and matched: 840 fixed-endpoint tests and 27,528 response queries, with both endpoint directions and individual-response residues modulo 2,3,4,5,6,8,11. These DP checks and the independent permutation checks use different counting methods.
