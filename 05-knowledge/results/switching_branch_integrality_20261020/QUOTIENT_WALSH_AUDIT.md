# Independent switching-quotient audit for the uniform 3/5 obstruction

Status: PROVED. This supplies exactly the product-of-odd-p-stabilizer-orders identity needed by the uniform obstruction, without assuming the stronger prime-valuation-distribution theorem.

## Quotient Walsh determinant identity

Let V = F_2^{E(K_n)}, let U be any S_n-invariant reversal subspace, and let Q be the affine V/U-torsor of tournaments modulo U. Put d = dim(V/U), R = Z[1/2]. The character functions phi_F on Q are indexed exactly by F in U^perp. They form an R-basis of R^Q because the Walsh inverse divides only by 2^d.

The S_n-action permutes these characters with the orientation-sign cocycle. If the stabilizer of F acts with nontrivial sign, its signed character orbit has no nonzero invariant vector: any such coefficient would equal its negative. If the sign is trivial, transporting coefficient 1 around that orbit yields a well-defined invariant signed sum Phi_F. Thus the Phi_F, indexed by sign-even graph orbits in U^perp, form an R-basis of (R^Q)^{S_n}.

Orbit indicator functions on Q form another R-basis of the same invariant module. Therefore the evaluation matrix

    M[q,F] = Phi_F(q)

is integral and invertible over R. In particular det(M) is a unit of R, hence a signed power of two. Let h_q = |Stab(q)| and k_F = |Aut(F)|. Grouping Walsh character orthogonality by vertex relabeling orbits gives

    M^t diag(n!/h_q) M = 2^d diag(n!/k_F).

Invert and rearrange (equivalently use the second Parseval relation) to get

    M diag(k_F) M^t = 2^d diag(h_q).

Taking determinants and valuations at any odd prime p gives

    sum_q v_p(h_q) = sum_F v_p(k_F).

No division by n! was used in the integral-basis claim. This is important: the denominators of the inverse matrix really are dyadic even when p divides n!.

## Consequence for any equivariant map bijective on types

Such a map matches each source orbit q to a distinct target orbit F and forces Stab(q) to be a subgroup of Aut(F) after compatible relabeling. Hence every difference v_p(k_F)-v_p(h_q) is nonnegative. The displayed aggregate equality forces every difference to vanish, separately for every odd prime p.

In particular the preimage orbit of the empty graph must have a stabilizer containing full ambient Sylow 3- and Sylow 5-subgroups.

## Switching specialization and the involution prohibition

For U equal to the cut space, U^perp is exactly the Eulerian edge space, and dim(V/U) = binomial(n-1,2). The target predicates are both Royle sign-even and Eulerian. The source stabilizer is the automorphism group of the switching class, not the automorphism group of a chosen tournament representative.

THM-479's level law says every permutation fixing a tournament switching class has all its cycle lengths of the same 2-adic valuation. Consequently every nontrivial involution in a switching-class stabilizer is fixed-point-free. In particular a double transposition supported on four points is forbidden whenever n >= 5.

This is precisely the hypothesis needed by the uniform 3/5 group obstruction. The support argument can be kept elementary:

- Two distinct single 3-cycles with supports overlapping in two points generate A_4 on their four-point union, hence a double transposition.
- If their supports overlap in one point, conjugate one by the other; that cycle and its original have a two-point overlap, again producing a double transposition.
- Thus all distinct supports of single 3-cycles in the stabilizer are disjoint.
- A full Sylow 3-subgroup of S_n contains floor(n/3) disjoint single 3-cycles. Their supports consequently constitute all such supports, a characteristic block family for the stabilizer.
- With k = floor(n/3), r = n-3k <= 2, the stabilizer embeds into S_3 wr S_k times S_r. Its 5-adic order is at most v_5(k!), which is strictly smaller than v_5(n!) for all n >= 5.

The last statement follows already from floor(n/5) > floor(floor(n/3)/5) for n >= 5. Thus the required full Sylow 5-subgroup cannot exist. This proves the switching-class half of HYP-9172(a) for every n >= 5, once the joint writeup includes the common group lemma.

## Inheritance

The affine character mechanism is THM-4531, `01-canon/theorems/THM-4531-twisted-mallows-sloane-principle-and-no-natural-bijection-for-switching-classes.md`, Theorem G. The full-tournament Walsh identity is THM-4721, `01-canon/theorems/THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md`, Theorem W. The already present broader quotient statement is `05-knowledge/results/prime_local_walsh_bijection_20261020/QUOTIENT_COROLLARY.md`. The argument above independently rederives just the needed specialization.
