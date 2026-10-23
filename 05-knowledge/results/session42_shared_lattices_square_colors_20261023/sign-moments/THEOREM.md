# Sharp universal bound for permutation edge-sign moments

**PROVED:** Every nonempty finite simple graph F has |c(F)|<=1/3. Equality holds exactly when deleting isolated vertices leaves P3 or C4. The conclusion is independent of vertex labeling. These elementary proofs are recorded without a literature-priority claim.

## Definition and covariance

Let F have vertex set [n]. Choose a uniformly random total order sigma, and put

c(F)=E_sigma product_({i,j} in E(F),i<j) sign(rank_sigma(j)-rank_sigma(i)).

The empty graph has moment1. Isolated vertices do not affect the moment, since the induced order on the remaining vertices is uniform. Relabeling may change its sign, but changes each reoriented reference edge by a deterministic minus sign, so |c(F)| is an ordinary graph invariant. Reversing the total order multiplies the integrand by (-1)^|E(F)|; hence every odd-edge graph has c(F)=0.

## First-vertex recurrence

For v in [n] define epsilon_v=(-1)^(number of neighbors of v smaller than v). Condition on v being first in sigma. Every edge from v to a smaller label then contributes -1, every edge to a larger label contributes+1, and the remaining order is uniform on [n] minus v. Therefore, with inherited labels,

c(F)=(1/n) sum_v epsilon_v c(F-v).                    (1)

Equivalently, writing S(F)=n! c(F),

S(F)=sum_v epsilon_v S(F-v).                         (2)

This recurrence uses an actual first position, not a freely chosen local sign. It works for every finite simple graph, connected or not.

## Star boundary

Suppose some vertex meets every edge. Then F is a star with d>=1 leaves, together with isolated vertices. Relabeling changes at most the global sign. Put the center below all leaf labels, and realize the random order by independent uniform real coordinates U_v in [0,1]. Conditional on center coordinate t, the expected product of leaf comparisons is (1-2t)^d. Thus

|c(K_(1,d))| = |integral_0^1 (1-2t)^d dt|
             = 0 for odd d, and 1/(d+1) for even d.   (3)

For d>=1 this is at most1/3; equality occurs exactly at d=2, the path P3. This handles every case in which a deletion could leave no edges.

## Universal bound

Induct on the number of vertices. Nonempty graphs on at most two vertices have moment0. For an arbitrary nonempty F, handle the star boundary above if a vertex meets all edges. Otherwise F-v is nonempty for every v, so by the inductive bound and (1),

|c(F)| <= (1/n) sum_v |c(F-v)| <= 1/3.

No degree restriction, enumeration or asymptotic passage is required. In particular the argument does not leave a maximum-degree-four residual class.

## Exact equality classification

We strengthen the induction to classify equality. Isolated vertices can be removed at the outset. Moments multiply over disjoint connected components: their induced relative orders are independent and uniform. Consequently two or more nontrivial components have absolute moment at most1/9, by the universal bound, so equality requires a single connected nontrivial component.

The star case is exactly P3. Assume F is connected and not a star and |c(F)|=1/3. Equality in the triangle-inequality chain forces |c(F-v)|=1/3 for every v. The inductive classification says each F-v consists of a P3 or C4 and isolated vertices. In particular its maximum degree is at most2.

Also |E(F)| is even, since the moment is nonzero. Every F-v has nonzero moment and therefore an even number of edges. Thus every degree d_F(v) is even. Because F is connected and nontrivial, its degrees are positive.

Fix any v. For u adjacent to v, d_(F-v)(u)=d_F(u)-1<=2; the positive even integer d_F(u) is therefore2. For u not adjacent to v and u!=v, d_(F-v)(u)=d_F(u)<=2, so it too is2. Vary v to conclude every vertex of F has degree2. Hence F is a cycle. Deleting one vertex leaves a path, which can be an extremal P3/C4-plus-isolates only when it is P3. Thus F=C4.

Conversely, P3 has magnitude1/3 by (3). For C4 labeled around a cycle with edges12,23,34,14, every deletion is a P3. Direct substitution in (1) gives c(C4)=-1/3 (the sign depends on edge-reference labeling); hence |c(C4)|=1/3. An independent integral check writes the C4 moment, up to global sign, as E[(1-2|U_1-U_3|)^2]=1/3 after integrating the two opposite vertices. Adding isolates does not change the value.

This completes the all-size classification.

## Fourier and coupling consequence: precise handoff

The moment is the Walsh Fourier coefficient of the uniform measure on the inversion masks of permutations. If an independently derived quotient-coupling eigenvalue is |c(F)|², its nontrivial values are at most1/9, with P3 and C4 attaining the bound whenever their sectors occur. Thus the corresponding nontrivial spectral gap is at least8/9, and exactly8/9 if an extremal sector is present. This note proves the moment bound and equality classification; the identification of a particular coupling operator's sectors is a separate theorem and must be checked there.

## Exact controls

The checker constructs all inversion masks for n=2..7 and computes every coefficient simultaneously by an integer Walsh–Hadamard transform. It verifies the universal bound, vanishing of all odd-edge coefficients, equality masks, and (2) exhaustively through n=6, plus a deterministic spread at n=7. Equality-mask counts are 3 binom(n,3)+3 binom(n,4), representing P3 and C4 with arbitrary isolated labels. These finite controls support, but do not replace, the induction.

A targeted literature search for graphical inversion/sign-imbalance and transitive-tournament Fourier coefficients did not supply a cited version of this exact bound during this session. No claim of novelty is inferred from that search.
