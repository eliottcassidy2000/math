# Independent observer-profile audit

Verdict: PASS for the submitted computational certificates and stated Reynolds normalization, with one non-mathematical certificate cleanup and one exposition caveat.

## Independent computations

The audit imports none of the submitted source scripts. `audit.py` generates ordinary graph and inversion-twisted tournament orbits using adjacent vertex transpositions and breadth-first search on the Boolean cube. It reconstructs the n=6 profile matrix from distinct tournament orbit members, restores permutation multiplicity via 6!/|orbit|, and compares every entry with the submitted certificate.

Verified exactly:
- 156 graph isomorphism classes and 56 tournament isomorphism classes.
- 78 complement-paired graph classes and 34 converse-paired tournament classes.
- Every entry of the submitted 34 by 78 integer profile matrix.
- Row sums 720, 34 pairwise distinct rows, and exact rational matrix rank 32.
- The submitted signed kernel vector has zero total mass and is annihilated by every profile coordinate.
- Its positive and negative parts have equal mass 122880.
- The submitted 9-support and 12-support probability measures are disjointly supported, each sums to one, and have equal observer profiles in all 78 coordinates.

`audit_walsh.py` propagates signs through adjacent-transposition graph orbit edges, independently of permutation enumeration. A sign conflict detects a forbidden stabilizer sign. For every submitted row for n=1,...,6 it verifies the representative, edge list, stabilizer order, orbit size, signed order sum, and derived rank.

Verified dimensions: 1,1,2,3,10,34.
Verified ranks: 1,1,2,3,10,32.
Verified n=6 zero representatives: 255 and 2015, and no others.

## Normalization

Let m=|Orb(X)|. With U_X=sum_Y chi_Y and S_X=sum_Y epsilon_Y chi_Y, both over distinct graph orbit elements, alpha_X=(sum_Y epsilon_Y)/m=(sum_p epsilon(p,X))/n!. Then

P_T chi_X=S_X/m,
P_T chi_Y=epsilon_Y S_X/m,
P_T U_X=alpha_X S_X,
P_G S_X=alpha_X U_X.

The same formulas hold when U_X and S_X are normalized by the same factor, including division by sqrt(m). They would need an orbit-size factor if one used a permutation sum for only one basis.

The invariant-function kernel consists of the S_X with alpha_X=0. When translating this to class-probability mixtures, use coefficients w(C)=|C| S_X(C), up to a common scalar. The submitted direct checker uses this weighting correctly. Unweighted S_X evaluations generally are not the row-matrix mixture kernel.

The complement-even condition is exactly |E(X)| even. Alpha changes sign under a change of labeled representative gauge, but its zero/nonzero status and squared singular value are intrinsic.

## Connected zero family

With centers a,b first, r common neighbors next, and the pendant last, conditioning on independent uniform coordinates x,y of a,b gives the exact average

integral_[0,1]^2 sign(y-x)(1-2x)(1-2|x-y|)^r dx dy
=2 integral_0^1 d(1-d)(1-2d)^r dd.

This is zero for odd r by d -> 1-d. For odd r>=3 the graph is connected, its two centers have distinct degrees r+2 and r+1, and its automorphism group is S_r on the degree-two common neighbors. Each such permutation acts evenly on edge orientations: swapping common neighbors affects the two identical center-neighbor edge blocks and contributes +1. The pendant is fixed. Adding isolated vertices preserves the signed average, since relative order on nonisolated vertices is uniform, and preserves the Royle-even property. In particular the r=3 example yields a nonzero invariant-function kernel for every n>=6.

## Distinctions and cleanup

Rank 32 in dimension 34 does not imply two individual tournament classes have equal profiles. The verified matrix has 34 distinct rows. The demonstrated loss is linear/mixture identifiability, not individual-class identifiability at n=6.

`analyze.py` currently emits `"graphs": 0` through an obvious placeholder expression. Remove this unused field or populate the correct graph counts 1,2,4,11,34,156. All other submitted Walsh fields match the independent computation.

The submitted `direct_check.py` imports `/tmp/session58-complement-tournaments/verify.py`. Vendor or otherwise package this dependency before publication, as already planned by the author. The independent audit scripts require only the Python standard library.

## Quadratic extension audit

Verdict: PASS. `audit_quadratic.py` reconstructs the selected 34 by 34 minor entrywise from the literal one-order matrix and the stored feature labels, then uses ordinary rational Gaussian elimination rather than the source Bareiss determinant implementation.

The selected 32 linear columns, together with squared columns indexed by graph masks 3 and 15, have determinant

211991592542375584355844096000000,

whose residue modulo 1000000007 is 967372386. Both values match exactly. The mask labels denote:
- 3: edges 01 and 02, hence P3 disjoint union 3K1 on six vertices.
- 15: edges 01,02,03,04, hence K1,4 disjoint union K1.

In each case the observable event is the graph-complement pair containing the named graph, not only that graph without its complement.

Let p_C(q)=A_C(q)/720. Drawing a single latent converse class C, and then two independent uniform vertex orders of that same tournament, gives conditional pair law p_C(q)p_C(r). The two added diagonal events are p_C(3)^2 and p_C(15)^2. One-order features are marginals of this pair law, so the verified minor proves that the two-order shared-latent observation identifies every probability mixture over the 34 n=6 converse classes. One order does not, by the earlier kernel certificate; therefore two orders are minimal for this mixture-identifiability property at n=6.

The determinant for the normalized 32-linear-plus-two-quadratic feature matrix is the displayed integer determinant divided by 720^36, so normalization cannot change rank.

For the submitted positive and negative mixtures, positive minus negative probabilities of the two diagonal events are exactly -1/57600 and +1/230400. These are obtained by dividing the raw weighted contrasts -1105920 and 276480 by 122880 * 720^2.

Important qualification: independently redrawing the latent class for each of the two observations gives the product of the one-order averaged laws and does not distinguish these mixtures. Squaring an already aggregated one-order probability likewise does not supply the additional information. The two orders must share the same latent class.
