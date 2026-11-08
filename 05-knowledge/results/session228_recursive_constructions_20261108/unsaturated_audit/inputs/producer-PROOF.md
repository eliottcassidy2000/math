# Unsaturated modules give exponential strong gcd-one fibers

Status: exact all-order construction; no claim of prime coverage or consecutive scalar values.

## The theorem
For every integer t>=1 there is a single physical rooted-tournament sidecar (E,U,d) with at least 2^t distinct Hamilton-path counts D. All its 2^t displayed witnesses are strong tournaments on 7t+2 vertices. The gcd of those displayed counts is 1. Explicitly,

 d=E=19^t,   U=(700*512^t-207*19^t)/493.

Take the six-vertex block B with outgoing adjacency masks [40,1,3,6,15,30]. Its Hamilton count H is 19; its ordered two-path-cover count Q, empty paths allowed, is 256. Its vertex-first two-path-cover weights are [21,20,29,35,73,59]. Root outneighbor subsets of masks 17={0,4} and 40={3,5} both have weight alpha=94. They define modules M0 and M1 on the SAME underlying block.

Delete the root. In forward block order put a singleton beaten by the root, t chosen modules separated by singletons beaten by the root, and finally a singleton beating the root. All arcs between different blocks point forward. Root arcs in each module are given by its chosen mask. Thus there are t+1 singletons and 6t module vertices, plus the root.

Every vertex is reachable from the root through the first singleton; every vertex reaches the root through the last singleton. These statements also cover the first and last singleton themselves, proving strongness, irrespective of any internal block reachability.

For choices epsilon_i in {0,1}, indexed i=0,...,t-1 in forward order, the scalar is

 D(epsilon)=D0(t)+sum_i epsilon_i delta_i(t),
 D0(t)=[19950*512^t+(13862*t-2289)*19^t]/15979,
 delta_i(t)=-[8400*512^i*19^(t-1-i)+1460*19^(t-1)]/493  (i<t-1),
 delta_(t-1)(t)=-[8400*512^(t-1)+474*19^(t-1)]/493.

All displayed divisions are exact because these expressions equal integer path counts or their differences, as proved below.

## Transfer convention and universal separator lemma
Use column vectors and six states (p,r)=(0,0),(0,1),(1,0),(1,1),(2,0),(2,1). The first coordinate records an empty P, or P ending at an inneighbor/outneighbor of the root. The second records an empty R or R already starting at an outneighbor. For a local ordered pair (P,R) covering a nonempty block, the eight endpoint categories (pp,rr) belong to {0,1,2}^2 minus (0,0). Here rr=1/2 means its first vertex is an outneighbor/inneighbor. From old state (p,r), discard a cover if r=0 and rr=2; otherwise update to (pp if pp!=0 else p, bool(r or rr)). The matrix entry is the number of covers making that transition. Singleton transfers are S0,S1, with bit1 meaning root beats the singleton.

Let v=(1,0,0,0,0,0)^T and ell=(1,1,1,1,0,0). Reading blocks forward applies their matrices right to left; D=ell M_last...M_first v.

Let W consist of vectors supported on coordinates 1,3,5 whose coordinate sum is zero. Let V consist of vectors satisfying z0=0, z2+z4=0, z1+z3+z5=0. Thus V is spanned by e1-e5,e3-e5,e2-e4, and W by the first two.

Every arbitrary physical nonempty block transfer M preserves V and acts on W as H times the identity. To see preservation, check each of the eight elementary category transfers on the three given basis vectors, so their weighted sum preserves V. For W, local covers with P nonempty erase the old p coordinate, and their contributions cancel because the input sums to zero. Covers with P empty preserve p and there are exactly H of these, giving H times the input. Each singleton S0 or S1 maps V into W: the P-nonempty covers erase the old p and cancel, and its P-empty contribution has R nonempty.

If two blocks or weighted blocks have equal H,Q,alpha, their difference F has category weights
 (0,1):x, (0,2):-x, (1,0):y, (2,0):-y,
 (1,1):a, (1,2):b, (2,1):-x-a, (2,2):x-b.
Equality of H gives the first two pairs, equality of alpha and Q gives the last two. Direct application of the update rule shows F maps into V and kills W. Consequently F_left C F_right=0 whenever C is any product of physical block transfers containing at least one singleton, at ANY position in that product. Before that singleton the image stays in V; afterward it stays in W, which the left defect kills. No binary-compatibility hypothesis is needed. Empty separation is not a valid universal claim: the checker includes a counterexample from the universal defect basis.

This kills all terms involving at least two chosen defects in our expansion, proving exact additivity. It also explains why arbitrary nonbinary baselines remain valid. It does NOT force the old product/lcm divisor: our baseline block has Q/H=256/19, not an integer binary power.

## Exact scalar certificates
The checker derives M0 and M1 by literal permutations split into two directed paths, not by the claimed formulas. Set A=S1 M0, F=M1-M0. It verifies the matrix identity

 A(A-19I)^2(A-512I)=0.

The baseline is ell S0 M0 A^(t-1) S1 v. For an interior defect, after applying S1 F its vector belongs to W. Every later module multiplies that vector by 19; singletons act as identity. Therefore its contribution is 19^(t-1-i) times g_i, where

 g_i=ell S1 F A^i S1 v=-(8400*512^i+1460*19^i)/493.

The final defect instead has

 k_i=ell S0 F A^i S1 v=-(8400*512^i+474*19^i)/493.

The checker verifies four initial terms for each of the three sequences and the displayed annihilating matrix polynomial. The proposed formulas are linear combinations of 512^i,19^i,i19^i, hence obey the same fourth-order recurrence. This proves them at every nonnegative index by induction, and proves the scalar formulas above.

## Sidecar
For any ordered physical blocks the root-started two-path sidecar obeys
 U/d=1+sum_i (alpha_i/H_i) product_(j>i)(Q_j/H_j), d=product_i H_i.
This follows by selecting the first block used by the path starting immediately after the root. Here putting R=512/19 gives U/d=1+(700/493)(R^t-1), as claimed. E counts full deletion paths ending at a root inneighbor. Every deletion Hamilton path ends in the final singleton, which beats the root, so E=d. Choices 17 and40 have equal alpha and the same physical deletion graph, hence preserve this sidecar exactly.

## Exponential distinctness and gcd one
Write w_i=-delta_i>0. Consecutive w_i have ratio greater than20. Cancelling their common factor 19^(t-1)/493 reduces the interior comparison to
 8400*(512/19-20)*(512/19)^(i-1)-19*1460>0.
For the final comparison replace the final constant1460 by474, subtracting an additional986; the inequality remains strict already at i=1. Thus each w_i exceeds the sum of all earlier w_j. Distinct binary subsets therefore have distinct sums. This proves 2^t distinct scalars at order7t+2, with logarithmic vertex order in fiber size.

For t=1 the counts are653 and635, with gcd1. For t>=2, delta_0=-20*19^(t-1). The final delta is not divisible by19, since its numerator modulo19 is8400*512^(t-1) and493 is invertible modulo19. It is not divisible by5, since its numerator modulo5 is474*19^(t-1) and493 is invertible modulo5. It is2 modulo4 up to sign, since8400 vanishes modulo4 and474*19^(t-1) is2 modulo4, with odd denominator493. Therefore gcd(delta_0,delta_last)=2. Baseline D0 is odd: its denominator15979 is odd and its numerator has an even first term plus an odd second term. The gcd of the fiber divides both D0 and all its one-bit differences, so is1.

## Scope and negative controls
The strong order-four H5,Q32 representative [8,1,3,6] has weights [5,6,9,7]. All16 subset sums are distinct, so no pair of root masks on it supplies this same-alpha mechanism. The checker exhausts the64 labeled order-four tournaments and verifies the same absence for every H5,Q32 block. This is a bounded negative result, not a classification of higher-order modules.

The new positive theorem removes the common-divisor defect of separated binary-compatible repeated modules. It does not imply that any displayed D is prime, that fibers contain intervals, or that arbitrary target Hamilton counts occur. The two root masks produce physically labeled witnesses; isomorphism is irrelevant because distinct D already distinguish the resulting rooted tournaments.
