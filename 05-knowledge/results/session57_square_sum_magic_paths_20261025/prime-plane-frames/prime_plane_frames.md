# Unique-sum affine-plane frames over a prime subfield

## Statement and scope
Let p be an odd prime, q=p^n with n>=2, G=(F_q,+), U an F_p-linear subspace of dimension two, and A=C+U. A frame is an unordered pair {B,D} of p-element subsets of G such that the addition map B x D -> A is bijective. No disjointness of the two vertex sets is imposed until expressly stated.

Every frame has an affine F_p-line as at least one factor. Its other factor is an arbitrary transversal to that line direction inside the appropriate affine U-plane. Consequently, with

C_p=(p+1)p^(p-1)-p(p+1)/2,

there are exactly C_p q frames, of which p(p+1)q/2 have two line factors and [(p+1)p^(p-1)-p(p+1)]q have exactly one line factor. Exactly C_p(q-p^2) frames have B and D disjoint. All these counts are for one fixed A, not summed over affine planes.

For p=3, C_3=30, with 6q line/line (AP/AP) and 24q line/nonline frames. The corresponding disjoint counts are 6(q-9) and 24(q-9). For p=5, C_5=3735, split as 15q and 3720q; the disjoint count is 3735(q-25).

## 1. Reduction to a direct factorization of U
Choose b_0 in B and d_0 in D. Differences b-b_0 and d-d_0 belong to U, since they are differences of elements of A with the other summand fixed. Thus X=B-b_0 and Y=D-d_0 are p-subsets of U. Since b_0+d_0 lies in A, bijectivity of addition is exactly X plus Y = U with unique representations.

Equivalently (X-X) intersect (Y-Y) = {0}. For p>3, one must NOT simply replace nonzero differences by their directions in this equivalence: a single difference pair supplies only plus/minus v, not all p-1 multiples. The following lemma repairs that gap using the full factorization.

## 2. Cyclotomic projection lemma
For every nonzero F_p-linear functional ell:U -> F_p, at least one of ell restricted to X and ell restricted to Y is bijective.

Proof. Put zeta=exp(2 pi i/p). Unique representations give

(sum over x in X of zeta^(ell(x))) (sum over y in Y of zeta^(ell(y)))
 = sum over u in U of zeta^(ell(u)) = 0.

Hence one factor sum vanishes. For that factor S, let n_j be the number of s with ell(s)=j, 0<=j<p. The integer polynomial P(T)=sum n_j T^j has degree at most p-1, P(zeta)=0, and P(1)=p. The minimal polynomial of zeta over Q is Phi_p(T)=1+T+...+T^(p-1). Thus P is a rational constant multiple of Phi_p, and evaluation at 1 makes that constant 1. Every n_j is 1, as required.

In geometric terms, for each direction L=ker ell, at least one factor meets every affine L-line exactly once, so does not determine direction L. Therefore the direction sets determined by X and Y are disjoint. This conclusion uses primality, exact p-point factors, and complete unique coverage of U.

## 3. A line factor is forced
The Redei-Megyesi direction theorem states that a noncollinear p-point subset of F_p^2 determines at least (p+3)/2 directions. A verified modern primary proof is Gabor Somlai, “A new proof of Redei's theorem on the number of directions,” Archiv der Mathematik 122 (2024), 575-580, Theorem 1.1, DOI https://doi.org/10.1007/s00013-024-01979-x . The publisher's introduction identifies the original result as Redei, Lacunary Polynomials over Finite Fields (English translation, 1973), Theorem 24', p.226. The proof here uses the stated theorem as its only finite-geometric external input; it does not claim this direction theorem as new.

If neither X nor Y were collinear, their disjoint direction sets would together contain at least p+3 directions, exceeding the p+1 directions of U. Hence one factor is collinear. A collinear p-point set is a full affine F_p-line. Translation transfers the conclusion to B or D.

## 4. Count all unordered frames
Count ordered frames (B,D) whose first factor is a line. Choose its direction L<=U in p+1 ways, and B=b+L among q/p affine L-lines of G. The second factor must be contained in A-b, an affine U-plane, and must select exactly one point from each of the p parallel L-lines there. Conversely every such selection gives a bijective sum with B. There are p^p selections. Hence the line-first count is

T=(p+1)(q/p)p^p=(p+1)p^(p-1)q.

For two line factors, their directions must be distinct. For each ordered direction pair (L,M), choose B among q/p lines, then D among the p M-lines inside A-b. Thus the ordered two-line count is p(p+1)q and the unordered two-line count is

H=p(p+1)q/2.

The line-first count represents a mixed unordered frame once and a two-line unordered frame twice. Since every frame has a line factor, the required total is T-H=C_p q. The mixed count is T-2H. No frame has B=D because its sum representations would fail uniqueness.

## 5. Gauge classes and disjointness
For an ordered frame define T_t(B,D)=(B+t,D-t), t in G. Addition and A are unchanged. This G-action is free: if B+t=B and D-t=D for t!=0, then t occurs as a nonzero difference in both B-B and D-D, contradicting uniqueness.

Swapping factors does not identify an ordered gauge orbit with itself. Indeed, if B+t=D then B and D are translates, so their nonzero difference sets agree, again impossible. Consequently the group generated by translations T_t and the swap acts freely on ordered frames, with orbits of size 2q. On unordered frames this defines gauge equivalence classes of size q. Strictly speaking T_t is not intrinsically a G-action on unordered pairs without choosing an orientation: reversing orientation replaces t by -t. The equivalence classes and their size are intrinsic, and this distinction avoids an unnecessary well-definedness error.

The map B x D -> G given by (b,d) -> d-b is injective: equal differences yield b+d'=b'+d, so addition uniqueness forces both pairs equal. Thus |D-B|=p^2. Now

(B+t) intersect (D-t) is nonempty iff 2t belongs to D-B.

Since p is odd, multiplication by 2 is bijective on G. Exactly p^2 of the q distinct representatives in each unordered gauge class intersect, and exactly q-p^2 are disjoint. There are C_p classes, proving the disjoint formula. The line/line property is invariant in a class, so each subtype's q-factor can also be replaced by q-p^2.

## 6. Characteristic two is different
If p=2, every two-point set is a line. The count formula for all frames still gives C_2=3 and 3q total frames. But 2t=0, so the odd-characteristic disjointness argument does not apply. In fact B and D intersect iff 0 is in A: in characteristic two a representation 0=b+d is exactly b=d. Thus all 3q frames are disjoint when 0 is not in A, and none are disjoint when 0 is in A. Each intersecting frame has exactly one common point. This is not 3(q-4) in general.

## 7. Verification and application guards
The accompanying standard-library Python script exhausts all C(25,5)=53,130 five-point subsets of F_5^2. Direction histogram: 30 sets determine 1 direction, 1,500 determine 4, 15,600 determine 5, and 36,000 determine 6. This verifies the sharp bound for p=5 independently; it is not a proof for all primes and not a full p=5 frame census.

It also exhausts all possible two-/three-point frame factors inside U-cosets for p=2,3 in dimensions 2 and 3, checks distinct differences and disjoint directions directly for every unique-sum frame, and verifies all frame and disjoint formulas. Results are in verification.json.

These are additive statements over F_q, equivalently over an elementary abelian p-group. Computation must not substitute Z/qZ when q is a nonprime prime power. The dimension-two U is over the prime subfield, not a plane over a larger subfield. The conclusion is not established here for larger factor sizes or nonprime subfield orders. No assumption concerning squares or a graph is needed for the additive theorem. To interpret disjoint frames as injective K_{p,p} graph embeddings, separately verify the graph's adjacency condition for every sum in A (and any required exclusion of zero); the present proof alone supplies no such graph condition.

## 8. Cross-skeleton Hamilton paths and a small-size obstruction
For all frames, the line/line fraction is

[p(p+1)/2] / C_p = 1/(2p^(p-2)-1).

It is the same fraction among disjoint frames when q>p^2; when q=p^2 there are no disjoint frames, so that conditional fraction is undefined. The values are 1/5 at p=3 and 1/249 at p=5.

A disjoint frame determines a complete bipartite cross-skeleton K_{p,p} on B union D, with edge sum labels bijectively equal to A. It has (p!)^2 unoriented Hamilton paths: order each side in p! ways, choose the initial side in two ways, then divide by reversal. Every such path has a unique unordered alternating bipartition, which recovers {B,D}. Thus no path is counted from two different frames. The number of unoriented Hamilton paths arising in these prescribed-plane complete cross-skeletons is exactly

(p!)^2 C_p (q-p^2).

The line/line contribution is (p!)^2 p(p+1)(q-p^2)/2, and the complementary contribution is (p!)^2 [(p+1)p^(p-1)-p(p+1)] * (q-p^2). These are statements about cross-skeleton paths. They do not assert that every path with sum labels in A extends to such a frame, or that an unrelated ambient graph has all these cross-edges. If the ambient graph includes every edge whose sum lies in A, each cross-skeleton is a subgraph; inducedness is not required.

Each Hamilton path uses 2p-1 edges. Two paths therefore use at most 4p-2 distinct edges/labels, strictly fewer than p^2 for every prime p>=5. Hence two Hamilton paths cannot cover all p^2 unique cross-edge labels for p>=5. The characteristic-three nine-label two-path phenomenon cannot be extrapolated to this larger p regime. In particular the general line-factor theorem should not be described as a general AP magic-square certificate.

## 9. Why triples have a special permutation symmetry
Fix an affine F_p-line L=b+F_p v, with v nonzero. Its bijective affine parameterizations F_p -> L are t -> b+(a t+c)v, where a is nonzero and c is arbitrary. Relative to one parameterization, these form AGL(1,p), a subgroup of S_p of order p(p-1), hence index (p-2)!.

Thus at p=3, AGL(1,3)=S_3: every permutation of a three-point affine line is an affine parameterization. At p=5 only 20 of the 120 permutations, a fraction 1/6, are affine parameterizations. For two fixed line factors the fraction of ordered parameterization pairs that are affine on both sides is 1/((p-2)!)^2. These parameterizations use the cyclic field index t in F_p; they are not claims about increasing real arithmetic-progression orders or real betweenness.

For triples this explains the familiar 72 ordered presentations of one unordered two-line frame: two choices of side order and (3!)^2 choices of point orders, all affine. More generally one unordered two-line frame has 2(p!)^2 fully ordered presentations but only 2[p(p-1)]^2 presentations affine on both sides. The equality of these two counts is special to p=2,3, and the nontrivial odd-prime case is p=3. This is an elementary permutation-group connection, not an additional representation-theoretic classification.
