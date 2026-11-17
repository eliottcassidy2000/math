# First births at large order: explicit near-extremal families and typical behavior

For any odd h, its birth order is the least n at which a tournament of order n has h directed Hamilton paths. This note concerns **first-birth classes** (a class whose own count first appears at its order) and **distinct first-birth counts** separately. In particular it does not claim that most odd integers are attained or that finite first-birth atlases classify large orders.

Write M(n) for the maximum Hamilton path count at order n and μ_n=n!/2^(n−1). Padding by an ordinal singleton preserves the count, so M is nondecreasing and all counts above M(n−1) at order n are certified births.

## A literal family and a constructive net of distinct births

Use Divoux and Wang, [The Maximum Number of Hamiltonian Cycles in Tournaments](https://arxiv.org/html/2610.23755v1), Theorem 1.4 and its proof in §10. Put c=cosh(1)/cos(1). Their theorem proves M(n)=(c+o(1))μ_n and its proof explicitly attains this asymptotic with R_n: for odd n, orient i→j when the positive circular difference j−i mod n is less than n/2; for even n, delete vertex 0 from the corresponding circular tournament on n+1 vertices. Thus **R_n itself is a first-birth constructor at every sufficiently large order**, since μ_n/μ_(n−1)=n/2. The preprint theorem supplies an eventual statement, not a numerical cutoff here.

One can generate many more distinct counts by walking from the transitive tournament to R_n in lexicographic order of the differing labelled arcs. Under one arc reversal x→y, the path sets not using that arc are unchanged. Any lost path using x→y splits, after deleting the consecutive pair x,y, into two ordered optional paths on the other n−2 vertices. Keeping the pair orientation and the two paths reconstructs the original. The same bound applies to the newly gained paths. Therefore the absolute change is at most P₂(T−{x,y}), **not** the count of paths on the deletion tournament times insertion slots: the seam between the two surviving pieces need not be a legal arc. For m=n−2 the universal upper bound is

    D_n = sum_(s=0..m) binom(m,s) M(s)M(m−s),  M(0)=1.

By the cited asymptotic, D_n=(c²+o(1)) m·m!/2^(m−2)=(4c+o(1)) M(n−1): for s and m−s both large, each summand normalized by m!/2^(m−2) tends to c² uniformly; the bounded number of small-endpoint terms have bounded normalized size and contribute o(m). Formally first exclude a fixed K terms at each end, then take K→∞ in the uniform error bound.

For any integer threshold t with M(n−1)<t≤H(R_n), choose the first walk vertex whose count is at least t. Its count belongs to [t,t+D_n−1], and is a birth. Thresholds spaced by D_n give disjoint output intervals. In particular the initial first-birth tail has an explicit one-sided O(M(n−1))-net, and the number of distinct first-birth counts b_n satisfies

    liminf_(n→∞) b_n/n ≥ 1/(8c).

This follows by using the thresholds from M(n−1)+1 through H(R_n) and H(R_n)/M(n−1)=(1+o(1))n/2. To avoid knowing M(n−1) computationally, an entirely explicit eventual prescription starts just above 3μ_(n−1), with spacing ceil(40μ_(n−1)); c<3 and 4c²<40. Compute the walk counts by subset DP and stop at each certified threshold. No exhaustive census of earlier orders is required asymptotically.

All these births are strong for large n. A singleton SCC can be removed without changing H. Otherwise a nontrivial SCC cut has parts a,b≥2, and the ordinal product identity gives H≤C² a!b!/2^(n−2)=O(M(n−1)/n), since a!b!≤2(n−2)! and some uniform constant C bounds M(k)/μ_k. This excludes decomposable tournaments above the previous-order maximum.

## A robust deterministic neighborhood of birth classes

Fix 0<η<1/(8c). Every tournament differing from R_n by at most floor(η n) arcs has H>M(n−1) eventually, by the same edge bound (total loss at most η n D_n). This ball contains exp(η n log n−O(n)) **isomorphism classes**, not merely labelings. Indeed if T,T' in the ball are isomorphic via σ, the permutation satisfies d(R_n,σR_n)≤2η n. Since 4η<1, this leaves a clean vertex incident to no disagreement. Its incoming and outgoing neighborhoods in R_n are transitive (also after the even-order deletion). Thus, upon choosing a clean anchor and its image in n² ways, σ must map the two neighborhoods onto the corresponding neighborhoods and their total inversions within the transitive orders are at most 2η n. Each permutation is specified by its inversion code; the number of pairs of such codes with the total bound is at most the number of weak compositions of ≤2η n into n−1 entries, which is exp(O(n)). Therefore each class meets the labelled edit ball in at most exp(O(n)) labelled tournaments, whereas the ball has at least binom(binom(n,2),floor(η n)) labelled tournaments. This gives a fully specified structured subset of the birth classes, regardless of the typical behavior below.

## Almost every tournament is already a first-birth representative

Here is an exact second-moment identity. In a uniformly random labelled tournament let H be its count, so E H=μ_n by Szele's expectation argument. Fix one ordered path, and let the other range uniformly over permutations. For any arc of the fixed path, the normalized joint expectation has factor 2 if the other path uses it in the same direction, 0 if the opposite direction, and 1 if it does not use the unordered pair. Expand the product using these factors in the form 1 plus a signed adjacency indicator. For a subset F of the n−1 ordered arcs, the selected consecutive runs must each appear as a consecutive block in the other path, forward or reversed. For a run of l arcs the two signs sum to 1+(-1)^l. After contracting the blocks there are (n−|F|)! possible permutations. Thus exactly

    E[H²]/μ_n² = sum_(F with all runs even) 2^{number of runs(F)} (n−|F|)!/n!.

A subset with 2k selected arcs is uniquely tiled by k disjoint dominoes on the n−1 positions, so there are binom(n−1−k,k) such F and at most k runs. The total term at size 2k is at most

    2^k binom(n−1−k,k) (n−2k)!/n! ≤ (4/n)^k/k!.

The k=1 term is exactly 2(n−2)/(n(n−1)). Hence Var(H)/μ_n²=2/n+O(n^-2), with all higher terms nonnegative and O(n^-2). In particular M(n−1)/μ_n=O(1/n) shows by Chebyshev that the fraction with H≤M(n−1) is O(1/n). Thus **all but O(1/n) of the labelled order-n tournaments are first-birth examples** (the over-threshold examples are even strong).

The same statement holds for uniform isomorphism classes, without assuming all orbits have maximum size. In Burnside's sum a permutation with an even cycle cannot fix a tournament: its half-power reverses an opposite pair. If only odd cycles occur and support size is s≥3, every moved unordered pair has orbit size at least three. There are s(n−s)+binom(s,2) moved pairs (a fixed pair would have two fixed vertices or a transposition). The proportion fixed by this permutation is at most 2^{-s(n−1)/3}. There are at most n^s permutations with support s, so the nonidentity part of Burnside's sum, normalized by the total labelled count, is O(n³2^{-(n−1)}). Apply the same sum to just the isomorphism-invariant bad set to obtain the O(1/n) result on classes. This says that almost all tournament structures are minimal witnesses for their own counts, despite little control on which odd counts those are or how many first-birth structures share one.

## Replay and scope

The earlier `birth-atlas-now` finite atlas describes exact small-order birth fibres and singleton padding. Here `replay.py` independently enumerates every labelled tournament through order five to check the moment identity, compares the circular endpoint to earlier maxima when the local atlas is available, and checks the corrected two-cover edge bound along literal walks through order ten. For example the order-nine circular tournament has 3267 paths versus the earlier maximum 661 from that atlas; the chosen lexicographic walk witnesses five distinct counts above 661. The asymptotic theorems do not rest on the replay range.

## What the typical count reads: the outdegree-score imbalance

The second-moment identity has a useful structural refinement, not just a concentration conclusion. Let s_v=2d⁺(v)−(n−1) be the signed row sum of the tournament's Seidel matrix. Then

    H/μ_n = 2 − (Σ_v s_v²)/(n(n−1)) + R_err,
    E[(R_err)²] = O(n^(−2)).

The symbol here denotes a random remainder, not the circular tournament.

To see it, let X_ij=+1 if i→j and −1 otherwise, so X_ji=−X_ij and the signs on unordered edges are independent. For each ordered path its normalized count contribution is the product of (1+X_ij) over the consecutive arcs. Average this product over permutations. The sum of all singleton Fourier terms vanishes by reversing the corresponding block. For degree two, a pair of disjoint edges likewise vanishes because each block can be reversed independently. A two-edge wedge i,j,k, on the other hand, has two block orientations with the same product, giving exactly

    2/[n(n−1)] · Σ_j Σ_{{i,k} excluding j} X_ij X_jk.

The sum on the right is [n(n−1)−Σ_j s_j²]/2 by expanding each square and using antisymmetry. Distinct edge-set monomials are orthogonal under independent signs; hence the variance of this degree-two term is exactly the k=1 summand of the preceding second-moment calculation. The orthogonal higher terms have total squared norm O(n^-2). Thus the imbalance in scores explains a proportion 1−O(1/n) of the variance in the typical count. It is an average statement, not valid pointwise for near-extremal circulants. The finite replay exhausts order 2 through 5 and checks both this squared residual and its orthogonality.

### A quantitative tail via the Fourier layers

In fact the exceptional fraction is at most exp(−Ω(sqrt(n))). Let f_k be the homogeneous Fourier layer of degree 2k in H/μ_n−1; every odd layer vanishes by the same signed block reversals. The exact second-moment computation by size gives `||f_k||₂²≤(4/n)^k/k!`. For completeness, the grading can be seen by taking two correlated random tournaments with correlation t on each unordered edge and independent between edges: their normalized joint expansion replaces each signed shared adjacency by t times that sign, so the size-2k contribution is exactly the coefficient of t^(2k), which is the squared layer norm by Fourier orthogonality.

Use the Boolean cube hypercontractive inequality of [Bonami](https://aif.centre-mersenne.org/articles/10.5802/aif.357/), in its degree-layer form `||f_k||_p≤(p−1)^k ||f_k||₂`. At real p=1+sqrt(n)/20 (at least 2 eventually) we get

    ||H/μ_n−1||_p ≤ Σ_(k≥1) (p−1)^k (2/√n)^k/√(k!) ≤ Σ_(k≥1) (1/10)^k = 1/9.

Therefore the probability of |H/μ_n−1|≥1/2 is at most (2/9)^p. The earlier threshold M(n−1)=o(μ_n) and the exponentially smaller Burnside correction transfer the bound to the fraction of all order-n isomorphism classes that fail to be first births. The elementary O(1/n) result above does not require hypercontractivity.
