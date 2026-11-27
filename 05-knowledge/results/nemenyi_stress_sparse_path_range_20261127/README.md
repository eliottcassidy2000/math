# Strong Moon geometry, sparse Hamiltonian-path heights

2026-11-27. Follow-up to the [Neményi–Moon cycle flags](../nemenyi_moon_blue_cuts_20261127/README.md) and the [third ternary chain](../chanpak_third_ternary_chain_shadow_20261127/README.md). Here a sample is the nondegenerate tournament obtained by comparing the determinants of actual planar deviators. This is a restriction on a realization model, not a restriction on the unrestricted Chan–Pak spectrum.

## 1. The sharp stress floor, with equality case

Every stress sample is coded by a binary word. Sort projective ray representatives in one open half-circle; a bit says whether the actual deviator points in that direction or the opposite direction. For i<j the arc is forward iff the two bits agree. Conversely every word is a stress sample: use Z=z²/2 on the right half-plane, λ=|z−2|²−1/4 and sites z_i=x_i(1+ε i√−1), with x_i=2 or 8 for negative or positive signs and ε>0 sufficiently small. Thus varying the sample sites in a single field suffices. Words with at most two runs are transitive; three successive runs supply a directed triangle and the angular insertion argument gives a spanning cycle.

[Grinberg–Stanley, Theorem 1.38 and its path-count specialization](https://arxiv.org/html/2307.05569) express H as a sum over collections of disjoint nontrivial directed odd cycles, with each collection of size j weighted 2^j; the empty collection contributes one. The exact Moon bound for canonical angular cycles gives at least binom(n−2,k−2) cycles at each odd length k≥3. Keeping one-cycle collections gives

    H ≥ 1 + 2 Σ_(k odd, 3≤k≤n) binom(n−2,k−2) = 1+2^(n−2).

It is attained by reversing the extreme arc of a total order: the paths using it are indexed by all subsets of middle vertices placed before the reversed arc; only the increasing path avoids it. This is itself a stress sample, the cyclic triple with two singleton blocks and one transitive block of size n−2.

In fact equality is unique up to isomorphism. It requires the triangle count to equal n−2, so the convex degree comparison in the preceding Moon note is tight for binom(d,2). For this function every nontrivial integer hinge has coefficient 1, and the hinge sums determine the degree multiset; any different degree multiset would have strict inequality at some hinge. For n≥4 the extremal multiset has two vertices v,w of outdegree n−2. Let v→w. The only in-neighbor of w is v, and the only in-neighbor of v is some u≠w. The vertices other than v,w are in the transitive out-neighborhood of w and hence are totally ordered. Both v and w beat all of them except that u beats v. If u were not last in that order, the tail after u would have no outgoing arc to its complement, contrary to strong connectivity. This gives the extreme-arc reversal. The n=3 case is the directed triangle.

## 2. A second cost: many runs force factorial growth

For a word with r runs, choose one site per run and, if necessary, drop the last to obtain an odd alternating subword of length k≥r−1. Each vertex has outdegree (k−1)/2: count later sites of its own bit and earlier sites of the other bit. Its intrinsic angular cycle therefore gives the regular carousel in which each vertex beats the next (k−1)/2 sites.

For odd k≥25, let m=floor(k/12) and take three intervals of m cyclic positions starting at 0, floor(k/3), floor(2k/3). A forward difference from one interval to the next lies between floor(k/3)−m+1 and ceil(k/3)+m−1, hence within 1 through (k−1)/2. All three sets of cross-arcs thus point cyclically. Two independent bijections from the first interval to the other two give (m!)² distinct collections of m disjoint directed triangles. Consequently H≥2^m(m!)², so log H≥c r log r for all sufficiently large r and some absolute c>0.

For H≤X the sharp floor forces n=O(log X), and the cycle collections force r=O(log X/log log X). If L and R are these upper bounds, the count of words is at most 2L Σ_(j≤R) binom(L,j), by recording the initial sign and the run boundaries. Taking logs bounds the number of distinct path-count values by

    exp(O(log X · log log log X / log log X)) = X^{o(1)}.

This is the obstruction hidden by the pancyclic result: abundant cycle collections make many candidate words too expensive to contribute *small* exact path totals.

## 3. Ordinal joins remain sparse, but have many distinct products

Write S for the strong sample heights above 1. Ordinal joins concatenate paths, so their heights are exactly products of elements of S. For every ε>0 the single-sample bound and a dyadic sum give Σ_(h∈S)h^(-ε)<∞. The product Π_(h∈S)(1−h^(-ε))^(-1) converges. It counts factor multisets with their weights, and thus majorizes the sum of n^(-ε) over the distinct product heights. Therefore the number of joined heights ≤X is O_ε(X^ε) for every ε>0.

There is also a nontrivial lower bound for joined heights, even using only the floor realizers. For m≥4 and m<k≤2m, [Bang–Zsigmondy](https://doi.org/10.1007/BF01692444) gives a prime of order 2k for base 2, hence a divisor of 2^k+1. It divides no 2^j+1 for other j in that interval: that would require k|j. The 2^m subset products of these floor heights are thus all distinct, and their maximum is ≤2^((2m+1)m). Choosing m on the scale sqrt(log X) shows that the joined range has at least exp(Ω(sqrt(log X))) values through X. This does not give positive density or reach all large odd numbers.

## 4. Chain 7 and the finite boundary

The exact saturated subset-and-last-vertex DP checks all locally transitive strong classes of orders 3–13 using the cyclic degree word modulo rotation, with one-word-per-complement. The sharp floor rules out n>14 for X=4097, but **order 14 itself contributes 4097**, absent in the orders through 13. This corrects omitting the endpoint in an order-13-only census. Including height 1 gives **76** single-sample heights and **326** product heights through 4097. Neither list contains any of 7,21,63,189,567,1701, the entire 7·3^j chain in this window. By the third-chain result, 63 is attained in the unrestricted spectrum. Of the two flanking triples (19,21,23) and (61,63,65), the single-sample model attains only 23 and 65. The exact list is in replay.out.

For the simpler palette of cyclic triple blow-ups with a singleton block, there is an all-exponent omission. With block sizes a,b,1, the count is (2^a−1)(2^b−1)+2: three block visits give 3 paths, four visits contribute (2^a−2)+(2^b−2), and five visits contribute (2^a−2)(2^b−2) by splitting both nonsingleton blocks. To hit 7·3^j at j≥1, reduction modulo 3 forces both a,b odd. Then modulo 8 the path count is 3 if both ≥3, 1 if exactly one is 1, and 3 if both are 1; the target is 5 or 7. At j=0 the product of the two Mersenne factors would have to be 5, impossible. Even ordinal joins of this palette omit the entire chain, since a factor carrying the single factor of 7 would itself be 7·3^t. No such all-exponent conclusion is asserted for the full stress family.

Run **node replay.js** for the finite order census, exact product closure, and the carousel interval checks.
