# Growing Hall deficits from two-base carry-free bipartitions

Status: PROVED, pending independent review. This quantitative addendum to the all-order HYP-9172(a) obstruction does not alter the frozen proof packet and does not address HYP-9172(b).

## 1. Results

Let D_T(n) be the number of tournament types minus the size of a maximum matching in the stabilizer-containment bipartite graph between tournament types and sign-even graph types. Let D_C(n) be the analogous deficit for tournament switching classes and sign-even Eulerian graphs. A positive deficit forbids a labelled S_n-equivariant map inducing a type bijection.

For all integers a,b>=1, put

    n=3^a 5^b-1,     M=min(3^a,5^b).

Then

    D_T(n) >= M,     D_C(n) >= (M+1)/2.

In particular, n_r=15^r-1 gives D_T(n_r)>=3^r and D_C(n_r)>=(3^r+1)/2.

There is a stronger explicit square-root subsequence. For each b>=1 choose a_b least with 3^(a_b)>=5^b, and put n_b=3^(a_b)5^b-1. Then

    D_T(n_b) > sqrt((n_b+1)/3),
    D_C(n_b) > (1/2) sqrt((n_b+1)/3).

Thus both deficits are unbounded, and grow at least on a square-root scale along an explicit infinite subsequence. These are lower bounds, not exact deficits or bounds asserted at every order n. The certificates use complete bipartite graphs and simultaneous carry-free binomial coefficients in bases 3 and 5.

## 2. Exact inherited input

The published source is `05-knowledge/results/uniform_equivariance_obstruction_20261020/UNIFORM_EQUIVARIANCE_OBSTRUCTION.md`, section “Corollaries: a joint-prime obstruction and an explicit Hall family.”

For either source/target pair, write A_p and B_p for the type subsets whose stabilizers contain a full Sylow-p subgroup of ambient S_n. That result proves, for every n>=5:

1. |A_3|=|B_3| and |A_5|=|B_5|, by the one-prime Walsh count identity.
2. A_3 and A_5 are disjoint, by the single-3-cycle support argument.
3. The neighborhood of A_3 union A_5 in the stabilizer-containment graph lies in B_3 union B_5.

Consequently

    D(n) >= |B_3 intersection B_5|.

The addendum constructs many distinct target types in this intersection. “Full Sylow” means p-adic stabilizer order v_p(n!), not merely the presence of an element of order p.

## 3. Complete bipartite targets

For 0<=k<=floor(n/2), let X_k=K_{k,n-k}, with X_0 the empty graph on n vertices. Distinct k give distinct graph types: positive k give different unordered bipartition sizes recoverable from the connected complete bipartite graph; X_0 alone is empty.

For 0<k<n-k, Aut(X_k)=S_k x S_(n-k). At k=n-k>0 its order is 2(k!)^2, because the parts can also be swapped. At k=0 its order is n!. Thus, for each odd p=3,5,

    v_p(n!) - v_p(|Aut(X_k)|) = v_p(binomial(n,k)).

The graph X_k therefore carries both full ambient Sylow subgroups exactly when gcd(binomial(n,k),15)=1.

Orient every bipartite edge from the first part to the second. Part-preserving automorphisms preserve every edge direction. When the parts have equal size k, their swap reverses k^2 edges, with sign (-1)^k. Hence X_k is sign-even precisely when the parts are unequal or their common size is even. The empty graph is sign-even.

For k>0, X_k is Eulerian exactly when both part sizes are even: their vertex degrees are n-k and k. The empty graph X_0 is Eulerian for every n, including odd n, and must be retained separately in a general formula.

## 4. General computable bounds

For every n>=5 define

    b(n) = #{k: 0<=k<=floor(n/2), gcd(binomial(n,k),15)=1,
                 and either 2k<n or k is even},

and

    e(n) = 1 + #{k: 1<=k<=floor(n/2), gcd(binomial(n,k),15)=1,
                     k is even and n-k is even}.

Then

    D_T(n) >= b(n),     D_C(n) >= e(n).

Proof. The graphs counted by b(n) are distinct sign-even targets with both full ambient Sylow predicates. Those counted by e(n), including its separate empty-graph term, additionally have even degrees. Apply the inherited Hall-deficit inequality. QED.

For odd n this particular bipartite family gives e(n)=1. This does not say that the full switching deficit is bounded at odd orders. Balanced odd parts are excluded from b(n); balanced even parts are allowed in both bounds.

Examples of these lower bounds, not actual matching deficits, are:

    n       5   6   7   8   14   44   224   674
    b(n)    1   1   2   4    6   18    36    27
    e(n)    1   1   1   2    3    9    18    14

The general index count can be strictly larger than the initial subfamily constructed next. For example at n=14 the initial family has M=3, while all k=0,...,5 qualify for b(14)=6.

## 5. Prime-power construction and square-root growth

Fix a,b>=1, put n+1=3^a5^b and M=min(3^a,5^b), and take 0<=k<M. Since the larger factor is at least 3, n+1>=3M and 2k<=2M-2<n. All these graphs therefore have unequal parts and are sign-even.

For p=3 put t=a, and for p=5 put t=b. Every 1<=j<=k satisfies j<M<=p^t, hence v_p(j)<t=v_p(n+1). Therefore

    v_p(n+1-j)=v_p(j).

The product formula

    binomial(n,k) = product_(j=1)^k (n+1-j)/j

now gives v_3(binomial(n,k))=v_5(binomial(n,k))=0. There are M indices, proving D_T(n)>=M.

The integer n is even and M is odd. Exactly (M+1)/2 of k=0,...,M-1 are even; their two part sizes are even, making the corresponding graphs Eulerian. This proves D_C(n)>=(M+1)/2.

For the square-root subsequence, choose a_b minimal with 3^(a_b)>=5^b. Then

    5^b <= 3^(a_b) < 3*5^b.

Here M=5^b, and n_b+1<3M^2, so M>sqrt((n_b+1)/3). The claimed bounds follow. The sequence n_b is strictly increasing and tends to infinity because 5^b increases strictly and a_b is nondecreasing. QED.

## 6. Two-base digit interpretation

The binomial unit test means every base-3 digit of k is at most the corresponding digit of n, and the same is true in base 5.

For a self-contained verification, in F_p[x] the identity (1+x)^p=1+x^p gives, when n has base-p digits n_i,

    (1+x)^n = product_i (1+x^(p^i))^(n_i).

Each selected exponent has digits at most n_i, has no carry, and is uniquely represented. Its coefficient is a product of binomial(n_i,k_i), which is nonzero modulo p because 0<=k_i<=n_i<p. Thus binomial(n,k) is nonzero modulo p exactly when each k_i<=n_i. Equivalently, adding k and n-k in base p has no carry.

For n=3^a5^b-1 the last a base-3 digits and last b base-5 digits are maximal. Every 0<=k<min(3^a,5^b) fits within those unrestricted digits in both bases. This independently explains the binomial claim and identifies the arithmetic coordinate producing the growing Hall family.

## 7. Exact controls and scope

Run `python3 -B hall_growth_controls.py --out hall_growth_controls.json`.

The independent standard-library-only verifier:

- enumerates every permutation of every X_k through n=8 and checks full automorphism order, edge-reversal sign, Eulerian parity, and the ambient Sylow test;
- compares exact binomial divisibility and valuations with independently computed base-3/base-5 digit conditions for all 0<=k<=n through n=250;
- checks every initial index 0<=k<M for a,b=1,...,5, all diagonal cases a=b=1,...,8, and the first eight balanced-exponent square-root cases;
- evaluates exact integer binomial coefficients as an additional path whenever n<=5000;
- records b(n), e(n) through n=250 and the displayed larger examples.

No census of all tournament or switching types at large orders is claimed. The inherited all-order Hall theorem and the construction prove the bounds; the finite computations are controls.
