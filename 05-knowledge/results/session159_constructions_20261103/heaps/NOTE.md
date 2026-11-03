# Strong square-sum Zeckendorf code

Status: proved construction; no claim of novelty or spectrum coverage. The count is that of directed Hamilton paths. The attachment lemma used below is established in Session 158, `square_overlap/CONTIGUOUS.md`.

Let F_2=1,F_3=2,F_(i+2)=F_(i+1)+F_i. Express N≥1 in its Zeckendorf expansion N=∑_(i=2)^m b_i F_i, with b_m=1 and no adjacent 1s. Start with the directed triangle on [1,3] whose only backward edge is 3→1. Put u=3 and k_prev=0. For i=2,...,m in order do the following:

1. Let s be the least positive integer with s²>2u+k_prev+1.
2. Put q=s+b_i, v=q²−u and k=v−u−1. Add labels u+1,...,v, orient all new edges forward except v→u, then replace (u,k_prev) by (v,k).

For N=0 use the initial triangle without additions. Every backward edge has square-sum endpoints: 3+1=4 initially, and u+v=q² for each addition. Because q²>2u+k_prev+1, k>k_prev, and the k interior labels are nonempty. Each addition glues the directed cycle u→u+1→...→v→u onto the preceding strong tournament at u; therefore all witnesses are strong.

By the attachment lemma the Hamilton-path count of the resulting tournament is H_N=3∏(1+2^k), with the product over all additions. Its hypothesis holds initially: on {1,2}, the only path is 1→2, so no path ends at p=1. At each subsequent step the newly added top vertex has exactly one backward outgoing edge, and deleting it leaves at least one new top-range interior vertex that cannot be followed by any old vertex. Thus the lemma's endpoint obstruction propagates exactly as in S158.

The mapping N→H_N is injective. Divide by 3; the strict increase of exponents implies that the least one is ν₂(H_N/3−1). Divide by its factor and iterate; the entire exponent sequence is recovered. From u and k_prev, the two possible values for k arise from distinct choices q=s and s+1, so each recovered k determines the digit and next u. The all-empty sequence is N=0. Hence the Zeckendorf digits and N are uniquely recovered from H_N.

Size: if m is the highest digit index, the number of vertices is O(m⁴), hence O((log(N+2))⁴). Indeed k_prev<u and q≤sqrt(2u+k_prev+2)+2, whence k≤k_prev+5+4sqrt(3u+2), while u_new=u+k+1. Induction gives k=O(j³), u=O(j⁴) after j steps, choosing sufficiently large absolute constants. This is a polynomial upper bound, not claimed sharp. Example N=1 gives b_2=1, k=9, 13 vertices and H=3(1+512)=1539. N=2 gives digits 0,1, exponents 2,12, 19 vertices and H=3·5·4097=61455.

`python3 check.py` verifies the empty word, words 0, 1, and 00 by independent endpoint-subset dynamic programming, and checks square sums and 2-adic exponent recovery over all binary words of length at most 8. Finite checks are controls, not proofs. Non-Zeckendorf words in the test are also allowed by the construction; their integer interpretation is not asserted to be unique under Fibonacci summation.

For one explicit bound, after j additions put U_j=10^6(j+1)^4 and K_j=10^4(j+1)^3. Inductively u≤U_j and k_prev≤K_j. As k_prev<u, the displayed bound gives k_new≤K_j+5+4sqrt(3U_j+2)≤10^4(j+2)^3. Then u_new≤U_j+10^4(j+2)^3+1≤10^6(j+2)^4. Both inequalities hold for j≥0. Therefore for N≥1, the number of vertices is at most 10^6 m^4, since the number of additions is m−1. This intentionally loose bound makes the asymptotic effective.

## Continued fraction, independence polynomial, and Riordan transform

This section packages a different family from S158, `square_overlap/ZECKENDORF.md` and `CONTIGUOUS.md`. It is a generating function over tournament witnesses, **not** the Hamilton-path count of a single tournament. For a fixed m and indices i=1,...,m put q_i=1+2^(i+3). Associate index i with the triangular-edge factor of index n=i+3, and for every independent set S of the path on 1,...,m take S158's explicit tournament T_S: its count is 265∏_(i∈S)q_i. Define

P_m(t)=∑_(S independent) t^|S| H(T_S)/265.

Then P_0=1, P_1=1+tq_1 and, by splitting according as m is absent or present,

P_m=P_(m−1)+t q_m P_(m−2) (m≥2),
P_m/P_(m−1)=1+t q_m/(1+t q_(m−1)/(1+.../(1+tq_1))).

The fractions are formal rational functions, legitimate because each P_j has constant term 1. This is the weighted hard-particle independence (clique) polynomial. For the Cartier–Foata heap monoid whose pieces are path vertices and commute exactly when nonadjacent, the standard Möbius polynomial is μ_m(t)=P_m(−t). The tournament count is encoded by the clique weights; the reciprocal heap series is a related formal object, not itself asserted to count tournament paths.

For an equal-weight illustration, choose m pairwise disjoint intervals, each with k interior labels and one reversed extreme edge. Retaining the reversal at precisely S yields H(T_S)=q^|S|, q=1+2^k. For nonadjacent selections S,

∑_(S independent) t^|S| H(T_S) = ∑_r binom(m−r+1,r)(qt)^r,
∑_(m≥0) z^m ∑_(S independent) t^|S| H(T_S) = (1+qt z)/(1−z−qt z²).

To justify the binomial coefficient, replace the selected positions s_1<...<s_r with s_j−(j−1), an r-subset of {1,...,m−r+1}. Equivalently, R_(n,r)=binom(n−r,r) is the ordinary Riordan array (1/(1−z), z²/(1−z)); the coefficient above is R_(m+1,r). These equal-weight tournaments are generally not strong because disjoint intervals concatenate. The preceding Zeckendorf construction is strong, but its adaptive weights depend on earlier digits and are not the fixed q_i in this continued fraction.

Primary dependencies: [S158 index](https://github.com/eliottcassidy2000/math/blob/6bb391f902df4a6fd34b947a7923f64a0d5e6057/00-navigation/S158.md), [attachment lemma](https://github.com/eliottcassidy2000/math/blob/6bb391f902df4a6fd34b947a7923f64a0d5e6057/05-knowledge/results/session158_creative_20261103/square_overlap/CONTIGUOUS.md), and [square-sum Zeckendorf family](https://github.com/eliottcassidy2000/math/blob/6bb391f902df4a6fd34b947a7923f64a0d5e6057/05-knowledge/results/session158_creative_20261103/square_overlap/ZECKENDORF.md).

### Audit targets

Check (1) the recursion’s q≥s and k>k_prev with no off-by-one error; (2) both hypotheses of the S158 attachment lemma initially and after every step; (3) 2-adic decoding even if factors share odd prime divisors; (4) validity of the effective polynomial bound; (5) the indexing q_i=1+2^(i+3) and the aggregate versus individual-count distinction; and (6) the Riordan row shift. `check.py` uses explicit checks, including under `python3 -O`, for all listed cases.
