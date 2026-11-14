# An exact arbitrary-Q source compiler and a least-resistance test

Write H(T) for Hamiltonian paths and Q(T)=Q_T(2) for ordered two-path covers with empty slots allowed. The [S325 digit bank](../session325_absolute_cut_overlap_20261114/FULL_99_DIGIT.md) requires the pair, not just H or a residue of Q. This note gives an exact constructive inverse for one-apex chain cores and their ordinal products, with all Q admitted. No claim of completeness outside that source language is made.

## 1. A cut decoder that uses every actual carry digit

A strong one-apex source is a forward transitive chain of length n with an apex whose outcut word is x_1…x_n. Strongness is exactly x_1=1,x_n=0: otherwise there is a source or sink, while these two arcs supply a spanning cycle. Let M_d count inversions 1…0 at distance d. The [S323 cycle-packing formulas](../session323_computability_and_native_descent_20261113/PAIRED_END_DESCENT.md) are

    H = 1 + Σ_(d=1)^(n−1) M_d(2^(d−1)+[d=1]),
    Q = 2^(n+1) + Σ_(d=1)^(n−1) M_d 2^(n−d)(3^(d−1)+1).

For H fixed, at most two n need consideration:

    ceil(log₂(H+1)) ≤ n ≤ 2+floor(log₂(H−1)).

Indeed the extreme pair supplies 2^(n−2), and every Hamiltonian path is determined by the subset of chain vertices before the apex, with the two sides traversed in chain order. Thus H≤2^n−1 by oddness. The Q formula also gives Q≥2·3^(n−2).

Read the word from both ends, a_i=x_(i+1), b_i=1−x_(n−i), with a_0=b_0=1. The k-th outer-span count is C_k=Σ a_i b_(k−i)=M_(n−1−k). Let R be the prescribed Q minus 2^(n+1); its exact contribution at layer k is

    c_k C_k,  c_k=2^(k+1)(3^(n−2−k)+1).

The valuations of these c_k are nondecreasing. After assigning the first k layers, compare the partial sum against R modulo 2 to the valuation of the next coefficient. Reject whenever they differ. Later layers cannot change those tested bits. Once the two ends meet, try either central bit if present and recompute the exact Q and H. This procedure accepts exactly the physical words at the prescribed pair. It generalizes the dyadic compiler of S323 (which insisted on a string of zero carry digits) to any specified digits of Q.

## 2. Uniform branching bound

For even n there is an initial singleton valuation layer k=0, and then the coefficient valuations tie in pairs (1,2),(3,4),…; for odd n they tie in pairs (0,1),(2,3),…. Successive paired valuations differ by two. For a complete interior pair (t,t+1), t≥2, factor out the common power of two. The coefficients are odd A=(3^e+1)/4 and B=(3^(e−1)+1)/2 for odd e=n−2−t. Once previous bits are fixed, the new four bits a_t,b_t,a_(t+1),b_(t+1) have coefficients modulo four

    A+B b_1, A+B a_1, B, B.

Thus the prescribed next base-four digit selects:
* at most six assignments when a_1=b_1=0 (four odd coefficients);
* at most six when a_1≠b_1 (three odd and one even);
* at most eight when a_1=b_1=1 and A+B=0 mod4, and exactly four in that case when A+B=2 mod4.

For odd coefficients complement bits to replace −1 by +1 at the cost of a fixed shift. For four plus-one bits the maximum multiplicity in a residue is six; for three it is three, doubled at most by the even bit. The last case uses two equal even coefficients and two equal odd coefficients, giving eight or four as stated. The expression A+B=(5·3^(e−1)+3)/4 alternates 0,2 mod4 on consecutive pairs because e−1 changes by two. The product of each consecutive two branching caps in the last case is at most 32<6². Boundary layers (including cases where the next unexposed span meets the center) add only a constant factor. Thus there are O(6^(n/4)) candidates or partial branches, with polynomial-in-n arithmetic per step. The worst-case search is O(Q^(log₃(6)/4+o(1))) arithmetic operations, exponent approximately 0.408, or O(H^(log₂(6)/4+o(1))) when the H depth bound is better. This is a guaranteed upper bound, not a lower bound or universal coverage proof.

## 3. Optimal ordinal source, not locally cheapest factor

Allow singleton SCCs and any number of strong one-apex chain SCCs. Both H and Q multiply on ordinal joins. For a candidate SCC of pair (d,e), retain one minimum-order witness found by the decoder; a singleton is the letter (1,2) of cost one. Then minimize the total order by the exact recurrence

    D(1,1)=0,
    D(h,q)=min [cost(d,e)+D(h/d,q/e)],

over available letters (d,e) dividing the respective coordinates. Every step strictly decreases q, and the divisibility constraints are exact on both axes. Recover a source by following a minimizing path and joining the literal SCCs. Factorization or a bounded physical catalogue can be chosen according to the input sizes; the carry bound above is only the per-letter decoder, not a promised complexity for factoring arbitrary h,q. Larger strong SCC types can be added to the alphabet without changing the recurrence.

There is a concrete failure of locally cheapest subtraction. The pair (9,336) occurs: use the near-transitive order-five strong core (H,Q)=(9,84) and two singleton pads. A directed triangle is a cheaper nontrivial step and its pair (3,12) divides both coordinates. But the residual is (3,28), impossible for **any** tournament. A count-three witness consists of exactly one directed triangle plus singletons, by the strong lower bound and the SCC product, hence its Q is 12·2^s. This counterexample passes the elementary requirement that residual Q be divisible by four. The least-resistance rule is valid only with the remaining exact two-coordinate feasibility as its state.
