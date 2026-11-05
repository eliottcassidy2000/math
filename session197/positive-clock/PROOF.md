# The sharp positive-unit ternary clock with a self-converse fixed root

Session 197. **The least positive amplitude of a root-fixed self-converse transitive-cloning family H_m=A*3^m+1 is exactly 72.** Eventual identities give the same conclusion. The attaining family contains explicit prime counts, but exactly one count of the form 7p with p prime. This is a theorem for one precise cloning mechanism, not unrestricted spectrum coverage. No literature-priority claim is made.

## 1. Exact physical source and identity

Take the nine-vertex tournament with outgoing row masks

    [414,260,8,50,294,327,415,318,12],    root r=4.

It is strong, and tau(v)=8-v reverses every arc and fixes r. Its exact root run coefficients are

    (a_1,a_2,a_3)=(217,216,72),

with no larger nonzero coefficient. Here a_j counts words using every nonroot vertex once and j indistinguishable root symbols, no two consecutive. The maximal run count is three because the root O-to-I matching number is two. The inherited run expansion therefore gives, for every positive clone multiplicity m,

    H_m = a_1 + a_2(2^m-2) + a_3(3^m-3*2^m+3)
        = 72*3^m+1.                                         (1)

The actual operation replaces r by a transitive m-vertex module with r's former exterior incidences. It preserves strongness: quotient strong connectivity lifts to every clone. Extend tau by reversing the clone chain. Every T_m is therefore a physically self-converse strong tournament of order m+8.

The first counts are 217,649,1945,5833,17497,52489,157465,472393. In particular the returning scalar recurrence is H_(m+1)=3H_m-2. This is a positive-unit family, not the negative-unit family A*3^m-1.

No deletion identity at m=0 is asserted. The run expansion and (1) are stated for positive multiplicity.

## 2. Reduction of any smaller amplitude to a complete finite universe

Suppose (T,r) has an anti-involution fixing r and transitive cloning satisfies H_m=A*3^m+1 eventually, with A>0. The nonnegative finite run expansion has maximal exponential base equal to its largest nonzero run index. Distinct exponential bases are linearly independent. Hence the maximum run index is three and

    a_3=A, a_2=3A, a_1=3A+1.                                (2)

Thus the identity holds for every positive m, A is an integer, and the root cross-matching number is two. Oddness of H_1 implies A is even. Extending the anti-involution to an even clone chain gives an even-order self-converse tournament, whose H is one modulo four. Consequently 9A+1 is one modulo four, so

    A=0 modulo 4.                                          (3)

The root cannot be a singleton SCC, since then transitive cloning would leave H constant. The other SCCs contribute a fixed odd factor c dividing both 3A+1 and 9A+1. But gcd(3A+1,9A+1)=gcd(3A+1,2)=1. Therefore those other components all have count one and are singletons. Deleting them preserves the entire clone-count family.

The remaining root component is strong and invariant under the anti-involution. Because that involution fixes r and can fix no second tournament vertex, this component has odd order.

For a smaller amplitude, (3) gives A<=68 and H_1<=205. Busch's strong-order inequality 5^(n-1)<=H_1^3 and the exact comparison

    5^10=9,765,625 > 205^3=8,615,125

force n<=10. Odd order leaves only 1,3,5,7,9. This applies after removing singleton padding, not to arbitrary weak order.

## 3. Complete marked enumeration through order nine

Fix the source anti-involution. Reversal followed by it acts as an involution on the odd set of Hamilton words, so at least one word is fixed. Relabel that word in order 0,...,n-1. The anti-involution becomes v->n-1-v and its fixed root becomes the middle vertex. All consecutive arcs point forward, and every remaining pair-orbit under reflection contributes exactly one independent bit.

Thus the complete normalized marked universe has 2^floor((n-1)^2/4) presentations at odd order n. The rank histograms are:

    n=1:  1 graph;     rank0:1
    n=3:  2 graphs;    rank0:1, rank1:1
    n=5:  16 graphs;   rank0:1, rank1:2, rank2:13
    n=7:  512 graphs;  rank0:2, rank1:6, rank2:154, rank3:350
    n=9:  65536 graphs; rank0:8, rank1:32, rank2:2544,
                       rank3:10016, rank4:52936.

The verifier regenerates all 66,067 presentations. For every one of the 2,711 rank-two sources it tests the exact positive-unit conditions

    H_2=3H_1-2,    H_3=9H_1-8.                              (4)

These are sufficient, not only necessary, because the maximum run index is already known to be three. The complete survivor inventory is:

    amplitude72: 11 normalized presentations, all at order9;
    amplitude200: 19 normalized presentations, all at order9;
    no other positive-unit amplitudes through order9.

In particular no smaller A is possible in the entire finite universe. Section2 makes this an all-order exclusion, while Section1 attains A=72. This proves the sharp theorem.

## 4. Prime members and an exact 7p obstruction inside the family

Direct trial division certifies the prime values

    H_5=17497,   H_6=52489,   H_8=472393.

The receipt includes outgoing masks and the physical anti-involution for each corresponding graph of orders 13,14,16. These are examples, not an assertion of infinitely many prime members.

The family's 7p members can be classified at all orders. Since 72=2 modulo7 and 3 has multiplicative order six modulo7,

    7 divides 72*3^m+1  iff  m=1 modulo6.

For such m, put t=(m+2)/3, a positive odd integer. Then

    72*3^m+1 = (2*3^t)^3+1
              = (2*3^t+1)(4*3^(2t)-2*3^t+1).                (5)

At m=1, t=1 and this is 7*31. At every later admissible m, t>=3, and both displayed factors are strictly greater than seven. A number 7p with p prime cannot factor into two integers both greater than seven. Therefore

    the only 7p value in this entire family is H_1=217=7*31.

This sharp limitation persists despite the family's exact ternary recurrence and prime examples. No limitation on other self-converse constructions or other clone roots is inferred.

## 5. Sources, scope and verification

The run/matching theorem and eventual-to-exact reduction are inherited from the S196 cloning framework and restated with a direct proof in `../../session197-prime-atoms/MINIMUM_TERNARY.md` (SHA256 `45d183b1df1f030aff1bf755bd2b5b30ce441e8a5c13e01b99b6ef8613957c54`). The marked self-converse normalization and negative-unit barrier were established in `../../session197-prime-atoms/SYMMETRIC_TERNARY.md` (SHA256 `2df7d079d49bb6e2d8287a610b275d827510f861a5793908619cea4b521408ce`). The even-order congruence is in `../EVEN_CONVERSE.md` (SHA256 `0fc0e7ca35a1d48626ffd1b49251b64dcc0564f66594bfd0e70307287c4c4cbb`). These earlier packets do not claim the present sharp positive-unit result.

Busch's published bound is Theorem1 of Arthur H. Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, EJC13(2006), N3, DOI10.37236/1141:
https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf

Run the standalone C++17 verifier in an optimized build and in an unoptimized build with library checks, comparing stdout. It uses exact forward Hamilton DP, subset-state bipartite matching, and physical clone adjacency construction. It checks every normalized source, the complete rank histograms and positive-unit inventory, source/clone strongness, every displayed anti-involution, exact clone counts, and trial-division prime certificates. It contains explicit exceptions rather than removable assertions.

The finite enumeration certifies the bounded universe. The strong-core cutoff proves all-order minimality, and the factorization proves the entire 7p claim independently of any finite scan. The strict root-fixed anti-involution hypothesis and the exact single-root cloning mechanism are essential scope restrictions.
