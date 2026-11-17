# An exponent ceiling for all one-apex sources in two-coordinate sevenfold readers

Let B consist of an apex and a transitive spine (L vertices), with bit w_i=1 if the apex beats spine vertex i. Its profile is h=H(B) and P=P₂(B), the latter counting ordered two-slot path covers with empty slots allowed. Write d=j−i. Then

    h = 1 + Σ_{i<j,w_i=1,w_j=0}(2^(d−1)+[d=1]),
    P = 2^(L+1) + Σ_{i<j,w_i=1,w_j=0} 2^(L−d)(3^(d−1)+1).

Indeed all directed odd cycles contain the apex, and each is specified by the inversion endpoints and an even subset of the d−1 intermediate spine vertices. There is no pair of disjoint such cycles. In the odd-cycle-pack expansion each cycle of size 3+k contributes 2 to h and 4·2^(L−2−k) to P. The sum of binomial coefficients at even k proves both formulas (d=1 gives a separate endpoint in h). Independent literal subset DP checks 255 words.

If a word admits an inversion, trim the leading 0's and trailing 1's, leaving an active interval of length m≥2 starting with 1 and ending with 0. The removed t=L−m vertices are ordinal singleton padding: h is unchanged and P is multiplied by 2^t. Equivalently this follows termwise from the formulas. The extreme inversion guarantees

    P ≥ 2^(t+1)(3^(m−2)+1).

There are at most 2^(m−2) words of active length m, and the division of t between the two ends makes no difference to the profile. If there is no inversion, B is transitive and contributes the sole profile (1,2^(L+1)) at each L.

Consequently every fixed actual two-coordinate reader q=aP+b h with a>0 and a+min(b,0)>0 has at most O(X^(log_3 2)) distinct values q≤X across **all** one-apex sources, allowing nonperiodic words, arbitrary length and arbitrary padding. Indeed q≥cP for c=a+min(b,0)>0 since h≤P. For each t the possible m satisfy 3^(m−2)≤X/(c·2^(t+1)); summing 2^(m−2) gives O((X/2^t)^(log_3 2)), and the sum over t is geometric. In particular this holds for the negative and positive sevenfold seam quotients 64P−h and 5P+h, and each fixed reader level in the opposite-sign bank. It rules out complete coverage of prime quotients in each fixed bank: the number of primes ≤X is of order X/log X, and hence a density-one fraction of primes are missed. The exponent is an upper bound, not a claim of a matching lower bound or best possible decay.

In the paired (n,7n) reader where both coordinates scale, this refutes using a single such fixed reader with a one-apex variable source to select all prime n. Arithmetic realization of 8n as a sum does not change the shared source bottleneck; 8n itself cannot be a tournament Hamiltonian-path count. The proof does not exclude variable readers whose coefficient of P tends to zero after normalization, or sources with two or more apices.

Exact inverse for a prescribed q within this family: enumerate legal (t,m) by the displayed inequality and active words of length m, read the two coefficients, and check equality. Padding may all go on one side. This searches O(q^(log_3 2)) words, not an increasing census in tournament order; the profiles above are polynomial time in word length.

Analogy: the LRC14 framework of THM-4160 (anchored Haar deletion cover) bounds the admissible repair labels by an analytic frequency cutoff before an exact candidate test. Its archived wedge-topology refinement (`continuing10_20260907_lrc_wedge_topology.md`) additionally retains edge incidence, rather than treating projected margins as freely compatible. Here the span bound is the analytic cutoff; the active first/last inversion is the required same-word incidence. Neither the modulus alone nor independently approximating h and P preserves it.

### All opposite-sign levels at once

For u_j=(8^j−1)/7, j≥1, the actual quotient readers are
`(1+4u_j)P+u_j h` and `(16+48u_j)P−u_j h`.
Both have q≥4u_j P, and u_j≥8^(j−1). The preceding constants are uniform after replacing X by X/(4u_j). Summing in j gives a convergent geometric sum, so the **union over every source word and every reader height** still has O(X^(log_3 2)) values ≤X. The transitive corner has O((log X)^2) values across all heights, negligible. Taking a subsequence of physical simultaneous prime pairs cannot evade this bound.

### Arbitrarily many ordinal one-apex sources still cannot cover the primes

Let B be an arbitrary finite ordinal sum of one-apex blocks, including transitive blocks, with no restriction on the number or sizes of the components. Under ordinal sum both h and P multiply, by partitioning each labelled path slot in block order. Normalize each individual block to its active word (plus external singleton factors), as above. The profile of B is determined by a **multiset** of active words and the total singleton padding; arbitrary order of these words cannot increase the set of profiles.

Give each active word the cost P of its normalized block. For any real s>log₃2, the span bound gives Σ_words P^(−s)<∞. All these costs are at least two, so the weighted multiset sum

    Σ_multisets (product of costs)^(−s) = ∏_words (1−P^(−s))^(−1)

is finite, since the logarithm is bounded by a constant times Σ_words P^(−s). Singleton padding adds the convergent factor (1−2^(−s))^(−1). Accordingly the number of distinct source profiles with P≤Y is O_s(Y^s): count possible multisets with padding, each contributing weight at least Y^(−s). Thus for every ε>0, even this unbounded composition closure reads only O_ε(X^(log₃2+ε)) values up to X for a fixed aP+bh with positive domination of P, or for the entire growing opposite-sign hierarchy by geometric summation of u_j^(−s). In particular it still misses a density-one fraction of primes for sufficiently small ε. No assertion that all profile products are distinct is needed.
