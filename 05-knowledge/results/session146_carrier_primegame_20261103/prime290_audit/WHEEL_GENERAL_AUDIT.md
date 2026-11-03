# Independent audit: complete +2 wheel path distribution

**Verdict: correct**, for distinct odd primes p<q. The smaller prime is essential to the specified parameterization; in unordered notation replace p and q by min(p,q) and max(p,q). This is a local coprimality theorem, not a primality or twin-prime theorem.

Modulo pq, write the +2 cycle as r=pq+2t, so deleting the nonunits is deleting multiples of p or q. The q deleted multiples of p divide the cycle into q gaps of p−1 surviving sites. The remaining p−1 q-multiples are cyclically separated by q>p, hence belong to distinct p-gaps, and their positions within those gaps are kq mod p for k=1,…,p−1, a permutation of 1,…,p−1. Deleting a position j yields lengths j−1 and p−1−j, omitting zeros. For each ℓ=1,…,p−2, exactly two deleted positions yield a component of length ℓ; the q−p+1 unaffected gaps yield components of length p−1. The total path count is 2(p−2)+q−p+1=p+q−3, and the vertex count is (p−1)(q−1). Therefore the edge count is (p−1)(q−1)−(p+q−3)=(p−2)(q−2). The histogram makes the upper bound sharp: maximum path length p−1.

The hole count p+q−1 and exactly two hole-to-hole edges follow separately by inclusion-exclusion and the two distinct CRT systems described in the audited proof. On the underlying pq-cycle, pq≥15 prevents a directed edge and its reversal from being the same physical edge. The actual wrap edge 2pq−1→1 exists exactly once.

Independent enumeration traversed the original odd-residue graph, rather than the producer's t-labeled run algorithm, and verified the entire histogram and all the other claims for 276 unordered distinct odd-prime pairs below 100. Examples: (3,5) gives {1:2,2:3}; (3,7) gives {1:2,2:5}; (5,29) gives {1:2,2:2,3:2,4:25}. `wheel_general_independent.py` uses explicit checks active under optimization; normal and `-O` outputs match. Producer's revised `wheel_general.py` likewise passes in both modes. These finite tests supplement the proof, not replace it.

