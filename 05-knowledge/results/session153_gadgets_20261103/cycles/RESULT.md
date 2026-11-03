# Sevenfold symmetry forces at least 175 paths

Let T be a finite tournament and H(T) its number of directed Hamiltonian paths. If Aut(T) has an element of order 7, then H(T) >= 175. The bound is attained, for example, by the 7-vertex circulant with connection set {1,2,3} modulo 7. Hence if H(T) < 175, 7 does not divide |Aut(T)|. In particular, independently of the gap theorem, a hypothetical H=7 tournament has trivial automorphism group; a hypothetical H=21 tournament has automorphism group of order 1 or 3. These are necessary conditions, not a proof of either gap.

Proof. Any order-7 permutation has a 7-orbit U. Its induced tournament is a Cayley tournament on Z/7 with connection set S containing one element of each pair {±1}, {±2}, {±3}. The eight possibilities have exactly 175 or 189 Hamiltonian paths, by the independent direct permutation and odd-cycle counts in check.py. Every directed Hamiltonian path in T[U] extends to one in T: insert each remaining vertex x in a prescribed order at the first position before a path vertex which x beats, or append it if none exists. The predecessor, if any, beats x by the first-position choice. This insertion is injective on the input paths, since restricting the output order to U recovers the input. Thus H(T) >= H(T[U]) >= 175. By Cauchy's theorem if 7 divides the automorphism-group order, an order-7 element exists. Finally, the automorphism group acts freely on Hamiltonian paths because a permutation fixing an ordered path fixes every vertex, so |Aut(T)| divides H(T). Combining the two facts gives the stated small-H consequences. The free-action divisibility already occurs in repository theorem THM-1745; novelty of the 175 threshold has not been established by full literature search.

The reproducible cross-check independently enumerates all 7! path orderings and oriented simple cycles (each cycle only once up to rotation), combines the latter by the odd-cycle formula, and verifies equality. Its output columns are orientation bits, connection set, counts c3,c5,c7 of directed cycles, number of vertex-disjoint triangle pairs, and H:

000 [4,5,6] 14 28 17 14 175
001 [3,5,6] 14 42 24  7 189
010 [2,4,6] 14 28 17 14 175
011 [2,3,6] 14 28 17 14 175
100 [1,4,5] 14 28 17 14 175
101 [1,3,5] 14 28 17 14 175
110 [1,2,4] 14 42 24  7 189
111 [1,2,3] 14 28 17 14 175

For n=7, H=1+2(c3+c5+c7)+4p, with p the number of pairs of vertex-disjoint directed triangles. No other nonempty odd cycle packing is possible.

Primary source: Grinberg and Stanley, Theorem 1.39, https://arxiv.org/abs/2307.05569 . It states U_T = sum_sigma 2^{psi(sigma)} p_type(sigma), summed over permutations whose cycles are odd and whose nonsingleton cycles are directed cycles of T; psi counts nonsingleton cycles. The Hamiltonian-path specialization maps every p_lambda to 1, yielding precisely the weighted cycle-packing expression. This is equivalently evaluation at x_1=1 and x_i=0 for i>1.
