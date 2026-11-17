# A density obstruction for the entire binary-contact suffix family

This closes the gap between the [triangle switch proof](README.md) and the question of whether an arbitrary binary suffix can place its center. It is **not** an obstruction for arbitrary transitive quotients and port incidence patterns from that proof.

Start with the two spine words 3 1^m 3 0 and 2 1^m 3 0, now including m=0. Write the spine as a, M, b, c, where M has m transitive vertices; the apices are x and y. Recall that the pair reverses just x–a. Then append any binary word of vertices in spine order, each an ordinary sink S or beating both x,a and losing to all other old vertices (N). The theorem in the README gives an exact +2 gap for every such word. Nevertheless, the number of *distinct* attainable values at most X, over both branches and all m and suffix words, is at most

    O(X^α),  α = log(2)/log(3) ≈ 0.63093.

In particular, infinitely many primes cannot be reached by this entire architecture: the standard prime-count lower bound is of order X/log X. This is stronger than a failure of a particular address rule or a finite target search.

## Three visits force the entropy loss

Discard all terminal sinks, which do not affect H. In a word with at least two N's let i be the position of the penultimate N and k the position of the last N, indexed from 1. Set r=k−i≥1. With i,k fixed there are at most 2^(i−1) words.

Use the following physical path trace, where repeated labels M_j and V_j denote distinct runs, with V the entire exterior binary suffix:

    M_1 V_1 a M_2 b c V_2 x M_3 y V_0 .

Check each joint: M and the nonport b,c,y beat all V; the first two V runs exit to a or x, both valid when their last vertex is N; the remaining joints are a→M, M→b, b→c, x→M and M→y. The arc x–a does not appear, so the trace works on both branches. Label the i-th N as the last member of V_1 and the k-th as the last of V_2. Those two runs may start with either type. If i≥2, reserve the first vertex for V_0; every other vertex before i has three free choices and every vertex strictly between i,k has two. Within each choice, the order in a run is forced. There are at least 3^(i−2)·2^(r−1) exterior decompositions. If i=1 and r≥2, reserve one middle vertex for V_0; this leaves at least 2^(r−2) decompositions, sufficient up to an absolute factor. The sole case (i,k)=(1,2) is the single word NN, treatable by the initial core lower bound.

When m≥3, reserve three distinct M vertices for distinct M runs and allocate the other m−3 independently, giving at least 3^(m−3) decompositions, independent of the suffix choices. When m=1,2 use instead V_1 a M b c V_2 x y V_0, with just one M run; for m=0 erase M (a→b). Both alternatives have the same exterior freedoms. In all cases including the exceptional NN word, we obtain with a uniform constant c>0

    H ≥ c·3^(m+i)·2^r .

The constants and all claims here come from a literal subcollection of Hamiltonian paths, not cancellations in a transfer polynomial.

## Summing the allowed addresses

For each m,r, the inequality H≤X allows only i with 3^i≤ X/(c 3^m 2^r). Sum 2^(i−1) over those integers i. The geometric sum is bounded by a constant times (X/(c 3^m 2^r))^α. The sums over all m≥0 and r≥1 are both convergent after factoring out X^α. Words with zero N have only the unextended outputs, O(log X) choices of m. Words with one N have shape S^rN: for m,r≥1 the README formula bounds each of m,r by O(log X), and the case r=0 retains the original core lower bound. At m=0 and r≥1 the same two-run exterior trace argument applies after omitting M: V_1 a b c x y V_0, where V_1 ends at N and receives any subset of the r sinks except a reserved sink for V_0. This gives a lower bound 2^(r−1), so the case contributes only O(log X) choices. Counting addresses is already enough for the stated upper bound on distinct outputs.

The important distinction is that the two latest N's create two return gates plus an unrestricted continuation. That three-way freedom increases path count faster than binary address choice grows. Changing the contact words alone, under this fixed placement of the cyclic triangle and the fixed incidences of M, does not resolve exact placement. One must change the quotient or the relation between a vertex and the ports, not just lengthen the tape.
