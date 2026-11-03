# A short finite certificate for the 7 and 21 tournament gaps

Let H(T) be the number of directed Hamiltonian paths in the finite tournament T. Then H(T) is neither 7 nor 21.

**Only external input.** Busch's theorem states that a strong tournament on n vertices has at least 5^((n−1)/3) Hamiltonian paths; see Arthur H. Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, Electronic Journal of Combinatorics 13 (2006), N3, https://doi.org/10.37236/1141. In particular, for n≥5 the bound is strictly above 7, and for n≥7 it is at least 25, strictly above 21. No global claim about other gaps follows.

**Two elementary reductions.** The strong components of any tournament are totally ordered, with every edge pointing forward between components. Consequently every Hamiltonian path visits the components contiguously in that order, and H(T) equals the product of their path counts. Singletons contribute 1. Since 7 is prime, H(T)=7 would require a strong component with H=7. Once 7 is excluded, a factorization of 21 forces a strong component with H=21: the only positive divisors of 21 are 1, 3, 7, and 21.

A strong tournament of order at least 3 has a directed Hamiltonian cycle. For completeness: take a longest directed cycle C. If a vertex v outside C has both an in-neighbor and an out-neighbor in C, some consecutive cycle edge c→d satisfies c→v→d, contradicting maximality. So each outside vertex either dominates all C (set A) or is dominated by all C (set B). Strongness forces A and B to be nonempty. It also forces an edge b→a with b in B and a in A; otherwise B has no outgoing edge to its complement. Replacing any edge c→d of C by c→b→a→d contradicts maximality. A strong tournament on at least three vertices has a cycle because an acyclic tournament has a source and a sink. Thus C is spanning.

**Finite certificate.** Fix the Hamiltonian cycle 0→1→⋯→n−1→0. Each remaining edge is freely oriented, giving 1, 4, 32, and 512 orientations for n=3,4,5,6 respectively. For every such orientation, count paths directly by testing each permutation. The complete distributions, obtained by running the accompanying Python 3 script, are:

| n | H : number of orientations |
|---|---|
| 3 | 3:1 |
| 4 | 5:4 |
| 5 | 9:10, 11:5, 13:10, 15:7 |
| 6 | 15:12, 17:12, 19:12, 23:36, 25:24, 27:12, 29:66, 31:36, 33:72, 37:126, 41:24, 43:60, 45:20 |

For an entirely hand-checkable exclusion of 7, note that the four possible orientations of the two chords of a 4-cycle are cyclic rotations of one another. In the representative with chords 0→2 and 1→3, the five paths are 0123, 1230, 1302, 2301, and 3012. The unique strong 3-vertex tournament has exactly three paths. Combined with Busch, this excludes 7 without a computer.

Every orientation has the fixed directed spanning cycle, so all are strong. Conversely every strong tournament in the relevant order appears after relabelling its Hamiltonian cycle. Thus the table excludes 7 for n≤4 and 21 for n≤6. Together with Busch's bound and the component argument, this proves the claim for all finite orders. The orientation counts in the table count labeled cycle-rooted orientations, not isomorphism classes.

The certificate is deliberately direct: the code uses no graph library, class generation, strong-connectivity routine, recurrence, or prior census. It checks at most 512×720 permutations, and independently asserts the distributions shown above.
