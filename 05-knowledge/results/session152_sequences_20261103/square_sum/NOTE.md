# Square-sum comparison of the two sets

Let A={15,17,19,23,25,27,29,33,37,41} and B={2,3,5,7,11,13,17,19,23,27}. They have intersection {17,19,23,27}. In THM-4114 §7, A is exactly the proper nonempty ear-response H-spectrum of a specific order-five strong tournament, not a definition of square-sum admissibility. No asserted map between the tournament statistic H and square-sum order follows from equal integers.

For Q_n, vertices are 1,...,n and distinct x,y are adjacent iff x+y is a square. The repo's finite-exact enumeration (procgen_smallgraph_20260926_small_graphs_arithmetic.md, §3.1 table) gives Hamiltonian-path orders n<=300 as {1,15,16,17,23} union [25,300]. Consequently A has 9 admissible members (all but 19), B has 3 (17,23,27). The broad range [25,300] makes the overlap comparatively unselective; it is not a theorem that A equals the admissible orders. The first admissible nontrivial order is 15, with unique path 8,1,15,10,6,3,13,12,4,5,11,14,2,7,9. The path uses targets 9,16,25, i.e. the (3,4,5) Pythagorean triple (THM-4505 and THM-4510). The number 15 is also the first member of A, but without a structure-preserving map this is a numerical coincidence, not a shared mechanism.

There is a direct test using the *same* square-sum graph relation on the two sets. Q[A] has precisely edges {17,19},{23,41},{27,37}, hence is 3K2 + 4K1. Q[B] has precisely {2,7},{2,23},{3,13},{5,11},{13,23},{17,19}, hence is P5 + 2K2 + K1, with P5=7-2-23-13-3. Different edge counts and degree profiles forbid any isomorphism. This comparison concerns the induced graphs, not Q_10 or the graphs Q_n whose order is a set member.

## Birth recurrence and generating function

Write E(n)=|E(Q_n)| and b(n)=E(n)-E(n-1), n>=1, E(0)=0. The new edges are exactly {n,k^2-n} for n<k^2<2n. Thus

b(n)=floor(sqrt(2n-1))-floor(sqrt(n)),   E(n)=sum_{m=1}^n b(m).

Equivalently, each k>=2 contributes a birth at precisely floor(k^2/2)+1 <= n <= k^2-1. Therefore, as formal power series,

sum_{n>=1} b(n)z^n = (sum_{k>=2}(z^(floor(k^2/2)+1)-z^(k^2)))/(1-z),
sum_{n>=0} E(n)z^n = (sum_{k>=2}(z^(floor(k^2/2)+1)-z^(k^2)))/(1-z)^2.

These are direct consequences of edge births, not claimed novel over Session 148. For the union A∪B, ordered n=2,3,5,7,11,13,15,17,19,23,25,27,29,33,37,41, the respective births are 0,1,1,1,1,2,2,1,2,2,2,2,2,3,2,3. At n=15, E=15, one more than a spanning path needs; the unused edge in the unique path is {1,3}. Births alone do not determine Hamiltonicity: b(17)=1 but Q17 has a path, while b(7)=1 and Q7 does not.

## Source caveats

THM-4510 gives a structural three-square construction from primitive Pythagorean triples; 15 is the first square-sum chain, not the first Pythagorean number in any general sense. It specifically notes the square-sum endpoint 8 derives from 16/2, whereas the Catalan identity 9-8=3²-2³ is a different mechanism. Session 148's shell bridge similarly proves a cardinality match is not by itself an isomorphism: its explicitly defined square-sum source and prime-side graphs have 12 and 14 edges, respectively. This is a useful methodological parallel, not a structural identification with A or B.
