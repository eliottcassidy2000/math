# Binary-coded strong tournament and a spectrum gate

Status: exact consequences of the S158 attachment lemma and S159 cycle-packing profile, not a proof of tournament path-count universality or novelty. No square-sum endpoint restriction is imposed.

## A quadratic-size injective count code

For an integer N>=1, write its standard base-two expansion b_1...b_m, starting with 1. Let N=0 correspond to the empty word. Begin with the directed triangle on 1,2,3, with sole backward edge 3->1. For i=1,...,m, if the current largest vertex is u, put k_i=2i+b_i and v=u+k_i+1, add u+1,...,v, and orient all new edges from smaller to larger labels except v->u.

The resulting tournament is strongly connected: the consecutive backward-edge intervals each induce a spanning directed cycle and consecutive intervals share a vertex. By the S158 contiguous attachment lemma, h=3 product_i (1+2^{k_i}). The initial hypothesis holds because the 2-vertex tournament after deleting 3 has no Hamiltonian path ending at 1. At every later stage, k_i>=1 and after deletion of the new top vertex no path can end at u, because the new interior vertices are dominated by all older vertices and hence no old vertex can follow an interior vertex. Hence the lemma continues to apply.

The factors are strictly increasing in their exponent. For a nonempty word, nu_2(h/3 - 1)=k_1 because k_1<k_j for every j>1; divide h/3 by 1+2^{k_1} and repeat. At step i, recover b_i=k_i-2i. The empty product represents N=0. Thus h itself injectively encodes N; it does not generally equal N, and all output counts are divisible by 3. The exact number of vertices is 3+m^2+2m+sum_i b_i = O(log(N+2)^2).

`python3 check.py` uses independent subset-and-endpoint dynamic programming on the literal tournaments to verify words of length at most two, including the leading-zero words 0, 00 and 01 as additional controls. The proof is the attachment and valuation argument, not finite checking.

## All integers as a triangle statistic, not as h

For N>=1 let X_{N+2} be the tournament on 1,...,N+2 in which all edges point upward except (N+2)->1. Every directed triangle uses that edge, and the third vertex can be any of 2,...,N+1, so X has exactly N directed triangles. The cycle 1->2->...->N+2->1 proves strong connectivity. For N=0 use a singleton. Therefore 2 times the triangle count realizes every nonnegative even integer, and 1 plus twice the triangle count realizes every positive odd integer. These are graph statistics, not Hamiltonian path counts. In fact S159's odd-cycle packing formula gives h(X_{N+2})=1+2^N: every nonempty odd-cycle packing has one odd cycle through both endpoints and one odd-sized choice of the N interior vertices. For N=3 the statistic 1+2N is 7 but h is 9; so the statistic cannot settle the 7 or 21 spectrum problem.

## A family-specific 7 and 21 exclusion

For nonempty X define G(X)=triangle(X,1,1). S159 proves h(G(X))=h(X)+F_X(2), while its positive cycle-packing formula gives F_X(2)>=2^{|X|}. If h(G)=7, |X|<=2, and the only cases give 3 or 5. If h(G)=21, |X|<=4. Inputs of orders 1,2,3 give respectively 3, 5, and 9 or 15. For order 4, the transitive tournament gives 17, while any nontransitive tournament has at least one triangle, so h(X)>=3 and F_X(2)>=16+8=24, giving h(G)>=27. Hence neither 7 nor 21 occurs as a one-wrapper output G(X), without implying that neither occurs in all tournaments.

## Self-converse versions

For the same bits and k_i, use a chain of overlapping backward intervals with interior-size list (1,k_1,...,k_m,k_m,...,k_1,1). The first interval is the initial directed triangle. The attachment lemma yields h=9 P^2 where P=product_i(1+2^{k_i}). Its order is 5+2m²+4m+2 sum_i b_i. The tournament is strong, since the interval cycles overlap consecutively. The interval list is palindromic, so reflection of the total vertex order fixes the set of backward edges and is an anti-automorphism. Thus the tournament is self-converse. The exact positive square root P=sqrt(h/9), followed by the previous 2-adic factor decoding, recovers N injectively from h. This is another count code, not a realization of h=N.

The family X_{N+2} in the triangle-statistic section is also self-converse: reflection i |-> N+3-i preserves its one backward edge. The singleton for N=0 is self-converse too. Therefore every N>=0 is the exact triangle count of a strong self-converse tournament; 1+2N is a derived statistic, not its number of Hamilton paths.

## Adaptive binary trajectory for exact triangle count

Start with a directed triangle, state x=1. For N>=1 read the remaining bits in its ordinary binary expansion. On the first remaining bit b, insert 1+b interior vertices into the single reversed-edge interval: that produces X_{4+b}, containing the original triangle as an induced subtournament, and state x=2+b. It also has the zero-endpoint property needed by the attachment lemma. If there are no remaining bits, stop at the original triangle.

For each later bit b, x>=2. Let u be the current maximum label, put k=x+b-1>=1, add k interior vertices and an extreme vertex v=u+k+1, and orient all new edges upward except v->u. The backward-edge intervals overlap in one endpoint, so the result stays strong. There are exactly k+1 new directed triangles: k consist of u,v and an interior vertex, and one is {p,u,v}, where u->p is u's sole outgoing backward edge. Old triangles persist and no other new triangles arise. Thus x'=x+k+1=2x+b. At termination x=N, and the vertex count is N+2. For N=0 use a singleton. The inserted or appended vertices preserve the previous tournament as an induced subgraph, after relabeling the original endpoint on the exceptional first step.

The path count has a different exact update. For the initial exceptional step it is 1+2^{2+b}; later h'=h(1+2^{x+b-1}) by the attachment lemma. The dynamic 2x+b and terminal value N concern the triangle statistic, not the Hamilton-path count. `python3 adaptive_check.py` checks literal triangle counts, strong connectivity and order for N=0,...,30 and separately checks the self-converse one-backedge family.

Dependencies: [S158 contiguous attachment lemma](https://github.com/eliottcassidy2000/math/blob/6bb391f902df4a6fd34b947a7923f64a0d5e6057/05-knowledge/results/session158_creative_20261103/square_overlap/CONTIGUOUS.md) and [S159 odd-cycle path-cover profile](https://github.com/eliottcassidy2000/math/blob/7589463e54e20585aff6a2e68019ee4a4af582cd/05-knowledge/results/session159_constructions_20261103/cycle_profiles/NOTE.md).
