# Sparse graph capacity and symmetry-preserving Euler repair

Status: PROVED elementary bounds. Numerical thresholds below are FINITE-EXACT checks, not optimal impossibility thresholds. No external priority claim is made.

## 1. Why an enlarged codec is essential

Let a(n) count tournament isomorphism classes, and b(n) count ordinary simple graph isomorphism classes on n vertices with n-1 edges, without a connectedness requirement. For n>=2,

    a(n) >= 2^binom(n,2) / n!,
    b(n) <= B(n) := binom(binom(n,2), n-1).

The first inequality divides the labelled universe by the largest possible relabelling orbit; the second simply forgets the quotient. Consequently an injective map from all tournament classes to ordered k-tuples of such same-size graph classes requires

    2^binom(n,2) <= n! B(n)^k.

There is an even smaller direct hostile: at n=3 the only two-edge graph class is P3, while the transitive tournament and directed triangle are distinct. Thus no fixed-length tuple of these unmarked same-size graph classes can encode both. This does not apply when marks, alignment, or varying output sizes carry additional information.

The exact integer test first violates this necessary condition at n=13,28,43 for k=1,2,3 respectively; the violation persists through n=500 in the supplied computation. These are sufficient obstructions from loose bounds, not the first actual class-count mismatches. They already refute any proposed universal same-size pair or triple encoding. For every fixed k, failure for all sufficiently large n follows from quadratic versus n log n growth, independently of the finite test.

If each of k output graphs instead has N vertices and N-1 edges, then

    k (N-1) log2(e N/2) >= binom(n,2) - log2(n!).

This uses binom(M,r) <= (e M/r)^r, with M=N(N-1)/2 and r=N-1. For fixed k, an injective family therefore needs

    N >= (1/(4k) + o(1)) n^2 / log2(n).

More precisely, the lower limit of N log2(n)/n^2 is at least 1/(4k). If that ratio is unbounded the assertion is automatic; otherwise log2 N <= 2 log2 n + O(1), giving the bound directly. This is a necessary information-capacity bound, not an existence theorem at that size.

When a tuple of fixed length k has variable sizes with total at most L, its number of possible size allocations contributes only a polynomial factor in L. The same upper bound O(L log L) on the logarithm of its capacity gives total size L=Omega(n^2/log n). No arbitrary unbounded auxiliary bit string is allowed for free. Named vertices, embeddings, external labels and mark data must be included in the output specification and its capacity.

The graph codec in ../trijection/ enlarges its carrier and therefore does not contradict these bounds. A sparse graph can still encode an arbitrary hard core; it cannot encode quadratically many independent input bits in only linearly many unnamed vertices.

## 2. The global Euler identity

For a graph G with v vertices, e edges, c connected components and total cycle rank beta,

    beta = e-v+c.

Thus e=v-1 says beta=c-1. It does not say beta=0. Write t for the number of tree components and let beta_i>=1 be the cycle ranks of the other components. Then

    t = 1 + sum_i (beta_i-1).

A connected retained core Q of cycle rank beta>=1 therefore forces at least beta tree components in any disjoint completion satisfying e=v-1. Adding further cyclic components cannot reduce this requirement.

This is the common mechanism behind the October20 tree-chart obstruction and the enlarged sparse codec: a cycle obstruction can be balanced in the scalar Euler characteristic by a disconnected component, but it is not removed as a relation.

## 3. Odd symmetry forces rigid, distinct tree components

If Aut(G) has odd order, then:

1. No two connected components are isomorphic. Exchanging two identical copies, using any fixed isomorphism and its inverse, gives an involution.
2. Every tree component is asymmetric. A nontrivial automorphism group of a finite tree contains an involution. To see this, use its center. If the center is an edge and an automorphism exchanges its endpoints, it can be chosen as an involution by matching the two rooted halves. Otherwise root at the fixed center (or subdivide the fixed central edge). Rooted-tree automorphisms are built from child-subtree automorphisms and permutations of identical child subtrees. A nontrivial child permutation gives a transposition, and descending a nontrivial child automorphism reaches such a transposition in a finite tree.

Thus an automorphism-exact disjoint repair of a connected core Q with odd automorphism group must supply at least beta pairwise nonisomorphic asymmetric trees. Merely adding beta isolated vertices preserves the class of Q but introduces an S_beta kernel when beta>=2, violating the symmetry contract.

## 4. Optimal order of the padding cost

Let t distinct unlabelled tree components have sizes s_1<=...<=s_t. The number of unlabelled trees of size j is at most the number of rooted plane trees of that size, Catalan_(j-1)<=4^(j-1). Therefore

    #{tree types of size <=s} <= (4^s-1)/3,
    s_i >= log_4(3i+1),
    sum_i s_i >= sum_(i=1)^t log_4(3i+1)
                = t log_4 t - O(t).

The bound holds for all distinct trees, so in particular for asymmetric ones. Combining with the previous section, a disjoint core-preserving repair with odd automorphism group needs Omega(beta log beta) additional vertices as beta tends to infinity.

The explicit distinct asymmetric caterpillar construction in ../trijection/ attains O(beta log(beta+1)). Hence its padding order is asymptotically optimal **within this disjoint core-preserving, odd-symmetry contract**. This is not a lower bound for arbitrary graph encodings or for constructions allowed to alter the core.

## 5. Transfer and limits

- Source: a connected graph carrier with an exact odd-order symmetry group. Target: a disconnected e=v-1 carrier with the same group.
- Preserved: the original component and all its automorphisms, provided the added components are pairwise distinct rigid trees and not isomorphic to the core.
- Additional cost: distinct component identities consume vertices. The excess rank is paid by real tree components, not by a formal cancellation of scalar terms.
- Connection to tournament memory: losing type identities can create spurious permutations just as losing named attachment alignments enlarges a stabilizer. This is an exact symmetry phenomenon, not an identification of trees with tournament ports.
- No claim: a new polynomial-time algorithm for graph/tournament isomorphism, a natural size-preserving Royle bijection, or a major open-problem solution.

Reproduce the exact integer capacity controls with `python3 check_capacity.py` and `python3 -O check_capacity.py`. The accompanying trijection checker tests the actual rigid-component codec.
