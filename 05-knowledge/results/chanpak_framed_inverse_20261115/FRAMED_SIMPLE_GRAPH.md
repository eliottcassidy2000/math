# What the smaller simple graph can and cannot select

Use the ascending convention for the fixed path 1→2→…→n. At every pair i<j<n, put an undirected edge ij in B exactly when j+1→i in the tournament. This is the staircase chart: every subset of edges of K_(n−1) is permitted; the orders of its vertices are indispensable. Fixing labels, the chart is the unique directed-path gauge within the switching class (THM-474), and forgetting the path yields the usual |Hamilton paths|/|Aut| fibre (THM-1430). An edge toggle in B is one free arc reversal, so cube geometry precedes the isomorphism quotient. Reflection i→n−i of B induces converse together with reversal of the frame; arbitrary permutations of the graph generally do not.

For any ordered vertex sequence P with distinct elements of [n], its admission indicator is a product of constraints on B:
- ascending consecutive pairs need no constraint;
- descending consecutive pairs reject;
- ascending nonconsecutive a,b need the nonedge {a,b−1};
- descending nonconsecutive b,a need the edge {a,b−1}.

Sum these products for Hamilton paths, or multiply the admission indicator for an even-order root-started sequence L by the indicator for a Hamilton path on its complement and sum to get E. Include one empty complementary path when appropriate. Thus H+2E is a partition sum over one- and two-path systems of the *ordered* tile graph, with both positive and negative literals. This is exact before projection to the merged metagraph.

## Same smaller graph, different strong physical outputs at every order

Let W_(n,b), 3≤b≤n, have only the backward arc b→1, and let u=1. The smaller graph has one edge {1,b−1}, so its unlabelled rooted isomorphism type (root 1) cannot see b. Its undirected Hamilton-path count is zero for n≥4. Its larger tournament has

    H(W_(n,b)) = 1+2^(b−2).
    E_u(W_(n,b)) = 2^(n−2) − (2^(n−b−1) if b<n, else 1).

Indeed, all odd directed cycles contain b→1, so no two coexist in the odd-cycle packing expansion. Each such cycle corresponds to an odd subset of {2,…,b−1}. A rooted path cannot use the only backward arc because it starts at 1. Its complement is transitive. Of the 2^(n−2) even vertex subsets containing 1, exactly the subsets whose second vertex in increasing order is b are forbidden. They are {1,b} plus an even subset from b+1,…,n. The first SCC is [b], so adding v with only out-neighbour 1 gives a genuinely strong tournament and exactly H+2E paths, not merely a formal count. For n=5 the outputs are 15,19,23. This family diagnoses a real obstruction: neither undirected Hamilton counts nor the full unlabelled rooted tile graph suffice to choose a prescribed target. A valid reduction must keep placement in the staircase or else actively choose it.

## A separate graph-to-tournament lens: fence weights

This is *not* the staircase chart. Let G be any simple graph on m labelled vertices. For each linear order ρ give every nonedge the forward direction in ρ and every edge the backwards direction; write T(G,ρ). This is also the remaining tournament of the source-star gauge after adding a dominating vertex, because any orientations are legal at pairs not incident with the source.

For a permutation P of the vertices, let β be the word of G-adjacencies along its successive pairs. Define e(β) to be the number of linear extensions of the fence poset with the respective forward and backward inequalities along P. In particular these are the numbers of permutations with a prescribed exact descent set. Then

    sum_ρ H(T(G,ρ)) = sum_β A_G(β) e(β),

where A_G counts the ordered vertex sequences of adjacency pattern β. Proof: exchange the sum over orders and over paths; each path requires precisely those adjacent comparisons, which never contain a directed cycle because their underlying graph is a path. For m≥2, the monochromatic A_G words are twice the ordinary undirected Hamilton-path counts of G and of its complement, with fence weight 1. The mixed words are unavoidable; e.g. for path words of length two the constant words have weight 1, while the alternating words each have weight 2.

The same exchange of order sums gives an exact averaged rooted target. Sum over even-order sequences L starting at u, and spanning sequences R of the complementary vertices. Each pair contributes

    binom(m,|L|) e(β_L) e(β_R)

to sum_ρ E_u(T(G,ρ)), since the respective fence orders are disjoint and can be interleaved arbitrarily. Set e=1 for an empty sequence; one-vertex sequences have e=1. The marked root stays at the start of L when choosing its sequence, but its position in ρ is unrestricted.

This explains exactly where the ordinary undirected Hamilton-path statistic sits inside a larger statistic with all adjacency words. It offers a moment identity for choosing an ordering; it is not a support theorem at any individual absolute target. In continuum language, it is a projection with loss of information: one must invert an ordered, boundary-marked observable, not just a scalar average. The repeated counterexample above pinpoints the information to keep.

## A sufficient condition under which a seam sign really is monotone

For one free arc a→b, paths using it are in bijection with Hamilton paths in the directed graph contracting a→b to one vertex q: set incoming arcs to q by incoming arcs to a and outgoing arcs from q by outgoing arcs from b; keep other arcs unchanged. This generally gives a digraph with empty and doubled pairs, not a tournament. Compare the analogous contraction after reversing the arc. The difference in H is precisely the difference of the two contracted counts; paths avoiding the arc cancel.

In particular if N^-(a)⊆N^-(b) and N^+(b)⊆N^+(a) on the other vertices, every original path using a→b survives as a path using b→a by interchanging the two consecutive vertices, so the reversal cannot lower H. The same assertion holds for the physical one-outlet objective, by applying the comparisons in the augmented tournament: they automatically agree at the new outlet if a,b are both different from u. This criterion detects when a seam edit is allowed to make a guaranteed sign, but does not imply a bound on the jump or coverage.

## The requested rooted extension is a single leaf in the tile graph

If u<n, the fixed Hamilton path of W remains fixed after appending v=n+1: its old last vertex n beats v. The sole backward entry in the new last column is v→u, so the new smaller graph is exactly B with a new *last-labelled leaf* n adjacent to u. A dominating sink, which leaves H unchanged, instead adds a last-labelled isolated vertex. Thus the exact target difference 2E_u is the increment from toggling that sole edge. For u=n the proposed append order is illegal; another frame is necessary. This distinction arises even though both operations have trivial unlabelled graph descriptions.

For G=B together with the last-labelled leaf, the ordinary undirected Hamilton path through that leaf, if it exists, must start or end there; it counts only paths of B starting at u. In contrast the directed Hamilton objective reads even-order rooted prefixes with arbitrary complement, via the exact weighted formula. This explains why two superficially parallel leaf operations calculate genuinely different statistics.

For n>b the family can be written 1+2^(b−2)+(2^(b−1)−1)2^(n−b). In particular b=5 realizes 9+15·2^t for every t≥1, including the recurring count 69 at t=2. Its concrete order-eight witness has all arcs ascending except 5→1 and 8→1; the final tile graph is just a two-leaf star and isolates. The old tournament is not strong (it has trailing singleton SCCs), but the added vertex gives a strong witness. None of this asserts minimum order or novelty for 69.
