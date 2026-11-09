# Hamiltonian path counts in semicomplete digraphs

Loops do not affect Hamiltonian paths, whose vertices are distinct. Removing the diagonal leaves the same count. The negation of an absent unordered pair is semicompleteness, not the tournament condition: both directions may occur.

## A self-contained generalized Camion lemma

Every strong semicomplete digraph on at least two vertices has a spanning directed cycle. Here a two-cycle is allowed at order two.

Proof. Strongness supplies a directed cycle. Choose a longest one C. If an outside vertex v can be inserted into a cycle edge, the cycle was not longest. Otherwise, if c_i -> v, noninsertability forbids v -> c_(i+1), so semicompleteness forces c_(i+1) -> v. Propagating around the cycle shows C -> v and that no reverse arcs occur. If no c_i -> v, semicompleteness gives v -> C, again with no reverse arcs. Thus every outside vertex lies in A (strictly dominates C) or B (strictly dominated by C). Strongness, if any outside vertices exist, forces both classes nonempty and a directed path from B to A. This path cannot first enter C from B, since there are no B -> C arcs. Its first arrival in A therefore supplies an arc b -> a with b in B and a in A. Replace a cycle edge c_i -> c_(i+1) with c_i -> b -> a -> c_(i+1), a contradiction. Thus C is spanning.

For order n >= 3, a Hamiltonian cycle uses different unordered pairs for all of its arcs. Choose its arcs and one direction of every other pair; this is a spanning strong tournament T inside D. Consequently H(D) >= H(T) >= f(n), the minimum strong-tournament count. The order-two exception has H=2.

## Finite-core theorem for counts 7 and 21

The SCCs of a semicomplete digraph are totally ordered, with every cross-component arc pointing forward and no reverse arcs. A Hamiltonian path traverses the components consecutively, so H(D) is the product of their Hamiltonian-path counts.

The tournament theorem in `01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md` supplies f(3)=3, f(4)=5, f(5)=9, f(6)=15, f(7)=25, and nondecreasing f(n) by Moon's strong-subtournament theorem plus path insertion.

It follows unconditionally from these ingredients that:
- If H(D)=7, there is exactly one nontrivial SCC, it has order four and count seven, and every other SCC is a singleton. Thus an exhaustive order-four classification is an all-order classification after arbitrary singleton padding before and after the core.
- If H(D)=21, either there is one count-21 strong core, necessarily of order at most six, or there are exactly two nontrivial cores with counts three and seven, in either order, and all remaining SCCs are singletons. The count-three core has order three and is either the directed triangle or the directed triangle with exactly one edge doubled; both have exactly three Hamiltonian paths. An exhaustive scan through order four excludes a count-21 strong core below order five, so only orders five and six remain.
- The finite classification of count-21 cores is complete only after the strong order-six census is finished and verified. No statement that all other counts occur follows.

## Verified seed and infinite families

On vertices 0,1,2,3 take adjacency rows
0111 / 0011 / 1001 / 1000.
The seven Hamiltonian paths are exactly
0123, 1203, 1230, 1302, 2013, 2301, 3012.
Its doubled pairs are 02 and 03. It is strong, has zero diagonal and no absent pair, and refutes the count-seven conjecture.

Prepending/appending any transitive chains gives count-seven counterexamples at every order n >= 4. Taking the ordered sum with a directed triangle gives count 21, and padding gives counterexamples at every order n >= 7. Smaller count-21 seeds may improve this threshold; the parent's scanner has reported order-five examples, which this note does not independently enumerate.

## Exact generators and path-count identities

For any nonempty subset S, let E_D(S,v) count Hamiltonian paths of D[S] ending at v, and B_D(S,v) those starting at v. Initialize E({v},v)=B({v},v)=1. The exact subset recurrence is
E(S,v) = sum over u in S minus {v} of A[u,v] E(S minus {v},u),
and H(D)=sum_v E(V,v). This generates/counts arbitrary binary adjacency matrices; loops never enter the recurrence.

To generate every labelled semicomplete digraph, add vertices in label order and, for each old vertex, choose exactly one of the three relations old->new, new->old, or both. There are 3^(n choose 2) such digraphs. For tournaments allow only the first two choices, giving 2^(n choose 2). Arbitrary off-diagonal binary matrices require a fourth state, neither direction, giving 4^(n choose 2); diagonal choices contribute an independent factor 2^n without changing H.

If an absent directed arc a->b is added, its exact increment is
Delta H = sum over S subset V with a in S and b not in S of E_D(S,a) B_D(V minus S,b).
The unique new arc splits each newly created path into these two parts. This applies in particular when adding the reverse of a tournament arc. Multiple reverse arcs can be added sequentially, recomputing the endpoint arrays each time. Interactions prevent replacing the increments by their values in the initial tournament.

If D has b doubled pairs and T ranges over its 2^b tournament resolutions, then
sum_T H(T) = sum over Hamiltonian paths P of D of 2^(b-k(P)),
where k(P) is the number of doubled pairs traversed by P. Thus resolution is a weighted identity, not preservation of H and not the unweighted sum of tournament counts. For the seven-path seed, the four resolution counts are 5,3,5,1; their sum is fourteen, not seven. No tournament resolution itself has seven paths.

## Scope

This note independently proves the generalized Camion reduction, SCC product and exact count identities. It uses the existing tournament floor theorem and its cited finite enumerations, rather than independently reproducing those censuses. Parent census results are marked as such. No novelty claim is made for classical structural lemmas.

## Parent's completed census and all-order corollary

The parent has now reported the completed exhaustive census:
- Three strong order-four count-seven isomorphism classes, each with exactly two reciprocal pairs.
- Sixteen strong order-five count-21 classes: six with two reciprocal pairs and ten with three.
- Sixteen strong order-six count-21 classes, all asymmetric: two with one reciprocal pair, eight with two and six with three.

Conditional only on the verified census artifacts, the structural proof above makes this a complete all-order catalogue of cores: three for seven and thirty-two for 21. Uniqueness follows from uniqueness of the SCC decomposition and its total ordering. The full count-seven family consists of a unique one of the three cores, with any number of singleton components before/after. The full count-21 family consists either of one of the thirty-two strong cores with arbitrary singleton prefix/suffix, or one of the two count-three cores (a directed triangle, or that triangle with one edge doubled) and one of the three count-seven cores in either order, with arbitrary singleton prefix, middle block and suffix. Relabel vertices arbitrarily. Independently choose every diagonal bit arbitrarily: the counts are unchanged.

For arbitrary digraphs (including loops), sharp necessary alternatives are therefore:
- H=7 implies at least one absent unordered pair OR at least two reciprocal pairs.
- H=21 implies at least one absent unordered pair OR at least one reciprocal pair.
The bounds on reciprocal pairs are attained by the semicomplete seeds. The diagonal is irrelevant and need not appear in either corrected statement.


## Correction: inert reverse arcs exist even in strong digraphs

The initial version incorrectly called the directed triangle the only count-three core. The corrected classification has two such cores. The independently inspected source is `/workspace/scratch/07d5532f2fce/semicomplete-scan-20261109/strong_h3.json`. Its matrices are 010/001/100 and 010/101/100. Both have exactly the three paths 012,120,201. In the second matrix the added arc 1->0 is inert: a three-vertex path using it would be either 102 (requiring absent arc 0->2) or 210 (requiring absent arc 2->1). Thus neither possible placement of the added arc completes a Hamiltonian path. The graph remains strong because the original triangle remains present.

The exact increment formula above explains inertness without any exception: all summands are nonnegative, and the increment vanishes exactly when every subset split has E_D(S,a) B_D(V minus S,b)=0. Strong connectivity supplies directed reachability, but does not require a Hamiltonian path to contain each available arc. Therefore one must never infer strict increase of H merely from adding a reverse arc, even to a strong tournament.

To check completeness of these two order-three cores directly, every strong semicomplete triple contains a directed triangle. With zero added reverse arcs it has H=3; one added reverse arc remains H=3 by cyclic symmetry and the calculation above; two added reverse arcs give H=4; all three give H=6. The strong-tournament floor excludes larger count-three cores, while the order-two strong core has count two.

There are consequently 2*3*2=12 ordered two-nontrivial-core types for count21, before singleton padding: choose one of two count-three cores, one of three count-seven cores, and either order. This corrects the six types that a triangle-only statement would have counted. The thirty-two strong count-21 cores, seven-core catalogue, reciprocal-pair lower bounds and diagonal independence remain unchanged.
