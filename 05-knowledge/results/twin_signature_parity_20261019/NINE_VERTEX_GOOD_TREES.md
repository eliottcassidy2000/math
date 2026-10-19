# One isolated vertex changes parity: the good trees on nine vertices

19 October 2026. **PROVED:** the eight-vertex padding theorem, the structural nine-vertex criterion, and the two explicit odd-parity families. **FINITE-EXACT:** the 36-of-47 census and orientation counts, with complete controls. Independently [reviewed](review/CHERRY_EXTENSION_REVIEW.md). No historical-priority claim is made. This is a separate addendum to the frozen [TWIN_SIGNATURE_PARITY.md](TWIN_SIGNATURE_PARITY.md) and [TWO_CHERRY_PARITY.md](TWO_CHERRY_PARITY.md) packets.

All embedding counts below are spanning: a graph on n vertices is embedded bijectively into an n-vertex tournament. A tree is *symmetric* if its undirected automorphism group is nontrivial. A graph is *good* when every orientation has constant embedding parity across tournaments of that order. “Odd” below means odd count in every such tournament, not an asserted exact integer count.

The inherited inputs from [THM-4533](../../../01-canon/theorems/THM-4533-redei-graphs-parity-of-shaved-tournaments-and-u9.md) and its [Rédei note](../procgen_shave4_20261001_redei_graphs.md) are:

1. Every tree on at most six vertices is good, since the smallest asymmetric tree has seven vertices and hereditary symmetry suffices.
2. Every orientation of P4 has odd embedding count (the Forcade path formula, or its power-of-two consequence).
3. Deletion–reversal: phi_S + phi_(reverse_e S) = phi_(S-e).
4. For disjoint rigid components of orders n_i and parity values c_i, the union has constant parity (n! / product n_i!) product c_i.
5. An even oriented automorphism group forces zero; an asymmetric underlying graph with at least one edge has nonconstant embedding parity, by its nonzero highest-degree Boolean monomial.

## 1. The eight-vertex padding theorem

**Theorem 1.** For an undirected tree T on eight vertices, all orientations of T disjoint union one isolated vertex have the SAME parity function on nine-tournaments. That function is constant if and only if T is symmetric.

**Proof of orientation independence.** Delete any edge of T, and write its component orders as a,8-a with 1<=a<=4. The extra isolated vertex remains.

- If a=1, the forest has two isolated vertices, whose swap forces zero.
- If 2<=a<=4, each component has at most six vertices and is good. Their union with the isolated vertex has parity multiplied by

      9! / (a!(8-a)!) = 9 binom(8,a),

  which is even: its values for a=2,3,4 are 252,504,630.

Thus every edge-deletion parity function is identically zero. Deletion–reversal makes every single orientation flip preserve the function. All orientations are connected by such flips.

**Proof of the constant/nonconstant alternative.** If T is asymmetric, T plus its unique isolated vertex is also asymmetric. Its seven-edge highest-degree Boolean coefficient is one, so the common function is nonconstant.

As established by the rooted-branch decomposition in the core packet, every nontrivial rooted-tree automorphism group contains an involution: either an isomorphic branch pair can be swapped, or the statement descends to a nontrivial child group. If T has an involution that fixes a vertex, that involution reverses no edge. Choose an invariant orientation and obtain zero by its free automorphism action. The common function is therefore zero.

The remaining possibility is that every nontrivial automorphism swaps the ends of the central edge. The two rooted halves are isomorphic, and each has four vertices and trivial rooted automorphism group; any nontrivial rooted automorphism would extend to a vertex-fixing one. Orient the two halves compatibly and direct the central edge either way. The half graph is good, with some constant parity c. Deleting the central edge yields two isomorphic oriented halves; swapping them interchanges the two possible directions of the central edge. Hence the count for the chosen central-edge orientation is exactly half the disconnected count.

Equivalently, sum over unordered partitions of the eight host vertices into two four-sets. On each partition, the product of the two half embedding counts has parity c^2=c. There are binom(7,3)=35 such partitions. Thus the chosen eight-tree orientation has constant parity 35c. After adding the isolated source vertex, sum over its nine possible images: the padded graph has constant parity 9*35c. This proves constancy in the remaining case. It does not assert an exact count of 35 or 315 in an arbitrary host. ∎

The full statement is about the FUNCTION on tournaments. An asymmetric eight-tree remains nonconstant after padding, even though changing its edge directions no longer changes that function.

**Exact values for the symmetric branch.** The only four-vertex un-oriented trees are P4 and the three-leaf star. Every rooting of the star has a nontrivial rooted automorphism, while P4 has exactly two rooted isomorphism types (endpoint and internal vertex), both rigid. Thus the only symmetric eight-trees with no vertex-fixing involution are P8 and the internally rooted double D described in Section 3. Their padded parity is one, since every orientation of each P4 half is odd. Every other symmetric eight-tree has padded parity zero. This is a structural two-object classification, not just a census observation.

## 2. A complete structural criterion at order nine

**Proposition 1b (universal orientation linearization at nine).** For EVERY nine-vertex tree U and every nine-tournament T, the embedding parity is an affine function of U's eight orientation bits. Goodness is not assumed.

To prove this, delete any two distinct edges of U. The resulting three component orders are positive and sum to nine. If the largest is at least seven, the orders must be 7,1,1, and the two isolated vertices force zero. Otherwise every component has at most six vertices and is good. The disjoint-union parity is its three component constants times the multinomial coefficient. That coefficient is even: an odd multinomial requires the positive component orders to have pairwise disjoint binary digits, whereas nine has only two nonzero binary digits and there are three positive parts. Thus EVERY orientation of every two-edge-deleted forest has zero parity.

Applying deletion–reversal twice identifies this zero function with the corresponding mixed second orientation derivative. Vanishing of all mixed second derivatives makes a Boolean-cube function affine. Equivalently, for any reference orientation S0,

    phi_S = phi_S0 + sum over flipped e of phi_(S0-e),

as FUNCTIONS on nine-tournaments. In particular, each first derivative is independent of the directions of the other edges. The coefficient functions can still be nonconstant; affineness is not goodness. This is why the first-cut tests below suffice after one constant anchor is found.


**Theorem 2.** A nine-vertex tree U is good if and only if all three conditions hold:

1. U is symmetric.
2. Deleting any leaf of U leaves a symmetric eight-vertex tree.
3. Whenever deleting an edge splits U into orders two and seven, the seven-vertex component is symmetric.

**Necessity of 1.** This is the highest-degree Boolean obstruction in the inherited result.

If U is good, every orientation of U-e is rigid: choose an orientation of the remaining forest and extend it in the two ways across e; their two constant parity functions sum to that forest's parity. For a leaf edge, U-e is an eight-tree plus one isolated vertex. Theorem 1 therefore proves condition 2.

For condition 3, suppose the seven-component A is asymmetric. The other component is P2. For ANY orientations, the forest A disjoint union P2 has nonconstant embedding parity. Here is the needed coefficient proof, rather than an unjustified appeal to the full underlying graph (which has an edge-swapping undirected automorphism).

Write tournament orientations in Boolean variables y_ij, with every directed-edge indicator a literal y_ij or 1+y_ij. Fix a monomial whose six unordered edge variables form a copy of the connected tree A on seven of the nine host labels. In the expansion of the embedding polynomial, the six selected source edges must be exactly the six edges of A: a selection using the P2 edge is disconnected and cannot map to that connected monomial. Since A is asymmetric, its un-oriented isomorphism to the chosen host copy is unique. The two remaining source vertices can map to the two remaining host labels in two orders. Their two P2 literals sum to one. The coefficient of this six-edge monomial is consequently one, independent of the chosen source orientation. This is a nonconstant function, contradicting rigidity of every U-e orientation. Thus condition 3 is necessary.

**Sufficiency.** Consider any edge deletion U-e, with component orders a,9-a and a<=4.

- For a=1, condition 2 and Theorem 1 make the forest rigid.
- For a=2, condition 3 makes the seven-tree hereditarily symmetric: all of its proper subtrees have at most six vertices and are symmetric. It is good. The P2 component is also good, so their union is rigid (in fact its shuffle coefficient binom(9,2)=36 is even).
- For a=3 or4, both components have at most six vertices and are good, so their union is rigid (the shuffle coefficients 84 and 126 are even).

Therefore an orientation flip of any one edge changes the embedding-parity function only by a constant. It remains to supply one rigid orientation. A nontrivial automorphism of an odd-order tree cannot interchange the two halves of a central edge, because that would make the total order even. Some rooted branch permutation therefore supplies a vertex-fixing involution, which reverses no edge. An invariant orientation has identically zero embedding parity. Since orientation flips connect all orientations, every orientation differs from this zero function by a constant. Thus U is good. ∎

**Corollary 2a (the complete orientation selector).** For a good nine-tree U, call a leaf edge critical if deleting its leaf leaves P8 or D. Choose an orientation S0 invariant under a vertex-fixing involution of U; it has parity zero. Then every orientation S satisfies

    phi_S = number of critical edges reversed from S0, modulo two.

Indeed all non-leaf edge-deletion forests have even shuffle coefficient, hence parity zero. By Theorem 1 and its exact-values refinement, a leaf edge-deletion forest has parity one exactly at a critical leaf. Deletion–reversal therefore gives the displayed formula along any sequence of flips. The invariant anchor fixes the constant term; changing that anchor does not change the result.

This criterion is a fixed-order theorem. It is not asserted at other orders: the cut sizes, the eight-tree padding theorem, and the even shuffle coefficients are load-bearing.

The exhaustive leaf-augmentation tree generator gives 47 unlabelled nine-trees. Exactly 36 satisfy the criterion, and exactly 24 of those are hereditarily symmetric. The 12 remaining good trees include three detected by the original A(U)=0 automorphism filter, three further two-cherry trees, four additional universally-even trees, and the two odd-orientation types in Section 3. The eleven excluded classes have explicit different-parity tournament witnesses. This is a complete census of GOOD nine-trees, unlike the earlier filter-only census.

## 3. Two explicit nonhereditary trees with Rédei orientations

There is an exact general operation converting a padding statement back to a tree. Let H have a pendant leaf x at v, and let U add a second pendant leaf y at v. Then

    U is good iff H disjoint union one isolated vertex is good.

Indeed equal-direction extensions have parity zero. For opposite directions, reversing the new edge produces an equal-direction zero term, and deleting it produces the corresponding orientation of H plus isolated y. Every orientation of H occurs among these opposite-direction extensions. Thus both implications hold, with no loss of the orientation quantifier. This is an embedding-parity equivalence at order |H|+1, not a claim about automorphism groups.


Let D have vertices 0,...,7 and edges

    01, 02, 03, 14, 15, 26, 47.

Its central edge is 01. Its rooted halves are the paths 3-0-2-6 and 5-1-4-7, rooted internally at 0 and 1. Let W be D together with an isolated vertex8.

**Proposition 3.** Every orientation of W has odd embedding count in every nine-tournament.

**Proof.** The orientation-independence part of Theorem 1 applies. For an anchor, orient the two rooted P4 halves compatibly and orient 01. Every oriented P4 has odd embedding count. On any eight-tournament the matching-half argument sums odd products over 35 unordered four-plus-four partitions, hence gives odd parity. The nine choices of the isolated-vertex image preserve odd parity. By orientation independence, all orientations of W are odd. ∎

Now form two nine-trees:

    U_c = D plus the edge 08;
    U_a = D plus the edge 28.

They are nonisomorphic: U_c has a degree-four vertex, while U_a has maximum degree three. Each contains D induced, and D contains the asymmetric S(1,2,3) obtained by removing vertex5. Thus neither is hereditarily symmetric.

In U_c the distinguished leaf pair is {3,8} at 0; in U_a it is {6,8} at 2.

**Theorem 4 (exact parity selector).** For either tree and any orientation S:

    phi_S = 0 if the distinguished two leaf edges point the same way;
    phi_S = 1 if they point opposite ways.

**Proof.** Equal directions admit the leaf transposition. For opposite directions, reverse the new edge incident to8. The reversed graph has equal directions and parity zero. Deletion–reversal leaves D plus the isolated vertex8, whose parity is one by Proposition3. ∎

The undirected automorphism group of each U is exactly C2, generated by that leaf transposition. For U_c, the unique degree-four vertex is fixed and its other two rooted branches have distinct shapes. For U_a, the unique vertex with two leaf neighbors is fixed, then the chain of degree-three vertices and all unequal outer branches are fixed. Only the displayed leaf pair remains interchangeable.

Each has 256 labelled orientation masks. Exactly 128 have opposite leaf directions; all of those have trivial oriented automorphism group. The C2 relabeling action is free on these masks, giving exactly 64 isomorphism classes of Rédei orientations. The other 128 masks are individually fixed by the leaf swap and give 128 zero-parity orientation classes. Thus each tree has 192 oriented isomorphism classes in total. These are exact class counts, not an inference from a count of labelled masks alone.

Consequently the possible repair “a good tree with some Rédei orientation must be hereditarily symmetric” is also false, at the same minimal order nine. This is a stronger failure than a universally-even counterexample.

## 4. Verification and scope

The companion `verify_nine_vertex.py` regenerates all 47 nine-trees, checks the structural criterion and its 36-good/11-bad split, and verifies the eleven supplied nonconstancy certificates by exact tree embedding counts. It also checks all 1,316 two-edge-deletion forests underlying the universal affine-orientation lemma. It verifies every orientation of W and both extensions on the transitive host by an independent linear-extension recurrence, the automorphism and class counts, all edge-cut sizes and shuffle factors, and the asymmetric subtree.

Run normally and with Python optimization; both receipts must agree. Exhaustive source-side orientation checks confirm the selector, but do not pretend to enumerate all nine-tournaments. The universal host quantifier comes from the proofs above. The earlier exact minimality witness packet supplies the through-eight boundary.

The general classification of good trees remains open. Theorem 2 is an order-nine classification, and Sections 1/3 isolate the padding and central-halving mechanisms that made this order tractable. No implication for planar Jacobian geometry, quadratic dynamics, or a general extremal Rédei-growth conjecture is claimed.
