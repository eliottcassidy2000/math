# Functional trijections, their first obstruction, and the exact 45/69 root switch

## 1. A trijection with disconnected carriers

The inherited [trijection](../trijection_four_nodes_three_edges_20261020.md) gives exact complementary-pair charts for an oriented spanning path, a directed cycle plus one isolated vertex, and an oriented in-star or out-star. The reason each chart has size n−1 is that an n−1-edge carrier leaves exactly binomial(n−1,2) free pairs. But having the correct number of pairs alone does not ensure the actions of carrier symmetries agree.

Here is an explicit extension. Let D be a spanning oriented carrier with one vertex r of outdegree zero and every other vertex v having exactly one outgoing arc v→f(v). Such f cannot contain loops or two-cycles. Every connected component is a tree oriented towards r, if it contains r, or consists of a directed cycle with in-trees attached otherwise. Conversely, any undirected pseudoforest with n−1 edges has exactly one tree component and can be so marked and oriented. Indeed, for a component with t vertices, e=t−1 if it is a tree and e=t if unicyclic, so the total number of tree components is one.

The complete free-pair chart onto pairs in V−r is:

* Keep each free pair not containing r, with its orientation.
* For a free pair {r,v}, the carrier edge v→f(v) is not incident to r (otherwise f(v)=r and the pair would already be marked). Read the bit r→v as v→f(v) at the previously missing pair {v,f(v)} in the smaller chart; read the opposite bit as its reverse.

Every marked edge not containing r arises once from its unique tail v, and every unaffected pair is kept. This proves bijectivity on free pairs and orientations. The map is vertex-equivariant, since each tail and its unique successor and the root are all part of the mark. Thus tournament completions of a fixed carrier D up to its marked isomorphisms are precisely the orbits of Tour(V−r) under Aut(D,r). This includes the older path, cycle-plus-solo and in-star charts by choosing the sink or solo as r. The out-star is obtained by converse, not by pretending its center has outdegree zero.

For a cover of all vertices by three directed runs, the mark consumes n−3 edges, leaving binomial(n−1,2)+2 free bits. Connecting the three runs by two available end-to-start arcs makes one functional in-tree (a path), but those same two boundary arcs cannot be dropped and then counted as free neutral data. In the two-apex setting the two connectors are replaced by two apex vertices and four prescribed arcs: the edge deficit still matches a single spanning path on the enlarged vertex set. Different endpoint choices give different marked paths, so the resulting complement charts are not one common unconstrained coordinate system. More generally one can point the unique root tree and allow directed cycles in other components. This yields a genuine chart on a disconnected carrier; it does not claim that merely counting arbitrary arms recovers a prescribed Hamiltonian count.

## 2. Exact converse criterion: Euler charge on fixed vertices

For any rooted oriented n−1-edge carrier D, not assumed functional, put G=Aut(D,r). The free-orientation space is isomorphic as a G-set to Tour(V−r) **if and only if** for every odd-order g in G,

    e(D restricted to Fix(g)) = |Fix(g)| − 1.

The condition is an Euler charge χ=1 on each odd fixed subgraph, not just on the entire carrier.

To prove it, compare numbers of orientations fixed by every subgroup K of G. If K has even order, it contains a nontrivial involution. A pair swapped by that element cannot be a marked arc of D, because the involution would reverse it; thus it obstructs every fixed completion. It also obstructs the smaller tournament orientation, since r is fixed. Both fixed sets are empty. If K has odd order, it cannot reverse an unordered pair and return it to itself, so the fixed orientation number is 2 to the power of the number of K-orbits on the relevant pair set. In the full pair set of V, pairs containing r correspond to V−r. Hence the two orbit exponents agree exactly when the carrier edge orbits and nonroot vertex orbits agree.

For an odd element g, each fixed unordered pair has two fixed endpoints, so the edge and nonroot vertex fixed counts are respectively e(D[Fix(g)]) and |Fix(g)|−1. Equality of these for all odd g yields equality of orbit numbers for every odd K by Burnside. Conversely equality of orbit numbers for every odd cyclic subgroup yields equality of fixed counts elementwise by induction on the divisor order (group the Burnside average by order of the cyclic powers). Finally, subgroup fixed counts determine a finite G-set: descending over possible stabilizer orders determines the multiplicity of each orbit type G/H after the larger stabilizers have been removed. This proves the criterion. On a functional carrier, whenever v is fixed so is its uniquely determined successor, giving the required edge count immediately.

**First non-pseudoforest hostile.** Take n=7 and D=K4 plus three isolated vertices. Orient the K4 as a source pointing to a directed triangle and mark one isolated root. Let g rotate the triangle. The free pairs of V have 15 elements, with six fixed (the pairs among the four fixed vertices, none marked), hence nine orbits. The pairs of V−r have 15 elements with only three fixed, hence seven orbits. The fixed orientation counts 512 and 128 rule out even a nonlinear equivariant chart for this rooted deletion action.

This is the first non-pseudoforest obstruction of that kind by order. Through n=5 all n−1-edge graphs are pseudoforests. For n=6 the only non-pseudoforest is a diamond plus two isolated vertices: a component with at least two cycles needs at least four vertices and five edges, using the whole edge budget. Its unrooted automorphism group is a 2-group, so every nonidentity odd order fixed-point test is vacuous for any rooted orientation of it. Abstract equivariance at that size does not imply a simple vertex-tail chart.

## 3. An exact Burnside enumerator for the whole functional family

For reference, the automorphism quotient can also be counted in one sum before specifying the shape of the carrier. Fix r among n labeled vertices, and let g on the nonroot vertices have only odd cycle lengths, with multiplicities m_l. A root-fixed g-equivariant f is determined by one allowed target per source cycle of length l; a target must be fixed by g^l. This gives X_l = sum over d dividing l of d·m_d choices per source cycle after adding the root and excluding a self-loop. A two-cycle can arise only between two distinct g-cycles of the same length: then one must map to the other by a cyclic shift, and vice versa by its inverse, giving l choices for each chosen pair. The forbidden pairs of source cycles must be disjoint. Inclusion/exclusion gives the factor

    F_l = sum over 0 ≤ j ≤ floor(m_l/2) of
          (−1)^j · m_l!/( (m_l−2j)! 2^j j! ) · l^j · X_l^(m_l−2j).

Multiply across lengths. For a list of odd lengths l_i the number of g-orbits on nonroot unordered pairs is

    b(g) = sum_i (l_i−1)/2 + sum_(i<j) gcd(l_i,l_j).

The first term classifies pair-distances inside a cycle and the second counts joint orbits between two cycles. Thus Burnside on pairs (f,T') gives the orbit count of all functional marked tournaments as the average of [product F_l]·2^b(g), summing only over odd-order g in S_(n−1); even-order g contribute zero. For n=2,…,7 the resulting exact orbit counts are 1,3,26,440,16684,1312512. The count 26 at order four already differs from the older three separate charts: rooting a tree at an interior path vertex or at a leaf of a star gives additional marks. These are counts of marks, not new tournament Hamiltonian counts.

## 4. Why the exact jump from 45 to 69 needs a new root

For [S328's one-way Fermat seam](../session328_fermat_gluing_and_phase_repair_20261114/FERMAT_RANGE.md), a rooted source has responses H and L, with L≥H, and the output is H+2^j L. A single fixed rooted source cannot reach both 45 and 69 even if j is any positive integer. Indeed 24=(2^k−2^j)L with k>j implies 2^j divides 24 and (2^(k−j)−1) divides 24. The only possible differences are one and two. Difference one forces H=45−24=21, forbidden. Difference two forces H=45−8=37, but L≤4 contradicts L≥H.

Changing the marked root, without changing source order or source clock, realizes both. The strong five-vertex source with outgoing bit rows (2,8,3,21,7) has H=9 and Q₂=84. Its full L menu over zero-based roots is (10,9,16,15,25). With the **same** genuine Fermat block R₂ (H=5), root 1 gives 9+4·9=45 while root 3 gives 9+4·15=69. Both literal outputs are strong of order eight, with the same source Q₂=84 but different output Q₂, namely 1416 and 1968. This distinction matters when trying to continue the construction. The first value 1416 is also the [inherited strong order-eight clock minimum](../session327_clock_orbits_and_activation_20261114/OPTIONAL_SLOT_FLOOR.md). The exact root swap reaches the 0-based odd coordinate 34 of the count 69, from coordinate 22 of the count 45. It does not treat even 34 itself as a Hamiltonian count.

Exhausting literal rooted sources through order four shows that neither target is produced there by a genuine Fermat block: at these levels one need only test j=1,2,4, since j≥8 exceeds both targets. Arbitrary near-transitive j=3 does produce 45 from H=L=5 in order four, but is not a genuine Fermat exponent and still cannot reach 69. Thus order five is smallest for the particular common source construction with genuine Fermat steps.

## 5. Historical square-sum 34 and 45: do not mistake multiplicity for transitions

The repo's [square-sum cutoff study](../crossroads_family_20260926_squares.md) separates Hamiltonian cycle existence from usable exchange structure. A direct independent search confirms 11 Hamiltonian cycles in the square-sum graph on [34], yet no four-cycle exists through size 45. Thus no two of those Hamilton cycles can be joined by a single degree-preserving two-edge switch. At sizes 46, 68, 69 the ambient graph has respectively 1,2,2 four-cycles; 69 is not a new local birth in that graph. In the tournament seam this same neighborhood is governed by the exact rooted profile, rather than a sparse switch graph. These are different predicates and no equivalence is claimed. The common methodological lesson is that a witness count or an unmarked isomorphism class can conceal a rigid continuation obstruction.

Reproduce with python replay_trijection_and_ports.py. It verifies the chart equivariantly through order five, validates the fixed-function factors by direct enumeration through order seven, and independently checks the numerical controls and two literal rooted outputs.

## 6. Exact absolute count slice in the two-apex three-run model

A separate use of three runs goes beyond how many lattices or orientations there are. Start with a transitive core and two fixed apices a,b. For the apex order a then b, the endpoint data are the last point q1 of the first run, the first and last p2,q2 of the middle run, and the first point p3 of the last run. A core vertex x is a ternary digit (eligible in all three runs) precisely if max(p2,p3) < x < min(q1,q2) for some valid four distinct anchors. The arc requirements are q1→a→p2 and q2→b→p3. In particular a ternary digit is a vertex eligible for all three choices, not a vertex used three times.

The seven vertices a,b,x,p2,p3,q1,q2 are in one strong component: the cycles through x are x→q1→a→p2→x and x→q2→b→p3→x. If two core vertices are ternary digits in one fiber, both lie there, with the same anchors, giving an eight-vertex induced subgraph. The replay tests the seven-vertex and eight-vertex marked forms exactly: at a possible center x, both apices must choose distinct out-neighbors below x and distinct in-neighbors above x. This requires the two lower eligibility sets nonempty with union of size at least two, and the same for the two higher sets. It is sufficient as well: select distinct anchors at each side. The first and last runs’ allowed ranges overlap, so every remaining vertex has at least one permitted run regardless of position. The relative apex direction is immaterial for this nonempty middle and is included in enumeration.

The induced seven-vertex cases have smallest count 41, and the induced eight-vertex cases with two ternary digits have smallest count 105. No full tournament with a ternary digit can have fewer than 41 directed Hamilton paths; a second ternary digit in the same profile costs at least 105. Monotonicity under adding vertices follows from the usual one-vertex insertion: every old Hamilton path extends at one of its endpoints or across a pair where its directions to the added vertex change, and deleting the new vertex recovers the old path injectively.

There is a stronger exact **absolute** statement: among all orders in this named model, the Hamiltonian counts at most 69 admitting one or more ternary digits are exactly

    {41,51,55,59,65,67,69}.

Indeed the active strong component alone has count at least 41, so all other components must contribute 1 (factorization across ordered components). The active component has order at most eight by Busch's bound H >= 5^((n−1)/3) for strong tournaments. It retains the apices, all anchors and a flexible vertex, and its deletion of the apices is still transitive. Complete checks of both core sizes 5 and 6 give precisely the displayed set. Thus even if one allows arbitrarily many outside vertices, the exact number 45 has no ternary digit in this model while 69 can; 67 also has one although it is excluded by the separate all-order one-way near-transitive seam of S328. This separates two fundamentally different ways to reach a prescribed absolute endpoint. It does not assert that such a digit is required to realize every odd count.

The literal seven-vertex masks are (62,28,56,112,32,66,23) at H=67 and (94,60,120,112,32,1,50) at H=69. Decomposing over both apex orders and the exact four-anchor fibers gives the following multiset counts, also recovered independently by assigning each core vertex its run and grouping the legal paths:

| Exact target | size-1 fibers | size-2 | size-3 | size-4 | size-8 |
| --- | --- | --- | --- | --- | --- |
| 67 | 19 | 9 | 2 | 4 | 1 |
| 69 | 19 | 12 | 2 | 5 | 0 |

The two graphs have the same number of ternary-choice fibers but different exact Hamilton counts. Absolute arithmetic and physical coverage cannot be read solely from how many triple-eligible positions are available.

There is also an actual unbounded ternary family, not just a list of small hits. In the strong seven-vertex frame (62,28,56,112,32,2,55), replace vertex 2 by a transitive block of size b≥1. The compressed counts for one, two, three block visits are 41,32,8, with no fourth visit; these are checked by all permutations of the six external vertices and all marker gaps. The number of assignments to j visits is the number of surjections to j ordered slots, giving

    H(b)=41 + 32(2^b−2) + 8(3^b−3·2^b+3) = 1+8·2^b+8·3^b.

The b flexible vertices all lie between the same four anchors, so b_ternary≥b. The family reaches the proved absolute minima for b=1 and b=2; no all-b minimum is claimed. And the family is very sparse, which is why physically realizing a clean 2/3-smooth formula is not yet arbitrary absolute targeting.
