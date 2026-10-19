# Seven-thin Gold: module reduction, a positional obstruction, and a triple-exchange rescue

Date: 2026-10-19. Status: PROVED AND FROZEN AFTER ROOT REVIEW. The general seven-thin weak-Gold problem and the good-pair-or-two-fifths dichotomy remain open.

## 1. Results and precise scope

This note does not prove or refute the proposed dichotomy

> H7: every non-total seven-thin poset has a good pair in the poset or its dual, or a comparison with probability in [2/5,3/5].

It gives three concrete results instead.

1. A minimum-cardinality counterexample to H7 has no proper non-total module. The same reduction holds for a minimum-cardinality counterexample to weak Gold. Consequently a sufficient minimal-counterexample-only version of H7 needs to be checked only for posets with no proper non-total module and with transitive strict majority order.
2. Every family of label involutions that swaps a specified second pair is subject to a universal positional matching cap. On one prescribed seven-thin, no-good-pair example, this cap excludes that entire method for an otherwise valid Gold comparison tree whose first heavy probability lies strictly between 3/5 and 2/3.
3. On that same tree, legal three-cycles of labels give an exact matching certificate large enough to prove Gold. The receipt includes the matching itself, its generating cycles, and an equal-cardinality vertex cover. This is a finite exact improvement of the exchange method, not a structural proof that every seven-thin poset admits such a certificate.

No poset census is used. The example is the ten-point perturbation Q′ already certified in the [conditional Gold packet](conditional_gold_exchange_and_minimal_obstruction_20261019.md). Its original poset has comparisons near one half, so it does **not** test H7's hard global probability gap.

## 2. Module transfer theorem

A nonempty set M of vertices is a **module** if every vertex outside M has the same relation to all vertices of M: below all of them, above all of them, or incomparable with all of them. Let P[M] be the induced poset.

**Theorem 1 (uniform restriction and transfer).** Restrict a uniformly random linear extension of P to M. Its induced order is uniform on the linear extensions of P[M]. More precisely, every extension of P[M] has exactly the same positive number k of extensions of P lying above it under restriction, so e(P)=k e(P[M]). Therefore:

- every comparison probability and every joint comparison-event probability involving only vertices of M is unchanged;
- a good pair in P[M] or its dual is a good pair in P or its dual;
- a weak-Gold two-comparison strategy in P[M] is a weak-Gold strategy in P.

**Proof.** Fix two internal extension words w,w′ of P[M]. In a full extension whose restriction is w, keep every outside label and every occupied position fixed, and replace the labels in the M positions, in order, by w′. Internal relations are respected because w′ is an extension. An outside vertex has a uniform relation to M, so changing the internal assignment at the same positions cannot violate an outside relation. Replacing w′ by w is the inverse. This proves equal restriction fibers.

For a,b in M, write D_P(a)=D_M(a)∪D_out and U_P(a)=U_M(a)∪U_out, where the outside sets are independent of a. The conditions D(a)⊆D(b) and “U(b)\U(a) is a chain” are therefore unchanged between the induced and ambient posets. The relevant probability is unchanged by equal fibers. This proves good-pair transfer, including the dual.

For a fixed strategy supported in M, every survivor count t_i in P equals k times the corresponding survivor count in P[M]. Thus t_0≥t_1+t_2 transfers verbatim. The same argument would also transfer an adaptive strategy because its branch histories involve only internal comparisons. ∎

**Corollary 1.** If H7 has a counterexample, one with the minimum number of vertices has every proper module a chain. The same statement holds for a minimum-cardinality seven-thin weak-Gold counterexample.

Indeed an induced subposet of a seven-thin poset remains seven-thin. A proper non-total module would have fewer vertices; either alternative of H7, or its Gold strategy, would transfer to the ambient counterexample.

In particular these minimum counterexamples have connected incomparability graphs. For completeness, distinct connected components of an incomparability graph are uniformly comparable: for an outside vertex z, the orientation of its comparisons with a component cannot change across an incomparability edge, because x<z<y would force x<y. Connectivity propagates the orientation. Thus the components are modules linearly ordered by the poset, giving an ordinal sum. A non-total component of a disconnected counterexample would be a proper non-total module.

This is not a reduction to completely indecomposable posets: proper **chain** modules are not ruled out. Contracting such a chain need not preserve the relevant comparison probabilities, so no contraction theorem is being assumed.

## 3. Why the existing ordinal-sum and blow-up families do not settle H7

The ordinal-sum families from [the ten-point witness packet](height_three_seven_thin_no_good_pair_20261019.md) and [its conditional extension](conditional_gold_exchange_and_minimal_obstruction_20261019.md) preserve every within-block comparison probability. Repeating either witness or placing chains between copies therefore preserves the already present [2/5,3/5] comparison. It cannot manufacture a counterexample to H7.

There is also a precise obstruction to the most immediate way of biasing the witness by vertex weights.

**Lemma 2 (chain substitution budget).** Replace a vertex v of P by a chain of m vertices, each having v's old relation to every other vertex. Each clone has incomparability degree d_P(v). An old vertex u incomparable with v acquires incomparability degree d_P(u)+m−1, and all other old degrees are unchanged. For an initially seven-thin P, this substitution remains seven-thin exactly when

    m−1 ≤ min_{u incomparable with v} (7−d_P(u)),

with an empty minimum imposing no restriction.

This follows by counting the replaced incomparable neighbor. The two prescribed seeds have predecessor masks

    Q  = (0,0,0,6,0,19,2,18,86,151),
    Q′ = (0,0,0,6,0,19,2,82,86,215).

In both seeds, vertices 0 and3 have degree7. Every vertex is incomparable with at least one of0 and3: vertex 0 is comparable only with5 and9, and each of5 and9 is incomparable with3. Hence every nontrivial chain substitution at **any** seed vertex violates seven-thinness. The receipt checks all ten substitution sites in both seeds.

This is only a seed-specific obstruction. It says nothing about other seven-thin constructions, about adding or deleting relations simultaneously with a substitution, or about the truth of H7.

## 4. Heavy-branch matching certificate

Use the weak Gold convention from [the conditional Gold packet](conditional_gold_exchange_and_minimal_obstruction_20261019.md): comparisons may be repeated or already determined. Let E=e(P). Choose a first comparison with a larger branch R containing b extensions, b≥E/2. Let u,v be a proposed second pair. Partition Ω=L(R) into

    Ω+ = {L: u precedes v},   Ω− = {L: v precedes u},

with cardinalities c and b−c. The exact heavy-branch criterion is

    max(c,b−c) ≤ E−b,
    equivalently |2c−b| ≤ 2E−3b.                 (1)

The smaller first branch is automatic, since two survivor counts there sum to at most2(E−b)≤E. Consequently the same second comparison can be used on both branches.

Let G be **any** bipartite graph with parts Ω+,Ω−, whose edges are certified allowable exchanges between full extension words. No group action or independence is required.

**Theorem 3 (matching and fractional matching certificate).** A matching with ν edges implies

    |2c−b| ≤ b−2ν.

More generally, assign nonnegative weights to edges so that the sum of incident edge weights at each vertex is at most1. If ν is the total edge weight, the same inequality holds. In either case, if

    ν ≥ 2b−E,                                  (2)

the first and second comparisons form a Gold strategy. Every heavy second cell has Gold margin at least E−2b+ν.

**Proof.** Every edge has one endpoint of each orientation. The capacity condition implies c≥ν and b−c≥ν. Hence max(c,b−c)≤b−ν and |2c−b|≤b−2ν. Inequality(2) makes b−ν≤E−b. The stated margin is E−b−(b−ν). ∎

A label involution σ with σ(u)=v, σ(v)=u gives disjoint edges L↔σL whenever both words are extensions of R. If M extensions lie in L(R)∩σL(R), it gives ν=M/2 and recovers precisely the single-overlap theorem from packet 136. Multiple involutions can be combined, but overlapping domains are not automatically additive: the per-extension capacity constraint is essential. Fractional weights or a verified matching are valid ways to enforce it.

This is elementary matching logic. No novelty claim is made for the graph or fractional-capacity argument.

## 5. A universal positional cap for pair-swapping involutions

Number word positions from0. For i<j, let

    A_ij = #{L∈Ω: position(u)=i, position(v)=j},
    B_ij = #{L∈Ω: position(v)=i, position(u)=j}.

**Theorem 4 (positional cap).** In the union of the exchange graphs of **any** family of label involutions swapping u and v, every matching or fractional matching has total weight at most

    ν_pos = Σ_{i<j} min(A_ij,B_ij)
          = (b−d_pos)/2,
    d_pos = Σ_{i<j} |A_ij−B_ij|.                (3)

Therefore this entire method cannot meet the Gold matching threshold if

    d_pos > 2E−3b.                              (4)

**Proof.** Applying any such involution exchanges the occupied positions of u and v exactly. Other labels may move, but the unordered position pair {i,j} is invariant. Thus every graph edge joins the A_ij cell to the B_ij cell for one fixed position pair. Capacity limits the weight in that block by its smaller cardinality. Summing gives(3), and(4) is equivalent to ν_pos<2b−E. ∎

The bound does not claim attainability: preserving the two positions is necessary, not sufficient, for two words to be connected by a permitted label involution. It also does not obstruct arbitrary relabelings, swaps involving only one of u,v, or three-cycles. Those can change the unordered pair of their positions.

The actual global imbalance is |Σ(A_ij−B_ij)|, which may be much smaller than Σ|A_ij−B_ij|. Cancellation between position blocks explains why an otherwise Gold pair can fail this method's certificate.

## 6. One hard-window tree: exact obstruction

Take Q′ from Section3. Its incomparability degrees are

    (7,3,6,7,5,6,5,5,5,3).

It has E=7249 extensions and no good pair in either orientation of the poset: exactly42 ordered candidates satisfy the structural conditions, and their smallest comparison numerator is3804>7249/2. These properties were proved in packet 136 and are independently rechecked here by the ideal recurrence.

Choose the first comparison2<6. Its heavy child R has predecessor masks

    R = (0,0,0,6,0,19,6,86,86,215)

and b=4710 extensions. Thus

    3/5 < 4710/7249 < 2/3,
    2E−3b = 368,
    2b−E = 2171.

Choose the second comparison3<7. The four joint counts, with first/second affirmative orientations as above, are

    [[2422,2288],
     [ 693,1846]].

The heavy cells have imbalance134≤368. All four Gold margins are positive:

    (117,251,4017,2864).

Nevertheless the position table gives

    d_pos = 2034,    ν_pos = 1338 < 2171.

Thus **no collection of pair-swapping label involutions**, including arbitrary additional disjoint transpositions of other labels, can certify this otherwise valid Gold tree via a matching or fractional matching.

There is a short structural explanation that already suffices without the complete table. In R,

    D_R(7)={1,2,4,6},   7<9.

Hence7 cannot occupy positions 0,1,2,3 or9. But3 occupies position 2 in208 extensions, position 3 in476, and position 9 in674. Any relabeling that literally swaps3 and7 would put7 in an impossible position on these disjoint sets. Therefore at least

    208+476+674 = 1358

extensions are isolated under **every** such exchange, already far above the allowed defect 368. The full position table sharpens1358 to2034; it is stored in the receipt. These are exact finite counts, checked by enumerating all 4710 child words, with total and comparison-cell counts independently checked by ideal dynamic programming.

This proves a limitation of a specified proof method. It is not a counterexample to Gold, to H7, or to the existence of a different pair with a good involution certificate.

## 7. Bounded-support exchanges rescue the same tree

For a label permutation σ and a child word L, retain the directed exchange L→σL only if σL is a legal extension of R and the second orientation flips from3<7 to7<3. The domains are explicitly restricted by these two conditions. A three-cycle is not being treated as an involution or as an automorphism of R.

On this one child, form the following finite graphs. The parts always have cardinalities2422 and2288.

| Allowed label permutations | Edges | Matching size | Equal-size cover, left+right |
|---|---:|---:|---:|
| All transpositions meeting {3,7} | 3210 | 1360 | 1120+240 |
| (3 7), plus three-cycles containing both 3 and7 | 6952 | 2086 | 1706+380 |
| All three-cycles meeting {3,7} | 12870 | 2206 | 1286+920 |

The two restricted families are included to locate the concrete improvement. The first permits swaps involving just one endpoint, so it is **not** subject to the pair-swapping positional cap. Even so, its best matching is too small. The second permits three-cycles but insists that both second-pair labels move; it too is insufficient. The successful family permits three-cycles involving only one endpoint as well.

For each row the receipt includes an explicit matching and a vertex cover of exactly the displayed size. A cover meets every graph edge; consequently any matching has at most its size. The equal-size certificates therefore prove the displayed optima without relying on an external matching theorem or on the search algorithm's asserted optimality.

The third row alone is enough for Gold: its ν=2206 exceeds the required2171 and gives

    certified defect b−2ν = 298 ≤ 368,
    certified heavy margin E−2b+ν = 35 > 0.

The actual heavy margins117 and251 are better. The exchange certificate need not recover the exact imbalance.

All three-cycles in this successful graph have support at most3 and change no positions outside those three occupied label positions. The supplied matching assigns disjoint source and target words, so it is a rigorously restricted collection of local injections. Its domains are finite certified sets rather than a uniform neighborhood formula valid for arbitrary posets. This is the remaining compression limitation of the example.

The mechanism is consistent with the older [one-relation twin Gold proof](one_relation_twin_gold_theorem_20261019.md): its injection from u<t<v to v<u<t is a three-cycle of labels on a restricted domain, rather than a label involution swapping a fixed pair. The [conditional twin composition theorem](disjoint_twin_gold_composition_20261019.md) exploits the same preserved outside skeleton. This note does not infer Gold from a general one-sided directed-defect bound; it supplies the full matching capacity required for this particular heavy branch.

## 8. A weaker minimal-counterexample-only target

The [conditional Gold packet](conditional_gold_exchange_and_minimal_obstruction_20261019.md) proves that a seven-thin weak-Gold counterexample minimal by number of incomparable pairs has:

- no good pair;
- no comparison in [2/5,3/5];
- no tied comparison and a transitive strict majority tournament.

The last assertion uses the classical cyclic-majority Gold criterion, stated by Peczarski in [Comments on the Golden Partition Conjecture, pages 108–109](https://cdm.ucalgary.ca/article/download/62349/46785/178291); no novelty is claimed for it. The 3/5 reduction uses Gold in a non-total conditioned child to obtain a 1/3–2/3 comparison there. A total larger child has one extension and forces the original E=2 case, which is already Gold.

Combining that reduction with Theorem1 gives a smaller sufficient target than H7:

> H7-restricted (OPEN): every non-total seven-thin poset with no proper non-total module and with transitive strict majority order has a good pair, or a comparison in [2/5,3/5].

**Proposition 5.** H7-restricted would imply weak Gold for all seven-thin posets.

**Proof.** If a counterexample exists, first minimize its number of vertices, then its number of incomparable pairs. Vertex minimality and Theorem1 exclude a proper non-total module. A good pair would imply Gold by the [good-pair theorem](good_pair_gold_synthesis_20261019.md). The classical cyclic criterion makes the strict majority order transitive. If H7-restricted supplies a [2/5,3/5] comparison, its larger child has b≤3E/5. Orienting the comparison preserves seven-thinness and strictly decreases the number of incomparable pairs. Unless total, the child has Gold by the second minimality choice, hence has a 1/3–2/3 comparison. Using that comparison next leaves each heavy cell at most2b/3≤E−b, proving Gold. If the larger child is total, b=1 and E=2, resolved by comparison and repetition. Every case contradicts the counterexample. ∎

This restricted hypothesis may also fail and is not supported as a universal conclusion by the selected examples. It is genuinely narrower in its stated assumptions than H7, but this note does not exhibit a poset separating their truth values.

If either dichotomy fails, Theorem3 still permits a direct branch-specific alternative. One needs a first heavy branch with b≤2E/3 and a second pair whose certified exchange matching has ν≥2b−E. For b/E in (3/5,2/3], a mere 1/3–2/3 comparison in the child is insufficient; its larger conditional cell must be at most(E−b)/b. The example shows that local triple exchanges can meet this sharper requirement when every pair-swapping-involution cover fails.

No theorem here forces the existence of such a root, second pair, or triple-exchange capacity in the global no-good-pair, no-two-fifths gap class. That existence problem, now also restricted to no proper non-total module and transitive majority order in a minimum Gold counterexample, is the unresolved structural obstruction.

## 9. Reproduction and status boundaries

Run the [self-contained verifier](../experiments/verify_gold_exchange_hierarchy_20261019.py) normally and with `python -O`. It uses only the Python standard library and emits a deterministic JSON receipt without modifying source files. The [normal](gold_exchange_hierarchy_controls_20261019.json) and [optimized](gold_exchange_hierarchy_controls_optimized_20261019.json) receipts agree.

The exact controls include:

- the specified ten-point seed and child, using ideal counts and full word enumeration;
- ambient seven-thinness and every structurally eligible good-pair orientation in the poset and dual;
- the full positional table and the shorter rank obstruction;
- every matching edge, its actual label cycle, its orientation flip, endpoint capacities, and a same-size cover of every graph edge;
- the two saturated seeds' chain-substitution budgets and small independent module-fiber controls.

The broad search question was deliberately replaced by a precise method obstruction and a concrete bounded-support rescue on one existing child. Neither the witnesses nor the controls prove H7, H7-restricted, or general seven-thin Gold. No historical novelty, minimum-order, or universal support-three-exchange claim is made.
