# What must a repair preserve?

Proof synthesis, 19 October 2026. Amended at 09:58 UTC with the height-three obstruction and its strict Gold certificate. This note separates established mechanisms from new research targets. It makes no historical-priority claim and does not settle the general Gold Partition, 1/3–2/3, Graceful Tree, Collatz, or planar Jacobian conjecture. The [companion frontier](liouville_20261019_overnight_frontier_after112.md) records exact scopes and remaining cases.

## 1. A precise transfer principle, and two different failure modes

Let a finite object set have an external-state map π. A π-preserving injection from A into B proves weighted(A)≤weighted(B) for every nonnegative weight constant on each π-fiber. A π-preserving bijection gives equality. More generally, a linear count inequality holds under every such weighting exactly when it holds on each fiber: necessity follows by concentrating the weight on one fiber, and sufficiency by summation.

This gives a useful proof audit, not a universal theorem connecting the conjectures:

1. Specify the admissible objects and the information retained by the repair.
2. Prove that every intermediate object remains admissible.
3. Prove the required inequality or invariant on the retained state.
4. Supply the input objects or trajectory guards independently.

Two opposite mistakes matter. Forgetting constraints enlarges the admissible measure or trajectory class and can invalidate a certificate. Freezing unnecessary data shrinks the move space and can disconnect it. A third operation, composing many steps into one, can preserve the endpoint formula while destroying the quantitative hypothesis needed for the proof.

## 2. Gold: retain the outside skeleton and recompute the measure

For the [one-relation twin gadget](one_relation_twin_gold_theorem_20261019.md), the swap and three-cycle change only three labels at their occupied positions. The outside word and marked positions are fixed. Thus A=C and B≤A hold on every skeleton, giving the exact weak Gold table [[A,B],[A,0]]. The [conditional composition theorem](disjoint_twin_gold_composition_20261019.md) consequently handles disjoint gadgets without assuming independent outcomes: each complete two-comparison block halves the surviving weight. Taking only one comparison from each gadget does not inherit Gold.

The [chain-surplus theorem](chain_surplus_gold_theorem_20261019.md) replaces three cells by insertion bins q₀=q₁≥q₂≥⋯. The [attributed good-pair synthesis](good_pair_gold_synthesis_20261019.md) needs only

- D(a)⊆D(b);
- U(b)∖U(a) is a chain;
- the current probability of a<b is at most 1/2.

Zaguia's insertion swaps give decreasing bins; adjacent cuts at the half-mass crossing give the two-comparison weak Gold strategy. This combines explicitly attributed arguments, without asserting priority. Under skeleton-invariant reweighting the bin inequalities survive, but the half-probability condition must be tested again. For a very good pair, equal downsets and chain differences in both directions allow orientation using the current measure.

The [initial height-two c=0 theorem](seven_thin_height_two_c0_very_good_pair_20261019.md) is now strengthened to [every non-total seven-thin height-at-most-two poset](seven_thin_height_two_very_good_gold_20261019.md), including isolated vertices. The proof uses actual all-vertex incomparability budgets: after choosing the larger level, a separated missing-incidence row family is too small, or forces a near column pair. Near means each directed set difference has size at most one, exactly the required very-good-pair condition. No critical-pair or common-neighbor restriction remains.

The [ten-point height-three obstruction](height_three_seven_thin_no_good_pair_20261019.md) now shows that the height-two structural conclusion cannot be extended even to height three: the example is seven-thin and has no good pair in either orientation or dual, hence no very good pair. Nevertheless, comparing (3,5) and then (2,4) gives a strict Gold certificate with minimum deficit 2,128 among 11,971 extensions. Good-pair injections are therefore sufficient, not necessary, for the two-comparison inequality. Ordinal sums of copies preserve both the absence of good pairs and the normalized Gold margin; their height grows. A repair program must admit certificates without a distinguished good pair, rather than mistake failure of that structural search for failure of Gold.

The endpoint boundary is concrete: two incomparable extra successors admit a six-point example with no Gold strategy confined to the four marked vertices, while an exterior endpoint gives a strict strategy. Keeping all witnesses inside the original gadget is too restrictive.

## 3. Square sums: retain positional parity, release a bounded edge set

In compatible interval-path pairs, every common vertex has the same positional parity. That is enough for either input member to represent the same labels under the affine lift. Internal ordering and traversal orientation may change; the certificates separately check coverage, connector squares, connectivity, and output parity.

Keeping every block intact was too restrictive. The complete construction permits specified 1–3 prefix cuts in intermediate states and 2–7 cuts in final states. The [exact-threshold theorem](seeded_optimal_threshold_6242_20261019/OPTIMAL_THRESHOLD_THEOREM.md) then gives a cycle containing 482–3362–6242 for every N≥6242, and no smaller N. Its finite gap certificate preserves one cyclic order, both seed edges, and the exact bad-edge set through each insertion/reversal; the eventual theorem supplies the infinite induction.

The [finite-cut retention theorem](square_sum_finite_cut_retention_20261019/FINITE_CUT_RETENTION.md) makes the cost explicit. For lift factor Q and retained parameters K,

    K′ = (union of [Qt−(Q−1)/2, Qt+(Q−1)/2] for t in K) minus D,
    |K′| ≥ Q|K|−|D|,

where D covers every prescribed copy whose actual edge a certified cut removes. Shared parity makes input-member choice harmless. This is a retention theorem, not a construction of missing compatible inputs.

The [all-size connected-path theorem](square_sum_all_size_connected_paths_20261019/ALL_SIZE_CONNECTED_PATHS.md) supplies those inputs for the prefix 1,3,6. It yields 2,303 fixed, pairwise vertex-disjoint triples in every interval N≥6,294,221, simultaneously with the seed. “All-size” means every sufficiently large interval size, not every prescribed path shape.

## 4. Graceful trees: retain closure, release some already-correct edges

The [rank-staircase theorem](graceful_upper_tail_two_color_connectivity_20261019/UPPER_TAIL_TWO_COLOR_CONNECTIVITY.md) works for disjoint matroid color classes with rank(Aᵢ)≥s−i+1. A recursive obstruction is repaired using a high-rank color. When two colors must change together, their selected set has the same closure before and after the move.

In the graphic case closure means the connected-component partition. Thus an unchanged connecting color remains valid after the two-color repair. For integer-difference colors, contracting any acyclic upper-tail choice leaves contact ranks at least r+1−d: distance-d edges have d residue chains. Every upper-tail completion fiber is therefore connected by at-most-two-color moves, with diameter at most 2^(r−1) for r free colors.

The [actual fixed-color counterexample](graceful_fixed_color_counterexample_20261019/ARITHMETIC_FIXED_COLOR_COUNTEREXAMPLE.md) shows the exact boundary. On labels 0,…,12, free colors 6,9,11 have five completions split into components of sizes three and two. The two contact assignments require all three colors to change. Global width-two connectivity still holds, but a connecting path must release at least one common, already-correct edge. The upper-tail hypothesis cannot be discarded. Tree shape may also change; this is not a prescribed-shape graceful-labeling theorem.

## 5. Poset balance: preserve exterior restrictions and let witnesses migrate

Nonnegative boundary increments express actual extension measures in a cone of core generators. A positive inequality on every generator transfers to every actual exterior. The reverse inference fails: a feasible relaxed mixture need not be an exterior, and failure to certify a core does not make a counterexample poset.

Three real constraints strengthen the transfer:

- [Saturation](c4_saturation_frontier_20261019/mixed12-c4-frontier-stage/seven_thin_saturation_increment_lemma_20261019.md): a core vertex already using all seven incomparable neighbors is comparable with the whole exterior, annihilating entire boundary-increment branches.
- [Directed defects](c4_analytic_reductions_20261019/README.md): exchange injections bound a comparison when its predecessor/successor defect is small; forced exterior comparabilities sharpen the bound.
- Maximum-span anchoring: no exterior vertex lies between the critical roots. If their global maximum span is at most three, the ambient poset is six-thin, so the explicitly cited [Peczarski theorem](https://doi.org/10.1007/s11083-008-9081-9) applies.

Even realizability is not enough for an overly narrow witness list. The [mixed-(2,2) realization](mixed22_c02_exact_exclusion_20261019/mixed22_core_realization_and_span_dependency_20261019.md) is an actual 16-point five-thin poset whose eight forward core comparisons all exceed 2/3; it has eight balanced exterior pairs. A core-only inequality excluding this actual measure would be false. The useful repair is to retain the global anchoring hypothesis or add exterior witnesses, not merely to demand product-form weights.

## 6. Collatz: retain guards and the expansion-then-reset factorization

The [bank reset theorem](collatz_bank_reset_20261019.md) uses expanding centers b/(a−1)≤C followed by contraction fixed points d/(1−r)≤P. Each bank's composition retains its center envelope. Once the total slope contracts, the normalized carry is less than C+P. Exact parity cylinders then turn the inequality into a uniform bound on excluded initial source lifts; independently guarded adjacent returns have the same cutoff at both endpoints. The [modular refinement](collatz_modular_bank_head_20261019.md) reduces the concrete bound from 534 to 107.

What must not be compressed away is demonstrated by [unbounded return centers](collatz_unbounded_return_centers_20261019.md). For genuine adjacent returns K^L F^floor(Lθ), the slopes lie between 1 and 4/3, every lower prefix stays above its source, yet a subsequence of compressed expanding centers tends to infinity as its slope tends to one. Thus no single bounded-center bank contains all these macroreturns. The original K/F factorization still has its small center and becomes contracting after one further lawful F block. Its enlarged guard cylinder must be checked anew. Neither presentation supplies arbitrary source entry.

## 7. Planar carriers: retain the divisor, then both faithful boundaries

The [single-carrier rigidity theorem](planar_carrier_factor_rigidity_20261019.md) retains the full factor support of an invariant polynomial, not merely one marked factor. An étale polynomial map preserving that polynomial permutes its irreducible factors. After an iterate, the compulsory pair t and 1+x²t is fixed up to scalars; determinant one then forces that iterate to be the identity.

A nonzero-time completed Hamiltonian flow with source c+p²(p³−y²)R has a nonzero leading source-adic jet, so has infinite order. It therefore cannot specialize to a polynomial source map. This is an exact incompatibility of preserved divisor and leading jet.

The [finite-composition theorem](planar_carrier_composition_rigidity_20261019.md) now removes the common-invariant restriction: every polynomial realization of a finite composition of these flows is the identity, allowing exact cancellation. Each flow fixes t=0 and a second boundary in the restricted common ring K[s][[τ]][τ⁻¹]. The faithful substitution s=xt,τ=t compares the actual source series; the continuous evaluation τ=−s² transfers fixedness to the original curve 1+x²t=0. Gwoździewicz's cited fixed-line theorem makes a hypothetical polynomial realization an automorphism; line fixedness reduces it to a shear, and the second boundary kills that shear. Independently, the completed cusp algebra intersects K[[τ]] only in constants, giving the same conclusion.

The ring restrictions are decisive: 1+x²t is a unit in the source t-adic completion, so apparent divisibility of approximants alone cannot justify the second boundary. General rational realizations and source-chart entry remain open here.

The common lesson is operational: preserve precisely what the next proof step uses, prove admissibility after each repair, and keep input supply and coverage as independent obligations.
