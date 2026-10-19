# Bounded seven-thin deformations of the K4 no-good-pair example

Date: 2026-10-19. Status: exact bounded finite result, with explicit witnesses. No general seven-thin good-pair theorem is claimed.

## 1. Question and precise result

The [twelve-point K4 example](k4_eight_thin_no_good_pair_gold_20261019.md) is eight-thin, has no good pair, but has weak Gold. Can adding a few relations lower its incomparability degree to seven while preserving the absence of good pairs?

For every set of one, two or three directed relations whose endpoints lie within either original six-element level, take its transitive closure and discard cyclic results. **Every resulting seven-thin poset has a very good pair**, hence weak Gold by the [attributed good-pair synthesis](good_pair_gold_synthesis_20261019.md).

There are 60 possible directed added relations, 30 within each level. The exact check covers 60+1,770+34,220=36,050 relation sets. Of these, 1,850 are cyclic and 240 repeat a closure already encountered. Among distinct acyclic closures, 33,336 still have maximum incomparability degree above seven. The remaining 624 are seven-thin, and every one has an explicitly recorded very good pair. None of the one- or two-relation patterns reaches seven-thinness; all 624 first occur with three additions.

This is a statement about those particular deformations, allowing chains and arbitrary closure consequences. It is not a census of seven-thin posets. Mixed patterns involving new cross-level relations, four or more added relations, vertex additions or deletions, and unrelated starting posets are outside scope. Cross-level additions alone stay height two and are already covered by the [height-two seven-thin theorem](seven_thin_height_two_very_good_gold_20261019.md) once they attain seven-thinness.

## 2. Why simple star additions create witnesses

Use the zero-based labels of the K4 example: lower star rows 0,1,2,3; lower matching rows 4,5; upper edge vertices 6,...,11. The original predecessor masks are

    (0,0,0,0,0,0,44,35,26,21,54,57).

The four stars correspond to vertices of K4 and the two matching rows to two chosen perfect matchings. Directed relations between star rows have two types under the automorphisms preserving the set of selected matchings: their common edge belongs either to a selected matching or to the unselected matching.

For the selected-edge representative, add 0<1. In the closure, lower vertices 0 and 4 have equal empty downsets. We have U(4)⊆U(0), and

    U(0)\U(4)={1,7}, with 1<7.

Thus (0,4) is a very good pair: the difference is a two-point chain, even though the original height-two incidence sets were separated.

For the unselected-edge representative, add 0<3. The maximal vertices 6 and 11 have equal empty upsets, while

    D(6)\D(11)={2}, D(11)\D(6)={4}.

They form a very good pair in the dual. Relabeling the K4 vertices covers every single directed star-to-star addition. These witnesses explain two mechanisms by which the no-good-pair incidence obstruction breaks: an upset difference becomes a chain, or transitivity makes two predecessor sets nearly equal.

Further additions can change these specific witnesses, so this single-relation observation is not used as an unproved monotonicity assertion. The bounded three-addition result independently tests the actual closed poset and records a valid witness for every seven-thin closure. The finite conclusion is exactly the 624-case result in Section 1.

## 3. Exact verification

Run `python -O ../experiments/verify_k4_three_relation_deformations_20261019.py` from the results folder. The standard-library program generates every bounded relation set in deterministic order, computes predecessor closure by bitset transitive closure, removes cycles and deduplicates exact predecessor tuples. For every seven-thin result it independently recomputes closure by graph reachability using sets.

The verifier constructs successor masks, checks actual incomparability degrees, and tests both P and its dual. A witness [dual,a,b] is accepted only when a,b are incomparable, their strict downsets in that orientation agree, and each directed upset difference is a chain under the actual transitive order. All 624 predecessor tuples, originating additions, degree lists and witnesses are stored in each receipt. The receipt also authenticates the deterministic record list by SHA256. Normal and optimized runs pass; guards remain active under `-O`.

No probability computation is required here: a very good pair can be oriented at probability at most one half, and the cited synthesis supplies weak Gold. The failure to find a no-good-pair example within this family does not establish the stronger conjecture that every seven-thin poset has a good pair. The original eight-thin example and the height-two theorem remain unchanged. No novelty claim is made for the finite controls or general methods.
