# Independent audit: no induced H13-to-H27 extension

**PASS.** The final theorem excludes an induced H13 tournament inside any H27 tournament, with no bound on the original number of vertices and no external strong-order estimate. It does not exclude either scalar count. It applies to any append-only graph operation, not only a port grammar or a recurring clone move.

## Final frozen binding

Source: `session197-hybrid-constructions/thirteen-to-twentyseven/`.

- PROOF.md: 10dd126f64a037589ce2dcf4387c02c4ff58e200b21db9c594e91d076320342c
- verify.py: 26cd39c6637aa7af4dada6048415c777bb820de25a8fa53627cea5fa2c76973b
- VALIDATION.json and VALIDATION_OPTIMIZED.json: 79333a563ebd7383c69bf596701ddccff18f9c5e7163d13bd2d62655f2f70163

Copied normal and optimized Python replays exactly match the frozen receipt. No frozen source was edited by the reviewer. The author incorporated the independently checked empty order-eight frontier and clarified the word “nontrivial” before this final binding. The previous external-bound argument is not needed by the final proof.

## Mathematical audit

The SCC reduction is exact. Hamilton-path counts multiply over the ordered strong components. A nontrivial strong component has at least three paths, so a count13 tournament has exactly one strong H13 component and only singleton components elsewhere. Its induced image remains strong and lies in one target SCC. Induced path monotonicity gives that target component count at least13, while the SCC product makes it a divisor of27. Hence its count is27. Arbitrarily many singleton components in either original graph do not evade this argument.

The strong-deletion lemma is sound for all orders at least four. Take a largest proper strong subset S. Unless it already has order n-1, every outside vertex is uniform toward S. Strongness forces an edge from the side dominated by S to the side dominating S. Adding this pair gives a larger strong subset. If precisely two outside vertices remain, deleting any one vertex of S instead leaves a strong n-1 subset: each retained S vertex lies on a directed triangle through the same outside pair, and these triangles form one strong component. This establishes the required predecessor at every order.

Every mixed birth from a strong graph is strong. Conversely the deletion lemma supplies such a mixed predecessor for every nontrivial strong graph of larger order. Hamilton-path monotonicity ensures a cap27 graph's whole chosen deletion chain remains within the cap. Thus generation under the cap is complete, and the source's degree-class canonicalization is exact because isomorphisms preserve outdegrees and every permutation within equal-degree classes is tried.

The empty order-eight frontier is decisive: a hypothetical larger nontrivial strong cap27 graph would have a strong cap27 induced deletion chain down to order eight, contradicting the complete empty frontier. This is closed induction, not extrapolation from a finite window. It removes any need for a Busch bound or a guessed maximum graph order.

The final count-only containment test is sufficient. The complete source H13 universe has one five-vertex core. Any induced copy of it in a target H27 core would be one of that target's five-subsets and have path count13. Absence of13 from all those subset counts excludes every such copy; the test does not assume a particular labelling of the source.

## Independent regeneration

I rebuilt the entire cap27 frontier from the triangle, using literal Hamilton vertex permutations rather than the source's subset DP. Accepted graphs are canonicalized by all actual Hamilton cycles and all rotations of each cycle, not by outdegree classes. Every isomorphism sends an actual Hamilton cycle to one of these cycles, so this is another exact canonicalizer.

All1,572 mixed extension candidates agree, with order rows:

- 4: 6 candidates,6 accepted,1 class.
- 5: 14 candidates,14 accepted,6 classes.
- 6: 180 candidates,48 accepted,14 classes.
- 7: 868 candidates,15 accepted,4 classes.
- 8: 504 candidates,0 accepted,0 classes.

Together with the initial triangle, this is26 nontrivial strong cap27 classes. After translating both banks into the independent cycle canonical form, every graph and exact path count matches the source.

The independent literal counter recovers one H13 core and five H27 cores. It checks all75 five-subsets of those H27 cores, finding no13. It also directly verifies the positive induced13-to25 clone example and the hostile induced15-to29 example. The latter confirms that the separate unit-clone obstruction at15 cannot be indiscriminately strengthened to every one-time induced extension.

Independent normal and optimized runs agree. Their artifacts are:

- closed_frontier.py: 5748e1483fa7dedf94fa2bd9a7dc40379049f899079fe04f248ed49f54e0bb0e
- closed-frontier.json and closed-frontier-optimized.json: 84f1415826654afc6916d2719ed6349e7286a6ab6e1796504925d97d3f77ff7f

The checker reads the frozen final-source/VALIDATION.json only for a final bank comparison; it regenerates all graphs, canonical forms and path counts independently. The frozen final-source snapshot is preserved beside it.

## Scope and a checked consequence

Changing to another H13 representative before an append-only construction does not help: all such representatives retain the same strong H13 core. A successful13-to27 operation must cease to preserve the old tournament as an induced graph.

The same closed cap27 inventory also independently recovers the familiar low-count gaps7 and21 after taking SCC products. This is a checked consequence of the finite closed proof, not a claim of priority or a solution of the spectrum beyond27. No stronger uniform binary-induction conclusion is asserted.

No unresolved issue remains in the stated theorem.
