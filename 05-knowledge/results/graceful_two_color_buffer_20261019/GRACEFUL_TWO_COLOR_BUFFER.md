# A width-two repair theorem for the unbounded frozen-tail hostile

## Theorem

Fix integers k≥2 and m≥2k+2, put b=m−1, and freeze the graceful upper-tail forest

    F = { {0,j}: k+1≤j≤m, j≠b } union { {1,m} }.

Its nontrivial component is H={0,1} union ({k+1,…,m} minus {b}); its other components are the singleton vertices2,…,k,b. A completion adds exactly one edge of each difference color1,…,k and forms a tree on0,…,m.

The completion fiber, with F fixed, is connected under moves that replace at most two low-color edges simultaneously, preserving their colors and the tree property. Any two completions are joined by at most2k−1 such moves. Since the previously proved one-color fiber has at least k components, permitting two colors is necessary and sufficient for connectivity of this family.

No assertion is made that2k−1 is the optimal diameter bound. The theorem changes the labeled tree's abstract shape. It does not prove global one-color connectivity GT1, general two-color connectivity for arbitrary graceful fibers, or the Graceful Tree Conjecture.

## Abstract sufficient buffer/contact condition

Let a finite edge-colored graph have a fixed rainbow forest F with r+1 components, and let r unused colors remain. Suppose one component is an isolated vertex b and another is H. Assume:

1. For each free color there is an allowed edge from b into H, and every allowed free-color edge incident with b goes into H.
2. For each free color d, delete b and contract the components of F. The graph of allowed color-d contacts between the remaining components is connected.

Then the rainbow spanning-tree completion fiber of F is connected under at-most-two-color moves, with diameter at most2r−1 whenever it is nonempty. The proof below applies with the attachment-color lemma's edge{z,z+d} replaced by a color-d contact across the non-H component cut. Condition2 guarantees such a contact. Condition1 makes b a leaf in every completion and supplies an attachment edge of whichever color is temporarily stored there.

These are sufficient conditions, not a characterization. They require actual contacts, not merely nonempty color classes. The integer-difference hostile has r=k, and its cut inequality verifies condition2 for every nonempty union of non-H fixed components. Thus its fiber connectivity follows from this abstract buffer theorem.

## 1. The isolated vertex is a color buffer

Every low-color neighbor of b lies in H. In particular b−d∈H for every1≤d≤k, because b−d≥k+1. A completion contains exactly one edge incident with b: no such edge disconnects b, and two would create a cycle through H. Thus b is a leaf carrying a uniquely determined free color d. Deleting b and its incident edge leaves a tree R on all other vertices.

### Attachment-color change lemma

Let e≠d be any other free color. Remove b's color-d edge and the color-e edge. The latter removal splits R into two components. Since F remains intact and connected on H, exactly one component contains H. Denote the other component by S; it is a nonempty subset of{2,…,k}.

Put z=max S. The edge{z,z+d} has color d and crosses the cut of R: z∈S but z+d>max S. It avoids b and remains in the label interval, since z+d≤2k<b≤m. Add{z,z+d} and the edge{b−e,b} of color e. The first reconnects R, the second attaches b to H. This is a graceful completion with the same F, and b now carries color e.

Only colors d,e changed. If a specified set of target edges is protected and the old color-e edge is not among them, the operation preserves that protected set: the only other deleted edge was incident with b.

## 2. Installing target edges without undoing previous ones

Fix a target completion T*. Let its attachment color at b be d*. The target has k−1 other low-color edges. Protect F and every target non-b edge already installed. Take any target non-b edge g of color e that is still absent.

If b does not currently carry e, apply the attachment-color lemma. Its removed color-e edge is not g, since g is absent, and cannot be another protected target edge, since a graceful tree has only one edge of color e. The protected set is therefore preserved. Now b carries e, and R contains every free color except e.

Insert g into R. It creates one fundamental cycle. That cycle contains an old edge f outside the protected target set: otherwise the cycle would be a subset of T*, contradicting that T* is a tree. In particular f is not in F, so its color q is low. Remove f from R and remove b's color-e edge; add g and{b−q,b}. The new R is a tree and b is attached to H. Exactly colors e,q changed. All protected edges remain, and g has been installed.

This takes at most two width-two moves per target non-b edge. After at most k−1 installations, R agrees with T* minus b. The only remaining color is d*, so b has the correct attachment color. If its actual attachment edge differs from T*'s, one single-color replacement puts it at the desired neighbor in H. Total moves are at most2(k−1)+1.

## 3. Scope and inheritance

The construction repairs precisely the hostile family in the frozen-fiber obstruction theorem. The original invariant was the color of b's unique attachment; the width-two operation explicitly changes it while retaining the entire upper tail. The proof then uses that attachment as a temporary color store during fundamental-cycle exchanges.

This uses elementary tree-cut and fundamental-cycle facts. The underlying exchange criterion and rainbow-tree viewpoint are established tools, not claimed as new operations. No historical-priority claim is made for this particular repair theorem without further literature comparison.

Unlike a generic square-sum analogy, the repair has a numerical certificate: z+d≤2k<b ensures a replacement color-d edge crosses the required cut without touching b. The argument depends on this special component geometry and cannot automatically be applied to arbitrary fixed upper tails.

## Exact controls and reproduction

Run `python verify_two_color_buffer.py` and `python -O verify_two_color_buffer.py`. Both receipts agree:891915 checks and84026 individually validated moves. The checker independently enumerates all completions for(k,m)=(2,6),(3,8),(4,10), obtaining6,26,150 states, and transforms every ordered source/target pair:36,676,22500 pairs respectively. Each intermediate state is checked for labels, one edge per difference, tree connectivity, unchanged F, protected target retention, and move width≤2. The maximum lengths produced by this particular algorithm are3,5,7; these are not computed shortest-path diameters.

The general proof covers every k≥2,m≥2k+2; finite controls do not establish that quantifier by extrapolation. Nonemptiness is explicit: attach b by its color-one edge{b−1,b}, and attach each j=2,…,k by{0,j}, together with F.

Source statements: [frozen-fiber obstruction](../graceful_fiber_obstructions_20261019.md) and [graceful exchange criteria](../graceful_square_exchange_20261019.md). The previous one-color obstruction remains valid; this result enlarges the permitted move width while preserving exactly the same frozen upper tail.
