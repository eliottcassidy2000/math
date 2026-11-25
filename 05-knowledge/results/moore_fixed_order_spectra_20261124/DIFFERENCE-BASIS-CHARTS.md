# Difference bases compress the order-11 Moore branch to 1,120 sites

This follows the [coupled action and voltage branches](COUPLED-AXES-AND-VOLTAGE-BRANCHES.md) and the [single-row erasure decoder](MOORE-ROUTE-MACHINE.md). The degree-57 construction is conditional both on existence and on the specified symmetry. Marked row/column labels and the coupled action are decoder metadata.

## Keep rows whose phase differences cover the group

Root on a fixed edge of an order-11 automorphism fixing a pentagon. The axis has one fixed row and five regular 11-cycles, and columns have the analogous action. To recover the moving-row grid it suffices to sample row phases \(D\) in each moving orbit, where \(D-D=\mathbb Z/11\). For any pair of moving rows with phases \(x,y\), choose \(u,v\in D\) with \(v-u=y-x\). Applying the same group element to both rows shifts the pair into the sample. Column labels move by the known same power on the other axis, so the entire matching for those rows can be recovered. This also works for two different rows of one orbit. Use \(D=\{0,1,2,5\}\). In contrast, simple overlap without coordinating phases needs six phases per orbit. Three phases cannot cover every nonzero difference because their ordered differences number at most six. Hence four per orbit is optimal **for covering all moving row-pair orbits by complete-row observations**; other types of decoders might take fewer sites.

The chart keeps \(5\cdot4=20\) rows and 56 columns, omitting the one fixed row. After recovering the moving rows, restore the fixed one by the erasure decoder, then restore the two heads' neighborhoods using the grid equation. Use the monotone orientation on the retained rows, reversing the true graph edges. Rotating or reordering the row interiors avoids the one forbidden graph neighbor at each spine boundary. The marked mate is decoded by undoing tile complementation.

The same method needs respectively \(2072,1848,1960\) sites for order 3 with Petersen core and order 5 with pentagon or Hoffman–Singleton core: choose two phases modulo 3 or three modulo 5, and keep all the fixed rows except one. In order 7 with a fixed star, one axis has no fixed rows; choose that as the row axis and retain three phases \(\{0,1,3\}\) in each of its eight cycles. This is 1344 sites and needs no fixed-row erasure. Bounds below are only for the order-11 chart.

## Explicit count window for the order-11 chart

For both the 1120-site tournament and marked mate,
\[
                          2^{5640}<H<2^{7566}.
\]
In the mate argument count its converse, so all backward graph edges are still present and the additional changes reverse the displayed spine.

**Lower.** Split the 20 rows after ten. The crossing graph is 10-regular on 560 sites per shore, so it has at least \((9^9/10^8)^{560}\) perfect matchings. Apply the [sparse deletion and insertion method](GRID-INSERTIONS-AND-UNMARKED-ROWS.md), with \(k=46\) deleted matched edges, \(q=209\) favorable positions, \(L=142\) eligible edges and \(B=141\) gaps. In the gap estimate, lose at most nine graph neighbors in the same shore, one same-shore mate change, \(k\) unmatched sites, one boundary gap, and at most eleven cross failures. The switching base is \(b=r-38\) and exponent \(9r+1\) for \(r=514\) retained blocks. The diameter-two graph structure supplies the regularity, not any randomness hypothesis. For deleted endpoint density at most four, Brégman bounds matching completions by \((4!)^{46/4}\). A fraction at least \(1-(46+9\cdot46\cdot45/141)/(4\cdot46+1)\) of deletion sets qualify. The path reads back the retained matching, the order, and the insertion slots; the possible deleted-edge completions are the only multiplicity divided out. The lower comparison takes fourth powers and uses integers.

**Upper.** Use the [affine weighted row product](AFFINE-THRESHOLD-AND-ARBITRARY-CONTACTS.md). The 20 rational weights and apex weight are in the replay. A site in row \(i\), local slot \(a\), has ordinary row contributions \(55-a\), one scaled \(w_i/w_j\) for each earlier row, 55 copies of \(w_i/w_j\) for each later row, and its apex arc. In the mate reverse the spine edges and correct the first and last slots of each row. Round the 1121 factors upward before comparing integer products. This is a new conditional chart, not a change to the original 3250-site window.

## Exact Hoffman–Singleton controls using the same decoder idea

In the standard two-sided affine presentation, use edges \(y=x+ij\) between \(A_{i,x}\) and \(B_{j,y}\), and within-cohort steps \(\pm1,\pm2\) on the respective sides. The automorphism sends \(A_{i,x}\mapsto A_{i,x-i}\), \(B_{j,y}\mapsto B_{j+1,y}\), fixes the \(A_{0,x}\) pentagon, and acts on both grid axes by a fixed point plus a 5-cycle. Root at \(A_{0,0}A_{0,1}\), index rows by \(A_{0,4},B_{0,0},\ldots,B_{4,0}\), and columns by \(A_{0,2},B_{0,1},\ldots,B_{4,1}\). Grid labels are \(6r+c\).

- **18 sites, full-row analogue.** Retain three consecutive moving phases (a difference basis modulo 5); the spine is \([6,7,8,9,10,11,12,13,14,15,16,17,19,18,20,21,22,23]\). The replay verifies that their observed graph edges reach exactly all moving-row graph edges under the symmetry. The fixed row follows by erasure.
- **14 sites, edge-orbit sampler.** The spine \([0,1,4,10,9,13,15,20,22,24,25,26,35,30]\) meets all 18 edge orbits, so its observed graph edges and their translates recover the entire labelled grid. This is an explicit sampler, not a minimum-order claim.

Subset endpoint DP gives the exact counts below. The lower bracket counts paths ending at one fixed displayed vertex; the upper bracket is the augmented permanent computed exactly by subset DP. For the mate take its converse throughout.

| Sites | Endpoint | Fixed last-site lower bound | Exact \(H\) | Permanent upper bound |
|---:|---|---:|---:|---:|
| 14 | source | 519,088 | 3,480,779 | 12,304,822 |
| 14 | mate | 752,938 | 5,796,107 | 20,737,356 |
| 18 | source | 596,935,418 | 5,571,834,411 | 24,670,179,110 |
| 18 | mate | 813,914,947 | 7,744,560,549 | 34,510,679,111 |

Both mates beat their sources in these controls. Testing whether that direction persists under more layouts is a concrete conjecture generator; two controls do not determine a degree-57 count or even a universal ordering of its endpoints. These are different encodings from the prior eight-site affine phase atlas.

Run `python3 difference-basis-charts.py` for exact comparisons, orbit coverage, and both Hoffman–Singleton controls.
