# Edge-rooted satellite tournaments and a narrower native window

All degree-57 results below are conditional on a degree-57 Moore graph. The bounds concern the named representations and the specified marked path, not arbitrary tournaments of either order.

| Representation | Sites | Uniform window for source and marked mate |
|---|---:|---:|
| Original native complement-spine chart | 3250 | \(2^{18997}<H<\tfrac9{10}2^{26853}\) |
| Edge-rooted satellite grid, monotone row order | 3136 | \(2^{19663}<H<2^{25627}\) |
| Edge-rooted satellite grid, two cyclic halves | 3136 | \(2^{19663}<H<2^{26795}\) |

The grid also gives a **selected original native spine** for which both 3250-site counts exceed \(2^{19663}\). This last floor is not uniform.

## 1. Remove both root neighborhoods; recover them from the grid

Choose an edge \(uv\) of a Moore graph of degree \(d\) and diameter two. Put \(s=d-1\). The \(s\) heads adjacent to \(u\) other than \(v\) supply row labels; the \(s\) heads adjacent to \(v\) other than \(u\) supply column labels. All heads are pairwise nonadjacent (an edge would make a triangle or a 4-cycle). Every row-head/column-head pair has a unique common neighbor, not among the roots or heads. Thus the other \(s^2\) vertices are the sites of an \(s\)-by-\(s\) grid, with coordinates recording their two incident heads.

Let \(M\) be adjacency on the satellite sites, and let \(R,C\) be the same-row and same-column indicator matrices **including diagonal**. The exact certificate is
\[
                         M^2+M=sI+J-R-C.                 \tag{1}
\]
Indeed a same-row or same-column pair already has a common head, while a pair differing in both coordinates is either adjacent or has exactly one common satellite. Each satellite loses just its two heads, so its satellite degree is \(s-1\). Conversely if a symmetric binary zero-diagonal \(M\) satisfies (1), there are no intra-row or intra-column edges. No site has two neighbors in any one row or column (that pair would then share a satellite despite having the same row or column). Degree \(s-1\) forces one neighbor in every other row and every other column. Add the heads and root edge back as above. Two heads in one part share their root; heads across parts share their grid intersection. Any nonincident satellite/head pair has a unique common satellite in the head's row or column; a satellite and either root have their specified head. The remaining root/head cases and satellite pairs follow from the added edges and (1). This gives a Moore graph, so (1) is a sufficient condition, not just a degree or eigenvalue filter.

For \(s=56\), row pairs carry permutation matchings \(p_{ij}\) on the 56 column labels. Between two different rows \(i,j\), the identity map, the direct matching \(p_{ij}\), and the 54 two-step permutations \(p_{tj}p_{it}\) for other rows \(t\) give **all 56 labels exactly once** at every input. This exact constraint is stronger than testing just a sample of routes. The matrix also has forced spectrum
\[
 55^1,\quad (-1)^{110},\quad 7^{1617},\quad (-8)^{1408}.
\]
On row or column contrasts it acts as \(-1\), on constants as 55, and on the remaining 3025 dimensions use (1) and \(\operatorname{tr}M=0\). In particular its 7-eigenspace projector in satellite coordinates must be
\[
  \frac1{15}\left(M+8I-\frac{R+C}{8}-\frac J{64}\right).
\]
Its diagonal is \(33/64\), and off-diagonal entries are \(-3/320\) within a row or column, \(21/320\) across an edge, and \(-1/960\) across a remaining pair. Thus a partial layout can be rejected if its prescribed Gram submatrices fail positive semidefiniteness or the corresponding rank tests.

### Two decodable tournament orientations on the grid

In the **monotone** template, rows are contiguous, all row-to-row arrows point forward, and each row is transitive. Reverse precisely the edges of \(M\). Order the sites within each row so that the displayed path can cross into the next row without using an edge of \(M\): at most one next-row site is excluded for the previous last site. In the **two cyclic halves** template, keep a forward orientation across the two halves of 28 rows each, but inside either half use a cyclic row tournament: from \(i\) go to \(i+1,\ldots,i+13\) modulo 28 and send antipodal pairs from lower to higher index. Consecutive rows are still aligned with the template. Again reverse precisely \(M\) and choose the spine boundaries by the same rule. The two source tournaments cannot be isomorphic: the monotone maximum score is 3080, whereas the cyclic-halves maximum is 2378.

The marked order and the specified template recover \(M\) by comparison, and the marked mate can first be undone using the same spine. Columns need not be tagged: away from the diagonal the zeros of \(M^2+M\) join rook neighbors; with the rows known, the unique zero in row 0 fixes each site's column. The degree sequence of the monotone source is especially rigid: at slot \(a\in\{0,\ldots,55\}\) of row \(i\) the score is \(3080-54i-a\). Only two tied scores occur at each row boundary; this optionally replaces the full order tag by one choice per tied pair. An unmarked count alone is not a decoder.

### Both grid lower bounds from an exactly regular cut

In either template the first/last 28 rows induce a 28-regular crossing bipartite graph on 1568 sites per side. Schrijver's regular bipartite bound gives \((27^{27}/28^{26})^{1568}\) perfect matchings. A matching pairs each right site to a left site along a backward arc. For consecutive blocks only left-to-right joins are needed. Between distinct blocks there are at most 27 forbidden successors or predecessors and \(1568\cdot27\) forbidden ordered pairs, since the graph edge of a matched block itself never appears between distinct blocks. Apply the [cyclic switching lemma](DEGREE57-FACTORS.md#2-count-compatible-block-orders-by-switching) with an unconstrained anchor. The one additional cross-boundary spine arc forbidden in the mate's converse can be removed last with the same switching budget: none of the four new adjacencies is the very pair being removed. The other constraints are unchanged. Expanding the distinct blocks is injective, and the actual path reads back the matching. Hence both counts are at least
\[
 (27^{27}/28^{26})^{1568}\ 1568!\ (1458/1459)^{42337}>2^{19663}.
\]
On the original native chart put the 114 omitted vertices before the grid. All heads are mutually nonadjacent, so order them starting and ending with left heads from the first 28 rows. Insert \(u\) between right heads and \(v\) between left heads. Choose the ending head in a different row from the first satellite. The whole displayed order is a complement spine. A forward prefix in the source, or the reverse prefix in the mate's converse, ends at a left head with no graph neighbors among the right-half satellites. It joins every counted alternating grid path, proving the selected-spine floor.

### Upper brackets from fixed row geometry

Use a universal bidirectional apex and the diagonal scaling of [the affine-layer proof](AFFINE-THRESHOLD-AND-ARBITRARY-CONTACTS.md). In the monotone grid, weights are constant within each row. A vertex in row \(i\), local slot \(a\), has an exact weighted row sum: \(55-a\) from its own row, \(w_i/w_j\) from each earlier row \(j\), \(55w_i/w_j\) from each later row \(j\), and \(w_i/w_{\rm apex}\) to the apex. Its threshold excess over one is the sum of the positive excesses from earlier rows and the apex. This uses the full per-row matching constraint and **no** spectral relaxation. For the mate's converse, add and subtract the two spine weights explicitly at row starts/ends. The affine threshold correction applies (at most 57 ordinary entries above threshold one even including the extra spine arc and apex). Multiplying upward-rounded row factors gives the exact ceiling \(2^{25627}\).

In the cyclic-halves chart take weights \(1,7\) on the halves and 3 at the apex. The regular internal row geometry and displacement across the cut give exact scaled sums. Within each half, the row tournament degrees are 13 and 14, each 14 times. A row of degree \(d'\) contributes \(54d'+27+55-a\) to the ordinary outgoing sum within its half; across the cut add \(28\cdot55/7\) on the left and \(28\cdot7\) on the right. Add the apex and correction \(5w_i\), then divide each upper factor by \(e\). The mate changes only four row weights: \(-1,1-1/7,7-1,+1\) at the path endpoints and two sites at the half boundary, in displayed order. Exact products on both endpoints give \(2^{26795}\).

Small controls for the monotone chart: the \(C_5\) grid has one site and both counts 1; the Petersen grid has four sites and both counts 5. Hoffman–Singleton has 36 satellites in a \(6\)-by-\(6\) grid. The same matching proof and exact integer Brégman bounds give \(2^{53}<H,H^\ast<2^{107}\) for the displayed control. None of these compressed charts asserts that the smaller grid's count equals the count of the original native tournament.

## 2. Uniform native lower bound: spread the auxiliary deficit

Reuse the [70-auxiliary repaired factor](WEIGHTED-AUGMENTATION.md) with \(r=22\), while starting the crossing fractional matching on *all* genuine edges. For crossing degrees \(d_i,d_j\), give their edge weight \(1/\max(34,d_i,d_j)\) in normalized units (or 22 times that in the factor flow). It is feasible at each site, bounded by \(1/22\) per edge, and its total is at most \(52812.5/34<1555\). The augmentation is then done in the original crossing flow, including the artificial edge if necessary, with each row/column capacity 22 and edge capacity 1 in unnormalized units. Existence of the repaired factor supplies a feasible value \(22\cdot1555\). Augmenting along residual source-to-sink paths never decreases the serving degree of an original vertex; internal residual exchanges conserve it. Thus final normalized deficits \(\delta_i\) satisfy \(\sum\delta_i=W=140\), \(\delta_i\le1-d_i/57\), and \(1-\delta_i\) is at least the prescribed initial weight at that site. Distribute deficits to the opposite auxiliaries as before. The genuine-edge entropy bound by the chord on \([0,1/22]\), fake-edge penalty, completion denominator and block-order bound remain valid.

Put \(D=\sum\delta_i(1-\delta_i)\), \(\mu=115/4\) and \(\sigma=365625/8\). The Moore shifted adjacency \(A+I/2\) has norm exactly \(7.5\) on the centered half-sign vector (its square acts by \(56.25\)), so \(\sum_i(d_i-\mu)^2=\sigma\). For the positive-deviation mass \(V=\sum_{d_i>\mu}(d_i-\mu)^2\), Cauchy–Schwarz gives for the weighted incidence \(S=\sum_i\delta_i d_i\)
\[
 S\ge W\mu-\sqrt{(\sigma-V)(W-D)} .
\]
Drop the positive deviations and note \(\sum\delta_i^2=W-D\). Now expand each initial weight as \(1/57\) plus the telescoping terms \(1/t-1/(t+1)\) for \(34\le t\le56\) on edges with **both** degrees at most \(t\). At most \(\lfloor V/(t+1-\mu)^2\rfloor\) sites have degree greater than \(t\). The loss from excluding each such site is at most \(57+(t+1)(56-t)/57\): incoming low-degree weighted incidences contribute at most 57; outgoing ones at most \(d_i(1-d_i/57)\). Consequently \(D\) must exceed or equal \(S/57\) plus the telescoping coefficients times these surviving nonnegative incidence lower bounds.

At a fixed \(V\), substituting the displayed lower bound on \(S\), the right side has Lipschitz constant at most \(\sqrt{\sigma}/(2\cdot34\sqrt{W-41.49})<1\) in \(0\le D\le41.49\). Thus testing its strict excess over 41.49 at that endpoint rules out all lower \(D\). Between the variance breakpoints where a degree bound changes, the right side is nondecreasing in \(V\). The **6,269** exact breakpoint tests give \(D>41.49\).

Write \(f(v)=(1-v)\log(1-v)-v\log v\). For \(q=70\) the auxiliary correction over the old chord is
\[
 qf(\delta/q)-q\delta f(1/q)
 \ge -\delta\log\delta-\frac{\delta(1-\delta)}{2(q-1)}.
\]
This follows by the second-derivative chord bound on \((1-x)\log(1-x)\) over \([0,1/q]\). Jensen with weights \(\delta/W\) gives total gain at least \(W\log(W/(W-D))-D/138\), an increasing function of \(D\). Using 41.49 and the rational alternating cubic lower bound for the final negative exponential, the exact product comparison gives \(H,H^\ast>2^{18997}\). This is a 70-bit gain over the previous uniform floor.

## 3. Uniform native upper bound: pay for within-bin edges

Reuse the [65-bin affine bound](AFFINE-THRESHOLD-AND-ARBITRARY-CONTACTS.md) but replace its trace-cost matrix \(G\) by \(G+\tau I\), \(\tau=1/4\). In the shifted adjacency compression \(K=B^{\mathsf T}(A+I/2)B\), \(\operatorname{tr}K=65/2+2e_{\rm in}/50\), where \(e_{\rm in}\) is the number of edges internal to these bins. The new controlled quantity is
\[
 \operatorname{tr}(KG)\le65S-455\tau+15\sigma_+-2\tau e_{\rm in}/50 ,
\]
where \(S=\sum G_{ij}/65\) and \(\sigma_+\) is the sum of positive eigenvalues of \(P(G+\tau I)P\) on the mean-zero space. This follows from the same trace norm inequality, using the shifted centered trace \(-S+64\tau\) and then subtracting \(\tau\operatorname{tr}K\). It both gives a uniform ceiling and forces fewer within-bin edges if an observed count approaches the ceiling of this certificate.

The top unshifted eigenvalue is below 104.7. At shifts \(1/20,2/20,\ldots,5/20\), exact LDL gives \(0,3,4,17,25\) additional positive eigenvalues. The right Riemann sum for their total positive parts bounds \(\sigma_+\) by \(2148/20\). Taking \(e_{\rm in}\ge0\) and inserting into the affine AM–GM certificate gives upward-rounded means 885.118 and 885.120 and \(H,H^\ast<\tfrac9{10}2^{26853}\). The two native ends exclude \(E_p\pm1\) for \(2\le p\le9499\) or \(p\ge13427\); the selected-spine floor extends its own small-exponent exclusion through 9832.

Reproduce:

    python3 moore-edge-grid.py
    python3 cyclic-half-grid.py
    python3 spread-native-deficits.py
    node diagonal-shift-upper.js

The exact spectral and factorial comparisons do not require a degree-57 graph to be instantiated.
