# Edge-rooted satellite tournaments and a narrower native window

All degree-57 results below are conditional on a degree-57 Moore graph. The bounds concern the named representations and the specified marked path, not arbitrary tournaments of either order.

| Representation | Sites | Uniform window for source and marked mate |
|---|---:|---:|
| Original native complement-spine chart | 3250 | \(2^{19025}<H<\tfrac9{10}2^{26853}\) |
| Edge-rooted satellite grid, monotone row order | 3136 | \(2^{19663}<H<2^{25627}\) |
| Edge-rooted satellite grid, two cyclic halves | 3136 | \(2^{19663}<H<2^{26795}\) |

The grid also gives a **selected original native spine** with both 3250-site counts in \((2^{19663},2^{26725})\). This last window is not uniform.

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

There is even a monotone representative whose **unmarked tournament** determines the order under this decoder. Choose uniform independent orders inside the rows. The probability of a legal spine is \((55/56)^{55}\), by exposing rows one at a time. Each site has one graph neighbor in each other row; membership of those neighbors in the four ambiguous edge slots is independent with probability \(1/14\). The chance that any site has more than 21 ambiguous neighbors is at most \(3136\binom{55}{22}/14^{22}<(55/56)^{55}\). Fix a layout that meets both conditions.

Any alternate layout passing the score and grid equation tests differs only inside disjoint four-site boundary slabs, since only the two tied scores per boundary can change assignments. Thus graph edits are inside the slabs, each site gains at most three edges, and every nonambiguous site's incident edges are fixed. If the two decoded graphs differ, their common regular degree ensures an old edge \(uv\) is removed. Among the at least 34 nonambiguous old neighbors \(z\) of \(u\), discard the at most two in possible old/new rows of \(v\) and the at most two in possible old/new columns of \(v\) (at most one per column in each layout). Also discard the at most three that are adjacent to an added neighbor of \(v\), since each such neighbor has at most one old common neighbor with \(u\), and \(z\)'s incidences are fixed. For a remaining \(z\), the old common-neighbor contribution \(u\) to \((M^2+M)_{vz}\) disappears with no compensation, and \(M_{vz}\) has not changed. This contradicts (1), as \(v,z\) differ in both coordinates in both layouts. Thus the decoded graph is fixed; reversing it inside the tournament recovers the transitive template and its unique order. This is a selected structural decoder, not a claim of a unique layout for every arbitrary local order.

### Both grid lower bounds from an exactly regular cut

In either template the first/last 28 rows induce a 28-regular crossing bipartite graph on 1568 sites per side. Schrijver's regular bipartite bound gives \((27^{27}/28^{26})^{1568}\) perfect matchings. A matching pairs each right site to a left site along a backward arc. For consecutive blocks only left-to-right joins are needed. Between distinct blocks there are at most 27 forbidden successors or predecessors and \(1568\cdot27\) forbidden ordered pairs, since the graph edge of a matched block itself never appears between distinct blocks. Apply the [cyclic switching lemma](DEGREE57-FACTORS.md#2-count-compatible-block-orders-by-switching) with an unconstrained anchor. The one additional cross-boundary spine arc forbidden in the mate's converse can be removed last with the same switching budget: none of the four new adjacencies is the very pair being removed. The other constraints are unchanged. Expanding the distinct blocks is injective, and the actual path reads back the matching. Hence both counts are at least
\[
 (27^{27}/28^{26})^{1568}\ 1568!\ (1458/1459)^{42337}>2^{19663}.
\]
On the original native chart put the 114 omitted vertices before the grid. All heads are mutually nonadjacent, so order them starting and ending with left heads from the first 28 rows. Insert \(u\) between right heads and \(v\) between left heads. Choose the ending head in a different row from the first satellite. The whole displayed order is a complement spine. A forward prefix in the source, or the reverse prefix in the mate's converse, ends at a left head with no graph neighbors among the right-half satellites. It joins every counted alternating grid path, proving the selected-spine floor.

### Upper brackets from fixed row geometry

Use a universal bidirectional apex and the diagonal scaling of [the affine-layer proof](AFFINE-THRESHOLD-AND-ARBITRARY-CONTACTS.md). In the monotone grid, weights are constant within each row. A vertex in row \(i\), local slot \(a\), has an exact weighted row sum: \(55-a\) from its own row, \(w_i/w_j\) from each earlier row \(j\), \(55w_i/w_j\) from each later row \(j\), and \(w_i/w_{\rm apex}\) to the apex. Its threshold excess over one is the sum of the positive excesses from earlier rows and the apex. This uses the full per-row matching constraint and **no** spectral relaxation. For the mate's converse, add and subtract the two spine weights explicitly at row starts/ends. The affine threshold correction applies (at most 57 entries above threshold one even including the extra spine arc and apex). Multiplying upward-rounded row factors gives the exact ceiling \(2^{25627}\). For the selected full spine give the 114 prefix sites weight 1. A satellite gains exactly two weighted backward head arcs, and the mate adds one prefix-boundary correction at the first satellite. In the prefix a head has at most one old backward edge to a root and one possible extra backwards spine edge, while the two roots may need 56 more each. Bound every prefix-to-grid row sum by \(56\sum_j1/w_j\), and use the threshold correction 5 since all prefix row weights are at most one. The replay's separate exact product gives the uniform-in-the-graph selected-spine upper bound \(2^{26725}\) for both endpoints.

In the cyclic-halves chart take weights \(1,7\) on the halves and 3 at the apex. The regular internal row geometry and displacement across the cut give exact scaled sums. Within each half, the row tournament degrees are 13 and 14, each 14 times. A row of degree \(d'\) contributes \(54d'+27+55-a\) to the ordinary outgoing sum within its half; across the cut add \(28\cdot55/7\) on the left and \(28\cdot7\) on the right. Add the apex and correction \(5w_i\), then divide each upper factor by \(e\). The mate changes only four row weights: \(-1,1-1/7,7-1,+1\) at the path endpoints and two sites at the half boundary, in displayed order. Exact products on both endpoints give \(2^{26795}\).

Small controls for the monotone chart: the \(C_5\) grid has one site and both counts 1; the Petersen grid has four sites and both counts 5. Hoffman–Singleton has 36 satellites in a \(6\)-by-\(6\) grid. The same matching proof and exact integer Brégman bounds give \(2^{53}<H,H^\ast<2^{107}\) for the displayed control. None of these compressed charts asserts that the smaller grid's count equals the count of the original native tournament.

## 2. Uniform native lower bound: spread the auxiliary deficit

Use the [auxiliary cut technique](WEIGHTED-AUGMENTATION.md) with **78 auxiliaries per side and factor degree 23**. The integer replay checks the flow demands: the shifted spectral cut bound leaves 99 size pairs; the even/odd second-moment test certifies all 99. No artificial repair edge is needed. Exactly 1547 true matching edges remain after the auxiliaries are removed.

Start the crossing fractional matching on all genuine edges with normalized weights \(1/\max(35,d_i,d_j)\), where \(d_i,d_j\) are crossing degrees (multiply by 23 for factor-flow units). The initial weights are feasible at each site and bounded by \(1/23\). The total is at most \(52812.5/35<1547\). Augment in the original crossing flow up to \(23\cdot1547\); the augmented 23-factor guarantees this value. Residual source-to-sink paths never decrease the serving degree of any original site. Thus the final normalized deficits \(\delta_i\) satisfy \(\sum\delta_i=W=156\), \(\delta_i\le1-d_i/57\), and \(1-\delta_i\) is at least the prescribed initial weight at each site.

Put \(D=\sum\delta_i(1-\delta_i)\), \(\mu=115/4\) and \(\sigma=365625/8\). The shifted adjacency \(A+I/2\) has norm exactly \(7.5\) on the centered half-sign vector, since its square acts there by \(56.25\). Thus \(\sum(d_i-\mu)^2=\sigma\). With positive deviation mass \(V=\sum_{d_i>\mu}(d_i-\mu)^2\), the weighted incidence satisfies
\[
 S=\sum_i\delta_i d_i \ \ge\ W\mu-\sqrt{(\sigma-V)(W-D)}.
\]
Drop the positive deviations and apply Cauchy–Schwarz using \(\sum\delta_i^2=W-D\).

Expand each initial weight as \(1/57+\sum_{t=35}^{56}(1/t-1/(t+1))\mathbf1[d_i,d_j\le t]\). Since \(\delta_i(1-\delta_i)\) dominates the initial degree at \(i\) times \(\delta_i\), it remains to lower-bound the weighted incidences at both-low edges. For threshold \(t\), at most \(k=\lfloor V/(t+1-\mu)^2\rfloor\) sites can have degree greater than \(t\). Outgoing incidences from these high sites cost at most \(k(t+1)(56-t)/57\), by \(\delta_i\le1-d_i/57\). The incoming incidence loss is at most
\[
 \min\left(57k,\quad W+\binom{k}{2},\quad\sqrt{(W-D)k(56+k)}\right).
\]
To see the last two bounds, if \(m_v\) counts high neighbors of \(v\), then \(\sum_v\binom{m_v}{2}\le\binom{k}{2}\) because distinct vertices have at most one common neighbor. Charge the first incoming hit to \(\sum\delta_v\), extras to the pair budget; or use \(\sum m_v^2\le57k+2\binom{k}{2}\) and Cauchy–Schwarz. All caps are nondecreasing in \(k\). Subtract the two losses from the lower bound on \(S\) and take the maximum with zero for the surviving incidence.

At fixed \(V\) the resulting lower bound on \(D\) has slope less than one for \(0\le D\le D_0=60.6\). Indeed the \(S\)-terms together have derivative at most \(214/(2\cdot35\cdot9)\). Whenever the square-root incoming cap is active and the surviving term is positive, the cap is at most \(S\le W\mu\) (use the substituted lower bound on \(S\)), so its derivative contributes at most \(W\mu/[2(W-D_0)]\) times the full telescoping weight \(1/35-1/57\). The sum is less than one. At \(D_0\) and between the variance breakpoints where a floor changes, the lower bound is nondecreasing in \(V\). Exactly **5,127** integer/rational breakpoint tests therefore prove \(D>D_0\).

For \(f(v)=(1-v)\log(1-v)-v\log v\) and \(q=78\), the auxiliary entropy correction over the chord is
\[
 qf(\delta/q)-q\delta f(1/q)
 \ge-\delta\log\delta-\frac{\delta(1-\delta)}{2(q-1)} .
\]
This follows by the second-derivative chord bound for \((1-x)\log(1-x)\) on \([0,1/q]\). Jensen under weights \(\delta/W\) gives gain at least \(W\log(W/(W-D))-D/154\), increasing in \(D\). Apply the same weighted permanent, exact \((78!)^2\) completion denominator and ordered-block switch as before, now to the sufficient genuine factor. The base product is
\[
 (22^{22}/23^{21})^{1547}(77^{77}/78^{76})^{156}
 \frac{1547!}{(78!)^2}(1309/1310)^{52929}.
\]
Multiply by the gain using \(D_0\), with a rational alternating cubic lower bound for the negative exponential term. The exact comparison gives \(H,H^\ast>2^{19025}\), a **98-bit** uniform improvement.

## 3. Uniform native upper bound: pay for within-bin edges

Reuse the [65-bin affine bound](AFFINE-THRESHOLD-AND-ARBITRARY-CONTACTS.md) but replace its trace-cost matrix \(G\) by \(G+\tau I\), \(\tau=1/4\). In the shifted adjacency compression \(K=B^{\mathsf T}(A+I/2)B\), \(\operatorname{tr}K=65/2+2e_{\rm in}/50\), where \(e_{\rm in}\) is the number of edges internal to these bins. The new controlled quantity is
\[
 \operatorname{tr}(KG)\le65S-455\tau+15\sigma_+-2\tau e_{\rm in}/50 ,
\]
where \(S=\sum G_{ij}/65\) and \(\sigma_+\) is the sum of positive eigenvalues of \(P(G+\tau I)P\) on the mean-zero space. This follows from the same trace norm inequality, using the shifted centered trace \(-S+64\tau\) and then subtracting \(\tau\operatorname{tr}K\). It both gives a uniform ceiling and forces fewer within-bin edges if an observed count approaches the ceiling of this certificate.

The top unshifted eigenvalue is below 104.7. At shifts \(1/20,2/20,\ldots,5/20\), exact LDL gives \(0,3,4,17,25\) additional positive eigenvalues. The right Riemann sum for their total positive parts bounds \(\sigma_+\) by \(2148/20\). Taking \(e_{\rm in}\ge0\) and inserting into the affine AM–GM certificate gives upward-rounded means 885.118 and 885.120 and \(H,H^\ast<\tfrac9{10}2^{26853}\). The two native ends exclude \(E_p\pm1\) for \(2\le p\le9513\) or \(p\ge13427\); the selected-spine floor extends its own small-exponent exclusion through 9832.

Reproduce:

    python3 moore-edge-grid.py
    python3 cyclic-half-grid.py
    python3 spread-native-deficits.py
    node diagonal-shift-upper.js

The exact spectral and factorial comparisons do not require a degree-57 graph to be instantiated.
