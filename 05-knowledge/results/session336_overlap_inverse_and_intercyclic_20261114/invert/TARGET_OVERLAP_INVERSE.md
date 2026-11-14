# A complete target-first overlap selector for two hubs

For a prescribed positive integer \(N\), there is an exact finite algorithm that either chooses a strong two-hub word with Hamilton count \(N\) and maximizes its overlap, or certifies that no word in the stated class has that target. We order feasible proper overlap frames first by the number \(\beta\) of triply eligible sites, then by the number \(\alpha\) of doubly eligible sites. Availability of any witness at each conjectured \(N\) is a separate problem, not proved here. The search is polynomial in \(N\), not in its bit length.

## An elementary all-order cutoff

Let the normalized hubs be \(x\to y\), the ascending spine have \(m\) vertices, and let \(b_h\) be the first spine vertex beaten by \(h\), \(a_h\) the last spine vertex beating \(h\). Take the formal values \(+\infty,-\infty\) respectively when absent. To each literal backwards bridge assign its interval on the spine; extreme intervals are
\[
 I_x=[b_x,a_x],\quad I_y=[b_y,a_y],\quad I_{xy}=[b_y,a_x].
\]
Intervals with reversed ends are empty. A spine cut separates a nonempty lower prefix from a nonempty higher suffix. If neither one-hub interval crosses it, a hub is either forced high (some higher vertex beats it), forced low (it beats some lower vertex), or flexible (neither). These conditions cannot both force one hub when no corresponding interval crosses. In the absence of a forced-high-to-forced-low hub arc, the two hubs can be assigned to the cut sides with all cross-arcs forward: assign flexible hubs to the side of the forced hub, or both to the same side if neither is forced. This contradicts strongness. The only other case is a forced-high \(x\) and forced-low \(y\) and their given arc, witnessed by \(I_{xy}\). Thus the three intervals cover all spine cuts.

The cross interval is contained in \(I_x\) if \(b_y\ge b_x\), in \(I_y\) if \(a_x\le a_y\), and otherwise contains both. Thus either one of the three spans the whole spine, or the two one-hub intervals together span it, with distinct left and right endpoints except possibly that the end of one is the start of the other. A spanning bridge gives at least \(2^{m-2}\) paths: follow the high endpoint, its connector hub(s), then the low endpoint; every intermediate spine vertex is independent before the high endpoint or after the low endpoint. For a one-hub bridge, insert the unused hub to extend all paths; extension is injective because deleting that hub recovers the path.

If the intervals overlap strictly, their four distinct anchors have the order \(i<j<k<l\), where the first two are the outgoing contacts and last two the incoming contacts. Both hub orders work, giving each as many as
\(2^{j-i-1+l-k-1}3^{k-j-1}\ge2^{m-4}\) paths, since \(i=1,l=m\).
If the intervals only meet at \(s\), use the high-end vertex, the high-interval hub, \(s\), the low-interval hub, and the low-end vertex. All four seams are exactly the bridge incidences: in particular the high-interval hub beats \(s\) and \(s\) beats the low-interval hub. The other \(m-3\) spine vertices choose before the high end or after the low end, giving \(2^{m-3}\) paths. Thus
\[
 H\ge2^{m-3},\qquad m\le\lfloor\log_2N\rfloor+3.
\]
The bound is vacuous for smaller \(m\) and needs no classification of general strong tournaments.

For the unrestricted (not necessarily strong) variant, delete all singleton spine SCCs, retaining the two hubs even if one or both are singleton SCCs. The path count is unchanged by ordinal factorization. Every nontrivial SCC contains a hub; there are at most two. If there are two, each is a strong one-hub tournament, whose single backwards bridge must span its entire spine and hence \(H_i\ge2^{m_i-2}\). The product gives total spine length at most \(\lfloor\log_2N\rfloor+4\). The one-component cases use the proved \(+3\) bound when both hubs lie in the component, or \(+2\) when just one does. A proper overlap has all its four anchors and both hubs in the same SCC: each bridge with its hub forms a directed triangle, and the two spine ranges overlap strictly, which is impossible in different linearly ordered SCCs. The multiply eligible sites in a proper overlap are inside this same SCC. Thus singleton stripping retains the optimum. Two different nontrivial SCCs can also be handled directly by trying divisor pairs \(H_1H_2=N\) and one-hub searches.

## Kernel and exact pruning

For ordered extreme anchors \(i<j<k<l\), let \(A(i,j)\) be the number of matchings assigning one distinct outgoing anchor to each hub, and \(B(k,l)\) the corresponding incoming matchings. These are 2-by-2 permanents on the appropriate incidence bits. For either hub ordering, each pair of such matchings gives a legitimate frame regardless of all other letters. The whole proper-overlap portion of the ordinary count is
\[
 \sum_{i<j<k<l} 2 A(i,j) B(k,l)\,2^{j-i-1+l-k-1}3^{k-j-1}.
\]
The factor 2 is the two hub orders. For each frame the five gaps along the spine have respectively 1,2,3,2,1 eligible runs. All other frames still contribute to the ordinary count; this is not a scalar shortcut that replaces them.

Along the prefix tree of normalized words, induced tournaments keep both hubs. \(H\) is nondecreasing with extension: insert a new vertex into each Hamilton path and retrieve the old path by deletion. So every prefix with \(H>N\) can be pruned. Evaluate \(H\) by the 36-state run transfer; at each exact hit test strongness, then try the ordered quadruples by the matching criterion. This chooses the complete global optimum over all strong realizations, not merely an arithmetic decomposition. The exhaustive candidate bound is \(O(N^2\operatorname{poly}(\log N))\). The replay uses signed long long transfer counts and accepts targets of machine int size; the mathematical algorithm has no such input limit.

## The two targets and the intercyclic restriction

| target | maximum \(\beta\) | max \(\alpha\) at that \(\beta\) | strong witness | ordered anchors, zero-based |
|---|---:|---:|---|---|
| 385 | 2 | 1 | 1030221 | 0,2,5,6 |
| 573 | 3 | 1 | 20302212 | 0,2,6,7 |

These strong witnesses are not subject to no-disjoint-triangles. In fact **any proper overlap**, including zero-length triple gap when \(j,k\) are consecutive, itself gives two vertex-disjoint directed triangles: append the two distinct spine anchors of each backwards bridge to its hub. Triangle-packing number one therefore prohibits all proper overlap. The constrained scan finds exactly 20 and 2 strong normalized labelled words respectively at the indicated targets with packing number one and no single-vertex triangle transversal; both targets occur at spine length 8 alone in this lane. These are word counts with named hubs, not tournament isomorphism-class counts.

Build **g++ -O3 inverse.cpp -o inverse** and run the two targets with and without an **intercyclic** argument to reproduce inverse.out. The cutoff and monotonic pruning make that finite scan a certificate for the target optima without any untested length range.

Also run **python replay_order_masks.py** for a separate subset-DP check of the exact path totals and anchor incidence, independent of the 36-state search recurrence.
