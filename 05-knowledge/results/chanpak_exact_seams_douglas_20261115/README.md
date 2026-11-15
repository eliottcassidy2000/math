# Exact seam rigidity and a near-floor gap for uniquely Hamiltonian tournaments

Let H count directed Hamiltonian paths, C count directed Hamiltonian cycles up to cyclic rotation, and f(n) be the Busch minimum for strong tournaments.

**Proved.**
1. Any strong Busch split with local surplus below n−2 has an exactly Ferrers cross matrix in the cycle orders of its factors. This strengthens the previous [bounded-edit theorem](../chanpak_ferrers_seams_20261115/README.md).
2. A tournament on n≥5 vertices with C=1 has either H=f(n) or H≥(47/45)f(n). A more precise bound when the sorted spine is a floor tournament is given below.
3. No strong tournament of order n≥8 has H=f(n)+2. If n≥9 and H=f(n)+4, it must have exactly two Hamiltonian cycles. In every Busch split the two factors are floor chains, the seam is Ferrers, and the second cycle has an invalid projection to a factor of order 4 or 5.

This excludes counts at a **fixed order of a strong tournament**, not from the unrestricted Chan–Pak spectrum: a count may be realized at smaller order or by several strong components. The two-cycle case and spectrum coverage remain open.

## 1. An adjacent defect costs n−2

Use X,Y,F,G from the [two-sided shuffle bound](../chanpak_ferrers_seams_20261115/README.md), and D[i,j]=1 for a_i→b_j in two cycles started at their shared root v. Write p=|A|−1, q=|B|−1, n=p+q+1. In the prefix calculation X never decreases down columns and Y never decreases along rows. At an interior cell D=0 forces X≥1, and D=1 forces Y≥1.

A violation has an adjacent horizontal 1,0 pair or vertical 0,1 pair. Write (i,j) for its 1, (u,w) for its 0, so u≤i, j≤w and (i−u)+(w−j)=1. At (i,w) both X and Y are positive. All cells below it in the same column, and all to its right in the same row, have F≥2. For a downward step, X persists if it is already at least two. Otherwise if X,Y are both positive above, bit 0 adds the positive Y to X and bit 1 regenerates Y from the positive total at the cell to the left, leaving X positive. The horizontal argument is dual. The southeast corner and its two rays contain p−i+q−w+1 distinct prefix bad cells. The converse calculation, with sides swapped and orders reversed, starts at the northwest corner (u,j) and gives u+j−1 suffix bad cells.

Both charges count towards L=H(T)−H(A)H(B), even if charged to the same pair of cycle rotations: each fiber contributes FG−1≥(F−1)+(G−1). Their sum is p+q−1=n−2. Thus L<n−2 forces D exactly Ferrers. The abstract shuffle bound is sharp: put the sole 0 of an all-1 matrix at (p,2). The prefix bad cells are exactly (p,j) for 2≤j≤q, each of multiplicity 2; in the converse the sole 0 lies at (q−1,1), and the bad cells are its last row with all positive columns, each of multiplicity 2. The two sets of fibers are disjoint, so the rotation surplus is q−1+p. This sharpness does not assert the absence of nonrotation or broken contributions to the tournament count.

## 2. Sort a Douglas tournament without losing cycle packings

Use [Douglas's unique-cycle characterization, Proposition](https://repository.lib.ncsu.edu/items/c636a073-5b23-48fe-bd24-ebd633b597b6), also [journal version](https://doi.org/10.1112/plms/s3-21.4.716). A unique Hamiltonian cycle has a presentation x,c_1,...,c_k,y,d_1,...,d_m,x with outdegree(x)=indegree(y)=1; the c-list has forward consecutive edges and backward nonconsecutive edges, and the d-list is transitive in its displayed order. The cross condition we use is: if i<j and r≤s, then c_i→d_r implies c_j→d_s.

Take the ascending spine y,c_k,...,c_1,x. Its only backward arcs are between consecutive marks. Each d-vertex belongs to a definite gap between consecutive spine marks: spine marks up to that gap beat it, and it beats the remainder. This follows from the cross condition with r=s and the two extreme degrees. The sequence of gap numbers in the d-list has no inversion of size ≥2. Otherwise two intermediate spine marks give c_i→d_r and d_s→c_j with i<j and r<s, contradicting the cross condition.

Sort the d-list by gap, retaining order inside gaps, and orient the d-vertices transitive according to this new order. The resulting tournament S has a natural total order (insert the d-vertices in their gaps), with all arcs forward except the consecutive spine arcs just described. It is strong: reverse down the whole spine, then pass forward through all the d-vertices and close, or close directly if the d-list is empty.

A directed cycle of S is specified uniquely by a spine interval of at least two marks, traversed backwards, and a subset of the unmarked vertices inside it, traversed in increasing order on the return. To verify this, each cut between marks has only one arc crossing backwards. Every cut between the least and greatest visited marks must use that arc. Unmarked vertices outside that interval cannot participate because their adjacent cuts have no backward crossing usable by the cycle. The return is therefore entirely unmarked between the endpoints. Degenerate selections of two vertices alone are not cycles and are omitted.

Embed such cycles into T by retaining the interval and subset and reordering the return in the original transitive d-list. The first and last spine marks have the required arcs to every selected unmarked vertex. This is an injection on cycles, preserves lengths and vertex supports, and thus embeds all vertex-disjoint odd-cycle packings. By [Grinberg–Stanley's Hamiltonian-path specialization](https://arxiv.org/abs/2307.05569), each packing of j nontrivial odd directed cycles contributes 2^j, with the empty packing contributing 1. Hence H(T)≥H(S).

Any inversion in the d-list gives a new triangle Δ on the inverted unmarked vertices and their separating spine mark. It has only one spine mark, so it is not in the image. Augment any packing from S−V(Δ) by Δ after embedding to obtain
```
H(T) >= H(S) + 2 H(S−V(Delta)).
```
If there are no inversions, T=S under this ordering. The cycle-packing argument is important: sorting is not an arc-by-arc monotonicity statement about Hamiltonian paths.

## 3. The gap word of the sorted spine has a three-coordinate exact count

Let a_1,...,a_t be the counts of unmarked vertices in successive gaps. An interval of r spine edges containing b unmarked vertices contributes weight 2^b if b>0 to the packing sum: the parity condition on the return subset selects half of them and each selected cycle carries weight 2. If b=0 the weight is 2 for r positive even, and 0 for r odd. Intervals are compatible exactly when their marked endpoints and interior marks are disjoint.

At the current rightmost mark let h be the packing total, s the sum of preceding packing totals over possible interval starts (weighted by 2 to the number of unmarked vertices up to the current mark), and z the analogous alternating sum over starts connected to the current mark by only zero gaps, including that mark. Start at (h,s,z)=(1,1,1). When appending gap a, put t=2^a and e=[a=0]. The update derived from the preceding interval weights is
```
(h,s,z) -> (h+t*s-e*z, h+t*s, h-e*z).
```
In particular after a positive gap s=h and 0≤z≤h/3; throughout 0≤z≤h≤s, as follows directly from the update.

One leading zero can be removed by moving the first spine mark into the next gap and increasing that gap by 1. The next spine mark becomes the natural top, and the original first mark is now its first unmarked successor. This does not change the tournament up to reordering. Similarly remove one trailing zero. For n≥3 this normalizes the two ends to positive gaps (the only single-gap all-zero graph has order 2).

Write F(d)=f(d+1) for d≥2. It satisfies F(d+3)=5F(d), and F(d_1)F(d_2)≥F(d_1+d_2) for d_1,d_2≥2 (by the floor formula, or by strong overlap of two floor chains). An all-positive normalized word has h=product(1+2^a). By the [Busch stability estimate](../chanpak_busch_stability_20261115/README.md), if this differs from F(sum(a+1)), it is at least (27/25)F(sum(a+1)).

For a normalized word with zeros, group its first zero run with the positive gap a immediately preceding it and the positive gap b immediately following it. Remaining gaps are either single positive gaps or a run of k zeros followed by a positive gap b. Before the first grouped pair there are only positive gaps, so the normalized state after a is (1,1,1/(1+2^a)) with a scalar factor 1+2^a pulled out. For the subsequent groups the normalized state is (1,1,z) with 0≤z≤1/3.

To check the comparisons, define V_k(t,z)=h_k+t*s_k after k zero updates from (1,1,z). It obeys the Tribonacci recurrence in k because the zero matrix has characteristic polynomial λ³−λ²−λ−1. Thus:
- V_k(t,z) is nonincreasing in z for t≥0, by its first three values (1+t), (2+2t−z), (3+4t−t*z).
- V_(k+3)≥5 V_k for t≥0, 0≤z≤1/3. For k=0,1,2 the respective differences are 1+2t−(1+t)z, 1+3t+(3−2t)z, and 5+4t+(t−3)z; all are nonnegative. The rest follows from the same recurrence.
- V_k(8t,z)≥5V_k(t,z) for t≥2 because h_k≤s_k.
- For the paired multiplier, write V_k(t,z)=A−Bz with A≥0,B≥0. Increasing a by 3 increases (1+2^a)A−B by at least a factor 5, since 1+8s≥5(1+s) for s=2^a≥2.

It is now enough to check k,a,b∈{1,2,3}. For a subsequent zero group the minimum of V_k(2^b,1/3)/F(k+b+1) is 29/27; for the first paired group the minimum of (1+2^a)V_k(2^b,1/(1+2^a))/F(k+a+b+2) is 49/45, at a=b=2,k=1. Each remaining positive gap costs at least its floor. Combining group costs by the product inequality gives H(S)≥(49/45)f(n) whenever an internal zero remains. This estimate is sharp at word (2,0,2), and also after appending any number of gaps 2. Overall a sorted spine below the ratio 27/25 must be a floor chain.

## 4. Consequences for the strong spectrum at fixed order

If S is exact floor, its normalized word is positive, all link orders are 3,4,5, with at most one order-5 link and no order-3 link if one of order 5 occurs. If a d-list inversion exists, its triangle touches only the two adjacent positive gaps participating in the inversion. Normalizing empty *end* gaps does not move either participating unmarked vertex or their separating mark out of the corresponding adjoining links. Hence all exterior links remain intact after deleting the triangle. Their counts multiply on the prefix and suffix chains; if both chains are present they have disjoint vertex sets. Their packings alone give H(S−V(Δ))≥the product of all exterior link counts. If h_1,h_2 are the two touched link counts, then
```
H(T)−f(n) >= 2 f(n)/(h_1 h_2) >= 2 f(n)/45.
```
If S is nonfloor the stronger 27/25 ratio already holds, so the claimed all-order 47/45 ratio follows for n≥5. More directly, a floor chain of order n≥9 has at least three links, since two have total order at most 8 under the budget. Thus any inversion costs at least 6. At n=8 the same packing argument guarantees only an excess of at least 2 if neither exterior chain remains. The boundary is checked exhaustively over Douglas's parameters: for r spine edges and 7−r unmarked vertices, enumerate all words of gap indices in {0,...,r−1} with no inversion of size at least 2, for r=2,...,7. This is 144 parameter words, not necessarily distinct isomorphism classes. The minimum H above 45 is 49. Counting each with the literal subset path recurrence gives the n=8 exclusion of 47. At n=7 the positive gap word (1,3) gives H=27=f(7)+2, so the cutoff n≥8 is sharp.

For n≥9, the [sharp cycle inequality](../chanpak_near_floor_crossing_20261115/GLOBAL.md) forces C=1 at excess 2, which is impossible by the preceding argument. At excess 4 it forces C≤2, while the same argument excludes C=1. At n≥9 one has 4<2f(n)/25, so both factors of any Busch split are exact floor. The cross matrix is Ferrers by §1, hence its one interleaving is the only good cycle over the unique factor cycles. The other cycle has an invalid projection. A factor of size m≥6 would force at least m−1≥5 broken rotations. An invalid projection on a triangle gives all n rotations broken. Thus the second cycle has an invalid projection on a factor of order 4 or 5. In particular where the floor budget forbids a link of order 5, it must be order 4.

[Replay](verify.py) uses literal subset Hamiltonian-path DP for the sorted-spine and Douglas checks, integer shuffle tables for the seam bound, and exact rationals for the finite transfer comparisons.
