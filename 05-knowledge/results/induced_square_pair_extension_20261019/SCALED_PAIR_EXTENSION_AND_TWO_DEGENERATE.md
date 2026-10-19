# Square-scaled pair extension and all finite two-degenerate graphs

2026-10-19. **PROVED**, with exact independent controls. No external priority claim. This note does not alter the earlier frozen two-subdivision packet.

## 1. A universal pair-extension theorem after square scaling

Let b_1<...<b_n be distinct positive integer labels, with H=b_n. Their existing square-sum graph can be arbitrary. Choose a power of four lambda satisfying

    32(H+1) <= lambda < 128(H+1).

Replace every original label by A_i=4lambda b_i. The factor 4lambda is an integer square, so this preserves the complete old adjacency relation, including nonedges.

For every selected unordered pair i<j, set

    d_ij=b_j-b_i,
    r_ij=lambda d_ij-1,
    X_ij=r_ij^2-4lambda b_i.

Then

    X_ij+A_i=r_ij^2,
    X_ij+A_j=(r_ij+2)^2.

All selected pair vertices may be added simultaneously. They are mutually nonadjacent, and each has exactly the two prescribed neighbors among the original vertices.

### Exact incidence

If X_ij+A_k=s^2, take s>=0 and put t=s-r_ij. The new label is 1 modulo four and every A_k is 0 modulo four. Therefore s and r_ij are odd, and t is even. The identity is

    4lambda(b_k-b_i)=2(lambda d_ij-1)t+t^2,

so lambda divides t(t-2). Writing lambda=2^(2a) and t=2u gives

    2^(2a-2) divides u(u-1).

The two factors are coprime. Hence t is congruent to 0 or 2 modulo lambda/2. Also

    |t| = 4lambda|b_k-b_i|/(s+r_ij)
         <= 4lambda(H-1)/(lambda-1) < 5H.

The modulus lambda/2 is at least 16(H+1). Thus only t=0 or t=2 is possible. Substitution gives k=i or k=j, respectively. There is no unintended third incidence.

### Positivity, collisions and new-new nonedges

Every X_ij is congruent to 1 modulo four, so the sum of two different new labels is 2 modulo four and is nonsquare.

Moreover

    X_ij >= (lambda-1)^2-4lambda H > 4lambda H.

The last inequality follows from lambda>=32(H+1), which is stronger than lambda>8H+2. All new labels exceed all original labels.

For a fixed difference d, distinct left endpoints b_i give distinct X_ij. For distinct positive integer differences d<e, the squared-root terms differ by at least

    (2d+1)lambda^2-2lambda >= 3lambda^2-2lambda > 4lambda H.

This dominates any change in the subtracted term 4lambda b_i. Therefore labels belonging to distinct difference classes are also distinct.

Finally,

    X_ij < lambda^2 H^2 < 16384 H^2(H+1)^2 <= 65536 H^4.

The scaled original labels satisfy this same deliberately loose upper bound. This proves the full simultaneous pair-extension theorem.

This theorem preserves the original graph after a specified common square scaling. It does not retain the original numerical labels unchanged.

## 2. All finite two-degenerate graphs occur induced

Every finite two-degenerate graph has an ordering in which each vertex has at most two earlier neighbors: repeatedly delete a vertex of degree at most two, then reverse the deletion order.

Realize the first vertex with label 1. Suppose an earlier induced realization has maximum H. Apply the common square scaling from Section 1.

- With two required earlier neighbors, use their pair vertex X_ij.
- With one required earlier neighbor of base label b_i, use the same formula with the fictitious second endpoint H+1. This endpoint is not included in the graph. Here d=H+1-b_i is between 1 and H. The incidence proof above still forces an original neighbor to have base label b_i or H+1, so exactly the desired one remains.
- With zero required earlier neighbors, use the new label 2. Every old label is 0 modulo four, so no old-new sum is a square. The label 2 is distinct from all scaled old labels.

These steps preserve every existing adjacency and introduce exactly the required earlier adjacencies. Induction proves that every finite two-degenerate graph is induced-realizable in the positive square-sum graph.

The same estimate H_new<=2^16 H^4 covers all three steps. Starting from H_1=1 gives the explicit, very loose height bound

    log_2 H_n <= (16/3)(4^(n-1)-1),       n>=1.

Thus the construction is effective; it does not supply a polynomial height bound for arbitrary two-degenerate graphs. The sharper quartic one-subdivision construction in Section 3 is simultaneous and avoids this iteration.

This proof does not establish a countably infinite two-degenerate extension: repeated rescaling need not stabilize any previously assigned label.

## 3. A sharper simultaneous one-subdivision construction

For a graph G on vertices i=0,...,n-1, n>=1, choose the least power of four lambda>=16n. Thus 16n<=lambda<64n. Set

    A_i=lambda(4i+1),
    X_ij=(lambda(j-i)-1)^2-A_i       for each edge i<j.

The original family is independent: lambda is a square, while the unscaled sums (4i+1)+(4j+1) are 2 modulo four. The new family is independent because every X_ij is 1 modulo four.

The intended incident sums are squares with roots lambda(j-i)-1 and lambda(j-i)+1. For any putative additional incidence with original k, writing the root as lambda(j-i)-1+t gives

    4lambda(k-i)=2(lambda(j-i)-1)t+t^2.

Again t is even and t is congruent to 0 or 2 modulo lambda/2. Rationalizing gives |t|<4n, whereas lambda/2>=8n. Only t=0,2 remain, forcing k=i,j. Positivity and collisions follow from the same gap argument: all new labels exceed H=lambda(4n-3), different gaps have squared-root separation at least 3lambda^2-2lambda>4lambda n, and equal gaps have different subtracted original labels.

Every new label is below lambda^2 n^2<4096n^4; original labels also obey this bound. Consequently the exact one-subdivision of every finite simple graph has an induced square-sum realization with all labels below 4096n^4. With K=n+|E(G)|, this is also O(K^4).

The theorem concerns the subdivision's entire adjacency relation. It does not embed G itself on the original labels, and one-subdivision loses odd-cycle parity and chromatic information.

## 4. Why three prescribed neighbors fail, even after square scaling

This boundary is inherited from [ap_triangle_common_neighbors_20261019.md](../ap_triangle_common_neighbors_20261019.md), Theorem 2 and Section 3.1: three AP anchors admit a rational common neighbor exactly when their gap is a congruent-number class. The gap-one exclusion is also the classical input in the nonsquare-gap corollary of [congruent_square_gap_classes_20261019.md](../congruent_square_gap_classes_20261019.md). The following is its specialization to the present extension question, not a new Fermat theorem or a new AP-anchor classification.

The seed labels 1,2,3 have no common positive square-sum neighbor after any common integer-square scaling. Indeed, if x+S^2, x+2S^2 and x+3S^2 were integer squares, division by S^2 would give positive rational roots u,v,w with

    v^2-u^2=1,       w^2-v^2=1.

The rational numbers w-u, w+u and 2v would be the sides of a right triangle, of area (w^2-u^2)/2=1. Fermat's theorem excludes a rational right triangle of area 1. A verified authoritative source is Keith Conrad, *The Congruent Number Problem*, [Theorem 2 and its descent proof](https://kconrad.math.uconn.edu/articles/congruentnumber.pdf).

Thus the universal square-scaled extension statement is sharp at two prescribed neighbors. This does not prove that any abstract three-degenerate graph fails to embed: the obstruction fixes this particular seed up to square scaling, and a different placement may behave differently.

## Independent controls

`verify_one_subdivision.py` checks complete pairwise adjacency for three separately derived constructions: square scaling, a dyadic/mod-seven construction, and an odd-power-of-three construction. It also checks simultaneous pair extension of arbitrary seeds, preserving all old adjacencies. `verify_two_degenerate.py` checks the incremental construction and the triple-seed obstruction over explicitly declared finite universes. All checks remain active under Python optimization. See the matching JSON receipts and independent review.
