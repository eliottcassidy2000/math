# Hierarchical ternary cells: exact kernels, shared frames and native path flags

23 October 2026. All-size proofs for an explicit family of regular tournaments and supplied named-cell attachment systems. This is not an unrestricted tournament wiring classification.

## 0. Inheritance and the extension

The incoming `odd_cell_completion_20261023/NETWORKS.md`, section 3, already proves the nine-cell example C3[C3,C3,C3], its common-outneighbor recognition, the target kernel C3^3, and a 3m-owner base in a connected network. THM-3134 records the iterated cyclic wreath subgroup and its order; THM-4743 derives full lexicographic automorphism groups via modules. THM-4754 identifies the corresponding ternary tower with a 3-adic first-difference model. These ingredients are inherited and source-pinned, not claimed new.

This note extends the nine-cell attachment to every ternary depth, computes every visibility-depth kernel and the exact joint named-cell group, constructs regular global tournaments with actual Hamiltonian paths, and obtains sharp indexed vertex/arc observation costs on their specified path orbits.

## 1. The internal tournament and its full automorphism group

Let T_k have vertices the words a_1...a_k in {0,1,2}^k. At the first differing coordinate j, orient a->b when b_j-a_j=1 modulo 3. Put T_0 equal to a singleton. This is the k-fold lexicographic C3 tournament, of order n=3^k and outdegree (n-1)/2.

The rooted ternary tree has one C3 rotation at every internal prefix. Its automorphism group is

    W_0=1, W_k=(W_(k-1))^3 semidirect C3,
    |W_k|=3^((3^k-1)/2).

These maps preserve T_k. In fact W_k is its FULL automorphism group. To see this without assuming a module algorithm, let the three top children have size r=3^(k-1). Two vertices in different children have exactly (r-1)/2 common outneighbors: if child A dominates child B, precisely the internal outneighbors of the B vertex are common. Two vertices in the same child have at least r common outneighbors, since both dominate the next whole child. Thus the top partition is intrinsically recovered by the threshold greater than (r-1)/2. An automorphism rotates its three parts cyclically and acts by automorphisms internally. Induction proves the claim. For k=1 the parts are singletons and the same conclusion is immediate.

The prefix rotations are compact coordinates, but their composition is semidirect: upper rotations permute the locations of lower rotations. For k>=2 this is not a cyclic group, nor an abelian product of independent chart residues under composition.

## 2. Attachment visibility and exact one-sided kernels

For 0<=l<=k define

    B_l(a,b)=1 iff the length-l prefixes of a and b agree.

For l=0 this is the uniform all-one block. For 1<=l<=k it is nonuniform, with every row and column sum s=3^(k-l).

Restriction to the depth-l tree gives a surjective homomorphism q_l:W_k->W_l. Its kernel is

    K_(k,l) = product_(3^l prefixes) W_(k-l),
    |K_(k,l)|=3^((3^k-3^l)/2).

A source-fixed automorphism h of the target preserves B_l exactly when q_l(h)=1: its columns must remain in their actual equal-prefix classes. Thus K_(k,l) is the exact one-sided attachment kernel inside the full internal tournament automorphism group.

For l<k, its minimum actual-vertex pointwise base size is exactly

    b(K_(k,l))=3^(k-1).

There are 3^(k-1) bottom triples, each with an independent C3 rotation in the kernel. Every bottom triple must contain an observed owner. Conversely one vertex from each bottom triple fixes every subtree rotation, by induction upward. At l=k the one-sided kernel is trivial, so its base is zero. This is a local statement; it does not erase a shared global frame.

The group order decreases at every increase of visibility depth, while the actual-owner base stays constant until the final depth. Counting group order and counting required marked owners measure different things.

## 3. A regular global tournament with complete informative support

Take m>=3 odd and any regular outer tournament R on m named cells, each containing T_k. For an outer arc i->j define the cross arcs by

    (i,a)->(j,b) iff B_l(a,b)=1,

and reverse them when B_l(a,b)=0. The reverse block is the complementary transpose, so every unordered vertex pair receives exactly one direction. Assume 1<=l<=k.

Each vertex has internal outdegree (n-1)/2, sends s arcs to each outer outneighbor cell and n-s to each outer inneighbor cell. Since both outer degrees equal (m-1)/2, its total outdegree is

    (n-1)/2 + (m-1)(s+n-s)/2 = (mn-1)/2.

Thus the whole mn-vertex tournament is regular. Every intercell block is nonuniform, so the coarse informative support is complete. Nevertheless, for l<k, independent local kernels retain a linear number of owner ambiguities.

## 4. Exact simultaneous chart group and the global-frame boundary

Retain the names of the m cells. Every cell-preserving automorphism restricts internally to W_k. A pair (g_i,g_j) preserves either B_l or its complement exactly when

    q_l(g_i)=q_l(g_j).

The complete support therefore gives the exact cell-preserving automorphism group

    H_(m,k,l) = W_k x_(W_l) W_k x_(W_l) ... x_(W_l) W_k,
    |H_(m,k,l)|=|W_l| |K_(k,l)|^m
      =3^((3^l-1 + m(3^k-3^l))/2).

This is a fiber product with ONE common depth-l frame, not independent choices on each attachment. A spanning tree propagates that common frame, but does not remove the lower subtree kernels.

Its exact actual-owner base is

    b(H_(m,k,l))=m·3^(k-1),  l<k;
    b(H_(m,k,k))=3^(k-1).

For l<k the independent bottom rotations in every cell prove the lower bound. One owner in every bottom triple also fixes every full internal chart W_k, proving sufficiency. At l=k the group is diagonal W_k: each bottom-triple address has a simultaneous rotation across all cells. At least one owner per such address is necessary, and selecting all of them in one cell suffices.

Hence zero one-sided kernel at full visibility does NOT mean zero global owner cost. It means all cells share one remaining chart. This distinction is precisely what an independent-edge alignment count would lose.

These are exact claims for the named-cell system. The full automorphism group of the uncolored global tournament may additionally move or mix cells; no classification of that larger group is asserted. It contains H, so the lower bounds still apply to it. The stated upper bounds rely on retaining the cell carrier.

## 5. Actual Hamiltonian paths and sharp indexed observation costs

The increasing lexicographic list of the ternary words is a directed Hamiltonian path in T_k: at each successor step the first changed digit increases by one, and that digit determines the arc even when later digits reset.

Choose an actual directed Hamiltonian path i_1->...->i_m of the outer tournament. Concatenate the internal lexicographic paths in reverse cell order i_m,...,i_1. At each cell boundary, the last word has all digits 2 and the next first word has all digits 0. Their prefixes differ for l>=1. Because the outer arc points from the next cell to the current one, the complement rule makes the actual boundary arc point forward along the concatenation. This constructs an actual Hamiltonian path P of the global regular tournament.

The group H acts freely on ordered spanning paths: fixing every position fixes every vertex. Consequently its orbit {gP:g in H} has exactly |H| distinct actual Hamiltonian paths.

An indexed vertex-owner observation specifies both a position in P and its actual vertex. Its surviving candidates correspond exactly to the pointwise stabilizer of those owners. Hence the minimum number of such observations identifying P within this orbit is the base size in section 4.

An indexed arc flag specifies a path position and both actual endpoints of that consecutive arc. Put

    t=m·3^(k-1) if l<k, and t=3^(k-1) if l=k.

Then the exact minimum number of indexed arc flags identifying P within the H orbit is

    ceil(t/2).

Necessity: a flag meets at most two independent bottom-triple rotation supports. If fewer than ceil(t/2) flags are observed, one such support is unmarked and its nontrivial rotation preserves every observed flag while changing the path.

Sufficiency: bottom triples occur consecutively in the lexicographic path. Select boundary arcs between the first and second bottom triples, third and fourth, and so on, then one internal arc of the last triple if a triple remains. This observes at least one owner in every required bottom triple. For l<k use the global concatenated list of triples; for l=k use the triples in one cell. Section 4 then kills the full group H.

These are sharp statements on the specified actual Hamiltonian-path orbit, not on every path outside that orbit. The flags retain their clock/position. Dropping the index is a different problem, as Session41's native arc examples already show. Observing the entire ordered path would of course name every vertex and trivially kill all symmetry; no claim of a nontrivial automorphism fixing that complete path is made.

## 6. Relation to lattice constraints and address systems

A cyclic chart action can make each attachment relation a coset encoded by two congruences. Here the full chart group W_k is nonabelian once k>=2. Its exact replacement is a prefix-tree group with a common quotient W_l and explicit lower kernels. Flattening it into unrelated residue choices discards the semidirect action and shared frame.

The construction has compact description length O(n): one ternary rotation coordinate for each internal prefix, plus actual owner/path data. Compact description does not imply a bounded owner budget; the example forces order mn/3 owners at every nonterminal visibility depth.

Fibonacci addresses and these ternary words both supply explicit prefix information, but their admissible languages and arithmetic differ. No identification of Zeckendorf carries with ternary rotations is asserted. Their useful common lesson is to retain actual addresses when a coarse color or prefix leaves a nontrivial kernel.

## 7. Finite checks and limits

The checker constructs internal tournaments and explicit regular three-cell hosts at k=1,2,3, verifies all arcs, common-outneighbor thresholds, regularity and the actual concatenated paths, and checks a generating set of every claimed local/joint group as actual automorphisms. It exhausts W_1 and W_2, not the much larger W_3, to verify small kernel orders, all one-owner-per-bottom-triple bases, and every compatible pair of charts. Exact proofs, not extrapolated census data, establish the general formulas.

No claim here settles arbitrary larger-cell kernels, all metagraph bridges, or unrestricted full wiring. It supplies an all-size hostile family, exact compact joint coordinates, and sharp native indexed-port observation bounds.
