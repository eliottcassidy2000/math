# The first increase above a product of disjoint strong tournaments

**PROVED.** Fix disjoint nonempty strong tournaments \(B_1,\ldots,B_k\), allowing singletons, and write \(h_i=H(B_i)\) and \(P=\prod_i h_i\). Retain their internal arcs. Allow every orientation of their cross edges, arbitrarily many additional vertices, and every orientation incident with those vertices. The smallest possible count is \(P\). There is an exact formula for the **next distinct possible count** in terms of orientations involving only two of the fixed blocks at a time.

For two fixed strong tournaments \(A,B\), let \(\delta(A,B)\) be the least \(H(U)-H(A)H(B)\) over tournaments on their union with edges crossing in both directions, preserving both internal orientations. Set it to infinity if the set is empty. With \(K_1\) denoting one vertex, the second distinct count over all extensions is \(P+\Delta\), where
\[
\Delta=\min\left\{2P,\
   \min_i \frac{P}{h_i}\delta(B_i,K_1),\
   \min_{i<j}\frac{P}{h_i h_j}\delta(B_i,B_j)\right\}.
\]
In particular, if \(k\ge2\) and \(a,b\) are the two largest \(h_i\) with multiplicity, the open interval
\[
                         (P,\ P+2P/(ab))
\]
contains no Hamilton-path count of any such extension.

**Proof.** Treat every additional vertex as its own singleton block. The odd-cycle packing identity first gives \(H\ge P\) by permitting only cycles internal to single blocks. If some pair of blocks has edges crossing in both directions, their union is strong. It contains a directed triangle using both blocks. Indeed, a vertex nonuniform on a strong block gives such a triangle along a path from the vertices it beats to the vertices beating it. If every vertex in each block is uniform on the other, all cross edges point in one direction. The crossing triangle adds a positive weight, so \(\delta\ge2\) whenever finite. Retaining the entire pair as one block in the same packing bound gives \(H\ge P+P\delta/(h_i h_j)\), with \(h_j=1\) for a singleton.

If all pairs of blocks are uniform, form the quotient tournament. A transitive quotient gives \(P\) by concatenation. An intransitive quotient contains a directed triangle and therefore has at least three Hamilton paths by the packing identity. Concatenate any Hamilton path of the quotient with one Hamilton path in each block. The resulting paths are distinct, giving at least \(3P\). Conversely, an ordered join of a minimizing mixed pair and all remaining blocks attains its displayed candidate. An ordered join of all blocks and a fresh directed triangle attains \(3P\). The transitive quotient attains \(P\). These constructions and the lower bounds prove the formula. [Grinberg–Stanley odd-cycle identity](https://arxiv.org/abs/2307.05569)

## Mixed edges as a cut on two cyclic orders

**PROVED.** Suppose the fixed strong blocks \(A,B\) have sizes \(n,m\). Whenever their edges cross in both directions, there are at least \(\lfloor nm/2\rfloor\) distinct directed odd cycles using both blocks. Consequently
\[
                       \delta(A,B)\ \ge\ 2\lfloor nm/2\rfloor .
\]
The bound allows a singleton on either side. It is attained by two directed triangles with a single exceptional cross edge, whose union has 17 Hamilton paths against the ordered product 9. For \(n,m\ge3\), equality requires a single exceptional cross edge.

Choose a directed Hamilton cycle in each strong block, using the unique one-site cyclic ordering for a singleton. On the \(N=nm\) cells of \(\mathbb Z_n\times\mathbb Z_m\), color \((i,j)\) by the bit \(X_{ij}=[a_i\to b_j]\). Use a cyclic interval of \(s\) vertices in the first block starting at \(i\), and one of \(t\) vertices in the second starting at \(j\). These intervals form an odd directed cycle precisely when \(s+t\) is odd, \(X_{i+s-1,j}=1\), and \(X_{i,j+t-1}=0\). Each such cycle has a unique cross edge of each direction, which recovers both runs. Singletons use no internal arc.

The admissible shifts \((-(s-1),t-1)\), for \(1\le s\le n,\,1\le t\le m\) and odd \(s+t\), are distinct nonzero elements of the torus. There are \(d=\lfloor N/2\rfloor\), including when both orders are odd. They form a simple directed Cayley graph of indegree and outdegree \(d\). The number of the displayed odd cycles is exactly the number of its arcs from the 1-class to the 0-class. The two cut directions have equal size by regularity. If the smaller class has size \(1\le r\le d\), it sends at least \(rd-r(r-1)=r(d-r+1)\ge d\) arcs outside. This proves the count. Each such cycle contributes at least two paths beyond the internal packing weight \(H(A)H(B)\).

Equality in this cut bound requires either \(r=1\), or \(r=d\) and the smaller class is a bidirected clique. For \(n,m\ge3\) the second case is impossible. If \(n\) is odd, two cells of the same column cannot be joined in both directions, since the two opposite nonzero vertical offsets have opposite parities. Such a clique has size at most \(m<d\). The same argument applies to odd \(m\). If both orders are even the parity of the sum of coordinates gives a bipartition, so a clique has size at most two. This proves the equality condition. The distinction between even and odd moduli enters the proof at the individual offsets. The directed-cut estimate holds in either case. [Two-triangle and singleton replays](factor_gap.py)

An additional weight makes the test sensitive to the particular blocks. For each directed simple path \(\alpha\) in \(A\) and \(\beta\) in \(B\), include the one-vertex paths, and use pairs with odd total vertex count. Give the contact arc
\[
  (\operatorname{end}\alpha,\operatorname{start}\beta)
  \ \longrightarrow\
  (\operatorname{start}\alpha,\operatorname{end}\beta)
\]
weight \(H(A-V(\alpha))H(B-V(\beta))\), with \(H(\varnothing)=1\), and add parallel weights. Twice the weight of the directed \(1\)-to-\(0\) cut is exactly the contribution from packings with **one** mixed odd cycle that crosses twice, and with all other cycles internal. Its two internal paths identify the mixed cycle uniquely. This is a lower bound on the full surplus for every mixed matrix, and is the exact surplus for a singleton color class, since then every mixed cycle must use the same exceptional arc and crosses precisely twice. A minimum-cut certificate and the cost of the cheapest singleton cut therefore bracket \(\delta(A,B)\). This extends the earlier [one-vertex weighted-cut grammar](../chanpak_minimal_core_grammar_20261122/README.md) and [one-reversed-edge two-block formula](../chanpak_arc_fiber_grammar_20261122/README.md).

For a single strong block on \(n\ge3\) sites, the second-value formula gives gain at least \(\min(2H(B),2\lfloor n/2\rfloor)\) even with an arbitrary exterior. For two or more nontrivial blocks a gain of two is also impossible, by the same two-block bound and the untouched factors. This recovers the inherited induced-extension \(+2\) row, whose only possibilities are \(1\to3\) and \(3\to5\). The two-block formula additionally allows changing every cross edge of the original product. [Prior full induced-gap theorem](../chanpak_log_induced_gap_cutoff_20261118/README.md)

## The exponent-61 exclusion with arbitrary cross edges

The archived lower realization at \(E=2^{60}(2^{61}-1)\) supplies disjoint strong factors of heights
\[
       3,\ 3,\ 175,\ 11,\ 13,\ 31,\ 41,\ 61,\ 151,\ 331,\ 1321,\ 2^{61}+1.
\]
Use the two triangles, the seven-vertex carousel, a strong core for each prime height, and the single-back-arc spike on 63 vertices. Their product is \(P=E-1\). Since the two largest heights are \(1321\) and \(2^{61}+1\), every tournament containing these disjoint factors has either count \(E-1\), or count at least
\[
                    E-1+\mathbf{1{,}745{,}528{,}394{,}559{,}950}.
\]
Thus even changing **every** edge between the factors and adjoining arbitrarily many vertices cannot give \(E+1\) while keeping all factor interiors. This quantifies the price of treating multiplication as disjoint physical storage. For the same archived construction at \(p=60\cdot2^s+1,\ s\ge1\), the two largest factors are \(2^p+1\) and \(2^{(p-1)/2}+1\). The analogous forced increase is at least \(2(2^{(p-1)/2}-1)\), by difference of squares. The lower bounds here need not be attained. [Factor construction and component rotations](../chanpak_fermat_solids_seams_20261122/README.md) · [Earlier single-arc cross-component theorem](../chanpak_arc_fiber_grammar_20261122/README.md)

For the specified 63-vertex spike, the torus cut strengthens the numerical bound to
\[
          H=P\qquad\text{or}\qquad H\ge P+\mathbf{493{,}984{,}535{,}660{,}465{,}850}.
\]
Indeed the block of height 1321 needs at least nine vertices since the maximum on eight is 661. Mixing it with the spike costs at least \(2\lfloor63\cdot9/2\rfloor=566\), times every untouched factor. Mixing the spike with another nonsingleton costs at least \(188\) times its untouched product, and the next possible height is at most 331; \(94\cdot1321/331>283\). Mixing it with an added singleton costs at least \(62P/(2^{61}+1)\). A mixed pair avoiding the spike retains its entire height as a factor in the gain, already larger than the stated bound. The uniform cyclic quotient costs \(2P\). This covers all terms in the second-value formula. [THM-4889, exact maxima by order](../../../01-canon/theorems/THM-4889-the-hamiltonian-path-spectrum-by-order-and-by-value-exact-censuses-to-order-11-resolved-by-cyclic-triangles-least-orders-strong-windows-level-floors-and-the-zone-law.md)

More generally, even with just three retained disjoint strong factors of heights \(2^{30}-1,\ 2^{30}+1,\ 2^{61}+1\), a positive increase would already cost at least \(2(2^{30}-1)\). The [125-site exact lower construction](SEAM.md) realizes this three-factor polynomial inside one strong tournament. Equality in the block bound shows that its three algebraic factors cannot be three disjoint subtournaments with those counts. This explains how a factorization can guide a construction while prescribing too much if read as a requirement on its physical parts.
