# Integral reflection splitting and arithmetic holonomy of the colored metagraph

**PROVED general graph statements; FINITE-EXACT order-five specialization.** The graph cycle spaces below are cycles of the **metagraph**, not directed cycles inside a tournament. They are not Hamiltonian-path counts.

## 1. Carrier and inherited pairing

Let G be the literal order-five merged metagraph: vertices are tournament isomorphism classes modulo full converse; its32 edges are complement pairs of64 fixed-Hamiltonian-path tilings. Edges retain multiplicity and loops. Reflection R fixes each merged vertex, fixes the8 blue edges pointwise, and exchanges the24 black edges in12 twin pairs. THM-4731 supplies this pairing; the present computation uses the Session129 independently reconstructed atlas and its owner display-node convention.

Let Q=G/R, keeping one edge for each black pair and each blue edge. Thus Q has20 edges and10 vertices. It is not the literal metagraph. In particular the4 literal black loops become2 loops, at nodes6 and8 in this convention. Orient each nonloop from its smaller node to its larger node. Orient each black loop and then transport that orientation to its R-twin. Orienting two loop twins independently by their smaller tiling masks could instead introduce an artificial sign; that is not the action used here.

The arguments generalize to any finite multigraph with a vertexwise-fixing involution, fixed edges pointwise fixed (never reversed), and all other edges exchanged in pairs. Write P for the number of exchanged pairs, B for the spanning fixed-edge subgraph, and c(-) for connected-component count including isolated vertices.

## 2. Exact integral cycle sequence

Write Z_1(G)=ker(partial:Z^E->Z^V), the integral cycle lattice. The quotient edge map q_* gives a split exact sequence

`0 -> Z^P -> Z_1(G) -> Z_1(Q) -> 0`.

The first map sends a basis vector to e-R(e). Every such difference is a cycle, since its two edges have the same boundary. A chain in the kernel has zero coefficient on every fixed edge and opposite coefficients on each exchanged pair, proving exactness. Choose one representative from each exchanged pair, and retain each fixed edge; lifting quotient edges to these representatives commutes with boundary, giving an integral section. The choice is not canonical.

For order5, G and Q are connected, so

`rank Z_1(G)=32-10+1=23`,
`rank Z_1(Q)=20-10+1=11`,
`P=12`.

The checker constructs both fundamental-cycle bases and independently verifies that quotient lifts together with12 twin differences form a unimodular basis of the full cycle lattice.

## 3. Rational symmetry splits; integral symmetry need not

Define the eigensublattices

`L_+={z:Rz=z}`, `L_-={z:Rz=-z}` inside Z_1(G).

Their intersection is zero. Over Q they span all cycles via `(z+Rz)/2` and `(z-Rz)/2`. Over Z these halves need not be integral.

**Theorem.** There is a canonical finite quotient

`Z_1(G)/(L_+ + L_-) ≅ (Z/2)^r`,

where

`r=beta_1(Q)-beta_1(B)=P+c(Q)-c(B)`.

**Proof.** The kernel of q_* is exactly L_-. The image q_*(L_+) consists precisely of quotient cycles whose coefficients on all nonfixed (black-orbit) edges are even. Necessity is immediate: equal twin coefficients add to an even number. For sufficiency, split each even black coefficient equally over its twins and leave the blue coefficient unchanged; this invariant lift has the same boundary as the quotient cycle.

Reduce integral quotient cycles modulo2 and retain their black coordinates. On the binary cycle space Z_1(Q;F2), the kernel of this black-coordinate restriction consists exactly of cycles supported on B. Graph fundamental cycles form integral bases whose reductions are bases of the binary cycle spaces, so the image dimension is beta_1(Q)-beta_1(B). The induced parity map identifies the displayed finite quotient. Its exponent is2, proving the full elementary-abelian group description rather than merely its order.

For order5 the blue spanning graph has three components: isolated nodes1 and2, and the other eight nodes. Its8 edges give beta_1(B)=1. Therefore

`r=11-1=10`, `[Z_1(G):L_+ + L_-]=2^10=1024`.

The two eigensublattices have ranks11 and12. The exact checker constructs an integral basis of L_+, adds the12 twin differences, and verifies a23-by23 determinant of absolute value1024. This is a genuine index obstruction. It does not identify these1024 cosets with the1024 labelled five-vertex tournaments merely because the counts match.

**Minimal local witness.** Take just one edge of a black loop pair, with coefficient1 on that loop and0 on its twin. It is an integral cycle. Its invariant and anti-invariant projections have coefficients1/2 and±1/2, hence neither is integral. The loop pair exposes one of the ten independent parity obstructions. This explains why dividing black multiplicities by2 can preserve a graph quotient while losing integral cycle information.

## 4. Four genuinely mixed-color cycles

Let K be the black-orbit spanning subgraph of Q. There is an exact sequence

`0 -> Z_1(B) ⊕ Z_1(K) -> Z_1(Q) -> Z_1(Gamma) -> 0`,

where Gamma is the bipartite multigraph whose vertices are the connected components of B and K, and whose edges are the original vertices of Q: a vertexv connects its blue component to its black component. Orient every Gamma edge from its blue component to its black component. For a cycle z=z_B+z_K in Q, assign Gamma edgev coefficient(partial z_B)(v). The sums at each blue component vanish, as do those at each black component because partial z_B=-partial z_K. The kernel is precisely the sum of monochromatic cycle lattices. Conversely any Gamma cycle gives integer demands summing to zero on every monochromatic component. Connected graph incidence realizes any zero-sum integer demand, so choose blue and black chains with opposite prescribed boundaries. Their sum is a lift cycle. This proves exactness integrally. The quotient is free, and the ranks give

`beta_1(Gamma)=|V(Q)|+c(Q)-c(B)-c(K)`.

Here c(B)=3,c(K)=4,c(Q)=1, so the mixed rank is4. The black-orbit graph has12 edges and rank6. Thus

`11 = 1 blue + 6 black + 4 mixed`.

The component-incidence graph makes this concrete. Its two central vertices are joined by five parallel edges, indexed by the five mixed nodes4,5,6,7,8. Its other five edges are leaves corresponding to the three pure-blue and two pure-black nodes. The five parallel edges yield four independent differences. Consequently the four mixed cycle classes cannot be generated by monochromatic cycles alone. The five here are actual mixed tournament nodes, not a deduction from Platonic solids or Fermat-prime counts.

## 5. Vertex arithmetic is flat; PRIMEGAME's control triangle is not

Let f:V(G)->Q_{>0} be any arithmetic node labeling, such as Session129's code7^a17^b55^c or the block-product code. Assign the oriented edge gain

`g(u->v)=f(v)/f(u)`.

For every closed walk its gain product is1, by telescoping. Loops have gain1 and all parallel edges with the same direction have equal gains. Thus vertex codes alone cannot distinguish black twins, detect their antisymmetric cycle directions, or transport nontrivial gain holonomy.

The three PRIMEGAME block products

`P=55/637`, `Q0=539/153`, `R0=585/34`

give the control triangle13->11,17->11,17->13. Its oriented triangular circulation has

`P*R0/Q0=2025/4802 != 1`.

Therefore no assignment of these gains to a closed metagraph walk using only endpoint ratios can preserve that circulation product. This is a precise obstruction to promoting the faithful static three-coordinate codec into an operation-preserving dynamics map.

**Repair with edge data.** For any connected graph, choose a spanning tree. A vertex gauge can make every tree-edge gain1, uniquely after fixing one vertex. The remaining chord gains are exactly independent cycle holonomies. Hence rational gains modulo vertex gauge are determined by one positive rational per integral fundamental cycle. Q requires11 such parameters; literalG requires23. This is the standard spanning-tree proof of graph gain classification, supplied explicitly here without a novelty claim. A single nontrivial PRIMEGAME holonomy needs at least one edge parameter, not an extra arbitrary node label.

The cycle lattice's integral ± obstruction explains an additional distinction: a rational eigenspace decomposition is always available, but integral cycle decompositions retain ten binary compatibility conditions. The arithmetic square test for binary phase decorations in the companion packet is another consequence of an involution acting on a different integral lattice. These are related mechanisms, not an asserted isomorphism between the two carriers.

## 6. Reproduction and source scope

Run `python3 check_lattice.py` and `python3 -O check_lattice.py`. Output includes the full20-edge colored atlas, quotient fundamental-cycle basis, even-lattice coordinate basis, ranks, determinant, component-incidence edges, and holonomy. Failure gates remain active under optimization; `--test-failing-gate` must exit1. The graph data are pinned, not silently downloaded at runtime. The independent session auditor reconstructs the atlas from all labelled tournaments rather than importing these constants.

Inherited sources: Session129 tournament census at commit5431278adcd6269b5eca7756acc7a842938cd8cf; THM-4731 black-twin pairing; the earlier `merged_metagraph_black_quotient_bridges_20261020/README.md` contains a binary black-doubling sequence. The general integral eigensublattice statement and its explicit order5 index are proved above, not inferred from that binary sequence. No universal Hamiltonian-count classification follows.
