# Native metagraph trees: exact leaf channels and a physical obstruction

Status: PROVED elementary statements below; FINITE-EXACT native order-four counterexample. The global logarithmic existence claim in Hajiaghayi, arXiv:2611.05189v1 (2026), Theorem 1.1 is CITED / CONDITIONAL here, not independently established by this packet. No priority claim.

## Objects, inheritance and losses

The labelled order-n graph has tournaments as vertices and an edge T--K_P(T) for each directed Hamilton path P, where K preserves all path arcs and reverses every other arc. This is not a single-arc-flip graph. For n>=3 each path yields a distinct neighbour, and degree is H(T). Work on one connected component G; no all-order connectivity is assumed. For quotient versions, first specify isomorphism/converse merging, retain parallel-line capacities, and delete loops for network purposes. The counts below concern only the labelled graph.

Closest mechanism: session207 universal metagraph transfer, section2 (actual padded complement neighbours) and section4 (incidence); hostile: replacing actual path-complement adjacency by arc flips; corrected near miss: interpreting a cut guarantee as an electrical one; least-used sidecar: session207's exact automorphism-stabilized neighbour degree. The anchor is tree transport of native metagraph potentials, niche is leaf/bridge exactness, wildcard is electrical energy. Their shared board is adjacency, fundamental cuts, count fibres, scalar potentials and dissipative energy.

The source-to-target map is G to a spanning tree S with capacities c_S(f)=sum of c(e) over original edges crossing the components of S-f. It preserves cut lower bounds and scalar total variation, loses cycles and electrical energy, and needs the original edge/path incidence operator as sidecar. Session207 degree labels alone do not reconstruct G or supply the missing incidence.

## 1. Exact bridge channels, regardless of degree amplification

**Theorem.** Every bridge e of a connected capacitated graph lies in every spanning tree, and its canonical capacity equals c(e). For any original lengths, the tree distance between its endpoints equals its original length. Thus every leaf channel is distortion-free, in every tree and every mixture.

**Proof.** Omitting e disconnects the graph. Deleting e leaves two components with e as their only original crossing edge, so the fundamental-cut sum is c(e). The sole tree path between its adjacent endpoints consists of e. This argument works with parallel edges too, provided bridge has its literal graph-theoretic meaning.

In the labelled tournament metagraph at order n>=3, every transitive tournament has H=1 and therefore is a leaf. By session207's universal transfer applied to a singleton core with two-sided transitive padding, its neighbour has unbounded tribonacci Hamilton count as order increases (equivalently use the older tribonacci-partner theorem THM-4731). Every such edge is nevertheless a bridge and is reproduced exactly by every tree. Consequently arbitrarily large adjacent Hamilton-count ratios do not force any distortion on those edges. This does NOT imply bounded distortion for their whole components. The statement is about actual path-complement neighbours, not padding moves as graph edges.

A useful exact preprocessing consequence: remove all bridges, choose trees independently inside the resulting connected blocks, and restore all bridges. This gives a spanning tree of G. A fundamental cut of an internal tree edge has the same crossing capacity as its cut inside that block: all attached blocks lie wholly on their attachment side. Therefore a common congestion bound on all blocks extends to G, with bridges having congestion1. For distances, the unique sequence of blocks between two vertices makes the same assertion true for expected stretch bounds at least1. This is the bridge-block reduction, not a claim that every graph decomposes into small blocks.

## 2. Scalar total variation, not electrical energy

For any real-valued potential phi on vertices define TV_G(phi)=sum_(uv in E)c(uv)|phi(u)-phi(v)|.

**Theorem.** Every canonical-capacity spanning tree satisfies TV_S(phi)>=TV_G(phi). If a distribution satisfies E[1_(f in S)c_S(f)]<=rho*c(f) for each original edge f, then also E TV_S(phi)<=rho TV_G(phi), simultaneously for every fixed phi.

**Proof.** Expand each tree capacity as its fundamental-cut sum and interchange sums:

TV_S(phi)=sum_(uv in E)c(uv) sum_(f in path_S(u,v)) |difference_f phi|.

The triangle inequality gives the lower bound. The assumed coordinatewise expected-capacity bound gives the upper bound directly. Equivalently integrate cut domination over thresholds of phi (discrete coarea). No distribution needs to be chosen after phi is supplied.

Taking phi=H or phi=log H yields rigorous count-landscape statements; taking phi=cyclic-triangle count links to the native flux/Kirchhoff framework of THM-4739. These are finite graph potentials, not monotone coordinates along every edge. Conditional on the paper's Theorem1.1(ii), rho=O(log N), where N is the size of the selected metagraph component, NOT the order n of its underlying tournaments. Since N<=2^(n choose2), the resulting bound is O(n^2); no O(log n) conclusion follows from that theorem alone.

## 3. A native physical counterexample at order four

A conductance network dissipates E_G(phi)=sum c(uv)(phi(u)-phi(v))^2. Canonical cut domination does not imply E_S(phi)>=E_G(phi), even for the actual tournament metagraph.

The reproducible certificate in RESULT.json specifies all64 labelled tournament vertices (bitmask orientation: lexicographic pairs i<j, set bit means i->j), the63 BFS-tree edges with their positive integer fundamental-cut capacities, and an integer potential. Its only nonzero values are phi(6)=phi(45)=-1 and phi(21)=phi(28)=phi(59)=1. The BFS starts at bitmask0 and visits neighbours in increasing numeric order. The native graph has96 edges, is connected, and has24 leaves. Direct integer sums give

    graph energy =31; canonical-tree energy =23;
    graph total variation =23; canonical-tree total variation =23.

All cuts are dominated by the path-crossing proof, not by an unreported exhaustive 2^64-cut computation. The potential and tree are a fully finite certificate; floating-point eigendecomposition in check.py is used only to discover an integer witness and is not evidence for its strict inequality. An auditor should check the saved certificate directly with integer arithmetic. Order3 is a positive control: two K_(1,3) components, hence every component already is a tree and all energies agree.

Physical meaning: a tree can reproduce or enlarge every domain-wall cost (binary cut / total variation) while underestimating Joule heating at a multilevel voltage. The preserved observable is L1 variation; quadratic energy and effective resistance require separate spectral control. In particular the paper's capacity result does not by itself supply a spectral sparsifier or preserve the metagraph Hodge/cycle currents.

## Scope and open next step

Unconditional: bridge exactness/reduction, coarea transport, and the integer native energy counterexample. Conditional: the logarithmic distribution imported from the new paper. Open: improve the O(n^2) generic bound using native metagraph structure, and identify which cycle blocks cause unavoidable distortion. The session207 tribonacci degree contrast is a hostile control for any proposed lower bound based only on endpoint-degree imbalance.

Reproduce: `OPENBLAS_NUM_THREADS=1 python3 session208/native-tree-physics/check.py`. Python/NumPy is used only for witness search; the numerical certificate is independently integer-checkable. The graph universe is all labelled tournaments of orders3 and4, with all Hamilton paths enumerated, no class or positivity filters.
