# Reflection pairs, colour tree polynomials, and parity obstructions to integer currents

Status: elementary all-order deductions for the precisely specified network, with finite-exact controls at orders 3 through 7. This is not a classification of the metagraph, a proof of black connectivity, or a claim of literature novelty. No reserved theorem is used.

## 1. Inheritance and the exact network

The source is the fixed-Hamiltonian-path tiling cube X={0,1}^m, m=binom(n-1,2), n>=3. C flips every free arc; R reflects positions without flipping bits. The owner map pi merges tournament isomorphism classes with converse. Thus pi R=pi, whereas C need not descend to a map on owners. An actual line is an unordered pair {t,Ct}. Loops and parallel lines are retained. A line is blue precisely when Rt=t; otherwise it is black.

Every R-orbit of black lines has two elements with the same unordered owner endpoints. Indeed R commutes with C. It could fix a black line only by Rt=Ct, but R has a fixed free coordinate, where this equality is impossible. Let B be the literal blue multigraph and Q the multigraph of R-pairs of black lines. The literal full network is therefore

    G = B + 2Q.

This does not say that a black support edge is a bridge or not a bridge: individual literal black edges cannot be bridges, but an edge of Q can be a bridge. Loops contribute two stubs to mass and zero to the electrical Laplacian.

For an ordinary tournament class T the tiling mass is H(T)/|Aut(T)|, an odd integer. A merged SC node has that mass; an NS node has twice that mass. At SC nodes the number of blue stubs is odd; NS nodes have no blue stubs. Consequently deg_Q(u) is odd at NS nodes. These are inherited results, not new claims here.

Closest proved sources:
- 01-canon/theorems/THM-643-blue-black-line-parity-structure.md, including its October 20 correction.
- 01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md, section 4.
- 05-knowledge/results/merged_metagraph_black_quotient_bridges_20261020/README.md, sections 1–3: paired-black cut charges and constructive Q bridges.
- 05-knowledge/results/merged_metagraph_endpoint_transport_20261020/endpoint-transport-and-complement-leakage.md, section 5: compression is not generally lumped dynamics.

The historical 07-reflections/the-smith-diagram-of-the-metagraph-opus-S307.md uses the unmerged d=1 wiggly-flip network. Its conductances and tree counts are not those of the merged d=m line network here; none are imported.

Canonical hostile: the order-five regular tournament has three pure-blue tilings and three distinct blue neighbors. Pure-blue does not mean leaf. Even-order pure-blue nodes, by contrast, are singleton blue leaves by the current purity theorem. Black-core connectivity and an all-order classification of Q bridges remain open. Full and blue-SC connectivity are also explicitly NOT PROVED all-order in THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md; its finite controls must not be used as a global premise. THM-4739 and THM-4740 were RESERVED empty stubs at the pinned input and are not dependencies.

## 2. A componentwise colour spanning-tree theorem

Fix one connected component Gamma of G. All statements below are componentwise, so they assume no unproved global connectivity. Write v=|V(Gamma)|. Let c_B and c_Q be the numbers of components of B and Q restricted to this vertex set, retaining isolated vertices. All edges have their literal integer multiplicities. Define

    Psi(z) = det((L_B + z L_Q)_red).

The same vertex is deleted in both Laplacians. Equivalently, Psi is the spanning-tree polynomial in which every paired-black edge has weight z and every literal blue edge has weight 1. It is a new symbol defined here, not an identification with any separately defined D(z).

### Theorem 1

Set a=c_B-1 and r=v+1-c_B-c_Q. Then r>=0 and

    Psi(z) = z^a sum_(j=0)^r p_j z^j,       p_j positive integers.

All r roots of the shifted polynomial are real and strictly negative. In particular every possible black-edge count between c_B-1 and v-c_Q occurs in a spanning tree: there are no gaps. The sequence p_j/binom(r,j) is log-concave (when r<=1 this is vacuous).

The endpoint coefficients are

    p_0 = product_i tau(B_i) * tau(Q/B),
    p_r = product_j tau(Q_j) * tau(B/Q).

Here Q/B contracts each blue component to a vertex, B/Q contracts each paired-black component, parallel edges retain multiplicities, and loops are discarded. A one-vertex graph has tree count 1. The actual literal network has

    tau(Gamma)=Psi(2).

In particular 2^(c_B-1) divides tau(Gamma).

### Proof

The ordinary matrix-tree identity follows directly by writing a reduced incidence matrix A and applying Cauchy–Binet to det(A diag(w) A^T): the nonzero maximal minors are precisely spanning trees, with squared determinant 1. This gives the stated weighted tree interpretation with no orientation dependence.

The matrix L=(L_B+L_Q)_red is positive definite, since B+Q has the same connected support as G. Put S=L^(-1/2)(L_Q)_red L^(-1/2). Both S and I-S are positive semidefinite, so its eigenvalues lie in [0,1]. Therefore

    Psi(z)=det(L) product_i (1-lambda_i+lambda_i z).

The multiplicity of lambda=0 is c_Q-1: a reduced Laplacian of a graph with c_Q components has nullity c_Q-1. The multiplicity of lambda=1 is c_B-1 by the same argument for I-S. The other r eigenvalues lie strictly between 0 and 1. Factoring out z^(c_B-1) leaves r positive linear factors, proving positivity, the full interval of coefficients, and strict negative real roots. The usual Newton inequalities can also be obtained by expanding this product of positive linear factors.

A spanning tree using exactly c_B-1 paired-black edges must contain a spanning tree inside each blue component and a tree connecting the blue components. Conversely any such choices give a tree. This proves p_0; exchange the colours for p_r. Substituting z=2 restores the two separate literal choices in every selected black pair. The divisibility follows. QED.

This is an application of standard incidence, matrix-tree, and real-rootedness mechanisms to the actual reflection-paired carrier. Its content for this problem is the forced substitution z=2, the correct treatment of owner multiplicities, and the resulting parity obstructions. No generic Laplacian fact is presented as a new theorem of tournament classification.

## 3. The exact mod-two obstruction and its electrical consequence

For a connected integer network define its critical group by

    Jac(Gamma)=Div^0(V)/L_G Z^V.

It has order tau(Gamma). Its 2-rank is dim_F2(Jac/2 Jac), equivalently the number of even nonzero Smith invariant factors. Define the blue bicycle dimension by

    b_B = dim_F2(ker A_B intersect im A_B^T),

where A_B is the full blue incidence matrix over F2. This definition includes every literal blue edge; loops give zero incidence columns, and parallel edges remain distinct columns.

### Theorem 2

For the actual literal network G=B+2Q on Gamma,

    dim_F2(Jac(Gamma)/2 Jac(Gamma)) = c_B-1+b_B.

Consequently

    2^(c_B-1+b_B) divides tau(Gamma).

There is a canonical surjection

    Jac(Gamma) -> {(epsilon_i) in F2^(c_B): sum_i epsilon_i=0},

sending a degree-zero divisor to its total chip count modulo 2 on each blue component. Its dimension is c_B-1. The extra b_B mod-two invariants are precisely the blue bicycle contribution, not additional black connectivity components.

If Gamma has r_NS non-self-converse nodes and at least one self-converse node, its 2-rank is at least r_NS. If Gamma consists entirely of NS nodes, the guaranteed rank is r_NS-1. These statements must not be silently interchanged.

For every NS node u and every different node w in the same component, the class of e_u-e_w has even order in Jac(Gamma). Equivalently, the grounded rational voltage vector solving

    L_G phi=e_u-e_w

has an even common denominator in lowest terms. The claim concerns the complete voltage vector, not just its single effective-resistance difference.

### Proof

Because literal black conductances are even, L_G reduces modulo 2 to L_B. The 2-rank of the critical group is therefore the nullity of a reduced L_B over F2. The full matrix L_B=A_B A_B^T has nullity c_B+b_B: map its kernel under A_B^T into the bicycle space; the image is exactly that space, and the kernel is ker A_B^T, of dimension c_B. Reduction at one vertex lowers the nullity by one. To see the last assertion explicitly, a function in the full kernel can be normalized to vanish at the root; the other rows are exactly the reduced equations, and the root equation follows from the zero sum of all rows. This proves the rank formula and divisibility.

Every integer Laplacian firing changes the total on a blue component by an even number: internal blue contributions cancel, and all crossing edges are black with even conductance. Thus the displayed parity map descends to Jac. Given any even-sum component parity vector, choose one vertex in each component and use differences from a fixed reference vertex to realize it. Hence the map is onto.

An NS node is a singleton blue component. Its individual chip parity is therefore invariant under every integer firing. The divisor e_u-e_w has odd parity at u, and any odd multiple still does, so no odd multiple is principal. The order is even. If d is the least positive integer making the rooted solution d phi integral, then d(e_u-e_w) is principal, and conversely any principal multiple gives an integral rooted solution; thus d equals the divisor order and is even. QED.

The blue bicycle space is genuine extra structure. Replacing it by component count alone already fails in the actual order-five network. Conversely the valuation of the tree count can exceed the 2-rank; order six is an actual hostile below.

## 4. Exact bounded controls

The supplied probe rebuilds every literal line from the atlas tiling-to-owner map, rather than treating the stored edge matrix as an oracle. For orders 3–5 it independently reconstructs every ordinary fibre by direct permutations and checks the complete owner classification. For orders 6–7 the frozen owner maps are inherited finite-exact inputs, not independently reclassified here. Source hashes are recorded in probe.json.

| n | merged nodes | NS nodes | blue components | blue bicycle dimension | critical-group 2-rank |
|---:|---:|---:|---:|---:|---:|
|3|2|0|1|0|0|
|4|3|1|2|0|1|
|5|10|2|3|1|3|
|6|34|22|23|3|25|
|7|272|184|185|11|195|

The full literal graph happens to be connected in these finite universes. That observation is not used as an all-order premise.

At order five,

    Psi(z)=20z^2+197z^3+599z^4+603z^5+165z^6,
    tau(G)=41096=2^3*5137,
    nonunit Smith invariants=(2,2,10274).

This is an especially useful hostile: the lowest-degree coefficient times 2^2 is 80, of valuation 4, but the entire evaluation Psi(2) has valuation 3. Thus the lowest-degree term alone does not determine the exact valuation. The minimum exponent 2 also misses the additional blue bicycle.

At order six,

    tau(G)=9672890651856613234344767407709495492608,
    v_2(tau(G))=26,              2-rank=25.

Thus the lower bound from the exact 2-rank is not always an equality. Neither the minimum exponent nor the 2-rank uniquely determines the full 2-primary group.

The probe computes integer Smith forms only through order six and mod-two ranks through order seven. It does not extrapolate a formula for the observed bicycle dimensions 0,0,1,3,11.

## 5. An actual Kron-reduction hostile: resistance forgets the parity witness

The distinction between the voltage vector and the terminal resistance is realized in the actual order-five metagraph, not merely in an abstract network. Use atlas node 2, represented by tiling mask 2, as the source, and node 0, the transitive class, as the grounded sink. The source is NS and has ordinary representatives of mass 5 each, so its merged mass is 10. Exact unit-current solution gives

    R_eff(2,0)=6092/5137,
    phi(1)=11327/10274,
    common denominator of the full rooted voltage vector=10274.

Thus terminal resistance has odd denominator, while the full voltage vector has the forced even denominator. Eliminating all internal nodes by exact Kron reduction gives the two-terminal Laplacian

    (5137/6092) [[1,-1],[-1,1]].

This preserves the exact Dirichlet-to-Neumann response and effective resistance, but its two terminal potentials alone no longer display the integral-current obstruction. The harmonic extension sidecar restores it: for an arbitrary source-terminal potential V with sink zero,

    phi(1)=(11327/12184) V.

Substituting V=6092/5137 recovers 11327/10274. Therefore an electrical quotient can be exact for boundary energy and still lose the parity data relevant to the integral critical group. One must keep the rational harmonic extension map, rather than treating a reduced rational-conductance network as the original integer-firing lattice.

The independent voltage checker solves and verifies L phi=e_u-e_w directly for every NS-source ordered pair through order six: 2, 18, and 726 pairs at orders 4, 5, and 6. It prints this and further actual odd-resistance-denominator hostiles. All common voltage denominators are even, as the theorem requires.

## 6. Wiring remains an independent coordinate

At order five, perform the following abstract switch in Q, using the frozen atlas node IDs:

    remove (1,4), (2,3);
    add    (1,3), (2,4).

Keep every blue edge unchanged and double the new Q edges to make the literal black network. This preserves each node's mass, blue count, black degree, reflection-pair degree and charge, the total number of edges, every loop count, and the entire blue graph. Therefore it preserves the critical-group 2-rank and all per-node parity laws above. But the physical spanning-tree count changes from 41096 to 39232.

This is an abstract parity- and degree-admissible rewire, not a realizable alternative tiling ownership assignment. Its purpose is precisely to show that the established masses and parity signatures do not determine the actual wiring. The complete tiling owner map or an equivalent incidence sidecar is still required. The electrical polynomial is a stronger diagnostic, not a complete invariant.

## 7. Operation and loss ledger

- Source: marked Hamiltonian tilings, exact isomorphism/converse owners, C-pairs, and R-pairs.
- Target: G=B+2Q, its colour spanning-tree polynomial, and its mod-two critical-group quotient.
- Preserved: exact literal conductance, loops as a separate stub ledger, colour counts on trees, and parity obstructions to integer currents.
- Forgotten: the actual Hamiltonian path realizing an endpoint, each individual R-pair label after collecting multiplicities, and all electrical effects of loops.
- Required sidecar: tiling-to-owner map, labelled C/R actions, and loop/stub counts. Tree counts and spectra alone do not restore them.
- Decisive hostiles: pure-blue order-five degree three; the order-five valuation mismatch; the order-six valuation/rank gap; the actual order-five Kron reduction; the degree-preserving order-five rewire.

No arbitrary orientation is used to turn the network into a tournament. Incidence orientations are computational gauges that cancel from the Laplacian and tree polynomial. There is no asserted transfer to LRC, a divisor-classification problem, or a trijection without an explicit operation-preserving map.

## 8. Reproduction

From this staging directory:

    python3 probe.py > probe.out
    python3 -O probe.py > probe.optimized.out
    cmp probe.out probe.optimized.out
    python3 hostile.py > hostile.out
    python3 voltage_audit.py > voltage_audit.out
    python3 -O voltage_audit.py > voltage_audit.optimized.out
    cmp voltage_audit.out voltage_audit.optimized.out

The source repository is read-only for this packet. All generated files are staged outside it. No Git operation or theorem-ID reservation is performed by this worker.
