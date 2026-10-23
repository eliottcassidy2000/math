# Session38: incidence information, complete order-eight wiring, and the dyadic electrical tower

## Main outcome

The full order-eight off-path-complement metagraph is now reconstructed and independently replayed, including every fixed-Hamiltonian-path presentation, every endpoint port, all loops, and the paired black lines. A twelve-symbol triangle-incidence signature identifies its ordinary tournament classes. With authenticated converse reflection, a simple **three-color rule** suffices to force the complete wiring. It retains only the sum of the head and edge cyclic-triangle counts modulo three. Explicit larger tournaments prove that these signatures cannot become all-size class classifiers. An exact residual-flow criterion replaces the inference that a stalled capacity-peeling procedure necessarily leaves an ambiguity.

The electrical computation supplies an unexpected finite result: the order-eight graph has a cyclic 2-primary critical-group factor of order **256**, not merely order four. The complete certificate terminates at that exponent. Electrical invariants are nevertheless not sufficient to authenticate wiring.

This session does not establish all-order black-core connectivity, a general classification of tournament metagraphs, or a solution to the Hamiltonian-count spectrum problem.

### Concurrent repository work and attribution

During this session, incoming commit `f2cac6db1bb8ef780b9211d1cf3f7adc109c5fca` independently completed the full order-eight port replay and the same Z/256 electrical computation. Its [frontier](../metagraph_order8_frontier_20261023.md) and linked packets are preserved, not replaced. Thus the full census and the dyadic factor are independently corroborated results, not claims of priority for this packet.

That incoming work also supplies a guarded typed-out-wedge bit repairing the first-round v4 signature, a strong order-13 hostile to vertex-only refinement, and an order-31 doubly-regular pair repaired by a transitive-five motif. The present residue alphabets are a different declared menu; their symbol counts must not be compared to that conditional one-bit repair as if both encoded the same initial information. The new emphasis here is converse-implied marginal repair, exact projected-flow uniqueness, the stronger order-15 edge-incidence obstruction and its strong lifts, and localization of the high dyadic factor in the bridge-free core.

## 1. The complete finite object

The object is the converse-merged **off-path-complement line multigraph**, not the single-arc-flip graph. The frozen ledger convention fixes the directed path 7→6→⋯→0. Reversing vertex labels gives the equivalent ascending convention used in some general proofs. Complement flips precisely the other arcs. Reflection pairs black lines; the blue fixed presentations are retained separately.

| Quantity | Exact value |
|---|---:|
| Ordinary isomorphism classes | 6,880 |
| Self-converse classes | 176 |
| Converse-merged nodes | 3,528 |
| Fixed-path presentations | 2,097,152 |
| Literal complement lines | 1,048,576 |
| Blue lines | 2,048 |
| Black lines | 1,046,528 |
| Loops | 622 = 4 blue + 618 black |
| Distinct nonloop supports | 437,981 |
| Ordinary endpoint-port ledger entries | 956,339 |
| Directed merged-ledger entries | 460,699 |

The full graph is connected. The black subgraph has one component on 3,525 nodes and three isolated nodes. The blue graph connects all 176 self-converse nodes in one component and has 3,352 other isolated nodes.

There are eight support bridges, all incident to a node of ordinary fibre h=H/|Aut|=1. Three are blue single-line bridges. Five are black double-line support bridges: they are not bridges in the literal multigraph, but become bridges when each black twin pair is counted as one edge. No interior support bridge occurs at order eight. This is a finite conclusion, not an all-order leaf-bridge theorem.

The producer classified presentations by normalizing Hamiltonian paths of the ordinary representatives. The independent audit instead classified every mask using directly recomputed triangle-edge incidence. These distinct routes agree on every one of the 2^21 presentations, and the independent reconstruction agrees on every ledger entry, not just summary counts.

Sources: [wiring packet](wiring/README.md), [independent audit](wiring-audit/).

## 2. What the smallest tested signature retains

Let c_v count cyclic triples through v and c_vw count cyclic triples through the pair {v,w}. Define

E(T) = multiset_v multiset_(v→w) (c_w, c_vw).

The root value c_v need not be stored separately: it equals the sum of the second coordinates in its outgoing row. Row length also recovers outdegree. Thus the useful information is the *association* of a neighbor's triangle count with the incidence of the edge leading to that neighbor, not merely independent histograms.

Reducing the two coordinates modulo 4 and modulo 3 gives E_(4,3), a twelve-symbol row alphabet. It separates all 6,880 ordinary order-eight classes. The exact signature E, the residue signature E_(4,3), and the inherited J4 signature with one additional scalar all force the native endpoint-port table when tested against its authenticated capacities.

The twelve-symbol statement is minimal for ordinary class separation within the specified two-coordinate modular family, not among every imaginable encoding. The exact-neighbor/edge-parity model and the neighbor-modulo-two and neighbor-modulo-three models with exact edge counts each have class collisions. Keeping only unassociated per-row marginals also fails. A lookup table or an entirely different signature family is not covered by this minimality claim.

The scalar repair of J4 uses

q_u = sum_(u→w) binom(c_uw,2),
M = sum_(u→v) q_u q_v^2.

It separates both surviving J4 collision pairs. In contrast, appending the unassociated global edge-incidence histogram does not separate those pairs.

Class collisions and wiring ambiguity are distinct claims. Reducing both coordinates modulo three leaves four collision pairs. However, the converse of every colliding class has a globally singleton signature. Reflection-implied margins therefore resolve each ambiguous target individually. This forces complete wiring with nine symbols in the reflection-constrained model, even though the signature alone does not identify every class. Exact mirrored counterfeits for coarser menus establish the nine-symbol minimum within this modular family; this is not a minimum over arbitrary encodings.

The residual tables and their allowed cells are retained so nonuniqueness is verified by actual circulation witnesses rather than inferred from unsuccessful peeling.

### A three-color rule beyond the product-residue menu

The smaller successful rule is

color(u→v) = (c_v + c_uv) mod 3.

Each vertex retains only its three outgoing color multiplicities; the ordinary signature is the multiset of these triples. This has 6,864 ordinary signatures, not 6,880, so it does not identify every class directly. Nevertheless its converse signature-incidence graph is a forest on 13,728 vertices and 6,880 class edges, with 6,848 components:

- 6,820 isolated single-edge components;
- 24 two-edge stars;
- four paths of length three.

Four class edges are not forced by the singleton-or-converse-singleton shortcut. A first parallel peeling pass forces 6,876 target classes; a second pass resolves the four middle edges of the three-edge paths. These four classes have Hamiltonian counts 191,197,191,197 and form two converse pairs. Thus the full forest criterion, rather than paired-signature injectivity or a singleton heuristic, forces every ordinary port count and the complete merged wiring.

The bounded search considers quotients of the nine E33 symbols, not arbitrary color functions. Three colors are the smallest alphabet in that family with a global forest certificate: 117 three-color partitions pass, whereas no two-color partition does.

The stronger wiring minimum is established separately, not inferred from failed forests. Every one of the 255 binary partitions has a literal reflection-equivariant four-port switch changing the merged wiring while preserving endpoint/signature inventories, blue wiring, loop multiplicities, black-twin parity and the cyclic-triangle current. Sixty-one distinct switches cover all 255 menus; no minimum witness-cover claim is made. Thus **three colors are necessary and sufficient within all quotients of this fixed E33 alphabet for the declared reflection-constrained aggregate-wiring model**. This does not assert minimality over arbitrary invariants, uniqueness of every individual port partner, or realization of the altered matching by off-path complement. See [complete binary lower bound](binary-counterfeits/README.md).

There is an exact all-size decoding identity underneath the successful sum rule. Let s_v be the sum modulo three of the outgoing colors at v, and s*_v the corresponding sum at the same aligned vertex in the converse. Then

s_v+s*_v = c_v mod 3,

because the corresponding unreduced incidence sum is c_v + sum_w c_w and sum_w c_w=3c3(T). Recovering c_uv modulo three additionally requires the actual arc/head association, not merely the aligned vertex histograms: subtract the recovered head value from that arc’s color. Precisely, color_T(u→v)−color_Top(v→u)=c_v−c_u, an exact gradient modulo three. This local statement requires vertex alignment; the outer multiset signature deliberately forgets it. It is not an all-order unlabelled classification theorem.

Sources: [incidence packet](incidence/README.md), [literal hostile examples](incidence/HOSTILES.json), [ternary quotient search](incidence/TERNARY_SIGNATURES.md), [aligned decoder](incidence/THREE_COLOR_DECODING.md), [peeling depth](incidence/CONVERSE_PEELING_DEPTH.md).

## 3. An exact all-size criterion for a specified wiring relaxation

Fix an authenticated feasible integer port table, its allowed cells, integer margins and any individual cell bounds. Let P record its merged colored edge multiplicities. On the residual bipartite digraph give a forward cell its projected contribution and a reverse cell its negative.

**All feasible tables have the same projected wiring if and only if every directed residual cycle has zero projected sum.** Equivalently, the vector costs are potential differences on each residual strongly connected component.

A nonzero cycle gives an explicit integer alternate table by unit augmentation. A potential assignment gives a compact uniqueness certificate. Table uniqueness is the special case P=identity; a merged projection can be unique even when the table is not.

The theorem applies at every order, but only to the specified transport model. If reflection compatibility is imposed, it must be built into genuine port-orbit channels or verified separately. Additional global equalities cannot be silently discarded. A different feasible table in a relaxation is not automatically realizable as another tournament metagraph.

A complementary all-size symmetry criterion is especially efficient here: form the bipartite graph whose class edges join S(c) and S(c*). If that multigraph is a forest, the two reflection-related margin families determine every class count. A singleton signature on at least one end of every class edge is an immediate sufficient case. Injectivity of the ordered pair alone is not sufficient, since separate margins can conceal a rectangular cycle. See [converse marginal repair](theory/CONVERSE_MARGINAL_REPAIR.md).

Proof and scope: [projected transport uniqueness](theory/PROJECTED_TRANSPORT_UNIQUENESS.md). Independent bounded-model tests include cases with nonunique tables but unique projections: [rule audit](rule-audit/AUDIT.md).

## 4. The finite signature has a proved infinite boundary

Two strong regular circulants on Z/15, with outgoing difference sets

S1={14,13,12,11,5,6,8},
S2={14,13,12,4,10,9,8},

have identical exact E. Every vertex has c_v=28 and outgoing edge-incidence multiset (3,3,3,4,4,5,6). Yet their Hamiltonian path counts are respectively 196,264,905 and 196,277,235. Consequently every residue version E_(p,q) also fails to distinguish them.

Adding a source or sink preserves E equality and H. Therefore such H-distinguishing collisions exist at every order n≥15, without a strongness assertion for this padding.

For a strong extension use L_r(T)=C3[T,{a},P_r], r≥1, where P_r is transitive and the block directions are T→a→P_r→T. Exact E equality survives this lift. If F2(T) counts ordered spanning covers by two nonempty directed paths, then

H(L_r(T))=(2^r+1)H(T)+(2^r−1)F2(T).

For the explicit seeds the lifted difference is 113,070−88,410·2^r, which never vanishes for r≥1. Thus **strong collisions occur at order 15 and every order n≥17**. The strong order-16 case is not resolved by this construction.

A stronger incoming-work test adds every transitive-subset count, not merely the new repository's transitive-five separator. Both seeds still have

F_T(z)=1+15z+105z²+315z³+375z⁴+150z⁵+15z⁶,

and every vertex in both has the same outneighborhood transitive polynomial

1+7z+21z²+25z³+10z⁴+z⁵.

These two polynomial equalities persist alongside E equality under the source and strong cyclic extensions. Thus the strengthened observer still misses the different Hamiltonian counts at every order n≥15, strongly at 15 and every n≥17. This does not mean all induced motifs agree: the seven-vertex outneighborhoods contain zero versus one induced regular five-tournament, and their Hamiltonian counts are 101 versus 89. The source-over-regular-five motif repairs this pair with global counts zero versus 15. See [transitive-profile obstruction](theory/TRANSITIVE_PROFILE_OBSTRUCTION.md).

The all-size information limit can also be quantified. With L colors, an ordered row-count array has nL bins and total mass N=n(n−1)/2, so there are at most binom(N+nL−1,nL−1) possible signatures even before forgetting row order. Comparing this with the strong-tournament isotype lower bound proves that ordinary class separation for all sufficiently large n requires

liminf L(n)/n ≥ 0.1469076866702077…,

where the exact constant is defined by f(2a)=1 and f(x)=(1+x)log₂(1+x)−x log₂x. In particular, a sublinear alphabet cannot classify all strong tournaments through this row-multiset format. A simple uniform finite bound already gives failure whenever n≥64 and L≤n/16. These are necessary information bounds, not constructions or sharp thresholds. They do not prove failure of reconstruction from a larger authenticated inventory.

The seeds also agree in exact E after taking converse. Adding a source preserves this ordered pair of signatures, because converse turns that addition into sink padding. Consequently the global signature-incidence graph for exact E, E33 and every E33 quotient has parallel distinct class edges at every order n≥15. The global forest certificate itself must therefore fail there. This is not a larger-order port counterfeit: additional capacities or authenticated constraints may still force the actual wiring. See [all-size forest obstruction](incidence/FOREST_OBSTRUCTION.md).

Proofs: [collision propagation](theory/INCIDENCE_COLLISION_PROPAGATION.md), [finite-alphabet boundary](theory/FIXED_ALPHABET_BOUNDARY.md), [sharper growing-alphabet bound](theory/GROWING_ALPHABET_BOUNDARY.md).

## 5. What the electrical hierarchy adds, and what it cannot add

Write the literal graph as G=B+2K, where B is blue and K counts paired black lines. Loops stay in the graph inventories but contribute zero to the Laplacian. The factor two is an exact orbit multiplicity, not a physical assumption.

For a grounded integer Laplacian A, let V_k consist of binary potentials admitting an integer lift x with Ax divisible by 2^k. The associated divided energy form q_k has radical V_(k+1); its rank counts Smith factors of exact 2-adic valuation k. This is standard integral linear algebra applied to the authenticated tournament network, not a new general Smith-normal-form theorem.

At order eight the complete, independently checked sequence gives

K(G_8)[2^∞] ≅ (Z/2)^3358 ⊕ (Z/4)^4 ⊕ Z/256,

and hence v_2(number of spanning trees)=3374. The last one-dimensional form vanishes at levels 3,4,5,6,7 and is nonzero at level 8. Integer lifts and their Laplacian images are provided at every stage, so 256 is certified rather than inferred from a precision cutoff.

The order-256 factor is not a leaf-bridge artifact. Removing the eight leaf bundles leaves a connected 3,520-node core with no support bridges and minimum support degree three. Its 2-primary group is (Z/2)^3353 ⊕ (Z/4)^4 ⊕ Z/256. The five removed black double bundles contribute exactly five additional Z/2 factors; the three blue singleton bundles contribute none. The all-size bridge decomposition and directly checked core certificates are in [bridge localization](electrical/BRIDGE_LOCALIZATION.md).

A paired-black switch changes the first divided form by a rank-at-most-two bilinear expression. This gives an all-size exact test for which such switches that form detects. But at order five a synthetic degree-preserving switch leaves the entire dyadic tower unchanged while changing the spanning-tree count from 41,096 to 36,680. Thus even the complete 2-primary electrical invariant cannot force the wiring. Moreover, any change of individual port pairing that leaves edge multiplicities unchanged is invisible to the whole Laplacian.

Proofs, exact lift certificates, and the hostile control: [electrical packet](electrical/README.md), [electrical proof](electrical/PROOF.md).

## 6. The structural lesson and next precise obligations

Three different levels must not be identified:

1. Local incidence signatures compress a tournament's cyclic interactions.
2. Port capacities and allowed-cell incidence constrain how those objects can be connected.
3. Electrical forms observe global responses of the resulting merged network.

The common algebraic theme is concrete: reconstruction depends on the kernel of the retained observations. Signature-incidence forests have independent columns, so reflected margins have no remaining class-count kernel. Residual cycles explicitly exhibit directions invisible to those margins. The electrical critical group instead records an integral Laplacian cokernel, and its dyadic tower recovers divisibility erased by reduction modulo two. The ternary color difference is an actual gradient modulo three. These are compatible tools for tracking information loss, not assertions that the underlying spaces or the two primes are interchangeable.

Each level supplies useful necessary constraints, but a compression at one level need not retain information required at another. The exact residual-cycle theorem specifies when the retained port information is sufficient; the explicit circulant and electrical hostiles specify where particular compressed invariants lose information.

The next targets are to strengthen the allowed-port signatures while retaining inexpensive updates, classify which constrained residual circulations survive reflection, and seek structural proofs for the observed black-core and support-bridge patterns. A universal fixed-alphabet outgoing-row-multiset class signature or a universal exponent-four law is no longer a defensible target. This does not exclude richer observation formats or reconstruction from a full authenticated inventory.

All finite census claims here are order-eight claims unless explicitly stated otherwise. No inference to Collatz, Hadwiger–Nelson, Yang–Mills, or the Riemann hypothesis is asserted by numerical analogy.
