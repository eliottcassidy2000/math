# The order-five cut aggregate is a distance-four graph, not a halved cube

2026-10-20. **PROVED, elementary with an explicit finite certificate.** This is a separate addendum; the frozen trijections/tree-gauge packet is unchanged. It concerns its cut-partition aggregate, not the regular-return icosahedra or the original tournament metagraph.

## 1. The proposed identification fails, but has an exact repair

The order-five cut aggregate has two components, each with 32 vertices and degree15. Each component is isomorphic to the graph

    V = {y in F2^6 : |y| is even},
    y adjacent to z iff |y+z|=4.                       (1)

It is **not** the halved six-cube, whose adjacency is Hamming distance2 on the same vertex set. Equivalently, (1) is the distance-two relation of the halved six-cube, not that graph itself.

Two decisive obstructions separate them:

| Graph | Adjacency spectrum | Local graph | Maximum clique |
|---|---|---|---:|
| Actual cut-aggregate component | 15^1, 3^10, (-1)^15, (-5)^6 | KG(6,2) | 4 |
| Halved six-cube | 15^1, 5^6, (-1)^15, (-3)^10 | triangular graph T(6) | 6 |

Superscripts are multiplicities. The aggregate weight is8 on every positive edge, so its weighted eigenvalues are120,24,-8,-40 with the corresponding multiplicities. The table above uses the unweighted support graph.

Thus the small coincidences K2 at order3 and the halved four-cube at order4 do not continue to a halved six-cube at order5. No all-order cube extrapolation follows.

## 2. Explicit linear certificate

Label the original vertices0,1,2,3,4 and use the star at0 to normalize cut classes. The six free coordinates are the pairs

    (12,13,14,23,24,34).

Write them as x=(a,b,c,d,e,f). At order5 every cut has even size, so edge parity is well-defined on the quotient. Its even component is

    a+b+c+d+e+f=0.

Define six output coordinates by

    y1=a+c+d,
    y2=b+c+d,
    y3=a+b+e,
    y4=b+c+e,
    y5=b+d+e,
    y6=c+d+e.                                         (2)

All arithmetic is in F2. These are parity rows with binary masks

    (13,14,19,22,26,28),

where the first pair12 is the least significant bit. Their sum is zero, so the image has even weight. The inverse on even output vectors is

    a=y1+y4+y5,
    b=y2+y4+y5,
    c=y1+y3+y5,
    d=y1+y3+y4,
    e=y1+y2+y3,
    f=y3+y4+y5.                                       (3)

Substitution, using sum_i yi=0, proves both inverse identities. Therefore (2) is a linear bijection of the two five-dimensional spaces.

For a Hamiltonian path P on the original five vertices, the aggregate generator is g_P=[E(K5)\P]. In the star coordinates its 15 distinct values are

    3,5,6,9,10,17,20,24,30,34,36,40,45,48,51.          (4)

Each has multiplicity8 when ordered Hamiltonian words are counted. Formula(2) sends these generators exactly to all fifteen weight-four vectors in F2^6. The resulting adjacency equivalence is

    x+z in {g_P} iff |L(x)+L(z)|=4.                    (5)

The other parity component is a translate of the first, so the same graph model applies after choosing a translation. The finite checker verifies every pair of component vertices, not merely the spectra.

## 3. Why the six coordinates are structural

The six coordinates in (2) are not arbitrary eigenvectors. The complete graph K5 has twelve undirected Hamiltonian cycles. The complement of each is another Hamiltonian five-cycle, so they form six complementary pairs.

For an even cut class x, attach to a complementary pair {C,E\C} the bit

    y_{ {C,E\C} } = |x intersect C| mod2.              (6)

This is well-defined in two ways:

- A cut meets any cycle in an even number of edges, so the bit is unchanged by changing the cut representative.
- Replacing C by E\C changes the bit by |x| mod2, which is zero on the even component.

Choose the following six representatives, each read as a cyclic word ending back at0:

    C1=(0,4,1,2,3,0),
    C2=(0,4,1,3,2,0),
    C3=(0,3,1,2,4,0),
    C4=(0,3,1,4,2,0),
    C5=(0,1,3,2,4,0),
    C6=(0,1,4,2,3,0).                                 (7)

Deleting the star edges at0 gives the six row masks above. Thus(6) is exactly(2), with its inverse already proved.

For the path P0=(0,1,2,3,4), the six intersection counts |P0 intersect Ci| are

    (2,1,1,0,2,2).

Since each Ci has five edges, the generator E\P0 evaluates to

    (1,0,0,1,1,1),                                   (8)

which has weight4.

### The exceptional six-point action

Every complementary cycle pair is also the pair of edge orbits of a cyclic subgroup of order5 acting regularly on the five original vertices. Conversely the rotation subgroup of either cycle recovers that subgroup. Hence these six pairs are exactly the six Sylow5-subgroups of S5.

The resulting six-point action of S5 is doubly transitive. It is transitive because all Hamiltonian cycles are relabellings of one another. A pair stabilizer has order120/6=20 and therefore a unique subgroup of order5. The order5 subgroup fixing one pair cannot fix another: otherwise it would be the unique order5 subgroup of both stabilizers and the pairs would coincide. Consequently it acts transitively on the other five pairs. This proves double transitivity.

This is not the action obtained by fixing one point among six and permuting the other five. That action splits the fifteen two-subsets into orbits5 and10. Here the six-point action is transitive on all fifteen two-subsets, exactly as the path-generator set requires.

The map(6) is equivariant for the **linear relabelling action on displacement edge fields modulo cuts**. An arbitrary relabelling of tournament state coordinates still has the affine observer correction from the main packet. No state-level linear equivariance has been smuggled into this statement.

Now all paths are relabellings of P0, so(8) and equivariance show that all generator images have weight4. Double transitivity gives every possible pair of zero coordinates, hence all fifteen weight-four vectors. This supplies a structural proof of(4)–(5) without inferring an isomorphism from spectrum alone. There are120 ordered paths and15 generators; transitivity makes every generator fibre equal, proving multiplicity8.

## 4. Spectrum, local structure, and the failed identification

Let A2 and A4 be the two adjacency matrices for Hamming distance2 and4 on even six-bit vectors, and let R send y to its antipode y+111111. Then

    A4=A2 R=R A2.                                    (9)

Indeed a difference has weight4 exactly when its complement has weight2. This changes one endpoint of the adjacency test; it is not conjugation by a vertex permutation and is not a graph isomorphism.

Characters of the even subspace are indexed by six-bit masks modulo the all-one mask. Representatives can be grouped by weights w=0,1,2,3, with multiplicities1,6,15,10. On such a character,

    eigenvalue(A2) = ((6-2w)^2-6)/2,
    eigenvalue(R) = (-1)^w.

Equation(9) gives the two spectra in Section 1. In particular their opposite-sign eigenvalues provide an immediate hostile to the halved-cube claim.

For the local graph of A4 at0, identify a weight-four neighbour with its missing two-subset of the six coordinates. Two such neighbours are adjacent exactly when their missing pairs are disjoint. Thus the local graph is KG(6,2), of degree6. A clique of pairwise disjoint pairs has at most three members, and a partition into three pairs attains it. Adding0 gives a clique of size4. Translation carries any clique to one containing0, so this is the global maximum.

By contrast the A2 local graph joins two-subsets that intersect in one point, giving T(6), of degree8. Its five edges through one fixed coordinate together with 0 make a six-clique. A pairwise-intersecting family of two-subsets has at most five members: if all pairs share one point this is immediate; otherwise pairs ab, ac, bc occur, and every further pair meeting all three must also lie in that triangle. Thus six is the maximum. The local graphs and clique sizes therefore separate the models independently of spectral calculations.

## 5. Reproduction and bounded scope

Run

    python3 check_distance_four.py > verification.json
    python3 -O check_distance_four.py > verification.optimized.json
    cmp verification.json verification.optimized.json

Both receipts are identical and report 5,977 explicit checks. The independent addendum checker imports no producer code. It enumerates the 120 Hamiltonian words, checks the inverse and every one of the 32^2 adjacency comparisons, computes both complete character spectra, checks the six cycle-pair action under all 120 permutations and all 32 even displacement vectors, and verifies the local maximum-clique obstruction.

The six-point action has 120 distinct permutations in that exact check, so it is faithful as well as doubly transitive. The action is derived from literal cycle pairs, not inferred only from its order.

The complete proof is restricted to order5. It preserves the cut-aggregate definition and its weight 8. It does not identify this graph with either regular-return icosahedron, does not lift quotient paths to fixed tournament representatives, and supplies no new original-metagraph connectivity claim.
