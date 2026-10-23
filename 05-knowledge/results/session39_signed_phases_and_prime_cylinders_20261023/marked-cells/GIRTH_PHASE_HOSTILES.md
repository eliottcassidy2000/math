# Hidden cycle phase in strong self-converse tournaments

**PROVED:** the scoped girth blindness, infinite tournament family, marked-motif thresholds, and one-mark repair below. **FINITE-EXACT:** the specified controls in the accompanying receipt. The general connection and its legal gauges are defined in AFFINE_CONNECTIONS.md.

## 1. Exact incidence is blind to phase when the support has girth at least five

Fix a signed matching support G of girth at least five, all block polarities, all uniform blocks, and the oriented cyclic triple cells. Then varying the phases delta leaves the exact incidence signature E unchanged in BOTH orientations. More strongly, each vertex's outgoing E row is unchanged as a multiset, and all vertices of a given cell have the same such row.

Here c_v counts cyclic triangles through v, c_uv counts those through {u,v}, and

    E(T) = multiset_u [ multiset_{u->v} (c_v,c_uv) ].

Proof. Any union of at most three cells has forest support. Its phases can be gauged away by cyclic cell rotations. This proves that c_v is constant on its cell and independent of phases, including with a specified root vertex fixed.

For c_uv, separate third vertices by their cell. If u,v lie in one cell, rotating other cells removes every phase without moving u or v. If their distinct cells i,j are support-adjacent, there is no common support neighbor, because that would make a triangle. Every contribution depends on at most the one direct relative coordinate y-pi_ij(x); all other incident phases can be removed. If i,j are not support-adjacent, at most one cell is a common support neighbor, because two would make a four-cycle. Only that cell can introduce a nonconstant contribution, and it uses just the one affine relative coordinate obtained by transport along its two-edge path. All other third-cell contributions are constant in the two root coordinates after forest gauges.

For a fixed root u and a target cell D_j, the head value c_v is constant. In the adjacent case, the outgoing restriction selects either the zero direct relative coordinate or its two nonzero values, and c_uv depends on that same coordinate. In the nonadjacent case, the block is uniform, so all three target vertices or none are used; the single relative coordinate ranges through all of F3. In either case the outgoing multiset contribution from D_j is phase-independent. The internal-cell contribution was already handled. Reversing all arcs merely exchanges block polarities and keeps cyclic-triangle incidences, so the proof also gives E(T^op).

This argument uses the absence of triangles and four-cycles. It does not assert that radius-two neighborhoods are trees, which would be false on a five-cycle.

## 2. An explicit infinite family with three invisible phases

For every g>=3, take cells D_0,...,D_(g-1), each labelled by F3 with x->x+1 internally. For every i<j:

- on the path pairs (i,i+1), use forward matching arcs y=x;
- on the closing pair (0,g-1), use forward matching arcs y=x+delta;
- on every other pair, direct all arcs from the lower cell index to the higher.

Write T_(g,delta), with delta in {0,1,2}. Its cell scores are

    s_0=3g-6,
    s_i=3g-2-3i       (1<=i<=g-2),
    s_(g-1)=5.                                                (1)

These are distinct: the first is 0 mod3, interior scores are 1 mod3 and distinct, and the last is 2 mod3. Consequently the cells are intrinsic ordinary-tournament invariants. Any isomorphism must rotate the cells individually. Path matchings force all rotations equal, so the closing phase cannot change. The three tournaments are pairwise nonisomorphic, and each has automorphism group exactly C3.

Each cell is strongly connected. Every support edge has actual arcs in both directions between its endpoint cells, so the connected cycle support makes the whole tournament strong.

Each tournament is self-converse. The literal map

    R(i,x) = (g-1-i, -x)

reverses every arc: it reverses each internal C3, reverses cell order on uniform blocks, preserves the selected path matchings, and preserves the closing matching equation after exchanging endpoints. Thus the three phases are three distinct converse-merged classes, not just converse aliases.

For every g>=5, Section1 now gives THREE distinct strong self-converse tournaments of order 3g with identical exact E, hence with identical E33 and all its quotient palettes. This is an infinite structural collision family; nonisomorphism does not require an unproved formula for Hamiltonian counts.

### Small exact controls

For delta=0,1,2 respectively, exact Hamiltonian-path counts are:

| g | order | H values |
|---|---:|---|
| 3 | 9 | 2697,2625,2721 |
| 4 | 12 | 401085,401397,401721 |
| 5 | 15 | 45132777,45150825,45102765 |

At g=3, all vertices have c_v=9, but exact E and the full ternary palette distinguish the phases. At g=4 they also distinguish them. At g=5, exact E in both orientations and the full ternary signature agree; the per-cell c_v values are 21,17,25,17,21, each repeated three times. The checker retains literal adjacency bitrows and both common exact signatures.

These collisions defeat the global converse-signature forest certificate, even on the self-converse inventory: distinct class edges are parallel. They do not establish nonuniqueness of any particular actual endpoint-capacitated port table.

## 3. Exact locality and marked detection thresholds

More generally, a specified collection of cells whose support is a forest has isomorphic induced unions for every phase assignment: solve the rotation gauges along that forest. Therefore every induced motif count with fewer than girth(G) vertices is phase-independent. Cell names may be retained. This is an exact finite statement about actual adjacency, not a physical effective-field-theory assertion.

For the cycle family, the two nonzero phases have identical cell-marked induced motif counts through g vertices. A subset using fewer than all g cells is handled by the forest argument. A subset using all g cells and at most g vertices is a transversal, one vertex per cell. Simultaneously replacing every chosen coordinate x by -x interchanges the two phases and preserves every induced arc, because no internal C3 arc is present.

The following marked observations show sharp thresholds, in the cell-marked sense.

1. With one vertex in every cell, require all selected support matchings. The count is three for delta=0 and zero for either nonzero phase. Thus zero versus nonzero is first detectable by some cell-marked induced motif at g vertices.
2. Use two vertices a,b in D_0 with a->b, and one vertex in each other cell. Require matching along the open path starting at a, and require the closing matching to join b to the last vertex. The path forces every other coordinate to equal a; b=a+1, so the closing constraint is a=a+1+delta. There are three such configurations exactly when delta=-1, and none when delta=+1. The remaining arcs are forced, so this is a particular cell-marked induced motif, not merely a partial compatibility statistic. Consequently the nonzero phases are first detectable by some cell-marked induced motif at g+1 vertices.

Do not promote these marked thresholds to an all-g statement about the first differing UNMARKED motif count. Forgetting cell names can merge several motif types and cause cancellations. The scoped marked theorem suffices to exhibit arbitrarily deep local indistinguishability on actual strong tournaments.

## 4. An actual mark repairs this family

Retain the intrinsic cell colors from (1), and distinguish one actual vertex. Directed neighbor refinement splits its C3 in one round: the other two vertices differ by whether they point to the mark. Once a triple is fully resolved, every support-neighbor triple is resolved at the next round, because its normalized matching gives three distinct outgoing records into the three known vertices. Complementary matching blocks work identically.

Resolution therefore propagates through connected support. Starting from any actual vertex of T_(g,delta), the partition becomes discrete within

    1 + eccentricity of its cell <= 1+floor(g/2)

rounds after cell colors have been supplied. A discrete colored outgoing record recovers every arc. If starting from a completely uniform seed, allow a preliminary degree-refinement round to obtain the intrinsic cells.

This is a scoped one-mark repair. It does not contradict the inherited DRT theorem that one mark can remain stuck at {r},N+(r),N-(r). It also demonstrates why absence of a global matching section is not an obstruction to reconstructing the actual tournament.

## 5. Compatibility with the Session38 ternary gradient

For an actual tournament arc u->v, the Session38 color chi_T(u->v)=(c_v+c_uv) mod3 satisfies

    chi_T(u->v)-chi_(T^op)(v->u) = c_v-c_u mod3.

This is an exact vertex gradient and sums to zero around every closed walk. In the cycle family the c_v values are independent of the attachment phase, so this gradient cannot recover the nonzero cell-connection holonomy. At g=3 it is identically zero because every c_v=9, although the three graphs have different H. The FULL ternary signature still distinguishes those order-nine controls. At g>=5 even exact E loses the attachment phase.

The two carriers must not be conflated: the gradient is a triangle-incidence observable on actual tournament vertices, while delta is transport of selected owners on the graph of triple cells. In the positive-sign case, delta is literally a graph F3-valued 1-cochain, rotations add a coboundary, and cycle phase is its cohomology class. This is a precise algebraic connection to flux and local invisibility, not a claim about physical electromagnetic fields or an unrestricted metagraph reconstruction rule.
