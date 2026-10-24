# Joint alignment of marked cyclic cells

**PROVED, scoped structural decoder.** This is a new tournament carrier for inherited affine holonomy, not a claim that affine fixed-point theory is new. The general fixed-section mechanism is already in `01-canon/theorems/THM-2622-affine-torsor-holonomy-fixed-section-spectrum-and-v4-c13-dictionary.md`; related inheritance is `01-canon/theorems/THM-2089-persistent-cut-affine-holonomy-normal-form.md` and `01-canon/theorems/THM-2637-derangement-character-fixed-branch-holotopy-principle.md`. The immediate input is `05-knowledge/results/incidence_wiring_20261023/MARKED-CELLS.md`, incoming commit e41a2c82b34fea4dffc6c1ee5dab07bfb8c56ded.

The new contribution is an exact joint decoder for actual cyclic-three-cell blocks, its automatic applicability to equitable triple partitions, and explicit tournament counterexamples when local chirality bits replace joint phase data.

## 1. The actual carrier, not a cosmetic tournament

Retain the names of a partition's cells. Each non-singleton cell D_i has three vertices and is a directed cycle. Choose a cyclic chart x in F3 so the internal arcs are x -> x+1. A different legal chart is x' = x+a_i, never an unrecorded reflection.

For each pair of triple cells, assume their binary adjacency block is either uniform, a permutation matrix, or the complement of a permutation matrix. Retain whether a nonuniform block has row sum one or two. Normalize the latter by complementing it. The normalized block selects a genuine bijection

    pi_ij : D_i -> D_j.

Every permutation of three elements is uniquely affine:

    pi_ij(x) = epsilon_ij x + delta_ij,
    epsilon_ij in {+1,-1}, delta_ij in F3.

The reverse transport is pi_ji = pi_ij^(-1). It is the SAME selected matching, not a second independent boundary decoder. The support graph G has triple cells as vertices and these nonuniform blocks as edges. Uniform blocks, all retained singleton attachments, and the entire tournament induced by singleton cells are separate fixed data.

The local chirality bit is epsilon_ij. It forgets the phase delta_ij. Rotating a cell chart changes phases by

    delta'_ij = delta_ij + a_j - epsilon_ij a_i.             (1)

This equation is obtained by substituting x'_i=x_i+a_i and x'_j=x_j+a_j into the actual selected matching. It is therefore a typed transport law for owned incidences.

### Automatic applicability to stable triple cells

Suppose a tournament partition is equitable for directed outgoing counts and every cell has size one or three. The internal outgoing degree in a triple is constant, so it is one; each triple is a directed cycle. For two triple cells the block has constant row sum q_ij. Its column sums are the constant 3-q_ji. Counting all nine entries gives q_ij=3-q_ji. A regular binary 3-by-3 block of degree zero or three is uniform; degree one is a permutation matrix; degree two is its complement. A triple-to-singleton block is uniform by equitability.

Thus **every such stable partition already has exactly the carrier above**. In particular the inherited seven-singleton/eight-triple two-mark partitions fall in this scope. Their actual refinement colors and ordered marks must remain retained. This theorem does not cover unresolved seven-vertex cells or arbitrary nonregular blocks.

## 2. Complete constructive classification

Fix cell identities, internal cyclic orientations, support, block polarities, uniform blocks, and singleton data. Then two remaining phase assignments give cell-preserving isomorphic tournaments if and only if they are related by (1).

Indeed, every isomorphism restricts on a directed C3 to a rotation a_i, and preservation of every nonuniform block is precisely (1). Conversely these rotations preserve all internal and uniform arcs, and (1) preserves every selected matching and its polarity. There are no forgotten arcs.

For each connected component choose a root and a spanning tree. Let g_i be the affine transport along the root-to-i tree path. For each oriented non-tree edge i->j retain the based holonomy

    h_ij = g_j^(-1) pi_ij g_i.                              (2)

Its linear sign is already determined by the retained epsilon values. Under a chart change (1), every based holonomy changes by the SAME root translation conjugacy. In coordinates,

    (E,D) -> (E, D+(1-E)a_root).                            (3)

The tuple of based chord holonomies, modulo simultaneous conjugacy (3), is a complete invariant. A constructive decoder first gauges tree phases to zero, reconstructs each chord phase from (2), and then fills every original block. The only remaining choice is the common root rotation. Comparing each chord independently up to its own conjugacy is insufficient: the same root chart must serve all chords.

Equivalently, (1) is the action of a linear map on phase vectors. For a connected support with v vertices and e edges, let kappa=1 if the product of signs around every cycle is +1, and kappa=0 otherwise. The kernel of this action has size 3^kappa, because a_j=epsilon_ij a_i propagates from one root and a negative cycle forces a_root=0. Consequently, with signs fixed, there are exactly

    3^(e-v+kappa)

phase orbits, and every phase assignment has rotational automorphism group of size 3^kappa. This is a cell-preserving count; ordinary uncolored isomorphisms can additionally permute cells unless their retained colors distinguish them.

## 3. Point sections, frames, and valid tournaments are different questions

A parallel owner section chooses one vertex x_i in each cell and requires x_j=pi_ij(x_i) on every support edge. Tree propagation determines all choices from x_root, so sections correspond exactly to the COMMON fixed points of every based holonomy in (2).

For one holonomy x -> E x+D over F3:

- E=+1,D=0 fixes all three points;
- E=+1,D!=0 fixes none;
- E=-1 fixes exactly one, solving 2x=D.

Several reflections may each fix one point but have no common fixed point. The explicit two-chord hostile h_1(x)=-x and h_2(x)=-x+1 fixes respectively 0 and 2, so no owner section exists.

A three-label parallel frame means THREE simultaneous disjoint sections, using all three points in every cell. Such a frame exists exactly when all based holonomies are identities. Frames here are unrestricted bijections with a root set. If each frame must also preserve the directed cyclic orientation, every individual edge sign must additionally be +1. A reflection on a tree does not obstruct a set-theoretic frame but does obstruct a rotation-only frame.

Most importantly, a tournament with nonzero holonomy is still a perfectly valid tournament. Lack of a section or frame obstructs a proposed synchronized-owner simplification, not the existence of the graph and not the uniqueness of its actual Hamiltonian-port wiring.

## 4. Exact small family classification

For three cells, make every i<j block a permutation matrix. There are 6^3=216 labelled block choices. The three cells have scores 3,4,5, so ordinary isomorphisms must preserve them. Among the eight edge-sign patterns, four have positive cycle sign and give three phase classes each; four have negative cycle sign and give one each. Thus the family has exactly sixteen ordinary isotypes: twelve have automorphism group C3 and orbit size nine; four have trivial automorphism group and orbit size twenty-seven.

The checker enumerates all 216 configurations, their 27 rotations, and exact Hamiltonian counts for the sixteen representatives. It also exhausts all 6^5=7,776 connections on a four-vertex, five-edge two-cycle support. The latter have 144 rotational classes, seventy-two with three automorphisms and seventy-two with one. The direct 3^4 owner-assignment search agrees with the common-fixed-point criterion.

## Connection and loss ledger

Source: a named equitable singleton/triple tournament partition. Target: fixed block data plus an affine F3 connection. Map: normalize each regular nonuniform block to its selected matching. Preserved: every actual arc, cell identity, and retained mark, modulo cyclic cell rotations. Lost by epsilon alone: the phase cocycle. Restoration sidecar: a common-root tuple of based chord holonomies, with only simultaneous translation conjugacy. Cheapest hostile: three cyclic cells with the same positive edge signs but three different cycle phases. The stronger exact-E hostile is proved in GIRTH_PHASE_HOSTILES.md.
