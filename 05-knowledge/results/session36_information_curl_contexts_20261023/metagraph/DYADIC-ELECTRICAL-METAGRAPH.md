# The metagraph has a dyadic electrical structure

Status: PROVED elementary graph/lattice statements below, specialized to the actual fixed-Hamiltonian-path line metagraph. FINITE-EXACT independent complete tiling censuses through order seven; integer Smith forms through order six; complete 2-primary Smith data through order seven. This is a new connection within this session's bounded source recovery, not a claim of literature priority. No all-order connectivity or wiring classification is asserted.

## 1. What was recovered, and what is genuinely being added

The source is the converse-merged **line** metagraph, not the single-arc-flip graph and not the one-tile-flip graph. A line pairs a fixed-path presentation with the presentation obtained by reversing all off-path arcs. Reflection of the fixed path is a second commuting involution. Its fixed presentations are blue; all other presentations are black.

The key older results are:

- THM-643, `01-canon/theorems/THM-643-blue-black-line-parity-structure.md`: odd ordinary presentation fibres, odd positive blue fibre on self-converse classes, and no blue presentation on a non-self-converse class.
- THM-4731, `01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md`: every black line has a distinct parallel reflected twin, including black loops. Thus the literal graph is Blue + twice PairedBlack.
- `05-knowledge/results/boundary_support_20261020_metagraph.md`, sections 2–4: the lost-fibre Gram defect, its rank, the exact bridge criterion after black contraction, and abstract commuting-involution realization. These are already proved, not rediscovered here.
- THM-4739, `01-canon/theorems/THM-4739-kirchhoff-and-hodge-laws-of-the-tournament-metagraph-intransitivity-currents-second-order-determination-and-the-consecutive-start-port-ledgers.md`: score/cyclic-triangle Hodge energies and intrinsic class currents. Its higher local wiring data do not determine the actual network: order six has 48,000 class-level solutions at the reported third level. This is a finite census, not an all-order impossibility theorem.
- THM-2467, `01-canon/theorems/THM-2467-bicycle-spaces-of-the-star-flip-split.md`: the binary cut/cycle intersection is a real obstruction, not a direct-sum slogan. Its graphs are K_n and K_n minus the fixed path, **not** the blue line graph studied below.
- Session33 `session33-metagraph-memory/NOTE.md` and `ADDENDUM_FULL_FLIPS.md`: literal line complement, fibre rerandomization, minimum linear memory and the full-labelled-flip obstruction are different mechanisms. They are not strengthened by renaming the graph Laplacian.

Two old reflection hazards matter. `the-blue-black-flip-lines-are-an-even-graph-plus-a-t-join-grid-reflection-is-the-complement-opus-20260701.md` says black loops are fixed by reflection times complement. That does not survive the current fixed-tile argument: this product has no fixed presentations for n>=3, and black loops also come in distinct reflected pairs. Likewise `the-half-tiling-is-an-abelian-square-complex-the-metagraph-is-an-anti-LTC-and-cut-cycle-is-alexander-duality-kps.md` claims an orthogonal direct cut/cycle split over F2 at even order; THM-2467 supplies the nonzero bicycle correction. The present result uses the repaired statements.

New scoped contribution: the twin involution controls the **integral Laplacian**, not just degrees. Its mod-two shadow is exactly the blue graph; the next dyadic layer is an explicit half-Laplacian pairing in which black wiring enters. This yields exact critical-group factors and a fake-wiring obstruction. A two-color spanning-tree polynomial then connects the old black quotient to a real conductance degeneration.

Concept board: involution multiplicity; blue harmonic modes; bicycle intersection; integral current solvability; black-component quotient; color-refined tree bases. Anchor: exact arithmetic electrical structure. Niche: the forgotten bicycle coordinate. Wildcard: real conductance collapse versus dyadic lifting.

## 2. Definitions and scope

Let G be any finite undirected multigraph whose blue subgraph is L and whose black part consists of two copies of a multigraph B:

    G = L + 2B,             Lap(G) = Lap(L) + 2 Lap(B).

All graphs retain the same vertex set, so isolated vertices count as components. Loops count twice for fibre degrees but contribute zero to electrical Laplacians, cuts, spanning trees and incidence columns. Parallel edges remain separate.

Write c_H for the number of connected components of H. Over F2 let C_H be a vertex-edge incidence matrix. The cycle space is ker C_H, the cut space im C_H^T, and

    bicycle(H) = ker C_H intersect im C_H^T,
    beta_H = dim_F2 bicycle(H).

For connected G, its critical group is the finite cokernel of a reduced integer Laplacian A_G, equivalently the group of zero-sum integer currents modulo Lap(G) applied to integer potentials. Its order is the spanning-tree count kappa(G). For disconnected G, take the product of the critical groups of its components.

An all-order actual tournament graph satisfies G=L+2B. The following theorems do not assume that the actual G is connected. Where a connected formula is used, apply it to each component.

## 3. The first dyadic theorem: blue topology is the exact 2-rank

For every graph G=L+2B,

    dim_F2 K(G)/2K(G) = c_L - c_G + beta_L.                 (1)

Equivalently, exactly c_L-c_G+beta_L of its nontrivial Smith invariant factors are even. Hence

    |K(G)[2]| = 2^(c_L-c_G+beta_L),
    v_2 |K(G)| >= c_L-c_G+beta_L.                           (2)

These are componentwise statements if G is disconnected.

Proof. Since Lap(G)=Lap(L) modulo two, their binary ranks agree. The map

    ker(C_L C_L^T) -> bicycle(L),        x -> C_L^T x

is onto: any bicycle vector is C_L^T x and lies in ker C_L. Its kernel is ker C_L^T, consisting of the functions constant on each blue component, of dimension c_L. Thus

    nullity_F2 Lap(L) = c_L + beta_L.

Deleting one row and column in each connected component of G reduces nullity by c_G without changing rank: the full Laplacian is integrally congruent to the reduced Laplacian plus one zero coordinate per component. Finally, the nullity modulo two of an integer full-rank matrix equals the number of even Smith invariant factors. This proves (1) and (2).

In particular, each non-self-converse merged tournament node is blue-isolated. If the full graph is connected and the self-converse blue vertices form one connected component, then

    2-rank K(G) = number of non-self-converse merged nodes + beta_L.

The connected-blue hypothesis holds in the present n<=7 census but is not promoted to an all-order theorem.

### Current interpretation

For an integer potential p, the current j=Lap(G)p satisfies

    u^T j = 0 mod 2                 for every u in ker_F2 Lap(L).    (3)

Indeed Lap(G) reduces to Lap(L). Thus every blue component receives an even total current, and the bicycle modes provide additional independent parity tests. Conversely, these tests characterize solvability of the potential equation **modulo two**, not over the integers. Higher divisibility layers remain.

The exact source-to-target map is G -> its integral critical group -> its reduction modulo two. It preserves parity obstructions to integral current transport; it loses odd torsion and every higher dyadic exponent. The bicycle space is the missing sidecar beyond component counts. It is not the whole cycle space or the ordinary first Betti number: at order seven the blue cycle space has dimension 169, whereas the blue bicycle space has dimension 11. Nor is it Session33's real observer-defect rank, which concerns different information and a different linear operator.

## 4. The second dyadic theorem: black wiring is a half-energy pairing

First assume G connected. Define the intrinsic binary potential space

    V = ker_F2 Lap(L) / <constant potential>.

A choice of root identifies V with ker(A_G mod 2), but the definition itself is unrooted. For u,v in V choose any integer vertex-potential lifts and define

    q_G(u,v) = (u^T Lap(G) v)/2 mod 2.                     (4)

This is a well-defined symmetric bilinear form over F2. Its radical has dimension

    dim rad(q_G) = number of Smith invariant factors of K(G)
                  divisible by 4.                        (5)

Moreover its exact color decomposition is

    q_G(u,v) = (u^T Lap(L) v)/2 + u^T Lap(B) v mod 2.      (6)

The first numerator is even because u and v are blue-harmonic modulo two. In edge language the terms are one half of the blue gradient-pairing sum plus the ordinary paired-black gradient-pairing sum.

Proof of well-definedness. Lap(G)v is even, so (4) is integral. Replacing u by u+2a changes (4) by a^T Lap(G)v, which is even. The other lift is analogous. Adding a constant changes nothing. The same calculation handles addition of binary vectors and proves bilinearity. Symmetry is inherited from the Laplacian.

Proof of (5). Use a reduced matrix A and a binary kernel vector x. It lifts to a vector x+2y satisfying A(x+2y)=0 modulo four precisely when

    A y = -A x/2 mod 2.

Over F2, symmetry gives im A=(ker A)^perp. Thus this is solvable precisely when w^T A x/2=0 for every w in ker A, or x lies in rad(q_G). In Smith coordinates, a kernel vector modulo two can lift to a kernel vector modulo four exactly in the coordinates whose invariant factor is divisible by four. Hence the radical dimension is the stated count. The argument is invariant under integral changes of potential basis, so the root and the chosen kernel basis do not matter.

For disconnected G, quotient constants separately on every full component and apply the same proof.

This is the first place in the dyadic hierarchy where **how the black edges are wired** matters. The mod-two operator erases every black twin; division by two exposes one representative of each twin. This is a standard modular-lifting mechanism applied to the actual involution graph, not a new general theorem about Smith forms.

### Why one should not replace this pairing by its diagonal

In characteristic two a symmetric bilinear form need not be determined by its diagonal. Its full matrix, including off-diagonal interactions between blue harmonic modes, is required. Nor do the separate ranks of its blue and black pieces determine the rank of their sum. More sharply, if the entire blue graph and every paired-black degree are held fixed, the diagonal of q_G is fixed: the black self-pairing u^T Lap(B)u reduces to the degree-parity dot product. Therefore the degree-preserving black hostile in section 6 isolates genuinely off-diagonal information, not a change of local self-energy.

## 5. Exact actual-source census

The checker starts with all 2^binom(n-1,2) tournaments containing 0->1->...->n-1. A one bit on (i,j), j>=i+2, means i->j. Complement toggles every bit. Reflection sends (i,j) to (n-1-j,n-1-i), without bit negation. Blue means equality with this reflection, not merely a self-converse unmarked class.

Canonicalization enumerates every permutation inside equal-outdegree blocks and minimizes both tournament and converse codes. Sorting the outdegrees only removes permutations that cannot preserve the degree-labelled canonical order. The graph is built from literal complement pairs once each. All black multiplicities are checked even, not imposed.

| order | vertices | c_L | beta_L | critical 2-rank | v2 kappa(G) | nonzero 2-primary exponents |
|---:|---:|---:|---:|---:|---:|---|
| 3 | 2 | 1 | 0 | 0 | 0 | none |
| 4 | 3 | 2 | 0 | 1 | 1 | 1 |
| 5 | 10 | 3 | 1 | 3 | 3 | 1,1,1 |
| 6 | 34 | 23 | 3 | 25 | 26 | twenty-four 1s; one 2 |
| 7 | 272 | 185 | 11 | 195 | 198 | 192 ones; three 2s |

For example, the 2-primary part is (Z/2)^24 plus Z/4 at order six, and (Z/2)^192 plus (Z/4)^3 at order seven. These are exhaustive finite statements, not a guess about all orders. They do not assert that the odd-primary part is trivial.

Integer Smith forms and exact tree determinants are independently compared through order six. At order seven, unit-pivot elimination modulo 2^32 determines every 2-primary exponent; the algorithm must end with no remaining coordinates before the precision expires. It does, at exponent two. Independent binary incidence ranks and the half-pairing give the same 2-ranks and numbers of exponents at least two. Normal and optimized Python use explicit exceptions for all checks.

The half-pairing radical dimensions are 0,0,0,1,3. The blue contribution in (6) has rank zero through order six and rank two at order seven; its separate rank is not enough to infer the full answer.

A useful hostile to monotonic torsion intuition is the blue graph itself. Its nonzero 2-primary exponents are:

- order five: 2;
- order six: 2,2,6;
- order seven: 1,1,3,3,4,4,5,5,6,9,10.

Adding doubled black edges therefore does not monotonically preserve the higher exponents. At order seven the blue network has a factor of order 1024 while the full network's 2-primary exponent is only 4. This is a concrete distinction between adding electrical conductance and preserving an integral quotient group.

### Uniform odd color rescaling, with a finite-order scope

For any fixed colored graph, multiplying all blue conductances by an odd positive integer u and all literal black conductances by an odd positive integer v preserves the first two dyadic counts: the mod-two operator is unchanged and the half-form becomes u q_L+v q_B=q_L+q_B modulo two.

For each actual graph of order 3 through 7, a stronger exact conclusion holds: **its entire 2-primary critical group is unchanged by every such odd pair (u,v)**. Multiply the reduced matrix u A_L+2v A_B by the odd unit v^(-1) modulo eight. It becomes t A_L+2 A_B with t in {1,3,5,7}. The checker exhausts these four cases at each order and finds the same valuations 0,1,2 with no invariant factor divisible by eight. Matrix congruence modulo eight determines these truncated Smith valuations, and their exhaustion makes them complete. Thus this is a uniform statement in the two conductance parameters for five fixed finite graphs, not an all-order tournament theorem.

## 6. A fake wiring invisible to node/color degrees

Using the local canonical node labels of the exact order-five census, replace blue edges

    0--5, 6--9     by     0--6, 5--9,

and leave every black edge unchanged. Every blue degree, black degree, ordinary fibre size, blue/black node category and global involution orbit count is unchanged. Both colored networks remain connected. But

    actual critical 2-rank = 3,      rewired critical 2-rank = 4;
    actual tree count = 41096,       rewired tree count = 37184.

This is a **synthetic colored-multigraph control built from actual data**, not a second tournament-realizable metagraph. The older abstract commuting-involution construction realizes the modified graph at that abstract level, but actual tournament wiring is fixed by the concrete complement operation and the modification changes it.

The first dyadic invariant already separates these two degree-identical networks. Thus the electrical arithmetic adds a genuine global obstruction beyond every local color inventory; it is not just another expression for fibre parity. The receipt persists the complete original and modified edge tables, so local node IDs are unambiguous.

### A stronger hostile invisible to the entire first dyadic layer

Keep the **whole actual order-five blue graph** unchanged. In the paired-black graph replace

    1--4, 3--7     by     1--3, 4--7.

Each displayed edge represents two literal black twins. Every paired-black degree remains unchanged, so again every color inventory and global involution count survives. Both graphs are connected and the full mod-two Laplacian is identical. Nevertheless,

    actual 2-primary group:  (Z/2)^3,
    rewired 2-primary group: (Z/2)^2 plus Z/16.

The exact tree count changes 41096 -> 37952. The half-energy pairing drops from rank 3 to rank 2 on the same three-dimensional blue harmonic space, exposing one 4-divisible factor; further dyadic lifting identifies its exact order 16. Thus all first-layer data and every local inventory can agree while the second layer detects a different integral network.

This is again an abstract colored-network hostile, not an actual tournament realization. It proves in particular that the abstract commuting involutions and actual node degrees cannot force the finite full-network exponent-four pattern.

## 7. Real conductance scales: a two-color tree theorem

Let G=L+2B be connected on r vertices. Give each blue edge conductance z and each literal black edge conductance1. Then

    P_G(z) = det(z A_L + 2 A_B)
           = sum_T 2^(r-1-b(T)) z^b(T),                  (7)

where T ranges over spanning trees of the graph with one representative per paired-black edge, and b(T) is its number of blue edges. Parallel edges remain distinct. Equation (7) is the matrix-tree theorem with the exact twin multiplicity retained.

The nonzero exponents form the complete interval

    c_B - 1 <= b(T) <= r - c_L.                           (8)

The endpoints are attained by taking a tree within every monochromatic component and a tree in the graph obtained by contracting those components. Every intermediate exponent occurs because spanning trees are connected under one-edge basis exchanges, which change blue count by at most one.

Let Q_L be the blue graph after contracting all B-components; let Q_B be the paired-black graph after contracting all L-components. Then

    [z^(c_B-1)]P_G
      = 2^(r-c_B) * product_i kappa(B_i) * kappa(Q_L),

    [z^(r-c_L)]P_G
      = 2^(c_L-1) * product_j kappa(L_j) * kappa(Q_B).      (9)

Proof. A tree attaining the lower bound must use a spanning tree in each black component, then precisely a quotient blue tree. The upper bound is the same argument with colors interchanged. The factors of two count choices of a black twin for every selected paired-black edge. Isolated components have tree count 1.

This recovers the old bridge quotient as the first nonzero coefficient of an exact electrical partition function. In particular a blue bridge is in every quotient tree, while a nonbridge has some quotient tree omitting it. No assumption that pure-blue nodes are leaves is made.

Exact polynomials include

    P_3(z)=z,
    P_4(z)=2z,
    P_5(z)=10560z^3+19296z^4+9584z^5+1576z^6+80z^7.

At order six all exponents 2 through 11 occur; the integer coefficients are in `tree_scales_receipt.json`. Coefficients are checked by exact determinant interpolation with the proved degree bound and independent test evaluations. Both endpoint product formulas are checked independently. In particular every coefficient at exponent k is divisible by 2^(r-1-k), as its literal twin choices demand.

### A conductance-limit interpretation

For the family Lap(G_epsilon)=2Lap(B)+epsilon Lap(L), contracting black components is the exact leading-order electrical limit. If u,v lie in different black components,

    epsilon R_G_epsilon(u,v) -> R_Q_L([u],[v]).             (10)

If they lie in the same black component C,

    R_G_epsilon(u,v) -> (1/2) R_C(u,v).                   (11)

These standard finite-network limits follow by minimizing flow energy: at scale epsilon, finite-cost intercomponent transport is entirely through quotient blue edges; within a black component all internal routing costs vanish after multiplying by epsilon. At unscaled finite energy, an intercomponent blue flow tends to zero, leaving the internal black network. Alternatively, split potentials into black-component constants and their orthogonal complement, invert the nonsingular black blocks and use the Schur complement. These limits are a corollary of the precise colored network, not an identification of the literal complement dynamics with a Markov walk.

## 8. Hostile extension: the blue form is not universally symplectic

The blue contribution to (6) happened to be alternating through order seven. This is not a general consequence of reflection or of the blue parity law. A separate complete **blue-only** census at order eight enumerates 4096 reflection-fixed presentations, finds 176 self-converse nodes and a connected blue graph, and gives

    beta_L=11,
    2-primary exponents of K(L) = 1,1,1,1,2,2,3,3,3,3,4.

Its half-energy form has rank 4 and is **nonalternating**. The receipt exhibits a subset S whose blue cut is simultaneously a cut and a cycle (every crossing degree is even), yet contains 914 literal blue edges, which is 2 modulo 4. Therefore

    q_L(1_S,1_S) = |delta_L(S)|/2 = 1 mod2.

All 2048 harmonic potentials modulo constants were checked; 914 is the minimum crossing size among those with nonzero self-pairing in this local canonical labeling. The full edge-multiplicity witness and source-independent producer are retained. This is an exact order-eight blue graph result, **not** a census of the entire order-eight blue-plus-black network.

Mechanism and strongest survivor: bicycle cuts are automatically even-sized, but need not be divisible by four. The correct object is a symmetric binary half-energy form, not automatically an alternating/symplectic form. The general theorems above deliberately retain its diagonal.

## 9. What this says most fundamentally

The metagraph is simultaneously an involution incidence object and an **integral electrical network**. The involution is not decorative symmetry: it forces a factor-of-two multiplicity on every black conductance. Consequently:

1. Over the reals, weakening blue conductance exposes black connected components and the old bridge quotient.
2. Modulo two, black conductances disappear and blue topology controls the whole first obstruction space.
3. At the next dyadic layer, dividing the energy form by two makes the detailed black wiring visible as an extension form on blue harmonic modes.

The same two-colored incidence data therefore support two honest, different quotients: a real low-conductance quotient and a binary harmonic quotient. Neither determines the original labelled tournaments. The relevant missing coordinates are component incidence, bicycle modes and then the half-energy pairing, not more unlabelled degree statistics.

The strongest immediate research direction is to derive the blue bicycle dimensions and the ranks of the half-energy pairing structurally, rather than guess them from the finite sequence0,0,1,3,11. An all-order exponent-four bound for the **full** network is **not** proposed from n<=7. The blue-only order-eight exponent 16 does not by itself refute a statement about the full network. Likewise these results neither settle all-order black-core connectivity nor construct actual wiring from abstract electrical data.

## Reproduction

From the enclosing workspace:

    python3 session36-metagraph/check_dyadic_network.py --max-n 7
    python3 session36-metagraph/check_dyadic_lift.py
    python3 session36-metagraph/check_bockstein.py
    python3 session36-metagraph/check_tree_scales.py
    python3 session36-metagraph/probe_blue_bicycle.py 8
    python3 session36-metagraph/check_odd_color_weights.py

Repeat each with `python3 -O`. The scripts use Python and SymPy only, import no repository census or producer, and store full finite receipts alongside the manuscript. The separate parent-requested audit is maintained outside this directory.
