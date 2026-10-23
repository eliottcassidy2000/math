# Phase, carries, and actual shared incidence

23 October 2026. **PROVED, scoped elementary statements**; **FINITE-EXACT** checks are in `receipt.json`. This note gives typed bridges, not a claim that a tournament is a quantum Hall fluid or that a square-sum graph is a tournament. General bundle classification and Bargmann invariants are classical. No external priority claim is made for the elementary identities below.

## Inheritance and live concepts

Closest mechanism: Session39's `AFFINE_CONNECTIONS.md`, published under `05-knowledge/results/session39_signed_phases_and_prime_cylinders_20261023/marked-cells/`, identifies every equitable singleton/triple tournament partition with retained block data and actual affine F3 transports. Its hostile is a chordless cycle of cyclic cells with a hidden phase. Its missing coordinate is joint alignment, not another marginal count.

The relevant older Sierpinski results are THM-4754, `01-canon/theorems/THM-4754-three-adic-sierpinski-tournament-collatz-arc-flip-steiner-midpoint.md`, and THM-4561, `01-canon/theorems/THM-4561-the-sierpinski-skew-hadamard-tower-in-the-littlewood-problem.md`. The former treats ternary self-similarity and Pascal reduction modulo3; the latter treats a dyadic skew-Hadamard tower. The binary carry lift below is a different observer, not a new proof of either theorem.

The corrected summand authority is THM-2422, `01-canon/theorems/THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md`. Its operation fibres must not be replaced by the old claims in `07-reflections/metagraph-summand-recursive.md`. In particular, distinct-summand closure from {2,3} omits exactly {1,4,6}; this is a hyperedge dependency fact, not three connected components. The old Fibonacci depth reset is a cutoff artifact. The synchronous depth from {1,2} is ceil(log2(n-1)) for n>=2. See the maintained warning atop `07-reflections/summand-graph-fermat-zeckendorf.md`.

Concept board: actual owner transport; cycle versus triangle tests; phase versus magnitude; dyadic carries; square-shell reflection domains.

## 1. An all-size sufficient synchronization rule using shared incidence

Let G be a finite connected simple graph. Every vertex i owns a finite set D_i of the same size. Every oriented edge i->j carries a specified bijection p_ij:D_i->D_j, with p_ji=p_ij^(-1). These must be actual shared-incidence maps, not independently matched histograms. A parallel owner section is a choice x_i in D_i satisfying x_j=p_ij(x_i).

**Theorem 1 (actual-triangle gluing).** If G is chordal and the composite transport around every triangle is the identity, then every closed walk has identity transport. Consequently any one marked owner at a root extends uniquely to a parallel owner section. In fact every root owner extends, giving a parallel frame.

**Proof.** An immediate backtrack cancels. A simple cycle of length at least four has a chord; splitting it at that chord expresses its transport as a product of conjugates of transports of shorter cycles. Induction reduces to triangles, whose transports are identity. Any closed walk reduces to simple cycles and backtracks. Thus transport along two paths with the same endpoints agrees. Transport the marked root owner along any path; path independence proves existence, and edge-by-edge propagation proves uniqueness. The same construction for all root owners gives the frame. This proof works for arbitrary bijections, not just an abelian phase group.

**Exact weaker criterion in the translation case.** If D_i=F3 and all transports are translations, let Z_1(G;F3) be the cycle space and let B_triangle be the span of boundaries of actual support triangles. Triangle tests force all cycle phases to vanish if and only if B_triangle=Z_1. Indeed phases define a linear functional on the cycle space. If the spans differ, choose a nonzero functional on Z_1/B_triangle and extend it to edge space; it gives an edge-phase assignment passing every triangle test but failing a cycle test. The remaining undetected phase dimension is exactly dim(Z_1/B_triangle). Root marking does not remove these untested cycle phases.

For a connected support, beta=e-v+1 independent cycle measurements are sufficient and necessary for unrestricted translation-phase recovery up to vertex gauge, when every measurement is one F3-valued linear cycle functional. There are 3^beta gauge classes, and fewer than beta such measurements have a nontrivial kernel. This is a restricted linear-measurement minimality statement, inherited from the Session39 phase-orbit mechanism, not minimality among arbitrary encodings.

**Hostile.** On a chordless C_g, g>=4, triangle tests are vacuous. Put identity on all but one edge and a nontrivial translation on the remaining edge. No parallel owner exists, despite every triangle test passing. Actual endpoint marking cannot cure inconsistent transport. Conversely, a reflection holonomy can fix one root owner without fixing a full frame: identity holonomy is sufficient, not necessary, for an individual section.

**Application boundary.** Equitable triple-cell tournament partitions supply precisely such actual maps after retaining block polarity, cyclic orientation, cell names, singleton attachments and uniform arcs. The theorem then certifies synchronization of that retained carrier. It does not authenticate Hamiltonian-port wiring unless the port ambiguity really is exhausted by these shared-incidence maps. Inventing chords or filling arbitrary triangles changes the problem and is not permitted.

## 2. Hopf: a precise phase model and a precise obstruction

The Hopf map sends a unit spinor z in C^2 to its rank-one projector P=zz*, equivalently a unit Bloch vector u in S^2 with P=(I+u.sigma)/2. Multiplying z by a common complex phase leaves P fixed. This is the circle-fibre loss that makes Hopf geometry a useful analogy for forgetting a local chart.

For three normalized spinors, the gauge-invariant Bargmann product is

    B(u,v,w)=<z_u,z_v><z_v,z_w><z_w,z_u>
            =Tr(P_u P_v P_w)
            =(1+u.v+v.w+w.u+i u.(v cross w))/4.

The identity follows directly from sigma_a sigma_b=delta_ab I+i epsilon_abc sigma_c. Pairwise transition probabilities are Tr(P_uP_v)=(1+u.v)/2. They do not determine the imaginary, oriented triple term.

**Exact hostile.** Take (u,v,w)=(e_z,e_x,e_y) and replace w by -e_y. All three pairwise probabilities remain 1/2, but B changes from (1+i)/4 to (1-i)/4. Thus independently recorded pairwise magnitudes discard a real shared triple-incidence coordinate. Complex conjugation reverses this phase; it is not automatically tournament converse.

This is a genuine model for why joint cycle products can repair marginal phase loss. It is not an embedding of arbitrary tournament affine transports into Hopf spinors. Such an embedding would require a positive-semidefinite rank-at-most-two Gram matrix satisfying all proposed overlaps, a condition absent from general tournament transports. Reflected affine edges act anti-linearly in the Session39 O(2) representation, not as ordinary U(1) edge phases. Even in an all-positive component, an arbitrary U(1) connection is not automatically a prescribed Hopf/Berry connection.

**Dimension obstruction.** A graph is a one-dimensional CW complex, so H^2(G;Z)=0. Every complex line bundle over G is topologically trivial. Equivalently, vertex choices of nonzero fibre vectors extend along each edge because C* is path connected; there are no two-cells imposing a Chern obstruction. Therefore pulling the Hopf line bundle back along a map G->S^2 gives first Chern class zero. Nontrivial connection holonomy can nevertheless remain on a trivial bundle. The degree-one cycle data used by tournament synchronization must not be called the degree-two Hopf Chern class.

Attaching two-cells can add relations among cycle transports, as in Section1, but a nonzero Chern class then needs additional transition/curvature data; merely drawing triangles does not manufacture it. This distinction also separates free graph cycle homology from the finite critical-group linking pairing used by the electrical construction.

Sources for these classical settings: [Hatcher, Vector Bundles and K-Theory, sections1.2 and3.1](https://pi.math.cornell.edu/~hatcher/VBKT/VB.pdf); [Rabei, Arvind, Simon and Mukunda, Bargmann Invariants and Geometric Phases, Physical Review A60 (1999)](https://arxiv.org/abs/quant-ph/9907080). The displayed Pauli calculation and hostile are derived here, rather than inferred from a visual analogy.

### 2.1 A concrete Hopf-root tournament encoder and its exact obstruction

Fix a Bloch root r. On distinct unit vectors u_i, define i->j when Im B(r,u_i,u_j)>0, and assume no zero values for distinct i,j. This is an intrinsic antisymmetric observable, so it really defines a tournament, rather than forcing a tournament onto tied data. By the formula above it is exactly the sign of det(r,u_i,u_j).

**Proposition (local-transitivity obstruction).** Every tournament realized by this fixed-root encoder has transitive in-neighborhood and transitive out-neighborhood at every vertex.

**Proof.** Project the u_i to the plane orthogonal to r. They are nonzero because otherwise all incident observables tie. Their angular directions are pairwise neither equal nor opposite, again by the no-tie assumption. Outneighbors of i occupy the open half-circle counterclockwise from i. Sorting those directions within this half-circle gives a total order; the difference of any two sorted angles lies between0 and pi, so every arc points forward. The inneighborhood is treated identically. Thus a directed3-cycle inside either neighborhood is impossible.

This immediately excludes a directed triangle plus a universal source or sink, already at order4. It also excludes the Paley tournament P7 with arcs i->j iff j-i belongs to{1,2,4} modulo7: the outneighbors1,2,4 of0 form the directed cycle1->2->4->1. In contrast, C3 and the cyclic C5 with steps{1,2} are realized by three or five equally spaced equatorial Bloch vectors and a north-pole root. Exact integer planar controls in the checker avoid trigonometric rounding.

In this particular encoder, simultaneously complex-conjugating every spinor and the root conjugates B and therefore induces tournament converse. This is an exact operation on the encoded family, not an identification with K->-K on the native electrical critical pairing; the Session39 no-antiisometry obstruction for that different construction remains intact.

This gives a concrete structural separation among C3, C5 and P7 under a declared Hopf-phase representation. It does not explain the missing Hamiltonian counts{7,21}, does not classify all Hopf-derived constructions, and does not prove the converse that every locally transitive tournament is realized. No physical realization is asserted. The preserved predicate is pairwise orientation; lost information includes latitude and overlap magnitudes. The cheapest failure witness is the four-vertex source-over-C3 tournament.

## 3. Sierpinski endpoint collapse has an exact mod8 carry repair

Let D=I+S, where S shifts an integer-valued sequence by one. Repeated adjacent summation at depth N has row polynomial (1+X)^N. Modulo2, for N=2^r, only the two endpoints remain. This is the familiar Pascal/Sierpinski collapse.

**Theorem 2 (five-tap mod8 row).** For every r>=2 and N=2^r,

    (1+X)^N = 1+4X^(N/4)+6X^(N/2)+4X^(3N/4)+X^N   (mod8).

Thus the aggregate at position i modulo8 is exactly

    f_i+4f_(i+N/4)+6f_(i+N/2)+4f_(i+3N/4)+f_(i+N).

**Proof.** The identity at N=4 is the binomial row (1,4,6,4,1). If P(Y)=1+4Y+6Y²+4Y³+Y⁴, then direct multiplication gives P(Y)²=P(Y²) modulo8. Squaring inducts on r.

For arbitrary integer inputs all five sample positions are indispensable in this coordinate-sampling model: changing only one retained input by1 changes the result by its nonzero coefficient modulo8. A delta input at N/4 has both endpoints zero but output4, the cheapest hostile to importing the parity-only collapse into a mod8 guard. For evaluation, the quarter samples only need their parity and the middle sample only needs its residue modulo4; the endpoints need residues modulo8. The two quarter parities occur only through their sum modulo2, so an aggregate sidecar can compress them. No unrestricted information-theoretic five-symbol minimality is claimed.

**General exact support.** For 0<j<2^r,

    v2(binomial(2^r,j))=r-v2(j).

To prove this, write binomial(2^r,j)=(2^r/j)binomial(2^r-1,j-1). Every coefficient of (1+X)^(2^r-1) is odd: modulo2 it is 1+X+...+X^(2^r-1), obtained by dividing 1+X^(2^r) by1+X. Hence modulo2^s the surviving positions, when r>=s-1, are exactly the multiples of 2^(r-s+1), including both endpoints. There are 2^(s-1)+1 such samples. For s>r+1 all positions survive.

**Two firewalls.** First, every row above is symmetric under j->N-j, so it cannot distinguish a signal from its reversal. The mod8 carry guard must accompany, not replace, a converse-sensitive observer. Second, every additive homomorphism C3->Z/(2^s) is zero, since its generator image a obeys3a=0 and3 is a unit modulo2^s. Therefore dyadic additive guards cannot replace the actual F3 cycle-phase sidecar. Nonlinear encodings are not excluded; this is a precise obstruction to a homomorphic identification of the two phase types.

The electrical mod8 guard and the F3 marked-cell phase should therefore be retained as separate typed coordinates. Nothing here claims that the native electrical guard is generated by this Pascal operator; that would require an explicit operator identification.

## 4. Square-sum edges are partial affine reflections (inherited, rechecked)

The shell reflection and birth-count mechanisms are inherited from `05-knowledge/results/collatz_mod6_20260922_w6_summand_closure_square_filter.md`, statements S10 and S16; the full word formula below is their elementary iterated form. That note also corrects the tempting false extrapolation that every even walk moves by a residue other than2 modulo4: the four-step walk1,3,6,10,15 has displacement14. Only a single two-step displacement, a difference of two squares, avoids2 modulo4.

For the undirected square-sum graph G_n on1,...,n, edge {x,y} means x!=y and x+y=s for a square s. Every square shell is a matching, because its partial map is

    R_s(x)=s-x,

with the domain restrictions 1<=x,s-x<=n and x!=s/2. Thus the graph is the union of actual partial involution matchings, not a tournament or an arbitrary three-colour complete graph.

**Theorem 3 (shell-word transport).** For a word s_1,...,s_l of square shells, its formal endpoint is

    R_(s_l)...R_(s_1)(x)=(-1)^l x + sum_(j=1)^l (-1)^(l-j)s_j.

This follows by induction; in particular R_t R_s(x)=x+t-s. An even closed walk forces the alternating shell sum to be zero. An odd closed walk forces twice its starting vertex to equal that alternating sum. All intermediate vertices must still belong to1,...,n, consecutive vertices must be different, and a simple cycle must have no earlier repetitions. These conditions are indispensable: the balanced word(4,9,9,4) from25 formally closes but immediately leaves the positive interval; the word(4,4) from1 is only a backtrack.

This is the same *affine composition grammar* as signed F3 transports, over a different carrier with essential partial domains. Reduction modulo3 forgets interval membership and injectivity, so it does not preserve Hamiltonian feasibility. The square condition is a restriction on allowed shell labels, not a graph-theoretical encoding of primality.

At birth of vertex n, every new edge has the form {n,k²-n}, where n<k²<2n. Hence its exact number is

    floor(sqrt(2n-1))-floor(sqrt(n)).

This counts new distinct-vertex edges only. Connectivity and Hamiltonian availability are global questions; this birth count alone cannot decide either.

## 5. What the bridges actually preserve

- Actual bijections -> cycle transport: preserves joint owner consistency; marginal histograms lose the maps.
- Hopf spinors -> pair probabilities: loses the oriented Bargmann triple; actual complex overlaps repair it, subject to Gram realizability.
- Pascal row -> mod2 endpoints: loses carry data; the mod8 five-tap formula repairs this observer only.
- F3 phases -> additive dyadic probes: necessarily loses all phase; retain a separate F3 sidecar.
- Square-shell word -> affine endpoint: preserves its formal composition but loses interval validity and distinct-vertex constraints; restore the full domain history.

These give useful all-size sufficient rules and decisive obstructions. They do not prove an all-size complete tournament metagraph wiring classification, a Hamiltonian-path existence theorem for every square-sum order, or a physical Hopf realization of a metagraph.

## Reproduction

Run `python3 check_connections.py` and `python3 -O check_connections.py` from this folder. The checker uses exact integers/Fractions only and explicit failure tests, not Python assertions. It checks every binomial coefficient through dyadic depth1024, all signed coordinate-axis Bargmann triples, all7,776 S3 connections on a four-vertex chordal support, all1,296 S3 connections on C4, and55,980 formal square-shell/start instances. These finite controls corroborate the proofs; they are not a native order-eight wiring census.
