# Hamilton-degree parity survives canonical tree transport

Status: elementary all-order identities; FINITE-EXACT sharp optimization at labelled tournament order four. No literature-priority claim. No use of the attached preprint's unverified logarithmic theorem.

## 1. Exact local excess, not just cut domination

Let G be a finite connected loopless undirected graph, with unit original edge capacities, and let T be a spanning tree. Give each tree edge f its fundamental-cut capacity L_T(f). For v in V(G), let R_T(v) be the number of original edges not incident to v whose tree path has v as an internal vertex. Finite path-incidence counting gives

    sum_(f in T incident to v) L_T(f) = deg_G(v)+2R_T(v).

An original edge incident to v contributes once, and every other original edge contributes either zero or twice. This refines Hajiaghayi arXiv:2611.05189v1, Propositions9.3/9.5, whose canonical capacities dominate cuts. It is an independent elementary identity.

For the actual labelled tournament path-complement metagraph at order n>=3, degree is the ordinary Hamilton-path count H(v), always odd by the inherited tournament parity theorem. Therefore every canonical tree weighted degree is odd, dominates H(v), and differs from it by an even transit count. This does not recover H from tree topology or capacities alone: the transit sidecar is required.

## 2. The odd-capacity forest and its exact optimization

Assume every vertex of G has odd degree. Handshaking shows that N=|V(G)| is even. For every A subset V,

    |delta_G(A)| = sum_(v in A) deg_G(v)-2|E(G[A])| = |A| (mod 2).

Consequently a tree edge f has odd L_T(f) if and only if either component of T-f has odd cardinality. Let J_T be the edges with odd L_T. The identity in section1 shows that every vertex has odd degree in J_T. Thus J_T is a spanning odd-degree forest, including every vertex, with at least N/2 edges. Every independent vertex set I supplies the stronger lower bound |J_T|>=|I|: each vertex of I requires at least one selected incident edge, and an edge cannot serve two vertices of I.

Call an edge set J with all vertex degrees odd a V-join. The minimum size tau(G) exists since E(G) itself is one. Then

    min_(T spanning tree of G) |J_T| = tau(G).

Proof. Every J_T is a V-join, giving >=. Choose a minimum-cardinality V-join J. It has no cycle, since deleting any cycle preserves every degree parity and decreases cardinality. Extend the forest J to a spanning tree T. The tree has a unique V-join: the symmetric difference of two would be an even-degree subgraph of a tree; a nonempty forest has a degree-one vertex and cannot be even-degree everywhere. Since J and J_T are both V-joins in T, they coincide. This proves <=. In particular the minimum is N/2 exactly when G has a perfect matching. These are parity/route identities, not expected stretch or congestion estimates.

## 3. Sharp native order-four answer: forty, not thirty-two

Define the labelled order-four metagraph directly. A vertex is a six-bit tournament code, with bit k for lexicographic pair i<j equal to1 if i beats j. For each directed Hamilton word P, preserve its three arcs and reverse the other three. This is the actual K_P move, not a single-arc reversal. The independent enumerator in check.py gives the following exact finite facts:

- There are64 vertices:24 with H=1,16 with H=3,24 with H=5.
- The40 vertices with H in {1,3} form an independent set.
- Each of the24 H=1 leaves has a distinct H=5 neighbour, covering all24 H=5 vertices.
- The16 H=3 vertices can be paired through eight distinct H=5 intermediate vertices. The emitted eight triples explicitly certify this; no cube-isomorphism recognition is required.

The independent set forces every V-join to contain at least40 edges. Include the24 leaf edges and the16 edges in the eight chosen two-edge paths. Every H=1 and H=3 vertex has selected degree1. Eight H=5 vertices have selected degree3, while the other sixteen H=5 vertices have selected degree1. The resulting40-edge set is a forest: its components are either a leaf-edge pair or a three-leaf star centred at an H=5 vertex. It is therefore a V-join that extends to a spanning tree. The checker supplies one63-edge extension and recomputes every fundamental cut directly, without using the parity formula to compute its capacity. Exactly40 of those63 capacities are odd.

Thus the exact native invariant is

    min_T #{f in T : L_T(f) is odd} = 40.

The graph also has maximum matching size exactly24: the24 H=5 vertices form a vertex cover, so any matching has at most24 edges, and the24 distinct leaf edges attain this bound. Its unmatched-vertex deficiency is64-2*24=16. This does not assert a general equality between minimum V-join size and N minus maximum matching size.

The general N/2 bound is32 here, so the count-stratum obstruction costs eight additional odd-capacity edges. In particular this native metagraph has no perfect matching despite all of its degrees being odd. At the expectation level, every distribution has expected number of odd-load edges at least40, and a deterministic tree attains40.

## 4. Verification and limits

Run `python3 check.py > check.out`. The verifier builds all64 tournaments from all24 vertex words, checks adjacency symmetry, classifies degree strata, verifies the independent-set lower bound, constructs and verifies the40-edge V-join and63-edge spanning tree, and computes every fundamental cut. Its output is a certificate with the actual tree and pairing triples.

The general parity statements apply to each connected labelled metagraph component; odd-degree handshaking forces each such component to have even size. The concrete order-four spanning-tree certificate proves its connectedness. Ordinary isomorphism classes and converse quotients require their own degree/loop conventions: H/|Aut| and loop deletion must not be silently substituted for the labelled graph. No assertion is made about order-five sharpness, global spectrum completeness, or a logarithmic distortion theorem.
