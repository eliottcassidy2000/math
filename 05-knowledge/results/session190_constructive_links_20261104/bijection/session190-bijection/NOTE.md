# Tournament/even-graph correspondence: probability, return response, and the path-complement boundary

Session 190, 2026-11-04. **PROVED elementary continuation results; FINITE-EXACT audits through order six. No global priority claim.**

## 1. Correct the premise before constructing anything

The cited equinumerosity concerns **automorphism-even graphs**, not graphs whose automorphism group has even order, not Eulerian graphs, and not graphs whose automorphisms induce even permutations of their edge sets.

Choose directions on every edge of a graph G. For an automorphism sigma define epsilon_G(sigma)=(-1) to the number of edges whose direction is reversed when mapped by sigma. This is independent of the chosen directions and is a character of Aut(G). Call G automorphism-even if this character is identically +1.

Royle, Praeger, Glasby, Freedman, and Devillers prove equality between unlabelled n-vertex tournaments and these graphs. Their paper explicitly leaves a *natural* type bijection open. This report does not claim to settle that public formulation, or claim the paper's open status has been comprehensively checked through 2026.

Primary source: [Tournaments and even graphs are equinumerous, J. Algebraic Combin. 57 (2023), 515–524](https://doi.org/10.1007/s10801-022-01197-0), especially Definition in §1, Theorems 1.1/3.3/4.4. [Author-deposited arXiv record](https://arxiv.org/abs/2204.01947).

Two smallest sign controls:

* K2: the endpoint swap reverses its one edge, so epsilon=-1, while its permutation of the one-element edge set has sign +1.
* P3: reflection exchanges its two edges and reverses both relative to a path orientation. Thus epsilon=+1, while the edge-permutation sign is -1.

At n=2 there is one tournament type and one automorphism-even graph type (the empty graph), but **two** graphs with even-order automorphism groups. The requested even-order interpretation is therefore false, already at order two. Our independent n=1..6 census gives tournament counts 1,1,2,4,12,56 versus even-order-group graph counts 0,2,4,11,34,148.

## 2. Inheritance: what is not new here

The repository already contains substantially more than the public equinumerosity statement:

* THM-4721, Theorem W: signed even-subgraph Walsh basis, its two orthogonality identities, its dyadic determinant, and symmetry refinements.
* `05-knowledge/results/tournament_even_bijection_20261011c.md`: an explicit all-order SCC/Lyndon/component type codec, with canonical rank choices at irreducible atoms. This is an effective bijection, not a rank-free geometric one.
* THM-4725 and linked continuations: constructive prime-circulant correspondences and precise obstructions to stronger naturality.
* Session34 `session34-representations/NOTE.md`, §§3–5: full-flip eigenvalue grading, converse parity, and free marked-Hamiltonian-path action.
* THM-848: complete Hamiltonian-path mean response on the full arc cube, and the distinction between mean closure and strong lumpability.
* Session33 `session33-metagraph-memory/NOTE.md`: literal fixed-path complement is an involution, whereas class rerandomization produces a different stochastic process.

The new use developed below is a **canonical stochastic type correspondence with an exact class-return interpretation**, an explicit information-loss witness for the actual path-complement operation, and an exact Walsh matrix for that latter operation. These are elementary consequences/continuations of the inherited basis, not a new enumeration proof.

## 3. Signed Walsh coordinates and normalization

Put N=binom(n,2). Label vertices 0,...,n-1 and choose reference orientation i->j for i<j. A tournament is a sign vector x in {+1,-1}^N. For an edge subset F write chi_F(T)=product_(e in F) x_e(T).

For an even graph G choose a representative, and transport its reference edge directions to each labelled copy G'. Because epsilon_G is trivial, the induced sign s(G') is independent of the transporter. Define

    Phi_G(T) = sum_(G' in orb(G)) s(G') chi_(G')(T).

Changing the representative or its starting directions multiplies the entire column Phi_G by a sign. No individual value or its sign is intrinsic without that convention, but its square is intrinsic.

Let r_T=n!/|Aut(T)|, r_G=n!/|Aut(G)|, and pi_T=r_T/2^N. The inherited orthogonality is

    sum_T pi_T Phi_G(T) Phi_H(T) = r_G [G=H].             (1)

The orbit-indicator functions and these signed graph sums are bases of the same invariant function space. Therefore the square matrix

    Q_(T,G) = sqrt(pi_T/r_G) Phi_G(T)                     (2)

is orthogonal. This can also be seen directly by Fourier-transforming the normalized indicator of each tournament orbit. Odd graph supports disappear because a negative stabilizer element cancels their signed orbit sum.

## 4. A canonical probabilistic type correspondence

**Theorem 1.** The matrix

    B_(T,G) = pi_T Phi_G(T)^2/r_G = Q_(T,G)^2             (3)

is independent of all label and orientation conventions and is doubly stochastic. Thus each tournament type determines a probability distribution on even-graph types, and a uniformly random tournament type maps to a uniformly random even-graph type.

*Proof.* Formula (3) is unaffected by column sign changes or relabelling. Orthogonality of Q makes the sum of squared entries in each row and each column one. Nonnegativity is immediate.

Such matrices are called orthostochastic. The term is standard, not new; see [Chterental–Djokovic, On orthostochastic, unistochastic and qustochastic matrices](https://arxiv.org/abs/0708.4051). All facts needed here follow from (1)–(3).

The column of the empty graph is particularly informative:

    B_(T,empty) = pi_T.

It remembers the labelled orbit size, and hence |Aut(T)|. Consequently equality of B rows forces equality of automorphism orders, though it need not force tournament isomorphism.

**Constructive support bijection.** The bipartite support graph B_(T,G)>0 has a perfect matching. Indeed, if a set S of tournament rows has neighbors N(S),

    |S| = sum_(T in S,G in N(S)) B_(T,G) <= |N(S)|.

Hall's condition follows. A deterministic augmenting-path algorithm, using canonical-mask row/column order, gives an explicit bijection whose matched signed correlation is nonzero. It still makes global ordering choices; it is not a local, relabelling-equivariant labelled map, and is not presented as resolving the paper's natural-bijection problem.

More strongly B is a finite rational mixture of these supported bijections. Find a support matching, subtract the minimum positive matched entry times its permutation matrix, and repeat on the common remaining row sum. At least one positive entry disappears at every step. This is an explicit Birkhoff decomposition and is implemented in the packet, not merely asserted to exist.

At n=3, order the tournament types as transitive, cyclic and the graphs as empty, P3. Then

    Phi = [[1,1],[1,-3]],
    B   = [[3/4,1/4],[1/4,3/4]].

So the probabilistic correspondence is three quarters the rank-preserving bijection plus one quarter the crossed bijection. Neither matching is forced by equality of counts alone.

## 5. Exact meaning: full-arc-flip return responses

This section uses the **full arc cube**, not the pinned-path tiling cube. Let A_(T,U) be the number of individual arc flips taking a representative of tournament type T to type U. Retain loops and multiplicity. A has row sum N and satisfies pi_T A_(T,U)=pi_U A_(U,T).

Flipping one coordinate multiplies chi_G by -1 exactly when that coordinate belongs to G, hence

    A Phi_G = (N-2|E(G)|) Phi_G.                         (4)

**Theorem 2.** For every integer m>=0,

    (A^m)_(T,U) = pi_U sum_G
       Phi_G(T) Phi_G(U) (N-2|E(G)|)^m / r_G.             (5)

The power m=0 uses the ordinary convention theta^0=1, including theta=0. For continuous time, where every labelled edge independently flips at rate one,

    P_t(T,U) = pi_U sum_G
       Phi_G(T) Phi_G(U) exp(-2|E(G)|t) / r_G.            (6)

In particular,

    P_t(T,T) = sum_G B_(T,G) exp(-2|E(G)|t),              (7)
    (A^m)_(T,T) = sum_G B_(T,G)(N-2|E(G)|)^m.            (8)

Thus B is a graph-resolved spectral measure of the isomorphism class. Its distribution of the number of graph edges exactly determines every full-flip return probability. Its first degree moment gives half the exit rate; explicitly, 2 sum_G |E(G)| B_(T,G) = N-A_(T,T).

For a discrete uniform single flip replace each eigenvalue in (5) by 1-2|E(G)|/N, when N>0. Formula (6) is a continuous semigroup; it is not obtained by substituting an exponential into the discrete one-step chain without Poissonization.

The full off-diagonal response needs the **signed** products Phi_G(T)Phi_G(U). B alone loses those signs and cannot recover (5) or (6) for general distinct states. This is an exact loss statement, not just a warning that the matrix is many-to-one.

## 6. A sharp loss witness and why a species bijection is the wrong target

Converse reverses every sign coordinate, so

    Phi_G(T^op)=(-1)^|E(G)| Phi_G(T),
    B_(T^op,G)=B_(T,G).                                  (9)

Even the complete probability row therefore loses converse orientation. There is also additional loss at n=4.

Use lexicographic edge order (01,02,03,12,13,23), with bit one meaning the lower endpoint wins. The canonical tournament masks are 0,2,4,5. Their even-graph masks are 0,3,7,30, respectively the empty graph, P3 plus an isolated vertex, K1,3, and C4. The exact B matrix is

    [3/8, 1/2, 0,   1/8]
    [1/8, 0,   1/2, 3/8]
    [3/8, 1/2, 0,   1/8]
    [1/8, 0,   1/2, 3/8].                                (10)

Mask 0 is transitive. Mask 4 is strongly connected with scores 1,1,2,2. Both are self-converse and have trivial automorphism group. They have identical full-flip return laws at every time, yet they have 1 and 5 Hamiltonian paths. Consequently B does **not** determine the actual path-complement degree or the tiling-line response. At n=5 and n=6 the exhaustive census finds that equal B rows are precisely converse pairs; that finite observation is not promoted to a general conjecture.

The distinctions between possible bijection goals are essential:

* Ordinary unlabelled type bijections exist, including the inherited SCC/Lyndon codec and the support matching above.
* A natural isomorphism of labelled species is impossible already at n=2: there are two labelled tournaments but only one labelled even graph.
* More generally, no relabelling-equivariant map from all labelled even graphs to tournaments can exist for n>=2. The empty graph is fixed by every permutation; no tournament is fixed by a transposition.
* Stronger obstructions to a tournament-to-graph orbit-bijective equivariant surjection are separately proved in the inherited repository. None follows merely from the preceding elementary argument, and none rules out all reasonable notions of an unlabelled natural correspondence.

## 7. Actual path-complement incidence in Walsh coordinates

For a directed Hamiltonian path p of T, let T_p keep the path arcs and reverse every other arc. Define the **unnormalized path-complement operator**

    (D f)(T) = sum_(p a directed Hamiltonian path of T) f(T_p).  (11)

Its row sum is H(T), so D is not the full-flip operator A. On labelled tournaments it is symmetric: the same path p is directed in T_p and (T_p)_p=T. The normalized transition D/H has stationary weights proportional to r_T H(T), rather than r_T.

### 7.1 Free incidence and exact quotient weights

The common directed arcs of T and T_p are exactly p. Hence an automorphism preserving both fixes the directed spanning path, and is the identity:

    Aut(T) intersection Aut(T_p) = 1.                    (12)

The map (T,p)->(T,T_p) is injective, since p is recoverable from the common arcs. This is a reformulation of the inherited free marked-path action, not a new free-action theorem.

For D_(T,U) counting paths from a fixed representative T whose complements have type U, the Aut(T)-action on these paths is free. Thus

    D_(T,U)/|Aut(T)| is an integer,
    D_(T,U)/|Aut(T)| = D_(U,T)/|Aut(U)|.                  (13)

These integers count ordered incidence orbits. For n>=3 the diagonal integers are even: exchanging the endpoints acts freely on ordered pair orbits. If a permutation exchanged the pair, it would preserve their common directed path and be the identity, which cannot exchange distinct endpoints. Equivalently, these are the two orientations of the actual fixed-path-complement lines.

### 7.2 Exact labelled Walsh matrix

Let U range over **undirected** spanning Hamiltonian paths of K_n, each taken once. Choose one of its two traversals. Write epsilon_U(F) for the product of the reference-orientation signs along a subset F of its edges. When |F| is even this is traversal-independent.

The indicator that either traversal of U is directed in T is

    I_U(T) = 2^(2-n) sum_(F subset E(U), |F| even)
                         epsilon_U(F) chi_F(T),          (14)

valid for n>=2. This follows by adding the two products of n-1 direction indicators: the odd terms cancel.

Since complementing off U changes chi_H by (-1)^|H\E(U)|, the uniform-labelled inner product gives the exact formula

    <chi_G, D chi_H>
      = 2^(2-n) sum_(U: G symmetric_difference H subset E(U),
                           |G symmetric_difference H| even)
             (-1)^|H\E(U)| epsilon_U(G symmetric_difference H).  (15)

This is a genuine link between the even-graph/Walsh correspondence and the literal path-complement metagraph. It shows:

* Walsh parity is preserved.
* A nonzero matrix entry requires G symmetric_difference H to be a linear forest with an even total number of edges.
* Equality of degree is not required. D mixes Fourier degrees and is not diagonal in the even-graph basis.
* Symmetry is visible in (15), because G\E(U)=H\E(U) for every contributing path.

At n=3, using the even-graph basis (1,Phi_P3), the matrix whose columns express D applied to a basis vector is

    [[ 3/2, -3/2],
     [-1/2, -3/2]].                                     (16)

Its constant column is H=3/2-Phi_P3/2. In contrast the full-flip operator is diagonal with eigenvalues 3,-1. The two operators cannot be identified even at the smallest nontrivial order.

## 8. A topological identity that prevents the sign confusion

Let sigma be an automorphism of any finite simple graph G. Write s_V, s_E, s_C for the signs of its permutations on vertices, edges, and connected components. Let d_1 be the determinant of its induced map on H_1(G;R). Then

    epsilon_G(sigma) = s_V s_E s_C d_1.                 (17)

*Proof.* On oriented edge chains C_1, the action is a signed permutation matrix, so det(sigma|C_1)=s_E epsilon_G(sigma). On C_0 the determinant is s_V, and on H_0 it is s_C. The equivariant exact sequence

    0 -> H_1 -> C_1 -> C_0 -> H_0 -> 0

gives det(C_1)det(H_0)=det(H_1)det(C_0). All factors have values ±1, yielding (17).

This keeps four different characters separate. It does not say that even graphs have even edge permutations or orientation-preserving cycle-space actions individually.

For a tree, H_1 vanishes, so epsilon=s_V s_E. More concretely, if the tree center is a vertex, orient every edge away from it: all automorphisms preserve every edge direction, so the tree is even. If its center is an edge, orient all other edges away from that edge. The only possible reversal is of the central edge. Hence a tree is odd exactly when its two rooted halves at the center edge are isomorphic. In particular every odd-order tree is even; P3 is even while P4 is odd. This is an elementary structural subclass, not a bijection for general tournaments.

## 9. Reproducibility and checked scope

Run from the enclosing scratch workspace:

    python3 session190-bijection/check.py --max-n 6 --out session190-bijection/results.json
    python3 -O session190-bijection/check.py --max-n 6 --out session190-bijection/results_optimized.json
    cmp session190-bijection/results.json session190-bijection/results_optimized.json
    python3 session190-bijection/path_complement.py

The scripts use only the Python standard library, import no repository producers, exhaust permutations without heuristic isomorphism pruning, use exact integers/rationals, and use explicit runtime checks rather than disabled-under-optimization asserts.

`check.py` exhausts every labelled tournament and graph through order six via complete permutation orbits. It checks signed-copy consistency, both stochastic marginals, weighted orthogonality, the full-flip eigenbasis, every full transition entry for powers 0..4, and the entire class-transition matrix at exact independent-flip correlation z=1/3. The latter is compared against a direct sum over all labelled targets, not against another spectral routine. It constructs and verifies a complete rational Birkhoff decomposition at every order. Homological signs are checked using a separately constructed spanning-forest fundamental-cycle basis and exact determinants, for all 3,635 graph automorphisms across the unlabelled graph representatives through order six.

`path_complement.py` directly enumerates paths, checks trivial common stabilizers, integer quotient weights and diagonal parity, and transforms the resulting class operator to the even-graph basis. For n<=4 it independently constructs D chi_H on every labelled tournament and Walsh-transforms it by butterfly sums, comparing every entry with (15). This checks 4,164 labelled Walsh matrix entries. Higher-order class checks do not masquerade as exhaustive full labelled matrix checks.

Observed B ranks for n=1..6 are 1,1,2,2,10,34. The numbers of distinct B rows and distinct edge-degree aggregated return measures are also 1,1,2,2,10,34 in this census. Converse-merged tournament counts are 1,1,2,3,10,34. The n=4 discrepancy is the explicit transitive/strong hostile above.

No efficient general graph-isomorphism algorithm, new asymptotic enumeration method, natural deterministic bijection, or all-order injectivity of these probability/return signatures is claimed. Publication and independent review belong to the coordinating session.
