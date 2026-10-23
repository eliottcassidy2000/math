# A positive replacement for the failed graph–tournament identification

23 October 2026. **PROVED at every size:** the weighted quotient coupling, its Walsh diagonalization, and its exact spectral gap 8/9 for n>=3. **FINITE-EXACT controls:** every labeled word and permutation through n=5. This is an elementary derivation with no literature-priority claim.

The antecedent is the Session42 correspondence packet: ordinary graph isomorphism and tournament isomorphism are different actions on the same edge-bit cube. Their complementary/converse quotients cannot be identified. Here they instead define an exact probabilistic coupling. The proof arose by combining the shared-cube action with signed graphical inversion averages; it does not identify the native Hamiltonian-path metagraph with this Markov chain.

## Inheritance: the signed basis was already known in this repository

[THM-4721, section 5, Theorem W](https://github.com/eliottcassidy2000/math/blob/5258d52ac95924d524ee816a7dffa7329ddcf095/01-canon/theorems/THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md), using [THM-4531, Theorem G](https://github.com/eliottcassidy2000/math/blob/5258d52ac95924d524ee816a7dffa7329ddcf095/01-canon/theorems/THM-4531-twisted-mallows-sloane-principle-and-no-natural-bijection-for-switching-classes.md), already constructs the signed Walsh orbit basis of tournament invariants. Its basis is indexed by automorphism-even graph types, meaning their automorphisms preserve orientation parity, not that their vertices have even degree. The source proves weighted orthogonality and a complete linear change of basis over Z[1/2]. These results are inherited, not discoveries of this session.

The present extension compares that SIGNED basis with the ordinary UNSIGNED graph-orbit sums. Their overlap is mu(H), and the alternating conditional-expectation operator squares that overlap. Thus the faithful signed linear correspondence survives while naive unsigned object identification and unsigned averaging lose information. The current claims are the explicit coupling spectrum and its sharp all-size gap, not novelty of the Walsh basis itself. Source hashes are retained in SOURCE_PROVENANCE.json.

## 1. Normalize the common probability space first

Fix n>=2, N=binom(n,2), X=F2^N, with uniform probability 2^-N and a fixed reference vertex order. Interpret x as an ordinary graph edge set, or as the backward-edge set of a tournament.

Let P_G be conditional expectation onto ordinary graph-isomorphism classes. It is the average of the ordinary coordinate-permutation action of S_n. Let P_T be conditional expectation onto tournament-isomorphism classes. It is the average of the affine action

    A_g(x)=P_g x+epsilon_g,

where epsilon_g is the inversion mask of the relabeled reference order. Each operator is an orthogonal projection in L2(X). On graph-invariant functions set

    K=P_G P_T P_G.

Equivalently, first sample uniformly inside the current graph class, then sample uniformly inside the tournament class observed at that word, then report its graph class. This is not a deterministic association of the two classes.

Write W_ij for the number of words belonging to graph class i and tournament class j, g_i=sum_j W_ij and t_j=sum_i W_ij. The exact transition matrix is

    K_ik=(1/g_i) sum_j W_ij W_kj/t_j.

It is row-stochastic, has positive diagonal, and is reversible with stationary measure pi_i=g_i/2^N. Moreover

    g_i K_ik = sum_j W_ij W_kj/t_j

is a symmetric positive-semidefinite matrix. Its quadratic form is a sum of squares. Thus all eigenvalues of K lie in [0,1].

There is an equally valid paired version: enlarge each action by x->x+1, which is graph complement on one side and tournament converse on the other. All statements below explicitly cover this version too.

## 2. Exact diagonalization by graph-shaped Walsh sums

For an edge set F, let

    chi_F(x)=(-1)^(sum_(e in F) x_e).

These characters form an orthonormal basis. Both group actions send chi_F to a signed character with an edge set isomorphic to F. Hence the linear span of the ordinary graph orbit of F is invariant under both projections.

For an ordinary graph shape H, let O_H be its orbit of edge subsets and put

    u_H=(1/|O_H|) sum_(F in O_H) chi_F.

The graph-invariant subspace in this orbit span is one-dimensional. Define the signed inversion average of a representative F by

    mu(F)=(1/n!) sum_(pi in S_n) (-1)^inv_F(pi),

where inv_F(pi) counts the graph edges whose endpoints reverse their reference order in pi. Relabeling F can change the sign of mu, but not its square.

**The diagonalization theorem:**

    K u_H = mu(H)^2 u_H.

Proof: the signed permutation representation on O_H has either no invariant vector, or one invariant line. If some stabilizer of F acts with sign -1, it has none and its average is zero. Otherwise choose its invariant vector with coordinates sigma_F'=+1 or -1 and sigma_F=1. Its mean coordinate is mu(F). Projection of the uniform orbit vector onto this signed line, and projection back, multiply by mu twice. Explicitly, with v=(1/|O_H|) sum sigma_F' chi_F', one has P_T u_H=mu v and P_G v=mu u_H. Thus K u_H=mu² u_H. QED.

The sign-zero stabilizer obstruction is the same orientation-sign issue appearing in the repository's automorphism-even graph codec, but an allowed stabilizer does not alone force mu to be nonzero: the signed orbit average can still cancel.

This diagonalization gives the full spectrum, with one eigenvalue per ordinary graph shape. It is not a numerical fit to the finite transition matrices.

## 3. The signed inversion bound, proved by deletion

**Lemma:** for every nonempty finite simple graph F,

    |mu(F)| <= 1/3.

Isolated vertices do not affect mu, because the relative order of the nonisolated vertices remains uniform. Global reversal changes the sign by (-1)^|E(F)|, so mu(F)=0 whenever the number of edges is odd.

For a direct all-size proof, condition on the first vertex v of a random order. Incident edges contribute the deterministic sign

    epsilon_v=(-1)^(number of neighbors of v preceding v in the reference order).

The remaining order is uniform on the remaining vertices, giving

    mu(F)=(1/n) sum_v epsilon_v mu(F-v).

Induct on the vertex count. If deleting every vertex leaves a nonempty graph, the triangle inequality and induction give |mu(F)|<=1/3.

The remaining case is that some vertex is incident to every edge. Then F is a star with d>=1 edges, together with isolates. The center's relative rank among itself and its d neighbors is uniform on 0,...,d, so, up to a fixed reference-order sign,

    mu(F)=(1/(d+1)) sum_(r=0)^d (-1)^r.

This is zero for odd d and has magnitude 1/(d+1)<=1/3 for positive even d. The single-edge base case also has value zero. This proves the lemma.

For the two-edge star P3, the absolute average is exactly 1/3. Adding isolates changes nothing. Thus this bound is attained at every n>=3.

## 4. Exact uniform spectral gap and convergence

The empty graph character is the constant function and has eigenvalue 1. Every other graph shape has eigenvalue at most 1/9 by section 3. The P3 shape supplies eigenvalue exactly 1/9 for every n>=3. Therefore

    spectral gap(K)=8/9,              n>=3.

At n=2 the unpaired operator projects directly to constants and has gap 1; the paired state space consists of a single state, so the mean-zero space is trivial.

The eigenvalue 1 is simple. Together with the positive diagonal and reversibility, this proves that K is irreducible and aperiodic. For every mean-zero graph-class statistic f,

    ||K^t f||_2 <= 9^(-t) ||f||_2,
    Var(K^t f) <= 81^(-t) Var(f).

For a starting graph state i, the total-variation distance from stationarity obeys

    TV(K^t(i,·),pi) <= (1/2) 9^(-t) sqrt(1/pi_i-1).

This follows by the L2 bound applied to the density of a point mass, reversibility, and Cauchy–Schwarz. A fixed spectral gap does not mean every starting state has an n-independent total-variation mixing time: pi_i can be exponentially small.

A preliminary bound obtainable only from adjacent-transposition edge flips and cube Poincare is gap>=4/(n n!). The signed-deletion argument replaces that weak lower bound by the exact constant 8/9.

## 5. Complement and converse do not spoil the theorem

The global complement involution acts on chi_F by (-1)^|F|. Thus paired graph classes retain exactly the even-edge Walsh orbit sums. Odd-edge shapes already have mu=0 by reversal, so pairing discards only zero-eigenvalue directions. The tournament converse projection likewise removes no nonzero coupled mode.

Consequently the paired and unpaired operators have the same nonzero spectrum with multiplicities. In particular both have exact gap 8/9 for n>=3, and P3 survives in the paired space. They generally have different dimensions and different multiplicities of eigenvalue zero.

This is a genuine structural transfer across the complementary graph/tournament-converse relation, but it is probabilistic and spectral. It is not a bijection of their isomorphism classes.

## 6. Fast mixing can coexist with exact information loss

The rank of P_T restricted to graph-invariant functions, and of K, is exactly the number of graph-shaped Walsh modes H with mu(H)!=0. Therefore the dimension of tournament-invariant observables invisible to every graph-class conditional expectation is

    number of tournament classes - rank(K).

This is an exact linear-information statement, not a count of tournaments recovered or lost. The full tournament-invariant dimension is A000568(n); use the corresponding converse-class count in the paired version.

The finite exact ranks are 1,2,3,10 for n=2,3,4,5. In the unpaired version the tournament dimensions are 1,2,4,12, leaving respectively 0,0,1,2 invisible dimensions. In the paired version the dimensions are 1,2,3,10, so none remain through n=5. This finite full-rank observation is not asserted at all sizes. On the graph side, zero-mode dimensions are 1,2,8,24 unpaired and 0,0,3,8 paired in the same range.

At order four an explicit invisible tournament observable is

    F(T)=sum_v product_(w!=v) a_vw,
    a_vw=+1 if v->w, and -1 otherwise.

It is invariant under tournament relabeling and changes sign under converse. Its values are -4 on C3 -> sink, +4 on source -> C3, and 0 on the transitive and strong four-tournaments. Nevertheless its conditional mean is zero on EVERY ordinary graph-isomorphism class in the fixed-order chart. It is the canceled cubic-star Walsh mode. This demonstrates concretely why rapid coupling mixing is not faithful reconstruction.

No enumeration speedup for A000568 follows merely from a spectral gap. Constructing all incidence weights already requires substantial information, and the chain can erase entire observable directions in one step. An algorithmic claim would need an efficient sampler, a quantified target estimator, and a comparison against existing enumeration methods.

## 7. Finite exact verification and limitations

The standard-library checker enumerates every word and every permutation at n=2,3,4,5. For both paired and unpaired versions it independently constructs W and K, checks stochasticity, detailed balance, positive diagonal and irreducibility, then verifies:

- rational LDL certificates for positivity and the sharp 8/9 gap inequality;
- every graph-shaped Walsh vector's exact eigenvector equation;
- exact trace and squared-trace identities;
- every signed inversion average and its predicted eigenvalue.

No floating-point eigenvalue is used. Normal and optimized runs agree byte-for-byte, with explicit failure checks still active under -O.

The chain is on graph/tournament isomorphism quotients of the whole edge cube under a fixed observer convention. It is not the native fixed-Hamiltonian-path metagraph, not a black-core walk, not a regular-tournament-only chain, and not a shared-port-preserving physical move. Restricting to such carriers changes the measure, action and spectrum; the 8/9 gap cannot be transferred without a new proof. Its value here is an exact, all-size positive replacement for a false direct identification.
