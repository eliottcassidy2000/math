# What the order-averaged graph-pair observer forgets

Status: PROVED elementary projection/kernel theorem and an all-size obstruction family; FINITE-EXACT ranks through n=6. Recovered input: THM-4721's signed Walsh basis and the earlier Session58 ordered observer. No novelty priority asserted. The main distinction is essential: **the observer loses linear information from n=6, but still separates every individual converse class through n=6.** No individual-profile collision or all-size separation theorem is claimed.

## 1. Objects and the observer

Let E=binom([n],2), m=|E|, and Omega={0,1}^E. A bit vector x is both a simple graph and the backward-edge encoding of a tournament relative to 0<...<n-1. The ordinary graph action of p in S_n permutes pair coordinates; the tournament action additionally toggles pairs whose endpoints p reverses. Global bit complement c commutes with both actions. Work over Q or R, with uniform inner product 2^{-m} sum_x f(x)g(x).

Let P_G and P_T be the respective ordinary and twisted Reynolds projections (average of pullbacks by all n! permutations), and P_C=(I+c)/2. All are self-adjoint orthogonal projections; P_C commutes with the other two. Let V_G^+ and V_T^+ be their respective invariant spaces with c acting as +1. They are functions on unlabelled graph-complement-pair types and tournament converse classes, respectively.

For a tournament x define A_x(q) as the number of orders producing graph-complement-pair type q. Equivalently,

    A_x(q)/n! = (P_T 1_q)(x),

where 1_q is the ordinary graph-pair orbit indicator. Thus the linear observer on test functions is L=P_T restricted to V_G^+ -> V_T^+; its adjoint is L*=P_G restricted to V_T^+ -> V_G^+. The row profile A_x is the collection of L1_q evaluations. A kernel of L* means lost linear/mixture information, not automatically two equal rows.

## 2. Walsh diagonalization with exact coefficients

For edge subset X define chi_X(x)=(-1)^{|X intersect x|}. For a permutation p set eps(p,X)=(-1)^{#{ij in X:i<j, p(i)>p(j)}. The twisted action sends chi_X to eps(p,X)chi_{pX} (with inverse convention if pullbacks are indexed oppositely; the group average is unchanged).

Call X **Royle-even** if eps(a,X)=+1 for every a in Aut(X). This is NOT the condition that all vertex degrees are even. If X is Royle-even, for Y=pX define e_X(Y)=eps(p,X), independent of the choice of p because the stabilizer character is trivial. Put

    U_X = sum_{Y in orb_G(X)} chi_Y,
    S_X = sum_{Y in orb_G(X)} e_X(Y) chi_Y,
    M_X = |orb_G(X)|,
    alpha_X = (1/M_X) sum_Y e_X(Y)
            = (1/n!) sum_p eps(p,X).

The ordinary sums U_X, one per graph type, are an orthogonal basis of graph invariants. On a graph orbit in Walsh space the twisted invariant subspace is one-dimensional with basis S_X if X is Royle-even, and zero otherwise: an odd stabilizer element forces a coefficient to equal its negative, whereas a trivial stabilizer character propagates one freely chosen coefficient consistently. This recovers the signed-Walsh mechanism of THM-4721, Theorem W; we only need it over Q/R and reprove the relevant part here.

Global complement acts on both orbit sums by (-1)^{|E(X)|}. Consequently retain only even-edge X for V_G^+, and only even-edge Royle-even X for V_T^+. (An even-edge graph need not be Royle-even.) Distinct graph orbits have disjoint Walsh support. Norms satisfy ||U_X||²=||S_X||²=M_X and inner product <U_X,S_X>=M_X alpha_X. Therefore

    P_T U_X = alpha_X S_X           if X is Royle-even,
    P_T U_X = 0                     otherwise,
    P_G S_X = alpha_X U_X.

In particular (P_T P_G)S_X=alpha_X² S_X. The sign of alpha_X depends on the labelled representative gauge, but zero/nonzero, its square, and all rank conclusions do not.

**Projection theorem.** The rank of the full graph-pair observer matrix is exactly the number of graph types X that are Royle-even, have an even number of edges, and satisfy alpha_X != 0. Its adjoint kernel has basis the S_X with those first two properties and alpha_X=0. This is an all-n exact criterion, not a census conjecture.

## 3. Connected obstruction beginning at n=6

For odd r>=3 define F_r with vertices a,b,z_1,...,z_r,w and edges

    ab; az_i and bz_i for every i; aw.

There are r+3 vertices and 2r+2 edges. Degrees of a,b,z_i,w are r+2,r+1,2,1, respectively. Hence Aut(F_r)=S_r permuting only the common neighbors. Label a,b before all z_i, and w after them. Every automorphism preserves every edge's chosen direction in this gauge, so F_r is Royle-even. Its edge count is even, and S_{F_r} is nonzero because its coefficient at its own support is +1.

Interpret alpha_F as the expectation of the product of edge order-signs when vertices receive independent continuous U[0,1] coordinates. Almost surely all coordinates differ and their ranks are a uniform random order. Write x,y for the coordinates of a,b. Up to the fixed global sign convention, integrating each common neighbor and the pendant gives

    alpha_F = integral_[0,1]^2 sign(x-y)(2x-1)(1-2|x-y|)^r dx dy.

The common-neighbor factor follows because sign(x-z)sign(y-z) equals -1 only when z lies strictly between x,y. The pendant's factor is 2x-1. Our bits use the opposite sign per upward edge, but F has an even number of edges so the product is unchanged.

Swap x and y in half the integral. The factor sign(x-y)(2x-1) becomes |x-y| after symmetrization. Using the density 2(1-d) for d=|x-y| gives

    alpha_F = 2 integral_0^1 d(1-d)(1-2d)^r dd = 0.

The last integrand changes sign under d -> 1-d because r is odd. Thus S_{F_r} is a genuine nonzero converse-invariant polynomial in ker(P_G). The restriction r>=3 matters: at r=1 degree classes merge, F_1 has an odd automorphism and its signed tournament invariant vanishes, so that smaller graph is not a valid kernel witness.

For r=3, labelled a=0,b=1,z_i=2,3,4,w=5, the edges are 01,02,03,04,05,12,13,14, i.e. mask 255 in lexicographic pair order. This is a connected six-vertex, eight-edge obstruction with Aut=S_3. Adjoining any number of isolated vertices preserves Royle-evenness, nonzero S_X, even edge count, and alpha=0: the relative order on the six nonisolated vertices remains uniform, and isolated-vertex permutations reverse no edges. Therefore **the observer is not linearly complete for every n>=6.**

### Independent disconnected mechanism

For disjoint components X_i, random relative orders on disjoint vertex sets are independent, so alpha_{disjoint union X_i}=product_i alpha_{X_i}. Reversing the order inside any odd-edge component negates its product, so that component has alpha=0. For example K1,3 is Royle-even with three edges; two disjoint copies give a nonzero converse-even signed invariant killed by ordinary averaging at n=8. The connected F_3 obstruction is earlier and stronger, but the disconnected mechanism explains a broad source of invisible modes.

## 4. Exact finite boundary and individual profiles

The following are independently checked over Q:

| n | dimension of converse-class function space | observer rank | distinct individual profiles |
|---|---:|---:|---:|
| 1 | 1 | 1 | 1 |
| 2 | 1 | 1 | 1 |
| 3 | 2 | 2 | 2 |
| 4 | 3 | 3 | 3 |
| 5 | 10 | 10 | 10 |
| 6 | 34 | 32 | 34 |

At n=6 the two invisible graph types have lexicographically minimal masks 255 and 2015. The first is F_3. The second's ten edges are 01,02,03,04,05,13,14,15,23,24 and its automorphism group has order 2. The verifier records their zero signed-order sums and all 34 basis entries; the infinite-family proof uses only the first graph.

The n=6 original profile matrix has 34 rows and 78 graph-complement-pair columns. Literal orbit enumeration, without relying on Walsh completeness, verifies 34 distinct rows and rational rank 32. A separate Walsh-Hadamard implementation obtains all signed invariant values on the 32,768 labelled tournaments and checks that removing the two zero-alpha coordinates still leaves 34 distinct vectors. Together these checks establish that the first linear failure does NOT itself create an individual profile collision.

The earlier n<=5 packet remains unchanged. Its exact matrices also have full row rank, which follows alternatively from the nonzero-alpha census in this packet. The least possible order of an individual collision is now n>=7; existence is OPEN here. We do not launch an unmotivated larger census.

## 5. A concrete ambiguity for mixtures, not individual tournaments

A nonzero mean-zero f in ker(P_G) does give two distinct probability laws on converse classes with the same averaged observer law. Let C range over converse classes, let b_C be the number of labelled tournaments in C, and choose f=S_{F_3}. Define signed weights w_C=b_C f(C). Orthogonality to constants gives sum_C w_C=0. For every profile column q,

    sum_C w_C A_C(q)
      = n! sum_x f(x)(P_T 1_q)(x)
      = n! sum_x (P_G f)(x)1_q(x) = 0,

using invariance of f and self-adjointness of the projections; the unnormalized sums and uniform inner products differ only by their common positive factor. Let M=sum_{w_C>0}w_C>0. Dividing positive and negative parts by M gives distinct probability laws with identical full graph-pair profile distributions.

At n=6 the certificate has M=122880. Represent each converse class by the smallest backward-edge mask under permutation and global complement, pair order 01,02,03,04,05,12,13,14,15,23,24,25,34,35,45. The positive law is

    {2:1/32, 8:3/32, 10:9/32, 80:1/32, 89:3/16,
     144:3/32, 153:3/32, 154:3/32, 156:3/32}.

The negative law is

    {5:1/32, 9:3/32, 18:3/16, 20:3/32, 72:1/16, 81:3/32,
     87:1/32, 88:1/32, 140:3/32, 145:3/32, 148:3/32, 152:3/32}.

These have disjoint supports, total mass one, and identical profiles in every one of the 78 columns. This is an explicit finite ambiguity of distributions, not an example of two individual tournaments with the same profile.

## 6. Two correlated observations repair the n=6 mixture loss

FINITE-EXACT: retain 32 linearly independent columns of the original profile matrix and append just the two columns

    A_T(3)^2 and A_T(15)^2.

Here graph-pair key 3 is the two-edge star K1,2 plus three isolated vertices, paired with its complement; key 15 is K1,4 plus one isolated vertex, paired with its complement. The exact selected 34-by-34 integer minor is stored in quadratic_certificate.json. Its Bareiss determinant is

    211991592542375584355844096000000,

which reduces to 967372386 modulo the prime 1000000007. Nonzero determinant certifies full rational rank 34; no numerical tolerance or floating rank is used. The selected linear columns and every minor entry are retained, so this is a replayable certificate rather than an unpinned rank report. All degree-two cross-products are unnecessary for this order-six result.

The statistical meaning requires shared provenance. Choose ONE latent converse class T from an unknown law mu, then choose two observer orders independently conditional on that same T. The probability both observed graph pairs equal q is

    sum_T mu(T) [A_T(q)/6!]^2.

Together with the one-observer marginal these two same-latent probabilities for q=3,15 identify every law on the 34 order-six converse classes. By contrast, choosing a fresh independent latent tournament for each order yields [sum_T mu(T) A_T(q)/6!]^2. That is merely the square of an already known marginal and cannot distinguish the two Section 5 laws. Conditional independence of observer orders does not mean independence of their shared latent tournament.

For the explicit positive-minus-negative laws in Section 5, the differences of unnormalized squared-count expectations are -9 at q=3 and 9/4 at q=15. Dividing by (6!)^2 gives differences -1/57600 and 1/230400 in the corresponding double-observation probabilities. These provide concrete separating consequences of the rank repair.

### General finite hierarchy, without an all-n quadratic assertion

For any fixed finite set of N latent classes, let v_T be their profile vectors. If two classes have equal vectors, no polynomial or arbitrary function of that vector separates them. Conversely, if the vectors are all distinct, choose a linear functional ell separating them. Such a functional exists outside finitely many hyperplanes; for integer nonnegative counts bounded by B-1, weights 1,B,B^2,... explicitly encode distinct vectors into distinct integers. Then the Lagrange polynomials

    p_T(v) = product_{U != T} [ell(v)-ell(v_U)]/[ell(v_T)-ell(v_U)]

have degree at most N-1 and equal the class indicators on this finite set. Hence same-latent observer moments through order N-1 recover every mixture law. Each monomial in normalized profile coordinates is a joint outcome probability of independent observer orders conditional on the SAME latent class. Repeating independent latent draws does not supply these moments. This is a finite interpolation guarantee, not a computational efficiency bound or a claim that degree two works at every n.

## 7. Recovery and exact relation to prior work

- THM-4721, slug `twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis`, Theorem W: signed Walsh graph-orbit basis of tournament invariants. We use the explicit monomial-action proof, not unrelated pre-audit subgroup extensions.
- Session58 first packet: fixed-order graph/tournament correspondence, smallest n=3 failure of unlabelled descent, and order-averaged graph-pair profile. This continuation studies the linear information loss of that specific profile, not a different graph encoding.
- Session56 labelled score/forest work: projected orientation cubes, complement-paired score fibers, and the warning that Eulerian-degree graphs are not Royle-even graphs. The score observer keeps labelled degree sums, whereas the present observer keeps a distribution of ordinary graph-pair types over orders. Neither observer is identified with the other.
- Session56 Cayley/forest source `HAMILTONIAN_FOREST_SIZE_BIAS.md`, Section 3: score-energy/Walsh-chaos covariance is an explicit inherited use of Fourier invariants. The common mechanism is retaining signed modes before projecting; its covariance formula does not imply any completeness theorem for this new observer.

The new conclusion relative to these recovered ingredients is the exact Reynolds coefficient alpha_X, the connected F_r kernel family, and the certified n=6 distinction between rank and row separation. Classical Walsh/representation principles and even-graph enumeration are not claimed new.

## 8. Replay and scope

Run analyze.py 6, profile_collision.py, direct_check.py, family_check.py, and quadratic_check.py under ordinary Python and python -O. They use explicit exceptions. analyze.py exhausts graph orbits and permutation signs; profile_collision.py uses independent Boolean Walsh transforms; direct_check.py constructs literal graph/tournament-converse orbits and the complete profile matrix, performs rational elimination, and verifies every mixture column. quadratic_check.py stores and checks the exact full-rank minor for the two extra squared features. family_check.py checks the integral polynomial exactly for a range of odd r and its nonzero even-r controls; the all-size proof is the symmetry argument above, not those finite integrations.

All output is confined to this task directory. No remote or repository edits are made. Source pins and artifact hashes accompany this file.
