# Simplex, sphere, cube: exact geometry and the missing observable

2026-10-20. **Status: PROVED elementary statements and derived consequences; CITED classical geometry; FINITE-EXACT controls.** The packet develops a reusable completion instrument, with sharp hostile examples. It does not claim a new classical Fisher–Rao theorem, Boolean moment theorem, or general convex-duality theorem. Collatz convergence, planar Jacobian, and LRC(14) are not resolved here.

## 1. What is inherited, and what is new here

The distinction that makes the proposed three-shape framework productive is between an **invertible change of coordinates**, an **information-losing observation**, and an **outer relaxation**. They can all look like replacing a polygon by something round, but their logical consequences differ.

### Recovered proved mechanisms

1. **Marked radial sections.** `05-knowledge/results/simplex_sphere_cube_geometry_20261009.md` already proves the three norm-section identifications, polarity facts, homogeneous affine-word matrices, and a conditional height bound from a marked coordinate. Its antecedent is **THM-2688**, `01-canon/theorems/THM-2688-affine-facet-normalization-vertical-sphere-diagonal-simplex-and-lens-quotient.md`. The useful, previously easy-to-drop datum is the scale/owner marker. Repeating the shape changes does not manufacture a height bound.
2. **The exact four-bit Fourier cone.** **THM-3396**, `01-canon/theorems/THM-3396-four-bit-pairwise-independent-fourier-cone.md`, gives a five-coordinate cubic/quartic completion after whitening, its sharp polar polytope, and a convolution law. The sign-symmetric consumer needs just the quartic in the four-root case. `05-knowledge/results/rooted_graph_costs_20261020/README.md` already proves this compression and even-star tomography. The parity example and star coefficient are inherited.
3. **Walsh orthogonality and stochastic coupling.** `05-knowledge/results/prime_local_walsh_bijection_20261020/SUPPORT_COUPLING_ADDENDUM.md` constructs an orthogonal class matrix U and the canonical doubly stochastic coupling B=U∘U. It retains nonzero-support matching and loses signs. Section 8 below identifies exactly what fails if this one-stage correspondence is composed as if squaring were a functor.
4. **Signed symmetry and positive routing.** `05-knowledge/results/sign_aware_augmentation_20261020/signed_root_augmentation.md` retains the automorphism sign character through a nonnegative stabilizer automaton. `05-knowledge/results/bounded_fitting_routing_20261020/BOUNDED_FITTING_ROUTING.md` separates nilpotent routing depth from invertible cycle length. Neither a group order nor a short scalar recurrence is the full transition carrier.
5. **Primitive Fourier transport.** `05-knowledge/results/primitive_hopf_lie_bijection_20261020/PRIMITIVE_HOPF_LIE_BIJECTION.md` detects connected graph primitives by signed cancellation across every nontrivial cut. Its Walsh map is a coalgebra map but not a graph-disjoint-union to tournament-ordered-join algebra map. This is a precise precedent for checking the operation, not only an equality of dimensions.

The canonical hostile here is the q=17 PSD correlation that satisfies every local equilateral-triangle test but violates a global Boolean cut inequality. The corrected near miss is “pairwise PSD means there is a joint sign law.” The least-used sidecar is the **complete multiplication table of even characters**: it turns the simplex of actual laws into a simultaneously diagonalizable positive-semidefinite matrix constraint.

The session contribution is the combined instrument below: a minimal-observer theorem; a sharp Fisher-diameter obstruction for incomplete character records; a same-scalar triangle lift; a missing-row hostile even when all scalar moments are present; quantitative q17 separation; and an exact defect formula for composing squared Walsh couplings. The individual tools use classical linear algebra and Fourier analysis; no historical priority is claimed.

Concept board: norm section; probability amplitude; character multiplication; consumer observability; covariance feasibility; signed composition. Anchor: define the proposed three-shape bridge precisely. Niche: recover unused moment/operation sidecars. Wildcard: quantify the loss in Fisher geometry and the Walsh composition defect.

## 2. Four distinct versions of the three shapes

### 2.1 Positive norm sections: an exact reparametrization

For z∈R^d_{≥0} other than zero, set N_p(z)=z/‖z‖_p, for p=1,2,∞. Then

    N_q(N_p(z)) = N_q(z).

The sections are respectively a (d−1)-simplex, a positive spherical orthant, and the union of d positive cube facets. In d=3 they are a filled triangle, a curved spherical triangle, and three square patches. They are not respectively an entire triangle boundary, an entire sphere, and a filled cube. Their zero-coordinate labels agree. Cube-facet ties must retain all maximizing coordinates.

Consequently alternating these exact section maps alone returns the same marked ray. This does not forbid a useful proof in better coordinates; it says that a genuine new estimate must come from a new inequality, a retained boundary datum, or an operation between changes of chart. In the earlier Collatz carrier z=(x,n,1), the last coordinate preserves x and the original source n by ratios. Removing it destroys integer scale and valuation guards.

### 2.2 Convex polarity: a different operation

With a fixed Euclidean pairing and origin,

    ([-a_1,a_1]×⋯×[-a_d,a_d])° = {y: Σ a_i|y_i|≤1},
    {x: xᵀQx≤1}° = {y: yᵀQ⁻¹y≤1}.

A cuboid's polar is a weighted crosspolytope. For d≥3 a cube is not combinatorially self-dual: a cube has 2^d vertices, whereas its polar has 2d vertices. Polarity exchanges vertices with facets. In two dimensions a square's polar is a rotated and scaled square. A centered regular d-simplex is similar to its polar after reflection and rescaling. Only the unit Euclidean ball is literally equal to its polar in the same Euclidean pairing: K=K° implies ‖x‖²≤1 for x∈K, whence K⊆B₂ and polarity gives the reverse inclusion.

These facts and their proofs are inherited from the October 9 note. Cubical subdivision, combinatorial duality, and polar duality are separate operations.

There is also a rigorous sense in which the cube **is** self-dual: as the coordinatewise ordered set [0,1]^d, the involution x↦1−x reverses the order and is its own inverse. On the vertices it is Boolean-lattice complementation, S↦[d]\S. Thus the user's cubical self-duality intuition has a precise order-theoretic realization. It must not be substituted for Euclidean polarity or face-lattice duality. A partial order, a convex body with pairing, and a face lattice are different typed objects.

### 2.3 Ellipsoid approximation: controlled loss

John's theorem provides a maximum-volume inscribed ellipsoid E of a full-dimensional convex body K. After translating its center to zero,

    E⊆K⊆dE;

for centrally symmetric K the factor improves to √d. The simplex and cube give the respective sharp constants. The contact-point certificate is a positive decomposition of identity, with zero weighted barycenter in the general case. [Keith Ball, Lecture 3, Theorem 3.1](https://math.uchicago.edu/~shmuel/AAT-readings/Combinatorial%20Geometry,%20Concentration,%20Real%20Algebraic%20Geometry/ball.pdf).

Replacing K by E preserves only the resulting certified inclusion bounds unless the original facets are retained. The same normalized inscribed ball can occur for very different bodies. This construction is neither the radial bijection nor the probability square-root map.

### 2.4 Probability amplitudes: an exact, nonlinear metric change

Let Ω have N atoms. A probability law is a point of the (N−1)-simplex,

    Δ(Ω)={p_ω≥0: Σp_ω=1}.

The map Φ(p)_ω=2√p_ω is a homeomorphism onto the positive orthant of the sphere of radius 2. In the interior it carries the Fisher–Rao metric

    g_p(u,v)=Σ u_ωv_ω/p_ω,       Σu_ω=Σv_ω=0,

to the induced Euclidean metric, because dΦ_p(u)_ω=u_ω/√p_ω. Its completed distance is

    d_FR(p,q)=2 arccos(Σ√(p_ωq_ω)).                         (1)

Thus 0≤d_FR≤π, with equality π precisely for disjoint supports. The formula and radius convention agree with the [official Manifolds.jl probability-simplex documentation](https://juliamanifolds.github.io/Manifolds.jl/stable/manifolds/probabilitysimplex/). They are also proved directly by the square-root calculation and the great-circle arc, which stays in the nonnegative orthant.

This spherical representation loses no probability information. The later projection to correlations is the lossy step. With three outcomes the simplex is a filled triangle, and its image is a spherical triangle on S², not a circle. With two outcomes its image is a quarter-circle on S¹. Separately, a planar triangle boundary is homeomorphic to a circle, and three noncollinear planar vertices lie on a circumcircle. Those statements use different maps and dimensions. In particular, a three-phase regular triangle on S¹ in a complex trace calculation is not the Fisher image of a three-outcome law.

### 2.5 No unqualified complexity hierarchy

At fixed affine dimension d≥2, simplex, ball and cube all have dimension d. Their extreme-point counts are d+1, infinite, and 2^d; their boundary regularity differs; their polar families differ; their John approximation factors differ. In dimension one, all three are intervals with two extreme points. A ball has more continuous rotational symmetry, but this is a specified symmetry invariant, not a universal notion of “highest dimension” or “greatest mathematical complexity.” A smooth sphere can encode the same finite probability law as a cornered simplex.

## 3. A minimal-observer theorem

This theorem tells us what enrichment a proposed bridge actually needs.

Let Ω be finite. Record expectations of functions o₁,…,o_k, including the constant 1. Let O be their row matrix on Ω. Suppose the desired consumers are functions f₁,…,f_s, with row matrix C.

**Theorem 1 (exact completion rank).** The following are equivalent:

1. Any two probability laws p,q with Op=Oq also satisfy Cp=Cq.
2. Every row of C lies in the row span of O.

The minimum number of additional arbitrary scalar expectations needed to determine all consumers for every probability law is

    rank([O;C]) − rank(O).                                  (2)

If the existing observations are selected distinct Fourier characters, and extra observations are restricted to individual Fourier moments, the minimum is the number of previously unobserved characters that occur with nonzero coefficient in at least one consumer.

**Proof.** Row-span inclusion gives C=AO and therefore sufficiency. If it fails, some h∈ker O has Ch≠0. Because 1 is recorded, Σh=0. Perturb the strictly positive uniform law to p±=uniform±εh, choosing ε>0 sufficiently small. These are genuine probability laws with the same observations and distinct consumers. Each new scalar expectation enlarges the row span by at most one, proving the lower bound in (2). A basis of the target rows modulo the old row span attains it. For a character-coordinate recorder, independence of the character basis proves the final assertion. ∎

This is a global statement over all laws, including a full-dimensional positive interior. On a separately proved restricted support or nonlinear family, fewer measurements might suffice. Its point is a computable, consumer-specific sidecar, rather than always carrying the full law.

**Rooted-cost application.** An unpinned graph-extension cost f(x) is invariant under x↦−x, so only even Fourier characters can occur. On at most three roots, these are constants and pairs. On four roots, the only missing even character is x₁x₂x₃x₄. Hence one quartic observation is necessary and sufficient for all such four-root costs, but it is unnecessary for a particular cost whose quartic coefficient vanishes. The inherited star has

    g₄=5/4−(1/4)Σ_{i<j}x_ix_j+(1/4)x₁x₂x₃x₄.             (3)

This is the target-specific content of the geometry.

## 4. The lost information can have the maximum possible spherical diameter

Work on Ω={±1}^m/{x∼−x}; its N=2^{m−1} atoms are global-sign orbits. Even subsets A⊆[m] index its characters χ_A(x)=∏_{i∈A}x_i. Let χ_A be a nonconstant even character omitted by a recorder that keeps some of the other characters.

For −1≤t≤1 define the law on Ω

    p_t(x)=N⁻¹(1+tχ_A(x)).                                 (4)

All recorded nonconstant moments vanish; the omitted moment is exactly t. Every p_t therefore lies in the same observation fiber.

**Theorem 2 (maximal invisible arc).** Along (4),

    g_FR = dt²/(1−t²),
    d_FR(p_s,p_t)=|arcsin(s)−arcsin(t)|.                    (5)

In particular p_−1 and p_1 have disjoint support and distance π. Every incomplete even-character recorder has a fiber of the full Fisher–Rao diameter. Interior laws in that fiber can have distances arbitrarily close to π.

For a consumer f=Σ_B f̂(B)χ_B,

    E_{p_t}f=f̂(∅)+t f̂(A).                                (6)

Any estimate using only the recorded data has worst-case error at least |f̂(A)| over this fiber; the midpoint estimate attains that error on this one-parameter family.

**Proof.** Orthogonality gives the moment assertions. Half the atoms have each sign of χ_A, so

    Σ (∂_t p_t)²/p_t = 1/(1−t²).

The square-root affinity of p_s and p_t is

    [√((1+s)(1+t))+√((1−s)(1−t))]/2.

Writing s=sin α and t=sin β, α,β∈[−π/2,π/2], makes this cos((α−β)/2). Equation (1) gives (5), including endpoints by continuity. Equation (6) is Fourier orthogonality. The two endpoint consumer values differ by 2|f̂(A)|, so one estimate cannot be closer than half that difference to both. ∎

For the four-root star, identical pair data can hide costs 1 and 3/2, at Fisher distance π. No algorithm based solely on pairs can uniformly estimate this cost to error less than 1/4. Adding the fourth moment resolves precisely that obstruction. The old parity counterexample supplies the consumer; the new geometric statement says the quotient can erase the largest possible information distance.

There is consequently no dimension-independent fixed Fourier-degree cutoff that determines every globally sign-symmetric consumer as the number of roots grows: choose an omitted even character, or the inherited even-star cost with that nonzero top coefficient.

## 5. A finite simplex is exactly a full character-positive matrix cone

Let μ_A be prescribed real even moments, μ_∅=1. Define the N×N matrix

    M(μ)_{A,B}=μ_{A△B},                                    (7)

with rows and columns indexed by **all** even subsets of [m]. The symmetric-difference index is the multiplication rule χ_Aχ_B=χ_{A△B}.

**Theorem 3 (complete even-character test).** The packet μ is the even-moment packet of a genuine sign law if and only if M(μ) is positive semidefinite. Its unique symmetrized orbit law is

    λ_x=N⁻¹Σ_A μ_Aχ_A(x),                                  (8)

and the eigenvalues of M are exactly Nλ_x. Thus full character positivity is an exact simplex description, not an elliptope relaxation.

**Proof.** Let H_{A,x}=χ_A(x), so HHᵀ=NI. Fourier inversion gives (8), and multiplication gives

    M=H diag(λ) Hᵀ.

The columns of H are orthogonal eigenvectors of M with eigenvalues Nλ_x. Since μ_∅=1 makes Σλ_x=1, PSD is exactly nonnegativity of every atom. Assigning λ_x/2 to x and −x produces a full-cube law. Conversely any law gives these symmetrized masses and the same matrix. ∎

This is the sign-quotient version of the classical finite Boolean moment criterion; compare [Laurent, *Sums of squares, moment matrices and optimization over polynomials*, updated February 6, 2010, §8.2.2, Lemma 8.15](https://homepages.cwi.nl/~monique/files/moment-ima-update-new.pdf). The full test is exponential in m; this is an exactness instrument, not an efficient universal cut algorithm.

### 5.1 The triangle: a minimal repair without new scalar moments

For three signs let a=μ₁₂, b=μ₁₃, c=μ₂₃. The usual 3×3 Gram matrix has unit diagonal and off-diagonal a,b,c. Instead use rows 1,x₁x₂,x₁x₃,x₂x₃:

    M = [[1,a,b,c], [a,1,c,b], [b,c,1,a], [c,b,a,1]].        (9)

Its four eigenvalues are

    1+a+b+c, 1+a−b−c, 1−a+b−c, 1−a−b+c.                  (10)

These are four times the four orbit masses. Their nonnegativity is exactly Boolean realizability of the three pair correlations. No new scalar moments were added: the extra information is the **multiplicative consistency of the same observables**.

The sharp cheap hostile is a=b=c=−1/2. The 3×3 Gram matrix is PSD, with eigenvalues 0,3/2,3/2. But 1+a+b+c=−1/2, so (9) rejects it. Pair PSD alone permits three unit vectors 120 degrees apart, which cannot be three jointly defined signs with those correlations.

### 5.2 All scalar moments are not enough if a character row is omitted

The next hostile separates scalar-coordinate completeness from matrix-test completeness.

**Theorem 4 (one-missing-row obstruction).** For m≥3 set N=2^{m−1} and prescribe

    μ_∅=1,       μ_A=−1/(N−2) for every nonempty even A.     (11)

Every proper principal submatrix of M(μ) is PSD. Nevertheless M(μ) is not PSD and μ is not a law. Its all-plus orbit has mass

    −1/[N(N−2)],

while each of the other N−1 orbits has mass

    (N−1)/[N(N−2)].                                       (12)

**Proof.** Every principal r×r matrix is equicorrelated, with diagonal 1 and off-diagonal c=−1/(N−2). Its eigenvalues are 1−c (multiplicity r−1) and 1+(r−1)c. For r≤N−1 both are nonnegative. For r=N the second is −1/(N−2). Formula (8) gives (12). ∎

For four roots, this is the concrete packet

    all six pairs = −1/6,       quartic = −1/6.

The 7×7 matrix on {1, six pairs} is PSD: eigenvalues 0 and 7/6 six times. The full eight-character matrix has eigenvalue −1/6 and a negative orbit mass −1/48. The quartic already occurs in the 7×7 matrix off-diagonal, so this failure is **not** omission of the quartic scalar. It is omission of the quartic row and its full algebraic consistency test.

The pair data by themselves are realizable. For example, mix equal weights of the uniform law and the uniform law on the six full-cube words with exactly two minus signs. The mixture has all pairs −1/6 but quartic +1/2. Thus this hostile rejects the proposed **full packet**, not its pair projection. It does not assert a new cut-polytope obstruction at four roots.

More generally, taking one orbit mass −ε and the other masses (1+ε)/(N−1), all (N−1)-row tests pass exactly when ε≤1/[N(N−2)]. This is a sharp threshold for this explicit family.

### 5.3 A genuine odd/even alternation in all dimensions

The regular-simplex Gram example generalizes exactly. This is classical complete-graph cut geometry, included as an explicit parity instrument rather than a priority claim.

**Proposition (equicorrelation boundary).** For m≥2 let R have diagonal one and every off-diagonal entry c. Then

    R is PSD  iff  −1/(m−1)≤c≤1.

There is a joint sign law with this correlation matrix if and only if

    −1/(m−1)≤c≤1,   when m is even;
    −1/m≤c≤1,       when m is odd.                         (17)

The law may moreover be chosen exchangeable and globally sign-symmetric.

**Proof.** R has eigenvalues 1−c, of multiplicity m−1, and 1+(m−1)c. For signs, S=Σ_iX_i has the parity of m, so

    E S² = m+m(m−1)c ≥ 0  (m even),
    E S² = m+m(m−1)c ≥ 1  (m odd).

These are the stated lower bounds. Uniformly distributing over all sign words with S=0, or with S=±1 in the odd case, attains them by permutation symmetry. A mixture with the globally sign-symmetric all-equal law attains every intervening c up to one. ∎

At the PSD endpoint c=−1/(m−1), the Gram vectors are the vertices of a centered regular simplex. For odd m their exact vector cancellation has no sign-valued counterpart: an odd number of ±1 values cannot sum to zero. The Boolean boundary lies higher by exactly 1/[m(m−1)]. For even m, balanced signs do attain the PSD endpoint. This is a precise parity alternation across dimensions, with the lost datum identified as integer parity, not a proposed universal ranking of the three body families.

## 6. Quantitative separation of the q=17 sphere-like relaxation

The proved finite q=17 packet supplies a graph with 289 vertices and 2601 norm-one edges. Every Boolean cut has at most 1730 edges. Its explicit radial PSD correlation R* has value −1/3 on every graph edge, predicting a cut value of 1734. Sources:

- `05-knowledge/results/plane_binary_transfer_20261020/proofs/finite_norm_triangle_bridge.md`, §4;
- the stronger bound in `05-knowledge/results/plane_binary_transfer_20261020/addendum-four-defects/README.md`.

Write r for the vector of correlations on these 2601 edges from any genuine joint sign law. Averaging the valid cut inequality gives

    Σ_e r_e ≥ 2601−2·1730 = −859,
    Σ_e R*_e = −867.

Therefore

    Σ_e(r_e−R*_e) ≥ 8.                                    (13)

**Corollary 5 (robust Boolean separation).** Every such r satisfies

    ‖r−R*‖₁ ≥ 8,
    ‖r−R*‖₂ ≥ 8/51,
    ‖r−R*‖∞ ≥ 8/2601.                                    (14)

**Proof.** Apply the triangle inequality and Cauchy–Schwarz to (13), using 2601=51². ∎

These are lower bounds, not proved exact distances. They remain valid if the true maximum cut is smaller than 1730. A numerical relaxation within less than 8/2601 uniformly of the witness's edge correlations is still non-Boolean. This upgrades an exact false positive to a tolerance-aware rejection test.

Here every equilateral triangle has a=b=c=−1/3. The four-character triangle eigenvalues (10) are 0,4/3,4/3,4/3, so the local exact triangle repair passes. Each triangle separately has the uniform law on its six nonmonochromatic assignments. The missing requirement is a common globally realizable law on all 289 vertices. One quartic moment does not universally repair this different problem; Theorem 3 gives an exact all-character completion, while finding a small-support separating lift remains a specific research question.

## 7. The triangle/circle alternation, in a useful precise form

There are now three legitimate arrows:

1. **Law simplex ↔ positive sphere:** p↔2√p, invertible and Fisher-isometric.
2. **Law simplex → pair-correlation polytope:** p↦E[xxᵀ], linear, generally non-injective.
3. **Pair-correlation polytope ⊂ elliptope:** retain only PSD and diagonal-one conditions, an outer relaxation.

To return from the round relaxation to a discrete law, supply either an actual distribution, a complete moment certificate, or a proved sufficient specialized condition. To determine only a consumer, use Theorem 1 and carry its minimal residual moments. To verify realizability, use Theorem 3 and do not drop the multiplication consistency or arbitrary character rows. These are different questions.

The numerical dimensions clarify the issue: an m-root symmetrized law has 2^{m−1}−1 free coordinates, while the pair recorder has at most binom(m,2). At m=3 both counts are three, and the exact simplex inequalities (10) suffice. At m=4 the law has seven coordinates, and the six pairs leave one quartic. Increasing spherical smoothness cannot erase this dimension count.

## 8. The new Walsh coupling needs an operation audit

Let U,V be real orthogonal N×N matrices and B(U)_{ij}=U_{ij}². Each B(U) is doubly stochastic. This is the form used by the inherited tournament/even-graph class coupling.

**Theorem 6 (exact composition defect).**

    B(VU)_{ij}−[B(V)B(U)]_{ij}
      = 2 Σ_{k<l} V_{ik}V_{il}U_{kj}U_{lj}.                (15)

Every row and every column of this defect sums to zero, but the defect need not vanish. If H is a symmetric Walsh matrix of order N=2^r, normalized to U=H/√N, then

    U²=I,       B(U)=J/N,
    B(U²)=I,    B(U)²=J/N.                                (16)

For any initial atom, the two resulting probability laws have total-variation distance 1−1/N.

**Proof.** Expand the square of (VU)_{ij}=Σ_kV_{ik}U_{kj}; the diagonal terms are the stochastic product and the cross terms give (15). Orthogonality gives the stochastic row/column sums, hence the defect sums. Walsh character orthogonality gives H²=NI, while all squared normalized entries are 1/N. The distance between a point mass and the uniform law is 1−1/N. ∎

The useful distinction is between three state types:

- A signed amplitude a on the sphere evolves as a↦Ua. Squaring a_i loses relative signs. For generic nonzero coordinates, restoring a up to global sign needs the relative-sign pattern.
- A diagonal probability state diag(p) evolves coherently to Q=U diag(p)Uᵀ. Taking its diagonal gives B(U)p and discards off-diagonal entries.
- Iterating stochastic kernels B(V)B(U) is valid classical Markov evolution, but it includes this intermediate diagonal projection. It generally differs from the single orthogonal transform VU followed by squaring.

Keeping Q, or the actual signed amplitudes when the state is pure, restores composition. The published one-stage support matching and doubly stochastic coupling remain correct. They never asserted that entrywise squaring respects matrix composition. This identifies the precise sidecar required by a proposed repeated triangle/sphere transform.

## 9. Connections and decisive next tests

| Source → target | Map and preserved predicate | Lost datum / required repair | Cheap decisive test |
|---|---|---|---|
| Marked integer ray → norm section | N_p; affine ratios and marked guards | Unmarked scale; retain last coordinate and arithmetic guards | Existing (3,5) versus (9,15) valuation hostile |
| Probability simplex → spherical orthant | 2√p; entire law and Fisher metric | None when labels retained | Inverse p_i=q_i²/4 |
| Law → prescribed consumers | O; only observed expectations | Quotient row space; add rank gap (2) | Nullspace perturbation of uniform law |
| Pair PSD → three-sign law | Four-character matrix (9) | Character products, not extra scalars | a=b=c=−1/2 |
| Full four-root moments → law | Full eight-character matrix | Dropping one row loses atom positivity | all nonconstant moments −1/6 |
| q17 radial PSD → Boolean field | No valid reverse map at this witness | Global joint law | tolerance bounds (14) |
| Walsh orthogonal map → stochastic coupling | Entrywise square; support, stochasticity | Composition cross terms | Walsh involution (16) |
| Signed symmetry → root-extension count | Signed Burnside / stabilizer transfer | Sign character and actual subgroup action | K_n versus edgeless graph, same S_n but different counts |

For Collatz, the sibling arithmetic lane's cyclic-word carry-difference ideal is a lattice/content datum. Euclidean charts and probability covariance do not automatically preserve it. A legitimate transfer should carry the exact integer ideal through the proposed coordinate map. For planar Jacobian, the sibling lane distinguishes complex bilinear trace data from Hermitian positivity; replacing the former by a PSD object can destroy the algebraic information needed to see poles. Neither application follows from the probability results alone.

Recommended next finite experiments:

1. Compute the Walsh support of an actual new rooted consumer before adding every moment. Use (2) to give the least linear sidecar.
2. On q17, search small-support higher-character consistency inequalities or symmetry-reduced moment completions; demand an exact dual certificate and test against (14). Local triangle completion is already exhausted at the displayed witness.
3. For a proposed repeated Walsh construction, evaluate (15) before identifying its stochastic and signed evolutions. If the defect is nonzero, decide which operation the target actually consumes.
4. For a purported simplex/spheroid/cuboid arithmetic bridge, test integrality, denominator content, owner/scale, and transition guards separately from geometry. An invertible real map is not automatically unimodular.

These tests are bounded, interpretable, and can kill a false bridge before a large census.

## 10. Reproducibility, scope, and review status

`verify_geometric_completion.py` is standard-library-only and makes no repository writes or repository imports. It runs **3,708 exact checks** using rational arithmetic, including independent Fourier inversion versus matrix eigenvector identities.

Declared controls:

- Every omitted nonconstant even character for m=2,…,6, with disjoint parity laws and all surviving moments checked;
- Five deterministic rational laws for each m=2,…,6, testing atom inversion, Gram construction, and every character eigenvector;
- The full 729-point grid a,b,c∈{−1,−3/4,…,1}: 249 Boolean-feasible packets, and 64 further packets passing ordinary 3×3 PSD but not Boolean feasibility;
- Missing-row hostiles for m=3,…,6, and a genuine four-root pair realization with an altered quartic;
- Even-star top coefficients and parity gaps through ten roots;
- Exact Fisher-speed identities at five interior parameters for m=2,…,6;
- Rational Householder composition defects in dimensions 2,…,6, and Walsh involutions of orders 2,…,32;
- Uniform balanced-sign controls for m=2,…,12, checking the exact equicorrelation boundary and parity gap;
- Arithmetic deduction of (14) from the inherited q17 cut bound.

Reproduce from this directory:

    python -B verify_geometric_completion.py
    python -B -O verify_geometric_completion.py

`verification.json` and `verification_optimized.json` agree byte-for-byte. The script does not re-prove the q17 finite certificate or independently rebuild the inherited repository theorems. The all-order claims above have displayed elementary proofs; finite tests are supporting controls. Independent external review of this new packet is still pending at initial handoff. The bounded startup reported a clean main checkout at cached HEAD `0ed0e3aa3ab2`; the parent owns remote freshness and publication. The coordination peek was attempted but unavailable because this checkout has no `.machine-id`; no registration or repository mutation was substituted. No scarce theorem identifier is reserved, and no Git mutation was performed by this worker.
