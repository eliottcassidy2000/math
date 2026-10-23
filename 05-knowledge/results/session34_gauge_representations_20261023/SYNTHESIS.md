# Session34 synthesis: fourth order becomes a map, not a numerology

Started 2026-10-23 03:05 UTC. This packet develops exact finite bridges among representation averaging, tournament Wilson loops, simplicial amplitudes, spectral sectors and a lattice vertex algebra. It does not claim a solution of four-dimensional Yang–Mills existence or mass gap. Standard mechanisms and inherited repository facts are credited; no literature-wide novelty claim is made.

## 1. Burnside in three representation-theoretic steps

For a finite group G acting on a finite set X, take the permutation representation V=C[X]. Its Reynolds operator is

    P=(1/|G|) sum_g rho(g).

P is a projector onto invariant vectors, one dimension per orbit. A permutation matrix has trace equal to its number of fixed basis elements. Therefore

    #orbits = trace(P) = (1/|G|) sum_g #Fix(g).

This is the short proof the user requested. It also identifies what scalar Burnside leaves out: every nontrivial representation channel.

The prime analogy has a precise domain. Over a characteristic-zero field with the appropriate splitting, finite-group representations decompose uniquely as sums of irreducibles. Tensor product instead obeys fusion rules. For S3, standard tensor standard equals trivial plus sign plus standard. For S2, sign tensor sign is trivial. In modular characteristic semisimplicity itself can fail. Irreducibles are additive atoms, not unique multiplicative primes.

The new tournament refinement keeps those channels. Multiplication by the product of all edge signs intertwines a permutation representation with its sign twist and reverses full single-arc-flip adjacency eigenvalues. Thus Specht labels lambda and lambda-transpose pair in precisely stated converse/spectral sectors, not necessarily inside the ordinary unlabelled quotient. The Hamiltonian-marked action is free under vertex relabelling; its complete character and converse decomposition is obtained explicitly. P7 provides a direct control:189 paths, nine automorphism-orbits, converse trace3 and eigenspace dimensions6 and3.

See [representation theorem and certificates](representations/NOTE.md). This is a refinement of counting, not yet a claimed faster all-order enumeration algorithm.

## 2. An exact fourth-order Wilson-loop expression for a tournament observable

A tournament supplies skew edge orientations. A Z2 gauge field instead uses symmetric edge transport signs after choosing a reference orientation; these are not interchangeable conventions. The gauge note fixes that choice and proves its switching, curvature and tetrahedral compatibility laws.

In the128-state fixed-Fano line family, let W_C(T)/W_C(T0) be the normalized Wilson product around a rainbow four-cycle, meaning its four edges lie on four different Fano lines. Exactly21 cycles have this property, three per nonzero fourth-order character. The inherited Fano Hamiltonian count satisfies the explicit identity

    H(T)=693/4+(3/4) sum_(21 rainbow C4) W_C(T)/W_C(T0).

This is a real edge-to-line map, not simply two expressions both called fourth order. It is scoped to the specified family. Globally H is not vertex-switching gauge invariant: a switch takes Paley189 to135 while leaving every Wilson loop unchanged. Averaging H over all vertex switches gives the constant n!/2^(n-1).

Extending the normalized rainbow action to all2^21 edge fields gives a rank-11 rainbow evaluation map, a1024-element ground kernel, and a complete finite partition histogram. Its kernel is the direct sum of the six-dimensional vertex-cut space and four-dimensional embedded Fano code. Ordinary vertex gauge transformations and the additional code symmetry remain distinct.

See [simplicial/Wilson bridge](gauge/SIMPLICIAL_WILSON_FANO.md) and [full ground locus](gauge/RAINBOW_GROUND_LOCUS_ADDENDUM.md).

## 3. Actual simplices and representation labels

For a finite group and a central inverse-symmetric face weight

    w(g)=sum_r d_r a_r chi_r(g),

the normalized tetrahedral-sphere partition function is exactly

    Z=sum_r d_r² a_r^4.

Root gauge fixing and character orthogonality prove this. For S3 with weight z at the identity and1 elsewhere, direct enumeration gives

    Z=(z^4+30z²+80z+105)/216.

Thus representation channels genuinely diagonalize a simplicial gluing computation. This supplies a small exact model for the requested “diagrams on simplices,” without asserting that diagrams, simplices and tournaments are equivalent structures. In a nonabelian model the Bianchi identities also need ordered multiplication and parallel transport; naive unsigned cancellation fails.

See [tetrahedral character addendum](gauge/TETRAHEDRON_CHARACTER_ADDENDUM.md).

## 4. Closed amplitudes can conceal the wrong boundary theory

A stronger hostile is constructed on S3. Two strictly positive central densities have representation coefficients

    (1,+1/4,1/8) and (1,−1/4,1/8).

For every connected closed orientable triangulated surface, they produce identical partition numbers. Indeed the character formula is sum_r d_r^(2−2g)a_r^F, and every closed triangular mesh has even F.

Yet the second model's reflected-boundary convolution operator has a negative eigenvalue−1/4. The first is positive definite and is the time-one operator of an explicit continuous-time Markov convolution. The second has negative determinant and cannot be a real matrix exponential on that six-state space.

So even complete agreement of this entire family of closed simplicial scalar amplitudes does not determine boundary positivity or time refinement. This is a concrete constructive-field-theory warning: retain the boundary operator and its representation signs, not only partition numbers. It is not a statement that every possible reflection or enlarged-state realization is impossible.

See [closed amplitudes and boundary positivity](CLOSED_AMPLITUDES_BOUNDARY_POSITIVITY.md).

## 5. A rigorous finite spectral comparison at fourth order

The first model uses seven Fano-line bits, transverse single-bit flips and a quartic code penalty. Its full128-dimensional Hamiltonian has16 symmetry sectors. The eight-dimensional invariant sector has an exact gap

    sqrt((J−6h)^2+28h²),

which grows like J. The FULL128-state gap instead shrinks as36h³/J² at fixed h>0, J/h→infinity. The slow modes live in sectors discarded by that projection. This does not identify which sector is physical in any Yang–Mills theory; the code symmetry is explicitly declared, not assumed to be physical gauge redundancy.

The second model uses the21 actual edge flips and the full rainbow potential, then projects only by vertex-cut gauge transformations. Its32768-dimensional gauge-invariant space has16 classical ground sectors. Their shortest logical moves have edge distances0,3,6,9 with multiplicities1,7,7,1. The seven distance-three moves are near-perfect matchings, not triangles.

The two models share a crown-graph third-order tunneling mechanism, hence the same cubic gap coefficient. Fourth order separates them:

- Line model, full128-state gap:36h³/J²+192h⁴/J³+O(h⁵/J⁴).
- Edge model, B-invariant32768-state gap:36h³/J²+O(h⁵/J⁴).

For the edge model all four-flip low-space terms are scalar; the independently enumerated518 return contribution and441 energy correction produce the common coefficient−77, which cancels from the gap. There is no weight-four ground-kernel word to create a new off-diagonal logical transition at that order.

This is the session's sharpest fourth-order result: a lower-order agreement is genuine, but it does not identify the dynamics. All formulas are finite-dimensional perturbation theorems with specified parameter limits, not continuum or thermodynamic mass-gap estimates.

See [finite spectral model](FINITE_CODE_SPECTRAL_BRIDGE.md), [logical distances](gauge/LOGICAL_EDGE_DISTANCE_ADDENDUM.md), [edge kinetic theorem](gauge/EDGE_KINETIC_FOURTH_ORDER_ADDENDUM.md), and [fourth-order comparison](FOURTH_ORDER_COMPARISON.md).

## 6. Vertex algebras: an exact obstruction and a positive realization

Linearize the tournament winner rule: e_i e_j is the winner, with e_i²=e_i, and adjoin the independently specified unit. Constant fields define a vertex algebra exactly when this product is associative, which occurs exactly for transitive tournaments. More strongly, its universal associative quotient collapses precisely each strongly connected component. Thus this tempting construction makes locality hold by erasing internal intransitivity.

A separate positive construction keeps coherence data instead of erasing it. On Z2, the3-cocycle omega(a,b,c)=(-1)^(abc) and braiding r(a,b)=i^(ab) satisfy the pentagon and hexagon identities. The same fusion labels as an untwisted Z2 model do not determine these associator/braiding data.

This is realized by the rank-one even lattice vertex algebra V_(sqrt2 Z). It has central charge1, two modules with weights0 and1/4, and twists1 and i. Its theta/eta characters begin

    q^(-1/24),       2q^(5/24).

The normalized character T matrix includes the common factor exp(-pi i/12); the ribbon twist matrix does not. The factor2 in the second leading term is necessary. The primary lattice theorem supports the construction and its semion category; this is not a new VOA or a VOA attached to every tournament.

The same standard lattice construction at rank r has central charge r and vacuum shift−r/24; at rank24 the shift is−1. This is a character/central-charge statement, not a claim that spacetime has24 dimensions or that every such lattice produces the Monster.

This gives an actual arithmetic/topology/vertex-algebra bridge through a lattice discriminant form and roots of unity. It also recovers the earlier1/24 theme as a precise character shift, not a spacetime-dimension analogy.

See [winner/coherence note](vertex/vertex_simplices_associator_bridge.md) and [lattice realization](vertex/LATTICE_VOA_ADDENDUM.md).

## 7. What has and has not advanced toward Yang–Mills

The session supplies explicit finite actions, representation formulas, a boundary-positivity obstruction, and rigorously distinguished low-energy expansions. These are useful test beds and constraints on a proposed proof architecture. They do not construct a nontrivial four-dimensional quantum field theory or prove its vacuum gap.

Fourth Walsh degree, quartic field-polynomial degree, fourth differential order and four spacetime dimensions are different notions. [The continuum-boundary note](FOURTH_ORDER_AND_CONTINUUM_BOUNDARY.md) explains their exact distinctions, classical scaling and quartic-gap hostiles. The [official Clay problem](https://www.claymath.org/millennium/yang-mills-the-maths-gap/) remains open.

A targeted [public constructive-frontier audit](frontier/constructive_ym_frontier.md) recovers the earlier Sussman–Alpoge project thread. Its Theorem0 inputs and current reflection-to-rotation manuscript remain source-qualified, not independently certified here. A concrete next obligation is to identify reflection-compatible Yang–Mills characteristic-functional matrices as positive semidefinite and prove their convergence; positivity passes to an established limit, but finite matrices do not establish that limit or the required uniform analytic bounds.

The attached-paper analyses are retained separately pending approval to publish their qualified critiques. They have not been silently included or replaced by public-source copies. The present packet consists of independently generated/public-source mathematical work.

See [next obligations](NEXT_PROOF_TARGETS.md) and [concept board](CONCEPT_BOARD.md). The main reusable lesson is concrete: representation averaging gives the correct invariant projection only after the physical or combinatorial meaning of that projection is specified, and closed scalar data do not determine the boundary/coherence/positivity structure needed for further gluing.
