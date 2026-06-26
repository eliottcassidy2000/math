# Tournament Matrices As A Universal Translation Layer

codex-2026-06-26-S210.  Prompt: explain how tournaments relate to matrices,
then generate many classic matrix results across domains and ask for the
corresponding useful tournament statement.

## Core Dictionary

For a tournament `T` on vertices `[n]`, the matrix language begins with:

- `A`: the `0/1` adjacency matrix, with `A_ij=1` iff `i -> j`.
- `A + A^T = J - I`: the complete-orientation constraint.
- `s = A 1`: the score vector.
- `S = A - A^T`: the skew signed matrix with entries in `{0,+1,-1}`.
- `iS`: a Hermitian matrix, so spectral theory applies honestly.
- `B = 2A - (J-I) = S`: the same signed relation off the diagonal.
- `L = diag(s) - A`: a directed Laplacian.
- `P`: a stochastic normalization of `A` or of a repaired positive kernel.
- `M(T)`: transfer matrices such as the repo's Hamiltonian-path transfer
  matrix, where `tr(M)=H` in odd order.
- `D S D`: switching by a diagonal sign matrix, the natural gauge action.

Thus a tournament is simultaneously a graph, a skew sign matrix, a Hermitian
operator, a Laplacian, a Markov kernel after normalization, a quiver, a
boundary operator source, a social-choice majority matrix, and a finite
comparison oracle.

## Assumption Challenge

The tempting assumption is that tournament matrix vertices are the original
tournament vertices.  For LRC and proof work, this is often false.  The useful
matrix may instead have rows/columns indexed by:

- runners,
- gaps,
- danger arcs,
- wall crossings,
- residues,
- fixed denominators,
- certificates,
- proof obligations,
- quotient fibers,
- cohomology classes,
- hidden coordinates,
- matrix entries or low-rank update directions.

The preserved predicate should be named each time.  Matrix compression often
preserves a scalar such as rank, trace, determinant, spectrum, or energy while
destroying the exact coordinate that decides the LRC predicate.  The core
guardrail is:

> A matrix invariant is a legal tournament quotient only when its fibers are
> route-pure, status-pure, reconstructible, dual-annihilated, descended by a
> family lemma, or explicitly routed to named residual debt.

## Immediate High-Yield Moves

1. Treat `S=A-A^T` as the primary tournament matrix whenever orientation
   reversal, switching, Pfaffians, determinants, and spectral symmetry matter.
2. Treat `A` as primary whenever reachability, Perron-Frobenius, Hamilton paths,
   Markov flow, or ranking dynamics matter.
3. Treat `L=diag(s)-A` as primary whenever cuts, bottlenecks, arborescences,
   Cheeger-type inequalities, and Hodge/TDA analogies matter.
4. Treat boundary matrices as primary when the obstruction is a cycle, cocycle,
   Cech nerve, path homology class, or owner-essential wall.
5. Treat low-rank updates as primary for edge flips, runner swaps, residual
   tooth repairs, and hidden-coordinate reattachments.

## Atlas Of Classic Matrix Results And Tournament Statements

Each line gives a classical matrix result or object, followed by the tournament
translation that seems worth retaining.

### Linear Algebra And Matrix Analysis

1. **Rank-nullity.**  A tournament quotient matrix should report both rank and
   nullspace.  Tournament statement: a scalar rank certificate is incomplete
   until the null vectors are interpreted as forgotten cycles, cuts, or owner
   coordinates.

2. **Gaussian elimination.**  Row reduction finds independent constraints.
   Tournament statement: eliminate proof obligations, not just vertices; pivots
   identify the first nonredundant hidden coordinate in a repair ladder.

3. **Determinant.**  Determinant measures signed volume and invertibility.
   Tournament statement: `det(I+S)` is a switching-invariant volume of the
   tournament simplex; extremal determinant points toward Hadamard/conference
   geometry rather than Hamilton-path ambiguity.

4. **Trace.**  Trace sums diagonal return data.  Tournament statement: traces
   of powers of `A` count closed directed walks, and traces of transfer
   matrices can encode Hamiltonian-path counts such as the repo's `tr(M)=H`.

5. **Characteristic polynomial.**  The spectrum packages all power traces.
   Tournament statement: the characteristic polynomial of `iS` is a compact
   signed-walk inventory, but it forgets labels and must be checked against
   route purity.

6. **Minimal polynomial.**  The least polynomial relation measures algebraic
   recurrence.  Tournament statement: a tournament with small minimal
   polynomial has few spectral modes; this is a candidate definition of a
   highly regular or dangerously over-compressed proof carrier.

7. **Cayley-Hamilton.**  A matrix satisfies its own characteristic polynomial.
   Tournament statement: high powers of `A`, `S`, or `L` collapse to finitely
   many lower modes; for LRC, this suggests finite recurrence certificates for
   repeated denominator or wall-crossing clocks.

8. **Jordan normal form.**  Non-diagonalizable matrices contain nilpotent
   memory.  Tournament statement: if a quotient operator has Jordan blocks,
   the quotient retained a transient direction; interpret it as delayed
   proof debt rather than equilibrium structure.

9. **Schur decomposition.**  Every complex matrix is unitarily triangular.
   Tournament statement: any tournament operator has an ordered spectral
   cascade; choosing the Schur order is like choosing a Hamiltonian path through
   modes, not vertices.

10. **Spectral theorem.**  Hermitian matrices diagonalize orthogonally.
    Tournament statement: `iS` gives a canonical orthogonal decomposition of
    the signed orientation field; this is the cleanest route from tournaments
    to Fourier/Hoffman-style bounds.

11. **Singular value decomposition.**  SVD separates domain modes, range modes,
    and strengths.  Tournament statement: the SVD of `A`, `S`, or a margin
    matrix splits a tournament into latent ranking, cyclic residual, and noise
    channels.

12. **Polar decomposition.**  `M=UP` separates rotation from stretch.
    Tournament statement: a tournament matrix can be decomposed into a pure
    orientation/phase part and a positive geometry part; LRC quotients should
    state which side preserves loneliness.

13. **QR/LU decompositions.**  Matrices factor into ordered triangular pieces.
    Tournament statement: elimination order is itself a tournament path; bad
    pivots are hidden coordinates surfacing too late.

14. **Cholesky decomposition.**  Positive definite matrices are Gram matrices.
    Tournament statement: whenever a tournament invariant becomes PSD, ask for
    the vectors whose inner products realize it; those vectors may be the
    missing proof witnesses.

15. **Schur complement.**  Eliminating a block leaves a corrected residual.
    Tournament statement: deleting a runner, proof tooth, or denominator must
    add the Schur complement correction; ignoring it is exactly how unsafe
    quotients lose hidden coordinates.

16. **Matrix determinant lemma.**  Low-rank updates have explicit determinant
    changes.  Tournament statement: a single edge flip or residual repair is
    low rank in `S`; determinant and Pfaffian updates should classify which
    flips change the obstruction.

17. **Sherman-Morrison-Woodbury.**  Inverses under low-rank updates are cheap.
    Tournament statement: maintain inverse/resolvent data for active LRC packet
    ledgers as swaps or tail changes are applied, instead of recomputing whole
    certificates.

18. **Cauchy-Binet.**  Minors of products are sums over matching minors.
    Tournament statement: tournament determinant and transfer identities should
    decompose into sums over subtournaments, paths, or owner subsets; this is
    the matrix form of "attach the missing address coordinate."

19. **Jacobi determinant formula.**  `d det(M)=det(M) tr(M^{-1} dM)`.
    Tournament statement: determinant products can be unrolled into derivative
    sums over edge directions; this is a direct bridge to sensitivity and
    hidden-coordinate repair.

20. **Adjugate/cofactor expansion.**  Cofactors locate determinant sensitivity.
    Tournament statement: the cofactor matrix of `I+S` ranks the edges or
    vertices whose flips most affect the signed volume.

21. **Eigenvalue interlacing.**  Principal submatrices have interlacing
    eigenvalues.  Tournament statement: subtournaments inherit constrained
    spectral windows; a forbidden LRC packet might be blocked by showing every
    deletion keeps the wrong interlacing profile.

22. **Weyl inequalities.**  Eigenvalues of sums are controlled by eigenvalues
    of summands.  Tournament statement: split a tournament matrix into AP core
    plus tail perturbation; bound the spectral damage of the tail.

23. **Horn inequalities.**  Possible spectra of sums obey polyhedral rules.
    Tournament statement: score, cyclic, and residual matrices may have
    compatible spectra only inside a Horn polytope; use violations as
    certificate impossibilities.

24. **Gershgorin disks.**  Eigenvalues lie near diagonal centers with row-sum
    radii.  Tournament statement: score imbalance gives crude spectral cages;
    useful as a cheap prefilter before exact tournament enumeration.

25. **Perron-Frobenius.**  Nonnegative irreducible matrices have a positive
    leading eigenvector.  Tournament statement: strong tournaments have a
    canonical positive centrality vector from `A` or a positive repair kernel;
    this vector is a ranking pressure, not necessarily an LRC certificate.

26. **Collatz-Wielandt.**  The Perron root is a min-max ratio.  Tournament
    statement: tournament ranking pressure can be certified by explicit
    positive test vectors; in LRC, such vectors may serve as dual witnesses for
    safe-mass bounds.

27. **Frobenius normal form.**  Reducible nonnegative matrices decompose into
    SCC blocks.  Tournament statement: SCCs of a carrier tournament are exactly
    the irreducible blocks of its adjacency operator.

28. **Matrix norms and condition numbers.**  Norms measure amplification.
    Tournament statement: a quotient with high condition number is a dangerous
    proof carrier; small hidden changes can flip route labels.

29. **Eckart-Young low-rank approximation.**  Truncated SVD is best in Frobenius
    norm.  Tournament statement: the best latent-ranking approximation of a
    tournament is low rank, but the residual cyclic matrix is where Condorcet,
    OCF, and LRC-like debt lives.

30. **Nuclear norm relaxation.**  Rank can be relaxed convexly.  Tournament
    statement: approximate transitivity could be found by nuclear-norm
    compression of the skew/margin matrix, with residual cycles flagged.

### Graph Theory And Combinatorics

31. **Adjacency powers count walks.**  `(A^k)_ij` counts directed walks.
    Tournament statement: powers of `A` give a hierarchy from pairwise wins to
    multi-step dominance and can detect hidden cyclic routing.

32. **Spectral moments.**  `tr(A^k)` counts closed walks.  Tournament statement:
    directed cycle inventories are trace moments; OCF is a refined independent
    odd-cycle support rather than only moment data.

33. **Laplacian matrix tree theorem.**  Cofactors of a Laplacian count trees.
    Tournament statement: directed Laplacian cofactors count arborescence
    redundancy; a robust ranking has many arborescences rooted at high-pressure
    vertices.

34. **Cheeger inequalities.**  Laplacian eigenvalues bound expansion.
    Tournament statement: a tournament carrier with small directed cut
    eigenvalue has a bottleneck proof split; LRC residual fibers should be
    searched there first.

35. **Expander mixing lemma.**  Spectral gap controls edge distribution.
    Tournament statement: quasirandom tournaments have uniform pairwise
    orientation against all large cuts; this is the tournament face of
    discrepancy theory.

36. **Hoffman bound.**  Graph eigenvalues bound independence/chromatic number.
    Tournament statement: apply Hoffman to symmetrized conflict graphs,
    Paley-derived graphs, or tournament comparability residues, not blindly to
    `A` itself.

37. **Lovasz theta.**  SDP gives sharp graph sandwich bounds.  Tournament
    statement: theta of cycle/conflict graphs may bound the maximum collection
    of independent obstruction cycles in a tournament.

38. **Ihara zeta determinant formula.**  A graph zeta is a determinant of a
    nonbacktracking operator.  Tournament statement: primitive directed cycles
    in a tournament should have a nonbacktracking determinant generating
    function; derivative coefficients expose cycle-address data.

39. **Hashimoto nonbacktracking matrix.**  Edges become states.  Tournament
    statement: when vertex quotients are too coarse, lift to directed-edge
    states; this is a matrix version of replacing runners by wall crossings.

40. **Tutte matrix.**  Matchings are detected by a skew symbolic matrix.
    Tournament statement: skew tournament matrices already live in this world;
    Pfaffian substructures can detect parity-rich cycle/path supports.

41. **Kasteleyn/Pfaffian method.**  Planar matchings become Pfaffians.
    Tournament statement: any planar or surface-embedded tournament residual
    should be checked for a Pfaffian orientation that turns path/cycle counts
    into determinant data.

42. **Seidel matrices and switching.**  Graph switching is diagonal conjugacy.
    Tournament statement: tournament switching is `S -> D S D`; any invariant
    meant for switching classes must be a function of this conjugacy class.

43. **Two-graphs.**  Switching classes correspond to parity data on triples.
    Tournament statement: triple-cycle parity is the first nontrivial matrix
    shadow of a tournament; it is a natural bridge to 3-state automata.

44. **Equitable partitions.**  Quotient matrices preserve selected eigenvalues.
    Tournament statement: a tournament quotient is safe only when the partition
    is equitable for the predicate-relevant operator; otherwise it hides route
    debt.

45. **Association schemes and Bose-Mesner algebras.**  Structured graphs have
    commuting adjacency bases.  Tournament statement: Paley and cyclotomic
    tournaments should be studied through their scheme algebra, where Fourier
    idempotents are exact.

46. **Strongly regular graph matrix equations.**  `A^2` collapses to a span of
    `I,A,J`.  Tournament statement: highly regular tournament shadows have a
    tiny matrix algebra; deviations from this algebra are useful obstruction
    coordinates.

47. **Hadamard and conference matrices.**  Orthogonal sign matrices saturate
    determinant bounds.  Tournament statement: skew conference matrices are
    determinant-extremal tournament signs; compare them against Hamilton-path
    and LRC-extremal behavior, which may be orthogonal.

48. **Incidence matrices of designs.**  `N N^T` encodes block intersections.
    Tournament statement: if tournament cycles, paths, or danger arcs form a
    design-like incidence system, matrix regularity can replace enumeration.

49. **Landau score sequence theorem.**  Tournament score sequences satisfy
    majorization inequalities.  Tournament statement: scores are row sums of
    `A`; they are necessary but not sufficient, so every score-based matrix
    quotient needs a cyclic residual sidecar.

50. **Redei Hamiltonian path theorem and parity.**  Every tournament has a
    Hamilton path, and the number is odd.  Matrix statement: Hamilton paths are
    permutation products of `A`, computable by transfer/dynamic matrices; parity
    is the mod-2 collapse that survives all tournaments.

51. **Permutation matrix conjugacy.**  Isomorphism is `A -> P A P^T`.
    Tournament statement: canonical tournament enumeration is orbit reduction
    under simultaneous row/column permutation; matrix invariants are filters,
    not complete labels unless proved.

52. **Weisfeiler-Leman refinement.**  Matrix/color refinement iterates local
    neighborhoods.  Tournament statement: WL on directed/margin matrices gives
    a fast pre-canonicalization layer for A000568-style enumeration.

53. **Graph complement matrices.**  Complement swaps `A` with `J-I-A`.
    Tournament statement: reversal is complement for tournament adjacency and
    negation for `S`; self-complementary tournaments are matrix fixed points up
    to permutation.

54. **Matrix product for graph counts.**  Multiplication composes relations.
    Tournament statement: `A^2` counts two-step dominance; compare it with the
    repo's "A squared" theme where second-order relations reveal hidden
    comparison geometry.

### Probability, Statistics, And Random Matrices

55. **Covariance matrices.**  PSD matrices encode second moments.  Tournament
    statement: random tournament observables have covariance matrices over
    edges, cycles, scores, and Hamilton paths; off-diagonal covariance reveals
    interacting proof teeth.

56. **PCA.**  Leading eigenvectors give principal variation.  Tournament
    statement: PCA of a tournament feature bank identifies dominant deformation
    directions, but rare residual classes may live in low-variance coordinates.

57. **Wishart matrices.**  Sample covariance has a known random matrix law.
    Tournament statement: random feature embeddings of tournaments should have
    Wishart baselines; deviations indicate structured classes like Paley,
    transitive, or AP-like objects.

58. **Wigner semicircle law.**  Random symmetric matrices have universal spectra.
    Tournament statement: `iS` for random tournaments is a Wigner-type Hermitian
    matrix with zero diagonal; spectral outliers measure nonrandom ranking or
    cyclic structure.

59. **Tracy-Widom edge laws.**  Extreme eigenvalues fluctuate universally.
    Tournament statement: spectral-extreme tournament tests should compare to
    random-tournament edge laws before declaring structure.

60. **Matrix Chernoff/Bernstein.**  Random matrix sums concentrate.  Tournament
    statement: random edge samples approximate tournament spectral/cut
    quantities; useful for active pairwise ranking and proof-certificate
    sampling.

61. **Martingale matrix inequalities.**  Adaptive matrix sums concentrate.
    Tournament statement: in active tournament querying, each comparison is a
    matrix martingale increment; concentration controls when the cyclic
    residual is real.

62. **Bradley-Terry-Luce models.**  Pairwise probabilities are logistic
    functions of latent utilities.  Tournament statement: decompose a noisy
    tournament into latent transitive matrix plus skew cyclic residual.

63. **Pairwise comparison eigenvector method.**  Positive reciprocal matrices
    produce rankings.  Tournament statement: majority margins can be exponented
    into reciprocal matrices; inconsistency eigenmodes correspond to Condorcet
    cycles.

64. **Confusion matrices.**  Prediction errors form matrices of class flow.
    Tournament statement: a classifier competition can be turned into a
    tournament; off-diagonal asymmetries locate cyclic model dominance.

65. **Kernel matrices.**  PSD kernels encode similarity.  Tournament statement:
    a tournament becomes kernel-friendly only after a symmetrization or feature
    lift; the lost skew part must remain as a separate obstruction channel.

66. **Matrix completion.**  Low-rank matrices can be recovered from samples.
    Tournament statement: sparse pairwise comparisons can recover a near-rank
    latent ranking, but exact tournament cycles are the non-low-rank residue.

67. **Sign rank and discrepancy.**  A sign matrix's rank controls communication
    complexity and discrepancy.  Tournament statement: sign rank of `S` is a
    powerful measure of how far a tournament is from low-dimensional order.

68. **Cut norm regularity.**  Large matrices have structured plus pseudorandom
    decompositions.  Tournament statement: tournament regularity should split
    into a low-complexity ranking kernel plus quasirandom cyclic dust; LRC
    proof carriers must keep the exceptional cells.

### Number Theory, Fourier, And Arithmetic Matrices

69. **Circulant matrices diagonalize by DFT.**  Fourier modes are exact.
    Tournament statement: circulant tournaments such as Paley/cyclotomic
    tournaments are the cleanest spectral testbed; each residue class is a
    frequency channel.

70. **Character tables.**  Finite group convolution operators diagonalize in
    irreducible representations.  Tournament statement: Cayley tournaments
    should be analyzed by group characters; nonabelian cases turn tournament
    spectra into character-ratio data.

71. **Gauss sum matrices.**  Quadratic characters give explicit spectra.
    Tournament statement: Paley tournament eigenvalues are Gauss-sum shadows of
    quadratic residue orientation.

72. **Ramanujan sums.**  Primitive roots of unity aggregate exact-period data.
    Tournament statement: exact-period tournament or LRC clocks can be stored
    as Ramanujan trace profiles rather than raw denominators.

73. **Smith normal form.**  Integer matrices decompose into invariant factors.
    Tournament statement: the Smith form of incidence, boundary, Laplacian, or
    skew matrices exposes p-adic hidden coordinates invisible over the reals.

74. **Hermite normal form and lattices.**  Integer row spaces have canonical
    bases.  Tournament statement: relation lattices among scores, cycles, and
    packet labels should be normalized before comparing families.

75. **LLL reduction.**  Short lattice vectors reveal hidden relations.
    Tournament statement: use LLL on tournament feature lattices to find small
    cycle/certificate identities, especially in LRC denominator packets.

76. **Vandermonde matrices.**  Evaluation at points is invertible for distinct
    nodes.  Tournament statement: if proof witnesses are evaluations over
    clocks, a Vandermonde condition says which clock samples reconstruct the
    hidden polynomial sidecar.

77. **Companion matrices.**  Polynomial roots become eigenvalues.
    Tournament statement: recurrence laws for tournament sequences or LRC
    clock families can be studied by companion matrices; root location becomes
    a stability guardrail.

78. **Sylvester/resultant matrices.**  Common roots are determinant zeros.
    Tournament statement: two proof clocks collide exactly when a resultant
    vanishes; this is a possible formalization of scalar/product collision.

79. **Hecke matrices.**  Arithmetic correspondences act linearly on modular
    spaces.  Tournament statement: prime-index transformations of structured
    tournaments may form Hecke-like operators on isomorphism classes.

80. **Transfer matrices in recurrence sequences.**  Powers encode evolution.
    Tournament statement: automaton/gap-language tournament carriers should
    keep their transition matrices, not only accepted words.

### Algebra, Representation Theory, And Quivers

81. **Representation matrices.**  Group actions become linear operators.
    Tournament statement: automorphisms of a tournament act on edge, cycle,
    path, and boundary spaces; decomposing these representations separates
    symmetric structure from residual debt.

82. **Maschke decomposition.**  Finite group representations split into irreps.
    Tournament statement: highly symmetric tournament packets should be reduced
    into irreducible orbit modules before enumeration.

83. **Schur's lemma.**  Intertwiners between irreps are scalar or zero.
    Tournament statement: if a tournament carrier respects enough symmetry,
    many possible hidden-coordinate maps are forced to be scalar, hence easy
    to certify or reject.

84. **Cartan matrices.**  Root systems encode reflection geometry.
    Tournament statement: local flip systems of tournaments may have Cartan-like
    matrices; mutations then become reflections in proof-coordinate space.

85. **Coxeter matrices.**  Products of reflections have spectral meaning.
    Tournament statement: repeated edge flips or switchings can be studied as a
    Coxeter-style dynamics on tournament chambers.

86. **Quiver adjacency matrices.**  Directed graphs encode representation
    categories.  Tournament statement: a tournament is the complete quiver with
    one orientation per pair; representation type may classify how complicated
    its cycle modules are.

87. **Cluster exchange matrices.**  Skew-symmetric matrices mutate.  Tournament
    statement: `S` is a skew exchange matrix; edge reversal and mutation ideas
    can generate controlled tournament deformations.

88. **Incidence algebras.**  Posets have triangular zeta and Mobius matrices.
    Tournament statement: transitive tournaments are total orders, so their
    matrices are triangular incidence algebras; cyclic tournaments measure
    deviation from Mobius inversion.

89. **Burnside character averaging.**  Orbit counts are traces of permutation
    matrices.  Tournament statement: A000568 enumeration is a matrix trace
    problem over the action of `S_n` on orientation assignments.

90. **Tensor/Kronecker products.**  Product systems multiply spectra.
    Tournament statement: product tournaments or blowups should track which
    invariants tensor cleanly and which cyclic residues create cross terms.

### Topology, Geometry, And Discrete Differential Operators

91. **Boundary matrices.**  Homology is kernel modulo image.
    Tournament statement: directed path homology of tournaments is matrix
    kernel/image data; owner-essential cycles are literal nullspace classes.

92. **Hodge Laplacian.**  Harmonic forms represent homology classes.
    Tournament statement: tournament obstruction cycles can be searched as low
    Hodge-Laplacian modes on path or Cech complexes.

93. **Persistent homology reduction.**  Column reduction pairs births/deaths.
    Tournament statement: edge flips, deletions, or threshold sweeps produce
    persistence barcodes for cyclic ambiguity and LRC danger-arc topology.

94. **Oriented matroid chirotopes.**  Sign matrices encode orientation data.
    Tournament statement: tournament signs are rank-2 shadows of richer
    oriented matroid data; if a proof needs realizability, the tournament alone
    may be too thin.

95. **Grassmann/Plucker matrices.**  Minors satisfy quadratic relations.
    Tournament statement: if tournament certificates come from subspace
    arrangements, Plucker relations are hidden constraints on allowable cycle
    signs.

96. **Gram matrices.**  Distances and angles are inner products.
    Tournament statement: convert geometric tournament rules, such as unit
    distance orientation, into Gram inequalities plus sign thresholds.

97. **Cayley-Menger determinants.**  Distances determine simplex volumes.
    Tournament statement: a geometric tournament on points has determinant
    feasibility constraints; not every abstract tournament is realizable by a
    distance threshold.

98. **Rigidity matrix.**  Infinitesimal motions are kernel vectors.
    Tournament statement: local rigidity of a tournament realization is a
    matrix-rank statement about which edge orientations can change together.

99. **Stress matrix.**  Self-stresses certify rigidity.
    Tournament statement: dual stresses are certificates that a tournament
    orientation pattern cannot be deformed without crossing a wall.

100. **Discrete exterior calculus.**  Incidence and Hodge-star matrices
     discretize differential operators.  Tournament statement: use incidence
     for orientation/cycle data and Hodge stars for weights such as margins,
     safe interval lengths, or certificate capacities.

### Analysis, PDE, And Operator Theory

101. **Matrix exponential.**  `exp(tM)` solves linear flow.
     Tournament statement: continuous-time ranking or proof-flow dynamics on a
     tournament use `exp(tA)`, `exp(tS)`, or `exp(-tL)`; transient behavior can
     reveal bottlenecks hidden by equilibrium.

102. **Lie bracket/commutator.**  `[X,Y]=XY-YX` measures noncommutativity.
     Tournament statement: two tournament carriers commute only when their
     proof coordinates are compatible; the commutator is a matrix form of
     hidden-coordinate debt.

103. **Baker-Campbell-Hausdorff.**  Products of exponentials have nested
     commutator corrections.  Tournament statement: composing quotient moves
     creates higher-order repair terms; ignoring them is unsafe.

104. **Trotter product formula.**  Alternating flows approximate a sum flow.
     Tournament statement: split a hard tournament operator into AP-core and
     tail operators, then alternate their flows to approximate combined
     dynamics.

105. **Toeplitz matrices.**  Constant-diagonal matrices encode translation
     invariance.  Tournament statement: interval and circulant tournaments have
     Toeplitz/circulant structure; exploit this before brute force.

106. **Hankel matrices.**  Moment sequences become constant anti-diagonal
     matrices.  Tournament statement: cycle-count or Hamilton-path moments may
     have Hankel constraints; rank drops reveal finite recurrence.

107. **Fredholm determinants.**  Infinite-dimensional operators have determinant
     zeta functions.  Tournament statement: graphons or infinite tournaments
     can be studied through determinant expansions whose coefficients are
     finite subtournament densities.

108. **Green's matrices.**  Inverses of differential operators are kernels.
     Tournament statement: inverse Laplacians on tournament carriers measure
     effective resistance between proof obligations.

109. **Finite element stiffness matrices.**  PDE energy discretizes as PSD
     matrices.  Tournament statement: turn LRC safe-mass energy into a stiffness
     form on intervals or danger arcs; zero modes become exact boundary debt.

110. **Lax pairs.**  `dL/dt=[B,L]` preserves spectrum.
     Tournament statement: switching or mutation dynamics may be organized as
     nearly isospectral flows, separating gauge motion from real obstruction
     changes.

### Optimization, Algorithms, And Computer Science

111. **Linear programming constraint matrices.**  Feasibility is matrix
     inequality.  Tournament statement: score-sequence realization, feedback
     arc constraints, and LRC packet sidecars can all be posed as polyhedral
     matrix systems.

112. **Farkas lemma.**  Infeasibility has a dual certificate.
     Tournament statement: every claimed impossible tournament/LRC packet should
     aim for a nonnegative dual certificate, not only failed search.

113. **Total unimodularity.**  Certain matrices make LP relaxations integral.
     Tournament statement: identify which tournament packet matrices are TU;
     those subproblems can be solved exactly by LP.

114. **Network matrices and max-flow/min-cut.**  Incidence matrices encode cuts.
     Tournament statement: residual proof ladders can be made into flow
     networks; min-cuts identify first necessary repair coordinates.

115. **Hungarian algorithm.**  Assignment is linear optimization on a matrix.
     Tournament statement: align two tournament structures by minimum-cost
     permutation before declaring a residual mismatch.

116. **Semidefinite programming.**  PSD matrix constraints relax hard discrete
     problems.  Tournament statement: feedback arc, Kemeny ranking, theta
     bounds, and LRC witness search may admit useful SDP relaxations, provided
     lost labels remain attached.

117. **Goemans-Williamson MaxCut.**  SDP vectors round to cuts.  Tournament
     statement: cyclic orientation debt can be relaxed into vectors, with
     rounding producing candidate edge-flip repair sets.

118. **Kemeny/Slater ranking.**  Ranking a tournament minimizes feedback arcs.
     Tournament statement: the optimal ranking is a matrix ordering problem;
     the residual feedback matrix is the cyclic core.

119. **PageRank.**  A stochastic matrix's stationary distribution ranks nodes.
     Tournament statement: ranking by random walks on `A` is useful for pressure
     but blind to exact cyclic obstruction unless paired with cycle sidecars.

120. **Boolean matrix multiplication.**  Relation composition gives reachability.
     Tournament statement: transitive closure of a tournament SCC is trivial
     once strong, so use Boolean products for carrier graphs rather than raw
     complete tournaments.

121. **Warshall/Floyd algorithms.**  Matrix dynamic programs compute closure and
     shortest paths.  Tournament statement: route-purity and repair ladders can
     be DP'd over carrier matrices of proof obligations.

122. **Min-plus algebra.**  Shortest paths are matrix products over a semiring.
     Tournament statement: minimum repair depth or witness cost is a tropical
     matrix problem on the proof-carrier tournament.

123. **Automata transition matrices.**  Regular languages are recognized by
     matrix products.  Tournament statement: automatic tournament/LRC shadow
     languages must keep transition matrices if state multiplicity matters.

124. **Communication matrices.**  Matrix rank lower-bounds protocols.
     Tournament statement: if the sign matrix of a tournament property has high
     rank, no small quotient can decide it without hidden coordinates.

125. **Randomized linear algebra.**  Sketches preserve matrix structure
     probabilistically.  Tournament statement: use sketches for large
     enumeration, but certify final residual fibers exactly.

### Control, Dynamics, Physics, And Information

126. **Kalman controllability matrix.**  Reachable states are span columns.
     Tournament statement: a set of allowed edge flips or proof moves controls
     a tournament family iff its update matrices span the needed obstruction
     directions.

127. **Observability matrix.**  Outputs determine hidden states.
     Tournament statement: a set of scalar invariants is safe iff its
     observability matrix recovers the hidden tournament coordinate on each
     relevant fiber.

128. **Lyapunov equation.**  Stability has a quadratic certificate.
     Tournament statement: ranking or proof-flow dynamics on a tournament can
     be certified stable by a PSD Lyapunov matrix.

129. **Riccati equation.**  Optimal control balances dynamics and cost.
     Tournament statement: active comparison querying is a control problem:
     spend edge observations to minimize posterior cyclic ambiguity.

130. **Unitary matrices.**  Quantum evolution preserves norm.
     Tournament statement: `exp(it iS)=exp(-tS)` or related Hermitian flows give
     spectral probes of signed orientation geometry.

131. **Density matrices.**  PSD trace-one matrices encode mixed states.
     Tournament statement: uncertainty over rankings or tournament classes can
     be represented as a density matrix; measurement corresponds to querying an
     edge or invariant.

132. **Pauli/Clifford matrices.**  Finite sign operators generate structured
     algebras.  Tournament statement: small tournament sign systems may embed
     into Clifford-like algebras, especially around parity/Pfaffian identities.

133. **Transfer matrices in statistical mechanics.**  Partition functions are
     traces of powers.  Tournament statement: families of tournaments generated
     layer by layer can be counted by transfer matrices, with Hamilton paths or
     OCF weights as fugacities.

134. **Scattering matrices.**  Incoming modes map to outgoing modes.  Tournament
     statement: a proof quotient is a scattering device: incoming local
     coordinates leave as route labels, residual debt, or annihilated modes.

135. **Entropy rate of Markov matrices.**  Stochastic matrices carry dynamical
     entropy.  Tournament statement: random walks on a tournament measure
     ambiguity/choice entropy; compare this with Hamilton-path ambiguity `H`.

## Candidate Theorem Schemas

The following are the most promising tournament-facing statements to test or
formalize next.

1. **Schur-complement quotient rule.**  If a tournament carrier quotient deletes
   a coordinate, the legal residual is the Schur complement, not the raw
   principal submatrix.  Violations predict hidden-coordinate collisions.

2. **Low-rank flip sensitivity.**  Every edge flip in `S` is a rank-2 skew
   update.  Determinant, resolvent, and Pfaffian update formulas should rank
   the edge flips that can change `H`, OCF, route, or LRC status.

3. **Observability criterion for safe tournament quotients.**  A list of
   matrix invariants is proof-safe exactly when its observability matrix
   separates all route/status labels inside the target fiber.

4. **Perron pressure plus cyclic residual split.**  Any ranking-like tournament
   matrix should be decomposed into a positive Perron pressure channel and an
   orthogonal skew/cyclic debt channel.  The first ranks; the second proves or
   refutes structural claims.

5. **Smith normal hidden-clock rule.**  When real/complex matrix invariants
   merge two tournament packets, compute Smith normal forms of the integer
   boundary/incidence/skew sidecar.  The missing coordinate is often p-adic.

6. **Path-Hodge owner-essential theorem.**  Boundary matrix ranks and harmonic
   representatives should certify when an AP/GW-like boundary cycle is
   owner-essential and when it is a deformable open packet.

7. **Cauchy-Binet packet expansion.**  Any product/determinant certificate over
   tournament matrices should be expanded into sums over subtournament minors;
   the surviving minor addresses are the proof coordinates.

8. **Circulant Fourier exactness.**  For circulant tournaments, every safe
   quotient should be expressible in Fourier/character coordinates.  Failures
   identify non-circulant residual debt.

9. **Matrix-tree route redundancy.**  Directed arborescence counts of carrier
   tournaments may measure how robust a proof route is under lost local labels.

10. **Tournament regularity with residual sidecars.**  Cut-norm or spectral
    regularity can compress bulk tournament behavior, but must leave a named
    exceptional matrix containing cycles, owners, endpoints, and p-adic clocks.

## Relation To The Current LRC Stack

The recent HYP-3039/HYP-3040/HYP-3041 sequence fits the matrix lesson cleanly.
The hidden-coordinate ledger says scalar quotients are illegal unless they
retain their missing coordinates.  The hidden-statement ledger names proof
obligations.  The AP-tail puncture atlas gives a miniature matrix theorem:

- the coarse owner-strip matrix merged two rows;
- the missing coordinate was `m mod 13`;
- adding the `q=13` puncture bit or reciprocal fixed-point witness made the
  fiber route-pure;
- this is exactly a Schur/observability story in finite form.

The next useful matrix experiment is to build a sidecar matrix for the target
automatic fiber whose columns are candidate hidden coordinates and whose rows
are residual packet pairs.  Then ask for the smallest column set with full
route/status observability.  That would turn the current ledger doctrine into
a finite linear algebra problem.

## Pointers

Local anchors: `05-knowledge/variables/d-tournament-determinant.md`,
`05-knowledge/variables/transfer-matrix.md`,
`05-knowledge/variables/residue-rank.md`,
`05-knowledge/variables/projection-defect.md`,
`07-reflections/a-squared-as-universal-principle.md`,
`07-reflections/address-coordinate-derivative-repair-s665.md`,
`00-navigation/LRC-CARRIER-PULLBACK-INDEX.md`.
