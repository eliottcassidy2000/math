# Historian Tangent Index

Purpose: preserve the odd doors this repo has opened, especially speculative
proper nouns, user-suggested investigation angles, and cross-domain theories
that can be used later as quick inspiration. This is not a proof ledger. For
status and details, check `00-navigation/TANGENTS.md`,
`00-navigation/INVESTIGATION-BACKLOG.md`, `05-knowledge/hypotheses/INDEX.md`,
and the named files in `07-reflections/`.

Coverage note: compiled from the existing tangent log, investigation backlog,
hypothesis index, explainer index, and the full `07-reflections/` title set as
of 2026-05-30.

## Core Tournament Geometry

- Isomorphism class graph `G_n`: the "space of all tournaments" after
  relabeling; often treated as the key object for dynamics, metagraph
  geometry, and quotient transport.
- Merged metagraph `G_n/Z_2`: complement symmetry factored out; source of
  spine/ribs/sea language, bucket parity, chromatic-number claims, and
  transport constraints.
- Tiling hypercube: fixed-Hamiltonian-path staircase encodes tournaments as
  binary tilings; most geometric metaphors in the repo pass through this cube.
- Wiggly/waggly layers: Hamming-distance layers of tiling moves; d=1 is
  wiggly, d=m is the old blue/black line layer, and all layers form the
  waggly system.
- Blue/black lines: grid-symmetric vs non-grid-symmetric tiling-complement
  pairs; important because older sessions overloaded "blue/black" with
  class-level SC/NS behavior.
- Spine/ribs/sea: class-level merged edge types, separating SC-SC, SC-NS, and
  NS-NS mass in the merged metagraph.
- Principal line of symmetry: recurring geometric axis in the staircase and
  metagraph pictures; tied to complement/transposition alignment.
- Two staircases / two sigma structures: the pin-grid staircase and tournament
  symmetries have subtly different sigma actions; several false starts came
  from conflating them.
- Two-sheeted cover: complement symmetry repeatedly appears as a double cover,
  with `(1-x)^(-1/2)` and Z/2 motifs.
- Unlabeling principle: many phenomena only become clean after quotienting by
  isomorphism; labeled data can hide the invariant.
- Orbifold `G_n`: tournament space as a finite quotient of a cube by `S_n`;
  useful language for singular fibers and Burnside counts.
- Tournaments at infinity: limit object as `n -> infinity`, with the cube
  becoming Cantor-like and quotients controlled by oligomorphic symmetry.

## H, Omega, and Forbidden Values

- Odd-Cycle Formula / OCF: `H(T)=I(Omega(T),2)` is the central bridge from
  Hamiltonian paths to independent sets of vertex-disjoint odd cycles.
- Conflict graph `Omega(T)`: odd directed cycles are vertices, disjointness is
  adjacency; its independence polynomial is the compressed H-count.
- Hard-core lattice gas at fugacity 2: `I(Omega,2)` is a statistical mechanics
  partition function; suggests transfer matrices, correlation inequalities,
  and polymer expansions.
- 2-adic tower / higher Redei: evaluating independence polynomials at
  `2,4,8,...` suggests higher parity layers beyond Redei's oddness theorem.
- Forbidden H values 7 and 21: permanent gaps; later work refuted 63 as
  universally forbidden.
- H=63 unlock: at n=8, exactly two classes realize H=63 via complete
  `Omega=K31`, not through the old disconnected-obstruction mechanism.
- Single-core complete-Omega signatures: H=63 examples have all odd cycles
  through one core vertex; binary insertion signatures `1001100` and `1100110`
  give weighted count 31.
- H=21 complete-Omega lens: reduce to absence of weighted signature counts
  `r=10` and non-core complete-Omega alternatives.
- Projection-kill / near-kill: H=63 is an exact deletion kill of all odd cycles;
  the n=9 real-root failure is a near-kill with a small surviving residue.
- Omega disjointness geometry: cycle-family core size, alpha-vector, disjoint
  pair count, support multiplicity, and independent-triple shape may explain
  H-gaps and root failures better than H alone.
- Real-rootedness of `I(Omega,x)`: true through n=8, false universally at n=9;
  remains a structural classification problem.
- Ultra-log-concavity / Newton / Turan: real-rootedness led to ULC strategies,
  Turan inequalities, and an unconditional k=1 style proof direction.
- Complete-Omega vs sparse-residue extremes: complete conflict graphs and
  concentrated non-real-root failures are opposite ends of one disjointness
  axis.
- Size-weighted independence polynomial `u_T`: explored as a refinement for
  H=21 and component reduction.
- Lee-Yang zeros / root gaps: independence-polynomial root locations were
  linked to physics-style phase boundaries and root-gap literature.

## Good Cuts, Buckets, and Transport

- Good-cut count `g`: now understood as strong-component defect
  `n - #SCC(T)` for a base path; descends to tournament isomorphism and merged
  classes.
- Good-cut bucket coordinate: early evidence showed merged classes pure in
  good-cut count; later SCC interpretation made this structural.
- No bucket 1 / top bucket gap: tilings cannot have exactly one good cut; top
  bucket behavior connects to strong connectivity.
- Good-cut interval gas: good-cut sets are interval unions on the cut path,
  with a 1D polymer recurrence.
- Bucket Morse coordinate: single-tile flips change interval covers one local
  piece at a time; possible cheap TDA feature.
- Quotient bucket balance: any finite quotient of a Boolean cube has half-line
  conservation and prescribed transport row sums.
- Boolean-cube mask theorem: xor masks are fixed-point-free involutions; Lean
  formalization gives reusable bucket-balance checks.
- Transport row checksum: row sums, even diagonal, cross-outflow parity, and
  Lucas-active layers constrain all quotient transport matrices.
- Endpoint transfer recursion: endpoint insertion gives quotient transfer
  matrices with row/column sums and parity boundary behavior.
- Projection defects: compare `Q_m -> G_n/Z_2` and `Q_m -> E_n`; defects split
  into tournament-only, even-only, joint, and silent.
- Structured projection probes: endpoint stars, strips, range flips, and
  Hamming shells act like different finite-difference operators, not
  interchangeable samples.
- Range parity polarity: even tile ranges tend even-graph biased; odd ranges
  tend tournament-class biased.
- Negative space / commutators: many errors cluster where two natural
  projections almost commute but leave a measurable residue.

## Even Graphs and Dual Quotients

- Even graph metagraph `E_n`: dual quotient of the tiling cube by even graph
  structure; first-class companion to `G_n`.
- Cut/cycle decomposition: tournament orientations and even graphs are two
  shadows of cube/cycle-space data.
- Bridge matrix between tournaments and even graphs: observed full rank in
  computed cases; suggests the two quotients retain complementary information.
- Even-graph projection for odd n: map from tournament orientation to
  degree-even cycle-space graph, separating Paley and interval examples.
- Royle-even / dark classes: graph-side "dark H" attempts overlapped too much;
  a combined cycle-space/order/entropy theory is still missing.
- Projection-defect features: proposed for Tournament TDA as a cheap dual-lens
  fingerprint beyond H, score, and Betti profile.

## Paley, QR, Codes, and Number Theory

- Paley tournaments: central maximizer family built from quadratic residues;
  verified maximizer at several sizes and linked to Gauss sums.
- Paley T_11 taxicab number: `H(T_11)=95095=55*1729`; whether the 1729
  structure is meaningful remains open.
- QR/NQR deletion symmetry: anti-automorphism such as `x -> 2x` in T_11 gives
  equal sub-tournament Hamiltonian-cycle counts.
- Paley deletion optimality: deleting a vertex from Paley appears to preserve
  maximizer behavior in known prime cases.
- Paley vs interval tournaments: same regular score can have different support
  multiplicity, even-graph projection, and alpha dominance regimes.
- QR resonance principle: congruences like `2^C(p,2) mod p^2` and Legendre
  symbols organize Paley orbit parity.
- Gauss sums / L-functions: Paley eigenvalues are direct Gauss-sum expressions,
  inviting Dirichlet L-function analogies.
- Krawtchouk coordinates: tournament tilings are band-limited in a coding
  theory basis; Paley T_7/T_23 resonate with Hamming/Golay codes.
- Paley gives dual codes: quadratic-residue tournament structure mirrors
  classical code dualities.
- Hamming [7,4,3] and Golay [24,12,8]: used as recurring perfect-code
  analogues for forbidden values, QR structure, and sphere packing.
- Fano plane / Steiner triples: T_7, triples, and alpha_2 structure connect to
  design theory.
- Skew-Hadamard / DRT / NDRT: doubly regular and nearly doubly regular
  tournament literature supplies maximizer candidates and families.
- Satake cyclotomic NDRTs: flagged as a possible new tournament family to test.
- A038375 / Paley maximizer sequence: extended computationally; still a
  sequence to preserve.
- `SC(n)=A000568` corrected identity: self-complementary class counting needed
  a q-deformed/counting correction.

## Path Homology, Topology, and TDA

- GLMY path homology: tournament path complexes are a major parallel theory to
  OCF, with Betti landscapes and vanishing questions.
- `beta_2=0`: central vanishing theorem target; multiple proof architectures
  via bad triangles, arc flips, and LES.
- Seesaw principle: adjacent or separated Betti dimensions trade off through
  image saturation; true in low ranges, corrected after n=8 counterexamples.
- Hidden cycles / ghost cycles: homology classes that vanish or survive under
  deletion/projection; tied to support-residue calculus.
- Path-homology projection residues: through-v cycles becoming boundaries after
  projection are the homological version of projection-kill.
- Paley homology eigenspace split: Fourier/circulant decomposition isolates
  trivial and nontrivial Betti sectors.
- Tang-Yau circulant Fourier homology: literature lead for symbol-matrix path
  homology of circulant digraphs.
- Magnitude-path spectral sequence: Hepworth/Roff/Asao leads for connecting
  path homology with magnitude/reachability homology.
- Poincare duality / Omega palindrome: tempting but not universal; T_7 was
  palindromic while T_11 was not.
- Pfaffian-Betti connection: even/odd determinant and path-count structures
  correlate with Betti split phenomena.
- Tournament TDA: proposed feature suite using H, scores, Betti profiles,
  Omega alpha-vector, projection defects, and bucket features.
- Protein folding / ligand binding: path Laplacian and persistent digraph TDA
  were noted as application bridges.

## Transfer Matrices, Walsh, and Algebra

- Transfer matrix symmetry `M[a,b]=M[b,a]`: a major proved identity with many
  routes (Hopf, Irving-Omar, Walsh).
- Irving-Omar det/per formula: external matrix-algebra bridge to endpoint
  conditioned Hamiltonian path counts.
- W-coefficient hierarchy: exact spectral decomposition stratifying odd-cycle
  complexity and transfer behavior.
- Walsh/Fourier OCF proof: constant-term and signed `S_n` symmetry approaches
  recast OCF on the tiling hypercube.
- H is band-limited: Walsh degree bounds imply tournament H only sees low-ish
  frequency structure.
- Spectral Legendre/ruler sequence: 2-adic valuation ladder in Walsh amplitudes;
  dual to 2-adic tournament-size column families.
- Signed adjacency polynomial identities: low-n polynomial identities, even
  cycle cancellation, and diagonal signed position theorems.
- Transfer chain metagraph: view metagraph moves as a Markov/transfer process.
- Hopf algebra / subset convolution: Redei-Berge and noncommutative analogues
  suggest a coproduct route to OCF.
- Noncommutative Redei-Berge deletion-contraction: Mitrovic lead with direct
  relevance to inductive proofs.
- Chromatic Redei bridge: Mitrovic-Stojadinovic connects Redei-Berge functions
  with chromatic symmetric functions and h/e-positivity questions.
- Stanley-Stembridge / Hikita: e-positivity for unit interval graphs may
  constrain adjacent tournament positivity frameworks.
- Worpitzky / Eulerian expansions: forward-edge polynomials, cumulants, mod-3
  universality, and signed variants.

## Recurrences, Constants, and Arithmetic Motifs

- k-nacci traces: `Tr(M_k^n)` creates Mersenne-like and tribonacci coincidences;
  causal links to forbidden H were later corrected.
- Tribonacci trace 7 and 21: numerical resonance with forbidden values, but
  not the proof mechanism.
- Golden shadow `f(n)=(n-2+sqrt(n^2+4))/2`: metallic-mean/continued-fraction
  family connected to transfer recurrences.
- Tau-phi clock: complex tribonacci eigenangle approximately tracks `ln(2)`;
  not exact but conceptually tied to information clocks.
- Cayley monad `D_4`: Cayley transform orbit of 2, fixed points +/-i, and
  Pauli/Bloch-sphere analogy.
- Cayley music / just intonation: superparticular ratios and arctanh/Cayley
  products connected to recurrence spectra.
- Prime telescope: prime distribution and forbidden-21 products seen through
  a telescoping recurrence lens.
- 42 -> 43 transition: 42 as Bernoulli/Hurwitz/tournament obstruction number,
  43 as the next prime transition.
- Bernoulli `B_6=1/42`: Von Staudt-Clausen puts primes 2,3,7 in one
  denominator; a recurring explanation for why 42 keeps appearing.
- Base-42 Erdős-Straus: Egyptian fraction decomposition of `4/p` organized by
  residues mod 42 and splitting in cyclotomic fields.
- Cyclotomic splitting for `k/N`: solvability of `k/N=1/b+1/c` connected to
  primes splitting in real cyclotomic subfields.
- Euler product analogy: odd cycles as "primes" of tournament structure,
  OCF as a finite Euler product.
- Riemann zeta tournament: zeta analogies were explored; several negative-zeta
  causal claims were corrected.
- Totient/divisor/harmonic arithmetic of 42: Euler-phi, divisor, harmonic, and
  Bernoulli patterns used as arithmetic metaphors.
- Wilson/double factorial fixed points: odd double factorials stabilize modulo
  squarefree even moduli; related to ladder ratios.
- Lucas / Zeckendorf / summand graph: non-consecutive representations and
  summand graphs were linked to recursive tournament geometry.
- Seven, twenty-one, forty-two, seventy: forbidden values multiplied by
  triangular or Bernoulli-related factors.
- Seventeen Vitali atoms: 17 appears as an atom count in certain decomposition
  pictures.
- Why h+1 is prime for exceptionals: Lie exceptional Coxeter data were mapped
  onto tournament numerology.

## Geometry, Polytopes, and Combinatorics

- Simplex inside cube: H geometry viewed as simplex/cuboid interaction in the
  tiling cube.
- Simplex-cuboid duality: doubling, gaps, forbidden H, and nesting
  obstructions framed by simplex-vs-cuboid volumes.
- Demicube corner theorem: user pattern that demicube corner pieces have
  volume `1/n!`; tied to simplex/cube geometry.
- Quotientope and clique complex: quotient polytope language for metagraph
  cells and Omega independence/clique structures.
- Polygon-simplex-staircase trinity: three geometric models of the same
  tournament data.
- Amplituhedron / 987 / chemistry: Fibonacci/Lucas and positive-geometry
  analogies to molecular structure.
- Barvinok-BIBD-QR triangle: counting, designs, and QR structure placed in one
  optimization triangle.
- TSSCPP / alternating sign matrices: self-evacuating SYT counts and staircase
  tableaux suggest ASM/TSSCPP links.
- Hook lengths / RSK / Young tableaux: staircase hook lengths are odd; sorting
  networks and tournament grids share the same staircase of comparisons.
- Viennot heaps: tournaments as comparison outcomes and sorting networks as
  comparison orders; heaps may give a Redei proof route.
- Ballot sequences / Dyck paths: clean binomial counts for Type-II positions
  suggest ballot/Dyck proofs.
- Lindstrom-Gessel-Viennot: explored as a possible bijective route to even/odd
  path splitting.
- Bags of sticks: decomposition metaphor from Redei-Berge literature; mostly
  a dead end for OCF but worth remembering.
- Bessis / dual braid monoid / noncrossing partitions: tiling model and path
  chambers have parallels with braid arrangements and cyclic sieving.
- Equidecomposability / Hilbert's Third Problem: `(H,beta_1)` as Dehn-like
  invariant at small n.

## Physics, Chemistry, and Systems Analogies

- Antiferromagnetic tournament: arcs as frustrated spin interactions; H
  correlates with frustration/ground-state structure.
- Fiber-bundle antiferromagnet: local arc choices as fibers over quotient
  geometry.
- Wick rotation / spin-1 Ising: `H(T)` as Wick-rotated partition function at
  complex temperature.
- Yang-Lee zeros: independence roots and fugacity axis interpreted through
  statistical physics.
- Fugacity axis / five points: special evaluations of `I(G,x)` organize
  vanishing, oblong numbers, and dimension transitions.
- Dimension axis: evaluation point x determines how much graph structure a
  theory sees; x=2 sees "all dimensions."
- Golden ratio / tribonacci boundary: phi as edge dimension, tau as triangle
  dimension, 2 as infinite-dimensional evaluation.
- Cayley-Dickson tower: R,C,H,O and beyond used as an analogy for transformer
  architecture and physics/chemistry unification.
- Cartan decomposition: tournament active/dark split compared with
  `gl(n)=so(n)+p+R` and gauge-theoretic decompositions.
- Ghosts and supersymmetry: cuboid ghost/partner cancellations interpreted
  with supersymmetry language.
- Super-orthogonality: repeated hidden orthogonality patterns in projections,
  residues, and spectral decompositions.
- Crystallization: when lossy invariants become rigid enough to reconstruct
  or classify structures.
- Protein folding tower: tournament/path-homology ideas as possible analogues
  for folding landscapes.
- Chemistry and independence polynomial: hard-core gas and molecular graph
  intuition for Omega independence.

## Computer Science, LLMs, and Engineering

- Tournament attention: threshold transformer attention matrices into
  tournaments; measure H, Omega, homology, and Cartan active/dark ratio.
- Soft OCF: differentiable sigmoid tournaments and soft independence
  polynomials as model regularizers or uncertainty signals.
- OCF for LLM uncertainty: layer-based OCF, MC dropout, and attention-head H
  proposed as correctness/uncertainty predictors.
- BLITZRANK connection: tournament graphs already used for zero-shot LLM
  ranking; this repo adds parity/H/sensitivity tools.
- Tournament-aware consensus: Paxos/Raft and distributed agreement interpreted
  with tournament inconsistency features.
- Active learning intelligence: spectral surprise and tournament uncertainty
  as a definition of where a model should query next.
- Polynomial agent: architecture whose model is literally a polynomial
  invariant/calculator.
- Code as compressed mathematics: repo scripts and Lean proofs as compressed
  theorem statements.
- mod_rank library: modular rank routines for sparse topology computations,
  intended as a practical package.
- circulant_homology module: reusable computation for circulant and Paley path
  homology.
- Tournament TDA feature extractor: practical feature block combining scores,
  H, Omega alpha, support residue, deletion loss, even graph, and buckets.
- Cryptographic vulnerabilities: Walsh puncturing, QR/cyclotomic splitting,
  QC-LDPC eigenspace attacks, and small-prime torsion checks.
- Genome contig ordering: post-assembly ordering as a near-transitive
  tournament Hamiltonian path problem.
- Ranking / sports / elections / manipulation: H and projection defects as
  inconsistency or match-fixing detectors.
- Photo compression wall / fractal codec: tournament compression and
  irreducible residuals as analogies for image/code compression.

## Literature Leads and Proper Nouns

- Grinberg-Stanley: OCF proof background and mod-4 theorem.
- Forcade: F2-invariance and parity of paths/circuits.
- Noga Alon: Paley maximizer conjecture and upper bounds for Hamiltonian paths.
- Szele: asymptotic lower/average Hamiltonian path context.
- Adler-Alon-Ross: random regular/asymptotic e behavior.
- Irving-Omar-Schur: determinant/permanent and positivity barriers.
- Chudnovsky-Seymour: claw-free graph real-rootedness theorem used for
  Omega real-rootedness through n=8.
- Heilmann-Lieb / line graph route: tempting real-rootedness route, later
  refuted by Beineke-style obstructions.
- Bezakova dichotomy / subdivided claws: computational complexity and
  claw-structure leads.
- Jerrum-Patel zero-free regions: polymer/independence-polynomial literature.
- Schweser-Stiebitz-Toft: Redei revisited, mixed graph extensions, and
  expository hierarchy.
- Mitrovic / Mitrovic-Stojadinovic: noncommutative deletion-contraction and
  chromatic-Redei bridges.
- Hikita / Stanley-Stembridge: proved e-positivity for unit interval graphs,
  potentially constraining related positivity questions.
- Tang-Yau: circulant digraph path homology via Fourier/symbol matrices.
- Ren: path independence complexes and embedding theorem leads.
- Hepworth-Roff / Asao / Chaplin: magnitude-path and random digraph homology.
- Awan-Bernardi B-polynomial, Sazdanovic-Yip categorification: digraph
  polynomial/categorification connections.
- Savchenko, Komarov-Mackey, Linial-Morgenstern: cycle-counting formulas and
  extremal regular tournament literature.
- Pantangi / Terwilliger algebras: distinguish Paley and other DRTs.
- Striker / Chapman: posets, alternating sign matrices, and tournament
  tableaux connections.
- Eplett / self-converse tournaments: classical self-complementary structure.

## Meta-Patterns and Process Lessons

- Support-residue calculus: choose support, project/forget, measure residue;
  now a unifying lens for good cuts, H=63, root failures, bucket transport,
  Paley/interval, and homology ghosts.
- Residue rank threshold: first obstructions may appear when surviving residue
  has just enough parity, rank, or independence to avoid cancellation.
- Projection commutators: where two natural projections nearly commute, the
  commutator often marks the next theorem.
- A squared principle: second-order structure (`A^2`, pair counts, two-step
  walks) repeatedly captures information invisible at first order.
- Everything is the triangle: 3-cycles are the primitive obstruction, but
  longer cycles encode layered corrections.
- Grand trichotomy / three types: fixed, boundary, free; lossless, lossy,
  crystallized; ramified, inert, split.
- Z/2 everywhere: complement, parity, two-sheeted covers, sign representations,
  and xor masks keep reappearing.
- Hidden orthogonality: cancellation is often exact because natural subspaces
  are orthogonal after the right transform.
- Computational irreducibility: some invariants resist compression; failures
  are data about what the quotient forgot.
- Territory and map: inert vs split behavior is a recurring metaphor for when
  observations preserve or fracture structure.
- Repo as tournament: multi-agent research itself forms a tournament of
  claims, mistakes, court cases, and convergence.
- Mistakes as boundary detectors: several major discoveries came from
  correcting overstrong analogies, especially blue/black terminology,
  zeta-forbidden claims, real-rootedness, and Mersenne-forbidden extrapolation.

## Dead-End Warnings Worth Preserving

- Naive arc-reversal cycle bijection fails; arc flips must be handled at the
  sum/algebra level.
- Tournament contraction is not generally a tournament because mixed cases
  leave arcs undefined.
- Naive Type-II/orphan involution fails because the swapped path can become
  invalid.
- Contiguous-block decomposition of cycles in Hamiltonian paths gives wrong
  counts; the correct decomposition is algebraic through Omega.
- Ratio extrapolation for Paley cycle counts was misleading.
- Line-graph/Heilmann-Lieb real-rootedness route does not characterize Omega.
- Universal real-rootedness/TRRT is false at n=9, though the subclass remains
  interesting.
- Mersenne numbers are not generally forbidden H values; 31, 63, and 127 are
  achievable in later data.
- "Blue/black" legacy terminology can be dangerously ambiguous; use strict
  line-color definitions or spine/ribs/sea for class edges.
- HYP number collisions exist; check theorem/hypothesis identity before
  trusting a bare number.

## Quick Source Map

- `00-navigation/TANGENTS.md`: dense chronological tangent entries, many with
  tags and certainty labels.
- `00-navigation/INVESTIGATION-BACKLOG.md`: actionable leads, literature
  pointers, and next steps.
- `05-knowledge/hypotheses/INDEX.md`: large hypothesis ledger including
  confirmed, refuted, corrected, and open speculative claims.
- `07-reflections/`: philosophical/meta-structural essays; file titles are
  often the best quick clue to a thread.
- `03-artifacts/explainers/`: plain-language summaries and prior-work
  positioning.
- `01-canon/MISTAKES.md`: essential guardrail for avoiding seductive but false
  repeats.

