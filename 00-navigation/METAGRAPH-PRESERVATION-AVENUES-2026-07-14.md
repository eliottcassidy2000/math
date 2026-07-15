# Metagraph Preservation Avenues — pull cards for future agents

**Created:** codex-2026-07-14-S4.  **Purpose:** a living backlog of old niche
threads that can refine the canonical tournament-tiling address or pull that
address back to the LRC(14) four-coordinate slope suspension.

The central rule is:

> Never ask only what a quotient remembers.  State the predicate it must
> preserve, the continuation operations it must survive, and the smallest
> sidecar that repairs its mixed fibres.

## Contribution protocol

Add the next `MPA-*` card with four fields:

- **Pull:** the proposed computation or lemma.
- **Old thread:** exact repository pointers being revived.
- **Must preserve:** theorem-facing predicate and continuation operations.
- **Deliverable:** one finite output, counterexample, lemma, visualization, or
  formally named residual.

Do not default vertices to runners or arcs.  Candidate vertices include
tilings, line orbits, merged nodes, gaps, fixed sections, boundaries, wall
events, residues, sheets, cover arcs, Fourier modes, matroid circuits,
certificate clauses, and proof obligations.

## Immediate metagraph pulls

### MPA-01 — Rooted blue/black WL frontier beyond `n=7`

- **Delivered at n=7:** exact enumeration of `32768` tilings, `456` classes,
  and `272` merged nodes found one connected line graph and 272 rooted line-WL
  cells after three rounds.  Raw incidence has only 159 cells; deletion-parent
  ancestry has 270 and leaves two twin pairs, both resolved by line recursion.
- **Pull:** push to `n=8` using existing nauty class machinery, or first prove
  a structural reason the refinement is complete on the `n<=7` line graphs.
- **Old thread:** HYP-6825; HYP-3808/3809; `merged_metagraph_*` scripts.
- **Must preserve:** line colour, multiplicity, loops, root, and converse
  orbit; do not collapse to unweighted adjacency.
- **Deliverable:** proof of completeness at `n=8`, or the first collision
  cells with canonical codes and exact tiling fibres.

### MPA-02 — Lift a node address to the line metagraph

- **Pull:** refine each blue/black edge by the simultaneous-isomorphism orbit
  of its endpoint pair and test whether line-orbit profiles separate any WL
  twins.
- **Old thread:** `04-computation/line_metagraph_klein_S163.py` and
  `the-line-metagraph-is-rigid-klein-S163.md`.
- **Must preserve:** a particular antipodal line, not only its unordered
  endpoint class pair.
- **Deliverable:** maps `tiling -> line orbit -> endpoint pair -> node pair`
  and all inverse fibres at `n<=7`.

### MPA-03 — Compare 1-WL, 2-WL, spectra, and canonical codes

- **Pull:** create the first honest hierarchy of metagraph fingerprints and
  report exactly which pairs each layer separates.
- **Old thread:** HYP-3811 spectral twins; `taxonomy-of-merged-metagraph.md`;
  HYP-3813 Gram/spectral bridge.
- **Must preserve:** edge channels independently; a cospectral or 1-WL twin
  is an object, not a failure to hide with a code.
- **Deliverable:** collision lattice and minimal separator for every twin.

### MPA-04 — Make deletion ancestry a genuine groupoid functor

- **Delivered endpoint case:** THM-796 identifies fixed-path tilings with
  marked Hamiltonian-path objects.  Terminal endpoint deletion is an honest
  functor; after forgetting the path it becomes a weighted span, and its node
  kernel is not Markov-compositional.
- **Pull:** replace “delete vertex zero from one representative” by all
  vertex-deletion orbits, their multiplicities, and transport under
  automorphisms/converse.
- **Old thread:** HYP-3106; perspective-depth/observer-cut controlled-
  forgetting reflections; explorer parent labels.
- **Must preserve:** which deletion perspective is being followed and whether
  different perspectives merge later.
- **Deliverable:** recursive address tower `G_n -> multiset(G_{n-1})` with
  invariance tests and first nontrivial ancestry twins.

### MPA-05 — Half-tiling recursion without the false dyadic tower

- **Delivered algebra:** THM-796 identifies lines with the quotient code
  `F_2^m/<1>` and blue lines with the zero fibre of
  `delta([t])=t+sigma(t)`.  Mode-B gives exact line/blue/defect kernels of
  dimensions `2n-5,n-2,n-3`; inherited and fresh blackness are distinct.
- **Pull:** express the SC blue subgraph through its anti-diagonal fixed
  subcube, then identify what recursion remains after the known “repeat the
  complement fold” correction.
- **Old thread:** HYP-3244; HYP-3811 quarter-tiling correction;
  `the-blue-subgraph-IS-the-half-tiling-metagraph-...md`.
- **Must preserve:** the actual involution and its fixed set; do not invent a
  second complement action after restriction.
- **Deliverable:** a corrected recursive address or a proof that recursion
  needs a different mirror/gauge.

### MPA-06 — Fibre size from reversible metagraph dynamics

- **Delivered complement-line kernel:** THM-796's coloured half-edge matrix is
  symmetric, has row sum equal to tiling-fibre size, and normalizes to a
  reversible kernel.  It is still recursively insufficient: repeated endpoint
  deletion fails the node Markov product in nearly every row by `n=7`.
- **Pull:** revisit the one-tile-flip Markov chain, its stationary distribution
  proportional to tiling fibres, and ask which weighted transition data
  reconstructs fibre size without circularly assuming it.
- **Old thread:** `tiling_count_from_metagraph_s279.py`,
  `tiling_to_metagraph_s20cu.py`.
- **Must preserve:** directed transition weights and self-loops, not only the
  simple underlying graph.
- **Deliverable:** a non-circular reconstruction theorem or an explicit pair
  of weighted shadows showing which datum is missing.

### MPA-07 — Strongly connected prime words as an address coordinate

- **Pull:** factor each tournament class into its ordered strongly connected
  components and test whether the free-monoid word predicts line-WL colour,
  depth, or fibre size.
- **Old thread:** `iso-classes-are-a-free-monoid-on-strongly-connected-primes-...md`.
- **Must preserve:** component order and converse reversal of the word.
- **Deliverable:** factor-word annotation of every `n<=7` node and collision
  savings before canonical fallback.

### MPA-08 — Depth versus width versus line address

- **Pull:** jointly plot MFAS depth, flip-rank/covering width, automorphism
  group, and rooted line-WL address.
- **Old thread:** `two-axes-of-the-tournament-metagraph-...md`;
  flip-rank/covering-code reflections; Paley/heptagon work.
- **Must preserve:** depth and covering width as distinct coordinates.
- **Deliverable:** determine whether the hardest covering classes occupy
  exceptional line-WL branches rather than merely maximum depth.

## Quotient and stalk pulls

### MPA-09 — The ultrafilter descent obstruction as address purity

- **Pull:** search for tile predicates or tile-family filters that are pure on
  rooted line-WL cells even though individual coordinate ultrafilters fail to
  descend to isomorphism classes.
- **Old thread:** HYP-2245 metagraph ultrafilter sheaf.
- **Must preserve:** which side of every antipodal tiling line is selected.
- **Deliverable:** first nontrivial quotient-pure filter, or a no-go bound on
  purity at `n<=6`.

### MPA-10 — Recover the metric-difference middle category

- **Pull:** decorate a node/tiling address with circular gap multiset and ask
  the minimum metric label that makes the LRC witness predicate fibre-pure.
- **Old thread:** `three-observer-categories-tiling-is-relative-tournament-is-blind.md`.
- **Must preserve:** translation-sensitive observer placement, scale-sensitive
  gap metric, and the distinction between both.
- **Deliverable:** smallest sampled pure quotient between raw tournament and
  full phase vector.

### MPA-11 — Mixed Haar square as the missing fibre cocycle

- **Pull:** attach the `2x2` mixed Haar sign to local metagraph squares and
  test whether it separates address fibres with identical row/column shadows.
- **Old thread:** HYP-2989 and the Haar-product/tournament-tiling scripts.
- **Must preserve:** the mixed sign
  `T00-T01-T10+T11`, not only margins.
- **Deliverable:** a cocycle-labelled square complex and its first quotient
  collision/repair theorem.

### MPA-12 — Continuation equivalence for controlled forgetting

- **Delivered prime-seven separator:** THM-773 gives two exact `c=7` states
  with the same node, least-path mask `27833`, six-moment vector, and
  owner-to-sheet assignment, but different next futures `(owner 10, sheet 6)`
  and `(owner 4, sheet 5)`.  Inverse steps and endpoint schedule are therefore
  mathematically necessary fields.
- **Pull:** formalize two states as equivalent only when every legal future
  wall, deletion, lift, and threshold continuation has the same terminal
  result.
- **Old thread:** HYP-3106; observer-cut orbit ledgers; HYP-3513 Nerode.
- **Must preserve:** future behaviour, not merely current witness truth.
- **Deliverable:** a minimized finite automaton for one bounded clock family,
  with counterexamples to any coarser state signature.

### MPA-13 — Automorphism stabilizers of the observer cut

- **Pull:** compute how `Aut(T)` acts on the fixed-path tiling fibre and on
  observer/gap markings; measure the label tax needed to select a cut.
- **Old thread:** observer-relative tiling hierarchy; iso-class/fibre formula;
  controlled-forgetting ledgers.
- **Must preserve:** stabilizer orbit, not just `|Aut|`.
- **Deliverable:** orbit decomposition and an exact cut-address compression.

## Pulling the atlas back to the LRC slope object

### MPA-14 — Constructible metagraph atlas on the affine suspension

- **Pull:** for named affine families `cP+R`, partition the `(u,t)` torus into
  phase-order chambers, attach node/tiling addresses, and record wall transport.
- **Old thread:** HYP-6815 affine-slope threshold suspension.
- **Must preserve:** nonemptiness of the integer-slope section at each
  threshold, exact clearance, owners, and transverse shear `R`.
- **Deliverable:** one exact chamber atlas with a verified round trip from
  slope section to LRC maximin.

### MPA-15 — Sector-carry movie rather than static tournament

- **Delivered prime-seven slice:** THM-773 gives the exact token update
  `k_a -> k_a-w_a^{-1}`, first-moment return holonomy, complete six-moment
  exact-cover test, and global carry `k_a(x+1)=k_a(x)-1`.  The general
  14-sector/metric movie remains open.
- **Pull:** implement the state `(14-sector potential, microphase weak order)`
  and annotate each wall crossing by its metagraph address change.
- **Old thread:** `lrc14-prime-grid-cover-sector-carry-and-threshold-tournament-...md`;
  sector-vector realizability.
- **Must preserve:** threshold walls at `1/14,13/14`, closed endpoint flag,
  and global microphase order.
- **Deliverable:** monodromy census: which address transitions occur, which
  carry laws forbid transitions, and whether good/bad futures become pure.

### MPA-16 — Prime-grid private-owner filtration

- **Pull:** treat a prime-grid translate cover as a filtration by lift depth;
  branch first on the translate with largest private fibre and record address
  of the blocker handoff graph.
- **Old thread:** persistent prime-grid cover reflection; published-sieve
  cover reformulation.
- **Must preserve:** private owner, common multiplication orbit, lift carry,
  and strict-versus-closed danger.
- **Deliverable:** exact lift tree for the first hard cover or an inverse
  theorem forcing AP/Goddyn-Wong structure.

### MPA-17 — THM-761 seven-exception cyclic tiling wall

- **Delivered at the prime lens:** THM-771 supplies the exact defect and wall
  tear; THM-773 identifies every `c=7` exact partition with a permutation of
  `F_7`, gives its polynomial/moment certificate, and supplies owner-event
  holonomy.  Ramified `c>7` normal forms remain open.
- **Pull:** classify covers of the sheet cycle `Z_c` by seven equal-length
  bad arcs, with gcd multiplicities and movable core-safe time `t_0`.
- **Old thread:** THM-761 multi-exception sheet theorem.
- **Must preserve:** shared sheet index, each exception's owner and gcd
  stratum, and freedom to move `t_0` through the core safe set.
- **Deliverable:** either a free-sheet strengthening or a complete tight-wall
  normal form plus handoff to a non-sheet witness.

### MPA-18 — Compare cyclic sheet tilings with staircase tilings

- **Delivered:** THM-773 constructs the map.  The cut gauge sends all 5040
  owner assignments to `n7-a000`; circular precedence sends all to
  `n7-a267`; a least-Hamiltonian-path gauge maps them onto exactly the 25
  masks in `a267`'s inverse fibre.  The sharp obstruction is also proved: the
  iso node is constant and even `(node,mask,assignment)` is not
  continuation-complete without inverse-step and endpoint transport.
- **Pull:** search for a functor, not an analogy, between seven-arc cyclic
  covers and fixed-Hamiltonian-path binary staircase fillings.
- **Old thread:** THM-761; HYP-6825; half-tiling fold.
- **Must preserve:** cover nonemptiness on the cyclic side and class/fibre
  membership on the staircase side.
- **Deliverable:** an explicit predicate-preserving map, or a sharp proof of
  which information blocks one.

### MPA-19 — Pair-sum obligations as metagraph vertices

- **Pull:** make HYP-6785 obligation/blocker pairs the tournament vertices;
  orient by forced ownership or discharge order and attach their state to the
  affine metagraph chamber.
- **Old thread:** HYP-6785; THM-668 pair-sum ruler.
- **Must preserve:** protected pair endpoints, nonunit multiplier, blocker
  owner, and endogenous denominator.
- **Deliverable:** a peel/descent rule that strictly decreases the obligation
  complex on the remaining scale-normal residual.

### MPA-20 — Endpoint sidecars as a minimal separating family

- **Pull:** combine divisor mask, cap ratio, signed endpoint blocks, peak
  witness, and metagraph address; compute the smallest subsets that determine
  each proof predicate on a named corpus.
- **Old thread:** `lrc14_endpoint_tournament_sidecar_audit_structure_S1.py/.out`.
- **Must preserve:** distinguish truth of covering, THM-755 cap, exact `M`,
  and endpoint measure rather than optimizing one scalar.
- **Deliverable:** a predicate-by-sidecar observability matrix and all minimal
  separating sets.

### MPA-21 — Restricted quotient stack revisited with exact addresses

- **Pull:** rerun source-deleted, observer-marked, gap-threshold, kinetic-gap,
  blocker-deficit, and apex-boundary quotients using `n-aXX` as the base label.
- **Old thread:** `lrc-restricted-tournament-quotient-stack-s535.md`.
- **Must preserve:** source/witness purity and transition purity separately.
- **Deliverable:** label-tax, purity-curvature, kinetic-torsion, and monodromy-
  defect tables over exact metagraph fibres.

### MPA-22 — Fixed sections and boundary events as vertices

- **Pull:** replace runner vertices by the 14 fixed circle sections or their
  boundaries; orient events by forced precedence/flux and fibre the result
  over metagraph addresses.
- **Old thread:** HYP-2024 and
  `lrc-section-boundary-functors-and-assumption-challenge-s539.md`.
- **Must preserve:** which observer-adjacent sections are clear, compactified
  wall events, and ownership of entry/exit.
- **Deliverable:** first section/event quotient whose fibres remain pure over
  a full arithmetic period, not merely open chambers.

### MPA-23 — Certificate sheaf over the address graph

- **Pull:** place proof obligations on nodes and edge transports, then test
  whether every source-avoiding section exports debt toward a terminal node.
- **Old thread:** HYP-2101 and `lrc-apex-lift-certificate-sheaf-s579.md`.
- **Must preserve:** local certificate, overlap compatibility, and residual
  debt under lift.
- **Deliverable:** a finite Cech/descent check on one hard family, with any
  non-gluing cocycle exposed explicitly.

### MPA-24 — Danger-Gram spectrum as a metric stalk

- **Pull:** attach the danger-overlap Gram spectrum to each metagraph/tiling
  state and test whether spectral changes localize at particular address
  branches or wall types.
- **Old thread:** HYP-3813 integrating merged metagraph with LRC; HYP-3571.
- **Must preserve:** the full relation or Gram spectrum, not trace/second
  moment alone.
- **Deliverable:** true LRC spectral twins (same elementary profile, different
  witness future) or a separator theorem for named tight/near-tight families.

### MPA-25 — Threshold persistence instead of a single `1/14` slice

- **Pull:** compute the filtration of safe-time components as `lambda` varies,
  with node/tiling/owner labels on births, mergers, and deaths.
- **Old thread:** affine suspension HYP-6815; endpoint/peak sidecar audit;
  persistence and cover lenses throughout the LRC frontier.
- **Must preserve:** closed endpoint convention and peak owner.
- **Deliverable:** labelled barcode/merge tree for AP, Goddyn-Wong, deep well,
  and the first scale-normal residual; identify which barcode features force a
  `1/14` section.

### MPA-26 — Learn the finite Nerode state, then prove its fields necessary

- **Pull:** start from the full exact packet, minimize it by continuation
  equivalence on bounded clock movies, and convert each surviving field into a
  lemma or an adversarial pair.
- **Old thread:** HYP-3513 private-firewall Nerode; certificate spigot/router
  work; HYP-3106.
- **Must preserve:** every future terminal witness/certificate class.
- **Deliverable:** a machine-minimal state schema plus a mathematical reason
  for each field.  This is the most direct answer to “precisely what
  information needs to be preserved?”

### MPA-27 — Fold the dynamics before quotienting it

- **Pull:** imitate THM-774's exact half-sum/half-difference conjugacy on the
  prime-seven stalk.  Test discrete Fourier modes of labelled occupancy,
  signed token chords, and elementary symmetric coefficients of the
  redundancy polynomial as coordinates for endpoint deletion and translation.
- **Old thread:** THM-774 folded diamond; THM-773 token polynomial and
  monodromy; Fourier/relation-line fields in HYP-6815.
- **Must preserve:** exact cover divisibility, the owner-specific update
  `k_a -> k_a-w_a^(-1)`, endpoint deletion, and the global carry.  A scalar
  norm, energy, or measure is not an acceptable substitute for the pointwise
  transformed predicate.
- **Deliverable:** either an invertible coordinate transform in which each
  wall event changes boundedly many fields, or an explicit pair proving that
  every proposed low-order mode truncation merges different continuations.

### MPA-28 — Continued-fraction substitutions on the labelled path stalk

- **Delivered base decoder:** THM-778 proves that centered Beatty ranks recover
  every simultaneous endpoint wall and that pairwise words carry an Euclidean
  parity cocycle.  The ordinary next-event tournament stays transitive; the
  eight-owner movie nevertheless visits 948 labelled Hamiltonian paths.
- **Delivered mask guardrail:** the ten covered walls are mask-labelled, but
  reflection does not descend to the 25-mask quotient: only 9/25 mask sources
  have unique images over all 5,040 assignments.
- **Delivered redundancy signal:** duplicate roots follow the five-step word
  `(1->0,4->6,2->4,0->2,6->5)` twice, a simpler transport coordinate than the
  nonpalindromic mask word.
- **Delivered finite base:** the nine between-wall block lengths are
  `(57,301,3,24,329,24,3,301,57)`, reducing the exact movie to five block types
  and their mirror.
- **Pull:** attach each Euclidean/centered-Christoffel substitution to its
  induced transition on owner tokens, redundancy roots, and the 25 masks over
  `n7-a267`.
- **Old thread:** THM-536 Sturmian cover word; THM-637 Farey roof; HYP-4078
  Ostrowski ladder; THM-745 Euclidean residual; THM-773 prime-sheet fibre.
- **Must preserve:** common gcd, midpoint phase bit, simultaneous blocks,
  owner-specific inverse steps, metric wall rank/time, core-safe component,
  mask, redundancy factor, and global carry.
- **Deliverable:** a centered Farey concatenation law plus a finite
  substitution automaton for the ten-wall r=8 word, minimized by continuation
  equivalence.

### MPA-29 — Hamiltonian-path orbit fibres instead of cube scans

- **Delivered theorem:** THM-781 gives the literal maps
  `tiling -> canonical class -> converse-merged node` and
  `node -> union HP(T)/Aut(T)` over its one or two classes.
- **Delivered size law:** each unmerged class fibre has
  `H(T)/|Aut(T)|` tilings; a non-self-converse merged pair has twice that
  number.  Every forward and inverse record through `n=7` agrees exactly.
- **Heptagon explanation:** `n7-a267` has `H=175`, `|Aut|=7`, hence its 25
  masks are intrinsic observer-cut orbits.  THM-773's 5,040 sheet assignments
  map onto these path orbits with nonuniform multiplicity.
- **Pull:** lift each `HP(T)/Aut(T)` element to its owner-labelled assignment
  fibre, then attach THM-778 endpoint transport to transitions between cut
  orbits.  Separate a mask transition from a complement-line orbit and from a
  mere endpoint-node incidence.
- **Concurrent join:** THM-783's wall `phi`-recurrence and visitor clusters are
  transport above this cut orbit, not node invariants.  Test whether `(path
  orbit, owner assignment, phi, minimal persistent stalk)` is sufficient for
  its flip-break step before adding more metagraph coordinates.  Raw wall count
  is excluded by MISTAKE-147.
- **Must preserve:** path orbit, owner assignment, simultaneous-isomorphism
  line orbit, inverse winding, centered phase/ties, wall rank, and carry.
- **Deliverable:** the first exact three-level incidence index
  `node <- path orbit <- owner-labelled LRC state`, with forward maps,
  inverse fibres, and continuation-purity audit at every quotient.

### MPA-30 — Disintegrate symmetric line flux over asymmetric node fibres

- **Delivered general laws:** THM-785 makes `C3` exact inverse score variance,
  proves complement flux `d0-dlast-1`, and decomposes it as
  `2(d0-n/2)-epsilon`.  Blue reflection sets `epsilon=0`; closed binomial
  formulas give every blue step, every total step, and every black
  `|epsilon|` count for arbitrary `n`.
- **Delivered joint law:** `(Delta C3,epsilon)` has Laurent kernel
  `v^(-p)(u+u^-1)(1+uv)^p(1+u^-1v)^p`, `p=n-3`.  Its coordinates have zero
  covariance and variances `(n-1)/2,(n-3)/2`; quotient drift is a
  disintegration of a centrally symmetric two-dimensional measure.
- **Delivered quotient signal:** line colour permits only the category path
  `pure_blue --blue-- mixed --black-- pure_black`.  At `n=7` the blue boundary
  is entirely outward, while increasing-`C3` black flow has 2,798
  pure-black-to-mixed instances against 1,254 mixed-to-pure-black.  This does
  not break complement symmetry: `(a,epsilon)` is sent to `(-a,-epsilon)`
  before quotienting.
- **Delivered first disintegration:** the `n=7` reverse raw current is a fibre-
  volume effect.  Pure-black has 25,460 black masks and mixed has 6,796; after
  source-mass normalization, reverse rate is `10.99%` while outward rate is
  `18.45%`.  Raw and density gauges point in opposite directions.
- **Pull:** condition the exact endpoint-bit measure `(R,X,Y)` successively on
  unmerged class, converse orbit, fibre size, automorphism stabilizer, line
  orbit, and pure/mixed/black phase.  Explain the residual density-normalized
  outward kernel and the one-sided `|epsilon|=3,4` tails, then predict both raw
  and normalized signs without enumerating masks.
- **Old thread:** score-majorization/Slater depth; half-tiling reflection;
  line-metagraph simultaneous isomorphism; HYP-3813 Gram/spectral bridge;
  MPA-02/05/07/13; THM-781/785/787.
- **Must preserve:** signed `epsilon`, line-instance multiplicity, both endpoint
  masks, simultaneous-isomorphism line orbit, converse action, and the `C3`
  orientation.  Absolute defect or simple projected adjacency alone cannot
  explain directional drift.  An unordered complement line must be oriented
  by an invariant such as `C3`; choosing the numerically smaller mask makes a
  “toward regular/transitive” ratio convention-dependent.
- **Deliverable:** an orbit-stabilizer formula for the conditional flux kernel,
  an `n=8` exact/semiexact census, and either a theorem predicting reverse
  black drift or the smallest counterexample.  Test separately whether every
  balanced node remains nondecreasing-`C3` reachable by a `B* K*` word.

### MPA-31 — Three-sorted recursive incidence and the missing line phase

- **Delivered theorem:** THM-796 proves
  `T_n=(T_(n-1) x_[T_(n-2)] T_(n-1)) x {0,1}_apex`; after choosing the
  apex-zero endpoint, `E_n` is the compatible-face pullback.  Passing to bare
  lower lines creates a uniform `C2` endpoint-phase torsor over
  `E_(n-1) x_[E_(n-2)] E_(n-1)` for every `n>=5`.
- **Delivered node law:** the coloured directed endpoint kernels are symmetric,
  have row sum equal to node fibre size, and have even diagonal; off-diagonal
  line multiplicity is `A_uv`, while loop multiplicity is `A_uu/2`.  The
  weighted lower-face row and its primitive normalization separate all nodes
  through `n=7`, but strong lumpability fails on `1206/1312` nonzero blocks.
- **Delivered joint coupling:** `Xi_n` records the upper endpoint-node pair,
  both face-lines' endpoint-node pairs, and the three colours.  All other
  matrices are its marginals.  It has `16031` cells for `16384` `n=7` lines;
  the actual lower-line pair plus one coherent phase bit is exactly invertible.
- **Delivered colour law:** upper/low/high blue are pairwise independent but
  not jointly independent.  Upper blue forces equal face colours; both-face
  blue is parity-constant on gap diagonals and all-three blue is Toeplitz.
- **Pull:** decide whether primitive face rows remain injective at `n=8`, and
  refine the first collision by rooted path-orbit lift distributions rather
  than by an opaque canonical code.  Quotient the line `C2` torsor only after
  attaching THM-793's leg current and THM-785's signed `C3` orientation.  Use
  HYP-6865's Smith harmonic potential as an independent horizontal coordinate:
  determine whether its rare `n=7` discordances localize exactly in primitive-
  face collision or non-lumpable blocks.
- **Old thread:** THM-345/477/643/781/793; HYP-6865 Smith diagram;
  HYP-3051..3054 observer-extension
  calculus; MPA-01/02/04/05/06/29/30; continued-fraction transport THM-778.
- **Must preserve:** both tiling endpoints of a line, coherent face phase,
  colour, line multiplicity, loop factor two, Hamiltonian-path orbit, and the
  individual—not merely aggregate—lower-tiling lift distribution.
- **Deliverable:** an `n=8` collision/lumpability atlas and a minimal rooted
  sidecar that turns the weighted node correspondence into a genuine recursive
  transition object.

### MPA-32 — Minimize the recursive three-sorted address

- **Delivered exact carrier and bounded minimization:** THM-796 separates
  merged nodes `M_n`, tiling
  half-edges `X_n`, and complement lines `L_n`; it adds reflection defect,
  same/converse-sheet loop holonomy, and the simultaneous high/low half-edge
  tensor.  Independently sorting endpoint pairs loses coupling at `n=4`, and
  bare node transfer fails Markov composition at `n=5`.  The S11 audit then
  minimizes literal lines over `3<=n<=7` under labelled or unordered-successor
  endpoint deletion plus reflection.  At `n=7`, node boundary recursively refines
  `6076->16359` cells, while `Xi` refines `16031->16382`; the latter leaves
  exactly the phase deck pairs `0x12ca/0x12cb` and `0x146c/0x146d`.
- **Finding:** in the labelled-face language, adding colour, class sheet, and
  exact defect to node boundary gives the same stable partition; adding sheet
  and defect to `Xi` also gives the same partition.  The exact phase address
  is singleton and agrees with the literal-line control.  Both residual pairs
  have the same defect because their deck move flips the reflection-fixed
  apex, so phase and defect are transverse.  Together the four residual lines
  form a commuting phase/reflection square of parallel black edges between a
  non-self-converse two-sheet node and a self-converse node.  On the rigid
  endpoint the unique class isomorphisms are reflection-conjugate interior
  5-cycles; the inverse endpoint-fixed choices on the self-converse endpoint
  give order-five relative holonomy.  Under unordered successor cells,
  node boundary and `Xi` retain `8310` and `16380` cells respectively.
- **Three-face join:** the pair is a chord-dual suspension of reflected black
  `n=5` core loops.  It is born as cross-lines at `n=7` and cannot directly
  face-lift to `n=8`.  THM-801's gap face detects it immediately: phase one
  gives a black gap loop `33--33`, phase zero a black gap cross-line `21--33`,
  so `Omega_7` resolves the square without `B2`.
- **Next pull:** join MPA-34/35 while extending the operation alphabet to
  inverse lifts, internal deletion/repair, and an explicit terminal proof-
  obligation predicate.  Run the `n=8` continuation minimizer and compare it
  with rooted weighted line-WL and THM-809's already-injective static `Lambda`
  filtration.  Classify any transport twins among the 52 doubletons surviving
  through `tau=6`; the static collision frontier has moved to `n=9`.
- **Old thread:** HYP-3513 Nerode; HYP-3106 controlled forgetting; MPA-02/04/12;
  THM-781/793/796/801; line-metagraph simultaneous-isomorphism work.
- **Must preserve:** the named line, simultaneous half-edge coupling, exact
  defect ancestry, class-sheet holonomy, path witness, chord-duality
  isomorphism torsor modulo endpoint automorphisms, its face restrictions and
  overlap holonomy, and every declared future observation.  Record which of
  these fields minimization proves redundant rather than dropping it in
  advance.  On an LRC pullback, keep
  class-sheet holonomy distinct from reduced token return and metric component
  translation; finite-state closure does not imply metric closure.
- **Delivered artifact:**
  `three_sorted_metagraph_continuation_minimization_codex_S11.py/.out` is the
  first minimal bounded address for the declared forward operation languages.
  Line-WL collision comparison, inverse operations, and a stabilization
  theorem or obstruction remain open.

### MPA-33 — Push the exact `B3` line descent beyond its first phase witness

- **Delivered theorem:** THM-801 lifts `A+B+C-D-E-F+G` to high-end,
  gap-contraction, and low-end restriction maps.  Their compatible tilings
  glue uniquely.  For `n>=6`, the nonempty triple overlap forces zero `F_2`
  phase holonomy, so pairwise-compatible lower complement lines glue to one
  upper line.  At `n=5`, the empty triple overlap leaves 32 balanced and 32
  unbalanced compatible triples.  MPA-32's only stable labelled `Xi_7`
  collision square is an exact local witness: its endpoint faces agree while
  its two gap-line pairs differ, and `Omega_7` separates all four lines.
- **Pull:** formulate the general cohomological line-gluing theorem for an
  arbitrary finite cover and identify exactly which deepest intersections
  kill which `H^1(nerve;F_2)` phases.  Test four-face covers of the staircase
  and covers adapted to Mode-B.
- **Must preserve:** literal overlap line, canonical endpoint phase, the nerve,
  and whether complement acts freely on each deepest intersection.
- **Deliverable:** an all-cover phase-obstruction theorem plus the smallest
  staircase cover with a genuine higher holonomy class.

### MPA-34 — Decide `Omega+S2` at `n=8`

- **Delivered finite codec:** THM-801's expanded `Omega_7` tensor has 16,308
  cells for 16,384 lines; the mirror crossing-layer `S2` vector from the B2
  chart separates all 76 remaining double collisions.  The coarser `(S2,S3)`
  address leaves only 16 collisions inside 6,126 coloured node-pair fibres.
- **Must preserve:** ordered upper/three-face endpoint-node pairs, four-role
  colour, raw crossing-layer counts, literal loop multiplicity, and witnesses
  for every collision.
- **Delivered by THM-809:** the strictly weaker lower-only key `Lambda` is injective
  on all `1,048,576` lines.  The exact excess ladder is
  `418,252,148,74,52,0,0` as `tau=3,...,7,fixed` are adjoined.  Therefore the
  new size-three `tau=7` layer, not the upper node pair or fixed layer, closes
  `n=8`.  Among the 418 initial doubletons, the first-separator counts are
  `166,104,74,22,52`; none has both endpoint faces equal or contains an S11
  residual line as a face.  The next pull is the same layer-by-layer test at
  `n=9`, with first positional moments retained as the guaranteed repair when
  counts cease to reconstruct a layer.

### MPA-35 — Minimize the corewise Boolean Möbius stalk

- **Delivered exact definition:** for literal core `c`, the function
  `I_(u,c)(p_L,p_H,a)` is the node-fibre indicator on the two legs and apex.
  Its full Boolean Möbius transform is invertible.  THM-796's `D`, THM-801's
  `K_u`, `Omega`, `S2`, and `S3` are smaller marginals of this stalk.
- **Pull:** census coefficient support by bidegree
  `(|S cap L|,|S cap H|,1_(a in S))`, then delete degrees greedily while
  checking both static injectivity and one-step continuation equivalence.
- **Must preserve:** the literal core, endpoint-path orbit, sign/phase action,
  and coefficient ownership; averaging over a core node is not legal.
- **Deliverable:** the smallest exact stalk through `n=7`, its first
  non-lumpable truncation, and an all-size finite-degree conjecture or refuter.

### MPA-36 — Disintegrate electrical Möbius curvature over the flow

- **Delivered identities:** for the interval-face Möbius operator `Omega`,
  `q=Omega C3` counts cyclic triples meeting both path endpoints.  The Smith
  modes recover the full THM-785 cloud:
  `D(Omega J_h-Omega J_v)=epsilon` and
  `D(Omega J_h+Omega J_v)=lambda+(2a-1)(n-2)`.  Every black projected-fibre
  signed-`epsilon` histogram and every black `(q0,q1)` coefficient is
  reflection-symmetric/even.
- **Pull:** refine MPA-30's raw-versus-normalized drift by
  `(C3,q0,q1,|epsilon|,lambda,Smith potential,path stabilizer)`.  Determine
  which unsigned curvature stratum creates the black boundary imbalance.
- **Must preserve:** apex orientation, raw-bit convention, both line endpoints,
  literal multiplicity, and reflection pairing.  Signed epsilon alone cannot
  explain quotient drift.
- **Deliverable:** a conditional-current table with an orbit-stabilizer formula
  that predicts the sign after node-fibre disintegration.

- **Delivered master law by THM-811:** the four-variable polynomial joining
  `(Delta C3,q0,q1,epsilon)` factors into two boundary bits and `n-4` internal
  four-state packets.  It proves
  `Cov(q_i,epsilon^2)=-(n-4)/8`, strict negative covariance after conditioning
  on black, and the support bound
  `|epsilon|<=n-4-q0+1_(q1-q0=1)`.  The node orbit polynomial
  `(C3,H_x)` separates every node through `n=7`.  Curvature is
  not an edge codec: it gives 7,248/8,064 black orbits at `n=7`, while the
  reflection-orbit `(B2,B3)` address gives all 8,064.  Source-normalized
  outward black flow dominates in every source-`q` stratum carrying a
  mixed--pure-black boundary edge at `n=6,7`; the nonempty `n=6,q=2` source
  populations carry no such edge and give `0=0`.  The `n=8` active-stratum
  table is now the decisive conjecture test.

### MPA-37 — Give the gap face an intrinsic tournament semantics

- **Delivered distinction:** the middle `d_B:(a,b)->(a-1,b)` face is a valid
  lower tournament tiling but not induced deletion of one vertex.  Its support
  row resolves five of the eight `n=7` endpoint-face support twins, while its
  primitive row is weaker than endpoint deletion.
- **Pull:** characterize `d_B` as an operation on a tournament with a marked
  Hamiltonian path—path-edge contraction, interval-root shortening, or a
  functor in a path-minor category—and state exactly how automorphisms act.
- **Must preserve:** marked path positions and interval-root ownership; an
  unmarked tournament node is insufficient.
- **Deliverable:** an intrinsic path-category definition and the induced
  weighted correspondence on isomorphism/converse classes.

### MPA-38 — Couple continued-fraction transport to Möbius stalk coefficients

- **Delivered analogy made precise enough to test:** THM-778 says a continued-
  fraction digit is safe only together with its action on the labelled token
  fibre.  THM-801 says a node/face address is safe only together with the
  corewise leg/apex interaction stalk.
- **Pull:** let centered Christoffel substitutions act on the low/high leg
  variables and compute their pullback on Möbius coefficients.  Search for a
  finite invariant coefficient sector under the ten-wall Euclidean word.
- **Must preserve:** common gcd, midpoint phase, tie blocks, owner inverse
  steps, literal core, and the substitution action—not only the digit or
  coefficient spectrum.
- **Deliverable:** one exact transported-stalk automaton, or a counterexample
  proving that no bounded-degree sector is closed.

- **Delivered first action/no-go by THM-808:** for `p+1` tokens covering
  `F_p` with duplicate root `d`, a centered-Christoffel block with owner counts
  `c_a` acts by `d'=d-sum c_a w_a^(-1)`.  On THM-778's ten-wall movie the five
  block types act by translations `4,3,3,4,3 mod 7`.  The repeated three-wall
  block begins twice at mask `31115` but reaches masks `14635` and `615`, so
  `(mask,CF block)` is not a continuation state; the two source roots `4,6`
  supply the missing transported phase.

- **Delivered direct stalk action by THM-812:** the centered schedule
  `(q,p,s)=(3,2,1)` with increment word `(1,2)` gives an explicit coordinate-
  copy embedding `X_5->X_6`.  It is complement/reflection equivariant and
  maps all 20 projected coloured edge cells injectively, although 10 bare
  source nodes spread over 23 target nodes.  Boolean Mobius coefficients obey
  the general subset-image pushforward law.  Degree at most three has one
  fixed-core collision, resolved by exactly three quartic cross-leg/apex
  coefficients.  This supplies the first exact CF-to-metagraph edge action;
  the owner/root stalk must still be coupled before any LRC predicate is
  transported.

- **Delivered sharp action boundary by THM-813:** the next centered copy
  `X_6->X_7`, with word `(1,2,1)`, is still injective and complement/reflection
  equivariant but splits 51/187 projected coloured edge cells.  The first
  failure is the blue cell `B(15,16)`, and it persists across all 40 natural
  equivariant core lifts.  The correct all-size action carrier is
  `Q_n=E_n/<sigma>`: all 272 source reflection orbits act without failure,
  although their target projections need not be distinct.  Degree-four node
  indicators remain exact after composing the first two CF replications;
  arbitrary bounded-degree sectors require saturation under coordinate
  collapse.

### MPA-39 — Tutte/Smith polynomial of the staircase cover

- **Delivered electrical bridge:** the canonical staircase Smith network is a
  series chain of parallel bundles, while THM-801's primal/dual face
  curvatures recover `(lambda,epsilon)`.  The colour law is a single top
  Boolean interaction, not an Eisenstein imbalance.
- **Concurrent delivery:** THM-805 gives the unmarked closed form
  `T(N_n;x,y)=prod_(k=1)^(n-2)(x-1+[k]_y)` and its classical specializations.
  Its unit-circle and harmonic-measure bridges are logged separately under
  HYP-6885.
- **Pull:** mark this product by `A/B/C` cover membership and test whether its deletion-
  contraction coefficients recover the `q` edge polynomial, the blue
  self-dual specialization, or the even-graph/Tutte duality left dormant by
  the old kps-S11 thread.
- **Must preserve:** planar cell network versus class-level wiggly network;
  do not identify their resistances or spanning trees.
- **Deliverable:** a proved marked Tutte specialization or a clean no-go
  separating the two electrical levels.

## Recommended next three pulls

1. Join `MPA-32/34/35`: compare `n=8` continuation cells with THM-809's
   injective lower `Lambda` and its 418-to-zero layer genealogy.  Include
   inverse lifts and internal deletion, then test static injectivity at `n=9`
   and explain the first pre-separation or genuine collision by an endpoint-
   fixed interior permutation or the smallest missing corewise Möbius
   coefficient.
2. Join `MPA-30/36`: use the exact Smith/Möbius coordinates to explain black
   drift after orbit/fibre disintegration, where signed symmetry itself cannot.
3. Join `MPA-28/38`: extend THM-812's exact one-block action through the
   ten-wall centered-CF word, auditing fibre-purity separately on tilings,
   lines, coloured edge cells, and nodes; then form the literal-witness fibre
   product with THM-808's transported owner/root stalk.
