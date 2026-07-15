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
- **Must preserve:** path orbit, owner assignment, simultaneous-isomorphism
  line orbit, inverse winding, centered phase/ties, wall rank, and carry.
- **Deliverable:** the first exact three-level incidence index
  `node <- path orbit <- owner-labelled LRC state`, with forward maps,
  inverse fibres, and continuation-purity audit at every quotient.

## Recommended next three pulls

1. Join `MPA-18/28/29`: attach redundancy roots and Hamiltonian path-orbit
   identities to THM-778's ten owner-and-mask-labelled r=8 first returns.
2. Continue `MPA-12/26`: minimize the bounded prime-sheet/Euclidean movie and prove
   whether inverse steps plus cyclic endpoint word are continuation-complete.
3. Join `MPA-27/28`: seek a folded coordinate system in which each Euclidean
   substitution changes only boundedly many fibre fields.
