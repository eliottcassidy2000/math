# LRC Tournament Technique Index

**Purpose:** A living catalog of tournament, metagraph, series, and quotient
techniques that have proved useful, or plausibly useful, for the LRC14 proof
program.  The immediate use is practical: future agents should be able to pick
a proof route, see which coordinates it must keep, and avoid rediscovering the
same quotient failures.

**Relation to the broader index:** `00-navigation/LRC-TECHNIQUE-INDEX.md` is
the main shared LTI atlas.  This file is a tournament-specific companion with
longer guardrail notes and pull prompts for agents working directly with
tournaments, metagraphs, A000568-style quotients, and fixed-path tilings.

**Current proof shape:** the LRC14 problem is no longer just "find a scalar
that distinguishes AP/Goddyn-Wong from everything else."  The best current
frame is a labelled packet proof.  A primitive row should route to a q-witness,
AP/GW equality atom, C27 petal/two-block discharge, K33/F7 state lift,
covering/boundary-moment positivity, or a rigorous Fejer/Toeplitz certificate.
Every technique below is a possible tooth in that zipper.

## Contribution Protocol

When adding a technique card, keep these fields:

Use `LTT-*` handles in this companion file.  The broader
`00-navigation/LRC-TECHNIQUE-INDEX.md` owns the `LTI-*` namespace.

- **Move:** the abstract trick or transformation.
- **LRC use:** how it could help LRC14 specifically.
- **Preserves:** the predicate or coordinate that survives the quotient.
- **Forgets / guardrail:** what is dangerous to discard.
- **Next pull:** one concrete thing another agent can do.
- **Pointers:** hypotheses, theorems, tangents, scripts, variables, or
  reflections.

Do not default tournament vertices to runners.  Before choosing vertices,
consider runners, gaps, endpoint owners, fixed circle sections, section
boundaries, wall-crossing events, residue packets, cover arcs, Fourier modes,
Haar rectangles, matroid topes/cocircuits, Fejer atom banks, state-lift
obligations, and proof-carrier interfaces.

## Fast Routing Map

- Need to separate strict-open mass from endpoint-only equality:
  use LTT-018 through LTT-021, plus LTT-069.
- Need to classify large pair-good decoy tails without trusting raw counts:
  use LTT-068.
- Need to repair failed automatic/fiber quotients without using route labels:
  use LTT-072.
- Need a rigorous positive-row certificate:
  use LTT-022, LTT-023, LTT-024, and LTT-026.
- Need to prevent an unsafe quotient:
  use LTT-001, LTT-025, LTT-039, and LTT-040.
- Need tournament enumeration speedups:
  use LTT-006 through LTT-012, plus LTT-033, LTT-035, and LTT-046
  through LTT-050.
- Need the AP/GW residue skeleton:
  use LTT-003, LTT-013, LTT-014, LTT-027, and LTT-028.
- Need a route for non-AP/GW zero-open residuals:
  use LTT-029, LTT-030, LTT-031, LTT-040, LTT-051, and LTT-052.
- Need cross-domain inspiration without scalarizing it away:
  use LTT-036 through LTT-044, plus LTT-053 through LTT-058.
- Need the broadest creative pullback menu across tournaments, metagraphs,
  series, automata, harmonic certificates, topology, and formal tooling:
  use LTT-073 and `00-navigation/LRC-CARRIER-PULLBACK-INDEX.md`.

## Core Guardrail

A quotient may forget a coordinate only if at least one of these is true:

1. The target LRC predicate is constant on each fiber.
2. The forgotten coordinate is reconstructible from retained data.
3. A dual certificate annihilates the coordinate.
4. The coordinate is routed to a named residual sector.

If none of these hold, the quotient is not a theorem; it is only a diagnostic.

## Technique Cards

### LTT-001: Controlled-Kernel / Zipper Theorem

- **Move:** Treat every proof route as a zipper of controlled quotients, each
  with a declared kernel and residual sector.
- **LRC use:** Assemble q-witness, Fejer, Ramanujan, endpoint, C27, K33, and
  covering lanes without forcing them into one scalar hierarchy.
- **Preserves:** LRC predicate, exact scale, packet family, endpoint ownership,
  dual certificate, and residual route.
- **Forgets / guardrail:** A quotient that forgets phase, scale, or ownership
  must prove fiber-constancy or emit a named residual.
- **Next pull:** Convert HYP-2987's F7 bucket into a formal named
  harmonic/state-lift residual sector.
- **Pointers:** HYP-2990, HYP-2987, HYP-2984, HYP-2985, HYP-2986, T1074.

### LTT-002: Proof-Carrier Tournament Analysis

- **Move:** Make proof carriers the tournament vertices rather than runners,
  then orient edges by retention, certifiability, or discharge strength.
- **LRC use:** Compare Fejer packets, Ramanujan packets, endpoint bridges,
  Kaczynski approaches, and state lifts as competing but composable carriers.
- **Preserves:** The proof obligation and the handoff relation.
- **Forgets / guardrail:** It can hide runner geometry unless the carrier card
  records which LRC predicate it certifies.
- **Next pull:** Add a standard fingerprint payload to each future LRC script:
  score histogram, directed 3-cycles, SCC sizes, and Hamiltonian path count.
- **Pointers:** AGENTS.md, HYP-2987, HYP-2990, CONCEPT-MAP "Tournament
  Analysis".

### LTT-003: Multi-Scale Winding Tournament Spectrum

- **Move:** Replace a single phase tournament by the finite weighted spectrum
  of isomorphism classes swept over all phase intervals.
- **LRC use:** Repairs magnitude-blind apex tournaments: AP and loose rows can
  share residues but differ in spectrum size, peak class, and binding scale.
- **Preserves:** Phase breakpoints, binding denominator, and class measure.
- **Forgets / guardrail:** A single apex class forgets magnitude and cannot
  distinguish AP from loose residue twins.
- **Next pull:** Compute spectrum fingerprints for the latest HYP-2963 packet
  bank and compare binding-scale migration to Fejer certificate degree.
- **Pointers:** HYP-2928, HYP-2927, HYP-2605, THM-568, OPEN-Q-108.

### LTT-004: Observer-Source Marked Tournament

- **Move:** Add the stationary observer as a distinguished vertex; orient
  observer-to-runner edges by whether a runner is outside the danger moat.
- **LRC use:** Loneliness is exactly the observer being a source; blocker count
  is observer indegree.
- **Preserves:** The LRC threshold predicate at a fixed moment.
- **Forgets / guardrail:** Unmarked tournament class is too coarse; it mixes
  safe and unsafe moments.
- **Next pull:** Attach endpoint-owner and exact-period labels to the
  observer-source cone for LRC14 hard packets.
- **Pointers:** THM-381, THM-385, HYP-1988, HYP-2486, S511, S517.

### LTT-005: A000568 Source Fiber / Perspective Lift

- **Move:** Project LRC states to source-deleted tournament classes and then to
  rooted perspective fibers over A000568.
- **LRC use:** Gives an exact small-n source slice and a warning that unrooted
  classes are not faithful.
- **Preserves:** Deleted runner class plus source-cone operation.
- **Forgets / guardrail:** Full threshold payload and incident word can be lost
  after unrooting.
- **Next pull:** Use source-fiber IDs as packet-family keys in LRC14 labelled
  classifier outputs.
- **Pointers:** HYP-2486, HYP-2120, HYP-2121, HYP-1977, A000568 scripts.

### LTT-006: Fixed Hamiltonian Path / Staircase Tiling

- **Move:** Fix a Hamiltonian path as the tie spine and encode all other arcs
  as tiles above or below it.
- **LRC use:** Supplies a canonical address system for local switches,
  good-cuts, metagraph moves, and half-tiling parity.
- **Preserves:** Path order, interval support, and tile range.
- **Forgets / guardrail:** The fixed path is arbitrary unless its role is
  declared as a tie Hamiltonian path, not a hidden theorem.
- **Next pull:** Build a direct dictionary from LRC endpoint-wall order to
  staircase tile intervals.
- **Pointers:** fixed-path tiling threads, good-cut variables, T1072, HYP-2992.

### LTT-007: Tournament Metagraph / SC-NS Spine-Ribs-Sea

- **Move:** Study the graph whose vertices are tournament isomorphism classes
  and whose edges are controlled flips, complement moves, or merged classes.
- **LRC use:** Gives a language for local deformation, regularity, and
  self-complementary symmetry; helps ask whether a packet is on a spine,
  rib, or sea.
- **Preserves:** Class adjacency and complement/merge structure.
- **Forgets / guardrail:** Raw blue/black visible edges are not automatically
  even/odd; parity usually lives in a boundary or line-addressed chain.
- **Next pull:** Place LRC14 packet-sign classes on the metagraph and measure
  whether hard rows sit near SC or high-good-cut regions.
- **Pointers:** COMPLETE_INVARIANT_CATALOG, metagraph parity S675, A000568
  census work.

### LTT-008: Even-Graph Projection / Projection Defect

- **Move:** Use the equinumerosity of tournaments and even graphs as a second
  projection lens, then measure where tournament and even-graph projections
  disagree.
- **LRC use:** Defects identify coordinates that one projection forgets, a
  close analogue of LRC quotient leakage.
- **Preserves:** Edge-flip parity and even-graph class.
- **Forgets / guardrail:** Even-only and tournament-only defects have different
  polarity by tile range; do not collapse them.
- **Next pull:** Test whether AP/GW boundary packets are pure under even-graph
  projection while K33/petal packets show projection defects.
- **Pointers:** projection-defect leads, HYP-2245, HYP-2187, S674b, S675.

### LTT-009: Good-Cut / SCC Support Coordinate

- **Move:** Count cuts of the fixed path crossed by at least one upward tile;
  structurally this equals `n - #SCC(T)`.
- **LRC use:** A compact support residue for how much strong connectivity a
  proof-carrier tournament has.
- **Preserves:** SCC condensation boundaries along the fixed path.
- **Forgets / guardrail:** Bucket 1 is impossible for interval reasons; other
  gaps must come from tournament structure, not interval support.
- **Next pull:** Add good-cut/SCC profiles to packet-sign and spectrum
  tournaments for named LRC14 rows.
- **Pointers:** THM-354, good-cut-count variable, INV-237.

### LTT-010: Hamiltonian Path Count H / OCF

- **Move:** Use Hamiltonian path count, OCF activities, and Moon strong atom
  products as a scalar shadow of tournament structure.
- **LRC use:** H can flag spread, regularity, and phase entropy, but only after
  threshold and observer labels are attached.
- **Preserves:** Directed path richness and some strong-component structure.
- **Forgets / guardrail:** H alone is not a max-gap meter; apex H is
  magnitude-blind.
- **Next pull:** Compare H distributions of LRC spectrum classes at the
  binding scale, not just at the apex phase.
- **Pointers:** THM-002, OCF variables, HYP-1970, HYP-2922 through HYP-2928.

### LTT-011: W-Polynomial / Walsh-Fourier Channels

- **Move:** Encode Hamiltonian path data in W-polynomial and Walsh/Fourier
  channel decompositions.
- **LRC use:** Separates even and odd channels, useful for midpoint/reversal
  gates and packet-sign symmetries.
- **Preserves:** Activity/Fourier channel data above raw H.
- **Forgets / guardrail:** Evaluation at one scalar can erase the channel that
  distinguishes AP/GW from nearby deformations.
- **Next pull:** Pair W-polynomial even/odd channel tests with Haar-product
  packet coefficients.
- **Pointers:** W-polynomial variable, OCF Fourier work, HYP-2992, HYP-2989.

### LTT-012: Burnside / CRT / A000568 Enumeration

- **Move:** Count tournament classes by group action, CRT factorization, and
  odd-divisor profiles.
- **LRC use:** Speeds exact enumeration of packet fibers and helps decide
  which quotient layers are feasible to exhaust.
- **Preserves:** Orbit structure under relabelling or complement.
- **Forgets / guardrail:** Orbit representatives need rooted/source/endpoint
  decorations before being used for LRC.
- **Next pull:** Build an LRC packet enumerator that canonicalizes decorated
  source-cones rather than raw tournaments.
- **Pointers:** A000568 scripts, Burnside/CRT reflections, HYP-2486.

### LTT-013: Three-State Tournament Automaton

- **Move:** Coarsen pair states into left, right, and middle, then propagate a
  small automaton instead of all phase data.
- **LRC use:** Models safe side, danger side, and boundary/tie as a minimal
  carrier for wall crossing.
- **Preserves:** Directional side and middle/boundary state.
- **Forgets / guardrail:** Middle can hide endpoint ownership and exact-period
  denominator unless labelled.
- **Next pull:** Implement a 3-state automaton over endpoint-wall intervals for
  AP/GW, K33, and petal rows.
- **Pointers:** tournament automata prompts, observer-source work, HYP-2986.

### LTT-014: Score-Sequence / Landau Prefix Majorization

- **Move:** Compare cumulative score prefixes rather than pointwise winners.
- **LRC use:** Mirrors the low-denominator prefix-majorization evidence: a row
  may beat AP in one bucket but pay mass earlier.
- **Preserves:** Aggregate dominance constraints.
- **Forgets / guardrail:** Pointwise q-bucket domination is false; use prefix
  ledgers and Farey tails.
- **Next pull:** Formalize low-q prefix majorization over residue patterns and
  combine with a three-gap/Farey tail packing lemma.
- **Pointers:** T1015, HYP-2866, THM-530, THM-531, HYP-2870.

### LTT-015: Packet-Sign Tournaments / K4 Orientations

- **Move:** Convert signs of packet interactions into small tournaments, often
  K4 or signed square cells.
- **LRC use:** Detects hidden orientation in four-term structures and
  distinguishes diagonal from anti-diagonal fixed-margin packets.
- **Preserves:** Local sign/cocircuit orientation.
- **Forgets / guardrail:** Three-term folds do not see every four-term hidden
  deformation.
- **Next pull:** Generate packet-sign K4 fingerprints for HYP-2963 packet
  families and compare with Fejer certificate margins.
- **Pointers:** HYP-2989, HYP-2992, HYP-2595, summand graph threads.

### LTT-016: Sector-State Transfer DP

- **Move:** Track runner insertion or packet migration by finite sector states
  rather than continuous phase.
- **LRC use:** A possible efficient checker for "does this orbit hit the safe
  box once" after endpoint-wall compression.
- **Preserves:** Section occupancy, boundary debt, and transition law.
- **Forgets / guardrail:** A sector DP must record exact endpoint ties or it
  will misclassify AP/GW equality.
- **Next pull:** Build a DP whose states are `(section occupancy, endpoint
  owners, exact-period class, source deficit)`.
- **Pointers:** LRC section functors, HYP-2570, HYP-2024, HYP-2949.

### LTT-017: Taut Bridge Graph

- **Move:** Treat positive safe intervals as directed bridges between endpoint
  owner labels; boundary-only equality points are zero-length taut vertices.
- **LRC use:** Separates AP/GW zero-open equality atoms from positive-open
  non-AP/GW rows.
- **Preserves:** Endpoint owners, bridge direction, and point-depth.
- **Forgets / guardrail:** Safe measure alone forgets owner-current and
  zero-sum boundary skeletons.
- **Next pull:** Prove boundary-owner skeleton rigidity: no positive bridge and
  AP/GW owner-current pattern implies AP/GW.
- **Pointers:** HYP-2975, HYP-2951, HYP-2949, T1043.

### LTT-018: Tope / Cocircuit Endpoint Arrangement

- **Move:** Cut the time circle by danger endpoints and regard open cells as
  topes, boundary equality atoms as cocircuits.
- **LRC use:** Rephrases a strict counterexample as a no-tope/no-cocircuit
  forbidden wall packet.
- **Preserves:** Exact endpoint arrangement and boundary owners.
- **Forgets / guardrail:** Residue-only tournaments miss cocircuit ownership.
- **Next pull:** Compare cocircuit owner-pair sums with Haar same-tile boundary
  atoms from HYP-2992.
- **Pointers:** HYP-2986, T1070, HYP-2975, HYP-2951.

### LTT-019: Haar-Baire Strict-Open Carrier

- **Move:** Split safe sets into regular-open Haar mass, Baire interior, and
  finite closed boundary support.
- **LRC use:** AP/GW are boundary-only; K33, petals, splices, and covering rows
  have positive strict-open mass in audited banks.
- **Preserves:** Open-vs-boundary topology and measure.
- **Forgets / guardrail:** Haar zero does not identify ownership unless
  boundary support is labelled.
- **Next pull:** Prove every strict-Haar-zero row after reductions has AP/GW
  boundary-owner skeleton.
- **Pointers:** HYP-2948, HYP-2949, HYP-2951, T1042.

### LTT-020: Haar Product Square / Fixed-Margin Switch

- **Move:** Use the elementary identity `h_I(x)h_J(y) = [[1,-1],[-1,1]]`
  on dyadic children.
- **LRC use:** It is the same local algebra as a 2-by-2 fixed-margin switch in
  tournament tilings.
- **Preserves:** Mixed product coefficient and row/column margins.
- **Forgets / guardrail:** Raw component count and margins alone hide the
  diagonal versus anti-diagonal sign.
- **Next pull:** Replace raw discrepancy component count by independent
  color-compatible mixed Haar switches.
- **Pointers:** HYP-2989, T1073, HYP-2594, HYP-2595.

### LTT-021: Haar Rectangle Tile Atlas

- **Move:** Enumerate products of dyadic Haar rectangles: orthogonal zero,
  same-tile atom, owner strip, cross handoff, and nested refinement.
- **LRC use:** Candidate vanishing lemma for zero-open residuals: some typed
  coefficient must survive unless the packet is AP/GW or a state-lift atom.
- **Preserves:** Two-coordinate interaction type and sign balance.
- **Forgets / guardrail:** Sign-balanced classes need endpoint and packet
  labels before they become proof data.
- **Next pull:** Compute typed Haar coefficients for named rows and cluster
  Fejer atoms by owner-strip/cross/nested class.
- **Pointers:** HYP-2992, T1072, HYP-2989, HYP-2988.

### LTT-022: Fejer / Toeplitz PSD Dual Certificate

- **Move:** Use trigonometric Fejer/Toeplitz forms as rigorous interval
  certificates on labelled packet fibers.
- **LRC use:** Floating audits certify every positive HYP-2963 row by degree
  at most 280; AP/GW remain the PSD-blind equality atoms.
- **Preserves:** Packet key, rational center, degree, atom bank, interval sign,
  and route handoff.
- **Forgets / guardrail:** A floating scalar evaluation is not enough; the
  interval certificate must be anchored to packet fibers.
- **Next pull:** Build a formal interval backend and lift selected row
  certificates to family templates.
- **Pointers:** HYP-2981, HYP-2974, T1065, Robbins/Robin guardrails.

### LTT-023: Multiplicity Moment Dual Barriers

- **Move:** Use danger-count moments, Bonferroni/Newton layers, and activation
  depth to bound covering pressure.
- **LRC use:** Handles covering/boundary-moment residuals after q-witness and
  packet routes have removed scalar-easy rows.
- **Preserves:** Depth parity and moment layer.
- **Forgets / guardrail:** Third-order truncation fails without activation and
  missing-depth parity ledgers.
- **Next pull:** Re-index moment barriers by first live Newton layer and packet
  family, then compare to Fejer degree.
- **Pointers:** HYP-2901, HYP-2903, T1016/T1017, HYP-2963.

### LTT-024: Ramanujan Exact-Period Projectors

- **Move:** Use Ramanujan sums and exact-period projectors to split
  denominator packets before analytic or Fejer certification.
- **LRC use:** Separates q=14, q=27, q=41, prime-power, and lcm-tail packets
  that scalar divisor functions mix.
- **Preserves:** Exact period and primitive root trace.
- **Forgets / guardrail:** Squarefree weights such as `mu^2/phi` erase
  prime-power packets like 27, 63, 84, 168, and 4312.
- **Next pull:** Add endpoint-owner Ramanujan profiles for q in `{14,27,41}`
  and the Fejer hard denominators.
- **Pointers:** HYP-2979, HYP-2978, T1062/T1063, HYP-2982.

### LTT-025: Divisor / Totient / Mobius Quotient Guardrails

- **Move:** Use divisor-function identities as packet laws, not final scalar
  invariants.
- **LRC use:** Tests which arithmetic quotients mix AP/GW, q-witness, K33,
  petal, splice, and covering routes.
- **Preserves:** Dirichlet convolution, unit orbit, Jordan capacity, and
  Ramanujan trace when labelled.
- **Forgets / guardrail:** qdiv, mod-14 residue, unit count, and lcm scalar
  quotients all mix routes in audits.
- **Next pull:** Turn the controlled-kernel criterion into a reusable quotient
  admissibility lemma for LRC scripts.
- **Pointers:** HYP-2978, HYP-2900, divisor-function web crawl, T1062.

### LTT-026: Analytic Sieve / Kaczynski Smoothing

- **Move:** Treat large sieve, Selberg weights, circle-method splits,
  exponential sums, and smoothing as middle certificates with explicit
  boundary approach classes.
- **LRC use:** Potentially handles late-denominator and true-wide packets that
  are beyond finite local atlases.
- **Preserves:** Kernel, smoothing transform, major/minor packet, exceptional
  set, and approach class.
- **Forgets / guardrail:** `Phi` capacity and `G=sum mu^2/phi` normalizers are
  scalar shadows unless exact-period and endpoint labels are restored.
- **Next pull:** Build an "explicit explicit formula" emitter for LRC atom
  banks with interval trig provenance and smoothing family tags.
- **Pointers:** HYP-2982, HYP-2983, HYP-2985, THM-548, T1066/T1069.

### LTT-027: Farey Mutation Scheduler / Kpq Wall

- **Move:** Use Farey neighbors and mutations of `p/q` as a branch scheduler:
  product, sum, denominator power, numerator power.
- **LRC use:** Splits `2/27` as the C27 petal branch and `3/41` as the first
  K33 wall beyond C27.
- **Preserves:** Exact M/Farey branch and product-incidence packet.
- **Forgets / guardrail:** Product value is not order-safe by itself; exact M
  and route labels come first.
- **Next pull:** Use Farey branch to pre-route packet families before Fejer or
  Ramanujan certificate generation.
- **Pointers:** HYP-2932 through HYP-2945, HYP-2984, T1030/T1036/T1041.

### LTT-028: C27 Shell Transfer / Unital Pair Completion

- **Move:** Work on antipodal shell pairs `{a,27-a}` and branch-local unital
  block completions.
- **LRC use:** Explains GW, K33 near-miss, and p=2 petal rows as typed shell
  transfers.
- **Preserves:** Hole/double ownership, unit/nonunit visibility, and branch
  local pair incidence.
- **Forgets / guardrail:** Global unital blocks can fail while branch-local
  charts remain useful; exact M/Farey branch must stay attached.
- **Next pull:** Formalize the p=2 petal/two-block discharge before invoking
  broader K33/F7 state lifts.
- **Pointers:** HYP-2937, HYP-2940, HYP-2942, HYP-2947, T1038.

### LTT-029: K33 / HYP-2908 State-Lift Residual

- **Move:** Treat K33 incidence and F7/HYP-2908 state lift as a named residual,
  not an error bucket.
- **LRC use:** The first p>=3 Farey child `12->36` and K33 splices route here.
- **Preserves:** Incidence wall, state-lift debt, and harmonic residual sector.
- **Forgets / guardrail:** Calling a row "positive" without explaining the
  state-lift address loses the main obstruction.
- **Next pull:** Define the F7 residual as a concrete finite tournament/state
  construction from zero-open non-AP/GW packets.
- **Pointers:** HYP-2908, THM-572, HYP-2987, HYP-2990, HYP-2945.

### LTT-030: Boundary-Moment Covering Residual

- **Move:** Route covering rows by boundary-moment charts and cap-cell pressure
  instead of brute force.
- **LRC use:** Discharges rows with low packet mass or late covering pressure,
  such as `12->168`-style families.
- **Preserves:** Covering cell, moment chart, endpoint pressure, and cap face.
- **Forgets / guardrail:** A single denominator chart is not an obstruction;
  covering needs multichart labels.
- **Next pull:** Tie cap-cell pressure to Fejer packet certificate degrees and
  Ramanujan exact-period projectors.
- **Pointers:** HYP-2963, HYP-2961, THM-398, S593, HYP-2976.

### LTT-031: Few-Apex / Exposure Poset

- **Move:** Order packet routes by exposure channels and check for hidden
  kernel rows not seen by Fejer or boundary atoms.
- **LRC use:** Current HYP-2988 audit has AP/GW as zero-safe atoms and no
  unexposed positive rows in the default bank.
- **Preserves:** Route, exposure channel, and named discharge.
- **Forgets / guardrail:** Exposure is only useful if every channel has a
  theorem-facing certificate, not only a label.
- **Next pull:** Prove no new normalized `p_0=0` wall exists in the Res_27
  carry/owner fiber.
- **Pointers:** HYP-2988, HYP-2963, HYP-2974, HYP-2981.

### LTT-032: Euler / Glaisher / Even-Graph Parity Carrier

- **Move:** Read even/odd partition and graph parity identities as quotient
  tests for which channel survives.
- **LRC use:** Helps distinguish reversal gates, midpoint gates, and parity
  address layers.
- **Preserves:** Parity channel or boundary over GF(2).
- **Forgets / guardrail:** Raw "black is even, blue is odd" visible-edge
  slogans are false; parity lives in the right chain complex.
- **Next pull:** Build parity-boundary tests for packet-sign metagraph layers.
- **Pointers:** S675 blue/black parity, HYP-2187, HYP-2245, Euler/Glaisher
  threads.

### LTT-033: Tiling Bucket Balance / Transport Matrix

- **Move:** Use quotient transport identities such as
  `2*self + incident_cross = |fiber|*|M|`.
- **LRC use:** Gives exact conservation laws for moving between packet fibers
  and metagraph classes.
- **Preserves:** Fiber size, incident transport, and balance row.
- **Forgets / guardrail:** Transport matrices need class labels; scalar totals
  can hide biased movement.
- **Next pull:** Make a transport matrix for HYP-2963 packet routes under
  one-swap and two-swap moves.
- **Pointers:** tiling-bucket-balance variable, metagraph transport results.

### LTT-034: Spectral Shadow / Trace / Zeta Gas

- **Move:** Use traces, cycle counts, Ihara/zeta analogues, and spectral gases
  as high-level shadows of tournament structure.
- **LRC use:** May find cheap prefilters for packet families before exact
  interval work.
- **Preserves:** Cycle/correlation statistics.
- **Forgets / guardrail:** Trace counts fail at overlap corrections; spectral
  shadows are not final certificates.
- **Next pull:** Test whether Fejer-hard packet families have distinctive
  trace-overlap or zeta signatures.
- **Pointers:** HYP-2498, tournament trace speedups, Ihara zeta scripts.

### LTT-035: Good-Cut Bucket Polynomial / Interval Gas

- **Move:** Count good-cut buckets by connected-run covers and recurrence
  rather than enumerating all tilings.
- **LRC use:** A model for replacing expensive packet enumeration with a gas
  recurrence over interval support.
- **Preserves:** Good-cut support count and run-cover structure.
- **Forgets / guardrail:** Support recurrence does not know endpoint-owner
  arithmetic by itself.
- **Next pull:** Derive a packet-support gas for endpoint walls or Haar
  rectangles and compare to exact bank scans.
- **Pointers:** good-cut-bucket-polynomial variable, INV-237, THM-349.

### LTT-036: Half-Tiling Parity Address

- **Move:** Split tiling moves by half-plane, anti-diagonal, or hypotenuse
  parity addresses.
- **LRC use:** Useful when AP/GW equality lives on a boundary and positive
  rows move into one side of a symmetric wall.
- **Preserves:** Side address and parity of crossing.
- **Forgets / guardrail:** Side labels are coordinate choices unless tied to an
  endpoint or cocircuit.
- **Next pull:** Pair half-tiling parity with AP/GW boundary owner skeletons.
- **Pointers:** half-tiling/parity tangents, HYP-2986, HYP-2992.

### LTT-037: Unit-Distance Spine / Ear Analogy

- **Move:** Translate point-set constraints into tournament edges or ears,
  asking whether mandatory Hamiltonian spines lie in unit or nonunit pairs.
- **LRC use:** Cross-domain analogy for when a local mandatory path changes
  carrier family, like AP/GW to K33/petal.
- **Preserves:** Spine/ear incidence and unit/nonunit edge status.
- **Forgets / guardrail:** Unit-distance geometry is not LRC unless the
  preserved predicate is explicitly mapped.
- **Next pull:** Use spine-flop language to name when an LRC packet leaves AP
  and first enters K33 or petal branch.
- **Pointers:** unit-distance H=21/u(21) threads, HYP-2298, HYP-2990.

### LTT-038: Octahedral Current / Hodge Curl Support

- **Move:** Decompose local support into divergence, curl, and harmonic
  current packets.
- **LRC use:** Candidate language for F7 residuals and for detecting hidden
  state-lift currents after scalar certificates fail.
- **Preserves:** Local conservation law and harmonic residue.
- **Forgets / guardrail:** A current without packet route is just a metaphor;
  it must attach to a wall or endpoint arrangement.
- **Next pull:** Express K33/C27 handoffs as divergence-free or curl-bearing
  support packets.
- **Pointers:** octahedral current threads, HYP-2990, HYP-2969.

### LTT-039: OCF Activity / Coimage / Noise Stability

- **Move:** Use activity polynomials, coimage decompositions, and noise
  stability as a higher-resolution shadow than H.
- **LRC use:** May detect when a packet family has the right support but the
  wrong activity address.
- **Preserves:** Activity sector and coimage map.
- **Forgets / guardrail:** Coimage scalarization can erase observer/source
  marks and endpoint ownership.
- **Next pull:** Compute OCF activity fingerprints for LRC spectrum classes at
  binding scales 14, 27, and 41.
- **Pointers:** OCF variables, THM-002, HYP-2990, HYP-2486.

### LTT-040: Labelled Packet Classifier / Source-Spectrum Unknown

- **Move:** Make every row classify into q-witness, AP/GW, petal, K33,
  covering moment, Fejer, or unknown source-spectrum route.
- **LRC use:** Turns the global proof goal into emptiness of a named unknown
  bucket.
- **Preserves:** Route, family, q-class, packet state, state lift, and
  threshold.
- **Forgets / guardrail:** The classifier is only theorem-facing if each label
  has a proof obligation attached.
- **Next pull:** Expand the classifier's packet key to include Haar tile class,
  Ramanujan projector, and spectrum binding scale.
- **Pointers:** HYP-2963, HYP-2976, HYP-2987, HYP-2988.

### LTT-041: Midpoint Scalar Gate / Odd-Channel Faulhaber

- **Move:** Center a balance equation at its midpoint so only odd Faulhaber
  moments survive; contrast with reversal gates where even channels survive.
- **LRC use:** Gives a clean model for deciding which channel a quotient is
  allowed to keep under symmetry.
- **Preserves:** Odd/even channel selected by the gate.
- **Forgets / guardrail:** The channel is gate-dependent; do not transfer an
  odd-channel proof through a reversal gate unchanged.
- **Next pull:** Recast AP/GW boundary equations in a midpoint/reversal channel
  table and compare to W-polynomial parity.
- **Pointers:** Faulhaber tower prompts, HYP-2992, W-polynomial variable.

### LTT-042: Prime / Irreducible Polynomial Tournament Bridge

- **Move:** Model reducibility as a convolution-lift tiling problem, or orient
  residue/sign/factor-capture choices as tournaments or hypertournaments.
- **LRC use:** Cross-trains the habit of keeping hidden lift coordinates:
  irreducibility means no coefficient convolution lift exists, like no hidden
  LRC packet lift may remain.
- **Preserves:** Factor allocation, residue local obstruction, and hidden lift.
- **Forgets / guardrail:** Coefficient signs or one witness value cannot decide
  global factor structure without the lift.
- **Next pull:** Borrow the convolution-lift infeasibility language for LRC
  packet handoffs: what hidden 2D lift would make a zero-open residual possible?
- **Pointers:** irreducibility/Bunyakovsky prompts, Singh witness prefilter,
  HYP-2990.

### LTT-043: CRT / p-Adic Residual Packet Tree

- **Move:** Split rows by prime-power carries, residue packets, and CRT
  recombination rather than by denominator alone.
- **LRC use:** The n=14 wall is `2*7`; q=27/C27 and late prime-power packets
  show that prime-power side channels matter.
- **Preserves:** Carry owner, shell height, prime-power period, and recombined
  CRT address.
- **Forgets / guardrail:** Squarefree or denominator-only quotients hide shell
  collapse and owner-carry fibers.
- **Next pull:** Build a p-adic packet tree for q in `{14,27,41,63,84,168}`
  with owner labels and route outcomes.
- **Pointers:** THM-568 correction, HYP-2929, HYP-2979, HYP-2982, Res_27 work.

### LTT-044: Robbins No-Bridge Assembly

- **Move:** Import the graph-theoretic Robbins lesson: a strong orientation is
  possible exactly when no bridge is forgotten.
- **LRC use:** A proof assembly is safe only if every bridge between atom bank,
  center, packet, interval certificate, and route handoff is retained or
  discharged.
- **Preserves:** Assembly bridges between local certificate pieces.
- **Forgets / guardrail:** A quotient that forgets the bridge can turn a true
  interval certificate into a floating numerical hint.
- **Next pull:** Add a "bridge checklist" to Fejer manifest generation and to
  every future packet-family certificate.
- **Pointers:** HYP-2981, T1065, Robbins theorem notes, HYP-2990.

### LTT-045: Tournament Fingerprint Payload

- **Move:** Standardize lightweight tournament diagnostics across research
  scripts.
- **LRC use:** Makes exploratory scripts comparable: even when the proof route
  fails, the fingerprint can reveal a reusable carrier.
- **Preserves:** Score histogram, SCC sizes, directed 3-cycles, Hamiltonian
  path count, complement behavior, and chosen tie path.
- **Forgets / guardrail:** A fingerprint is not a proof unless its preserved
  predicate is declared.
- **Next pull:** Create a small shared helper for score histogram, SCCs,
  3-cycles, H count for small carriers, and Hamiltonian tie path reporting.
- **Pointers:** AGENTS.md Tournament Analysis default, HYP-2987, HYP-2990,
  T1074.

### LTT-046: Deletion-Contraction Summand Depth

- **Move:** Use the Redei/Mitrovic deletion-contraction identity
  `H(D)=H(D\e)+H(D/e)` as a metagraph recursion, then record the depth at
  which a Hamiltonian-path value appears in the summand graph.
- **LRC use:** Gives a proof-debt metric for packet reduction: if a labelled
  packet can only be discharged after several contractions, those lost boundary
  labels must become explicit coordinates.
- **Preserves:** Deletion edge, contraction endpoint, summand pair, SCC/product
  factor, and irreducible versus product source.
- **Forgets / guardrail:** One summand is usually a near-tournament, not a
  tournament.  Dropping the missing-arc or contracted-boundary label destroys
  the theorem-facing predicate.
- **Next pull:** Build a deletion-contraction depth ledger for the F0-F7
  labelled packet classifier and mark which branches are genuine tournaments.
- **Pointers:** THM-082, `07-reflections/metagraph-summand-recursive.md`,
  LTT-007, LTT-040.

### LTT-047: Delta-H Flip Energy Ledger

- **Move:** Treat an arc flip as a local energy move with
  `Delta H=H(T/e)-H(T'/e')`; old scans show these deltas are even and carry
  rich distribution data.
- **LRC use:** Packet handoffs, Haar switches, and endpoint-owner transfers
  can be given an "energy" comparable across rows instead of being collapsed
  to a safe/unsafe bit.
- **Preserves:** Flipped pair, contracted residue, `2`-adic valuation of the
  delta, level edge, and before/after fingerprint.
- **Forgets / guardrail:** `Delta H` is only an H-shadow.  It does not know the
  LRC endpoint owner, Farey branch, or exact-period label unless these are
  stored beside the flip.
- **Next pull:** Compute the flip-energy analogue for one-swap LRC packet
  transitions: AP to GW, GW to K33, C27 petals, and Fejer-hard rows.
- **Pointers:** `07-reflections/metagraph-summand-recursive.md`, LTT-018,
  LTT-033, LTT-045.

### LTT-048: Tiling Count / Hamiltonian-Path Overlap Partition

- **Move:** Use the fixed-path theorem
  `tiling_count([T])=H(T)/|Aut(T)|`, the master equation
  `2^m=sum_[T] H(T)/|Aut(T)|`, and Hamiltonian-path overlap partitions as
  normalization checks.
- **LRC use:** Prevents overcounting when packet fibers are represented by
  fixed-path staircase tilings, and gives an automorphism-aware denominator for
  packet census claims.
- **Preserves:** Automorphism size, fixed Hamiltonian base path, tiling fiber,
  and path-overlap partition block.
- **Forgets / guardrail:** Tiling count is not an LRC certificate by itself;
  observer/source marks and packet labels must survive the quotient.
- **Next pull:** Add `aut_size`, `H`, `H/aut`, and fixed-path-overlap fields to
  the HYP-2963 representative packet schema.
- **Pointers:** `05-knowledge/results/tiling_count_theorem.md`,
  `05-knowledge/results/master_identities.md`, LTT-006, LTT-011.

### LTT-049: Burnside Perturbation / Orbit-Cost Expansion

- **Move:** Read Burnside enumeration as a statistical-mechanics expansion:
  the identity orbit is the vacuum, non-identity cycle types are excitations,
  and each symmetry has positive cost at inverse temperature `log(2)`.
- **LRC use:** Separates generic labelled packet fibers from symmetric
  exceptional atoms.  AP/GW-like equality should be stabilizer-supported, not
  an average over a large anonymous orbit.
- **Preserves:** Stabilizer, conjugacy/cycle type, orbit cost, identity
  dominance term, and non-identity correction.
- **Forgets / guardrail:** Asymptotic orbit dominance is not a local proof.
  Small symmetric packets can be exactly the dangerous rows.
- **Next pull:** Add a Burnside cost column to packet fibers: identity mass,
  stabilizer size, first non-identity cost, and whether the row is exceptional
  because of symmetry.
- **Pointers:** `07-reflections/burnside-perturbation-theory.md`, LTT-012,
  LTT-048.

### LTT-050: Score-Class H-Spread / Magic Measure

- **Move:** For a score sequence `s`, record
  `spread(s)=max H-min H` over tournaments with that score.  Spread zero means
  H is score-determined; positive spread means hidden cycle-space content.
- **LRC use:** Score or ranker quotients can be used as cheap prefilters only
  after checking whether the relevant score fiber is a stabilizer class or a
  magic class.
- **Preserves:** Score sequence, within-score H range, stabilizer/magic flag,
  and extremal witnesses.
- **Forgets / guardrail:** Score sequence alone misses the LRC-active
  cycle-space exactly when spread is positive.
- **Next pull:** Compute H-spread for tournament shadows of AP, GW, K33,
  petals, and weakest Fejer-margin packets.
- **Pointers:** `05-knowledge/results/tournament_magic_measure_Hspread_kps.md`,
  LTT-014, LTT-039.

### LTT-051: Path-Homology Betti Carrier

- **Move:** Use GLMY path homology of tournaments as a residual-topology
  detector: `Omega_2` is transitive triples, tested tournaments have
  `beta_2=0`, and `beta_1`/`beta_3` phases appear mutually exclusive in
  audits.
- **LRC use:** A zero-open non-AP/GW packet might be forced into a directed
  topological hole; Betti payloads can distinguish a harmless scalar shadow
  from a genuine residual class.
- **Preserves:** Path-homology phase, transitive-triple structure, complement
  behavior, and OCF/topology correlation.
- **Forgets / guardrail:** Betti values are not an LRC proof unless attached
  to labelled packet routes and state-lift obligations.
- **Next pull:** Run path-homology fingerprints on the carrier tournaments
  created by HYP-2963 packet classes and HYP-2995 cocycle carriers.
- **Pointers:** `05-knowledge/results/path_homology_synthesis.md`, HYP-301,
  HYP-302, HYP-303, LTT-038, LTT-039.

### LTT-052: Permanent H-Gap Obstruction Grammar

- **Move:** Treat the permanent H-gaps `{7,21}` as forbidden atoms, with the
  H=21 proof pattern reducing to strong components plus lower bounds on odd
  cycles from pancyclicity.
- **LRC use:** Converts a residual packet target into a concrete endpoint:
  construct a tournament-conflict or OCF packet whose H-value must be `7` or
  `21`, then invoke impossibility.
- **Preserves:** Strong-component factor, OCF odd-cycle count, connected
  forbidden atom, and state-lift equality.
- **Forgets / guardrail:** Ordinary digraphs and partial orientations can
  realize forbidden-looking values.  The lift must land in the tournament or
  connected OCF category, not merely a binary relation.
- **Next pull:** Make a verifier that checks whether each proposed F7 state
  lift really constructs a complete tournament/OCF packet rather than a loose
  digraph.
- **Pointers:** THM-200, THM-572, HYP-2908,
  `07-reflections/h21-proof-complete-s680.md`, LTT-029, LTT-040.

### LTT-053: Metric Comparator / Trienerment Dichotomy

- **Move:** Read a tournament as a comparator on a metric with a threshold and
  possible tie layer; separate geometric comparators that collapse to the
  circular runner picture from arithmetic comparators that retain residue,
  exact-period, or p-adic channels.
- **LRC use:** A sanity check for new analogies: ask whether the construction
  is only a geometric shadow of existing circular order, or whether it carries
  a new arithmetic coordinate needed by LRC14.
- **Preserves:** Metric, threshold/tie convention, arithmetic channel, and
  observer/source mark.
- **Forgets / guardrail:** Monotone geometric metrics usually add no new LRC
  data.  Arithmetic metrics need explicit residue/carry labels before they are
  proof carriers.
- **Next pull:** Add a metric-comparator audit line to future technique cards:
  geometric collapse, arithmetic channel, or mixed.
- **Pointers:**
  `07-reflections/lrc-a-tournament-is-a-comparator-on-a-metric-the-geometric-arithmetic-dichotomy-s541o.md`,
  LTT-004, LTT-043, LTT-057.

### LTT-054: Nonabelian Character-Ratio / Alternating-Group Carrier

- **Move:** Treat highly symmetric nonabelian carriers through character
  ratios, not only abelian Fourier modes; alternating-group and icosahedral
  lanes are candidates for residual parity and path-action structure.
- **LRC use:** Provides a possible nonabelian Fourier transform for packet
  sectors where cyclic Ramanujan projectors are too abelian to see the
  obstruction.
- **Preserves:** Group action, conjugacy class, representation character
  ratio, parity sector, and orbit stabilizer.
- **Forgets / guardrail:** Matching a character ratio or Platonic count is
  not a proof unless the LRC packet predicate is transported through the group
  action.
- **Next pull:** Test whether F7/Johnson-harmonic residuals admit a
  nonabelian character-ratio formulation over a small action group.
- **Pointers:** alternating-group graph reflections, icosahedral tangents,
  LTT-030, LTT-049.

### LTT-055: Converse-Z2 / Half-Arc-Transitive Orientation

- **Move:** Use self-converse, half-arc-transitive, and orientation-reversal
  phenomena as tests for whether a quotient has accidentally identified a
  direction that the proof needs.
- **LRC use:** AP/GW boundary rows and symmetric packet fibers often have a
  hidden reversal or antipodal ambiguity; this card asks whether the ambiguity
  is harmless, reconstructible, or a live residual.
- **Preserves:** Orientation bit, complement/reversal action, half-arc
  orbit, and fixed path under the symmetry.
- **Forgets / guardrail:** A self-converse shadow is not a complete
  tournament proof unless endpoint owners and tie directions are retained.
- **Next pull:** Mark every AP/GW and C27 transfer by its reversal/complement
  action and record where the action changes packet route.
- **Pointers:** half-arc-transitivity/self-converse notes, LTT-036, LTT-044.

### LTT-056: Round-Tournament Realizability Filter

- **Move:** Before trusting a tournament shadow of a runner movie, filter it
  through circular/round-tournament realizability constraints and record which
  SCC patterns are impossible for clock movies.
- **LRC use:** Blocks false tournament counterexamples whose abstract
  orientation cannot come from a circular threshold relation.
- **Preserves:** Circular order, interval neighborhood, SCC end structure, and
  tie Hamiltonian path.
- **Forgets / guardrail:** Abstract tournament enumeration produces many
  shadows that are not realizable by LRC clock images.  Do not feed these into
  a state-lift as if they were packet rows.
- **Next pull:** Add a round-realizability flag to tournament shadows used by
  the HYP-2963 classifier and the HYP-2995 carrier tournament.
- **Pointers:** round tournament threads, THM-354, HYP-2924, LTT-003,
  LTT-004, LTT-045.

### LTT-057: Paley / Frobenius Arithmetic Tournament Carrier

- **Move:** Use Paley and Frobenius tournaments as arithmetic difference
  carriers: quadratic residues orient edges, automorphism groups normalize
  path counts, and Gauss/Ramanujan sums expose exact-period structure.
- **LRC use:** Gives a prime-modular model for exact-period channels near the
  `14=2*7` wall, especially when a row behaves like a residue-difference
  packet rather than a circular metric packet.
- **Preserves:** Prime modulus, residue class, Frobenius action, automorphism
  size, and Gauss/Ramanujan side channel.
- **Forgets / guardrail:** Paley smoothness or high H is often coincidence.
  Use Paley as an arithmetic carrier, not as a scalar maximizer proof.
- **Next pull:** Compare Ramanujan exact-period projectors on q=7, 14, 27, and
  41 packets against Paley/Frobenius residue tournaments.
- **Pointers:** Paley path-homology notes, Paley ratio tangents, HYP-2979,
  LTT-030, LTT-043, LTT-051.

### LTT-058: Pascal-Slope / Additive-Basis Farey Packet Schema

- **Move:** Treat Pascal-slope row vectors, additive-basis representation
  fibers, Zeckendorf carry normal forms, and Farey operator lanes as vertices
  in a proof-carrier tournament.
- **LRC use:** Adds a labelled way to ask whether a hard row is being explained
  by additive abundance, ternary smoothing, bounded polygonal residue cover,
  Fibonacci/Zeckendorf carry normal form, product-incidence debt, or only a
  power/magnitude shadow.
- **Preserves:** Exact `p/q`, Farey excess, `p+q` additive lane, `p*q`/`Kpq`
  factor fiber, Pascal row vector, representation entropy, local residue rank,
  and carry width.
- **Forgets / guardrail:** The scalar Fibonacci term, raw Goldbach count,
  polygonal cover count, product value, or power value is unsafe unless the row
  fiber, carry rule, exact Farey root, and endpoint packet labels are retained
  or discharged.
- **Next pull:** Add `additive_basis_regime`, `representation_entropy`,
  `local_residue_rank`, `carry_width`, `pascal_slope_row_id`, and
  `farey_operator_lane` to the HYP-2963 packet classifier.
- **Pointers:** HYP-2999, HYP-2998, HYP-2940, HYP-2934, HYP-2932,
  HYP-2931, HYP-2523, S501 additive-basis reflections, LTI-149.

### LTT-059: Curried Packet Functional Tower

- **Move:** Treat every LRC proof route as a curried evaluator
  `S -> packet -> root -> lane -> fiber -> certificate -> verdict`.
  Tournament vertices are the partial-evaluation carriers, not runners.
- **LRC use:** Makes every quotient auditable: once a coordinate is fixed,
  summed, or forgotten, emit the lost-coordinate function on the remaining
  fiber and discharge it before scalarizing.
- **Preserves:** Argument order, exact packet root, Farey lane, residual
  section, certificate family, and named verdict.
- **Forgets / guardrail:** A raw runner movie, scalar safe mass, H-value,
  product, or representation count is illegal unless all earlier partial
  evaluations have closed by zero/reconstruction/coboundary/dual/family/
  AP-GW-boundary/residual debt.
- **Next pull:** Add `curried_call_signature` and
  `lost_coordinate_function` to HYP-2963 records and Fejer/Ramanujan manifests.
- **Pointers:** HYP-3002, HYP-3000, HYP-2999, HYP-2997, HYP-2996,
  HYP-2995, HYP-2974, HYP-2963, LTI-152.

### LTT-060: Poincare Worldline Frame Ledger

- **Move:** Treat an LRC row as worldlines `x_i(t)=v_i*t` in a time/phase
  cylinder with a danger tube around the observer worldline.  Use the
  Poincare/Galilean analogy only after declaring whether the observer velocity,
  tube metric, integer lattice, and sign orientation are retained.
- **LRC use:** True anchored-LRC automorphisms include runner permutations,
  independent sign flips, reflection/time reversal, and integer dilation with
  primitive scaling.  Common boosts are exact only in observer-coupled packets;
  stationary speed translation is a warning signal for a forgotten observer
  label.
- **Preserves:** Observer-relative predicate when frame labels are retained,
  sign-kernel status, relative speed normal form, primitive scale, tube metric,
  and boost/recentering cocycle.
- **Forgets / guardrail:** Bare winding tournament order forgets sign-kernel
  debt, metric gaps, observer placement, and exact Farey scale.  Lorentz-like
  velocity addition forgets the integer-speed lattice and fixed circle metric
  unless they are carried as deformed-tube cochains.
- **Next pull:** Add worldline-frame fields to HYP-2963 packet records and test
  sign-kernel/boost-admissibility lemmas before any Poincare-flavored proof
  route scalarizes speeds.
- **Pointers:** HYP-3007, HYP-3006, HYP-3002, HYP-2997, HYP-2963,
  HYP-2486, HYP-2291, THM-381, THM-385, LTI-157.

### LTT-061: Automatic-Gap / Power-Lift Packet Ledger

- **Move:** Treat Moser-de Bruijn words, fibbinary/Zeckendorf carry words,
  Hadamard-lacunary support ratios, Fermat-Catalan power guards, Hurwitz
  doubling-CF states, and visibility approximations as proof-carrier vertices.
- **LRC use:** Before a sequence, gap, or perfect-power shadow is used in an
  LRC14 proof, test whether it is route-pure on labelled packets or whether it
  mixes AP/GW, C27, K33, Res_27, and covering residual families.
- **Preserves:** Automatic-language state, carry status, lacunary address,
  power-lift obstruction label, dyadic doubling clock, visibility guard, and
  residual route.
- **Forgets / guardrail:** Raw scalar words collapse important distinctions:
  AP13 and GW `12->24` share the word `MFCMMCCFFFCCC`, so automatic shadows are
  not certificates unless packet labels remain attached.
- **Next pull:** Add `automatic_language_class`,
  `fibbinary_carry_status`, `moser_even_bit_status`,
  `ostrowski_digit_system`, `lacunary_gap_ratio`, `power_lift_guard`,
  `fermat_catalan_residual`, `hurwitz_doubling_cf_state`, and
  `visibility_potato_approx_guard` to HYP-2963 records, then compute
  route-purity and first mixed fibers.
- **Pointers:** HYP-3009, HYP-3008, HYP-3007, HYP-3003, HYP-3000, HYP-2998,
  HYP-2963, HYP-2950, HYP-2944, HYP-2937, HYP-2702, HYP-2698, HYP-1920,
  HYP-1902, LTI-159.

### LTT-062: Automatic / Lacunary Safe-Component Filter

- **Move:** Treat finite automata and lacunary sequence languages as packet
  filters: a row may carry DFA state, 2-adic window state, gap-block profile,
  and exception-ledger labels, but those labels are not final scalar
  invariants.
- **LRC use:** Imports the user's Fermat-Catalan, 2-adic/Hurwitz,
  Ostrowski-Hadamard, fibbinary, and Moser-de Bruijn prompts as proof roles.
  This is the exact safe-component companion to HYP-3008's automatic
  gap-language membership audit.  S171 keeps AP and GW `12->24` as zero-open
  denominator-14 boundary atoms; the first-13 fibbinary and Moser-de Bruijn
  rows instead have positive safe mass `66077/399840` and
  `4264747/40348854`.
- **Preserves:** LRC predicate, exact `M=p/q`, endpoint owners,
  exact-period labels, finite-state word state, gap-block/lacunarity profile,
  first safe component, and named certificate/residual route.
- **Forgets / guardrail:** The raw sequence name, growth rate, gap ratio, or
  automatic-language inclusion forgets real circle gaps.  It is proof-safe only
  when paired with packet labels or discharged by a strict-safe component,
  AP/GW boundary equality, Fejer/Ramanujan/endpoint annihilation, family
  descent, or named F7/THM-572 residual debt.
- **Next pull:** Pair HYP-3008's `automatic_gap_carrier` with exact
  safe-component fields (`largest_component`, `safe_measure`,
  `boundary_units`, and `automatic_filter_exit`) in HYP-2963 packet manifests,
  then audit whether any zero-open non-AP/GW packet survives after these
  labels and existing dual exits are attached.
- **Pointers:** HYP-3011, HYP-3008, HYP-3002, HYP-3000, HYP-2997, HYP-2996,
  HYP-2963, HYP-1902, THM-572, LTI-160.

### LTT-063: Gap Automaton Carrier Tournament

- **Move:** Treat finite automata, lacunary support, valuation budgets,
  visibility cores, and induced tournament class ledgers as proof-carrier
  vertices.  Orient edges by retained LRC predicate dimensions or by declared
  rank-priority gauges.
- **LRC use:** A sequence-shadow quotient is admissible only after retaining its
  automaton state, native transition (`x->2x` versus `x->4x`), gap-boundary
  label, finite-exception budget, packet route, and induced tournament class
  census.
- **Preserves:** Automaton language, state word, support ratio, doubling
  transition, base-4 digit mask, Zeckendorf carry state, valuation budget,
  visibility-core label, and induced isomorphism-class word.
- **Forgets / guardrail:** Raw fibbinary/Moser membership mixes all `14` residue
  classes and cannot distinguish AP/GW equality atoms or hard residual routes.
- **Next pull:** Build product automata on the unit-excess lane `q=14p-1` and
  compare hard non-AP/GW packet tournaments against the S173 `n=4,5,6`
  canonical class words.
- **Pointers:** HYP-3012, HYP-3011, HYP-3010, HYP-3009, HYP-3008, HYP-3007, HYP-3006, HYP-2998,
  HYP-2997, HYP-2983, HYP-2982, HYP-2963, THM-572, LTI-161.

### LTT-064: Perfect-Number Divisor Packet Tournament

- **Move:** Treat perfect-number controls, divisor-lattice packets, Farey
  unit-excess address, Kpq product incidence, and automatic power states as
  proof-carrier vertices.  Orient edges by retained LRC predicate dimensions,
  not by raw product size.
- **LRC use:** The Euclid-Euler `n=2` chain is an exact calibration lane, while
  LRC14 `q=14a-1` rows are deficient only under a prime-q side condition.
  Composite `q14` rows can flip abundant, so factorization and abundancy defect
  are load-bearing packet fields.
- **Preserves:** Exact `M`, unit-excess apex, prime/composite `q` flag,
  divisor factorization, abundancy defect, product/Kpq route, and automaton
  transition state.
- **Forgets / guardrail:** Raw product scalar, power-of-two address, or
  fibbinary/Moser membership erases the difference between perfect controls,
  deficient prime-q LRC14 shadows, and abundant composite-q rows.
- **Next pull:** Add `unit_excess_apex`, `perfect_control_status`,
  `abundancy_defect`, `divisor_lattice_factorization`, `prime_q_flag`,
  `product_incidence_rank`, and `automaton_transition_state` to HYP-2963
  sidecars, then rerun unit-excess route purity on `q=14p-1`.
- **Pointers:** HYP-3013, HYP-3012, HYP-3009, HYP-3008, HYP-2946, HYP-2945,
  HYP-2941, HYP-2221, HYP-2220, HYP-2963, THM-572, LTI-162.

### LTT-065: Creative Packet-Lens Tournament

- **Move:** Treat speculative proof ideas as exact packet lenses: Cech nerve
  classes, tropical slack potentials, CRT/solenoid charts, endpoint
  chip-firing currents, danger-count distributions, matroid tope/cocircuit
  walls, and automaton/divisor sidecars.
- **LRC use:** Separates AP/GW zero-open boundary controls from K33/petal
  positive-open rows that still have denominator-14 boundary witnesses, and
  from covering rows whose first useful chart can move to denominator `41`.
- **Preserves:** Safe-component topology, exact slack, CRT witness chart,
  endpoint current, danger-count distribution, tope/cocircuit state, and typed
  arithmetic/automatic sidecars.
- **Forgets / guardrail:** Any one creative lens forgets another coordinate:
  Cech loses arithmetic period, tropical slack loses owner identity, CRT loses
  real interval topology, and danger-count duals lose route labels.
- **Next pull:** Add `cech_nerve_class`, `positive_component_count`,
  `tropical_slack_margin`, `crt_solenoid_first_chart`,
  `endpoint_current_word`, `danger_count_distribution`,
  `tope_cocircuit_wall_state`, and `automaton_divisor_sidecar` to HYP-2963
  sidecars, then test hard covering families.
- **Pointers:** HYP-3014, HYP-3013, HYP-3012, HYP-3008, HYP-2974, HYP-2973,
  HYP-2970, HYP-2969, HYP-2965, HYP-2963, HYP-2949, HYP-2948, THM-572,
  LTI-163.

### LTT-066: Automaton Fiber-Mixing Quotient Tournament

- **Move:** Treat quotient candidates as tournament vertices: exact labelled
  packets, Farey magnitude height, residue+automaton words, terminal DFA
  words, perfect-power words, gap-ratio buckets, and raw counts.  Orient by
  fiber purity, finite-state checkability, magnitude retention,
  residue/endpoint retention, route compatibility, anti-scalar guard, and
  proof cost.
- **LRC use:** A finite-state automaton quotient is admissible only after exact
  fiber-mixing stress.  S175/HYP-3016 finds that residue+Moser/fibbinary
  terminal fields mix AP/Goddyn-Wong boundary atoms with strictly open
  single-swap rows.  S187/HYP-3023 extends the stress to the full HYP-2963 bank
  and finds exact magnitude cocycle is the first tested non-route coordinate
  with zero mixed theorem-route fibers.
- **Preserves:** The tested quotient surface and exact safe-component status,
  including mixed boundary/open fibers and the need for magnitude/Farey
  side-channel data.
- **Forgets / guardrail:** Terminal automaton state and residue data forget the
  magnitude coordinate distinguishing AP from `12->26`/`12->96` and GW from
  later one-dipole tails.
- **Next pull:** Prove the family magnitude-cocycle lemma inside
  automatic/residue fibers, starting with `MFCMMCCFFFCCC`; use persistence
  barcode, Fejer, Ramanujan, Haar, and packet zippers as certificate anchors
  when magnitude alone needs a handoff.
- **Pointers:** HYP-3023, HYP-3022, HYP-3021, HYP-3020, HYP-3019, HYP-3018, HYP-3017, HYP-3016, HYP-3015,
  HYP-3014, HYP-3013, HYP-3012, HYP-3011, HYP-3009, HYP-3008, HYP-3002,
  HYP-2997, HYP-2963, HYP-2928, THM-572, LTI-170, LTI-169, LTI-168, LTI-167,
  LTI-166, LTI-165.

### LTT-067: Discrepancy-Height Trident Tournament

- **Move:** Treat proof carriers as tournament vertices: exact labelled packet,
  safe topology / barcode, discrepancy-height-Hensel trident, residue-height
  pair, Hensel singular-lift guard, Erdos-Turan residue discrepancy,
  Mahler-height proxy, automatic word sidecar, and raw scalar family name.
  Orient edges by retained predicate purity, route purity, magnitude
  retention, discrepancy retention, local lift stability, certificate handoff,
  finite cost, and anti-scalar guard.
- **LRC use:** The trident is the first bounded carrier after HYP-3016/HYP-3017
  that splits AP/Goddyn-Wong boundary rows from open rows inside mixed
  residue/automatic fibers without using exact speed tuple identity.  On the
  S184 named plus single-swap bank it has `0` mixed boundary/open fibers,
  while automaton words, residue+MFC pairs, residue discrepancy alone, Hensel
  alone, and height alone all still leak.
- **Preserves:** Boundary-vs-open predicate in the bounded scout, route-purity
  telemetry, small-denominator discrepancy, lost magnitude/Farey scale,
  Hensel singular-root status, and a Beck-Fiala-style bounded feature
  incidence surface.
- **Forgets / guardrail:** The full trident is nearly exact (`2167` fibers for
  `2173` rows), so it is not yet a compressed theorem coordinate.  It also
  forgets endpoint-owner geometry, exact Fejer atom banks, and barcode
  topology unless HYP-3015/HYP-2981 fields are reattached.
- **Next pull:** Add trident fields to the full HYP-2963 packet bank; then
  coarsen residue denominators, Erdos-Turan bins, height buckets, and Hensel
  flags until the smallest route-pure signature is found for the large mixed
  automatic fibers, especially `MFCMMCCFFFCCC`.
- **Pointers:** HYP-3020, HYP-3017, HYP-3016, HYP-3015, HYP-3014, HYP-3009,
  HYP-3008, HYP-2997, HYP-2995, HYP-2991, HYP-2989, HYP-2963, THM-572,
  LTI-167.

### LTT-068: Pair-Good Decoy Generator / Barcode Classifier Tournament

- **Move:** Treat source pair lanes, blocker residue teeth, active blocker
  decks, generator-cover sizes, barcode relations, normal-fan support
  relations, and raw pair-good counts as tournament vertices.  Orient by
  retained blocker mechanism, modular tooth inequality, source-lane
  specificity, support/barcode relation, and whether the carrier separates
  AP/GW boundary atoms from positive-open rows.
- **LRC use:** HYP-3021 replaces raw pair-good decoy counts by the exact
  generator rule `14*min(c*p mod q,q-c*p mod q)<q`; HYP-3022 then records how
  each generated false switch sits relative to HYP-3015 barcode bars and
  HYP-3018 normal-fan peak supports.  Large tails such as `12->200`,
  drop6-add180, and covering `12->84` mostly recycle AP-core blocker teeth
  rather than producing new obstruction families.
- **Preserves:** Source pair lane, source shell, blocker role, blocker-depth
  bucket, zero-tooth flag, generator-cover size, top generator key, barcode
  relation, and normal-fan support relation.
- **Forgets / guardrail:** A raw count of pair-good times forgets that many
  times can land on the same modular tooth cover, and can also forget that the
  same false switch is outside all strict bars or already controlled by a
  normal-fan support.  It should not be used as an obstruction-size scalar
  without the generator and barcode/normal-fan sidecars.
- **Next pull:** Add the decoy-generator and barcode/normal-fan fields beside
  active-pair and arc-Cech sidecars in HYP-2963 packet records, then test
  whether residual packet families have bounded generator-cover and
  blocker-deck templates.
- **Pointers:** HYP-3022, HYP-3021, HYP-3019, HYP-3018, HYP-3015, THM-524,
  HYP-2990, HYP-2963, LTI-169, LTI-168.

### LTT-069: Closed Arc-Cech Nerve Carrier Tournament

- **Move:** Treat proof carriers around the exact danger cover as tournament
  vertices: endpoint tope/cocircuit wall, individual-arc closed Cech nerve,
  boundary owner current, safe interval measure, runner quotient nerve, Fejer
  dual certificate, automaton sidecar, and raw speed/sequence scalar.  Orient
  by preserved circle-cover predicate, endpoint equality retention, quotient
  defect visibility, and certificate handoff strength.
- **LRC use:** Makes the exact closed cover topology primary.  AP/Goddyn-Wong
  are full-cover cycles in the endpoint-completed individual-arc Cech nerve
  with open pieces glued by six boundary cocircuits; K33, petal, covering,
  fibbinary, and Moser controls have closed arc `beta1=0` and positive safe
  mass.
- **Preserves:** Closed arc Betti numbers, open arc component count, boundary
  cocircuit facet word, owner sums mod `14`, runner quotient Betti defect,
  private arc/runner counts, safe tope count, and exit route.
- **Forgets / guardrail:** Runner-level nerves and sequence shadows can merge
  several disjoint danger arcs owned by one speed, erasing either disconnected
  cover pieces or the full-cover cycle unless `runner_quotient_betti_defect`
  remains attached.
- **Next pull:** Add the arc-Cech fields to HYP-2963 packet records or a
  sidecar, run the full packet bank, prove K33/petal `beta1=0` exits
  familywise, and define F7 as a good-cover quotient-defect class before using
  more scalar filters.
- **Pointers:** HYP-3025, HYP-3024, HYP-3023, HYP-3022, HYP-3021, HYP-3020, HYP-3018,
  HYP-3016, HYP-3015, HYP-3014, HYP-3013, HYP-3012, HYP-3011, HYP-3010,
  HYP-3009, HYP-3008, HYP-2997, HYP-2986, HYP-2975, HYP-2974, HYP-2970,
  HYP-2963, HYP-2990, THM-572, LTI-172.
### LTT-070: Fiber-Zipper Convergence Tournament

- **Move:** Treat zipper gates as tournament vertices: automatic word,
  residue-terminal fiber, exact ET clocks, Henselian unit-root rule, coarse
  ET+unit gate, magnitude cocycle, magnitude+ET+unit, and barcode/packet
  zipper.  Orient by route purity, max mixed-fiber size, ET discrepancy
  retention, unit-lift stability, magnitude retention, topology, packet
  labels, finite-state checkability, and proof cost.
- **LRC use:** HYP-3024 shows the full HYP-2963 bank admits a cleaner
  two-stage proof target.  Exact ET clocks at `14,27,41` split to singleton
  fibers and are too close to an address coordinate.  The coarse ET+unit gate
  has `21702` fibers, `15` mixed theorem-route fibers, max mixed `4`, and
  `0` mixed boundary/open fibers.  Thus status convergence can be attacked
  before full route purity.
- **Preserves:** Boundary/open status, coarse residue discrepancy, p-adic
  unit-root clocks, zero-root scale debt, and the residual route-purity
  telemetry needed to schedule Fejer/Ramanujan/Haar, magnitude, K33/F7, or
  covering certificates.
- **Forgets / guardrail:** Exact ET clocks can become a packet address, and
  unit-Hensel counts alone still mix AP/GW boundary with open rows.  The
  quotient is admissible only as a status-preserving gate plus a named
  certificate scheduler.
- **Next pull:** Prove the coarse ET+unit gate cannot mix AP/GW boundary
  equality with strict-open packets inside automatic/residue fibers; then
  route the remaining open-route collisions through magnitude formulas and
  the existing certificate stack.
- **Pointers:** HYP-3024, HYP-3023, HYP-3020, HYP-3017, HYP-3016, HYP-3015,
  HYP-2963, THM-572, LTI-171, LTI-170, LTI-167.

### LTT-071: Carrier-Fusion Switchboard Tournament

- **Move:** Treat fused proof sidecars as tournament vertices: labelled packet
  fusion signatures, safe-stick/potato bodies, lonely-profile barcodes,
  magnitude-cocycle fibers, endpoint currents, CRT charts, danger-count duals,
  Hurwitz doubling states, automatic words, and raw row names.  Orient by
  retained predicate, exactness, topology, arithmetic, endpoint data, route
  split, family transfer, computability, and anti-scalar discipline.
- **LRC use:** Converts the HYP-3014 creative lenses, HYP-3015 barcode,
  HYP-3016 magnitude-cocycle guardrail, HYP-3017 sidecar route-purity
  failure, HYP-3018 active normal-fan sidecar, HYP-3020
  discrepancy-height trident, HYP-3021 pair-good decoy generator classifier,
  HYP-3022 pair-good barcode/normal-fan refinement, HYP-3023 automatic
  fiber zipper, HYP-3024 fiber-zipper convergence audit, and HYP-3025 closed
  arc-Cech carrier into one exact packet switchboard.  In the S189
  named bank,
  automatic words and chart denominators still mix AP/Goddyn-Wong boundary
  atoms with open rows, while barcode shape, magnitude cocycle, and the full
  fusion signature remove boundary/open leakage.
- **Preserves:** Strict-open status, exact safe-component topology, largest
  safe-stick and safe-body mass, CRT first chart, endpoint current, magnitude
  cocycle, automatic/doubling sidecars, ET/Henselian convergence sidecar,
  danger-count duals, and route labels.
- **Forgets / guardrail:** Any proper subcarrier can forget a load-bearing
  coordinate: automatic state forgets magnitude, chart denominators forget safe
  topology, endpoint current forgets open bars, and safe-stick geometry forgets
  arithmetic period unless packetized.
- **Next pull:** Add `fusion_signature`, `largest_safe_stick`,
  `safe_body_mass`, `barcode_shape`, `magnitude_cocycle`,
  `endpoint_current_word`, `crt_first_chart`, `danger_distribution_word`,
  `et_henselian_unit_zipper`, and `doubling_transition_word` to HYP-2963
  packet sidecars, then rerun the full bank to test whether
  `automatic_plus_barcode_shape` remains boundary/open pure.
- **Pointers:** HYP-3026, HYP-3025, HYP-3024, HYP-3023, HYP-3022, HYP-3021,
  HYP-3020, HYP-3018, HYP-3017, HYP-3016, HYP-3015, HYP-3014, HYP-3013,
  HYP-3009, HYP-2963, HYP-2974, HYP-2969, THM-572, LTI-173, LTI-172,
  LTI-171, LTI-170, LTI-169, LTI-168, LTT-070, T1106, T1105, T1104.

### LTT-072: Side-Channel Repair Ladder Tournament

- **Move:** When a quotient fails fiber-purity, make the possible repairs the
  tournament vertices: automatic word, q-threshold, q-factor/power guard, tail
  magnitude, exact `M`, exact `M+q`, boundary topology, HYP-3020 trident
  coordinates, HYP-3021 decoy-generator coordinates, HYP-3022 barcode/normal-fan
  coordinates, HYP-3023 automatic zipper coordinates, HYP-3024 convergence
  coordinates, HYP-3025 arc-Cech topology, packet labels, guarded non-route
  signature, and circular route/exit diagnostics.
- **LRC use:** HYP-3027 turns the HYP-3017/HYP-3023 automatic-word failure and
  the HYP-3020/HYP-3021/HYP-3022/HYP-3024/HYP-3025 side-channel carriers into
  an ordered proof target.  Exact `M` repairs open/boundary status but leaves
  `366` mixed theorem-route fibers;
  `M+q` and boundary topology each leave one mixed route pair; packet labels and
  the guarded non-route signature are route-pure on the audited HYP-2963 bank.
- **Preserves:** Fiber-purity tests for open/boundary status and theorem route,
  the declared side-channel retained by each repair, and the noncircularity
  status of the repair coordinate.
- **Forgets / guardrail:** Route labels and automatic-filter exits are circular
  diagnostics, not admissible proof quotients.  Exact scale alone forgets
  packet route; tail magnitude alone remains route-mixed.
- **Next pull:** Prove local zipper repairs for the two residual mixed pairs:
  `two drop(10,13)->add(20,26)` versus `two drop(8,12)->add(16,24)`, and
  `two drop(12,13)->add(26,36)` versus `single swap 12->72`.
- **Pointers:** HYP-3027, HYP-3026, HYP-3025, HYP-3024, HYP-3023, HYP-3022, HYP-3021, HYP-3020, HYP-3018, HYP-3017, HYP-3016,
  HYP-3015, HYP-3014, HYP-2963, HYP-2997, HYP-2995, HYP-2992, THM-572,
  LTI-174, LTI-173, LTI-172, LTI-171, LTI-170, LTI-169, LTI-168, LTI-167, LTI-166.

### LTT-073: Carrier Pullback Mega-Index Tournament

- **Move:** Treat project techniques themselves as proof-carrier pullback
  vertices.  Each vertex is a row in
  `00-navigation/LRC-CARRIER-PULLBACK-INDEX.md`: a tournament/metagraph,
  series, topology, harmonic, arithmetic, automaton, geometric, formal, or
  computation technique with named retained LRC fields and a quotient
  guardrail.
- **LRC use:** Gives future agents a large operational menu before they reach
  for a new scalar.  The index currently has `90` `CPI-*` rows that pull back
  A000568/source lifts, OCF/Walsh/deletion-contraction, path homology,
  metagraph Laplacians, exact interval topology, Cech/topes/cocircuits,
  Fejer/Toeplitz/Ramanujan/Haar, Mobius/totient/Farey/Fibonacci/Zeckendorf,
  Moser/fibbinary/Fermat-Catalan, perfect-number divisor packets,
  Poincare/worldline frames, K33/C27/state-lift geometry, and formal/checkpoint
  tooling into LRC packet fields.
- **Preserves:** The proof-interface discipline: boundary/open status,
  theorem route, exact scale, endpoint/topology, arithmetic period, harmonic
  certificate, state-lift residual, family-transfer promise, formal/checkable
  payload, and proof cost.
- **Forgets / guardrail:** A big index can become a passive bibliography.  A
  row is only useful after another agent instantiates its retained fields,
  runs a fiber-mixing audit, and records a carrier tournament fingerprint or a
  negative quotient guardrail.
- **Next pull:** Pick one row from each bundle in the pullback index, implement
  its packet sidecar on HYP-2963 or a named stress family, and promote any
  successful row into a theorem-facing `LTI-*` or `LTT-*` card.
- **Pointers:** T1108, LTI-175,
  `00-navigation/LRC-CARRIER-PULLBACK-INDEX.md`, LTI-173, LTI-156,
  LTI-147, LTI-021, LTM-071, LTM-072.

## Immediate Pull List

1. Expand the HYP-2963 labelled packet classifier with Haar tile class,
   Ramanujan exact-period projector, spectrum binding scale, additive-basis
   regime, Pascal-slope carry width, curried call signature,
   lost-coordinate function, observer velocity, sign-kernel status, boost
   cocycle, tube metric, automatic language class, fibbinary/Moser carry
   status, lacunary gap ratio, gap-block profile, safe-component readout,
   power-lift guard, Hurwitz doubling-CF state, automaton state word, native
   doubling transition, finite-exception budget, induced tournament class word,
   unit-excess apex, perfect-control status, abundancy defect,
   divisor-lattice factorization, prime-q flag, product-incidence rank,
   automaton transition state,
   Cech nerve class, positive-component count, tropical slack margin, CRT
   solenoid first chart, endpoint current word, danger-count distribution,
   tope/cocircuit wall state, automaton-divisor sidecar,
   lonely-profile barcode, bar count, persistence margin,
   residue-automaton fiber ID, magnitude cocycle, Farey magnitude height,
   fiber anchor row, fiber-mixing exit, automatic fiber zipper,
   fiber-convergence delta, sidecar redundancy flag, closed arc-Cech beta,
   open arc component count, boundary cocircuit facet word,
   runner quotient Betti defect,
   route-purity split report, discrepancy-height trident signature,
   residue discrepancy bins, Erdos-Turan proxy bins,
   Henselian unit-root rule, zero-root scale debt,
   coarse ET+unit status-convergence gate,
   Mahler/Farey height bucket, Hensel singular counts,
   fusion signature, largest safe stick, safe-body mass, barcode shape,
   danger-distribution word, doubling-transition word, carrier-fusion exit,
   repair-ladder stage, first nonzero repair cochain, route-purity status,
   status-purity status, guarded non-route signature,
   carrier-pullback row ID, destroyed-coordinate ledger,
   Burnside cost, score-class H-spread, and round-realizability flag.
2. Make a Fejer certificate manifest bridge checklist based on LTT-044, then
   add interval-arithmetic proof anchors for the floating Fejer evaluations.
3. Compute multi-scale tournament spectra for AP, GW, K33, petals, splices,
   covering rows, and the weakest Fejer-margin rows.
4. Build a decorated source-cone canonicalizer: observer-source mark,
   endpoint owners, exact-period class, Haar tile class, and p-adic carry
   owner.
5. Prove or refute the HYP-2992 vanishing lemma: no owner-strip,
   no cross-handoff, no nested-refinement coefficient implies AP/GW boundary or
   F7 state-lift.
6. Build a deletion-contraction packet metagraph and record DC depth,
   near-tournament branches, and Delta-H flip energy for named rows.
7. Run path-homology and OCF-activity fingerprints on packet-carrier
   tournaments, especially K33, P10+GW, C27 petals, and hard Fejer rows.
8. Add an H-gap state-lift verifier that rejects loose digraphs and accepts
   only complete tournament/connected OCF packets.
9. Audit each new analogy with the metric-comparator dichotomy: geometric
   collapse, arithmetic channel, or mixed carrier.
10. Turn the controlled-kernel rule into a short reusable theorem template in
    the style of HYP-2990.
