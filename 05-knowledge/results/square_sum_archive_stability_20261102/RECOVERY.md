# Square-sum recovery: mechanisms, corrections and less-used neighbours

## Scope of the search

Baseline: `66cf347030dc`. `CATALOG.json` indexes **936 tracked Markdown files,
665 distinct Git blobs**, found by square-sum/squareful/Gerbicz/nice-pair and
sequence-name variants. There are 299 direct-subject hits and 241 after the
path-based copy/audit filter. These are files, **not independent sessions**:
some messages are replicated, several packets are successive checkpoints,
and some 'square sums' concern a different operation entirely.

A separate all-reachable-ref filename-history scan examined 2,641 paths.
The six absent from the baseline concern square-prefix algebra, Danielewski
surfaces and elliptic splitting, rather than missing Hamiltonian square-sum
sessions. This supplements current-tree discovery; it does not exhaust
historical text that never used any searched name, untracked chats or files
outside the repository. The index preserves all hits so later searches need
not repeat this discovery pass.

Below is the mechanism-level review. Statements marked inherited PROVED are
read from their theorem/proof packets or explicit repaired statements; this
session does not claim to rerun every old certificate. Hand-offs were used to
find stronger successors, not to overrule corrections or proof documents.

## The actual strongest construction frontier

Ordinary interval paths are not the remaining conjecture: they exist exactly
for n in {1,15,16,17,23} union [25,infinity). Ordinary cycles exist exactly
for n>=32. The repository has reconstructed Gerbicz's complete construction.
The useful stronger results are **unchanged prescribed labels, synchronized
phase, whole boundary relations, exact closed subgraphs, and multiplicity**.

| Recovered route | Mechanism and scope | Tangential concepts / retained sidecar |
|---|---|---|
| [September22 Hamiltonicity](../collatz_mod6_20260922_w6_square_sum_hamiltonicity.md) and [square filter](../collatz_mod6_20260922_w6_summand_closure_square_filter.md) | Exact early census; connectivity from14; reflection x->s²-x; the 25-fold lift. Corrected September25 endpoint forcing supersedes the original certificate claim. | Distinct-summand closure, Collatz admission filters, self-loop holes, born-leaf windows. Three components are not the closure holes {1,4,6}. |
| [September25 decoder](../decoder_prime_square_20260925.md) | Full endpoint-aware obstruction at24; every path at25 uses 11--25--24. Halving toggles square classes1 and2, quadrupling restores. | Odd-square brackets give {2,3,11}; valuation-layer seams need an arithmetic sheet AND attachments. |
| [September26 square exchanges](../crossroads_family_20260926_squares.md) | Unicyclic criterion explains15,16,17; the first four-cycle and fixed-endpoint two-edge cycle switch occur at46. | Signless incidence is the exact lift operator; even-cycle equations and odd-cycle pins, not a generic 'fractal recurrence'. |
| [Optimal seeded threshold](../seeded_optimal_threshold_6242_20261019/OPTIMAL_THRESHOLD_THEOREM.md) | Unchanged 482--3362--6242 extends iff N>=6242. Exact trace covers6242..128452, followed by the universal induction. | Source preservation, sequence reversal certificates, no disconnected 2-factor accepted as a cycle. |
| [Simultaneous matching theorem](../square_sum_simultaneous_matching_20261019/SIMULTANEOUS_MATCHING_THEOREM.md) | At T_k=(107·49^k-1)/2, all 48·49^(k-1) specified edges of sum4·49^k coexist for k>=3, along with the seed. | Positional parity makes lifted internal edges agree; different existential witnesses cannot be multiplied. |
| [2303 connected paths](../square_sum_all_size_connected_paths_20261019/ALL_SIZE_CONNECTED_PATHS.md) | Every N>=6294221 has all2303 fixed three-vertex paths plus the seed. Some seed/path edges overlap compatibly. | Shared ports and simultaneous versus disjoint predicates. Later threshold and recursive-ray work is stronger in different directions. |
| [THM-4850](../../../01-canon/theorems/THM-4850-ninja-pairs-of-the-square-sum-graph-start-at-42-and-the-base-49-numeration.md) | Weak-B pairs start sharply at42; Gerbicz-normalized nice pairs start41. These types must not be conflated. | Canonical radix49 digits24..72; a reducing address does not imply vertex insertion. |
| [THM-4851](../../../01-canon/theorems/THM-4851-the-prehistory-of-the-square-sum-numbers-unique-cycles-at-32-and-33-terminal-spectrum-with-the-hole-at-39.md) | Unique cycles32/33; terminal cycles at34..38 and>=40, with real hole39. | Terminal predicates can fail inside an interval of ordinary existence. |
| [Session124 ternary schedules](../session124_mixed_radix_square_rays_20261031/TERNARY_SCHEDULE_THEOREM.md) | Three actual prefix operators A,B,C; every coded prefix extends on all N>=(85R-1)/2. | Mixed radix25/49; exact word decoding and prefix antichains. Sparse infinite rays are not Hamiltonian rays of all positive integers. |
| [Overnight a](../square_sum_spectrum_overnight_20261031a/README.md) | The unchanged eleven-vertex path from the old24-obstruction has exact eventual threshold72. | A former hostile becomes a protected positive seed after its boundary deficit is paid. |
| [Overnight e](../square_sum_spectrum_overnight_20261031e/README.md), proof sections1--3 | Same-phase two-choice bases112..2811 make every signed offset independent: 2^R cycles at N>=(225R-1)/2. Uniform exponential bound and liminf log2(count)/N>=2/225. | Full phase functions, not only parity within one pair. This supersedes polynomial and stretched-exponential multiplicity bounds. |
| [Overnight f](../square_sum_spectrum_overnight_20261031f/README.md), proof sections1--3 | Fully refine the switch code: a subcubic spanning subgraph has exactly2^b cycles. Edge constraints become consistent bit literals. | At N70312,25 tied bits and625 free bits have the SAME union graph. A code need not be closed under its own graph's legal operations. |
| [Overnight g port orders](../square_sum_spectrum_overnight_20261031g/PORT-ORDERS.md), scoped result recovered through its handoff | Rooted reference-port order constraints cost at most floor(k/2) bits in that carrier; shared edge/order demands must agree. | Order, weights and incidence must use one bit assignment. This is not an arbitrary full-graph selector theorem. |
| [Overnight i](../square_sum_spectrum_overnight_20261031i/README.md) and [j](../square_sum_spectrum_overnight_20261031j/README.md) | One-copy same-square-class joins use a marked nontorsion elliptic point and finitely many forbidden scale ratios; tree gluing preserves jointly conditioned counts. | Denominator clocks, square-class descent, Farey rectangles, arithmetic collision avoidance. Sources are explicitly square-scaled, not numerically unchanged. |
| [Overnight l correction](../square_sum_spectrum_overnight_20261031l/README.md) | Cyclic products of local Hamiltonian cycles may miss global cycles. Prism+C4 has3 global cycles but2 such products. | Retain every spanning linear-forest boundary state, including states that reclose locally to multiple cycles. |

The all-N exponential and closed-subgraph theorems are the starting point for
future square-sum work. Reproving a long finite interval or merely embedding
a sparse path would be weaker. The unresolved extension is arbitrary
compatible prescriptions / interval absorption under authenticated boundary
relations, not ordinary Hamiltonicity.

## Older and less frequently used structural neighbours

1. **July28 discriminants and graceful labels.**
   [THM-2761, graph-edge-sum-discriminant-codegree-factorization-and-graceful-sign-gauge](../../../01-canon/theorems/THM-2761-graph-edge-sum-discriminant-codegree-factorization-and-graceful-sign-gauge.md)
   splits edge-sum collisions into adjacent-edge codegree factors and
   disjoint-edge additive-energy factors. A bipartite sign gauge converts
   sums to oriented differences; a mirror factor is still required to detect
   equal absolute differences. The path labels0,1,2 are the minimal hostile.
   This is directly relevant to the present distinction between incident
   wedges P21 and disjoint matching terms, but it is not the same polynomial.

2. **July28 universal rulers versus selected edges.**
   [THM-2786, binary-golomb-universal-edge-difference-separation-and-graceful-compression-boundary](../../../01-canon/theorems/THM-2786-binary-golomb-universal-edge-difference-separation-and-graceful-compression-boundary.md)
   decodes pairs by valuation and gives quadratic Sidon height. Universal
   pair separation must pay for all binomial(n,2) differences; graceful
   labeling asks only for n-1 selected edges. [THM-2765, rooted-nullstellensatz-linear-range-distinct-edge-labeling](../../../01-canon/theorems/THM-2765-rooted-nullstellensatz-linear-range-distinct-edge-labeling.md) already has a stronger
   tree-specific bound, so the general construction must not replace it.
   The analogous square-sum distinction is all native continuations versus
   one fixed transition catalog.

3. **August12 Pell selectors and clock positions.**
   [THM-3335, square-triangular-pell-markov-pythagorean-selector](../../../01-canon/theorems/THM-3335-square-triangular-pell-markov-pythagorean-selector.md)
   retains an actual Pell/Markov/Pythagorean parameter map and exhibits LRC
   rows with identical coarse norm/owner data but different extrema. Its
   tournament arc-count and design-order maps are different. The useful
   bridge is a selector retaining position, not the recurrence of a number.
   No LRC result follows from the new graph certificate.

4. **August14 branch averages versus branch order.**
   THM-3357, `../../../01-canon/theorems/THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit.md`,
   recovered through its August reflection, distinguishes Walsh level sums,
   signed sibling circuits and the actual ordered word. This is a provenance
   lead, not a dependency audited afresh here. It motivates retaining the
   ordered root word alongside a compressed spectrum.

5. **October19 graceful exchanges.**
   [Graceful/square exchange note](../graceful_square_exchange_20261019.md)
   proves midpoint rotations preserve long edges, a double-end splice escapes
   that obstruction, and fixed-endpoint graceful two-edge switches are
   impossible at EVERY size. Square-sum switches instead first appear at46.
   Thus the superficially same exchange has a permanently different arithmetic
   obstruction. The three-edge graceful connectivity conjecture remains open.

6. **October19 induced embeddings.**
   [Induced forest theorem](../induced_square_forests_subdivisions_20261019/INDUCED_SQUARE_FORESTS.md)
   uses square-gap factorization for codegree and gives near-quadratic forest
   height. Its separate [mod7 construction](../induced_square_forests_subdivisions_20261019/SUBDIVISION_MOD7_CONSTRUCTION.md)
   proves quartic-height induced two-subdivisions of every finite graph.
   The former theorem's additional candidate sections must be distinguished
   from the separately reviewed final subdivision theorem. The elementary
   fixed-seed bridge is re-proved in our NATIVE_CONTEXTS.md, so no unpromoted
   candidate is a hidden dependency.

7. **October19 AP triangles, elliptic curves and odd zeta.**
   [AP square/zeta note](../ap_square_odd_zeta_20261019.md) gives an exact
   Dirichlet-series bridge between signed AP square-sum triangles and
   primitive Pythagorean hypotenuses. Positive graph labels require an angular
   inequality that the height alone forgets. The first positive seed is
   482,3362,6242, explaining why this same seed later carries interval
   constructions. The formulas involving zeta(5),zeta(7),zeta(9) give no new
   irrationality theorem; reinstating the positivity cone changes the series.

8. **October23 forests and the golden cat map.**
   [Session32 cat-map proof](../session32_square_monster_20261023/CAT_MAP_PATH_FORESTS.md)
   gives the integral conjugacy of [[3,-1],[1,0]] and [[2,1],[1,1]]. Rooted
   path-forest counts are F_(2m), whereas traces are L_(2m). The values7 and21
   are different observables, not explanations of forbidden tournament H.
   Replacing the branched five-vertex block by a path changes52 to55.

9. **October23 continuation algebra.**
   [Session33 PORT_CONNECTION_ALGEBRA](../session33_continuation_states_20261023/ports/PORT_CONNECTION_ALGEBRA.md)
   ALREADY proves every different port pairing distinguishable by an abstract
   context. Deterministic state count differs from linear connection rank;
   F2 compression can erase a positive even count. The native arithmetic
   strengthening in this packet combines that theorem with item6 rather than
   reproving its weaker abstract form.

10. **October23 covariance and actual endpoints.**
    [Session40 path-scatter theorem](../session40_marked_ports_square_sum_20261023/PATH_SCATTER_AND_MARKING.md)
    identifies log-partition Hessians with edge covariance. At Q32 its global
    affine rank45 falls to12 within endpoint sectors; two extra conditional
    conserved edge costs remain beyond degree potentials. These are exact
    finite ranks, not a universal charge law or a minimum bit count.

11. **October23 square clocks and shared integer lifts.**
    [Session42 synthesis](../session42_shared_lattices_square_colors_20261023/SYNTHESIS.md)
    and [shared-clock theorem](../attachment_lattices_square_sum_20261023/square_sum/SHARED_CLOCKS.md)
    retain x_i=(-1)^i t+c_i, interval bounds, collision holes, and one common
    integer root per attached component. Odd/even words have one/two possible
    Hamiltonian realizations; first square-clock twins occur36, palindromic
    twins48. Marginal feasibility and residue agreement fail to synchronize
    actual ports. An odd loop can pin a half-integer or an illegal repeated
    vertex, so closure alone is insufficient.

12. **October23 ancestry and Fibonacci addition.**
    The same Session42 proof synthesis separates founders1,2,4 from later
    components: every n>=71 has all three ancestries;70 is the last exception.
    Parent-selector forests are not uniformly distributed colorings. Canonical
    Fibonacci colors give exact rotation-interval convolution laws, but
    shared target selectors cannot be multiplied as independent marginals.

13. **October24 Fortunate saturation.**
    [Session52 theorem](../session52_saturation_residues_metagraph_20261024/fortunate-square-pattern/THEOREM.md)
    identifies the tempting17-bit match precisely and refutes it at the next
    index: F18=109>67 while Q27 is Hamiltonian. Its surviving graph quantity
    is path-cover deficit above component count. This prevents treating a
    visually matching run word as a universal arithmetic law.

14. **October25 square-magic incidence.**
    [Session57 synthesis](../session57_square_sum_magic_paths_20261025/SYNTHESIS.md)
    distinguishes pair-square edges from fixed-sum triple hyperedges. A
    decorated K4,4 preserves sixteen square forms; a bare biclique does not
    retain the affine magic constraints. The mod47 Lagrange-family failure
    excludes one construction, not all magic squares. The compatible doubled
    elliptic-coordinate condition is another shared-incidence obligation.

15. **October31 parity seams and November1 machines.**
    [Session131 square proofs](../session131_reflection_numeration_moonshine_20261031/square_sum/PROOFS.md)
    use mod8 channels to bound even-root seams; ordered reflections retain
    more than the root multiset. [Session135 machine theory](../session135_fractran_tm_structure_20261101/square_machines/THEORY.md)
    gives a literal square predicate and permutation verifier. The repeated
    (1,8) word passes local squares and moments but not vertex coverage.
    Verification clocks do not themselves produce a universal constructor.

The July20 shear-catalog messages were inspected as a provenance lead.
Their 'square sums' are array sums, not adjacent integer squares. They are
retained in the discovery index but not counted as an independent Hamiltonian
theorem or evidence of a shared recurrence.

## Connections that survived exact typing

| Source -> target | Map / preserved predicate | Information lost without a sidecar | Cheapest decisive test / result |
|---|---|---|---|
| Abstract port separator -> native square-sum extension | Replace each connector by its explicit fixed-label length-three bridge; selected-state extendability is preserved | Whole initial-interval coverage; not supplied | Four-cycle labels1,3,35,46: one8-cycle versus two4-cycles; PROVED all-size in NATIVE_CONTEXTS.md |
| Cycle expansion -> ordinary endpoint graphs | Pay signed directed-cycle couplings, then retain q_e=AijAji and its matching/wedge polynomial | Orientations before payment; distinguished-cycle multiplicities | (2,4,4) detects the missing factor; all12,005,168 endpoint graphs pass the repaired reduction |
| Square-clock network -> signed incidence / graceful gauge | Edge equations propagate an integer affine root with signs | Absolute values, integer pins, forbidden collisions and global coverage | 4,9,16 closes at11/2; 4,9,9 closes at an illegal loop; inherited exact tests |
| Phase-fixed recursive paths -> independent Boolean choices | Fix every label's positional parity; each offset becomes a disjoint readable carrier | Synchronization ties and source position | N70312 same graph,25 versus625 bits; primitive refinement restores completeness |
| AP square triangle -> hypotenuse-height series | Explicit invertible Pythagorean formula, with gcd/scale retained | Angular positivity and Gaussian factor allocation | (3,4,5) gives a negative label; old odd-zeta formulas cannot be used unfiltered |
| Graph/tournament observer charts -> averaged coupling | Retain ordering signs and orbit weights | Actual source direction and path information | Session42's exact8/9 coupling gap is not a spectral theorem for the native metagraph |

These are maps and failure boundaries, not a claim that all listed structures
are isomorphic. The two new useful outcomes are order-ten stability and
native source-preserving separation; arbitrary-interval prescription,
universal odd tournament-H completeness, and Collatz convergence remain open.

The clean-checkpoint fetch then supplied additional square-sum work, including
the prime-label construction and exact guarded tournament response at Q25.
[INCOMING.md](INCOMING.md) records their separate verification, the singleton
scope correction, and how the new triangle surplus changes the next target.
