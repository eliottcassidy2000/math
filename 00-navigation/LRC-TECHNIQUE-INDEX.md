# LRC Technique Index

This is the shared pull-from / contribute-to index for LRC and LRC14 proof
technology.  It is deliberately broad: tournament analysis, metagraphs,
sequence transforms, harmonic carriers, endpoint topology, arithmetic packets,
and proof-zipper templates all belong here when they preserve a usable LRC
predicate.

Use this file as a technique menu before starting a new pass.  Add entries when
you introduce a genuinely new carrier, repair a failed scalar quotient, or show
that an old tournament/metagraph idea gives a named LRC packet, obstruction, or
certificate.

## Contribution Protocol

Every new entry should answer these questions.

```text
ID: LTI-###
Name:
Carrier / vertex set:
Pairwise observable:
Binary relation or gauge:
Preserved LRC predicate:
Destroyed information:
Best use:
Failure mode:
Anchors:
Next agent hook:
```

Do not default to runners as tournament vertices.  Consider proof carriers,
safe components, gaps, sections, endpoints, residues, wall crossings, exact
period packets, Fourier modes, Haar tiles, quotient fibers, state-lift sectors,
and proof obligations first.

## Current Assembly Principle

The live LRC14 proof object is a labelled packet over several compatible
coordinates:

```text
exact M/Farey node
qdiv / denominator shell
regular-open safe component
closed boundary atom
endpoint owners
C27 / unital / K33 state labels
Fejer / Toeplitz certificate data
Ramanujan / exact-period denominator labels
Haar / Baire / discrepancy tile coefficients
boundary-moment chart data
Kaczynski / analytic smoothing approach class
state-lift or F7 residual sector
```

A quotient may forget a coordinate only when the LRC predicate is constant on
fibers, the coordinate is reconstructible, a dual certificate annihilates the
forgotten direction, or the loss is routed to a named residual sector.  This is
the HYP-2990 no-free-slider rule.

## Agent Workflow

1. Pick an `LTI-*` technique and one target packet family.
2. State the vertex set before building a tournament or metagraph.
3. State the observable and switch/gauge that turns pair data into a relation.
4. Record tournament fingerprints when feasible: score histogram, SCCs,
   directed cycles, edge flips, Hamiltonian paths, and tie path.
5. State what the technique preserves and what it destroys.
6. If it fails, add the failure as a quotient guardrail instead of deleting it.

## Technique Atlas

### A. Tournament And Metagraph Carriers

| ID | Technique | Carrier / vertex set | LRC pull | Anchors | Next agent hook |
|----|-----------|----------------------|----------|---------|-----------------|
| LTI-001 | Proof-carrier tournament | Proof obligations / certificates | Rank proof routes by retained LRC predicate rather than by row statistics | HYP-1941, HYP-2987, HYP-2990 | Build a carrier tournament for any new proof pass before making scalar claims. |
| LTI-002 | Gauge / switchboard analysis | Binary gauges over carriers | Turn pairwise observables into relation switches and edge-flip ledgers | HYP-1940, HYP-1941 | Declare the switch and tie Hamiltonian path in every computational note. |
| LTI-003 | Hamiltonian path count `H` as shadow | Tournament iso classes | Use `H` as a spread/complexity shadow only after marking observer, threshold, or packet labels | HYP-1970, HYP-1971, THM-370 | Audit when `H` predicts LRC gap data and when endpoint labels break it. |
| LTI-004 | OCF independence polynomial | Odd-cycle compatibility graph | Transfer tournament cycle packing to LRC packet packing without scalarizing cycles | THM-081, HYP-2618, HYP-2990 | Test whether an LRC residual is an activity/coimage sector. |
| LTI-005 | Walsh / Fourier tournament spectrum | Arc monomials / cycle indicators | Identify which low-degree Fourier modes are lost by a quotient | HYP-2426, `H_fourier_structure.py`, `walsh_legendre_proof.py` | Attach Walsh support to Fejer or Haar packet certificates. |
| LTI-006 | Deletion-contraction for `H` | Edges and contractions | Model LRC packet removal / state contraction as a recurrence with boundary debt | `dc_ocf_tracking.py`, THM deletion-contraction notes | Try a packet deletion-contraction identity for boundary atoms. |
| LTI-007 | Trace / closed-walk overlap correction | Cycle traces plus support overlaps | Separate true cycle structure from trace overcount, analogously separating component count `K` from Haar switch labels | HYP-2498 | Use overlap corrections when raw moment or trace counts look promising. |
| LTI-008 | A000568 quotient and source lift | Unmarked, rooted, source-deleted tournaments | Raw iso class is too lossy; observer-source lift is exact for loneliness predicate | HYP-1977, HYP-2486 | Pair source-deleted class with exact endpoint/Farey packet data. |
| LTI-009 | Observer-source tournament | Observer plus runner vertices | LRC witness iff observer is a source; blocker count equals observer indegree | THM-381, THM-385, HYP-1988 | Lift LRC14 packets to observer-source deficit layers. |
| LTI-010 | Threshold-decorated class fiber | Tournament class plus threshold colors | Preserve equality walls, gap lengths, and observer-adjacent threshold data | HYP-1982, THM-382, THM-383, THM-384 | Use for q=14 wall packets and compactified source-gap forcing. |
| LTI-011 | Pair-cell operation grid | Unordered runner pairs | Detect additive, dyadic, product, and same-chain branch structure at second order | HYP-1976 | Use pair-cell vertices for C27/K33/Farey branch audits. |
| LTI-012 | Section-boundary functors | Circle sections / boundaries / void sectors | Replace runner vertices by fixed regions, events, or proof obligations | HYP-2024, `lrc_section_boundary_functors_s539.py` | Re-run difficult row families with section vertices. |
| LTI-013 | Metagraph walk | Structural classes as nodes, flips as edges | Treat proof progress as movement in a class graph with boundary labels | `metagraph-as-transfer-chain.md`, `directed-graph-metagraph.md` | Define a metagraph for F0-F7 packet transitions. |
| LTI-014 | Blue/black boundary parity | Metagraph color layers | Keep GF(2) boundary data when raw color subgraphs fail parity claims | HYP-2250 | Use color-boundary parity for packet flip ledgers. |
| LTI-015 | Good-cut / SCC support | SCC support intervals | Convert support structure into zipper coordinates, not raw connectivity | THM-354, HYP-2990 | Test good-cut support as a state-lift coordinate. |
| LTI-016 | Rado / monodromy tournament | Infinite or universal tournament models | Use as branch-locus language for families of local states | Rado reflections, HYP-2127 | Only apply after exact finite packet coordinates are declared. |
| LTI-017 | VT / cyclic rigidity split | Vertex-transitive classes | Separate cyclic polygon rigidity from dihedral bracelet and nonabelian Cayley mesh | HYP-2125 | Use for AP uniqueness claims; avoid calling all VT structure polygonal. |
| LTI-018 | Tournament spectrum `Sigma(S)` | Isomorphism-class path over phase | Restore magnitude awareness lost by a single apex tournament | HYP-2928, HYP-2927 | Compute binding-scale migration for new near-miss rows. |
| LTI-019 | Apex tournament census | Apex winding tournaments | Gives residue-layer necessary conditions, not tightness by itself | HYP-2924, HYP-2927 | Pair apex iso class with Farey denominator scale. |
| LTI-020 | Metatournament / repo-as-tournament | Research claims / routes as vertices | Use route conflicts to decide which proof carrier dominates | `the-meta-tournament.md`, `repo-as-tournament.md` | Add meta-edge reasons when two carriers disagree. |

### B. LRC Packet And Reduction Carriers

| ID | Technique | Carrier / vertex set | LRC pull | Anchors | Next agent hook |
|----|-----------|----------------------|----------|---------|-----------------|
| LTI-021 | Labelled packet classifier | Primitive rows with route labels | Classify into q-witness, AP/GW, petal, K33, covering, unknown | HYP-2963 | Extend row banks only with retained packet labels. |
| LTI-022 | Family/sporadic decision tree | F0-F7 packet families | Converts "counterexample" into named residual cases | HYP-2961, HYP-2956 | Give F7 a formal harmonic/state-lift definition. |
| LTI-023 | Fixed-margin packet signature | Fixed row/column margins and swaps | Import local-swap proof architecture while preserving packet fibers | HYP-2962, arXiv-pattern notes | Build swap connectivity inside a labelled LRC family. |
| LTI-024 | Packet migration gauntlet | qdiv>=14 exact packets | Show residuals migrate to positive Haar/open fronts except AP/GW | HYP-2955 | Re-run on larger or structured source families. |
| LTI-025 | Source-spectrum pullback | Farey-indexed rooted packet spectrum | Pull global proof data back to source packets before quotienting | HYP-2954 | Attach source-spectrum fingerprints to F6/F7 candidates. |
| LTI-026 | Skeleton gate | qdiv=14 boundary skeletons | Separate AP/GW skeleton from impostor shells and loose rows | HYP-2960 | Use as first filter before K33/state-lift claims. |
| LTI-027 | Forbidden H7 state lift | State-lift obstruction sector | Route K33/nonunit residuals to THM-572-style forbidden lift | HYP-2908, THM-572 | State exact lift data emitted by each failed carrier. |
| LTI-028 | Few-apex lift packets | Few changed apex runners | Keep small lifted packet structure before wide analytic estimates | HYP-2968 | Check whether few-apex residuals already have strict open mass. |
| LTI-029 | Nork pinch templates | Pinched endpoint packet templates | Name local pinches rather than storing them as anonymous exceptions | HYP-2966 | Add Nork template labels to boundary ledgers. |
| LTI-030 | Boundary-gap packet bridge | Safe interval endpoints and owner transitions | Convert positive safe intervals into exact boundary bridges | HYP-2965 | Prove zero bridge packets have moment or state-lift discharge. |
| LTI-031 | Boundary-moment packet ledger | Exact-period denominator charts | Avoid false obstruction from one all-covered chart; use multi-chart ledger | HYP-2969 | Couple with Haar tile coefficient detection. |
| LTI-032 | Danger-count moment dual | Moment duals over danger counts | Use count moments as labelled dual certificates, not raw scalar counts | HYP-2971, HYP-2970 | Attach moment duals to packet signatures. |
| LTI-033 | Endpoint-credit winding cycle | Endpoint owner graph | Preserve active owner cycle data at boundary equality | HYP-2970 | Compare with tope/cocircuit boundary atoms. |
| LTI-034 | Taut bridge graph | Positive bridges vs zero-length taut vertices | AP/GW are taut equality atoms; positive rows have bridges | HYP-2975 | Combine with kernel support radii and Fejer packets. |
| LTI-035 | Kernel homotopy ledger | Safe components / boundary defects | Strict components persist under smoothing; zero-open emits defect | HYP-2984 | Use support radius as admissibility precondition for smoothing. |
| LTI-036 | Certificate handoff atlas | Proof carriers as vertices | Glue q-witness, Fejer, Ramanujan, endpoint, twist, moment, analytic, state-lift certificates | HYP-2987 | Prove the six atlas arrows O1-O6. |
| LTI-037 | Abstract zipper theorem | Two toothed carriers over packet base | No-free-slider rule for all quotients | HYP-2990 | Turn every technique in this index into a zipper tooth. |
| LTI-038 | Alternating gauntlet | Ordered proof teeth | Iterate carrier tests until certificate, equality atom, or forbidden lift | HYP-2990, HYP-2956 | Define a well-founded debt order for F0-F7. |
| LTI-039 | Harmonic residual sector | Named final failure bucket | Replace anonymous F7 by representation/homology/cocircuit/state-lift sector | HYP-2990 | First priority: define F7 predicate and destroyed data. |
| LTI-040 | Holistic proof-object sheaf | All retained packet coordinates | Keep exact scale, topology, endpoint, state, harmonic, and dual labels together | HYP-2976 | Use as checklist before declaring any reduction complete. |

### C. Geometry, Topology, And Event Algebra

| ID | Technique | Carrier / vertex set | LRC pull | Anchors | Next agent hook |
|----|-----------|----------------------|----------|---------|-----------------|
| LTI-041 | Borel/Baire/Haar event carrier | Danger/safe finite arc unions | Safe sets are Borel, Baire, Haar measurable; boundary debt is finite | HYP-2948, HYP-2949 | Separate strict-open witnesses from closed equality atoms. |
| LTI-042 | Haar-Baire taut wave | Regular-open set, Haar mass, boundary debt | Any-angle path-planning analogy for propagating safe fronts | HYP-2949 | Add endpoint owner labels to every wavefront. |
| LTI-043 | Boundary-owner skeleton | Active endpoint-owner pairs | AP/GW share skeleton; Haar alone cannot distinguish C27 transfer | HYP-2951 | Prove skeleton rigidity for zero-open packets. |
| LTI-044 | Tope/cocircuit wall | Endpoint arrangement cells | Open all-safe cells are strict witnesses; all-safe endpoints are boundary cocircuits | HYP-2986 | Classify no-tope/no-cocircuit packets. |
| LTI-045 | Oriented-matroid wall language | Topes, cocircuits, circuits | Use endpoint arrangement combinatorics before analytic smoothing | HYP-2986 | Search for minimal forbidden wall packet. |
| LTI-046 | Haar-product square | 2-by-2 Haar checkerboard | Fixed-margin switch changes mixed coefficient, so row/column margins are lossy | HYP-2989 | Count independent color-compatible mixed switches. |
| LTI-047 | Haar rectangle product atlas | Dyadic rectangles and interaction classes | Owner strips, cross handoffs, nested refinements carry discrepancy signal | HYP-2988 | Test non-AP/GW zero-open residuals for nonzero coefficients. |
| LTI-048 | Colored discrepancy reservoir | Color grids mod 14 | Replace raw component count by resonance-aware discrepancy bound | HYP-2593, HYP-2594, HYP-2595 | Fuse with Haar mixed-product switches. |
| LTI-049 | Safe component support radius | Exact open intervals | Gives local smoothing radius and familywise stability target | HYP-2984, HYP-2975 | Compute support radii for new packet families. |
| LTI-050 | Section check-off Hall graph | Section supports and matchings | Translate one-runner-per-section into Hall/wall-switch problem | HYP-2570 | Apply to C27 and covering residual packets. |
| LTI-051 | Compactified source-gap forcing | Observer-adjacent threshold gaps | Observer lonely iff adjacent gaps are threshold-long in compactified model | THM-384, HYP-1986 | Try compactified forcing on q=14 wall packets. |
| LTI-052 | Void-pressure / boundary-flux | Fixed sectors and events | Measures pressure through sections rather than runners | HYP-2024 | Combine with tope/cocircuit wall scan. |
| LTI-053 | Baire rank recombination | Measurable rank packets | Recombine measurable packet fronts after C27/K33 labels survive | HYP-2947 | Link with Haar tile discrepancy. |
| LTI-054 | Any-angle algorithm analogy | Visibility / lazy verification / finite atlas | Converts path-planning carriers into proof roles | HYP-2948, HYP-2949 | Use only when it preserves event labels. |
| LTI-055 | Unit-distance endpoint ears | Unit-distance graph ears and endpoints | Cross-domain zipper for endpoint extension blockers | HYP-2620, THM-408, HYP-2990 | See if ears model LRC endpoint owner extension debt. |

### D. Harmonic, Fourier, Analytic, And Certificate Carriers

| ID | Technique | Carrier / vertex set | LRC pull | Anchors | Next agent hook |
|----|-----------|----------------------|----------|---------|-----------------|
| LTI-056 | Fourier-Toeplitz PSD dual | Fejer / Toeplitz matrices | Positive rows have PSD-dual violations; AP/GW are equality atoms | HYP-2974 | Formalize familywise PSD certificates. |
| LTI-057 | Fejer interval certificate | Packet-keyed rational interval data | Convert floating negatives into checkable certificates | HYP-2981 | Replace prototype trig intervals with formal backend. |
| LTI-058 | Fejer packet manifest | Stable certificate records | Retain packet key, center, degree, atom formula, interval sign, route | HYP-2981 | Extend manifest to K33/petal/covering families. |
| LTI-059 | Robbins no-bridge guardrail | Graph bridges / certificate assembly | Forgotten certificate fields are bridges unless reconstructed | HYP-2981 | Audit new quotients for bridge loss. |
| LTI-060 | Robin / divisor scalar warning | Sigma/tau/divisor shadows | Divisor scalars help pre-split but cannot replace packet labels | HYP-2981, HYP-2978 | Use divisor features only with lost-label repair. |
| LTI-061 | Ramanujan exact-period projectors | Exact-period denominators | Restore prime-power and repeated-prime labels erased by squarefree weights | HYP-2979, HYP-2978 | Build endpoint-aware exact-period projectors for late denominators. |
| LTI-062 | Ramanujan-divisor quotient guardrail | Divisor / unit / exact-period fibers | Quotient admissible only if predicate constant, reconstructible, annihilated, or residualized | HYP-2978 | Run fiber-mixing checks on every new scalar feature. |
| LTI-063 | Analytic sieve packet weights | `Phi`, `G`, `mu`, prime sums | Large-sieve/Selberg weights are middle certificates with squarefree blindness | HYP-2982 | Attach exact-period labels to any analytic estimate. |
| LTI-064 | Kaczynski boundary approach | Boundary functions / approach classes | Treat exceptional sets as labelled approach classes | HYP-2983, HYP-2679 | Define approach class for true-wide residuals. |
| LTI-065 | Admissible smoothing dispatcher | Smoothing policies | Route packet families to the proof clock that retains labels | HYP-2985 | Prove dispatcher lemma over HYP-2963 fibers. |
| LTI-066 | Major/minor arc labelled template | Arc regimes and smoothing kernels | Imports circle-method architecture without losing endpoint/state labels | HYP-2982, HYP-2983 | Split late denominator walls into labelled arc regimes. |
| LTI-067 | Abel / signed decay | Off-resonance analytic tails | Discharge true-wide generic packets with approach-labelled decay | HYP-2985 | Quantify signed decay after packet labels are fixed. |
| LTI-068 | Selberg / large-sieve precondition | Squarefree upper-bound weights | Useful only as labelled precondition, not final equality carrier | HYP-2982, HYP-2985 | State which exact-period packets survive the quotient. |
| LTI-069 | Riesz product witness route | Positive trigonometric products | Potentially certificate lonely points through weighted test functions | T827, Riesz reflections | Revisit with packet labels and Fejer interval backend. |
| LTI-070 | Bonferroni / activation-depth ledger | Inclusion-exclusion layers | Use live Newton layer and missing-depth parity, not universal third-order truncation | T1016, T1017, HYP-2901 | Add activation depth to wide packet classifier. |
| LTI-071 | Relation-lattice theta form | Integer relation lattice | Interpret lonely measure as positive-definite theta sum / signed gas | T826, HYP-2540 | Couple relation height with product-kernel tail bounds. |
| LTI-072 | Subtorus relation-height split | Affine relations in far packet | Low relation height finite check plus high-height decay | HYP-2599 | Apply to true-wide Kaczynski/Abel packets. |
| LTI-073 | Resonance-channel Fourier split | Congruence-compatible frequencies | Kills false count terms and exposes colored resonances | HYP-2595, HYP-2867 | Merge with Haar product coefficient signs. |
| LTI-074 | Spectral shadow dual | Low-frequency shadows | Useful as pre-splitter; must retain endpoint/Farey labels | `lrc14_spectral_shadow_dual_codex_20260624.py` | Compare with Ramanujan exact-period projectors. |
| LTI-075 | Krawtchouk / orthogonal shadow | Orthogonal basis coordinates | Detect which shadows anonymize cycle or packet structure | Krawtchouk reflections | Use as anti-scalar guardrail for new bases. |

### E. Arithmetic, Farey, Sequence, And Number-Theory Carriers

| ID | Technique | Carrier / vertex set | LRC pull | Anchors | Next agent hook |
|----|-----------|----------------------|----------|---------|-----------------|
| LTI-076 | Exact `M=p/q` / Farey node | Binding rational | Always keep exact scale before derived products or powers | HYP-2931, HYP-2984 Farey scheduler | Require `M`, `q`, and `e=14p-q` in packet records. |
| LTI-077 | Farey mutation scheduler | Mutations of `p/q` | Product mutation is route index only after unit-excess gate | HYP-2984 Farey scheduler | Use `e=0,1,>1` to dispatch route. |
| LTI-078 | Farey product / K33 incidence | Factor fibers and bipartite minors | `p>=3` unit-excess branch signals K33/state-lift/Fejer route | HYP-2932, HYP-2934 | Preserve factor fiber, not just product value. |
| LTI-079 | Summand/multiplicand bridge | Additive and product graphs | Separate Farey summand shell from product incidence wall | HYP-2934 | Compare C27 p=2 branch with K33 p>=3 branch. |
| LTI-080 | C27 shell-transfer spectrum | Antipodal summand shells | Mark hole/double transfer labels before exact Farey branch | HYP-2937 | Prove second-gap transfer rigidity. |
| LTI-081 | C27 two-swap splice | Two-hole second-gap rows | Known splices combine petal with GW/K33 branches | HYP-2940 | Classify all second-gap two-hole residuals. |
| LTI-082 | q=3 unital block lift | Pair-completion blocks | Use branch-local unital charts, not global scalar design claims | HYP-2942 | Attach unital chart IDs to petal/K33 packets. |
| LTI-083 | Pi / flower / unital alias guardrail | Normalization labels | Prevent accidental equivalence of radians, turns, and unital meanings | HYP-2938 | Add normalization field to analogy-driven passes. |
| LTI-084 | Perfect-product Farey guardrail | Even-perfect product chains | Product value is meaningful only with Farey determinant and minor labels | HYP-2945, HYP-2946 | Use perfect-product analogies as carrier roles, not averages. |
| LTI-085 | Totient curvature | Euler factor residuals | Tournament recurrences fail on `phi` by labelled CRT curvature | HYP-2900, HYP-2899 | Feed curvature labels into exact-period packet ledgers. |
| LTI-086 | Mobius / product packet ledger | Multiplicative packet data | Track squarefree, unit, and product support separately | HYP-2899, HYP-2982 | Avoid using `mu^2/phi` without prime-power repair. |
| LTI-087 | Zeta `-1/12` finite core | One-tail AP core | Regularized infinite signal becomes finite divisibility and binding margin | HYP-2896 | Use as warning against regularized scalar shortcuts. |
| LTI-088 | Pisano / Fibonacci band ladder | Modular recurrence quotients | Sequence periodicity can mark fibers but may be partition-bound | Pisano reflections, HYP-2523 | Require band and fiber labels before using recurrence data. |
| LTI-089 | Zeckendorf boundary normal form | Zeckendorf packets | Potential normal-form carrier for boundary decompositions | HYP-1902 | Revisit only with qdiv/Farey labels attached. |
| LTI-090 | Euler / Glaisher / Witt parity | Partition and parity transforms | Parity identities can expose doubled/forgotten boundary labels | HYP-2685, Glaisher/Witt reflections | Test parity identities on exact packet fibers. |
| LTI-091 | Pollock defect-pair descent | Tetrahedral defect pairs | Defect-pair carry ledgers model forbidden residual descent | HYP-2491, HYP-2497 | Use as analogy for F7 defect-pair closure. |
| LTI-092 | Heegner / prime boundary | Prime boundary tournaments | Number-theory horizons can serve as labelled carrier boundaries | HYP-2225, HYP-2226, HYP-2227 | Use only after identifying preserved LRC predicate. |
| LTI-093 | Adelic / 7-adic split | Local primes and archimedean data | Separate local congruence gates from real gap geometry | 7-adic reflections, `lrc14_7adic_archimedean_split_s6.out` | Attach local-prime gate to denominator wall packets. |
| LTI-094 | Singular series / circle method | Arithmetic density factors | Good for middle-layer estimates after packet and endpoint labels exist | singular-series scripts, HYP-2982 | State major/minor packet labels explicitly. |
| LTI-095 | Freiman / additive energy | Small-excess finite cores | Finite additive structure can compress small-excess packets | HYP-2638, HYP-2635 | Pair with Fejer family certificates. |

### F. Formalization, Computation, And Verification Carriers

| ID | Technique | Carrier / vertex set | LRC pull | Anchors | Next agent hook |
|----|-----------|----------------------|----------|---------|-----------------|
| LTI-096 | Lean divisibility core | Formal arithmetic statements | Separate proved local arithmetic from shell-collapse conjecture | HYP-2929, `LRCApexShell.lean` | Formalize shell collapse or state-lift alternative. |
| LTI-097 | Exact rational interval unions | Rational endpoints | Prevent floating topology errors in safe/danger component audits | HYP-2948, HYP-2975, HYP-2984 | Use exact endpoints for every zero-open claim. |
| LTI-098 | Full-bank Fejer audit | Large finite row banks | Turns finite positive rows into PSD-certificate candidates | HYP-2974, HYP-2981 | Convert floating certificates to interval payloads. |
| LTI-099 | Family compression prototype | Parameterized packet families | Replace named-row certificates with reusable templates | HYP-2987 O3, HYP-2990 | Build one template for K33 or petal family. |
| LTI-100 | Stale-quotient fiber-mixing test | Feature fibers | Prove a scalar is unsafe by finding mixed LRC routes in one fiber | HYP-2978 | Add a fiber-mixing table to every scalar proposal. |
| LTI-101 | Edge-flip / order-flip stress | Pairwise row order under scalar | Refute order scalar candidates before they enter proof route | Farey scheduler, HYP-2984 | Store flip counts for any new ranker. |
| LTI-102 | Bounded atlas plus escape route | Finite stress bank and theorem debt | Make bounded computation useful by naming the infinite escape lemma | HYP-2963, HYP-2981 | Never leave a finite audit without its family theorem target. |
| LTI-103 | Result fingerprint convention | Scores, SCCs, cycles, paths | Standardizes Tournament Analysis output across agents | AGENTS.md, HYP-2987, HYP-2990 | Add this fingerprint to new scripts by default. |
| LTI-104 | Session checkpoint discipline | Small coherent commits | Prevent ID collisions and preserve incoming signal | CONCURRENT-SESSIONS.md | Claim IDs with honest stubs before large runs. |
| LTI-105 | POKE/forum synthesis | Shared narrative artifacts | Keeps active theorem targets visible to parallel agents | comms/POKE-FORUM.md, comms/POKE-COORDINATION.md | Post a short technique delta when updating this index. |

## Cross-Cutting Guardrails

- A scalar that works on named rows is not a proof carrier until its fibers are
  checked for mixed LRC routes.
- A tournament invariant on runners is usually a shadow.  Ask what marks,
  thresholds, endpoints, or packet labels were forgotten.
- AP/Goddyn-Wong are boundary equality atoms.  K33 is a state-lift / residual
  route, not the same kind of atom.
- Exact `M=p/q`, `qdiv`, and Farey excess should be retained before any
  product, summand, power, divisor, or spectral quotient.
- Haar/Baire/Haar-product data separates strict-open witnesses from closed
  boundary debt; it does not replace C27, endpoint owner, or K33 labels.
- Fejer, Ramanujan, Toeplitz, Kaczynski, and analytic-sieve carriers should be
  used as packet certificates or admissible arrows, not as global scalar
  monotonicities.
- Every final residual must be named: harmonic sector, cocircuit, state lift,
  curl, homology class, representation component, or exact packet family.

## Open Contribution Slots

| Slot | Needed artifact | Suggested starting techniques |
|------|-----------------|-------------------------------|
| LTI-TODO-1 | Formal `F7` definition | LTI-036, LTI-037, LTI-039, LTI-073 |
| LTI-TODO-2 | Family Fejer certificate template for K33 | LTI-056, LTI-057, LTI-058, LTI-099 |
| LTI-TODO-3 | Haar tile coefficient test on labelled LRC packets | LTI-046, LTI-047, LTI-048, LTI-073 |
| LTI-TODO-4 | Section/tope wall theorem for no-tope/no-cocircuit packets | LTI-012, LTI-044, LTI-045, LTI-050 |
| LTI-TODO-5 | Exact-period Ramanujan repair for squarefree blindness | LTI-061, LTI-062, LTI-063, LTI-085 |
| LTI-TODO-6 | Kaczynski true-wide approach-class classifier | LTI-064, LTI-067, LTI-071, LTI-072 |
| LTI-TODO-7 | Source-spectrum binding-scale migration audit | LTI-018, LTI-019, LTI-025, LTI-076 |
| LTI-TODO-8 | Metagraph of F0-F7 transitions | LTI-013, LTI-022, LTI-038, LTI-103 |
| LTI-TODO-9 | Lean shell-collapse or state-lift formal route | LTI-027, LTI-096, LTI-102 |
| LTI-TODO-10 | Contributor script that emits a technique-card stub | LTI-103, LTI-104, LTI-105 |

## Minimal Entry Template

Copy this block when adding a new technique.

```text
### LTI-NEW: name

Carrier / vertex set:
Pairwise observable:
Binary relation or gauge:
Preserved LRC predicate:
Destroyed information:
Best use:
Failure mode:
Anchors:
Next agent hook:
```
