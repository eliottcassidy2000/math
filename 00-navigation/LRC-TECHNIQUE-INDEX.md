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
closed arc-Cech beta / open tope component count
boundary cocircuit facets
runner quotient Betti defect
endpoint owners
C27 / unital / K33 state labels
Fejer / Toeplitz certificate data
Ramanujan / exact-period denominator labels
Haar / Baire / discrepancy tile coefficients
boundary-moment chart data
Kaczynski / analytic smoothing approach class
state-lift or F7 residual sector
cochain / cocycle obstruction class
additive-basis regime / Pascal-slope row / carry-width labels
Farey operator lane (`root`, `p+q`, `p*q`, `q^p`, `p^q`)
curried call signature / lost-coordinate function
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
| LTI-047 | Haar rectangle product atlas | Dyadic rectangles and interaction classes | Owner strips, cross handoffs, nested refinements carry discrepancy signal | HYP-2992 | Test non-AP/GW zero-open residuals for nonzero coefficients. |
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
| LTI-106 | Haar zipper cocycle | `2 x 2` Haar/fixed-margin packet tables | Retain the mixed cocycle `zeta=T00-T01-T10+T11` when row/column margins collide | HYP-2991, HYP-2990, HYP-2989, HYP-2595 | Compute packet-level `zeta` signatures and route unpaired cocycles to F7/THM-572 debt. |
| LTI-107 | Haar product tiling synthesis | Haar rectangles and fixed-path staircase tilings | Treat product-rule descent and Walsh/Haar multiplicativity as quotient guardrails before using tiling shadows | HYP-2989, HYP-2988, HYP-2452, THM-351 | Build a packet grid whose product-rule classes feed color-compatible discrepancy and Fejer handoff arrows. |
| LTI-108 | Cocycle-sheaf exactness | Cocycle carriers / cochain complex `C0 -> C1 -> C2` | Treat every quotient as retaining a cocycle, cancelling it by exactness/duality, restricting it to a smaller packet, or routing it to named residual cohomology | HYP-3006, HYP-2991, HYP-2990 | Build the finite emitted-cocycle matrix on HYP-2963 packet banks and test exactness at `C1`. |
| LTI-109 | Packet cocycle atlas | Labelled packet fibers / quotient fiber graphs | Measures the signed coordinate lost by a quotient and demands exact, dual-annihilated, descended, boundary, reconstructed, or residual exits | HYP-2995, HYP-2991, HYP-2990, HYP-2978 | Define `omega_Q` for every HYP-2963 quotient and test the cocycle exits familywise. |
| LTI-110 | Cocycle obstruction atlas | Cochain carriers / proof obligations | Classify each forgotten coordinate as exact coboundary, closed cocycle, torsion/period class, or named residual | HYP-2994, HYP-2993, HYP-2991, HYP-2988, THM-572 | Emit `C0/C1/C2/H2` fields and sparse preserved-label warnings for every new quotient. |

### G. Recovered Tournament, Metagraph, Series, And Cohomology Carriers

These rows promote high-value methods from the long-form bank, older tangent
threads, and the S167 cocycle pass into compact handles.  They are intentionally
carrier-first: choose the vertex set, declare what the quotient preserves, and
record what it destroys before asking it to prove anything.

| ID | Technique | Carrier / vertex set | LRC pull | Anchors | Next agent hook |
|----|-----------|----------------------|----------|---------|-----------------|
| LTI-111 | Deck-derivative reconstruction | Deleted cards, child fibers, derivative traces | Repair global scalar quotients by asking what the deck of local deletions reconstructs | HYP-2936, HYP-2237, HYP-2246 | Build an LRC packet deck: delete one runner/obligation and test whether endpoint and route labels reconstruct. |
| LTI-112 | Burnside orbit-tax filter | Orbit representatives, stabilizers, fixed-point taxes | Use automorphism size and fixed-point side channels to detect when A000568-style quotienting hides LRC labels | HYP-1978, HYP-2228, THM-409, THM-417 | Add orbit-size/stabilizer tax to any apex or winding-tournament census. |
| LTI-113 | Merged metagraph quotient transport | Metagraph nodes, complement-merged edges, transport heights | Track what leaks across a merged quotient instead of treating the merge as lossless | HYP-1835, HYP-2245, `quotient-transport-and-good-cut-gas.md` | Build the F0-F7 transition metagraph with explicit lost-label transport arrows. |
| LTI-114 | Good-cut interval gas | SCC condensation intervals / support blocks | THM-354 turns base-path good cuts into SCC defect, giving a real support coordinate | THM-354, HYP-1779, HYP-1853 | Recast LRC endpoint protection as an SCC/good-cut support theorem over packet carriers. |
| LTI-115 | OCF activity/coimage sector | Odd-cycle conflict graph, activity packets, coimage classes | Transfer `H=I(Omega,2)` from scalar path count to compatibility/activity ledgers for LRC residuals | THM-002, HYP-2618, HYP-2619 | Compare F7 residuals with OCF activity/coimage sectors before inventing a new bucket. |
| LTI-116 | Noncommutative Redei-Berge recurrence | Ordered paths, edge deletion, contraction terms | Use noncommutative deletion-contraction as a packet-removal grammar with boundary debt | Redei/Berge notes, `dc_ocf_tracking.py`, HYP-2990 | Try a labelled packet deletion-contraction identity for Fejer/Haar handoff debt. |
| LTI-117 | GLMY path-homology witness | Directed path chains, Betti residues, boundary images | Treat path-homology classes as topological residuals when cycle/packet obstructions survive scalar tests | HYP-1815, HYP-2171, `beta2_zero_proof_v2.md` | Compute Betti/path-boundary features for cocycle-carrier SCCs from HYP-2992. |
| LTI-118 | Heap/tableau/sorting-network carrier | Heaps, staircase cells, hook/swap conflicts | Translate fixed-path tournament tiles into heap conflicts and sorting-network moves with visible local dependencies | `heap-tournament-tableau-triangle.md`, staircase reflections | Attach heap conflict labels to Haar staircase tilings and pair-cell operation grids. |
| LTI-119 | Transfer-matrix / Hopf convolution | Transfer states, subset convolutions, strong-component products | Replace brute-force enumeration with composable state algebras while keeping packet boundary states | `transfer-matrix.md`, HYP-1835, HYP-2228 | Build a transfer matrix whose states are labelled LRC packets, not scalar rows. |
| LTI-120 | Walsh puncture support gate | Low-degree Walsh support, punctures, Reed-Muller-like layers | Detect which low-degree cancellation gates a quotient preserves or destroys | HYP-2426, HYP-2459, HYP-2991 | Add Walsh support vectors to Haar zipper cocycle outputs. |
| LTI-121 | Krawtchouk association shadow | Orthogonal polynomial coordinates / association schemes | Use Krawtchouk coordinates as a warning system for basis shadows that anonymize packet identity | HYP-2723, LTI-075 | Test whether Krawtchouk shadows separate AP/GW, K33, C27, and F7 packets. |
| LTI-122 | Paley / QR circulant orbit carrier | Quadratic-residue orientations, circulant paths, character orbits | Separate finite-field character symmetry from real Farey/magnitude information | Paley path-homology notes, HYP-2227, QR reflections | Use QR/NQR orbit labels only with exact `M`, endpoint, and Ramanujan period fields attached. |
| LTI-123 | Terwilliger / DRT algebra carrier | Distance-regular layers, association modules, idempotents | Look for hidden module blocks behind packet residuals rather than one scalar obstruction | DRT/Satake notes, HYP-2990 | Try a module decomposition of F7 after endpoint/Fejer/Ramanujan coordinates are fixed. |
| LTI-124 | Matroid / gammoid obstruction test | Circuits, cocircuits, gammoid reachability | Decide when packet independence behaves like a matroid and when odd-cycle compatibility refutes it | gammoid tests, HYP-2986, HYP-2618 | Run matroid-circuit tests on zero-open residual packet systems. |
| LTI-125 | Hypergraph blocker allocation game | Blocker hyperedges, failed twists, factor-capture packets | Model failed certificates as a hypertournament/allocation game instead of a pair graph | HYP-2972, THM-566, HYP-2447 | Build blocker hypergraphs whose vertices are proof obligations and whose edges are simultaneous covering failures. |
| LTI-126 | SCC condensation / source-deleted deck | Source-deleted classes, condensation DAGs, observer roots | Lift A000568 source tests by retaining source-deleted deck and adjacent gap data | THM-381, THM-385, HYP-1981, HYP-1988 | For each LRC packet, store observer-source deck plus blocker indegree and adjacent gap labels. |
| LTI-127 | Operation-state vector atlas | Add/drop/double/product mutations, pair-cell states | Keep route history and second-order pair cells before comparing rows | HYP-1835, HYP-1976, HYP-2931 | Make a compact operation vector for HYP-2963 rows and test edge flips under each scalar. |
| LTI-128 | Quotientope / clique-complex carrier | Quotient cells, topes, clique complexes | Check whether a quotient has a cell complex whose faces preserve lonely/non-lonely predicates | HYP-2986, HYP-2947, HYP-2248 | Build the quotient cell complex for no-tope/no-cocircuit candidates. |
| LTI-129 | Zeta / Ihara / path-torsion carrier | Prime cycles, closed walks, torsion factors | Use zeta products as cycle-packet ledgers, not as naked analytic numerology | Ihara reflections, HYP-2896, THM-002 | Compare emitted cocycles with prime-cycle factors and record which factors are annihilated. |
| LTI-130 | Sequence-shadow companion bank | Fixed, merged, bisection, q-shadow, skinny companion sequences | Treat sequences as weather reports over carriers, useful only with the labels that generate them | HYP-2209, sequence-shadow lab, HYP-2618 | Add companion sequence shadows to packet families and note when shadows mix routes. |
| LTI-131 | Lucas / Pisano active-layer checksum | Modular recurrence layers, CRT periods, active clocks | Use recurrence periods as checksums for local p-adic or shell clocks, not final proof quotients | HYP-2437, HYP-2444, THM-486 | Attach Pisano/CRT layer labels to denominator-wall packets where period data survives. |
| LTI-132 | Legendre / ruler 2-adic amplitude ladder | Character signs, ruler valuations, Walsh amplitudes | Separate 2-adic parity ladders from real endpoint geometry | HYP-2459, HYP-2426, signed-LRC notes | Test whether Legendre/ruler signs predict which Haar `zeta` switches are color-compatible. |
| LTI-133 | Redei parity / higher 2-adic tower | Hamiltonian-path parity, mod `2^k` residues | Parity is a guaranteed floor; higher residues can expose obstruction layers but cannot replace packet labels | THM-001, THM-002, residue feedback loop | Audit LRC packet tournaments for mod-4/mod-8 path-count residue changes at route boundaries. |
| LTI-134 | Hamiltonian-path sheaf / Yoneda sections | Local source sections, restriction maps, global path sections | Phrase "observer can be lonely" as existence of compatible local sections over source decks | HYP-2161, HYP-2101, THM-381 | Build a small sheaf whose sections are source witnesses and whose failed gluings are F7 debts. |
| LTI-135 | Automorphism tax / orbit-size divisibility | Stabilizers, anti-cosets, transporter records | A quotient may forget labels only after stabilizer/anti-coset side channels are accounted for | THM-409, HYP-2208, HYP-2228 | Add `(Aut,Anti,orbit-size)` records to apex, self-converse, and signed LRC quotient tables. |
| LTI-136 | Product-rule / irreducibility no-lift | Coefficient convolutions, hidden factor grids, boundary totals | Reducibility and LRC quotienting share the same danger: boundary totals may hide a 2D lift | HYP-2452, HYP-2455, HYP-2447, HYP-2993 | Port convolution no-lift infeasibility to Haar/Fejer packet grid compression. |
| LTI-137 | ASM / tableau / alternating-sign bridge | ASMs, monotone triangles, local switch signs | Use alternating-sign local conservation as a model for fixed-margin Haar switch ledgers | ASM reflections, HYP-2991 | Test whether packet `zeta` ledgers have ASM-like local conservation laws. |
| LTI-138 | Visibility taut-path proof search | Intervals, line-of-sight, taut obstacles, lazy checks | Turn any-angle path planning into a finite proof search over visible safe components and certified obstacles | HYP-2948, HYP-2949, LTI-054 | Implement a lazy verifier that expands only boundary atoms not killed by Fejer/Haar labels. |
| LTI-139 | Relation-lattice / subtorus packet complex | Relation lattices, subtori, height strata | Split low-height exact checks from high-height Weyl/decay while preserving packet labels | HYP-2540, HYP-2599, T826 | Tensor relation-height with exact-period and endpoint-owner fields in true-wide packets. |
| LTI-140 | Faulhaber odd-moment bridge | Odd moments, midpoint anchors, boundary-moment charts | The odd-moment survivor is the moment analogue of odd-cycle/tournament residuals | HYP-2457, HYP-2459, HYP-2530 | Compare danger-count moment duals with Faulhaber odd-moment Hankel packets. |
| LTI-141 | A000568 iso-class deck lift | Unlabelled classes plus marks/decks/magnitude fields | A000568 is useful only after restoring marks, thresholds, source, and magnitude/Farey data | HYP-1977, HYP-1978, HYP-2924 | Create an "allowed-forgetting" checklist for every unlabelled tournament class used in LRC. |
| LTI-142 | Packet-residual metagraph Laplacian | Residual packets as nodes, handoff edges, Laplacian eigenvectors | Measure bottleneck clusters and second axes of failure in the proof metagraph | metagraph reflections, HYP-2987, HYP-2990 | Compute SCCs, cuts, and low eigenvectors of the F0-F7 packet metagraph. |
| LTI-143 | Winding-tournament phase entropy | Phase order movies, loneliness masks, entropy over `H` | Use winding tournaments to locate structural bridges while admitting magnitude blindness | HYP-2605, HYP-2920, HYP-2924 | Pair every winding class with Farey-neighbor scale data before making tightness claims. |
| LTI-144 | Binding-pair switch tournament | Pair crossings, clearing switches, active blocker decks | THM-524 says the hard gap problem lives at binding-pair switches, but only after other-runner clearance is retained | THM-524, HYP-3019, HYP-2571, HYP-2621 | Add binding-switch fields to HYP-2963 packet records and test family fibers for mixed clearance routes. |
| LTI-145 | Coimage wall-address atlas | Projective coimage classes, signed wall mass, low-height relations | Route signed residual tails by class before taking absolute bounds | HYP-2617, HYP-2618, HYP-2619, HYP-2626 | Attach coimage class IDs to Ramanujan/exact-period and Haar residual packets. |
| LTI-146 | Namespace collision audit | Hypothesis/tangent IDs, filenames, anchors, claimed packets | Prevent research pointers from identifying the wrong quotient or theorem after concurrent pushes | CONCURRENT-SESSIONS.md, duplicate HYP-2992 signal | Add a small audit script/report for duplicate frontmatter IDs and reused tangent numbers. |
| LTI-147 | Cocycle normal form | Cocycle channels / proof obligations in the LRC packet complex | Every forgotten coordinate must descend, become a coboundary, be dual-annihilated, or route to a named residual class | HYP-2997, HYP-2996, HYP-2995, HYP-2994, HYP-2992, HYP-2991, HYP-2990, THM-572 | Build HYP-2963 packet-level cocycle ledgers and tag every low-frontier row by first nonzero class. |
| LTI-148 | Residual-section packet grid | Residual sections and Haar-product exits | Make F7 executable by requiring hard non-AP/GW packets to have owner-strip, cross-handoff, or nested-refinement exits | HYP-2996, HYP-2963, HYP-2989, HYP-2995 | Lift the bounded section map into family templates and Fejer/Ramanujan certificate manifests. |
| LTI-149 | Pascal-slope additive-basis Farey packet schema | Representation hypergraph rows, additive fibers, Farey operator lanes | Complement HYP-2998 by keeping Pascal-row, carry-width, and Farey `p+q,p*q,q^p,p^q` lane fields as labelled packet data | HYP-2999, HYP-2998, HYP-2940, HYP-2934, HYP-2931, HYP-2523, S501 | Add `additive_basis_regime`, `representation_entropy`, `local_residue_rank`, `carry_width`, `pascal_slope_row_id`, and `farey_operator_lane` fields to packet records. |
| LTI-150 | Fibonacci path-rank additive-basis bridge | Path independent-set ranks, Zeckendorf carry graph, additive-basis hypergraphs, Farey unit-excess payloads | Reads `F_n=sum_k binom(n-k-1,k)` as `I(P_{n-2};1)`, separating smoothing, bounded-arity, normal-form, additive, product, and power-stress proof currencies | HYP-3000, HYP-2999, HYP-2998, HYP-1902, HYP-2984, HYP-2934, HYP-2982 | Classify HYP-2963 packet debts as smoothing, bounded-invoice, normal-form carry, additive Farey scale, product incidence, or stress-only magnitude. |
| LTI-151 | Labelled smoothing switchboard | Packet routes plus admissible kernels | Choose Fejer/Ramanujan/boundary-moment/large-sieve/Kaczynski only after the packet route is known | HYP-2984, HYP-2981, HYP-2979, HYP-2982, HYP-2983 | Scale the `16`-packet route matrix to HYP-2963 families and count any unrouted packet class. |
| LTI-152 | Curried packet functional tower | Higher-order packet evaluators, quotient derivatives, partial applications | Treat each LRC route as `S -> packet -> root -> lane -> fiber -> certificate -> verdict`; a quotient may forget a coordinate only after its lost-coordinate function is discharged | HYP-3002, HYP-3000, HYP-2999, HYP-2997, HYP-2996, HYP-2995, HYP-2974, HYP-2963, HYP-2954 | Add `curried_call_signature` and `lost_coordinate_function` fields to packet records and certificate manifests. |
| LTI-153 | Summand/multiplicand Farey basis merge | Farey packets with summand antidiagonal and multiplicand hyperbola fibers | Keep `p+q` as additive/pinch fiber and `p*q` as factor/Kpq incidence fiber before comparing Goldbach, polygonal, Zeckendorf, or power shadows | HYP-3003, HYP-3000, HYP-2999, HYP-2998, HYP-2935, HYP-2934, HYP-2083, HYP-2067 | Add `summand_fiber_id`, `multiplicand_factor_fiber`, `zeckendorf_carry_width`, and `farey_shadow_lane` to sequence-shadow packet records. |
| LTI-154 | Dichotomy recursion mode atlas | Binary/operator mode switches as proof carriers | Treat odd/even, sign, addition/multiplication, `+1`/`/2`, `*2`/`+2`, Farey sum/product, Zeckendorf, Collatz, and triune recursion as allowed-forgetting switchboards; S171's summand/multiplicand LTI-153 is the add/multiply subcase | HYP-3004, HYP-3003, HYP-3002, HYP-3001, HYP-3000, HYP-2999, HYP-2998, HYP-2262, HYP-2238, HYP-2134, HYP-2091 | Add `parity_block`, `sign_cut_status`, `additive_pair_sum_owner`, `multiplicative_unit_orbit`, `recursion_boundary_state`, and `smoothing_route` to HYP-2963 packet records. |
| LTI-155 | Representation-economy carrier | Goldbach smoothing, polygonal arity, Zeckendorf carries, Farey address clocks | A sequence shadow is usable only after declaring the representation economy whose fibers preserve the LRC predicate | HYP-2998, HYP-3004, HYP-3003, HYP-3000, HYP-2999, HYP-2523, HYP-2416, HYP-2931, HYP-1920 | Add `representation_economy` fields to packet classifiers: smoothing, bounded arity, normal form, or Farey address. |
| LTI-156 | Technique multiverse annex | Technique families as vertices; retained LRC payload dimensions as pairwise observable | Keeps the broad tournament/metagraph/series/sieve/cocycle toolbox as a cyclic contribution surface rather than a scalar ranking | HYP-3005, T1090, `00-navigation/LRC-TECHNIQUE-MULTIVERSE-INDEX.md` | Promote theorem-facing `LTM-*` cards into `LTI-*`/`LTT-*` entries and run packet-level fiber-mixing tests before trusting scalar quotients. |
| LTI-157 | Poincare worldline frame ledger | Time/phase cylinder worldlines, observer velocity, danger tube metric, sign kernel, boost cocycle | Separates true anchored-LRC automorphisms from observer-coupled boosts and Lorentz-like deformed-tube shadows | HYP-3007, HYP-3006, HYP-3002, HYP-2997, HYP-2963, HYP-2486, HYP-2291, THM-381, THM-385 | Add `observer_velocity_label`, `relative_speed_normal_form`, `sign_kernel_status`, `primitive_scale_gcd`, `tube_metric_label`, `worldline_frame_label`, `boost_cocycle_status`, and `orientation_debt_for_winding` to packet records. |
| LTI-158 | Automatic gap-language carrier | Moser-de Bruijn DFA states, fibbinary no-adjacent carry states, product phase automata, valuation/SML gates, exact safe-component side channels | Sequence shadows are safe only after retaining `2`-adic phase, carry-boundary, exact `M=p/q`, endpoint/safe-component geometry, or eventual-periodic recurrence coordinates; Moser is stable under `n -> 4n`, fibbinary under `n -> 2n` | HYP-3008, HYP-3010, HYP-3004, HYP-3003, HYP-3001, HYP-3000, LTI-132, LTI-150, LTI-154 | Add `automatic_gap_carrier` to packet records: `none`, `fibbinary_2adic_normal_form`, `moser_even_phase_gap`, `product_phase_automaton`, `sml_eventual_periodic_gate`, or `valuation_gate_only`; pair it with exact `M`, Farey excess, residue owner, and safe-component fields before using sequence membership as evidence. |
| LTI-159 | Fermat-Catalan automatic-gap / power-lift extension | Moser-de Bruijn even-bit words, fibbinary/Zeckendorf no-adjacent carries, lacunary gap support, Fermat-Catalan perfect-power guards, Hurwitz doubling states, exact safe-component side channels | Extends LTI-158 by keeping sequence/gap/power shadows as packet fields before using them to split LRC14 residual routes; HYP-3017 has now run the full HYP-2963 route-purity audit and found `143` mixed-route automatic-word fibers, downstream of HYP-3016's single-swap magnitude-cocycle guardrail | HYP-3017, HYP-3016, HYP-3015, HYP-3014, HYP-3009, HYP-3010, HYP-3008, HYP-3007, HYP-3003, HYP-3000, HYP-2998, HYP-2963, HYP-2950, HYP-2944, HYP-2937, HYP-2702, HYP-2698, HYP-1902 | Split the largest mixed automatic fibers, especially `MFCMMCCFFFCCC`, by exact `M/q`, magnitude cocycle, endpoint geometry, persistence barcode, C27/K33, Haar/Fejer, covering, valuation, and residual labels before using any automatic word as a theorem coordinate. |
| LTI-160 | Automatic lacunary safe-component filter | DFA states, 2-adic windows, gap-block labels, exact safe components, exact-period packet fields | Use fibbinary/Moser/2-adic window automata as packet tags before scalarization; HYP-3011 shows first-13 automatic rows are positive-open while AP/GW stay zero-open boundary atoms | HYP-3011, HYP-3008, HYP-3002, HYP-2997, HYP-2963, HYP-1902, THM-572 | Add safe-component readouts to the HYP-3008 `automatic_gap_carrier` field before using automatic membership as LRC evidence. |
| LTI-161 | Gap automaton induced-class carrier | Fibbinary/Moser automata, lacunary support, 2-adic doubling, Fermat-Catalan exponent budgets, visibility cores, induced tournament classes | Extends HYP-3011's automatic lacunary safe-component filter and HYP-3010's exact maximin/safe-component audit by adding induced tournament isomorphism-class ledgers; sequence/gap analogies are theorem-safe only after retaining automaton state, native transition, gap-boundary label, valuation budget, exact packet route, safe-component side channel, and induced class census | HYP-3012, HYP-3011, HYP-3010, HYP-3009, HYP-3008, HYP-3007, HYP-3006, HYP-2998, HYP-2997, HYP-2983, HYP-2982, HYP-2963, THM-572 | Add `automaton_language_id`, `automaton_state_word`, `gap_support_ratio_label`, `hadamard_boundary_warning`, `doubling_transition_state`, `base4_digit_mask`, `zeckendorf_carry_state`, `valuation_exponent_budget`, `finite_exception_budget`, `visibility_core_label`, `safe_component_label`, and `induced_tournament_class_word` to packet records or sidecar manifests. |
| LTI-162 | Perfect-number divisor packet merge | Euclid-Euler controls, LRC14 `q=14a-1` shadows, abundancy defect, divisor factorization, Kpq product incidence, automaton state | Imports the perfect-number carrier as a typed side channel: the `n=2` chain is exact, prime `q=14a-1` rows are deficient by `12/q`, and composite `q14` rows can flip abundant, so product analogies require exact `M`, prime/composite `q`, factorization, defect, Kpq route, and automaton labels | HYP-3013, HYP-3012, HYP-3009, HYP-3008, HYP-2946, HYP-2945, HYP-2941, HYP-2221, HYP-2220, HYP-2963, THM-572 | Add `unit_excess_apex`, `perfect_control_status`, `abundancy_defect`, `divisor_lattice_factorization`, `prime_q_flag`, `product_incidence_rank`, and `automaton_transition_state` to unit-excess/product sidecars before using perfect-product analogies. |
| LTI-163 | Creative exact packet-lens atlas | Cech nerves, tropical slack, CRT solenoid charts, endpoint currents, danger-count distributions, tope/cocircuit walls, automaton-divisor sidecars | Turns speculative LRC angles into exact side-channel fields; AP/GW are zero-open boundary controls, K33/petal rows are positive-open despite denominator-14 boundary witnesses, and covering rows can require a first CRT chart such as `17/41` | HYP-3014, HYP-3013, HYP-3012, HYP-3008, HYP-2974, HYP-2973, HYP-2970, HYP-2969, HYP-2965, HYP-2963, HYP-2949, HYP-2948, THM-572 | Add `cech_nerve_class`, `positive_component_count`, `tropical_slack_margin`, `crt_solenoid_first_chart`, `endpoint_current_word`, `danger_count_distribution`, `tope_cocircuit_wall_state`, and `automaton_divisor_sidecar` to named-row and covering-family packet sidecars. |
| LTI-164 | Lonely-profile persistence barcode | Superlevel bars of `m_S(t)=min_v ||v t||`, exact peak heights, bar lengths, stability margins | Promotes exact safe components into barcode sidecars: AP/GW are zero-bar atoms, while positive rows retain bar count, peak height, length, and persistence before scalarization | HYP-3015, HYP-3014, HYP-3013, HYP-3012, HYP-3011, HYP-3010, HYP-2996, HYP-2975, HYP-2974, HYP-2963, THM-572 | Add barcode sidecar fields to HYP-2963 packets and route lowest-persistence families to endpoint-owner / Fejer interval certificates. |
| LTI-165 | Automaton fiber-mixing stress | Residue-language fibers, Moser/fibbinary terminal states, exact safe components, magnitude/Farey cocycles | Tests whether finite-state sequence quotients preserve the LRC14 boundary-vs-open predicate; HYP-3016 shows residue+automaton terminal fields mix AP/GW boundary atoms with strictly open rows, so automaton state is only admissible with a retained or discharged magnitude cocycle | HYP-3016, HYP-3015, HYP-3014, HYP-3013, HYP-3012, HYP-3011, HYP-3009, HYP-3008, HYP-3002, HYP-2997, HYP-2963, HYP-2928, THM-572 | Add `residue_automaton_fiber_id`, `automaton_terminal_word`, `magnitude_cocycle`, `farey_magnitude_height`, `fiber_anchor_row`, `safe_component_measure`, `safe_component_count`, `largest_safe_component`, and `fiber_mixing_exit` to automaton packet sidecars. |
| LTI-166 | Active-bottleneck normal fan | Bar endpoint owner sets, peak bottleneck owner sets, support speeds, residue sums, local equation counts | Strengthens barcode and automaton-fiber carriers by keeping the local active inequalities: AP/GW are zero-bar boundary owner-pair atoms, while open rows in mixed automaton fibers emit explicit peak supports such as `(5,7)`, `(5,96)`, `(5,36)`, `(7,20)`, `(1,26)`, or `(13,84)` | HYP-3018, HYP-3017, HYP-3016, HYP-3015, HYP-3014, HYP-3013, HYP-3012, HYP-3011, HYP-3010, HYP-3009, HYP-3008, HYP-2996, HYP-2963, THM-572 | Add normal-fan sidecars to HYP-2963 packets and prove mixed automatic fibers become pure after retaining peak support and endpoint residue sums, or route the remaining support classes to Fejer/state-lift certificates. |
| LTI-167 | Discrepancy-height trident carrier | Residue discrepancy, Erdos-Turan proxy bins, Mahler/Farey height, Hensel singular-lift counts, bounded feature incidence | Splits the bounded AP/GW automaton-fiber leaks by restoring three independent proof clocks: noisy denominator, lost magnitude, and local p-adic lift status; HYP-3020 has `0` mixed boundary/open fibers for the full trident on the named plus single-swap bank, but warns that the signature is nearly exact and must be compressed | HYP-3020, HYP-3017, HYP-3016, HYP-3015, HYP-3014, HYP-3009, HYP-3008, HYP-2997, HYP-2995, HYP-2991, HYP-2989, HYP-2963, THM-572 | Add trident sidecars to HYP-2963 packets, then coarsen them until the full bank keeps `mixed_status=0` while routing nonzero coordinates to magnitude-cocycle, Ramanujan/Haar/Fejer, Hensel, or K33/F7/THM-572 exits. |
| LTI-168 | Pair-good decoy generator classifier | Source pair lanes, blocker residue teeth, active blocker decks | Replaces raw decoy counts by the exact modular tooth inequality `14*min(c*p mod q,q-c*p mod q)<q`; large counts mostly mean many source crossings hit the same AP-core tooth cover | HYP-3021, HYP-3019, HYP-3018, THM-524, HYP-2990, HYP-2963 | Add decoy generator sidecars to HYP-2963 packet records and test whether each packet family has a bounded generator-cover template. |
| LTI-169 | Pair-good decoy barcode/normal-fan refinement | Pair-good false switches, sum/diff lanes, pair shells, denominator lanes, pair gaps, blocker decks, barcode/normal-fan relations | Extends LTI-168 by recording whether each generated false switch is outside a strict barcode bar, in an AP/GW zero-bar row, or overlaps/disjoins HYP-3018 active peak supports through the good pair or blocker deck | HYP-3022, HYP-3021, HYP-3019, HYP-3018, HYP-3015, HYP-2963, THM-524, THM-572 | Add pair-good barcode/normal-fan sidecars to HYP-2963 packets; prove blocker-deck lemmas for common singleton blockers and route overlap classes by normal-fan support. |
| LTI-170 | Automatic fiber zipper splitter | Automatic-word fibers, residue-terminal automaton fibers, exact magnitude cocycles, lonely-profile barcode shadows, packet zippers, route-labelled packets | Converts the HYP-3017/HYP-3016 mixed-fiber warning into a proof ladder; HYP-3023 shows automatic words and residue-terminal fields still mix theorem routes on the full HYP-2963 bank, while exact magnitude cocycle is the first tested non-route coordinate with zero mixed theorem-route fibers | HYP-3023, HYP-3022, HYP-3021, HYP-3020, HYP-3019, HYP-3018, HYP-3017, HYP-3016, HYP-3015, HYP-3014, HYP-3009, HYP-3008, HYP-2963, THM-572 | Prove the family magnitude-cocycle lemma inside automatic/residue fibers, starting with `MFCMMCCFFFCCC`; use barcode and packet zippers as certificate anchors. |
| LTI-171 | Fiber-zipper convergence / Henselian unit rule | Residue-terminal fibers, coarse Erdos-Turan clocks, p-adic unit roots, zero-root scale debt, magnitude/barcode zipper fields | Extends HYP-3023 and HYP-3020 on the full HYP-2963 bank: exact ET clocks are too sharp, but coarse ET+unit data has `0` mixed boundary/open fibers with only `15` mixed open-route fibers and max mixed `4` | HYP-3024, HYP-3023, HYP-3020, HYP-3017, HYP-3016, HYP-3015, HYP-2963, THM-572 | Prove coarse status convergence before route purity: every automatic/residue fiber either has bounded ET discrepancy, Hensel unit/zero-root debt, family magnitude formula, or named K33/F7/THM-572 residual. |
| LTI-172 | Closed arc-Cech nerve carrier | Individual threshold danger arcs, endpoint-completed closed Cech nerves, open tope nerves, boundary cocircuit facets, runner quotient defects | Makes the actual circle-cover topology the proof carrier: AP/GW are closed arc full-cover cycles with open pieces glued only by endpoint cocircuits, while K33, petals, covering, fibbinary, and Moser controls have closed arc `beta1=0` and positive safe mass | HYP-3025, HYP-3024, HYP-3023, HYP-3022, HYP-3021, HYP-3020, HYP-3018, HYP-3016, HYP-3015, HYP-3014, HYP-3013, HYP-3012, HYP-3011, HYP-3010, HYP-3009, HYP-3008, HYP-2997, HYP-2986, HYP-2975, HYP-2974, HYP-2970, HYP-2963, HYP-2990, THM-572 | Add `closed_arc_cech_beta`, `open_arc_component_count`, `boundary_cocircuit_facet_word`, `boundary_owner_sum_word_mod_14`, `runner_quotient_betti_defect`, `private_arc_count`, `private_runner_count`, `safe_tope_count`, and `arc_cech_exit_route` to HYP-2963 packets or a sidecar before trusting runner-level nerves. |
| LTI-173 | Carrier-fusion switchboard | Fusion signatures, safe-stick/potato bodies, barcode shapes, CRT charts, endpoint currents, magnitude cocycles, automatic/doubling and ET/Henselian sidecars, danger-count duals | Combines the HYP-3014..HYP-3025 sidecars and tests quotient leakage: automatic words and chart-denominator sidecars still mix AP/GW boundary atoms with open rows, while barcode shape, magnitude cocycle, and the full fusion signature remove boundary/open mixing in the named bank; HYP-3024's coarse ET/Henselian gate is the cheapest known full-bank boundary/open status carrier inside the fusion packet | HYP-3026, HYP-3025, HYP-3024, HYP-3023, HYP-3022, HYP-3021, HYP-3020, HYP-3018, HYP-3017, HYP-3016, HYP-3015, HYP-3014, HYP-3013, HYP-3009, HYP-2963, HYP-2974, HYP-2969, THM-572 | Add `fusion_signature`, `largest_safe_stick`, `safe_body_mass`, `barcode_shape`, `magnitude_cocycle`, `endpoint_current_word`, `crt_first_chart`, `danger_distribution_word`, `et_henselian_unit_zipper`, and `doubling_transition_word` to HYP-2963 packet sidecars, then rerun route-purity on automatic fibers. |
| LTI-174 | Side-channel repair ladder | Automatic-word fibers plus exact M/q, boundary topology, packet labels, tail magnitude, trident/decoy sidecars, fiber-convergence data, arc-Cech topology, automatic zipper data, and guarded non-route signatures | Turns failed quotients into an ordered repair theorem after HYP-3023/HYP-3024/HYP-3025: exact `M` repairs open/boundary status but not theorem route, `M+q` and boundary topology leave one mixed route pair each, HYP-3020/HYP-3021/HYP-3022 supply adjacent side-channel clocks, and C27/K33/transfer labels or the guarded non-route signature are route-pure on the HYP-2963 bank | HYP-3027, HYP-3026, HYP-3025, HYP-3024, HYP-3023, HYP-3022, HYP-3021, HYP-3020, HYP-3018, HYP-3017, HYP-3016, HYP-3015, HYP-3014, HYP-2963, HYP-2997, HYP-2995, HYP-2992, THM-572 | Add `repair_ladder_stage`, `first_nonzero_repair_cochain`, `route_purity_status`, `status_purity_status`, and `guarded_nonroute_signature` to automatic sidecar manifests; prove local zippers for the two remaining mixed pairs against the HYP-3023 magnitude-cocycle splitter, HYP-3024 convergence audit, and HYP-3025 arc-Cech topology. |
| LTI-175 | Carrier pullback mega-index | Cross-surface pullbacks from tournament, metagraph, series, automaton, topology, harmonic, arithmetic, geometry, computation, and formalization techniques | Converts the project-wide technique archive into `90` operational LRC carrier rows, each naming source handles, retained packet fields, quotient guardrail, and a next-agent pull; vertices are proof-carrier pullbacks, not runners | T1108, LTT-073, `00-navigation/LRC-CARRIER-PULLBACK-INDEX.md`, LTI-173, LTI-156, LTM-* | Pick a `CPI-*` row, instantiate its retained fields on HYP-2963 packets or a named stress family, run status/route fiber mixing, then promote successful rows into theorem-facing LTI/LTT cards. |
| LTI-176 | Residual status-gate switchboard | Residual proof obligations after coarse ET+Henselian status convergence | Reframes the `15` remaining full-bank coarse-gate mixed-route fibers as certificate scheduling debt because the gate has `0` mixed boundary/open fibers; the first theorem should prove status convergence, not full route purity | HYP-3028, HYP-3026, HYP-3024, HYP-3023, HYP-3020, HYP-2963, THM-572 | Add a cached HYP-2963 packet-ledger mode, list the 15 residual fibers, and attach the first successful tooth: q-witness, safe-stick/barcode/Fejer/Haar, unit-petal, K33/F7/THM-572, covering, or magnitude-cocycle formula. |
| LTI-177 | Safe-component stalk descent | Largest strict safe component stalks, endpoint owner residues, peak bottleneck residues, local length/height, open/boundary status | Reverses HYP-3026 fusion: on target word `MFCMMCCFFFCCC`, residue-terminal fibers have `27` mixed routes, owner-only stalks reduce to `7`, coarse stalks to `2` size-2 open-route residuals, and exact largest-component stalks to `0`, matching exact magnitude while retaining local topology | HYP-3029, HYP-3026, HYP-3025, HYP-3024, HYP-3023, HYP-3018, HYP-3015, HYP-2963, THM-572 | Prove largest-stalk descent inside the target automatic fiber, discharge the two coarse residual pairs `13->159/117` and `13->118/104`, then stress exact stalk keys on the full HYP-2963 bank against Cech and normal-fan sidecars. |
| LTI-178 | Status-topology gate bridge | Arc-boundary cycle, coarse ET+unit status gate, magnitude route splitter, barcode packet scheduler, residue-terminal word | Orders HYP-3025 topology before HYP-3024 arithmetic inside the HYP-3026/HYP-3028/HYP-3029 switchboard stack: residue-terminal fibers have exactly `2` boundary/open collisions, AP and GW are the only boundary rows and carry the arc-Cech `(1,1)` cycle, while every open cohabitant and every coarse route-mixed survivor has closed arc `beta1=0` | HYP-3030, HYP-3029, HYP-3028, HYP-3027, HYP-3026, HYP-3025, HYP-3024, HYP-3023, HYP-3020, HYP-3018, HYP-3016, HYP-3015, HYP-2963, THM-572 | Prove zero-open implies AP/GW arc-boundary cycle or named F7/THM-572 residual, then allow coarse ET+unit to forget route labels while preserving boundary/open status. |
| LTI-179 | Haar-tile repair ladder dictionary | Mixed Haar cocycles, fixed-margin switches, fixed-path tournament tiles, automatic-fiber repair cochains | Identifies HYP-2989/HYP-2991/HYP-2992 and HYP-3023..HYP-3030 as one product-rule guardrail: row/column shadows, automatic words, raw component counts, and scalar tournament classes are unsafe unless the mixed coordinate `zeta=T00-T01-T10+T11` is reconstructed, annihilated, descended, boundary-stopped, or routed to F7/THM-572 | HYP-3031, HYP-3030, HYP-3029, HYP-3028, HYP-3027, HYP-3026, HYP-3025, HYP-3024, HYP-3023, HYP-3020, HYP-2992, HYP-2991, HYP-2989, HYP-2997, HYP-2995, HYP-2990, HYP-2963, THM-572 | Add `zeta_repair_class` to packet manifests with values `orthogonal_zero`, `same_tile_boundary`, `owner_strip`, `cross_handoff`, `nested_refinement`, `residual`; build two-coordinate packet grids for the two HYP-3027 residual mixed pairs and use HYP-3030/HYP-3028/HYP-3029 to schedule open-route residual teeth. |
| LTI-180 | Analytic sieve-clock bridge | Mobius tails, squarefree inverse-unit capacity, large-sieve minor-arc budgets, circle-method splits, exponential-sum checksums, smoothing defects, Kaczynski approach classes, and labelled repair packets | Converts HYP-2982/HYP-2983 analytic-number-theory imports into HYP-3027 repair-ladder clocks. HYP-3032/S196 shows `mu^2/phi` is a capacity meter with a blindness certificate: C27 prime-power petals and fibbinary `q=25` vanish under it, while the HYP-3027 `q=23` petal/covering pair remains mixed even after exact denominator and non-route analytic packet fields. | HYP-3032, HYP-3031, HYP-3027, HYP-3026, HYP-3024, HYP-3023, HYP-3020, HYP-2985, HYP-2984, HYP-2983, HYP-2982, HYP-2979, HYP-2978, HYP-2963, HYP-2997, HYP-2995, HYP-2992, THM-572, T1113, LTT-078 | Add `mobius_tail_clock`, `mu2_phi_capacity`, `squarefree_blindness_report`, `large_sieve_budget`, `exponential_sum_checksum`, `smoothing_defect`, and `kaczynski_approach_class` to repair-ladder manifests; use the `q=23` petal/covering pair as the first analytic-clock zipper test. |
| LTI-181 | Residual certificate teeth | S194 residual coarse ET+unit route-mixed ledger, arc topology compact signature, safe-tope/quotient-defect topology bucket, unit-scale tooth | Turns HYP-3030's `15` open route-mixed coarse fibers into a small route scheduler: topology alone leaves `3` mixed classes, unit-scale alone leaves one mixed class, exact `M` alone leaves `2` mixed classes, but full topology or topology bucket plus unit-scale tooth yields `21` residual fibers with `0` route mixing | HYP-3033, HYP-3031, HYP-3030, HYP-3029, HYP-3028, HYP-3026, HYP-3024, HYP-2963, THM-572 | Promote `residual_topology_bucket`, `unit_scale_tooth`, and `residual_certificate_tooth` into packet sidecars, rerun on the full bank instead of stored S194 text, then prove the route scheduler familywise and connect teeth to HYP-3031 zeta repair classes. |
| LTI-182 | Arc-boundary path-lift carrier | Closed danger-arc Cech boundary ranks, GF(2) H1 representative, owner-deletion persistence | Pulls older path-homology/deletion machinery back to the HYP-3030 topology gate without using runner tournaments.  On the HYP-2963 residue-terminal status-collision surface, AP and GW are the only closed-H1 rows; each has a `58`-edge representative and deleting any owner speed kills H1, while all open cohabitants have closed H1 `0` | HYP-3034, HYP-3031, HYP-3030, HYP-3029, HYP-3028, HYP-3025, HYP-3024, HYP-3018, HYP-2963, THM-572 | Add `closed_arc_boundary_rank_d1`, `closed_arc_boundary_rank_d2`, `closed_arc_h1_rep_edge_count`, `closed_arc_h1_owner_support`, and `owner_deletion_beta1_word` sidecars; prove zero-open implies owner-essential AP/GW path-cycle or named residual debt. |
| LTI-183 | Residual tooth atlas | Coarse ET+unit residual fibers, arc-topology compacts, safe-component stalks, magnitude cocycles, q/covering certificate labels | Turns the `15` coarse ET+unit route-mixed open fibers into a first-tooth manifest layered after HYP-3034's topology-front path lift and HYP-3033's certificate scheduler: arc topology separates `13`, coarse safe-component stalk separates `2`, and all first repairs are HYP-3031 `owner_strip` teeth; exact stalk, magnitude, and q/covering labels are nested backups rather than first proof objects | HYP-3035, HYP-3034, HYP-3033, HYP-3032, HYP-3031, HYP-3030, HYP-3029, HYP-3028, HYP-3027, HYP-3024, HYP-3023, HYP-2963, THM-572 | Add `residual_tooth_class` and `first_tooth` to HYP-2963 sidecars, then prove arc-topology and coarse-stalk owner-strip descent lemmas instead of re-counting residual fibers. |
| LTI-184 | Ramanujan primitive-period route scheduler | Primitive denominator phases `a/q`, safe-residue count decks, first primitive safe denominator, Ramanujan trace profiles, coarse ET+unit residual fibers | Turns HYP-3030's `15` strict-open coarse route-mixed residual fibers into route-pure scheduler fibers without exact `M`: on the `38` S194 residual packets, `coarse+first_primitive_safe_q_2_13` and `coarse+primitive_safe_deck_2_13` both have `0` mixed route and `0` mixed status; all residual covering rows have zero `q<=13` primitive safe mass | HYP-3036, HYP-3033, HYP-3032, HYP-3031, HYP-3030, HYP-3028, HYP-3027, HYP-3024, HYP-3023, HYP-2963, LTT-082, LTI-181, LTI-180, LTI-179, LTI-178, LTI-176, CPI-043, OPEN-Q-108 | Add `primitive_safe_deck_2_13` and `first_primitive_safe_q_2_13` to packet sidecars, run a cached full-bank ledger, and prove that post-status Q-witness residuals are exactly the positive `q<=13` primitive-deck rows while zero-deck rows route to covering/q=14/boundary-moment certificates. |
| LTI-185 | Residual capacitor flow cuts | Mixed strict-open route pairs after status protection, first-cut stages, exact scale, topology, stalk, fusion, primitive-period, first-tooth, and packet-label sidecars | Treats the two HYP-3027 residual mixed-route pairs as capacitors in a proof-flow network downstream of HYP-3036's primitive-period scheduler, HYP-3035's residual tooth atlas, HYP-3034's arc-boundary path lift, HYP-3033's topology-bucket-plus-unit-scale scheduler, and HYP-3032's analytic clock bridge: the petal/covering exact-scale collision is cut by boundary topology (`nested_refinement`), while the K33/covering topology collision is cut by exact `M+q` (`cross_handoff`); later closed arc, stalk, fusion, primitive-period, first-tooth, and packet-label carriers split both | HYP-3037, HYP-3036, HYP-3035, HYP-3034, HYP-3033, HYP-3032, HYP-3031, HYP-3030, HYP-3029, HYP-3028, HYP-3027, HYP-3026, HYP-3024, HYP-3023, HYP-2992, HYP-2991, HYP-2990, HYP-2963, THM-572 | Compare the two capacitor min-cuts with HYP-3036's primitive-period split, HYP-3035's first-tooth owner-strip atlas, HYP-3034's path-lift topology, and HYP-3033's full residual certificate teeth, then add `residual_capacitor_id`, `first_cut_stage`, and `zeta_exit_class` to a direct HYP-2963 packet ledger. |
| LTI-186 | q=23 drop/add Haar square | Two-coordinate drop-pair/add-pair fixed-margin squares, exact-M zeta, safe-body geometry, endpoint-owner strips | Turns the HYP-3032 q=23 petal/covering residual into a HYP-3031 two-coordinate grid downstream of HYP-3037/HYP-3036: diagonal doubling keeps `M=2/23`, cross-swaps open as q-witness rows, exact-M zeta is nonzero, and endpoint-owner strips split the diagonal routes | HYP-3038, HYP-3037, HYP-3036, HYP-3035, HYP-3032, HYP-3031, HYP-3027, HYP-3026, HYP-2991, HYP-2989, HYP-2963, THM-572, T1119, LTT-084, LTI-185, LTI-184 | Add `drop_add_square_id`, `diagonal_doubling_match`, `exact_M_zeta`, `safe_body_zeta`, `endpoint_owner_strip`, and `zeta_repair_class` to residual-pair manifests; test whether more double-pair squares either open off diagonal, descend through a q-diagonal, or expose owner-strip routing data. |
| LTI-187 | Hidden-coordinate ledger | Controlled-forgetting stages, visible/hidden relation types, primitive decks, capacitor IDs, zeta classes, endpoint-owner strips, anti-wedge debt | Reads HYP-3024..HYP-3038 as a proof-order ladder rather than a scalar chain: every quotient is legal only after the next hidden coordinate is exposed, dual-annihilated, cut, or routed to named debt.  Connects address-retention, visible/hidden folds, anti-wedge transitivity, cochannel transfer, and pair-good generator classification to the current residual stack | HYP-3039, HYP-3038, HYP-3037, HYP-3036, HYP-3035, HYP-3034, HYP-3033, HYP-3032, HYP-3031, HYP-3030, HYP-3028, HYP-3027, HYP-3024, HYP-3023, HYP-3022, HYP-3021, HYP-3018, HYP-2963, THM-572, T1120, LTT-085 | Add `hidden_coordinate_stage`, `visible_hidden_relation_type`, `residual_capacitor_id`, `first_cut_stage`, `drop_add_square_id`, `endpoint_owner_strip`, and `anti_wedge_debt_count` to a cached packet ledger; reject scalar quotients that skip this ledger. |
| LTI-188 | Hidden statement ledger | Micro-statements / proof obligations connecting owner-essential H1, first-tooth owner strips, primitive period decks, analytic blindness, residual capacitors, q=23 zeta squares, automaton shadows, and decoy generator teeth | Extends HYP-3039 by naming `11` proof-facing micro-statements across HYP-3034..HYP-3038 and older HYP-3021/HYP-3022/HYP-3023/HYP-3029/HYP-3031 work.  Field histogram shows `status`, `route`, `owner`, `residual`, and `topology` as the recurring missing coordinates; the strongest hidden edge is residual first-tooth atlas to safe-stalk/barcode/normal-fan geometry. | HYP-3040, HYP-3039, HYP-3038, HYP-3037, HYP-3036, HYP-3035, HYP-3034, HYP-3033, HYP-3032, HYP-3031, HYP-3030, HYP-3029, HYP-3028, HYP-3027, HYP-3024, HYP-3023, HYP-3022, HYP-3021, HYP-3018, HYP-2963, THM-572, T1121, LTT-086 | Add `boundary_h1_owner_support`, `first_tooth`, `primitive_safe_deck_2_13`, `residual_capacitor_id`, `first_cut_stage`, `exact_M_zeta`, `endpoint_owner_strip`, `analytic_blindness_report`, and `automaton_shadow_class` to packet sidecars; prove the candidate layered packet principle or find the first leaking sidecar field. |
| LTI-189 | AP-tail q13 puncture / reciprocal fixed-point repair | AP-core one-tail rows, exact-period `q=13` punctures, reciprocal fixed points, mod-14 owner strips, safe-component stalk residuals | Proves the family `S_m={1,...,12,m}` by a two-clock split: `13 does not divide m` gives `t=1/13`, while `m=13s>=26` gives `t=s/(13s+1)`; the two HYP-3029 coarse-stalk residual pairs are exactly mod-14 owner-strip collisions that forgot this `m mod 13` coordinate, complementing the analytic sieve-clock warning that capacity clocks need blindness certificates and HYP-3033's residual-tooth scheduler | HYP-3041, HYP-3040, HYP-3039, HYP-3033, HYP-3032, HYP-3031, HYP-3029, HYP-3028, HYP-3027, HYP-3024, HYP-3023, HYP-3017, HYP-2963 | Before sending an AP-core tail residual to Fejer or THM-572, attach `ap_tail_certificate_kind`, `q13_puncture_bit`, and reciprocal fixed-point witness; then search for two-tail analogues whose hidden prime clock repairs coarse owner-strip collisions. |
| LTI-190 | Owner-strip filtration | Primitive safe decks, drop/add Haar zeta, endpoint-owner strip currents, residual proof pages | Synthesizes HYP-2997/HYP-3018 with HYP-3035..HYP-3038 and absorbs HYP-3041 as the concrete AP-tail example: post-status residuals should pass through primitive-period deck, Haar/drop-add zeta, and endpoint-owner strip current before route labels; a surviving counterexample must be invisible on all three pages or become named F7/THM-572 debt | HYP-3042, HYP-3041, HYP-3038, HYP-3037, HYP-3036, HYP-3035, HYP-3031, HYP-3018, HYP-2997, HYP-2963, THM-572, T1123, LTT-088, LTI-189, LTI-186, LTI-184, LTI-183, LTI-179 | Add `primitive_safe_deck_2_13`, `ap_tail_certificate_kind`, `q13_puncture_bit`, `drop_add_square_id`, `exact_M_zeta`, `endpoint_owner_strip_current`, `owner_strip_page`, and `first_surviving_filtration_page` to cached packet sidecars. |
| LTI-191 | Comprehensive LRC lens map | Lens families, preserved predicates, destroyed coordinates, required sidecars, handoff targets, proof-obligation tournament vertices | Creates `00-navigation/LRC-LENS-MAP.md` as the coordination layer above LTI-175/LTI-187/LTI-188: packet, topology, owner/barcode/normal-fan, arithmetic/Farey/period, harmonic/analytic, automaton/sequence, pair-decoy, residual/state-lift, formal, and external-analogy lenses all point to controlled forgetting over a labelled packet sheaf; HYP-3041's AP-tail puncture atlas and HYP-3042's owner-strip filtration are concrete examples of coarse owner-stalk data repaired by hidden period/current pages | HYP-3043, HYP-3042, HYP-3041, HYP-3040, HYP-3039, HYP-3038, HYP-3037, HYP-3036, HYP-3035, HYP-3034, HYP-3032, HYP-3024, HYP-3023, HYP-3022, HYP-3021, HYP-3018, HYP-3015, HYP-3012, HYP-3009, HYP-2997, HYP-2995, HYP-2963, THM-572, T1124, T1123, T1122, LTT-089, LTT-088, LTT-087, LTI-190, LTI-189, LTI-188, LTI-187, LTI-175 | Add `lens_family`, `preserved_lrc_predicate`, `destroyed_coordinate`, `required_sidecar`, `handoff_target`, `status_mixing_result`, `route_mixing_result`, `tournament_vertex_choice`, and `challenged_assumption` to new lens notes and packet-ledger experiments before promoting another scalar or analogy. |
| LTI-192 | Residual topology-exception teeth | Compact arc-topology failures, single-swap collar drops, coarse safe-component stalks, primitive safe-denominator decks | Zooms into the two HYP-3035 topology failures and shows they are exactly four strict-open single-swap collar rows at drops `9` and `11`; compact topology mixes `9->99/155` and `11->121/163`, but coarse stalk and primitive `q<=13` deck split both, with covering rows zero-deck and Q-witness rows first-safe at the drop speed | HYP-3044, HYP-3041, HYP-3040, HYP-3039, HYP-3038, HYP-3037, HYP-3036, HYP-3035, HYP-3034, HYP-3033, HYP-3031, HYP-3030, HYP-3029, HYP-3028, HYP-2963, THM-572, LTT-090, LTT-087, LTT-086, LTT-085, LTT-084, LTT-083, LTT-082, LTT-081 | Add `residual_topology_exception_id`, `topology_exception_drop`, `topology_exception_stalk_key`, and `topology_exception_first_primitive_q` to packet sidecars, then prove topology-failure implies owner-stalk collar plus primitive-deck split. |
| LTI-193 | Endpoint-owner transfer carrier | External boundary-owner strips, endpoint-owner transfer deltas, residue deltas, and largest-safe-component owner stalk joins inside coarse endpoint-current fibers | Refines HYP-3042's owner-strip filtration, complements HYP-3044's topology-exception collar audit, HYP-3039/HYP-3040, and HYP-3041 by making the hidden `B18Z6` owner coordinate explicit: all four audited residual capacitor packets have coarse endpoint word `B18Z6`, but external owner strips split the q=23 diagonal and both residual capacitors (`12:26x6,6:20x4`, `2:16x6`, `12:26x6,8:36x4`, `2:72x6`), refining both exact-M and topology first cuts without using route labels | HYP-3045, HYP-3044, HYP-3042, HYP-3041, HYP-3040, HYP-3039, HYP-3038, HYP-3037, HYP-3036, HYP-3035, HYP-3032, HYP-3031, HYP-3027, HYP-3026, HYP-3018, HYP-2963, THM-572, T1126, T1125, T1123, LTT-091, LTT-090, LTT-088, LTI-193, LTI-192, LTI-190, LTI-189, LTI-188, LTI-187, LTI-186, LTI-185 | Add `endpoint_owner_strip`, `endpoint_owner_transfer_delta`, `endpoint_owner_residue_delta`, `safe_component_owner_stalk`, and `owner_transfer_carrier` to residual manifests; prove that protected `B18Z6` residual packets route by owner transfer, reconstruct owners, are dual-annihilated, or emit named F7/THM-572 debt. |
| LTI-194 | Hidden connection accelerator index | Connector lemmas, source-marker audit, residual-section/Haar exit alphabet, endpoint-owner transfer, topology-exception teeth, exact-period decks, owner-strip filtration page, repair cochains, cocycle exactness | Merges hidden overlaps across HYP-3045/HYP-3044 and HYP-3038..HYP-3034 with older carriers: capacitor exits are residual-section/Haar exit codes, endpoint-owner transfer is the `B18Z6` address lift, topology exceptions are owner-stalk collars plus primitive q<=13 deck splits, the q=23 square is the first nested-refinement normal form, primitive decks are exact-period packet atlases, owner strips route through HYP-3042 plus normal-fan/path-lift support, q=14 and q13 punctures belong to exact-period/puncture guardrails, decoys are blocker decks, and squarefree blindness needs divisor-lattice fields | HYP-3046, HYP-3045, HYP-3044, HYP-3043, HYP-3042, HYP-3041, HYP-3040, HYP-3039, HYP-3038, HYP-3037, HYP-3036, HYP-3035, HYP-3034, HYP-3032, HYP-3031, HYP-3027, HYP-3022, HYP-3018, HYP-3013, HYP-3006, HYP-2996, HYP-2995, HYP-2992, HYP-2886, THM-523, THM-566, T1127, T1126, LTT-092, LTT-091, LTI-193, LTI-192, LTT-090, LTI-191, LTT-089 | Add the HYP-3046 sidecar merge set to the HYP-2963 packet ledger, beginning with `zeta_exit_class`, `residual_section_exit`, `coarse_endpoint_word`, `external_endpoint_owner_strip`, `endpoint_owner_transfer_delta`, `endpoint_owner_residue_delta`, `primitive_safe_deck_2_13`, topology-exception fields, AP-tail puncture/fixed-point fields, owner-support fields, `first_surviving_filtration_page`, `first_cut_stage`, `drop_add_square_id`, `omega_Q_class`, exact-period chart fields, `divisibility_threshold_qS`, divisor-lattice fields, and blocker-deck fields. |
| LTI-195 | A000568 k-depth perspective ladder | Rooted/node perspectives, k-depth in/out role colors, directed-edge sectors, directed-cycle conflicts, transitive-clique insertion cuts, rootless Burnside sidecars, endpoint-owner packet sheaves | Recasts the old A000568 perspective curiosity as a controlled-forgetting alarm: `P(5)=48` while `U(6)=56`, and k-depth node refinement reaches exact rooted type by depth `2`, so the missing eight classes are incident-word/cross-coupling payload rather than deeper node memory.  The companion Burnside audit shows a fixed-point-free `[3,3]` term for `U(6)` with `32` fixed tournaments and `0` fixed vertices, so the defect is rootless/cyclic and matches T1128 matrix observability: a leaked coordinate must survive as an observability, kernel, or Schur-complement sidecar. | HYP-3047, HYP-2120, HYP-2121, HYP-3046, HYP-3045, HYP-3043, HYP-3042, HYP-3041, HYP-3040, HYP-3039, HYP-3018, HYP-3015, HYP-1824, HYP-1825, THM-381, THM-385, T1129, T1128, T1127, LTT-093, LTT-092, LTT-091, LTT-089, LTT-088, LTT-087, OPEN-Q-108 | Build the exact map from 5-edge perspectives to 6-tournament classes, isolate the eight-class defect by edge-sector, cycle-chirality, and `[3,3]` rootless sidecars, then add `perspective_root_type`, `perspective_depth_k`, `observer_cut_position_word`, `incident_sector_deck`, `edge_zone_profile`, `cycle_relation_word`, `clique_root_shape`, and `cross_sector_orientation_word` sidecars to LRC threshold-packet experiments. |
| LTI-196 | Expanded tournament-matrix atlas | Matrix-result domains, skew signs, Hermitian `iS`, Laplacians, incidence/boundary matrices, transfer matrices, Smith sidecars, KKT/Farkas/SOS duals, observability matrices | Extends S210 with `165` additional classic matrix hooks (`300` total with S210) and reframes matrix invariants as proof carriers whose rows/columns may be edges, cycles, endpoint owners, denominator clocks, proof obligations, quotient fibers, cohomology classes, sidecar fields, or low-rank update directions rather than runners | HYP-3048, HYP-3047, HYP-3046, HYP-3043, HYP-3042, HYP-3040, HYP-3039, HYP-2121, HYP-2120, THM-381, THM-385, T1130, LTT-094, LTT-093, OPEN-Q-108 | Build the sidecar observability matrix: rows are coarse-fiber packet pairs, columns are hidden coordinates, and proof-safe sidecar sets separate, reconstruct, annihilate, descend, or name debt for every route/status-changing pair. |
| LTI-197 | Exact non-node perspective-depth carrier extension | Exact arc/triple rooted orbits, ordered-pair extension states, local edge/triple depth colors, edge tail/tip sectors, four-sector ordered-pair decks, cross-sector orientation, transitive/cyclic triple splits, cycle-conflict roots, matrix sidecar observability | Extends HYP-3047 and HYP-3049 while using HYP-3048's matrix sidecar language. HYP-3049 first shows rooted 5-perspective plus incident word is exact ordered-pair perspective on `U(6)` (`1408=1408`), with directed-edge quotient `704` and cross-sector orientation separating sector decks from `55/56` to `56/56`. HYP-3050 then computes broader non-node carrier totals through the first failure: at `m=5`, directed-edge perspectives and triple perspectives both total `88`, splitting triples as `64` transitive and `24` cyclic. Under the stricter directed-WL convention, node depth reaches exact rooted node orbits by depth `3`, while local edge/triple depth `2` recover exact first-failure carriers. The missing `8` A000568 states remain observer-extension cut data; edge sectors and cyclic triples are sidecars for pair-good switches, residual capacitors, and Omega/conflict payloads. | HYP-3054, HYP-3050, HYP-3049, HYP-3048, HYP-3047, HYP-3040, HYP-3039, HYP-2210, HYP-2120, HYP-1978, HYP-1977, HYP-1824, HYP-1825, THM-381, THM-385, THM-260, THM-409, T1136, T1132, T1131, T1130, T1129, LTT-099, LTT-095, LTT-094, LTT-093, LTI-196, LTI-195, LTI-194, LTI-188, LTI-187 | Add `observer_endpoint_role`, `incident_word`, `ordered_pair_sector_deck`, `cross_sector_orientation_word`, `edge_tail_tip_sector_word`, `triple_perspective_kind`, `cycle_conflict_pair_id`, `observer_extension_cut_signature`, and sidecar-observability matrix columns; extend exact carrier counts to `m=6` and compare edge-sector signatures against pair-good blocker teeth and residual capacitor IDs. |
| LTI-198 | Rooted layer-extension flow | Fixed-path layer boundary words, rank-one/coboundary cross-sheets, parent automorphism cycle indices, incident-word orbits, rooted children, unrooting collision fibers | Adds HYP-3051's enumeration-side controlled-forgetting rule: the complete bipartite sheet between layer sizes `k` and `k+1` has `k^2+k` apparent line bits but only `2k` boundary-potential bits when `e_ij=x_i XOR y_j`; one-vertex tournament extension similarly counts incident-word orbits under `Aut(parent)`, giving `E(n->n+1)=R(n+1)` before the final unrooting quotient.  The old perspective coincidence `R(3)=A(4)` and `R(4)=A(5)` fails at `R(5)=48 != A(6)=56`, so rooted extension state is the theorem carrier, not raw A000568. | HYP-3051, HYP-3050, HYP-3049, HYP-3048, HYP-3047, HYP-3046, HYP-3043, HYP-3039, T1133, T1132, T1131, T1130, T1129, T1128, LTT-096, LTT-095, LTT-094, LTT-093, LTT-092, LTT-089, LTT-085, LTT-048, LTT-036, LTT-006, A000568, half-tilings | Add `parent_class`, `root_orbit`, `incident_word_orbit`, `layer_boundary_word`, `rank_one_sheet_id`, `rectangle_defect_rank`, and `unrooting_collision_fiber` to tournament/LRC packet sidecars before promoting an unrooted class quotient. |
| LTI-199 | Diagonal-layer transport law | Parent isomorphism classes, diagonal word orbits under parent automorphisms, `K_{k,k+1}` position-line profiles, aligned triangle-flow sidecars, deletion-parent fibers, half-tiling branch data | Complements HYP-3051's rooted layer-extension flow with the exact unrooted sink and half-tiling guardrail: `parent class + diagonal word orbit under Aut(parent) -> rooted child -> unrooted child sink`.  At `5 -> 6`, raw labelled extensions are `384`, word orbits/rooted children are `296`, and all `56` child sinks are reached.  The `k^2+k` line carrier is geometric and rank-one-labelled by two words (`N00*N11=N01*N10`), while HYP-3049 identifies the same incident word as an ordered-pair sector deck whose first missing column is `cross_sector_orientation_word`; half-tilings mark the converse-fold branch locus rather than enumerating isomorphism classes | HYP-3052, HYP-3051, HYP-3050, HYP-3049, HYP-3048, HYP-3047, HYP-3046, HYP-3043, HYP-3039, HYP-3031, HYP-2685, HYP-2690, HYP-2120, HYP-2121, THM-549, THM-550, THM-381, THM-385, T1134, T1133, T1132, T1131, T1130, LTT-097, LTT-096, LTT-095, LTT-094, LTI-198, LTI-197, LTI-196, LTI-195, LTT-093, OPEN-Q-108 | Build the displayed `5 -> 6` child-class ledger with deletion-parent profiles, rooted orbit counts, self-converse status, score sequences, and aligned triangle-flow summaries; then add `diagonal_word_orbit`, `K_position_line_profile`, `aligned_pair_counts`, `newest_link_bit`, `cross_sector_orientation_word`, and `deletion_parent_profile` sidecars to threshold-packet experiments; emit the low-rank matrix as a HYP-3048/HYP-3051 sidecar observability row. |
| LTI-200 | Fixed-path diagonal-layer flow | Binary half-tiling diagonal layers, `K_{k,k+1}` coboundary line flow, rectangle and hourglass cycle laws, fixed Hamiltonian-path fibers, path-reversal/converse diagonal quotient | Complements HYP-3052's diagonal transport and HYP-3051's rooted layer-extension law with the fixed-path/global fiber audit: each `k^2+k` bridge is a `GF(2)` cut-space presentation with `k(k+1)` lines but rank `2k`, with `k(k-1)` 4-cycle redundancies; globally adjacent-layer flow has rank `#tiles-1` and redundancy `2*C(n-1,3)+C(n-2,2)`; A000568 is the full `S_n` orbit count, while fixed-path half-tilings cover it with fiber `H(T)/|Aut(T)|` and diagonal `Z2` symmetry is only a sidecar quotient | HYP-3053, HYP-3052, HYP-3051, HYP-3050, HYP-3049, HYP-3048, HYP-3047, HYP-3043, HYP-3039, HYP-3031, HYP-2991, HYP-2989, HYP-2120, HYP-2121, THM-381, THM-385, T1135, T1134, T1133, T1132, T1131, LTT-098, LTT-097, LTT-096, LTT-095, LTT-094, OPEN-Q-108 | Emit rectangle/hourglass cycle bases for layer bridges, join them to HYP-3052 diagonal transport, HYP-3051 rooted extension fibers, HYP-3049 cross-sector orientation, HYP-3050 exact edge/triple carriers, and matrix sidecar observability, and add `tile_layer_id`, `interlayer_bridge_id`, `rectangle_cycle_basis_id`, `line_potential_word`, `cycle_defect_word`, and owner/barcode support fields to packet experiments. |
| LTI-201 | Observer-extension cut payload law | Quotient, next operation, observer roots, incident words, ordered-pair/edge sectors, cross-sector orientation, deletion-parent fibers, rectangle/hourglass residues, endpoint-owner payloads, proof-obligation sidecars, sidecar observability columns | Abstracts the HYP-3050 observer-extension/cut defect into the general controlled-forgetting test and corrects the first shifted A000568 arithmetic to `48+8=56`, not `48+12=56`.  The recurring `12` is a fold/parent/fixed-locus count (`R(4)`, `U(5)`, both `5->6` source/sink deletion slices, and `SC(6)`), while the additive defect is the first observer-extension payload.  A quotient is legal only after the payload for its next outside operation is retained, reconstructed, dual/cut/cocycle-annihilated, descended, boundary-stopped, or routed to named residual debt.  This reclassifies pair-good decoys by blocker generator and active-owner/barcode relation, residual capacitors by first cut/zeta/endpoint-owner exit, AP-tail collisions by q13 clock, automaton shadows by magnitude/topology/owner handoff, diagonal-layer lines by rectangle/hourglass defect, and matrix invariants by observability columns. | HYP-3054, HYP-3053, HYP-3052, HYP-3051, HYP-3050, HYP-3049, HYP-3048, HYP-3047, HYP-3046, HYP-3045, HYP-3043, HYP-3040, HYP-3039, HYP-3038, HYP-3037, HYP-3034, HYP-3027, HYP-3024, HYP-3022, HYP-3021, HYP-3018, HYP-2997, HYP-2995, HYP-2991, HYP-2989, HYP-2963, THM-381, THM-385, THM-572, T1136, T1135, T1134, T1133, T1132, T1131, LTT-099, LTT-098, LTT-097, LTT-096, LTT-095, LTT-094, LTT-093, LTI-200, LTI-199, LTI-198, LTI-197, LTI-196, LTI-187, OPEN-Q-108 | Add manifest fields `quotient`, `next_operation`, `observer_extension_payload`, `observer_payload_stage`, `incident_word_orbit`, `edge_sector_cross_orientation`, `deletion_parent_fiber`, `rectangle_hourglass_residue`, `endpoint_owner_payload`, `barcode_active_owner_support`, `extension_address`, `cut_or_cycle_defect`, `route_owner_certificate`, and `payload_exit/discharge_rule`; build HYP-3048 observability matrices whose rows are coarse-fiber packet pairs and whose columns are these payloads before trusting scalar counts. |
| LTI-202 | Duodecimal observer-extension cut payload | Dozen control slice, first-failure inclusion-exclusion, observer-extension cut payload, incident-word orbit, cross-sector orientation, deletion-parent fibers, rectangle/hourglass residues | HYP-3055 specializes HYP-3054's calculus to the first A000568 defect: `48+12=60`, while `U(6)=P(5)+U(5)-U(4)=48+12-4=56`.  The dozen remains structural because `P(4)=U(5)=SC(6)=12`, and the missing `8` is `SC(6)-U(4)`.  The reusable rule is that rooted/node caches, apex iso classes, fixed-path tilings, and raw line counts may forget observer-extension coordinates only after incident-word orbit, endpoint role, cross-sector orientation, deletion-parent fiber, and rectangle/hourglass residues are reconstructed, annihilated, descended, or routed to named residual debt | HYP-3055, HYP-3054, HYP-3053, HYP-3052, HYP-3051, HYP-3050, HYP-3049, HYP-3048, HYP-3047, HYP-3043, HYP-3039, HYP-3031, HYP-2991, HYP-2989, HYP-2928, HYP-2120, HYP-2121, THM-381, THM-385, T1137, T1136, T1135, T1134, T1133, T1132, T1131, LTT-100, LTT-099, LTT-098, LTT-097, LTT-096, LTT-095, LTT-094, LTI-201, LTI-200, LTI-199, LTI-198, LTI-197, LTI-196, OPEN-Q-108 | Build the sidecar observability matrix promised by HYP-3055: rows are coarse-merged tournament or packet pairs, columns are `incident_word_orbit`, `observer_endpoint_role`, `ordered_pair_sector_deck`, `cross_sector_orientation_word`, `deletion_parent_profile`, `rectangle_residue`, `hourglass_residue`, `self_converse_status`, `endpoint_owner_packet`, and `binding_scale`; then transfer the same columns to LRC endpoint-owner packets, Haar zeta grids, pair-good decoys, residual capacitors, and tournament spectra. |
| LTI-203 | Observer-cut orbit ledger | Visible-fiber automorphism groups, cut-payload orbits, observer kinds, changed LRC predicates, sidecar columns, discharge modes | Refines LTI-201 from "retain the next payload" to an auditable orbit object `C_q(x,o)=orbit_Aut_q(x)(boundary slice, incidence word, extended shadow)`.  The ledger records `base_quotient`, `fiber_id`, `observer_kind`, `visible_automorphism_group`, `cut_payload_orbit_id`, `changed_lrc_predicate`, `separating_sidecar`, `discharge_mode`, and `residual_debt_name`.  It applies the same test to A000568 incident/cross-sector orientation, AP/GW closed-H1 owner support, q=23 Haar zeta plus endpoint-owner strip, K33 state lifts, automaton exact-packet handoffs, analytic blindness reports, diagonal-layer cycle residues, and matrix observability columns. | HYP-3056, HYP-3055, HYP-3054, HYP-3053, HYP-3052, HYP-3051, HYP-3050, HYP-3049, HYP-3048, HYP-3047, HYP-3043, HYP-3039, HYP-3037, HYP-3034, HYP-3032, HYP-3031, HYP-3024, HYP-3018, HYP-2997, HYP-2995, HYP-2963, THM-572, T1138, T1137, T1136, LTT-101, LTT-100, LTT-099, LTI-202, LTI-201, OPEN-Q-108 | Build an observer-cut ledger over HYP-2963 coarse fibers; for each mixed route/status pair enumerate admissible next observers, quotient automorphisms, payload orbit ids, separating sidecars, and discharge modes before promoting any scalar or automaton quotient. |
| LTI-204 | Tournament value-origin ledger | Value-origin types for small integer shadows: class counts, rooted perspectives, self-converse fixed branches, incident-word orbits, deletion fibers, edge-sector sidecars, fixed-path fibers, rectangle/hourglass residues | Corrects the first shifted A000568/rooted failure to `48+8=56` and reframes the `12` signal as a typed diagonal alignment `R(4)=U(5)=SC(6)=12`.  The theorem carrier is a ledger of where each number comes from before it is used: `U(5)=12`, raw incident extensions `384`, rooted word orbits `R(6)=296`, unrooted sinks `U(6)=56`, ordered-pair perspectives `1408`, and directed-edge perspectives `704`.  Cross-sector orientation repairs the only U(6) edge-sector collision, while S217 rectangle/hourglass residues are the fixed-path duplication analogue | HYP-3057, HYP-3054, HYP-3053, HYP-3052, HYP-3051, HYP-3050, HYP-3049, HYP-3048, HYP-3047, HYP-3039, HYP-2991, HYP-2989, THM-381, THM-385, T1139, T1136, T1135, T1134, T1133, T1132, T1131, LTT-102, LTT-099, LTT-098, LTT-097, LTT-096, LTT-095, LTT-094, OPEN-Q-108 | Add `value_origin_type`, `parent_class`, `incident_word_orbit`, `root_orbit_count`, `child_sink_class`, `deletion_parent_profile`, `self_converse_status`, `fixed_path_fiber`, `edge_sector_cross_orientation_word`, `rectangle_residue_class`, `hourglass_residue_class`, and `lost_coordinate_exit` to tournament/LRC packet experiments before trusting small numerical coincidences. |
| LTI-205 | Hyperbolic reciprocal packet sidecar | Three-lane packet orders, reciprocal sums, orbifold Euler sign, triangle signatures, Fermat-Catalan power guards, hyperbolic debt discharge routes | Imports `1/a+1/b+1/c<1` as the negative-orbifold-Euler / hyperbolic triangle threshold rather than a proof shortcut: define `hyperbolic_triple_signature`, `reciprocal_sum`, `curvature_margin=1-sum`, `orbifold_euler=sum-1`, and `hyperbolic_debt_discharge_route`; use the first hyperbolic `(2,3,7)` signature and its `41/42` margin only as route pressure near `14=2*7`, `q=27=3^3`, `3/41`, C27, K33, and AP/GW | HYP-3058, HYP-3055, HYP-3054, HYP-3043, HYP-3039, HYP-3012, HYP-3009, HYP-3003, HYP-3002, HYP-2998, HYP-2963, HYP-2945, HYP-2937, HYP-2934, HYP-2928, THM-572, T1140, T1137, T1136, LTT-103, LTT-100, LTT-099, LTI-202, LTI-201, OPEN-Q-108 | Pick HYP-2963 packet families, define honest triples from exact scale, route incidence, automaton depth, observer-cut depth, primitive-period deck, Fejer degree, or state-lift obligation, then test spherical/euclidean/hyperbolic sign against discharge route without dropping exact `M`, endpoint-owner, topology, deletion-fiber, rectangle/hourglass, or certificate sidecars. |
| LTI-206 | Observer-extension cut payload | Incident words, source/sink slices, deletion fibers, self-converse branch loci, ordered-pair sector decks, rectangle/hourglass residues, endpoint owners, sidecar observability matrices | Refines HYP-3056's observer-cut orbit ledger, HYP-3055's duodecimal bridge, and HYP-3054's general calculus with the exact first-defect arithmetic: `48+12=60` but `U(6)=P(5)+U(5)-U(4)=48+12-4=56`; audits the nearby `12` carriers (`P(4)`, `U(5)`, source/sink `5->6` slices, `SC(6)`, and `K_{4,5}` rectangle redundancy), source/sink overlap `4`, `SC(6)` deletion decks touching all `12` parents, and cross-sector decks separating `56/56`; abstracts these into `base quotient + observer/cut word + stabilizer + sidecar + sink map + legal exit` as the controlled-forgetting protocol | HYP-3059, HYP-3056, HYP-3055, HYP-3054, HYP-3053, HYP-3052, HYP-3051, HYP-3050, HYP-3049, HYP-3048, HYP-3047, HYP-3046, HYP-3043, HYP-3039, HYP-3031, HYP-3013, HYP-3008, HYP-2120, HYP-2121, THM-381, THM-385, T1141, T1138, T1137, T1136, T1135, T1134, T1133, T1132, T1131, LTT-104, LTT-101, LTT-100, LTT-099, LTT-098, LTT-097, LTT-096, LTT-095, LTT-094, LTI-203, LTI-202, OPEN-Q-108 | Build the observer/cut-payload theorem: for each quotient, retain `observer_cut_word`, `stabilizer_orbit`, `source_sink_slice_id`, `deletion_fiber_payload`, `self_converse_branch_bit`, `cross_sector_orientation_word`, `rectangle_hourglass_residue`, `endpoint_owner_packet`, and `legality_exit`, then prove the payload is retained, reconstructed, annihilated, descended, AP/GW equality, or named residual debt. |
| LTI-207 | Desargues/Beal finalizer carrier | Desargues girth-six incidence residue, Beal common-owner gate, labelled packets, observer-cut orbits, rectangle/hourglass residues | Adds the post-rectangle finalizer test for LRC14: after S217 rectangle and hourglass `4`-cycle/coboundary residues vanish, a survivor may still carry a girth-six incidence address.  The Desargues graph is the canonical scout model (`20` vertices, `30` edges, cubic bipartite, girth `6`, automorphism count `240`, cycle counts `{6:20,8:30,10:132}`).  The Beal cue is used only as a proof-owner rule: a primitive three-channel equality/collision should retain a common owner/prime/packet coordinate before quotienting, or else be named F7/THM-572 debt. | HYP-3060, HYP-3057, HYP-3056, HYP-3055, HYP-3054, HYP-3053, HYP-3052, HYP-3051, HYP-3048, HYP-3037, HYP-3034, HYP-3031, HYP-2991, HYP-2963, THM-572, T1142, LTT-105, LTI-204, LTI-203, LTI-201, OPEN-Q-108 | Add `desargues_girth6_residue` and `beal_common_owner_gate` to packet/fiber ledgers after `rectangle_residue_class` and `hourglass_residue_class`; prove a residual with both gates closed discharges through existing sidecars, and a residual with an open gate descends to a family or becomes named F7/THM-572 debt. |

## Cross-Cutting Guardrails


- A scalar that works on named rows is not a proof carrier until its fibers are
  checked for mixed LRC routes.
- A tournament invariant on runners is usually a shadow.  Ask what marks,
  thresholds, endpoints, or packet labels were forgotten.
- AP/Goddyn-Wong are boundary equality atoms.  K33 is a state-lift / residual
  route, not the same kind of atom.
- Exact `M=p/q`, `qdiv`, and Farey excess should be retained before any
  product, summand, power, divisor, or spectral quotient.
- Representation counts must name their economy: smoothing, bounded arity,
  normal form, or Farey address.  Otherwise they are raw scalar shadows.
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
| LTI-TODO-11 | Packet-level Haar zipper cocycle signatures | LTI-046, LTI-047, LTI-048, LTI-106 |
| LTI-TODO-12 | Product-rule packet grid for fixed-path tilings | LTI-047, LTI-048, LTI-106, LTI-107 |
| LTI-TODO-13 | Packet cocycle theorem formalization | LTI-037, LTI-040, LTI-056, LTI-061, LTI-109 |
| LTI-TODO-14 | Executable F7 cocycle residual record | LTI-038, LTI-106, LTI-110 |
| LTI-TODO-15 | Cocycle-sheaf obstruction matrix over HYP-2963 packets | LTI-108, LTI-109, LTI-110, LTI-117, LTI-145 |
| LTI-TODO-16 | Namespace collision audit and repair plan for reused IDs | LTI-104, LTI-105, LTI-146 |
| LTI-TODO-17 | Deck-derivative endpoint reconstruction table | LTI-111, LTI-126, LTI-134, LTI-141 |
| LTI-TODO-18 | Burnside/A000568 marked quotient tax report | LTI-112, LTI-135, LTI-141, LTI-143 |
| LTI-TODO-19 | F0-F7 residual metagraph Laplacian and SCC atlas | LTI-013, LTI-022, LTI-113, LTI-142 |
| LTI-TODO-20 | Path-homology/Pfaffian residual classifier | LTI-039, LTI-110, LTI-117, LTI-129 |
| LTI-TODO-21 | Binding-pair switch tournament for covering rows | LTI-031, LTI-050, LTI-144, LTI-145 |
| LTI-TODO-22 | Product-rule irreducibility no-lift compression test | LTI-107, LTI-136, LTI-137, LTI-139 |
| LTI-TODO-23 | Packet-level first-nonzero cocycle ledger | LTI-044, LTI-046, LTI-056, LTI-147 |
| LTI-TODO-24 | Familywise residual-section templates for HYP-2996 | LTI-047, LTI-107, LTI-109, LTI-148 |
| LTI-TODO-25 | Additive-basis / Farey-operator packet schema | LTI-076, LTI-077, LTI-078, LTI-089, LTI-130, LTI-149 |
| LTI-TODO-26 | Additive-basis proof-currency packet classifier | LTI-063, LTI-090, LTI-094, LTI-149, LTI-150 |
| LTI-TODO-27 | Family-scale smoothing switchboard for HYP-2963 packets | LTI-035, LTI-049, LTI-057, LTI-061, LTI-065, LTI-151 |
| LTI-TODO-28 | Curried packet evaluator manifest | LTI-147, LTI-148, LTI-149, LTI-150, LTI-152 |
| LTI-TODO-29 | HYP-2963 dichotomy-mode packet fields | LTI-011, LTI-063, LTI-076, LTI-150, LTI-151, LTI-152, LTI-153, LTI-154 |
| LTI-TODO-30 | Representation-economy field for sequence-shadow packet classifiers | LTI-076, LTI-089, LTI-109, LTI-149, LTI-150, LTI-153, LTI-154, LTI-155 |
| LTI-TODO-31 | Multiverse annex promotion queue | LTI-207, LTI-206, LTI-205, LTI-204, LTI-203, LTI-202, LTI-201, LTI-200, LTI-199, LTI-198, LTI-197, LTI-196, LTI-195, LTI-194, LTI-193, LTI-192, LTI-191, LTI-190, LTI-189, LTI-188, LTI-187, LTI-186, LTI-185, LTI-184, LTI-183, LTI-182, LTI-181, LTI-180, LTI-179, LTI-178, LTI-177, LTI-176, LTI-175, LTI-174, LTI-173, LTI-172, LTI-171, LTI-170, LTI-169, LTI-168, LTI-167, LTI-166, LTI-165, LTI-164, LTI-163, LTI-162, LTI-161, LTI-160, LTI-159, LTI-158, LTI-157, LTI-156, LTI-155, LTI-154, LTI-153, LTI-152, LTT-105, LTT-104, LTT-103, LTT-102, LTT-101, LTT-100, LTT-099, LTT-098, LTT-097, LTT-096, LTT-095, LTT-094, LTT-093, LTT-092, LTT-091, LTT-090, LTT-089, LTT-088, LTT-087, LTT-086, LTT-085, LTT-084, LTT-083, LTT-082, LTT-081, LTT-080, LTT-079, LTT-078, LTT-077, LTT-076, LTT-075, LTT-074, LTT-073, LTT-072, LTT-071, LTT-070, LTT-069, LTT-068, LTT-067, LTT-066, LTT-065, LTT-064, LTT-063, LTT-062, LTT-061, LTT-060, LTT-059, LTM-*, CPI-* |
| LTI-TODO-32 | Mixed automatic-word family templates after HYP-3023 zipper audit | LTI-177, LTI-176, LTI-175, LTI-174, LTI-173, LTI-172, LTI-171, LTI-170, LTI-169, LTI-168, LTI-167, LTI-166, LTI-165, LTI-164, LTI-163, LTI-162, LTI-161, LTI-160, LTI-159, LTI-158, LTI-155, LTI-150, LTI-151, LTI-152, LTT-076, LTT-075, LTT-074, LTT-073, LTT-072, LTT-071, LTT-070, LTT-069, LTT-068, LTT-067, LTT-066, LTT-065, LTT-064, LTT-063, LTT-062, LTT-061 |
| LTI-TODO-33 | Automaton-state exact-gap side-channel audit | LTI-089, LTI-130, LTI-149, LTI-150, LTI-158, LTI-159, LTI-160, LTI-161, LTI-162, LTI-163, LTI-164, LTI-165, LTI-166, LTI-167, LTI-168, LTI-169, LTI-170, LTI-177, LTI-176, LTI-175, LTI-174, LTI-173, LTT-076, LTT-075, LTT-074, LTT-073, LTT-072, LTT-071, LTT-070, LTT-069, LTT-068, LTT-067, LTT-066, LTT-065, LTT-064, LTT-063, LTT-062 |
| LTI-TODO-34 | DFA/gap-block safe-component fields for HYP-2963 packets beyond the HYP-3017 sidecars | LTI-089, LTI-130, LTI-131, LTI-158, LTI-160, LTI-161, LTI-162, LTI-163, LTI-164, LTI-165, LTI-166, LTI-167, LTI-168, LTI-169, LTI-170, LTI-177, LTI-176, LTI-175, LTI-174, LTI-173, LTT-076, LTT-075, LTT-074, LTT-073, LTT-072, LTT-071, LTT-070, LTT-069, LTT-068, LTT-067, LTT-066, LTT-065, LTT-064, LTT-063 |
| LTI-TODO-35 | Perfect-product/divisor packet merge audit | LTI-084, LTI-155, LTI-158, LTI-159, LTI-161, LTI-162, LTI-163, LTI-164, LTI-165, LTI-166, LTI-167, LTI-168, LTI-169, LTI-170, LTI-177, LTI-176, LTI-175, LTI-174, LTI-173, LTT-076, LTT-075, LTT-074, LTT-073, LTT-072, LTT-071, LTT-070, LTT-069, LTT-068, LTT-067, LTT-066, LTT-065, LTT-064 |
| LTI-TODO-36 | Creative exact packet-lens sidecar audit | LTI-177, LTI-176, LTI-175, LTI-174, LTI-173, LTI-172, LTI-171, LTI-170, LTI-169, LTI-168, LTI-167, LTI-166, LTI-163, LTI-147, LTI-148, LTI-151, LTI-152, LTI-160, LTI-161, LTT-076, LTT-075, LTT-074, LTT-073, LTT-072, LTT-071, LTT-070, LTT-069, LTT-068, LTT-067, LTT-065 |
| LTI-TODO-37 | Full-bank lonely-profile barcode sidecar audit | LTI-177, LTI-176, LTI-175, LTI-174, LTI-173, LTI-172, LTI-171, LTI-170, LTI-169, LTI-168, LTI-167, LTI-166, LTI-164, LTI-138, LTI-099, LTI-089, LTI-158, LTI-160, LTI-161, LTI-162 |
| LTI-TODO-38 | Magnitude-cocycle fiber-purity theorem | LTI-177, LTI-176, LTI-175, LTI-174, LTI-173, LTI-172, LTI-171, LTI-170, LTI-169, LTI-168, LTI-167, LTI-166, LTI-165, LTI-164, LTI-163, LTI-161, LTI-158, LTI-147, LTT-076, LTT-075, LTT-074, LTT-073, LTT-072, LTT-071, LTT-070, LTT-069, LTT-068, LTT-067, LTT-066, LTT-065, LTT-063 |

| LTI-TODO-39 | Normal-fan packet sidecar purity theorem | LTI-166, LTI-165, LTI-164, LTI-099, LTI-044, LTI-057, LTI-147, LTI-151 |
| LTI-TODO-40 | Discrepancy-height trident compression on the full HYP-2963 bank | LTI-171, LTI-167, LTI-166, LTI-165, LTI-164, LTI-159, LTI-147, LTI-151, LTT-070, LTT-067, LTT-066, LTT-044 |
| LTI-TODO-41 | Pair-good decoy barcode/normal-fan blocker-deck grammar theorem | LTI-169, LTI-168, LTI-166, LTI-144, LTI-044, LTI-057, LTI-099, LTI-151 |
| LTI-TODO-42 | Full-bank closed arc-Cech nerve and quotient-defect audit | LTI-172, LTI-041, LTI-044, LTI-045, LTI-097, LTI-147, LTT-069 |
| LTI-TODO-43 | Coarse ET + Henselian-unit status convergence theorem | LTI-181, LTI-176, LTI-171, LTI-170, LTI-167, LTI-165, LTI-164, LTI-151, LTT-079, LTT-074, LTT-070, LTT-067, LTT-066, LTT-044 |
| LTI-TODO-44 | Carrier-fusion sidecar purity theorem | LTI-177, LTI-176, LTI-175, LTI-174, LTI-173, LTI-172, LTI-171, LTI-170, LTI-169, LTI-168, LTI-167, LTI-166, LTI-165, LTI-164, LTI-163, LTI-159, LTI-158, LTI-076, LTI-097, LTT-076, LTT-075, LTT-074, LTT-073, LTT-072, LTT-071, LTT-070, LTT-069, LTT-068, LTT-067, LTT-066, LTT-065, LTT-061 |
| LTI-TODO-45 | Side-channel repair zipper theorem | LTI-174, LTI-173, LTI-172, LTI-171, LTI-170, LTI-169, LTI-168, LTI-167, LTI-166, LTI-165, LTI-164, LTI-163, LTI-147, LTI-109, LTT-072, LTT-071, LTT-070, LTT-069, LTT-068, LTT-067, LTT-066, LTT-065 |
| LTI-TODO-46 | Carrier-pullback row instantiation queue | LTI-175, LTI-174, LTI-173, LTI-156, LTI-147, LTI-021, LTT-073, LTT-072, LTT-071, LTT-070, LTT-069, LTM-071, LTM-072, CPI-* |
| LTI-TODO-47 | Residual status-gate packet-ledger teeth | LTI-181, LTI-176, LTI-171, LTI-170, LTI-167, LTI-165, LTI-151, LTT-079, LTT-074, LTT-070, LTT-067, LTT-066, HYP-3033, HYP-3028 |
| LTI-TODO-48 | Full-bank largest safe-component stalk descent stress | LTI-177, LTI-176, LTI-175, LTI-174, LTI-173, LTI-172, LTI-171, LTI-170, LTI-166, LTI-164, LTI-057, LTI-044, LTI-041, LTT-076, LTT-075, LTT-074, LTT-073, LTT-072, LTT-071, LTT-070, LTT-069, OPEN-Q-108 |
| LTI-TODO-49 | Status-topology gate / zero-open arc-boundary cycle theorem | LTI-182, LTI-181, LTI-178, LTI-177, LTI-176, LTI-172, LTI-171, LTI-166, LTI-164, LTI-147, LTT-080, LTT-079, LTT-076, LTT-075, LTT-074, LTT-070, LTT-069, THM-572, HYP-3034, HYP-3033, HYP-3030 |
| LTI-TODO-50 | Haar-tile repair class for residual mixed pairs | LTI-193, LTI-189, LTI-188, LTI-187, LTI-186, LTI-185, LTI-184, LTI-183, LTI-182, LTI-181, LTI-180, LTI-179, LTI-178, LTI-177, LTI-176, LTI-174, LTI-148, LTI-147, LTI-109, LTI-108, LTI-107, LTI-106, LTI-046, LTI-047, LTT-091, LTT-087, LTT-086, LTT-085, LTT-084, LTT-083, LTT-082, LTT-081, LTT-080, LTT-079, LTT-078, LTT-077, LTT-076, LTT-075, LTT-074, LTT-072, LTT-006, HYP-3045, HYP-3041, HYP-3040, HYP-3039, HYP-3038, HYP-3037, HYP-3036, HYP-3035, HYP-3034, HYP-3033 |
| LTI-TODO-51 | Analytic repair-clock zipper theorem | LTI-180, LTI-174, LTI-171, LTI-151, LTI-065, LTI-064, LTI-061, LTI-057, LTI-049, LTI-035, LTT-078, LTT-072, LTT-070, LTT-044, HYP-3032 |
| LTI-TODO-52 | Residual certificate teeth theorem | LTI-181, LTI-179, LTI-178, LTI-176, LTI-171, LTT-079, LTT-077, LTT-076, LTT-074, HYP-3033 |
| LTI-TODO-53 | Residual first-tooth sidecar manifest | LTI-183, LTI-182, LTI-181, LTI-179, LTI-178, LTI-177, LTI-176, LTI-171, LTI-170, LTT-081, LTT-080, LTT-079, LTT-077, LTT-076, LTT-075, LTT-074, HYP-3035, HYP-3034, HYP-3033 |
| LTI-TODO-54 | Primitive-period deck full-bank packet sidecar | LTI-184, LTI-181, LTI-180, LTI-179, LTI-178, LTI-176, LTI-171, LTI-061, LTI-062, LTT-082, LTT-079, LTT-078, LTT-077, LTT-074, CPI-043, HYP-3036 |
| LTI-TODO-55 | Residual capacitor family cut audit | LTI-193, LTI-188, LTI-187, LTI-185, LTI-184, LTI-183, LTI-182, LTI-181, LTI-180, LTI-179, LTI-178, LTI-177, LTI-176, LTI-174, LTI-172, LTI-171, LTI-151, LTT-091, LTT-086, LTT-085, LTT-083, LTT-082, LTT-081, LTT-080, LTT-079, LTT-078, LTT-077, LTT-076, LTT-075, HYP-3045, HYP-3040, HYP-3039, HYP-3037, HYP-3036, HYP-3035, HYP-3034, HYP-3033, HYP-3032, HYP-3028 |
| LTI-TODO-56 | Double-pair drop/add square family audit | LTI-193, LTI-188, LTI-187, LTI-186, LTI-185, LTI-184, LTI-180, LTI-179, LTI-174, LTI-166, LTI-164, LTI-047, LTI-046, LTT-091, LTT-086, LTT-085, LTT-084, LTT-083, LTT-082, LTT-078, LTT-077, HYP-3045, HYP-3040, HYP-3039, HYP-3038, HYP-3037, HYP-3036 |
| LTI-TODO-57 | Hidden-coordinate ledger and anti-wedge debt audit | LTI-187, LTI-186, LTI-185, LTI-184, LTI-183, LTI-182, LTI-181, LTI-180, LTI-179, LTI-178, LTI-176, LTI-169, LTI-168, LTI-166, LTT-085, LTT-084, LTT-083, LTT-082, LTT-081, LTT-080, LTT-079, LTT-078, LTT-077, LTT-076, HYP-3039, HYP-3038, HYP-3037, HYP-3036, HYP-3035, HYP-3034, HYP-3033, HYP-3022, HYP-3021 |
| LTI-TODO-58 | Hidden-statement packet-sidecar synthesis | LTI-188, LTI-187, LTI-186, LTI-185, LTI-184, LTI-183, LTI-182, LTI-181, LTI-180, LTI-179, LTI-178, LTI-177, LTT-086, LTT-085, LTT-084, LTT-083, LTT-082, LTT-081, LTT-080, LTT-079, LTT-078, HYP-3040, HYP-3039, HYP-3038, HYP-3037, HYP-3036, HYP-3035, HYP-3034, HYP-3033, HYP-3032, HYP-3029, HYP-3023, HYP-3022, HYP-3021, HYP-3018 |
| LTI-TODO-59 | Owner-strip filtration full-bank ledger | LTI-190, LTI-189, LTI-187, LTI-186, LTI-185, LTI-184, LTI-183, LTI-179, LTI-166, LTI-147, LTT-088, LTT-087, LTT-085, LTT-084, LTT-083, LTT-082, LTT-081, HYP-3042, HYP-3041, HYP-3039, HYP-3038, HYP-3037, HYP-3036, HYP-3035 |
| LTI-TODO-60 | Lens-family packet manifest and status/route audit | LTI-191, LTI-190, LTI-189, LTI-188, LTI-187, LTI-175, LTI-174, LTI-173, LTI-172, LTI-171, LTI-170, LTI-169, LTI-168, LTI-166, LTI-164, LTI-147, LTT-089, LTT-088, LTT-087, LTT-086, LTT-085, LTT-073, LTT-072, LTT-071, HYP-3043, HYP-3042, HYP-3041, HYP-3040, HYP-3039, HYP-3026, HYP-3024, HYP-3023, HYP-3022, HYP-3021, HYP-3018, HYP-2963 |
| LTI-TODO-61 | Residual topology-exception collar lemma | LTI-192, LTI-189, LTI-188, LTI-187, LTI-184, LTI-183, LTI-177, LTI-176, LTI-171, LTT-090, LTT-087, LTT-086, LTT-085, LTT-082, LTT-081, LTT-075, HYP-3044, HYP-3041, HYP-3040, HYP-3039, HYP-3036, HYP-3035, HYP-3029, HYP-3028 |
| LTI-TODO-62 | Endpoint-owner transfer sidecar theorem | LTI-193, LTI-192, LTI-190, LTI-189, LTI-188, LTI-187, LTI-186, LTI-185, LTI-183, LTI-182, LTI-177, LTI-174, LTI-033, LTT-091, LTT-090, LTT-088, LTT-087, LTT-086, LTT-085, LTT-084, LTT-083, LTT-081, LTT-080, LTT-075, HYP-3045, HYP-3044, HYP-3042, HYP-3041, HYP-3040, HYP-3039, HYP-3038, HYP-3037, HYP-3035, HYP-3034, HYP-3029 | Add owner-strip and owner-transfer fields to the residual packet ledger, test the full `B18Z6` residual surface, and prove that owner names are retained, reconstructible, dual-annihilated, or routed to named residual debt. |
| LTI-TODO-63 | Hidden connector sidecar merger | LTI-194, LTI-193, LTI-192, LTI-191, LTI-190, LTI-189, LTI-188, LTI-187, LTI-186, LTI-185, LTI-184, LTI-183, LTI-182, LTI-180, LTI-179, LTI-174, LTI-169, LTI-166, LTI-162, LTI-061, LTT-092, LTT-091, LTT-090, LTT-089, LTT-088, LTT-087, LTT-086, LTT-085, LTT-084, LTT-083, LTT-082, LTT-081, LTT-080, LTT-078, LTT-077, LTT-072, HYP-3046, HYP-3045, HYP-3044, HYP-3043, HYP-3042, HYP-3041, HYP-3040, HYP-3039, HYP-3038, HYP-3037, HYP-3036, HYP-3035, HYP-3034, HYP-3032, HYP-3027, HYP-3022, HYP-3018, HYP-3013, HYP-2996, HYP-2886 | Merge the HYP-3046 sidecar list into packet manifests after preserving HYP-3045 endpoint-owner transfer fields; every hidden connector must state whether the forgotten coordinate is retained, reconstructed, annihilated, descended, or emitted as named residual debt. |
| LTI-TODO-64 | Observer-extension perspective carrier theorem | LTI-203, LTI-202, LTI-201, LTI-197, LTI-196, LTI-195, LTI-194, LTI-191, LTI-187, LTI-175, LTI-169, LTI-168, LTI-144, LTI-141, LTI-112, LTI-008, LTT-101, LTT-100, LTT-099, LTT-095, LTT-094, LTT-093, LTT-092, LTT-089, LTT-085, LTT-073, HYP-3056, HYP-3055, HYP-3054, HYP-3050, HYP-3049, HYP-3048, HYP-3047, HYP-3046, HYP-3043, HYP-3040, HYP-3039, HYP-2210, HYP-2120, HYP-1978, HYP-1977, THM-381, THM-385, THM-260, THM-409, T1138, T1137, T1136, T1132, T1131, T1130, T1129 | Define observer-extension cut perspectives as payload orbits under visible automorphisms, extend exact node/edge/triple/cycle/conflict carrier counts to `m=6`, and test whether ordered-pair sector decks, cross-sector orientation, and edge tail/tip sector words classify pair-good decoy generators and residual capacitor cuts before raw counts matter. Every new tournament-analysis script should declare whether vertices are nodes, ordered pairs, directed edges, directed cycles, proof obligations, payload columns, or conflict fibers before adding scalar sidecars. |
| LTI-TODO-65 | Rooted extension quotient theorem | LTI-198, LTI-197, LTI-196, LTI-195, LTT-096, LTT-095, LTT-094, LTT-093, HYP-3051, HYP-3050, HYP-3049, HYP-3048, HYP-3047, T1133, T1132, T1131, T1130, T1129 | Prove when a tournament or tiling quotient may pass from parent-rooted extension state to an unrooted class without losing LRC payload; every experiment should retain `parent_class`, `root_orbit`, `incident_word_orbit`, `layer_boundary_word`, `rank_one_sheet_id`, `rectangle_defect_rank`, and `unrooting_collision_fiber` until the lost coordinate is reconstructed, annihilated, descended, or named residual debt. |
| LTI-TODO-66 | Diagonal transport deletion-fiber theorem | LTI-199, LTI-198, LTI-197, LTI-196, LTI-195, LTT-097, LTT-096, LTT-095, LTT-094, LTT-093, HYP-3052, HYP-3051, HYP-3050, HYP-3049, HYP-3048, HYP-3047, T1134, T1133, T1132, T1131, T1130, T1129 | Prove that diagonal word-orbit transport plus deletion-parent fiber data is the correct recursive carrier for tournament/tiling growth; before using a half-tiling or unrooted A000568 class, retain `diagonal_word_orbit`, `K_position_line_profile`, `aligned_pair_counts`, `newest_link_bit`, `cross_sector_orientation_word`, and `deletion_parent_profile`, or certify the loss by reconstruction, dual annihilation, descent, or named residual debt. |
| LTI-TODO-67 | Diagonal-layer flow sidecar theorem | LTI-203, LTI-202, LTI-201, LTI-200, LTI-199, LTI-198, LTI-197, LTI-196, LTI-195, LTI-194, LTI-187, LTT-101, LTT-100, LTT-099, LTT-098, LTT-097, LTT-096, LTT-095, LTT-094, LTT-093, HYP-3056, HYP-3055, HYP-3054, HYP-3053, HYP-3052, HYP-3051, HYP-3050, HYP-3049, HYP-3048, HYP-3047, HYP-3043, HYP-3039, THM-381, THM-385, T1138, T1137, T1136, T1135, T1134, T1133, T1132, T1131, T1130 | Emit rectangle and hourglass cycle bases for adjacent diagonal-layer bridges, attach line-potential and cycle-defect words to owner/barcode/route packet experiments, and prove when the `K_{k,k+1}` coboundary redundancy descends, reconstructs, annihilates, or names a hidden LRC sidecar coordinate. |
| LTI-TODO-68 | HYP-2963 observer-cut orbit ledger | LTI-203, LTI-202, LTI-201, LTI-196, LTI-187, LTT-101, LTT-100, LTT-099, LTT-094, HYP-3056, HYP-3055, HYP-3054, HYP-3048, HYP-3039, HYP-2963, T1138, T1137, T1136 | Build the ledger schema from HYP-3056 over HYP-2963 coarse fibers: enumerate admissible observers, visible automorphism groups, payload orbit ids, changed LRC predicates, separating sidecars, and discharge modes; then emit the induced payload-column tournament and flag any noncommuting discharge cycles. |
| LTI-TODO-69 | Value-origin sidecar theorem | LTI-204, LTI-201, LTI-200, LTI-199, LTI-198, LTI-197, LTI-196, LTI-195, LTT-102, LTT-099, LTT-098, LTT-097, LTT-096, LTT-095, LTT-094, LTT-093, HYP-3057, HYP-3054, HYP-3053, HYP-3052, HYP-3051, HYP-3050, HYP-3049, HYP-3048, HYP-3047, HYP-3039, HYP-2991, HYP-2989, T1139, T1136, T1135, T1134, T1133, T1132, T1131 | Add `value_origin_type` and `lost_coordinate_exit` to every tournament/LRC quotient experiment that cites a small count; prove that the origin tag plus sidecars separate, reconstruct, annihilate, descend, or name debt for all route/status-changing fibers. |
| LTI-TODO-70 | Hyperbolic reciprocal sidecar audit | LTI-205, LTI-202, LTI-201, LTI-159, LTI-076, LTT-103, LTT-100, LTT-099, HYP-3058, HYP-3055, HYP-3054, HYP-3009, HYP-2963, T1140, T1137, T1136 | Define honest triples `(a,b,c)` on HYP-2963 packet families from exact scale, route incidence, automaton depth, observer-cut depth, primitive-period deck, Fejer degree, or state-lift obligation; test whether spherical/euclidean/hyperbolic sign predicts discharge route without dropping exact `M`, endpoint-owner, topology, deletion-fiber, rectangle/hourglass, or certificate payloads. |
| LTI-TODO-71 | Desargues/Beal finalizer ledger | LTI-207, LTI-204, LTI-203, LTI-201, LTI-200, LTI-187, LTI-185, LTI-182, LTI-179, LTT-105, LTT-102, LTT-101, LTT-099, LTT-098, LTT-083, LTT-080, LTT-077, HYP-3060, HYP-3057, HYP-3056, HYP-3054, HYP-3053, HYP-3037, HYP-3034, HYP-3031, HYP-2991, HYP-2963, THM-572, T1142 | Insert `desargues_girth6_residue` after rectangle/hourglass checks and `beal_common_owner_gate` after three-channel route collisions; prove each remaining residual either has zero gate and discharges through existing sidecars, descends through a girth-six incidence family, shares a common owner/factor coordinate, or becomes named F7/THM-572 debt. |

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

## Detailed S166 Technique Bank

The compact `LTI-*` atlas above is the canonical registry.  The bank below is
the long-form creative technique menu from the parallel S166 pass.  Future
agents can either promote one of these entries into a new `LTI-*` record or use
it directly as a reminder of what labels, observables, and quotient guardrails
must survive.

**Purpose:** A shared pull-from / contribute-to index of techniques for LRC14
work, especially techniques imported from tournament analysis, metagraphs,
sequence shadows, packet carriers, harmonic certificates, and cross-problem
analogies.

**Use this file when starting an LRC session.**  Pick a technique family, keep
its required side labels attached, and add a new entry when a session discovers
a reusable carrier, quotient guardrail, theorem template, or failed route that
should save future agents time.

## Contribution Template

Add new entries in the same shape:

```text
### Technique Name

Sources:
Preserves:
Destroys if used alone:
Best vertex sets:
Pairwise observable / gauge:
LRC use:
Next contribution:
```

For LRC or Tournament Analysis work, do not assume vertices are runners.  First
try at least a few of:

```text
runners, gaps, fixed circle sections, section boundaries, wall-crossing events,
endpoint owners, residues, cover arcs, Fourier modes, Haar rectangles,
Ramanujan exact-period packets, Fejer atoms, packet families, metagraph nodes,
state-lift obligations, matroid circuits, proof obligations.
```

## Current Proof-Use Rule

The strongest common lesson across the repo is:

```text
Keep labels until they have done their proof job.
```

A quotient is allowed to forget a coordinate only when:

```text
1. the target LRC predicate is constant on quotient fibers;
2. the forgotten coordinate is reconstructible from retained labels;
3. a dual/orthogonality certificate annihilates the forgotten direction; or
4. the forgotten direction is routed to a named residual sector.
```

This is the HYP-2978 / HYP-2990 controlled-kernel rule.  Use scalar invariants
only after this check.

## Fast Route Map

| If the session wants... | Start with... | Keep attached |
|---|---|---|
| a near-term LRC14 proof route | source-kernel / labelled-packet zipper | qdiv, exact `M`, Haar-open/boundary, endpoint owners, C27/K33, Fejer/Ramanujan handoff |
| a finite classifier | HYP-2963 packet families plus HYP-2956 F0-F7 | route, family, exact safe mass, q-threshold, state-lift debt |
| a computational pruning engine | observer-source reachability, good-cut/SCC, threshold quotient stack | marked observer, source target, SCC defect, endpoint-owner labels |
| a harmonic certificate | Fejer/Toeplitz, Ramanujan, boundary moments, Haar products | packet key, exact period, center, degree, atom bank, chart labels |
| a metagraph route | transfer matrices / operation metagraph / deck derivatives | predecessor class, extension strip, deleted-card derivative, torsion label |
| a weird cross-topic route | unit-distance ears, octahedral currents, OCF coimages, sequence shadows | endpoint sections, curl/divergence labels, compatibility relation, companion sequences |

## Technique Bank

### 1. Qdiv And Direct Rational Witnesses

Sources: THM-366, THM-523, HYP-2980.

Preserves: exact direct lonely witness at a denominator `q <= 14`.

Destroys if used alone: AP/GW boundary structure, off-apex escape size, endpoint
owner labels.

Best vertex sets: denominators, q-witness obligations, endpoint owners.

Pairwise observable / gauge: which packet has the smaller first missing
denominator, or which exact denominator supplies a strict witness.

LRC use: immediate discharge for rows that omit a multiple of some
`q in {2,...,14}`.

Next contribution: standardize qdiv as the first field in every new packet
record.

### 2. Endpoint-Owner Protection Cycles

Sources: THM-358, THM-360, THM-365, THM-379, HYP-2970, HYP-2975.

Preserves: which speed protects which threshold endpoint and whether endpoint
credit winds or cancels.

Destroys if used alone: continuous safe interval geometry and late exact-period
behavior.

Best vertex sets: endpoint owners, boundary points, taut bridge vertices, owner
pairs.

Pairwise observable / gauge: endpoint credit, owner pair sum, bridge sign, or
taut/positive bridge relation.

LRC use: separates AP/GW equality atoms from positive-open rows; feeds
state-lift and boundary-cocircuit exits.

Next contribution: attach endpoint-owner pairs to every Fejer/Ramanujan packet
certificate.

### 3. Haar/Baire Open-Boundary Split

Sources: HYP-2948, HYP-2949, HYP-2951, HYP-2984.

Preserves: strict-open safe mass versus closed threshold boundary support.

Destroys if used alone: hidden C27/GW transfer, endpoint owner identity.

Best vertex sets: open components, boundary cocircuits, safe interval fronts.

Pairwise observable / gauge: strict Haar mass, Baire interior, boundary-only
support, smallest open component.

LRC use: AP/GW are zero-open equality atoms; non-AP/GW tested packets have
positive strict-open mass.

Next contribution: prove boundary-owner skeleton rigidity beyond the current
finite AP-neighborhood.

### 4. Tope-Wall And Boundary-Cocircuit Arrangement

Sources: HYP-2986, HYP-2984, HYP-2951.

Preserves: the exact endpoint arrangement cut by `(14k +/- 1)/(14v)`.

Destroys if used alone: smoothing policy and formal dual certificate payload.

Best vertex sets: topes, boundary cocircuits, wall crossings, owner-zero walls.

Pairwise observable / gauge: existence of all-safe open tope, all-safe boundary
cocircuit, or no-tope/no-cocircuit wall candidate.

LRC use: turns "covered or not" into an oriented-matroid style proof object.

Next contribution: couple the tope/cocircuit output to Fejer interval family
templates.

### 5. Fejer/Toeplitz PSD Dual Certificates

Sources: HYP-2974, HYP-2981, HYP-2987.

Preserves: a dual certificate that the danger multiplicity cannot cover.

Destroys if used alone: packet route, exact safe component, endpoint-owner
labels, and interval rigor.

Best vertex sets: Fejer atoms, Fourier modes, Toeplitz moment rows, packet
fibers.

Pairwise observable / gauge: certificate degree, interval upper bound, atom-bank
size, negative margin.

LRC use: selected hard rows such as K33 `12->36` and `P10+GW` are
Fejer-certified; AP/GW are PSD-blind equality atoms.

Next contribution: compress named-row certificates into packet-family
templates with formal interval payloads.

### 6. Ramanujan Exact-Period Projectors

Sources: HYP-2978, HYP-2979, HYP-2981.

Preserves: primitive exact-period phase content via `c_q(n)`.

Destroys if used alone: endpoint ownership and strict-open safe measure.

Best vertex sets: exact periods, primitive unit residues, Ramanujan autocorrelation
packets.

Pairwise observable / gauge: weak/strict primitive witness by `q`, `c_q(a-b)`,
or Ramanujan danger-energy defect.

LRC use: pre-splits late denominators and explains why scalar divisor channels
mix routes.

Next contribution: add endpoint-owner Ramanujan profiles for q=25,27,41,4312.

### 7. Twist-Ladder And Blocker Hypergraphs

Sources: HYP-2972, THM-566, HYP-2901.

Preserves: rational twist witnesses and the finite hypergraph of failed twists.

Destroys if used alone: continuous endpoint geometry and adaptive denominator
growth.

Best vertex sets: twists `a/q`, failed twist blockers, denominator-ladder rungs.

Pairwise observable / gauge: first successful q, blocker inclusion, dynamic
ladder extension.

LRC use: q<=42 certifies the HYP-2963 bank, while THM-566 blocks any fixed
finite-ladder theorem.

Next contribution: turn failed ladders into recursive packet reducers rather
than treating them as search failure.

### 8. Boundary-Moment Multi-Chart Ledger

Sources: HYP-2969, HYP-2954, HYP-2987.

Preserves: missed-depth sector vectors and exact-period chart labels.

Destroys if used alone: continuous torus geometry between selected charts.

Best vertex sets: packet sectors, chart denominators, missed-depth histograms,
moment obligations.

Pairwise observable / gauge: which layer preserves strict-counterexample status
and owner labels before scalarization.

LRC use: warns that one all-covered denominator chart is not a contradiction;
the bridge must be multi-chart or feasible-region.

Next contribution: replace finite chart proxies by a true gK8/L_y feasible
region theorem.

### 9. Delsarte / gK8 / Linear Programming Moments

Sources: HYP-2823, HYP-2969, Lean `LRCFactorialAtom`.

Preserves: finite binding-row moment inequalities and explicit dual vectors.

Destroys if used alone: packet ownership and chart identity.

Best vertex sets: binding rows, moment coordinates, dual LP constraints.

Pairwise observable / gauge: stronger bound, smaller residual chart, or
compatible dual vector.

LRC use: useful terminal inequality layer after packet labels are fixed.

Next contribution: connect boundary-moment sectors to the Lean gK8 all-row
checks.

### 10. Source-Kernel Exposure Posets

Sources: HYP-2988 exposure-poset work, HYP-2980, HYP-2987.

Preserves: whether a labelled packet has any remaining unexposed source kernel
after all known certificates.

Destroys if used alone: family compression and formal interval details.

Best vertex sets: certificate channels, exposure obligations, source kernels.

Pairwise observable / gauge: which certificate channel exposes a row first.

LRC use: turns a counterexample into a no-hidden-kernel object rather than a
raw row.

Next contribution: define the first genuine `UNEXPOSED_SOURCE_KERNEL` row as
the F7 residual if it appears.

### 11. Zipper / Controlled-Kernel Handoff

Sources: HYP-2990, HYP-2987, HYP-2978.

Preserves: compatibility of two certificate carriers over a shared labelled
packet fiber.

Destroys if used alone: nothing by design, if the kernel checklist is obeyed.

Best vertex sets: proof obligations, packet fibers, residual sectors, carrier
interfaces.

Pairwise observable / gauge: predicate retention, fiber labels, declared
kernel, formal checkability, anti-scalar guardrail.

LRC use: global glue theorem target for Fejer, Ramanujan, endpoint, tope,
moment, smoothing, and state-lift certificates.

Next contribution: make F7 a named harmonic/state-lift sector rather than an
anonymous failure bucket.

### 12. Haar Product / Fixed-Margin Switches

Sources: HYP-2988 Haar tile, HYP-2989 Haar-product discrepancy, HYP-2990.

Preserves: mixed product signs erased by row/column margins.

Destroys if used alone: endpoint owners and packet family identity.

Best vertex sets: Haar rectangles, dyadic children, fixed-margin switch squares,
cocircuit signs.

Pairwise observable / gauge: mixed Haar coefficient, switch sign, row/column
margin collision.

LRC use: imports discrepancy/fixed-margin switch structure into packet
tilings; explains why margins alone are unsafe.

Next contribution: count independent color-compatible mixed switches for
structured LRC banks.

### 13. Admissible Smoothing Dispatcher

Sources: HYP-2985, HYP-2984, HYP-2982, HYP-2983.

Preserves: which proof clock is allowed to certify which packet.

Destroys if used alone: endpoint-owner clock or exact-period clock if a policy
is applied globally.

Best vertex sets: smoothing policies, packet exits, kernel clocks, approach
classes.

Pairwise observable / gauge: label retention across endpoint, exact-period,
smoothing, and far-approach clocks.

LRC use: AP/GW route to endpoint/Kaczynski, K33/covering to Fejer, q=27 petals
to Fejer plus Ramanujan, wide packets to Abel/Freiman/state-lift exits.

Next contribution: prove the admissible-smoothing lemma over HYP-2963 fibers.

### 14. Analytic Sieve / Kaczynski Boundary Packets

Sources: HYP-2982, HYP-2983, HYP-2679.

Preserves: smoothing kernel, exceptional-set boundary, major/minor arc split,
and approach class.

Destroys if used alone: prime-power exact-period packets and endpoint owners.

Best vertex sets: analytic proof modules, smoothing kernels, exceptional
approach classes.

Pairwise observable / gauge: labelled LRC packet versus large-sieve/Selberg
weight versus raw Mobius/prime count.

LRC use: late-denominator and true-wide packets can use analytic estimates only
as middle certificates.

Next contribution: formalize a Kaczynski/Abel off-resonance relation-height
bound.

### 15. Farey Branch And Unit-Excess Scheduling

Sources: HYP-2930, HYP-2931, HYP-2932, HYP-2934, HYP-2984 scheduler.

Preserves: exact optimum `M=p/q`, Farey excess `e=14p-q`, and route index `p`
on the unit-excess branch.

Destroys if used alone: endpoint geometry and state-lift labels.

Best vertex sets: exact fractions, Farey pairs, branch obligations.

Pairwise observable / gauge: theorem safety of exact scale, affine sum
mutation, product scheduler, or power stress test.

LRC use: after `e=1`, `p=1` is q-parent/AP-GW, `p=2` is C27 petal/two-block,
`p>=3` is K33/state-lift/Fejer.

Next contribution: prove every remaining q=14 non-AP/GW atom enters C27 or
K33 packet lanes.

### 16. Complete-Bipartite Product / K33 Incidence Wall

Sources: HYP-2932, HYP-2934, HYP-2945.

Preserves: the incidence depth of the Farey product `p*q` as `K_{p,q}`.

Destroys if used alone: denominator shell, summand shell, minor label.

Best vertex sets: factor fibers, complete-bipartite sides, K33 owner packets.

Pairwise observable / gauge: `p=2` planar two-block versus `p>=3` K33
incidence.

LRC use: isolates `3/41` / `12->36` as the first three-owner incidence wall.

Next contribution: convert K33 incidence owner packets into an explicit
THM-572 state-lift construction.

### 17. C27 / Unital / Goddyn-Wong Hidden Transfer

Sources: HYP-2942, HYP-2937, HYP-2920.

Preserves: q=3 unital branch labels, AP/GW calibration, petal/two-block route.

Destroys if used alone: global realizability and K33 branch distinction.

Best vertex sets: C27 shell slots, pair-completion blocks, branch-local
unital packets.

Pairwise observable / gauge: unique local completion, AP/GW transfer, or
petal/K33 route.

LRC use: explains why GW is invisible to Haar/Baire boundary owners yet
arithmetically distinct from AP.

Next contribution: lift branch-local completions into a global state-lift or
petal discharge theorem.

### 18. Summand / Multiplicand Bigraded Signatures

Sources: HYP-2934, HYP-2161, HYP-2083.

Preserves: denominator shell, Farey additive pair, factor-fiber ledger.

Destroys if used alone: exact LRC safe interval data.

Best vertex sets: summand nodes, factor fibers, denominator shells.

Pairwise observable / gauge: branch retention across additive shell and
multiplicand product.

LRC use: separates `2/27` as C=27 second-gap shell from `3/41` as K33 blow-up.

Next contribution: use the C=27 summand-unit toolkit to close the p=2 branch.

### 19. Perfect-Product / Abundancy Calibration

Sources: HYP-3013, HYP-2945, HYP-2220, HYP-2221.

Preserves: exact unit-excess product chain as an arithmetic control family.

Destroys if used alone: LRC14 denominator and graph-minor labels.

Best vertex sets: perfect-product fractions, divisor abundance shadows,
incidence walls.

Pairwise observable / gauge: perfect fixed point at n=2 versus deficient
LRC14 shadow.

LRC use: calibrates product-chain intuition without treating product as a
proof scalar.

Next contribution: compare LRC14 power-of-two/prime-q shadows against K33 rank
and divisor-lattice features.

S174 addendum: HYP-3013 computes the comparison directly.  Prime `q=14a-1`
rows at `a=1,16,256` have defect `12/q`, while composite `q14` rows in the
scan are abundant, starting with `q=27=3^3`.  This promotes the perfect-product
lane from a reminder to a packet-field requirement: keep `prime_q_flag`,
factorization, abundancy defect, product/Kpq route, and automaton state.

### 20. Divisor / Multiplicative Function Guardrails

Sources: HYP-2978, divisor/Ramanujan quotient audit.

Preserves: controlled arithmetic packet laws such as Dirichlet convolution,
Jordan capacity, and Ramanujan primitive traces.

Destroys if used alone: route labels; scalar divisor signatures mix AP/GW,
K33, petal, and covering rows.

Best vertex sets: divisor fibers, exact-period packets, unitary divisor
features.

Pairwise observable / gauge: which quotient avoids route-mixing collisions.

LRC use: falsifies unsafe scalar quotients and supplies admissibility rules.

Next contribution: extend quotient-collision tests to the full HYP-2963 bank
with endpoint-owner labels.

### 21. Singular-Series Relation Lattices

Sources: THM-501, THM-515, THM-518, HYP-2980.

Preserves: exact additive relation lattice and archimedean sinc-weighted
covering-deficit density.

Destroys if used alone: closed-threshold endpoint debt and finite resonant
stranger behavior.

Best vertex sets: relation packets, additive lattices, stranger tails.

Pairwise observable / gauge: relation-lattice suppression strength, finite
resonance, AP-like energy.

LRC use: explains why AP-like cores are hard and why generic/Sidon-like rows
are easy.

Next contribution: control finite resonant strangers around the AP-tail
families.

### 22. Riesz / Prime Point-Mass Measures

Sources: OPEN-Q-104, THM-518.

Preserves: point-mass certificates and relation-tuned Fourier weights.

Destroys if used alone: exact finite resonance; Riesz products fail AP
extremizers.

Best vertex sets: point masses, Fourier relations, AP-tail stranger packets.

Pairwise observable / gauge: certificate ratio to `1/14`, relation-count
gain, finite resonance gap.

LRC use: near-miss analytic route; prime point mass reaches `2/29`, not enough
but informative.

Next contribution: test composite/relation-tuned point masses against AP-tail
finite resonant strangers.

### 23. P0 / Bonferroni / Witness-Floor Route

Sources: HYP-2832, Lean `LRCWitnessBonferroni`, `LRCDenseCovers`.

Preserves: event-level lower bound from safe set and cover-set complement.

Destroys if used alone: final finite-ruler Part A realization and exact
goodSet bridge if not attached.

Best vertex sets: events (`GOOD`, `safeSet`, `coverSet`, `denseSet`), proof
obligations.

Pairwise observable / gauge: inclusion, Bonferroni handoff, p0 margin, cap
duality.

LRC use: converts wide p0 bounds into positive witness floors.

Next contribution: keep this formal route aligned with the packet classifier
rather than duplicating analytic floors.

### 24. Finite Part A / GoodSet Formal Bridge

Sources: THM-527, `LRCGoodSet`, `LRCGapReach`, `LRCWitnessPartA`.

Preserves: implication from a finite phase gap to a lonely margin.

Destroys if used alone: analytic floor establishing that the gap event has
positive measure.

Best vertex sets: phase gaps, arcs, finite-ruler approximations, Part A
obligations.

Pairwise observable / gauge: arc-count error budget, positive finite witness,
Mreach handoff.

LRC use: formal end of the p0/Bonferroni route.

Next contribution: connect packet-family safe components to the existing
goodSet witness-margin theorems.

### 25. Observer-Source Marked Reachability

Sources: THM-381, S511, S583.

Preserves: exact LRC predicate in a marked tournament: observer lonely iff the
marked observer is a source.

Destroys if used alone: continuous time geometry after cell projection.

Best vertex sets: marked observers, source targets, reachability cells,
threshold-decorated classes.

Pairwise observable / gauge: source reachability, marked H-loneliness mask,
observer-source pressure.

LRC use: turns LRC into a finite marked-tournament target for speedup engines.

Next contribution: precompute marked-class transition automata for residual
packet families.

### 26. Apex / Phase Winding Tournament Atlases

Sources: HYP-2924, HYP-2923, HYP-2913.

Preserves: half-clock order class and denominator-14 phase winding.

Destroys if used alone: q-threshold divisibility and off-apex escape magnitude.

Best vertex sets: apex points on `Z/14Z`, terminal runner phases, wall crossings.

Pairwise observable / gauge: clockwise displacement, half-clock cutoff, tie
Hamiltonian path by residue or speed.

LRC use: exact finite atlas for AP/GW/near-miss phase classes; magnitude-blind
guardrail is explicit.

Next contribution: use phase classes only as state-lift interfaces with q data
retained.

### 27. Packet-Sign K4 / Six-Sector Tournaments

Sources: HYP-2677, HYP-2676, THM-546.

Preserves: chamber topology of one-missed-sector packet signs and weighted
opposite-pair balances.

Destroys if used alone: exact magnitudes and additive profile.

Best vertex sets: six missed sectors, K4 edges, opposite edge-pair balances.

Pairwise observable / gauge: cyclic pair-sum, opposite-compensated pressure,
K4 edge sign orientation.

LRC use: separates same-sign pockets, KPS third-pocket cancellation, boundary
collars, and true-wide branches.

Next contribution: attach packet-sign fingerprints to HYP-2963 route rows.

### 28. Good-Cut / SCC Support Residue

Sources: THM-354, INV-237, HYP-2990.

Preserves: support-residue coordinate `g_P(T)=n-#SCC(T)`.

Destroys if used alone: arithmetic labels on protector intervals.

Best vertex sets: good cuts, SCC condensation components, support residues.

Pairwise observable / gauge: SCC defect, top-bucket strong connectivity,
component-boundary transport.

LRC use: possible state-lift coordinate for endpoint-protection graphs and
quotient transport.

Next contribution: attach SCC/private-leaf peeling to endpoint-owner pressure
graphs.

### 29. Metagraph Transfer Matrix Chain

Sources: `metagraph-as-transfer-chain.md`, A000568 metagraph work.

Preserves: incremental class extension from `n` to `n+1` via strip variables.

Destroys if used alone: arithmetic endpoint labels and LRC threshold data.

Best vertex sets: metagraph classes, transfer-matrix states, extension strips.

Pairwise observable / gauge: successor fan-out, H-range transfer, hereditary
maximizer, strip tile contribution.

LRC use: model LRC recursion as a transfer chain rather than a single invariant.

Next contribution: build an LRC packet transfer matrix with qdiv/Farey labels.

### 30. Operation-Metagraph State Vectors

Sources: HYP-1835, S378.

Preserves: repair deficit, torsion order, endpoint core, `phi(n)`, `tau(n)`,
and operation closure across denominator changes.

Destroys if used alone: exact safe interval endpoints.

Best vertex sets: residue vectors, one-coordinate moves, endpoint feature
states, operation closures.

Pairwise observable / gauge: repair balance, selected torsion order, exposed
cell creation.

LRC use: predicts prime/composite seams and n=14 torsion leakage better than
raw missed-cell count.

Next contribution: run n=14 support-shell metagraphs ordered by repair deficit
and torsion labels.

### 31. Burnside / Restricted Orbit Counts

Sources: HYP-2097, HYP-2099, S583, A000568 quotient scripts.

Preserves: symmetry quotient representatives and exact orbit counts.

Destroys if used alone: endpoint/pair certificates unless reattached.

Best vertex sets: orbits, fixed classes, round/merged/fixed-boundary classes.

Pairwise observable / gauge: orbit size, stabilizer, converse/complement
transport, fixed-boundary purity.

LRC use: compresses the `n=14` search surface from raw rows toward fixed
round/certified classes.

Next contribution: refine round316/merged190/fixed64 with endpoint-owner
labels.

### 32. A^2 Isomorphism Fingerprint Lookup

Sources: S339, S583, fractal-isomorphism work.

Preserves: fast tournament isomorphism fingerprints through verified ranges.

Destroys if used alone: proof certainty beyond range unless canonical fallback
exists.

Best vertex sets: adjacency-square row fingerprints, quotient-class caches.

Pairwise observable / gauge: sorted `A^2` rows, collision status, canonical
fallback hit.

LRC use: speed up repeated quotient walks and state-lift searches.

Next contribution: add A^2 cache keys to LRC packet JSON fingerprints.

### 33. Cheap TDA / Feature Prefilters

Sources: `tournament_tda.py`, S355, S583.

Preserves: score, c3, SCC defect, residue rank, Omega packet summaries.

Destroys if used alone: LRC threshold predicate.

Best vertex sets: quotient cells, TDA features, residue summaries.

Pairwise observable / gauge: feature dominance, residual rank, deletion residue.

LRC use: ranks residuals for exact search; never certifies alone.

Next contribution: emit standard TDA feature blocks for every HYP-2963 packet.

### 34. Near-Transitive Perturbation DP

Sources: S90 applications, S583.

Preserves: local updates around AP/Vstar/tie-wall neighborhoods with small
upset sets.

Destroys if used alone: global class exactness when upset set grows.

Best vertex sets: upset arcs, perturbation neighborhoods, tie-wall flips.

Pairwise observable / gauge: upset-set inclusion, local update cost, AP
distance.

LRC use: efficient finite audits near AP/GW and near-miss rows.

Next contribution: parameterize packet families by upset count and exact
route labels.

### 35. Tournament Deck Derivatives

Sources: HYP-2237, S660 deck atlas.

Preserves: paired derivative `(card type, deleted vertex outdegree)` rather
than scalar deck data.

Destroys if used alone: LRC carry/owner derivative unless paired to the card.

Best vertex sets: deleted cards, carry coordinates, owner derivatives.

Pairwise observable / gauge: which paired derivative resolves quotient
collisions.

LRC use: motivates carry-deletion derivative theorem for the Res_27 / C=27 seam.

Next contribution: define LRC card + carry/owner derivative records for floor
atoms.

### 36. Carry-Deletion Derivative

Sources: HYP-2237, HYP-2167, HYP-2230.

Preserves: exact changes in `M`, odd wall pairs, gcd-shell mass, n-clock word,
and owner route under local deletion/toggle.

Destroys if used alone: global lift/CRT seam.

Best vertex sets: carry cocycles, deletion derivatives, owner route cards.

Pairwise observable / gauge: zero-tax derivative versus strict derivative.

LRC use: possible no-leak theorem around C=27 / odd-wall seam.

Next contribution: compute derivative vectors for all HYP-2164 floor/near-floor
atoms.

### 37. Owner-Concurrency Jackknife

Sources: HYP-2237, finite-field Kakeya/Falconer prompts.

Preserves: pinned owner fibers and hidden packets deciding floor versus strict.

Destroys if used alone: exact endpoint geometry if owner labels are not tied
to intervals.

Best vertex sets: owner fibers, pinned directions, concurrency packets.

Pairwise observable / gauge: which owner packet decides certificate loss under
jackknife deletion.

LRC use: search for decisive hidden owner channels after visible supports
saturate.

Next contribution: hold odd-wall/C27 support fixed and vary owner/carry packets.

### 38. Apex Sheaf Gluing

Sources: HYP-2099, HYP-2100, HYP-2101, HYP-2237.

Preserves: local cheap-pair certificates as sections over unit/nonunit/apex
charts.

Destroys if used alone: global restriction residuals unless charts are glued.

Best vertex sets: sheaf charts, local sections, restriction maps, carry moves.

Pairwise observable / gauge: section compatibility, failed gluing residual,
endpoint-owner positivity.

LRC use: turns local cheap certificates into a finite sheaf-gluing problem.

Next contribution: include carry derivatives as restriction maps in the apex
sheaf.

### 39. Two-Block Helly Extractor

Sources: HYP-2237, S599, HYP-2107.

Preserves: small determinant/CRT contradictions in multiple/Cprime residuals.

Destroys if used alone: full owner packet context.

Best vertex sets: component languages, pair constraints, CRT residue states.

Pairwise observable / gauge: minimal incompatible component certificate.

LRC use: discharge residual rows after cheap-pair and local carry derivatives
fail.

Next contribution: build a library of two-block Helly contradictions keyed by
packet route.

### 40. Three-State L/M/R Middle Automata

Sources: HYP-2109, S582, S583.

Preserves: left/right owner exits and middle wall survival.

Destroys if used alone: metric margin unless middle cells keep length/residue
labels.

Best vertex sets: endpoint-owner cells, residue cells, wall-crossing events,
proof obligations.

Pairwise observable / gauge: event words over L/M/R and closed-middle circuits.

LRC use: encodes "side cannot flip without middle witness" for owner/carry
obligations.

Next contribution: compile HYP-2963 route rows into L/M/R event words.

### 41. Threshold-Decorated Quotient Stack

Sources: S535, S539, S583.

Preserves: good-only class fibers after observer/source/gap/threshold labels.

Destroys if used alone: purity if labels are dropped.

Best vertex sets: threshold-colored classes, section boundaries, quotient
fibers.

Pairwise observable / gauge: good-only purity, threshold label retention,
section boundary behavior.

LRC use: pre-sieve quotient cells before exact interval search.

Next contribution: define a standard threshold-color schema for packet-bank
rows.

### 42. Support-Residue Calculus

Sources: `support-residue-calculus.md`, THM-354, THM-025/H=63 work.

Preserves: surviving support geometry after projection/cancellation.

Destroys if used alone: exact arithmetic coordinates.

Best vertex sets: supports, projection residues, deletion residues, Omega
packets.

Pairwise observable / gauge: which projection kills supports and what residue
independence profile remains.

LRC use: asks what residue survives after qdiv/Farey/Haar/Fejer projections.

Next contribution: define an LRC source-kernel residue rank analogous to
Omega deletion-residue rank.

### 43. OCF Hard-Core Activity / Coimage

Sources: HYP-2618, HYP-2936, HYP-2990.

Preserves: compatible packet address and activity-two partition evaluation.

Destroys if used alone: signed/compatible distinction if read as raw noise
stability.

Best vertex sets: odd-cycle packets, compatibility classes, coimage fibers.

Pairwise observable / gauge: activity value, forbidden H-value signal,
compatibility-preserving map.

LRC use: terminal state-lift category after LRC residual is lifted into a
tournament conflict packet.

Next contribution: define the exact LRC packet coimage evaluation feeding
THM-572.

### 44. Anti-Poisson / Coimage Atlases

Sources: HYP-2156, HYP-2936.

Preserves: visible/hidden side channels and coimage maps under quotient.

Destroys if used alone: exact LRC threshold scale.

Best vertex sets: coimage fibers, anti-Poisson packets, visibility states.

Pairwise observable / gauge: which hidden coordinate becomes visible under
carrier refinement.

LRC use: language for residual packets that are invisible to scalar quotient
but visible after state lift.

Next contribution: connect anti-Poisson visibility to F7 residual definition.

### 45. Path Homology / Boundary-Image Residues

Sources: support-residue calculus, HYP-1815, path-homology threads.

Preserves: boundary image, old-cycle residue, rank/torsion of packet incidence.

Destroys if used alone: metric threshold information.

Best vertex sets: packet incidence matrices, boundaries, cokernels, old path
residues.

Pairwise observable / gauge: incidence rank, torsion, ghost versus genuine
cycle residue.

LRC use: possible formal language for source-kernel / root-packet residuals.

Next contribution: compute packet incidence rank for AP/GW, K33, petal, and
covering packet banks.

### 46. Root-Packet Incidence Rank

Sources: HYP-1815, HYP-2663, HYP-2664.

Preserves: transfer rank/torsion after support compatibility is fixed.

Destroys if used alone: shell/carry gate information.

Best vertex sets: root packets, packet columns, holes/tails, endpoint fibers.

Pairwise observable / gauge: rank collapse, private packet columns, support
versus incidence distinction.

LRC use: after shell-1 gate, classify shell-full AP-tail packets by root rank
and mouth owner.

Next contribution: build rank matrices for HYP-2664 shell-full residues.

### 47. Shell-1 Carry Gate / Mouth Owner

Sources: HYP-2661, HYP-2664, HYP-2990.

Preserves: whether the tower `{1,2,4,8}` survives and which old mouth owns the
residue.

Destroys if used alone: root-packet incidence rank.

Best vertex sets: shell bits, AP holes/tails, mouth-owner labels.

Pairwise observable / gauge: shell-damaged versus shell-full, old-mouth survival,
carry conservation.

LRC use: drastically reduces AP-tail finite residue burden before comb
enumeration.

Next contribution: prove the shell-1 deletion theorem independently of bounded
scans.

### 48. Octahedral Hodge Current

Sources: HYP-2887, HYP-2892, HYP-2990.

Preserves: divergence, face curl, and tail defects on the octahedral `L(K4)`
carrier.

Destroys if used alone: endpoint/Farey packet labels unless attached.

Best vertex sets: octahedral vertices, faces, currents, wall deletions.

Pairwise observable / gauge: divergence source, curl obstruction, Abel/Cauchy
tail.

LRC use: repeated-packet lift realizability and support-six current obstruction.

Next contribution: prove the octahedral Hodge estimate after low-height wall
deletion.

### 49. Clebsch / Half-Cube / Unital Design Frames

Sources: HYP-2891, HYP-2892, HYP-2942.

Preserves: covariance/frame structure in finite carrier quotients.

Destroys if used alone: exact route labels and sparse-design exact-tiler claims.

Best vertex sets: Clebsch classes, unital points/blocks, Bruhat faces.

Pairwise observable / gauge: pair-design incidence, curl frame, covariance
carrier.

LRC use: finite side-channel frames for residual leak, not scalar quotient
proofs.

Next contribution: align Clebsch/Bruhat faces with packet route labels in the
K33/C27 frontier.

### 50. Walsh / Fourier Degree Towers

Sources: support-residue calculus, THM-164, HYP-2977.

Preserves: which degree first sees a non-score cycle residue.

Destroys if used alone: endpoint ownership and exact M/Farey scale.

Best vertex sets: Fourier modes, Walsh degrees, cycle supports.

Pairwise observable / gauge: first separating degree, support tower level,
spectral shadow.

LRC use: explains why low-degree PSD/moment tests can miss AP/GW and where
Fejer regularization helps.

Next contribution: record the first separating harmonic degree for each
HYP-2963 route family.

### 51. Paley / Interval Gradient Amplifiers

Sources: THM-143, support-residue calculus, HYP-2936.

Preserves: co-occurrence gradients and multiplicity/disjointness residue.

Destroys if used alone: LRC-specific endpoint arithmetic.

Best vertex sets: cycle supports, distance shells, Fourier magnitudes.

Pairwise observable / gauge: gradient slope, alpha_j amplification,
disjointness excess.

LRC use: analogy for how a small structured residue can dominate higher
packet layers.

Next contribution: seek LRC analogues of interval-gradient amplifiers in
endpoint-owner pressure graphs.

### 52. Sequence-Shadow Weather Reports

Sources: S633, sequence-shadow recursions, A000568/SC/Burnside work.

Preserves: companion shadows of a hard count: fixed, merged, nonfixed,
q-shadow, bisection, skinny, transporter.

Destroys if used alone: theorem predicate unless the shadow is tagged.

Best vertex sets: sequence layers, bisections, transporter orbits.

Pairwise observable / gauge: which shadow changes when the proof problem
changes.

LRC use: build companion sequences for packet counts and route families instead
of chasing one hard next term.

Next contribution: maintain route-family companion sequences for HYP-2963
under qdiv, exact period, and state-lift labels.

### 53. Anti-Coset Transporters

Sources: HYP-2208, THM-409, S632.

Preserves: `Anti(X)=Iso(X,JX)` transporter records for complement/converse,
reflection, conjugation.

Destroys if used alone: threshold and endpoint labels.

Best vertex sets: rooted perspectives, automorphism orbits, anti-cosets.

Pairwise observable / gauge: transporter size, stabilizer size, orbit action.

LRC use: disciplined complement/reflection quotienting in round/self-converse
lanes.

Next contribution: add transporter records to fixed-boundary LRC quotient
tables.

### 54. Fixed-Margin Packet Classifiers

Sources: HYP-2962, HYP-2963, HYP-2956.

Preserves: qdiv, q-cover vector, U14 profiles, C27 states, exact M, Haar-open
mass, route.

Destroys if used alone: family proof templates; it classifies but does not
prove.

Best vertex sets: packet signatures, route families, fixed-margin fibers.

Pairwise observable / gauge: signature equality, route mixing, family/sporadic
split.

LRC use: current backbone for "all possible counterexamples into families and
sporadics."

Next contribution: attach Fejer/Ramanujan/Haar-product certificate IDs to every
fixed-margin packet signature.

### 55. NORK / Pinch Templates

Sources: HYP-2966, HYP-2964, HYP-2965.

Preserves: pinch-template local obstruction and moon-core/few-apex structure.

Destroys if used alone: global packet route if not tied to boundary-moment
ledger.

Best vertex sets: pinch sites, apex apertures, local packet obligations.

Pairwise observable / gauge: which pinch reduces a packet or creates boundary
debt.

LRC use: local structural reductions in covering / few-apex lanes.

Next contribution: add NORK pinch tags to HYP-2963 packet-family records.

### 56. Pincer / Pressure DAGs

Sources: HYP-2950, HYP-2952, HYP-2953.

Preserves: endpoint pressure, derived relative profiles, transitive apex
pressure fronts.

Destroys if used alone: C27/K33 hidden state labels.

Best vertex sets: pressure sources, endpoint cells, derived profile classes.

Pairwise observable / gauge: pressure dominance, source-spectrum pullback,
profile hash.

LRC use: routes qdiv=14 skeleton-gate subclasses and identifies where
state-lift labels are needed.

Next contribution: connect pressure DAG features to the no-hidden-kernel
exposure poset.

### 57. Cauldron / Residue Stream Games

Sources: HYP-2190, HYP-2936.

Preserves: correlated residue stream behavior that scalar density misses.

Destroys if used alone: exact LRC endpoint geometry.

Best vertex sets: residue streams, two-block game states, cauldron buckets.

Pairwise observable / gauge: stream correlation, block transition, residue
leak.

LRC use: stress-test density-style proofs against correlated residue packets.

Next contribution: apply cauldron-style stream diagnostics to p=2 C27 and
p>=3 K33 branches.

### 58. Unit-Distance Endpoint-Ear Recursion

Sources: HYP-2620, THM-408, HYP-2990.

Preserves: endpoint-compatible Hamiltonian unit spine through deletion/extension.

Destroys if used alone: geometric unit-cocyclic realizability.

Best vertex sets: endpoint masks, deletion ears, unit-cocyclic neighbor sets.

Pairwise observable / gauge: endpoint-universal deletion, traceability, spine
extension.

LRC use: cross-domain model for "graph-only certificate survives; geometry is
the real obstruction."

Next contribution: import endpoint-ear language into LRC state-lift handoffs:
which certificate endpoint remains extendable?

### 59. Moser / Unit-Distance Lattice Ladders

Sources: unit-distance toolkit leads, HYP-2298, THM-432.

Preserves: exact lattice unit vectors and angle-ladder fields.

Destroys if used alone: tournament/LRC packet predicates.

Best vertex sets: unit vectors, norm layers, CM-field rungs.

Pairwise observable / gauge: unit-vector capacity, endpoint universality,
lattice extension.

LRC use: analogy for retaining norm/field side channels before scalar edge
counting.

Next contribution: compare exact-period LRC packets with unit-distance
field/norm side-channel discipline.

### 60. Graph Minor / Kuratowski Guardrails

Sources: HYP-2932, HYP-2945, graph prompts.

Preserves: minor/subdivision closure and irreducible obstruction cores.

Destroys if used alone: Farey denominator and LRC packet labels.

Best vertex sets: minor cores, Kpq incidence sides, obstruction packets.

Pairwise observable / gauge: minor containment, K33 incidence wall, anti-scalar
edge-count guardrail.

LRC use: product-edge counts are safe only with the complete-bipartite and minor
labels retained.

Next contribution: define exactly what LRC state-lift obstruction corresponds
to "contains K33" in the packet category.

### 61. Exposure / Haar / Zipper Integration

Sources: HYP-2988, HYP-2989, HYP-2990.

Preserves: no-hidden-exposure, mixed Haar switch labels, and controlled-kernel
handoff at once.

Destroys if used alone: none if it uses the full controlled-kernel checklist.

Best vertex sets: exposure channels, Haar switches, packet zippers, residual
sectors.

Pairwise observable / gauge: which label survives the quotient and which
certificate consumes it.

LRC use: current "summit" integration lane after the S165 rebase.

Next contribution: define an integrated packet record with exposure channel,
Haar-switch signature, and zipper residual status.

### 62. Haar Zipper Cocycle

Sources: HYP-2991, HYP-2990, HYP-2989, HYP-2988, HYP-2595.

Preserves: the local mixed Haar / fixed-margin switch coordinate
`zeta(T)=T00-T01-T10+T11` on labelled `2 x 2` packet tables.

Destroys if used alone: global endpoint owners, exact-period labels, and
packet-family context; `zeta` is a local tooth, not the whole proof object.

Best vertex sets: local packet tables, Haar rectangles, fixed-margin fibers,
boundary cocircuits, and unpaired cocycle residuals.

Pairwise observable / gauge: whether two packets share margins but differ in
mixed cocycle, and whether the cocycle cancels, hands off, descends, or remains
as state-lift debt.

LRC use: make row/column margin quotients admissible only after `zeta` is
retained, reconstructed, annihilated by a discrepancy certificate, or routed to
named F7/THM-572 residual debt.

Next contribution: compute packet-level `zeta` signatures for the HYP-2963
bank and attach each nonzero owner-strip/cross/nested coefficient to a
HYP-2987 handoff arrow.

### 63. Haar Product Tiling Synthesis

Sources: HYP-2989, HYP-2988, HYP-2452, HYP-2450, THM-351, THM-346, THM-345.

Preserves: the product-rule classes behind Haar rectangles, fixed-margin
switches, and fixed-path staircase Walsh tilings.

Destroys if used alone: endpoint owners and packet-family labels; product rules
explain local propagation but do not certify an LRC packet without its retained
coordinates.

Best vertex sets: Haar rectangle products, staircase tiles, product-rule
classes, fixed-margin fibers, and discrepancy-compatible color packets.

Pairwise observable / gauge: whether the product of two local carriers is
orthogonal zero, same-tile atom, owner strip, cross handoff, nested descent, or
residual product debt.

LRC use: refuse raw component-count, row/column-margin, scalar-discrepancy, or
tournament-isomorphism quotients until product-rule descent is homogeneous,
owner coordinates are reconstructed, orthogonality/dual annihilates the lost
mode, or the mode is sent to a named residual.

Next contribution: instantiate the product-rule classes on HYP-2963 packet
families and compare independent color-compatible product modes against the
HYP-2595 `k+c_GP` discrepancy target.

### 64. Packet Cocycle Atlas

Sources: HYP-2995, HYP-2991, HYP-2990, HYP-2989, HYP-2978, HYP-2970,
HYP-2230, HYP-2241, HYP-2618.

Preserves: quotient fibers, lost-coordinate cochains, local closedness/gluing
laws, and the legal exits that make forgetting a coordinate theorem-safe.

Destroys if used alone: direct numerical proof content.  This is a proof
interface, not a certificate; each cocycle still needs exact values on packet
families and a verified exit.

Best vertex sets: packet fibers, proof carriers, endpoint arrows, Haar squares,
carry derivatives, exact-period modes, product-rule classes, transfer states,
octahedral currents, OCF coimage sectors, and residual proof obligations.

Pairwise observable / gauge: which carrier better retains the LRC predicate,
fiber labels, local defect, closedness law, coboundary/dual/boundary/descent
exit, formal checkability, and named residual routing.

LRC use: define a cochain `omega_Q` for every quotient `Q:P(S)->Y` measuring
the coordinate forgotten by `Q`; admit the quotient only when the predicate is
fiber-constant, the coordinate is reconstructed, `omega_Q` is exact, a dual
certificate annihilates it, it descends to a smaller packet family, it is
AP/Goddyn-Wong boundary equality, or it becomes named F7/THM-572 state-lift
debt.

Next contribution: implement an `omega_Q` packet schema over the HYP-2963 bank
and fill it for Haar `zeta`, endpoint-credit, Ramanujan, Fejer, Farey/K33,
C27/unital, carry-owner, product-rule, and boundary-moment quotients.

## Pull Lists

### Good First Contributions

- Add source/result links to any technique entry whose source list is too terse.
- Add a small "used by" note when a new session applies an entry.
- Add a route-family packet example for any technique that currently has only
  analogy evidence.
- Add a failed-use warning when a scalar quotient mixes proof routes.

### High-Leverage Contributions

- Build a standard HYP-2963 packet JSON schema containing qdiv, exact M/Farey,
  Haar/Baire state, endpoint owners, Fejer certificate ID, Ramanujan profile,
  Haar-product switch signature, exposure channel, state-lift debt, and zipper
  residual.
- Define `F7` as a concrete harmonic/state-lift residual sector.
- Compress Fejer certificates familywise for K33, petal/two-block, covering,
  and few-apex packets.
- Add endpoint-owner Ramanujan profiles to exact-period projectors.
- Build an LRC packet transfer matrix, inspired by the tournament metagraph
  chain, with qdiv/Farey labels retained.
- Implement the HYP-2995 `omega_Q` packet-cocycle schema and verify one legal
  exit for each quotient in the HYP-2963 bank.

### Red-Flag Patterns

- Raw tournament isomorphism class without qdiv/off-apex data.
- Product edge count without `K_{p,q}` and minor labels.
- Divisor scalar signature without exact-period and endpoint labels.
- One denominator chart treated as a global obstruction.
- Fejer floating value without packet key, center, degree, and interval backend.
- TDA feature used as a certificate rather than a residual ranker.
- Sequence next-term chasing without fixed/merged/nonfixed/q-shadow companions.

## Current Best LRC14 Assembly Sketch

```text
primitive row
  -> qdiv / direct witness gate
  -> exact M/Farey branch and unit-excess scheduler
  -> AP/GW boundary atom or positive Haar/Baire/tope exposure
  -> C27 p=2 shell or K33 p>=3 incidence packet
  -> labelled packet family classifier
  -> Fejer/Ramanujan/twist/boundary-moment/smoothing handoff
  -> exposure + Haar-switch + zipper controlled-kernel check
  -> strict certificate, AP/GW equality, or named F7 state-lift residual
```

The challenged assumption is that one object type, such as runners or
tournament isomorphism classes, should be the universal vertex set.  The index
above records the opposite lesson: choose vertices from the load-bearing proof
obligation, and state exactly what the quotient preserves and destroys.
