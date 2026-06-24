# Investigation Backlog

**Purpose:** Systematic catalog of every lead, reference, connection, and unexplored direction extracted from the repo. Claude agents should consult this before choosing what to work on, and add new leads as they emerge. Prioritized by potential impact on proving OCF (Claim A).

**Last full repo scour:** opus-2026-03-06-S10
**Last web research:** kind-pasteur-2026-03-10-S50 (Tang-Yau circulant path homology, Schweser-Stiebitz-Toft Redei revisited, Mitrovic NC Redei, Satake DRT)
**Last engineering update:** kind-pasteur-2026-03-10-S54 (sparse T_19 Omega computation, CLAUDE.md engineering mandate)
**Last packing framework update:** opus-2026-03-14-S71f (nesting obstruction, (z-2)(z-3) recurrence, 2-Bridge)
**Last S90 update:** opus-2026-03-15-S90 (simplicial Rédei, Cayley monad, τ-φ clock, equidecomposability)
**Last gauge theory analysis:** kind-pasteur-2026-03-21-S12 (Napolitano paper, Cartan decomposition bridge, TournamentProbe)
**Last literature sweep:** opus-2026-04-05-S24 (Stanley-Stembridge proved, Mitrovic NC deletion-contraction, Tang-Yau circulant Fourier)
**Last LRC web trawl:** mac-mini-2026-06-21-S20 (Rosenfeld 2025 proves LRC n=8,9,10; Huffer-Shepp Schur-convex coverage; Cusick AP-tight; see leads below)
**Last full repo scour:** opus-2026-03-05-S9
**Last web research:** opus-2026-03-05-S9 (Paley maximizer, n=8 anomaly)

---

## Lead codex-2026-06-24-S166: Haar zipper cocycle for LRC14

**Status:** PROOF-INTERFACE / local cocycle audit and theorem target complete (HYP-2991/T1075).  User asked to synthesize recent agent work through discrepancy theory, the `2D` Haar product rule, tournament tiling, and more zipper theorems.  This pass rebased over HYP-2990's abstract zipper atlas and makes the local no-free-slider obstruction explicit.  The script `04-computation/lrc14_haar_zipper_cocycle_codex_s166.py` stores output in `05-knowledge/results/lrc14_haar_zipper_cocycle_codex_s166.out`.
**Core readout:** For a `2 x 2` table, the mixed Haar / fixed-margin switch coordinate is `zeta(T)=T00-T01-T10+T11`.  Auditing all nonnegative tables through total `10` gives `1001` tables, `506` margin fibers, `285` nontrivial fibers, `0` duplicate `(margins,zeta)` keys, and fixed-margin zipper-step gcd `4`.  Thus row/column margins are an unsafe quotient unless `zeta` is retained, reconstructed, annihilated, or routed.
**Zipper theorem target:** On every labelled LRC14 packet fiber, each local mixed Haar cocycle is either cancelled by color-compatible discrepancy, stopped at a boundary cocircuit, handed to endpoint/period/certificate/smoothing clocks, descended to a smaller packet family, or converted into a named F7/THM-572 residual.  The depth-4 dyadic product census supplies the tooth taxonomy: orthogonal zero, same-tile boundary atom, owner strip, cross handoff, nested descent, residual cocycle.
**Next:** (1) build the actual HYP-2963 packet grid and compute packet-level `zeta` signatures; (2) count independent color-compatible cocycles and compare with HYP-2595's `k+c_GP`; (3) attach each nonzero owner-strip/cross/nested coefficient to a HYP-2987 handoff arrow; (4) define F7 as unpaired mixed-cocycle state-lift debt.  Namespace: HYP-2991 / T1075.
## Lead codex-2026-06-24-S166: LRC technique index for tournament/metagraph/series methods

**Status:** NAVIGATION / shared contribution index created and merged with the parallel S165 `LTI-*` registry.  User asked to look back through the many tournament, metagraph, sequence, and carrier methods and make a big index other agents can pull from and add to.  Added `00-navigation/LRC-TECHNIQUE-INDEX.md`.
**Contents:** the merged artifact has `107` compact `LTI-*` rows plus `63` long-form technique entries after adding the incoming HYP-2991 Haar zipper cocycle lane and HYP-2989 product-rule tiling synthesis.  Together they cover qdiv and endpoint geometry, Haar/Baire/tope/frontier methods, Fejer/Ramanujan/twist/moment duals, Farey/C27/K33 arithmetic carriers, metagraph transfer and operation-state methods, tournament speedup engines, deck derivatives, source-kernel/zipper/exposure methods, OCF/coimage/path-homology residues, sequence shadows, unit-distance analogies, graph-minor guardrails, and formal verification workflows.
**Contribution rule:** every new method should state carrier / vertex set, pairwise observable/gauge, preserved LRC predicate, what it destroys if scalarized, anchors, and next contribution.  This explicitly challenges the assumption that LRC Tournament Analysis vertices must be runners.
**Next:** keep this file as the living first-stop index for LRC sessions; update entries when scripts add certificate IDs, packet JSON schemas, endpoint-owner Ramanujan profiles, or F7 residual definitions.
## Lead codex-2026-06-24-S166: LRC14 zipper theorem pattern atlas

**Status:** SYNTHESIS / proof-interface atlas complete (HYP-2993/T1076), concrete LRC14/Haar-Fejer extension of HYP-2990's abstract zipper/no-free-slider atlas.  User asked to pull recent agents together, think discrepancy theory, use the 2D Haar product rule as the same structure as tournament tiling, and search for more zipper theorems.  The script `04-computation/lrc14_zipper_theorem_pattern_atlas_codex_s166.py` stores its output in `05-knowledge/results/lrc14_zipper_theorem_pattern_atlas_codex_s166.out`.
**Core schema:** A zipper theorem has two labelled local certificate sides, a labelled interface, declared stops, and named residuals.  A quotient may forget a coordinate only when the opposite tooth reconstructs it, orthogonality/boundary atoms annihilate it, or the coordinate is emitted as residual data strong enough for the next theorem.
**Computed atlas:** Ten patterns were scored: Haar-Fourier product, Fejer interval packet, tope/cocircuit wall, exposure-poset kernel, Ramanujan exact-period, smoothing/Kaczynski policy, fixed-margin/Johnson sector, apex sheaf gluing, convolution irreducibility lift, and unit-distance cyclotomic norm.  Tournament Analysis is transitive: `score_hist={0:1,...,9:1}`, `directed_3cycles=0`, `SCC_sizes` all singleton, and one Hamiltonian path.  Retention spine: `haar_fourier_product > tope_cocircuit_wall > exposure_poset_kernel > fejer_interval_packet > convolution_irreducibility_lift > ramanujan_exact_period > fixed_margin_johnson > smoothing_kaczynski_policy > apex_sheaf_gluing > unit_distance_cyclotomic_norm`.
**Next:** (1) build a Haar-Fejer compression engine over HYP-2963 packet rows by grouping mixed Haar switch signatures before interval certificate generation; (2) attach primitive-period `c_q` labels to endpoint-wall/Haar rectangle cells; (3) prove the HYP-2988 no-hidden-kernel target after HYP-2992 typed coefficients and HYP-2981 interval certificates are attached; (4) define F7 as a preserved fixed-margin/Johnson harmonic residual; (5) port HYP-2452 convolution no-lift certificates to LRC blocker ledgers.  Namespace: HYP-2991 / T1075.

## Lead codex-2026-06-24-S165: Haar product discrepancy and tournament tiling square

**Status:** PROOF-INTERFACE / synthesis scout complete (HYP-2989/T1073).  This pass integrates the recent HYP-2985/HYP-2986/HYP-2987/HYP-2988 packet, handoff, and exposure work with the older HYP-2594/HYP-2595 colored discrepancy program.  The script `04-computation/lrc14_haar_product_discrepancy_tiling_codex_s165.py` stores its output in `05-knowledge/results/lrc14_haar_product_discrepancy_tiling_codex_s165.out`.
**Core identity:** On dyadic children, the 2D Haar product `h_I(x)h_J(y)` is `[[1,-1],[-1,1]]`, which is exactly the 2-by-2 fixed-margin switch.  Diagonal and anti-diagonal packets have the same row/column margins but mixed Haar coefficients `+2` and `-2`; applying the switch jumps the coefficient by `4`.  Therefore row/column margins and raw continuous component counts are quotient shadows, while the mixed Haar switch is the packet that must be kept.
**Synthesis:** HYP-2594's `K` bound counts micro-boundaries before product cancellation; HYP-2595's colored resonance identity says only color-compatible mixed modes survive.  HYP-2986's tope/cocircuit wall language and HYP-2985's smoothing-dispatcher labels are the geometric and analytic forms of the same rule.
**Next:** (1) express the HYP-2595 resonance condition in a dyadic/Walsh-Haar packet basis; (2) count independent mixed switches for structured banks and compare to `k+c_GP`; (3) route switch families through HYP-2987's handoff atlas as O3 family-compression data and HYP-2988 exposure audits; (4) test whether Fejer packet certificates can be precompressed by shared mixed-Haar switch signatures.  Namespace: HYP-2989 / T1073.
## Lead codex-2026-06-24-S165: Haar-product tile discrepancy for LRC14

**Status:** PROOF-INTERFACE / exact finite product-table scout complete and folded into HYP-2989/T1073, alongside the minimal Haar-square switch and the broader Haar/Walsh product-algebra scout.  User prompt asked to synthesize recent agents through discrepancy theory and the two-dimensional Haar product rule, specifically the way it creates the same structure as the tournament tiling model.  The script `04-computation/lrc14_haar_product_tile_discrepancy_codex_s165.py` enumerates dyadic Haar rectangle products and stores the readout in `05-knowledge/results/lrc14_haar_product_tile_discrepancy_codex_s165.out`.
**Status:** PROOF-INTERFACE / exact finite product-table scout complete (HYP-2992/T1072), companion to HYP-2989/T1073's minimal Haar-square switch.  User prompt asked to synthesize recent agents through discrepancy theory and the two-dimensional Haar product rule, specifically the way it creates the same structure as the tournament tiling model.  The script `04-computation/lrc14_haar_product_tile_discrepancy_codex_s165.py` enumerates dyadic Haar rectangle products and stores the readout in `05-knowledge/results/lrc14_haar_product_tile_discrepancy_codex_s165.out`.
**Core readout:** Through depth `3`, `225` rectangles give `50625` ordered products: `43736` orthogonal zeros, `225` same-tile boundary atoms, `1020+1020` one-coordinate owner strips, `2312` cross handoffs, and `2312` nested refinements.  Every nonzero non-atom class is sign-balanced.  This suggests a discrepancy proof by typed coefficient detection rather than scalar positivity.
**Transfer to LRC14:** Use a rectangular packet grid whose axes are endpoint/tope walls and proof clocks (exact period, Fejer scale, Ramanujan packet, K33/state-lift route).  Then prove a vanishing lemma: if a primitive zero-open packet has no owner-strip, cross-handoff, or nested-refinement Haar coefficient, it collapses to AP/GW boundary skeleton or creates the missing THM-572/F7 state-lift atom.
**Next:** (1) build the actual packet grid for named rows AP, GW, `12->36`, `10->20`, `13->26`, `P10+GW`, `12->168`, and lcm-tail rows; (2) compute typed Haar coefficients by packet family; (3) compare same-tile boundary atoms with the HYP-2986 boundary cocircuit owner-pair sums; (4) test whether HYP-2981 Fejer atoms cluster in owner-strip or nested-refinement classes; (5) formalize the dictionary from fixed-path tournament staircase tiles to Haar rectangles.  Namespace: HYP-2989 / T1073.
## Lead codex-2026-06-24-S165: abstract zipper theorem atlas across past topics

**Status:** EXPLORATORY PROOF-TECHNOLOGY / HYP-2990/T1074 created after rebasing over incoming HYP-2988 exposure and HYP-2989 Haar-product discrepancy lanes.  User asked for more zipper theorems and free creative exploration of past topics.  Added `04-computation/abstract_zipper_theorem_atlas_codex_s165.py` and stored `05-knowledge/results/abstract_zipper_theorem_atlas_codex_s165.out`.
**Core rule:** A zipper theorem is a controlled-kernel theorem.  A quotient may forget a coordinate only if the target predicate is fiber-constant, the coordinate is reconstructible, a dual certificate annihilates it, or it is routed to a named residual sector.
**Carriers compared:** LRC14 certificate handoff, kernel/tope/smoothing, octahedral Hodge current, C27/unital pair completion, Farey/K33 incidence, boundary-moment multi-chart, shell-1/root packet, unit-distance endpoint ear, OCF activity/coimage, good-cut/SCC support, and raw scalar shadow.
**Tournament fingerprint:** `score_hist={0:1,1:1,2:1,4:2,5:1,6:1,7:2,9:1,10:1}`, `directed_3cycles=4`, `SCC_sizes=[1,1,1,6,1,1]`, `Hamiltonian_path_count=15`.  The six-carrier middle SCC says several past topics form typed handoff teeth rather than a scalar hierarchy.
**Next:** (1) define HYP-2987's `F7` as a named harmonic/state-lift residual sector; (2) turn the no-free-slider rule into an LRC14 quotient checklist; (3) test family compression on K33/petal Fejer packets; (4) revisit octahedral divergence/curl and good-cut/SCC support as state-lift coordinates.
**Next:** (1) build the actual packet grid for named rows AP, GW, `12->36`, `10->20`, `13->26`, `P10+GW`, `12->168`, and lcm-tail rows; (2) compute typed Haar coefficients by packet family; (3) compare same-tile boundary atoms with the HYP-2986 boundary cocircuit owner-pair sums; (4) test whether HYP-2981 Fejer atoms cluster in owner-strip or nested-refinement classes; (5) formalize the dictionary from fixed-path tournament staircase tiles to Haar rectangles.  Namespace: HYP-2992 / T1072.

## Lead codex-2026-06-24-S164: admissible smoothing dispatcher for LRC14

**Status:** PROOF-INTERFACE / routing theorem target complete (HYP-2985/T1069), complementary to the incoming HYP-2984/T1068 kernel-homotopy stub.  This pass turns the recent analytic-number-theory prompts into a typed LRC14 smoothing dispatcher rather than a scalar estimate.  The script `04-computation/lrc14_smoothing_dispatcher_codex_20260624.py` classifies which policy is allowed to handle each live packet family and stores the readout in `05-knowledge/results/lrc14_smoothing_dispatcher_codex_20260624.out`.
**Core split:** AP/GW boundary atoms route to Kaczynski/endpoint labels; K33 and covering fronts route to Fejer/Toeplitz interval certificates; q=27 petals and P10+GW splices require Fejer plus Ramanujan prime-power side channels; late prime-power denominator walls can use Selberg/large-sieve only as labelled preconditions; true-wide off-resonance packets route to Kaczynski/Abel signed decay; true-wide resonant packets route to Freiman finite reduction or HYP-2908/THM-572 state lift.
**Clocks to keep:** endpoint-owner clock, exact-period denominator clock, smoothing/certificate clock, and far-approach boundary clock.  Raw prime counts, `sum mu`, and `sum mu/n` are diagnostics or tails unless reattached to packet labels.
**Next:** (1) add endpoint-owner Ramanujan profiles for q=25,27,41,4312; (2) tag HYP-2981 Fejer certificates with smoothing-policy labels; (3) prove a Kaczynski/Abel off-resonance relation-height bound; (4) build a Freiman finite-atlas handoff for resonant true-wide packets; (5) encode the admissible-smoothing lemma over HYP-2963 fibers.  Namespace: HYP-2985 / T1069.
## Lead codex-2026-06-24-S164b: LRC14 certificate handoff atlas and zipper theorem target

**Status:** PROOF-INTERFACE / handoff atlas complete, six theorem arrows open (HYP-2987/T1071).  User prompt asked for another creative LRC proof pass while keeping the labelled-packet, tournament-analysis, boundary-moment, quotient guardrail, Fejer PSD, Ramanujan/divisor, and fixed-margin themes in play.  S164 adds `04-computation/lrc14_certificate_handoff_atlas_codex_s164.py` and stores `05-knowledge/results/lrc14_certificate_handoff_atlas_codex_s164.out`.  The script makes proof carriers the Tournament Analysis vertices and scores retained LRC predicate, exact scale, phase/period, topology, endpoint owners, packet family, dual certificate, formal checkability, and residual routing.
**Core readout:** The carrier tournament has `score_hist={0:1,1:1,2:1,4:3,6:1,7:1,8:1}`, one directed 3-cycle, SCC sizes `[1,1,3,1,1,1,1]`, and `3` Hamiltonian paths.  The nontrivial SCC is the live middle layer: Ramanujan exact-period packets, endpoint bridge graphs, and twist ladders can all be correct depending on the predicate preserved, so handoffs must declare their retained labels.
**Zipper theorem target:** If O1 source-kernel exclusion, O2 formal interval backend, O3 family compression, O4 admissible smoothing, O5 state-lift construction, and O6 F7 definition hold, then every primitive LRC14 row either has a strict witness/dual certificate, is the AP/Goddyn-Wong equality atom, or constructs the forbidden HYP-2908/THM-572 tournament state lift.
**Fixed-margin import:** `arXiv:2606.22636` proves a fixed-margin binary swap-chain spectral-gap theorem by preserving margin fibers, reducing to a low-row heat-bath core, and decomposing by count plus Johnson harmonic sectors.  The LRC14 transfer is to preserve packet fibers, reduce arbitrary primitive rows to a finite packet core, and define F7 as a named harmonic residual sector instead of an untyped failure bucket.
**Rebase signal:** HYP-2984/T1068 is now the sibling kernel-homotopy and boundary-defect lane, HYP-2985/T1069 is the admissible-smoothing dispatcher lane, and HYP-2986/T1070 is the tope-wall/cocircuit lane.  Treat all three as admissible-arrow proofs inside this HYP-2987 handoff atlas: smoothing or Fourier kernel deformations must preserve packet certificates or emit named boundary defects, smoothing policies must keep endpoint/exact-period/kernel/far-approach clocks separate, and endpoint arrangements must separate open topes from boundary cocircuits.
**Next:** (1) make O3 concrete by compressing selected Fejer certificates into packet-family templates; (2) make O5 concrete by constructing the THM-572 state lift from zero-open non-AP/GW residuals; (3) formalize O6 as a Johnson-harmonic residual predicate with count-sector separation; (4) write the O1 source-kernel exclusion as a typed reduction theorem over qdiv/Farey/Haar/endpoint labels.  Namespace: HYP-2987 / T1071.

## Lead codex-2026-06-24-S163: analytic sieve packet weights and Kaczynski boundary guardrails

**Status:** EVIDENCE / finite arithmetic atlas and proof-carrier guardrail complete (HYP-2982/T1066).  User prompt asks to merge explicit sums over primes, `sum mu(n)`, `sum mu(n)/n`, `sum mu(n)^2/phi(n)`, large-sieve/circle-method improvements, upper-bound quadratic/Selberg sieve ideas, exponential sums, ternary-Goldbach parabolic-cylinder/saddle-point/explicit-formula machinery, smoothing choices, and the repo's Kaczynski/Kaczorowski threads.  S163 computes a finite arithmetic atlas through `N=200000`: `Phi(z)=sum_{q<=z}phi(q)` is quadratic primitive packet capacity, while `G(z)=sum_{d<=z}mu(d)^2/phi(d)` is logarithmic inverse-unit normalizer (`G(200000)-log(200000) ~= 1.332518`).  Therefore analytic sieve estimates are useful as middle certificates only when their kernel, smoothing transform, exceptional-set boundary, and retained LRC labels are declared.
**Denominator warning:** `mu^2/phi` is valuable but squarefree-blind.  It keeps `q=14` and prime `q=41`, while zeroing live prime-power or repeated-prime exact-period packet denominators `25`, `27`, `36`, `63`, `84`, `98`, `168`, `280`, and `4312`.  Thus a large-sieve or quadratic/Selberg-sieve quotient can be an upper-bound/minor-arc component, not a final equality atom, unless prime-power, endpoint-owner, Ramanujan, interval-Fejer, or Kaczynski approach-class labels are restored.
**Helfgott/circle-method import:** Treat the ternary Goldbach architecture as a proof-design pattern: local arithmetic labels, sieve weights, exponential sums, major/minor arc split, smoothing choices, explicit transform/special-function backend, and rigorous inequality.  LRC14 analogue: q/Farey/endpoint labels, Ramanujan/divisor packets, Fejer/Toeplitz exponential sums, strict-open/boundary packet split, familywise smoothing, interval trig backend, and labelled certificate or state lift.
**Kaczynski import:** HYP-2679/THM-548 make approach classes proof data: bounded core is boundary datum, far packet is approach path, and resonance corrections are ambiguous-boundary labels.  Do not forget the boundary approach class when applying smoothing or decorrelation.
**Rebase signal:** HYP-2983/S162 now supplies the companion Kaczynski/exponential-sum synthesis and motif census.  Treat HYP-2982 as the finite `Phi`/`G` weight atlas and squarefree-blindness audit; treat HYP-2983 as the coupled proof-template lane for labelled source kernels, exponential sums, smoothing, and boundary resonance.
**Next:** (1) add prime-power side channels to HYP-2979 for `q=25,27` and Fejer denominator `4312`; (2) split HYP-2981 interval certificates by smoothing family; (3) test `mu^2/phi` capacity only after reattaching exact-period labels; (4) record Kaczynski boundary approach class in true-wide/far-speed packet certificates; (5) build the LRC14 "explicit explicit formula" emitter with atom-bank provenance, trig intervals, endpoint owners, and packet handoff.  Namespace: HYP-2982 / T1066.

## Lead codex-2026-06-24-S162: Robbins/Robin bridge guardrails and interval Fejer packet certificates

**Status:** ACTIVE PROOF-INTERFACE / packet-anchored scaffold, precision blueprint, and named-row interval prototype complete; production proof backend pending (HYP-2981/T1065).  User prompt asks to use Robbins theorem in graph theory, Robin/Robbins number-theory divisor-function readings, quotient guardrails, Ramanujan sums, and the current Fejer result.  HYP-2981 now combines `04-computation/lrc14_packet_fejer_interval_scaffold_codex_s162.py`, `04-computation/lrc14_fejer_interval_packet_certificates_codex_s162.py`, and `04-computation/lrc14_fejer_interval_named_packet_certificates_codex_s162.py`.  S157 already found floating Fejer PSD-vector violations for every positive HYP-2963 packet-bank row (`21911/21911`) by degree `<=280`, with AP/Goddyn-Wong as the only zero-safe atoms.  The scaffold attaches selected hard Fejer forms to packet keys `P(S)` and gives rational interval upper bounds `<0` for `12->36`, `P10+GW`, `12->168`, `drop(12,13)->add(14,29)`, and `6->63`; the budget script expands hard certificates into divisor-curried atom banks and shows the interval burden looks finite: `P10+GW` is the high-degree/large-bank case (`d=280`, `862` atoms), K33 `12->36` needs `d=159`, and the weakest full-bank margin row has `122` atoms and about `27` conservative precision bits.  The direct `mpmath.iv` named-row prototype anchors exact safe components and certifies all nine positive named packets, including `K33 12->36` and `P10+GW`; the least negative interval upper is `-1.30202901999956e-05` on covering `12->168`.
**Robbins/Robin guardrail:** Graph Robbins says a connected graph can be strongly oriented iff no bridge remains; Robin's divisor-function theorem reduces the RH-equivalent sigma inequality to the extremal divisor-density fiber; Neville Robbins partition/cyclotomic side readings reinforce that divisor functions are fibers before they are scalars.  LRC14 translation: a quotient is admissible only when every bridge/fiber it forgets is named and discharged, so a Fejer certificate tuple must retain family/route, exact safe component, rational center, degree, coefficient-fiber formula, signed interval `Q` with `upper(Q)<0`, packet fiber, and route handoff.
**Next:** (1) replace the scaffold's hard-coded `pi` interval/Taylor backend and the named prototype's `mpmath.iv` backend with formally sourced Lean/arb-compatible interval payloads; (2) lift exact sign certificates from rows to packet-family templates (AP/GW, K33, petal/two-block, covering, few-apex); (3) use HYP-2979 Ramanujan projectors to pre-split late-q packet families; (4) prove a no-bridge quotient lemma: forgotten coordinates are reconstructible, annihilated by a dual certificate, or routed to a named residual/state-lift bucket.  Namespace: HYP-2981 / T1065.
**Status:** ACTIVE PROOF-INTERFACE / packet-anchored scaffold, manifest, and precision blueprint complete, production interval backend pending (HYP-2981/T1065).  User prompt asks to use Robbins theorem in graph theory, Robin/Robbins number-theory divisor-function readings, quotient guardrails, Ramanujan sums, and the current Fejer result.  HYP-2981 now combines `04-computation/lrc14_packet_fejer_interval_scaffold_codex_s162.py`, `04-computation/lrc14_fejer_interval_packet_certificates_codex_s162.py`, and the S163 manifest `04-computation/lrc14_fejer_packet_certificate_manifest_codex_s163.py`.  S157 already found floating Fejer PSD-vector violations for every positive HYP-2963 packet-bank row (`21911/21911`) by degree `<=280`, with AP/Goddyn-Wong as the only zero-safe atoms.  The scaffold attaches selected hard Fejer forms to packet keys `P(S)` and gives rational interval upper bounds `<0` for `12->36`, `P10+GW`, `12->168`, `drop(12,13)->add(14,29)`, and `6->63`; the budget script expands hard certificates into divisor-curried atom banks and shows the interval burden looks finite: `P10+GW` is the high-degree/large-bank case (`d=280`, `862` atoms), K33 `12->36` needs `d=159`, and the weakest full-bank margin row has `122` atoms and about `27` conservative precision bits.  S163 stores these five selected hard rows as theorem-facing certificate records with `certified_negative=True`, exact packet keys, rational centers, degree, interval sign/digit sizes, and the Robbins bridge contract from center to route handoff.
**Robbins/Robin guardrail:** Graph Robbins says a connected graph can be strongly oriented iff no bridge remains; Robin's divisor-function theorem reduces the RH-equivalent sigma inequality to the extremal divisor-density fiber; Neville Robbins partition/cyclotomic side readings reinforce that divisor functions are fibers before they are scalars.  LRC14 translation: a quotient is admissible only when every bridge/fiber it forgets is named and discharged, so a Fejer certificate tuple must retain family/route, exact safe component, rational center, degree, coefficient-fiber formula, signed interval `Q` with `upper(Q)<0`, packet fiber, and route handoff.
**Next:** (1) replace the scaffold's hard-coded `pi` interval/Taylor backend with a formally sourced Lean/arb-compatible interval backend; (2) lift the S163 selected-row records to packet-family templates (AP/GW, K33, petal/two-block, covering, few-apex); (3) use HYP-2979 Ramanujan projectors to pre-split late-q packet families; (4) prove a no-bridge quotient lemma: forgotten coordinates are reconstructible, annihilated by a dual certificate, or routed to a named residual/state-lift bucket.  Namespace: HYP-2981 / T1065.

## Lead codex-2026-06-24-S161: Ramanujan-divisor quotient guardrails for LRC14

**Status:** ACTIVE INQUIRY / web-crawled, named-row collision audit complete, packet-family extension pending.  User prompt: read the divisor-function neighborhood and use it to formalize quotient guardrails across irreducibility, unital designs, Faulhaber moments, Pollock defects, unit-distance carriers, tiling/solid analogies, multiplicative functions, and Ramanujan sums.  Core question: what may a quotient forget without breaking the proof predicate?  Candidate proof object: a divisor/cyclotomic packet ledger where coarse multiplicative signatures are explicitly stress-tested against phase-sensitive Ramanujan signatures and LRC14 route labels.
**S161 collision-audit update:** named-row audit confirms qdiv, open-state, mod-14 residue, `c_14`, unit-count, and lcm-scalar quotient collisions across AP/GW, q-witness, K33, petal, and covering routes; only the over-labelled guarded packet signature avoids route mixing.
**Web/source update:** divisor functions supply the scalar pushforward; Dirichlet convolution supplies packet laws; Jordan totients supply primitive tuple capacity; divisor summatory functions expose boundary defects after a product-lattice quotient; Lambert/Eisenstein series and Ramanujan expansions connect divisor fibers to harmonic coefficients; Ramanujan-sum orthogonality and supercharacters justify unit-orbit quotients only when the primitive trace is retained.
**Next:** (1) extend quotient-collision audit to packet families in the HYP-2963 bank; (2) test endpoint-owner Ramanujan profiles `R_q^+`, `R_q^-` for `q in {14,27,41}`; (3) turn the controlled-kernel criterion into a theorem-facing admissibility lemma; (4) route any mixed fibers to endpoint-owner, C27/K33, Toeplitz, moment-dual, or THM-572 state-lift labels.  Namespace: HYP-2978 / T1062.
**Artifacts:** `04-computation/lrc14_ramanujan_divisor_quotient_guardrails_codex_s161.py`; `05-knowledge/results/lrc14_ramanujan_divisor_quotient_guardrails_codex_s161.out`; `07-reflections/lrc14-ramanujan-divisor-quotient-guardrails-codex-s161.md`.

## Lead codex-2026-06-24: Ramanujan exact-period projectors for LRC14

**Status:** ACTIVE PROOF-INTERFACE (HYP-2979/T1063, companion to HYP-2978). Ramanujan sums turn primitive
**Status:** ACTIVE STUB (HYP-2979/T1063, companion to HYP-2978). Ramanujan sums turn primitive
**Status:** EVIDENCE / quotient-admissibility guardrail, not an LRC14 proof.  User prompt: read the divisor-function neighborhood and use it to formalize quotient guardrails across irreducibility, unital designs, Faulhaber moments, Pollock defects, unit-distance carriers, tiling/solid analogies, multiplicative functions, and Ramanujan sums.  HYP-2978 now states the rule: a quotient may be a proof carrier only when the LRC predicate is constant on its fibers, or when the quotient carries a named lost-label certificate.  Exact audit `lrc14_ramanujan_divisor_quotient_guardrails_codex_s161.py` checks `2694` rows.  Scalar divisor signatures have `138` mixed qdiv/safe-route fibers and `239` bad pair-collisions; unitary divisors reduce to `12/18`; exact-period packets reduce to `14` one-pair misses but still identify AP with positive `12->96`, proving endpoint/safe-measure labels or K33 state-lift debt must be reattached.  Tournament Analysis over quotient carriers is transitive with path `endpoint_measure > full_row > exact_period_packet > unitary_divisor > scalar_divisor > ramanujan_pair > ramanujan_speed > qcover`.
**Next:** (1) add endpoint-owner labels to the HYP-2979 exact-period projector and rerun the fiber test on the HYP-2963 bank; (2) compare `c_14(v_i+v_j)` directly to HYP-2970 endpoint credit `K=14(rm-sn)+r+s`; (3) test shifted Carmichael/Ramanujan autocorrelation of danger multiplicity against HYP-2973 moment-duals and HYP-2974 Toeplitz PSD certificates; (4) keep multiplicative functions as irreducibility ledgers, not proof-ending scalars, unless the lost-label certificate is explicit.  Namespace: HYP-2978 / T1062.
## Lead codex-2026-06-24: Ramanujan exact-period projectors for LRC14

**Status:** ACTIVE STUB (HYP-2979/T1062), complementary to HYP-2978's quotient guardrail. Ramanujan sums turn primitive
**Status:** INITIAL AUDIT COMPLETE.  User prompt: read the divisor-function neighborhood and use it to formalize quotient guardrails across irreducibility, unital designs, Faulhaber moments, Pollock defects, unit-distance carriers, tiling/solid analogies, multiplicative functions, and Ramanujan sums.  Core question: what may a quotient forget without breaking the proof predicate?  S161 verifies the arithmetic identities and shows named-row route collisions for every coarse quotient tested; only the guarded packet signature avoids route mixing.
**Next:** (1) extend the quotient-collision audit from named rows to HYP-2963-bank rows; (2) compare primitive phase profiles with Toeplitz/Fejer centers and danger-count moments; (3) formalize the admissibility theorem target.  Namespace: HYP-2978 / T1062.
## Lead codex-2026-06-24: Ramanujan exact-period projectors for LRC14

**Status:** ACTIVE STUB (HYP-2979/T1063, child route of HYP-2978). Ramanujan sums turn primitive
q-th roots into an integer projector:
`c_q(n)=sum_{(a,q)=1}exp(2*pi*i*a*n/q)=sum_{d|gcd(q,n)}d*mu(q/d)`.
This looks tailor-made for the current LRC14 stack because the hard rows are
already organized by exact denominators, primitive unit witnesses, endpoint
taut equalities, and Fourier/Toeplitz shadows. For `q=14`, `c_14` splits
differences into the four parity/seven-adic classes `6,-6,-1,1`; use
`c_14(r+s)` for AP/Goddyn-Wong zero-credit traces and `c_14(r-r')` for residue
coincidence packets. Next tasks: compute primitive phase witness profiles for
the HYP-2963 bank; test shifted Carmichael autocorrelation of `N_S(t)` against
HYP-2973 moment-duals; compare Ramanujan modes to HYP-2974 Toeplitz Fejer
certificates; and determine whether a Ramanujan-subspace packet catches any
qdiv>14 zero-open SOURCE-SPECTRUM-UNKNOWN residual. External source trail:
Ramanujan sums/roots of unity, Carmichael orthogonality, finite-duration
Ramanujan subspaces, cyclotomic primitive-root traces, and supercharacter
Ramanujan sums over `Z/n`. HYP-2978 now supplies the guardrail: exact-period
projectors must retain or reattach q/Farey, Haar, endpoint, C27/K33, and
state-lift labels before they are admissible proof quotients.

**2026-06-24 update:** `04-computation/lrc14_ramanujan_exact_period_projector_codex_20260624.py`
audits named rows plus `21906` AP-neighborhood rows.  Every row has a weak
primitive witness by `q<=42`; only AP and GW lack a strict primitive witness by
`q<=42`.  q=14 primitive phase packets still mix routes, so the next task is to
scale the strict primitive-witness audit to the full HYP-2963 bank and
interval-enclose late q packets `{25,27,34,40,41}`.

## Lead mac-mini-2026-06-21-S20: LRC LAYER-3 sharpened routes + key literature (from the creative-lead trawl)
**Source:** S20 trawl + drill workflows (HYP-2760..2763). **Status:** ACTIVE — LAYER 3 (consec Schur-maximizes measS7) is the last wall; LAYERS 1-2 PROVED+Lean.
- **Conductance invariant (HYP-2760):** measS7 = EXACT function of c_r = sum_{e=r mod7} 1/|e| (0 collisions); Foster sum_r g(r)=112; bottleneck/minimax on C_7. NEXT: derive the exact F(c) and prove consec maximizes the bottleneck via Chebyshev equalization.
- **Windows extremality (HYP-2761):** consec UNIQUELY maximizes WIN=harmonic sum of binding speeds (0 ties, sharper than the wall). OBSTRUCTION: WIN/DISCONNECTED split does NOT compose (they anti-correlate). NEXT: harmonic-sum exchange argument.
- **AP-unification (HYP-2762):** BOTH extremalities = 'the additive AP Schur-maximizes' (consec / cyclic Interval); Paley is small-p only (crossover p=19, HYP-479/THM-135). Wiener-Khinchin: argmin sum|lambda|^4 = QR difference set. NEXT: the elementary-symmetric c_k crossover mechanism; does it transfer to LRC?
- **ROUTE 5 RESOLVED (mac-mini-S20, the c_k crossover + transfer test):** The unification is NOT one c_k theorem — the two extremalities are driven by OPPOSITE level + OPPOSITE-sign fugacity. (i) The spectral "H = c_0 + sum c_k e_k(|lambda|^2)" of THM-134 is UNDERDETERMINED by circulant data (only 2 classes at p=7, 4 at p=11 vs m unknowns) and NON-canonical; THM-505 already shows H is not a low-degree spectral polynomial. The GENUINE canonical decomposition on both sides is a SIGNED SUBSET-COMPLEX LEVEL-SUM. (ii) TOURNAMENT: H = sum_j (+2)^j alpha_j (odd-cycle independence poly, fugacity +2). Paley wins level 1 (max cycle counts c_ell, every ell), Interval/AP wins levels j>=2 (disjoint packings). Crossover fugacity x*(p): 3.00 (p=7) -> 2.61 (p=11) -> below 2 by p=19 (= THM-135). The AP wins because +2 AMPLIFIES the high-order packing levels where it dominates. (iii) LRC: measS7 = sum_k (-1)^k MISS_k (sector-miss inclusion-exclusion, fugacity -1; MISS_k = sum_{|S|=k} P(all of S missed)). consec UNIQUELY MINIMIZES MISS_1 (best single-sector coverage, rank 1/319) and wins via the LOW-order level that the alternating -1 fugacity rewards. Sweep G_x = sum_k x^k MISS_k over all 319 shapes: at x=-1,-1/2 consec is unique argmax (rank 1); at x=+1/2,+1,+2 consec drops to rank 319/319, 316/319, 220/319 (worst/near-worst). So consec-optimality is SPECIFIC to x=-1. PRECISE OBSTRUCTION (refines HYP-2758 no-Parseval-simplex): the tournament AP wins HIGH-order packing at fugacity +2; the LRC AP wins LOW-order coverage at fugacity -1 — same "additive AP extremal" surface, reversed level and reversed fugacity sign. Scripts: 04-computation/lrc14_route5_*.py; results: 05-knowledge/results/lrc14_route5_*.out.
- **KEY LITERATURE (HYP-2763):** **Rosenfeld 2025 (arXiv:2509.14111, arXiv:2511.22427) PROVES LRC n=8,9,10** via finite-checking + the Tao 2018 -> Malikiosis-Santos-Schymura 2025 bound n^{2n} (arXiv:2411.06903). STRATEGIC Q: is the finite check feasible at n=14 (n^{2n}=14^28 brute is not, but our sector-route + stratum-localization + sharp L7 atlas p<=14 reduces it massively)? Huffer-Shepp 1987 reflection lemma (port to per-cell W_a). Cusick view-obstruction (AP {1..n} proven tight). Chung-Graham 3-distance (AP gap-rigidity => LAYER 3 as a rigidity statement).

---

## Lead kind-pasteur-2026-06-17-S2: close LRC(14) covering case via a BOUNDED-SPEED reduction (THM-525, OPEN-Q-108)

**Status:** 🔴 HIGH — the single most valuable next step toward LRC(14). THM-525 reduces the covering case to OPEN-Q-108 (Q=uniform meas(G_C)≥c) and has ~105k covering 13-sets with ZERO counterexamples, but the scan is EVIDENCE, not a certificate. **The move that would CLOSE it:** a bounded-speed reduction (Tao "Some remarks on the lonely runner" Thm 1.3; Malikiosis–Sgall–Somer / Dubickas finite-checking) proving LRC(14) on covering sets need only be verified for v_max ≤ an explicit V₀. Then the exact scan becomes a finite proof up to V₀ and only the V₀-bound remains. SECOND target: the named sub-gap **G2** (transversality) — prove G_C (a union of arcs whose endpoints are binding-pair crossings k/(v_a±v_b), all v ≤ bound) cannot concentrate in O(1/w)-neighborhoods of {a/w}, so w's danger comb (measure 1/7) cannot CONTAIN G_C. Either one moves OPEN-Q-108 from "localized + strongly evidenced" to "closed."
**Next (ranked):** (1) find/derive the explicit V₀ for N=14 covering sets from Tao Thm 1.3 / MSS; (2) run the exact covering-set scan to v_max≤V₀ (the THM-525 tooling already does covering enumeration); (3) attack G2 directly via the binding-crossing-endpoint structure of G_C; (4) test whether the circle/local-tournament structure (HYP-2576) pins G_C away from the parked comb.
**Sources:** THM-525, OPEN-Q-108 (GAP A + new GAP G2), HYP-2573..2576; `04-computation/lrc14_{easy_dominates_hard,parked_centering,verify_reduction_*}_kps-S2-wf.py`; reflection `the-perfect-middle-is-the-wrong-fixpoint-and-the-lrc-tournament-is-local-not-lonely-kps-S2.md`.

## Lead mac-mini-2026-06-18-S7: ANGLE C — three-gap / dilation-orbit extremality of meas(S7) (THM-532 closure)

**Status:** LITERATURE-SCOUTED + COMPUTE. The THM-532 crux ("consecutive AP maximizes meas(S7) for dangerous rows k=8..11") is now grounded against the three-gap / dilation-orbit literature, and the deployment is HONESTLY MAPPED. KEY EXTERNAL CONFIRMATIONS: (a) The LRC state of the art (Pham–Sawhney–Tidor–Trakulthongchai? "Eleven, twelve, and thirteen lonely runners", arXiv 2604.23906) proves LRC(k) for k≤12 and **identifies (1,2,…,k) — the consecutive AP — as the unique bottleneck / hardest tuple** that forces the polynomial-method step; k=13 (=LRC(14)) is genuinely the first open case. This is INDEPENDENT literature confirmation that "AP is the extremal configuration" for the dangerous rows. (b) Bedert, "Riesz products and the LRC" (arXiv 2511.16636): `ML(1,2,…,n)=1/(n+1)` exactly — consecutive integers are the *unique* extremal case — and the proof dichotomy is **low additive dimension (AP-like) vs. high additive dimension (dissociated → Riesz product)**, the SAME split as THM-532's relation-height (low W = AP-rich vs. high height = dissociated). (c) Tao "Some remarks…": union bound is tight to `(1+O(1/n))` exactly for AP-like sets. (d) The orbit `{0,x,…,(k-1)x}` of the AP is EXACTLY a Steinhaus three-gap configuration (≤3 distinct gaps — verified k=8..13); a general E gives a union of d dilated APs with up to **3d gaps** (Liang 1979, "A short proof of the 3d distance theorem", Discrete Math 28; generalization arXiv 1910.00865) — so non-AP orbits are LESS uniform, the mechanism behind AP-extremality. COMPUTE (this session): AP is the UNIQUE maximizer of BOTH meas(S7) AND the difference-energy `D(E)=Σ_d mult(d)²` for k=8..11 (exact small-box census, 0 beaters); but Pearson corr(meas(S7),D)≈0.53–0.64 only, so **meas(S7) is NOT a monotone function of any single scalar** — AP-extremality is a JOINT extremum, not reducible to a one-line three-gap inequality. HONEST VERDICT: the literature CONFIRMS and CONTEXTUALIZES AP-extremality (strong external corroboration) but does NOT supply a ready-made theorem that PROVES `meas(S7(E))≤meas(S7(AP_k))` for k=8..11; the three-gap law bounds the *number* of gaps, not the *cover measure*. No clean three-gap closed form for `meas(S7(AP_k))` either (denominators 1470,5880,17640,5880,194040,21560 = 7·lcm-type, not a tidy product).
**Next (ranked):** (1) Try the difference-energy as a TWO-piece certificate: prove (i) meas(S7) ≤ f(D(E)) for an explicit increasing envelope f, then (ii) the finite residual {E : D(E) ≥ D₀} is AP-near and small — this mirrors THM-532's W-split but with D (a genuine additive-energy, literature-standard) instead of the ad-hoc W. (2) Pursue the Liang 3d-distance bound as a hard cap: if orbit(E) has ≥ 4 distinct gaps on a positive-measure set, meas(S7) loses a quantifiable amount vs the 3-gap AP — make this quantitative (the missing inequality). (3) Read the FULL 2604.23906 / 2511.16636 PDFs (WebFetch failed on binary; use the HTML or ar5iv) for the exact AP-bottleneck lemma — it may directly bound the AP cover and transfer. (4) Check whether Bedert's Riesz-product certificate, applied to the dangerous *clusters* (not all of V), gives meas(S7)≤cap directly for high-additive-dimension E, leaving ONLY the AP-near finite residual.
**Sources:** THM-532; `04-computation/lrc14_angleC_threegap_{extremality,mechanism}_macmini_0618s7.py`; `05-knowledge/results/lrc14_angleC_threegap_*_macmini_0618s7.out`; arXiv 2604.23906, 2511.16636, 1910.00865; Liang Discrete Math 28 (1979) 325; Tao "Some remarks on the lonely runner conjecture" (2017); Balog–Granville–Solymosi arXiv 1410.8404.

## Lead codex-2026-06-18: HYP-2608 empty-window moment lower certificate for LRC(14)

**Status:** NEW complement-side certificate.  Instead of upper-bounding the bad seven-sector cover `S7(E)`, this route lower-bounds the good anchored empty-window union `EWLB_A(E)=meas(union_j W_j(E))`, where `W_j` is the event that the open arc `(j/14,j/14+1/7)` is empty.  Let `R(x)=# {j:x in W_j(E)}` and `T_s=E[C(R,s)]`.  Polynomial minorants `h_R(r)<=1[r>0]` give proved per-set lower bounds `EWLB_A(E)>=Phi_z(E)=sum z_s T_s(E)`.  Degrees `4,3,3,2,1` clear the AP rows `k=8..12`; exact bounded primitive banks (`5797` rows total) have `0` threshold failures, with AP minimizers for dangerous rows `k=8..11`.  This matches the user's region-first prompt: the natural tournament vertices are fixed loop regions/proof obligations, not runners.  AP region-load tournaments are transitive with score histogram `0..6`.
**Next (ranked):** (1) Prove the scalar rearrangement `Phi_k(E)>=Phi_k(AP_k)` for `k=8..11`, likely by a three-gap/Sturmian or relation-height split for the empty-window count distribution. (2) For `k=12`, prove a coarse first-moment bound; AP extremality is unnecessary because the degree-1 certificate has large slack. (3) Compare HYP-2608's minorant polynomials to THM-534's majorants and look for a shared Krawtchouk/Bonferroni algebra that proves both scalar extremalities. (4) Build a fast integer common-refinement engine for `T_s(E)` and widen the exact banks to match HYP-2604's AP-frontier boxes.
**Sources:** `04-computation/lrc14_empty_window_moment_lower_codex_0618.py`, `05-knowledge/results/lrc14_empty_window_moment_lower_codex_0618.out`, `05-knowledge/hypotheses/HYP-2608-lrc14-empty-window-moment-lower-certificate.md`; HYP-2603, HYP-2604, THM-531, THM-532, THM-534.

## Lead codex-2026-06-17-S4: unit-distance n=22 Mathieu residue (HYP-2572, T840)

**Status:** NEW exact finite scout.  The unresolved `60 <= u(22) <= 61` frontier is a one-point extension problem: any 61-edge graph has a degree-4 or degree-5 deletion ear over a `57/56`-edge 21-core.  The Mathieu chain suggests retaining the `M22` point-stabilizer side-channel: fixing a point exposes `M21 = L_3(4) = PSL(3,4)` on `PG(2,4)`, whose 21 lines are exactly the fixed-point residual hexads of `S(3,6,22)`.  The new script verifies the plane and classifies ear-neighbor sets: degree-5 has `line_5` (`21`), `near_line_4_plus_1` (`1680`), `arc_5_no_three` (`1008`), and `three_collinear_5` (`17640`); degree-4 has `line_4` (`105`), `arc_4_no_three` (`2520`), and `three_collinear_4` (`3360`).  Unit-circle caps give at most `4`/`3` internal unit chords for degree `5`/`4` ears.  Tournament Analysis over proof carriers is transitive with leader `line_5_hexad_ear`.
**Next (ranked):** (1) Attach `PG(2,4)` ear types to stored `n=21` core and `n=22` extension candidates; record which extensions are line, punctured-line, near-line, or scattered. (2) Prove/certify coherent-ear exclusions using the circle-cap counts plus the existing Moser cap-endpoint ledger. (3) Build a scattered-ear obstruction library indexed by secant profile and compare it against the totally-unfaithful killers from the graph-only 62-edge coimage. (4) Test whether the `M21` residue aligns with `P_2^- / P_2^+` spine-ladder rows and the `57=20+37` centered-hex carrier split.
**Sources:** `04-computation/unit_distance_n22_mathieu_residue_codex.py`, `05-knowledge/results/unit_distance_n22_mathieu_residue_codex.out`, `05-knowledge/hypotheses/HYP-2572-unit-distance-n22-mathieu-residue.md`, `07-reflections/unit-distance-n22-mathieu-residue-codex.md`; HYP-2176, HYP-2188, HYP-2203, HYP-2467.

## Lead mac-mini-2026-06-17-S2: LRC regions/sections reframe → the binding-pair switch reduction (THM-524, HYP-2571)

**Status:** PROVED (sawtooth lower-envelope lemma): M(S)=max(g(½), max over pairs (a,b),k of g(k/(v_a±v_b))) — the LRC gap is a pairwise SWITCH; LRC(N) ⟺ some pair-crossing gives gap≥1/N with others clear (~78 switches for 13 runners, polynomial). Regions/SDR = the on-grid q=14 witness (blind off-grid). Complement=reversal is the one exact tournament bridge; overtaking tournament is trivially transitive (Rédei link dead). Covering hard core: M=7m/(84m+5), min 7/89.
**Next (ranked):** (1) the THM-524 NEXT step — prove "grid-sharpness (gridM=M) ⟺ a binding complement pair (v,N−v)" and that it fails EXACTLY for covering sets, merging the section lens with THM-523 (off-grid criterion). (2) Merge with codex's Hall/wall-switch program (T838) and dihedral mouth-exchange (T837/HYP-2569) — three takes on the same regions prompt; the binding-pair (mine), the Hall packets (codex), the dihedral mouths (codex) should unify on the covering hard core. (3) Bound inf M over covering sets ≥ 7/89 (the gap-side compactness frontier = THM-522's measure-side). (4) Is there a NON-trivial tournament built from the SWITCHES (crossing parities) rather than the snapshot order? — the only place a Rédei/OCF link could survive.
**Sources:** 04-computation/lrc14_{regions_sections_grounding,binding_pair_reduction,gap_M_exact}_mac-mini-2026-06-17-S2.py + workflow angle scripts (sdr_hall_symmetry, angleD_overtaking_tournament, angle_c_exotic_switches, angle_e_cross_modulus) + .out; THM-524, HYP-2571; reflection regions-not-runners-the-lonely-gap-is-a-pairwise-switch-mac-mini-S2.md; THM-523/522, codex T837/T838/HYP-2568/2569/2570.

## Lead mac-mini-2026-06-16-S3: LRC(14) GAP-SIDE reduction to covering sets (THM-523, HYP-2566/2567)

**Status:** PROVED the q-witness reduction: LRC(14) ⟺ M(S)=max_τ min_v||vτ|| ≥1/14 for all primitive COVERING SETS (a multiple of every q∈{2,…,14}); everything else is M≥1/14 via τ=1/q. No counterexample over extensive search; residual covering-set min M = 7/89≈0.0787 (10% margin). The disproof search produced the proof reduction.
**Next (ranked, to finish LRC(14)):** (1) **HYP-2566 — bound inf M over covering sets > 1/14.** Covering sets are bounded/structured (need multiples of 8,9,11,13,14 → spread out); combine scale-invariance + kind-pasteur's bounded-lcm compactness (THM-522, measure-side) on the GAP side. Is 7/89 the true inf? Search covering sets with larger speeds + 2 multiples of 14. (2) Prove the residue-refinement δ-existence lemma rigorously (the τ=a/14+δ construction for C'(14) with multiple large multiples of 14). (3) The finite-witness-cover endgame: a fixed finite family of τ's (1/q's + a/14+δ's) certifying M≥1/14 for ALL covering sets except a finite checkable residual (the destroy-lonely-points agent: 12 witnesses cover all but 12/2000). (4) Connect M-side covering reduction to L-side inf=1/1260 (HYP-2561) — both fence the same hard core.
**Sources:** 04-computation/lrc14_{gap_M_exact,constructive_dichotomy,residue_criterion,residual_hardcore_minM,strong_necessary_condition}_mac-mini-2026-06-16-S3.py + ~30 workflow angle scripts (lrc14_{disprove_*,covering_obstruction_proof,M_quantization_floor,destroy_lonely_points,family_reduction,…}) + .out; THM-523, HYP-2564/2565/2566/2567; reflection the-disproof-search-builds-the-proof-lrc14-covering-reduction-mac-mini-S3.md; THM-360/398/501/522.

## Lead mac-mini-2026-06-16-S2: LRC-14 inf L>0 — corroboration + the four closed doors (HYP-2562/2563, T834, complements kind-pasteur THM-522)

**Status:** kind-pasteur-S7 owns the headline (THM-522 scale-invariance+quantization; HYP-2561 inf=1/1260; MISTAKE-075). I independently CORROBORATED inf=1/1260 (single perturbations w≤600, doubles w≤140: nothing in (0,1/1260)) and CLOSED four positivity routes (LLL/Shearer, Selberg–Beurling, Abel/Cesàro, OCF-bridge), all blocked by ONE obstruction: L is an archimedean SIGNED singular integral with no termwise floor + positively-correlated danger events.
**Next (ranked, toward inf L>0):** (1) **The crux per THM-522/HYP-2561:** classify the (conjecturally finite, bounded-lcm) TIGHT LOCUS = LRC(14) extremal configs (Goddyn–Wong); then quantization+compactness give inf L≥1/1260. Extend my `lrc14_perturbation_inf_search` to 3-element perturbations + non-AP bases to map the full tight locus. (2) Prove the UNIFORM (core-independent) large-speed decoupling constant C (Erdős–Turán discrepancy of stranger 14m against the fixed bounded core's lonely set; constant via #core-runners + min core-gap) — the one missing ingredient for the m→∞ tail (my workflow Angle 8). (3) The Bedert level-bound |E_k∩P|≤(C log|P|)^k → Abel-controlled alternating Λ_k (OPEN-Q-097, the complementary analytic route). DO NOT re-attempt the 4 closed doors.
**Sources:** 04-computation/{lrc14_exact_rational_measure,lrc14_resonant_removal_7adic,lrc14_family_completeness_search,lrc14_perturbation_inf_search}_mac-mini-2026-06-16-S2.py + workflow angle scripts (lrc14_{7adic_archimedean_split,lll_shearer,abel_cesaro_levelmass,selberg_minorant,largespeed_decoupling,resonant_stranger_certificate,ocf_bridge_probe}) + .out; 07-reflections/why-every-lrc-positivity-method-fails-the-signed-integral-with-no-floor-mac-mini-S2.md; THM-522, HYP-2561/2562/2563, MISTAKE-075, OPEN-Q-104/097.

## Lead mac-mini-2026-06-16-S1: the octal lens on the H-spectrum + the triangular tile-count (T832, HYP-2558..2560)

**Status:** NEW probe past the kind-pasteur boat session (T830/HYP-2557, which mapped Fib/tri/square/prime + Heron — all re-verified, 3rd confirmation). The octal identity odd²=8T_k+1 aimed at H(T) (always odd, Rédei): at n≤6 the H-spectrum gaps {7,21,35,39} all AVOID residue 1 mod 8; the odd-square residue class is gap-free. Candidate necessary sieve for the forbidden-H frontier.
**Next (ranked):** (1) HYP-2558 — one Held-Karp sweep at **n=7 (and n=8 if feasible)**: is residue 1 mod 8 still gap-free? do 35,39 get realized? extend `H_mod8_octal_probe`. (2) HYP-2560 — extend the G_n/E_n invariant table to n=8 (m=21, first nontrivial Fibonacci tile-count); does anything special happen at n∈{3,4,8,12}? (3) HYP-2559 — align THM-067 c₁-vanishing Mersenne indices with the perfect-number tile-counts m=T_{2^k-1} at n=5,9,33. (4) reproduce {1,3,21,55} inside the V_4 (complement×transpose) Burnside engine (`burnside_enum_v2.c`) — a positive run upgrades the Klein-four link from suggestive to structural.
**Sources:** 04-computation/{fib_tri_square_prime_heron,H_mod8_octal_probe}_mac-mini-2026-06-16-S1.py + .out; 07-reflections/the-octal-H-spectrum-and-the-triangular-tilecount-mac-mini-S1.md; T830/HYP-2557 + triangular_fibonacci_heron_boat_kps.py (prior art); THM-462/463 (H gap-freeness), THM-485/486/067/224, T232. (arXiv:1809.09936 = uniqueness Diophantine, per T830.)

## Lead mac-mini-2026-06-15-S4: LRC-14 via the lonely measure + the Riesz-product key (THM-515, T826, HYP-2540..2543, OPEN-Q-104)

**Status:** THM-515 (theta/lonely-measure form PROVED; additive-energy predictor verified; Riesz-product route set up + feasibility-probed; closed doors confirmed). inf L>0 OPEN. Adversarial check running.
**Next (ranked):** (1) OPEN-Q-104 — build the optimized dissociated Riesz product for the j=6 interior-drop extremizer cores and push ∫M·R below 1 (the certificate); port arXiv:2511.16636's construction. (2) Bonami hypercontractive level-k bound on the cores' AP relation lattice. (3) Gaussian-subordinated Selberg minorant (finitize + theta-positivity). (4) Linnik/ternary-form reduction (Pollock-style, THM-501).
**Sources:** 04-computation/{lrc14_theta_lattice,lrc14_riesz_product}_macmini_0615s4.py + .out; THM-515; reflection lrc14-is-the-lonely-measure-and-the-key-is-a-riesz-product.md; arXiv:2511.16636, Tao 1701.02048.

## Lead mac-mini-2026-06-15-S3: FKN ⟹ Arrow — the tournament cube IS the social-choice cube (THM-512, T823, HYP-2534..2536, OPEN-Q-102)

**Status:** THM-512 proved/verified (Arrow/Condorcet bridge exact; c3=Guilbaud quadratic; Möbius sieve verified c3 n≤6). Builds on THM-511. Adversarial check running.
**Next:** (1) OPEN-Q-102 — is the OCF H=I(Ω,2) a noise-stability Stab_ρ functional (mirroring Kalai ρ=-1/3)? forbidden H ⟷ forbidden Condorcet spectra. (2) the H multi-vertex Möbius sieve (HYP-2534, depth=band-limit D). (3) Friedgut/KKL: which invariants are juntas / which arc is decisive. (4) connect Guilbaud (arcsine/level-2) to the project's pi/Wallis constants (everything-is-the-triangle).
**Sources:** 04-computation/{fkn_tiling_cube,mobius_sieve_arrow}_macmini_0615s3.py + .out; THM-512; reflection fkn-arrow-and-the-tournament-as-a-vote.md.

## Lead mac-mini-2026-06-15-S1: baby Hodge = the det/permanent (P/#P) wall through the moment problem (THM-509, HYP-2526..2528, OPEN-Q-099, T820)

**Status:** THM-509 proved (Layer-1 det-side spectral-blindness rigorous) + all holes certified moment-interior (n=6,7); holes shown to be pure integrality gaps (flag-cut refuted). Adversarially verified.
**Next steps:**
1. **OPEN-Q-099 positivstellensatz:** prove no polynomial moment inequality (any degree, spectral or overlap) cuts a baby-Hodge hole for all n — making "hole = integer lattice point interior to the flag-feasible body, skipped by #P" a theorem.
2. **Generating function for the non-spectral defect:** does the necklace/zeta moment-cumulant pair (HYP-2526) give a clean handle on the non-spectral defect dimension A000009(n)-3 (THM-505)? The Witt/zeta Π(1-u^k)^{-W_k} and the restricted-partition Π_{k odd≥3}1/(1-x^k) are both Euler products — connect them.
3. **The c3-fiber hole structure:** the regular/near-regular score class is the richest hole source (c3=8 at n=6 → 3 holes; c3=14 at n=7 → 12). Characterize the holes as a function of the score sequence / the extremal (densest-cycle) fiber.
4. **Cross-link to d⊥H and the BSD/Hodge self-dual-middle reflection:** the det/per wall unifies them (d det-side, H per-side).
**Sources:** 04-computation/{baby_hodge_moment_region,baby_hodge_convex_certificate}_macmini_0615s1.py + .out; THM-509; reflection baby-hodge-is-the-det-permanent-wall-read-through-the-moment-problem.md.

## Lead mac-mini-2026-06-14-S1: dual-tower (self-dual) + LRC(14) singular-series structure (THM-503) + the Pascal-slope-d Pisot tower (T818/T819, HYP-2520..2525, OPEN-Q-097/098)

**Status:** session complete; THM-503 adversarially verified + adopted by the mesh (THM-504 builds on it); dual-tower answered; Pascal family identified+extended.
**Highest-value next steps:**
1. **OPEN-Q-097 (the LRC(14) prize):** prove inf L>0 over the dilated-AP cores. Now reframed (THM-503 + the concurrent THM-504) as a CROSS-LEVEL Abel-summation bound on the conditionally-convergent alternating-in-|T| series L=(6/7)^13+(6/7)^11 C_2−(6/7)^10 Σ_3+… with a Pólya–Vinogradov/Weil bound on each convergent signed sinc-lattice level Σ_k. Joint target with kind-pasteur's PZ thread.
2. **HYP-2409 indecomposability:** prove the skew tower's row code stays the INDECOMPOSABLE d+ (not e8⊕e8) for all k — the weight-4 support graph stays connected under doubling [[H,H],[−Hᵀ,Hᵀ]]. Self-dual Type II half is provable + verified to order 64; indecomposability is the deep part (rides SO(32) forever).
3. **OPEN-Q-098 (gap-d tournament realization):** does a_d(n) count a natural gap-d tournament/staircase family for d≥3? Prove the realized 2·Fib(m-2) circular-tournament count (d=2, S518). Define the "d-graded metagraph" whose H-level sizes are the slope-d ridge.
4. **The plastic-number bridge:** d=5 of the family (x⁵=x⁴+1) shares its root with Padovan (x³=x+1, the monad/free-factorial THM-438 thread). Is there a slope-5 ↔ monad tournament bridge, or pure algebraic coincidence?
**Sources:** 04-computation/{lrc14_singular_series_adelic,skew_tower_selfdual,pascal_slope_d_family}_macmini_0614s1.py + .out; THM-503; reflections the-dual-of-the-skew-sylvester-tower-is-itself.md, the-pascal-slope-d-family-and-its-pisot-towers.md.

## Lead codex-2026-06-13: tournament trace speedups and first overlap correction (HYP-2498/T817/OPEN-Q-093)

**Done:** added `04-computation/tournament_structure_speedup_patterns_codex.py` and stored output. The script turns the "efficiency becomes proof" theme into a trace-correction boundary. It verifies `c_k=tr(A^k)/k` for `k=3,4,5`, then proves/tests the first correction `tr(A^6)=6*c6+3*c3+6*p33_meet`, where `p33_meet` is the number of unordered pairs of distinct directed triangles with nonempty intersection. The naive midpoint-return correction `sum_v(A^3_vv)^2` is explicitly shown insufficient. Validation: exhaustive `n=3..6` and random `n=7..9` samples have zero mismatches; fixed random `n=14` corrected `c6` matches brute and is about `106x` faster. Exhaustive `n=6` information tournament over `score,c3,c4,c5,c6,H` is transitive with champion order `H>score>c4>c5>c3>c6`; bucket audit shows `(c3,c4,c5,c6)` determines `H` at `n=6`, while `score+c5+c6` still leaves one mixed H bucket. Incoming S5/THM-499 sharpens the meaning: `H=1+2(c3+c5)+4D` uses the disjoint-triangle-pair count `D=alpha_2`, and this trace correction uses its complement `p33_meet=C(c3,2)-D`.
**Next:** enumerate support types for `tr(A^7)` and `tr(A^8)`, separating scalar corrections from placement-sensitive corrections. Then run the same information-tournament and bucket-mixing audit at `n=7`, ideally over isomorphism classes or an incremental labelled sweep, to test whether corrected trace vectors continue to compress `H`.

## Lead codex-2026-06-13: Pollock Sierpinski carry-pair lift (HYP-2497/T816/OPEN-Q-092)

**Done:** added `04-computation/pollock_sierpinski_carry_scout_codex.py` and stored output. The script extends HYP-2491 by testing the Sierpinski/Waring dyadic analogy. Result: pure single-residue dyadic obstruction is unlikely, because tetrahedral atoms hit every residue class modulo `2^e` for `1<=e<=12` in the scan. Lucas parity remains exact (`Te_k` odd iff `k=1 mod 4`), but it is not a missing residue. After lifting to defect pairs `r,r+tri(k) in D_4`, dyadic compression becomes informative: for `k>=100`, observed pair classes stabilize at `168` by `2^8` while possible pair classes grow as `4^e`; the dyadic compression tournament is transitive `12>11>...>3`, no directed 3-cycles, one Hamiltonian path. Carry-window data shows `85/241` four-defects are within `100` below a tetrahedral number and `240/241` within `5000`, but the largest defect `343867` is `5637` below `Te_127`, so carry windows are a diagnostic, not the whole proof.
**Next:** prove the 2-adic surjectivity lemma `{Te_k mod 2^e}=Z/2^eZ` as a formal anti-obstruction. Then build pair-address feasibility tests that combine `(r mod 2^e, k mod 2^e)`, triangular gap carries, and four-tetrahedron convolution lifts. The target remains HYP-2491's no-long-pair theorem for `k>825`, paired with the finite width-3 stencil through `k<=825`.

## Lead codex-2026-06-13: Pollock tetrahedral defect-pair descent (HYP-2491/T815/OPEN-Q-091)

**Done:** added `04-computation/pollock_tetrahedral_defect_descent_codex.py` and stored output. The scout treats Pollock's tetrahedral conjecture as a four-defect problem. Let `D_4` be numbers not in the at-most-four tetrahedral sumset. Through `10^6`, the computation finds exactly `241` defects, largest `343867`, and verifies no misses for at most five tetrahedra. In shell `[Te_k,Te_{k+1})`, the one-back anchor fails exactly when `r` and `r+tri(k)` are both in `D_4`; among known defects the last such triangular pair is `3142 -> 343867 = 3142 + tri(825)`. Shell stencil audit through `k<=1200` shows offsets `0..3` cover every shell, offsets `0..1` suffice after `k=825`, and the anchor-offset tournament is transitive `3>2>1>0`.
**Next:** prove the strong tail `D_4 subset [1,343867]`, or the weaker no-long-triangular-self-correlation lemma for `D_4`. Build a pair-residue scout: for moduli `m`, study observed pairs `(d,d+tri(k))` rather than single defects, because single four-defects have no obvious small local obstruction. Then convert the width-3 shell stencil up to `825` into a compact finite certificate.

## Lead codex-2026-06-13: LRC A000568 source fiber (HYP-2486/T814/OPEN-Q-090)

**Done:** added `04-computation/lrc_a000568_source_fiber_codex.py` and stored output. The script isolates the exact A000568 layer hidden in LRC: add the observer by threshold edges `0 -> i iff ||v_i t|| >= 1/N`; then LRC-good is exactly observer-source. At good states, the rooted class is `source_cone(deleted runner class)`. Canonical enumeration through `m=6` verifies the source-cone bijection with A000568 (zero collisions/deletion failures). Exact small-clock audits show moving-runner classes are mixed but rooted source fibers are pure (`rooted_mixed=0`, `cone_exact=True`). LRC14 snapshots for AP13, one-stranger-611, the two HYP-2470 exceptions, and the THM-497 band-2 refuter all have observer outdegree `13` at first witness and leave 13-vertex runner fingerprints.
**Next:** attach deleted A000568 fingerprints to the existing Q27/Q31/band-cover ledgers. Test whether hard rows with the same deleted fingerprint share active Q31 constraints, or whether the multiplicative cover fiber separates them. The proof target is a dichotomy: the blocked-band walk enters a source-cone class, or its avoidance forces balanced-cover congruences feeding HYP-2471/HYP-2480.

## Lead codex-2026-06-13: LRC14 blocking-height dominance atlas (HYP-2481/T811/OPEN-Q-089)

**Done:** added `04-computation/lrc14_blocking_height_dominance_codex.py` and stored output. The script imports the dilated-band viewpoint from the recent Q31/KPS work: `h(S)` is the first shell with an uncovered unit, and pre-height shells `14<=q<h(S)` are all fully blocked. It orients a speed tournament by cumulative cover-mask dominance across those pre-height shells. The prompt's dominance question is now nuanced: raw cumulative mean pair margin grows with blocking height (`0.779` Pearson in the one-stranger family, `0.942` in random primitive rows), but per-shell/per-unit normalized dominance falls (`-0.711` and `-0.729` for normalized margin). Every named hard packet has a transitive speed tournament (`score_hist=0..12`, no directed 3-cycles, one Hamiltonian path), so the binary speed quotient is a pipeline shadow rather than the final proof object.
**Next:** prove or refute the dichotomy. Branch A: a high cumulative/private-load carrier is peelable, yielding a witness, Bprime opening, or lower-core descent. Branch B: no peelable carrier exists, so balanced-cover congruences should force the Q31/band-2 ramified portal already isolated by HYP-2471/HYP-2480. Add leave-one-out support-criticality, first-leak deficit motion, and unit-obligation/shell-vertex tournaments; then feed the resulting typed budget into the below-eight-core `e=5` search.

## Lead codex-2026-06-13: irreducibility tricks as LRC14 ramified local gates (HYP-2480/T810/OPEN-Q-088)

**Done:** added `04-computation/irreducibility_tricks_lrc_transfer_codex.py` and stored output. The script turns polynomial irreducibility tactics into LRC14 proof carriers: primitive normal form, mod-p residue blockers, Eisenstein/Newton valuation repairs, Singh/Cohn factor-capture budgets, Hensel/recombination ledgers, and Cohn/Perron dominance. The concrete diagnostic is the two HYP-2470 Q27-feasible four-deletion packets. Both have `12/13` speeds divisible by `7`; each has exactly one non-7 primitive escape; both primitive escapes are divisible by `13` (`936=2^3*3^2*13` and `1066=2*13*41`). They open at the missing plain shells `q=33` and `q=31`, with Bprime(any) and positive exact safe measure. Proof-carrier Tournament Analysis is transitive with leader `integral_convolution_lift_ilp`, reflecting a pipeline rather than competing heuristics.
**Next:** prove the ramified 7-ideal/13-clock portal lemma: any four-deletion Q27 packet with this valuation shape opens at `q in {31,33,41}` or by Bprime. Extract dual/Farkas certificates from the HYP-2465/HYP-2470 MILPs so Q27 infeasibility becomes a human-readable local-prime certificate. Then build a below-eight-core survivor ledger with shell-27 class, 13-clock debt, divisor fiber, support-load, and owner/Bprime channels; add a Cohn/Perron outside-window normalizer for speeds beyond `1092`.

## Lead codex-2026-06-13: LRC14 Q31 fiber repair for the eight-core exceptions (HYP-2471/T812/OPEN-Q-087)

**Done:** preserved the HYP-2471 addendum alongside HYP-2470/HYP-2480. The full delete-four/budget-five carry-window scan over `495` deletion shapes gives `489` Q27-infeasible rows, `2` Q27-feasible rows, and `4` short-timeout rows later resolved Q27-infeasible. The two Q27-feasible rows are exactly HYP-2470's exceptions, `(28,42,56,84)` and `(42,56,70,84)`. The direct theorem uses plain shells through `41`, but the fibered repair is sharper structurally: both exceptions become infeasible over `Q31={d*m:d in {1,2,7,14},m<=31}\{1}`.
**Next:** extract dual/Farkas or hitting-set certificates for the two Q31 infeasibilities and compare their active constraints with HYP-2480's 7-ideal/13-clock diagnosis. Then test whether the same Q31 certificate pattern survives when five or more core speeds are deleted, or whether the proof must switch to owner/Bprime and outside-window normalization.

## Lead codex-2026-06-13: LRC14 dilated-band covering correction (THM-497/T813/OPEN-Q-087)

**Done:** integrated kps1's covering reformulation and correction. THM-497 writes a plain shell witness as an uncovered unit in `(Z/q)^*` outside the `13` dilated danger bands of the row. The band-0 divisibility lemma is exact, but the cardinality side goes the wrong way for a pure counting proof: `26k` unit residues can in principle be covered by roughly `14k` band positions. Companion scouts show the danger of scalar ceilings: a primitive non-Bprime-dominant row blocks all plain shells `q<=41`, first opens at `q=43`, and has first Q27/Q31 witness `(44,7)`; resource climbing finds non-dominant full blockers through `K=55` with first witnesses `29,43,56`, then no greedy full blocker at `K=69,83`.
**Next:** turn "structure forbids core" into a theorem. Candidate retained channels are 7-adic occupancy, 13-clock primitive escape, divisor-fiber address, owner/Bprime support, and support-load dual certificates. Treat the covering model as a search/pruning oracle, not as a standalone scalar proof.

## Lead codex-2026-06-13: LRC14 eight-core shell-41 exception gate (HYP-2470/T809/OPEN-Q-087)

**Done:** pushed the first below-nine-core finite boundary past HYP-2465. Added `04-computation/lrc14_eight_core_q27_setcover_codex.py` plus stored output. Q27-only set-cover over all `binom(12,4)=495` four-core deletion sets gives `493` infeasible deletion addresses and exactly two Q27-feasible addresses after repairing `12` sparse short-cap unknowns. The two addresses are `(28,42,56,84)` and `(42,56,70,84)`. Example Q27 packets open at plain `q=33` and `q=31`, respectively, and both have Bprime(any) plus positive exact safe-measure certificates. Adding the missing plain shells through `41` makes both exceptional addresses infeasible. Therefore every primitive carry-window row retaining at least eight core speeds has either a Q27 witness or a plain witness `q<=41`.
**Next:** convert the MILP certificate into a dual/combinatorial proof: high-obligation cases should yield a uniform coverage deficit, while sparse-address repairs need their own typed certificate. Then attack the true remaining portal: no-Q27/no-plain-41 rows must delete at least five core speeds or exit the carry window. Build a below-eight-core typed budget with separate 13-clock, divisor-fiber, shell-27, owner/Bprime, and low-clock resources.

## Lead codex-2026-06-13: unit-distance small-factor resonance capacity gap (HYP-2467/T807/OPEN-Q-085)

**Done:** turned THM-493's resonance bonus into an exact connected-factor capacity atlas. Added `04-computation/unit_distance_resonance_capacity_atlas_codex.py` and stored output. The script enumerates every connected triangular-lattice patch through size `9` up to translation and `D6`, computes edge counts and norm-`t` displacement spectra, and maximizes the non-degenerate `t>1` resonance bonus over every relative `D6` orientation. Exact result: `27=3*9` maxes at `75<81`; `28=4*7` reaches `85>84`; `30=5*6` ties; `32=4*8` crosses. The size-3 stress test is the useful proof object: `K3` has generic `75` but zero bonus, while the resonance-bearing 3-point paths reach only `69/70`.
**Next:** prove the size-3 capacity inequality without enumeration, then build a Moser-patch compression verifier: given a dense rank-4 patch, extract displacement packets and try to factor them through small connected triangular shadows. The target theorem is that any 27-point 82-edge patch either compresses to an impossible capacity lane or exposes a genuinely irreducible obstruction worth classifying directly.
## Lead codex-2026-06-13: LRC14 Church-Frobenius descent upgrade (HYP-2469/T808/OPEN-Q-086)

**Done:** read Church's arXiv:2508.14876 through the existing HYP-2445 bridge and reprocessed it against the newer HYP-2463/HYP-2464/HYP-2465 LRC14 finite atlases. Added `04-computation/lrc14_church_frobenius_descent_codex.py` plus stored output, HYP-2469, and reflection `lrc14-church-frobenius-descent-upgrade.md`. The key upgrade: import Church's proof grammar, not just its `1092/91/78` arithmetic. In the paper, scalar Shioda supersingularity is too weak; diagonal forms on every asymmetric partial Frobenius twist survive; non-exceptional curves descend by projection-degree drop. LRC14 analogue: raw plain-shell blocking is the scalar shadow; Q27 obligations/resource coordinates are the retained support channel; a no-Q27 row must either enter a certified finite atlas, open a named exception, or descend. Certified blocks now cover one-stranger (`936/936`), hard replacement hull (`77520/77520`), two-stranger plain residuals (`877/877`), and near-core set-cover through `|D|<=3` (`299/299`). Tournament Analysis ranks `lrc14_Q27_obligation_setcover` first, Church's route second, with `4/28` edge flips versus scalar-only ranking.
**Next:** integrate HYP-2470's correction into the descent theorem: the `|D|=4` portal is closed once plain shells through `41` are admitted, so a no-Q27/no-plain-41 row must either delete at least five core speeds or leave the carry window. First, build the below-eight-core typed set-cover scout with separate budgets for 13-clock, divisor-fiber, shell-27, owner, and low-clock roles. Second, prove/test an outside-window normalizer: speeds beyond `1092` should open Bprime(any), dominate/transport an existing core speed, or reduce through a divisor/carry fiber without losing blockedness. Also formalize the exception catalogue: AP, Vstar, nonprimitive 2AP, `q=91`, `q=161`, owner-private/Bprime, and low-clock exits.

## Lead codex-2026-06-13: LRC14 two-stranger compression stress (HYP-2464/T805/OPEN-Q-083)

**Done:** extended the HYP-2463 Q27 resource-independence route beyond the old hard-residue list. Added `04-computation/lrc14_two_stranger_compression_stress_codex.py` and stored output. The script deletes one runner from `7*{1,...,12}` and adds any two distinct non-core speeds up to `13*84`, scanning `6,868,368` primitive true two-stranger rows by bitset safe-twist masks. Only `877` rows block every plain shell `q<=27`; all `877` have a Q27 witness, so there are `0` Q27 misses. The important correction to HYP-2463 is that `636/877` residuals use zero old hard residues, but every residual still has at least one added speed divisible by `13`, no residual deletes `7,21,49`, and the late rescues are divisor fibers (`70,84,91` plus one `161=7*23`). Compression-map Tournament Analysis is nontransitive with one directed 3-cycle and leader `divisor_fiber_Q27`.
**Next:** prove the resource-coordinate compression lemma: arbitrary primitive rows that block the plain shell should compress to `13`-clock debt + deleted-core address + shell-27 pair class + divisor fiber, or else open a low clock, AP/Vstar/2AP descent, or odd owner/Bprime channel. Upgrade the bounded `r<=13*84` statement using a fast Bprime/large-speed escape, and build a faster Bprime certificate engine for the `877` residuals.
## Lead codex-2026-06-13: LRC14 near-core Q27 set-cover compression (HYP-2465/T806/OPEN-Q-084)

**Done:** strengthened HYP-2463 from the eight named hard residues to arbitrary bounded replacements. Added `04-computation/lrc14_near_core_q27_setcover_codex.py` plus stored output. For `B=CORE\D`, the script builds the Q27 obligation hypergraph of twists safe for `B`; added speeds cover obligations where they are dangerous. In the HYP-2444 carry window `1..1092`, primitive binary set-cover MILPs with add budget `|D|+1` are infeasible for all deletion sets through `|D|=3`: `1/1`, `12/12`, `66/66`, and `220/220`, zero feasible/unknown. Thus every primitive bounded replacement row retaining at least 9 of 12 core speeds has a Q27 witness. A direct one-deletion/two-add scan finds `877` plain `q<=27` blockers but `0` Q27 misses, so plain-shell residue names are noisy and Q27 set-cover number is the stable object.
**Next:** prove a global compression/descent theorem. Any LRC14 Q27 blocker must either normalize into the carry window and delete at least four core speeds, or leave the carry window and trigger divisor-fiber/owner/Bprime opening. Then analyze below-nine-core rows: do four core deletions force low clocks, AP/Vstar/2AP descent, support-load contradiction, or a reusable MILP lower bound with a new carrier?

## Lead codex-2026-06-13: LRC14 hard resources do not stack (HYP-2463/T804/OPEN-Q-082)

**Done:** implemented the parity-typed Q27 hard-resource ledger for LRC14. Added `04-computation/lrc14_parity_typed_q27_ledger_codex.py` and stored output. The eight HYP-2444 one-stranger plain-shell blockers are treated as resource atoms inside the `7*{1,...,12}` core. With bitset safe-twist masks, the script scans the complete hard replacement hull `sum_k binom(8,k)binom(12,k-1)=77520`: every row has a Q27 witness, no Q27 misses. Only ten rows miss plain `q<=27`: the original eight one-stranger rows plus `delete (28), add (351,1053)` caught by `q=30`, and `delete (28,63), add (351,962,1053)` caught by `q=34`. The only `q=91` rows are the original `611,702` packets. Proof-obligation Tournament Analysis is transitive with leader `typed_Q27_ledger`.
**Next:** prove the compression/resource-independence theorem. Show any primitive LRC14 row with no Q27 witness can be parity-typed and compressed to the hard replacement hull without losing blockedness, unless it opens a low clock, divisor-fiber witness, AP/Vstar/2AP descent, or odd owner/Bprime deletion. First technical sublemma: any row blocking all plain `q<=27` shells has a marked subpacket projecting to shell-27 class `0` or `+-10` plus the 13-clock; second sublemma: copying that packet forces enough 7-core support loss to expose Q27.

## Lead codex-2026-06-13: parity projector channel gate (HYP-2459/T803/OPEN-Q-081)

**Done:** turned the prompt's midpoint/reversal slogan into an exact projector calculus. Added `04-computation/parity_projector_channel_atlas_codex.py` and stored output. The scalar side records that pair-differences around a midpoint keep odd offset powers, while the Faulhaber interval balance also has one fixed central atom `c^p`. The tournament side encodes arcs by signs; converse is global sign reversal, so invariant functions have even Walsh support and anti-invariant functions have odd support. Exact labelled audits for `n=3..5` verify `H,c3` even-Walsh, writhe odd-Walsh, rooted `start0` mixed but `start0+end0` even and `start0-end0` odd, raw `H` flip delta even, and oriented `H` gradient odd. Proof-carrier TA is transitive with leader `lrc14_q27_owner_carry_ledger`.
**Next:** implement a parity-typed LRC14 Q27 ledger. Fields should declare `even_scalar`, `odd_marked`, `transported`, or `compatibility_packet`; AP/Vstar/2AP and HYP-2444 one-stranger rows are the first regression cases. Quotient on even scalar clocks only after splitting transported source/sink or start/end fields into sum/difference; then use odd owner/carry/deletion channels to force strict witnesses, descent to known wall atoms, or owner-private openings.

## Lead codex-2026-06-13: Faulhaber anchor expansion (HYP-2457/T801/OPEN-Q-079)

**Done:** sharpened the HYP-2454 Bernoulli/Faulhaber route for the user's power-balance anchor. Added `04-computation/triangular_faulhaber_anchor_expansion_codex.py` and stored output. With midpoint `c=a+n`, the exact defect is `D_p(c,n)=c^p-2*sum_{r odd} binom(p,r)c^(p-r)S_r(n)`, so only odd Faulhaber moments survive. Writing `u=n(n+1)`, the formal root expansion is `c=p*u+alpha_p+beta_p/u+gamma_p/u^2+...`; `alpha_p`, the user's `beta_p`, and the factored `gamma_p` all carry `(p-1)(p-2)`, explaining exact p=1/p=2 towers. The p=2 face is the square-pyramidal cuboid identity `6*sum_{j<=n}j^2=n(n+1)(2n+1)=2*S1`. Tournament Analysis over proof carriers is transitive with leader `odd_moment_projection`.
**Next:** prove the fixed-`p` asymptotic with a uniform remainder, then use it or direct odd-moment inequalities to prove HYP-2454's bracket `D_p(p*n(n+1),n)<0<D_p(p*n(n+1)+1,n)` for all `p>=3`. Transfer target: build an LRC14 odd-wall ledger whose fields include blocked twists, owner support, shell-27 class, divisor fiber, carry residue, endpoint atom, and moment/resource defect rather than scalar "q blocked" alone.

## Lead codex-2026-06-13: boundary-lift irreducibility transfer (HYP-2455/T799/OPEN-Q-077)

**Done:** merged recent agent work into a scalar-shadow/hidden-lift synthesis. Added `04-computation/boundary_lift_analogy_atlas_codex.py` and stored output. The atlas treats polynomial convolution factor grids, LRC Q27/Pisano support ledgers, unit-distance non-product Moser fibers, triangular moment/fractional addresses, p-curvature operator ledgers, product-quotient diagonal gates, and `[72,36,16]` support-design incidence as instances of one proof interface. Tournament Analysis over carrier/proof-obligation vertices is nontransitive with `3` directed 3-cycles, SCC sizes `[5,1,1]`, `9` Hamiltonian paths, and leader `polynomial_convolution_lift`.
**Next:** implement a common lift-feasibility schema; extend HYP-2452 to bounded degree-6 ILP/SAT convolution lifts; turn LRC multi-stranger rows into allocation ledgers; split unit-distance `N=27/28` candidates into product-reducible and Moser-irreducible fibers; encode `[72,36,16]` minimum-word support incidence over the `78/90` address.

## Lead codex-2026-06-13: Faulhaber odd-moment OCF bridge (HYP-2458/T802/OPEN-Q-080)

**Done:** added an OCF compatibility-packet addendum to the newly landed HYP-2457 anchor expansion. Added `04-computation/faulhaber_odd_moment_ocf_bridge_codex.py` and stored output. The script rechecks the odd-moment identity and p=1/p=2 anchors, then separates the p=2 square-pyramid cuboid `6*S_2(n)=n(n+1)(2n+1)` from the antisymmetric balance layer, which uses only odd `S_1`. Main new point: OCF supplies the lift template, because odd atom inventory is not enough without compatibility packets `alpha_k`. Carrier Tournament Analysis over proof carriers has `8` directed 3-cycles, SCC sizes `[6,1,1]`, and `45` Hamiltonian paths.
**Next:** build a finite odd-moment compatibility object whose one-particle shadow is `S_1,S_3,...`; test it on HYP-2456 Beatty/Pell boundary atoms, HYP-2444 Q27 owner-support ledgers, and `[72,36,16]` minimum-word support packets; compare whether `alpha_k`-style compatibility predicts LRC strict/wall/open status better than raw moment or residue scalars.
## Lead codex-2026-06-13: Beatty-Pell crossover word (HYP-2456/T800/OPEN-Q-078)

**Done:** turned the HYP-2453 side-containment word into an explicit Beatty-Pell normal form, as a concrete exact instance of HYP-2455's boundary-lift grammar. Added `04-computation/triangular_tower_beatty_pell_decomposition_codex.py` and stored output. For `B_m.L`, the square-shell address is exactly `n_m=floor(sqrt(2m^2+m))=floor(sqrt(2)(m+1/4))`, with `d_m=n_m-m` and Pell/carry remainder `r_m=m^2-2m*d_m-d_m^2`; integer inequalities in `(d_m,r_m)` determine both `B_m.L` and `B_m.R` states and match the direct floor-sqrt classifier through `m=250000`. The `d_m` increment word is Sturmian, while the visible two-side token word is a six-interval rotation/carry decoration with zero-density Pell atoms `LR` and `RL`. Tournament Analysis ranks `exact_carry_normal_form` first.
**Next:** prove the six-interval coding from the exact inequalities and equidistribution of `{sqrt(2)(m+1/4)}`; build the LRC14 Q27 analogue with fields `(q, shell class, unit quotient class, divisor fiber, owner support, carry residue, endpoint atom, deletion/opening target)`; compare AP/Vstar/2AP and HYP-2241 owner-private deletion rows to the `LR/RL` wall-atom grammar.

## Lead codex-2026-06-12: triangular power-balance towers (HYP-2454/T798/OPEN-Q-076)

**Done:** integrated the user's two triangular towers as exact interval power balances. Added `04-computation/triangular_power_balance_towers_codex.py` and stored output. The ordinary tower is `D_1(C,n)=0` at `C=2T_n`; the square tower is `D_2(C,n)=0` at `C=4T_n`. The first tower partitions square shells `[n^2,(n+1)^2-1]`; the second tower has square equality but ordinary defect `L2-R2=2T_n` and total `L2+R2=4S1`. Crossovers reproduce the user's examples and find the unique checked full side equality `Q_L(3)=[21,22,23,24]=F_R(4)`. At the same row, ordinary shadows are `90=S1(4)` and `78=C(13,2)`, matching the existing `[72,36,16]` `lambda_5` beacon. Tournament Analysis is nontransitive and ranks the `78/90` support shadow first.
**Next:** prove the general p>=3 bracket `D_p(2pT_n,n)<0<D_p(2pT_n+1,n)` or find an exception; solve the Pell-style endpoint alignment families; test whether the `78/90` shadow becomes a concrete LRC14/code support-incidence constraint rather than scalar numerology.

## Lead codex-2026-06-12: convolution-lift split survivors (HYP-2451/T795/OPEN-Q-073)

**Done:** extended HYP-2449 by turning reducibility into the hidden 2D lift `a_k=sum_{i+j=k}b_i c_j`. Added `04-computation/convolution_lift_irreducibility_carrier_codex.py` and stored output. Exact brute-force residue lifts agree with symbolic mod-p split survivors in small examples. In the degree-4 coefficient scout, least mod-p convolution blockers through `31` certify `3058/3096` irreducibles (`98.77%`) with zero false positives and cut sign-bucket mixing from `16` to `8` mixed `signs+least_blocker` buckets. Newton examples show the complementary face: Eisenstein-style rows can look reducible mod p while a one-edge p-adic lower hull proves irreducibility.
**Next:** add Newton/valuation certificates to the `38` no-small-blocker irreducibles; extend split-survivor signatures to degrees `5/6`; attach Singh-depth and Cohn-depth only after cheap residue/valuation gates; translate survivor ledgers to LRC14 Q27 denominator/resource fibers.
## Lead codex-2026-06-12: convolution factor-capture tiling (HYP-2452/T796/OPEN-Q-074)

**Done:** extended HYP-2450/HYP-2451 from coefficient and residue/valuation split-survivor quotients to an integer factorization lift with value pruning. Added `04-computation/convolution_factor_capture_tiling_codex.py` and stored output. Reducibility is encoded as an integral convolution grid with diagonal sums `a_k=sum_{i+j=k}b_i c_j`; irreducibility means no nontrivial lift. The exact primitive degree-4/5 oracle agrees with Sympy on all tested rows: `3856` degree-4 rows (`792` reducible) and `2016` degree-5 rows (`488` reducible), with zero mismatches. Added factor-capture witness scores from `Omega(f(m))`, residue-class tournaments, sign-cube chamber summaries, and a proof-route tournament whose nontrivial SCC puts `convolution_lift_disprover`, `factor_capture_hypertournament`, and `sign_cube_chamber_tournament` in live tension.
**Next:** build bounded ILP/SAT/SMT convolution feasibility for degree `>=6`, using the exact degree-4/5 solver as a regression oracle; add Newton-polytope boundary-layer constraints for sparse/multivariate rows; transfer the boundary-total/hidden-lift grammar to LRC14 blocker ledgers and `[72,36,16]` support/design incidence.

## Lead codex-2026-06-12: coefficient-tiling prime/irreducible carrier (HYP-2449/T793/OPEN-Q-071)

**Done:** extended HYP-2447/HYP-2448 with the user's coefficient tiling model.  The degree-5 row sizes are literally fixed-path skip rows: `a5` apex down to `a1`, and the stronger `constant_spine` model places `a0` on the Hamiltonian-path row.  Added a degree-4 exhaustive sign+magnitude scout (`3888` rows): every bare unmarked coefficient-tournament quotient is mixed for irreducibility, marked signs remain mixed, and `marked_signs+local_zero_primes` separates fixed-divisor obstruction in the scout.  Cohn repunit rows show sign-only transitive tournaments can be reducible or irreducible depending on place-value address.
**Next:** build exact p-adic Newton-row tournaments for Eisenstein/Dumas/Perron irreducibility criteria; extend the coefficient sweep to degree 5/6; turn Cohn base-`b` addresses into weighted skip-row tournaments; transfer fixed-divisor row detection to LRC14 Q27 resource rows.

## Lead codex-2026-06-12: irreducibility-prime tournament prism (HYP-2447/T791)

**Done:** merged Bunyakovsky/Cohn/Singh/Iravanian with the repo's Heegner/THM-410 prime-polynomial horizon. The new atlas treats prime values as finite irreducibility witnesses, fixed divisors as local obstructions, Cohn digit polynomials as prime+address -> irreducible, and real-factor recombination as subset-sum over looser atoms.
**Next:** build an LRC14 fixed-divisor analogue: identify runner/denominator atoms that block every lift the way `x^2+x+2` is always even. Then run a recombination audit over HYP-2443 blocker atoms across maps `q->2q`, `q->7q`, and `27->9->3`.

## Lead codex-2026-06-12: Grothendieck-Katz p-curvature ledger for LRC14 (HYP-2446/T790)

**Done:** merged the Grothendieck-Katz p-curvature conjecture as an operator/carry local-global template for LRC14. Added a toy local p-jet atlas showing scalar mod-p shadows can lie both ways: `1/(1-z)` has zero operator p-curvature but full naive scalar rank, while `z/(1-z)` has full operator p-curvature but zero naive scalar rank. The LRC14 transfer is the denominator-curvature ledger `K_q(S)=(blocked twists,tau_q,Pisano class,13-clock,divisor fiber,Bprime/owner target)`.
**Next:** compute support-transport defects across denominator maps (`q -> 7q`, `q -> 2q`, `q -> d*m`, `27 -> 9 -> 3`) for the HYP-2443 high-pressure rows. Prove that zero defect over a long ladder forces AP/Vstar/2AP or descent, while positive defect gives a finite witness or owner-private opening.
## Lead codex-2026-06-12: irreducible-prime certificate-state addendum (HYP-2448/T792/OPEN-Q-070)

**Done:** merged the user's Bunyakovsky/Buniakowski bidirectional atom lens with Singh arXiv:2411.18366 and Iravanian arXiv:2410.15880, then rebased over the newly landed HYP-2447 prism and kept this as an addendum. Added `04-computation/irreducible_prime_carrier_tournament_codex.py` and stored output. The atlas separates forward prime production, reverse Singh/Murty/Cohn certificates, fixed-divisor obstructions, and real-factor recombination survivor tests. Key small facts: `x^2+x+2` is irreducible but fixed-divisor-blocked; `9841` in base `3` gives `1+x+...+x^8` with factor degrees `[2,6]` and `Omega(9841)=2`; `x^4-10x^2+1` is irreducible but has two false integer-trace subset candidates. Proof-carrier TA is nontransitive with `3` directed cycles and SCC sizes `[5,1,1,1]`.
**Next:** replace the floating real-root trace scout by exact algebraic trace lattices; build certificate states `C(f;X)` for a larger polynomial family and measure edge flips as `X` grows; transfer the same retained-certificate pattern to LRC14 Q27 resource vectors and to `[72,36,16]` support/matroid/design construction moves.

## Lead codex-2026-06-12: coefficient tiling prime bridge (HYP-2450/T794/OPEN-Q-072)

**Done:** turned the user's triangular coefficient/tile picture into a diagonal quotient of fixed-Hamiltonian-path tournament tilings. Added `04-computation/coefficient_tiling_prime_bridge_codex.py` and stored output. Counts `c_d` on gap-`d` diagonals are Cohn digits; centered magnetizations `A_d=2c_d-(N-d)` are coefficient signs/magnitudes. Full grids through `N=7` have zero positive-degree Cohn-prime irreducibility mismatches. Fixed-path `N=6` has 120 coefficient profiles over 1024 tilings, 57 positive-degree Cohn-prime profiles, 96 centered-irreducible profiles, and 91 profiles with hidden `H` variation, max spread 34.
**Next:** prove the exact fixed-path Cohn-diagonal lemma; classify magnetization magnitude slices, especially the `N=6` parity-minimum slice `(1,0,1,0,1)` where all 8 distinct polynomials are irreducible in the pilot; add SCC/cycle/Hamiltonian-path fiber fingerprints; translate the same profile/fiber split to LRC14 resource vectors and `[72,36,16]` support/matroid/design moves.

## Lead codex-2026-06-12: LRC14 marked ladder support gate (HYP-2443/T787)

**Done:** added a marked blocker-hypergraph atlas for HYP-2438. For each denominator `q`, compute which runners block each unit twist and the minimum support `tau_q(S)`. This separates pure shell band-1 failures from fibered ladder failures: the single-stranger `S(r)` rows fail pure `q<=27` exactly in the `13|r`, `r mod 27 in {0,+/-10}` pattern, while fibered addresses such as `q=91=7*13` or rung-up shell `q=40/41` catch them.
**Next:** prove the support-pressure dichotomy: either a finite ladder witness exists, or the blocker ledger has a universal/apex blocker reducible to HYP-2256, or repeated cover-load marks a runner for Bprime(any runner)/owner-private deletion. Exclude AP/Vstar/2AP as normalized wall atoms before applying this to primitive loose targets.

## Lead kind-pasteur-2026-06-11-S3: the pentagonal product is a hub — random-sign Lyapunov γ_pent, Euler-sign rigidity, the η²⁴ code-discriminant bridge, [72,36,16] localized (THM-487/488/489, HYP-2417..2423)

**Done:** THM-488 (γ_pent≈0.206 new Lyapunov constant, validated vs Viswanath; Euler's signs the UNIQUE subexponential pattern = analytic shadow of the pentagonal #thm; IVT half proved, hard half certified on 1585 sets via argument principle). THM-489 (code discriminant P₂₄ = 16η²⁴ exactly; extremal correction c₁(m)=−42m proved; HYP-2420 MOS-mechanism CORRECTED to secular, n=3696 reproduced). THM-487 ([72,36,16] obstruction is code-combinatorial — W₇₂ all-positive, Γ₇₂ exists, Paley gauge stalls at d=12). Renumbered from 485/486 per claudebox-S5 first-come.
**Next:**
1. Prove the rigidity (HYP-2417): the uniform min-modulus bound on the sparse pentagonal polynomial P_S over a boundary circle ⟹ Rouché ⟹ interior zero for all finite flip sets.
2. The η^{−b} Lyapunov family: γ for general b interpolating partitions (b=1) and codes (b=24).
3. The ternary Gleason–Pierce analog — extremal ternary is already negative at n=72 (Pierce); is there a ternary γ-rigidity?
4. Sparse-lag extension of Goldsheid–Zeitouni (arXiv:2505.00377) ⟹ rigorous existence/positivity of γ_pent.
5. HYP-2422: probe the random-extremal never-negative prefix tail at larger m (does deterministic stay the strict maximum-delay outlier?).

## Lead kind-pasteur-2026-06-11-S2: Gleason × tournaments closed (THM-481 joint w/ claudebox-S3), zigzag law opened (THM-483) — Erdős–Moser #1216 corrected

**Done:** THM-481 (merged per MSG-870): both Gleason Type II generators tournament-generated
(g₂₄ = C(I+S(Paley₂₃)), ê₈ = Paley₇); eQR identity proved then attributed — Kim–Solé, Des.
Codes Cryptogr. 49 (2008) Prop. 3 (found by the adversarial round's literature mandate);
Gleason-generation framing novel. THM-483: zigzag law trans(D(T)) = z(T); HYP-2360 +2
sandwich REFUTED by the alternating family A_l (δ = l unbounded, counterexample at n = 7,
one past the n=6 census); THM-455 addendum; Erdős–Moser bounded-increment route closed.
**Next:**
1. HYP-2413/HYP-2442: tower support gates - `trans(T127)=15` and `trans(T255)=23` are now exact, so the frontier is `T511`. Do not brute-force first; mine marked packets `X subset T127` with `trans(X)=15` and `q(X)=trans(D(D(X)))-30`. A packet with `q>=1` proves the predicted `trans(T511)>=31`; chain-only packets have `q=0`.
2. HYP-2409 claims 2–3 (dictionary glue map; check claudebox THM-482's even-row hypothesis covers all tower levels).
3. Is RM(2,5) a tournament gauge of ANY order-32 skew-Hadamard? (finite sharp question, THM-481 remark).
4. GF(p) Kim–Solé: p-ary Paley tournament codes vs the Gleason–Pierce self-dual rings.

## Lead mac-mini-2026-06-11-S1: REED-MULLER ON THE TILING CUBE — the blue code + the digit ladder (THM-477/478, HYP-2406..2408, T779)

**Status:** session complete; both theorems proved + adversarially verified (fresh conventions, zero corrections).
**What landed:** blue tilings = dual-containing linear code ker(1+sigma), self-dual defect = hypotenuse f, Plotkin=Mode-B recursion b(n)=b(n-2)+(n-2); digit ladder of H on the tiling cube: digit_0 = RM(0,m) (Redei), deg(digit_1) = D = 2*floor((n-1)/2) (bound PROVED via cycle-reversal pairing, equality n=4..7), digits >= 3 saturate; (D+1)-flat XOR annihilation = the literal RM(r)^perp = RM(m-r-1) statement for H mod 4. Digit-degree question for counting functions confirmed UNPUBLISHED (model case: Hamming weight digits = e_{2^i}, Canteaut-Videau).
**Next steps:**
1. **Prove HYP-2406 equality** (deg digit_1 = D for all n): exhibit a single ANF monomial of degree D with odd coefficient — candidate: the top-cycle monomial through the base path; likely a short argument on top of THM-478's pairing.
2. **HYP-2408 Ax shadow:** apply Ax divisibility (zeros of deg-D functions divisible by 2^(ceil(m/D)-1)) to joint digit level sets; target a coding-theoretic re-proof or strengthening of H ∉ {7,21}; compute achievable-residue tables of H mod 8, 16 per n.
3. **d_2 sequence 3,6,7,11** (n=4..7): extend to n=8 (2^21 tilings x H — feasible in C), find the law (HYP-2407).
4. **The extended-Hamming slot:** dual of the skew-tower code (RM(1,m)=Sylvester analog exits RM at order 16, THM-451) — what tournament object is RM(m-2,m)? Unexamined.
5. **Blue slice:** restrict digit functions to B_n; compute degrees in the sigma-invariant algebra; connect to SC/transpose-self classes (canon: transpose-self never pure black).
6. **Cite/absorb Sole-Zaslavsky 1994** (cosets of cocycle code = switching classes of signed graphs) into the THM-474 orbit — the tournament/base-path-gauge version remains the project's.
**Sources:** 04-computation/rm_digit_ladder_macmini_0611s1.py + thm477/478 adversarial rechecks; 05-knowledge/results/*0611*.out; 07-reflections/reed-muller-on-the-tiling-cube.md.

## Lead mac-mini-2026-06-10-S2: THE DETERMINANT LENS — d(T) = det(I+S)/2^(n-1) on the metagraph; tilings = switching classes (THM-468/472/473/474, HYP-2383..2389+2398..2400, T777, MISTAKE-071)

**Status:** session complete; theorems proved + adversarially verified; census exact n≤9; conjecture witnesses at n=13.
**What landed:** floor (d=1 ⟺ local order, attribution Knuth/Babai-Cameron/Boussairi-D₁), ceiling (DRT switching classes ⟺ skew-Hadamard, verified zero-corrections), average (involutions/Hermite, attribution KMPRS LAA 707 (2025)), gauge theorem (tilings = switching classes), tournament Barba conjecture (OPEN-Q-058), mod-4 score law (HYP-2398), carousel takeover (HYP-2388/2399), d ⊥ H + d-not-second-eigenvector (new metagraph coordinate), HYP-2312 universal form refuted (MISTAKE-071).
**Next steps (ranked):**
1. **Prove the tournament Barba bound** (OPEN-Q-058) — the genuinely new extremal problem; publishable with Klanderman et al. as the companion.
2. **Build the switching metagraph S_n** (A049313: 1,1,2,2,6,12,79): nodes = switching classes up to iso, per-node d + H-spectrum + SC composition; Babai–Cameron Thm 3.2 gives per-class tournament counts. New quotient layer above G_n/Z_2.
3. **Prove the mod-4 score law** (HYP-2398) via Pfaffian-minor parity expansion — clean THM candidate.
4. **OEIS submissions:** max det(I+S) 2,4,16,32,160,512,4096,8192 and d-version 1,1,2,2,5,8,32,32 (absent, checked); extend A334123 (labeled max-H counts: a(7)=240 confirmed; a(8),a(9) computable from census + |Aut| of the 6/1 max classes); extend A038375 to n=12+ via Pf=±1-existential pruning (HYP-2312 amended) or circulant-first search.
5. **Klanderman arc-flip update on the metagraph:** det S' = det S(1+2S⁻¹_ij)² assigns every wiggly line an exact determinant weight; flip-silent ⟺ S⁻¹_ij ∈ {0,−1}. Compute the det-silent wiggly fraction per class — the d-analog of silent mutations.
6. **Grinberg–Stanley mod-4 Rédei (H ≡ 1 + 2·#odd-cycles mod 4, published 2023) = depth-1 OCF truncation** — formalize against THM-062/OCF; their U_D machinery vs kind-pasteur's THM-466 2-adic digit tower (flagged to kp in MSG).
7. **GLMY path homology per iso class** (literature: uncomputed for general tournaments; Ω₂ spanned by transitive triangles + square differences) — new metagraph invariant; cross with THM-120's |Pf| phase separation.
8. **Schweser–Stiebitz–Toft arXiv:2510.10659 (Rédei revisited, Feb 2026)** — read for the Dirac/Berge strengthenings; Satake 2025 (nearly-DRTs ↔ almost difference sets, conditional on H-L Conjecture F) for the q ≡ 1 mod 4 side.
**Sources:** 05-knowledge/results/{hadamard_det_census,hadamard_det_n9,circulant_odd_det,tournament_simplex,metagraph_det_gradient}_macmini_s2.out + thm472/thm473_adversarial_check.out + verify_hyp2312_and_mod4law_subagent_2026-06-11.out; reflection 07-reflections/the-determinant-lens-sgn-vs-chi-and-the-three-geometries.md.

## Lead mac-mini-2026-06-11-S2: BSD + HODGE merged into the self-dual-middle frame (T786, HYP-2420..2424, OPEN-Q-095, reflection the-self-dual-middle-where-parity-lives)

**Status:** session complete; reflection + 5 HYPs; rigorous spine adversarially verified; one duplicate (THM-490) caught & withdrawn.
**What landed:** (1) RIGOROUS SPINE — alternating pairing ⟹ square (Pfaffian) is the shared symplectic mechanism behind |Ш(E)|=□ (Cassels), det(I+2A)=Pf² (THM-174), Weil/Hodge alternating forms; the square-defect is 2-ADIC in both worlds (Poonen-Stoll c∈Ш[2], the 8=2³ in THM-442's 8Q) — the precise home of the user's "doubling." (2) FRAME — the self-dual-middle table (RH s↔1−s / BSD root number / Hodge p↔q / Goldbach swap / project complement-σ-SC), parity lives at each fixed middle; the one rigorous ½/doubling link is λ(2n)=−λ(n). (3) MERGE — additive Goldbach (primes, segment) vs multiplicative Euler (strong-tournament H-primes); the s=2 polygonal segment bridges them; the project's PROVEN genus-2 H-semigroup R=odds\{7,21} (HYP-2271) is the tournament ζ-side.
**Next steps:**
1. **OPEN-Q-095** (genuinely open): is there a finite abelian group / pairing attached to a tournament (Ш-analog) whose order = det(I+2A) and whose alternating-vs-antisymmetric type (hence square-vs-2×square) is controlled by a tournament parity (SC / transpose-self / blue-black)? Compute the Smith normal form of I+2A and S; see whether SNF squares track det=Pf² and whether the 2-part carries a c∈[2]-style defect. Cf. Klanderman et al. (SNF of skew D-optimal designs).
2. **HYP-2421** (polygon side = additive arity): make the s↔arity alignment a precise statement on the project's polygonal↔Goldbach ladder (HYP-1962); does s=2 segment / s=3 triangle / s=4 square map to a tournament/figurate hierarchy with a clean transition law?
3. **The root-number analog:** is there a tournament "root number" (a global ±1 = product of local ±1's over spine/ribs/sea or over strong components) whose sign predicts a parity (e.g. of some rank-like Q-grading)? The honest version must avoid the decorative local-global product; look for a genuine multiplicative ±1 over strong components paralleling w(E)=∏wᵥ.
4. **Cross-link** HYP-2424 to the existing band-gap/numerical-semigroup reflection (polarized-delta-fields-...-s699) — the genus-2 semigroup already has a physics (band gap) and NT (Frobenius) reading; add the Euler-product/Goldbach-duality reading.
**Sources:** 04-computation/h_realizability_goldbach_macmini_0611s2.py + .out; verify_moon_realizability_fresh_opus_0614.py/.out; reflection the-self-dual-middle-where-parity-lives.md.

## Lead kind-pasteur-2026-06-11-S1: Erdős 592 session 4 — the seam EXPLAINED (sum-free gradings), the t=7 wall, first Schipperus cutoffs, m=3 probe (THM-469/470/471, HYP-2390..2397)

**Done:** Picked up mac-mini-2026-06-10-S1's Next items 2/3/4. (3) SOLVED → THM-469 (PROVED):
v_p level sets sum-free ⟺ p=2 (parity closure); odd-p algebras die by 3-term-AP mono-collapse
(14-clause minimized UNSAT core at (3,4)); b=3 branching control ⟹ the seam follows the SCHUR
ARITY, not branching (HYP-2390 CONFIRMED; THM-464 D slogan corrected); leading-digit refinement
rescues every odd p (HYP-2391 CONFIRMED: (sign,v₃,lead) SAT (3,4)/(3,5) verified). (2) ANSWERED
NEGATIVELY SO FAR → THM-470: coarsening collapse (rungs must REFINE F2); dyadic 1-jet F2J dies
PER-SIZE at (3,7) ⟹ F2 itself dies at t=7 (refinement) ⟹ gap-determined algebras antitone in t:
decided alive t≤6/dead t≥7 — master run on the FULL invariant algebra at (3,7) [in flight at
write-up]; HYP-2396 candidate: R(n,2) = 2n+1 (3, 5, 7?; lower bounds Q(n,2n) SAT match at
n=1,2,3). (4) BUILT → THM-471: B3 general-shape tuple grammar implemented + validated (m=1
reproduces THM-460 D); j1-march size explosion documented (BT(3,2) j1 = 3.5M leaves, vacuous
guard); first Schipperus-forced cutoffs ever computed: m=2 M=1 UNSAT at (2,2) and (2,3); m=3
probes at (2,2) in flight.
**Next:**
1. Harvest the in-flight runs (invariant wall at (3,7); F2X t=6,7; m=2 M=2-j0 sweep; m=3 M=1/M=2-j0 at (2,2); crossval).
2. If invariant dies at (3,7): prove/refute R(n,2) = 2n+1 (HYP-2396) — try the "+2 per level" induction via THM-459's L1-style row reduction; if it lives: extract the invariant witness's non-dyadic structure = the true next rung.
3. NON-invariant rungs (value-dependent features, e.g. v₂ of coordinate values themselves) — the only route left to a t-uniform table if the wall is real.
4. m=3 vs m=2 differential at matched (s,c) with M=2-j0 and M=3 shapes; first Chang number (m=1 (3,3)) still open.
5. THM-459 hand closure; f_grid(4) optimality (inherited).

## Lead mac-mini-2026-06-10-S1: Erdős 592 session 3 — Q(3,5)/Q(3,6) SETTLED, the algebra ladder is strict, the seam is 2-adic at n=3; cubic lens placed (THM-464/465, T770, HYP-2373..2376)

**Done:** Q(3,5) = SAT and Q(3,6) = SAT, settled by explicit bi-dyadic (sign+v₂ per coordinate) witnesses, independently verified — the 80k-clause timeout instance falls in 2.8s in the feature quotient (POKE Task 2 progressed, Task 1.3 advanced: R(3,2) > 6). Uniform-table SAT over t=4..7 is feature-UNSAT (0.3s) ⟹ NO infinite (sign,v₂)-measurable strong witness: the ladder signs ⊊ signs+v₂ ⊊ (open) is strict. Triadic control: sign+v₃ (equal-size algebra) is feature-UNSAT at (3,4) instantly ⟹ the 2-adic seam is real at n=3; at n=2 ALL gradings share cutoff 5 (v₃ control). R_b(1) = R(3,b) proved+verified (classical Ramsey = height-1 row). Cubic lens (HYP-2376 census): cubes sum-free forever / non-Sidon from 1729 = first C4 of the cubic summand graph; signed from (3,4,5,6); cites kind-pasteur THM-462/463 + HYP-2370.
**Live at close:** batched (3,3) Chang (first Chang number) and ternary free sweep t≥7 — outputs stream to 05-knowledge/results/erdos592_{chang33_batched,ternary_seam}_macmini_s3.out.
**Next:**
1. Harvest the two live runs; R_3(2) exact; first Chang number.
2. CLIMB THE ALGEBRA LADDER: find the next rung past sign+v₂ that admits a t-uniform table (candidates: unbounded-v₂ with tail conditions, Larson partial-sum/scheme features, mixed row-grading × column-algebra) — a uniform rung + König = constructive strong Specker at ω³.
3. Explain the v₂/v₃ asymmetry algebraically (parity closure odd+odd→even vs the mod-3 escape; THM-464 D open note).
4. m≥2 tower enumerator (THM-460 B3 grammar) — still the open-case probe (POKE Task 1.2).
5. THM-459 hand closure; f_grid(4) optimality.

## Lead mac-mini-2026-06-09-S2: Erdős 592 session 2 — Chang towers, dyadic witnesses, R(2,2)=5 lemma layer (THM-459/460, HYP-2363..2366)

**Source:** second 592 session; continues the S1 lead below.
**Done this session:** THM-459 (R(2,2)=5: lemma layer L1–L5 — doubly-dark clique, C₅/triple dichotomy, trace cliques, composition-freeness — + machine closure; hand closure open); THM-460 (tower miniature for ω^(ω^m): full-type ⟺ stacked towers PROVED, recursive binary shape grammar, König bridge ⟹ Chang/Schipperus FORCE finite tower-Ramsey cutoffs; m=1 sweep partial); invariant cutoff = free cutoff = 5 at n=2; DYADIC sufficiency (B_g through v₂(g)) verified through (3,4) — the 2-adic seam shows up inside 592 witnesses; f_grid = 1, 7, ≤30; witness-tournament H flier (odd, transitive-leaning).
**HARVESTED post-close (both runs completed; results in 05-knowledge/results/):** Q(2,5) UNSAT independently RE-CONFIRMED (R(2,2)=5 double-certified by two complete-verifier implementations); Q(3,5) TIMEOUT at 3605s/80,111 clauses (free) and 1801s/62,854 (invariant) — UNDECIDED; Chang M=2: **(2,5) SAT with an 83-edge hard-fought witness (916s)** — the tower game survives at the exact ambient where the grid game dies, so Chang numbers sit strictly deeper than Specker numbers; (3,3) and all larger instances TIMEOUT — the first Chang number remains uncomputed (C2 guarantees it exists); M=3 vacuity guard: s ≥ 3 required (height-3 grids need 3 split positions); Q(4,3) SAT (115 edges).
**Next steps:**
1. Decide Q(3,5) and the (3,3)/(2,6) Chang instances: needs symmetry breaking (row/value permutations), incremental cardinality, possibly parallel CEGAR — the instances are at the edge of naive encodings.
2. Implement the THM-460 B3 general-shape enumerator (towers of towers) ⟹ run m=2 (Schipperus-forced cutoffs) and the **m=3 probe of the open $1000 case** — the designated centerpiece.
3. Extract a closed-form dyadic rule from the invQ(3,t)/dyadic witnesses (B_g via v₂(g) + within-row R) and try to prove it works for all t ⟹ via THM-453 D1, a constructive strong form of Specker's ω³ negative.
4. Hand-close THM-459 (the two regimes are bounded case analyses; SAT core is small).
5. Witness-ensemble statistics for the H-flier before claiming anything (HYP-2364 note).

## Lead mac-mini-2026-06-09-S1: Erdős 592 — the R(n,t) tree-grid ladder and the rank-3 graded-relation frontier (THM-453, T768, HYP-2344..2346)

**Source:** Erdős 592 session (ω^β → (ω^β,3)²; open case α=ω^(ω³) = exponent with three CNF summands).
**Status:** ACTIVE — exact finite results at n≤2, n=3 partially computed, ω^(ω^γ) lift not started.
**What is established (THM-453):** witness frame (triangle-free graphs, no full-type independent set); grid characterization of full-type subsets; the compactness bridge Q(n,t)-SAT-∀t ⟹ ω^n ↛ (ω^n,3)² (strong witness), hence positive relations force finite cutoffs R(n,2); computed **R(1,2)=3, R(2,2)=5 (exact)**, Q(3,3) SAT (35 edges), Q(3,4) SAT (346 edges, complete verifier, 692s); no pure order-pattern witness exists at ω³ (all 13 maximal triangle-free pattern sets grid-avoidable); the unique triangle-free pattern class at ω² is the SHIFT GRAPH.
**Next steps, in order of leverage:**
1. Decide Q(3,5) (was running at session close; CEGAR + complete verifier; needs either patience or a smarter verifier — bitset adjacency, subgrid search ordered by solver phase).
2. Hand-prove R(2,2)=5 (UNSAT certificate is small: 532 CEGAR clauses; extract a human Zarankiewicz/Mantel argument — would make THM-453 D+E a self-contained finite theorem with a clean constant under Specker's theorem).
3. Prove Q(3,t) SAT for all t by a UNIFORM rule (extract from the t=3,4 witnesses; the invariant quotient (R,{B_g}) format is the right ansatz) ⟹ via THM-453 D1 a NEW constructive proof of Specker's ω³ ↛ (ω³,3)², SAT-discovered.
4. Lift F1 to ω^(ω^γ): write the graded-relation algebra over the gap semigroup (CNF of γ; summand count = grading rank); recover Schipperus's ≤2-summand impossibility / ≥4-summand constructions as the two regimes; instantiate rank 3 = the open $1000 case.
5. Import literature exacts: Darby (JCTB, "Negative partition relations for ordinals ω^(ω^α)"), Larson "pentagram" (JSL 2000), Schipperus APAL 2010 constructions; EHMR's presentation of Specker's witness — compare with the SAT-discovered structure.
**Files:** 04-computation/erdos592_*.py (6 scripts), 05-knowledge/results/erdos592_*.out, 03-artifacts/drafts/erdos-592-survey-and-reframes-macmini-s1.md, 07-reflections/the-additive-ladder-reaches-the-ordinals-macmini-s1.md. Mistakes: MISTAKE-066 (bridge direction), MISTAKE-067 (incomplete CEGAR verifier — read the witness).

## Lead kind-pasteur-2026-06-09-S2: Erdős problems through the doubling lens (THM-455/456/457, T769, HYP-2356..2361, MISTAKE-068/069)

**Status:** ACTIVE — all branches adversarially verified.
**Erdős 64 (#64, Erdős–Gyárfás):** THM-456 (blowup spectrum law spec(G[K₂]) = gap-free [3, 2s(G)]; a single edge plants a twin C₄ ⟹ blowups never counterexamples; Turán corridor closed n≤9; all 71 C₄-free δ≥3 graphs n=10–12 killed by forced C₈) + THM-457 (dyadic gate ladder: girth ladder 24 / 28 (NEW EXACT) / ≤32 / >46 / 58; ladder principle: closing each dyadic gate inflates the next; dihedral 3-reflection Cayley family with dyadic spectrum exactly {4,32}; Exoo G78 reconstructed — contains C32, new beyond the 2014 paper). MISTAKE-069: McGee HAS 34 C8s (S710 enumeration-order artifact corrected).
**Erdős–Moser (#1216):** THM-455 — tower trans = 3,5,7,11 (extremality window 7–31, closes by 63); sandwich trans(D(T)) − trans(T) ∈ {1,2} over all 32768 n=6 tournaments (HYP-2360); Reid–Parker FORCES the Paley₇ +2 exception (verified trans(D(Paley₇)) = 5); T31 ties published trans(QR₃₁) = 7 (solver externally validated vs Momihara–Suda Table 2).
**Next steps:** half-life ladder prediction at order 128 (HYP-2361); 2-adic character condition for the {4,32} dihedral family; classify alternating-chain +2 exceptions; hand trans data to mac-mini's 592 thread (THM-453 Part A bridge).

## Lead kind-pasteur-2026-06-09-S1: the skew-Sylvester doubling D(T) — Walsh/skew-Hadamard recursion on tournaments (THM-447/448/450..452/454, T767, HYP-2332..2337+2350..2355)

**Source:** user directive (n*2 doubling, three copies + one negated, skew-Hadamard normalization = tiling frame).
**Status:** ACTIVE — THM-447 proved + verified n=3..5; D(border(C3)) core ≅ Paley T_7 (iso verified); branch fan-out in progress (H-formula, Mersenne tower, tiling σ-eigenspaces, Hadamard equivalence, Clifford family classification).
**Key objects:** D(T) = [[M, M+I],[M−I, −M]]; spectral law M'² = I₂⊗(2M²−I) (Chebyshev T₂, λ²+1 doubles = Hadamard order); canonical Ham path (p, twin, reversed p); twin arcs saturate the anti-diagonal (hypotenuse) of the doubled staircase; tiling splits into σ-symmetric copy blocks + σ-antisymmetric cross block.
**Next steps:** (1) H(D(T)) / H(T[K₂]) functional in I(Ω(T),x) — would answer OPEN-Q-045 Q1; (2) identify tower DRT_15 and DRT_31 vs Paley; (3) Ω(D(T)) structure theorem; (4) engineering: skew-Walsh butterfly transform O(m log m).

---

## Lead monad-explorer-2026-06-07-S18: the two humps are one resurgent series — Watson bridge with explicit coefficients `bⱼ` (THM-438 ADD-16)

**Status:** PARTIALLY DONE (this session). The free factorial law's spectral tail `ρ(x)~e^{1−x}Σbⱼx⁻ʲ`
and moment ratio `A088368(k)/k!~e·Σbⱼ/(k)ⱼ` are Watson-lemma images sharing the EXACT rationals
`b=1,2,10,178/3,1178/3,42494/15,…` (`bⱼ=[tʲ]exp(R(σ(t))−1)`, `t=σ/Q(σ)`). VERIFIED: tail vs
parametric density `1e−16`; moments `k≤60`. Refines OEIS Kotesovec `a(n)~e·n!`. The `bⱼ` are
Gevrey-1 (divergent) ⟹ the `e`-overshoot on both sides is optimal truncation of one resurgent series.
**Next steps:** (1) **prove the Watson identity coefficient-by-coefficient** (leading order rigorous via
Carleman-determinacy + `e^{−x}` tail; full series = standard "moment asymptotics from tail expansion"
uniform-remainder theorem, resurgent ⟹ asymptotic only). (2) **closed form / GF for `bⱼ`** — the `bⱼ`
are built from the same Gompertz/`E₁` `g`; is there an `E₁`-type closed form, or an OEIS-negative
recurrence? (3) feed `bⱼ` back to the EDGE: the analogous `lnln` edge expansion (ADD-15 NEXT (3)) and
whether an analogous "edge bridge" links the `x→0` density to small-`k`/negative-`k` data.
**Files:** `04-computation/paley_starstar_tail_moment_watson_monad.py` (+`.out`); reflection
`07-reflections/the-two-humps-are-one-resurgent-series.md`; THM-438 ADDENDUM-16; HYP-2308/INDEX.

---

## Lead monad-explorer-2026-06-07-S17: the free factorial law's density is a closed-form parametric curve — push the crossing-`q` family `μ_q` (THM-438 ADD-15)
**Source:** THM-438 ADDENDUM-15, reflection `the-two-singularities-of-the-exponential-integral-shape-the-density.md`, HYP-2308.
**Status:** OPEN (new). ADD-15 gave the whole density as a parametric curve `x(u)=−u²g(u)`, `ρ(u)=−Im(u)/(π|u|²)` on
`Im(u²g(u))=0` (`u=−1/G`, `g=eᵘE₁(u)`), verified vs root-found ρ to `1e−12` and moments A088368 to `<0.3%`. The two
ends of the support = the two singularities of `E₁` (log at `0` → edge √log; cut `(−∞,0]` → tail). The tail constant `e`
is DERIVED as `e^{R(0)}=e^{κ₁}` via the Stokes term, with `ρe^x=e^{R(G(x))}` overshooting `e` (resurgent hump = the
`m_k/k!→e` hump of MISTAKE-063). DONE from the S16 lead: the edge structure is now fully transparent from the
parametrization (`arg u→−π/2`, `x~ε²(ln(1/ε)−γ)`); full `lnln` resummation still a clean finite task.
**Next steps:** (1) **`K_q(z)` / the parametric density of the crossing-`q` family `μ_q`** (the priority now). The
parametric machine `{x=K(w), ρ=−Im w/π on Im K=0}` works for ANY `μ` with explicit `K`. If `K_q` has an `E₁`-type
closed form, the two singularities should slide continuously across `0≤q≤1`: the `E₁` log (free edge √log) weakening
into the classical bounded edge `→e⁻¹` + atom `e⁻¹δ₀`. WHERE (in `q`) does the log→atom transition happen — `q=0⁺` or a
critical `q_c`? (2) **Belinschi–Nica `B_t(μ_free)`** — `B_t` in closed `K` form; any named image. (3) Resum the full
edge `lnln` expansion (now a finite calc from the parametrization). (4) Off-diagonal `t(k,m)` THIRD deformation; `t(7,5)`;
HYP-2308 remainder (non-circulant DRT n=15).

## Lead monad-explorer-2026-06-07-S16: the free factorial law's edge is `√log`, not a constant — resum it, and find where the atom dissolves in `μ_q` (THM-438 ADD-14)
**Source:** THM-438 ADDENDUM-14, reflection `eulers-divergent-series-is-the-free-factorial-laws-r-transform.md`, HYP-2308.
**Status:** OPEN (new). ADD-14 gave the free factorial law's K-transform a CLOSED FORM `K(z)=−(1/z²)e^{−1/z}E₁(−1/z)=−(1/z²)g(−1/z)`
(`g`=Gompertz fn; Euler's divergent series resummed; verified `<7e−17`), proved `K(−1)=−δ` (Gompertz constant), and
showed the `x→0` edge has NO finite constant: `π ρ√x → √(ln|G|−γ) ~ √(½ln(1/x))` (corrects ADD-12 `1/π` / ADD-13 `0.4–0.6`).
**Next steps:** (1) **Resum the full edge asymptotic.** `π ρ√x=√(ln|G|−γ)` is exact; solve `r²x=ln r−γ` with `2θ→π` to higher
order for the complete `x→0` expansion of the density (the `lnln`/`γ` corrections). Closed-form `K` makes this tractable.
(2) **`K_q(z)` for the crossing-`q` family.** Does `K_q` have its own exp-integral closed form? The classical end (`q=1`)
has a BOUNDED density (`→e⁻¹`) plus an atom `e⁻¹δ₀`; the free end (`q=0`) has no atom and a `√log·x^{−1/2}` divergence.
Where in `0<q<1` does the atom dissolve / the divergence switch on — at `q=0⁺` or a critical `q_c`? (Ties to S15.)
(3) **Belinschi–Nica `B_t(μ_free)`** is now computable in closed form (`K` explicit) — any named `B_t` image gives the
wild end a second analytic handle.

## Lead monad-explorer-2026-06-07-S15: the Belinschi–Nica `B_t` semigroup vs the crossing-`q` family `μ_q` (THM-438 ADD-13)
**Source:** THM-438 ADDENDUM-13, reflection `the-two-endpoints-are-bercovici-pata-partners.md`, HYP-2308.
**Status:** OPEN (new). ADD-13 established the two named endpoints A000262/A088368 are Bercovici–Pata partners
(shared cumulants `n!`, classical↔free CP of `ν=e^{-x}dx`), joined by a positive-definite measure family `μ_q`
(`m_k(q)=Σ_π q^{cr(π)}∏|B|!`, Hankel `D_n(q)` all-nonneg `q`-coeffs, `q=0` free, `q=1` classical).
**Next steps:** (1) The Belinschi–Nica `B_t` transform interpolates `∗`/`⊞` analytically with `B_1=Λ` — is the
crossing-`q` family `μ_q` the SAME interpolation? If so it pins `μ_q`'s Cauchy transform to a known closed form
and may yield the free density's exact edge constant (numerics give `≈0.4–0.6`, NOT `1/π`; see ADD-13 part 6).
(2) The off-diagonal `t(k,m)` columns are now ruled out as BOTH the crossing-`q` triangle AND the rate-marked
`N(k,j)` (ADD-12) — they need a THIRD deformation, transverse to both named axes. Find it.
(3) `μ_q`'s own cumulants (q-deformed) — does `μ_q` stay compound-Poisson-like for `0<q<1`?

## Lead monad-explorer-2026-06-07-S12: PROVE the diagonal `t(k,k)=A088368(k)=Σ_{NC(k)}∏|B|!` (the cleanest THM-438 sub-target)
**Source:** THM-438 ADDENDUM-10, MISTAKE-063, reflection `the-cancellation-runs-between-two-named-endpoints.md`.
**Status:** VERIFIED k≤7 + OEIS match (A088368 = Callan "partitions of [n] into sets of noncrossing lists", arXiv:0711.4841). Closed form `t(k,k)=Σ_{π∈NC(k)}∏_B|B|!` VERIFIED k≤7. Asymptotic `~e·k!` (Kotesovec) — ADD-9's "`≁e·m!`" RETRACTED (the ratio overshoots e, peaks at m=8, descends back; MISTAKE-063).
**Next steps:** (a) **Prove the bijection** diagonal even-series patterns (= max-cycle-rank = doubled plane trees, Euler tour visits v `deg(v)` times, weight `∏(deg−1)!`) ↔ NC(k) with block-factorial weight `∏|B|!`. Finite, number-theory-free; would upgrade `t(k,k)=A088368(k)` to PROVED and give `h_m(m)→e` as a corollary. (b) The two open handoffs (#1 `(k)_m|t(k,m)`; #2 `g_m(−1)=(−1)^m(2^m−1)`) live at the TAME end — the wild diagonal is now closed-form, so attack the tame end via the `(1+s)|Q_m` line-parity involution. (c) Off-diagonal columns/residues are all OEIS-NEGATIVE (P_m(1)=1,3,20,181; subdiag 9,72,580,4845; col3 13,72,230,560; unsigned rowsum 1,4,23,160,1262,10944; `Σ_{NC}∏(|B|−1)!=1,2,6,23,105,553,3311`) — do NOT re-hunt them in OEIS. (d) Still need `t(7,5)` (the k=7 background run DIED at k=6): a core-aware enumerator validated vs k≤6.

## Lead monad-explorer-2026-06-07-S5: prove `(★★)` via walk-induced ribbon genus + A215257 bijection
**Source:** THM-438 ADDENDUM-3, HYP-2308, reflection `the-drt-engine-is-S-squared-equals-J-minus-nI-the-catalan-is-genus-zero.md`.
**Status:** OPEN. `(★★) Σ_{even-series σ}μ(0̂,σ)=(−1)^kC_k` VERIFIED k≤5; value = free cumulant of the two-point law. Leading-order engine now `S²=J−nI` (DRT-universal, number-theory-free, VERIFIED p≤43).
**Next steps:** (a) build the rotation system the Euler walk `x_0..x_{2k}` induces on `G_σ`, compute ribbon genus per pattern (k≤5 enumerated on disk), test `Σ_{genus 0}μ=(−1)^kC_k` & `Σ_{genus>0}μ=0` (naive non-crossing-PARTITION version is REFUTED — do not retry). (b) write the explicit bijection even-series patterns ↔ indecomposable deque-sortable permutations (OEIS A215257); its sign structure may prove `(★★)`. (c) HYP-2308 remainder: expander-mixing `o(n^{k+1})` bound, tested on a verified non-circulant DRT n=15. (d) extend cycle-rank triangle to k=6,7 (needs enumerator smarter than all-set-partitions) to pin its recursion.

## Priority opus-2026-05-16-S1: TRRT and All-0 Staircase

### INV-191: H=63 Unlocks at n=8 via Complete Conflict Graph
**Source:** opus-2026-05-29-S8; exact census opus-2026-05-29-S10
**Status:** EXACT at n=8; single-core mechanism identified in S11
**What:** HYP-1754 ("H=63 is universally forbidden") is refuted. A concrete n=8 tournament has H(T)=63 by both DP and direct permutation enumeration. Its odd-cycle conflict graph Ω(T) has 31 directed odd cycles and is complete, so OCF gives H=I(K31,2)=1+2·31=63. S10 upgraded this to a finite theorem: among all 6880 n=8 isomorphism classes, exactly two have H=63; both have |Aut|=1, score sequences (1,2,2,3,3,5,6,6) and (1,1,2,4,4,5,5,6), and Ω(T)=K31. This explains how 63 bypasses the old disconnected K3-factor obstruction: it realizes 63 through a complete Ω, not through K3⊔2K1.
**S11 update (opus-2026-05-29-S11):** Both H=63 classes are **single-core**: every odd directed cycle contains one vertex, and deleting that core vertex leaves a transitive 7-vertex tournament. Their transitive-insertion signatures are `1001100` and `1100110`; the weighted signature count
`r(s)=Σ_{i<j, s_i=1, s_j=0} 2^{max(j-i-2,0)}` equals 31, matching the number of odd cycles. A complete-Ω census over isomorphism classes n=3..8 shows cycle counts r=3 and r=10 are absent whenever Ω is complete; a signature target search finds r=3 and r=10 absent for all single-core signatures up to m=16, while r=31 first appears at m=7. This gives a new focused H=21 lens: prove the r=10 gap in the single-core family, then handle non-core complete Ω and non-complete α-tuples separately.
**S12 update (opus-2026-05-29-S12):** Reframed the H=63 mechanism as a projection defect. The two H=63 classes are exact old-projection kills: deleting the core vertex loses 31/31 directed odd cycles and leaves `H(T-v)=1`, `alpha(T-v)=[1,0]`. A core-stratified complete-Ω census through n=8 confirms r=3 and r=10 are absent in every core stratum, and r=31 occurs only in core-size-1 classes. The single-core target search was extended to m≤40: r=3 and r=10 remain absent; r=31,42,63 appear with simple linear count laws after their first occurrence.
**S5 update (opus-2026-05-30):** The applied residue/phase contrast table classifies H=63 as the clean exact-kill endpoint of the standard anomaly set: `rho=1`, max-loss residue `(0,0,0)`, residue rank 0. This sharpens the next step: treat the single-core signature count as a finite-state target problem for absent values `r=3,10`, separate from THM-025-style near-kill failures.
**Next:**
  1. Prove the single-core signature formula and the persistent gaps r=3,10.
  2. Classify complete-Ω tournaments with empty cycle-family core; compare their r-support to the single-core support.
  3. Revisit H=21 by separating the complete-Ω case Ω=K10 from the non-complete α-vector cases.
  4. Locate the H=63 classes in the merged metagraph/principal-line coordinates.
  5. Prove or refute the projection-kill/near-kill hypothesis by scanning more real-root failures and complete-Ω classes.
**Files:** `04-computation/h63_counterexample_audit_s8.py`, `04-computation/h63_n8_isoclass_census_s10.py`, `04-computation/omega_extreme_fingerprints_s11.py`, `04-computation/projection_defect_bridge_s12.py`, `05-knowledge/results/h63_counterexample_audit_s8.out`, `05-knowledge/results/h63_n8_isoclass_census_s10.out`, `05-knowledge/results/omega_extreme_fingerprints_s11.out`, `05-knowledge/results/single_core_signature_targets_s11.out`, `05-knowledge/results/projection_defect_bridge_s12.out`, THM-344, MISTAKE-050, HYP-1757, HYP-1758, HYP-1760, HYP-1761, HYP-1762.

### INV-189: Real-Rootedness of I(Ω(T), x) for All Tournaments (TRRT)
**Source:** opus-2026-05-16-S1
**Status:** STALE/REFUTED AS UNIVERSAL — see THM-025
**What:** The universal conjecture "I(Ω(T), x) has all real, negative roots for every tournament T" is false. Canon THM-025 gives an n=9 counterexample with score sequence [1,1,3,4,4,4,6,6,7] and I(Ω,x)=1+94x+10x²+x³, violating Newton's inequality. What remains interesting is to characterize the large real-rooted subclass and understand why failures are rare.
**Key implication:** Conditional ultra-log-concavity and product-formula ideas remain useful only for the real-rooted subclass, not universally.
**S5 update (opus-2026-05-30):** The contrast table identifies THM-025 as a dangerous near-kill rather than an exact kill: one vertex supports 92/94 odd cycles, but the max-loss residue is `(alpha1,alpha2,alpha3)=(2,1,0)`, residue rank 2, and `I(R,2)=9`. New working target: failures may require both a tiny rank-2 deletion residue and a root/phase imbalance.
**Next:**
  1. Reconcile the newer Hermite-Biehler/interlacing notes with THM-025.
  2. Characterize the THM-025 failure class and the generic real-rooted class.
  3. Search for a stability/negative-dependence theorem with explicit hypotheses.
  4. Run an n=9 census stratified by `omega_near_kill_rank2_vertices` and Newton/root margins.
**Reflection:** `real-rootedness-omega-conjecture.md`

### INV-190: All-0 Staircase H Sequence and Algebraic Structure
**Source:** opus-2026-05-16-S1
**Status:** EXTENDED to k=12 (monad-researcher-2026-06-02-S578); c3 formula formalized in Lean (monad-formalizer-2026-06-04-S1)
**What:** H values of all-0 interleaved staircase at n=2k. 5,29,233 are Markov numbers (breaks at k=5: 2489=19×131). c3=k(k-1) confirmed through k=12.
**Full sequence (k=2..12):** 5, 29, 233, 2489, 33773, 562685, 11222321, 262755369, 7110764837, **219612027389, 7658921303353**
**Proved:** # directed 3-cycles = k(k-1) exactly.  The staircase specialization was formalized sorry-free in `eliott-monad/math-lean` commit `b5ffcde` as `Math.Tournaments.allZeroStaircase_threeCycleCount`; THM-410 records the more general interval-reversal counting theorem.
**Growth ratios:** r(k)=H(k)/H(k-1): 5.80, 8.03, 10.68, 13.57, 16.66, 19.94, 23.41, 27.06, 30.88, 34.87. Differences ≈ +3.8/step (slowly increasing). Empirical: r(k) → 3k asymptotically; deficit 3k−r(k) peaked at k≈6 then decreases toward 0.
**Ruled out:** No order-2/3 linear recurrence. Markov fails all consecutive triples. OEIS search inconclusive (oeis.org blocked all fetches; web search for specific terms found no match — likely novel sequence).
**Next:** (1) Submit to OEIS. (2) Formalize THM-410's general interval-reversal theorem in Lean. (3) Find algebraic structure (product formula, generating function). (4) Try k=13 if feasible (~5 min runtime with array.array). (5) Investigate why r(k)→3k — is there a combinatorial explanation?

---

## Priority S24: New Leads from Session S24

### INV-184: Stanley-Stembridge PROVED (Hikita, arXiv:2410.12758)
**Source:** opus-2026-04-05-S24 web search
**Status:** NEW — needs integration
**What:** e-positivity of chromatic symmetric functions for unit interval graphs proved by Tatsuyuki Hikita (Oct 2024, revised Dec 2025). Since Mitrovic-Stojadinovic (arXiv:2506.08841) connects Redei-Berge to chromatic functions, and Stanley-Stembridge concerns e-positivity, this potentially constrains the h-positivity of U_T.
**Next:** Read Hikita's proof. Check implications for our U_T framework (THM-062/063). Does h-positivity of U_P follow?

### INV-185: Mitrovic NC Deletion-Contraction (arXiv:2504.20968)
**Source:** opus-2026-04-05-S24 web search
**Status:** NEW
**What:** The original Redei-Berge function does NOT satisfy deletion-contraction. Mitrovic introduces a noncommutative analogue that DOES. This enables inductive proofs — exactly what we need for Claim A-type results.
**Next:** Read the NC deletion-contraction construction. Can it give a new proof of Claim A? Does it yield new tournament invariants?

### INV-186: Real-Rootedness of I(Ω(T), x)
**Source:** opus-2026-04-05-S24 computation
**Status:** PROVED for n ≤ 8, DISPROVED universally at n=9 (THM-025)
**What:** All zeros of the independence polynomial of the odd-cycle conflict graph are real and negative for n≤8 by claw-freeness + Chudnovsky-Seymour. At n=9, THM-025 gives a counterexample. The right problem is now structural characterization, not universal proof.
**Next:** Characterize the n=9 failure mechanism and identify sufficient tournament conditions for real-rootedness.

### INV-187: E[H(T)] = n!/2^{n-1} 
**Source:** opus-2026-04-05-S24 computation
**Status:** PROVED (linearity of expectation)
**What:** Clean closed form for average Hamiltonian paths. W(n) = n! × 2^{C(n-1,2)}. Probably known but not explicitly in our bibliography.
**Next:** Check literature for this result. Is it in Moon's book?
**EXTENSION (monad-explorer-2026-06-07, HYP-2307):** the Paley/circulant maximizer's RATIO over this average, `R(p)=H(T_p)/(p!/2^{p−1})`, → **e** — PROVEN via a character-sum cluster expansion (`R=E_σ[∏(1+χ(d_k))]→exp(Σa_L)`, only the cherry cluster survives, `a_2=−χ(−1)=1`). Universal across circulant tournaments (cherry weight = tournament condition `g` odd). **SUB-LEMMA CLOSED → THM-438 (monad-explorer-2026-06-07, same day):** `a_{2k}=0 ∀k≥2` PROVEN UNIFORMLY (no per-k Weil): `B_L=0` ⟹ `A_L=−Σ`coincidence-patterns; no-leaf forces `V≤2k`; only the single `V=2k` pattern (one even cycle) needs Weil, all others `o(p^{2k})` trivially. **Plus the CATALAN LAW** `A_{2k}=C_k p^{k+1}+O(p^{k+1/2})` (C_k=Euler tours of plane trees=moment-method tree-walks). See `04-computation/paley_cluster_{sharp_order,catalan}_monad.py`, reflection `the-paley-cluster-integrals-are-catalan-numbers-tree-walks-and-the-moment-method.md`, OPEN-Q-013.
**RATE DIAGNOSTIC (monad-explorer, for handoff #2 / the p=31,43,47 compute node):** exact `R(p)` (p=3..23) give `(e−R)·p = 2.15,2.23,3.07,3.63,3.71` (climbing, decelerating — consistent with `e−R~K/p`, `K→≈3.8`, i.e. `R=e(1−C/p)`, `C≈1.4`) BUT `(e−R)·√p = 1.24,0.84,0.92,0.83,0.77` is actually FLATTER — so the 5 points do NOT cleanly distinguish `O(1/p)` from `O(1/√p)`. The a_4-sector cluster argument predicts a `+O(1/p)` term (`~p·A_4/p^5 = 2/p+…`), favoring `1/p`; settle + pin `C` at p≥31. This is HYP-2307 handoff #2 (the smooth analytic Paley signature) — now a *prediction to test*, not a blind extrapolation.
**RATE RESOLVED → 1/p (monad-explorer 3rd session, MISTAKE-060 + THM-438 ADDENDUM):** the cluster-integral error term is `A_{2k}=C_k p^{k+1}+O(p^k)` — i.e. `O(p^k)`, NOT the `O(p^{k+1/2})` originally claimed. Verified: `(A_4−2p^3)/p^2` is STABLE (≈−7.1…−7.8→≈−8) while `/p^{2.5}` drifts to 0. So the analytic correction to `R` is relative `O(1/p)`, confirming `R=e(1−C/p+…)` over `O(1/√p)`. (The √p-vs-1/p ambiguity above was a finite-data artifact; the Catalan-skeleton error term settles it analytically.) `C` still to be pinned at p≥31 — the `R=e(1−C/p)` ansatz is now justified, `C≈1.4` from current 5 points.

### INV-188: I(Ω(T), x) as Tournament Invariant
**Source:** opus-2026-04-05-S24 computation
**Status:** ESTABLISHED
**What:** The full polynomial I(Ω(T), x) is a strictly finer tournament invariant than H(T). Distinguishes tournaments at n=6 that share the same H. Real-rooted decomposition: H = Π(1 + 2/r_i) where -1/r_i are the zeros.
**Next:** Classify which iso classes are distinguished. Connection to Redei-Berge function?

---

## Priority A: Key structural questions (OCF PROVED by Grinberg-Stanley)

### INV-150: Simplicial Rédei Theorem (THM-220) — PROVED for all n ≥ 4
**Source:** opus-2026-03-15-S90
**Status:** COMPLETE. sim_H ∈ {0,1} for all tournaments on n ≥ 4 vertices.
**Key results:** Algebraic proof (Key Lemma + case analysis), near-transitive construction (H=2^{n-2}+1), OCF proof (Ω = K_{2^{n-3}}), β₁=1 for all near-transitive.
**Char poly:** p(λ) = λ^{n-2}(λ²+1) - (1+λ)^{n-2} (binomial coefficients!).
**Next:** Write as standalone paper. This is publishable.

### INV-151: The Cayley Monad D₄ Framework
**Source:** opus-2026-03-15-S90c/d/f
**Status:** ESTABLISHED. Q⁴=Id generates D₄ of order 8 acting on ℂP¹.
**Key results:** Cross-ratio of Q-orbit of 2 = 2 (the fugacity). Fixed points ±i. Bloch sphere connection: D₄ = Pauli group mod center. Q⁴(M)=M on transfer matrix. Spectral zeta ζ_M(-3)=7, ζ_M(-5)=21 (forbidden H!).
**Next:** Investigate D₄ action on tournament invariants. Possible paper.

### INV-152: The τ-φ Clock and Quasicrystal Structure
**Source:** opus-2026-03-15-S90h
**Status:** DISCOVERED. arg(λ_c)/π ≈ ln(2) (4 sig figs, NOT exact).
**Key results:** Gear ratio ≈ ln(2), period ≈ 2/ln(2), Tr mod 8 period = 8 exactly. Bott periodicity connection. Information interpretation: 1 bit per half-turn.
**Next:** Investigate why the approximation is so good. Connection to Rauzy fractal.

### INV-153: Tournament Equidecomposability (Hilbert's Third Problem analog)
**Source:** opus-2026-03-15-S90k
**Status:** PARTIALLY PROVED. β₁ as Dehn invariant. (H,β₁) complete at n=5.
**Key results:** 8 equidecomposability classes at n=5. Near-transitive = regular tetrahedron (non-decomposable). Within each (H,β₁) class, I(Ω₃,x) constant.
**Next:** Verify at n=6,7. Prove or disprove completeness for general n.

### INV-154: The Golden Shadow f(n) = (n-2+√(n²+4))/2
**Source:** opus-2026-03-15-S90o/p, kind-pasteur-S112
**Status:** EXPLORED. CF = [n-2; n,n,n,...]. Pell norm f·f'=-n. n-tribonacci family.
**Key results:** Unifies golden ratio, √2, silver ratio in one formula. f(n) satisfies λ²=(n-2)λ+n, transfer matrix [[n-2,n],[1,0]]. Memory correction T_n-M_{n-2} maximal at n=3.
**Next:** Investigate the n-tribonacci family more. Connection to Pisot substitution conjecture.

### INV-155: Tournament Wick Rotation and Spin-1 Ising
**Source:** opus-2026-03-15-S90g/i
**Status:** DISCOVERED. arctanh(2) = log(3)/2+iπ/2 → complex temperature.
**Key results:** H(T) = Wick-rotated partition function. ln(2) = renormalized criticality. Arrow's theorem = arctanh pole. Discriminant Δ(x) = 4x(x²-11x-1).
**Next:** Formalize the Ising correspondence. Yang-Lee zero connections.

### INV-141: T_19 Full Omega Dims — Degrees 9-18 Pending (Engineering Priority)
**Source:** kind-pasteur-2026-03-10-S54
**Status:** PARTIALLY COMPUTED. Degrees 0-8 done: [1,9,72,540,3753,23832,136260,688266,2987622].
**What:** Complete the T_19 Omega dim sequence to degrees 0-18. Degree 9+ requires C/C++ or LinBox.
**Key constraint:** chi per eigenspace = 1 (expected). Partial chi through m=8 is 2415061. Remaining (m=9..18) must sum to -2415060 alternating.
**Blocking issue:** Python dict-based pivot storage needs ~TB RAM for degree 9 (pivot density grows: max 10→26→95→455 at degrees 5→6→7→8).
**Next steps:**
  1. Implement in C/C++ (or use SageMath/LinBox) for degree 9+
  2. Or prove palindrome conjecture: Omega_{18-m} = Omega_m (gives free half via symmetry). T_7 IS palindromic; T_11 is NOT. T_19 unknown.
  3. Use chi constraint to check partial palindrome: if Omega_9 ≈ 12.95M (extrapolating ratio 4.34→3.6), check if sum through m=18 = 1.
**Priority:** MEDIUM (engineering interest, pattern understanding)
**Scripts:** t19_omega_dims_sparse.py, results at 05-knowledge/results/t19_omega_dims_sparse.out

### INV-142: Engineering Product: mod_rank Library (PyPI Target)
**Source:** kind-pasteur-2026-03-10-S53 (synthesis), S54 (update)
**Status:** MOSTLY DONE. mod_rank_library.py has core functions. 24 pytest tests written (test_mod_rank_library.py, all pass). Needs: README, benchmark vs dense/scipy, create setup.py.
**What:** General-purpose small-prime modular rank library. Key functions: gauss_rank_uint8, gauss_rank_nullbasis_uint8, certified_rank, betti_number_from_boundary_ranks. Memory table: T_11 deg 9 (52550×15745) = 6.6 GB int64 → 827 MB uint8.
**Engineering value:** Useful for any combinatorics/topology computation over finite fields with large sparse matrices.
**Next steps:** Write README, benchmark vs dense/scipy, create setup.py for PyPI

### INV-143: Engineering Product: circulant_homology Python Module
**Source:** kind-pasteur-2026-03-10-S54
**Status:** COMPLETE. circulant_homology.py has CirculantHomology and PaleyHomology classes. Fixed and verified.
**What:** Clean API for computing Omega dims and Betti numbers of circulant tournaments. Uses sparse column reduction (THM-125 eigenspace identity gives n× speedup). Verified for T_3, T_7, T_11.
**S55 fixes:** (1) betti_numbers() had wrong formula — fixed to use correct boundary map ranks via eigenspace computation. Verified T_3=[1,1,0], T_7=[1,0,0,0,6,0,0] exactly. (2) Added 27 pytest tests (test_circulant_homology.py, all pass). (3) Added caching for omega_basis_k and face_data for performance.
**New finding (HYP-453):** T_7 eigenspace structure: k=0 = H_0 only; k=1..6 = one H_4 generator each. T_11 eigenspace: all non-trivial homology at k=0.
**monad-compute (2026-06-03) from-scratch re-verification (`04-computation/verify_t11_betti_s_monad.py`, results `..._s_monad.out`/`.._NOTES.md`):**
  - Ω dims `[1,5,20,70,205,460,700,690,450,180,30]` RE-CONFIRMED from scratch (use_cache=False), χ=1. Raw |A_m| path counts = `[1,5,25,110,430,1430,3970,8735,14395,15745,8645]`. root-field prime=23.
  - Boundary ranks rank(d_m^(k)): **k=0** `[0,0,5,15,55,150,305,390,300,150,30,0]` (388s); **k=1** `[0,1,4,16,54,151,309,390,300,150,30,0]` (382s).
  - Per-eigenspace Betti: k=0 → `[1,0,0,0,0,5,5,..]`; k=1 → `[..,0,1,0..]` (+1 to β_6 only).
  - **Structural finding (REFINES HYP-453):** β_5=5 lives entirely at k=0; β_6 is distributed — k=0 carries 5, each non-principal eigenspace carries +1. Predicts β_6 = 5+10·1 = **15** (matches cached). HYP-453 "all T_11 homology at k=0" holds for β_5 but NOT β_6.
**monad-compute-2026-06-03-S2 — FULLY COMPLETE from scratch (all 11 eigenspaces):**
  - **Cache-bloat fixed:** clearing `_omega_basis_cache` per eigenspace held timing flat at 384–413 s/eigenspace (no more >12 min slowdown). Script now persists each eigenspace to `verify_t11_betti_s_monad_ranks.json` (resumable) + auto-commit.
  - **All k=1..10 have IDENTICAL boundary ranks** `[0,1,4,16,54,151,309,390,300,150,30,0]`; k=0 `[0,0,5,15,55,150,305,390,300,150,30,0]`.
  - **β = `[1,0,0,0,0,5,15,0,0,0,0]` CONFIRMED.** Per-eigenspace β_6 = `[5,1,1,1,1,1,1,1,1,1,1]` → 5+10·1=15; β_5 = `[5,0,...]` → 5 (k=0 only). The "+1 to β_6 per non-principal eigenspace" pattern holds for ALL 10.
  - **Euler check clarification:** the script's first "OVERALL MISMATCH" was a FALSE ALARM — it compared χ_Betti=11 against single-copy χ_Omega=1. By THM-125 each of the n=11 eigenspaces carries a FULL copy of Ω_m, so the correct identity is χ_Betti = n·χ_Omega = 11·1 = 11. ✓ Fixed in the script. The Betti numbers were always correct (reproduce the library's official `betti_numbers()` formula exactly).
**STATUS: β_5=5, β_6=15 mechanically re-verified from scratch across all 11 eigenspaces. INV-143 T_11 Betti re-verification CLOSED.** Remaining open engineering item: C/LinBox reimplementation for routine degree-9+ re-verification (cf. INV-141).

### INV-192: Engineering Product: Odd-Cycle Disjointness Features for Tournament TDA
**Source:** opus-2026-05-29-S11
**Status:** PARTIALLY IMPLEMENTED (opus-2026-05-29-S15)
**What:** S11's H=63 / THM-025 comparison suggests that practical tournament fingerprints should expose the disjointness geometry of odd cycles, not only H, scores, and Betti numbers. Candidate features: `alpha_vector(Omega)`, cycle-family core size, complete-Ω flag, disjoint-pair count, independent-triple supports, single-core signature, `r_core(s)`, and deletion profile `(H(T-v), |Omega(T-v)|, complete?)`.
**S12 update:** Extend this feature layer to projection defects: support count vs support excess, max support multiplicity, deletion loss profile `(lost, kept, loss_frac, alpha(T-v))`, exact-kill/near-kill flags, and even-graph projection weight/degree sequence for odd n.
**S15 update:** `04-computation/tournament_tda.py` now computes exact `omega_*` features for n≤9: `omega_alpha_1`, `omega_alpha_2`, `omega_alpha_3`, complete-Ω flag, disjoint-pair density, support count/excess, max support multiplicity, cycle-family core size, projection-kill vertex count, and max deletion-loss fraction. These are included in the flat ML feature vector and demo output is saved in `05-knowledge/results/tournament_tda_omega_features_s15.out`.
**Why it matters:** These features separate two phenomena that H alone compresses: (1) H=63 unlocks through a complete-core Ω=K31; (2) real-rootedness fails through a no-core, highly concentrated independent triple. For ranking data, this distinguishes "all inconsistency localized at one pivot" from "three disjoint inconsistency groups with lopsided coupling."
**Next steps:**
  1. Benchmark the implemented features on synthetic rankings, sports/election data, and attention-derived tournaments.
  2. Use them as prefilters before expensive full Ω or path-homology computations.
  3. Add odd-n even-graph projection signatures and deletion `alpha(T-v)` profiles as a second feature layer.
**Files:** `04-computation/tournament_tda.py`, `04-computation/omega_extreme_fingerprints_s11.py`, `04-computation/projection_defect_bridge_s12.py`, `05-knowledge/results/tournament_tda_omega_features_s15.out`, `05-knowledge/results/omega_extreme_fingerprints_s11.out`, `05-knowledge/results/projection_defect_bridge_s12.out`, `07-reflections/omega-extremes-as-cycle-disjointness-axis.md`, `07-reflections/projection-defect-as-common-residue.md`.

### INV-193: Projection-Defect Axis Across Ω, Even Graphs, and Path Homology
**Source:** opus-2026-05-29-S12
**Status:** NEW META-STRUCTURAL LEAD
**What:** Several mature threads may be instances of the same question: what residue survives a forgetful projection? S12 compared three projections:
  1. Directed odd cycles → vertex supports (multiplicity defect).
  2. Vertex deletion / old-coordinate projection (cycle loss, projection kill).
  3. Tournament orientation → degree-even cycle-space graph, well-defined only at odd n.
**Key examples:** H=63 is an exact old-projection kill: deleting the core removes all 31 odd cycles. THM-025 is a near-kill: one vertex supports 92/94 odd cycles, but the two surviving old cycles carry `alpha=[1,2,1]`, enough to keep the real-root failure alive. Paley T7 and interval T7 have the same 36 odd-cycle supports but different support-excess and even-graph projection weights. Path-homology HYP-408/ghost cycles asks the same structural question in chain-complex language: when do through-v-only cycles become boundaries after old projection?
**Engineering angle:** Add projection-defect features to Tournament TDA: max deletion cycle-loss fraction, old-projection kill vertices, support multiplicity defect, and odd-n even-graph projection signature. These are cheap fingerprints that may prefilter root failures, localized inconsistency, and homology anomalies.
**Next steps:**
  1. Scan known THM-025-like non-real-root examples for high max deletion loss.
  2. Compare projection-defect statistics against beta_3/beta_4 path-homology anomalies.
  3. For odd n, correlate even-graph projection fibers with Ω support-multiplicity defect.
  4. Extend the S15 `tournament_tda.py` `omega_*` feature block with explicit deletion-residue alpha profiles and even-graph projection signatures.
**Files:** `04-computation/projection_defect_bridge_s12.py`, `05-knowledge/results/projection_defect_bridge_s12.out`, `07-reflections/projection-defect-as-common-residue.md`, HYP-1760, HYP-1763, T282, T283.

### INV-194: Merged Tiling Bucket Constraints
**Source:** kind-pasteur-2026-05-29-S5
**Status:** THEOREM + BOOLEAN-CUBE LEAN SPECIALIZATION COMPLETE + OPEN STRUCTURAL EXTENSION
**What:** Treat quotient maps out of the tiling cube as bucket maps. THM-346 proves the general half-line balance law for any quotient `q: Q_m -> B` and any mask family `M`: `2*self_b + incident_cross_b = |q^{-1}(b)|*|M|`. THM-345 specializes this to `pi: Q_m -> G_n/Z_2`, proving that bucket size parity detects SC/NS type exactly: SC buckets are odd, NS buckets are `2 mod 4`, and `sum_M B_M=2^m`. For every Hamming layer `d`, the ordered transport matrix `W_d(M,N)` is symmetric, has row sums `B_M*C(m,d)`, has even diagonal, and therefore has forced cross-outflow parity. Lucas' theorem says the active parity layers are exactly the binary submasks of `m=C(n-1,2)`.
**S6 Lean update (kind-pasteur-2026-05-29-S6):** The local good-cut bucket gap was strengthened in `TournamentH7.GoodCuts`: nonempty good-cut support is equivalent to the existence of an upward tile, any upward tile forces at least two distinct good cuts, and every tiling satisfies `goodCutCount = 0 ∨ 2 ≤ goodCutCount`. The `TournamentH7.Verify` audit confirms this core depends only on Lean foundations.
**S1 Lean update (kind-pasteur-2026-05-30-S1):** Added THM-348 / `TournamentH7.BucketBalance`, a generic axiom-free Lean proof of the oriented finite-set balance `|selfHalf|+|crossHalf|=|fiber|*|moves|`, plus the zero-cross closure criterion. This isolates the remaining Lean bridge for THM-346: prove the fixed-point-free involution pairing that turns internal half-lines into unordered self-lines counted twice.
**S2 Lean update (kind-pasteur-2026-05-30-S2):** Added THM-350 / the unordered abstract layer in `TournamentH7.BucketBalance`: the partner map `(x,u)->(step u x,u)` preserves internal half-lines for involutive moves, fixed-point-free moves forbid self-partners, and the unordered balance follows from `Even selfHalf.card`. Also added the concrete staircase-to-tournament bridge and top-good-cut/SC audit via `StaircaseConnectivity`.
**Codex 2026-05-30 orbit update:** The generic orbit-parity bridge is now Lean-proved. `BucketBalance.even_card_of_fixedPointFree_involutiveOn` proves that any finite fixed-point-free involution has even cardinality, and `BucketBalance.unordered_balance_of_involutive_fixedPointFree` removes the separate evenness assumption from the abstract unordered balance.
**Opus 2026-05-30 Boolean-cube update:** THM-351 closes the remaining Boolean-mask bridge in Lean. `BoolCube`, `xorMask`, `xorMask_involutive`, `xorMask_fixedPointFree_of_nonzero`, and `unordered_balance_boolCube_masks` specialize the abstract THM-350 layer to finite Boolean cubes with nonzero xor masks. This gives a reusable formal model for tiling-cube mask families before quotient-specific tournament maps are attached.
**Why it matters:** This turns the merged metagraph into a constrained reversible transport system, not just an unweighted graph. It also corrects the old S202 "merged tiling excess" narration: merged buckets still partition the fixed-base cube exactly.
**Next steps:**
  1. Attach the Boolean-cube theorem to the concrete tiling-coordinate type used by the staircase explorer, if a more semantic wrapper is wanted.
  2. Determine whether generic NS-sea transport is approximable from bucket sizes alone, with SC/rib corrections.
  3. Seek a Burnside formula for the bucket-size distribution, not just for the number of buckets.
  4. Promote the Boolean cube/mask API into the stable Lean hierarchy if later modules need it.
  5. Add bucket parity and normalized `W_d` features to the future `tournament_tda.py` extractor.
**Files:** `01-canon/theorems/THM-336-good-cuts-structure.md`, `01-canon/theorems/THM-345-merged-bucket-parity.md`, `01-canon/theorems/THM-346-tiling-quotient-bucket-balance.md`, `01-canon/theorems/THM-348-finite-bucket-halfline-balance.md`, `01-canon/theorems/THM-350-finite-unordered-bucket-balance-layer.md`, `01-canon/theorems/THM-351-boolean-cube-mask-bucket-balance.md`, `04-computation/merged_bucket_constraints_s5.py`, `04-computation/tiling_quotient_bucket_balance_s5.py`, `04-computation/lean/TournamentH7/TournamentH7/GoodCuts.lean`, `04-computation/lean/TournamentH7/TournamentH7/BucketBalance.lean`, `04-computation/lean/TournamentH7/TournamentH7/StaircaseConnectivity.lean`, `05-knowledge/results/merged_bucket_constraints_s5.out`, `05-knowledge/results/tiling_quotient_bucket_balance_s5.out`, `05-knowledge/results/lean_goodcuts_bucket_strengthening_kind_pasteur_s6.out`, `05-knowledge/results/lean_bucket_balance_kind_pasteur_2026-05-30-S1.out`, `05-knowledge/results/lean_verify_unordered_kind_pasteur_2026-05-30-S2.out`, `05-knowledge/results/lean_boolcube_bucket_balance_opus_2026-05-30-S1.out`, `05-knowledge/results/lean_verify_boolcube_bucket_balance_opus_2026-05-30-S1.out`, `05-knowledge/variables/merged-bucket-size.md`, `05-knowledge/variables/tiling-bucket-balance.md`, `07-reflections/merged-tiling-bucket-constraints.md`, `07-reflections/unordered-bucket-balance-orbits.md`, `07-reflections/boolean-cube-balance-as-checksum.md`.

### INV-195: Lean Import Narrowing for Fast Distributed Formalization
**Source:** kind-pasteur-2026-05-29-S6
**Status:** NEW ENGINEERING LEAD
**What:** A cold build of `TournamentH7.GoodCuts` pulled in the broad `Mathlib.Tactic` umbrella and compiled 2956 targets before the small local theorem file built. The module itself only uses elementary finset/cardinality facts plus arithmetic automation (`omega`, `simp`, `rw`, `Finset.card_bij`, `Finset.card_pos`, `Finset.card_eq_zero`, `Finset.card_le_card`). Narrowing imports in proof-heavy files could materially improve agent turnaround time and reduce wasted compute in fresh workspaces.
**Why it matters:** Formalization is now an active project deliverable. Faster cold builds make multi-agent Lean work more practical and reduce the cost of CI or new-machine onboarding.
**Next steps:**
  1. Replace broad `import Mathlib.Tactic` in `GoodCuts.lean` with the minimal tactic/data imports and verify a cold-ish build.
  2. Scan other Lean modules for broad imports that can be narrowed without churn.
  3. Document the minimal import pattern in `04-computation/lean/TournamentH7/README.md`.
**Files:** `04-computation/lean/TournamentH7/TournamentH7/GoodCuts.lean`, `05-knowledge/results/lean_goodcuts_bucket_strengthening_kind_pasteur_s6.out`.

### INV-135: Tang-Yau (arXiv:2602.04140): Path Homology of Circulant Digraphs via Fourier
**Source:** Web research kind-pasteur-2026-03-10-S50
**Status:** INVESTIGATED (kind-pasteur-2026-03-10-S50). February 2026.
**What:** Tang-Yau develop a "symbol-matrix" approach to GLMY path homology for circulant digraphs using Fourier decomposition on the shift automorphism. Provides computable Betti number formulas, with dependence on whether n is prime or composite.
**Key findings (from paper fetch):**
- Main theorems: C→_5^{1,2} and C→_n^{1,s} with s≠2 have β_1=1,β_2=1, otherwise β_1=1.
- Symbol matrix M_m(t) reduces rank to evaluating Laurent polynomials at roots of unity.
- Stability theorem (Thm 1.4): for large primes p∉Q+(S), Betti numbers stabilize.
- DOES NOT compute Betti numbers for Paley T_p (connection sets of size (p-1)/2).
- DOES NOT have Omega palindrome or top vanishing results for arbitrary connection sets.
**Impact:** The symbol-matrix framework is the right tool for T_p, but they leave Paley application open. The Stability Theorem (Thm 1.4) implies our pattern β_6(T_p)=p-1 holds for ALL large p∉Q+(QR_p). Computing Q+(QR_p) would be the key step.
**Next step (UPDATED S52):** COMPLETED. Tang-Yau symbol matrix applied to T_7 (deg 2-5) and T_11 (deg 2-5). Key discovery: Q+(QR_p) is EMPTY (not just "avoids p-th roots" — ALL t values work equally). Symbol matrix M_m(t) is CONSTANT (THM-125: proven algebraically). Eigenspace identity is trivial consequence. New open: prove Q+(QR_p) empty algebraically for all p.
**S71c UPDATE:** Full HTML paper read. Confirms our eigenspace approach matches their Prop 4.1. Their Conjecture 4.8 (H_m=0 for m≥3, no-wrap-around S) does NOT apply to Paley (QR has wrap-around). Our HYP-824 Betti concentration at d=m,m+1 is genuinely new. Key open: prove β_d^(0)=0 for 3≤d≤m-1 (the intermediate vanishing).
**Priority:** CLOSED for computational part. Proving intermediate vanishing is Priority A.

### INV-136: Schweser-Stiebitz-Toft (arXiv:2510.10659): Redei's Theorem Revisited (Oct 2025)
**Source:** Web research kind-pasteur-2026-03-10-S50
**Status:** INVESTIGATED (kind-pasteur-2026-03-12-S55).
**What:** 6-page paper unifying three equivalent stronger forms of Rédei's theorem: (i) Rédei's stronger theorem (adding undirected vertices creates even path count), (ii) Berge's theorem (complementary mixed graphs have same parity), (iii) Dirac's theorem (Hamiltonian paths through non-oriented edges are even). All three are interconnected.
**Relevance:** Directly relevant to H(T) parity. The mixed-graph generalization provides a structural lens complementary to OCF. THM-016/017 (our even-odd split) and DC (THM-082/083) could potentially be re-derived via Berge/Dirac. Possible route to extend OCF to mixed digraphs.
**Action:** Add to TANGENTS as "Mixed graph OCF extension via Schweser-Stiebitz-Toft stronger Rédei".

### INV-137: Satake (arXiv:2502.12090): Cyclotomic Nearly-Doubly-Regular Tournaments (Feb 2025)
**Source:** Web research kind-pasteur-2026-03-10-S50
**Status:** DEFINITIVELY RESOLVED (kind-pasteur-2026-03-12-S56). Satake NDRTs do NOT maximize H.
**What:** For prime powers q ≡ 5 (mod 8), cyclotomic tournament CT_q is NDR iff q = s² + 4. When true: full adjacency spectrum computed explicitly (eigenvalues (q-1)/2 and (-1 ± i√(q ∓ 2√q))/2). Under Hardy-Littlewood conjecture F: infinitely many such q. This is the n ≡ 1 (mod 4) analog of Paley (n ≡ 3 mod 4).
**RESOLUTION (S56 satake_analysis_ext.py, exhaustive at n=13):**
- q=5: H_sat=15 = H_max (trivially tied — ALL 4 circulants tie at n=5)
- q=13: H_sat=3,703,011, rank 40/64 — FAR from maximum (gap=8,164, HYP-456 REFUTED)
- Maximizer at n=13: cyclic interval S={7,...,12}, H=3,711,175 (unique, rank 1/64) — HYP-455 CONFIRMED
- NDR property does NOT predict H-optimality for q≡5 mod 8
**Pattern discovered (HYP-455):** At Paley primes p≡3 mod 4 → Paley maximizes. At q≡5 mod 8 primes → cyclic interval maximizes.
**Next step:** Verify cyclic interval pattern at q=29 (needs sparse/C++ Held-Karp). q=29 currently hits MemoryError in Python.

### INV-138: Ren (arXiv:2504.15126): Path Independence Complexes of Digraphs (Apr 2025)
**Source:** Web research kind-pasteur-2026-03-10-S50
**Status:** INVESTIGATED (kind-pasteur-2026-03-12-S55). HIGH PRIORITY — closely adjacent to our work.
**What:** Studies "path independence complexes" of digraphs — simplicial complexes whose faces are vertex sets with no directed path between any two. Uses GLMY infrastructure. Main results: (i) canonical embeddings from independence complex (undirected) into path independence complex (digraph), (ii) these are Σ_k-equivariant and isometric giving double-parametrized persistent homology, (iii) Shannon capacity consequences.
**NOT the same as Omega(T):** Omega(T) is the conflict GRAPH of directed cycles (edges = two cycles sharing a vertex). The path independence complex is a SIMPLICIAL COMPLEX on vertices with different independence condition. However, both are topological structures related to the tournament and use GLMY machinery.
**Key connection:** The embedding theorem (independence complex → path independence complex) might explain the tight relationship between Omega(T) and tournament path homology. Our beta_2=0 result and seesaw mechanism (THM-095) may follow from an embedding theorem structure.
**Next step:** Read sections on Betti number bounds and explicit computations. Compare beta numbers from Ren's framework to our beta_2=0, beta_3 onset results. This is potentially a theoretical unification.

### INV-139: Tang-Yau (arXiv:2402.05682): Cellular Homology of Digraphs (Feb 2024)
**Source:** Web research kind-pasteur-2026-03-10-S50
**Status:** NEW — background.
**What:** Introduces cellular chain complexes for digraphs using admissible minimal paths. Connects to GLMY under strongly regular conditions.

### INV-140: Hepworth-Roff (arXiv:2404.06689): Magnitude-Path Spectral Sequence (Apr 2024)
**Source:** Web research kind-pasteur-2026-03-10-S50
**Status:** NEW — background.
**What:** The MPSS connects magnitude homology (page 1) to path homology (page 2). Bigraded path homology satisfies Kunneth, Mayer-Vietoris, excision. Could provide extra structure for tournament path homology via the bigraded version.

### INV-144: Awan-Bernardi B-polynomial for Digraphs (arXiv:1610.01839)
**Source:** Web research kind-pasteur-2026-03-14-S69
**Status:** NEW — high priority connection
**What:** Awan and Bernardi defined a Tutte polynomial for directed graphs satisfying deletion-contraction B(D) expressed via B(D\a) and B(D/a). Published JCTB 2020. Our F_T(x) = F_{T\e}(x) + (x-1)*F(T/e,x) is structurally the same type of relation. The B-polynomial is a 3-variable polynomial detecting acyclicity, strong connectivity, and directed paths.
**Next step:** Read paper, check if F(T,x) is a specialization of their B-polynomial.

### INV-145: Sazdanovic-Yip Categorification of Chromatic Function (arXiv:1506.03133)
**Source:** Web research kind-pasteur-2026-03-14-S69
**Status:** NEW — highest priority for categorification
**What:** Sazdanovic and Yip categorified the chromatic symmetric function using Khovanov-style techniques. Combined with Mitrovic-Stojadinovic chromatic-Redei-Berge bridge (arXiv:2506.08841), this suggests categorifying F(T,x) via long exact sequences in homology, producing "Khovanov homology for tournaments."
**Next step:** Read and check if the deletion-contraction for F(T,x) lifts to a long exact sequence in GLMY path homology.

### INV-146: Asao Magnitude-Path Spectral Sequence (arXiv:2201.08047)
**Source:** Web research kind-pasteur-2026-03-14-S69
**Status:** NEW — connects magnitude and path homology
**What:** Asao proved magnitude homology and GLMY path homology are pages of the same spectral sequence. Published Bull. London Math. Soc. 2023. Magnitude homology categorifies magnitude (Hepworth-Willerton, arXiv:1505.04125), analogous to Khovanov categorifying Jones. This creates an indirect chain: Khovanov ↔ magnitude homology ↔ path homology.
**Next step:** Understand what the spectral sequence looks like for tournaments.

### INV-147: Hepworth Reachability Homology (arXiv:2312.01378)
**Source:** Web research kind-pasteur-2026-03-14-S69
**Status:** NEW — unifies magnitude and path homology
**What:** Hepworth defined reachability homology of digraphs, unifying magnitude and path homology. Satisfies homotopy invariance, Kunneth, excision, Mayer-Vietoris. Published IMRN 2025. Potentially the "right" homology theory for tournaments.
**Next step:** Compute reachability homology for small tournaments, compare to GLMY.

### INV-050: Fourier Decomposition Proof of OCF — OCF PROVED AT n=5 AND n=7
**Source:** opus-2026-03-06-S11b (continued^7, ^8)
**Status:** OCF PROVED AT n=5 AND n=7 via Fourier decomposition. All identities at both n proved.
**What:** OCF decomposes into independent degree-homogeneous Fourier identities:
- **Fourier Homogeneity Theorem:** w_{n-1-2k} is a homogeneous polynomial of degree 2k in centered edge variables s_e = A_e - 1/2.
- **Degree 0:** Trivial (expected values). PROVED for all n.
- **Degree 2:** Proved via proportionality constants c_{2j+1} = C(n-3,2j-2)*(2j-2)!/2^{2j-2}. PROVED for all n.
- **Degree n-1:** Proved via path-cycle bijection: w_0 = 2*[deg-(n-1) of t_n]. PROVED for all n.
- **Degree 4 at n=7:** PROVED via counting lemmas. The degree-4 Fourier space is 2-dimensional:
  - Type P: 5-vertex spanning paths (coefficients ±1 in t5_d4)
  - Type Q: 6-vertex disjoint P₂ pairs (coefficients ±1 in α₂_d4)
  - [deg-4 of t₇] = (1/2)·[deg-4 of t₅] + [deg-4 of α₂] (EXACT, verified all 5985 monomials)
  - w₂/4 = 3·[deg-4 of t₅] + 6·[deg-4 of α₂] (EXACT, verified all 5985 monomials)
  - Counting: c₅=2, c₇=4, paths=12 (type P); c₇=8, paths=24, a₂=1 (type Q)
**Key insight:** Fourier supports of t₅_d4 and α₂_d4 are DISJOINT (spanning paths vs P₂ pairs), reducing the identity to independent counting arguments.
**Scripts:** `04-computation/fourier_homogeneity.py`, `fourier_degree2_identity.py`, `ocf_fourier_proof_framework.py`, `degree4_identity_n7.py`, `degree4_proof_n7.py`
**CRITICAL FINDING (opus-2026-03-07-S37):** At n=9, the degree-4 Fourier space has dimension **>> 200** (no saturation at 200 random tournaments with 300 probe monomials). Within the P4 type alone, coefficients are NOT proportional (117/126 vertex sets have non-constant |coeff|). The |coeff| values range from 1 to 27. This means the n=7 "two type" decomposition DOES NOT generalize. The Fourier approach is **infeasible at n=9** for middle degrees.
**Scripts:** `04-computation/degree4_n9_rank.py`, `degree4_n9_rank2.py`, `degree4_n9_saturation.py`
**Next step:** (1) The Fourier proof cannot extend to n>=9 for middle degrees. Focus on algebraic approaches (OCF already proved by Grinberg-Stanley). (2) The degree-0, degree-2, and degree-(n-1) identities still hold for all n and have clean proofs. Can they be combined differently?

### INV-123: THM-086 Universal Taylor Zeros mod 3 — PROOF SKETCH COMPLETE
**Source:** kind-pasteur-2026-03-07-S37
**Status:** PROOF SKETCH COMPLETE. Verified n=5-10, inductive structure identified.
**What:** c_j(T) = 0 mod 3 for all tournaments T on n vertices and all j < val(n), where val(n) = 2*floor((n-1)/2). This means (x-1)^{val(n)} | F(T,x) mod 3. For n odd, F(T,x) mod 3 is determined by a SINGLE parameter alpha = c_{n-1}(T) mod 3.
**Proved cases:** j=0,1,2 (THM-085, algebraic). j=3 (palindrome + THM-085). j>=4 (DC induction + palindrome, verified computationally).
**Key corollary:** Eulerian conjecture: 3|A(n,k) => 3|F_k(T) for all T. Follows from (x-1)-adic valuation of A_n(x) mod 3 being exactly val(n).
**What remains:** The "almost-tournament claim" — c_j(T\e) = 0 mod 3 for j < val(n)-1 — needs formal proof, likely via nested DC induction. Verified exhaustively at n=5, sampled at n=6-8.
**Scripts:** `04-computation/thm086_verify.py`, `dc_induction_proof.py`, `c4_induction_test.py`, `taylor_cj_mod3_analysis.py`, `eulerian_zeros_from_palindrome.py`
**Next step:** (1) Prove almost-tournament claim algebraically (N_uv formula reduces it to Taylor zeros of the "adjacent pair" polynomial). (2) Extend to mod 9. (3) Mod p for p>=5 INVESTIGATED (S38): universal zeros match Eulerian val for n >= p+2 but Eulerian conjecture FAILS for p>=5 (multiple free parameters).

### INV-124: THM-094 F_k mod 2 Tournament-Independent — PROOF SKETCH COMPLETE
**Source:** kind-pasteur-2026-03-07-S38
**Status:** PROOF SKETCH COMPLETE. Verified exhaustively n<=6, sampled n=7,8.
**What:** F_k(T) = A(n,k) = C(n-1, k) mod 2 for ALL tournaments T. F(T,x) = (1+x)^{n-1} mod 2 is COMPLETELY tournament-independent. Proof via universal Taylor zeros mod 2 (c_j = 0 for j < n-1) + Redei's theorem (F_{n-1} = Hamiltonian path count is always odd). The mod-2 result is the strongest possible: individual F_k are determined, not just linear combinations.
**Key insight:** p=2 is special because (1) val_2(A_n(x)) = n-1 (maximal), giving a single free parameter, and (2) Redei pins that parameter to 1.
**Mod-p generalization (S38):** For p >= 5, universal Taylor zeros match Eulerian valuation only for n >= p+2. The Eulerian conjecture (p|A(n,k) => p|F_k(T)) FAILS for p=5 at n=7 because multiple free parameters in F(T,x) mod 5 allow different zero patterns.
**Scripts:** `04-computation/fk_mod2_proof.py`, `taylor_zeros_mod_p.py`, `mod_p_general_conjecture.py`
**Next step:** (1) Prove universal Taylor zeros mod 2 algebraically (c_j = 0 for j < n-1). (2) Is there an elementary proof not using THM-086 machinery?

### INV-032: Omega(T) structural properties — PARTIALLY DISPROVED
**Source:** Web research opus-S5, opus-S7 (disproof), opus-S9 (line graph disproof), opus-S10 (structure analysis)
**Status:** DISPROVED: Omega(T) is NOT always claw-free (fails n=9, 90%) or perfect (fails n=8, 53.8%). NOT a line graph (K_5-e found at n=6, 45%). S_{1,1,1}-free through n=11, fails n=12.
**DISPROVED (THM-025, opus-S18):** Real-rootedness of I(Omega(T),x) FAILS at n=9. Counterexample: score [1,1,3,4,4,4,6,6,7], I=[1,94,10,1], Newton k=2 fails (100 < 141), complex roots confirmed.
**Structural characterization (opus-S19):** Failure requires (a) three vertex-disjoint 3-cycles partitioning V, AND (b) near-total inter-group domination (9-0, 9-0, 7-2 arc counts), creating a transitivity bottleneck with hub vertex in 92/94 cycles. The extreme tournament (full domination) gives I=(1+x)^3 with disc=0 exactly. One arc flip creates disc<0. Failure is MAXIMALLY RARE: 0 in 10000 random samples at n=9, 0 at n=10,11.
**What remains true:** Real-rootedness holds at n<=8 (Chudnovsky-Seymour, claw-free) and for "generic" tournaments at all n. But it is NOT a universal structural property.
**Next step:** (1) Characterize exact tournament class where failure occurs. (2) Check if the clique-deletion interlacing approach (INV-038) can prove real-rootedness under a claw-free assumption.
**What remained true:** All-real-roots of I(Omega(T), x) appears to hold even for imperfect/non-claw-free Omega (tested n<=10, 0 failures). This is a deep structural conjecture NOT explained by any forbidden subgraph property.
**Note:** OCF is now proved by Grinberg-Stanley, so this is no longer a proof strategy — it's a structural question. Real-rootedness explanation must be algebraic (Irving-Omar/Grinberg-Stanley symmetric function framework).
**Extended testing (opus-S18):** Real-rootedness tested for I(Omega_3(T), x) at n=9-21 with 0 failures across 1470+ samples (degrees up to 5). Log-concavity and Newton's inequalities hold in all cases. The "Omega_3 complement = matching" structure holds exhaustively at n≤6 (31088/31088) but fails at n≥7 (75.3%).
**Turán-based proof for n≤11:** At n=9-11, alpha(Omega_3) = 3, so the disjoint-pair graph is triangle-free. Turán gives a2 ≤ c3²/4, proving Newton's first inequality a1² ≥ 3a2. Combined with the degree-3 discriminant bound, this could give a complete proof at n≤11. For n≥12, Turán alone fails.
**Next step:** (1) Complete Turán+discriminant proof for n=9-11. (2) Find tournament-specific bounds on a2 for n≥12. (3) Investigate Irving-Omar determinantal formula for algebraic proof.

### INV-038: Clique-deletion interlacing for Omega(T)
**Source:** opus-2026-03-06-S17, T100, interlacing-clique-deletion.md
**Status:** STRUCTURAL INSIGHT. Proof sketch for n<=8.
**What:** Through-v cycles always form a CLIQUE in Omega(T) (proved: sharing vertex v). Deleting vertex v = deleting this clique from Omega(T). The sequential deletion recurrence I(G,x) = I(G-u,x) + x*I(G\N[u],x) can be applied step-by-step. At n=5, 100% of remaining cycles are adjacent to some through-v cycle (Omega is very dense).
**Key insight:** For n<=8 (claw-free), Chudnovsky-Seymour guarantees each step preserves real roots. For n>=9, real-rootedness can FAIL (THM-025), so the interlacing approach cannot extend to all tournaments. However, 84 claws at n=9 counterexample all share the same 3 leaves — the claw structure is very specific.
**Verification:** 0 failures: n=5 (5120 exhaustive), n=6 (196608 exhaustive), n=7-8 (random).
**Impact:** Proves real-rootedness for n<=8. For n>=9, would need a tournament-specific claw-free condition.
**Next step:** (1) Characterize when Omega(T) is claw-free at n>=9. (2) Check if "generically claw-free" suffices for applications.
**Scripts:** `04-computation/interlacing_verify.py`, `04-computation/interlacing_structure.py`
**Writeup:** `03-artifacts/drafts/interlacing-clique-deletion.md`

### INV-039: Blueself odd-n obstruction — PROVED for ALL odd n
**Source:** opus-2026-03-06-S17, THM-022 Theorem 5 (upgraded)
**Status:** PROVED. Pure algebraic proof, no exhaustive search needed.
**What:** No blueself tilings exist at any odd n. Grid-symmetry forces k_0+k_{n-1}=n-2 (endpoint constraint). Flip changes endpoint multisets: {1+k_0, n-2-k_0} -> {n-1-k_0, k_0}. For these to be equal as multisets, need k_0=(n-2)/2 (non-integer at odd n) or 1=0 (impossible). Therefore sorted scores always differ, so flip(T) is never isomorphic to T.
**Script:** `04-computation/blueself_odd_n_proof.py`
**Impact:** Upgrades THM-022 Theorem 5 from "proved n<=7" to "proved all n". Completes the odd half of the blueself existence dichotomy.

### INV-040: Blueself vs SC maximizer — DISPROVED
**Source:** opus-2026-03-06-S17, T099
**Status:** DISPROVED at n=6.
**What:** At n=6, blueself class with H=41 is NOT the SC maximizer in score class (3,3,3,2,2,2) (SC max is H=45, also blueself). Blueself classes are always SC and have regular scores, but not always max-H. The blueself with higher disjoint pair count (alpha_2=4) beats the one with more total cycles (alpha_1=16).
**Script:** `04-computation/blueself_sc_maximizer_connection.py`

### INV-041: Quasi-regularity of Omega(T) — EXPLAINED
**Source:** opus-2026-03-06-S17 (T101), opus-2026-03-06-S18 (T103, proof)
**Status:** EXPLAINED. Theoretical argument + verified n=5-20.
**What:** Omega_3(T) is quasi-regular because adjacency depends on vertex-set intersection (sharing ≥1 vertex), not arc orientations. This makes Omega_3 an induced subgraph of J(n,3) (Johnson graph), inheriting its regularity. All 3-element subsets have identical intersection statistics, so degree of each 3-cycle concentrates around E[deg] = (C(n,3)−C(n−3,3))/4−1. The coefficient of variation CV = O(1/√m) → 0 as n→∞, giving λ_max/avg_deg ≈ 1+CV² → 1. Verified: CV drops from 0.05 (n=6) to 0.03 (n=20). This does not directly explain real-rootedness but constrains spectral structure.
**Scripts:** `04-computation/omega_spectral_fast.py`, `04-computation/omega_quasireg_proof.py`

### INV-042: Paley deletion maximizer — VERIFIED p=3,7,11
**Source:** kind-pasteur-2026-03-06-S18e (T097), opus-2026-03-06-S18
**Status:** VERIFIED at p=3,7,11. Conjecture for all Paley primes.
**What:** H(T_p − v) = a(p−1) (OEIS A038375 max H at n=p−1). Verified: T_3−v: H=1=a(2), T_7−v: H=45=a(6), T_11−v: H=15745=a(10). By vertex-transitivity all deletions equivalent. Combined with T053 (T_p achieves a(p)), the maximizer chain is "hereditary" via vertex deletion. Claim A decomposition for T_7: diff=144=2×72, sum_mu=6×3+30+24=72 (all 3-cycle complements have a 3-cycle in Paley).
**Next step:** Verify at p=19 (need H(T_19−v) = a(18)). Investigate: does the n=p−1 maximizer always come from Paley deletion, or can it be achieved by non-Paley tournaments too?
**Scripts:** `04-computation/paley_deletion_test.py`

### INV-043: Anti-aut involution existence — PROVED (THM-024)
**Source:** opus-2026-03-06-S18 (T102, THM-024), correcting kind-pasteur S18e (T095)
**Status:** PROVED. Clean group theory argument.
**What:** Every SC tournament has ≥1 involution anti-automorphism. Proof: (1) Moon's theorem: |Aut(T)| is odd. (2) H = ⟨Aut(T), σ₀⟩ has order 2|Aut(T)| (even). (3) By Cauchy, H has order-2 element. (4) Can't be in Aut(T) (odd order group). (5) Must be in σ₀·Aut(T) = set of anti-auts. NOT all anti-auts are involutions (counterexamples at n=6 with |Aut|>1), but at least one always is.
**Scripts:** `04-computation/anti_aut_involution_test.py`, THM-024

### INV-044: Hereditary Maximizer Chain — CORRECTED (regular-only at odd n)
**Source:** kind-pasteur-2026-03-06-S18f, S18g (correction), T104, T105, MISTAKE-010
**Status:** CORRECTED. Only REGULAR maximizers at odd n are hereditary.
**What:** Previous claim "all maximizers at odd n hereditary" was WRONG. Exhaustive check:
- n=3: 2/2 hereditary (all regular)
- n=4: 0/24 hereditary (non-regular)
- n=5: 24/64 hereditary (only 24 regular, NOT 40 with score (1,2,2,2,3))
- n=6: 0/480 hereditary (non-regular)
- n=7: 240/240 hereditary (all regular)

**R-minimization conjecture (NEW, S18g, REFUTED at n=7 — see OPEN-Q-017):** The H-maximizer was conjectured to minimize R(T) = sum_v H(T-v)/H(T). Formula R(T) = n - E_weighted[|U(S)|] is PROVED. But R-minimization FAILS at n=7: a tournament with H=123 has R=1.585 < R(max)=5/3. Valid only at n=3-6.

**Key insight:** Being hereditary (R = n*H_{n-1}/H_n for regular maximizers) is NOT the same as minimizing R. The non-regular n=5 maximizers have LOWER R (1.4) than regular ones (5/3 ≈ 1.667) despite not being hereditary.
**Next step:** (1) Verify R-minimization at n=7 (running). (2) Prove R-minimization from OCF. (3) Test if regular n=9 maximizers are hereditary.
**Scripts:** `04-computation/hereditary_maximizer.py`, `04-computation/hereditary_correction.py`, `04-computation/R_minimization_proof.py`, `04-computation/R_min_n7_check.py`
### INV-032: Omega(T) is always claw-free AND perfect — Dyer-Jerrum decomposition
**Source:** Web research opus-S5, arXiv:1909.03414 (Dyer-Jerrum-Müller-Vušković)
**Status:** CLAW-FREENESS PROVED for n<=8 (THM-020, vertex counting: claw needs 9+ vertices). PERFECTNESS FAILS at n=8 (THM-019, 53.8% of random n=8 tournaments have C5 in Omega).
**What:** Omega(T) is always claw-free for n<=8 but NOT always perfect for n>=8. The Dyer-Jerrum decomposition framework still applies to claw-free graphs (without perfectness), but the structure is less constrained.
**Critical consequence:** Chudnovsky-Seymour (2007) proved that I(G,x) has ALL REAL ROOTS for claw-free G. Since Omega(T) is claw-free for n<=8, ALL roots of I(Omega(T), x) are negative real (THM-020). At n>=9, claw-freeness fails (90% of random tournaments have a claw in Omega). Real-rootedness may still hold by other mechanisms.
**Next step:** (1) Test claw-freeness at n=9 exhaustively. (2) Test line graph hypothesis (Beineke forbidden subgraphs). (3) Test subdivided-claw-freeness at n=9.

### INV-033: Redei-Berge Hopf algebra formalization of OCF
**Source:** Web research opus-S5, arXiv:2402.07606 (Grinberg), arXiv:2506.08841 (Mitrovic-Stojadinovic)
**Status:** CONNECTION IDENTIFIED. Key bridge found (S36).
**What:** The Redei-Berge symmetric function U_X for digraphs has comultiplication Delta([X]) = sum_S [X|S] tensor [X|V\S] — this IS our subset convolution. The character zeta counts Hamiltonian paths. The antipode S(U_X) = (-1)^|V| U(X-bar) encodes Berge's theorem.
**NEW (Mitrovic-Stojadinovic, arXiv:2506.08841, June 2025):**
  - X_{inc(P)} = omega(U_P): Chromatic function of incomparability graph = omega of Redei-Berge
  - "Converse of Redei": if poset is not a chain, quasi-linear extensions are even
  - Bags-of-sticks decomposition: U_X = sum of simpler digraphs via inclusion-exclusion on edges
  - Stanley-Stembridge connection: e-positivity of X_{inc(P)} <=> h-positivity of U_P
  - Noncommutative deletion-contraction: W_X = W_{X\e} - W_{X/e}^up
  - Mitrovic-Stojadinovic phi(pi) = sum_{gamma X-cycle} (len(gamma)-1) is EXACTLY our S = sum(l_i-1)!
**Verified (S36):** OCF specialization p_1->1, p_{odd>=3}->2, p_{even}->0 gives H(T) from U_T.
**NEW (Mitrovic, arXiv:2504.20968, April 2025):** Noncommutative Redei-Berge function W_X has deletion-contraction: W_X = W_{X\e} - W_{X/e}↑. Thm 3.16: cycle decomposition via inclusion-exclusion over cycle edges. Cor 3.12: tournament formula W_X = Σ(2^{ψ(σ)} p_{Type(σ)}) for odd-cycle permutations = exactly OCF.
**h-POSITIVITY TEST (kind-pasteur-S39b):** h-positivity of U_T FAILS for all non-transitive tournaments. At n=3: 1/2 h-positive, n=4: 1/4, n=5: 1/11. Only the transitive tournament (H=1) is h-positive. The h(2,1) and h(2,2,1) coefficients are always negative for non-transitive. This is expected since tournament posets are NOT (3+1)-free in general, and Stanley-Stembridge conjecture requires (3+1)-freeness.
**Next step:** (1) Express OCF via bags-of-sticks decomposition. (2) Check if deletion-contraction on W_T gives a direct proof of Claim A. (3) Explore chromatic function connection for imperfect Omega(T). (4) Study Thm 3.16 cycle decomposition for odd cycles.

### INV-034: Björklund cycle cover reduction adapted for OCF — TESTED (NEGATIVE for new identities)
**Source:** Web research opus-S5, arXiv:1008.0541, arXiv:1301.7250
**Status:** TESTED (kind-pasteur-S39b). No new identity beyond OCF.
**What:** Björklund reduces Hamiltonian cycle counting to cycle cover counting via inclusion-exclusion and determinants. Tested 6 formulations at n=3-6:
1. Full-vertex all-odd CC weighted: FAILS (0%)
2. Partial odd CC weighted by 2^k: MATCHES H(T) 100% — but this IS OCF
3. Inclusion-exclusion sum (-1)^{n-|S|} perm(T[S]): FAILS (for paths)
4. Irving-Omar odd traces: exploratory only
5. perm(I + x*A): FAILS
6. Odd permanent (unweighted): FAILS (0%)
**Conclusion:** OCF = partial odd cycle cover polynomial at weight 2. This is a restatement, not a new identity. The Björklund approach doesn't give a new route to proving OCF.
**Script:** `04-computation/bjorklund_cycle_cover.py`

### INV-035: Tribonacci structure — OCF for T_full family via interval graphs
**Source:** opus-2026-03-05-S6 (Tribonacci web research), kind-pasteur-S11 (Tribonacci discovery)
**Status:** VERIFIED n=3,...,8. Both sides match Tribonacci(n) = A000213 independently.
**What:** T_full_n (full tiling tournament) has H(T_full_n) = Tribonacci(n) (proved via run decompositions). INDEPENDENTLY, Omega(T_full_n) is an INTERVAL GRAPH on odd-length consecutive intervals [k, k+2j], and I(Omega, 2) satisfies the same Tribonacci recurrence via a weighted interval packing DP that telescopes: f(n) = f(n-1) + 2f(n-3) + 2f(n-5) + ... = f(n-1) + f(n-2) + f(n-3).
**Key structural insight:** All directed odd cycles of T_full_n are consecutive intervals. The clique-cutset decomposition of this interval graph mirrors the DP structure computing H(T_full_n). Both sides produce Tribonacci by the same algebraic mechanism (telescoping) through different combinatorial objects.
**Why this matters:** Shows OCF's "both sides match" emerges from parallel decomposition structures. If this parallelism generalizes (clique-cutset of Omega mirrors Ham path DP), it could prove OCF.
**Extended results (opus-S13):**
- **Transitive+flip(i,j):** H = 1 + 2^(j-i-1). All odd cycles form a clique in Omega, so I(Omega,2) = 1 + 2·(#cycles) = 1 + 2^(j-i-1). Clean OCF-based proof.
- **Cone theorem:** H(source_cone(T')) = H(sink_cone(T')) = H(T') for ALL T'. Proved: source must be first in every Ham path. Verified exhaustively through n'=6.
- **Partial cones palindromic:** H(k) = H(n'-k) where k = out-degree of cone vertex. From self-converse symmetry.
- **Circulant S={1}:** H(T_{n,{1}}) = n, order-2 recurrence. No circulant with |S|>=2 has low-order recurrence.
- **CORRECTED (S41):** Circulant S={1,5,6,7} at n=9 DOES give H=3357 = max. Non-circulant maximizers also exist.
**Next step:** (1) Find direct bijection between run decompositions and weighted interval packings. (2) Check if the transfer matrix for T_full has Tribonacci characteristic polynomial factor.

### INV-001: Prove transfer matrix symmetry for all n — PROVED (THM-030)
**Source:** T214, T103 (tangents), symmetry_check.py, symbolic_symmetry_proof.py
**Status:** PROVED FOR ALL n by induction (kind-pasteur-2026-03-06-S25). Verified computationally m=2,...,6 all (a,b) pairs.
**What:** M[a,b] = sum_S (-1)^|S| E_a(S)*B_b(M\S) is always symmetric. This is STRONGER than the even-odd split.
**BREAKTHROUGH (opus-S4):** M[a,b]-M[b,a] = 0 as a polynomial in the arc variables t_{ij} AFTER applying the tournament constraint T[j,i]=1-T[i,j]. With independent arc variables the difference is NONZERO (12 terms at n=4, 48 at n=5). The tournament constraint is essential and sufficient.
**Equivalent formulation:** M_{T^op} = (-1)^{n-2} M_T (converse identity). Combined with path reversal M_{T^op}[i,j]=(-1)^{n-2}M_T[j,i], gives symmetry.
**Key insight:** Connects to Feng's Dual Burnside (Q=AB symmetric under detailed balance). The tournament constraint T[x,y]+T[y,x]=1 plays the role of the "detailed balance" condition.
**c-TOURNAMENT GENERALIZATION (opus-S19b):** Symmetry holds for ALL c-tournaments where t_ij + t_ji = c (any constant c, not just c=1). Verified symbolically n=3,4,5; numerically n=6,7. The constraint must be UNIFORM across pairs (non-uniform c_ij gives 100% failure) and ALL pairs need it. In skew coordinates t_ij = c/2 + s_ij: M is EVEN in s for n even, ODD for n odd. The c^{n-2} coefficient is (n-2)!/2^{n-2} for even n, 0 for odd n. This c-generalization SIMPLIFIES the proof problem: we can work with skew-symmetric part S = A-A^T and ignore the specific value c=1.
**New findings (opus-S6):**
- **THM-027 PROVED:** Trace formula tr(M) = H(T) for odd n, 0 for even n. Clean bijection proof via (-1)^{pos(a,P)} formula for diagonal entries.
- **MISTAKE-011:** Old claim M = [[1,0],[0,-1]] always is FALSE (2199/2500 failures at n=4). M entries range from -3 to +3.
- **Off-diagonal sum:** sum_{a≠b} M[a,b] = 0 (odd n), 2*H(T) (even n). Verified n=3,...,7 but NOT yet proved.
- **Complement pairing D(S)+D(U\S) is constant at n=4 but NOT at n=5**, ruling out the simplest telescoping argument.
- **Cauchy-Binet decomposition:** M = E^T * Lambda * B where E[S,v]=E_v(S), B[S,v]=B_v(U\S), Lambda=diag((-1)^|S|). Symmetry equivalent to E^T*Lambda*B = B^T*Lambda*E.
**PATH REVERSAL PROOF AT c=0 (kind-pasteur-S23):** COMPLETE proof when c=0 (pure skew weights). Path reversal: B_v(S+v) = (-1)^|S| E_v(S+v). This gives M[a,b] = (-1)^{n-2} sum_S E_a(S+a) E_b(R+b) — unsigned, manifestly symmetric by S<->R relabeling. Verified n=3,4,5,6.
**EVEN r-POWERS CONJECTURE (kind-pasteur-S23):** At general c, M(r,s) where r=c/2 has ONLY even r-powers. Equivalent to symmetry. Verified n=3,4,5,6. Path reversal gives B_v(c,s) = E_v(c,-s), which yields M[b,a](c,s) = (-1)^{n-2} M[a,b](c,-s). So symmetry reduces to M having definite s-parity (-1)^{n-2}, i.e., only even r-powers.
**ALGEBRAIC PROOF (kind-pasteur-S23b): M[a,b](-r) = M[b,a](r)** — 5-step proof: T(-r)=-T^T, path reversal under negated transpose, sign bookkeeping, S↔R relabeling. Verified n=4,5,6. This proves the EQUIVALENCE between (i) symmetry, (ii) even-r-powers, (iii) s-parity.
**TOGGLE ANALYSIS (S23b):** At n=4, r^1 monomials cancel pairwise between different S-subsets. At n>=5, cancellation is multi-way (not simple pairwise). No clean single-vertex toggle involution found on whole subsets.
**H(U) MATRIX (S23b, from Kogan/Hamiltonian cycle polynomial):** H(U)_{i,j} = sum of Ham path weights from i to j. Identity: U*H(U)^T = H(U)*U^T. For c-tournaments U+U^T = c(J-I), this gives UH^T = H(cJ-cI-U), but does NOT directly imply H=H^T. Also note: M[a,b] is NOT the same as H(T)_{a,b} (M has inclusion-exclusion signs, H is a direct sum).
**PROVED INDEPENDENTLY by both opus-S25 and kind-pasteur-S25 (THM-030).** Key Identity: odd_r(B_b(W)) = r * col_sum_W(b), equivalently B_b(W)+(-1)^m E_b(W) = 2r*col_sum(b). Inductive proof using column recurrence + first-edge decomposition. The Sigma identity (r*Sigma = odd(T)) follows from summing the inductive hypothesis. The proof closes because odd(sum s*even(B_v)) = 0. Verified computationally m=2..6 all (a,b) pairs. See complete_even_r_proof.py, key_identity_complete_proof.py.
**Scripts:** `04-computation/symbolic_symmetry_proof.py`, `04-computation/transfer_symmetry_analysis.py`
**COEFFICIENT STRUCTURE (opus-S22 continuation):**
- **[r^{n-2}] = (n-2)!** when n even, **0** when n odd. Proof: counting argument, sum_k C(n-2,k)(-1)^k k!(n-2-k)! = (n-2)! * sum_k(-1)^k. Verified n=3,...,6.
- **[r^2] for n=5 = 2·sum_{u∈U}(s_{au}+s_{bu})**. For n=6: degree-2 in s with all coefficients ±2. For n=4: just 2 (constant).
- **[r^1] telescoping (n=4):** Each s_{uv} (u∈U, v∈{a,b}) appears exactly once with + and once with - across subsets. Moving vertex u between S and R flips the sign contribution.
- **M is NOT a cofactor** of A=rJ'+S (exhaustive test n=3,4). Cofactors have degree n-1; M has degree n-2.
- **M is NOT a permanent minor** of A either. The fundamental identity A(-r)=-A^T is clean but M does not decompose as a simple matrix function of A.
- **Key algebraic identity:** A(-r,s) = -A(r,s)^T (since J' symmetric, S skew). Any expression built from A·A^T+A^T·A (which is even in r) could explain the property, but no such expression matching M has been found.
- **Literature update:** Irving-Omar (arXiv:2412.10572), Mitrovic noncommuting Redei-Berge (arXiv:2504.20968) with deletion-contraction W_X = W_{X\e} - W_{X/e}^up, El Sahili-Ghazo Hanna proving T and T^op have same Hamiltonian path type distribution.
**APPROACH RULING (opus-S22 continuation):**
- ❌ Simple cofactor/minor of A (degree mismatch)
- ❌ Permanent of A minor (doesn't match)
- ❌ Adjugate entries of A, I±A, J-A (all fail)
- ❓ Deletion-contraction via Mitrovic noncommuting Redei-Berge (unexplored, most promising NEW lead)
- ❓ Irving-Omar walk generating function det(I+zXĀ)/det(I-zXA) (connection to M unclear)
- ❓ Direct r^1=0 proof via telescoping + induction on n (promising for base case)
**Next step:** (1) Try Mitrovic deletion-contraction approach — express M[a,b] recursively and prove even-r by induction. (2) Understand Irving-Omar matrix formula and whether it encodes M[a,b]. (3) Prove [r^1]=0 directly via the telescoping structure observed at n=4,5. (4) Previous approaches (Feng, Hopf, involution) remain viable but untested.
**Scripts:** `04-computation/symbolic_symmetry_proof.py`, `04-computation/transfer_symmetry_analysis.py`, `04-computation/determinantal_identity_test.py`, `04-computation/det_compare_explicit.py`, `04-computation/r1_coefficient_analysis.py`, `04-computation/r_coefficient_structure.py`

### INV-002: Subset convolution identity — the core algebraic challenge
**Source:** proof-landscape-for-general-ocf.md (Approach B), T047
**Status:** Correct framework identified. No simplification found.
**What:** sum_S [f_i(S)*g_j(R) - f_j(S)*g_i(R)] = sum_{k>=1} 2^k Delta(alpha_k). Both sides are multilinear polynomials. The bracket B(u,w) has a 4-way type structure where Z0 rows and Z1 columns vanish.
**Next step:** Use the bracket table to decompose the convolution into 6 nonzero bracket types. Try to show the resulting expression telescopes via induction on the number of M-/M+ vertices.

### INV-003: Sign-reversing involution on the subset convolution
**Source:** proof-landscape-for-general-ocf.md (bottom), signed-adjacency-identity.md
**Status:** Idea only. Not attempted.
**What:** Find an involution on the terms of sum_S [f_i(S)*g_j(R) - f_j(S)*g_i(R)] that cancels everything except the cycle terms. This is the "bijective" approach to the algebraic identity.
**Connection:** The sigma-invariance (s -> -s, B is even in s-variables) reduces OCF to proving all s-degree-1 terms vanish: C_w + D_w = 0 for each w. This is a per-vertex condition that might be provable.

### INV-004: Flip-class proof strategy (prove for R-cones, extend via cut-flip)
**Source:** T046, paper-connections.md (CONNECTION 1, "Flip Class + OCF")
**Status:** Strategy identified. Not attempted.
**What:** Rajkumar et al. show every tournament is in the flip class of an R-cone. For R-cones (vertex beating/losing to everyone), Ham paths all start or end at universal vertex, simplifying both H(T) and I(Omega,2). Prove OCF for R-cones, then show cut-flip phi_S preserves E(T) = H(T) - I(Omega(T),2) = 0.
**Key gap:** Need to track how both H and I change under phi_S (reversing all arcs across a cut). This is a MULTI-arc flip, not a single arc flip.

### INV-005: Induction on mu(T) (flip-feedback dimension)
**Source:** paper-connections.md (CONNECTION 4)
**Status:** Idea only. Not formalized.
**What:** Rajkumar et al.'s mu(T) measures minimum flip-feedback node set. Dim bound: 2(mu(T)+1). For mu=0, T is transitive (OCF trivial). Could induct on mu(T), with each step being a cut-flip. Need: does cut-flip increase or preserve mu?

### INV-006: n=8 exhaustive proof completion
**Source:** OPEN-Q-009, ocf_n8_full.c
**Status:** PROVED by opus-S4 (2^27 configs, 57min, all passing). Independent verification by opus-S4b C implementation (3M+ configs, 0 fails through partial run).
**Next:** Close this out. Focus on n=9 strategy or general proof.

### INV-053: Even Cycle Vanishing Theorem — PROVED
**Source:** opus-2026-03-06-S10, T148
**Status:** PROVED. Clean involution argument.
**What:** For any tournament T on [n], p_mu(U_T) = 0 whenever mu has an even part. The proof pairs each permutation sigma with even k-cycle c with sigma' (c reversed); the sign flips because (-1)^{k-1} = -1 for even k. Verified computationally n=3 through n=7.
**Consequences:** (1) U_T lives in the subspace spanned by p_mu with all odd parts — drastically fewer terms. (2) At n=4, only types (1^4) and (3,1) contribute; at n=5, only (1^5), (3,1,1), (5); at n=6, only (1^6), (3,1,1,1), (3,3), (5,1). (3) The Schur expansion simplifies: [s_lambda]U_T = sum over odd-part mu only. (4) This is the SAME T<->T^op involution as in the path reversal proof (T147).
**Connection to INV-001:** The even-r-powers conjecture (kind-pasteur-S23) is the transfer matrix version of this same phenomenon. Both arise from the perpendicular grid symmetry.

### INV-054: Hook Schur Positivity for Tournaments — PARTIAL (fails at n=7)
**Source:** opus-2026-03-06-S10, T149
**Status:** PROVED at n=4 (clean sign argument). VERIFIED at n=5 (11/11), n=6 (40/40). **FAILS at n=7** (231/242, 11 failures all for middle hook (4,1,1,1)).
**What:** [s_{(k,1^{n-k})}]U_T >= 0 holds at n=4,5,6 but fails at n=7.

**n=4 proof:** Only p-types (1^4) and (3,1) matter (by INV-053). All hook characters non-negative at both → sum of non-negative terms.

**n=7 failure mechanism:** chi^{(4,1,1,1)}((7)) = -1 and chi^{(4,1,1,1)}((5,1,1)) = 0. Regular tournament T_7 has 48 directed 7-cycles contributing -48/7, overwhelming positive 3-cycle terms. Result: [s_{(4,1,1,1)}]U_{T_7} = -83/28 ≈ -2.96.

**Which hooks always hold:** Hooks (n) and (1^n) have all-positive characters at odd types → always positive. Hooks (n-2,1,1) and (3,1^{n-3}) also have all-positive chars (at n=7). Only hooks with j odd and j near n/2 can fail.

**Non-hook negativity:** Non-hook chars at (3,1,...,1) are always negative → non-hook coefficients negative for all non-transitive tournaments (verified n=4,5,6).

**Refined question:** For which hooks is positivity universal? Is there a simple characterization (e.g., j even, or |j - n/2| > threshold)?
**Scripts:** `04-computation/schur_hook_analysis.py`, `04-computation/tournament_cycle_structure.py`, `04-computation/even_cycle_vanishing_proof.py`, `04-computation/hook_positivity_n6.py`, `04-computation/hook_positivity_n7.py`

---

### INV-200: Theorem number collisions — 10 numbers have multiple files (HOUSEKEEPING)
**Source:** opus-2026-04-01-S1 (systematic audit)
**Status:** DOCUMENTED. Needs human decision on which file keeps the original number.
**Collisions (22 files):**
- THM-130: `c5-paley-closed-form` vs `paley-betti-formula`
- THM-133: `h-trace-formula-p7` vs `spectral-ocf-chain`
- THM-134: `paley-local-max-parseval` vs `schur-concavity-dichotomy`
- THM-137: `paley-crossover-mechanism` vs `paley-eigenvector-interaction`
- THM-201: `grand-energy-formula` (CORRECTED) vs `k3-component-impossibility`
- THM-217: `binomial-truncation` vs `transfer-matrix-gk`
- THM-224: `golden-exceptional-points` vs `simplicial-up-laplacian`
- THM-260: `bandlimited-all-n` vs `edge-count-structure` vs `rooted-tournament-layer-decomposition` (THREE files!)
- THM-261: `petersen-root-orthogonality` vs `sc-orbit-pairing`
- THM-262: `dual-lie-embedding` vs `seesaw-identity`
**Next step:** Human should decide which file keeps the number. Duplicates should be renumbered to THM-266+. Check for references before renaming.

---

## Priority B: Important structural understanding

### INV-007: Odd-cycle bijection (Open Problem 3 in paper)
**Source:** oq:bijection in tex, bijection_search.py, T046
**Status:** Searched computationally at n=3,4. No natural bijection found.
**What:** Construct Phi: Ham(T) -> {(C, f) : C vertex-disjoint odd-cycle collection, f: C -> {0,1}}. The "2-colored cycles" interpretation means each independent set of k cycles contributes 2^k paths. At n=3 with one 3-cycle: 3 paths = 1 + 2. Need to identify which paths correspond to which colored cycle sets.
**Key obstacle:** The correspondence is NOT local/contiguous-block based (T035, confirmed dead end). Must be global.

### INV-008: Striker-Chapman S3-equivariance (Open Problem 5 in paper)
**Source:** oq:striker in tex, \cite{striker2011,chapman2001}
**Status:** INVESTIGATED (kind-pasteur-S22 agent). Question imprecise as stated.
**What:** Is the Striker-Chapman bijection between ASMs and tournaments S3-equivariant under the barycentric identification? Striker (2011) gives a unifying poset perspective connecting ASMs, plane partitions, Catalan objects, tournaments, and tableaux. Chapman (2001) connects ASMs to tournaments directly.
**Finding:** Chapman's bijection maps oriented monotone triangles (not ASMs directly) to tournaments. Striker's tournament bijection uses the disjoint three-color poset {b,r,(g)}, while ASMs use the four-color poset {b,y,o,g} — different subposets. No S3 action is defined in either paper. The question needs precise formulation: (a) define "the bijection" (Chapman's Phi? Striker's poset? composition?), (b) define how S3 acts on both sides. Only then can it be tested at n=3,4.
**Why it matters:** If S3-equivariant, then the S3 orbit counts in our Section 6 could be computed via ASM symmetries, potentially giving new structural constraints on H(T).
**Next step:** Read Striker (2011) and Chapman (2001) carefully. Test equivariance computationally at small n.

### INV-009: Self-evacuating SYT bijection (Open Problem 6 in paper)
**Source:** oq:se_bijection in tex, Section 8 (tetrahedral geometry)
**Status:** Count verified (2^{m^2} for n=5,7). Bijection NOT constructed.
**What:** Natural bijection between 2^{m^2} self-evacuating SYT of delta_{n-2} and 2^{m^2} sigma-fixed tilings of Grid(2m+1). Both count the same thing but via very different combinatorial objects.
**Connection:** TSSCPPs (Totally Symmetric Self-Complementary Plane Partitions) and the ASM conjecture. The TSSCPP count for order m is known to equal the ASM count.

### INV-010: Mixed graphs extension (Open Problem 4 in paper)
**Source:** oq:mixed in tex, \cite{schweser2025}
**Status:** NOT INVESTIGATED.
**What:** Extend the Q-Lemma to complete mixed graphs, recovering Schweser-Stiebitz-Toft (2025) strengthening of Redei.
**Next step:** Read arXiv:2510.10659 carefully. Understand what "complete mixed graph" means and how the Q-Lemma generalizes.

### INV-011: Mod-4 score-sequence criterion (Open Problem 9 in paper) — RESOLVED (NO)
**Source:** oq:mod4_struct in tex
**Status:** RESOLVED NEGATIVELY (kind-pasteur-2026-03-07-S39b).
**What:** Does score sequence determine H(T) mod 4?
**Answer:** NO for n >= 5. The formula H mod 4 = (1 + 2·c3) mod 4 holds only for n <= 4.
**Results:**
- n=3,4: YES (exhaustive, H mod 4 constant within every score class)
- n=5: NO — score (1,2,2,2,3) has H mod 4 in {1,3} (c5 varies within this class)
- n=6: NO — 5/22 score classes have varying H mod 4
- n=7: NO — 27/59 sampled score classes have varying H mod 4
**Key insight:** H mod 4 = (1 + 2·alpha_0) mod 4 where alpha_0 = total odd cycle count.
For n <= 4, alpha_0 = c3 which is score-determined. For n >= 5, 5-cycles contribute to alpha_0 but c5 is NOT determined by score sequence (confirmed independently).
**Also found:** c5 is NOT determined by (score, sum_d², edge_score_sum). Even (score, common_out_neighbor) pairs vary. c5 requires genuine graph structure beyond all local/pairwise statistics.
**Scripts:** `04-computation/mod4_score_test.py`, `04-computation/c5_score_determination.py`

### INV-012: BlackSelf(8) exceptional class (Open Problem 7 in paper) — PARTIALLY RESOLVED
**Source:** oq:n8 in tex
**Status:** DEEP INVESTIGATION by opus-2026-03-05-S8. See `03-artifacts/drafts/n8-anomaly-deep-dive.md`.
**What:** Unique isomorphism class at n=8 that is self-converse, has |Aut|>1, |Fix(beta)| odd, but H(T)/|Fix(beta)| is even. Is it related to a Hadamard matrix of order 8 or a skew conference matrix?
**Key findings:**
- Exhaustive search over ALL 65536 SC tournaments (alpha = reversal): ZERO have Fix(beta) | H with H/Fix even.
- The definition likely means (H-Fix)/2 is even (number of beta-orbit pairs is even).
- Under this interpretation: T_657 (H=657, Fix=33, |Aut|=3) is the best candidate: (657-33)/2 = 312 (even).
- T_657 CONTAINS P(7) (Paley tournament on 7 vertices) as vertex-deletion. P(7) ↔ unique skew Hadamard matrix of order 8. THIS IS THE HADAMARD CONNECTION.
- T_657 has perfectly uniform D_v = 54 (mu-weighted 3-cycle count) for all 8 vertices.
- Full survey: 10 distinct (H, Fix, Aut) combinations among 2560 SC+Aut>1 tournaments.
**Next step:** Confirm T_657 is isomorphic to a known Paley extension. Resolve definition ambiguity with paper author.

### INV-013: Realizable odd-cycle conflict graphs (Open Problem 8 in paper) — INVESTIGATED
**Source:** oq:realizable in tex
**Status:** INVESTIGATED (opus-S13 background agent). Key structural findings.
**What:** Which graphs G arise as Omega_3(T) for some tournament T?
**Results:**
- Realizable isomorphism classes: 2 (n=3), 3 (n=4), 6 (n=5), 18 (n=6), ~97+ (n=7 sampled)
- n≤5: Omega_3 is ALWAYS a complete graph (two 3-cycles on ≤5 vertices must share a vertex)
- n=6: Always "complete minus matching" — complement is disjoint union of edges
- n=7: First non-perfect graphs (11/97 classes have chi=omega+1)
- alpha(Omega_3) ≤ floor(n/3) — proved by vertex counting (3k vertices for k disjoint 3-cycles)
- 100% real-rooted across all realizable classes (exhaustive n≤6, sampled n=7)
**Key insight:** The low independence number (alpha ≤ 2 for n≤8) means I(Omega_3, x) has degree ≤ 2, so real-rootedness is "easy" for Omega_3. The full Omega (including 5,7-cycles) has alpha ≤ floor(n/3) similarly, keeping degree low.
**Next step:** Characterize which "complete minus matching" graphs at n=6 are NOT realizable. Extend to n=8 where alpha=2 still but degree may increase.

### INV-014: 2-adic tower / higher Redei theorems — PARTIALLY RESOLVED
**Source:** OPEN-Q-008, T007, tex Section 5.5
**Status:** COMPUTED (opus-S13). v_2(H(T)) = 0 ALWAYS (= Redei's theorem).
**What:** I(Omega(T), x) at x=4,8,... gives mod-4, mod-8 invariants of H(T). v_2(H(T)) = 0 universally.
**Results:** H mod 4 ≡ 1+2*alpha_1 (mod 4) via OCF. At n=3,4 this equals 1+2*c3 (mod 4) exactly. At n≥5 the c3 formula breaks (5-cycles contribute to alpha_1). H mod 2^k approaches uniform on odd residues as n grows.
**Impact:** OPEN-Q-008 partially resolved. No deeper 2-adic structure at level of H(T). The mod-4 structure is fully explained by alpha_1 parity via OCF.

### INV-136: Dimensional Meta-Patterns (Tournament = Simplex Orientation)
**Source:** kind-pasteur-2026-03-09-S44
**Status:** IN PROGRESS. Core data collected n=3-10.
**What:** Framing tournaments as binary relations on simplices (T_n = orientation of Delta_{n-1}). Systematic study of how homological invariants scale with dimension d=n-1.
**Key findings:**
1. Transitive tournament = contractible simplex: dim(Omega_p) = C(n,p+1) (HYP-302)
2. Filling ratio f_p = dim(Omega_p)/C(n,p+1) inflates at high p for n>=6 (HYP-303)
3. H(T_4) = 2*c3+1 for ALL 4-vertex tournaments (HYP-304)
4. excess_4 = 2*c3*(n-3) universally (HYP-305)
5. P(beta_1>0) decays exponentially ~exp(-0.755n) (HYP-310)
6. P(beta_3>0) grows: 0->0.4->7.2->19->23% for n=5-9, may saturate ~25%
7. Beta_5 NOT YET observed at n=7,8,9 — onset unknown
8. Chi(T) in {0,1} for n<=7; chi up to 6 at n=8 (HYP-312)
9. dim(Omega_2) NOT determined by (c3, score) (HYP-308/309)
10. |A_p| mod 2 = C(n,p+1) mod 2 via local Redei
11. Poincare polynomial P(T,1)/(2^n-1) grows with n — path complex exceeds simplex
12. Surplus = excess_paths - rank(constraints) exactly (HYP-314)
**Scripts:** dimensional_crossover.py, filling_ratio_formula.py, local_redei_investigation.py, euler_char_scaling.py, chi_A_identity.py, omega2_formula.py, betti_rate_scaling.py, poincare_polynomial.py, omega_parity_structure.py, beta5_onset_search.py
**Next steps:** ~~(1) Find beta_5 onset~~ FOUND: beta_5 at n=8. ~~(3) Prove beta_1*beta_3=0~~ PROVED: THM-095 seesaw. (2) Formula for filling ratio. (4) Explain defect rate U-shape. (5) PROVE beta_2=0 algebraically (critical). (6) Investigate n=9 exotic profiles.

### INV-137: Seesaw Mechanism and Tournament Homology Structure
**Source:** kind-pasteur-2026-03-09-S45
**Status:** MAJOR RESULTS. THM-095 and THM-096 established.
**Key findings:**
1. THM-095: beta_1*beta_3=0 via seesaw through beta_2=0. im(d_2) mediates: 2 values only.
2. beta_2=0 for ALL tournaments (0/1000 at n=8). DEEPEST structural invariant.
3. beta_4>0 at n=8 (~1.1%), values {1, 5}. Even Betti CAN appear, just not beta_2.
4. beta_5 onset at n=8 (with beta_1=1). Profile [1,1,0,0,0,1,0,0], chi=-1.
5. beta_3+beta_4 coexist at n=8 (~0.15%), chi=1.
6. chi ranges over {-1, 0, 1, 2, 6} at n=8.
7. Constraint ratio NA_faces/|A_p| < 1 for p=2 always (may explain beta_2=0).
8. THM-096 corrected: simplicity holds for n<=7 only.
**S45 continuation updates:**
9. THM-119 (was THM-097) (PROVED): Disjoint support at Omega_2 — each 2-path has at most 1 non-allowed face. Constraint matrix always full rank. dim(Omega_2) = |A_2| - #NA_faces exactly.
10. Completeness is SHARP: removing 1 edge from tournament creates beta_2>0 (13/500 at n=6).
11. beta_2=0 confirmed at n=9 (0/500) and n=10 (0/100). Disjoint support verified at n=9.
12. H_1(A-complex) is NOT always 0 (only for transitive tournaments). Omega restriction is essential.
13. rank(d_3) is NOT a simple function of c3 (multiple values per c3 class).
14. 3-path NA face distribution: exactly 25%/50%/25% for 0/1/2 faces (very clean, universal).
15. At level 4: paths can have multiple NA faces => overlapping constraint rows => rank deficit => beta_4 possible.
16. Complete bidirectional graph: beta_2=0, homology at top dimension only.
**S45 defect rate analysis (omega2_exact_formula.py, defect_rate_ushape.py, ecyc_formula.py):**
17. dim(Omega_2) = C(n,3) + 2*c3 - e_cyc EXACTLY (exhaustive n=4,5,6). e_cyc = #{directed edges in ≥1 three-cycle}.
18. e_cyc NOT determined by c3 alone — depends on cycle arrangement (edge sharing). Constant for most score seqs, varies near-regular.
19. Defect rate is WAVE PROPAGATION, not U-shape: beta_1 rate decreasing (29.7%→1%), beta_3 increasing (0%→21%), beta_4 appears at n=8 (2%).
20. beta_4 CAN be nonzero (values 1,2 at n=8) — "all even beta vanish" is FALSE. Only beta_2=0 always.
21. Only 3-5 distinct Betti profiles per n. Extremely constrained.
22. beta_3*beta_5=0 at n<=8 (trivially since beta_5 rarely nonzero). Generalized seesaw needs testing at n>=10.
**Open:** (1) Beta_2=0 is PROVED (THM-108+109). (2) Does generalized seesaw beta_{2k-1}*beta_{2k+1}=0 hold generally? (3) Onset of beta_6 (not seen at n<=8 in 600+ samples). (4) Why exactly 3-5 Betti profiles per n?
**Scripts:** beta1_beta3_mediator.py, even_betti_quick_v2.py, beta4_investigation.py, beta_coexistence_analysis.py, beta2_algebraic_analysis.py, beta2_disjoint_support_proof.py, beta2_completeness_argument.py, beta2_exactness_proof.py, beta2_only_n9.py, omega2_exact_formula.py, defect_rate_ushape.py, ecyc_formula.py

### INV-138: Beta_3 ≤ 1 Proof Architecture (LES Induction)
**Source:** kind-pasteur-2026-03-09-S46
**Status:** PROOF ARCHITECTURE COMPLETE. Both algebraic ingredients computationally verified. No algebraic proof yet.
**What:** THM-123 (was THM-110): beta_3(T) ≤ 1 for all tournaments T. Equivalent to rank near-saturation: rank(d_4) ≥ ker(d_3) - 1.
**Proof strategy:** LES induction on n using pair (T, T\v):
  ... → H_3(T\v) → H_3(T) → H_3(T,T\v) → H_2(T\v) = 0
  Since H_2(T\v) = 0 (THM-108), map H_3(T) → H_3(T,T\v) is surjective.
  Find v with beta_3(T\v) = 0. Then beta_3(T) = dim H_3(T,T\v).
**Key ingredients verified:**
1. Good vertex existence for beta_3: ∃v with beta_3(T\v)=0 when beta_3(T)>0.
   - n=6: 320/320 exhaustive. beta_3 COMPLETELY fragile (ALL 6 deletions give 0).
   - n=7: 34/34 sampled. 5-7 good vertices per tournament.
   - n=8: 31/31 sampled.
2. Relative H_3 bound: dim H_3(T,T\v) ≤ 1 ALWAYS.
   - n=6: ALL 1920 pairs give dim=1 (exhaustive for beta_3>0).
   - n=7: dim ∈ {0,1}. Max=1.
3. LES isomorphism: beta_3(T\v)=0 ⟹ beta_3(T) = dim H_3(T,T\v). Perfect n=6 (1920/1920).
**Additional findings:**
4. Quotient proportionality: ALL ker(d_3) basis vectors project proportionally to H_3 (240/240 Type B at n=6).
5. Cokernel direction varies by tournament (NOT universal).
6. Two beta_3=1 types at n=6: Type A (scores 1,1,1,4,4,4, Omega_4=0, 80 tours) and Type B (2,2,2,3,3,3, ker_d3=7, 240 tours).
7. H_3 generator: Type A uses 9 paths/9 vertex sets; Type B uses 36 paths/all 15=C(6,4) vertex sets.
8. At n=7 (2000 samples): max beta_3 = 1. ker(d_3) ranges 10-46. When beta_3=1, rank(d_4) = ker(d_3)-1 always.
9. Relative complex dims at n=6: two H_3=1 profiles. Type A: (d2,d3,d4)=(9,6,0). Type B specific: (12,14,8).
10. Filling ratio f_2 nearly linear in c3 (from 1.0 to 1.08 at n=6). Higher f_p grow rapidly.
**Scripts:** rank_near_saturation.py, beta3_homology_structure.py, beta3_les_analysis.py, beta3_good_vertex_and_relative_h3.py, beta3_proportionality_proof.py, relative_h3_structure.py, defect_ushape_filling_ratio.py
**Next steps:** (1) PROVE good vertex existence algebraically (key open). (2) PROVE relative H_3 bound algebraically. (3) Investigate whether quotient proportionality can be proved directly. (4) Extend LES approach to beta_5.
**NOTE:** HYP-342 (Boolean odd Betti) needs correction: TRUE for k=1,2 (beta_1,beta_3 ∈ {0,1}), but FALSE for k≥3 (beta_5=10 at n=9 Paley maximizer). The "Boolean" property is specific to beta_1 and beta_3.

---

## Priority C: References to investigate

### INV-015: Rajkumar et al. (arXiv:2110.05188) — tournament representations
**Source:** paper-connections.md, T046, paper-deep-connections.md Section 2
**Status:** INVESTIGATED (opus-S4b, opus-S5). Key theorems extracted. Connections documented.
**What:** Flip classes, locally transitive = rank 2, R-cones, sign rank. Key results: every T in flip class of R-cone (Prop 1, distance exactly 1 cut-flip); mu(T) dimension bound <= 2(mu(T)+1) (Thm 11); sign-rank bound (Thm 12).
**Key finding:** Proposition 1 is constructive: T' = phi_{i ∪ T_i^-}(T) is R-coned by i. Cut-flip distance to R-cone is exactly 1. This directly enables INV-004 strategy.
**Assessment:** mu(T) induction (INV-005) less promising than FAS induction — mu may change unpredictably under cut-flips.
**Action needed:** ADD to bibliography. Concretely develop INV-004 (R-cone + cut-flip proof).
**Tested:** Locally transitive tournaments DO have 5/7-cycles (T046). OCF passes 100% for LT, R-cones, automorphism-symmetric tournaments.

### INV-016: Feng (arXiv:2510.25202) — dual Burnside process
**Source:** paper-connections.md, paper-deep-connections.md Section 3, tex line 1600/2151
**Status:** INVESTIGATED (opus-S4b, opus-S5). Key theorems extracted. Deep connection found.
**What:** Q=AB factorization, primal-dual spectral correspondence, lumping theory.
**Key findings:** (1) Q=AB is REVERSIBLE with pi(g)=|X_g|/(|G|*z) — detailed balance gives symmetry (Thm 3.3). (2) Block-flip M=[[0,A],[B,0]], M^2=[[Q,0],[0,K]] — bipartite structure with period 2. (3) Eigenvector intertwining (Thm 3.10): A maps K-eigenvectors to Q-eigenvectors, B maps back. (4) Our transfer matrix has EXACTLY this AB structure: A maps subsets to "path ends at vertex", B maps "path starts at vertex" to complement subsets. Transfer matrix symmetry (INV-001) = hidden detailed balance condition.
**Action needed:** Try to formalize the "hidden detailed balance" — identify the group action and show it satisfies Feng's reversibility conditions. This could prove INV-001.

### INV-017: El Sahili & Abi Aad (2020) — parity of paths in tournaments
**Source:** tex bibliography, \cite{elsahili2020}
**Status:** Referenced in tex for decisive/concordant classification. NOT deeply investigated for connections to OCF.
**What:** Discrete Math 343 (2020), Art. 111695. Mod-4 congruences.
**Action needed:** Read the paper. Check if their mod-4 results constrain or relate to our alpha_1 = |C_3| (mod 2) conjecture (INV-011).

### INV-018: El Sahili & Ghazo Hanna (2023) — number of Ham paths/cycles
**Source:** tex bibliography, \cite{elsahili2023}
**Status:** Referenced. NOT investigated for OCF connections.
**What:** J. Graph Theory 102 (2023), 684-701. About the number of oriented Hamiltonian paths and cycles in tournaments.
**Action needed:** Read the paper. They study H(T) directly — any bounds or structural results could inform OCF.

### INV-019: Schweser-Stiebitz-Toft (arXiv:2510.10659) — Redei revisited
**Source:** tex bibliography, \cite{schweser2025}, paper-deep-connections.md Section 1
**Status:** INVESTIGATED (opus-S5). Key theorems extracted.
**What:** Redei's Stronger Theorem (Thm 1.1): add non-oriented edges to tournament, #Ham paths beginning AND ending in tournament vertices is EVEN. Berge's Stronger Theorem (Thm 1.2): G and G-bar have same Ham path parity. Dirac's Stronger Theorem (Thm 2.1): inclusion-exclusion on edge subsets.
**Key findings:** (1) Direct connection to Open Problem 4 (mixed graphs). Non-oriented edges DOUBLE insertion opportunities in Q-Lemma. (2) Berge's theorem gives H(T) ≡ H(T^op) (mod 2) for tournaments, and constrains I(Omega(G),2) under complementation. (3) Strategy for extending Q-Lemma: verify computationally that inshat remains odd for mixed graphs.
**Action needed:** Test inshat parity for mixed graphs computationally. If confirmed, Q-Lemma proof extends directly.

### INV-020: Striker (2011) — unifying poset perspective
**Source:** tex bibliography, \cite{striker2011}
**Status:** Referenced in Open Problem 5. NOT investigated.
**What:** Adv. Appl. Math. 46 (2011), 583-609. Connects ASMs, plane partitions, Catalan objects, tournaments, tableaux via posets.
**Action needed:** Read the paper. Check S3-equivariance (INV-008). The poset perspective may give new structural insights for OCF.

### INV-021: Chapman (2001) — alternating sign matrices and tournaments
**Source:** tex bibliography, \cite{chapman2001}
**Status:** Referenced in Open Problem 5. NOT investigated.
**What:** Adv. Appl. Math. 27 (2001), 290-298. Direct connection between ASMs and tournaments.
**Action needed:** Read the paper. The ASM connection could provide algebraic tools (determinantal formulas, etc.) for H(T).

### INV-022: Eplett (1979) — self-converse tournaments
**Source:** tex bibliography, \cite{eplett1979}
**Status:** Referenced briefly. NOT investigated for OCF.
**What:** Canad. Math. Bull. 22 (1979), 23-27. Self-converse tournament counts.
**Action needed:** Check if self-converse tournaments have special OCF properties. The BlackSelf(8) class (INV-012) is self-converse.

### INV-023: Forcade (1973) — parity of paths and circuits
**Source:** tex bibliography, \cite{forcade1973}
**Status:** Referenced heavily. Our paper gives a NEW combinatorial proof of his F2-invariance.
**What:** Discrete Math 6 (1973), 115-118. Original F2-invariance proof via generating functions.
**Action needed:** Compare his generating function approach to our subset convolution. His GF machinery may contain seeds of an OCF proof.

---

## Priority D: Computational targets

### INV-024: H(T_p) for Paley primes — EXTENDED
**Source:** OPEN-Q-013, opus-S9
**Status:** COMPUTED through p=23.
**Results:**
- H(P(3))=3, H(P(7))=189, H(P(11))=95095, H(P(19))=1,172,695,746,915, H(P(23))=15,760,206,976,379,349
- All match OEIS A038375 where known (a(3)=3, a(7)=189, a(11)=95095)
- P(7) confirmed as GLOBAL maximizer by exhaustive check of all 2^21 n=7 tournaments
- |Aut(T_19)| = 171, H/|Aut| = 6,857,869,865
- Ratio H(P(p))/(p!/2^{p-1}): 2.000, 2.400, 2.440, 2.527, 2.557 — converging toward e=2.718
**Sequence:** H/|Aut| = 1 (p=3), 9 (p=7), 1729 (p=11), 6857869865 (p=19).
**Next step:** Submit H(P(p)) values and H(P(p))/|Aut| sequence to OEIS. Compute H(P(31)) if feasible (2^31*31 DP — ~66B ops, might take hours).

### INV-025: Integrality conjecture C(p,k) | c_k(T_p) for k >= (p+3)/2
**Source:** T036/T153 (tangents), OPEN-Q-013 table
**Status:** VERIFIED at p=7 (kind-pasteur-S39b). Previously observed at p=11.
**What:** For Paley primes p = 3 mod 4, the cycle count c_k(T_p) is divisible by C(p,k) when k >= (p+3)/2. (CORRECTED from (p+1)/2: at p=11, c_6=1595 is NOT divisible by C(11,6)=462, but c_7=3960 IS divisible by C(11,7)=330.)
**Results at p=7:** c_3=14 (C(7,3)=35 does NOT divide, but k=3 < 4 = (p+1)/2), c_5=42 (C(7,5)=21 DIVIDES, quotient=2), c_7=24 (C(7,7)=1 trivially divides). Conjecture HOLDS.
**Explanation:** Aut(T_p) = Z_p acts on k-subsets, partitioning them into orbits of size p (except the full vertex set which is fixed). Each k-subset orbit has the same cycle count by symmetry. So c_k = p * (cycle count per orbit) when k < p, giving p | c_k. For C(p,k) divisibility: the orbit structure under the full Aut group (which has order p*(p-1)/2 for Paley) should give the stronger divisibility.
**Next step:** Verify at p=19. Prove C(p,k) divisibility from Aut(T_p) orbit counting.

### INV-026: Alpha_1 vs |C_3| mod 2 — systematic test
**Source:** INV-011, oq:mod4_struct
**Status:** TESTED (opus-S5). CONJECTURE IS FALSE.
**Result:** Counterexamples at every n tested:
  - n=3: 2/8 counterexamples (the 3-cycle tournaments have c3=1 odd, alpha_1=0 even)
  - n=4: 16/64 counterexamples (R-cone and near-R-cone tournaments)
  - n=5: 384/1024 counterexamples
  All counterexamples have alpha_1=0 but c3 odd. The conjecture fails because alpha_1 counts independent PAIRS of odd cycles in Omega(T), which is 0 whenever #cycles <= 1.
**Impact:** Open Problem 9 needs reformulation. Alpha_1 ≠ c_3 mod 2 in general.

### INV-027: Realizable conflict graphs catalog
**Source:** INV-013, conflict_graph_catalog.py
**Status:** DONE (opus-S5). Major structural finding.
**Results:**
  - n=3: 2 distinct Omega structures. n=4: 3. n=5: 6. n=6: 24. n=7 (sampled): 172.
  - **Omega(T) is ALWAYS PERFECT** (exhaustive n<=6, 2000 random n=7). This is a significant constraint — independence number = clique cover number.
  - Omega(T) is NOT always chordal (14% non-chordal at n=6, 12% at n=7).
  - For n<=5, Omega is always complete (pigeonhole: any two subsets of size>=3 share a vertex in a 5-element set).
  - At n=6, non-edges correspond exclusively to complementary 3-cycles (vertex sets partition {0,...,5}).
  - Omega can be disconnected at n=6 (80/32768 tournaments, always exactly 2 complementary 3-cycles).
**Impact:** Perfectness of Omega(T) constrains the independence polynomial and could simplify OCF proof strategies. For perfect graphs, Lovasz theta = clique number, and the fractional chromatic number = chromatic number.

### INV-028b: Fix DR mod-4 proof (Thm 7.4 in tex)
**Source:** tex-deep-analysis.md (ISSUE-1)
**Status:** Proof is BROKEN (arithmetic produces v_2 = -2). Result verified for n=3,7,11 only.
**What:** The proof attempts Moon's formula arithmetic but fails. Need proper v_2 analysis using Kummer's theorem, or prove via alpha_1 parity directly (not just |C_3|).
**Next step:** Compute alpha_1 mod 2 for DR_n using OCF. Possibly downgrade to "Verified Conjecture" in tex.

### INV-029b: Fix SE-SYT formula (Thm 7.3 in tex)
**Source:** tex-deep-analysis.md (ISSUE-2)
**Status:** Classical formula cited gives non-integer (2^{3/2} for m=2). Result verified n=5,7.
**What:** Find correct classical reference for SE-SYT count on 2-core shapes. Likely Stembridge (1996) or similar.
**Next step:** Look up Stembridge's "Canonical bases and self-evacuating tableaux." Give clean proof or correct citation.

### INV-030b: Pin grid S_3 symmetry for OCF
**Source:** tex-deep-analysis.md (Section E)
**Status:** NOT explored.
**What:** The S_3 action on barycentric coordinates constrains the polynomial identity. Can it reduce the proof of delta_H = delta_I by exploiting the 6-fold symmetry? The subset convolution lives on Boolean lattice 2^{others} which is a sublattice of the pin grid.
**Next step:** Check if delta_H = delta_I as polynomial has S_3 symmetry. If so, proving it on a fundamental domain suffices.

---

## Priority E: Tangents needing investigation

### INV-028: Hard-core lattice gas at fugacity 2
**Source:** T006, hard_core_lattice_gas.py, hard_core_fast.py
**Status:** INVESTIGATED (opus-S5). Key finding: non-perturbative regime.
**What:** H(T) = I(Omega(T), 2) = Z(Omega(T), lambda=2). Lambda=2 is ABOVE all cluster expansion convergence thresholds for any max degree Delta >= 2:
  - Shearer bound: 1/(Delta-1) << 2 for Delta >= 2
  - LLL/tree bound: (Delta-1)^{Delta-1}/Delta^Delta << 2 for Delta >= 2
  - Kotecky-Preiss: 1/(e*(Delta+1)) << 2 for Delta >= 1
  This means OCF is a non-perturbative identity — standard polymer expansion / cluster expansion methods CANNOT prove it.
**Omega(T) structure (n=4,5):** #cycles dist ranges from 0 to 6 (n=5). Max degree of Omega grows with n. Density is moderate. Independence number = 1 for all n=4 tournaments with cycles (all cycles share vertices).
**Impact:** Rules out perturbative approaches. OCF requires exact cancellations, not convergence arguments.

### INV-029: Ballot sequence / Dyck path connection
**Source:** T001, OPEN-Q-005, ballot_sequence_test.py
**Status:** RESOLVED (opus-S5). Bijective proof FOUND.
**What:** C(L-2, 2k-1) counts signatures with exactly k Type-II positions in an L-cycle window.
**Bijective proof:** The L-cycle through v has L-1 non-v vertices, giving L-1 signature values (s_1=1 forced, s_{L-1}=0 forced, L-3 free). There are L-2 consecutive pairs. Define transition indicators t_j = (s_j != s_{j+1}). Since s starts at 1 and ends at 0, total transitions must be odd. Transitions alternate fall-rise-fall...fall, so k Type-II = (2k-1 transitions + 1)/2. Choosing which 2k-1 of the L-2 positions are transitions gives C(L-2, 2k-1). QED.
**Convention note:** Initial attempt with wrong convention (L-4 free vars, sig length L-2) gave C(L-3, 2k-1). Correct convention: L-1 non-v vertices, sig length L-1, L-3 free vars, L-2 pairs.

### INV-030: Tower hypothesis (L-cycle corrections from (L+2)-cycles)
**Source:** T012, OPEN-Q-012
**Status:** Hypothesis only. NOT tested.
**What:** At n=2k, the first cycle with mu>1 has length 2k-1. Excess from shorter cycles may be compensated by (L+2)-cycle contributions. Is there a recursive tower structure?

### INV-042: Fano-Paley design structure and alpha_2 — MAJOR PROGRESS
**Source:** T102, opus-2026-03-06-S4; T114, kind-pasteur-2026-03-06-S18h
**Status:** PROVED at n=7. BIBD arrangement MINIMIZES alpha_2 but MAXIMIZES H.
**What:** The cyclic triples of Paley T_p form a 2-(p, 3, (p+1)/4) BIBD. At p=7: lambda=2, 7 disjoint pairs = MINIMUM among all regular tournaments. The BIBD is Aut-transitive.
**CRITICAL CORRECTION (S18h):** Previous hypothesis was that BIBD maximizes alpha_2, driving H. This is WRONG. The BIBD actually MINIMIZES alpha_2 (7 vs 10 or 14 for other regular tournaments). But H-maximization is driven by alpha_1 (total DIRECTED odd cycles), not alpha_2. The BIBD forces every 5-vertex subtournament to be regular T_5 (2 directed Ham cycles each), giving 42 directed 5-cycles vs 28-36 for non-BIBD. Combined: alpha_1=80, alpha_2=7, H=189 vs alpha_1=59, alpha_2=14, H=175. Three rigid classes at n=7 (THM-027).
**Formula:** For regular tournaments, D = C(b,2) - p*C(r,2) + sum C(lambda_e, 2). BIBD minimizes the convex sum by Jensen's inequality.
**Next step:** (1) Verify at p=11: does BIBD also maximize directed 5-cycles? (2) Prove that BIBD forces subtournament regularity. (3) Can we prove alpha_1 maximization from BIBD structure at general p?

### INV-043: Paley deletion extended to p=19
**Source:** T104, opus-2026-03-06-S4
**Status:** COMPUTED. Consistent with conjecture.
**What:** H(T_19 - v) = 117,266,659,317 for all vertices v. Scores: (8^9, 9^9), self-complementary. H(T_19)-H(T_19-v) = 1,055,429,087,598 = 2*527,714,543,799.
**Conjecture:** a(18) = 117,266,659,317 in OEIS A038375. Verified chain: a(2)=1=H(T_3-v), a(6)=45=H(T_7-v), a(10)=15745=H(T_11-v).
**Cannot verify** against OEIS (only goes to a(11)=95095).
**Next step:** If someone computes a(18), compare. Or prove Paley deletion gives maximizer using design theory + SC maximizer mechanism.

### INV-031: Lindstrom-Gessel-Viennot (LGV) approach to bijection
**Source:** T046
**Status:** Idea only. NOT attempted.
**What:** The bijection between Ham paths and 2-colored cycle sets, if it exists, might require a global construction like LGV lattice path counting. The non-local nature of the correspondence (T035 dead end) suggests a determinantal approach.

### INV-036: Tiling grid geometry and class structure
**Source:** opus-2026-03-06-S1 (deep tiling investigation)
**Status:** INVESTIGATED. Key structural findings.
**What:** How does the {0,1}^m tiling space geometry relate to tournament isomorphism classes?
**Results:**
- **Sigma (converse) acts cleanly on classes:** sigma permutes bits (no complement), preserves weight. Self-converse: 2,2,8,12 classes at n=3,4,5,6. Sigma-fixed tilings = 2^floor((n-1)^2/4).
- **Complement does NOT respect classes:** unlike sigma, flipping all non-path arcs does not map classes to classes.
- **Standard invariants almost distinguish:** At n=6, (score, c3, c5, omega_deg, H) fails only for sigma pairs (converse-paired classes) plus occasional self-converse coincidences.
- **Triangle 3-cycle probability:** P=1/2 for consecutive triples (path arcs), P=1/4 for all others. E[c3] = (C(n,3) + n-2)/4.
- **Strong H~c3 correlation:** r=0.956 at n=5,6. H = 1+2c3 exact at n<=4, breaks at n>=5.
- **Bit-position variance:** Longest arc (gap=n-1) most predictive of class. Middle arcs vary most.
- **Class transition graph:** Always connected. ΔH always even. E[ΔH]=0 for every arc position.
- **Weight distributions distinguish:** Can separate classes sharing all tournament invariants.
**Full writeup:** `03-artifacts/drafts/tiling-symmetry-analysis.md`
**Scripts:** `04-computation/tiling_*.py` (5 files)
**Next step:** (1) Investigate which tiling properties predict H beyond c3. (2) Connect sigma reduction to arc-flip proof strategy. (3) Look for grid-local rules that determine class.

### INV-037: Pin-grid sigma vs tournament sigma — two-sigma structure
**Source:** opus-2026-03-06-S2 (sigma structure investigation)
**Status:** INVESTIGATED. Clean structural results, but no proof path yet.
**What:** The pin-grid sigma (r,c)->(c,r) and tournament sigma (i,j)->(n-1-j,n-1-i) are DIFFERENT symmetries. Pin sigma acts within strips; tournament sigma acts across strips. They agree only on diagonal r=c.
**Key results:**
- **POS-free identity:** free(strip k) = cumul_POS(k) = floor(k/2). Growth rate: delta_free(k) = POS(k) = [k even].
- **n->n+2 structure:** Adds strips n and n+1 with exactly n sigma-free bits and exactly 1 POS (midpoint arc).
- **Tournament sigma always preserves H** (converse operation, verified n=3,...,7).
- **Pin-grid sigma does NOT preserve H** in general (only 5% at n=7).
- **Two sigmas don't commute;** composition has order 3; generate S_3-like group.
- **Mod-4 structure:** Neither sigma preserves H mod 4 reliably.
**Scripts:** `04-computation/sigma_structure.py`
**Next step:** (1) Understand algebraic significance of the S_3 group. (2) Can the n->n+2 POS structure be used differently (not through H preservation)? (3) Relate to transfer matrix symmetry (INV-001).

### INV-039: SC Maximizer Theorem and sigma* structure
**Source:** kind-pasteur-2026-03-06-S18/S18e, opus-2026-03-06-S4, OPEN-Q-016
**Status:** VERIFIED exhaustive n=4,5,6,7. Mechanism deeply analyzed. NOT proved.
**What:** Within each self-complementary score class, max H is always achieved by SC tournament. The mechanism: involutory anti-automorphism sigma induces sigma* on directed odd cycles, which is an involutory automorphism of Omega(T). At even n, sigma* is fixed-point-free, pairing all cycles. Some pairs are vertex-disjoint (giving alpha_2 contributions). At even n, sigma is fixed-point-free on vertices (proved: fixed point implies score=(n-1)/2, non-integer).
**Key results (opus-S4 deepened):**
- **PROVED: sigma always induces Omega automorphism** (clean proof: sigma maps directed C to reverse of sigma(C), preserving vertex-sharing)
- **PROVED: At even n, sigma* has NO fixed cycles** (3-cycles: can't fix set of 3 with fpf involution; 5-cycles at n=6: can't fix set of 5 with 3 two-cycles)
- SC and NSC have the SAME number of 3-cycles within score class — the difference is in ARRANGEMENT (disjoint pairs vs all overlapping)
- SC max alpha_2 >= NSC max alpha_2 within every score class at n=6
- **Path reversal identity:** M_{T^op}[i,j] = (-1)^{n-2} M_T[j,i] (proved)
- **At odd n=5:** alpha_2 = 0 for ALL SC tournaments (fixed point forces all cycles to overlap)
- **At odd n=7 (Paley):** 21 anti-auts, 7 involutions (one per fixed point), each finds 1 disjoint 3-cycle pair
**Scripts:** sc_maximizer_n7_fast.py, anti_aut_analysis.py, anti_aut_exhaustive.py, clique_antiaut_connection.py
**Draft:** sc-maximizer-mechanism.md
**Next step:** (1) Test at n=8 (even n). (2) Prove SC always achieves max alpha_2 within score class. (3) Formalize the "arrangement advantage" into an algebraic proof.

### INV-040: Paley deletion gives H-maximizer
**Source:** kind-pasteur-2026-03-06-S18e, opus-2026-03-06-S4
**Status:** VERIFIED at p=3,7,11. Conjecture.
**What:** Deleting any vertex from Paley tournament T_p gives a tournament with H = max H at n=p-1 (= OEIS A038375(p-1)).
**Results:**
- T_3 → T_2: H=3 → H=1 = a(2) ✓
- T_7 → T_6: H=189 → H=45 = a(6) ✓
- T_11 → T_10: H=95095 → H=15745 = a(10) ✓
- All vertex deletions give the same H (by Aut(T_p) transitivity)
- T_11 - v has self-complementary scores (4,4,4,4,4,5,5,5,5,5)
- H(T_p) - H(T_p-v) = 2 * (sum of mu-weighted cycles through v) (Claim A)
**Conjecture:** T_p - v is the GLOBAL H-maximizer at n = p-1 for all Paley primes p ≡ 3 mod 4.
**Next step:** (1) Verify at p=19 (need H(T_19 - v), n=18 DP ~2^18*18 ~ 5M, feasible). (2) Test whether T_p - v is SC. (3) Relate to lattice theory or QR structure.

### INV-038: Blueself parity theorem and census structure
**Source:** opus-2026-03-06-S3 (deep census investigation)
**Status:** THM-023 PROVED. Census in progress through n=8.
**What:** Blueself (GS + self-flip) exists if and only if n is even. Proved algebraically: flip changes endpoint scores by score'(0) = n - score(0), so same-score requires score(0) = n/2 (integer only at even n).
**Census results (exhaustive n=3,...,6, in progress n=7,8):**
- POS orientation is perfectly UNIFORM: each pattern gets exactly 2^(m-#POS) tilings
- GS POS is also perfectly UNIFORM
- SC always maximizes H within each score sequence class (confirmed with kind-pasteur findings)
- Blueself at n=4: H=5 (rank 1/4), n=6: H=41,45 (ranks 5,1/56) — near or at global maximum
- Blackself at odd n is in SC classes; at even n exclusively in NSC (paired) classes
- SF tilings come in flip-pairs; SF count per class is 2 at n=6, 4 at n=5
- Self-flip fraction decreases: 25%, 12.5%, 1.56% at n=4,5,6
**Scripts:** `04-computation/deep_census_analysis.py`, `04-computation/pos_tiling_census.py`, `04-computation/census_n8.py`
**Theorem:** `01-canon/theorems/THM-023-blueself-parity.md`
**Next step:** (1) Complete n=7 and n=8 census. (2) Investigate why blueself achieves max H. (3) Count blueself at n=8 (1280 eligible GS tilings, need canonicalization).

---

## Completed / Closed investigations

- [DONE] OCF verified n<=8 exhaustive (opus-S3/S4)
- [DONE] Transfer matrix symmetry discovered and verified (opus-S4b)
- [DONE] Locally transitive tournaments tested — DO have 5/7-cycles (opus-S4b)
- [DONE] Feng + Rajkumar connections documented in paper-connections.md (opus-S4b)
- [DONE] T_11 cycle table complete, H(T_11)=95095 confirmed (kind-pasteur-S2/S5)
- [DONE] Per-path identity failure characterized (THM-009)
- [DONE] Even-odd split is consequence not equivalent to OCF (MISTAKE-008)
- [DONE] Bracket structure B(u,w) analyzed (T047, bracket_structure.py)
- [DONE] H(T_19) computed: 1,172,695,746,915; H/|Aut|=6,857,869,865 (opus-S5)
- [DONE] Deep paper analysis: SST, Rajkumar, Feng — all key theorems extracted (opus-S5)
- [DONE] Ballot sequence bijective proof for C(L-2, 2k-1) (opus-S5)
- [DONE] Hard-core lattice gas: lambda=2 is non-perturbative regime (opus-S5)
- [DONE] Alpha_1 ≡ c_3 (mod 2) conjecture DISPROVED (opus-S5)
- [DONE] Conflict graph catalog: Omega(T) is PERFECT for n<=7, FAILS at n>=8 (53.8% have C5 in Omega_3). See OPEN-Q-014.
- [DONE] Omega(T) is CLAW-FREE for n<=8, FAILS at n>=9 (90%). See OPEN-Q-014.
- [DONE] Web research: 9 new connections documented in web-research-connections.md (opus-S5)
- [DEAD] Per-vertex decomposition of unmatched counts (T045)
- [DEAD] Cycle bijection under arc reversal (MISTAKE-005)
- [DEAD] Contiguous block decomposition (T035)
- [DEAD] Contraction approach (T017)
- [DEAD] Alpha_1 = c_3 (mod 2) conjecture — counterexamples at all n (opus-S5)
- [DONE] Paley maximizer conjecture verified at p=3,7,11 (exhaustive at p=7), extended with H(P(19)), H(P(23)) (opus-S9)
- [DONE] n=8 H-maximizer identified: H=661=a(8), self-converse, |Aut|=1, does NOT contain P(7) (opus-S9)
- [DONE] Full Omega structure at n=8: 76-78 vertices, density 0.98, 70-75% of H from 5/7-cycles (opus-S9)
- [DONE] Ratio H(P(p))/(p!/2^{p-1}) converges toward e: 2.00, 2.40, 2.44, 2.53, 2.56 for p=3,7,11,19,23 (opus-S9)
- [DONE] Deep web synthesis: Irving-Omar, Grinberg-Stanley, Grujić-Stojadinović, Feng, DRT theory (opus-2026-03-06-S5)

---

## Priority F: New leads from web synthesis (opus-2026-03-06-S5)

### INV-045: Hopf algebra route to transfer matrix symmetry — SUPERSEDED by THM-030
**Source:** T114, Grujić-Stojadinović arXiv:2402.07606, Feng arXiv:2510.25202
**Status:** INVESTIGATED (kind-pasteur-S22 agent). No direct proof obtained.
**What:** The Hopf comultiplication Δ([T]) = Σ_S [T|_S]⊗[T|_{V\S}] encodes our subset convolution. Feng's dual Burnside process proves Q=AB is symmetric under detailed balance. The tournament constraint T[x,y]+T[y,x]=1 is the detailed balance condition.
**Finding:** Feng's framework FAILS because it requires positivity (probability transitions), but our M = E^T * Lambda * B has Lambda = diag((-1)^|S|) with negative entries. Grujic-Stojadinovic gives U_X = U_{X^op} (proved) but this is a GLOBAL identity, not per-entry M[a,b]=M[b,a]. The Hopf comultiplication IS cocommutative but this doesn't directly imply transfer matrix symmetry because E_a and B_b are different types of objects. **Most promising remaining paths:** (1) inductive s-variable approach, (2) Irving-Omar det/per formula (INV-046), (3) "signed Feng" extension (new theorem needed), (4) pointed Hopf algebra tracking distinguished vertices.
**Why it matters:** Transfer matrix symmetry ⟹ OCF (via even-odd split + s-coefficient identity). A Hopf algebra proof would be clean and publication-worthy.
**Next step:** (1) Express our transfer matrix M in Feng's AB framework precisely. (2) Verify detailed balance condition algebraically. (3) Apply Feng Thm 3.3 to conclude symmetry.

### INV-046: Irving-Omar det/per formula → transfer matrix — SUPERSEDED by THM-030
**Source:** T118, Irving-Omar arXiv:2412.10572 Proposition 2
**Status:** SUPERSEDED. Transfer matrix symmetry now proved directly (THM-030, opus-S25) without needing Irving-Omar det/per.
**What:** ham(D) = Σ_S det(Ā[S])·per(A[S^c]). The walk generating function W_D(z)=det(I+zXĀ)/det(I-zXA) encodes Hamiltonian path structure.
**Remaining interest:** Irving-Omar's framework may still provide insight into WHY the Key Identity works — e.g., is there a matrix-algebraic interpretation of B_b + (-1)^m E_b = 2r·col_sum?

### INV-047: Paley maximizer via DRT theory
**Source:** T116, Reid-Brown 1972, Nozaki-Suda arXiv:1202.5374
**Status:** CLASSICAL EQUIVALENCE + our computational evidence.
**What:** DRTs ↔ skew Hadamard matrices. Paley T_p is the canonical DRT. Nozaki-Suda characterize skew Hadamard via spectra of tournaments of size n-2. Our spectral regularity finding (corr(H,λ₁)=-0.97) + DRT theory could explain why Paley maximizes H.
**Key question:** Among all DRTs on p vertices, does Paley ALWAYS maximize H? Or is this specific to Paley among all tournaments?
**Next step:** (1) Check if non-Paley DRTs exist at small p and compare H values. (2) Relate DRT cycle balance to alpha_k maximization.

### INV-048: Asymptotic convergence H(T_p)/(p!/2^{p-1}) → e
**Source:** T117, Adler-Alon-Ross 2001
**Status:** COMPUTATIONAL EVIDENCE. Not proved.
**What:** Adler-Alon-Ross proved max H(T) ≥ (e-o(1))·n!/2^{n-1} using random regular tournaments. Our Paley ratios 2.00→2.56 suggest convergence to e. Paley tournaments are quasi-random in Chung-Graham-Wilson sense.
**Next step:** (1) Compute H(T_p)/(p!/2^{p-1}) for p=31 if feasible. (2) Check if quasi-randomness implies near-optimal H. (3) Try Stirling approximation on the cycle-count formula.

### INV-049: El Sahili-Ghazo Hanna type-preserving converse symmetry
**Source:** arXiv:2101.00713 (2023), J. Graph Theory 102
**Status:** PUBLISHED RESULT, connection to our work identified.
**What:** T and T^op have the same number of oriented Hamiltonian paths of EVERY type. Our transfer matrix identity M_{T^op} = (-1)^{n-2} M_T is a STRONGER result that implies this. The El Sahili result follows from transfer matrix symmetry + path reversal.
**Impact:** Our transfer matrix results strengthen known literature results.
**Next step:** Note in paper draft. Check if Ai (2025) "New Digraph Polynomials" extends further.

### INV-050: Satake's cyclotomic NDRTs (arXiv:2502.12090, Feb 2025)
**Source:** T116 area
**Status:** IDENTIFIED. Not investigated.
**What:** Nearly-doubly-regular tournaments from almost difference sets. Savchenko's conjecture on canonical spectrum. Under Hardy-Littlewood conjecture F, infinitely many NDRTs with canonical spectrum exist.
**Next step:** Read paper. Check if NDRTs approach Paley's H-maximization. Compare spectra.

---

## Priority G: New leads from web search (kind-pasteur-2026-03-06-S19)

### INV-051: Mitrovic noncommuting Rédei-Berge function (arXiv:2504.20968, Apr 2025)
**Source:** kind-pasteur-2026-03-06-S19 web search; opus-2026-03-06-S9 (detailed paper read)
**Status:** DEEPLY READ — HIGH PRIORITY
**What:** Stefan Mitrovic introduces the Rédei-Berge function in NONCOMMUTING variables, which satisfies deletion-contraction (W_X = W_{X\e} - W_{X/e}↑). The commutative version does NOT have deletion-contraction. Key properties: W_X = W_{X^op}, product rule W_{X·Y} = W_X·W_Y. For tournaments: W_X = sum over permutations with all odd cycles of 2^{psi(sigma)} p_{Type(sigma)} with positive integer coefficients.
**Why it matters:** Deletion-contraction enables INDUCTIVE PROOFS. This could provide an inductive framework for OCF or transfer matrix symmetry. The noncommutative structure preserves more information than the commutative version.
**TESTED (kind-pasteur S19):** Direct deletion-contraction does NOT preserve OCF. At n=4:
  - H(T) = H(T\e) - H(T/e): only 18.8% match (DC is for W_X not H)
  - OCF for T\e (non-tournament): only 39.3% hold
  - OCF for T/e (contracted): only 60.7% hold
  OCF is TOURNAMENT-SPECIFIC and does not hold for general digraphs from deletion/contraction.
  The noncommuting framework operates at a different level than H(T).
**DETAILED READING (opus-S9):**
  - Definition 3.1: W_X = sum_{f:V->P} sum_{sigma in Sigma_V(f,X)} delta_f(sigma) x_{f(v1)}...x_{f(vn)} where x_i are NONCOMMUTING. Depends on vertex labeling (unlike commutative U_X).
  - Theorem 3.7 (Deletion-Contraction): W_X = W_{X\e} - W_{X/e}↑ where e=(v_{n-1},v_n). The ↑ operation doubles the last variable: (x_{i1}...x_{in-1})↑ = x_{i1}...x_{in-2} x_{in-1}^2. This is the KEY technical device.
  - Theorem 3.10: W_X = sum_{sigma in S_V(X,Xbar)} (-1)^{phi(sigma)} p_{Type(sigma)} where Type is now a SET PARTITION (not integer partition). The noncommutative p_pi tracks which vertices are in which cycle.
  - Corollary 3.12: For tournaments, W_X = sum_{sigma in S_V(X), all odd cycles} 2^{psi(sigma)} p_{Type(sigma)}. This is the NONCOMMUTATIVE OCF — same sum but tracking set-partition cycle types.
  - Key insight: The ↑ operation has no obvious tiling interpretation, but it corresponds to "merging the last two vertices while remembering which was which." This is precisely what our contraction approach (T017) failed to handle in commutative setting.
**Next step:** (1) Instead of naive DC -> OCF induction, investigate whether DC can be used to relate the SYMMETRIC FUNCTION U_T across tournaments (e.g., U_T = U_{T'} + correction for single arc reversal). (2) Check if the noncommutative framework gives a new proof of even-odd split or transfer matrix symmetry. (3) CRITICAL: Investigate whether Theorem 3.10 (noncommutative power-sum expansion over SET partitions) implies transfer matrix symmetry when projected to integer partitions.

### INV-052: Mitrovic-Stojadinovic chromatic↔Rédei-Berge connection (arXiv:2506.08841, Jun 2025)
**Source:** kind-pasteur-2026-03-06-S19 web search; opus-2026-03-06-S9 (FULL PAPER READ)
**Status:** DEEPLY READ — HIGH PRIORITY. Contains multiple directly relevant results.
**What:** Proves the chromatic function of a graph and the Rédei-Berge function of a digraph are "almost identical" at the poset level. For any poset P: X_{inc(P)} = omega(U_P) (Theorem 3.2). This extends to noncommutative versions: Y_{inc(P)} = omega(W_P) (Theorem 5.7). Commutative diagram of four Hopf algebra morphisms (p.8, Remark 3.4).
**KEY RESULTS FROM FULL READ (opus-S9):**
  1. **Theorem 4.1 (Converse of Rédei):** If P is a poset that is NOT a chain, then the number of quasi-linear extensions is EVEN. Proof: chi_{inc(P)}(1) and (-1)^|P| u_P(-1) have the same parity. u_P(1) counts quasi-linear extensions; chi_{inc(P)}(1)=0 unless inc(P) is discrete (P is a chain). THIS GENERALIZES RÉDEI'S THEOREM TO POSETS.
  2. **Theorem 4.8:** For any poset P and partition lambda: #{broken-cycle-free subsets of inc(P) with component sizes lambda} = #{permutations in S_V(D_P-bar) with cycle type lambda}. This is a CYCLE-TYPE-PRESERVING bijection between broken cycles and permutations — potentially a route to bijective OCF proof.
  3. **Corollary 3.3:** chi_{inc(P)}(m) = (-1)^|P| u_P(-m). The chromatic polynomial at m colors = Rédei-Berge polynomial at -m. Our H(T) = u_T(1) = (-1)^n chi_{inc(P)}(-1) for the associated poset.
  4. **Theorem 4.6:** D_P-bar contains a Hamiltonian cycle iff P is irreducible. Combined with Corollary 4.5: U_{P_n} forms an algebraic basis of Sym for any sequence of irreducible posets.
  5. **Section 6 (Positivity):** Conjecture 6.3: If P is (3+1)-free, U_P is h-positive. Theorem 6.4: Already proved s-positive. Theorem 6.11: h-positivity propagates through deletion-contraction if there's a "sink-like" vertex.
  6. **Section 7 (Bags of sticks):** Decomposition into bags of sticks (disjoint unions of directed paths) gives explicit formulas. Triple deletion property generalized.
**Why it matters for us:**
  - The bridge X_{inc(P)} = omega(U_P) means we can use 30 years of chromatic symmetric function theory on tournament problems
  - Theorem 4.1 is a clean generalization of the parity theorem we're studying
  - Theorem 4.8 gives a TYPE-PRESERVING bijection — this is exactly what a bijective OCF proof would need
  - The h-positivity results (Section 6) may apply to tournament-associated posets
**Next step:** (1) For a tournament T, identify the associated poset P such that D_P = T. (2) Check if Theorem 4.8 gives a new proof of OCF when specialized to tournaments. (3) Investigate whether tournament posets are (3+1)-free (this would give h-positivity of U_T). (4) Test computationally: does the broken-cycle bijection from Thm 4.8 match our OCF terms?

### INV-053: Savchenko cycle counting formulas for regular tournaments — VERIFIED AT n=7
**Source:** kind-pasteur-2026-03-06-S19 web search; Savchenko J. Graph Theory 83 (2016), Discrete Math (2017), arXiv:2403.07629 (2024)
**Status:** VERIFIED at n=7. Cycle counts are class invariants. DRT vs LTT classification matches.
**What:** Savchenko has a series of papers giving EXACT polynomial formulas for c_k(T) (number of k-cycles) in regular tournaments:
- c5, c6 formulas (2016, J. Graph Theory 83)
- c7 formula (2017, Discrete Math)
- c8 for DRTs vs locally transitive tournaments (2024, arXiv:2403.07629)
Key finding: c8(DRT_n) is INDEPENDENT of which DRT is chosen. Phase transition at n=39: DRTs have more 8-cycles than locally-transitive for n≤35 but FEWER for n≥39.
**Why it matters:** These exact formulas could determine whether Paley tournaments maximize cycle counts at EVERY length, or only for short cycles. The phase transition at n=39 suggests our cycle-maximization mechanism may reverse at larger n. Also, the spectral methods used (eigenvalue-based cycle counting) could connect to our transfer matrix work.
**VERIFIED (kind-pasteur S19):** At n=7, the three regular tournament classes are EXACTLY:
  - DRT (Paley): 240 tours, dc={3:14, 5:42, 7:24}, H=189
  - Locally Transitive: 720 tours, dc={3:14, 5:28, 7:17}, H=175
  - Other Regular: 1680 tours, dc={3:14, 5:36, 7:15}, H=171
Cycle counts are CLASS INVARIANTS (exactly one vector per class). DRT maximizes directed 5-cycles and 7-cycles. LTT has "diametrically opposite" properties per Savchenko.
**EXTENDED (kind-pasteur S21):** Savchenko (2024) proves c_m(DR_n) > c_m(RLT_n) for ALL m = 1,2,3 mod 4 (including all odd m). Only m = 0 mod 4 has the phase transition. This directly explains DRT's H-maximization via OCF.
**DRT n=11 ANALYSIS (kind-pasteur S21, CORRECTED S39b per MISTAKE-017):** The connection set {1,2,3,5,8} is NOT a valid tournament (S∩(-S)={3,8}≠∅). The ONLY valid circulant DRT at n=11 is the Paley tournament QR={1,3,4,5,9} (H=95095, c3=55, c5=594, |Aut|=55). All claims about "non-Paley DRT" H=69311, c3=44 were computed on an invalid digraph. Whether a non-circulant DRT exists at n=11 remains open.
**Next step:** (1) Obtain Savchenko's exact polynomial c_k formulas. (2) Test at n=19 or n=23 (multiple DRT classes). (3) Prove Paley maximizes H among ALL DRTs.

### INV-054: Komarov-Mackey exact 5-cycle formula (arXiv:1410.6828, JGT 2017) — PARTIALLY INVESTIGATED
**Source:** kind-pasteur-2026-03-06-S19 web search
**Status:** PARTIALLY INVESTIGATED (kind-pasteur-S39b). c5 is NOT score-determined.
**What:** Exact formula for number of directed 5-cycles in any tournament in terms of edge score sequence. Maximum c5 ≈ (3/4)*C(n,5), achieved by almost all random tournaments. Lower bounds also proved.
**NEW FINDING (S39b):** c5 is NOT determined by score sequence, even combined with sum_d², edge_score, or common_out_neighbor statistics. At n=5, score (1,2,2,2,3) has c5 in {1,2,3}; at n=6, 9/22 score sequences have varying c5. The Komarov-Mackey formula likely involves CUBIC or higher-order graph statistics (e.g., directed walks of length 3+). For regular tournaments, c5 IS a class invariant (Savchenko, verified n=5,7).
**Why it matters:** This rules out O(n²) c5 computation from scores alone. Cycle enumeration (O(n^5) for 5-cycles, or O(2^n) bitmask DP) remains necessary.
**RESOLVED (S39b, THM-118):** c_5 = tr(A^5)/5 gives O(n^3) computation via matrix multiplication. This IS the "cubic invariant" — tr(A^5) is a sum over all length-5 closed walks, and THM-118 proves all such walks in tournaments are simple cycles (no vertex repetition possible for length <= 5).
**Next step:** Read Komarov-Mackey formula to see if it matches tr(A^5)/5.

### INV-055: Linial-Morgenstern cycle density conjecture and extremal tournaments
**Source:** kind-pasteur-2026-03-06-S19 web search; arXiv:2011.14142 (Ma-Tang), arXiv:1902.00572
**Status:** NEW LEAD — MEDIUM PRIORITY
**What:** Linial-Morgenstern conjecture: among tournaments with fixed c3 density d, the c4 density is minimized by random blowups of transitive tournaments. Proved for d ≥ 1/36 using spectral methods. Ma-Tang extend to c_ℓ for ℓ ≢ 2 mod 4 when d is near 1.
**Why it matters:** This is the "dual" to our maximization question. We show Paley maximizes total directed cycles; this literature characterizes minimizers. The spectral methods used here (eigenvalue-based cycle density bounds) could provide tools for our Paley maximizer proof.
**Next step:** Check if the extremal results constrain H(T) via OCF.

### INV-056: Jerrum-Patel zero-free regions for H-free graphs (JLMS 2026)
**Source:** kind-pasteur-2026-03-06-S19 web search
**Status:** NEW LEAD — MEDIUM PRIORITY (for real-rootedness question)
**What:** Jerrum & Patel (2026, JLMS) prove zero-free regions for the independence polynomial of H-free graphs for various H. For claw-free: all zeros on negative real line (= Chudnovsky-Seymour). For subdivided claws: related zero-free regions. KEY: for H NOT a subdivided claw or path, there exist H-free graphs of max degree 3 with zeros NOT on the negative real line.
**Why it matters:** Our Omega_3(T) has all real roots for n≤20 but is NOT always claw-free (fails n≥9). Jerrum-Patel's results on subdivided claw avoidance may explain why real roots persist beyond n=8. The tournament-specific constraint on Omega_3 structure may ensure avoidance of exactly the "bad" subgraphs.
**Next step:** (1) Check what specific subdivided claws appear in Omega_3(T) at n≥9. (2) Apply Jerrum-Patel to determine if their zero-free regions explain our observations.

### INV-057: Herman's Terwilliger algebras of DRTs (arXiv:2404.11560, 2024)
**Source:** kind-pasteur-2026-03-06-S19 web search
**Status:** NEW LEAD — LOW-MEDIUM PRIORITY
**What:** Allen Herman computes Terwilliger algebras for DRTs (asymmetric rank-3 association schemes). Thin irreducible modules, dimension 4k+9. Key: Terwilliger algebras distinguish non-isomorphic DRTs up to n=23, but FAIL at n=27 (need rational Terwilliger algebras). There are 237 non-isomorphic DRTs at n=27.
**Why it matters:** (1) If all DRTs at small n have the same H(T), that would be a DRT invariant. (2) If Terwilliger algebra structure constrains H(T), this gives an algebraic route to Paley maximization. (3) The n=27 DRT classification gives test cases for our conjectures beyond Paley primes.
**Next step:** (1) Check if all DRTs at n=7 (there's only one: Paley) or n=11 have the same H. (2) At n=27, compare H across different DRT isomorphism classes.

### INV-058: Pantangi critical groups distinguish Paley from other DRTs
**Source:** kind-pasteur-2026-03-06-S19 web search; Pantangi arXiv:1905.08568 (2019)
**Status:** CONNECTION IDENTIFIED
**What:** Pantangi shows critical groups (sandpile groups) distinguish Paley from non-Paley DRTs. Chandler-Sin-Xiang computed Smith/critical groups of Paley GRAPHS. Different DRT constructions (Szekeres-Whiteman 2-block, Wallis-Whiteman 4-block) are distinguished by their critical groups.
**Why it matters:** If H(T) is a DRT invariant AND different DRTs have different critical groups, then H could be read off the critical group. This would give a purely algebraic characterization of the H-maximizer.
**Next step:** Compute critical groups for DRTs at n=11,19 and check correlation with H values.

### INV-062: Universal Master Polynomial and Central Factorial Numbers — VERIFIED (THM-059)
**Source:** opus-2026-03-06-S30
**Status:** VERIFIED computationally (22/22 cases, n=4..9). Algebraic proof pending.
**What:** The per-invariant r-polynomial C_I(r,n) = 2^{parts(I)} * F_f(r) where f = free position count and F_j is determined by the central factorial number triangle (OEIS A036969) via b_{k,j} = b_{k-1,j-1} + j^2 * b_{k-1,j}. This completely determines the entire W-coefficient hierarchy for all tournaments at all n. The shift principle is a corollary.
**Key findings:** (1) F_j(1/2) = 1 for all j. (2) Leading coefficient = (j+1)!. (3) Predictions made for n=11 without computation. (4) Complete n=8 even-n table computed.
**Next step:** (1) Algebraic proof of the central factorial recurrence from position pattern analysis. (2) Verify F_8 prediction at n=11 computationally. (3) Investigate C_0(r) (constant/background polynomial).
**Scripts:** `04-computation/universal_master_polynomial.py`, `04-computation/w1_n8_complete.py`

### INV-059: Cyclic subsets of tournaments (arXiv:2508.03634, Aug 2025)
**Source:** kind-pasteur-2026-03-06-S19 web search; Hunter-Liu-Milojević-Sudakov
**Status:** NEW LEAD — LOW PRIORITY
**What:** Optimal bounds on probability that random induced subtournament of high min-degree tournament is Hamiltonian. Extends to p-biased measure. Proves that high min-degree forces high Hamiltonicity probability.
**Why it matters:** Paley T_p has min-degree (p-1)/2 (doubly regular). This paper could give explicit bounds on the fraction of induced subtournaments that are Hamiltonian, which connects to our cycle counting.
**Next step:** Apply their bounds to Paley tournaments. Check if this gives lower bounds on c_k counts.

### INV-060: Eulerian cycle trace formula (arXiv:2502.02915, Feb 2025) — CLOSED
**Source:** kind-pasteur-2026-03-06-S19 web search; Ye Luo
**Status:** CLOSED (too remote, kind-pasteur-S22 agent investigation)
**What:** Trace formula counting Eulerian cycles via "twisted" vertex and edge adjacency matrices. Uses homological spectral graph theory.
**Finding:** Connection too remote to be actionable. Luo counts Eulerian cycles (edge traversals) via H_1 characters; we count Hamiltonian paths (vertex traversals) via inclusion-exclusion. Different domain, different algebraic structure. Closest classical analogue to our transfer matrix is Ryser's permanent formula, not Luo's trace formula.
**Why it matters:** Our transfer matrix tr(M) = H(T) (THM-027) is also a trace formula. This paper's approach—using twisted adjacency matrices with spectral antisymmetry—could provide a template for proving our trace formula properties (symmetry, off-diagonal sum) at general n.
**Next step:** Read the paper and check if "twisted adjacency" techniques apply to tournament transfer matrices.

### INV-061: Hamilton transversals in tournaments (Combinatorica 2024)
**Source:** kind-pasteur-2026-03-06-S19 web search; Chakraborti-Kim-Lee-Seo arXiv:2307.00912
**Status:** NEW LEAD — LOW PRIORITY
**What:** For collections of sufficiently many tournaments on the same vertex set, transversal Hamilton paths/cycles exist. For m=|V|-1 tournaments, there's a transversal Ham path; for m=|V| with m-1 strongly connected, transversal Ham cycle.
**Why it matters:** The "transversal" perspective could give a new way to relate Ham paths across different tournaments, potentially connecting to how H(T) changes under arc reversals.

### INV-062: Forward arc maximization in tournaments (arXiv:2602.10713, Feb 2026)
**Source:** kind-pasteur-2026-03-06-S19 web search; Guo-Gutin-Lan-Shao-Yeo-Zhou
**Status:** NEW LEAD — LOW PRIORITY
**What:** Characterizes maximum forward arcs in Hamilton cycles/paths for semicomplete and locally semicomplete digraphs. Polynomial-time algorithms.
**Why it matters:** Forward arcs in Hamilton paths relate to our "position-based" analysis (pos(a,P) in THM-027 trace formula). The maximum forward arc structure could inform transfer matrix properties.

### INV-063: Spectral pseudorandomness and Paley clique bounds (Exp. Math. 2024)
**Source:** kind-pasteur-2026-03-06-S19 web search; Kunisky-Yu arXiv:2303.16475
**Status:** NEW LEAD — LOW PRIORITY
**What:** Studies spectral pseudorandomness of Paley graphs via subgraph eigenvalue distributions. Conjecturally, minimum eigenvalue convergence would improve clique number bounds beyond √p.
**Why it matters:** Spectral properties of Paley graphs/tournaments are central to our theory. If Paley tournaments have stronger spectral pseudorandomness than other DRTs, this could explain H-maximization via eigenvalue-based cycle counting formulas.

### INV-064: Mitrovic Hopf algebra new bases (arXiv:2407.18608v3, Mar 2025)
**Source:** kind-pasteur-2026-03-06-S19 web search
**Status:** CONNECTION IDENTIFIED — supplements INV-033
**What:** Introduces two new combinatorial Hopf algebras of posets and permutations with Rédei-Berge functions. Constructs new bases for symmetric functions whose generators are Rédei-Berge functions. Investigates which digraph invariants are extractable from the Rédei-Berge function.
**Why it matters:** If H(T) can be expressed as a coefficient in one of these new bases, it gives an algebraic handle on Hamiltonian path counting.
**Next step:** Check which digraph invariants the paper extracts. Is H(T) among them?

### INV-065: Independence polynomial root gap (arXiv:2510.09197, FSTTCS 2025)
**Source:** kind-pasteur-2026-03-06-S19 web search; Om Prakash & Vikram Sharma
**Status:** NEW LEAD — LOW PRIORITY
**What:** Quantifies the gap between the smallest real root β(G) of I(G,x) and all other roots. For connected graphs, β(G) is a simple real root smaller than 1, but previous proofs gave no gap bound. This paper provides explicit bounds.
**Why it matters:** For our Omega(T) real-rootedness question, having a gap bound could help prove that all roots are real by showing they're well-separated from the complex plane.

### INV-066: Low-rank matrices from tournaments and symmetric designs (arXiv:2401.14015, 2024)
**Source:** kind-pasteur-2026-03-06-S19 web search; Balachandran-Sankarnarayanan
**Status:** NEW LEAD — LOW-MEDIUM PRIORITY
**What:** Constructs symmetric matrices from tournament structures where rank depends on design-theoretic properties. Symmetric designs (BIBDs) give matrices with rank near n/2. The rank-topology relationship involves bipartite graph eigenvalues.
**Why it matters:** Our transfer matrix M is constructed from a tournament and is symmetric. This paper's framework connecting tournament-derived matrices with design theory could explain structural properties of M (e.g., why symmetry holds, what the rank structure is).
**Next step:** Check if our M fits their M_T(f,a) framework.

### INV-067: Alpha_1 gap theorem and converse of Redei — CORRECTED S22
**Source:** kind-pasteur-2026-03-06-S21 (computation), CORRECTED S22
**Status:** PARTIALLY PROVED (THM-029 corrected)
**What:** alpha_1=3 is impossible at n<=6 but ACHIEVABLE at n>=7 (~9.2% of c3=3 at n=7). Common-vertex property fails at n>=7. HOWEVER, H=7 remains impossible for ALL n by refined argument: H=7 requires (alpha_1=3, i_2=0), but i_2=0 forces common vertex => c5>=1 => alpha_1>=4; while alpha_1=3 implies i_2>=1 => H>=11. H=21 absent through n=7.
**Achievable H values:** n=5: {1,3,5,9,11,13,15}. n=6: {1,3,5,9,11,...,45}\{7,21,35,39}. At n=7: 35 and 39 become achievable but 7 and 21 remain gaps.
**Connection:** Relates to Mitrovic-Stojadinovic "converse of Redei" (INV-052).
**Next step:** (1) Prove alpha_1=3 impossibility for ALL n (not just n<=6). (2) Find all permanent H-gaps. (3) Check OEIS for the sequence of achievable alpha_1 values.
**Scripts:** h7_impossibility.py, alpha1_gaps.py, alpha1_gap3_proof.py, c3_forces_c5.py, redei_converse_fast.py
**Theorem:** THM-029

### INV-068: DRT non-uniqueness and Paley dominance at n=11 — CORRECTED (MISTAKE-017)
**Source:** kind-pasteur-2026-03-06-S21 (computation), CORRECTED kind-pasteur-2026-03-07-S39b
**Status:** CORRECTED — previous "non-Paley DRT" was INVALID
**What:** Previous claim of "2 DRT classes at n=11" was based on connection set {1,2,3,5,8} which is NOT a valid tournament (S ∩ (-S) = {3,8} ≠ ∅, creating bidirectional edges). ALL claims about "non-Paley DRT" cycle counts (c3=44, c5=407, H=69311) are INVALID (MISTAKE-017).
**Corrected facts:** The only valid tournament (11,5,2)-difference sets in Z_11 are {1,3,4,5,9} (QR) and {2,6,7,8,10} (NQR), which give isomorphic Paley tournaments. There is no non-Paley circulant DRT at n=11. Whether a non-circulant DRT exists at n=11 remains open.
**Paley T_11 correct data:** H=95095, c3=55, c5=594, c7=3960, c9=11055, c11=5505, |Aut|=55.
**Next step:** (1) Check literature for non-circulant DRT existence at n=11 (all groups of order 11 are Z_11, so non-circulant DRT must be non-Cayley). (2) Test DRT uniqueness at n=23 where multiple constructions are known.
**Scripts:** drt_n11_analysis.py (CONTAINS BUG — uses invalid connection set), drt_n11_verify.py (correction script)

### INV-069: Scalar M characterization — M=(H/n)*I ↔ H-maximizer at n=5, ↔ VT at n=7
**Source:** kind-pasteur-2026-03-06-S25c (T156), opus-2026-03-06-S26 (T158)
**Status:** COMPUTED, OPEN CONJECTURE
**What:** Transfer matrix M[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b) is scalar (= (H/n)*I) for certain tournaments. At n=5: ALL 64 H-maximizers (H=15) have M=3*I, including 40 non-regular (scores 1,2,2,2,3) with |Aut|=3 (Z/3Z). These non-VT tournaments still have uniform position distribution. At n=7: M scalar ⟺ VT (100 sampled regular, 0 counterexamples). All circulant tournaments have scalar M. The characterization seems to be "uniform position distribution" (every vertex appears equally often at each position in Ham paths), which is weaker than VT at n=5 but equivalent at n=7.
**Key results:**
- n=5: M scalar ⟹ H=15 (max), but NOT ⟹ VT or regular
- n=7: M scalar ⟺ VT (conjecture, verified 100 samples)
- All circulant ⟹ M scalar (verified n=3,5,7)
- Paley T_7: |Aut|=21 (Frobenius Z/7Z⋊Z/3Z), sigma(x)=2x mod 7 is extra aut
- opus finding: non-regular n=5 scalar M has Z/3Z aut, uniform endpoint counts, uniform 3-cycle participation
**NONHAM VANISHING (kind-pasteur S25c):** The split pair decomposition M[a,b] = HAM(a,b) + NONHAM(a,b) shows:
- NONHAM(a,b) = 0 trivially when T[a,b]=1 (all split pairs are Ham paths)
- NONHAM(a,b) = M[a,b] when T[a,b]=0 (junction always fails)
- So NONHAM=0 for all (a,b) ⟺ M[a,b]=0 whenever T[a,b]=0
- VERIFIED: NONHAM=0 for ALL position-uniform n=5 (64/64 exhaustive) and ALL circulant n=7 (8/8)
- NONHAM ≠ 0 for general (non-uniform) tournaments at n=3,4,5
**PROOF CHAIN (verified n=3,5,7):**
1. Position-uniform ⟹ NONHAM=0 ⟹ M[a,b]=0 when T[a,b]=0
2. THM-030: M[a,b]=M[b,a]
3. For T[a,b]=1: T[b,a]=0 ⟹ M[b,a]=0 ⟹ M[a,b]=0 by symmetry
4. M[a,a]=H/n from position uniformity
5. M=(H/n)*I. QED (modulo proving step 1 at general n)
**Algebraic proof of step 1:** OPEN. The cancellation mechanism: for uniform T, nonzero E*B products at adjacent subset sizes pair up and cancel. For non-uniform T, "orphan" terms remain.
**Connection to Björklund et al.:** M[a,b] is a subset convolution in their "Fourier meets Möbius" framework.
**Next step:** (1) Prove NONHAM=0 for position-uniform tournaments at general n. (2) Test at n=9 (circulant only, n=9 exhaustive infeasible). (3) Find algebraic proof of the pairing mechanism.
**Scripts:** `04-computation/transfer_matrix_maximizers.py`, `04-computation/circulant_scalar_m_proof.py`, `04-computation/scalar_m_characterization.py`, `04-computation/scalar_m_n5_analysis.py`, `04-computation/circulant_m_scalar_proof.py`, `04-computation/nonham_vanish_general.py`, `04-computation/nonham_vanish_uniform.py`, `04-computation/nonham_proof_analysis.py`

### INV-070: Fibonacci determinant formula for transitive tournament
**Source:** kind-pasteur-2026-03-06-S25c (T157)
**Status:** VERIFIED n=2,...,11
**What:** For the transitive tournament T_n (i beats j iff i<j), det(M) = (-1)^{n(n-1)/2} * F(n+1) where F is Fibonacci. The matrix D*M = I+U-L (tridiagonal: 1s on diagonal/superdiagonal, -1s on subdiagonal) satisfies the Fibonacci recurrence. The Chebyshev eigenvalue conjecture (eigenvalues = 2cos(kπ/(n+1))) is FALSE.
**Connection:** The transitive tournament is the unique acyclic tournament, so Omega(T_n) is empty and I(Omega,2)=1=H(T_n). The transfer matrix has a clean tridiagonal structure reflecting the total order.
**Next step:** (1) Does the Fibonacci structure extend to near-transitive tournaments? (2) What is det(M) for other named tournament families?
**Scripts:** `04-computation/fibonacci_determinant_proof.py`

### INV-071: det(M) as tournament invariant — exhaustive n=5
**Source:** kind-pasteur-2026-03-06-S25c (T157)
**Status:** COMPUTED
**What:** At n=5, det(M) takes values {-27,-9,1,8,9,16,17,243} across 10 isomorphism classes (up to converse). The 10 eigenvalue patterns of M classify tournaments into exactly the 10 isomorphism classes. det(M) is NOT determined by (H,c3) — e.g., H=5,c3=2 gives both det=9 and det=17.
**Key finding:** det(M) is not related to det(A), per(A), or any simple matrix expression of the adjacency matrix A. Tested det(I±S), det(A+A^T), det(A*A^T), det(I+A*A^T) — none match. The transfer matrix is fundamentally different from the adjacency matrix.
**opus finding:** char poly det(λI-M) encodes hierarchy of path correlations. Scalar M ⟹ char poly = (λ-H/n)^n. PSD threshold at n=5: H≥13 ⟹ M is PSD.
**Next step:** (1) Compute det(M) exhaustively at n=7 (sample). (2) Find if any spectral graph invariant matches. (3) Investigate PSD threshold at general odd n.
**Scripts:** `04-computation/det_m_general_formula.py`, `04-computation/det_m_vs_adjacency_spectrum.py`, `04-computation/char_poly_M_analysis.py`, `04-computation/pos_skeleton_connection.py`

### INV-072: IO walk GF vs transfer matrix bridge
**Source:** opus-2026-03-06-S26, kind-pasteur-2026-03-06-S25c
**Status:** STRUCTURAL COMPARISON, NO BRIDGE FOUND
**What:** Irving-Omar walk GF W_D(z) = det(I+zXA^T)/det(I-zXA) uses cycle covers (det/per), while M[a,b] uses path decomposition (E_a*B_b). M[a,b] ≠ H(a→b) (direct endpoint-conditioned count). W(-z,-r) = W(z,r) at commutative level (opus finding). No simple matrix expression of A gives M.
**Key structural difference:** IO is multiplicative (det/per = products over cycles), M is additive (sum over subsets with inclusion-exclusion signs). Bridge might exist through Hopf algebra coproduct structure.
**Next step:** (1) Express M[a,b] in terms of IO's det/per framework. (2) Check if deletion-contraction on M matches Mitrovic's W_X = W_{X\e} - W_{X/e}^up.

### INV-073: Palindromic N => Scalar M for circulant tournaments — PROVED (THM-052)
**Source:** kind-pasteur-2026-03-06-S25d, building on THM-050 (opus-S26), THM-051 (opus-S26)
**Status:** PROVED for all circulant tournaments at odd n.
**What:** For circulant T on Z/nZ, the consecutive-position count N(a,b,j) = f(b-a mod n, j) is palindromic: f(d,j) = f(d,n-2-j). The proof uses three ingredients: (1) translation symmetry gives f depends only on d, (2) N symmetry gives f(d)=f(n-d), (3) self-complementarity via sigma: i->-i gives f(d,j)=f(n-d,n-2-j). Combining (2)+(3): f(d,j)=f(d,n-2-j). At odd n, palindromic N forces alternating sum = 0, so M[a,b]=0 for a!=b. Combined with M[a,a]=H/n, gives M=(H/n)*I.
**Verification:** n=5 (64/64 exhaustive pos-uniform), n=7 (8/8 circulant), n=9 (16/16), n=11 Paley, n=13 Paley.
**Key finding:** ALL position-uniform n=5 tournaments are self-complementary (64/64). No position-uniform tournaments exist at even n (n=4,6).
**Open extension:** Prove for general vertex-transitive (non-circulant) tournaments. The proof uses circulant-specific translation symmetry. At n=15, non-circulant VT tournaments exist (Babai-Kantor doubly-regular tournaments).
**Scripts:** `04-computation/palindromic_N_proof.py`, `04-computation/palindromic_N_posuniform.py`, `04-computation/palindromic_N_n9.py`, `04-computation/palindromic_N_n11.py`, `04-computation/selfcomp_posuniform_n7.py`

### INV-074: Diagonal signed position theorem — VERIFIED n=5
**Source:** opus-2026-03-06-S11b (continued³)
**Status:** VERIFIED computationally at n=5 (all 12 iso classes).
**What:** M[v,v] = sum_P (-1)^{pos(v,P)} where the sum is over all Hamiltonian paths P of T and pos(v,P) is the 0-indexed position of v in P. This means M[v,v] counts even-position appearances minus odd-position appearances. M[v,v] can be NEGATIVE (not a path count). For VT tournaments at odd n: M[v,v] = H/n. "Defect vertex" = vertex whose position distribution is biased relative to average.
**Connection to THM-027:** The trace formula tr(M) = H(T) = sum_v sum_P (-1)^{pos(v,P)} = sum_P sum_v (-1)^{pos(v,P)} = sum_P 1 (at odd n, since alternating sum of 0..n-1 positions = 1). This reproves the trace formula.
**Next step:** Prove from IE formula definition. Connect to defect vertex characterization.
**Scripts:** `04-computation/diagonal_signed_position_theorem.py`

### INV-075: Perpendicularity of M-directions across H-classes — CONFIRMED n=7
**Source:** opus-2026-03-06-S11b (continued³)
**Status:** CONFIRMED computationally at n=7 (790 iso classes).
**What:** The non-scalar part of M (i.e., M - (H/n)*I) has a "direction" in matrix space. Measuring cosine similarity between these directions across iso classes shows:
  - Low H (H<85): positive cosine (~0.5-0.8, aligned)
  - Mid H (H≈85-105): near zero (perpendicular!)
  - High H (H>105): negative cosine (anti-aligned)
  - Overall mean cosine = -0.0485 (near perpendicular)
The crossover (true perpendicularity) occurs near the MEDIAN H value. This means the non-scalar perturbation "rotates" continuously through eigenspace as H varies.
**Connection:** This is the "perpendicularity" the user hypothesized in earlier sessions. It connects to the inverted-U of position variance and the grid symmetry structure.
**Next step:** (1) Prove analytically why crossover occurs at median H. (2) Check at n=9. (3) Connect to the even cycle vanishing theorem (INV-053) which uses the same T↔T^op involution.
**Scripts:** `04-computation/perpendicularity_cosine_n7.py`

### INV-076: All H-maximizers have scalar M — VERIFIED n=5,7
**Source:** opus-2026-03-06-S11b (continued³)
**Status:** VERIFIED exhaustively at n=5 (both max-H classes) and n=7 (all 43 maximizers).
**What:** Every tournament achieving max H has M = (H/n)*I (scalar transfer matrix). At n=7: all 43 maximizers have M = 27*I. At n=5: both H=15 classes (VT circulant and non-regular VT) have M = 3*I.
**Conjecture:** For ALL odd n, the H-maximizer has scalar M. This is equivalent to saying H-maximizers are always vertex-transitive (or at least "position-uniform").
**Connection:** Combines with THM-052 (circulant => scalar) and the Paley maximizer conjecture. If Paley maximizes H and Paley is circulant, then Paley gives scalar M. The deeper question is whether scalar M is NECESSARY for H-maximization.
**Next step:** Verify at n=9 (need to find maximizers first).

### INV-077: VT tournament NOT self-converse at n=21 — THM-052 DISPROVED for non-SC VT
**Source:** kind-pasteur-2026-03-06-S25e
**Status:** RESOLVED. M is NOT scalar for non-SC VT tournaments.
**What:** ALL 22 non-circulant VT tournaments at n=21 (from McKay's database) are NOT self-converse. These are Cayley tournaments on F_21 = Z/7 x| Z/3. All 88 circulant VT tournaments at n=21 ARE self-converse.
**Computation (n=21, 1075s):**
- H(T) = 123,522,430,238,361 (divisible by 21)
- N(0,1,j) is NOT palindromic: N[0]=581,223,220,317 vs N[19]=581,314,958,778
- Alternating sum = M[0,1] = 45,478,409 != 0
- Therefore M != (H/n)*I for this VT tournament
**Conclusion:** THM-052 is PROVED for self-converse VT (including all circulants) but DISPROVED for non-SC VT. Self-converse is the exact boundary.
**Scripts:** `04-computation/frobenius21_palindromic_N.py`, `04-computation/mcKay_vt21_selfconverse.py`

### INV-078: Aut(T) union Anti(T) transitivity characterizes scalar M
**Source:** opus-2026-03-06-S26 (scalar_m_aut_anti_characterization.py)
**Status:** VERIFIED at n=5 (exhaustive). CONFIRMED at n=21: F_21 non-SC has Anti=empty, M not scalar.
**What:** Scalar M (M = (H/n)*I) holds iff Aut(T) union Anti(T) acts transitively on V. For the F_21 non-normal tournament: Aut=F_21 (transitive) but Anti=empty, so Aut+Anti = Aut alone. The conjecture predicts M not scalar, which is CONFIRMED by computation.
**Next step:** Verify at n=7 (exhaustive).

### INV-079: W(r) coefficient stratification by odd-cycle complexity
**Source:** kind-pasteur-2026-03-06-S25f, opus-2026-03-06-S27 (THM-055)
**Status:** PROVED (k=0,1), VERIFIED (k=2). Connected to Hopf algebra coproduct.
**What:** W(r) = sum_P prod(r + s_e) has coefficients w_{n-1-2k} = sum_P e_{2k}(s_P).
  - w_{n-1} = n! (universal)
  - w_{n-3} = 2*(n-2)!*t_3 - const (depends on t_3 only — PROVED)
  - w_{n-5} depends on t_3 AND the 4th moment of f_P (opus THM-055)
  - At n=5: w_0 = -t_3 + 2*t_5 + 1 (EXACT, kind-pasteur verified exhaustive)
**Key identity:** H = 1 + 2*(t_3 + t_5) at n=5 (OCF simplification since a_2=0)
**Recursive Hopf structure:** overlap=3 contribution to w_{n-5} at n=7 uses OCF at n=5 on each 5-element sub-tournament. This is the Hopf algebra coproduct Delta([T]) evaluated on fibers.
**Connection to THM-055:** e_{2k}(s_P) is a polynomial of degree 2k in f_P. All power sums p_{2l} are constant; only p_1 = f - (n-1)/2 varies. So everything reduces to moments of f_P.
**Next step:** (1) Find explicit formula for w_0 at n=7 in terms of tournament invariants. (2) Prove the Hopf algebra recursion algebraically. (3) Determine whether the 4th moment has a cycle-theoretic interpretation.

### INV-080: Pfaffian-path duality at even/odd n
**Source:** kind-pasteur-2026-03-06-S25f (pfaffian_path_duality.py)
**Status:** COMPUTATIONAL. Interesting correlations found.
**What:** At even n: det(S) = Pf(S)^2 is a nonzero odd square. At odd n: det(S) = 0 but tr(M) = H > 0.
  - n=4: det(S) is EXACTLY determined by t_3: det=1 if t_3 even, det=9 if t_3 odd
  - n=6: det(S) is NOT determined by t_3 alone (needs finer invariants)
  - |Pf| always odd (Fisher-Ryan), values in {1,3,5,...,n-1} at even n
  - The Pfaffian and H are NOT functionally related but are both determined by S
**Connections:** D_k classification (Zeng-You-Zhao 2025), Seidel tournament matrices (determinant gap phenomena)
**Speculative:** Is there a formal duality: paths (odd n) ↔ cycle covers (even n)?
**Next step:** (1) Check if det(S) at n=6 is determined by (t_3, t_5). (2) Investigate the eigenvalue connection between Pf and H.

### INV-081: Paley tournament W(r) structure
**Source:** kind-pasteur-2026-03-06-S25f (eigenvalue_W_connection.py)
**Status:** COMPUTED at p=7.
**What:** Paley T_7 has W(r)/7! = [1/320, 0, 1/80, 0, 1/4, 0, 1].
  - Eigenvalues of A_Paley: (p-1)/2 = 3 (mult 1), (-1±sqrt(p))/2 (mult (p-1)/2 each)
  - All non-trivial eigenvalues degenerate → scalar M → W coefficients have maximal symmetry
  - Paley requires p ≡ 3 (mod 4) (so -1 is not a QR)
**Next step:** Compute W(r) for T_11 and T_3. Check if W/p! ratios have a closed form involving p.

### INV-082: EXACT W-coefficient hierarchy — spectral decomposition of tournaments
**Source:** kind-pasteur-S25g
**Status:** VERIFIED at n=5 (exhaustive) and n=7 (20 random samples, 0 error)
**What:** W(r) coefficients form a hierarchy: w_{n-1-2k} depends on cycle data of complexity k.
  - w_{n-1} = n! (universal)
  - w_{n-3} = (n-2)! * [2*t_3 - C(n,3)/2] (depends on t_3 only; CENTERED: zero mean over random T)
  - w_{n-5} at n=7: -60*t_3 + 12*t_5 + 24*alpha_2 + 231
  - w_0 at n=7: 2*t_3 - t_5 + 2*t_7 - 2*alpha_2 - 17/4
  - Each level adds ONE new cycle complexity (t_{2k+1} or alpha_k)
  - H - w_0 penalty SHIFTS to higher-order cycles at larger n
  - n=5: H - w_0 = 3*t_3; n=7: H - w_0 = 3*t_5 + 6*alpha_2 + 21/4
  - Analogous to spectral decomposition / renormalization group / characteristic classes
**Next step:** Verify w_{n-3} = (n-2)! * [2*t_3 - C(n,3)/2] at n=9. Compute full hierarchy at n=9.

### INV-083: Rooted tournaments = OEIS A093934
**Source:** kind-pasteur-S25g
**Status:** VERIFIED n=2 through 6
**What:** P(n) = sum over iso classes of (# vertex orbits) = # rooted tournament iso classes.
  - P(2)=2, P(3)=4, P(4)=12, P(5)=48, P(6)=296
  - Matches OEIS A093934 (with offset)
  - P(n) = 2*(n-1)! for n<=5 (coincidence), FAILS at n=6
  - Orbit distributions: n=5 {1:1,3:4,5:7}, n=6 {2:5,4:10,6:41}
**Next step:** Check A093934 description more carefully. Compute P(7) if feasible.

### INV-084: W-coefficient hierarchy as Mobius inversion on cycle complex
**Source:** kind-pasteur-S25g (creative synthesis)
**Status:** CONJECTURED
**What:** The W-coefficients can be viewed as evaluations of I(Omega, x) at different points.
  - H = I(Omega, 2) (OCF, x=2)
  - chi = I(Omega, -1) (Euler characteristic, x=-1)
  - w_0 is an "intermediate" evaluation (not simply I(Omega, c) for some c)
  - The hierarchy parallels Fourier decomposition, but REVERSED: high "frequencies" are simple
  - Connects to renormalization: each level "integrates out" one cycle scale
**Next step:** Determine if w_0 = I(Omega, c) for some specific c, or if it's genuinely different.

### INV-085: Bipartite skeleton and t3 parity (THM-060)
**Source:** kind-pasteur-2026-03-06-S25h
**Status:** PROVED at n=3,5,7,9. Structural argument + exhaustive verification.
**What:** Blue line skeleton (GS flip graph on SC classes) is bipartite at odd n, with bipartition determined by t3 parity. At even n, skeleton is NOT bipartite (has 3-cycles).
**Key mechanism:** Consecutive triples each contribute 1 to t3(T)+t3(flip(T)), total n-2 (odd at odd n). Non-consecutive triples contribute even total for GS tilings.
**Open:** Algebraic proof of Type B evenness; spectral structure at large n; connection between skeleton eigenvalues and tournament invariants.
**Scripts:** 04-computation/t3_parity_proof_complete.py, bipartition_invariant.py, bipartition_n7_verify.py
**Writeup:** 01-canon/theorems/THM-060-bipartite-skeleton.md, 03-artifacts/drafts/bipartite-skeleton-synthesis-S25h.md

### INV-086: Silver ratio in skeleton eigenvalues
**Source:** kind-pasteur-2026-03-06-S25h, skeleton_spectral.py
**Status:** OBSERVED at n=5. Eigenvalues {±(1+√2), ±1, ±1, ±(√2-1)}.
**What:** The skeleton adjacency matrix at n=5 has eigenvalues involving the silver ratio 1+√2. K^2 diagonal = GS class sizes. Is this coincidence or does it persist at n=7?
**Next step:** Compute skeleton eigenvalues at n=7 (88×88 matrix). Check if silver ratio generalizes.

### INV-087: Antiferromagnetic interpretation of skeleton
**Source:** kind-pasteur-2026-03-06-S25h
**Status:** CONCEPTUAL. Skeleton = Ising model with antiferromagnetic coupling.
**What:** SC classes have "spin" = (-1)^{t3}. GS flip edges connect opposite spins. At odd n: perfect Neel order (unfrustrated). At even n: frustrated (odd cycles). Connects tournament theory to statistical mechanics.
**Next step:** Compute partition function Z(beta) = sum over SC of H(T)^beta. Check for phase transitions.

### INV-088: Schweser-Stiebitz-Toft — Rédei revisited (expository)
**Source:** arXiv:2510.10659 (Oct 2025, revised Feb 2026), found by web search opus-2026-03-07-S36
**Status:** CATALOGED. Expository paper, likely low priority.
**What:** Revisits the classical theorems of Rédei, Dirac, and Berge on Hamiltonian paths in tournaments. Exhibits the stronger theorems and explains connections between them. Does NOT mention independence polynomials, odd cycles, or conflict graphs.
**Next step:** Skim for any novel structural insight about H(T) parity not already in our framework. Low priority.

### INV-089: Irving-Omar authorship correction
**Source:** opus-2026-03-07-S37 (this session)
**Status:** CORRECTED in THM-002, CONJ-001, THM-070.
**What:** arXiv:2412.10572 ("Revisiting The Rédei-Berge Symmetric Functions via Matrix Algebra") is by **Irving & Omar**, NOT Grinberg & Stanley. Their Corollary 20 restates Grinberg-Stanley's Theorem 1.39 + Lemma 6.5 from arXiv:2307.05569. The OCF result itself is correctly attributed to Grinberg-Stanley; only the paper authorship was wrong.
**Remaining:** Some computation scripts and broadcast messages still reference "Grinberg-Stanley" for arXiv:2412.10572. These are historical and low priority to fix.

### INV-090: Three equivalent H(T) formulas and the even-cycle cancellation
**Source:** opus-2026-03-07-S37
**Status:** VERIFIED computationally. Structural explanation needed.
**What:** Three equivalent ways to compute H(T):
1. **Direct**: count Hamiltonian paths (Held-Karp)
2. **OCF = I(Ω(T), 2)**: sum over T-only disjoint odd cycle collections, weight 2^k
3. **ps_1(U_T)(1)**: sum of ALL Rédei-Berge p-coefficients (uses T + T^op cycles, signed by (-1)^φ)

KEY FINDINGS at n=7 (Paley):
- All p-coefficients of U_T have ALL-ODD partitions (no even-part lambdas appear!)
- Sum of coefficients = 189 = H(T) [ps_1 at m=1]
- OCF specialization of FULL U_T = 433 ≠ H(T)
- OCF specialization of T-ONLY part = 189 = H(T)
- The "OCF specialization" (p_1→1, p_odd→2, p_even→0) is NOT how GS proves OCF
- GS uses ps_1(U_T)(1) which uses ALL cycles, not just T-direction ones

OPEN: Why are the p-coefficients supported only on all-odd partitions at n=7? Is this true for all tournaments? Omega constraint: (-1)^{n-l(λ)} symmetry forces sum over even-l terms = 0, but here they're individually zero. Is the all-odd support a coincidence of n=7 (Mersenne prime) or universal?

**Mixed-direction findings at n=9:** Mixed T/T^op cycle pairs DO exist at n=9 (100+ per tournament). So U_T at n=9 WILL have even-part lambda contributions, but they must cancel in ps_1(1).
**RESOLVED (opus-S38 agent):** The all-odd support IS universal for tournaments (proved by Grinberg-Stanley Theorem 1.39). The factor 2^k in OCF is the **orientation multiplicity**: each odd cycle has exactly 2 directed orientations (T and T^op), both contributing sign +1 (since (-1)^{k-1}=+1 for odd k). An independent set of k cycles thus contributes 2^k copies. Even-part lambdas vanish because even cycles contribute opposite signs from T and T^op directions.
**Scripts:** `04-computation/omega_constraint.py`, `mixed_sum_n9.py`, `ut_specialization_n9.py`

### INV-091: H=21 permanent gap — PROVED through n=7
**Source:** opus-2026-03-07-S38, kind-pasteur-2026-03-07-S28
**Status:** PROVED through n=7 (exhaustive). Strong conjecture for all n.
**What:** H(T)=21 is never achieved by any tournament on n<=7 vertices. Exhaustive computation: 2,097,152 tournaments at n=7, zero with H=21. The gap 19→23 appears at both n=6 and n=7. OCF analysis: none of the valid (alpha_1,alpha_2) decompositions for H=21 are achievable.
**Key structural insight:** At n=6, the achievable alpha_1 values jump in a way that skips all decompositions summing to H=21. The constraint is that certain cycle counts force additional cycles or disjoint pairs, pushing H past 21.
**Scripts:** `04-computation/h21_n7_fast.py`, `04-computation/h21_theory_fixed.py`
**Theorem:** THM-075
**Next step:** (1) Prove H=21 impossibility at general n. (2) Characterize ALL permanent gaps. (3) Is {7, 21} the complete list, or are there more?

### INV-092: Type-count sequence = A000009 (partitions into distinct parts)
**Source:** opus-2026-03-07-S38 agent, kind-pasteur-S29
**Status:** PROVED (bijective argument).
**What:** The number of OCF cycle types at size n (multisets of odd parts >=3 summing to <=n, plus empty) equals A000009(n) = number of partitions into distinct parts. Bijection: remove all 1's from partition into odd parts. Null-dim sequence 0,0,1,3,6,11,19,29,44,65 is NOT in OEIS (likely novel).

### INV-093: Tangent number connection proved
**Source:** opus-2026-03-07-S38 agent, kind-pasteur-S29
**Status:** PROVED.
**What:** P_n(0,0) = 2^{(n-1)/2} * T_n where T_n is the n-th tangent number. Proof: P_n(u,0) = A_n(t)/t^{(n-1)/2} with u=t+1/t. Setting u=0 gives t=i, and A_n(i)/i^{(n-1)/2} = 2^{(n-1)/2} * T_n by the EGF of Eulerian polynomials evaluated at t=i.
**Connection:** Hetyei (2017, arXiv:1704.07245) on "alternation acyclic tournaments" connects tournament counts to median Genocchi numbers (same family as tangent numbers).

### INV-094: Mitrovic noncommutative deletion-contraction explored
**Source:** opus-2026-03-07-S38, arXiv:2504.20968 (Mitrovic, Apr 2025)
**Status:** EXPLORED. Not useful for OCF.
**What:** W_X = W_{X\e} - W_{X/e}^up for the noncommutative Redei-Berge function. Edge contraction T/e is NOT a tournament (bidirectional edges possible). OCF fails for T\e and T/e. The deletion-contraction approach is useful for algebraic properties of U_T but not for proving OCF.
**Scripts:** `04-computation/edge_deletion_contraction.py`

### INV-095: Bags-of-sticks for OCF — DEAD END
**Source:** opus-2026-03-07-S38 agent
**Status:** CLOSED (dead end).
**What:** The bags-of-sticks decomposition (Mitrovic-Stojadinovic Theorem 4.2) expresses U_X via inclusion-exclusion over edge deletions. Under OCF specialization, every bag of sticks contributes 1 (acyclic digraphs have empty Omega). So the decomposition reduces to: H(T) = sum of inclusion-exclusion coefficients, which is trivially true. No new information for OCF.

### INV-096: H=21 Component Reduction (THM-079) — PROVED FOR ALL n
**Source:** opus-2026-03-07-S39 (partial), kind-pasteur-2026-03-07-S33 (completion)
**Status:** PROVED. H(T) ≠ 21 for ALL tournaments on ALL n vertices.
**What:** For H(T)=21, the OCF requires I(Omega(T),2)=21. Component factorization gives:
  - Disconnected case: IMPOSSIBLE. 21=3*7, but I(component)=7 impossible by THM-029 argument.
  - P_4 case: IMPOSSIBLE. P_4 realization blocked because sharing 3-cycles force extra cycles.
  - K_6-2e case: SUPERSEDED by Dichotomy Theorem proof.
**Dichotomy Theorem (Part R):** For cycle-rich T on n≥9, either (a) 3 disjoint 3-cycles exist (⟹ H≥27), or (b) safe deletion to cycle-rich T−v exists (⟹ H≥H(T−v)+2≥27). Proved via poisoning graph DAG argument. Combined with base case n≤8 (exhaustive) and Part J (non-cyclic vertex removal), gives H(T)≠21 for ALL n.
**Key ingredients:** Lemma Q (cycle-rich ⟹ no source/sink), poisoning graph has out-degree ≤1 and is acyclic, DAG source deletion preserves cycle-rich.
**Scripts:** `04-computation/h21_gap_mechanism.py`, `h21_dichotomy_proof.py`, `h21_poisoning_graph.py`, `h21_cycle_rich_auto_no_ss.py`
**Writeup:** `03-artifacts/drafts/dichotomy-proof-formal.md`

### INV-097: u_T Size-Weighted Independence Polynomial (THM-078) — PROVED
**Source:** opus-2026-03-07-S39
**Status:** PROVED. u_T(m) = sum_j sw(j)*m^{n-2j} where sw(j) = sum over j-element independent sets of 2^|S|.
**What:** The size-weighted independence polynomial identity connects u_T(m) to the Omega(T) independence structure. Q_T(w) = u_T(sqrt(w))/sqrt(w) has all real non-positive roots for n<=8 (Leake-Ryder/Chudnovsky-Seymour stability for claw-free graphs). Fails at n>=9 when claws appear.
**Next step:** Check if Q_T root structure has implications for achievable H values.

### INV-098: Lichiardopol's Conjecture and Disjoint Cycle Forcing — EXPLORED
**Source:** kind-pasteur-2026-03-07-S33 (web research)
**Status:** EXPLORED. Used in H=21 proof context but not directly needed.
**What:** Lichiardopol's conjecture (proved for q=3 by Bang-Jensen, Bessy, Thomassé): tournaments with min out-degree ≥ (q-1)k-1 have k vertex-disjoint q-cycles. For 3-cycles with k=3: min outdeg ≥ 5. However, at n=9 cycle-rich, min outdeg is ALWAYS ≤ 4 (100% of 106,424 tested), so Lichiardopol doesn't directly fire. The poisoning graph argument covers ALL cases including those below the Lichiardopol threshold.
**Papers:** [Bang-Jensen-Bessy-Thomassé](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v27i2p52)
**Next step:** Could be useful for other permanent H-gap proofs where k≥4 disjoint cycles are needed.

### INV-099: Chen-Chang 2024 Disjoint Cycles in Tournaments — CATALOGED
**Source:** kind-pasteur-2026-03-07-S33 (web research)
**Status:** CATALOGED. Not yet deeply investigated.
**What:** Chen-Chang (2024, J. Graph Theory) prove results on disjoint cycles in tournaments. Extends cycle-matching theory. Could provide tools for proving existence of k≥4 disjoint 3-cycles under weaker conditions than Lichiardopol.
**Paper:** [Chen-Chang 2024](https://onlinelibrary.wiley.com/doi/10.1002/jgt.23038)
**Next step:** Read the paper for potentially stronger theorems applicable to H-gap proofs.

### INV-100: Frankl's Proof of Erdos Matching Conjecture (k=3) — CATALOGED
**Source:** kind-pasteur-2026-03-07-S33 (web research)
**Status:** CATALOGED. Provides context for 3-uniform hypergraph matching.
**What:** Frankl proved the Erdos matching conjecture for k=3: bounds the max number of 3-element sets with no matching of size s+1. The cycle vertex sets in Omega(T) form a 3-uniform hypergraph (for 3-cycles). Frankl's bound could constrain the maximum number of 3-cycles with bounded matching number.
**Connection:** If mm(T) ≤ 2, Frankl's bound limits |Omega_3(T)| ≤ max(C(5,3), 3*3-3+1) depending on exact formulation. This could give an independent route to the dichotomy.
**Next step:** Check exact Frankl bound for our setting (n vertices, 3-uniform, matching ≤ 2).

### INV-101: Other Permanent H-Gaps Beyond 7 and 21 — CONFIRMED THROUGH n=8 EXHAUSTIVE
**Source:** kind-pasteur-2026-03-07-S33, opus-2026-03-07-S43
**Status:** STRONG CONJECTURE that H=7 and H=21 are the ONLY permanent gaps.
**What:** With H=7 and H=21 both proved as permanent gaps (never achieved for ANY n), the natural question is: are there other permanent gaps?
**Computational evidence (opus-S43):**
  - ALL n=7 gaps (63, 107, 119, 149, 161-169, 173) fill at n=8 (sampling, very quickly)
  - n=8 exhaustive computation running (268M tournaments)
  - For H≥27 (w≥13): 20+ graph-feasible decompositions, blocking all seems impossible
  - Decomposition analysis: for w≥4, the (w-2,1,0,...) decomposition is available; it fails at w=10 due to cascade forcing (Part N) but works at all other w
**Algebraic argument (opus-S43):**
  - For w≥13: alpha_3≥1 decompositions become available (3 disjoint cycles feasible)
  - The number of decompositions grows rapidly with w (14 at w=10, 20 at w=13, 60 at w=20)
  - Each decomposition needs an independent tournament obstruction to block it
  - Only w=3 (1 feasible decomp) and w=10 (4 feasible decomps, all blocked) have this property
**Mod-4 result (Grinberg-Stanley Theorem 7.1):** H(T) ≡ 1 + 2·(# nontrivial odd cycles) mod 4. Does not directly rule out any odd H.
**Conjecture: H=7 and H=21 are the ONLY permanent gaps in the H-spectrum.**
**n=8 EXHAUSTIVE RESULT (opus-S45):** All 268,435,456 tournaments enumerated. Max H=661. Only missing odd values in [1,300]: H=7 and H=21. This CONFIRMS the conjecture through n=8. No new gaps appear. All n=7 gaps (63, 107, 119, etc.) fill at n=8.
**Status:** STRONG EVIDENCE. Conjecture holds through n=8 exhaustive enumeration.

### INV-102: Grinberg-Stanley Mod-4 Theorem (Theorem 7.1) — CATALOGED
**Source:** opus-2026-03-07-S43 (web research)
**Status:** CATALOGED. Read and extracted from arXiv:2307.05569.
**What:** Theorem 7.1: H(T) ≡ 1 + 2·(# nontrivial odd D-cycles) (mod 4). This is the OCF mod-4 reduction: since H = 1 + 2·alpha_1 + 4·(...), H mod 4 = 1 + 2·alpha_1 mod 4. The proof uses the power-sum expansion and the specialization map zeta. Not directly useful for gaps but confirms the algebraic structure.
**Next step:** Check if higher modular refinements (mod 8, mod 16) exist in the Grinberg-Stanley framework.

### INV-103: Non-Separating Vertices in Tournaments — CATALOGED
**Source:** opus-2026-03-07-S43 (web research)
**Status:** CATALOGED. Related but weaker than cycle-rich deletion.
**What:** A vertex in a strongly connected tournament is "non-separating" if its removal preserves strong connectivity. For min in/out-degree ≥ p, at least min{|V|, 4p-2} non-separating vertices exist. Our "good deletion" requirement is stronger: preserve cycle-richness (every vertex in a 3-cycle), not just strong connectivity.
**Next step:** Could the non-separating vertex techniques be adapted to our stronger requirement?

### INV-104: "Cycle-Rich" as Novel Concept — NOTED
**Source:** opus-2026-03-07-S43 (web research)
**Status:** The term "cycle-rich" (every vertex in a directed 3-cycle, no source/sink) does not appear in the literature. This is a novel concept from our project. The poisoning graph argument (Part R) may be publishable as a standalone result about cycle-rich tournaments.

### INV-105: Deletion-Contraction for H(T) — VERIFIED
**Source:** kind-pasteur-2026-03-07-S34
**Status:** VERIFIED COMPUTATIONALLY. H(T) = H(T\e) + H(T/e) holds 100% at n=4,5. Commutative specialization of Mitrovic's W_X = W_{X\e} - W_{X/e}↑ (arXiv:2504.20968).
**Convention:** Contraction merges tail/head: w inherits IN from tail, OUT from head.
**Next step:** Prove algebraically (should follow from Mitrovic by specialization). Use for inductive proof of Redei/OCF/H-gaps.

### INV-106: GLMY Path Homology of Tournaments — IDENTIFIED
**Source:** kind-pasteur-2026-03-07-S34 (web research)
**Status:** IDENTIFIED. Tang-Yau (2026, arXiv:2602.04140) compute path homology of circulant digraphs using Fourier decomposition — directly applicable to Paley tournaments.
**Key connection:** H_1(T) should relate to cycle space of T and thus to Omega(T) and alpha_1.
**Next step:** Compute GLMY path homology for small tournaments (n=4,5). Check if rank(H_1) = alpha_1 or c_3.

### INV-107: Extended Root Polytope Deletion-Contraction — IDENTIFIED
**Source:** kind-pasteur-2026-03-07-S34 (web research). Kalman-Tothmeresz 2024, arXiv:2409.18902.
**Status:** IDENTIFIED. h*-polynomial of extended root polytope is monotone under deletion-contraction.
**Next step:** Compute root polytope for small tournament digraphs. Check if h* relates to I(Omega,x).

### INV-108: Lee-Yang Zeros of F(T,x) — DISCOVERED (opus-S44)
**Source:** opus-2026-03-07-S44
**Status:** MAJOR DISCOVERY. F(T,x) zeros come in reciprocal pairs (palindrome). Cluster at ±2pi/3 on unit circle. H=9 at n=5: ALL zeros on unit circle.
**Key findings:** F(T,omega) real at n=7, F(T,i) pure imaginary, universal divisibilities F(T,omega) = 0 mod 9 and F(T,i) = 0 mod 16i.
**Next step:** Prove Lee-Yang property for specific tournament classes. Connect to phase transitions in statistical mechanics.

### INV-109: Walsh/Fourier Spectral OCF — DISCOVERED (opus-S35c)
**Source:** opus-2026-03-07-S35c9
**Status:** MAJOR DISCOVERY. THM-081: hat{t_k}[S] = (1/2^k) sum (-1)^{asc(S,C)}. Counting identity provides new proof path for OCF.
**Next step:** Prove counting identity algebraically for d=1 (single edge). Extend to n=7 (need hat{t7} and hat{bc35}).

### INV-110: Ihara Zeta Function of Tournaments — TESTED
**Source:** kind-pasteur-2026-03-07-S34
**Status:** TESTED. z_inv(1/2) = det(I-A/2+(D-I)/4) strongly correlated with H (r=-0.95 at n=5) but NOT uniquely determined.
**Conclusion:** Ihara zeta constrains H but doesn't determine it. Consistent with cycles constraining but not determining independence structure.

### INV-111: p-adic Structure Beyond p=2 — TESTED
**Source:** kind-pasteur-2026-03-07-S34
**Status:** TESTED. H mod 3 = (1+2*alpha_1+alpha_2) mod 3 from OCF. At n=4: H mod 3 uniquely determined by c3 via (1+2c3) mod 3. H mod 7 = 0 impossible at n<=6, first achievable at n=7.
**Next step:** Investigate H mod p for larger p. Is there a p-adic tower for p=3?

### INV-112: Converse Invariant Digraph Polynomials — IDENTIFIED
**Source:** kind-pasteur-2026-03-07-S34 (web research). Ai-Gutin-Lei-Yeo-Zhou 2024, arXiv:2407.17051.
**Status:** IDENTIFIED. New digraph polynomial for converse invariance testing. H(T)=H(T^op) means Ham paths are converse invariant.
**Next step:** Read the paper. Check if their polynomial gives new information about H(T).

### INV-113: Stanley-Stembridge Resolution Implications — IDENTIFIED
**Source:** kind-pasteur-2026-03-07-S34 (web research). Hikita 2024.
**Status:** IDENTIFIED. Stanley-Stembridge (e-positivity of chromatic SF for 3+1-free posets) proved by Hikita 2024.
**Connection:** Via Mitrovic-Stojadinovic, Redei-Berge U_T connects to chromatic SF. If tournament poset is (3+1)-free, U_T inherits e-positivity.
**Next step:** Check if tournament arc ordering posets are (3+1)-free. Investigate Hessenberg varieties approach.

### INV-114: Flip Formula F(T,x) - F(T',x) = (x-1)*D(x) — PROVED (THM-083)
**Source:** opus-2026-03-07-S45 (computational discovery), kind-pasteur-2026-03-07-S35 (algebraic proof)
**Status:** PROVED algebraically (THM-083). Verified at n=4,5 exhaustive.
**What:** For arc u->v in T, flip gives T'. The difference F(T)-F(T') factors as (x-1)*D(x), where D = F(T/e) - F(T'/e') is anti-palindromic.
**Key identities (THM-083):**
  - F_T(x) = F_{T\e}(x) + (x-1) * F(T/e, x)  (polynomial deletion-contraction)
  - G_{u,v}(x) = F(T/e, x)  (contraction = conditional path polynomial)
  - D(x) = -x^{n-2} D(1/x)  (anti-palindromicity from tournament palindrome)
  - H(T) - H(T') = D_{n-2}  (leading coefficient of D)
**CORRECTION (kind-pasteur-S35):** H(T) ≠ H(T') under arc flip in general (deltas up to ±12 at n=5). The opus claim "H(T)=H(T')" was WRONG. The correct statement: F(T,1)=n!=F(T',1) (total permutation count).
**CORRECTION (kind-pasteur-S35):** G_uv + G_vu = 2*F(T/e) only when T/e is a tournament (requires u,v to have identical profiles to other vertices). In general, G_uv + G_vu = F(T/e) + F(T'/e') which is palindromic but ≠ 2*F(T/e).
**Scripts:** `04-computation/f_poly_flip_formula.py`, `04-computation/flip_formula_D_analysis.py`, `04-computation/poly_deletion_contraction.py`, `04-computation/flip_reduction_via_contraction.py`

### INV-115: Matroid Structure of Vertex-Disjoint Odd Cycles — BOUNDARY at n=5
**Source:** opus-2026-03-07-S45 (computational discovery)
**Status:** VERIFIED. Exchange axiom holds at n=5 (1024/1024) but FAILS at n>=6 (15360/32768 at n=6).
**What:** Collections of vertex-disjoint odd directed cycles in T form the independent sets of a matroid if and only if n<=5. At n>=6, maximal independent sets can have different sizes.
**Script:** `04-computation/gammoid_matroid_test.py`
**Next step:** (1) Prove n<=5 case (small enough for case analysis). (2) Characterize failure at n=6 — which exchange pairs fail? (3) Relationship to Omega(T) perfectness boundary (perfect through n=7, fails n=8).

### INV-116: Transfer Matrix W(x) and per(W) — EXPLORED
**Source:** opus-2026-03-07-S45
**Status:** EXPLORED. per(W(1)) = D_n (subfactorial) universally. F(T,x) = Hamiltonian path sum over W entries. per(W(x)) palindromic for certain tournament classes.
**Scripts:** `04-computation/transfer_matrix_F_connection.py`, `04-computation/per_W_analysis.py`
**Next step:** (1) Explore eigenvalues of W(x) at specific x values. (2) Connection to Irving-Omar det formula (INV-046). (3) Can det(I-zW) generating function extract F?

### INV-117: Archer-Gessel-Graves-Liang Strong Tournament Descent Polynomial — RESEARCHED
**Source:** opus-2026-03-07-S45 (background agent), Discrete Math 343 (2020)
**Status:** RESEARCHED. Paper fully reviewed. t_n(u) = descent poly for strong tournaments. Palindromic, divisible by (1+u)^{floor(n/2)}. GF: U(x) = 1/(1-T(x)).
**Connection:** Different statistic from F(T,x) (global descent vs path-local forward edges), but same palindromic structure. The (1+u)^{floor(n/2)} divisibility may connect to OCF factor-2 structures. The "Eulerian graphic GF" framework (q=(1+uy)/(1+y)) is the natural algebraic home for descent statistics on tournaments.
**Notes:** `04-computation/gessel_strong_tournament_notes.md`
**Next step:** (1) Compute t_n(u) at small n and compare to F(T,x) aggregates. (2) Check if the strong component decomposition gives new structural insights for H(T). (3) Test the (1+u)^{floor(n/2)} divisibility analogue for F(T,x).

### INV-118: F(T,omega) mod 9 Universality — CONFIRMED NOVEL
**Source:** opus-2026-03-07-S44/S45 (computational discovery + background agent literature search)
**Status:** CONFIRMED NOVEL. Extensive literature search found NO prior work on roots-of-unity evaluations of F(T,x). Not a consequence of Grinberg-Stanley mod-4. Chebikin et al. studied cyclotomic factors of descent set polynomial Q_n but Phi_3 doesn't appear for n<=23.
**What:** F(T,omega) ≡ 0 mod 9 at n=7 (all 5040 tournaments). Equivalently S_0 = sum_{k≡0 mod 3} F_k ≡ 0 mod 6.
**Next step:** (1) Prove algebraically using OCF or Fourier decomposition. (2) Check at n=9,10. (3) Generalize: are there universal congruences for F(T,zeta_k) mod k^2?

### INV-119: Deletion-Contraction for Hamiltonian Paths — PROVED (THM-082)
**Source:** kind-pasteur-2026-03-07-S35
**Status:** PROVED by clean bijection argument. Verified exhaustive n=4,5.
**What:** For any digraph D with directed edge e=(u→v):
  H(D) = H(D\e) + H(D/e)
where D\e = deletion, D/e = contraction (w inherits IN from u, OUT from v).
**Proof:** Ham paths not using e = H(D\e). Ham paths using e biject with Ham paths of D/e via collapsing ...→u→v→... to ...→w→...
**Corollary:** Arc-flip H-difference reduces to contraction: H(T)-H(T') = H(T/e) - H(T'/e').
**Key structural insight:** T/e and T'/e' differ ONLY in how w connects to other vertices (profile swap). If u,v have identical profiles, T/e = T'/e' and H(T) = H(T').
**Connection to Mitrovic:** Commutative specialization of W_X = W_{X\e} - W_{X/e}↑ (arXiv:2504.20968).
**Scripts:** `04-computation/deletion_contraction_test.py`, `04-computation/flip_reduction_via_contraction.py`

### INV-120: Polynomial Deletion-Contraction for F(T,x) — PROVED (THM-083)
**Source:** kind-pasteur-2026-03-07-S35
**Status:** PROVED algebraically. Verified exhaustive n=4,5.
**What:** F_T(x) = F_{T\e}(x) + (x-1) * F(T/e, x). Generalizes THM-082 to polynomial level.
**Key identification:** G_{u,v}(x) = F(T/e, x) — the "conditional path polynomial" summing over permutations with u immediately before v equals the forward-edge polynomial of the contraction.
**Flip formula as corollary:** F_T - F_{T'} = (x-1) * [F(T/e) - F(T'/e')], with D anti-palindromic.
**Anti-palindromicity proof:** D(x) = -x^{n-2}D(1/x) follows from palindromicity of F_T, F_{T'}.
**Scripts:** `04-computation/poly_deletion_contraction.py`

### INV-121: F(T,omega) mod 9 universality — PROVED (THM-085)
**Source:** kind-pasteur-2026-03-07-S36 (extending S35 analysis)
**Status:** PROVED algebraically. Complete proof via Taylor expansion.
**What:** 9 | F(T,omega) for ALL tournaments on n >= 6 vertices. Proof:
1. Taylor expansion F(T,x) = sum c_k (x-1)^k. Over F_3: x^3-1 = (x-1)^3.
2. c_0 = n! (tournament-indep), c_1 = n!(n-1)/2 (tournament-INDEPENDENT!), both divisible by 3.
3. c_2 = A_non + (n-2)!*dp(T), where A_non is tournament-independent and dp(T) = directed 2-path count. Both A_non and (n-2)! divisible by 3 for n >= 5, so c_2 = 0 mod 3 regardless of T.
4. Therefore (x-1)^3 | F(T,x) mod 3 for n >= 5, giving S_r = 0 mod 3.
5. Combined with v_3(n!) >= 2 for n >= 6: 9 | F(T,omega).
**Additional:** Eulerian conjecture: 3|A(n,k) => 3|F_k(T) for all T (verified n=5-8). But this alone doesn't explain n=9,10 (where all A(n,k) = 1 mod 3); the Taylor proof covers all n.
**Sharp:** n=5 has S_r=0 mod 3 but v_3(5!)=1 blocks mod 9. n=4 has c_2 NOT forced 0.
**Scripts:** `04-computation/c2_mod3_proof.py`, `fk_mod3_conjecture.py`, `sr_mod3_n9_check.py`, `f_omega_mod27_analysis.py`

### INV-122: THM-084 naming fix + Corollary 2 error
**Source:** kind-pasteur-2026-03-07-S36
**Status:** FIXED.
**What:** opus-S46 created THM-082-flip-factorization-anti-palindrome.md, colliding with kind-pasteur's THM-082-deletion-contraction-ham-paths.md. Renamed opus's to THM-084. Also fixed Corollary 2 which incorrectly claimed H(T)=H(T') under arc flip (FALSE: H(T) != H(T') in general, deltas up to +-12 at n=5). Correct: F(T,1)=n!=F(T',1) trivially.

### INV-123: Worpitzky Expansion of F(T,x) — PROVED (THM-084)
**Source:** opus-2026-03-07-S46b
**Status:** PROVED algebraically, verified n=3..7
**What:** F(T,x)/(1-x)^n = sum a_m x^m where a_m is polynomial in m of degree n-1.
  - Top 2 coefficients are UNIVERSAL: n and C(n,2)
  - For transitive tournament: a_m = (m+1)^n - m^n (binomial coefficients)
  - Deviation from binomial: delta_2 = 2(n-2)*t3, delta_3 = (n-2)(n-3)*t3
  - At n=6, deeper coefficients need invariants beyond t3
  - Spectral connection: delta_2 = 2(n-2)/3 * tr(A^3)
**Analogy:** F(T,x) is an h*-vector; a_m is an Ehrhart-like polynomial. Transitive tournament corresponds to unit cube h*-vector.
**Scripts:** `04-computation/worpitzky_coefficients.py`, `04-computation/worpitzky_deeper.py`

### INV-124: Signed Forward-Edge Polynomial SF(T,x) — PROVED (THM-085b)
**Source:** opus-2026-03-07-S46b
**Status:** PROVED algebraically
**What:** SF(T,x) = sum sgn(sigma) x^{fwd_T(sigma)} is palindromic with parity (-1)^{C(n,2)}.
  - SF(T,1) = 0 always, so (x-1) | SF(T,x)
  - Quotient SF/(x-1) is anti-palindromic
  - At n=4: SF = c(T) * (x-1)^2(x+1) (since anti-palindromic of even degree has (x+1) factor)
  - SF determines F at n<=5 but NOT at n>=6 (coarser invariant)
**Connection:** SF is a "path immanant" for the sign character. F is the "path permanent."
**Scripts:** `04-computation/signed_F_analysis.py`

### INV-125: Forward-Edge Variance Formula — PROVED (THM-086)
**Source:** opus-2026-03-07-S46b
**Status:** PROVED
**What:** Var[fwd] = (n+1)/12 + 4*t3/(n(n-1)). Exact formula.
  - Non-adjacent forward indicators are UNCORRELATED (tournament completeness: C(n-2,2) edges always)
  - Adjacent covariance = -1/12 + 2*t3/(n(n-1)(n-2))
  - Directed 2-path count = C(n,3) + 2*t3
  - At t3=0 (transitive): Var = (n+1)/12 = Eulerian variance
**Scripts:** `04-computation/worpitzky_w_connection.py`

### INV-126: Cross-Domain Connections for F(T,x)
**Source:** opus-2026-03-07-S46b
**Status:** EXPLORED
**What:** Multiple connections between F(T,x) and other mathematical structures:
  1. **q-analogue F(T,x,q):** q-marginal (inv distribution) is UNIVERSAL = [n]_q! for all T
  2. **det(W(x)):** det((J-I)+(x-1)A) at x=1 = (-1)^{n-1}(n-1) for all T
  3. **Descent algebra:** SF is the sign-character evaluation of the "tournament descent" element
  4. **Worpitzky = Ehrhart:** a_m analogous to Ehrhart polynomial, F to h*-vector
**Next step:** (1) Explore F(T,x,q) as bivariate polynomial. (2) Find polytope whose h*-vector is F(T,x). (3) Connect Worpitzky coefficients to W-polynomial hierarchy (INV-082).

### INV-127: GLMY Path Homology of Tournaments — EVEN BETTI VANISHING
**Source:** opus-2026-03-07-S46e (path_homology_phase2.py output + new analysis)
**Status:** PARTIALLY CONFIRMED — β₂=0 exhaustive n<=6, sampled n<=9 (0 failures in ~50k tests). β₄=0 exhaustive n<=6, sampled n<=7 (0/3000 random). BUT Paley T_7 has β₄=6! And β₄=1 found in 0.6% of random n=8 tournaments. So only β₂=0 appears truly universal.
**What:** β₂(T) = 0 for ALL tournaments T (HYP-207). β₄ can be nonzero starting at n=7 (Paley) and n=8 (random). β₁ and β₃ NOT mutually exclusive at n=8 (need to check). χ(Ω) ∈ {-11,...,7} at n=7, NOT {0,1} — HYP-267 REFUTED.
**Evidence (S42):** Exhaustive n=3,4,5,6; sampled n=7 (5000+), n=8 (1000+), n=9 (100). β₂=0 in ALL cases.
**Only 3 Betti profiles at n=5,6:** (1,0,...), (1,1,0,...), (1,0,...,β₃=1,...). At n=7: same 3 + Paley's (1,0,0,0,6,0). At n=8: adds (1,0,0,0,1,0,0).
**Key formulas (proved n=5):** dim(A₂) = C(n,3) + 2c₃; dim(Ω₂) = dim(A₂) - #{non-allowed pairs with mediators}; rank(d₂) = C(n-1,2) - β₁.
**Algebraic mechanism:** ker(d₂|Ω₂) = im(d₃|Ω₃) always. The "swap cycle" characterization: pure v-chain 2-cycles have form Σ M_{ab}[(a,b,v)-(v,a,b)] with zero row/col sums. ALL swap cycles are boundaries (confirmed exhaustive n=5,6).
**Literature (S42):** Tang-Yau (2026): H_m=0 for m>=2 for circulant tournaments. Burfitt-Cutler: Ω₂ generated by transitive triples only. No paper addresses β₂=0 for general tournaments — this is genuinely open.
**Not in literature:** Confirmed via comprehensive search of Caputi-Menara, Burfitt-Cutler, Fu-Ivanov, Tang-Yau, Chaplin, all GLMY papers.
**Scripts:** `beta_parity_pattern.py`, `beta2_algebraic_mechanism.py`, `beta2_deformation_retract.py`, `beta2_large_n_sample.py`, `beta_paley_verify.py`, and many more
**Cone-from-T' construction (S42):** For vertex v with swap cycle z, the filling w = Σ α_{abc} [(v,a,b,c)+(a,b,c,v)] over T'=T\{v} 2-paths always works. The T'-internal faces cancel: d₃(v,a,b,c)+d₃(a,b,c,v) has zero (a,b,c) component. The resulting B·α=z system is always solvable.
- **Filtered** (only T' paths with v→a AND c→v): Works exhaustive n=5,6 (32768/32768). Fails 1/1000 at n=8.
- **Unfiltered** (ALL T' paths): Works 500/500 at n=7, 500/500 at n=8, 200/200 at n=9. Zero failures.
- **Multi-vertex** (combine all vertices): ALWAYS works, including the n=8 filtered failure.
- **Ω₃ auto-membership**: Cone filling automatically in Ω₃ at n=5,6 (100%). Breaks at n≥7 (~98% at n=7, ~93% at n=8).
- **Rank surplus grows**: rank(B)-swap_dim min is 2(n=5), 2(n=6), 6(n=7), 11(n=8), 15(n=9). System increasingly overdetermined.
- **β₂=0 confirmed through n=10**: Exhaustive n≤6, sampled n=7(500), n=8(1000), n=9(200), n=10(50). Zero failures.
- **Dimension formula**: rank(d₃) = ker(d₂) EXACTLY for every tournament tested (n=5-9).
**Next step:** (1) Prove B·α=z always solvable algebraically (rank argument). (2) Prove Ω₃ membership of filling. (3) Try inductive proof using LES of pair (T, T\v). (4) Check if result follows from multisquare-free property (Fu-Ivanov).

### INV-128: Universal Coefficient Theorem — PROVED (THM-117)
**Source:** opus-2026-03-07-S46e
**Status:** PROVED
**What:** coeff(t_{2k+1} in κ_{2k}) = 2/C(n, 2k) for all k. Proved via forward path formula + OCF + multinomial expansion. Resolves OPEN-Q-023.
**Scripts:** `universal_coeff_proof.py`

### INV-129: Celano-Sieger-Spiro A_T(t) — NOT same as F(T,x)
**Source:** opus-2026-03-07-S46e (web research)
**Status:** CLARIFIED (dead end for direct application)
**What:** arXiv:2309.07240 defines A_T(t) = sum over labelings of t^{des_T(sigma)} where des counts descents across ALL arcs. This has degree C(n,2), not n-1. The (1+t)^{floor(n/2)} divisibility applies to A_T(t), not to our F(T,x). The two polynomials encode different statistics (all-arc descents vs Hamiltonian-path forward edges).
**Impact:** The Celano-Sieger-Spiro result cannot be directly applied to F(T,x). However, it establishes that tournaments have a universal structural constraint on A_T(t) depending only on n, which is analogous to our universal constraint on F(T,x) mod 3 (THM-086).

### INV-130: Pfaffian-Betti Connection — EXHAUSTIVE at n=6, extended n=7,8
**Source:** opus-2026-03-07-S46e, kind-pasteur-2026-03-08-S40
**Status:** VERIFIED EXHAUSTIVE n=6. Sampled n=7,8. THM-120 (was THM-098) + THM-099 documented.
**What:** The Pfaffian of the skew-adjacency matrix constrains path homology Betti numbers. At n=6 (exhaustive): β₁>0 ⟹ |Pf(S)| ∈ {1,3}; β₃>0 ⟹ |Pf(S)| ∈ {7,9}. Perfect separation. At n=7 (odd): spectral gap separates phases. At n=8: |Pf| NOT perfect separator but strongly correlated.
**CORRECTED (S40):** H-maximizers at n=6 are NOT all S-phase. 480 maximizers split 240 C-phase (|Pf|=1) + 240 S-phase (|Pf|=7), both with score (2,2,2,3,3,3) and c3=8. Complementation preserves phase.
**Scripts:** `pfaffian_betti_check.py`, `pfaffian_betti_n7.py`, `pfaffian_topology_deep.py`, `pfaffian_betti_mechanism.py`, `spectral_betti_gap.py`, `spectral_topology_n8.py`, `s_phase_structure.py`, `s_phase_maximizer_n7.py`, `maximizer_betti_deep.py`
**Next step:** (1) Prove Pfaffian separation algebraically at n=6. (2) Why does spectral gap separate at n=7?

### INV-135: H-Maximizer Betti Dimension Shift — THM-099
**Source:** kind-pasteur-2026-03-08-S40
**Status:** VERIFIED EXHAUSTIVE n=4,5,6. Sampled n=7.
**What:** H-maximizers always have nontrivial GLMY path homology, with the topological dimension increasing:
- n=4: ALL 24 max have β₁=1 (C-phase)
- n=5: ALL 64 max have β₁=1 (C-phase)
- n=6: 480 max split 240 β₁=1 + 240 β₃=1
- n=7: ALL 240 max have β₄=6 (beyond S-phase classification)
At n=7, all maximizers are conference-matrix (gap=0, eigenvalues all √7). Second-highest H=175 has β₁=1 (C-phase). Third H=171 is contractible. Topology stratifies H values.
**Scripts:** `betti_dimension_shift.py`, `betti_dimension_shift_v2.py`, `maximizer_betti_deep.py`, `maximizer_betti_n8.py`
**Next step:** (1) Check n=8 maximizers (H=661). (2) Why β₄=6 specifically at n=7? (3) Algebraic mechanism connecting eigenvalue uniformity to high Betti.

### INV-131: Path Homology Hidden Invariant — Cycle Overlap Pattern
**Source:** opus-2026-03-07-S46e
**Status:** VERIFIED EXHAUSTIVE n=5
**What:** Path homology β₁ is NOT determined by (F-polynomial, t₃, score sequence, strong connectivity). At n=5, the distinguishing invariant is the 3-cycle overlap pattern: β₁=1 iff all 3-cycles form a "star" (share a common edge); β₁=0 when cycle overlaps are heterogeneous. Two non-isomorphic tournament types with identical F-poly and t₃ have different β₁.
**Scripts:** `betti_hidden_invariant.py`, `f_poly_betti_deep.py`
**Next step:** (1) Formalize the overlap pattern criterion. (2) Check at n=6,7.

### INV-132: Mod-2 Cumulant Collapse via OCF
**Source:** opus-2026-03-07-S46e
**Status:** VERIFIED
**What:** THM-094 (F mod 2 universal) implies all integer moment sums n!·E[fwd^r] ≡ 0 mod 2 for n≥3. Root cause: OCF evaluates at x=2 ≡ 0 mod 2, so I(Ω,2) mod 2 = I(Ω,0) = 1, giving H ≡ 1 = Rédei. Combined with THM-086 (mod 3 zeros), all moments vanish mod 6 for n≥3.
**Scripts:** `mod2_cumulant_connection.py`

### INV-133: Tang-Yau Circulant Path Homology (arXiv:2602.04140)
**Source:** opus-2026-03-07-S46e (web research)
**Status:** NEW LEAD
**What:** Tang-Yau (Feb 2026) compute GLMY Betti numbers of circulant digraphs using Fourier decomposition via shift automorphism. Additive structure of connection set S governs low-dimensional chains. Could be adapted to circulant tournaments (odd n, S={1,...,(n-1)/2}).
**Next step:** Read paper. Check if their Fourier method can prove β₂=0 for tournament subclass.

### INV-134: GLMY Torsion-Zeta Bridge
**Source:** opus-2026-03-07-S46e (web research)
**Status:** NEW LEAD
**What:** Grigoryan-Muranov define Reidemeister/analytic torsion on digraphs via path homology Hodge Laplacian, proving Cheeger-Müller analogue. Torsion invariants are classically expressible through zeta-type regularized determinants. This could connect the Ihara/Bartholdi zeta of a tournament to its path homology.
**References:** arXiv:2012.07302, grigor/torsion.pdf
**Next step:** Investigate whether Bartholdi zeta of a tournament encodes path homology info.

### INV-135: Asao — Transitive Tournaments Have Trivial Path Homology
**Source:** opus-2026-03-07-S46e (web research, arXiv:2503.06722)
**Status:** CONFIRMED (known)
**What:** Asao (March 2025) proves transitive tournaments are "regularly diagonal" digraphs with vanishing reduced regular path homology. Confirms our computational finding that transitive tournaments have β=(1,0,...,0).
**Impact:** Provides a published reference for one endpoint of the topological trichotomy.

### INV-136: Chaplin — Random Digraph β₁ Phase Transitions
**Source:** opus-2026-03-07-S46e (web research, arXiv:2111.13493)
**Status:** NEW LEAD
**What:** Chaplin (2022) shows β₁ of random Erdős-Rényi digraphs has two phase transitions. Since tournaments are "density 1/2" digraphs, this places them in a specific regime. Could explain why ~30% of tournaments at n=5 have β₁>0.
**Next step:** Check if their density threshold matches tournament β₁ fraction.

### INV-137: THM-118 Trace-Cycle Identity — PROVED (extended to k=3,4,5)
**Source:** kind-pasteur-2026-03-07-S39b
**Status:** PROVED
**What:** tr(A^k) = k * c_k(T) for k=3,4,5 in any tournament. Extended to k=4 (no bidirectional edges => length-4 closed walks must be simple 4-cycles). Gives O(n^3) c_4 and c_5 computation via matrix multiplication. Sharp: fails at k>=6 (compound (3,3) walks at k=6). Correction for k=6 is NOT a simple polynomial in global cycle counts (tested and failed).
**Impact:** c4_fast() and c5_fast() in tournament_fast.py. Speedups: 3.8x for c4, 5.4x for c5 at n=8.
**Scripts:** `trace_cycle_k4.py`, `c6_correction_formula.py`, `c6_from_trace.py`

### INV-138: Björklund Cycle Cover — OCF Connection Tested (NEGATIVE)
**Source:** kind-pasteur-2026-03-07-S39b
**Status:** TESTED (NEGATIVE for new identities)
**What:** Tested 6 Björklund-style cycle cover formulations for OCF connections. Only Test 2 (partial odd cycle cover weighted by 2^{num_cycles}) matches OCF — but this IS OCF restated. No new route to proving OCF found. Permanent of A+I counts cycle covers but doesn't simplify OCF.
**Scripts:** `bjorklund_cycle_cover.py`

### INV-139: h-Positivity of U_T — CLOSED (fails for all non-transitive)
**Source:** kind-pasteur-2026-03-07-S39b
**Status:** CLOSED (dead end)
**What:** U_T is NOT h-positive for any non-transitive tournament. Only the transitive tournament (H=1) has h-positive U_T. This closes the Stanley-Stembridge connection for tournament Rédei-Berge functions. The e-positivity question from INV-051/052 is also resolved negatively.
**Scripts:** `bjorklund_cycle_cover.py` (h-positivity test section)

### INV-140: THM-097 Alpha_2 Trace Formula — PROVED
**Source:** kind-pasteur-2026-03-07-S39b
**Status:** PROVED. O(n^3) computation of vertex-disjoint 3-cycle pairs.
**What:** alpha_2(Omega_3) = C(c3,2) - sum_v C(t3(v),2) + s2, where t3(v) = (A^3)[v][v] and s2 = sum_{edges a->b} C((A^2)[b][a], 2). Proof via inclusion-exclusion on pair overlap counts. Valid for full Omega at n<=7 (since 5+3=8>7 prevents cross terms). Implemented as alpha2_from_trace() in tournament_fast.py.
**Scripts:** `trace_ocf_bridge.py`, `alpha2_formula.py`

### INV-141: H(T) Polynomial Trace Formula — VERIFIED n<=9
**Source:** kind-pasteur-2026-03-07-S39b
**Status:** VERIFIED. 100% match at n=5,6 (exhaustive), n=7 (500), n=8 (100), n=9 (200).
**What:** H(T) = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3, all computable from matrix trace data. At n<=7: O(n^3). At n=8: O(n^5) (cross terms). At n=9: O(n^5) with additional O(2^7*C(n,7)) for c7 and O(C(n,3)^3/6) for alpha_3. Timing: trace method 7x slower than DP at n=9 but POLYNOMIAL.
**Key n=9 findings:** alpha_3 nonzero 86% of time (3+3+3=9). alpha_2(3,5) cross dominates alpha_2(3,3) by 2:1. H contribution: 56% alpha_1, 41% alpha_2, 2.3% alpha_3.
**Scripts:** `h_from_trace_n8.py`, `h_polynomial_n9.py`, `alpha_structure_n9.py`

### INV-142: Spectral Characterization of H-Maximizers — NEW FINDING
**Source:** kind-pasteur-2026-03-07-S39b
**Status:** COMPUTED. Key structural finding.
**What:** Paley T_p has the CONFERENCE MATRIX property: S^2 = -pI + J where S = A-A^T. This means ALL nonzero skew eigenvalues equal ±i*sqrt(p) (zero spectral gap). This property CHARACTERIZES Paley among DRTs at n=11 (non-Paley would fail). At n=7: Paley has zero spectral gap while other regular tournaments have gap 3.46-3.90. The general spectral correlation (|skew_max| vs H) is weak (-0.03) but among REGULAR tournaments, zero spectral gap → max H.
**Also found:** tr(S^k) = 0 for ALL odd k in ALL tournaments (skew-symmetry). The adjacency spectrum does NOT distinguish DRTs (all have same eigenvalues), only the skew spectrum does. Paley's conference matrix property is a VERY strong constraint.
**Scripts:** `spectral_h_maximizer.py`, `spectral_cycle_density.py`

### INV-143: MISTAKE-017 — Invalid DRT at n=11 Corrected
**Source:** kind-pasteur-2026-03-07-S39b
**Status:** CORRECTED
**What:** The "non-Paley DRT" from {1,2,3,5,8} was NOT a valid tournament (S∩(-S)={3,8}≠∅). All claims about c3=44, c5=407, H=69311 were computed on a non-tournament digraph. The only valid circulant DRT at n=11 is Paley. Exhaustive search found exactly 2 valid (11,5,2)-difference sets in Z_11: QR and NQR, which give isomorphic tournaments.
**Impact:** INV-068 corrected. MEMORY.md and TANGENTS.md updated.

### INV-144: Circulant Digraph Path Homology (arXiv 2602.04140, Feb 2026) — CONJ 4.8 DISPROVED
**Source:** opus-2026-03-08-S40 (web search), kind-pasteur-S41 (counterexample)
**Status:** CONJECTURE 4.8 DISPROVED. New characterization found.
**What:** Uses exactly our Fourier eigenspace decomposition approach for circulant digraphs. Key results:
- Strong Stability (Thm 4.5): Betti numbers constant for large primes
- ~~Conjecture 4.8: H_m = 0 for m >= 3 under "no-wrap-around" condition~~ **FALSE**
- S={1,s} with s!=2 gives H_2 = K (nonzero!)
- No results on tournaments or Paley specifically
**COUNTEREXAMPLES to Conj 4.8 (kind-pasteur-S41):**
- C_8^{1,5}: |S|=2, S cap (-S) = empty, but beta_3=1, beta_4=1
- C_8^{3,7}: same structure, also beta_3=1, beta_4=1
- P_7 = C_7^{1,2,4}: tournament with beta_4=6
- Z_9 = C_9^{1,5,6,7}: tournament with beta_5=10
- Their conjecture may hold for |S|=1 only (directed cycles have beta=[1,1,0,...])
**NEW FINDING (HYP-213):** For |S|=2, beta_2=0 iff {s1,s2} is "doubling-closed" (2s1=s2 or 2s2=s1 mod n). Perfect correlation at n=5,7,9,11,13. One exception at n=8 (s2-s1=n/2).
**Relevance:** Their Fourier decomposition matches our per-eigenspace approach. For tournament beta_2=0, the mechanism is tournament completeness, NOT the Fourier structure.
**Scripts:** tang_yau_counterexample.py, beta2_nonzero_analysis.py
**Next step:** (1) Notify Tang-Yau of counterexamples. (2) Investigate whether their techniques prove beta_2=0 for tournaments specifically. (3) Generalize doubling-closure to larger |S|.

### INV-145: Ω_2 Structure — Cancellation Chains in Tournaments
**Source:** opus-2026-03-08-S40
**Status:** DISCOVERED
**What:** Ω_2 ≠ span(transitive triples). Non-transitive 2-paths with shared non-allowed faces form "cancellation chains" in Ω_2. Gap dim(Ω_2) - |TT| ranges 0-5 at n=5. Cancellation chains never individually in ker(∂_2), but mixed elements (TT + cancellation) can be 2-cycles.
**Impact:** Previous β_2 analysis assumed Ω_2 = TT, which was incomplete. Corrected computation still gives β_2 = 0 through n=6 (exhaustive).

### INV-146: P_11 Path Homology — Non-palindromic Ω Dims
**Source:** opus-2026-03-08-S40
**Status:** COMPUTING (dims 8-10 in progress)
**What:** P_11 per-eigenspace Ω dims: [1, 5, 20, 70, 205, 460, 700, 690, ?, ?, ?].
Inner sequence NOT palindromic: 460≠700, 700≠690. Contrasts with P_7's palindromic [3,6,9,9,6,3].
Using J^H J + eigvalsh method for memory-efficient rank computation.
**Next step:** Complete Ω_8, Ω_9, Ω_10 to determine Betti concentration dimension.

### INV-147: Eigenspace Decomposition of β_top — Trivial vs Non-trivial Split
**Source:** kind-pasteur-2026-03-08-S41
**Status:** VERIFIED (P_7 and Z_9)
**What:** For circulant maximizers, β_top decomposes across Z/nZ eigenspaces as:
- P_7: trivial (k=0) gives β_4=0, each non-trivial (k=1..6) gives β_4=1, total = 6
- Z_9 S={1,5,6,7}: trivial (k=0) gives β_5=2, each non-trivial (k=1..8) gives β_5=1, total = 2+8=10
- Formula: β_top = (n-1) + δ, where δ=0 for prime n (P_7), δ=2 for n=9 (9=3²)
- All eigenspaces have IDENTICAL Om_5 dim=74 and Om_6 dim=63
- The difference: trivial has ker(∂_5)=39, non-trivial has ker(∂_5)=38
**Key question:** Why does trivial eigenspace contribute extra at n=9? Is δ=2 because 9=3²? What is δ for n=11 (Paley, prime)?
**Conjecture (HYP-212):** δ=0 for prime n, δ>0 for composite. CONFIRMED for P_11: β_8 = 10 = p-1 + 0 (opus-S42 + kind-pasteur-S41 independent confirmation).
**P_11 data (opus-S42):** Om dims (k≠0): [1, 5, 20, 70, 205, 460, 700, 690, 450, 180, 30]. β_8^(triv)=0 (kind-pasteur-S41 confirmed), β_8^(k≠0)=1 per eigenspace.
**Scripts:** `04-computation/n9_beta5_eigenspace.py`, `04-computation/p7_eigenspace_verify.py`, `04-computation/p11_beta8_v5.py`
**Next step:** (1) Test another composite (n=15?). (2) Find algebraic reason for δ. (3) Extend to P_13?

### INV-148: Arc-Flip Induction Proof for β₂=0 — STRONGEST LEAD
**Source:** kind-pasteur-2026-03-08-S41, opus-2026-03-08-S43 (arc-flip local invariance)
**Status:** VERIFIED EXHAUSTIVELY n=5,6. Key structural mechanism identified.
**What:** β₂=0 can potentially be proved by arc-flip induction:
1. Base: transitive tournament (β₂=0 trivially)
2. Step: flipping any arc preserves β₂=0
The "surplus" = dim(Ω₃) - dim(Z₂) satisfies surplus ≥ |drop| ALWAYS.
**Key findings (kind-pasteur-S41):**
- n=5: 10240 arc flips, 0 violations. Surplus=0 cases: max_drop=0
- n=6: 491520 arc flips, 0 violations. Surplus=1 cases (tightest): max_drop=0
- Surplus=0 stability mechanism: joint (δΩ₃, δZ₂) = {(0,0), (1,0), (2,1), (4,2)} — always δΩ₃ ≥ δZ₂
- "Every new Z₂ cycle comes with at least one new Ω₃ chain to fill it"
- This is the 2-for-1 principle: tournament completeness ensures enough Ω₃ for every Z₂
**Key lemma needed:** For any tournament T with arc (u,v), and T' = flip(T, u, v):
  dim(Ω₃(T')) - dim(Z₂(T')) ≥ dim(Ω₃(T)) - dim(Z₂(T))  when starting surplus ≤ 1
  (or more generally: surplus(T') ≥ 0)
**Why completeness matters:** Non-tournament arc flips CAN create β₂>0 (seen in circulant digraphs).
The tournament constraint ensures every pair of vertices has an arc, providing the intermediary
vertices needed for Ω₃ chains.
**NEW FINDINGS (kind-pasteur-S41 continued):**
- THM-121 (was THM-100) PROVED: delta_|A_3| = (n-3)*delta_|A_2| exactly, for ALL tournaments, ALL arcs
- delta_|A_2| = 2*(d_u - d_v - 1) depends ONLY on out-degrees
- n=7 sampling (10k): 0 violations, min surplus = 9
- n=8 sampling (20k): 0 violations, min surplus <= 25
- Min surplus floor: 0, 1, 9, <=25 for n=5,6,7,8 — grows super-linearly
- Transitive tournament: surplus = C(n-1,4). One-flip delta = -(n-2-gap)
- DT paths: |DT| >= dim(Z_2) for 100% (n=5), 97.1% (n=6). Rest filled by cancellation
- Omega_2 NOT just TT paths: dim(O2) > |TT| for 76.6% (n=5), 94.6% (n=6)
- Tang-Yau Cor 3.15: H_m=0 for m>=2 when S={1,...,d} — applies to circulant tournaments
- Algebraic identity: surplus = beta_3 + rk(d_4) - beta_2
**Scripts:** beta2_arcflip_proof.py, beta2_surplus_zero_stability.py, beta2_arcflip_n7_sample.py,
  beta2_arcflip_mechanism.py, beta2_arcflip_counting.py, beta2_delta_ratio_*.py,
  beta2_min_surplus*.py, beta2_omega_ratio.py, beta2_injectivity_analysis.py, beta2_surplus_formula.py
**Next step:** (1) Prove the key lemma: surplus(T') >= 0 for all arc flips, using THM-121
  (2) Generalize Tang-Yau deformation retract to non-circulant tournaments
  (3) Prove beta_2 = 0 by induction on number of flips from transitive

### INV-149: β₂=0 Density Threshold for Circulant Digraphs
**Source:** kind-pasteur-2026-03-08-S41
**Status:** CHARACTERIZED. New conjecture (HYP-219).
**What:** For C_n^S with S∩(-S)=∅, β₂=0 when |S| is large enough:
- n=7: |S|≥3 (all β₂=0)
- n=9: |S|≥4 (all β₂=0)
- n=11: |S|≥4 (all β₂=0)
- n=13: |S|≥5 (all β₂=0, 96/96 tested)
- n=15: |S|=5 still has 15/201 failures; threshold at |S|≥6?
- |S|=2 perfect characterization: β₂=0 iff has-doubling-pair (HYP-217)
- Exceptions without doubling pair: coset structure (S = a + H for subgroup H)
**Scripts:** beta2_doubling_closure_general.py, beta2_threshold_analysis.py
**Next step:** (1) Find exact threshold formula (2) Prove for tournaments (|S|=(n-1)/2)

### INV-148: Nesting Obstruction Theory — Why H=7 is Simplex-in-Cuboid
**Source:** opus-2026-03-14-S71f
**Status:** DISCOVERED. Algebraic framework established, geometric interpretation given.
**What:** H=7 = (1+2(1+x))|_{x=2}: composing simplex (1+x) into cuboid (1+2y) gives 3+2x, which has constant term 3 ≠ 1 → NOT a valid independence polynomial. Tournament geometry prevents nesting/composition of brick structures, leaving only multiplicative (disjoint union) decomposition.
**Key results:**
- H=7 impossible = tesseract impossibility (HYP-1230)
- H=21 = 3×7 inherits obstruction multiplicatively (THM-079)
- α₁=3 forces α₂≥2 at n=7,8 (verified 500k samples each) → H=15 not 7
- Structural proof for all n sketched (pigeonhole on 9 vertex slots in n≥7 vertices)
**Scripts:** nesting_obstruction.py, alpha1_3_n8_verify.py
**Next steps:** (1) Complete rigorous proof that α₁=3 → α₂≥2 for ALL n. (2) Publish nesting obstruction as appendix to main paper.

### INV-149: The (z-2)(z-3) = 0 Recurrence — Simplex-OCF Bridge
**Source:** opus-2026-03-14-S71f
**Status:** DISCOVERED. Complete algebraic analysis.
**What:** The characteristic polynomial (z-2)(z-3) with roots 2 (OCF point) and 3 (simplex value) generates the forbidden sequence via the pure z=3 orbit from seed 7: {7, 21, 63, 189, ...}. Only first two are permanently forbidden; 63+ are achievable at n≥8. Mixed orbit H = 3^{k+1} - 2^{k+1} gives all achievable values.
**Scripts:** knacci_tournament_recurrence.py
**Next steps:** (1) Does the recurrence connect to deletion-contraction? (2) Higher-order recurrences (z-2)(z-3)(z-5)=0 for cuboid inclusion?

### INV-150: The 2-Bridge — Unified Origin of the Factor 2
**Source:** opus-2026-03-14-S71f (connecting kind-pasteur S72 Degree Drop Theorem)
**Status:** DISCOVERED. Three manifestations verified.
**What:** The number 2 appears in: (a) I(Ω,2)=H (OCF evaluation), (b) top-degree coefficients ±2 (Degree Drop Theorem), (c) ΔH=2·Δα₁ (arc flip derivative). All arise from the binary arc choice / path reversal involution. 1+(-1)^{n-1} = 2 for odd n.
**Scripts:** degree_drop_packing.py, delta_h_gap.py
**Next steps:** (1) Prove ΔH=2·Δα₁+4·Δα₂ for general n. (2) Connect to Vassiliev theory formally.

### INV-151: Simplicial Selection Interpretation
**Source:** opus-2026-03-14-S71f
**Status:** DISCOVERED. Geometric framework.
**What:** Each tournament T selects H(T) simplices from the standard n!-simplex triangulation of [0,1]^n. H/n! = 1/2^{n-1} on average. The f-polynomial of Δ^{n-1} at x=2 gives 3^n (simplex brick), and the f-polynomial of □^n at x=2 gives 5^n (cuboid brick).
**Scripts:** simplex_cuboid_geometry.py
**Next steps:** (1) Is the selection pattern random-looking or structured? (2) Does the geometric constraint explain the ΔH gap?

---

## Engineering & Cross-Domain Leads

### INV-180: Tournament Structure of Transformer Attention Patterns
**Source:** kind-pasteur-2026-03-21-S12 (analysis of Napolitano "Mathematics Is All You Need")
**Status:** INITIAL EXPLORATION. Computational proof of concept complete.
**What:** Threshold attention matrices to get tournaments on token sequences, then apply our full OCF / path homology / spectral machinery. Key findings:
1. OCF verified on ALL attention-derived tournaments (200/200 at n=3..6)
2. Cartan decomposition of gl(n,R) = so(n) + p + R decomposes attention into "tournament" (antisymmetric) and "similarity" (symmetric) parts
3. Dark/active ratio = (n+1)/(n-1) exactly. For n=4 (Napolitano): 10/6 = 5/3
4. Random softmax attention puts ~72% energy in symmetric sector (beyond dimensional prediction)
5. Soft tournament (differentiable thresholding) converges to hard tournament as tau->0, making invariants differentiable
6. n=11 is the UNIQUE order where regular tournament transitivity = EXACTLY 2/3
**Papers analyzed:**
- Napolitano (Zenodo 19120857): LOW rigor, metaphorical physics. Empirical finding (dark modes carry correctness) possibly interesting
- van Nierop (arXiv:2412.14543): MODERATE rigor. SO(d-1) gauge symmetry from LayerNorm is genuine
- NeurReps 2025 (OpenReview YC9O7OyLFK): HIGH rigor. Principal fiber bundle on transformer params
- GET (NeurIPS 2021): HIGH rigor. Gauge equivariant transformer
**Scripts:** tournament_attention_analysis.py, cartan_attention_theorem.py, phase_transition_universality.py, tournament_probe_design.py
**Reflection:** 07-reflections/tournament-gauge-bridge.md
**Full analysis:** 03-artifacts/drafts/napolitano-gauge-theory-analysis-S12.md
**Next steps:**
1. Run TournamentProbe on actual LLM (GPT-2 or similar) to test if trained attention has different tournament profile than random
2. Test if H(T_attention) correlates with model correctness
3. Test if training shifts Cartan energy from symmetric to antisymmetric (making attention more tournament-like)
4. Build TournamentProbe as PyPI package (parameter-free LLM analyzer)
5. Investigate soft OCF: does I(Omega(T_soft), 2) approximate H(T_hard) for soft tournaments?

### INV-181: The n=11 Transitivity Theorem
**Source:** kind-pasteur-2026-03-21-S12
**Status:** PROVED (trivial algebra).
**What:** n=11 is the UNIQUE order where the transitivity of a regular tournament equals exactly 2/3. Formula: transitivity = 3(n-3)/(4(n-2)), setting = 2/3 gives n=11 uniquely. Asymptotic: transitivity -> 3/4 as n -> inf.
**Connection:** Coincidence with Napolitano's 67% depth phase transition. Probably NOT deep, but n=11 being special in tournament theory (it's our Paley T_11) adds to its mystique.
**Next steps:** Investigate other "magic" values of transitivity. Does transitivity = 1/2 have special meaning? (Solve: n = 5.)

### INV-182: Soft OCF — Differentiable Independence Polynomial
**Source:** kind-pasteur-2026-03-21-S12
**Status:** PROOF OF CONCEPT. Soft H converges to hard H.
**What:** Replace hard tournament T[i,j] in {0,1} with soft tournament T_soft[i,j] = sigmoid((A[i,j]-A[j,i])/tau). The soft Hamiltonian path count H_soft(T) = sum_sigma prod T_soft[sigma(k)][sigma(k+1)] is a differentiable function of the attention matrix A. As tau->0, H_soft -> H_hard (verified computationally).
**Application:** Use H_soft as a DIFFERENTIABLE LOSS TERM to encourage Paley-like attention structure during training. If Paley attention is optimal (as our maximizer results suggest), this could improve model performance.
**Scripts:** tournament_probe_design.py
**Next steps:**
1. Define soft conflict graph and soft independence polynomial
2. Prove soft OCF: does I_soft(Omega_soft, 2) approximate H_soft?
3. Test as regularizer in a small transformer training run

### INV-183: The Antiferromagnetic Tournament Framework
**Source:** opus-2026-04-04-S16
**Status:** ESTABLISHED. Dictionary validated exhaustively at n=4,5,6.
**What:** Tournaments are antiferromagnets on the staircase lattice. The frustration index c₃/C(n,3) correlates with H at r≈0.97. Score variance = Néel order parameter. Regular tournaments = AFM ground state. Magnon spectrum is flat over labeled ensemble (S_n isotropy). Boltzmann-weighted correlations show true AFM order at β>0. Phase transition at β_c≈0.7 (n=5). H=7 gap explained by frustration propagation (α₁ gap: 3 three-cycles force a 5-cycle).
**Key files:** antiferromagnetic_tournament_s15.py, afm_deep_analysis_s15.py, afm_remaining_s15.py, the-antiferromagnetic-tournament.md
**Open directions:**
1. Per-class magnon dispersion (condition on iso class for anisotropy)
2. β_c(n) scaling as n→∞ (thermodynamic limit)
3. Yang-Lee zeros of Z(β) in complex plane
4. Connection between seesaw β₁·β₃=0 and AFM selection rules
5. α₁ gap mechanism for H=21 (needs n≥7 analysis)
6. The staircase lattice as Brillouin zone of the tournament AFM
**Priority:** HIGH (provides physics intuition for all open questions)

### INV-052: Paley maximization vs Interval — complete status (opus-2026-04-16-S1)
**Source:** opus-2026-04-16-S1 (this session)
**Status:** RESOLVED for p=7,11,19,23; OPEN for p≥31
**What:** Does Paley tournament T_p maximize H among all circulant tournaments? Among all tournaments?
**Results:**
- p=7:  H(Paley)=189 > H(Interval)=175 → Paley wins. EXHAUSTIVE: Paley = global max at n=7.
- p=11: H(Paley)=95,095 > H(Interval)=93,027 → Paley wins (exhaustive among circulants).
- p=19: H(Paley)=1,172,695,746,915 < H(Interval)=1,184,212,824,763 → Interval wins.
- p=23: H(Paley)=15,760,206,976,379,349 < H(Interval)=16,011,537,490,557,279 → Interval wins by 1.57%.
**Mechanism:** At p=23, Paley has MORE α₁ and α₂ than Interval but loses on α₃+.
  The term 8(α₃_I - α₃_P) = +152T swamps the combined 2Δα₁+4Δα₂ = +86T advantage.
**Dominant term crossover:** k=1→2 at n≈11; k=2→3 predicted at n≈29. The Paley→Interval crossover 
  happens between p=11 and p=19, likely because k=2 becomes the dominant term around n=11.
**Open questions:**
  (1) Does INTERVAL maximize H among ALL tournaments (not just circulants) for n≥15?
  (2) Is the Paley→Interval crossover at p=13? (p=13 ≡ 1 mod 4, no Paley, skip to p=19)
  (3) Why does Interval win on α₃,α₄,...? Is this related to arithmetic spacing allowing more vertex-disjoint packings?
**Scripts:** alpha_from_cc_bin.py, alpha_crossover_analysis.py, alpha_ratio_trends.py

---

## 2-Adic Column Family Investigations (oracle-2026-05-15)

### INV-187: SC Blowup — $\Omega(T_{\mathrm{SC}})$ Structure and Cross-Lane Cycles
**Source:** oracle-2026-05-15-S2 (sc-blowup-and-twin-gaining.md)
**Status:** OPEN. Universal Score Theorem proved. H concentration data at n=3,4,5.
**Key results (verified exhaustively n=3,4,5):**
- Universal Score: all v₀ have out-degree n, all v₁ have out-degree n-1 (PROVED).
- H_SC varies only 4.2% at n=5 (14937–15565). Minimized by transitive, maximized by regular.
- H_SC = H_Lex iff T is regular (at n=3,5; first departure at Paley(7): H_SC=24453597 ≠ H_Lex=24589929).
- SC preservation: T_SC is SC iff T is SC (PROVED via anti-automorphism τ(v_i) = σ(v)_{1-i}).
- Kronecker: A(T_SC) = A(T)⊗I₂ + A(T)ᵀ⊗Φ + I_n⊗e₀₁.
- Eigenvalue splitting (circulant): λ_{k,0}(T_SC) = 2Re(λ_k(T))+1, λ_{k,1} = 2i·Im(λ_k(T))−1.
**Open questions:**
1. What is Ω(T_SC) in terms of Ω(T)? Same-lane cycles = copies of Ω(T); cross-lane cycles = ?
2. Is there a formula H(T_SC) = f(I(Ω(T), x)) for some x or some operation on the polynomial?
3. Why does T_Lex ≅ T_SC for regular n=3,5 but NOT for Paley(7)? What breaks at n=7?
4. Tower: H(Trans_n^SC) for n=3..?: 1, 41, 530293, ... No recurrence found. OEIS?
5. HYP-SC-1 (H monotonicity): H(T₁) ≤ H(T₂) ⟹ H_SC(T₁) ≤ H_SC(T₂)? (Checked n=3,4; needs n=5 verify.)
**Scripts:** blowup_study3.py, results: blowup_study.out
**Next steps:** (1) Compute α₁(T_SC) for C3 and Trans-3. (2) Find Ω(T_SC) for n=3 cases by direct cycle enumeration. (3) Look for H_SC formula via independence polynomial at different evaluation point.

### INV-184: Tournament Blowup $T[K_2]$ — H Formula and Family Inheritance
**Source:** oracle-2026-05-15 (2-adic column family analysis; see `07-reflections/adic-column-families.md`)
**Status:** OPEN. Computationally accessible immediately.
**What:** The tournament blowup $T[K_2]$ — each vertex $v$ splits into $v_0 \to v_1$,
each original arc $u \to v$ expanded to all four arcs $u_i \to v_j$ — is the concrete
realization of the "row step" in the 2-adic column family grid: $(r, k) \to (r+1, k)$.
It doubles $n$ and stays within the same column family.
**Key questions:**
1. Is there a formula $H(T[K_2]) = f(H(T), n)$?
2. Does blowup preserve SC status? SF status? If yes, the column family inherits SC/SF structure.
3. The pairs anomaly: $\lfloor n/2 \rfloor$ gains +1 extra at the odd→even ($r=0 \to r=1$)
   transition only. Does H have a similar anomalous first-blowup jump?
4. What is the isomorphism class of $T[K_2]$? Does it stay in the same G_n/Z_2 sector?
**Expected behavior:** Blowup creates a natural Z_2 action (swap $v_0 \leftrightarrow v_1$
is NOT an automorphism, but the structure is highly symmetric). This may force the blowup
into SC or near-SC classes.
**Connection:** Linial-Morgenstern conjecture (INV-013) uses "random blowup of transitive
tournament" — the row step applied to the transitive class. Our framework says this is the
column-1 family at any depth.
**Computed results (oracle-2026-05-15, exhaustive n=3,4,5):**

| $T$ | $n$ | $H(T)$ | score$(T)$ | $H(T[K_2])$ | score$(T[K_2])$ |
|-----|-----|--------|------------|-------------|-----------------|
| Transitive | 3 | 1 | (0,1,2) | 1 | (0,1,2,3,4,5) |
| Cyclic C3 | 3 | 3 | (1,1,1) | **45** | (2,2,2,3,3,3) |
| Transitive | 4 | 1 | (0,1,2,3) | 1 | (0,1,2,3,4,5,6,7) |
| H=3 class | 4 | 3 | (0,2,2,2) | **45** | (0,1,4,4,4,5,5,5) |
| H=5 class | 4 | 5 | (1,1,2,2) | 393 | (2,2,3,3,4,4,5,5) |
| C5^{1,2} (interval) | 5 | 15 | (2,2,2,2,2) | 15565 | (4,4,4,4,4,5,5,5,5,5) |
| Paley(5) C5^{1,4} | 5 | 10 | (2,2,2,2,2) | 7910 | (4,4,4,4,4,5,5,5,5,5) |

**STRIKING FINDING:**
- Blowup of cyclic C3 (max-H at n=3) → H=45 = **max H at n=6** ✓
- Blowup is near-regular: scores $(n/2-1)^{n/2} \cup (n/2)^{n/2}$ — exactly the SC∩SF score.
- Blowup of max-H regular tournament → max-H SC∩SF class at the doubled size.
- This matches the max-H pattern: "even $n$ → max H by SC∩SF class."

**CONJECTURE (HYP-new):** Blowup of the max-H regular tournament at odd $n$ IS the
max-H SC∩SF class at even $2n$. Equivalently: the column-family row step maps
max-H-regular to max-H-SC∩SF.

**Next steps:**
1. Verify: is H=15565 = max H at n=10? (Run exhaustive or sampling-based search.)
2. Verify: is H=393 = max H at n=8? Compare to known SC∩SF n=8 classes.
3. Prove: blowup of regular tournament → near-regular tournament → SC∩SF class.
4. Connect to the SC∩SF = SC(n-2) identity: does the blowup construction provide the bijection?

### INV-185: HYP-217 Proof via 2-Adic Orbit Structure
**Source:** oracle-2026-05-15 (2-adic column family analysis; see `07-reflections/adic-column-families.md`)
**Status:** OPEN. Theory suggests a proof route.
**What:** HYP-217 says: for circulant tournament $C_n^{\{s, 2s \bmod n\}}$, $\beta_2 = 0$.
A "doubling pair" $\{s, 2s \bmod n\}$ is a 2-adic orbit of size 2 in $\mathbb{Z}/n\mathbb{Z}$.
**Proposed proof route:** The column family structure of $n = 2^r(2k-1)$ induces a splitting
of the chain complex $\Omega_*(C_n^S)$ along 2-adic eigenspaces. When $S = \{s, 2s\}$, the
doubling pair selects a single eigenspace in this splitting, and the chain complex restricted
to that eigenspace is contractible (hence $\beta_2 = 0$).
**Predicted generalization (stronger than HYP-219):** $\beta_2(C_n^S) = 0$ whenever $S$
contains at least one complete 2-adic orbit $\{s, 2s, 4s, \ldots, 2^{d-1}s\} \bmod n$
where $d$ is the multiplicative order of 2 modulo $n / \gcd(s, n)$.
**Connection to density threshold (INV-149):** The threshold $|S|$ needed for $\beta_2 = 0$
should equal the size of the smallest complete 2-adic orbit in $\mathbb{Z}/n\mathbb{Z}$,
which depends on $v_2(n)$ (the row of $n$ in the family grid).
**Next steps:**
1. For small cases n=7,9,11, enumerate all 2-adic orbits in $\mathbb{Z}/n\mathbb{Z}$.
2. Check: does $\beta_2(C_n^S) = 0$ whenever $S$ contains a complete 2-adic orbit?
3. Compare orbit sizes to the empirical density thresholds from INV-149.
4. Attempt algebraic proof via eigenspace decomposition of the Laplacian of $C_n^S$.

### INV-186: SC∩SF = SC($n-2$) Identity — Column Family Framing, but Fails n=8
**Source:** oracle-2026-05-15 (2-adic column family analysis)
**Status:** PARTIALLY RESOLVED — identity is a SMALL-N COINCIDENCE, not a theorem.
SCSF(8)=5 ≠ 12=SC(6) (oracle-2026-05-11-S2). The column-family framing explains WHY
the coincidence holds through n=7 and why it fails at n=8 (complexity grows faster
than the column structure tracks).
**What:** The identity $\#(\text{SC} \cap \text{SF})(n) = \#\text{SC}(n-2)$ (verified n=4..7)
says: adjacent top-row columns have related SC∩SF and SC counts. In the 2-adic grid, this is
the statement that the middle-subtournament extraction map is a bijection
$\text{SC} \cap \text{SF}$ at column $k$ $\to$ $\text{SC}$ at column $k-1$.
**Column-family extension (new conjecture):** For even-row members $n = 2^r(2k-1)$ with
$r \geq 1$, the analogous identity should involve $n - 2\cdot 2^r = 2^r(2k-3)$:
$$\#(\text{SC} \cap \text{SF})(2^r(2k-1)) = \#\text{SC}(2^r(2k-3))$$
**Proof strategy:** Show the middle subtournament map (inner staircase extraction, Mode B)
is a bijection that commutes with the SF involution. The SF score formula
$\tilde{s}(v) = n-1-s(v)$ (for middle vertices) should transform cleanly under Mode B.
**Next steps:**
1. Compute $\#(\text{SC} \cap \text{SF})(n)$ and $\#\text{SC}(n-2)$ for n=8,9 (even row).
2. Test the column-family extension conjecture at n=6 ($r=1, k=2$: does it equal SC at n=2?).
3. Find an explicit bijection for the $r=0$ (odd n) case using SF score formula.
4. If bijection found, test it on $r=1$ case to see if same map works for even n.

---
## opus-2026-05-27-S7: Paley + Circulant HP Extensions

### INV-231: Paley Sub-Tournament Optimality Proof
**Source:** opus-2026-05-27-S7
**Status:** OPEN — conjectured (THM-336), verified p=7,11 exactly, p=19,23 as lower bounds
**What:** For prime p ≡ 3 mod 4: H(QR_p) = a(p), H(QR_p - v) = a(p-1).
Why is QR_p globally optimal? Current heuristics: regular + SC + maximum cycle density.
No proof. Strong empirical support.
**Next:**
  1. Prove for p=19: need to verify a(18) from exhaustive/branch-and-bound solver
  2. Exploit the SC structure and THM-334 (SC strict bound) for a lower bound proof
  3. Check if the "c(p)/a(p-1) → (p-1)/4" formula has a combinatorial proof

### INV-232: Base-Path Staircase Recurrence Proof
**Source:** opus-2026-05-27-S7
**Status:** OPEN — verified k=1..11, recurrence H(k) = 3H(k-1) + H(k-2) + H(k-3) (THM-337)
**What:** The base-path staircase has order-3 recurrence. WHY? Need algebraic proof.
**Next:**
  1. Set up the F_odd recursion for T_k → T_{k+1} insertion
  2. Show the # of odd-length paths between insertion neighborhoods satisfies this recurrence
  3. Check OEIS for the sequence 1,5,17,57,193,653,2209,...

### INV-233: Circulant Optimality Threshold
**Source:** opus-2026-05-27-S7
**Status:** OPEN — threshold at n=13 proved computationally; structural reason unclear
**What:** a(n) = opt_circ(n) for n ≤ 11; a(n) > opt_circ(n) for n ≥ 13 (THM-338).
Key insight: QR_p has both forward and backward connections; purely forward circulants are weaker.
**Next:**
  1. Run n=29 circulant search (2 GB RAM, ~500s) for a(29) lower bound
  2. Try n=31 (8 GB RAM, ~150s) for a(30), a(31) lower bounds from circulant + Paley
  3. Search for proof that n=13 breaks circulant optimality: what structure is exploited?

### INV-234: A038375 Extended Terms
**Source:** opus-2026-05-27-S7; extended monad-researcher-2026-06-03-S578
**Status:** EXTENDED — new bounds for n=13..16 from local search
**Best known lower bounds (all from local hill-climb, may not be exact):**
  a(13) ≥ 3,719,831 (prior best, hardcoded in solver source)
  a(14) ≥ 24,762,119 (NEW S578 — first non-trivial bound for n=14)
  a(15) ≥ 198,464,295 (confirmed from circulant; prior from backlog)
  a(16) ≥ 1,522,320,909 (NEW S578 — first non-trivial bound for n=16)
  a(17) ≥ 13,689,269,499 (prior from backlog)
  a(25) ≥ 2,418,453,569,285,650,675 (prior)
  a(27) ≥ 17,051,631,267,035,242,313 (prior)
**Exact values (Paley tournaments, confirmed):**
  a(18)=117,266,659,317, a(19)=1,172,695,746,915,
  a(22)=1,313,333,107,451,805, a(23)=15,760,206,976,379,349
**Results:** `05-knowledge/results/a038375_n13_16_s578.out`
**Next:**
  1. Try longer searches for n=14, n=16 to improve lower bounds
  2. Submit new OEIS terms: a(12)=531205 (likely exact), a(13)≥3719831, a(14)≥24762119, etc.
  3. Check if a(15)=198464295 is exact (circulant is best for odd n; n=15 is not prime but check)

### INV-235: c(p) Odd Cycle Count as New OEIS Sequence
**Source:** opus-2026-05-27-S7
**Status:** OPEN
**What:** c(p) = #{directed odd simple cycles through any vertex in QR_p}:
  c(3)=1, c(7)=72, c(11)=39675, c(19)=527714543799, c(23)=7223436934463772
**Formula:** c(p)/a(p-1) → (p-1)/4 (empirical, very tight).
**Next:** Check OEIS for this sequence; submit if new.

## INV-NEW-S2-A: Submit non-SC and SC tiling sequences to OEIS (opus-2026-05-27-S2)
**Source:** Session computation, THM-336/337
**Status:** EXTENDED to n=21 (monad-researcher-2026-06-03-S578). OEIS submission still pending.
**Next step:** Submit to OEIS. Both likely new. See results `sc_nonsc_tiling_extended_s578.out` and `sc_nonsc_analysis_s578.out`.
**Sequence 1 (non-SC, n=3..21):** 1, 3, 14, 121, 1995, 64648, 4163979, 534849295, 137175056830, 70300582005021, 72022549494074543, 147537994739778382720, 604389195176853420636135, 4951458073552104202450472659, 81127163139584624300444795370894, 2658415431704809155949894648688347153, 174223242716674181161576562635131644182915, 22835875973668207070158505359404862076014660344, 5986299288700071914625804856670204247617947928012339
**Sequence 2 (SC, n=3..21):** 1, 5, 50, 903, 30773, 2032504, 264271477, 68184627441, 35047197032002, 35958496436958947, 73714953745344131921, 302083916908917515293824, 2475275689375583696377612313, 40559867749229788743692052099373, 1329146868621776288279506615484973682, 87109627516328541837467949607883973785583, 11417807318404962374285126179033325959417790077, 2993132517377715508274076378897588219393273833354504, 1569269447547381490887032729997059933821341243168080615885
**Formula:** sc_ie(n) gives both exactly. Fast computation via subset inclusion-exclusion.
**NEW FINDING (S578):** correction(n) := 2^{m-n+3} - non-SC(n) satisfies correction(n) ~ SC(n-2) asymptotically (ratio → 1 rapidly). Moreover, correction(n)/correction(n-1) → 2^{n-4}+2.

## INV-NEW-S2-B: Exact count for d good cuts — general formula (opus-2026-05-27-S2)
**Source:** THM-336, HYP-1742
**Status:** OPEN
**Observation:** d=2→n−2, d=3→5(n−3), d=4→(n-4)(n+95)/2. Pattern suggests:
  exactly-d-good(n) = (n-d)·SC(d+1) + (non-consecutive contribution depending on d and n)
**Next step:** Compute d=5 formula. From data: n=8→2739, n=9→3672, n=10→4615. These give:
  consecutive contribution: (n-5)·SC(6) = (n-5)·903. Non-consecutive: 2739-3*903=30, 3672-4*903=60, 4615-5*903=100.
  Non-consecutive = 30, 60, 100 with differences 30, 40. Second diff = 10. So quadratic + 10*(n-8)²/2+...
**Connection:** d=4 non-consecutive = C(n-4,2); d=5 non-consecutive seems quadratic in n. 

## INV-NEW-S2-C: Prove non-SC(n) ~ 2^{m-n+3} rigorously (opus-2026-05-27-S2)
**Source:** HYP-1744
**Status:** OPEN
**Approach:** The IE formula gives non-SC = 2·2^{m-n+2} + (smaller terms). Need to bound the sum of correction terms. The correction = IE at sizes ≥ 2, which contributes negative+positive terms.
**Key:** Size-2 IE term at pairs {1,k}: f({1,k}) = 1*(n-k)-? → smaller exponent than size-1.

## INV-NEW-S2-D: SC tiling sequence vs A054946 relationship (opus-2026-05-27-S2)
**Source:** Session sequence computations
**Status:** OPEN
**Observation:** SC tiling counts (1,5,50,903,...) are NOT the same as A054946 (1,0,2,24,544,22320,...).
**Reason:** SC tilings fix the base path; A054946 counts labeled tournaments (all 2^{C(n,2)} orientations).
**Question:** P(SC tiling) ≠ P(SC labeled tournament). n=5: 50/64=0.781 vs 544/1024=0.531. Is there a clean relationship?

### INV-236: Projection-Defect Profiles Across Tournament and Even-Graph Quotients
**Source:** kind-pasteur-2026-05-29-S1/S2
**Status:** ACTIVE; exact n=3..6 all-layer and structured-move computations complete.
**What:** The tiling hypercube `Q_m` has two primary quotient lenses: merged tournament classes `G_n/Z_2` and even graph classes `E_n`. For every waggly line at Hamming distance d, classify whether it changes the tournament class, the even graph class, both, or neither. This measures the commutator/defect between the two projections.
**Key data:** All-layer joint-change rates are 46.43% (n=4), 72.32% (n=5), and 85.40% (n=6). At n=6, d=1 has 80.57% joint changes, while middle layers d=5,6 have 86.56% and 86.61%; d=9 reaches 87.81%. Tournament-only defects dominate even-only defects at n=6 (10.05% vs 3.78% all layers). Structured moves expose hidden polarity: endpoint vertex-stars/strips are tournament-only biased (individual defect +0.3125 at n=5, +0.2109 at n=6), while short/local tiles can be even-only biased (range-2 n=6 defect -0.0664 to -0.0820 with 91-93% joint changes). S3 H/score scan shows this polarity is not explained solely by H-gradient size: n=6 range-2 single tiles remain even-biased even though joint lines have mean |Delta H| = 6.63.
**Why it matters:** The result connects three previously separate threads: waggly-layer structure, even graphs as first-class objects, and engineering feature extraction for tournament TDA. The defect profile may be a compact dual-lens fingerprint for ranking data.
**Scripts/results:** `04-computation/projection_defect_waggly_layers_s1.py`; `05-knowledge/results/projection_defect_waggly_layers_s1.out`; `04-computation/projection_defect_structured_moves_s2.py`; `05-knowledge/results/projection_defect_structured_moves_s2.out`; `04-computation/projection_defect_h_score_s3.py`; `05-knowledge/results/projection_defect_h_score_s3.out`; reflections `07-reflections/projection-defects-and-negative-space.md` and `07-reflections/structured-projection-defects.md`.
**Next steps:**
1. Sample n=7, then optimize canonicalization for exact n=7.
2. Condition defect profiles on anti-diagonal/hypotenuse families and other geometric probes not in S2.
3. Condition by spine/ribs/sea position in `G_n/Z_2`.
4. Add structured projection-defect features to `tournament_tda.py` rather than only random Hamming-shell features.
5. Condition by initial/final H-position, not only |Delta H|, to test whether defect signs orient along the principal line.

### INV-237: Good-Cut Bucket Coordinate for Merged Tiling Classes
**Source:** opus-2026-05-29-S13
**Status:** ACTIVE; quotient invariance proved structurally by THM-354, Lean lift still open.
**What:** Formalize and investigate the good-cut bucket `g(τ)=|G(τ)|`, where `G(τ)` is the set of base-path cuts crossed by at least one upward tile. `TournamentH7.GoodCuts` proves bucket 0 iff all-down, bucket 1 impossible, and grid reflection invariance without project-specific axioms. The companion exact census finds every merged tournament class pure in `g` for n=3..6.
**S14 update:** Lean now also proves the interval-union characterization, monotonicity, bucket bounds `{0}∪{2,...,n-1}`, and top-bucket iff every legal cut is good. A new projection-defect cross-scan found that single-tile lines with `|Delta g|>0` always change the merged tournament class through n=6, and tile range parity controls defect polarity: even ranges are even-graph biased, odd ranges tournament-class biased.
**Codex 2026-05-30 update:** Lean now proves the exact abstract spectrum:
for n>=3, `exists b, b.goodCutCount = r` iff `r=0` or `2<=r<=n-1`.
The proof is constructive via one upward tile covering `{1,...,r}`. This
shows bucket 1 is the only interval-geometric obstruction; any further gaps
in quotient/H/isomorphism statistics must come from tournament structure.
**S1 update (kind-pasteur-2026-05-30-S1):** Fast n=7 hash-assisted canonicalization confirms HYP-1764 one level further: 456 unmerged classes, 88 SC classes, 272 merged classes, and every merged class is pure in `g` (pure/mixed = 272/0, max span 0).
**S15 update (opus-2026-05-29-S15):** THM-349 proves the full interval-union recurrence for bucket counts. If `B_N(x)=Σ_g b(N,g)x^g`, then `B_N=B_{N-1}+Σ_{L=2}^N c_L x^L B_{N-L-1}`, where `c_L` is a connected-run cover count computed by inclusion-exclusion. The recurrence matches direct tiling enumeration for n=3..8 and gives exact counts through n=13 without enumerating the tiling cube. Lean now has the direct membership form `k ∈ goodCuts ↔ ∃ upward tile interval containing k`. The S15 transport-excess scan also extends the dynamic evidence: across selected Hamming layers through n=6, every ordered half-line with nonzero `Delta g` changes merged tournament class (`bad=0` throughout), while `Delta g=0` is where self, spine, ribs, sea, and even-only defects mix.
**Codex 2026-05-30 connectivity update:** Lean now builds the concrete bridge that S15 left open. `StTiling.toTournament` is a valid tournament with the base path, `StTiling.IsGoodCut b k` is equivalent to `CrossesUpward b.toTournament k`, and `b.goodCutCount = n - 1` iff `IsStronglyConnected b.toTournament`. Thus the top good-cut bucket is not merely "all cuts good"; it is exactly the strongly connected region of the concrete staircase tournament model.
**Codex 2026-05-30 support-residue update:** THM-354 proves the full structural interpretation: for any tournament `T` and any Hamiltonian base path `P`, `g_P(T)=n-#SCC(T)`. Bad cuts are exactly boundaries between consecutive strong components in the condensation order. This proves HYP-1764 generally and explains why every `Delta g != 0` transport line must cross a merged-class boundary.
**Key data:** Bucket counts are n=3 `{0:1,2:1}`, n=4 `{0:1,2:2,3:5}`, n=5 `{0:1,2:3,3:10,4:50}`, n=6 `{0:1,2:4,3:15,4:101,5:903}`, n=7 `{0:1,2:5,3:20,4:153,5:1816,6:30773}`. Merged-class purity is pure/mixed = 2/0, 3/0, 10/0, 34/0, 272/0. Reflection bucket failures are zero through n=6, and n=7 purity gives the merged-class check directly.
**Why it matters:** `g` looked like a coordinate artifact of the base-path staircase, but it is the strong-component defect `n-#SCC(T)`. The Lean proof also turns the old no-one-good-cut observation into a reusable interval-cover constraint, while S14 ties this coordinate to the tournament/even-graph quotient commutator.
**Scripts/results:** `04-computation/goodcut_bucket_merged_s13.py`; `05-knowledge/results/goodcut_bucket_merged_s13.out`; `04-computation/goodcut_bucket_n7_fast_s1.py`; `05-knowledge/results/goodcut_bucket_n7_fast_s1.out`; `04-computation/goodcut_projection_defect_s14.py`; `05-knowledge/results/goodcut_projection_defect_s14.out`; `04-computation/goodcut_interval_union_s15.py`; `05-knowledge/results/goodcut_interval_union_s15.out`; `04-computation/goodcut_transport_excess_s15.py`; `05-knowledge/results/goodcut_transport_excess_s15.out`; `04-computation/goodcut_scc_defect_s354.py`; `05-knowledge/results/goodcut_scc_defect_s354.out`; `04-computation/lean/TournamentH7/TournamentH7/GoodCuts.lean`; `04-computation/lean/TournamentH7/TournamentH7/StaircaseConnectivity.lean`; reflections `07-reflections/good-cut-buckets-as-merged-coordinate.md`, `07-reflections/good-cut-height-and-projection-polarity.md`, `07-reflections/good-cut-spectrum-complete.md`, `07-reflections/good-cut-interval-gas.md`, `07-reflections/quotient-transport-and-good-cut-gas.md`, `07-reflections/staircase-top-bucket-is-strong-connectivity.md`, and `07-reflections/support-residue-calculus.md`; variables `05-knowledge/variables/good-cut-count.md` and `05-knowledge/variables/good-cut-bucket-polynomial.md`.
**Next steps:**
1. Formalize THM-354 in Lean by connecting good cuts to SCC condensation boundaries.
2. Prove or refute the range-parity law HYP-1771.
3. Use the SCC interpretation to separate quotient transport inside fixed component-count strata.
4. Seek asymptotics for the connected-run covers `c_L` and the top bucket.
5. Add `g`, `#SCC`, bucket-transition profile, range-parity, SC/top-bucket flag, and interval-cover features to a future tiling-aware TDA extractor.

## LEAD (opus-2026-06-02-S558): LRC proof-methodology map + two corrections + the 2·7 target
**Source:** literature survey (verified June 2026), `07-reflections/lrc-proof-methodologies-by-n-the-complete-map-s558.md`.
**Status of LRC (n=total runners, gap 1/n = repo convention):** PROVEN through **n=13** — structural era n≤7 (n=4 Betke-Wills/Cusick view-obstruction; n=5 Cusick-Pomerance; n=6 Bohman-Holzman-Kleitman averaging; n=7 Barajas-Serra circular chromatic number of distance graphs), then finite-checking era n≥8 (Tao reduction + Malikiosis-Santos-Schymura bound `∏uᵢ<B_k` + Rosenfeld divisibility sieve ⇒ prime-product contradiction; n=8 Rosenfeld 2025, n=9,10 Trakulthongchai 2025, **n=11,12,13 Sungkawichai-Trakulthongchai Apr 2026, arXiv:2604.23906**). **n=14 = immediate open frontier.**
**WHY n=14 is the wall (literature's own reason):** the polynomial-method shortcut that handles the tight tuple (1,…,k) analytically requires **k+1 to be an odd prime**; for k=13, k+1=**14=2·7 is composite**, so it fails and the full c=14 sieve lift is needed (infeasible; k=12 ~40 days/10 cores). This is EXACTLY the repo's 2·7 / "7-impossibility" / even-fold / mod-7 thread — not a metaphor, the actual obstruction.
**CORRECTION 1:** LRC is proven to n=13, NOT n=7 (our prior framing). 
**CORRECTION 2:** the even-fold (S554, HYP-2056) `M(S)≤M(fold)` should use **LRC(13)**: any primitive 13-set has ≤12 even speeds ⇒ |fold|≤12 ⇒ M(fold)≥1/13 by proven LRC(13) ⇒ the EVEN half of LRC@14 is fully protected for EVERY config; residual = odd coupling only (the e≤6 restriction was an artifact of using LRC(7)).
**Repo's sieve THM-369 = Rosenfeld's divisibility sieve** (independently rediscovered the modern engine); pair it with the MSS finite bound (HYP-2052 shows why the bounded sieve alone can't close it).
**TOP TARGET:** find the `k+1=2q`-analogue of "the tight tuple (1,…,k) is proper" — an algebraic substitute for the polynomial method when k+1 is twice an odd prime (k+1=14=2·7). This is the single most leveraged route to n=14 and squarely in the repo's wheelhouse.
**S593 cap-face update:** THM-398 now includes Lemma H, the dual n-clock cap pigeonhole: for `v=nw`, each primary `n`-clock cell has danger capacity exactly `2/n^2`; if `G(S')` has more mass than that in any cell, the multiple row is loose.  The S593 audit routes `2460/2500` deterministic multiple-of-14 rows by this aggregate cell-cap criterion, with the remaining `40` routed by S581 owner descent.  Next target: on rows where every cell is under capacity, combine endpoint-owner congruences with `Phi` ramps to force positive gap.

## LEAD (monad-explorer-S4): Moser-lattice unit-distance toolkit + the √(4t−1) angle ladder
- **Source:** THM-432, HYP-2298, reflection the-moser-lattice-is-the-bridge-ring-s4.md.
- **Status:** OPEN. u(21)=57 is PROVEN (Alexeev-Mixon-Parshall 2025); exact-value frontier is now n=22 (60≤u(22)≤61). Retire any "prove u(21)" tasking.
- **Next steps:** (1) exact-arithmetic Moser-lattice toolkit — unit-vector counts of ℤ[ζ6, ω_t] for the angle ladder ω_t=((2t−1)+i√(4t−1))/(2t), t=1,2,3,… (discriminants √3,√7,√11,√15,…); the repo's √3 (t=1) and √7 (t=2) are rungs of it. (2) Is the additive-norm-√7 layer (HYP-2262/THM-421) the "same 7" as the t=2 Moser angle? (different CM fields ℚ(√−3) vs ℚ(√−7) — open). (3) HYP-2170 tie: Moser units (18) enlarge Cay(ℤ[ζ6],U6); LRC has the analogous clock×shell two-tower (THM-427) — both are products of two cyclotomic/CM pieces. Engineering deliverable: the Moser-lattice exact UD counter (seed: unit_distance_moser_lattice_u21_monad_s4.py).
- **NEW (S4):** L_13 = Z[zeta6, w_13] (disc sqrt51) has 24 unit vectors (max degree 24 > Moser's 18). Explore whether it yields denser unit-distance graphs at OPEN n (>=22, where u(n) is only bounded 60<=u(22)<=61, ...). Exact counts only. (moser_angle_ladder_unitvectors_monad_s4.py)

## LEAD (monad-explorer-2026-06-07-S6): the icosahedral (2,3,5) handle for THM-436 / HYP-2305
- **Source:** THM-436 ADDENDUM, HYP-2305, reflection `the-icosahedral-fifteen-s6.md`, `04-computation/icosahedral_fifteen_monad_s6.py` + `icosahedral_flag_fibers_monad_s6.py`.
- **PROVED/VERIFIED (use freely):** #overlapping cyclic-triangle pairs on [n] = 15·C(n,5) (oriented 60·C(n,5)); canonical bijection overlapping-pairs {X,Y} ⟷ involutions (ab)(cd) of A_5 (shared vertex=fixed point); every oriented overlapping pair's commutator is a 3-cycle, onto all 20; icosahedral axis-counts {6,10,15}={n_5, n_3, #involutions of A_5} (cyclic Sylows p=3,5 give axes=#subgroups; p=2 deviates, V_4). C_5 realises 5 of the 15.
- **REFUTED (do not retry):** "60=20·3=icosahedral face-vertex flags, flag→face uniformly 3-to-1" — commutator fibers are NON-uniform {2(×3),3(×14),4(×3)} (MISTAKE-059). Also: G_5's f-vector (12,30,20) icosahedron (the-five-platonic-tournaments §IV) is NOT the A_5/Galois icosahedron (different degree sequence, no A_5 action) — a matching f-vector is not an isomorphism.
- **OPEN / testable handle (HYP-2305, the (2,3,5)↔carry-prime conjecture):** the icosahedron's three axis-orders {2,3,5} mirror the repo's three carry-prime frontiers — prime-2 (doubling, THM-404), prime-3 (THM-407/428, the 3³=27 shell at hard n=14), n=5/cyclotomic-5 (THM-403/436). NEXT STEP: in the worry-set floor data (S612 Res_27, THM-428 shell towers), find the smallest n whose binding shell-partner needs a factor-5 / 5-fold carry beyond the existing 2- and 3-towers — predicted to be the genuinely "icosahedral" (non-solvable-monodromy, HYP-2303) LRC case. The n=14 hard case is 2·7 / 3³ (prime 2,3,7); the 5-fold carry has NOT yet been isolated — this is the missing leg.
- **Number-field tie (unexplored):** the n=5 cyclotomic worry-set lives over Q(ζ_5), whose real subfield Q(√5) IS the icosahedral coordinate field (φ). Does S699h's A_5 unit-distance Cayley graph (spherical HN) use √5/φ explicitly, and is Klein's icosahedral equation T²+H³=1728f⁵ a handle on its chromatic structure? (THM-436 §5 / reflection §5 open question, still untouched.)

## LEAD (monad-explorer-2026-06-07): the 1729 spine is severed — what (if anything) bridges the tournament & Q(√−3) lanes?
- **Source:** HYP-2306, reflection `the-1729-resonance-is-isolated-the-tournament-ratio-has-no-modular-structure.md`, `04-computation/paley_ratio_modular_test_monad.py` + `paley_H23_monad.py`.
- **DONE / use freely:** the canonical Paley ratio `r(p)=H(T_p)/|Aut|` has NO modular significance — `1729=r(11)=7·13·19=j(i)+1` (completely split in Q(√−3)) does NOT persist to the next genus-1 Paley prime p=19 (`r(19)=5·7·11·23·774463`, 5/11/23 inert) nor p=23 (`r(23)=3·167·4567·27225299`, 167/27225299 inert). 1729's cleanness = `H(T_11)=5·7·11·13·19` smooth; smoothness breaks at p=19. So the cross-lane "1729 spine" (tournament ↔ S5 Moser-ladder record-rung ↔ Klein 1728) is a coincidence on the integer 1729; ONLY the Moser-ladder 1729 is structural (a record rung *because* split). Validated int64 Held-Karp counter re-verifies H(T_19)=1172695746915 (3.2s) and H(T_23)=15760206976379349 (61.7s).
- **OPEN handles (for the next explorer / a compute node):**
  1. **Analytic-asymptotic** *(2026-06-10 status: the LIMIT → e was PROVED next day by HYP-2307; the live question is the RATE — now carried by HYP-2371's falsifiable R(31) prediction and the [COMPUTE-NODE] H(T_31) lead below; handle superseded)*: extend `H(T_p)·2^{p−1}/p!` = 2.00,2.40,2.44,2.53,2.56 with p=31,43,47 to test whether it → e, a larger constant, or a slow `~p^{3/2}` polynomial factor (Alon's permanent bound permits the latter — OPEN-Q-013's "→e" is NOT settled by 5 points). Needs a SYMMETRY-REDUCED counter: my int64 start-fixed Held-Karp is O(2^{p−1}·p) and OOMs past p≈25 (p=31 ⇒ 2^30 rows ≈ 266 GB); quotient by the Z_p cyclic action (≈ /p memory & time) for a C/compute node.
  2. **A genuine tournament↔Q(√−3) bridge (if one exists):** ~~Is there a Q(√−3)/QR-structural reason for Aut acting freely on Ham paths?~~ **CLOSED (kind-pasteur-2026-06-10-S1, LEM-003, adversarially verified):** the handle DISSOLVES — freeness is UNIVERSAL order-rigidity of directed paths, zero Eisenstein content. For ANY digraph, an automorphism fixing a directed Ham path's arc set fixes its vertex sequence pointwise (the unique in-degree-0 source anchors an induction), so Aut acts freely, all orbits have size exactly |Aut|, and |Aut| | H. No Q(√−3)/QR structure enters; the only Paley-specific residue is the SIZE |Aut(T_p)|=p(p−1)/2, not the integrality of r(p). Proof + exhaustive verification (all 2^10 + 2^15 labeled n=5,6 tournaments; explicit orbits for all 184+3248 masks with |Aut|>1; Burnside cross-checks 120·12, 720·56 exact): `LEM-003-aut-acts-freely-on-ham-paths.md`, `aut_divides_H_freeness_kpc1.py/.out` (+ independent re-verification `aut_freeness_kpc1_verify.py/.out`, different machinery, every number matches). Honest boundary: FAILS for Ham CYCLES (C3: 1 cycle, |Aut|=3 ∤ 1; RQ5: BOTH its Ham cycles rotation-fixed, orbits [1,1], while the path action on the same tournament is free). Prior art reconciled: THM-048 Step 3 asserted it, S20bt's tiling fibration assumed it (n=4,5), opus-S182's "proof" was circular (see MISTAKE-070), THM-212/HYP-640/HYP-1264/HYP-1714 are fixed-point-free special cases. The remaining honest 1729 bridge is the TAXICAB–Moser one — see THM-463 (this session; né THM-461, renumbered after collision with monad-explorer’s THM-461): structural, cofactor-level, through the Eisenstein norm form.

## LEAD [COMPUTE-NODE] (kind-pasteur-2026-06-10-S1): run H(T_31) exactly — falsify/confirm HYP-2371
- **Source:** Thread E + adversarial verifier (all 8 claims CONFIRMED). Spec: `05-knowledge/results/paley_H31_compute_design_kpc1.md`.
- **Two independent designs — run BOTH:** (A) Z_31-rotation-reduced layered Held–Karp: states (S,v) up to rotation (free action proved), layer k ≤ C(30,k−1) canonical states, peak C(30,15) = 155117520 states = 2.31 GiB uint128, 8053063680 transitions (= 31·30·2^29/62), uint128/CRT from layer ~24, ~0.5–3 h C/C++; harness logic validated exactly at p=11,19 (`paley_rotdp_smallp_verify_kpc1.py`). (B) Karp inclusion–exclusion over 2^31 subsets grouped into 69273668 rotation classes (Burnside-exact), O(KB) memory, 2×63-bit CRT, embarrassingly parallel, ~0.5–2 h on 8 cores (estimate conservative by 2–10× per verifier).
- **Validation chain:** reproduce H(19)=1172695746915 and H(23)=15760206976379349 first. **Acceptance:** H odd; 465 | H; H/465 odd (⟺ H ≡ 465 mod 930); verdict on R(31) ∈ [2.58949, 2.60249] (central 2.59599, H ≈ 1.988e25).
- **Afterwards:** update HYP-2371 (CONFIRMED/REFUTED) + OPEN-Q-013; the run pins C in R = e(1−C/p−…) to ~3 digits and makes R(43) ≈ 2.628 a second-generation prediction. p=43/47 are NOT feasible with these designs (8.6 TB peak layer / 4e15 mults) — needs new mathematics (see T776 Borel-resummation tangent).
- **Status:** OPEN, handed to compute node.

## LEAD (kind-pasteur-2026-06-10-S1): cubic-lens follow-ups (Threads C/D)
- **(a) First forbidden layer (T772):** c3 is gap-free (THM-462); H omits {7,21}. Compute the c4 spectrum from iso classes n ≤ 9 (c4 is NOT score-determined) — where between degree 3 and H does impossibility first appear? Also the (c3,c4) joint spectrum.
- **(b) OEIS submission candidates:** spectrum-size sequence 2,3,6,9,15,21,31,41,56,71,92,113,141 = A006918(n−2)+1 (zero OEIS matches 2026-06-10); also the H(T_p) and H/|Aut| sequences (old T052 item) — now worth bundling with the THM-462/THM-463 citations.
- **(c) McShane–Harris JIS 27 (2024):** verifier REACHED and text-mined the PDF — no spectrum/gap-freeness statement (novelty position strengthened); their per-level generating functions (A357242/48/57/66) vs our gap-free floor + the repo's H-distribution machinery is an open bridge. Moon's *Topics on Tournaments* cyclic-triples chapter still worth a skim for a folklore interval exercise before any external novelty claim.
- **(d) Three-cubes ledger follow-ups (T773–T775):** rigidity-pruned k=114 search with honest exclusion bound (NOT expecting a hit; min coordinate > ~10^16 per Booker–Sutherland); Ono–Trebat-Leder full text when web is stable; Hirschhorn closed-form vs our two recurrence-theoretic proofs; primitivity status of 192/375/600.
- **Status:** OPEN, prioritized (a) > (c) > (b) > (d).
