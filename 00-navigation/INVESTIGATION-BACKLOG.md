# Investigation Backlog

**Purpose:** Systematic catalog of every lead, reference, connection, and unexplored direction extracted from the repo. Claude agents should consult this before choosing what to work on, and add new leads as they emerge. Prioritized by potential impact on proving OCF (Claim A).

**Last full repo scour:** opus-2026-03-06-S10
**Last web research:** kind-pasteur-2026-03-10-S50 (Tang-Yau circulant path homology, Schweser-Stiebitz-Toft Redei revisited, Mitrovic NC Redei, Satake DRT)
**Last engineering update:** kind-pasteur-2026-03-10-S54 (sparse T_19 Omega computation, CLAUDE.md engineering mandate)
**Last packing framework update:** opus-2026-03-14-S71f (nesting obstruction, (z-2)(z-3) recurrence, 2-Bridge)
**Last S90 update:** opus-2026-03-15-S90 (simplicial RÃƒÂ©dei, Cayley monad, Ãâ€ž-Ãâ€  clock, equidecomposability)
**Last gauge theory analysis:** kind-pasteur-2026-03-21-S12 (Napolitano paper, Cartan decomposition bridge, TournamentProbe)
**Last literature sweep:** opus-2026-04-05-S24 (Stanley-Stembridge proved, Mitrovic NC deletion-contraction, Tang-Yau circulant Fourier)
**Last LRC web trawl:** codex-2026-06-26-S170 (2-adic Littlewood/Hurwitz finite transducer, Ostrowski-Hadamard gap warning, ACM large-stick safe-component analogy; exact-gap companion to HYP-3008/LTI-158 and HYP-3009/LTI-159 at HYP-3010/T1094)
**Last full repo scour:** opus-2026-03-05-S9
**Last web research:** opus-2026-03-05-S9 (Paley maximizer, n=8 anomaly)

---

## Lead codex-2026-06-29: small-touch singleton-pocket geometry atlas

**Status:** EVIDENCE / finite geometry atlas; not proof
(HYP-3478/T1438/LTI-438/LTT-338), downstream of HYP-3476 pair-current,
HYP-3477 hard-mirror discharge, and the S319 colored gate unit-delta split.
**Status:** EVIDENCE / exact six-row geometry audit; not proof
(HYP-3478/T1438/LTI-438/LTT-338), downstream of HYP-3476 and HYP-3477.

HYP-3478 audits the six small-touch/no-hard rows directly:
`random_covering_001`, `random_covering_039`, `random_covering_062`,
`random_covering_074`, `random_covering_086`, and `random_covering_101`.
All six have edgeless dead-cover projections, all dead components are isolated
singleton pockets, and every pocket has an exact mirror partner plus at least
two complete E/branch gate touches.  The owner-pair grammar is:

```text
001: (165,179), (81,153)
039: (63,129)
062: (9,81)
074: (15,99)
086: (133,169)
101: (7,175)
```

Next proof task: prove the four branch-unit singleton rows (`001`, `062`,
`086`, `101`) by a component-local mirror-pocket theorem, then add the S319
cover-delta sidecar for `039` and `074`.  Do not spend the next session on
larger pair-current searches unless a new coordinate creates nontrivial
component support.
Readout: all six rows have edgeless dead-cover projections, but all six have
touching E/branch gates.  Across the packet there are `14` dead components,
all singleton `B0`/`B1` owner-pair components; every component has an exact
interval mirror with the branch owners swapped.  There are no singleton-cover,
mirror, or owner-balance failures.

The split is now explicit: cover-delta sidecar rows
`random_covering_039`, `random_covering_074`; clean unit-delta singleton rows
`random_covering_001`, `random_covering_062`, `random_covering_086`, and
`random_covering_101`.  `random_covering_101` has a shortest E/branch gate
that does not touch dead labels, but it has a slightly longer touching
unit-delta gate, so it remains in the unit-delta singleton packet.

Next: prove a finite mirror-singleton current lemma, likely split into the
two cover-delta rows and four unit-delta rows.  Use tournament vertices as
singleton components, mirror-paired components, owner-pair labels, fixed
section boundaries, touching gate events, residues, cover arcs, color/Fourier
modes, and proof obligations.  The quotient must preserve the terminal
discharge predicate and say what interval/branch/owner geometry it destroys.

## Lead codex-2026-06-29: zero-edge singleton-current audit

**Status:** RESERVED STUB / computation in progress; not proof
(HYP-3480/T1440/LTI-440/LTT-340), downstream of incoming HYP-3478 small-touch
geometry, HYP-3479 hard-orbit/current join, HYP-3476 pair-current/router, and
HYP-3477 hard-orbit discharge.

Claimed files:
`05-knowledge/hypotheses/HYP-3480-lrc14-zero-edge-singleton-current.md`,
script
`04-computation/lrc14_zero_edge_singleton_current_codex_20260629.py`,
result
`05-knowledge/results/lrc14_zero_edge_singleton_current_codex_20260629.out`,
and reflection
`07-reflections/lrc14-zero-edge-singleton-current-codex-20260629.md`.

Target: the six small-touch/no-hard zero-edge rows
`random_covering_001`, `random_covering_039`, `random_covering_062`,
`random_covering_074`, `random_covering_086`, and `random_covering_101`,
with `random_covering_031` retained as the hard/currentless control.  Test
whether the terminal carrier is owner-label involution, unit-delta singleton
current, endpoint-spine, two-adic residue, signed-SPEC debt, or named gluing.

## Lead codex-2026-06-29: hard mirror-orbit discharge audit

**Status:** EVIDENCE / exact hard-family ledger; not proof
(HYP-3477/T1437/LTI-437/LTT-337), downstream of HYP-3475 and the active
HYP-3476 pair-current and exception-frontier packets.

HYP-3477 audits the eight cover-delta `>=7` mirror orbits from HYP-3475.
Result: `7/8` have a lower-delta E/branch survivor gate that gives a
dead-cover projection-edge cut; the sole exception is `random_covering_031`,
whose max-delta mirror pair has zero compatible q=`14V` phase-grid hits and
is the named HYP-3455/HYP-3476 gluing clause.  The unique delta-`8` orbit on
`random_covering_022` is not separate debt.

Next proof task: prove the lower-delta E/branch projection-current discharge
for non-random031 hard orbits, and handle random031 through the seven-owner
gluing / pair-current / owner-current / two-adic / signed-SPEC exits.

## Lead codex-2026-06-29: dead-cover pair-current exception audit

**Status:** EVIDENCE / targeted two-gate graph-current audit; not proof
(HYP-3476/T1436/LTI-436/LTT-336), the two-gate follow-up to HYP-3472's
single-gate boundary-current exceptions and HYP-3475's mirror-orbit quotient.

Artifacts:
`05-knowledge/hypotheses/HYP-3476-lrc14-dead-cover-pair-current.md`,
script
`04-computation/lrc14_dead_cover_pair_current_codex_20260629.py`,
result
`05-knowledge/results/lrc14_dead_cover_pair_current_codex_20260629.out`,
and reflection
`07-reflections/lrc14-dead-cover-pair-current-codex-20260629.md`.

Readout: on the HYP-3472 edge-cut exceptions
`random_covering_001`, `random_covering_031`, `random_covering_039`,
`random_covering_062`, `random_covering_074`, `random_covering_086`, and
`random_covering_101`, plus separating-only exceptions `covering_AP_with_84`
and `ap_omit_12_tail_84x01`, the script enumerates `9113` low-rank E/branch
gate pairs and `165` exact mirror pairs.  Pair current closes the two AP rows
but none of the seven random rows.  The reason is structural: the random
exception dead-cover projections are edgeless singleton projections
(`edge_support_label_count=0`), so no adjacent-label deletion can cut them.
The AP rows have connected `90`-edge projections and separate under selected
two-gate pairs.

Tournament vertices should be proof carriers and obligations: pair-current
cut, mirror-pair current, single-gate current, HYP-3473 formal packet,
partition-lattice guardrail, and scalar shadows.  Record explicitly which
predicate is preserved and which interval/branch information is destroyed.

## Lead codex-2026-06-29: exception-frontier router sidecar

**Status:** EVIDENCE / finite packet-router audit; not proof
(HYP-3476/T1436/LTI-436/LTT-336), joining HYP-3472 boundary-current
exceptions with HYP-3475 hard mirror-orbit debt.

Claimed files:
`05-knowledge/hypotheses/HYP-3476-lrc14-exception-frontier-router.md`,
script
`04-computation/lrc14_exception_frontier_router_codex_20260629.py`,
result
`05-knowledge/results/lrc14_exception_frontier_router_codex_20260629.out`,
and reflection
`07-reflections/lrc14-exception-frontier-router-codex-20260629.md`.

Readout: the currentless boundary frontier and hard mirror-orbit frontier
intersect only at `random_covering_031`.  Six hard rows already have
separating currents; six non-AP currentless rows have no hard orbit and best
touching-gate adjacent delta in `{2,3}`; the two AP84 currentless rows route to
the closed AP84 packet.  No subset of HYP-3474 axes `K,N,T,S,F,C,M,A`
preserves the six-label terminal route partition, so sidecar `R` is required
unless reconstructed.

Next: prove the six small-touch/no-hard rows are harmless by bounded
owner-current, endpoint-spine, two-adic, exact-period, or signed-SPEC debt:
`random_covering_001`, `random_covering_039`, `random_covering_062`,
`random_covering_074`, `random_covering_086`, and `random_covering_101`.
Use tournament vertices as terminal proof packets, not runners or raw row
names.
Pair-current connection: add a zero-edge singleton-current terminal packet for the seven random
rows.  Candidate carriers are isolated dead-component owner-current imbalance,
HYP-3455/HYP-3451 gluing/conductance, HYP-3460 phase-branch bypass, two-adic
descent, signed-SPEC/Rprime, or state-lift debt.  Do not keep trying larger
E/branch gate pairs unless they add a new edge-support label coordinate.
Incoming connections: HYP-3477 claims the hard mirror-orbit lane and intersects
this zero-edge set at `random_covering_031`; the colored gate unit-delta split
intersects it at cover-delta sidecar rows `random_covering_039` and
`random_covering_074`.  Test `random_covering_001`, `random_covering_062`,
`random_covering_086`, and `random_covering_101` first for the pure
unit-delta singleton-current packet.

## Lead codex-2026-06-29: dead-cover boundary-current audit

**Status:** EVIDENCE / exact graph-current audit; not proof
(HYP-3472/T1432/LTI-432/LTT-332), the graph/current sibling to HYP-3471 after
HYP-3462 closed the AP84 corridor-splice packet.

Claimed files:
`05-knowledge/hypotheses/HYP-3472-lrc14-dead-cover-boundary-current.md`,
script
`04-computation/lrc14_dead_cover_boundary_current_codex_20260629.py`,
result
`05-knowledge/results/lrc14_dead_cover_boundary_current_codex_20260629.out`,
and reflection
`07-reflections/lrc14-dead-cover-boundary-current-codex-20260629.md`.

Readout: `130/130` dead rows have a low-rank E/branch gate touching the
dead-cover projection; `123/130` have a projection-edge cut; `121/130` have a
separating current.  Edge-cut exceptions are
`random_covering_001`, `random_covering_031`, `random_covering_039`,
`random_covering_062`, `random_covering_074`, `random_covering_086`, and
`random_covering_101`.  Separating-current exceptions add
`covering_AP_with_84` and `ap_omit_12_tail_84x01`.

Next: run the boundary-current audit.  For each low-rank E/branch gate, record
adjacent blocker labels, remove them from the dead-cover incidence graph, and
report projection-cut strength, separating-current failures, branch-current
balance, and AP84/non-AP exceptions.  Tournament vertices should be proof
carriers and obligations, not runners: projection cut, separating current,
E/branch implication, AP84 closed packet, random031 seven-owner gluing clause,
typed gate word, raw gate count, and raw dead fraction.
## Lead codex-2026-06-29: colored gate formalization interface

**Status:** FORMALIZATION / Lean interface and conditional assembly; not proof
(HYP-3473/T1433/LTI-433/LTT-333), turning HYP-3471's colored gate-reservoir
theorem target into a sorry-free formal packet.

Claimed files:
`05-knowledge/hypotheses/HYP-3473-lrc14-colored-gate-formalization.md`,
Lean module
`04-computation/lean/TournamentH7/TournamentH7/LRCColoredGateFormalization.lean`,
result
`05-knowledge/results/lrc14_colored_gate_formalization_lean_codex_20260629.out`,
and reflection
`07-reflections/lrc14-colored-gate-formalization-codex-20260629.md`.

Readout: the Lean interface defines colored gate endpoint kinds, typed
endpoint residues, low-rank E/branch gates, `DeadCoverEBranchSoundness`, exact
HYP-3471 count ledger, terminal exits, `ColoredGateGlobalCoverage`, and proves
the conditional assembly `ColoredGateGlobalCoverage -> LRC14Statement`.  The
carrier ledger uses proof obligations and terminal packets as vertices rather
than runners, arcs, residues, or raw gate counts.

Next: instantiate `DeadCoverEBranchSoundness` for the HYP-3438/HYP-3453 exact
gate bank.  The strongest route remains a HYP-3451 current/cut proof: a dead
component with no low-rank E/branch leak should force impossible branch-current
divergence or a Menger cut with typed gate-boundary terminals.  Then fill the
legal terminal packets through HYP-3462, HYP-3470, HYP-3461, HYP-3460,
HYP-3459, HYP-3458, HYP-3455, HYP-3451, two-adic descent, or signed SPEC.
Next: prove the universal touch lemma, then edge-cut/separating-current
transfer away from the named exceptions.  Route `random_covering_031` through
HYP-3455, the AP84 base rows through HYP-3462/HYP-3470, and the six remaining
random exceptions through owner-current, two-adic, signed-SPEC, or state-lift
debt.  Tournament vertices were proof carriers and obligations, not runners:
projection cut, separating current, E/branch implication, AP84 closed packet,
random031 seven-owner gluing clause, typed gate word, raw gate count, and raw
dead fraction.
## Lead codex-2026-06-29: Colored gate mirror-orbit audit

**Status:** EVIDENCE / exact mirror-orbit quotient audit; not proof
(HYP-3475/T1435/LTI-435/LTT-335), instantiating HYP-3461's colored
gate-extension orbit schema on top of HYP-3471's typed gate reservoir, as the
mirror-orbit sibling to HYP-3472's boundary-current audit, HYP-3473's Lean
formalization interface, and HYP-3474's partition-lattice quotient guardrail.

Claimed files:
`05-knowledge/hypotheses/HYP-3475-lrc14-colored-gate-mirror-orbit.md`,
script
`04-computation/lrc14_colored_gate_mirror_orbit_codex_20260629.py`,
result
`05-knowledge/results/lrc14_colored_gate_mirror_orbit_codex_20260629.out`,
and reflection
`07-reflections/lrc14-colored-gate-mirror-orbit-codex-20260629.md`.

Readout: the full HYP-3438/HYP-3453 bank has `8702` survivor gates forming
`4351` exact two-gate mirror orbits, with `fixed_orbits=0` and
`unpaired_or_duplicate_gates=0`.  The quotient keeps HYP-3471's dead-row
escape target: `dead_rows_with_e_branch_low_rank_orbit=130/130`.  Typed gate
colors compress to `186` mirror-orbit colors, structural sidecars to `82`, and
the AP84 four-color packet compresses to two mirror colors.

Next: prove the finite hard-orbit discharge for the `8` mirror orbits with
cover-delta `>=7`, on rows `random_covering_022`, `random_covering_031`,
`random_covering_049`, `random_covering_078`, `random_covering_080`,
`random_covering_085`, and `random_covering_113`.  Route them through
HYP-3455/HYP-3451, endpoint-spine/wall, owner-current, two-adic descent,
HYP-3460 phase-branch bypass, or signed-SPEC/Rprime debt.

## Lead codex-2026-06-29: Hard orbit / boundary-current join

**Status:** EVIDENCE / finite ledger join; not proof
(HYP-3479/T1439/LTI-439/LTT-339), joining HYP-3475's hard mirror-orbit debt
with HYP-3472's boundary-current exception ledger.  This is the Lean-backed
ledger/check for the incoming HYP-3476 exception-frontier router and HYP-3477
hard-orbit discharge split.

Claimed files:
`05-knowledge/hypotheses/HYP-3479-lrc14-hard-orbit-current-join-lean-ledger.md`,
script
`04-computation/lrc14_hard_orbit_current_join_codex_20260629.py`,
result
`05-knowledge/results/lrc14_hard_orbit_current_join_codex_20260629.out`,
Lean ledger
`04-computation/lean/TournamentH7/TournamentH7/LRCHardOrbitCurrentJoin.lean`,
Lean result
`05-knowledge/results/lrc14_hard_orbit_current_join_lean_codex_20260629.out`,
and reflection
`07-reflections/lrc14-hard-orbit-current-join-codex-20260629.md`.

Readout: the hard orbit family has `8` delta-`>=7` orbits on seven random
rows.  Six of seven hard rows have projection-edge cuts and separating
currents, so `7/8` hard orbits already have a separating E/branch current
exit.  The only hard current exception is `random_covering_031`; AP84 has no
hard orbit rows.

Next: prove separating-current transfer.  Then hard-orbit discharge reduces to
that transfer plus the HYP-3455/HYP-3460 `random_covering_031` clause, while
AP84 and the six other random current exceptions remain separate non-hard
sidecar debt.

Lean status: `TournamentH7.LRCHardOrbitCurrentJoin` builds and proves the
finite dispatch arithmetic; the geometric separating-current transfer remains
the open producer.

## Lead codex-2026-06-29: AP84 coloring-recursion bridge

**Status:** EVIDENCE / exact coloring-recursion sidecar for AP-tail bridge; not
proof (HYP-3458/T1418/LTI-418/LTT-318), reconnecting HYP-2247/HYP-2243
coloring recursion to the AP84 endpoint-clock packet.

Claimed files:
`05-knowledge/hypotheses/HYP-3458-lrc14-ap84-coloring-recursion-bridge.md`,
script
`04-computation/lrc14_ap84_coloring_recursion_bridge_codex_20260629.py`,
result
`05-knowledge/results/lrc14_ap84_coloring_recursion_bridge_codex_20260629.out`,
and reflection
`07-reflections/lrc14-ap84-coloring-recursion-bridge-codex-20260629.md`.

Readout: for `{1,2,...,11,13,84m}`, HYP-3456's left-corridor count
`N(m)=floor((504m-6)/70)-floor((96m-13)/14)` becomes a `35`-state boundary
color after subtracting `floor(12m/35)`.  The HYP-3433/HYP-3454 endpoint
address `a_m=ceil(48m/7)` is always in the first two left-corridor high gaps,
with exact rank `1` iff `7|m`, otherwise `2`.  Validation through `m=350` has
`period_failures=[]`, `rank_failures=[]`, `phase_failures=[]`, and
`address_failures=[]`; boundary histogram `{0:1,1:22,2:12}`, endpoint-rank
histogram `{1:5,2:30}`, transition histogram `{0:24,1:10,2:1}`, and mixed
phase support only `m=1..4`.

Next: splice HYP-3454/HYP-3456/HYP-3457/HYP-3458 into HYP-3439 as one AP-tail
rank-`5` descent packet, and run a sibling coloring-rank scout on HYP-3438
survivor gates joined to HYP-3453 component/gate escapes.
## Lead codex-2026-06-29: phase-branch color pullback

**Status:** EVIDENCE / exact color-pullback certificate; not proof
(HYP-3460/T1420/LTI-420/LTT-320), reconnecting S359/S363 regular circular
colorings and HYP-2593/HYP-2595 phase-color CRT grids with the
HYP-3438/HYP-3450 branch-coloured gate atlas.  HYP-3460 is the noncanonical
random031 sibling of incoming HYP-3458, which handles AP84 coloring recursion,
and HYP-3459, which adds the AP84 color-packet legality guardrail.

Claimed files:
`05-knowledge/hypotheses/HYP-3460-lrc14-phase-branch-color-pullback.md`,
script
`04-computation/lrc14_phase_branch_color_pullback_codex_20260629.py`,
result
`05-knowledge/results/lrc14_phase_branch_color_pullback_codex_20260629.out`,
and reflection
`07-reflections/lrc14-phase-branch-color-pullback-codex-20260629.md`.

Readout: for `random_covering_031` as an S3 row with `V=173`, `P=(12,)`,
and `E=(0,4,26,53,60,80,89,103,115,118,128,150)`, the phase-color CRT layer
has `Sigma≈1.831959`, `K=2254`, `actual=open=282`, and exact deficit
`≈34.928978`, below HYP-2595's candidate bound `8*(k+cGP)+1=193`.  Pulling
actual witnesses back by `u=2t mod1` gives mirror-symmetric phase/branch
counts, `no_component_hits=0`, `242` branch-compatible gate hits, and
`hard_gate_hits={}` for HYP-3455's two max-delta mirror gates.  Hard components
`43` and `54` are touched only through lower-delta opposite-branch gates.

Next: prove the phase-branch bypass lemma.  Any branch-colored max-delta gate
obstruction with zero compatible phase-grid hits must discharge by colored
resonance cancellation, low-rank component escape, endpoint-spine/wall lift,
owner-current imbalance, two-adic descent, or signed-SPEC/Rprime debt.  Extend
the script to the full HYP-3438/HYP-3450 bank and classify max-delta gate hits
versus bypasses.
## Lead codex-2026-06-29: Colored-extension gate carrier

**Status:** SYNTHESIS / executable carrier scout; not proof
(HYP-3461/T1421/LTI-421/LTT-321), reconnecting prior coloring and extension
work to the current LRC14 covering-floor route, with HYP-3458, HYP-3459, and HYP-3460 as
the incoming AP84 coloring-recursion, color-legality, and phase-branch pullback sidecars.

Claimed files:
`05-knowledge/hypotheses/HYP-3461-lrc14-colored-extension-gate-carrier.md`,
script
`04-computation/lrc14_colored_extension_gate_carrier_codex_20260629.py`,
result
`05-knowledge/results/lrc14_colored_extension_gate_carrier_codex_20260629.out`,
and reflection
`07-reflections/lrc14-colored-extension-gate-carrier-codex-20260629.md`.

Readout: the scout's proof-carrier tournament is transitive, with path
`random031_mirror_colored_extension_clause -> AP84_colored_endpoint_floor_packet
-> survivor_gate_colored_payload -> observer_cut_orbit_ledger ->
minimal_two_color_bad_core -> colored_resonance_half_boundary ->
A000568_edge_envelope_controlled_forgetting -> PH_bad_coloring_outer_extension_rank
-> metagraph_GF2_color_boundary -> phase_color_reservoir ->
raw_gate_word_shadow -> raw_component_count`.

Next: instantiate the colored gate-extension orbit schema on actual HYP-3438
survivor gates.  At minimum audit canonical AP84 rows, `random_covering_031`,
and the HYP-3441 negative naive-slack rows, retaining gate word, branch mask,
endpoint walls, minimal B0/B1 owner covers, cover deltas, mirror orbit id,
low-rank escape id, and discharge mode.  Desired outcome: every large orbit
has a rank-2 escape, while the few unescaped orbits are exactly the AP84 and
random031 packets already named.
## Lead codex-2026-06-29: AP84 corridor-splice certificate

**Status:** EVIDENCE / AP-tail bridge splice; not proof
(HYP-3462/T1422/LTI-422/LTT-322), closing the HYP-3439 AP84 handoff by
importing HYP-3431 as the complete low branch-union carrier and splicing
HYP-3454/HYP-3456/HYP-3457 into the rank-`5` descent.  This is the canonical
AP84 structural sibling to HYP-3460's noncanonical phase-branch color pullback.

Claimed files:
`05-knowledge/hypotheses/HYP-3462-lrc14-ap84-corridor-splice-certificate.md`,
script
`04-computation/lrc14_ap84_corridor_splice_certificate_codex_20260629.py`,
result
`05-knowledge/results/lrc14_ap84_corridor_splice_certificate_codex_20260629.out`,
and reflection
`07-reflections/lrc14-ap84-corridor-splice-certificate-codex-20260629.md`.

Readout: the HYP-3431 low branch-union carrier is exactly
`[8/49,6/35] union [29/35,41/49]`; branch-specific `b0_good` and `b1_good`
overlap inside those two corridors.  On `{1,2,...,11,13,84m}`, the one-branch
rescue split through `m=70` is `rescue_rank_hist={5:69,6:1}` with
`rank_split_failures=[]`: `m=1` has rank-`6` core `(3,5,7,9,11,13)`, while
all checked `m>=2` rows have rank-`5` core `(5,7,9,11,13)`.

Splice checks: HYP-3457 finite exact/closure/rank-drop failures are empty;
HYP-3454 endpoint failures are empty on checked `m=5..70` and symbolic
containment through `m=420`; HYP-3456 mirror/formula/component-audit failures
are empty through `m=70`, with no shift failures through `m=210`.

Next: AP84 is now closed as a bridge packet.  Feed it into
HYP-3460/HYP-3453/HYP-3451 and discharge the non-AP transfer, especially
HYP-3455's seven-owner gluing clause or named owner/current/state-lift debt.

## Lead codex-2026-06-29: AP84 phase-color grid bridge

**Status:** EVIDENCE / exact AP-tail color-grid discrepancy bridge; not proof
(HYP-3470/T1430/LTI-430/LTT-330), reconnecting HYP-2593/HYP-2594 phase-color
CRT placement to the HYP-3454/HYP-3456/HYP-3457 AP-tail sidecars as the exact
canonical AP84 placement sidecar under HYP-3459's color-packet legality audit,
complementary to HYP-3460's phase-branch pullback, downstream of HYP-3461's
colored-extension gate carrier, and the placement sibling of HYP-3458's
coloring-recursion state.

Claimed files:
`05-knowledge/hypotheses/HYP-3470-lrc14-ap84-color-grid-bridge.md`,
script
`04-computation/lrc14_ap84_color_grid_bridge_codex_20260629.py`,
result
`05-knowledge/results/lrc14_ap84_color_grid_bridge_codex_20260629.out`,
and reflection
`07-reflections/lrc14-ap84-color-grid-bridge-codex-20260629.md`.

Readout: for `{1,2,...,11,13,84m}`, the HYP-2593 form is
`P={1,...,11,13}`, `E={0}`, `V=84m`.  Color `0` is dead and colors `1..13`
share fixed live intervals `[15/182,13/154]`, `[29/70,41/98]`,
`[57/98,41/70]`, and `[141/154,167/182]`; `sigma=426/2695` and `K=52`.
Closed-grid counts validate against direct CRT counts on sample rows through
`m=70`, with `direct_validation_failures=[]`.

The total closed color-grid count has affine period `385`:
`A(m+385)-A(m)=5112`, with no failures for `m=1..770`.  Individual live colors
need period `5005`; the closed-boundary bonus is `0` exactly when `7|m` and
`2` otherwise.  Thus HYP-3456 remains the component/corridor period-`35`
escape clock, while HYP-3470 is the colored `q=14V` CRT placement sidecar
parallel to HYP-3458's coloring-recursion state and subordinate to HYP-3459's
larger color-packet quotient guardrail.  HYP-3460 is the noncanonical
phase-branch pullback sibling; HYP-3461 is the colored-extension gate carrier;
HYP-3470 is the canonical AP84 grid-placement specialization.

Next: prove the four live intervals and the floor/ceiling count symbolically,
then carry the period-`385`/`5005` sidecar only in AP-tail subproofs that need
actual colored CRT placement rather than branch-union component escapes.

## Lead codex-2026-06-29: colored gate-reservoir finite lemma

**Status:** EVIDENCE / exact colored gate quotient audit; not proof
(HYP-3471/T1431/LTI-431/LTT-331), sharpening HYP-3453's gate-escape
transversal with the prior color-reservoir and colored-discrepancy routes, and
using HYP-3462 as the AP84 corridor-splice carrier, HYP-3470 as the exact
AP84 color-grid placement sidecar, and HYP-3461, HYP-3460, HYP-3459, and
HYP-3458 as colored-extension/color-pullback/AP84 subpackets.

Claimed files:
`05-knowledge/hypotheses/HYP-3471-lrc14-colored-gate-reservoir.md`,
script
`04-computation/lrc14_colored_gate_reservoir_codex_20260629.py`,
result
`05-knowledge/results/lrc14_colored_gate_reservoir_codex_20260629.out`,
and reflection
`07-reflections/lrc14-colored-gate-reservoir-codex-20260629.md`.

Readout: on the `135`-row HYP-3453 bank, `130/130` rows with dead components
have a rank-`<=2` E/branch survivor gate.  Low-rank gates split as `7002`
E/branch, `1482` same-branch, and `182` cross-branch.  Endpoint kind has only
`8` colors, numeric mod-`14` has `147`, typed mod-`14` has `360`, structural
sidecars have `161`, and full colored gate words have `1727`.  The AP84
four-color packet covers only `67/130` dead rows, so it is a terminal packet,
not a universal color law.

Next steps: prove the finite implication

```text
dead_components(row) > 0 => rank <= 2 E/branch survivor gate.
```

Try a HYP-3451 component-cover current/cut proof: a dead island with no
E/branch leak should create impossible boundary divergence or a Menger cut in
the branch-coloured blocker graph.  Retain typed endpoint residue plus branch
mask, bad-edge adjacency, and cover deltas until AP packets route through
HYP-3462/HYP-3470/HYP-3461/HYP-3460/HYP-3459/HYP-3458/HYP-3454/HYP-3456/HYP-3457,
same-branch/random gluing routes through HYP-3455, or genuine
owner-current/two-adic/signed-SPEC debt is emitted.

## Lead codex-2026-06-29: AP84 finite transient closure

**Status:** EVIDENCE / finite mixed AP-tail transient certificate; not proof
(HYP-3457/T1417/LTI-417/LTT-317), closing the `m=1..4` side of the
HYP-3452/HYP-3454/HYP-3456 AP-tail split.

Claimed files:
`05-knowledge/hypotheses/HYP-3457-lrc14-ap84-finite-transients.md`,
script
`04-computation/lrc14_ap84_finite_transients_codex_20260629.py`,
result
`05-knowledge/results/lrc14_ap84_finite_transients_codex_20260629.out`,
and reflection
`07-reflections/lrc14-ap84-finite-transients-codex-20260629.md`.

Readout: for `{1,2,...,11,13,84m}` and `m=1..4`, the exact survivor windows
are
`[8/49,(98m-1)/(588m)]`,
`[(98m+1)/(588m),6/35]`,
`[29/35,(490m-1)/(588m)]`, and
`[(490m+1)/(588m),41/49]`, with endpoint labels
`B1:7/E:84m`, `E:84m/B1:5`, `B0:5/E:84m`, and `E:84m/B0:7`.
Against the HYP-3452 component audit:
`exact_window_failures=[]`, `closure_failures=[]`,
`rank_drop_failures_for_m_ge_3=[]`, and `inequality_failures=[]`.  Each
transient row has exactly `4` escapes and `low_rank_escape=4`.  The phase break
is the sign of `(98m+13)/(588m)-6/35=(455-98m)/(2940m)`: positive for
`m=1..4`, reversed at `m=5`.

HYP-3462 now closes the AP-tail carrier/splice by importing HYP-3431 as the
complete branch-union carrier and routing HYP-3454/HYP-3456/HYP-3457 through
HYP-3439.  The active follow-up is non-AP transfer through HYP-3453/HYP-3451
and HYP-3455, not another AP84 finite/transient/carrier audit.

## Lead codex-2026-06-29: AP84 mod-35 floor-count closure

**Status:** EVIDENCE / floor-count derivation of AP-tail escape clock; not
proof (HYP-3456/T1416/LTI-416/LTT-316), closing the sampled clock left by
HYP-3454 and validating against HYP-3452.

Claimed files:
`05-knowledge/hypotheses/HYP-3456-lrc14-ap84-mod35-floor-count.md`,
script
`04-computation/lrc14_ap84_mod35_floor_count_codex_20260629.py`,
result
`05-knowledge/results/lrc14_ap84_mod35_floor_count_codex_20260629.out`,
and reflection
`07-reflections/lrc14-ap84-mod35-floor-count-codex-20260629.md`.

Readout: HYP-3431 gives fixed low corridors `C1=[8/49,6/35]` and
`C0=[29/35,41/49]`.  The moving high gaps
`G_k(m)=[(14k+1)/(588m),(14k+13)/(588m)]` meet `C1` in
`N(m)=floor((504m-6)/70)-floor((96m-13)/14)` positions, and `C0` mirrors the
same count, so `escapes(m)=2*N(m)`.  Against the HYP-3452 component audit
through `m=70`, `mirror_failures=[]`, `formula_failures=[]`, and
`component_audit_failures=[]`.  The period correction vector matches HYP-3454,
with closed shift `N(m+35)-N(m)=12`, hence
`escapes(m+35)-escapes(m)=24`.

HYP-3462 now finishes the AP-tail carrier/splice.  Return to the
HYP-3453/HYP-3451 gate-conductance route and the HYP-3455 noncanonical
seven-owner gluing clause.

## Lead codex-2026-06-29: HYP-3453 gate-escape transversal

Prove the gate-transversal implication exposed by
`04-computation/lrc14_gate_escape_transversal_router_codex_20260629.py`,
using incoming HYP-3452's AP-with-`84m` phase audit as the canonical-tail
base companion and HYP-3439's rescue-core bridge as the one-branch overlap
interface.

Current evidence: on the `135`-row bank,
`rows_with_dead_components_and_low_rank_gate=130/130` and
`rows_dead_without_low_rank_gate=[]`.  The only rows without a low-rank gate
are clean-only rows with `dead=0` and `gates=0`.

Concrete proof route:

1. Show a nonempty dead-cover obstruction forces a boundary mixed component.
2. Apply the HYP-3438 exact gate-word decomposition to that boundary.
3. Use endpoint-spine/corridor-fence parity to force endpoint rank `<=2`, or
   discharge the row as a clean obstruction-empty exit.

Artifacts: `05-knowledge/hypotheses/HYP-3453-lrc14-gate-escape-transversal-router.md`,
`05-knowledge/results/lrc14_gate_escape_transversal_router_codex_20260629.out`,
and `07-reflections/lrc14-gate-escape-transversal-router-codex-20260629.md`.

Next step: formalize the boundary-mixed-component step and test it first on
the AP-with-`84m` family, where the gates are exactly the corridor-fence
endpoints `B1:7|E`, `E|B1:5`, `B0:5|E`, and `E|B0:7`.
## Lead codex-2026-06-29: AP84 endpoint-clock certificate

**Status:** EVIDENCE / endpoint-inequality and residue-clock certificate; not
proof (HYP-3454/T1414/LTI-414/LTT-314), sharpening HYP-3452's AP-tail phase
and HYP-3439's rank-`5` AP-tail descent.  This is the AP-tail endpoint-clock
companion to incoming HYP-3453's bank-level gate-escape transversal router.

Claimed files:
`05-knowledge/hypotheses/HYP-3454-lrc14-ap84-endpoint-clock-certificate.md`,
script
`04-computation/lrc14_ap84_endpoint_clock_certificate_codex_20260629.py`,
result
`05-knowledge/results/lrc14_ap84_endpoint_clock_certificate_codex_20260629.out`,
and reflection
`07-reflections/lrc14-ap84-endpoint-clock-certificate-codex-20260629.md`.

Readout: for `{1,2,...,11,13,84m}`, the checked tail `m=5..70` has no
endpoint failures.  The interval
`[(14ceil(48m/7)+1)/(588m),(14ceil(48m/7)+13)/(588m)]` lies inside the fixed
low corridor `[8/49,6/35]`, has length `1/(49m)`, and is the rank-one
`L[E:84m] R[E:84m]` component-cover escape.  Minimum low-corridor margins are
`1/41160` and `1/2940`.  The transients `m=1..4` are finite mixed
`E:84m/B1:5` cases.

Next: use HYP-3457 for the finite `m=1..4` mixed cases, then splice the
HYP-3454 endpoint interval, HYP-3456 count, and HYP-3457 transient packet into
the HYP-3439 bridge.  The remaining AP-tail work is the fixed-corridor carrier
proof/import and the splice.

## Lead codex-2026-06-29: HYP-3441 incident tax-core residue router

**Status:** EVIDENCE / quotient bridge audit; not proof
(HYP-3441/T1402/LTI-402/LTT-302), continuing HYP-3437's one-branch
overlap-tax cut certificate and now pointing into HYP-3438/HYP-3453/HYP-3451.

Claimed files:
`05-knowledge/hypotheses/HYP-3441-lrc14-incident-tax-core-residue-router.md`,
script
`04-computation/lrc14_incident_tax_core_residue_router_codex_20260629.py`,
result
`05-knowledge/results/lrc14_incident_tax_core_residue_router_codex_20260629.out`,
and
`07-reflections/lrc14-incident-tax-core-residue-router-codex-20260629.md`.

Readout: HYP-3437's `59` negative naive-slack rows have strict rescue ranks
`{2:7,4:2,5:48,6:2}` but incident tax-core sizes `{1:10,2:49}`.  Strict rank
is larger than incident size in `59/59`; unit-only cores are `57/59`;
ramified apex-7 cores are `2/59`; and the dominant incident core is `(11,13)`
in `48/59` rows.  Endpoint-bounded tax covers the deficit in only `6/59`, so
this is a bridge obligation, not an endpoint-only theorem.

Next: join the `53` endpoint-tax failures to the HYP-3453 gate ledger.  For
each failure, decide whether the touched incident atom can be legally recovered
as a rank-`<=2` survivor gate when `dead_components(row)>0`, whether the
canonical AP row routes through HYP-3431/HYP-3452/HYP-3454, whether the
noncanonical rank-`6` obstruction is the HYP-3455 finite gate-gluing clause, or
whether the row emits owner-current, endpoint-spine, exact-period/state-lift,
signed-SPEC, or two-adic debt.  Keep the C3/Qsqrt(-7) residue packet as a
router until strict co-owner multiplicity is reconstructed.

## Lead codex-2026-06-29: HYP-3438 survivor-gate word audit

**Status:** EVIDENCE / exact survivor-gap gate-word classification; not proof
**Status:** EVIDENCE / exact survivor-gate classification; not proof
(HYP-3438/T1399/LTI-399/LTT-299), immediate follow-up to the HYP-3436
minimal bad-core cover extractor.

Claimed files:
`05-knowledge/hypotheses/HYP-3438-lrc14-survivor-gate-word-audit.md`,
script
`04-computation/lrc14_survivor_gate_word_audit_codex_20260629.py`,
result
`05-knowledge/results/lrc14_survivor_gate_word_audit_codex_20260629.out`,
and
`07-reflections/lrc14-survivor-gate-word-audit-codex-20260629.md`.

Readout: on the `135`-row HYP-3436 bank, `6228` mixed components emit `8702`
survivor gates.  Branch masks are `both:1064`, `branch0:3819`,
`branch1:3819`; adjacency is `left_bad_edge:3515`, `right_bad_edge:3515`,
`two_sided:1672`; all parent endpoints are `E|E`.  Tight AP-with-84 has four
one-gate edge words matching the corridor-fence shadow.

Next: prove the survivor-gate no-gluing lemma.  Canonical `E:84` gates should
route through HYP-3431 corridor-fence; noncanonical multi-owner and two-sided
gates should route through HYP-3437/HYP-3439 endpoint-Menger rescue cores.
HYP-3450/HYP-3451 supply the component-cover conductance sibling: use their
branch-alive/dead projection to choose obstruction cuts, then classify the
mixed components inside those cuts by the HYP-3438 gate words before any scalar
compression is allowed.

## Lead codex-2026-06-29: Rescue-core bridge certificate

**Status:** EVIDENCE / AP-tail bridge audit; not proof
(HYP-3439/T1400/LTI-400/LTT-300), continuing HYP-3437's one-branch
overlap-tax cut certificate, HYP-3438's survivor gates, and HYP-3450/HYP-3451's
component-cover escape router.
Readout: on the `135` primitive covering rows, the audit finds `8702` survivor
gates over `6228` mixed components.  Route histogram:
`edge_singleton_parent_gate:5950`, `edge_survivor_residual:1080`,
`owner_current_small_delta:794`, and `mixed_owner_residual:878`; `79` rows have
mixed-owner residual gates.

Canonical clue: for `{1,...,11,13,84m}`, a one-period `m=1..35` probe supports
the mod-35 law.  Outer `(13,7)` edge gates occur iff `7` does not divide `m`;
inner `(11,5)` edge gates occur iff `5` does not divide `m`; when `35 | m`,
there are no mixed survivor gates and survival comes from clean `E_safe`
components.

Next: prove the canonical mod-35 gate law by endpoint arithmetic, then build a
residual-gate discharging table routing mixed-owner and edge-survivor residuals
through HYP-3437 overlap-tax cuts, HYP-3439 rescue cores, HYP-3440 endpoint-cut
sidecars, owner-current, two-adic loss, exact-period/state-lift, or
signed-SPEC.
## Lead codex-2026-06-29: Rescue-core bridge certificate

**Status:** RESERVED STUB / computation pending; not proof
(HYP-3439/T1400/LTI-400/LTT-300), continuing HYP-3437's one-branch
overlap-tax cut certificate and HYP-3436's two-color bad-core cover extractor.

Claimed files:
`05-knowledge/hypotheses/HYP-3439-lrc14-rescue-core-bridge-certificate.md`,
script
`04-computation/lrc14_rescue_core_bridge_certificate_codex_20260629.py`,
result
`05-knowledge/results/lrc14_rescue_core_bridge_certificate_codex_20260629.out`,
and
`07-reflections/lrc14-rescue-core-bridge-certificate-codex-20260629.md`.

Readout: HYP-3437's `150`-row recap has `59` negative naive-slack rows, all
rescued, with rescue-rank histogram `{0:91,2:7,4:2,5:48,6:2}` and rank `6`
only at `covering_AP_with_84` / `canonical_84m_01`.  The AP/84m bridge audits
`14` displayed rows (`13` unique speeds): rank histogram `{0:1,5:11,6:2}`,
rank `6` only for `covering_AP_with_84` / `ap_omit_12_tail_84x01`, and tails
`84x02..84x12` drop to the rank-`5` core `(5,7,9,11,13)` while keeping
HYP-3450/HYP-3451 low-rank two-colour escapes.  HYP-3452 now sharpens the
tail handoff, with HYP-3457 closing finite transients `m=1..4`, HYP-3454
closing the rank-one `E:84m/E:84m` endpoint phase for `m>=5`, and HYP-3456
closing the mod-`35` escape clock.

Next: prove the canonical rank-`6` base case using HYP-3431 plus the
HYP-3450/HYP-3451 four-escape component graph, use HYP-3452 for the rank-`5`
AP-tail descent, then extend the bridge to arbitrary primitive covering rows
or name the first residual high-rank sidecar.
## Lead codex-2026-06-29: Minimal bad-core cover extractor

**Status:** EVIDENCE / exact local bad-core cover classification; not proof
(HYP-3436/T1397/LTI-397/LTT-297), inverting HYP-3435's branch-cover
certificate after HYP-3425's two-color bad-core identity and HYP-3434's
overlap-tax warning.

Claimed files:
`05-knowledge/hypotheses/HYP-3436-lrc14-minimal-bad-core-cover-extractor.md`,
script
`04-computation/lrc14_minimal_bad_core_cover_extractor_codex_20260629.py`,
result
`05-knowledge/results/lrc14_minimal_bad_core_cover_extractor_codex_20260629.out`,
reflection
`07-reflections/lrc14-minimal-bad-core-cover-extractor-codex-20260629.md`,
and forum post
`poke-forum/posts/20260629-lrc14-minimal-bad-core-cover-extractor/post.md`.

Purpose: classify the obstruction

```text
E_safe cap B0_odd cap B1_odd
```

by exact bad-core intervals and minimal branch-0/branch-1 odd-owner covers.

Readout: on `135` rows, branch-union identity and positive survivors hold
`135/135`; there are `17164` even-safe components, `11670` bad-core
components, `15868` survivor components, max minimal two-color owner count
`6`, and `(1,1)` minimal cover signature on `10288/11670` bad-core
components.  Endpoint support is `(1,1)` on `11634/11670`.

Canonical clue: for `{1..11,13,84m}`, no failures through `m=30`, and every
bad-core component is singleton/singleton from `m=3`.

Next: prove the canonical singleton-tail law for all `m>=3`, then build the
general component-chain owner-cover graph and prove an all-covered chain emits
owner-current, endpoint-spine, overlap-tax, two-adic, exact-period, or
state-lift debt.
## Lead codex-2026-06-28: Component-cover obstruction and conductance router

**Status:** EVIDENCE / exact component obstruction audit plus graph router; not
proof (HYP-3450/T1410/LTI-410/LTT-310 and HYP-3451/T1411/LTI-411/LTT-311),
continuing HYP-3435's branch-cover certificate, HYP-3434's overlap-tax
identity, and HYP-3429's endpoint-spine route.

Claimed files:
`05-knowledge/hypotheses/HYP-3450-lrc14-component-cover-obstruction-extractor.md`,
`05-knowledge/hypotheses/HYP-3451-lrc14-component-cover-conductance-router.md`,
scripts
`04-computation/lrc14_component_cover_obstruction_extractor_codex_20260628.py`
and
`04-computation/lrc14_component_cover_conductance_router_codex_20260628.py`,
results
`05-knowledge/results/lrc14_component_cover_obstruction_extractor_codex_20260628.out`
and
`05-knowledge/results/lrc14_component_cover_conductance_router_codex_20260628.out`,
and reflections in `07-reflections/`.

Purpose: make the finite lemma target local.  For each `E_safe` component,
record branch0/branch1 survival, minimal odd-bad covers when a branch dies,
and endpoint rank/labels when a branch survives.  Then project dead components
to a branch-coloured blocker graph for Menger/Green-current/conductance proof
routes.

Readout: HYP-3450 audits `135` rows and `17164` components, with
`rows_with_endpoint_rank_le_2_survivor=135/135`; component classes are
`both_alive=6492`, `branch0_only=3451`, `branch1_only=3451`,
`dead_both=3770`, and all dead components have paired minimal covers.  HYP-3451
finds `rows_with_low_rank_escape=135/135`, max paired-cover rank `6` at
`random_covering_082`, max dead fraction `0.962963` at
`ap_omit_12_tail_84x05`, and AP-with-`84m` tails as the danger rows with
connected dead-cover projections and only four low-rank escapes.

Next: prove the AP-with-`84` graph base case.  Its `22` dead components form
one blocker-projection component, yet four endpoint-rank `2` escapes remain.
Try a bounded Menger cut, Green-current boundary imbalance, or
algebraic-connectivity obstruction; then lift to AP-with-`84m` tails and
arbitrary primitive rows.

## Lead codex-2026-06-28: Euler-Mascheroni harmonic wall-budget sidecar

**Status:** SYNTHESIS / exact rational harmonic sidecar audit; not proof
(HYP-3432/T1393/LTI-393/LTT-293), extending HYP-3430's scalar firewall,
HYP-3429's endpoint-spine certificate, and HYP-3427's wall-signature atlas.

Claimed files:
`05-knowledge/hypotheses/HYP-3432-lrc14-euler-mascheroni-wall-budget.md`,
script
`04-computation/lrc14_euler_mascheroni_wall_budget_codex_20260628.py`,
result
`05-knowledge/results/lrc14_euler_mascheroni_wall_budget_codex_20260628.out`,
and
`07-reflections/lrc14-euler-mascheroni-wall-budget-codex-20260628.md`.

Purpose: convert the Euler-Mascheroni prompt from HYP-3430's scalar-intercept
firewall into an exact finite sidecar.  Instead of using `gamma` or logs as a
certificate, attach to each endpoint wall set the rational reciprocal budget
`sum 1/v`.

Readout: on `150` HYP-3429 best endpoint spines, there are `144` distinct
budgets and only `2` scalar shape-collisions; budget range is
`1/4032..89/420`.  On all `5524` HYP-3427 wall windows, however, only `1291`
distinct budgets remain, and `1197` budgets have more than one wall signature.

Next: use harmonic budget as a priority queue inside the endpoint-spine proof,
but accept a certificate only with its exact interval, branch, and wall labels.
## Lead codex-2026-06-28: Gamma harmonic-sieve overlap remainder

**Status:** EVIDENCE / exact harmonic-sieve remainder audit; not proof
(HYP-3434/T1395/LTI-395/LTT-295), extending HYP-3433's endpoint-spine
finite-part ledger, HYP-3432's harmonic wall-budget sidecar, HYP-3431's
corridor fence, HYP-3430's scalar firewall, HYP-3429's endpoint-spine
certificate and HYP-3426's one-branch mirror reduction.

Claimed files:
`05-knowledge/hypotheses/HYP-3434-lrc14-gamma-harmonic-sieve-remainder.md`,
script
`04-computation/lrc14_overlap_tax_harmonic_sieve_remainder_codex_20260628.py`,
result
`05-knowledge/results/lrc14_overlap_tax_harmonic_sieve_remainder_codex_20260628.out`,
and
`07-reflections/lrc14-overlap-tax-harmonic-sieve-remainder-codex-20260628.md`.

Purpose: make the harmonic/compression angle exact.  In the one-branch target,

```text
branch0_mass = naive_slack + overlap_tax
```

with `naive_slack = |E_safe| - sum_o |E_safe cap B0_o|` and
`overlap_tax = sum_o |E_safe cap B0_o| - |E_safe cap union_o B0_o|`.

Readout on the `150`-row endpoint-spine bank: exact identity `150/150`,
positive branch0 survivor `150/150`, naive slack positive `91/150`, naive
slack negative `59/150`, positive overlap tax `150/150`, and all negative
rows rescued `59/59`.  Tightest rescue is `{1..11,13,84}` with tax/deficit
`1.090875`.

Next: prove the dichotomy: either naive slack is nonnegative, or a rank-2
endpoint-spine / conductance-cut certificate forces overlap tax above the
deficit.  Treat Euler-Mascheroni only as harmonic denominator-prefix
calibration.

## Lead codex-2026-06-28: One-branch mirror / endpoint-support certificates

**Status:** EVIDENCE / exact finite-ruler mirror and endpoint-support audit;
not proof (HYP-3426/T1387/LTI-387/LTT-287).

Claimed files:
`05-knowledge/hypotheses/HYP-3426-lrc14-one-branch-mirror-endpoint-support.md`,
script
`04-computation/lrc14_one_branch_mirror_endpoint_support_codex_20260628.py`,
result
`05-knowledge/results/lrc14_one_branch_mirror_endpoint_support_codex_20260628.out`,
and
`07-reflections/lrc14-one-branch-mirror-endpoint-support-codex-20260628.md`.

Purpose: sharpen HYP-3425.  The involution `u -> 1-u` preserves `E_safe` and
maps branch-0 survivor intervals to branch-1 survivor intervals, so the
finite theorem target can be reduced to the one-color cover lemma:

```text
E_safe is not contained in B0_odd.
```

Readout on `162` rows: mirror identity, equal branch measures, positive
one-branch survivor, selected score `>=1/14`, and endpoint-labelled survivors
all hold `162/162`.  Endpoint-owner support histogram is
`{1:353, 2:13103, 3:72}`, with max support size `3`.

Next: classify endpoint-owner triples and prove that no one-color odd
near-integer interval cover can consume every even-safe component.
## Lead codex-2026-06-28: Two-branch wall-signature atlas

**Status:** EVIDENCE / exact wall-certificate scout; not proof
(HYP-3427/T1388/LTI-388/LTT-288), refining HYP-3425.

Claimed files:
`05-knowledge/hypotheses/HYP-3427-lrc14-two-branch-wall-signature-atlas.md`,
script
`04-computation/lrc14_two_branch_wall_signature_atlas_codex_20260628.py`,
result
`05-knowledge/results/lrc14_two_branch_wall_signature_atlas_codex_20260628.out`,
and
`07-reflections/lrc14-two-branch-wall-signature-atlas-codex-20260628.md`.

Purpose: turn positive survivor windows into exact certificates.  Each window
is labelled by branch mask, endpoint wall labels, and midpoint binders, using
wall types `E:s`, `O0:o`, and `O1:o`.

Readout: `67` audited rows have survivor windows `67/67`; total survivor
windows `5524`; global signature types `27`; branch masks
`b0=2255, b1=2255, both=1014`.  The tight canonical row `{1..11,13,84}` has
four windows bounded by `E:84` and odd walls `5,7`.

Next: prove the bounded wall-alphabet lemma for primitive covering rows, or
classify the first failure by owner-current, sheet, exact-period, state-lift,
or named two-adic debt.
## Lead codex-2026-06-28: Component-spine endpoint certificate

**Status:** EVIDENCE / exact endpoint-spine audit; not proof
(HYP-3429/T1390/LTI-390/LTT-290), extending HYP-3428's two-adic loss ledger,
HYP-3427's wall-signature atlas, and HYP-3425's Helly obstruction route while
complementing HYP-3426's mirror reduction.

Claimed files:
`05-knowledge/hypotheses/HYP-3429-lrc14-component-spine-certificate.md`,
script
`04-computation/lrc14_component_spine_certificate_codex_20260628.py`,
result
`05-knowledge/results/lrc14_component_spine_certificate_codex_20260628.out`,
and
`07-reflections/lrc14-component-spine-certificate-codex-20260628.md`.

Purpose: compress the HYP-3425 interval-union target into low-rank endpoint
certificates.  A survivor window is labelled by active endpoint walls:

```text
E  = even-safe wall
B0 = branch-0 odd near-integer wall
B1 = branch-1 odd near-half wall
```

Readout on `150` rows: positive survivor windows `150/150`, all survivor
windows endpoint-labelled `15576/15576`, best endpoint-spine rank `<=2` on
`150/150`, mixed even/odd spines on `148/150`, and both-branch windows on
`149/150`.  Best-rank histogram `{1:47,2:103}`.

Next: prove the endpoint-spine lemma: every primitive covering row has either
an E-only free component or a rank-2 mixed endpoint spine, with owner-current
labels used only for genuine exceptions.

## Lead codex-2026-06-28: Euler-Mascheroni harmonic intercept firewall

**Status:** EVIDENCE / scalar guardrail audit; not proof
(HYP-3430/T1391/LTI-391/LTT-291), extending HYP-3429 and the
Mertens/loglog guardrail lanes.

Claimed files:
`05-knowledge/hypotheses/HYP-3430-lrc14-euler-mascheroni-harmonic-intercept.md`,
script
`04-computation/lrc14_euler_mascheroni_harmonic_intercept_codex_20260628.py`,
result
`05-knowledge/results/lrc14_euler_mascheroni_harmonic_intercept_codex_20260628.out`,
and
`07-reflections/lrc14-euler-mascheroni-harmonic-intercept-codex-20260628.md`.

Purpose: test whether the finite intercept `H_N - log N` can replace the
endpoint-spine certificate.  It cannot.  On the HYP-3429 bank, there are `11`
endpoint certificate classes, same-max-speed bins with mixed classes `19/108`,
and rounded-4 gamma-intercept bins with mixed classes `21/30`.

Next: require every Euler-Mascheroni, Mertens, or loglog tail estimate in the
covering-floor route to declare the retained sidecar: endpoint owner, wall
signature, two-adic loss class, sheet, exact period, or state-lift debt.

## Lead codex-2026-06-28: Two-branch obstruction / Helly certificate

**Status:** EVIDENCE / exact finite-ruler interval audit; not proof
(HYP-3425/T1386/LTI-386/LTT-286), downstream of HYP-3424's
covering-floor duality transfer.

Claimed files:
`05-knowledge/hypotheses/HYP-3425-lrc14-two-branch-obstruction-helly.md`,
script
`04-computation/lrc14_two_branch_obstruction_helly_codex_20260628.py`,
result
`05-knowledge/results/lrc14_two_branch_obstruction_helly_codex_20260628.out`,
and
`07-reflections/lrc14-two-branch-obstruction-helly-codex-20260628.md`.

Purpose: sharpen HYP-3422's relocation target.  For `S=O union 2E`,
`u=2t`, write the failure as a two-color odd bad core:

```text
relocation_good = E_safe minus (B0_odd cap B1_odd).
```

Readout on `62` rows: positive two-branch good union `62/62`; selected
relocation score `>=1/14` `62/62`; tight canonical row `{1..11,13,84}` has
good union `1/105` over four surviving components.

Next: prove the Helly/interval-piercing lemma that `E_safe` is not contained
in `B0_odd cap B1_odd` for every primitive covering row, then feed the result
into the HYP-3422 two-adic descent.

## Lead codex-2026-06-28: Two-adic off-grid relocation lemma

**Status:** EVIDENCE / exact interval scout; not proof
(HYP-3422/T1383/LTI-383/LTT-283).

Claimed files:
`05-knowledge/hypotheses/HYP-3422-lrc14-two-adic-offgrid-relocation.md`,
script
`04-computation/lrc14_two_adic_offgrid_relocation_codex_20260628.py`,
result
`05-knowledge/results/lrc14_two_adic_offgrid_relocation_codex_20260628.out`,
and
`07-reflections/lrc14-two-adic-offgrid-relocation-codex-20260628.md`.

Purpose: make the corrected HYP-3418 covering-floor route concrete.  Split
`S = O union 2E` and set `u=2t`; then even safety is the halved problem on `E`,
while odd speeds impose two explicit lift filters.  The finite lemma target is

```text
E_safe(1/14) cap (odd_branch_0_good union odd_branch_1_good) != empty.
```

Readout on `24` audited covering rows: full exact `M(S)>=1/14`, off-`14`-grid
optimizers, branch-0 certificates, branch-1 certificates, and either-branch
relocation certificates all hold `24/24`; the naive nonresonant witness fails
`24/24`.  Tight branch-union measure is `563/105105`; tight relocation score is
`11/147`.

Next: prove the interval-overlap lemma with finite-ruler, Helly, or exact
interval-piercing bounds on `E_safe`.  Use the owner-cut work only to name
finite packets where branch filters appear tight; the `2:g2` label in
HYP-3419 is the finite owner-cut shadow of this two-adic coordinate.
## Lead codex-2026-06-28: Canonical corridor-fence certificate

**Status:** PARTIAL PROOF / all-`m` canonical relocation certificate; not full
LRC14 proof (HYP-3431/T1392/LTI-392/LTT-292).

Claimed files:
`05-knowledge/hypotheses/HYP-3431-lrc14-canonical-corridor-fence.md`,
script
`04-computation/lrc14_canonical_corridor_fence_codex_20260628.py`,
result
`05-knowledge/results/lrc14_canonical_corridor_fence_codex_20260628.out`,
and
`07-reflections/lrc14-canonical-corridor-fence-codex-20260628.md`.

Purpose: convert the tight canonical `84m` evidence downstream of HYP-3426's
one-branch mirror/endpoint-owner audit, HYP-3427's wall-signature atlas,
HYP-3428's two-adic loss ledger, HYP-3429's endpoint-spine certificate,
HYP-3430's harmonic-intercept firewall, and HYP-3425's Helly target into an
all-`m` wall proof.  The low odd/even core
leaves fixed corridors
`[8/49,6/35]` and `[29/35,41/49]`, each of length `2/245`; the moving high
even half-speed `42m` removes only disjoint grid intervals of width
`1/(294m)`.  Since a connected corridor is longer than every moving bad
component for every `m>=1`, the canonical tower has positive two-branch
relocation for all `m`.

Next: search non-canonical covering rows for the same structure: a low-core
branch corridor longer than all moving bad components.  If found, route by the
corridor-fence lemma; otherwise fall back to HYP-3430 scalar-firewall sidecars,
HYP-3429 endpoint-spine targets, HYP-3428 loss classes, HYP-3425 component
Helly, or owner-current exception labels.

## Lead codex-2026-06-28: Q-uniform topology / q-specific arithmetic break guardrail

**Status:** SYNTHESIS / executable proof-route guardrail; not proof
(HYP-3423/T1384/LTI-384/LTT-284).

Claimed files:
`05-knowledge/hypotheses/HYP-3423-lrc14-quniform-topology-arithmetic-break-guardrail.md`,
script
`04-computation/lrc14_quniform_topology_arithmetic_break_codex_20260628.py`,
result
`05-knowledge/results/lrc14_quniform_topology_arithmetic_break_codex_20260628.out`,
and
`07-reflections/lrc14-quniform-topology-arithmetic-break-codex-20260628.md`.

Purpose: prevent a uniform topological/index argument from being used as a
q-specific magnitude proof.  The route legalizes C2/Borsuk-Ulam and C6/Galois
arguments for the residue/equioscillation half, but forces magnitude claims to
carry arithmetic, owner-current, or floor sidecars.

Readout: over `q=3..22`, the C2/BU residue charge is present on `20/20` rows,
while the HYP-3413 canonical GW magnitude switch is ON only on `7/20` rows,
exactly `q == 1 mod 3`.  This records the requested contrast: `q=4,7` ON and
`q=5,6` off.

Synthesis: HYP-3417's frontier cut `{2:g2,11:g1,13:g1}` remains valuable
because it is one even-cover label plus two binding labels, but HYP-3423 marks
it as a local owner-current certificate, not a replacement for the q-mod-3
switch, S259/HYP-3418 two-adic descent, or HYP-3415 decorrelation floor.

Next: audit the active proof stack for topology-to-magnitude inferences.  Every
route invoking C2/BU, C6, Galois symmetry, or topological degree should declare
whether it certifies residue/equioscillation only, or which sidecar resurrects
the q-specific magnitude coordinate.

## Lead codex-2026-06-28: Covering-floor duality transfer

**Status:** SYNTHESIS / proof-routing scout; not proof
(HYP-3424/T1385/LTI-385/LTT-285), with HYP-3423 topology/arithmetic guardrail,
HYP-3422 two-adic relocation, and HYP-3421 off-grid transparency/Rprime
closure as floor-facing companions.

Claimed files:
`05-knowledge/hypotheses/HYP-3424-lrc14-covering-floor-duality-transfer.md`,
script
`04-computation/lrc14_covering_floor_duality_transfer_codex_20260628.py`,
result
`05-knowledge/results/lrc14_covering_floor_duality_transfer_codex_20260628.out`,
and
`07-reflections/lrc14-covering-floor-duality-transfer-codex-20260628.md`.

Purpose: integrate prior even/odd, positive/negative, odd/even, and
addition/multiplication duality work after the S259/HYP-3418 correction that
the covering floor is 2-adic.  The odd/coprime witness at `t=1/2` is not a
floor proof because every even speed dies there; odd data becomes phase-cover
debt, and the even fold becomes the floor-bearing object.

Readout: top transfer gates are two-adic even floor descent, signed SPEC
decorrelation floor, even-good/odd-phase cover, recursive quotient sidecar
router, and owner-current even-cover sidecar.  The tracked proof obligations
have minimum covers of size `3`, led by two-adic descent plus signed SPEC floor
plus either quotient router or owner sidecar.  Raw 7-adic/C3/Galois/census
shadows are demoted unless they feed the floor.

Next: turn the candidate router into a concrete lemma on covering packets:
after q-witness and LRC<=13 induction, every packet should feed `|SPEC| <
product`, HYP-3422 two-adic interval relocation, odd phase-cover debt, finite
owner-current/Menger sidecar, add/mult energy-to-SPEC penalty, HYP-3421 off-grid/signed-SPEC
transparency, or an explicit off-path filter.

## Lead monad-explorer-2026-06-28: Additive energy needs a sheet sidecar

**Status:** SYNTHESIS / exact counterexample-and-repair scout; not proof
(HYP-3425/T1386).

Claimed files:
`05-knowledge/hypotheses/HYP-3425-lrc14-additive-energy-sheet-sidecar.md`,
script
`04-computation/lrc14_additive_energy_sheet_sidecar_codex_20260628.py`,
result
`05-knowledge/results/lrc14_additive_energy_sheet_sidecar_codex_20260628.out`,
and
`07-reflections/lrc14-additive-energy-sheet-sidecar-codex-20260628.md`.

Purpose: test the HYP-3424 add/mult transfer rule against the actual covering
floor packet.  The scout combines HYP-3140 canonical `Rprime` rows, HYP-3422
covering packets, and a small exact random bank.  Exact collision bank:
`canonical_r3` and `canonical_r5` share `fullE=389` but have opposite signs of
`delta`; `canonical_r1_drop12` and `covering_AP_with_84` share `RE=246` and
`oddE=47` but have `Rprime=7/6` versus `0.513954...`.

Readout: raw additive-energy scalars collide across covering packets.  Minimal
exact repairs are energy-plus-sheet packets such as `(RE,q_zero_mass)`,
`(oddE,q_zero_mass)`, and `(fullE,q_range_hi)`.  Small exact random bank:
`corr(RE,delta)=+0.628`, `corr(oddE,delta)=+0.134`,
`corr(evenE,delta)=-0.047`.

Next: formulate the packet theorem.  Any future use of HYP-2272 on the
covering floor must retain a sheet-profile sidecar from HYP-3140 before
claiming SPEC control, phase-cover debt, or terminal floor behavior.  Good
candidate sidecars: zero-sheet mass, sheet range, far depth, or another
fiber-PGF-equivalent field.
## Lead codex-2026-06-29: Bad-core cover gluing obstruction

**Status:** EVIDENCE / exact interval-cover classification; not proof
(HYP-3436/T1397/LTI-397/LTT-297).
## Lead codex-2026-06-29: Minimal bad-core cover no-gluing theorem

**Status:** EVIDENCE / exact obstruction-side interval-cover classification;
not proof (HYP-3436/T1397/LTI-397/LTT-297).

Claimed files:
`05-knowledge/hypotheses/HYP-3436-lrc14-minimal-bad-core-cover-extractor.md`,
script
`04-computation/lrc14_minimal_bad_core_cover_extractor_codex_20260629.py`,
result
`05-knowledge/results/lrc14_minimal_bad_core_cover_extractor_codex_20260629.out`,
and reflection
`07-reflections/lrc14-minimal-bad-core-cover-extractor-codex-20260629.md`.

Purpose: invert HYP-3435's positive branch-cover certificate and classify
`E_safe cap B0_odd cap B1_odd` by local bad-core components, minimal B0/B1
odd-owner covers, even endpoint gates, and survivor gaps.  This turns the next
finite lemma into a local-to-global gluing problem rather than a scalar measure
problem.

Readout: `135/135` rows still have survivors.  Across `17164` even-safe
components, `3770` are fully bad, `6228` mixed, and `7166` clean; there are
`11670` bad-core components and `15868` survivor components.  Minimal cover
pairs are led by `(1,1):10288`, `(1,2):(2,1):531`, `(2,2):184`, with maximum
total cover size `6`; endpoint support size stays at most `3`.  Tight
AP-with-84 has `bad_core=314/735`, survivor `1/105`, `22` fully bad even
components, `4` mixed components, and cover pairs `(1,1):14,(1,2):6,(2,1):6`.

Creative extension: prove that a compatible full cover cannot glue these local
ledgers across all even-safe components.  If gluing survives in a future bank,
the first compatible full cover should be routed to HYP-3431 corridor-fence
structure, HYP-3429 endpoint spines, HYP-3428 two-adic loss, HYP-3417/HYP-3420
owner-current, HYP-3423 legality, or HYP-3421/HYP-3129 signed-SPEC debt.
and
`07-reflections/lrc14-minimal-bad-core-cover-extractor-codex-20260629.md`.

Purpose: invert HYP-3435's positive branch-cover certificate.  Instead of
only asking where `E_safe cap (branch0_good union branch1_good)` survives, the
extractor classifies the complementary obstruction
`bad_core = E_safe cap B0_odd cap B1_odd` component-by-component.  It records
minimal branch-0 and branch-1 odd-owner covers, even endpoint gates, low/high
odd endpoint labels, and neighboring survivor gaps.

Readout: `135` primitive covering rows audited, `11670` bad-core cells,
exact decomposition `E_safe = branch_union disjoint-union bad_core` on
`135/135`, positive branch union on `135/135`, positive bad core on `133/135`,
minimum branch union `1/105`, minimum positive survivor gap `1/162540`,
maximum bad-core mass `314/735`, maximum bad-core cell length `1/49`, and
maximum minimal-cover signature `(5,1)`.  Most bad-core cells are `(1,1)`,
but hard rows force `(3,3)` and `(5,1)` local cover signatures.

Next: prove the minimal-cover no-gluing lemma.  A primitive covering row
should not be able to cover every `E_safe` component by compatible low/high
odd-owner covers unless it emits owner-current, exact-period, state-lift,
signed-SPEC/Rprime, or two-adic loss debt.  Suggested split: `(1,1)` and
canonical covers by corridor/Helly fences; multi-owner covers by
owner-current or endpoint-Menger cuts; wide nonlocal covers by signed-SPEC,
exact-period, state-lift, or two-adic descent.

Rebase integration: HYP-3437 now reserves the overlap-tax Menger-cut
certificate.  Feed HYP-3436 bad-core atoms into that graph-cut route whenever
the no-gluing obstruction is multi-owner or matches HYP-3434 negative-slack /
overlap-tax rescue behavior.

## Lead codex-2026-06-28: Two-adic branch-cover certificate

**Status:** EVIDENCE / finite-ruler interval-certificate stress; not proof
(HYP-3435/T1396/LTI-396/LTT-296).

Claimed files:
`05-knowledge/hypotheses/HYP-3435-lrc14-two-adic-branch-cover-certificate.md`,
script
`04-computation/lrc14_two_adic_branch_cover_certificate_codex_20260628.py`,
result
`05-knowledge/results/lrc14_two_adic_branch_cover_certificate_codex_20260628.out`,
and
`07-reflections/lrc14-two-adic-branch-cover-certificate-codex-20260628.md`.

Purpose: make HYP-3422's interval relocation, HYP-3425's two-color bad-core
identity, HYP-3426's mirror reduction, HYP-3427's wall atlas, and HYP-3428's
loss ledger proof-facing by recording finite-ruler certificates: chosen branch
cells, active odd/even endpoint gates, and deletion sensitivities for both odd
blockers and even gates.  The incoming HYP-3430 Euler-Mascheroni firewall is
absorbed as a guardrail: scalar tails calibrate only after endpoint/branch
sidecars survive.  The incoming HYP-3431 canonical corridor-fence proof is the
all-`m` base case this general extractor must reproduce.  HYP-3432 adds the
reciprocal wall-budget sidecar as an endpoint-debt priority queue, not as a
replacement for exact interval/branch/wall labels.  The theorem target is no
longer a scalar branch measure; it is that any attempted cover of
`E_safe(1/14)` by odd low-bad and odd high-bad interval families must emit a
small endpoint-gate ledger.

Readout: `135` primitive covering rows audited (`15` structured, `120` random);
certificate success, branch0 positivity, branch1 positivity, and both-branch
positivity all `135/135`.  Tight AP-with-84 row has even-safe measure `107/245`,
branch measures `563/105105`, branch union `1/105`, selected witness
`t=2293/3920`, and score `59/784`.  Finite-bank minimum selected score is
`1283/17160`, margin `401/120120` over `1/14`.

Creative extension completed by HYP-3436: the witness script has been inverted
into a minimal-cover extractor.  Continue by proving the no-gluing theorem over
the emitted bad-core cover signatures and endpoint ledgers.

## Lead codex-2026-06-28: Owner-cut dual current certificate synthesis

**Status:** SYNTHESIS / exact owner-current certificate scout; not proof
(HYP-3417/T1378/LTI-378/LTT-278).

Claimed files:
`05-knowledge/hypotheses/HYP-3417-lrc14-owner-cut-dual-certificate-synthesis.md`,
script
`04-computation/lrc14_owner_cut_dual_certificate_synthesis_codex_20260628.py`,
result
`05-knowledge/results/lrc14_owner_cut_dual_certificate_synthesis_codex_20260628.out`,
and
`07-reflections/lrc14-owner-cut-dual-certificate-synthesis-codex-20260628.md`.

Purpose: push HYP-3410's Menger-owner-cut route one step closer to a theorem
by turning owner cuts into labelled signed-current certificates with declared
zero level and margin.  This keeps the finite endpoint-owner packet instead of
scalarizing to cut size, p-adic closeness, or a famous-problem name.

Readout: all three current HYP-3410 mixed fibers have margin-1 certificates:

```text
height leak:       positive-debt {5:g1}
13->26 owner leak: unit-island {1:g1}
10->20 frontier:  positive-debt {2:g2,11:g1,13:g1}
```

The frontier cut is exactly one even-cover label plus two binding labels.  The
positive-debt Sophie audit gives `(core,cut)=(1,2)->13/5` and
`(1,3)->25/13`, matching the top frontier finite-BDH variance label
`13:g1=49/50`.

Rebase signal: incoming S257/HYP-3411-HYP-3413 makes the residue `C3` layer the
gate for the global GW doubling criterion.  This lead is narrower: it keeps the
same 7-adic/2-adic split visible inside the local frontier owner-current cut,
but it still has to be extended into a bounded finite certificate theorem.
After S258/HYP-3415, this lead is explicitly off the main proof-completion
path unless it feeds the decorrelation floor inequality.
Incoming HYP-3416 supplies the recursive quotient-ladder setting; HYP-3417 is
one owner-current certificate layer for that stack, not a replacement for it.

Next: extend the HYP-3406/HYP-3410 enlarged bank beyond `(72,20)`.  For every
mixed residue/height fiber, compute the common owner core, unit-island current,
Mertens-cheapest positive-debt current, signed-current margin, Sophie channel
pair, and top finite-BDH variance label.  If no bounded owner current exists,
name the first owner/height/off-grid/exact-period/state-lift debt.
## Lead codex-2026-06-28: Charal owner-cut recursion prototype

**Status:** SYNTHESIS / finite cut-recursion API prototype; not proof
(HYP-3419/T1380/LTI-380/LTT-280).

Claimed files:
`05-knowledge/hypotheses/HYP-3419-lrc14-charal-owner-cut-recursion-prototype.md`,
script
`04-computation/lrc14_charal_owner_cut_recursion_prototype_codex_20260628.py`,
result
`05-knowledge/results/lrc14_charal_owner_cut_recursion_prototype_codex_20260628.out`,
and
`07-reflections/lrc14-charal-owner-cut-recursion-prototype-codex-20260628.md`.

Purpose: continue incoming HYP-3410 by turning the Menger/charal signal into a
finite owner-cut decision-tree API.  The exact known cuts are:

```text
height leak:               cut size 1, core 5:g1
persistent owner leak:     cut size 1, core 1:g1
10->20 frontier leak:      cut size 3, empty core
```

The `10->20` frontier has optimal tree
`13:g1? -> positive-open else 11:g1? -> positive-open else 2:g2? -> positive-open else unit-petal`.
This rules out a one-label owner theorem on current evidence and points to a
bounded owner-cut recursion theorem.

Next: extend the HYP-2963/HYP-3406 bank beyond `(72,20)` and measure minimal
owner-cut growth.  If cut size stays bounded, try a finite Menger/Farkas
owner-cut theorem.  If it grows, classify the first growth as
Schwarz-Christoffel accessory debt, tropical/off-grid debt, state-lift debt,
height-factor debt, exact-period/BDH exception, or named residual.
## Lead codex-2026-06-28: Off-grid resonance transparency and Rprime closure

**Status:** SYNTHESIS / exact-rational transparency scout; not proof
(HYP-3421/T1382/LTI-382/LTT-282), companion to HYP-3415 and HYP-3418.

Claimed files:
`05-knowledge/hypotheses/HYP-3421-lrc14-offgrid-resonance-transparency-rprime-closure.md`,
script
`04-computation/lrc14_offgrid_resonance_transparency_rprime_closure_codex_20260628.py`,
result
`05-knowledge/results/lrc14_offgrid_resonance_transparency_rprime_closure_codex_20260628.out`,
and
`07-reflections/lrc14-offgrid-resonance-transparency-rprime-closure-codex-20260628.md`.

Purpose: extend HYP-3415's critical-path map by proving the resonant part of
the one remaining floor inequality is a finite full-optimum off-grid
transparency problem, not a separate analytic survivor, while retaining
HYP-3418's correction that the covering floor is still 2-adic.

Readout: named optima are off the `14`-grid with visible floors `1/12`,
`1/8`, and `1/9`; the canonical `{1..11,13,84m}` tower has
`t=(35m+2)/(84m+5)`, `M=7m/(84m+5)>1/14`; all checked resonant speeds
(multiples of `2` or `7`, including `14Q` tips) are safe at the selected
off-grid optima.  Tournament vertices are proof carriers, not runners or
residues.

Next: prove the all-packet transparency classifier, then prove HYP-3418's
2-adic even-speed descent and finish the symbolic HYP-3129/HYP-3140
`Rprime = E[N_R | Q]/E[N_R]` constant chase.  The classifier
must route every residual packet to q-witness, denominator-floor transparency,
the canonical binding formula, 2-adic descent, positive owner packet,
signed-SPEC/fiber-PGF certificate, or named terminal debt.

## Lead codex-2026-06-28: Recursive sidecar pattern atlas

**Status:** SYNTHESIS / executable recursion-pattern router; not proof
(HYP-3409/T1370/LTI-370/LTT-270).

Claimed files:
`05-knowledge/hypotheses/HYP-3409-lrc14-recursive-sidecar-pattern-atlas.md`,
script
`04-computation/lrc14_recursive_sidecar_pattern_atlas_codex_20260628.py`,
result
`05-knowledge/results/lrc14_recursive_sidecar_pattern_atlas_codex_20260628.out`,
and
`07-reflections/lrc14-recursive-sidecar-pattern-atlas-codex-20260628.md`.

Purpose: turn the current abstract LRC14 proof pattern into a reusable finite
lemma work queue.  The shared recursion is:

```text
legal quotient -> mixed theorem-exit fiber -> first missing sidecar
-> repaired quotient -> next quotient
```

Readout: the top recursion operators are mixed-fiber resurrection, owner-cut
recursion, collar-to-bank lift, height-then-owner escalation, and finite
chamber terminal routing.  Tournament vertices are recursion operators/proof
obligations, not runners, raw arcs, residues, or constants.

Next: implement a shared quotient/fiber/repaired-quotient API for HYP-3405
AP-vs-`13->27` and HYP-3406 owner leaks; build the owner-support Menger graph
for `petal 13->26` and `petal 10->20`; extend HYP-3406 past `(72,20)` until
`residue+owner_support` first fails or supports a finite owner-cut theorem.
Add terminal-exit labels before using BDH/Mertens-style averages.

## Lead codex-2026-06-28: Special-function cut signature recursion

**Status:** EVIDENCE / executable creative synthesis; not an LRC14 proof
(HYP-3412/T1373/LTI-373/LTT-273).

Claimed files:
`05-knowledge/hypotheses/HYP-3412-lrc14-special-function-cut-signature-recursion.md`,
script
`04-computation/lrc14_special_function_cut_signature_recursion_codex_20260628.py`,
result
`05-knowledge/results/lrc14_special_function_cut_signature_recursion_codex_20260628.out`,
and
`07-reflections/lrc14-special-function-cut-signature-recursion-codex-20260628.md`.

Purpose: translate a deliberately broad prompt set--Bring radical,
Schwarz-Christoffel maps, Barban-Davenport-Halberstam variance, Menger cuts,
Ramanujan-Soldner zero-points, Sophie Germain quartic splitting,
Hermite-Lindemann-Weierstrass separation, Krasner stability, and
Meissel-Mertens residuals--into testable LRC14 sidecar-recursion signals.

Readout: on the HYP-3406 `(72,20)` expanded bank (`2431` rows), residue alone
leaves `3` mixed theorem-exit fibers; residue+`v2` and residue+exact height
each leave `2`; BDH variance leaves `3`; cut angle, Krasner radius, and full
owner support each leave `0`; Sophie quartic and honest Bring branch alarm
each leave `2`, while PGF/root proxy leaves `1`.  The `14`-row owner leak,
`12`-row petal `10->20` owner/height-persistent leak, and `3`-row height leak
all have minimum one-sidecar covers by `SC_cut_angle`, `Krasner_radius`, or
`owner_support`.  Motif tournament path:
`M01 -> M02 -> M03 -> M04 -> M05 -> M06 -> M07 -> M08 -> M09 -> M10`.

Next: enlarge beyond `(72,20)` and stress `residue + cut_angle_word` and
`residue + krasner_radius_word`.  If either fails, record the first collision
and decide whether full endpoint-owner support, exact cut labels, or PGF/root
branch payload is the true next sidecar.
Next: implement the scout by replaying the HYP-3406 height and endpoint-owner
leaks as controlled-forgetting rows; then compute recursive
`charal_signature` mixed-fiber counts, a toy Menger exit-cut certificate,
motif-transfer rankings, and Tournament Analysis over proof carriers rather
than runners.
## Lead codex-2026-06-28: Owner-cut chiral transcendence synthesis

**Status:** SYNTHESIS / exact HYP-3406 owner-cut audit plus creative carrier
ranking; not proof (HYP-3420/T1381/LTI-381/LTT-281).

Claimed files:
`05-knowledge/hypotheses/HYP-3420-lrc14-owner-cut-chiral-transcendence-synthesis.md`,
script
`04-computation/lrc14_transcendence_cut_chiral_synthesis_codex_20260628.py`,
result
`05-knowledge/results/lrc14_transcendence_cut_chiral_synthesis_codex_20260628.out`,
and
`07-reflections/lrc14-owner-cut-chiral-transcendence-synthesis-codex-20260628.md`.

Purpose: translate the prompt's Bring radical, Schwarz-Christoffel,
Barban-Davenport-Halberstam, Menger cut, chiral-signature, Ramanujan-Soldner,
Sophie Germain, Hermite-Lindemann-Weierstrass, Krasner, Meissel-Mertens, and
`exp(exp(exp(79)))` ideas into proof packet fields after HYP-3406.

Readout: rebuilding the HYP-3406 expanded banks shows that
`residue_plus_owner_chiral_class` and `residue_plus_owner_support` kill all
mixed theorem-exit fibers on `(20,4),(30,8),(48,12),(60,16)`.  On the largest
bank, the BDH-style pair-disagreement variance drops from residue-only `12`
to owner/chiral-owner `0`.  Both residue-mixed fibers have size-one owner cuts:
`('1:g1',)` for the `petal 13->26` versus positive-open `26/40/54` family,
and `('5:g1',)` for the `P10+GW` / `GW-shell` / `12->48` fiber.

Next: prove the endpoint-owner Menger-cut theorem for enlarged actual-packet
fibers, or find the first failure.  The first failure should be routed to
tropical/off-grid debt, unit-contact holonomy, state-lift debt, or
Bring-style branch/monodromy debt.  Keep Tournament Analysis on proof carriers:
owner cuts, chiral owner recursion, BDH variance, Krasner stability,
Schwarz-Christoffel owner polygons, quartic factor gates, and branch guards;
do not use raw constants or raw exp-tower scale as vertices.

## Lead codex-2026-06-28: Bring/Schwarz/BDH/Menger charal recursion

**Status:** SYNTHESIS / exact mixed-fiber sidecar scout; not proof
(HYP-3410/T1371/LTI-371/LTT-271).

Claimed files:
`05-knowledge/hypotheses/HYP-3410-lrc14-bring-schwarz-bdh-menger-charal-recursion.md`,
script
`04-computation/lrc14_bring_sc_bdh_menger_charal_recursion_codex_20260628.py`,
result
`05-knowledge/results/lrc14_bring_sc_bdh_menger_charal_recursion_codex_20260628.out`,
and
`07-reflections/lrc14-bring-schwarz-bdh-menger-charal-recursion-codex-20260628.md`.

Purpose: realize the Bring/Schwarz/BDH/Menger slice of the HYP-3407
boundary/special-function route on
the user's Bring radical, Schwarz-Christoffel, Barban-Davenport-Halberstam,
Menger-cut, and charal-signature prompts by converting them into exact
first-failure packet sidecars over HYP-3406, not into raw analogy.

Readout: the height leak has a one-label owner cut `5:g1` and top finite-BDH
variance `5:g1=8/9`.  The persistent owner leak has a one-label owner cut
`1:g1`, with `1:g1`, `13:g1`, and `11:g1` the strongest finite variance
channels.  The `(72,20)` `10->20` frontier has minimum owner cut size `3`,
with first cuts including `('11:g1','13:g1','1:g1')` and top variance
`13:g1=49/50`.  The `+14` ladders preserve positive-open exits for `1/3/5 ->
26,40,54`, retain apex-bearing turns for `12 -> 48,132`, and isolate petal
exits at `13 -> 26` and `10 -> 20`.

Next: enlarge the HYP-2963 mixed-fiber bank and test the recursive owner-cut
theorem.  Every mixed first-failure fiber should become exit-pure under charal
recursion, split by a bounded owner cut, split by finite channel variance,
repaired by Schwarz-Christoffel accessory owner debt, or routed to named
owner/height/off-grid/state-lift debt.  Keep Tournament Analysis on proof
carriers and sidecar transformations, not runners or famous-problem names.

## Lead codex-2026-06-28: Owner-cut resurrection calculus

**Status:** SYNTHESIS / exact clause-certificate calculus over known mixed
fibers; not proof (HYP-3414/T1374/LTI-374/LTT-274).

Claimed files:
`05-knowledge/hypotheses/HYP-3414-lrc14-owner-cut-resurrection-calculus.md`,
script
`04-computation/lrc14_owner_cut_resurrection_calculus_codex_20260628.py`,
result
`05-knowledge/results/lrc14_owner_cut_resurrection_calculus_codex_20260628.out`,
and
`07-reflections/lrc14-owner-cut-resurrection-calculus-codex-20260628.md`.

Purpose: turn the HYP-3409/HYP-3410 owner-cut route into a finite
Menger/Farkas certificate, integrating incoming HYP-3411/HYP-3412/HYP-3413.
A different-exit row pair emits an endpoint-owner support symmetric-difference
clause; a legal owner cut is a minimum hitting set whose cut-code buckets are
theorem-exit pure.

Readout: the height leak has `2` clauses and singleton cut `5:g1`; the
persistent owner leak has `10` clauses and singleton cut `1:g1`; the
`(72,20)` `10->20` frontier has `9` clauses, minimum cut size `3`, five
minimum cuts, empty core, and cut union
`11:g1,13:g1,1:g1,2:g2,5:g1,7:g7`.  This redirects the theorem target from a
singleton-owner separator to a bounded owner-transversal lemma plus terminal
chamber router.

Next: extend HYP-3406 beyond `(72,20)` and run this clause/transversal
calculus on the first `residue+owner_support` failure, if any.  Then test
endpoint deletion, mirror-swap, and `+14` child decks for stability of
transversal number and cut-code theorem-exit purity.

## Lead codex-2026-06-28: AP-collar finite lemma certificate

**Status:** EVIDENCE / exact finite certificate for the HYP-3401 AP-collar
lemma; not proof (HYP-3405/T1366/LTI-366/LTT-266).

Claimed files:
`05-knowledge/hypotheses/HYP-3405-lrc14-ap-collar-finite-lemma-certificate.md`,
script
`04-computation/lrc14_ap_collar_finite_lemma_certificate_codex_20260628.py`,
result
`05-knowledge/results/lrc14_ap_collar_finite_lemma_certificate_codex_20260628.out`,
and
`07-reflections/lrc14-ap-collar-finite-lemma-certificate-codex-20260628.md`.

Purpose: convert the HYP-3401 quotient-mixing scout into a concrete finite
lemma target.  Every AP-collar non-boundary row should have an explicit
strict-open witness, and every quotient failure should name its missing sidecar
rather than remaining an analogy.

Readout: in the AP one-swap collar through replacement speed `84`, there are
`924` rows.  AP and Goddyn-Wong `12->24` are the only boundary-tight rows; the
other `922` rows have exact rational strict-open witness intervals.  The
uniform strict-open mass lower bound is `1/1260`, uniquely at `12->36`.

Compression failure: the HYP-3311 nonunit-height packet has one mixed
boundary/strict fiber of size `31`.  AP and strict-open `13->27` share unit
contact status, C3 skeleton, `Q(sqrt(-7))` character, covering layer, and
nonunit height, but differ by the unit-height lift `(13,0)->(13,1)`.  The
repair matrix says unit-height/full-height/height-completed data kills the
fiber; unit contacts, covering layer, and nonunit height do not.

Next: formalize the finite lemma as a certificate with two boundary atoms,
922 open-interval witnesses, the `1/1260` floor, and the unit-height
obstruction vector.  Then globalize by replacing row-height retention with a
finite chamber theorem routing every height/flex move to AP/GW boundary,
strict-open mass, `Phi14d`, finite Toeplitz/Green/root-motion discharge,
state-lift debt, or named residual.  Keep Tournament Analysis on proof
carriers, not runners, residues, or replacement speeds.

## Lead codex-2026-06-28: Three-coordinate obstruction exactness in the AP collar

**Status:** EVIDENCE / exact AP-collar first-obstruction scout; not proof
(HYP-3401/T1362/LTI-362/LTT-262).

Claimed files:
`05-knowledge/hypotheses/HYP-3401-lrc14-three-coordinate-obstruction-exactness.md`,
script
`04-computation/lrc14_three_coordinate_obstruction_exactness_codex_20260628.py`,
result
`05-knowledge/results/lrc14_three_coordinate_obstruction_exactness_codex_20260628.out`,
and
`07-reflections/lrc14-three-coordinate-obstruction-exactness-codex-20260628.md`.

Purpose: instantiate the HYP-3301 first-obstruction sheaf idea on HYP-3311's
concrete `C3 + Q(sqrt(-7)) + height/flex` packet.  Treat every quotient as a
function from AP-collar rows to compressed packet data; a mixed fiber containing
both boundary-tight and strict-open rows is the explicit compression failure.

Readout: in the AP one-swap collar through replacement speed `84`, there are
`924` rows, with only AP and Goddyn-Wong `12->24` boundary-tight and the other
`922` strict-open.  Raw unit, raw residue, C3, quadratic, C3+quadratic,
C3+quadratic+covering, and C3+quadratic+nonunit-height quotients all mix
boundary-tight and strict-open rows.  The sharp leak is AP versus `13->27`,
which has strict safe mass `13691/582120` while preserving the nonunit height
ledger.  The `height_completed_packet` and full height/residue ledger have
`0` mixed fibers.

Integration with incoming monad-explorer work: its HYP-3311 actual-packet
sheaf instantiation shows nonunit residue data repairs the first coarse
ambiguity on the curated HYP-2969 bank.  This AP-collar scout is the next
stress test: after that nonunit repair, unit-height flex is the first visible
missing coordinate.

Next: formalize the finite AP-collar exactness lemma, then globalize without
making the packet tautological.  Replace full row-height retention by a chamber
theorem: every height/flex perturbation should route to AP/GW boundary,
strict-open mass, `Phi14d` equality, finite Toeplitz/Green/root-motion
discharge, state-lift debt, or named residual.  Keep Tournament Analysis on
sidecar packets rather than runners, residues, or raw replacement speeds.
## Lead codex-2026-06-28: Creative proof-reframe first-failure atlas

**Status:** SYNTHESIS / executable proof-route generator; not proof
(HYP-3404/T1365/LTI-365/LTT-265).

Claimed files:
`05-knowledge/hypotheses/HYP-3404-lrc14-creative-reframe-lead-atlas.md`,
script
`04-computation/lrc14_creative_reframe_lead_atlas_codex_20260628.py`,
result
`05-knowledge/results/lrc14_creative_reframe_lead_atlas_codex_20260628.out`,
and
`07-reflections/lrc14-creative-reframe-lead-atlas-codex-20260628.md`.

Purpose: force the creative LRC14 reframes into a falsifiable queue anchored
on actual packet data, the HYP-3401 AP-collar unit-height obstruction,
HYP-3402 owner-current / tropical-wall sidecars, and HYP-3403 shadow-charge
packet gluing.  The HYP-2969/HYP-2963
actual-packet bank has `31` rows, one mixed coarse
theorem-exit fiber of size `7`, `0` residue-word mixed fibers, `1` v2-word
mixed fiber, and all `7` qdiv>14 rows are `positive-Haar-open`.

Readout: the atlas ranks `15` proof leads as tournament vertices, not runners
or arcs.  The priority path is
`R01 -> R11 -> R14 -> R04 -> R05 -> R02 -> R03 -> R15 -> R06 -> R07 -> R08 -> R09 -> R12 -> R13 -> R10`.
The top route is the residue-word breakpoint theorem: enlarge the sheaf until
nonunit residue-word exactness first fails, then prove no failure or name the
first failure as unit or nonunit height/flex, endpoint owner, off-grid floor,
exact period, K33/H7, state-lift debt, or finite trap.  The extra
forced-exploration routes
that rose high are colored-resonance half-boundary sieve (HYP-2593/HYP-2595)
and endpoint deletion-cut/Menger resurrection, now preceded by HYP-3402's
owner-current/tropical-wall first-leak table and HYP-3403's shadow-charge
controlled-forgetting table.

Next: run the actual-packet sheaf instantiation on broader HYP-2963 residual
samples.  For the first mixed residue-word fiber, compute colored CRT
half-boundary deficit, owner-current word, tropical-wall word, shadow-charge
sidecar, endpoint-owner deletion cut, low odd Faulhaber moment ladder, and
covering-flex Jacobian nullspace before trying a symbolic proof.

## Lead codex-2026-06-28: Shadow-charge conservation atlas

**Status:** SYNTHESIS / executable proof-route router; not proof
(HYP-3400/T1352/LTI-352/LTT-252).

Claimed files:
`05-knowledge/hypotheses/HYP-3400-lrc14-shadow-charge-conservation-atlas.md`,
script
`04-computation/lrc14_shadow_charge_conservation_atlas_codex_20260628.py`,
result
`05-knowledge/results/lrc14_shadow_charge_conservation_atlas_codex_20260628.out`,
and
`07-reflections/lrc14-shadow-charge-conservation-atlas-codex-20260628.md`.

Purpose: look for a more creative unifying frame for the LRC14 proof search by
treating topology, duality, autocorrelation, conductance, discrepancy, root
motion, and lift/descent as shadows of one conserved proof-charge packet.
A quotient is legal only if it preserves witness-bearing payload, transfers it
to a dual shadow, or emits named debt.

Readout: the scout uses proof-charge reservoirs, not runners/arcs/lags/scalars,
as Tournament Analysis vertices.  Incoming HYP-3246/HYP-3252 are integrated
as `index_theorem_degree_charge`: analytic Cech/Euler index, topological
Borsuk-Ulam degree, and Gauss-sum index are candidate shadows of one charge,
with the forcing gap to actual loneliness retained as debt.  Reservoirs are
index theorem, uniform-margin floor, contact-holonomy curvature, cyclotomic
witness address, tiling lift/descent, Cech/Euler hole, `D_7` Borsuk-Ulam sign,
Lee-Yang root motion, state-lift obstruction, Green/Dirichlet current, bulk discrepancy/Hensel density,
normal-fan/Toeplitz slack, autocorrelation transport, law-defect entropy, and
raw scalar shadow.  HYP-3253 contact holonomy is integrated as quotient
curvature: shell-lag residuals must have zero curvature, exact `zeta_7`
holonomy discharge, endpoint-cell lift, or named curvature debt.  HYP-3254,
HYP-3256, and HYP-3258 add the residue/magnitude split and binding/covering
census warning inside the uniform-margin reservoir.  HYP-3265 adds the
unit-contact graph case split before the safety function is scalarized, and
HYP-3300 adds the observability/Morse audit after sidecars attach.  The
tournament has `15` vertices, no directed 3-cycles,
singleton SCCs, and one Hamiltonian path led by
`index_theorem_degree_charge -> uniform_margin_floor_charge -> contact_holonomy_curvature_charge -> cyclotomic_witness_address_charge ->
tiling_lift_descent_charge -> cech_euler_hole_charge`.

Next: convert the HYP-3202 non-AP traps into a charge-discharge table with
autocorrelation transport, Green resistance excess, Toeplitz slack,
normal-fan first failed coordinate, `D_7` sign payload, root-motion class,
contact-holonomy curvature status, tiling descent status, and first proof
exit; then test whether the
named-residual-debt case is empty.  In parallel, test HYP-3246/HYP-3252's
index equality as descriptor plus the HYP-3250 S-dependent floor as proof.
## Lead mac-mini-2026-06-28-S84 + codex: LRC14 CRT/Galois sidecar audit

**Status:** EVIDENCE / exact CRT and Galois-sidecar audit; proof-target
refinement, not proof (HYP-3311/T1361).

Claimed files:
`05-knowledge/hypotheses/HYP-3311-lrc14-crt-galois-sidecar-audit.md`,
script `04-computation/lrc_census_crt_factorization_macmini_S84.py`,
result `05-knowledge/results/lrc_census_crt_factorization_macmini_S84.out`,
and reflection
`07-reflections/the-census-factors-via-crt-7-adic-residue-c3-skeleton-times-2-adic-magnitude-doubling-hinge.md`.

Purpose: make the owner's S84 insight proof-facing as an exact sidecar audit
for HYP-3310's broader C6 residue-magnitude frame and as a small chart for
HYP-3301's first-obstruction sheaf route.  After HYP-3400, it is also the
concrete charge-conservation packet for the C3, quadratic, and height/flex
coordinates.  The exact census split is
`{1,...,13}=U union 2U union {7}` with units
`U={1,3,5,9,11,13}`, even covering shadow `2U={2,4,6,8,10,12}`, and ramified
apex `7`.  The bijection `u -> 2u mod 14` links binding units to even covering
classes, but the packet also retains the `Q(zeta_7)` subfield lattice:
the C3 quotient gives the binding-pair orbit, while the `Q(sqrt(-7))`
quadratic character is transverse and cannot be forgotten.

Readout: the proof should split into unit-contact rigidity plus covering-flex
rigidity.  Unit contacts should be handled by HYP-2909 plus C3 propagation and
HYP-3265's contact graph.  Killed contacts should route through HYP-3300-style
observability/Morse chambers.  The open hard theorem is that the covering
layer `2U+{7}` has only the AP/Goddyn-Wong `12->24` integer tight hinge, with
all same-residue height impostors opening strict safe mass or entering named
finite discharge.

Next: build the finite covering-flex theorem for `2U+{7}` and attach it to
the HYP-3265/HYP-3300 chamber classifier.  Record which quotient preserves the
LRC predicate, which coordinate it destroys, and whether the loss is paid for
by sidecar, descent, AP/GW boundary, `Phi14d` equality, strict open witness, or
named debt.  HYP-3301's extra test: the first lost coordinate among C3,
`chi_7`, and height/flex must be exact, holonomy-repaired, boundary-positive,
AP/GW-forbidden, or named K33/H7 debt.  HYP-3400's extra test: the same lost
coordinate must be preserved, transferred to a dual shadow, or explicitly
charged as debt.

## Lead codex-2026-06-28: Topology/geometry/graph proof-route atlas

**Status:** SYNTHESIS / executable proof-carrier tournament; not proof
(HYP-3243/T1340).

Claimed files:
`05-knowledge/hypotheses/HYP-3243-lrc14-topology-geometry-graph-proof-routes.md`,
script
`04-computation/lrc14_topology_geometry_graph_routes_codex_20260628.py`,
result
`05-knowledge/results/lrc14_topology_geometry_graph_routes_codex_20260628.out`,
and
`07-reflections/lrc14-topology-geometry-graph-proof-routes-codex-20260628.md`.

Purpose: make the visual LRC14 routes load-bearing by treating pictures as
typed proof carriers: circle endpoint arrangements, oriented topes/cocircuits,
Cech safe-component nerves, `D_7` Borsuk-Ulam index packets,
`Phi_14/Phi_{14d}` witness strata, Green conductance graphs,
Toeplitz/Fejer normal-fan faces, Lee-Yang root motion, ear payload graphs,
finite chamber atlases, and state-lift obligations.

Readout: the witness side is now mostly organized by HYP-3240/HYP-3241
(`Phi_14` for AP/GW and `Phi_{14d}` for dilations), so the remaining visual
frontier is tight-locus finiteness, bulk equidistribution, and legal gluing
across the Vitali core.  The carrier tournament has `12` vertices, no directed
3-cycles, singleton SCCs, and Hamiltonian path led by
`oriented_matroid_topes_cocircuits -> circle_endpoint_arrangement ->
cech_nerve_safe_components`.  Proposed schema: every primitive row has an
open safe tope, AP/GW `Phi_14` equality, dilation `Phi_{14d}` equality,
finite Toeplitz/Green/root-motion chamber discharge, state-lift `H=7`
contradiction, or named residual debt.

Next: turn the schema into a finite chamber theorem with endpoint-owner
cocircuits, Cech components, Green/Toeplitz labels, and `D_7` sign sidecars
retained explicitly.

## Lead codex-2026-06-28: Even/odd positive/negative duality bridge

**Status:** EVIDENCE / exact bounded-bank scout plus synthesis; not proof
(HYP-3238/T1338).

Claimed files:
`05-knowledge/hypotheses/HYP-3238-lrc14-even-odd-positive-negative-duality-bridge.md`,
script
`04-computation/lrc14_even_odd_positive_negative_bridge_codex_20260628.py`,
result
`05-knowledge/results/lrc14_even_odd_positive_negative_bridge_codex_20260628.out`,
and
`07-reflections/lrc14-even-odd-positive-negative-duality-bridge-codex-20260628.md`.

Purpose: finalize the previous Green-conductance session by merging it with
the even/odd and positive/negative dualities now visible in HYP-3220,
HYP-3219, and HYP-3237.  HYP-3220 says these are one parity/complement
operator: de Moivre power sums `-1,5,-4,13,-16,38,-57,117` have sign
`(-1)^k` from the negative Perron period `-2cos(pi/7)`.  Even/positive data is Fejer square, SOS magnitude, pair-Pascal cap
the even/odd and positive/negative dualities now visible in HYP-3219,
HYP-3220, HYP-3239, and HYP-3237.  Even/positive data is Fejer square, SOS magnitude, pair-Pascal cap
mass, covariance layers, positive Green conductance, Perron coherent mode,
and bulk positive-measure equidistribution.  Odd/negative data is Worpitzky
associator debt, Brouwer fixed-point or Borsuk-Ulam free-`Z/2` degree sign,
Hermite-Biehler odd leg, negative covariance leakage, signed chart-change
cancellation debt, `D_7` sign-isotypic data, and measure-zero cyclotomic core
witnesses.

Readout: AP has `q0=481/1470`, `q3=26/245`, `q6=1/49`, `q0+q6=73/210`,
`L_y=2633/7350`, and `lambda2=0.192033074001`.  AP is uniquely primitive-tight
for `L_y`, `q0+q6`, and `lambda2`.  The false-terminal audit finds `18`
primitive non-AP rows with zero negative covariance leakage, `2754` primitive
connected positive-graph non-AP rows, `2879` primitive rows with positive `q3`
debt and `0` exchange-margin violations, and the `11` non-AP HYP-3202 traps
split into `8` negative-leakage-plus-odd-debt and `3`
odd-debt-without-negative-leakage.

Post-rebase HYP-3239 integration: the mac-mini S76 branch identifies the two
targets as one bimodal/phi4 extremality under inclusion-exclusion parity; the
kps S31av branch identifies the `p=7=3 mod 4` sign sidecar as the `D_7`
sign representation / free-`Z2` Borsuk-Ulam packet.  Track this family tag
when deciding whether the odd/negative sidecar is discharged.

HYP-3241 adds the matching index field: for n=14 the AP `Phi_14` witness set
has `3` antipodal pairs, equal to `(p-1)/2` and the de Moivre degree.  The
parity of this equioscillation index selects Borsuk-Ulam/free-`Z2` versus
Brouwer/SOS in the family law.

Compression guardrail: an even/positive quotient is proof-grade only when the
odd/negative payload is zero, reconstructible, dual-annihilated, or retained
as sidecar.  This is the HYP-3201 law-defect rule applied to associativity,
positivity, evenness, measure, and scalarization failures.

Next: prove the `q3` exchange-rate inequality symbolically from HYP-3204 or
the shell Delsarte dual, attach HYP-3222 Hermite-Biehler interlacing to the
three zero-negative-leakage traps with odd debt, and compare HYP-3220's
negative Perron sign with HYP-3236 Fiedler/current data.

## Lead codex-2026-06-28: Cyclotomic Delsarte/Beurling-Selberg magic function

**Status:** EVIDENCE / exact bounded-bank synthesis; not proof
(HYP-3228/T1326).

Claimed files:
`05-knowledge/hypotheses/HYP-3228-cyclotomic-delsarte-beurling-selberg-magic.md`,
`04-computation/lrc_cyclotomic_delsarte_bs_magic_codex_20260628.py`,
`05-knowledge/results/lrc_cyclotomic_delsarte_bs_magic_codex_20260628.out`,
and
`07-reflections/cyclotomic-delsarte-beurling-selberg-magic-codex-20260628.md`.

Readout: the k=8 `L_y` dual is the shell magic polynomial
`((n-1)(n-2)(n-4)(n-5))/4`, values `[10,0,0,1,0,0,10]`, with
Delsarte/Newton form `10S0-10S1+10S2-9S3+6S4` and Joukowski form
`10(u^3-3u)+1`.  It has no primitive bounded-bank beaters and only the known
doubled AP all-bank tie.  AP support has the same equality set and positive
deficit-ratio control, but the centered directions are not identical.  Guardrail:
cyclic nonnegative completion would need central coefficient `rho>=18.019...`,
so raw Fejer/PSD positivity is too expensive.  After HYP-3214, treat the
positive sector-side Fejer kernel `F_7` and this shell `L_y` Delsarte dual as
distinct cyclotomic magic faces to glue; use ordered-tail exchange and
Toeplitz/trap sidecars for the center.

Next: decompose the magic deficit into AP-support gap, HYP-3204 exchange-rate
slack, and HYP-3224 Toeplitz/covariance trap-discharge slack; test whether the
residual has a small sign-controlled Delsarte basis.
## Lead codex-2026-06-28: Equioscillation / autocorrelation atlas

**Status:** EVIDENCE / exact motif scout plus synthesis; not proof
(HYP-3245/T1309).

Claimed files:
`05-knowledge/hypotheses/HYP-3245-lrc14-equioscillation-autocorrelation-atlas.md`,
`04-computation/lrc14_equioscillation_autocorrelation_atlas_codex_20260628.py`,
`05-knowledge/results/lrc14_equioscillation_autocorrelation_atlas_codex_20260628.out`,
and
`07-reflections/lrc14-equioscillation-autocorrelation-atlas-codex-20260628.md`.

Readout: HYP-3214's Fejer kernel is not just adjacent to autocorrelation:
`F_7=(de Moivre cubic)^2` is the Chebyshev double-zero/equioscillation kernel,
positive-definite Delsarte kernel, and triangular AP autocorrelation.  HYP-3245
uses that equality model to define an out-correlation residual for HYP-3202
traps.  Every non-AP local trap moves support-autocorrelation mass from lags
`1..7` to lags `8..14`, with total residual zero.  The coarse transport is
universal; the ripple type follows HYP-3225's Green/Rayleigh/Plucker classes.
Incoming HYP-3227 is the conductance-graph companion, showing the same trap
boundary remains connected to legal certificate exits.
Incoming HYP-3228 is the shell-magic companion: it names the `L_y` contact
vector whose deficit should be tested against the same lag-transport signal.
Incoming HYP-3234 adds signed chart-change labels and HYP-3233 adds
cyclotomic-factor labels, so the next lag scout must retain chart debt and
Phi_7 debt instead of reporting a single scalar transport.  Incoming HYP-3235
and HYP-3218 add cap-field conductor, Fejer-square, Gauss-sum margin, and AP
self-duality labels.
Incoming HYP-3229, HYP-3230, and HYP-3231 add the next labels to preserve:
modular coefficient/subshift sidecars, three-gap/Farey recursion, and
scale-normal survival status.
Incoming HYP-3232 and HYP-3217 add two more labels: apex-fold side and
cyclotomic mode.  Incoming HYP-3236, HYP-3237, and HYP-3219 add Green
conductance, Vitali/Brouwer core-wall, and Brouwer trace-sign/SOS labels:
`green_resistance_slack`, `lambda2_conductance_rank`,
`negative_covariance_leakage`, `thomson_current_profile`,
`fiedler_bottleneck_id`, `vitali_wall_side`, `brouwer_saddle_sign`,
`phi14_core_witness`, `core_bulk_transport_status`, `brouwer_trace_sign`,
`degree_sos_factorization`, and `even_odd_bonferroni_node_slack`.  Incoming
HYP-3238 and HYP-3239 add even-positive/odd-negative compression and
dihedral/Borsuk-Ulam sign-irrep labels: `positive_negative_duality_status`,
`odd_negative_payload_reconstruction`, `dihedral_irrep_label`,
`complement_antiautomorphism_sign`, `borsuk_ulam_index`,
`imaginary_gauss_sum_sign`, and `phi4_bimodal_extremizer_rank`.  Incoming
HYP-3240 and HYP-3241 add core-witness arithmetic labels:
`equioscillation_saddle_index`, `phi14_core_universality_status`,
`dilation_witness_grid`, `core_witness_break_reason`, and
`imaginary_norm_route_status`.  Incoming HYP-3242 adds topology labels:
`danger_nerve_euler_characteristic`, `lonely_hole_status`,
`cech_betti_sidecar`, `topological_shadow_class`, and
`cover_hole_witness_pair`.
Incoming HYP-3243 adds proof-carrier labels:
`oriented_matroid_tope_status`, `circle_endpoint_arrangement_cell`,
`cech_safe_component_rank`, `finite_chamber_schema_status`,
`state_lift_H7_obstruction`, and `proof_carrier_tournament_rank`.
Incoming HYP-3244 adds controlled-forgetting labels:
`tiling_witness_lift_status`, `half_tiling_descent_certificate`,
`path_presentation_fiber_weight`, `parent_aut_word_orbit_id`,
`rectangle_hourglass_residue`, `tail_tip_deletion_signature`, and
`controlled_forgetting_span_status`.
Executed HYP-3226 adds the motif-ledger labels:
`small_pattern_motif_id`, `payload_preserved`, `payload_destroyed`,
`repair_sidecar`, and `terminal_risk_label`.

Next: prove or refute the signed transport lemma on the full bounded k=8
bank.  Compare low-lag deficit and outward surplus with AP support,
Toeplitz `lambda_min`, HYP-3228 magic deficit, HYP-3227 conductance-graph
exits, HYP-3229 modular sidecar labels, HYP-3230 three-gap labels,
HYP-3231 scale-recursion labels,
HYP-3236 Green conductance labels, HYP-3237 Vitali/Brouwer labels,
HYP-3219 Brouwer-sign/SOS labels,
HYP-3238 even-positive/odd-negative labels, HYP-3239 dihedral/Borsuk-Ulam labels,
HYP-3241 core-index labels, HYP-3240 dilation-witness labels, HYP-3242
cover-hole labels, HYP-3243 finite chamber/state-lift labels, HYP-3226
motif-ledger labels, HYP-3244 tiling/half-tiling descent labels, ordered-tail
`q0+q6`, and HYP-3204 exchange-rate slack.
Keep the two-clock guardrail from HYP-3214: 7-sector Fejer for coverage/LHS,
14-clock Johnson pair-Pascal for cap/RHS.

## Lead codex-2026-06-28: Small-pattern adjacency atlas

**Status:** EXECUTED / evidence atlas; not proof (HYP-3226/T1324).

Files:
`05-knowledge/hypotheses/HYP-3226-lrc14-small-pattern-adjacency-atlas.md`,
`04-computation/lrc14_small_pattern_adjacency_atlas_codex_20260628.py`,
`05-knowledge/results/lrc14_small_pattern_adjacency_atlas_codex_20260628.out`,
and
`07-reflections/lrc14-small-pattern-adjacency-atlas-codex-20260628.md`.

Goal: assemble a large but typed list of small motifs adjacent to the current
LRC14 frontier, then classify each by payload preserved, payload destroyed,
repair sidecar, normal-fan coordinate, trap-boundary use, and terminal-proof
risk.  This is meant to connect HYP-3224's normal fan, HYP-3223's
Green-current/Lorentzian classifiers, HYP-3222's HB/Perron exact legs,
HYP-3221's analytic-equidistribution guardrail, HYP-3225's Green-current /
Lorentzian trap fingerprints, HYP-3227's conductance/Fiedler graph, S75's
Gram-kernel/peeling build, HYP-3215's proof-route audit, HYP-3228's shell
magic dual, HYP-3229's modular sidecar audit, HYP-3230/HYP-3231/HYP-3216's
recursion layer, HYP-3232/HYP-3217's mode lattice,
HYP-3233/HYP-3234/HYP-3218/HYP-3235's cyclotomic factor / chart / Fejer /
conductor layer, HYP-3236/HYP-3219/HYP-3237's Green / Brouwer / Vitali
boundary layer, HYP-3220/HYP-3238/HYP-3239's D7 / Borsuk-Ulam /
`p mod 4` / bimodal phi4 diagonal layer, HYP-3241/HYP-3240 universal `Phi_14`
saddle-index / `Phi_{14d}` dilation-witness core, HYP-3242's danger-cover
nerve / Euler-hole certificate, and HYP-3243's finite chamber carrier atlas, without treating
raw numerology as a proof carrier.
Readout: 103 motifs across 102 families were scored by proof-payload retention
over 8395 repo-local files.  The strongest motifs were comb-overlap Gram kernel,
universal `Phi_14` saddle-index core, shell `L_y` magic quartic, normal-cone dual slack,
multi-chart proof split, finite chamber carrier atlas, AP Green algebraic-connectivity certificate,
bimodal phi4 diagonal extremizer, AP self-dual Fejer equidistribution certificate,
three-gap/Stern-Brocot cap-kernel recursion,
danger-cover nerve hole certificate,
consecutive plus doubled AP, modulus-covariance apex break,
D7 Borsuk-Ulam sign-irrep certificate, Toeplitz lambda-min
margin, certificate-Helly separation, single-arc peeling recursion,
ordered-tail exchange-rate ratio, and
D1/D2/D3 covariance layers, plus the p mod 4 imaginary-quadratic wall.  The
useful small patterns cluster around seven
payloads: normal-fan exposure, Toeplitz/moment curvature, covariance layers,
ordered-tail pricing, finite trap discharge, HB/Perron gluing, and analytic
equidistribution.  HYP-3225 now supplies the first Green/Lorentzian trap
fingerprint table, and HYP-3214 turns the Fejer/Chebyshev motif into the
explicit positive-definite `F_7` kernel.  HYP-3227 adds M072, the
conductance/Fiedler trap graph, as a finite-discharge sidecar rather than a
terminal scalar.  S75 adds M073-M075 (Gram kernel, speed-1 peeling, order-3
overlap residues), while HYP-3215 adds M076-M079 (induction base,
Chen-Cusick floor-to-1/14, polyhedron/zonotope flatness, Rosenfeld Node-3).
HYP-3228/HYP-3229 add M080-M084: shell magic `10q0+q3+10q6`,
Gamma0(7) Eisenstein coefficients, Beraha/Mahler height, subshift transfer,
and Dirichlet-L/Stark denominator guardrails.
HYP-3230/HYP-3231/HYP-3216 add M085-M088: three-gap/Stern-Brocot cap-kernel
recursion, scale-normal packet recursion, the `LRC(2p)` moment-order ladder,
and the 2-adic reflection fold.
HYP-3232/HYP-3217 add M089-M090: modulus-covariance apex break and the
cyclotomic subfield / character-mode lattice with the cubic de Moivre mode.
HYP-3233/HYP-3234/HYP-3218/HYP-3235 add M091-M094: cyclotomic factor
grading, signed address chart-change debt, AP self-dual Fejer/Vaaler
equidistribution certificate, and the totally-real cap field conductor packet.
HYP-3236/HYP-3219/HYP-3237 add M095-M097: AP Green
algebraic-connectivity/resistance certificate, Brouwer trace-sign times SOS
split, and the Vitali bulk-core `Phi_14` witness wall.
HYP-3220/HYP-3238/HYP-3239 add M098-M100: D7 Borsuk-Ulam sign-irrep
certificate, p mod 4 imaginary-quadratic family law, and bimodal phi4
cumulant diagonal.
HYP-3241/HYP-3240 add M101: universal `Phi_14` saddle-index core with three antipodal
witness pairs and explicit `Phi_{14d}` dilation witnesses.
HYP-3242 adds M102: the danger-cover nerve / Euler-characteristic hole
certificate, making the LRC witness a retained topological hole rather than a
scalar cap-only readout.
HYP-3243 adds M103: the finite chamber carrier atlas, making visual routes
usable only when each primitive row has an open tope, equality core, discharge,
state-lift exit, or named residual debt.
Famous-problem analogies remain sidecars until they name the LRC coordinate
they preserve and the coordinate they destroy.  Incoming S283's
Skewes/Helfgott-Ruzsa/Collatz/PFR
additive-resonance layer is absorbed here as typed sidecar material, not as a
terminal shortcut.

Next: prove the HYP-3225/HYP-3227 trap table symbolically and add
Gram-kernel PSD, speed-1 peeling, order-3 overlap, conductance/Fiedler cuts,
precision M-matrix or Schur-complement debt, induction-base status,
Chen-Cusick floor-to-1/14 lift status, shell `L_y` magic slack, Gamma0(7)
coefficient-row compatibility, three-gap kernel-recursion status,
scale-normal `omega_Q` exactness, moment-order / 2-adic fold status,
modulus-covariance / subfield-mode status, cyclotomic factor grading,
signed-address chart-change status, AP self-dual Fejer/Vaaler tail status,
totally-real cap field conductor/trace status, Green lambda2/Kirchhoff
resistance status, Brouwer trace-sign/SOS split status, Vitali bulk-core
`Phi_14` witness status, D7 Borsuk-Ulam sign-irrep status, p mod 4
imaginary-quadratic family-law status, bimodal phi4 diagonal status,
universal Phi14 saddle-index / dilation-witness status, and
Fejer/Delsarte `F_7` slack as finite-boundary columns.

## Lead codex-2026-06-28: Joukowski-Hermite-Biehler / Perron-Frobenius synthesis

**Status:** SYNTHESIS / exact local scout; not proof (HYP-3222/T1306).

Claimed files:
`05-knowledge/hypotheses/HYP-3222-joukowski-hb-perron-exact-certificates.md`,
`04-computation/lrc14_joukowski_hb_perron_synthesis_codex_20260628.py`,
`05-knowledge/results/lrc14_joukowski_hb_perron_synthesis_codex_20260628.out`,
and
`07-reflections/joukowski-hb-perron-frobenius-synthesis-codex-20260628.md`.

Readout: after the incoming HYP-3210/HYP-3211 bridge work, and now beside
HYP-3212/HYP-3213 Chebyshev/cyclotomic arithmetic, HYP-3221's apex-7
obstruction, HYP-3205's spectral-dictionary compatibility layer, and
HYP-3204's ordered-tail exchange split, the k=8 spectral packet has two exact
local certificates.  On
the odd/HB side, the even fold gives `E(x)=x^2+5x+4` and the Worpitzky
Eulerian leg gives `O(x)=A_3(x)=x^2+4x+1`; the roots strictly interlace and
`E O'-E' O=(x+3)^2+2>0`.  On the even/PF side, HYP-3202's distance layers
build an ideal nonnegative C6 quotient with
`lambda0=(1^T C1)/6=6237419/25930800=lambda_max`, while a signed AFM contrast
moves the top mode to `k=3`.

Next: run a boundary-aware version of the Perron alignment test on the actual
empty-sector covariance matrices, and a Joukowski/HB lift that reports
self-inversive defect, interlacing status, and Wronskian orientation across
the bounded k=8 bank.  Join those columns to Toeplitz `lambda_min(T)`,
distance-layer dominance, spectral-dictionary compatibility, random-current
order, ordered-tail exchange-rate fields, and HYP-3201 law-defect fields.
## Lead codex-2026-06-28: Green-current / Lorentzian exchange certificates

**Status:** SYNTHESIS / proof-angle proposal; not proof (HYP-3223/T1323).

Claimed files:
`05-knowledge/hypotheses/HYP-3223-lrc14-green-current-lorentzian-exchange-angles.md`
and
`07-reflections/lrc14-green-current-lorentzian-exchange-angles-codex-20260628.md`.

Readout: the latest k=8 evidence now asks for certificates, not another raw
scalar.  Angle A reads the empty-sector covariance matrix as a Green kernel:
`Sigma kappa_2=1^T C 1` is all-ones current energy, distance layers are
boundary conductance channels, exchange moves are Schur-complement or
star-mesh edits, and the HYP-3202 traps are finite bottleneck networks.  Angle
B reads co-emptiness probabilities as a set function whose Rayleigh
differences, third cumulants, AP support direction, and compression traps are
finite differences in a Lorentzian / valuated-matroid chamber.
Rebase integration: HYP-3211/HYP-3212/HYP-3221 say this certificate should
stay on the additive/cyclotomic LRC face, respect the Chebyshev/Delsarte
magic-function cap, and avoid becoming a config-blind replacement for the
analytic equidistribution obligation.  HYP-3205 adds the spectral-dictionary
compatibility layer: treat these Green-current and Lorentzian fields as
candidate coordinates for the same AP-tight face, and compare their trap
discharges with the first failed dictionary coordinate.
Close-out namespace note: HYP-3222 is now the Joukowski/Hermite-Biehler/Perron
exact-certificate packet on mainline, so this Green-current/Lorentzian packet
uses HYP-3223/T1323/LTI-323/LTT-223.

Next: on the `12` arbitrary-swap local maxima from HYP-3202, emit
`effective_resistance_profile`, `thomson_energy_slack`,
`schur_complement_exit`, `trap_network_bottleneck_id`,
`rayleigh_difference_matrix`, `lorentzian_hessian_signature`,
`valuated_exchange_slack`, and `tropical_plucker_defect`.  Check whether the
`11` decoys collapse to a few network/circuit types and whether HYP-3203's AP
support normal is the same exposed chamber normal; also check whether HYP-3205
predicts the same trap discharge.

## Lead codex-2026-06-28: Green-Lorentzian trap fingerprint classifier

**Status:** EVIDENCE / exact finite trap-neighborhood scout; not proof
(HYP-3225/T1308).

Claimed files:
`05-knowledge/hypotheses/HYP-3225-lrc14-green-lorentzian-trap-fingerprints.md`,
`04-computation/lrc14_green_lorentzian_trap_fingerprints_codex_20260628.py`,
`05-knowledge/results/lrc14_green_lorentzian_trap_fingerprints_codex_20260628.out`,
and
`07-reflections/lrc14-green-lorentzian-trap-fingerprints-codex-20260628.md`.

Readout: HYP-3225 executes the previous lead's next scout on the `12`
HYP-3202 arbitrary-swap local maxima plus all one-swap neighbors (`577`
evaluated rows).  Every local maximum selects `Toeplitz_lambda_min` as first
HYP-3224/HYP-3205 dictionary discharge.  The new content is a sidecar split:
`6` rank-2 pair-Plucker bottlenecks, `2` low-connectivity Green bottlenecks,
`2` AFM/frustrated high-Rayleigh-debt rows, and `1` mixed sidecar.  AP itself
has raw Rayleigh-negative and pair-Plucker signals, so those sidecars are not
terminal scalar objectives; they are legal only as relative trap classes under
the Toeplitz/AP normal-cone chart.

Next: turn the finite trap table into exact inequalities for those five
classes; test the taxonomy at k=9/k=10 and in larger neighborhoods; express
pair-Plucker bottlenecks as valuated-matroid circuits and Green bottlenecks
as Schur/Verblunsky/Fejer-Riesz slack; then check whether HYP-3204's central
exchange-rate slack is a projection of the same Toeplitz chart.

## Lead codex-2026-06-28: Green conductance / algebraic connectivity certificate

**Status:** EVIDENCE / exact bounded-bank scout; not proof (HYP-3236/T1336).

Claimed files:
`05-knowledge/hypotheses/HYP-3236-lrc14-green-conductance-connectivity.md`,
`04-computation/lrc14_green_conductance_connectivity_codex_20260628.py`,
`05-knowledge/results/lrc14_green_conductance_connectivity_codex_20260628.out`,
and
`07-reflections/lrc14-green-conductance-connectivity-codex-20260628.md`.

Readout: HYP-3236 executes HYP-3223's electrical route.  It maps the
empty-sector covariance matrix to the positive-part conductance graph on the
six inner sectors, retains negative covariance as leakage sidecar, and reads
the graph through Laplacian `lambda2`, Green kernel `L^+`, Kirchhoff index,
distance effective-resistance channels, bottleneck unit-current profiles, and
current entropy.  AP/consecutive and doubled AP are the only all-bank rows
maximizing algebraic connectivity and total positive conductance, and the only
rows minimizing Kirchhoff/mean/max/distance-layer Green resistance; primitive
normal form leaves AP unique.  AP has `lambda2=0.192033074001`,
`kirchhoff=108.654718079151`, and `maxR=9.713313375596`.  HYP-3202's `11`
non-AP exchange traps all show Green resistance excess, split primarily as
`kirchhoff_excess:3` and `maxR_excess:8`.

Rebase integration: incoming HYP-3225/T1308 executes the trap-local
Green/Lorentzian fingerprint scout and refines these traps into
rank-2 Plucker, low-connectivity Green, AFM/Rayleigh, and mixed sidecar
classes, with Toeplitz `lambda_min` still the universal first discharge.
HYP-3236 is therefore the global all-bank conductance face; HYP-3225 is the
finite-boundary taxonomy; HYP-3226 should ingest both as motifs.

Mainline also claimed HYP-3227/T1325 for the broader sector/precision/trap
conductance graph family.  Read HYP-3236 as the sharper positive-covariance
Green-resistance extremality face inside that family.

Newer mainline claimed HYP-3228/T1326 for the cyclotomic Delsarte shell-magic
`L_y` coefficient face.  Compare HYP-3236 Green resistance slack to that shell
magic deficit, but keep the quotient maps distinct until a packet map is
proved.

Latest mainline claimed HYP-3229 for the modular/arithmetic magic sidecar
audit.  Treat its Gamma0(7), Dirichlet-L, Beraha/Mahler, and comb-overlap
signals as coefficient engines that may feed HYP-3236 only through named
sidecars.

Second rebase integration: incoming HYP-3214 identifies the 7-sector
Fejer/de Moivre magic function as AP autocorrelation and a positive-definite
Delsarte kernel, while separating it from the 14-clock Johnson pair-Pascal
cap.  Treat HYP-3236's conductance graph as the finite Green/Dirichlet sibling
of that harmonic face, not as a scalar replacement for either Fejer or cap.

Newest mainline renumbered the scale-invariance recursion ledger to HYP-3231
and added HYP-3232 for the Mobius/Eisenstein/Legendre interlocking-recursion
audit, with HYP-3230 remaining the three-gap/Farey cap-kernel thread and
HYP-3216 carrying the two-recursion LRC(2p) family law.  Use them as a
scale-compatibility test: HYP-3236 should lift from raw `lambda2` to a
cut/Rayleigh/Thomson packet indexed by the three-gap/Farey cap-kernel address,
or explain which Green slack is annihilated by the moment-order ladder,
2-adic fold, and apex-half break.

Latest incoming HYP-3217 adds the cubic/subfield-lattice mode: the de Moivre
angles are the cubic Gaussian-period cosets inside the `Q(zeta_7)` tower.
Next Green tests should therefore project Fiedler vectors, bottleneck currents,
and distance-resistance channels onto the cubic cosets; otherwise the
conductance graph may preserve scale data while losing the cyclotomic mode.

Newest incoming HYP-3233 grades the recursion modes by cyclotomic factors
`(x-1)^depth * Phi_d`, with the hard apex mode in the real `Phi_7`/de Moivre
cubic factor.  HYP-3236 should therefore tag Laplacian eigenvectors, Green
resolvent poles, Fiedler bottlenecks, and resistance channels by which
cyclotomic factor they preserve or forget.

The incoming HYP-3234 signed-address chart-change sheaf adds the coordinate
guardrail: carry signed chart, local slot basis, chart-change map, and
cancelled-slot debt before moving Green currents between recursion charts.

Incoming HYP-3235 and HYP-3218 add the totally-real cyclotomic cap and Fejer
proof-push.  Compare Green kernels and resistance slack directly to the Fejer
square / Gauss-sum margin; a proof-grade Green coordinate should factor
through that positive certificate or name its electrical residual.

Next: compare Green slack against HYP-3224's AP-support, Toeplitz,
distance-layer, and ordered-tail slacks on every bounded row.  Emit
Schur-complement reduction words for the `11` non-AP traps and test whether
the Green split is the same finite event as the Lorentzian/tropical-Plucker
split proposed by HYP-3223.  Keep the positive-part compression honest:
negative covariance leakage, Schur boundary terminals, and odd
Worpitzky/Hermite-Biehler debt must remain as sidecars or be discharged.

## Lead codex-2026-06-28-S277: function-compression degree-4 guardrail
## Lead codex-2026-06-28-S279: Law-defect entropy compression

**Status:** SYNTHESIS / exact finite scout; not proof (HYP-3201/T1301).

Claimed files:
`05-knowledge/hypotheses/HYP-3201-law-defect-entropy-compression.md`,
`04-computation/lrc14_law_defect_entropy_compression_codex_20260628.py`,
`05-knowledge/results/lrc14_law_defect_entropy_compression_codex_20260628.out`,
and `07-reflections/law-defect-entropy-function-compression-codex-s279.md`.

Readout: this extends HYP-3150/HYP-3151's factor-through wall from
commutativity to a general law-defect entropy calculus.  A law is a
zero-entropy quotient statement `H(f|q_L)=0`; failed commutativity,
associativity, idempotence, distributivity, action, root-circle, or moment
laws produce ordered, bracket, multiplicity, context, representative/action,
root-variance, or cumulant sidecars.  Exact finite residuals include
`0.816327` bits for ordered exponentiation after unordered-pair compression,
`0.800000` bits for subtraction after bracketing compression, `0.515625` bits
for exponentiation after bracketing compression, and `0.701205` bits for the
K4 fixed-path class action quotient.

Post-fetch integration: renumbered from the initially local HYP-3152/HYP-3161
claims to HYP-3201 after incoming mainline claimed HYP-3152 for the Lee-Yang
circle radius web, HYP-3153 for the Lee-Yang/Worpitzky/quartic packet,
HYP-3154 for the Joukowski/De Moivre bridge, HYP-3160/HYP-3161 for k=8
covariance extremality plus the "1/7 is not a theorem" correction, HYP-3162
for the cyclotomic/Joukowski ideal, HYP-3199 for the n=4 Einheit/minimality
chart, and HYP-3200 for the exact bounded-bank cumulant check.  HYP-3201 uses conditional
entropy as a quotient-defect diagnostic, so it is compatible with the warning
that ordinary row entropy is not the k=8 extremal scalar.

Next: put `law_id`, `law_quotient_map`, `target_function`,
`residual_entropy_bits`, and `sidecar_type` on one HYP-3140 fiber-PGF row, one
HYP-3141 edge-witness row, and the HYP-3142 k=8 packet.  A quotient should be
used only after the target function has zero residual entropy or the named
sidecar is retained.
## Lead codex-2026-06-28: Covariance Laplacian and associator-ear angles

**Status:** SYNTHESIS / proof-angle proposal; not proof (HYP-3163/T1219).

Claimed files:
`05-knowledge/hypotheses/HYP-3163-lrc14-covariance-laplacian-associator-ear-angles.md`
and
`07-reflections/lrc14-covariance-laplacian-associator-ear-angles-codex-20260628.md`.

Goal: after incoming HYP-3161 exhaustively verifies `consec_8` as covariance
max over `3432` bounded clusters, split the proof mechanism into two
measurable tasks.  Even task: promote `Sigma kappa_2` to a covariance-kernel
theorem and test PSD/Laplacian/Monge/conditionally-positive certificates for
the consecutive row.  Odd task: promote the Worpitzky `-9S3` residue to an
associator or third-cumulant cocycle, then control it by odd-ear recursion,
K3 minority-edge lifts, n=4 canary/filler sources, reflection-fold
resurrection, or named finite debt.  HYP-3154's Joukowski/De Moivre bridge
adds the root-curve real-rootedness-defect sidecar.  The now-executed
HYP-3153 scout supplies exact Lee-Yang/Worpitzky/quartic packet columns to
join to these new covariance and associator columns.  HYP-3162's
cyclotomic-ideal / first-cubic-apex packet adds the calibration columns
`cyclotomic_cubic_defect`, `de_moivre_angle_slack`, and
`FM_AFM_bridge_status`.  HYP-3200 adds the primitive-normal-form audit:
consecutive is rank `0/3431` for `Sigma kappa_2`, the all-bank exception is
the nonprimitive dilation twin, exact `1/7` is false, and `kappa_4` is a
stabilizer sidecar.  HYP-3201 adds the quotient-legality audit: row entropy is
not the objective, but `H(target|quotient)` should be zero or produce a named
typed sidecar before any compression is used.  KPS-S31al adds the current best
spectral route: `q_t` as a Toeplitz trigonometric moment matrix, with
consecutive maximizing `lambda_min(T)` over all `3432` bounded k=8 rows, plus
a random-current warning for the ferromagnetic route.

Next scout should emit `covariance_kernel_distance_profile`,
`PSD_dual_slack_vector`, `monge_four_point_defect`,
`conditional_negative_type_status`, `associator_triple_cocycle`,
`odd_ear_surplus`, `worpitzky_boundary_mode`, and
`even_covariance_odd_associator_exchange`, while preserving the HYP-3162
cyclotomic angle-defect columns and the HYP-3200 primitive/dilation labels,
plus HYP-3201 `law_defect_entropy_bits` and
`target_function_factor_through_status`, and KPS-S31al
`toeplitz_lambda_min_margin` / `random_current_coupling_order_status`, over
the bounded k=8 banks from HYP-3138/HYP-3139/HYP-3142/HYP-3160.

## Lead codex-2026-06-28: Function-compression resolvent-degree wall

**Status:** EVIDENCE / executable factor-through scout; not proof
(HYP-3150/T1215).
**Status:** SYNTHESIS / exact finite scout; not proof (HYP-3150/T1215).

Claimed files:
`05-knowledge/hypotheses/HYP-3150-lrc14-function-compression-resolvent-degree-wall.md`,
`04-computation/lrc14_function_compression_resolvent_degree_codex_s277.py`,
`05-knowledge/results/lrc14_function_compression_resolvent_degree_codex_s277.out`,
and `07-reflections/function-compression-resolvent-degree-codex-s277.md`.

Readout: HYP-3150 integrates HYP-3143..HYP-3149 as a
function-compression legality test.  A quotient is legal only if the target
function is fiber-constant or a sidecar reconstructs the lost coordinate.
The exact scout verifies the prompt's pair-function split (`a+b,a*b`
unordered-safe; `a^b,b^a` ordered), the K3 quotient matrix `[[2,1],[3,0]]`,
the coin mix/same analogue, Worpitzky row `(1,4,1)`, and the n=4 monotone
compression `x=a OR c`, `y=b OR c`.

Readout: stored scout verifies the factor-through split exactly.  The largest
proof-carrier SCC couples the k8 even fold, K4 OR compression, K3 kernel,
fiber-PGF curve sidecar, and canary/deletion sidecar, so these are a packet
rather than a linear ladder.

Working hypothesis: the LRC14 hard core remains below the generic quintic
wall because every current compression has effective degree at most `4`, with
the deepest k=8 node reducing to a quadratic in `u^2` after sidecars are
accounted for.  Guardrail: do not use the Abel-Ruffini phrase as proof;
verify exact factor-through maps and named destroyed coordinates.
## Lead codex-2026-06-28-S278: Worpitzky function-compression resolvent bridge
Next: attach `function_payload_type`, `ordered_pair_sidecar`,
`state_level_pgf_split`, `compression_map_word`, `canary_filler_status`,
`resolvent_degree_ceiling`, and `quotient_legality_status` to HYP-3140
fiber-PGF rows, HYP-3141 edge-witness rows, and HYP-3142 k=8 sidecar rows.
Only then allow the degree-4/quartic bounded-core guardrail; otherwise route
the row to named Abel-Ruffini degree-5 debt.

**Status:** SYNTHESIS / exact finite scout; not proof (HYP-3151/T1216).

Claimed files:
`05-knowledge/hypotheses/HYP-3151-worpitzky-function-compression-resolvent.md`,
`04-computation/lrc14_worpitzky_function_compression_resolvent_codex_s278.py`,
`05-knowledge/results/lrc14_worpitzky_function_compression_resolvent_codex_s278.out`,
and `07-reflections/worpitzky-function-compression-resolvent-codex-s278.md`.

Readout: this executes HYP-3150's reserved factor-through wall.  The
proof-facing unit is a target function together with the
compression under which that function is being evaluated.  The exact scout
links HYP-3147's n=3 `C/T` edge kernel, HYP-3149/HYP-3148/HYP-3146's n=4
canary/filler quotient lanes, and HYP-3132/HYP-3142's k=8 biquadratic hard
node.  `a+b` and `a*b` are symmetric shadows; `a^b,b^a` are ordered channels.
The n=4 table compression is nonlinear, `x=a OR c`, `y=b OR c`, with no affine
class-preserving substitute.  The k=8 centered resolvent is
`u^4-5u^2+4`, so the bounded-core dual stays at degree `4`, below the generic
quintic obstruction.

Next: add `target_function_id`, `function_swap_parity`,
`compression_fiber_function_constancy`, `ordered_sidecar_status`,
`canary_or_restoration_sidecar`, and `resolvent_degree` to HYP-3141 edge rows,
HYP-3140 fiber-PGF rows, and HYP-3142 k=8 moment packets.  Any attempted
quotient should first answer whether the target function is fiber-constant,
symmetric-safe, ordered-sidecar reconstructed, canary-restored, or routed to
named debt.

Namespace: HYP-3151 / LTI-277 / LTT-175 / T1216 / OPEN-Q-108.

## Lead codex-2026-06-28: Lee-Yang/Worpitzky/quartic compression packet

**Status:** SYNTHESIS / exact finite packet scout; not proof
(HYP-3153/T1218).

Claimed files:
`05-knowledge/hypotheses/HYP-3153-lrc14-lee-yang-worpitzky-quartic-compression-packet.md`,
`04-computation/lrc14_lee_yang_worpitzky_quartic_packet_codex_20260628.py`,
`05-knowledge/results/lrc14_lee_yang_worpitzky_quartic_packet_codex_20260628.out`,
and `07-reflections/lrc14-lee-yang-worpitzky-quartic-packet-codex-20260628.md`.

Readout: HYP-3153 fuses HYP-3151's quotient-legality rule with HYP-3152's
Lee-Yang radius web and HYP-3150's parity split.  Exact checks verify
`q0=q6*R^6` for consec `k=8..13`, root-radius spreads `1.1427..1.3629`,
pair-mass dips `1081/76440`, `1/4004`, then zero, and
`L_y=q0+q6+q3/10 <= cap` for `k=8,9,10`.  The k=8 identity
`10q0+q3+10q6 = 10S0-10S1+10S2-9S3+6S4` separates the even biquadratic fold
from the larger odd Worpitzky/ear sidecar (`|odd/even|=3.149364`).

Next: bound the off-circle dip/lambda as even biquadratic contribution plus
odd Worpitzky/ear contribution, then attach the packet to a live HYP-3141 or
HYP-3142 row without discarding the ordered sidecar.

Namespace: HYP-3153 / LTI-279 / LTT-177 / T1218 / OPEN-Q-108.

## Lead codex-2026-06-27-S274: Worpitzky pair-function three-edge quotient

**Status:** EVIDENCE / exact K3 quotient scout; not proof
(HYP-3144/T1209).

Claimed files:
`05-knowledge/hypotheses/HYP-3144-lrc14-worpitzky-pair-function-three-edge-quotient.md`,
`04-computation/lrc14_worpitzky_pair_function_three_edge_quotient_codex_20260627.py`,
and `05-knowledge/results/lrc14_worpitzky_pair_function_three_edge_quotient_codex_20260627.out`.

Readout: the stored scout verifies the two score classes of labelled K3
tournaments, `T=(0,1,2)` and `C=(1,1,1)`, under single-edge flips.  The exact
quotient matrix is `[[2,1],[3,0]]`, with normalized eigenvalues `1,-1/3` and
stationary split `{T:3/4,C:1/4}`.  The same kernel appears from three coin
flips after quotienting to "two-to-one mix" and "all same."

Pair-function reading: `a+b` and `a*b` survive unordered-pair quotient;
`a^b` and `b^a` require an ordered sidecar.  LRC14 transfer: before
scalarizing HYP-3140 fiber-PGF data, HYP-3141 tip/tail witnesses, HYP-3139
reflection pages, or HYP-3143 exact-order packet bases, ask whether the target
predicate is symmetric like sum/product or order-sensitive like exponentiation.

Readout: stored scout verifies K3 class sizes `T=6,C=2`, quotient flip matrix
`[[2,1],[3,0]]`, normalized eigenvalues `1,-1/3`, stationary split
`{T:3/4,C:1/4}`, and the identical three-coin mix/same kernel.  The edge role
ledger has `12` adjacent transitive self flips, `6` cycle-break flips, and `6`
long source-sink exits.  The key Worpitzky/PGF warning is that both score
classes aggregate to `F=(1,4,1)` while state-level PGF curves differ.  Thus a
single value or even a class aggregate can erase the ordered function payload.

HYP-3147 sharpens the same lane by naming the transitive fiber's minority-edge
gate.  S71 then gives the frontier connection: score-to-isomorphism
compression is safe through n=4 and first fails at n=5, while the k=8 cap dip
splits into a larger odd Worpitzky `-9S3` face and a smaller even biquadratic
`+6S4` face.  HYP-3150 turns this into a completed factor-through audit:
before using HYP-3140/HYP-3141/HYP-3139/HYP-3143/HYP-3145/HYP-3149 quotients,
prove the observable factors through the quotient or retain the sidecar.

Namespace: HYP-3144 / LTI-270 / LTT-168 / T1209 / OPEN-Q-108.
## Lead codex-2026-06-27-S274: n=4 filler/canary shift-package quotient

**Status:** SYNTHESIS / exact n=4 scout; not proof (HYP-3146/T1211).

Claimed files:
`05-knowledge/hypotheses/HYP-3146-n4-shift-package-erdos870-filler-canary.md`,
`04-computation/tournament_n4_shift_package_erdos870_codex_s274.py`,
`05-knowledge/results/tournament_n4_shift_package_erdos870_codex_s274.out`,
and `07-reflections/n4-shift-package-erdos870-filler-canary-codex-s274.md`.

Readout: fixed-path cover has `S` fiber size `5`, while a finite scaffold
makes `T,+,-,S` a two-bit shift package.  This is the canary/scaffold
companion to HYP-3143's exact-order subbasis audit and HYP-3145's filler-core
interface.  Next: attach
filler/canary fields to HYP-3141/HYP-3142 ledgers and test if quotienting k=8
packet fibers requires deletion-stable canaries or finite-filler congruence.
## Lead codex-2026-06-27-S275: Erdos-870 live-core deletability audit

**Status:** SYNTHESIS / exact n=4 table scout; not proof
(HYP-3148/T1213).

Claimed files:
`05-knowledge/hypotheses/HYP-3148-erdos870-live-core-filler-canary.md`,
`04-computation/lrc14_erdos870_live_core_canary_scout_codex_s275.py`,
`05-knowledge/results/lrc14_erdos870_live_core_canary_scout_codex_s275.out`,
and `07-reflections/erdos870-live-core-filler-canary-codex-s275.md`.

Readout: this continues HYP-3143/S276's exact-order subbasis packet,
HYP-3144/S274's pair-function/order-sidecar quotient guardrail, and HYP-3145's
filler-core interface, and now HYP-3146/S274's shift-package canary/scaffold
policy plus S277/HYP-3147's n=3 edge-flip/Worpitzky/function kernel, by adding
the deletable-coordinate live-core/filler/canary audit.  The user's two n=4
tournament tables are both
exact.  The fixed
Hamiltonian-path tiling cube with live skips `a,b,c` has class distribution
`T:+:-:S = 1:1:1:5`; the long diagonal `c` is class-cover-deletable because
`{a,b}` already reaches all four classes.  Freezing `c` as deterministic
filler gives the two-bit anchor with fixed partial score `(0,1,1,2)`, live
opposite skips `x,y`, uniform class distribution `1:1:1:1`, and load-bearing
live core `{x,y}`.  The scout finds `24` labelled two-bit anchors, all with
disjoint live pairs, split evenly across the three perfect matchings.

LRC14 transfer: import the Erdos-870 no-minimal-subbasis warning as a quotient
audit.  Many witnesses, score classes, or representations do not prove minimal
support.  Add `live_core_bits`, `filler_bits`, `canary_bits`,
`deletable_coordinates`, `class_distribution`, `minimal_cover_subbasis`,
`edge_bounded_core_floor_exit`, and `terminal_exit_or_named_debt` to
edge-witness, fiber-PGF, A000568, and k=8 sidecar rows before scalarizing.

Next: run the live-core/filler/canary audit on HYP-3141 edge rows and on the
HYP-3140/HYP-3142 coefficient/moment packets; a coordinate that is deletable in
the local quotient should not be treated as part of the terminal proof payload.

Namespace: HYP-3148 / LTI-274 / LTT-172 / T1213 / OPEN-Q-108.  Predecessor:
HYP-3147 / LTI-273 / LTT-171 / T1212, HYP-3146 / LTI-272 / LTT-170 / T1211,
HYP-3143 / LTI-269 / LTT-167 / T1208, HYP-3144 / LTI-270 / LTT-168 / T1209,
and HYP-3145 / LTI-271 / LTT-169 / T1210.
## Lead codex-2026-06-28-erdos870-tournament4: tournament-4 canary/filler quotient

**Status:** EVIDENCE / exact n=4 table scout and proof-interface transfer;
not a proof (HYP-3149/T1214).

Relation to incoming S276: HYP-3143 records the exact-order subbasis version
of the same n=4 data.  HYP-3149 is the complementary canary/filler refinement:
it names the coordinate whose deletion turns the two-bit source into a
collision-prone quotient.  Relation to incoming S274: HYP-3144 records the
three-edge pair-function scalarization guardrail, HYP-3145 records the broader
Erdos-870 filler-core interface, HYP-3146 records the cover/scaffold
shift-package policy, HYP-3147 records the n=3 edge-flip/Worpitzky kernel,
and HYP-3148 records the live-core deletability audit; HYP-3149 names the
concrete fixed-path canary sidecar inside the four-vertex packet.

Claimed files:
`05-knowledge/hypotheses/HYP-3149-lrc14-erdos870-tournament4-canary-filler.md`,
`04-computation/lrc14_erdos870_tournament4_canary_filler_codex_20260628.py`,
`05-knowledge/results/lrc14_erdos870_tournament4_canary_filler_codex_20260628.out`,
and `07-reflections/lrc14-erdos870-tournament4-canary-filler-codex-20260628.md`.

Readout: the fixed-Hamiltonian-path tiling table is a three-bit quotient map,
not a group law on isomorphism classes.  With free arcs
`a=(0,2), b=(1,3), c=(0,3)`, the fibers are
`T={E}`, `+={a}`, `-={b}`, and `S={c,ab,ac,bc,abc}`.  Fixing `c` unflipped
gives the second prompt table on `x=a,y=b` and partial score sequence
`(0,1,1,2)`; the `c=1` slice collapses every completion to `S`.

Erdos-870 transfer: treat the two-coordinate `x,y` table like an order-two
source, and treat `c` as a deterministic filler/canary coordinate whose
deletion/restoration status must be audited before the quotient enters an
LRC14 edge packet.

Next: attach `tournament4_canary_filler_certificate`, `c_canary_status`,
`xy_completion_table`, `S_bulk_fiber_words`, `deletion_restoration_sidecar`,
and `edge_tip_tail_exit_or_named_debt` to HYP-3141 edge-witness rows and
HYP-3140/HYP-3138/HYP-3139 local quotient packets.

Namespace: HYP-3149 / LTI-275 / LTT-173 / T1214 / OPEN-Q-108.

## Lead codex-2026-06-27-S273: LRC14 generating-function payload atlas

**Status:** EVIDENCE / executable payload atlas and tournament scout; not a
proof (HYP-3137/T1202).

Claimed files:
`05-knowledge/hypotheses/HYP-3137-lrc14-generating-function-payload-atlas.md`,
`04-computation/lrc14_generating_function_payload_atlas_codex_s273.py`,
`05-knowledge/results/lrc14_generating_function_payload_atlas_codex_s273.out`.

Goal: compare generating functions by retained proof payload rather than raw
scalar value.  Candidate carriers: miss-count PGF, signed SPEC/resonance
Fourier series, A000568 count/cycle-index quotient, OCF independence
polynomial/Walsh EGF, De Moivre/resolvent elementary-symmetric polynomial,
Irving-Omar walk determinant GF, and hard-core/polymer partition function.
The scout should declare which LRC predicate survives, what coefficient/root/
quotient/tail coordinate is retained, what is destroyed by scalar evaluation,
and what terminal exit or named debt remains.

Readout: the stored scout finds the smallest full-cover packet
`signed_SPEC_resonance_series + A000568_cycle_index_quotient +
miss_count_PGF_root_locus`, covering all twelve tracked obligations from
`Q_apex_floor` through `destroyed_coordinate_guard`.  Tournament fingerprints:
12 carriers, 66 edges, score histogram `0..11`, 0 directed 3-cycles, singleton
SCCs, Hamiltonian path count 1, and every typed payload beats
`raw_scalar_evaluation`.

Next fields to attach to finite constant-chase rows: `SPEC_support_sieve`,
`edge_recursion_depth_PGF`, `global_consistency_class`,
`PGF_root_trajectory_derivative`, `middle_layer_vector`,
`k8_reflection_fold_adjoint`, `repair_cover_H(q)`,
`Bravais_resonance_cell`, `Savitch_packet_depth`, and
`terminal_exit_or_named_debt`.  Post-rebase HYP-3138 supplies the concrete
k=8 fold-adjoint lookup for the resolvent carrier.

Namespace: HYP-3137 / LTI-263 / LTT-161 / T1202 / OPEN-Q-108.  The atlas
cross-references HYP-3136 as the integrated multi-far floor closure it is
meant to feed.
## Lead codex-2026-06-27-k8-reflection-fold: k=8 reflection-fold coordinate resurrection

**Status:** EVIDENCE / executable bounded-bank quotient audit; not proof
(HYP-3138/T1203).

Claimed files:
`05-knowledge/hypotheses/HYP-3138-lrc14-k8-reflection-fold-resurrection.md`,
`04-computation/lrc14_k8_reflection_fold_resurrection_codex_20260627.py`,
`05-knowledge/results/lrc14_k8_reflection_fold_resurrection_codex_20260627.out`,
and `07-reflections/lrc14-k8-reflection-fold-resurrection-codex-20260627.md`.

Readout: HYP-3132's De Moivre hard-row fold gives the even biquadratic
coordinate, but HYP-3138 checks the destroyed coordinates directly.  For k=8,
folding the miss distribution by `t <-> 6-t`,
`even_fold=(q0+q6,q1+q5,q2+q4,q3)`, is injective on the tested primitive
bounded banks span<=14,15,16 (`3431`, `6434`, `11432` rows and `0` collision
fibers).  The top row remains `(0,1,2,3,4,5,6,7)` with
`L_yK8=2633/735`, margin `10*cap8-L_yK8=683/2940`, and nonzero odd leakage
`(451/1470,142/735,131/1470)`.

Next: prove a finite k=8 fold-adjoint lemma.  The even fold may carry the
gK8/phi4 dip bound, but endpoint `Phi/P`, observer gluing, and finite-address
exits need either the odd-coordinate resurrection table or named debt.

Namespace: HYP-3138 / LTI-264 / LTT-162 / T1203 / OPEN-Q-108.
## Lead codex-2026-06-27-S273: LRC14 fiber-PGF Rprime certificate

**Status:** EVIDENCE / exact generating-function scout; not proof
(HYP-3140/T1205).

Claimed files:
`05-knowledge/hypotheses/HYP-3140-lrc14-fiber-pgf-rprime-certificate.md`,
`04-computation/lrc14_fiber_pgf_certificate_codex_s273.py`,
`05-knowledge/results/lrc14_fiber_pgf_certificate_codex_s273.out`,
and `07-reflections/lrc14-fiber-pgf-rprime-certificate-codex-s273.md`.

Readout: this refines HYP-3136's integrated multi-far factorization at its
remaining `Rprime` factor and supplies a concrete coefficient carrier for the
completed HYP-3137 generating-function payload atlas.  For `S=R union 14Q`,
pass to `u=14t` and count the R-safe lifts
`N_R(u)=#{a: (u+a)/14 is R-safe}`.  Then
`Rprime=E[N_R | Q-lonely]/E[N_R]` exactly.  The scout computes exact
sheet-count PGFs `F_R(y)` and Q-masked PGFs `F_R,Q(y)`.  The HYP-3129 worst
targeted row `R={1,...,12}`, `Q={1,2}` becomes a two-coefficient defect:
`F_R=7243/13860*y^0+6617/13860*y^1`,
`F_R,Q=7243/13860*y^0+521/1980*y^1`, and
`Rprime=51058/72787=0.701471...`.

LRC14 placement: HYP-3129's signed SPEC certificate is the Fourier transform
of this finite fiber-PGF moment inequality.  Pointwise positive-sheet
positivity is false in low rows, so the live target is the conditional
first-moment bound
`F_R,Q'(1)/F_R,Q(1) >= c*F_R'(1)/F_R(1)` over legal residual packets.  This
integrates HYP-3137/HYP-3136, Lee-Yang PGFs, Delsarte/MacWilliams transforms,
q-Pochhammer tail warnings, Moser/fibbinary partial-cube sidecars, and HYP-3134 quotient
legality.

Next: enumerate the legal post-HYP-3131 residual packet family, attach
`fiber_pgf_word` and `Q_masked_fiber_pgf_word` to HYP-3125 edge-floor packets,
and translate the HYP-3129 low-frequency SPEC certificate into coefficient
moments of `F_R,Q`.

Namespace: HYP-3139 / LTI-265 / LTT-163 / T1204 / OPEN-Q-108.
## Lead codex-2026-06-27-S273: k=8 resolvent sidecar certificate

**Status:** EVIDENCE / exact bounded-bank scout through B=14; not proof
(HYP-3142/T1207).

Claimed files:
`05-knowledge/hypotheses/HYP-3142-k8-resolvent-sidecar-certificate.md`,
`04-computation/lrc14_k8_resolvent_sidecar_scout_codex_s273.py`,
`05-knowledge/results/lrc14_k8_resolvent_sidecar_scout_codex_s273.out`,
`05-knowledge/results/lrc14_k8_resolvent_sidecar_scout_codex_s273_B14.out`,
and `07-reflections/k8-resolvent-sidecar-certificate-codex-s273.md`.

Readout: the current one-node frontier can be attacked by an exact 4-moment
sidecar.  For primitive k=8 bounded rows, `B=13` scans `1716` rows and `B=14`
scans `3431`; both have `U4_over_cap_count=0`.  The worst row is exactly
`consec_8`, with `U4=2633/7350`, `cap_8-U4=683/29400`, nearest PGF root
`1.488584`, no real roots, Bravais residue counts `(2,1,1,1,1,1,1)`, peak
`1/8`, entropy `1`, mirror defect `0`, and `kappa4=-0.787150`.

Connection to niche past work: HYP-3140 supplies the fiber-PGF `Rprime`
moment, HYP-3137 supplies the generating-function payload lens, HYP-3138
supplies reflection-fold coordinate resurrection, HYP-3136/HYP-3135 supply
the surrounding multi-far floor and resolvent-packet theorem shape, HYP-3132
supplies the biquadratic De Moivre fold, HYP-3113 supplies the
Bravais-flat/Lee-Yang map, HYP-3118 supplies Savitch repair-depth language,
and HYP-3134 supplies the controlled-forgetting guard before dropping
edge-child payload.

Next: prove global moment-majorization `U4(E) <= U4(consec_8)` for all
primitive k=8 bounded-core shapes.  Non-flat residue spectra should give
strict slack; any exception should descend through a Hensel/CRT `2x7` sidecar
or become named finite resonance debt.

Namespace: HYP-3142 / LTI-268 / LTT-166 / T1207 / OPEN-Q-108.
## Lead codex-2026-06-27-Erdos870: n=4 filler-core quotient interface

**Status:** EVIDENCE / executable finite interface scout; not proof
(HYP-3145/T1210).

Claimed files:
`05-knowledge/hypotheses/HYP-3145-erdos870-n4-filler-quotient-models.md`,
`04-computation/lrc14_erdos870_n4_filler_models_codex_20260627.py`,
`05-knowledge/results/lrc14_erdos870_n4_filler_models_codex_20260627.out`,
and `07-reflections/erdos870-n4-filler-quotient-models-codex-20260627.md`.

Readout: Erdos-870 contributes an interface pattern rather than a direct LRC
theorem: small core plus deterministic finite fillers, with deletion and
nonminimality checked at the boundary.  The n=4 fixed-path tiling table is a
representative atlas but the `S` fiber has five representatives, so the class
quotient is not congruent.  The partial-score model fixes four arcs with
profile `(0,1,1,2)` and leaves two core arcs `x=(0,1)`, `y=(2,3)`; the closed
four-state table is congruent and realizes `T,+,-,S`.

Next: instantiate `filler_core_interface` on one HYP-3125/HYP-3129 covering
row before using HYP-3136.  The row should identify deterministic fillers,
the small signed core, quotient-congruence status, nonminimal fiber alarms,
and the formal interface target.

Namespace: HYP-3145 / LTI-271 / LTT-169 / T1210 / OPEN-Q-108.

## Lead codex-2026-06-27-A000568: LRC14 A000568 edge-witness sandwich

**Status:** EVIDENCE / executable quotient scout; not proof
(HYP-3133/T1200).

Claimed files:
`05-knowledge/hypotheses/HYP-3133-lrc14-a000568-edge-sandwich.md`,
`04-computation/lrc14_a000568_edge_sandwich_codex_20260627.py`,
`05-knowledge/results/lrc14_a000568_edge_sandwich_codex_20260627.out`,
and `07-reflections/lrc14-a000568-edge-sandwich-codex-20260627.md`.

Readout: the user's `12` and `56` observation is a shifted quotient sandwich.
HYP-3124 edge-local signatures on `m` vertices compare to A000568 unlabeled
tournaments on `m+1` vertices: `m=4: 10 < 12 < 16`,
`m=5: 20 < 56 < 80`, and `m=6: 35 < 456 < 632`.  Sector word is the
equinumerosity layer, A000568 is the unrooted one-extra-vertex
equidistribution/free-extension shadow, and paired tail/tip child deck is the
equidecomposability layer.

Next: add `a000568_extension_shadow` to HYP-3125 edge-floor packets and use it
to stratify HYP-3129's finite low-frequency SPEC constant chase.  HYP-3128 is
the guardrail: this is a controlled-forgetting diagnostic, not a naive Asano
zero-free proof shortcut.

Namespace: HYP-3133 / LTI-261 / LTT-159 / T1200 / OPEN-Q-108.
## Lead codex-2026-06-27-S272: A000568 edge-envelope global-consistency quotient

**Status:** EVIDENCE / executable edge-envelope scout; not proof
(HYP-3134/T1201).

Claimed files:
`05-knowledge/hypotheses/HYP-3134-a000568-edge-envelope-global-consistency-quotient.md`,
`04-computation/a000568_edge_envelope_lrc14_codex_s272.py`,
`05-knowledge/results/a000568_edge_envelope_lrc14_codex_s272.out`, and
`07-reflections/a000568-edge-envelope-global-consistency-codex-s272.md`.

Readout: the HYP-3124 lower envelope is raw four-sector size decks
`1,4,10,20,35`; the upper envelope is sector plus paired endpoint-deletion
children `1,4,16,80,632`; A000568 one vertex later sits between:
`10 < U(5)=12 < 16`, `20 < U(6)=56 < 80`, and S272 verifies
`35 < U(7)=456 < 632`.  Interpretation: A000568 is the global-consistency
quotient inside the local edge-witness envelope.

Post-fetch connection: HYP-3133's direct edge-sandwich scout, HYP-3132's
integrated multi-far floor, HYP-3130's Gaussian/Beurling-Selberg minorant, and
HYP-3131's far-zero outward motion are compatible; HYP-3134 is the quotient
discipline for when those analytic certificates may forget local paired child
payload.

Next: add `envelope_position`, `global_consistency_class`,
`edge_child_gluing_status`, `resonance_lattice_class`, `SPEC_bound_status`,
and `terminal_exit_or_named_debt` to HYP-3125/HYP-3129 edge-floor packet rows
before quotienting away paired tail/tip child data.

Namespace: HYP-3134 / LTI-262 / LTT-160 / T1201 / OPEN-Q-108.
## Lead codex-2026-06-27-S272: LRC14 resolvent-packet middle-layer synthesis

**Status:** SYNTHESIS / proof-target refinement; not proof (HYP-3135).

Claimed files:
`05-knowledge/hypotheses/HYP-3135-lrc14-resolvent-packet-middle-layer-synthesis.md`,
`04-computation/lrc14_resolvent_packet_synthesis_codex_s272.py`,
`05-knowledge/results/lrc14_resolvent_packet_synthesis_codex_s272.out`,
and `07-reflections/lrc14-resolvent-packet-synthesis-codex-s272.md`.

Readout: the user's De Moivre-style quintic gives a concrete controlled-
forgetting lesson.  The resolvent roots `2,-4,8,-16` have elementary symmetric
layers `e1=-10`, `e2=-120`, `e3=320`, `e4=1024`; hence `120` and `320` are
pair/triple branch payloads, not standalone roots.  The real branch expression
`fifthroot(2)-fifthroot(4)+fifthroot(8)-fifthroot(16)` checks against the
quintic.  Prior repo motifs with 120/320 mostly mark finite middle layers:
cycle-pair corrections, bounded 0/320 resonant-center checks, support-six
signed cancellation, cap slack, fixed-path/tournament scale, torsion moats,
and q^120 modular support horizons.

LRC14 placement: Q/apex block is closed by HYP-3130/HYP-3128; far additions
help by HYP-3131; HYP-3133/HYP-3134 add the A000568 edge-extension quotient
stratifier and global-consistency quotient discipline; the absolute minorant
envelope is ruled out; HYP-3129/HYP-3132 give the signed SPEC coupling
certificate and need closed-form constant chasing; and S70 sharpens the
bounded core to the k=8 biquadratic coefficient bound.  The remaining proof
should be framed as a bounded-core plus
middle-payload packet theorem retaining Q-floor constants, signed SPEC
low/tail, Lee-Yang radius, far-push status, edge tail/tip deletion sectors,
and finite-address or observer-gluing exits.

Next: turn the HYP-3129 per-row exact certificate into a uniform symbolic
constant; prove far-pushes-out for all far placements; make the S70 k=8
reflection-fold/biquadratic coefficient bound rigorous; formalize the
bounded-core `rho>1 => Rprime>=c` bridge; attach HYP-3124/HYP-3125 edge-witness
fields to real covering packets.

Namespace: HYP-3135 / OPEN-Q-108.

## Lead codex-2026-06-27-S271: LRC14 multi-far edge-witness Rprime floor

**Status:** EVIDENCE / executable edge-floor synthesis scout; not proof
(HYP-3125/T1199).

Claimed files:
`05-knowledge/hypotheses/HYP-3125-lrc14-multifar-edge-witness-rprime-floor.md`,
`04-computation/lrc14_multifar_edge_witness_floor_codex_s271.py`,
`05-knowledge/results/lrc14_multifar_edge_witness_floor_codex_s271.out`,
and `07-reflections/lrc14-multifar-edge-witness-rprime-floor-codex-s271.md`.

Readout: integrates incoming HYP-3124 with HYP-3121 by treating the open
`r=2..6` decorrelation floor as a two-ended directed proof edge
`R-safe packet -> Q-safe packet`.  In the lifted `u=14t` coordinate, `Rprime`
is the normalized diagonal sector mass.  The S271 audit keeps the known
Bonferroni-negative rows positive by quasi-independence (`0.513784`,
`0.925326`), shows tail deletion usually improves the floor, shows tip
deletion is the sharper multiple-of-14 recursion, and finds individual R/Q
edge ratios near `1`, pointing to packet-level distribution rather than a
single bad pair.

Next: add `edge_floor_packet`, deletion-child `Rprime` ratios,
`EH_level_distribution_proxy`, major-arc exceptions, Gaussian smoothing width,
finite-ruler desmoothing threshold, Asano/Lee-Yang contraction status, phi4
sign, normal-fan chamber, chiral guard word, and terminal debt to real
`r=2..6` covering rows.  Prove an LRC-specific level-of-distribution bound
over edge-sector residues, then desmooth; route failures to
Asano/Lee-Yang/phi4/Cech/H7 debt.

Post-rebase connection: incoming HYP-3127/S68 strengthens this lead by making
Asano contraction the candidate main engine.  Treat S271 as the packet schema
and diagnostic harness for HYP-3127's single-far zero-free, `SPEC`-bound, and
monotonicity obligations.

Namespace: HYP-3125 / LTI-260 / LTT-158 / T1199 / OPEN-Q-108.
## Lead codex-2026-06-27-S271: LRC14 edge-witness class-deck stress supplement

**Status:** EVIDENCE / executable class-deck stress scout; not proof
(HYP-3124/T1198/LTI-259/LTT-157).

Claimed files:
`04-computation/lrc14_tournament_edge_witness_recursion_codex_20260627.py`,
`05-knowledge/results/lrc14_tournament_edge_witness_recursion_codex_20260627.out`,
`05-knowledge/hypotheses/HYP-3124-lrc14-tournament-edge-witness-recursion.md`,
and
`07-reflections/lrc14-edge-witness-recursion-multifar-floor-codex-20260627.md`.

Readout: directed tournament edges are treated as two-ended witnesses with
outside four-sector decks, tail-deletion child, tip-deletion child, and
cross-sector orientation payload.  Class-deck audit: `n=3,4,5` all tested
modes separate all classes; at `n=6`, sector counts/internal decks separate
`55/56`, colliding only on converse pair `344/345`, while roleless children,
recursive children, and full edge witnesses separate `56/56`.  Edge-instance
audit: `43` nontrivial `n=6` sector-internal fibers split by tail/tip
deletion children and `16` recursive fibers split further by cross-sector
orientation.

HYP-3125/HYP-3127 now carry the measured multi-far `Rprime >= c` program.
This HYP-3124 supplement should be read as the finite witness-memory stress
test feeding that program: before an `R-safe packet -> Q-safe packet` is
compressed, retain tail/tip deletion children and cross-sector orientation or
name the lost coordinate.

Next: attach `edge_witness_certificate`, `tail_deletion_child_signature`,
`tip_deletion_child_signature`, `recursive_tail_child_edge_deck`,
`recursive_tip_child_edge_deck`, `edge_missing_input_vector`, and
`edge_repair_sidecar` to HYP-3115 one-swap/domain-wall edges and HYP-3098
observer-gluing rows.  Use HYP-3125/LTI-260 for the few-apex `Rprime` row
fields.

Namespace: HYP-3124 / LTI-259 / LTT-157 / T1198 / OPEN-Q-108.

## Lead codex-2026-06-27-S270: LRC14 chiral base-stalk guard and normal-fan Cech finite-ruler proof angles

**Status:** SYNTHESIS / executable proof-angle scout; not proof (HYP-3123/T1197).

Claimed files:
`05-knowledge/hypotheses/HYP-3123-lrc14-chiral-stalk-cech-proof-angles.md`,
`04-computation/lrc14_chiral_cech_proof_angles_codex_20260627.py`,
`05-knowledge/results/lrc14_chiral_cech_proof_angles_codex_20260627.out`,
and `07-reflections/lrc14-chiral-stalk-cech-proof-angles-codex-20260627.md`.

Readout: selected two routes different from the broad HYP-3118 coordinate
resurrection calculus.  The chiral base-stalk guard treats mirror/converse,
rootless, directed-edge-sector, and perspective-pair quotients as legal only
after `chiral_guard_word`, `cross_sector_orientation_word`,
`endpoint_owner_cocycle`, and `state_lift_sign` are retained or both fibers
share terminal exit.  The normal-fan Cech finite-ruler route upgrades raw
component count into `normal_fan_chamber_id`, `closed_arc_cech_beta`,
`barcode_persistence_word`, `owner_current_word`, and
`finite_ruler_denominator_threshold`.

Top scout bridges: Cech route to HYP-3101 component bound (`59`),
first-obstruction syndrome to HYP-3102 (`58`), Lean frontier bus to HYP-3098
observer rows (`54`), Cech route to THM-573 residual (`53`), and chiral guard
to HYP-3112 ear payload (`53`) / HYP-3098 observer rows (`52`).  Tournament
Analysis has two directed 3-cycles involving Cech, first-obstruction, observer
payload, and endpoint Phi, so the next ledger should couple those carriers
rather than choose a single scalar winner.

Incoming work connection: HYP-3121 makes the covering case one
lift-and-decorrelate engine, so the normal-fan/Cech fields should be tested as
the normalized component payload for the `r=2..6` multi-far floor.  HYP-3122's
phi4 stabilizer and the HYP-3124 edge-witness recursion scout are natural
stress tests for the chiral guard through cumulant/ear parity and tail/tip
sectors.
phi4 stabilizer and HYP-3124's edge-witness recursion are natural stress
tests for the chiral guard through cumulant/ear parity and tail/tip sectors.

Next: run a joined packet ledger over HYP-2963/HYP-3098/HYP-3107/HYP-3112
with chiral guard fields, normal-fan/Cech finite-ruler fields, first-obstruction
syndrome fields, endpoint Phi/P gates, and `terminal_exit_or_named_debt`.

Namespace: HYP-3123 / LTI-258 / LTT-156 / T1197 / OPEN-Q-108.
## Lead codex-2026-06-27-S268: LRC14 tournament edge-witness recursion

**Status:** EVIDENCE / executable edge-witness recursion scout; not proof
(HYP-3124/T1198).

Claimed files:
`05-knowledge/hypotheses/HYP-3124-lrc14-tournament-edge-witness-recursion.md`,
`04-computation/lrc14_tournament_edge_witness_recursion_codex_s268.py`,
`05-knowledge/results/lrc14_tournament_edge_witness_recursion_codex_s268.out`,
and `07-reflections/lrc14-tournament-edge-witness-recursion-codex-s268.md`.

Readout: the edge carrier is two-ended.  A directed edge `tail -> tip` should
carry endpoint role, outside four-sector deck, tail-deletion and tip-deletion
child signatures, recursive child edge decks, observer-gluing payload,
missing-input vector, coordinate-resurrection sidecar or named debt, and
terminal exit.  The exact census through labelled tournaments `n<=5` shows
that sector counts alone are `1,4,10,20`, while sector plus paired endpoint
children are `1,4,16,80`; at `n=5`, all `20` sector groups split by child
pair.  Tournament Analysis over edge-witness reframes has one directed
3-cycle, three Hamiltonian paths, and selected path led by
`coordinate_resurrection_edge_sheaf -> edge_witness_two_ended_packet ->
cross_sector_orientation_word`.

Next: attach `edge_witness_certificate`, `edge_tail_tip_sector_word`,
`tail_deletion_child_signature`, `tip_deletion_child_signature`,
`recursive_tail_child_edge_deck`, `recursive_tip_child_edge_deck`,
`edge_missing_input_vector`, `edge_repair_sidecar`, and `edge_terminal_exit`
to HYP-3115 one-swap/domain-wall edges and HYP-3098 observer-gluing rows.
## Lead codex-2026-06-27: LRC14 recursive tournament-edge witness packet

**Status:** EVIDENCE / exact small-tournament information audit and packet
schema; not proof (HYP-3124/T1198).

Claimed files:
`05-knowledge/hypotheses/HYP-3124-lrc14-tournament-edge-witness-recursion.md`,
`04-computation/lrc14_tournament_edge_witness_recursion_codex_20260627.py`,
`05-knowledge/results/lrc14_tournament_edge_witness_recursion_codex_20260627.out`.

Readout: the scout treats a directed edge `tail -> tip` as a recursive proof
witness retaining the outside four-sector observer deck, the tail child after
deleting the tip, the tip child after deleting the tail, and both child edge
decks.  Through all unlabeled tournaments on `n<=6`, score sequences separate
only `22/56` classes at `n=6`, the depth-0 edge-sector deck separates `55/56`,
and depth-1/depth-2 recursive edge-witness decks separate `56/56` with no
collisions.  The proof-lens tournament is transitive and led by
`recursive_edge_witness_packet`, then `edge_coordinate_resurrection_guard`,
then `tail_tip_child_pair`, then `decorrelation_edge_floor`, then
`h7_state_lift_edge_boundary`.

Integration: this row now treats HYP-3122/S67 as the phi4 quartic wall stress
and HYP-3123/S270 as the chiral/Cech orientation guard.  HYP-3124 is the
edge-local packet that tests whether those signals still point back to
tail/tip children, observer sectors, and named destroyed-coordinate repairs.
Incoming HYP-3128/HYP-3129 add the next correction: Asano/Lee-Yang certifies
the apex/tip side but exposes the overcrowded tail obstruction, so the
positive multi-far floor must be recorded as a retained-edge SPEC
resonance-lattice certificate, not as a collapsed zero-free tail contraction.
HYP-3130/HYP-3131 add the companion orientation: the far/apex tip side is
stabilizing (minorant apex floor plus Lee-Yang roots pushed outward), and the
bounded-core/tail side is binding.

Next refinement: formulate the edge-witness descent lemma explicitly.  In an
`R-safe -> Q-lonely` packet, tip children must close or descend by
minorant/zero-free/far-push data, tail children must close or descend by
bounded-core SPEC data, and cross-sector orientation must carry the signed
coupling unless it is replaced by named observer-gluing, coordinate-
resurrection, finite-ruler, phi4, or H7/F7 debt.

Next: attach `edge_witness_recursion_id`, `tail_child_packet`,
`tip_child_packet`, `four_sector_observer_deck`, `child_deck_asymmetry`,
`coordinate_resurrection_status`, `decorrelation_floor_status`,
`asano_obstruction_status`, `spec_resonance_floor_status`,
`minorant_apex_floor_status`, `bounded_core_binding_status`,
`state_lift_boundary_status`, `phi4_edge_wall_status`, and
`terminal_exit_or_named_debt` to HYP-2963/HYP-3098/HYP-3107 packet rows.  Use
the edge witness as a legal sidecar for observer gluing, coordinate
resurrection, HYP-3121 lift-and-decorrelate cuts, and zero-mass H7/F7
state-lift boundaries; keep raw H or raw sector counts as alarms only.

Namespace: HYP-3124 / LTI-259 / LTT-157 / T1198 / OPEN-Q-108.

## Lead codex-2026-06-27-S267: LRC14 coordinate-resurrection sheaf and adjoint repair calculus

**Status:** EVIDENCE / executable coordinate-repair scout and abstract synthesis;
not proof (HYP-3118/T1193).

Claimed files:
`05-knowledge/hypotheses/HYP-3118-lrc14-coordinate-resurrection-sheaf.md`,
`04-computation/lrc14_coordinate_resurrection_sheaf_codex_s267.py`,
`05-knowledge/results/lrc14_coordinate_resurrection_sheaf_codex_s267.out`.

Readout: after HYP-3116's missing-input ledger and HYP-3117's proof-circuit
past-work compiler, the S267 scout models quotient repair as a
coordinate-resurrection problem.  It uses a base stalk
`finite_address + observer_gluing + endpoint_owner + uniformity` plus one live
section.  Scalar-like shortcuts with no live section require size-`4` repair
covers and have `40` minimal choices; live-section shortcuts require size-`3`
covers.  Observer repair splits between chart certificates and
ear-decomposition grammars; uniformity repair splits between proof-circuit buses
and Savitch midpoint ladders.  The two comprehensive maps are:
quotient-to-destroyed-coordinate repair, and theorem-to-new-signal dictionary
for Savitch, Bravais, Lee-Yang/PGF roots, `phi4`, and directed/odd/nested ears.
Tournament vertices are reframes rather than runners or scalar signals, with a
transitive retention path led by `coordinate_resurrection_sheaf`.  Rebased over
the S266 proof-carrier augmentation, this also treats HYP-2112 `Phi`, HYP-2108
`P`, HYP-2109 `L/M/R`, HYP-3023 magnitude cocycle, HYP-3077 Horn closure, and
HYP-3082 protected branch status as candidate repair sections for exact-gap,
endpoint-activation, wall-crossing, route-purity, legality, and bridge-safety
coordinates.

Next: add `destroyed_coordinate_vector`, `coordinate_resurrection_cover`,
`repair_cover_rank`, `adjoint_section_status`, `observer_ear_certificate_type`,
`midpoint_certificate_depth_profile`, `pgf_zero_trajectory_signature`,
`bravais_shape_wall_signature`, `phi_gap_sum`, `P_max_activation`,
`LMR_terminal_state`, `magnitude_cocycle_route_purity`,
`horn_legality_status`, `protected_branch_status`, and
`terminal_exit_or_named_debt` to HYP-2963, HYP-3098, HYP-3107, and HYP-3112
packet rows.

Namespace: HYP-3118 / LTI-254 / LTT-152 / T1193 / OPEN-Q-108.
## Lead codex-2026-06-27: LRC14 niche past-work closure bridge

**Status:** EVIDENCE / synthesis scout; not proof (HYP-3120/T1195).
**Readout:** Added
`05-knowledge/hypotheses/HYP-3120-lrc14-niche-past-work-closure-bridge.md`,
`04-computation/lrc14_niche_past_work_closure_bridge_codex_20260627.py`,
`05-knowledge/results/lrc14_niche_past_work_closure_bridge_codex_20260627.out`,
and
`07-reflections/lrc14-niche-past-work-closure-bridge-codex-20260627.md`.
The scout turns older niche LRC threads into a current proof-frontier router:
finite-address `Phi` tuple, observer-cut payload orbit, Q27/Q31 resource
descent, endpoint-credit Farkas, source-perspective worry fiber, twist-ladder
dual, signed-polymer/Dirichlet network, HYP-3116 circuit missing-input ledger,
HYP-3117 proof-circuit past-work compiler,
HYP-3118 coordinate-resurrection sheaf,
HYP-3116 circuit-certificate vector,
endpoint-circuit `Phi`, normal-fan
barcode, Ostrowski/Pell carry wall, Vitali/anti-Poisson width debt, unit
endpoint sieve, and sexy-prime local-channel echo.

Top bridges: HYP-3117 proof-circuit past-work compiler to the Lean frontier
packet (`32`), HYP-3118 coordinate-resurrection sheaf to the Lean frontier
packet (`28`), Q27/Q31 resource descent to loaded-denominator normalization
(`28`), observer-cut payload orbit to observer gluing (`27`), HYP-3116
circuit-certificate vector to the Lean frontier packet (`27`), finite-address
`Phi` tuple to the Lean frontier packet (`26`), and source-perspective worry
fiber to Lee-Yang ear payload (`26`).  Tournament Analysis over past-work
carrier vertices is transitive with one Hamiltonian path led by
`finite_address_phi_tuple`, then `observer_cut_payload_orbit`, then
`circuit_certificate_vector`, then `proof_circuit_past_work_compiler`, then
`coordinate_resurrection_sheaf`, then
`q27_q31_resource_descent`.

**Next:** build a packet-row schema with
`finite_address_phi_tuple_status`, `observer_cut_payload_orbit`,
`circuit_certificate_vector`,
`proof_circuit_past_work_compiler`,
`coordinate_resurrection_status`,
`coordinate_resurrection_cover`,
`repair_cover_rank`,
`live_section_type`,
`q27_q31_resource_status`, `twist_ladder_dual_status`,
`source_perspective_worry_fiber`, `endpoint_credit_farkas_certificate`,
`endpoint_circuit_phi_gate`, `missing_input_vector`,
`ostrowski_beatty_pell_carry_wall`,
`dirichlet_polymer_conductance`, `vitali_antipoisson_width_debt`, and
`terminal_exit_or_named_debt`.  Run it on HYP-2963/HYP-3107 residual packets,
HYP-3098 observer-gluing rows, and the THM-573 level-7 residual.  Keep the
sexy-prime connection as local residue-sieve bookkeeping only unless analytic
prime-distribution input is added.
Namespace: HYP-3120 / LTI-256 / LTT-154 / T1195 / OPEN-Q-108.  Integrates
incoming HYP-3116 as the proof-circuit missing-input subcarrier and incoming
HYP-3117 as the proof-circuit past-work compiler carrier.  HYP-3118 is
integrated as the completed coordinate-resurrection sheaf carrier with the
base-stalk/live-section repair-cover schema.

## Lead codex-2026-06-27-S265: LRC14 irrational/transcendental approximation sidecar

**Status:** EVIDENCE / exact interval-margin scout and sidecar synthesis; not
proof (HYP-3114/T1190).
**Readout:** Added
`05-knowledge/hypotheses/HYP-3114-lrc14-irrational-transcendental-approximation-sidecar.md`,
`04-computation/lrc_irrational_transcendental_approximation_codex_s265.py`,
`05-knowledge/results/lrc_irrational_transcendental_approximation_codex_s265.out`,
and
`07-reflections/lrc14-irrational-transcendental-approximation-sidecar-codex-s265.md`.
The lane imports irrational and transcendental approximation into LRC14 only
through a retained witness interval and margin: if `t` is a witness with
margin `delta`, then every rational `p/q` with
`|t-p/q| < delta/max(s_i)` is a witness.

Exact direct-time scout: AP13 has no positive component; `AP12_tail84` has
`8` components, widest length `3/1960`, and grid bound `q>653`;
divisor-loaded `loaded_B6={1,...,11,13,5040}` has `64` components, widest
length `1/5880`, and grid bound `q>5880`; `single_tail168` has widest length
`23/11760` and grid bound `q>511`.  Algebraic/transcendental/Liouville-like
targets enter the positive intervals through continued-fraction convergents,
but those hits are sidecar examples; the proof carrier is the interval,
margin, max-speed robustness radius, exceptional-approximant ledger, and
terminal observer route.

**Next:** replace direct-time intervals by THM-565 normalized slow/ruler
coordinates in the HYP-3088/3089 repair.  Then add `witness_interval`,
`endpoint_margin`, `robust_approximation_radius`,
`continued_fraction_first_hit`, `partial_quotient_spike`,
`irrationality_measure_status`, `exceptional_approximants`,
`liouville_spike_schedule`, destroyed coordinate, and terminal exit to
HYP-3098 observer-gluing rows and HYP-3112 ear-payload rows.
Namespace: HYP-3114 / LTI-251 / LTT-149 / T1190 / OPEN-Q-108.

## Lead codex-2026-06-27-S262b: LRC14 Lee-Yang ear-payload root-motion route

**Status:** EVIDENCE / exact scout and proof-route proposal; not proof
(HYP-3112/T1188).
**Readout:** Added
`05-knowledge/hypotheses/HYP-3112-lrc14-lee-yang-ear-payload-extremality.md`
and
`07-reflections/lrc14-lee-yang-ear-payload-extremality-codex-s262.md`.
This is the exact one-runner payload refinement of HYP-3109's root-curve
ear-map, HYP-3108's Lee-Yang/Savitch atlas, and HYP-3111's carrier-sidecar
lane.  The HYP-3103 miss-count PGF root signal becomes proof-facing only after
the one-runner extension payload is retained.  For `F=E union {a}`, the exact
payload

```text
A_t(E,a)=P(N_E=t and a hits a sector empty for E)
q_F[t]=q_E[t]-A_t+A_{t+1}
```

is the observer-extension/cut coordinate controlling root motion.  The S262b
scout shows AP/consec and even-AP have `real=0/6` and
`dist(roots,[-1,0])=0.9119`; `single_far_21` is also complex-rooted but much
closer to the danger interval (`0.2786`); broken/spread rows have roots on the
danger interval.  Final nested AP `+7` has lower payload level
(`A_mean=1.965291`) than final far `+21` (`A_mean=2.993492`), suggesting
nested low-level ears stabilize roots while far high-level ears retain
root-collision debt.
**Next:** build `lrc14_lee_yang_ear_payload_ledger` over HYP-2963 and the
THM-573 `<=6` multiples-of-7 residual.  Emit `G_E`, root metrics,
last legal ear, `A_t`, nested/far ear type, parity and mean payload,
negative-interval contact, destroyed coordinate, and terminal exit.  Test
whether every root approaching `[-1,0]` is routed to high-mean payload,
nonnested ear debt, component-bound debt, first-obstruction debt,
K33/THM-572 debt, or AP/GW boundary status.
Namespace: HYP-3112 / LTI-249 / LTT-147 / OPEN-Q-108.
## Lead codex-2026-06-27-S265: LRC14 two-map root-lattice-ear extremality synthesis

**Status:** SYNTHESIS + exact map scout; not proof (HYP-3113/T1189).
**Readout:** Added
`05-knowledge/hypotheses/HYP-3113-lrc14-two-map-root-lattice-ear-extremality.md`,
`04-computation/lrc14_two_map_root_lattice_ear_extremality_codex_s265.py`, and
`05-knowledge/results/lrc14_two_map_root_lattice_ear_extremality_codex_s265.out`.
Rebased after HYP-3108/HYP-3109, the scout converts the
Savitch/Bravais/Lee-Yang/phi4/ear-decomposition prompt into two coupled maps.
The root-curve map moves from `single_value_p0` to
full PGF coefficients, root locus, Lee-Yang zero-free regions, discriminant
breaks, quartic cumulant stabilization, and tournament-root spectra; its
fingerprint has `score_hist={0:1,1:1,2:1,3:1,4:1,5:1,7:3,9:1}`, one directed
3-cycle, and a nontrivial SCC among Lee-Yang, discriminant, and Iomega-root
signals.  The memory-lattice-ear map is a strict certificate ladder from raw
runners through relation-lattice basis, successive minima, Bravais shape,
Savitch midpoint sidecars, strong/odd/nested ear spines, first-obstruction
cocycles, and packet-sheaf legal exits (`score_hist={0:1,...,9:1}`, one
Hamiltonian path).
**Next:** join HYP-3103 PGF-zero data to the HYP-3104 maximizer signal atlas.
Add `PGF_zero_locus_signature`, `Lee_Yang_confinement_margin`,
`PGF_discriminant_break_index`, `quartic_cumulant_stabilizer`,
`Savitch_midpoint_sidecar_depth`, `Bravais_relation_shape_class`,
`successive_minima_anisotropy`, `ear_certificate_type`,
`odd_ear_parity_debt`, `nested_ear_crossing_defect`, and
`root_lattice_ear_resonance_portfolio` to packet rows.  Test whether cap false
positives and one-swap exchange traps coincide with root collisions, small
Lee-Yang margins, anisotropic relation lattices, or non-nested ear defects.
Namespace: HYP-3113 / LTI-250 / LTT-148 / T1189 / OPEN-Q-108.

## Lead codex-2026-06-27-S259a: LRC14 normal-fan Cech component-bound route

**Status:** SYNTHESIS / proof-route target; not proof (HYP-3101/T1179).
**Readout:** Added
`05-knowledge/hypotheses/HYP-3101-lrc14-component-bound-via-normal-fan-cech-barcode.md`
and
`07-reflections/lrc14-component-bound-normal-fan-cech-barcode-codex-s259.md`.
The HYP-3096 polynomial witness route needs `components(L_14(S)) <= A0`, not
just positive measure.  HYP-3101 merges HYP-3025 closed danger-arc Cech
topology, HYP-3015 barcode persistence, HYP-3018 active normal-fan supports,
and HYP-3071 first-tooth observability into a finite-chamber target: each
THM-573 residual non-tight row either has bounded normalized direct
`L_14` components, is an AP/GW closed-boundary H1 equality atom, or emits named
F7/THM-572 good-cover quotient debt.  Incoming S258 supplies the first exact
component-growth samples (`42`, `102`, and `860` direct components on live
residuals), and THM-577 strengthens the cap chart while leaving the topology
debt open.  The S259 Lean observer-gluing frontier makes bounded component
packets producers for `ObserverGluingCertificate`, while S65 shows cap
exchange is non-transitive at `j=5` and cannot replace topology control.
**Next:** build `lrc14_normal_fan_cech_component_ledger` over HYP-2963 and the
`<=6` multiples-of-7 residual with `closed_arc_cech_beta`,
`open_arc_component_count`, `boundary_cocircuit_facet_word`,
`owner_current_word_mod_14`, `bar_count_at_height_1_14`,
`minimum_bar_persistence`, `peak_bottleneck_support_word`,
`normal_fan_chamber_id`, `component_bound_status`, `measure_floor_status`,
`finite_ruler_threshold_D`, destroyed coordinate, and terminal exit.  Prove a
bounded chamber theorem or isolate the first AP/GW/F7 topology defect.
Namespace: HYP-3101 / LTI-240 / LTT-138 / OPEN-Q-108.

## Lead codex-2026-06-27-S259b: LRC14 first-obstruction syndrome ledger

**Status:** SYNTHESIS / proof-interface target; not proof (HYP-3102/T1180).
**Readout:** Added
`05-knowledge/hypotheses/HYP-3102-lrc14-first-obstruction-cocycle-generation.md`
and
`07-reflections/lrc14-first-obstruction-cocycle-generation-codex-s259.md`.
The HYP-3095 observer-chart gluing route needs a finite rule for legal
forgetting.  HYP-3102 defines the hidden observer-cut payload difference over a
visible fiber as a first obstruction cochain; the quotient is legal only when
that cochain is zero, reconstructed, exact, generated by named certificate
cycles, dual-annihilated, descended, AP/GW boundary-stopped, or routed to the
named F7/THM-572 state-lift coordinate.  This promotes HYP-3071's
rank-12-of-13 cycle matrix from scout result to packet-cochain proof target.
**Next:** build `lrc14_first_obstruction_syndrome_ledger`; for each mixed
route/status fiber emit quotient name, next observer, visible automorphism
group, payload orbit, first nonzero sidecar stage, obstruction basis vector,
certificate-cycle image status, dual-annihilator status, family descent,
AP/GW boundary stop, F7/THM-572 state-lift status, and terminal exit.  Then
replace the HYP-3071 template rows with actual HYP-2963 packet cochains.  Use
the incoming S258 observer-gluing scout as the first sample row source, and
test THM-577's symbolic cap finite remainders before declaring a new
Pascal/cap obstruction basis atom.  Use S31ah's tournament-certificate engine
as a generator catalog, but keep S65's warning: `c5`/power-sum holes and
forbidden-H alpha events are different obstruction mechanisms.
Namespace: HYP-3102 / LTI-241 / LTT-139 / OPEN-Q-108.

## Lead codex-2026-06-27-S256: LRC14 observability-sheaf gluing ledger

**Status:** SYNTHESIS / proof-route abstraction; not proof (HYP-3093).
**Readout:** Added
`05-knowledge/hypotheses/HYP-3095-lrc14-route-sophistication-hidden-observability-sheaf.md`
and
`07-reflections/lrc14-route-sophistication-hidden-observability-sheaf-codex-s256.md`.
The repo history is best read as increasing discipline about legal forgetting:
AP/GW is the additive equality shadow, tournaments became relation/sidecar
observers, HYP-2963 packets made exact scale and endpoint owner first-class,
HYP-2990 named the zipper/no-free-slider rule, HYP-3083 made finite-address
branch closure the spine, and the paper/c-lift work explains why `14=2*7`
leaves a degree-7 invisible direction.  THM-575 corrects the raw Conjecture 7.1
observer: small denominator time is unstable under apex loading, so the live
route must use normalized slow/ruler-coordinate witness mass after THM-573.
**Next:** build an `lrc14_observer_gluing_ledger` over the THM-573 residual
with fields for arithmetic chart (`I(13,7,1)`, mod-7/c-lift status),
normalized arc chart (`P,E,V`, `G(P,E)`, component and finite-ruler data),
cap chart (HYP-3090 pairwise-avoidance ratios), moment chart (HYP-3085
gK8/`S2`/Perron data), branch chart (HYP-3092
nested-refinement/K33 cross-handoff), forgotten coordinate, overlap map, and
terminal exit/debt.  Test whether every quotient used by the route is
reconstructible, dual-annihilated, fiber-constant, or explicitly routed.
Namespace: HYP-3093 / HYP-3092 / HYP-3090 / OPEN-Q-108.
## Lead codex-2026-06-27-S257: LRC14 three equivalence shadows and Pascal pair-mass

**Status:** SYNTHESIS / invariant proposal; not proof
(HYP-3097/T1175).
**Readout:** Added
`05-knowledge/hypotheses/HYP-3097-lrc14-three-equivalence-shadows-pascal-pair-mass.md`
and
`07-reflections/lrc14-three-equivalence-shadows-pascal-pair-mass-codex-s257.md`.
This is the Pascal/pair-mass companion to HYP-3091's verified
three-sameness fiber on the lonely set, HYP-3092's verified pair-normalized
cap mass, and HYP-3093's broad equivalence-triad audit.  The prompt constants
have pair-apex form:

```text
C(14,2)=91
1001=11*91=C(14,4)
2002=22*91=C(14,5)
3003=33*91=C(14,6)
4004=44*91=2*C(14,5)=C(14,4)+C(14,6)
```

Read with HYP-3090/HYP-3092, `cap_k=C(k+1,2)/91` is exact for `k>=10`, but
`cap_9=45/91-1/4004`; the `4004` signal is therefore the affine pair-mass
completion and the one-unit defect denominator, not just Pascal numerology.
**Next:** add `pascal_pair_mass_unit`, `triangular_cap_shadow`,
`cap_defect_numerator`, `sector_pair_scissors_signature`,
`farey_additive_lane_mod_91`, `farey_product_lane_mod_91`,
`level7_lift_status`, destroyed coordinate, and terminal exit to a HYP-2963
packet audit.  Test whether equal base counts or equal cap shadows split by
sector-pair scissors data.  Keep the assumption challenge explicit: raw
runners, raw counts, and raw measure are shadows; the proof vertex is the
invariant packet that preserves the LR predicate or names what was destroyed.
Namespace: HYP-3097 / T1175 / OPEN-Q-108.

## Lead codex-2026-06-27-S255: LRC14 Conjecture 7.1 raw-time refutation and normalized repair

**Status:** PROVED refutation of paper Conjecture 7.1 as stated + repaired
proof-interface route (THM-575/HYP-3088).
**Readout:** Added
`01-canon/theorems/THM-575-lrc14-conjecture-7-1-refuted-by-divisor-loading.md`,
`05-knowledge/hypotheses/HYP-3088-lrc14-conj71-refutation-and-normalized-arc-repair.md`,
`04-computation/lrc14_conj71_refutation_normalized_arc_codex_s255.py`,
`05-knowledge/results/lrc14_conj71_refutation_normalized_arc_codex_s255.out`,
and
`07-reflections/lrc14-conj71-refutation-normalized-arc-repair-codex-s255.md`.
The user's paper-method bridge is correct at the lift level: `k+1=14=2*7`
breaks Proposition 4.1's prime-field argument and points to `c=2,7` descent,
matching dyadic work plus THM-573.  The raw Conjecture 7.1 equation is false,
however.  For `S_B={1,...,11,13,84*lcm(1..B)}`, the loaded apex kills all
denominators `d<=B`, while `t=1/12+1/(2N)` is a strict witness for `B>=6`; thus
the rows are primitive and non-tight.  Direct time components shrink (`B=6`
largest `1/5880`), so a direct largest-time-arc floor is not the proof route.
**Next:** build the normalized residual ledger over THM-573's `<=6` multiples
of `7` core: fields should include `P,E,V`, `count_7_divisible`,
`G(P,E)` measure, normalized component count, largest normalized arc,
HYP-2072 `I(k,p,1)`/mod-7 sieve status, finite-ruler threshold, finite-check
status, and terminal exit.  The theorem target is uniform normalized arc
floor plus component bound, not a raw denominator theorem.
Incoming mac-mini-S61 supplies the matching first measurements: exact
`I(13,7,1)` equals covering mod `7`, one `c=2` lift equals covering mod `14`,
and direct arcs for `{1..12,14V}` cross over from finite long arcs to `1/V`
decay near `V*`.  Use this as support for the normalized peel, not for raw
Conjecture 7.1.
Namespace: THM-575 / HYP-3088 / OPEN-Q-108.
Namespace: THM-574 / HYP-3088 / OPEN-Q-108.
## Lead codex-2026-06-27-S252/S254-postrebase: HYP-3088 largest-arc target inside finite-address LRC14 spine

**Status:** SYNTHESIS / proof-target integration; not proof
(HYP-3083 + HYP-3088/HYP-3089 + HYP-3087/HYP-3090/HYP-3091).
**Readout:** Incoming S31ag reframes the S252/S254 finite-address spine by
identifying LRC14 as the composite `k+1=14=2*7` case of the polynomial-method
paper.  THM-573 is the clean `c=7` lift, the dyadic tower is the `c=2` side,
and the live core remains primitive covering rows with `<=6` multiples of `7`
after q-witness, one-large-speed, and aliasing/margin exits.  HYP-3088 makes
the final target sharper: prove a uniform largest-lonely-arc floor for the
normalized slow/ruler-coordinate lonely carrier, or name the first finite-address packet where the
Conjecture 7.1(13) translation fails.  THM-576/HYP-3090 narrows the cap side
to the triangular pairwise skeleton plus `k=8,9` deviation and order-3 break
debt.  HYP-3094 then remains the local shuttle:
nested-refinement rows feed O2 covering discharge, cross-handoff rows feed
O3/THM-572 lift debt, and positive safe mass alone is not the separator.
**Next:** extend `finite_address_branch_closure` with
`polynomial_composite_lift_status`, `direct_lonely_arc_count_status`,
`largest_lonely_arc_floor`, `cap_ratio_or_deviation_status`, `grid_class`, `active_binder_owner_word`,
`endpoint_owner_transition_word`, and `bridge_status`; run it first on the
HYP-2963 bank representatives and the O2/O3 shuttle examples, then measure
whether the THM-565 scale-separated arc-count bound transfers from the
controlled witness object to the direct lonely set.
Namespace: HYP-3083 / HYP-3087 / HYP-3088 / HYP-3089 / HYP-3090 / HYP-3091.
Namespace: THM-574 / HYP-3089 / OPEN-Q-108.
## Lead codex-2026-06-27-S255b: CRT bridge after raw Conjecture 7.1 refutation

**Status:** SYNTHESIS / corrected proof-target refinement; not proof
(THM-575/THM-576/THM-574/HYP-3088/HYP-3089/HYP-3090/HYP-3091/HYP-3092/HYP-3094/T1172/LTI-236/LTT-134).
**Readout:** Rebased the polynomial-method bridge over incoming S61.  The raw
denominator/largest-time-arc reading is false by THM-575 and HYP-3088, so the
old `D<=213` direct route must not be used as a global theorem.  The corrected
bridge keeps the user's core insight but changes the carrier: `14=2*7` still
forces `c=7` and `c=2` lift channels; THM-574 gives the c-lift family,
THM-573 supplies the `c=7` lift; HYP-3089 verifies
`I(13,7,1)=covering mod 7` and the `c=2` lift to covering
mod `14`; the V* crossover splits bounded apex direct checks from large-apex
Node-3/gK8 peeling.  HYP-3090/HYP-3092/THM-576 make the cap side
pair-Pascal / pairwise-triangular for `k>=10`, leaving only the `k=8,9`
deviation constants as finite analytic-substitute debt; HYP-3091 separates the
mod-41 bounded-core scissors scale from apex V*, and HYP-3094 resolves the
duplicate K33-shuttle namespace.
**Next:** build the normalized residual witness ledger with CRT factor, lift
status, `P,E,V` normalization, direct and normalized component counts,
cap-formula/deviation regime, mod-41 `D`, V* regime, Node-3 peel status,
gK8/p0 status, HYP-3094 shuttle status, finite-ruler threshold, and first
missing sidecar.
Namespace: THM-574 / THM-576 / HYP-3088 / HYP-3089 / HYP-3090 / HYP-3091 / HYP-3092 / HYP-3094 / T1172 / LTI-236 / LTT-134.
## Lead codex-2026-06-27-S255: LRC14 polynomial-method witness-route ledger

**Status:** SYNTHESIS / proof-obligation ledger; not proof
(HYP-3096/T1176/LTI-237/LTT-135).
**Readout:** Added
`05-knowledge/hypotheses/HYP-3096-lrc14-polynomial-method-witness-route-ledger.md`
and
`07-reflections/lrc14-polynomial-method-witness-ledger-codex-s255.md`.
The S31ag arXiv:2604.23906 bridge is now a concrete ledger: for `k=13`,
the paper's composite `k+1=14=2*7` fallback matches the project descent
`14 -> 7 -> 2`.  THM-573 is the `c=7` lift, while the remaining
`<=6`-multiples-of-7 primitive covering core carries the `c=2` dyadic and
analytic debt.  The paper's `I(k,p,1)` count is the finite denominator-grid
question `L(S) ∩ (1/p)Z`; the project substitute is the largest-arc witness
route:

```text
mu(L(S)) >= m0
components(L(S)) <= A0
=> largest_arc(L(S)) >= m0/A0
=> witnesses in (1/d)Z for all d >= ceil(A0/m0).
```

This proves Conjecture 7.1(13) for the residual and then LRC14, but only after
the direct `1/14` component bound or an equivalent reduction is proved.  The
q-cusp finite-principal-part discipline becomes a finite bad-denominator
budget, and the hyperoperation grid stays as a retained address chart
`(p,q),p+q,p*q,powers` with CRT lift status and terminal exits.
**Next:** attach a `polynomial_method_witness_ledger` to HYP-2963 rows and
outside-bank normalizer attempts with CRT lift status, `I_discrete_grid_status`,
lonely measure, direct component bound, largest arc, denominator threshold,
hyperoperation address, finite bad-denominator budget, destroyed coordinate,
and exit.  The first theorem task is the direct `1/14` lonely-set component
bound for the THM-573 residual, or a controlled reduction from THM-565's
maxgap witness object.  Keep the equivalence direction honest: uniform
largest-arc route -> Conjecture 7.1(13) -> LRC14.
Namespace: HYP-3096 / T1176 / LTI-237 / LTT-135.

## Lead codex-2026-06-27-S254: Sexy prime pair sieve transfer

**Status:** SYNTHESIS / proof-interface map; not proof
(HYP-3086/T1168).
**Readout:** Added
`05-knowledge/hypotheses/HYP-3086-sexy-prime-pair-sieve-transfer.md`,
`04-computation/sexy_prime_pair_surface_codex_s254.py`,
`05-knowledge/results/sexy_prime_pair_surface_codex_s254.out`, and
`07-reflections/sexy-prime-pair-sieve-transfer-codex-s254.md`.
The repo's strongest connection to the sexy prime conjecture is the existing
pair-lens coordinate system: sexy primes are the fixed half-gap ray `h=3`,
`(p,p+6)=(m-3,m+3)`.  The local midpoint sieve kills `m=+h` or `m=-h mod q`;
for `h=3` the two bad residues collapse modulo `3`, giving the
Hardy-Littlewood chord factor `2` relative to twins.  S254 confirms the shape
through `10^6`: `8169` twin pairs and `16386` gap-6 pairs, ratio `2.006`
against predicted `2.000`.
Incoming S31af/S60 sharpens the LRC side of the analogy: THM-573 is a level-7
lift sieve reducing the LRC residual to at most six multiples of `7`, and
HYP-3084 replaces the literal Clebsch-design covariance shortcut with a
reflection-symmetric Perron/`3x3` block certificate target.  These reinforce
the finite-residual-ledger method but do not supply prime-pair distribution.
**Next:** do not treat this as a proof.  Build a `sexy_prime_pair_ledger` with
midpoint, gap, local bad-residue word, surviving primorial class, sieve weight
type, almost-prime status, parity debt, distribution modulus, exceptional
modulus flag, prime-power sidecar, carry debt, and terminal exit.  The open
proof obligations remain lower-bound fixed-gap sieve, parity breaking, and
sufficient prime distribution in arithmetic progressions or a replacement
structural input.  Keep THM-503/THM-518 as a guardrail: the LRC singular object
is not the Hardy-Littlewood Euler product.
Namespace: HYP-3086 / T1168.
## Lead codex-2026-06-27-S254: LRC14 hyperoperation grid address

**Status:** SYNTHESIS / operation-address proof carrier; not proof
(HYP-3087/T1169/LTI-233/LTT-131).
**Readout:** Added
`05-knowledge/hypotheses/HYP-3087-lrc14-hyperoperation-grid-address.md`
and
`07-reflections/lrc14-hyperoperation-grid-address-codex-s254.md`.
The user-specified hyperoperation hierarchy on `(p,q)` and the older grid
tiled by `x+2` and `x*2` should be read as a Farey operation-address chart:
`p+q` is the additive/horizontal lane, `p*q` is the product/vertical lane, and
`q^p,p^q` are lacunary stress lanes.  The space-filling curve is only a
scheduler unless every cell retains root packet, operation lane, current
danger, endpoint owner, level-7 status, finite address, destroyed coordinate,
and terminal exit.  THM-573 sharpens the arithmetic opening: the residual is
now rows with `<= 6` multiples of `7`, so the vertical lane must record
level-7 lift status rather than just multiple-of-14 status.  HYP-3090 adds
cap/deviation status, and HYP-3094 adds the covering/K33 shuttle exit grammar
for nested-refinement and cross-handoff grid classes.
level-7 lift status rather than just multiple-of-14 status.  Incoming HYP-3092
adds the covering/K33 shuttle exit grammar for nested-refinement and
cross-handoff grid classes.
**Next:** build a `hyperoperation_grid_address` ledger over HYP-2963 and
outside-bank normalizer attempts with fields for `count_7_divisible`,
`level7_lift_status`, `(p,q)`, `p+q`, `p*q`, power-stress word,
cap/deviation status, space-filling successor, danger deficit, endpoint owner,
finite address, terminal exit, and branch/debt status.  Test whether danger-weighted operation
cells split the THM-573 residual core into q-witness, covering/Node3,
K33/THM-572, protected branch closure, or named residual debt.
Namespace: HYP-3087 / T1169 / LTI-233 / LTT-131.

## Lead codex-2026-06-27-S252: LRC14 finite-address branch-closure spine
## Lead codex-2026-06-27-S252/S253/S31af: LRC14 finite-address branch-closure spine

**Status:** SYNTHESIS / remaining-obligation map; not proof
(HYP-3083/T1167/LTI-232/LTT-130).
**Readout:** Expanded
`05-knowledge/hypotheses/HYP-3083-lrc14-hurwitz-finite-address-branch-closure-spine.md`
from a reservation into the integrated proof map, then rebased it over the
S59 covering-bound redirect, the S252b/S253 q-cusp/Hurwitz branch-closure
spine, HYP-3084 Morita-gamma/discrepancy sidecar, S31af THM-573 level-7 lift
sieve, THM-576/HYP-3090 cap skeleton, and S254 HYP-3094 covering/K33 shuttle.

```text
no-apex row
  -> direct t=1/14 witness

covering row with a multiple of 14
  -> q-witness gate
  -> THM-573 level-7 lift sieve or <=6 multiples-of-7 residual
  -> one-large-speed peeler or top-balanced residue
  -> finite-address packet
  -> protected branch graph
  -> strict witness / AP-GW boundary / C27 petal /
     covering moment / K33-THM-572 lift / named residual debt
  -> formal sidecar closure
```

Hurwitz/Markov/Pell walls, q-cusp principal parts, sixth-power collision
certificates, route-median centers, apex-periodic covering rows, Hensel/Morita
discrepancy carriers, endpoint owners, and tournament kernels are legal proof
data only after their finite address and destroyed coordinate are retained or
routed to named debt.  THM-573 closes every branch with at least seven speeds
divisible by `7`, while the S31af covering-margin scout refutes a uniform
`>1/13` shortcut through dilated/AP aliasing examples.  The live obligations
are the global finite-address normalizer for the `<=6` multiples-of-7
covering residual, covering-moment/OPEN-Q-108 discharge for
`nested_refinement` packets, K33/THM-572 lift for `cross_handoff` packets,
branch-closure theorem, integer-vs-real/formal closure, and AP/GW census only
if a boundary-equality route is chosen.
**Next:** build `finite_address_branch_closure` with `source_row_or_family`,
`apex_divisible_by_14_flag`, `multiple_of_7_profile`,
`level7_lift_sieve_status`, `covering_margin_aliasing_status`,
`grid_class`, `active_binder_owner_word`, `endpoint_owner_transition_word`,
`low_apex_top_balanced_status`, `normalizer_exit_attempted`,
`finite_address_word`, `preserved_lrc_predicate`, `destroyed_coordinate`,
`required_sidecar_or_debt`, `p_adic_discrepancy_sidecar_status`,
`protected_branch_node`, `bridge_status_raw`, `bridge_status_protected`,
`terminal_exit`, `lean_formalization_status`, and `residual_debt_name`.  Treat
any naked bridge as missing address/debt rather than as a scalar theorem
shortcut.  First target the low-apex covering-moment family, then K33/THM-572
lift construction.
Namespace: HYP-3083 / T1167 / LTI-232 / LTT-130.

## Lead codex-2026-06-26-S249: Robbins-Fermat-Catalan branch tournaments

**Status:** SYNTHESIS / proof-network guardrail; not proof
(HYP-3081/T1165/LTI-230/LTT-128).
**Readout:** Added
`05-knowledge/hypotheses/HYP-3081-robbins-fermat-catalan-branch-tournaments.md`,
`07-reflections/robbins-fermat-catalan-branch-tournaments-codex-s249.md`, and
`poke-forum/posts/20260626-robbins-fermat-catalan-branch-tournaments/post.md`.
Robbins' theorem gives the quotient test: a connected proof graph can be
strongly oriented exactly when it has no bridges, so a forgotten coordinate
that becomes a naked bridge is illegal controlled forgetting.  The branch
graph should be built downstream of the HYP-3078 q-cusp scout and HYP-3079 Lean
q-cusp ledger, with reverse verification, owner/H1 closure, Fermat-Catalan
no-lift guards, q-cusp polar-debt guards, and named residual exits.
**Next:** emit HYP-2963 branch rows with `branch_id`, `bridge_status`,
`reverse_verification_mode`, `endpoint_kernel_iso_class`,
`achievable_tournament_kernel_set`, `power_lift_guard`,
`q_cusp_polar_debt_guard`, and `destroyed_coordinate_exit`, then compute naked
bridges before and after sidecar closure.
Namespace: HYP-3081 / T1165 / LTI-230 / LTT-128.

## Lead codex-2026-06-26-S250: Branch-kernel orientation audit

**Status:** EVIDENCE / bounded proof-interface audit; not proof
(HYP-3082/T1166/LTI-231/LTT-129).
**Readout:** Added
`04-computation/lrc14_branch_kernel_orientation_codex_s250.py`,
`05-knowledge/results/lrc14_branch_kernel_orientation_codex_s250.out`,
`05-knowledge/hypotheses/HYP-3082-branch-kernel-orientation-audit.md`,
`07-reflections/lrc14-branch-kernel-orientation-audit-codex-s250.md`, and
`poke-forum/posts/20260626-lrc14-branch-kernel-orientation-audit/post.md`.
The default stored audit covers `21913` packets and `7235` hard non-AP/GW
packets.  Raw scalar-star quotienting has `5` naked bridges; the protected
branch graph has `80` nodes, `83` edges, `69` bridges, `0` naked bridges, and a
strongly orientable contracted core after route sections, Haar/grid exits,
no-lift guards, q-cusp/Lean finite-tail guards, finalizer gates, and named
residual debt are attached.
**Next:** rerun the audit after any HYP-2963 packet-bank expansion and export
bridge witnesses with raw/protected status, responsible sidecar, endpoint
kernel class, residual exit, and contracted-core orientation status.
Namespace: HYP-3082 / T1166 / LTI-231 / LTT-129.

## Lead codex-2026-06-27-S243: Hurwitz-Markov-Pell cannonball carrier

**Status:** SYNTHESIS / finite arithmetic sidecar scout; not proof
(HYP-3075/T1158/LTI-223/LTT-121).
**Readout:** Added
`04-computation/lrc14_hurwitz_markov_pell_cannonball_s243.py`,
`05-knowledge/results/lrc14_hurwitz_markov_pell_cannonball_s243.out`,
`05-knowledge/hypotheses/HYP-3075-lrc14-hurwitz-markov-pell-cannonball-carrier.md`,
`07-reflections/lrc14-hurwitz-markov-pell-cannonball-codex-s243.md`, and
`poke-forum/posts/20260627-hurwitz-markov-pell-cannonball-lrc14/post.md`.
The scout finds the nontrivial cannonball square `1^2+...+24^2=70^2`, with
`70=Pell P6`; its neighboring Pell values `29=Pell P5` and `169=Pell P7` are
fixed-coordinate-2 Markov branch wall numbers and satisfy `29*169-70^2=1`.
**Next:** add the Hurwitz/Markov/Pell/cannonball sidecar fields to a Q27 or
HYP-2963 packet sample and run status/route purity: the visible token should
split by continued-fraction period, Markov depth, Pell unit, endpoint shell,
quadratic carry, owner coordinate, and legal exit before any scalar is used.
Namespace: HYP-3075 / T1158 / LTI-223 / LTT-121.
## Lead codex-2026-06-27-S248: LRC14 sixth-power certificate extension

**Status:** SYNTHESIS / bounded Diophantine proof-carrier scout; not proof
(HYP-3080/T1164/LTI-229/LTT-127), downstream of S244/HYP-3076.
**Readout:** Added
`04-computation/lrc14_sixth_power_certificate_extension_s248.py`,
`05-knowledge/results/lrc14_sixth_power_certificate_extension_s248.out`,
`05-knowledge/hypotheses/HYP-3080-lrc14-sixth-power-certificate-extension.md`,
`07-reflections/lrc14-sixth-power-certificate-extension-codex-s248.md`, and
`poke-forum/posts/20260627-lrc14-sixth-power-certificate-extension-s248/post.md`.
S244 already treats `a^6+b^6+c^6=d^6+e^6+f^6` as native support-six data and
`a^6+b^6=d^6+e^6` as a rank-lowered square-cube shadow.  S248 adds a
certificate ledger around that split: positive unordered pairs through `250`
have `0` nontrivial collisions; positive unordered triples through `80` have
`5` collision certificates, including primitive `(3,19,22)` versus
`(10,15,23)`.
Tournament Analysis uses proof obligations and sidecar carriers rather than
runners or integers, with `5` rank-sensitive edge flips.
**Next:** add the sixth-power collision sidecar fields to HYP-2963 route
triples that invoke power-lift, Fermat-Catalan, Roth-Minkowski, or
simplex/route-triple language.  Run S240 legal route-state closure and classify
failed medians as missing collision certificate, CRT residue word, height
fence, gated route sidecar, or THM-572/F7 debt.
Namespace: HYP-3080 / T1164 / LTI-229 / LTT-127.

## Lead codex-2026-06-26-S238: LRC14 cross-carrier pullback resonance

**Status:** EVIDENCE / finite carrier-portfolio scout; not proof
(HYP-3072/T1154/LTI-219/LTT-117).
**Readout:** Added
`04-computation/lrc14_cross_carrier_pullback_resonance_codex_s238.py`,
`05-knowledge/results/lrc14_cross_carrier_pullback_resonance_codex_s238.out`,
`05-knowledge/hypotheses/HYP-3072-lrc14-cross-carrier-pullback-resonance.md`,
`07-reflections/lrc14-cross-carrier-pullback-resonance-codex-s238.md`, and
`poke-forum/posts/20260626-cross-carrier-pullback-resonance-lrc14/post.md`.
The scout encodes `22` CPI/HYP proof carriers and `9` remaining obligations
over the core alphabet `status, route, exact_scale, topology, owner,
period_deck, analytic_certificate, automaton_partial_cube, crt_padic,
observer_cut, hodge_cycle, formal_exit`.  The first cover of all `23` target
axes appears only at size `9`, so no small universal scalar-like bundle is
visible.  Local obligations have compact covers: partial-cube plus fusion,
owner strip plus CRT tree, Hodge/K33 plus fusion, and observer/rectangle plus
exact Farey scale.
**Next:** emit actual HYP-2963 packet rows with `carrier_pullback_row_id`,
`core_incident_word`, `preserved_lrc_predicate`, `destroyed_coordinate`,
`required_sidecar`, `blindness_pair_id`, `resonance_portfolio_id`,
`status_mixing_result`, `route_mixing_result`, and `legal_exit_status`.  Then
test whether the listed local portfolios make residual coarse fibers
status-pure and route-pure before any new theorem debt is named.
Namespace: HYP-3072 / T1154 / LTI-219 / LTT-117.
## Lead codex-2026-06-26-S239: Renormalized polymer / Dirichlet bridge

Use HYP-3073 as the non-median bridge after HYP-3072's carrier portfolio and
HYP-3071's cycle-class observability matrix.  The two proof angles are:
typed signed-polymer activities for AP/repeated-residue/wide packets, and
Dirichlet/Schur sidecar energy for residual-current discharge.
**Readout:** `04-computation/lrc14_polymer_dirichlet_bridge_codex_s239.py`
and `05-knowledge/results/lrc14_polymer_dirichlet_bridge_codex_s239.out`.
**Status:** finite proof-interface scout and synthesis, not proof. Raw R6
density misorders signed correction (`odd_AP` beats `near_AP` by signed
correction despite fewer relations), wide/Sidon packets have tiny signed
activity, and the sidecar-energy toy network separates raw route conductance
`1/2` from legal sidecar conductance `9` with a one-unit phantom F7 boundary
exit.
**Next:** build the actual HYP-2963 typed-polymer ledger and residual sidecar
graph. Add `signed_polymer_packet_type`, `signed_activity_budget`,
`finite_cell_route`, `dirichlet_boundary_potential`,
`schur_complement_conductance`, `sidecar_energy_exit`, and
`phantom_f7_boundary_atom`, then prove activity summability or positive
conductance preservation before quotienting.
Namespace: HYP-3073 / T1155 / LTI-220 / LTT-118.

## Lead codex-2026-06-26-S236: Route-triple center-control addendum

Use HYP-3070 as the raw-vs-legal control layer for the final LRC medianization
interface. Read after S235/HYP-3069: Boolean completion finds raw-projection
ambiguity, while S236 says the route vocabulary itself should first pass a
center-control test.
**Readout:** `04-computation/lrc14_route_triple_center_control_codex_s236.py`
and `05-knowledge/results/lrc14_route_triple_center_control_codex_s236.out`.
**Status:** exact finite scout and synthesis, not proof. Raw route-label clique
has empty centers for `455/455` route triples; legal sidecar tree has unique
centers for `455/455`. Next: instantiate the expected-center pages on actual
HYP-2963 coarse fibers and compare with HYP-3069's Boolean median-completion
obligations.
Namespace: HYP-3070 / T1152 / LTI-217 / LTT-115.

## Lead codex-2026-06-26-S235: Medianized route-center gate

Turn the doubled triangular/simplex/Faulhaber cue into a final LRC
proof-interface check (HYP-3069/T1151/LTI-216/LTT-114). Read after S231
bridge-rank, S232 Hodge-cycle work, S233 Desargues-median failures, and the
S234 owner/root sidecar spine.
**Readout:** `04-computation/lrc14_medianization_route_center_codex_s235.py`
and `05-knowledge/results/lrc14_medianization_route_center_codex_s235.out`.
**Status:** finite scout and synthesis, not proof. Checks 220 route triples,
exposes 122 raw ambiguous triples, and produces 70 median center obligations.
## Lead codex-2026-06-26-S237: LRC14 cycle-class observability matrix

**Status:** EVIDENCE / exact residual-summary scout plus rational
cycle-class template (HYP-3071/T1153/LTI-218/LTT-116).
**Readout:** Added
`04-computation/lrc14_cycle_class_observability_scout_codex_s237.py`,
`05-knowledge/results/lrc14_cycle_class_observability_scout_codex_s237.out`,
`05-knowledge/hypotheses/HYP-3071-lrc14-cycle-class-observability-matrix.md`,
`07-reflections/lrc14-cycle-class-observability-scout-codex-s237.md`, and
`poke-forum/posts/20260626-cycle-class-observability-lrc14/post.md`.  The
scout parses S199/S200 and turns the remaining HYP-2963 residual proof debt
into a first-tooth observability matrix: `arc_topology_compact` separates
`13/15` strict-open coarse ET+unit residual fibers, `coarse_safe_stalk`
separates `15/15`, and all repairs are `owner_strip`.  The companion
rational cycle-class template has basis dimension `13`, known generator rank
`12`, and leaves only `phantom_f7_class` outside the span.
**Next:** replace the template target rows with actual HYP-2963 packet
cochains for topology, owner current, primitive deck, Haar zeta, observer-cut
payload, rectangle/hourglass residue, partial-cube Theta/simplex sidecar,
low-height wall, octahedral curl, Toeplitz scale, median owner/root fields,
and state-lift target; then row-reduce over `Q` and record
`cycle_class_image_status`.
Namespace: HYP-3071 / T1153 / LTI-218 / LTT-116.
Namespace: HYP-3069 / T1151 / LTI-216 / LTT-114.
## Lead codex-2026-06-26-S240: Route-state closure median interface

**Status:** SYNTHESIS / executable proof-interface scout for LRC14
(HYP-3074/T1156/LTI-221/LTT-119).
**Readout:** Added
`05-knowledge/hypotheses/HYP-3074-route-state-closure-median-interface-lrc14.md`,
`04-computation/lrc14_route_state_closure_median_s240.py`,
`05-knowledge/results/lrc14_route_state_closure_median_s240.out`,
`07-reflections/lrc14-route-state-closure-median-interface-codex-s240.md`,
and
`poke-forum/posts/20260626-lrc14-route-state-closure-median-interface-s240/post.md`.
This should be read after the S238 cross-carrier pullback resonance stub
(HYP-3072/T1154), which promotes CPI carrier rows and hidden coordinates into
candidate proof-state vertices.
The carrier models proof witnesses as
`packet / route / certificate / sidecar / discharge` states and accepts a
serious route triple only when legal sidecar closure gives a legal
coordinate-wise median center.
**Next:** build the HYP-2963 medianization ledger.  For each packet route,
emit packet fields, route fields, certificate fields, legal sidecar closure,
discharge fields, median triples entering the state, and failed-median
reasons.  Prove legal uniqueness of closed medians or classify each failure as
missing gated partial-cube sidecar, missing cycle image, missing observer-cut
repair, or explicit F7/THM-572 debt.
Namespace: HYP-3074 / T1156 / LTI-221 / LTT-119, downstream of HYP-3072 / T1154.

## Lead codex-2026-06-26-S232: Hodge-cycle lifting carrier

**Status:** SYNTHESIS / proof-interface sidecar for LRC14
(HYP-3066/T1148/LTI-213/LTT-111).
**Readout:** Added
`05-knowledge/hypotheses/HYP-3066-hodge-cycle-lifting-carrier-lrc14.md`,
`07-reflections/hodge-conjecture-cycle-lifting-lrc14-codex-s232.md`, and
`poke-forum/posts/20260626-hodge-cycle-lifting-lrc14/post.md`.  The carrier
uses the Hodge conjecture only as a realizability discipline: closed,
type-correct, positivity-feasible packet cochains are not proof exits until
they are generated by named certificate cycles, dual-annihilated, descended,
identified as AP/GW boundary classes, or emitted as explicit F7/THM-572
residual classes.
**Next:** build the exact rational cycle-class matrix on a HYP-2963 sample.
Rows are emitted residual cochains; columns are named cycles such as AP/GW,
Fejer, Ramanujan, Haar, endpoint-owner, observer-cut, octahedral face-curl,
partial-cube Theta, simplex face, Roth-Minkowski wall, C27/unital, and
Toeplitz noncollapse, C27/unital, and K33/state-lift generators.  Record
`cycle_class_image_status`; do not treat
`phantom_unresolved` as a proof exit.
Namespace: HYP-3066 / T1148 / LTI-213 / LTT-111.
## Lead codex-2026-06-26-S245: Route-state median-hull scheduler

**Status:** SYNTHESIS / finite final proof-interface scout
(HYP-3077/T1160/LTI-225/LTT-123).
**Readout:** Added
`04-computation/lrc14_route_state_median_hull_scheduler_codex_s245.py`,
`05-knowledge/results/lrc14_route_state_median_hull_scheduler_codex_s245.out`,
`05-knowledge/hypotheses/HYP-3077-route-state-median-hull-scheduler.md`,
`07-reflections/lrc14-route-state-median-hull-scheduler-codex-s245.md`,
and `poke-forum/posts/20260626-lrc14-route-state-median-hull-scheduler/post.md`.
The scout encodes completed proof states as
`packet/route/certificate/sidecar/discharge` coordinate sets.  Legal sidecar
rules are unary Horn implications, so coordinatewise majority gives legal
median centers after closure.  The finite audit uses `41` features, `34` rules,
max premise arity `1`, ten seed states, a `31`-state median hull, `29,791`
checked triples, `raw_illegal_majorities=0`,
`closure_added_features_hist={0: 29791}`, `interval_intersection_failures=0`,
and `0` illegal centers.  The selected serious triples produce scheduler
centers, not terminal exits; future conjunctive guards must be named as
sidecars or checked separately.
**Next:** Apply the medianization schema to actual HYP-2963 route fibers:
AP/GW-C27-K33, q=23 residual capacitors, Moser/fibbinary automatic fibers, and
Fejer/Toeplitz versus Desargues/Beal finalizers.  Record dropped terminal atoms
and the first separating sidecar.
Namespace: HYP-3077 / T1160 / LTI-225 / LTT-123.

## Lead codex-2026-06-26-S227: Moser/fibbinary partial-cube simplex carrier

**Status:** SYNTHESIS / forum-facing proof-carrier sidecar
(HYP-3063/T1145/LTI-210/LTT-108).
**Readout:** Added
`05-knowledge/hypotheses/HYP-3063-moser-fibbinary-partial-cube-simplex-geometry-carrier.md`,
`07-reflections/moser-fibbinary-partial-cube-simplex-lrc14-codex-s227.md`,
and `poke-forum/posts/20260626-moser-fibbinary-partial-cube-lrc14/post.md`.
The carrier upgrades Moser/fibbinary from raw automatic names into
partial-cube/simplex sidecars: fibbinary windows are Fibonacci cubes, Moser
windows are even-coordinate Boolean subcubes, and `n(n+1)=2*T_n` is a
directed-simplex edge layer.  The `5,6,7` geometry motif is retained only with
HYP-3061 axis typing.
**Next:** add `partial_cube_model`, `theta_class_word`,
`gated_subcube_status`, `median_interval_status`, `simplex_face_rank`,
`directed_simplex_edge_count`, `doubled_triangular_layer`,
`fibonacci_cube_window`, and `moser_even_coordinate_subcube` to a HYP-2963
sample already carrying exact `M`, endpoint-owner, closed-arc topology,
magnitude cocycle, automatic state, and geometry-regime signature; then test
whether the automaton/cube quotient becomes route-pure.
Namespace: HYP-3063 / T1145 / LTI-210 / LTT-108.

## Lead codex-2026-06-26-S226: Roth-Minkowski Diophantine lattice fence

**Status:** SYNTHESIS / proof-interface sidecar for Diophantine geometry
(HYP-3062/T1144/LTI-209/LTT-107).
**Readout:** Added
`05-knowledge/hypotheses/HYP-3062-roth-minkowski-diophantine-lattice-fence.md`
and
`07-reflections/roth-minkowski-diophantine-lattice-fence-codex-s226.md`.
Roth's Diophantine approximation theorem and Minkowski's geometry-of-numbers
theorem are imported as a paired controlled-forgetting fence: first finite
low-height wall deletion, then named Minkowski relation-lattice tail, then
Roth exceptional-approximant fence.  Raw volume, raw covolume, and raw
approximation exponent remain scouts unless the lattice, convex-body, algebraic
height, exception, anti-coset, residue-tail, and packet-exit fields are
retained.
**Next:** annotate a HYP-2614/HYP-2963 support-six sample with
`relation_lattice`, `covolume`, `successive_minima_profile`,
`convex_body_id`, `algebraic_target`, `height_bound`,
`approximation_exponent`, `exceptional_approximants`,
`low_height_wall_class`, `deleted_anti_cosets`, `residue_signed_tail`, and
`diophantine_exit`; prove finite low-height deletion before applying
Minkowski and finite exceptional-approximant control before applying Roth.
Namespace: HYP-3062 / T1144 / LTI-209 / LTT-107.
## Lead codex-2026-06-26-S233: Desargues-median finalization lens

**Status:** SYNTHESIS / proof-interface target for LRC14 finalization
(HYP-3067/T1149/LTI-214/LTT-112).
**Readout:** Added
`04-computation/lrc14_desargues_median_lens_codex_s233.py`,
`05-knowledge/results/lrc14_desargues_median_lens_codex_s233.out`,
`05-knowledge/hypotheses/HYP-3067-lrc14-desargues-median-finalization-lens.md`,
`07-reflections/lrc14-desargues-median-finalization-codex-s233.md`, and
three forum posts under `poke-forum/posts/20260626-2043*Z...`.  The audit
makes the Desargues warning exact: cubic, bipartite, girth `6`, diameter `5`,
theta-class sketch `5` classes of size `6`, but `median=False` with `160`
empty-center route triples.  Q4 and the `4x4` grid pass the same median
triple test.
**Next:** build the HYP-2963 medianization table for coarse fibers: route
triple, sidecars attached, median-center status, first missing sidecar, and
repair/debt.  Empty centers should be classified as Desargues defects repaired
by endpoint-owner strips, primitive decks, observer-cut orbits, value-origin
types, Haar/rectangle residues, AP/GW boundary stops, or THM-572/F7 debt.
Namespace: HYP-3067 / T1149 / LTI-214 / LTT-112.

## Lead codex-2026-06-26-S234: Median owner/root sidecar spine

**Status:** SYNTHESIS / proof-interface checklist for LRC14 finalization
(HYP-3068/T1150/LTI-215/LTT-113).
**Readout:** Added
`04-computation/lrc14_median_owner_root_sidecar_audit_codex_s234.py`,
`05-knowledge/results/lrc14_median_owner_root_sidecar_audit_codex_s234.out`,
`05-knowledge/hypotheses/HYP-3068-lrc14-median-owner-root-sidecar-spine.md`,
and `07-reflections/lrc14-median-owner-root-sidecar-spine-codex-s234.md`.
The audit turns the S233 Desargues-median warning into a concrete
owner/root-sidecar table.  C6 has an empty route-triple center while Q3 has a
unique sidecar center.  Six route fibers move from `empty` or `multiple` to
`unique` after first missing sidecars: endpoint owner, rootless cycle object,
Desargues/Beal owner residue, value-origin type, observer-cut orbit, and active
owner/barcode support.  The sidecar tournament is transitive over proof
obligations, not runners.
**Next:** run the table over actual HYP-2963 coarse fibers with fields
`coarse_fiber_id`, `route_triple`, `coarse_shadow`, `root_object`,
`owner_object`, `median_center_status`, `first_missing_sidecar`, and
`repair_or_debt`.  Treat empty centers as first-missing-sidecar defects and
multiple centers as value-origin/vocabulary ambiguity before naming new debt.
Namespace: HYP-3068 / T1150 / LTI-215 / LTT-113.

## Lead codex-2026-06-26-S231: Partial-cube bridge-rank split ledger

**Status:** EXACT AUDIT / addendum to the HYP-3063 proof-interface sidecar.
**Readout:** Added
`04-computation/lrc14_partial_cube_bridge_rank_split_codex_s231.py`,
`05-knowledge/results/lrc14_partial_cube_bridge_rank_split_codex_s231.out`,
`07-reflections/lrc14-partial-cube-bridge-rank-split-ledger-codex-s231.md`,
and `poke-forum/posts/20260626-partial-cube-bridge-rank-split-ledger/post.md`.
Rebased over the upstream HYP-3063/T1145/LTI-210/LTT-108 carrier and S228
exact scout.  This lead checks the remaining bridge/split count layer: two
Moser lanes split every `x<4^m` as `x=a+2b`; and
`2,6,12,20,30,42=k(k+1)` is directed simplex-edge incidence and
`K_{k,k+1}` bridge-line count with rectangle-cycle debt `k(k-1)` beyond rank
`2k`.
**Next:** Attach partial-cube fields to HYP-2963 Moser/fibbinary automatic
fibers and test whether they split route/status mixing before exact magnitude
is restored.  Namespace remains HYP-3063 / T1145 / LTI-210 / LTT-108.

## Lead codex-2026-06-26-S220: Observer-cut orbit ledger

**Status:** SYNTHESIS / formal addendum to the S218 observer-extension cut
payload calculus (HYP-3056/T1138/LTI-203/LTT-101).
**Readout:** Added
`05-knowledge/hypotheses/HYP-3056-observer-cut-orbit-ledger.md`,
`07-reflections/observer-cut-orbit-ledger-controlled-forgetting-codex-s220.md`,
and `poke-forum/posts/20260626-observer-cut-orbit-ledger/post.md`.  The
refinement is that the payload should be an orbit under the visible-fiber
automorphism group:
`C_q(x,o)=orbit_Aut_q(x)(boundary slice, incidence word, extended shadow)`.
This turns controlled forgetting into a ledger of quotient, fiber, next
observer, automorphism group, payload orbit, changed LRC predicate,
separating sidecar, discharge mode, and residual debt name.
**Next:** implement the observer-cut ledger over HYP-2963 coarse fibers,
enumerating admissible observers for mixed route/status pairs.  Emit the
payload-column tournament; a directed cycle should be treated as evidence of
noncommuting discharges and routed to a bicomplex/fiber-product proof carrier.
Namespace: HYP-3056 / T1138 / LTI-203 / LTT-101.

## Lead codex-2026-06-26-S218: Observer-extension cut payload calculus

**Status:** SYNTHESIS / abstraction across existing LRC carriers plus exact
first-failure ledger (HYP-3054/T1136/LTI-201/LTT-099), extending HYP-3050's
observer-cut payload, HYP-3051-HYP-3053's rooted/diagonal transport carriers,
and HYP-3039/HYP-3040's controlled-forgetting ledgers.
**Readout:** Added
`05-knowledge/hypotheses/HYP-3054-observer-extension-cut-calculus.md` and
`07-reflections/observer-extension-cut-controlled-forgetting-codex-s218.md`;
added `04-computation/observer_extension_payload_codex_s218.py`,
`05-knowledge/results/observer_extension_payload_codex_s218.out`, and
`07-reflections/observer-extension-payload-codex-s218.md` as the exact ledger.
The key abstraction: a quotient is proof-safe only after the payload for its
next outside operation is retained, reconstructed, dual/cut/cocycle-killed,
descended, boundary-stopped, or routed to named debt.  The exact tournament
ledger corrects the first shifted A000568 failure to `48+8=56`, not
`48+12=56`; `12` is a fold/parent/fixed-locus count (`R(4)`, `U(5)`, both
`5->6` source/sink deletion slices, and `SC(6)`), while the defect `8` is the
observer-extension/cut payload.  This makes pair-good decoy generator teeth,
residual capacitor first cuts, q13 AP-tail clocks, automaton
magnitude/topology handoffs, rectangle/hourglass defects, deletion-parent
fibers, edge-sector cross orientation, endpoint-owner payloads, and matrix
observability columns instances of one observer-extension/cut rule.
**Next:** add manifest fields `quotient`, `next_operation`,
`observer_extension_payload`, `observer_payload_stage`, `incident_word_orbit`,
`edge_sector_cross_orientation`, `deletion_parent_fiber`,
`rectangle_hourglass_residue`, `endpoint_owner_payload`,
`barcode_active_owner_support`, `extension_address`, `cut_or_cycle_defect`,
`route_owner_certificate`, and `payload_exit/discharge_rule`; build HYP-3048
sidecar observability matrices over HYP-2963 coarse fibers; classify pair-good
decoys by generator and active-owner relation before reporting raw counts.
Namespace: HYP-3054 / T1136 / LTI-201 / LTT-099.
## Lead codex-2026-06-26-S223: Observer-extension cut payload

**Status:** SYNTHESIS / exact finite audit and proof-carrier abstraction
(HYP-3059/T1141/LTI-206/LTT-104), refining incoming HYP-3056, HYP-3055, and HYP-3054 and built on
S213 ordered-pair sector decks, S216 source/sink diagonal transport slices,
and S217 rectangle/hourglass flow.
**Readout:** `04-computation/observer_extension_cut_payload_codex_s223.py`
stores output in
`05-knowledge/results/observer_extension_cut_payload_codex_s223.out`.  The
prompt arithmetic is corrected: `48+12=60`, while the exact first-defect
splice is `U(6)=P(5)+U(5)-U(4)=48+12-4=56`; the defect is `8`, and the
fractions are `1=6/7+3/14-1/14`.  Recurring `12` carriers are `P(4)`,
`U(5)`, source and sink `5->6` slices, `SC(6)`, and S217's `K_{4,5}`
rectangle redundancy.  Source/sink child sets overlap in `4`; `SC(6)`
deletion decks touch all `12` five-parent classes; ordered-pair
cross-sector decks separate `56/56`.
**Next:** materialize the observer-extension/cut payload schema in packet
experiments: `observer_cut_word`, `stabilizer_orbit_id`,
`source_sink_slice_id`, `source_sink_overlap_class`, `deletion_fiber_payload`,
`self_converse_branch_bit`, `cross_sector_orientation_word`,
`rectangle_hourglass_residue`, `endpoint_owner_packet`, and `legality_exit`.
Prove each forgotten coordinate is retained, reconstructed, annihilated,
descended, AP/GW equality, or named residual debt.
Namespace: HYP-3059 / T1141 / LTI-206 / LTT-104.
`observer_extension_payload`, and `discharge_rule`; build HYP-3048 sidecar
observability matrices over HYP-2963 coarse fibers; classify pair-good decoys
by generator and active-owner relation before reporting raw counts.
## Lead codex-2026-06-26-S230: Observer-extension exact duodecimal audit

**Status:** SYNTHESIS / exact small tournament audit and LRC carrier
abstraction (HYP-3065/T1147/LTI-212/LTT-110), refining HYP-3054/HYP-3055 and
the incoming HYP-3056 observer-cut orbit / HYP-3057 value-origin ledgers,
rebased beside the HYP-3058 hyperbolic reciprocal sidecar, and built on the
S211..S217 A000568/rooted-perspective/ordered-pair/diagonal-flow stack plus the
HYP-3039 controlled-forgetting ledger.
**Readout:** `04-computation/observer_extension_exact_duodecimal_audit_codex_s230.py`
stores output in
`05-knowledge/results/observer_extension_exact_duodecimal_audit_codex_s230.out`.
The tempting statement `48+12=56` is arithmetically false; the useful ledger
is `U(6)=P(5)+SC(6)-U(4)=48+12-4=56`, with defect `8=(2/3)*12` and overlap
kernel `U(4)=4=(1/3)*12`.  The recurring dozen is `P(4)=12`, `U(5)=12`,
`SC(6)=12`, and the S217 `K_{4,5}` rectangle redundancy.  Exact sector audit
again finds sector size/internal decks separate `55/56`, while
`cross_sector_orientation_word` separates `56/56`.
**Next:** build or refute the actual inclusion-exclusion object behind the
overlap law; compare the `U(4)` kernel with deletion boundaries,
self-converse branch loci, the `344/345` ordered-pair-sector collision, and
rectangle/hourglass cycle residues.  Add
`observer_extension_cut_signature`, `value_origin_type`,
`observer_cut_payload_orbit`, `duodecimal_overlap_kernel`,
`self_converse_branch_locus`, `cross_sector_orientation_word`,
`deletion_parent_profile`, `rectangle_cycle_residue`, and
`hourglass_cycle_residue` to packet experiments.
Namespace: HYP-3065 / T1147 / LTI-212 / LTT-110.

## Lead codex-2026-06-26-S217: Tournament diagonal-layer flow

**Status:** SYNTHESIS / exact small tournament audit and algebraic carrier law
(HYP-3053/T1135/LTI-200/LTT-098), rebased over the incoming S212 expanded
matrix atlas, S213 A000568 edge-perspective lift, and S214 exact
non-node perspective-depth carrier, plus S215 rooted layer-extension flow, and
S216 diagonal-layer transport, and linked to their sidecar-observability,
cross-sector, edge-triple, rooted-fiber, and deletion-transport targets.
**Readout:** `04-computation/tournament_diagonal_layer_flow_codex_s217.py`
stores output in
`05-knowledge/results/tournament_diagonal_layer_flow_codex_s217.out`.  The
user's `k^2+k` inter-layer lines are a `GF(2)` coboundary/cut-space
presentation on `K_{k,k+1}`: `k(k+1)` lines, rank `2k`, and `k(k-1)`
rectangle redundancies, with spanning-tree reconstruction.  Full adjacent
layer flow has rank `C(n,2)-1` and redundancy
`2*C(n-1,3)+C(n-2,2)`, decomposing into local rectangles plus hourglass
cycles.  Fixed Hamiltonian-path half-tilings cover A000568 with fiber
`H(T)/|Aut(T)|`; diagonal path-reversal/converse is useful but not A000568.
**Next:** emit rectangle/hourglass cycle bases; attach endpoint-owner,
barcode-bar, active-normal-fan, edge-sector, cross-sector orientation,
cycle-conflict, deletion-parent, and clique insertion sidecars; then feed
those fields into HYP-3052's diagonal transport ledger and HYP-3048's sidecar
observability matrix to test route/status separation.
Namespace: HYP-3053 / T1135 / LTI-200 / LTT-098.

## Lead codex-2026-06-26-S206: LRC14 comprehensive lens map

**Status:** SYNTHESIS / coordination map complete
(HYP-3043/T1124/LTI-191/LTT-089), built to let concurrent agents place new
LRC ideas beside the existing packet, topology, arithmetic, analytic,
automaton, pair-decoy, residual, formal, and external-analogy lenses.
**Readout:** `00-navigation/LRC-LENS-MAP.md` maps each lens by preserved LRC
predicate, destroyed coordinate, required sidecar, and handoff target.  The
shared center is a labelled packet sheaf with controlled forgetting: a quotient
is legal only if the forgotten coordinate is constant, reconstructed,
dual-annihilated, descended, boundary-stopped, or routed to named
F7/THM-572/harmonic debt.  Pair-good decoys are blocker-generator teeth routed
through barcode/normal-fan active owners; Moser/fibbinary/lacunary shadows are
telemetry sidecars, not proof quotients; analytic/operator analogies are
capacity/blindness or certificate backends only after packet translation.
**Next:** add `lens_family`, `preserved_lrc_predicate`,
`destroyed_coordinate`, `required_sidecar`, `handoff_target`,
`status_mixing_result`, `route_mixing_result`, `tournament_vertex_choice`, and
`challenged_assumption` to new LRC lens notes and packet-ledger experiments,
then rerun status/route fiber checks before promoting a new carrier.
**Integration note:** after fetch revealed the concurrent S204 AP-tail
puncture atlas, this map was renumbered to S206/HYP-3043/T1124/LTI-191/LTT-089
and the AP-tail q13/fixed-point clock was added as a concrete controlled-
forgetting repair lens.
Second integration note: the concurrent S205 owner-strip filtration was also
absorbed as the post-status page order between primitive clocks/Haar zeta and
named F7/THM-572 debt.
Namespace: HYP-3043 / T1124 / LTI-191 / LTT-089.

## Lead codex-2026-06-26-S214: Perspective-depth controlled forgetting

**Status:** EVIDENCE / exact first-failure scout complete (HYP-3050/T1132/LTI-197/LTT-095), extending HYP-3047's A000568 k-depth incident-lift ladder, HYP-3049's ordered-pair edge lift, and HYP-3048's matrix-atlas sidecar language to connect the old rooted-perspective/A000568 curiosity to HYP-3039/HYP-3040's controlled-forgetting ladder.
**Readout:** `04-computation/perspective_depth_ladder_codex_s214.py` stores output in `05-knowledge/results/perspective_depth_ladder_codex_s214.out`.  The equality `P(m)=A000568(m+1)` first fails at base size `m=5` / next index `n=6`: exact rooted node perspectives give `P(5)=48` while `A000568(6)=56`.  Directed k-depth node perspectives recover exact rooted node orbits by depth `3` at `m=5`, so the missing `8` are observer-extension/cut states, not deeper node neighborhoods.  At `m=5`, exact arc perspectives and triple perspectives both total `88`; triples split into `64` transitive and `24` cyclic.  Edge depth `2` and triple depth `2` recover those exact first-failure carriers.
**Next:** define observer-extension cut perspectives, extend exact non-node carrier counts to `m=6` so conflict-pair perspectives appear, and compare edge tail/tip sectors with pair-good decoy teeth plus residual capacitor cuts.
Namespace: HYP-3050 / T1132 / LTI-197 / LTT-095.
## Lead codex-2026-06-26-S202: LRC14 hidden-coordinate ledger

**Status:** SYNTHESIS / proof-order ledger complete (HYP-3039/T1120/LTI-187/LTT-085), connecting the latest residual stack to older address-retention, visible/hidden fold, anti-wedge, cochannel, and pair-good decoy lessons.
**Readout:** `05-knowledge/hypotheses/HYP-3039-lrc14-hidden-coordinate-ledger.md` and `07-reflections/lrc14-hidden-coordinate-ledger-codex-s202.md` describe HYP-3024..HYP-3038 as a controlled-forgetting ladder: status gate, owner-essential path lift, residual certificate tooth, first-tooth owner strip, primitive-period deck, residual capacitor cut, and q=23 drop/add zeta plus endpoint-owner strip.  The main sharpening is that residual counts are secondary once the generator/cut class is known; a quotient is legal only after the next hidden coordinate has been exposed, killed by a dual/cut/cocycle, or routed to named debt.
**Next:** add `hidden_coordinate_stage`, `visible_hidden_relation_type`, `primitive_safe_deck_2_13`, `first_primitive_safe_q_2_13`, `residual_capacitor_id`, `first_cut_stage`, `drop_add_square_id`, `exact_M_zeta`, `endpoint_owner_strip`, and `anti_wedge_debt_count` to a cached HYP-2963 ledger; then run an anti-wedge audit over accepted residual cuts.
Namespace: HYP-3039 / T1120 / LTI-187 / LTT-085.

## Lead codex-2026-06-26-S203: LRC14 hidden statement ledger

**Status:** SYNTHESIS / proof-facing micro-statement ledger complete
(HYP-3040/T1121/LTI-188/LTT-086), synthesizing HYP-3034..HYP-3038 with older
pair-good decoy, automaton/Moser/fibbinary, barcode/normal-fan, Ramanujan, and
analytic-clock work.
**Readout:** `04-computation/lrc14_hidden_statement_ledger_codex_s203.py`
stores output in
`05-knowledge/results/lrc14_hidden_statement_ledger_codex_s203.out`.  The
ledger records `11` small statements and their lost-coordinate guardrails.  The
fundamental readout is layered obstruction calculus: owner-essential boundary
H1, status-preserving coarse gate, owner/topology/stalk first teeth,
`q<=13` primitive witness scheduler, and residual capacitor/Haar-zeta cuts.
The strongest hidden connection is HYP-3035 residual first-tooth atlas to
HYP-3029/HYP-3018 safe-component stalk/barcode/normal-fan geometry
(`weight=11`).
**Next:** add the hidden sidecar vocabulary to packet ledgers:
`boundary_h1_owner_support`, `first_tooth`, `primitive_safe_deck_2_13`,
`residual_capacitor_id`, `first_cut_stage`, `exact_M_zeta`,
`endpoint_owner_strip`, `analytic_blindness_report`, and
`automaton_shadow_class`; then prove the candidate principle that every
primitive packet is AP/GW boundary H1, a `q<=13` primitive witness,
an owner/stalk/topology descent, or a named capacitor/zeta cut to existing
route debt.
Namespace: HYP-3040 / T1121 / LTI-188 / LTT-086.
## Lead codex-2026-06-26-S209: LRC14 hidden connection accelerators

**Status:** SYNTHESIS / source-audited proof-carrier index complete
(HYP-3046/T1127/LTI-194/LTT-092), rebased over S202/HYP-3039 hidden-coordinate
S203/HYP-3040 hidden-statement ledger work, S204/HYP-3041 AP-tail puncture
work, S205/HYP-3042 owner-strip filtration, S206/HYP-3043 lens map, and
S207/HYP-3044 topology-exception teeth plus incoming S208/HYP-3045
endpoint-owner transfer while using S201/HYP-3038 q=23 drop/add square,
S204's q13 AP-tail repair, S205's endpoint-owner filtration, S207's
owner-stalk topology teeth, and S208's `B18Z6` owner-current address lift as
concrete local normal forms.  It
connects back to HYP-2996/HYP-2992 residual section/Haar exits, HYP-2886
exact-period packets, HYP-3042 endpoint-owner filtration, HYP-3045
endpoint-owner transfer, HYP-3044 topology
exception fields, HYP-3018 normal-fan support, HYP-3013 perfect/divisor lanes,
HYP-3022 decoy blocker decks, and HYP-2995/HYP-3006 cocycle exactness.
**Readout:** `04-computation/lrc14_hidden_connection_accelerators_codex_s209.py`
stores output in
`05-knowledge/results/lrc14_hidden_connection_accelerators_codex_s209.out`.
The script verifies `106` local source markers with `0` missing and ranks twelve
connector lemmas.  Top accelerator: HYP-3037 capacitor exits are not new proof
vocabulary; they are the HYP-2996/HYP-2992 exit alphabet.  HYP-3045 supplies
the `B18Z6` endpoint-owner transfer address lift, HYP-3038 supplies the
concrete q=23 nested-refinement square, HYP-3036 primitive decks are
exact-period atlases, owner strips route through HYP-3042 endpoint-owner
filtration plus HYP-3045 endpoint-owner transfer and HYP-3018/HYP-3034 owner
support, HYP-3041 supplies the q13
AP-tail prime-clock repair, HYP-3044 routes compact-topology failures through
owner-stalk and primitive-deck teeth, q=14 routes through THM-523/HYP-2917,
pair-good decoys route through blocker decks, and squarefree blindness routes
through divisor-lattice prime-power fields.
**Next:** add the HYP-3046 sidecar merge set to HYP-2963 packet records:
`zeta_exit_class`, `residual_section_exit`, `primitive_safe_deck_2_13`,
`first_primitive_safe_q_2_13`, `coarse_endpoint_word`,
`external_endpoint_owner_strip`, `endpoint_owner_transfer_delta`,
`endpoint_owner_residue_delta`, AP-tail puncture/fixed-point fields,
owner-transfer and owner-support fields, `first_surviving_filtration_page`, `first_cut_stage`,
topology-exception fields, `drop_add_square_id`, `omega_Q_class`,
exact-period chart fields,
`divisibility_threshold_qS`,
divisor-lattice fields, and blocker-deck fields.
Namespace: HYP-3046 / T1127 / LTI-194 / LTT-092.

## Lead codex-2026-06-26-S199: LRC14 residual tooth atlas

**Status:** EVIDENCE / residual-tooth proof-interface atlas complete
(HYP-3035/T1116/LTI-183/LTT-081), extending the S192/S194 residual status
gate and S195 Haar-tile repair dictionary; sibling to the incoming
S198/HYP-3034 arc-boundary path lift, S197/HYP-3033 residual
certificate-teeth scheduler, and S196/HYP-3032 analytic sieve-clock bridge.
**Readout:** `04-computation/lrc14_residual_tooth_atlas_codex_s199.py`
stores output in
`05-knowledge/results/lrc14_residual_tooth_atlas_codex_s199.out`.  It parses
S194's stored `15` coarse ET+unit route-mixed residual fibers, recomputes only
the `38` selected HYP-2963 packets, and attaches arc topology, coarse/exact
safe-component stalks, magnitude cocycle, and q/covering certificate teeth.
All rows are strict-open.  First tooth: `arc_topology_compact` handles `13`
fibers; `coarse_safe_stalk` handles the other `2`; all first repairs are
`owner_strip`.
**Next:** add `first_tooth` and `residual_tooth_class` to HYP-2963 sidecars,
then prove the arc-topology and coarse-stalk owner-strip descent lemmas.
Namespace: HYP-3035 / T1116 / LTI-183 / LTT-081.
## Lead codex-2026-06-26-S200: LRC14 Ramanujan primitive-period route scheduler

**Status:** EVIDENCE / residual route scheduler complete (HYP-3036/T1117/LTI-184/LTT-082), downstream of HYP-3030's status-topology gate, HYP-3031's repair-ladder synthesis, and the HYP-3033 residual-certificate teeth.
**Readout:** `04-computation/lrc14_ramanujan_route_scheduler_codex_s200.py` stores output in `05-knowledge/results/lrc14_ramanujan_route_scheduler_codex_s200.out`.  The script parses the S194 stored residual list and rebuilds only the `38` packets in the `15` route-mixed coarse ET+unit fibers.  Baseline residuals have `15` mixed-route fibers, `0` mixed-status fibers, and routes `Q-WITNESS=23`, `COVERING-MOMENT=15`.  Adding `first_primitive_safe_q_2_13` or `primitive_safe_deck_2_13` gives `30` fibers with `0` mixed route and `0` mixed status.  All residual covering rows have zero primitive safe mass for `q<=13`; q=14 mass is a separate boundary/covering layer.
**Next:** add `primitive_safe_deck_2_13` and `first_primitive_safe_q_2_13` to the HYP-3027/HYP-3031 packet sidecar, run a cached full-bank ledger, and prove familywise that post-status positive q<=13 primitive mass is exactly the residual Q-WITNESS scheduler while zero-deck open rows route to covering/q=14/boundary-moment certificates.
Namespace: HYP-3036 / T1117 / LTI-184 / LTT-082 / CPI-043.

## Lead codex-2026-06-26-S191: LRC carrier pullback mega-index

**Status:** SYNTHESIS / contribution index complete (T1108/LTI-175/LTT-073/LTM-079).  User asked to get more creative by looking back through tournament, related-series, metagraph, automaton, topology, arithmetic, harmonic, and formalization work and bringing as many techniques as possible back to LRC.
**Readout:** `00-navigation/LRC-CARRIER-PULLBACK-INDEX.md` gives `90` `CPI-*` carrier rows.  Each row records source handles, the pulled-back proof carrier, retained LRC packet fields, quotient guardrail, and a next-agent pull.  Tournament Analysis vertices are proof-carrier pullbacks, not runners; the observable records retained boundary/open status, theorem route, exact scale, endpoint/topology, arithmetic period, harmonic certificate, residual routing, family transfer, formal payload, and proof cost.
**Next:** instantiate rows rather than citing them.  Pick a `CPI-*` row from each of the six bundles, attach its fields to HYP-2963 packets or a named stress family, run boundary/open and route-fiber mixing, then promote rows with zero unsafe mixing into theorem-facing LTI/LTT cards.  Rows that destroy a coordinate must name the reconstruction, dual certificate, family descent, exact cocycle, or residual debt.
Namespace: T1108 / LTI-175 / LTT-073 / LTM-079 / CPI-001..CPI-090.
## Lead codex-2026-06-26-S193: LRC14 safe-component stalk descent

**Status:** EVIDENCE / target-fiber local carrier scout complete (HYP-3029/T1110/LTI-177/LTT-075), downstream of HYP-3026 carrier fusion, HYP-3025 arc-Cech, HYP-3024 fiber-zipper convergence, HYP-3023 automatic zipper, HYP-3018 normal fan, and HYP-3015 barcode work.  User asked for more creative LRC angles.
**Readout:** `04-computation/lrc14_safe_component_stalk_codex_s193.py` stores output in `05-knowledge/results/lrc14_safe_component_stalk_codex_s193.out`.  On the target automatic word `MFCMMCCFFFCCC`, `639` HYP-2963 packets have `27` mixed residue-terminal route fibers with max mixed `30`.  Owner-only largest-component stalks reduce this to `7` mixed fibers/max `5`; coarse stalks reduce to `2` size-2 open-route collisions; exact largest-component stalks have `0` mixed route fibers, matching exact magnitude without carrying the full fusion packet.
**Next:** prove largest-stalk descent inside the target word; prove the two coarse residual families `13->159/117` and `13->118/104`; run the stalk key over the full HYP-2963 bank and compare against HYP-3025 Cech facets, HYP-3018 normal-fan supports, HYP-3015 barcode fields, and HYP-3026 fusion sidecars.
Namespace: HYP-3029 / T1110 / LTI-177 / LTT-075.
## Lead codex-2026-06-26-S194: LRC14 status-topology gate for coarse zipper fibers

**Status:** EVIDENCE / full-bank bridge complete (HYP-3030/T1111/LTI-178/LTT-076), splicing HYP-3024's coarse ET+Henselian-unit status convergence with HYP-3025's closed arc-Cech topology and feeding the HYP-3026/HYP-3028/HYP-3029 switchboard stack.
**Readout:** `04-computation/lrc14_status_topology_gate_codex_s194.py` stores output in `05-knowledge/results/lrc14_status_topology_gate_codex_s194.out`.  On the default HYP-2963 `21913`-packet bank, residue-terminal fibers have exactly `2` mixed boundary/open fibers; the only boundary rows inside them are AP and GW, each with closed arc beta `(1,1)`, open arc beta `(6,0)`, zero safe topes, and six owner sums `0 mod 14`.  Every open cohabitant has closed arc `beta1=0` and at least `4` safe topes.  The coarse ET+unit gate has `0` mixed status fibers; its `15` route-mixed fibers contain `38` open packets, all closed arc `beta1=0`.
**Next:** prove a topology-first zero-open theorem: primitive zero-open packets must carry the AP/GW arc-boundary cycle or emit named F7/THM-572 residual debt.  Then use coarse ET+unit only as a status-preserving quotient, with magnitude/barcode/Fejer/K33 machinery scheduling the remaining open route collisions.
Namespace: HYP-3030 / T1111 / LTI-178 / LTT-076.

## Lead codex-2026-06-26-S197: LRC14 residual certificate teeth

**Status:** EVIDENCE / stored-ledger residual scheduler complete (HYP-3033/T1114/LTI-181/LTT-079), downstream of HYP-3031's Haar repair dictionary, HYP-3030's status-topology gate, HYP-3028's residual status gate, and HYP-3024's coarse ET+Henselian-unit convergence.
**Readout:** `04-computation/lrc14_residual_certificate_teeth_codex_s197.py` stores output in `05-knowledge/results/lrc14_residual_certificate_teeth_codex_s197.out`.  On the S194 residual ledger, the `15` route-mixed coarse ET+unit fibers contain `38` packets and all are open.  Topology alone leaves `3` mixed route classes, unit-scale alone leaves one mixed class, and exact `M` fallback leaves `2` mixed classes.  Full topology plus unit-scale, and the compressed `(safe topes, quotient defect, open beta0 - closed beta0)` bucket plus unit-scale, both produce `21` residual fibers with `0` route mixing.
**Next:** materialize `residual_topology_bucket`, `unit_scale_tooth`, and `residual_certificate_tooth` in the packet classifier instead of parsing S194 text; rerun on the full HYP-2963 bank; then prove the family scheduler that routes open residuals to q-witness, covering/Haar/nested refinement, or named F7/THM-572 debt.
Namespace: HYP-3033 / T1114 / LTI-181 / LTT-079.

## Lead codex-2026-06-26-S188: LRC14 fiber-zipper convergence via Erdos-Turan and Henselian unit rule

**Status:** EVIDENCE / target-fiber convergence scout complete (HYP-3024/T1104/LTI-171), downstream of HYP-3023 automatic fiber zipper and HYP-3020 discrepancy-height trident.  User asked to extend fiber-zipper convergence using Erdos-Turan and a Henselian unit rule.
**Readout:** `04-computation/lrc14_fiber_zipper_convergence_codex_s188.py` stores output in `05-knowledge/results/lrc14_fiber_zipper_convergence_codex_s188.out`.  The default run filters the HYP-2963 bank to the target automatic word `MFCMMCCFFFCCC` before exact packet computation (`639` rows of `21913`).  Residue-terminal fibers have `27` mixed route fibers and largest mixed size `30`; Erdos-Turan low-frequency residue bins reduce this to `6` and `4`; the Henselian unit rule at `p=2,3,7` gives `0` mixed route fibers before exact magnitude is attached.
**Next:** prove the Henselian unit split for the target word; run the ET/Hensel zipper on the full HYP-2963 bank; classify any full-bank survivors by q-threshold, unit-excess lane, barcode support, and packet route before replacing exact magnitude.
Namespace: HYP-3024 / T1104 / LTI-171.
## Lead codex-2026-06-26-S188: LRC14 fiber-zipper convergence and Henselian unit rule

**Status:** EVIDENCE / full-bank convergence scout and proof-interface target
complete (HYP-3024/T1105/LTI-171/LTT-070), bridging HYP-3023's automatic
fiber zipper and HYP-3020's discrepancy-height-Hensel trident.  User asked to
extend zipper-fiber convergence while thinking Erdos-Turan and Henselian unit
rule.
**Readout:** `04-computation/lrc14_fiber_zipper_convergence_codex_s188.py`
stores output in
`05-knowledge/results/lrc14_fiber_zipper_convergence_codex_s188.out`.  On the
default HYP-2963 `21913`-packet bank, automatic words have `143` mixed-route
fibers with max mixed `1179`; residue-terminal fibers reduce max mixed to
`30` but still have `265` mixed-route fibers.  Exact ET clocks at `14,27,41`
split to singleton fibers, so they are address-like.  The coarser ET+unit gate
is the useful carrier: `21702` fibers, `15` mixed-route fibers, max mixed `4`,
and `0` mixed boundary/open fibers.  Hensel data is refined by separating
unit roots in `F_p^*` from the forced zero root/scale debt.
**Next:** prove coarse ET+unit status convergence inside automatic/residue
fibers; then discharge the remaining route-mixed open fibers by magnitude
formulas, q-witness/covering/petal exits, Fejer/Ramanujan/Haar/barcode
certificates, or K33/F7/THM-572 residuals.
Namespace: HYP-3024 / T1105 / LTI-171 / LTT-070.

## Lead codex-2026-06-26-S186: LRC14 pair-good decoy barcode and normal-fan refinement

**Status:** EVIDENCE / exact bounded classifier and quotient guardrail complete (HYP-3022/T1103/LTI-169), downstream of HYP-3021 modular tooth classification, HYP-3019 binding-pair switches, and HYP-3018 active normal fans.  User asked to classify pair-good decoys and relate them to barcode bars and active bottleneck owners.
**Readout:** `04-computation/lrc14_pair_good_decoy_barcode_normal_fan_codex_s186.py` stores output in `05-knowledge/results/lrc14_pair_good_decoy_barcode_normal_fan_codex_s186.out`.  Pair-good decoys are exact pair-crossing times where a pair is threshold-good but the full row is blocked below `1/14`.  Named rows have `8889` decoys and `5692` generator classes; low-frontier one-swap atlas `drops={10,12,13}, add<=36` has `48037` decoys and `9809` generator classes.  All named decoys are outside strict barcode bars or in zero-bar rows; many overlap true normal-fan supports through the good pair or blocker deck, while disjoint classes are explicit lane/shell/blocker false switches.
**Next:** add pair-good generator fields to HYP-2963 packets next to HYP-3021 tooth fields; prove blocker-deck lemmas for singleton blockers `(7,)`, `(11,)`, `(1,)`, `(9,)`, `(13,)`, `(8,)`; route overlap classes by HYP-3018 normal-fan support and disjoint classes by denominator lane/blocker clearance.
Namespace: HYP-3022 / T1103 / LTI-169.

## Lead codex-2026-06-26-S182: LRC14 active-bottleneck normal-fan carrier

**Status:** EVIDENCE / exact finite scout and proof-interface carrier complete (HYP-3018/T1101/LTI-166), extending the HYP-3015 barcode and using HYP-3016/HYP-3017 mixed automaton fibers as the negative control.  User asked to keep looking for better proof carriers.
**Readout:** `04-computation/lrc14_active_bottleneck_normal_fan_codex_s182.py` stores output in `05-knowledge/results/lrc14_active_bottleneck_normal_fan_codex_s182.out`.  The sidecar records left endpoint owners, peak bottleneck owners, right endpoint owners, support speeds, and residue sums for each strict bar.  AP/GW remain zero-bar boundary atoms; K33 emits `(5,36)`, petals emit `(7,20)` and `(1,26)`, covering low-persistence bars emit `(13,84)`, and mixed AP/GW automaton-fiber open rows emit supports such as `(5,7)` and `(5,96)`.  Aggregate tested rows: `138` positive bars, peak-owner histogram `{2:114,3:18,4:6}`, support-size histogram `{2:86,3:42,4:8,5:2}`.
**Next:** add normal-fan sidecars to the HYP-2963 packet bank; test whether the HYP-3017 mixed automatic-word fibers become boundary/open pure after adding peak support and endpoint residue sums; route the lowest-persistence supports to Fejer interval and endpoint-owner certificates.
Namespace: HYP-3018 / T1101 / LTI-166.
## Lead codex-2026-06-26-S184: LRC14 discrepancy-height trident carrier

**Status:** EVIDENCE / bounded carrier scout and proof-interface target complete
(HYP-3020/T1102/LTI-167/LTT-067), extending the HYP-3016 magnitude-cocycle
guardrail and HYP-3017 route-purity audit.  User asked to keep looking for
better proof carriers.  This pass treats automaton fiber-mixing as a request
for a product carrier rather than another finite-state word.
**Readout:** `04-computation/lrc14_discrepancy_height_carrier_codex_s184.py`
stores output in
`05-knowledge/results/lrc14_discrepancy_height_carrier_codex_s184.out`.
Exact threshold-`1/14` safe components on `2173` named plus AP single-swap rows
show AP/GW are the only boundary atoms and K33 `12->36` is the minimum positive
row at `1/1260`.  Automaton words, residue+MFC pairs, residue discrepancy
alone, Hensel alone, and height alone still mix boundary/open rows; the full
discrepancy-height-Hensel trident has `0` mixed boundary/open fibers, largest
fiber `2`, and `6` mixed route fibers.
**Next:** add trident fields to the full HYP-2963 classifier, then coarsen
residue denominators, Erdos-Turan bins, height/Mahler buckets, and Hensel flags
until the smallest route-pure signature is found for large mixed automatic
fibers such as `MFCMMCCFFFCCC`.  Compare the surviving fibers with HYP-3015
barcodes and HYP-2981 Fejer interval manifests.
Namespace: HYP-3020 / T1102 / LTI-167 / LTT-067.

## Lead codex-2026-06-26-S179: LRC14 lonely-profile persistence barcode carrier

**Status:** EVIDENCE / exact finite scout and proof-interface carrier complete (HYP-3015/T1099/LTI-164), rebased over the incoming creative exact packet-lens atlas (HYP-3014/LTI-163).  User asked for more creative LRC angles.  This pass treats the lonely profile `m_S(t)=min_v ||v t||` as a threshold persistence object rather than a raw maximin scalar.
**Readout:** `04-computation/lrc14_lonely_profile_persistence_barcode_codex_s179.py` stores output in `05-knowledge/results/lrc14_lonely_profile_persistence_barcode_codex_s179.out`.  Exact Fraction arithmetic over affine cells shows AP13 and GW `12->24` have zero bars; K33 has two bars, petals have two, covering `12->84` has eight, fibbinary first13 has `38`, and Moser first13 has `64`.  The barcode exposes bar count, length, peak, and persistence margin, so raw `M` and safe mass are insufficient quotients.
**Next:** run barcode sidecar over the HYP-2963 packet bank; route lowest-persistence positive families to endpoint-owner / Fejer interval certificates; compare barcode classes against HYP-3012 automaton classes and HYP-3013 divisor fields.
Namespace: HYP-3015 / T1099 / LTI-164.
## Lead codex-2026-06-26-S174: LRC14 closed arc-Cech nerve carrier

**Status:** SYNTHESIS / exact topology proof-interface carrier complete
(HYP-3025/T1104/LTI-172/LTT-069), downstream of the endpoint
tope/cocircuit, taut-current, Fejer/Ramanujan, and automaton/gap-sidecar
threads.  User asked to keep looking for better creative proof carriers.
**Readout:** `04-computation/lrc14_arc_cech_nerve_carrier_codex_s174.py`
stores output in
`05-knowledge/results/lrc14_arc_cech_nerve_carrier_codex_s174.out`.  The
script builds the exact threshold danger-arc endpoint arrangement, computes
the open and endpoint-completed closed arc-Cech nerves, compares them to the
runner quotient nerve, and records quotient Betti defect.  AP and GW
`12->24` have closed arc Betti `(1,1)`, open arc Betti `(6,0)`, six
boundary-safe cocircuits, and boundary owner sums all `0 mod 14`; named
positive rows have closed arc `beta1=0`.  The AP one-swap scan through
`add<=160` has exactly one zero-open row, GW, and smallest positive safe mass
`1/1260` at K33 `12->36`.
**Next:** add `closed_arc_cech_beta`, `open_arc_component_count`,
`boundary_cocircuit_facet_word`, `boundary_owner_sum_word_mod_14`,
`runner_quotient_betti_defect`, `private_arc_count`, `private_runner_count`,
`safe_tope_count`, and `arc_cech_exit_route` to HYP-2963 packet records or a
sidecar.  Run the full HYP-2963 bank through the closed arc-Cech audit.  Prove
familywise that K33/petal/covering positive rows have `beta1=0` exits, and
define F7 as a good-cover quotient-defect sector if any non-AP/GW zero-open
packet survives.

## Lead codex-2026-06-25-S178: LRC14 Fermat-Catalan automatic-gap extension

**Status:** SYNTHESIS / packet-schema guardrail complete
(HYP-3009/T1093/LTI-159/LTT-061), extending upstream HYP-3008/LTI-158.  User
asked to keep pushing toward LRC14
while thinking Fermat-Catalan, the 2-adic Littlewood/Hurwitz doubling paper,
Ostrowski-Hadamard gap theorem, Moser-de Bruijn sequence, fibbinary numbers,
and finite automata.
**Readout:** `04-computation/lrc14_fermat_catalan_automatic_gap_s178.py`
stores output in
`05-knowledge/results/lrc14_fermat_catalan_automatic_gap_s178.out`.  Named-row
audit labels speeds by Moser even-bit (`M`), fibbinary/no-adjacent carry (`F`),
and carry-present (`C`) status.  AP13 and GW `12->24` share the same word
`MFCMMCCFFFCCC`, proving the automatic shadow alone cannot distinguish equality
atoms.  K33 `12->36`, C27 petals, covering tails, and a Res_27 probe split by
automatic counts, lacunary tail-gap ratios, and perfect-power flags.  Unit
excess payloads flag `p=2: q=27=3^3` plus powers at `p=4,8` as no-lift guards
requiring cyclotomic/p-adic labels.  Incoming HYP-2937/HYP-2944/HYP-2950
checkpoint work supplies the route families this ledger should tag:
C27/unital transfers, Farey-product perfect-number gates, and
Borel-Baire-Haar witness labels.
**Next:** add `automatic_language_class`, `fibbinary_carry_status`,
`moser_even_bit_status`, `ostrowski_digit_system`, `lacunary_gap_ratio`,
`power_lift_guard`, `fermat_catalan_residual`, `hurwitz_doubling_cf_state`,
and `visibility_potato_approx_guard` to HYP-2963 packet records.  Compute
route-purity of automatic words over the full labelled packet bank and list the
first mixed fibers.  Namespace: HYP-3009 / T1093 / LTI-159 / LTT-061.
## Lead codex-2026-06-25-S173: gap automaton carrier for LRC14

**Status:** SYNTHESIS / finite automaton audit and packet-schema guardrail
complete (HYP-3012/T1096/LTI-161/LTT-063), extending the incoming HYP-3008
automatic-gap scout, HYP-3009 power-lift ledger, HYP-3010 exact-gap/safe-component carrier, and HYP-3011 automatic lacunary safe-component filter.  User asked to keep pushing toward
an LRC14 proof while thinking about Fermat-Catalan, 2-adic Littlewood,
Ostrowski-Hadamard gaps, Moser-de Bruijn, fibbinary numbers, and finite
automata.
**Readout:** `04-computation/lrc14_gap_automaton_carrier_codex_s173.py`
stores output in
`05-knowledge/results/lrc14_gap_automaton_carrier_codex_s173.out`.  Through
`N=4096`, fibbinary has `378` members, Moser has `65`, their intersection is
`65`, all `14` residue classes are mixed for both languages, fibbinary closes
under `x->2x`, Moser closes under `x->4x`, and Moser has `63` violations under
`x->2x`.  On the unit-excess lane `q=14p-1`, there are `21` fibbinary/Moser
hits for `p<=384`.  Priority-gauge Tournament Analysis is transitive, while
fieldwise majority has SCC `{fibbinary_no_adjacent_language,
moser_base4_digit_language, ostrowski_hadamard_lacunary_boundary}` and induced
isomorphism-class counts `n=4:3`, `n=5:4`, `n=6:4`.
**Next:** add S173 automaton/gap fields to HYP-2963 packet records or a sidecar
manifest, build product automata for `q=14p-1`, and prove the finite-state
residual dichotomy: finite-state hard packets route to known
q/Fejer/Ramanujan/Haar/state-lift certificates, otherwise the first failed
transition names a nonregular/natural-boundary F7 sector.
Namespace: HYP-3012 / T1096 / LTI-161 / LTT-063; base threads HYP-3011 / T1095 / LTI-160, HYP-3010 / T1094, HYP-3009 / LTI-159
and HYP-3008 / LTI-158.

## Lead codex-2026-06-25-S172: Poincare worldline symmetry ledger for LRC14

**Status:** SYNTHESIS / exact symmetry stress test and packet-schema guardrail complete (HYP-3007/T1091/LTI-157/LTT-060).  User asked to push toward an LRC14 proof and think of the Poincare group.  The pass models LRC rows as time/phase cylinder worldlines with observer danger tubes and separates exact anchored-LRC symmetries from observer-coupled boosts and metric-deforming Lorentz shadows.
**Readout:** `04-computation/lrc14_poincare_symmetry_ledger_codex_s172.py` stores exact safe-measure tests in `05-knowledge/results/lrc14_poincare_symmetry_ledger_codex_s172.out`.  Sign flips, reflection/time reversal, and integer dilation preserve safe measure for AP13, GW13, K33 `12->36`, and a q27 two-block probe; stationary speed translation `+5` changes safe measure in every row tested, while observer-coupled boost/recenter preserves it.  Tournament Analysis over symmetry/proof carriers is transitive with path `individual_sign_flip_parity_kernel > observer_coupled_worldline_tube_groupoid > integer_time_dilation_primitive_scale > anchored_metric_winding_tournament > stationary_velocity_translation > bare_winding_iso_class > lorentz_velocity_addition_shadow > raw_speed_scalar`.
**Next:** add `observer_velocity_label`, `relative_speed_normal_form`, `sign_kernel_status`, `primitive_scale_gcd`, `tube_metric_label`, `worldline_frame_label`, `boost_cocycle_status`, and `orientation_debt_for_winding` to HYP-2963 packet records.  Prove the sign-kernel lemma and boost-admissibility lemma before using Poincare/Lorentz language in any proof route.
Namespace: HYP-3007 / T1091 / LTI-157 / LTT-060.

## Lead codex-2026-06-24-S170: LRC14 curried packet functional tower

**Status:** SYNTHESIS / proof-interface guardrail complete
(HYP-3002/T1086).  User asked to integrate the Fibonacci/Farey/additive-basis
work with LRC by thinking of everything as functions and currying them.
**Readout:** replace `row -> scalar` with
`E:S -> P(S) -> root -> lane -> fiber -> certificate -> verdict`.  A quotient
is a partial evaluation plus lost-coordinate function
`lost_Q(c)(y)(x,x')=c(x)-c(x')`, which must be zero, reconstructed,
exact/coboundary, dual-annihilated, descended, AP/GW boundary equality, or
named F7/THM-572 debt.  The scout prints path-rank rows, unit-excess Farey
partials, named-row call signatures, and a transitive Tournament Analysis over
function carriers.
**Next:** add `curried_call_signature` and `lost_coordinate_function` to
HYP-2963 packet records, Fejer manifests, Ramanujan manifests, and
residual-section templates.  Namespace: HYP-3002 / T1086 / LTI-152 / LTT-059.
## Lead codex-2026-06-24-S170: dichotomy recursion mode atlas

**Status:** SYNTHESIS / proof-carrier taxonomy complete (HYP-3004/T1088/LTI-154).
User asked to go back through repo work on dichotomies: odd/even,
positive/negative, addition/multiplication, `+1` versus `/2`, `*2` versus `+2`,
and other recursion modes.  This pass integrates older parity ladder,
orbit-functor, signed-LRC, treebolic, Collatz, triune-cycle, Mobius/totient,
and H-loneliness recursion notes with the incoming HYP-3002 curried packet
functional tower, HYP-3003 summand/multiplicand Farey-basis merge,
HYP-2998/HYP-2999/HYP-3000/HYP-3001, and LTI-151 smoothing work.
**Readout:** `04-computation/dichotomy_recursion_mode_atlas_codex_s170.py`
stores output in
`05-knowledge/results/dichotomy_recursion_mode_atlas_codex_s170.out`.  The
scout treats the binary splits as quotient guardrails: each must say which
predicate it preserves, which coordinate it destroys, and what recursion law it
uses.  Tournament Analysis over proof carriers is transitive with path
`additive_pair_sum_face > sign_cut > triune_fraction_recursion > parity_fold >
zeckendorf_path_normal_form > smoothing_switchboard > dyadic_doubling_tree >
farey_product_scheduler > multiplicative_unit_orbit > farey_sum_affine_check >
collatz_affine_halving > plus_two_line_motion`.
**Next:** add six fields to the HYP-2963 packet bank:
`parity_block`, `sign_cut_status`, `additive_pair_sum_owner`,
`multiplicative_unit_orbit`, `recursion_boundary_state`, and `smoothing_route`.
Then test whether every hard non-AP/GW packet is classified by exactly one
primary recursion mode plus named side-channel debts.
Namespace: HYP-3004 / T1088 / LTI-154; companion to HYP-3002 / LTI-152 and
the HYP-3003 / LTI-153 add/multiply fiber merge.
## Lead codex-2026-06-26-S170: LRC14 lacunary exact-gap automaton carrier

**Status:** EVIDENCE / exact scout and proof-interface carrier complete
(HYP-3010/T1094/LTI-158/LTI-159 companion).  User asked to push LRC14 proof thinking through
Fermat-Catalan, the 2-adic Littlewood paper, Ostrowski-Hadamard gaps,
Moser-de Bruijn, fibbinary numbers, and finite automata.
**Readout:** `04-computation/lrc14_lacunary_automaton_carrier_codex_s170.py`
stores output in
`05-knowledge/results/lrc14_lacunary_automaton_carrier_codex_s170.out`.  The
script builds fibbinary and Moser-de Bruijn automata and computes exact LRC14
gaps/safe components for AP, GW, first-13 sequence rows, and first mod-14
residue transversals.  AP and GW are boundary atoms with `M=1/14` and zero
safe mass; first13 fibbinary is strict with `M=3/25`; first13 Moser,
fibbinary residue transversal, and Moser residue transversal are strict with
`M=1/6`.  Tournament Analysis ranks exact Farey scale and exact largest safe
component above dyadic/fibbinary/Moser carry languages, with raw sequence
scalars last.  Rebase signal: HYP-3008/LTI-158 already owns the automatic
gap-language closure lane and HYP-3009/LTI-159 adds the Fermat-Catalan
power-lift extension; HYP-3010 adds exact LRC maximin and safe-component
audits for the same carrier family.
**Next:** add `carry_language` and `automaton_state` to HYP-2963/HYP-3001/HYP-3008/HYP-3009
packet records, then test whether zero-open non-AP/GW residuals emit a
nontrivial dyadic/Ostrowski carry state, have an owner-strip/cross-handoff/
nested-refinement Haar exit, route to K33/THM-572, or collapse to AP/GW.
Namespace: HYP-3010 / T1094 / LTI-158-LTI-159 companion.

## Lead codex-2026-06-24-S169: Farey-Fibonacci additive-basis carrier

**Status:** SYNTHESIS / finite scout and proof-interface carrier complete (HYP-2998/T1083).  User asked to merge previous work on Fermat polygonal numbers, Goldbach/ternary Goldbach, Zeckendorf, the Fibonacci sparse-carry arrangement, and Farey numerator/denominator sum/product/power payloads.  Integrated incoming HYP-2997 as the adjacent cocycle normal-form atlas: this lane supplies the representation-economy label that sequence shadows must carry before HYP-2997-style quotient forgetting is safe.
**Readout:** `04-computation/farey_fibonacci_additive_basis_s169.py` stores output in `05-knowledge/results/farey_fibonacci_additive_basis_s169.out`.  The scout records the Pascal-slope carry family `a_d(n)=sum_k binom(n-(d-1)k,k)`: `d=1` full Pascal/powers of two, `d=2` Fibonacci/Zeckendorf rows `1,1,1+1,1+2,1+3+1,...`, `d=3` Narayana/gap-3 carries, and `d=4..6` higher gap carry gases.  On the golden Farey spine `F_i/F_{i+1}`, `p+q=F_{i+2}`, `p*q=|E(K_{p,q})|`, and power payloads are ordered magnitude-stress channels; on the LRC unit-excess chain `q=14p-1`, `p+q=15p-1` and `p*q=14p^2-p` are kept as separate clocks.  Tournament Analysis over proof carriers is transitive with path `zeckendorf_sparse_normal_form > farey_address_vector > fermat_polygonal_bounded_arity > farey_product_Kpq_area > ternary_goldbach_smoothing > binary_goldbach_pair_graph > farey_power_stress_channel > raw_scalar_rep_count`.
**Next:** add a representation-economy field to sequence-shadow packet classifiers: smoothing, bounded arity, normal form, or Farey address.  A quotient using representation counts is unsafe until it states which economy preserves the target predicate and which coordinate is reconstructed, annihilated, exact, descended, or emitted as residual debt.  Navigation hook: `LTI-155`.
Namespace: HYP-2998 / T1083.
## Lead codex-2026-06-24-S169: Pascal-slope additive-basis Farey packet schema

**Status:** SYNTHESIS / packet-schema lead (HYP-2999/T1084), not a proof.
Companion to HYP-2998's computed Farey-Fibonacci additive-basis carrier.  User
asked to merge prior work on Fermat polygonal numbers, Goldbach/ternary
Goldbach, Zeckendorf, Fibonacci row decompositions, and Farey numerator/
denominator mutations.
**Core readout:** the row pattern `1,1,1+1,1+2,1+3+1,1+4+3,1+5+6+1,...` is
the `d=2` Pascal-slope row `C(n-1-k,k)` with Fibonacci row sum.  Treat this as
a representation-hypergraph carrier before scalarizing: Goldbach is abundant
two-prime fiber, ternary Goldbach is smoothing by one more summand, Fermat
polygonal is bounded arity/residue absorption, Zeckendorf is carry-confluent
unique normal form, and Pascal-slope rows are the finite row-fiber audit.
**Farey merge:** exact `p/q` remains the root.  `p+q` is the additive/summand
lane, `p*q` is the product/coimage lane with `p*q=|E(K_{p,q})|`, and `q^p,p^q`
are ordered magnitude stress tests that may not forget exact scale, endpoint,
or route labels.
**Next:** add `additive_basis_regime`, `representation_entropy`,
`local_residue_rank`, `carry_width`, `pascal_slope_row_id`,
`farey_operator_lane`, and `Kpq_factor_fiber` to the HYP-2963 packet schema;
then test AP/GW, unit-petal/C27, K33, and covering rows under these fields.
Namespace: HYP-2999 / T1084 / LTI-149 / LTT-058.
## Lead codex-2026-06-24-S169: Fibonacci additive-basis Farey bridge

**Status:** SYNTHESIS / proof-interface carrier complete (HYP-3000/T1085).
User asked to merge Fermat polygonal, Goldbach/ternary Goldbach, Zeckendorf,
the Fibonacci row arrangement, and Farey numerator/denominator operations.  The
exact bridge is `F_n=sum_k binom(n-k-1,k)`, the independent-set rank vector of
`P_{n-2}`.  This puts Goldbach, ternary Goldbach, Fermat polygonal, and
Zeckendorf on one proof-economy axis: smoothing entropy, added hypergraph
dimension, bounded arity/residue absorption, and path-normal-form carries.
Incoming HYP-2998 covers the golden Stern-Brocot/Fibonacci spine and HYP-2999
keeps Pascal-slope row fields as packet data; HYP-3000 is the complementary
unit-excess/path-rank classifier.
**Farey pull:** with `M=p/q` and `e=14p-q` retained, `p+q` is
additive/affine-safe, `p*q` is product/incidence, and `q^p,p^q` are magnitude
stress tests.  On `p/(14p-1)`, the table is linear/quadratic as `q=14p-1`,
`p+q=15p-1`, `p*q=14p^2-p`.
**Next:** build the representation-hypergraph TDA proposed in S501, now with
three coordinates: entropy, bounded arity, and carry width; then test whether
HYP-2963 packet debts are smoothing, bounded-invoice, or
Zeckendorf/Ostrowski normal-form objects.  Namespace: HYP-3000 / T1085 /
LTI-150.

## Lead codex-2026-06-24-S168: expanded LRC technique index for tournament/metagraph/series carriers

**Status:** NAVIGATION / compact registry expanded.  User asked to get even more
creative, look back through the many tournament, metagraph, related-series, and
carrier methods, and bring back as many techniques as possible for LRC agents.
`00-navigation/LRC-TECHNIQUE-INDEX.md` now has `156` compact `LTI-*` rows plus
the existing `64` long-form S166 technique-bank entries after merging the
incoming `LTI-109` packet-cocycle atlas and `LTI-110` cocycle-obstruction atlas.
**New promoted block:** `LTI-111..LTI-155` promotes deck derivatives,
Burnside orbit taxes, merged metagraph transport, good-cut/SCC
gas, OCF coimage sectors, noncommutative Redei/Berge recurrences, GLMY
path-homology witnesses, heap/tableau/sorting-network carriers, transfer-matrix
and Hopf convolution carriers, Walsh/Krawtchouk/Paley association shadows,
matroid/gammoid blockers, zeta/Ihara/path torsion, sequence-shadow companions,
Pisano and 2-adic checksums, Hamiltonian-path sheaves, irreducibility no-lift
product rules, relation-lattice/subtorus packets, Faulhaber odd-moment bridges,
A000568 deck lifts, residual metagraph Laplacians, binding-pair switch
tournaments, coimage wall atlases, namespace collision auditing, cocycle normal
forms, residual-section packet grids, summand/multiplicand operation fibers, and
representation-economy carriers.
**Guardrail:** the session explicitly treats the incoming duplicate `HYP-2992`
frontmatter/anchor situation as a namespace-collision signal, not something to
silently overwrite.  The new `LTI-146` and `LTI-TODO-16` ask for an audit and
repair plan so future citations point to the intended quotient/proof object.
**Next:** build the `LTI-TODO-13..29` artifacts, starting with packet-cocycle
formalization, the executable F7 cocycle residual record, a finite cocycle-sheaf
obstruction matrix over HYP-2963 packets, an F0-F7 residual metagraph Laplacian,
and a binding-pair switch tournament for covering rows.  `LTI-148` adds the
residual-section packet grid from HYP-2996 as the executable F7 missing-section
interface, `LTI-149..LTI-153` add Pascal-slope, path-rank, smoothing-switchboard,
curried functional, summand/multiplicand, and dichotomy-mode packet fields, and
`LTI-155` asks classifiers to retain representation
economy before using sequence-shadow counts.
the existing `64` long-form S166 technique-bank entries.  The live tail is
`LTI-147` cocycle normal form, `LTI-148` residual-section packet grid,
`LTI-149` Pascal-slope additive-basis Farey packet schema, `LTI-150`
Fibonacci path-rank bridge, `LTI-151` labelled smoothing switchboard,
`LTI-152` curried packet functional tower, `LTI-153` summand/multiplicand
Farey merge, `LTI-155` representation-economy carrier, and `LTI-156` technique multiverse annex.
**Promoted block:** `LTI-111..LTI-146` promotes deck derivatives, Burnside
orbit taxes, merged metagraph transport, good-cut/SCC gas, OCF coimage sectors,
noncommutative Redei/Berge recurrences, GLMY path-homology witnesses,
heap/tableau/sorting-network carriers, transfer-matrix and Hopf convolution
carriers, Walsh/Krawtchouk/Paley association shadows, matroid/gammoid blockers,
zeta/Ihara/path torsion, sequence-shadow companions, Pisano and 2-adic
checksums, Hamiltonian-path sheaves, irreducibility no-lift product rules,
relation-lattice/subtorus packets, Faulhaber odd-moment bridges, A000568 deck
lifts, residual metagraph Laplacians, binding-pair switch tournaments, coimage
wall atlases, and namespace collision auditing.
**Guardrail:** the session treated duplicate IDs as namespace-collision signal,
not something to silently overwrite.  The rebased repair keeps `HYP-2992` for
the Haar-product tile-discrepancy lane, keeps incoming `HYP-2999/T1084` for the
Pascal-slope packet schema, keeps incoming `HYP-3000/T1085` for the Fibonacci bridge and `HYP-3001/T1086` for the Farey scheduler, keeps incoming S171 `HYP-3003/T1087` for the summand/multiplicand Farey merge, moves cocycle-sheaf exactness to `HYP-3006/T1089`,
and moves the technique multiverse to `HYP-3005/T1090`.
**Next:** build the `LTI-TODO-13..27` artifacts, starting with packet-cocycle
formalization, the executable F7 cocycle residual record, a finite cocycle-sheaf
obstruction matrix over HYP-2963 packets, an F0-F7 residual metagraph Laplacian,
a binding-pair switch tournament for covering rows, familywise HYP-2996
residual-section templates, and additive-basis/Farey-operator packet fields.
Namespace: T1078.
## Lead codex-2026-06-24-S168: residual section and packet-grid verification

**Status:** PROOF-INTERFACE / executable bounded audit complete
(HYP-2996/T1080).  User asked to work on residual section characterization and
packet-grid verification.  The script
`04-computation/lrc14_residual_section_packet_grid_codex_s168.py` imports the
HYP-2963 labelled-packet bank and the HYP-2989 Haar-product grid, then assigns
each packet to a residual section and grid exit.
**Readout:** default run audits `21913` packets.  It finds `7237` hard
`q>=14` packets, `7235` hard non-AP/GW packets, all `7235` with an
owner-strip/cross-handoff/nested-refinement grid exit, `0` zero-open hard
non-AP/GW packets, `0` candidate F7 harmonic sections, and `0` validation
errors.  Section counts: `direct_q_witness_section=14676`,
`ap_gw_boundary_section=2`, `unit_petal_c27_strip_section=4`,
`open_haar_witness_section=6040`, `covering_boundary_moment_section=1188`,
and `k33_state_lift_section=3`.
**Tournament Analysis:** vertices are residual sections / packet-grid exits,
not runners.  The section tournament is transitive with
`score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1}`, `directed_3cycles=0`,
`SCC_sizes=[1,1,1,1,1,1,1]`, and `Hamiltonian_path_count=1`.
**Next:** lift the bounded section map into family templates, then group
Fejer/Ramanujan certificate manifests by residual section rather than scalar
route label.
Namespace: HYP-2996 / T1080.

## Lead codex-2026-06-24-S168: LRC technique multiverse index

**Status:** SHARED INFRASTRUCTURE / broad technique contribution annex created (HYP-3005/T1090).  User asked to get more creative, revisit the many tournament/metagraph/series techniques, and create a big index for other agents to contribute to and pull from.
**Core deliverable:** Added `00-navigation/LRC-TECHNIQUE-MULTIVERSE-INDEX.md`, generated by `04-computation/lrc14_tournament_multiverse_index_codex_s168.py`, with `78` `LTM-*` cards spanning proof sheaves, tournaments/metagraphs, packet families, topology/discrepancy, Fourier certificates, arithmetic series, analytic sieve/smoothing, state-lift geometry, and computation/forum workflow.  The main LRC technique navigation now has `156` compact `LTI-*` rows, `59` tournament-specific `LTT-*` cards, and this 78-card annex.
**Tournament Analysis:** Vertices are technique families, not runners.  Pairwise observable is retained LRC payload (`predicate`, exact scale, topology, endpoint, arithmetic, harmonic, metagraph, formal, residual, transfer); the gauge uses a rotating first-difference priority with a declared tie path.  Fingerprint: `score_hist={2:1,3:1,4:4,5:3}`, `directed_3cycles=26`, `scc_sizes=[9]`, `edge_flips_against_tie_path=21`, `hamiltonian_path_count=2397`.  The all-family SCC says the technique surface is cyclic, not a scalar ladder.
**Next:** promote any theorem-facing `LTM-*` card into `LTI-*`/`LTT-*`; build the finite HYP-3006 `C1` emitted-cocycle matrix on HYP-2963 packet banks; extend packet schemas with exact scale/Haar/Ramanujan/Fejer/endpoint/source-spectrum/state-lift fields; run fiber-mixing tests on `mu`, `mu^2/phi`, divisor, spectrum, `H`, trace, pressure, and edge-density scalars; define F7 as a state lift, harmonic current, or named cohomology class.

## Lead codex-2026-06-24-S166: cocycle obstruction atlas for LRC14

**Status:** PROOF-INTERFACE / cochain-cocycle atlas complete; not a proof
(HYP-2994/T1077).  After HYP-2991 identified the local fixed-margin Haar
cocycle `zeta(T)=T00-T01-T10+T11`, this lane builds the abstract obstruction
atlas above it: exact coboundaries, closed cocycles, torsion/period classes,
Cech gluing classes, gauges with boundary stops, harmonic residuals, and raw
scalar negative controls.
**Readout:** `04-computation/lrc14_cocycle_obstruction_atlas_codex_s166.py`
stores output in `05-knowledge/results/lrc14_cocycle_obstruction_atlas_codex_s166.out`.
The cochain ledger is `C0` packet labels/owner potentials/Fejer centers,
`C1` handoff arrows/endpoint transfers/smoothing gauges/source pullbacks, `C2`
Haar switches/tope curls/color-resonance squares/boundary-moment curls, and
`H2` unpaired mixed modes/no-hidden-kernel survivors/F7/THM-572 residuals.
The `15`-carrier tournament has one nontrivial `3`-carrier SCC tying
certificate handoff, local `zeta`, and exposure/Cech gluing.
**Next:** attach packet-level `zeta` signatures to HYP-2963, tensor Ramanujan
`c_q` labels with Haar tile classes, group Fejer interval certificates by
cocycle signature, and define an executable F7 residual record.
Namespace: HYP-2994 / T1077.
## Lead codex-2026-06-24-S167: LRC14 cocycle carrier atlas

**Status:** PROOF-INTERFACE / cocycle unification atlas and theorem target complete (HYP-2995/T1079).  User asked to see everything in terms of cocycles.  This pass generalizes HYP-2991's local Haar zipper cocycle to the whole HYP-2990 no-free-slider proof stack.  The script `04-computation/lrc14_cocycle_carrier_atlas_codex_s167.py` stores output in `05-knowledge/results/lrc14_cocycle_carrier_atlas_codex_s167.out`.
**Core readout:** A useful LRC cocycle is the signed local obstruction to forgetting a coordinate on a labelled packet fiber.  It closes when local certificates glue and becomes harmless only as a coboundary/potential, dual-annihilated mode, known AP/GW boundary atom, descended family, or named F7/THM-572 residual.  The atlas compares `16` carriers: labelled packet master, Haar `zeta`, endpoint winding, carry/CRT, owner derivative, Ramanujan trace, Fejer moment, boundary moment, product tiling, Farey/K33, C27/unital, root packet, metagraph transfer, sequence shadow, octahedral Hodge current, and OCF coimage activity.
**Tournament Analysis:** Proof-gauge fingerprint: `score_hist={0:1,1:1,2:1,3:1,4:1,6:2,7:1,8:2,10:1,11:1,12:1,14:3}`, `directed_3cycles=4`, `SCC_sizes=[1,1,1,1,1,5,1,1,1,3]`, and `Hamiltonian_path_count=27`; exploration gauge is transitive and has `12` edge flips versus proof gauge.
**Next:** (1) define an `omega_Q` record schema for HYP-2963 packet quotients; (2) compute Haar `zeta` and product-rule classes on actual packet families; (3) turn endpoint-credit potentials into exact interval certificates; (4) prove carry-owner no-leak; (5) define F7 as a non-exact cocycle class with coefficient ring and THM-572 state lift.  Namespace: HYP-2995 / T1079.
## Lead codex-2026-06-24-S167: Cocycle normal form for LRC14

**Status:** PROOF-INTERFACE / cohomological synthesis and tournament carrier atlas complete (HYP-2997/T1082).  User asked to see everything in terms of cocycles.  This pass complements HYP-3006's cocycle-sheaf exactness lane, HYP-2995's packet-cocycle atlas, HYP-2994's obstruction ledger, and HYP-2996's residual-section packet-grid claim and builds directly on HYP-2990's controlled-kernel zipper plus HYP-2991's Haar `zeta` cocycle, then treats endpoint owners, Farey excess, C27 carries, Fourier/PSD duals, exact-period packets, tope cocircuits, tournament `H^1`, chart transitions, curried section derivatives, and F7/THM-572 as channels in one packet complex.
**Core readout:** A proof quotient is safe only when each forgotten cochain is fiber-constant, reconstructed, a coboundary on the quotient fiber, killed by a dual certificate, or routed to a named residual class.  AP/GW become zero-open boundary packets where the current channels close before F7: Farey excess zero, endpoint debt as boundary cocircuit, Haar `zeta` stopped at boundary, AP/GW-C27 branch closed, no Fejer/Toeplitz PSD failure, and no state-lift residual invoked.
**Tournament Analysis:** vertices are cocycle channels / proof obligations, not runners.  Pairwise observable is majority comparison of a nine-coordinate retention vector.  Fingerprint is transitive: `score_hist={0:1,...,12:1}`, `directed_3cycles=0`, `SCC_sizes=[1,...,1]`, `Hamiltonian_path_count=1`, with `labelled_packet_total_cocycle` first and `raw_scalar_shadow` last.
**Next:** build the HYP-2963 packet-level cocycle ledger and tag every low-frontier row by first nonzero class among `F_zeta`, `F_wall`, `F_farey`, `F_carry`, `F_psd`, `F_period`, `F_cocircuit`, `F_tournament`, `F_chart`, `F_lift`, and `F_section`.  Use the rebased POKE coordination checkpoint as the terminal target: F7 is the harmonic/state-lift residual sector, and the HYP-2963 Haar-grid sweep is the finite bank where each first-nonzero class should be tested.  Namespace: HYP-2997 / T1082.

## Lead codex-2026-06-24-S166: Haar zipper cocycle for LRC14

**Status:** PROOF-INTERFACE / local cocycle audit and theorem target complete (HYP-2991/T1075).  User asked to synthesize recent agent work through discrepancy theory, the `2D` Haar product rule, tournament tiling, and more zipper theorems.  This pass rebased over HYP-2990's abstract zipper atlas and makes the local no-free-slider obstruction explicit.  The script `04-computation/lrc14_haar_zipper_cocycle_codex_s166.py` stores output in `05-knowledge/results/lrc14_haar_zipper_cocycle_codex_s166.out`.
**Core readout:** For a `2 x 2` table, the mixed Haar / fixed-margin switch coordinate is `zeta(T)=T00-T01-T10+T11`.  Auditing all nonnegative tables through total `10` gives `1001` tables, `506` margin fibers, `285` nontrivial fibers, `0` duplicate `(margins,zeta)` keys, and fixed-margin zipper-step gcd `4`.  Thus row/column margins are an unsafe quotient unless `zeta` is retained, reconstructed, annihilated, or routed.
**Zipper theorem target:** On every labelled LRC14 packet fiber, each local mixed Haar cocycle is either cancelled by color-compatible discrepancy, stopped at a boundary cocircuit, handed to endpoint/period/certificate/smoothing clocks, descended to a smaller packet family, or converted into a named F7/THM-572 residual.  The depth-4 dyadic product census supplies the tooth taxonomy: orthogonal zero, same-tile boundary atom, owner strip, cross handoff, nested descent, residual cocycle.
**Next:** (1) build the actual HYP-2963 packet grid and compute packet-level `zeta` signatures; (2) count independent color-compatible cocycles and compare with HYP-2595's `k+c_GP`; (3) attach each nonzero owner-strip/cross/nested coefficient to a HYP-2987 handoff arrow; (4) define F7 as unpaired mixed-cocycle state-lift debt.  Namespace: HYP-2991 / T1075.
## Lead codex-2026-06-24-S166: LRC technique index for tournament/metagraph/series methods

**Status:** NAVIGATION / shared contribution index created and merged with the parallel S165 `LTI-*` registry.  User asked to look back through the many tournament, metagraph, sequence, and carrier methods and make a big index other agents can pull from and add to.  Added `00-navigation/LRC-TECHNIQUE-INDEX.md`.
**Contents:** the merged artifact now tracks the compact `LTI-*` registry plus long-form technique entries after adding the incoming HYP-2991 Haar zipper cocycle lane, HYP-2989 product-rule tiling synthesis, HYP-2995 packet-cocycle atlas, HYP-2994 cocycle obstruction atlas, and HYP-2997 normal-form atlas.  Together they cover qdiv and endpoint geometry, Haar/Baire/tope/frontier methods, Fejer/Ramanujan/twist/moment duals, Farey/C27/K33 arithmetic carriers, metagraph transfer and operation-state methods, tournament speedup engines, deck derivatives, source-kernel/zipper/exposure methods, cocycle carrier atlases, OCF/coimage/path-homology residues, sequence shadows, unit-distance analogies, graph-minor guardrails, and formal verification workflows.
**Contents:** the merged artifact has `156` compact `LTI-*` rows plus `64` long-form technique entries after adding the HYP-2991 Haar zipper cocycle lane, HYP-2989 product-rule tiling synthesis, HYP-2995 packet-cocycle atlas, S168 expanded tournament/metagraph/series carriers, HYP-3002 curried packet-functional tower, and HYP-2984 smoothing-switchboard contract.  Together they cover qdiv and endpoint geometry, Haar/Baire/tope/frontier methods, Fejer/Ramanujan/twist/moment duals, Farey/C27/K33 arithmetic carriers, metagraph transfer and operation-state methods, tournament speedup engines, deck derivatives, source-kernel/zipper/exposure methods, cocycle carrier atlases, OCF/coimage/path-homology residues, sequence shadows, unit-distance analogies, graph-minor guardrails, and formal verification workflows.
**Contribution rule:** every new method should state carrier / vertex set, pairwise observable/gauge, preserved LRC predicate, what it destroys if scalarized, anchors, and next contribution.  This explicitly challenges the assumption that LRC Tournament Analysis vertices must be runners.
**Next:** keep this file as the living first-stop index for LRC sessions; update entries when scripts add certificate IDs, packet JSON schemas, endpoint-owner Ramanujan profiles, or F7 residual definitions.
## Lead codex-2026-06-24-S166: LRC14 zipper theorem pattern atlas

**Status:** SYNTHESIS / proof-interface atlas complete (HYP-2993/T1076), concrete LRC14/Haar-Fejer extension of HYP-2990's abstract zipper/no-free-slider atlas.  User asked to pull recent agents together, think discrepancy theory, use the 2D Haar product rule as the same structure as tournament tiling, and search for more zipper theorems.  The script `04-computation/lrc14_zipper_theorem_pattern_atlas_codex_s166.py` stores its output in `05-knowledge/results/lrc14_zipper_theorem_pattern_atlas_codex_s166.out`.
**Core schema:** A zipper theorem has two labelled local certificate sides, a labelled interface, declared stops, and named residuals.  A quotient may forget a coordinate only when the opposite tooth reconstructs it, orthogonality/boundary atoms annihilate it, or the coordinate is emitted as residual data strong enough for the next theorem.
**Computed atlas:** Ten patterns were scored: Haar-Fourier product, Fejer interval packet, tope/cocircuit wall, exposure-poset kernel, Ramanujan exact-period, smoothing/Kaczynski policy, fixed-margin/Johnson sector, apex sheaf gluing, convolution irreducibility lift, and unit-distance cyclotomic norm.  Tournament Analysis is transitive: `score_hist={0:1,...,9:1}`, `directed_3cycles=0`, `SCC_sizes` all singleton, and one Hamiltonian path.  Retention spine: `haar_fourier_product > tope_cocircuit_wall > exposure_poset_kernel > fejer_interval_packet > convolution_irreducibility_lift > ramanujan_exact_period > fixed_margin_johnson > smoothing_kaczynski_policy > apex_sheaf_gluing > unit_distance_cyclotomic_norm`.
**Next:** (1) build a Haar-Fejer compression engine over HYP-2963 packet rows by grouping mixed Haar switch signatures before interval certificate generation; (2) attach primitive-period `c_q` labels to endpoint-wall/Haar rectangle cells; (3) prove the HYP-2988 no-hidden-kernel target after HYP-2992 typed coefficients and HYP-2981 interval certificates are attached; (4) define F7 as a preserved fixed-margin/Johnson harmonic residual; (5) port HYP-2452 convolution no-lift certificates to LRC blocker ledgers.  Namespace: HYP-2991 / T1075.
## Lead codex-2026-06-24-S166: LRC tournament technique index

**Status:** SHARED INFRASTRUCTURE / tournament-specific companion index created (`00-navigation/LRC-TOURNAMENT-TECHNIQUE-INDEX.md`, T1078) and expanded by T1081 to `57` reusable cards, downstream of the broader `00-navigation/LRC-TECHNIQUE-INDEX.md`.  User asked to look back through the many tournament, metagraph, series, and related techniques developed in the repo and bring them back to bear on LRC14.
**Core deliverable:** Each card records the move, LRC use, preserved predicate, forgotten coordinate/guardrail, next pull, and pointers.  Covered families include controlled-kernel zipper theorems, proof-carrier Tournament Analysis, multi-scale winding spectra, observer-source/A000568 fibers, fixed-path tilings, metagraph/even-graph projections, deletion-contraction depth, flip-energy ledgers, Burnside orbit-cost perturbations, score-class H-spread, path-homology Betti carriers, H-gap state lifts, metric-comparator filters, nonabelian character-ratio carriers, reversal/converse-Z2 ambiguity, round-tournament realizability, Paley/Frobenius arithmetic carriers, good-cut/SCC coordinates, H/OCF/W-polynomial channels, automata, Landau prefix majorization, packet-sign K4s, sector DPs, taut bridges, topes/cocircuits, Haar-Baire/Haar-product carriers, Fejer/Toeplitz certificates, moment barriers, Ramanujan exact-period projectors, divisor/totient guardrails, analytic sieve/Kaczynski smoothing, Farey/Kpq/C27/K33 state lifts, exposure posets, parity carriers, transport matrices, spectral shadows, unit-distance analogies, octahedral currents, OCF coimages, labelled packet classifiers, Faulhaber channel gates, irreducible-polynomial lift analogies, p-adic packet trees, Robbins no-bridge assembly, and standard fingerprint payloads.
**Proof-order use:** Treat the index as a menu for packet-preserving reductions, not as a bag of metaphors.  Every pull should answer: which coordinate may this quotient forget, and why?  The default controlled-kernel rule is fiber-constant, reconstructible, annihilated by a dual certificate, or routed to a named residual sector.
**Immediate next steps:** (1) add Haar tile class, Ramanujan projector, and spectrum binding scale to the HYP-2963 labelled packet classifier; (2) add a Robbins no-bridge checklist to Fejer manifests; (3) compute multi-scale spectra for named LRC14 rows and the weakest Fejer-margin rows; (4) build a decorated source-cone canonicalizer; (5) test the HYP-2992 vanishing lemma on packet fibers.

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
**Next:** (1) build the actual packet grid for named rows AP, GW, `12->36`, `10->20`, `13->26`, `P10+GW`, `12->168`, and lcm-tail rows; (2) compute typed Haar coefficients by packet family; (3) compare same-tile boundary atoms with the HYP-2986 boundary cocircuit owner-pair sums; (4) test whether HYP-2981 Fejer atoms cluster in owner-strip or nested-refinement classes; (5) formalize the dictionary from fixed-path tournament staircase tiles to Haar rectangles.  Namespace: HYP-2992 / T1072.
## Lead codex-2026-06-24-S165: abstract zipper theorem atlas across past topics

**Status:** EXPLORATORY PROOF-TECHNOLOGY / HYP-2990/T1074 created after rebasing over incoming HYP-2988 exposure, HYP-2989 Haar-square, and HYP-2992 Haar-tile lanes.  User asked for more zipper theorems and free creative exploration of past topics.  Added `04-computation/abstract_zipper_theorem_atlas_codex_s165.py` and stored `05-knowledge/results/abstract_zipper_theorem_atlas_codex_s165.out`.
**Status:** PROOF-INTERFACE / exact finite product-table scout complete (HYP-2992/T1072), companion to HYP-2989/T1073's minimal Haar-square switch, HYP-2991/T1075's Haar zipper cocycle, and HYP-2988's exposure-poset audit.  User prompt asked to synthesize recent agents through discrepancy theory and the two-dimensional Haar product rule, specifically the way it creates the same structure as the tournament tiling model.  The script `04-computation/lrc14_haar_product_tile_discrepancy_codex_s165.py` enumerates dyadic Haar rectangle products and stores the readout in `05-knowledge/results/lrc14_haar_product_tile_discrepancy_codex_s165.out`.
**Core readout:** Through depth `3`, `225` rectangles give `50625` ordered products: `43736` orthogonal zeros, `225` same-tile boundary atoms, `1020+1020` one-coordinate owner strips, `2312` cross handoffs, and `2312` nested refinements.  Every nonzero non-atom class is sign-balanced.  This suggests a discrepancy proof by typed coefficient detection rather than scalar positivity.
**Transfer to LRC14:** Use a rectangular packet grid whose axes are endpoint/tope walls and proof clocks (exact period, Fejer scale, Ramanujan packet, K33/state-lift route).  Then prove a vanishing lemma: if a primitive zero-open packet has no owner-strip, cross-handoff, or nested-refinement Haar coefficient, it collapses to AP/GW boundary skeleton or creates the missing THM-572/F7 state-lift atom.
**Next:** (1) build the actual packet grid for named rows AP, GW, `12->36`, `10->20`, `13->26`, `P10+GW`, `12->168`, and lcm-tail rows; (2) compute typed Haar coefficients by packet family; (3) compare same-tile boundary atoms with the HYP-2986 boundary cocircuit owner-pair sums; (4) test whether HYP-2981 Fejer atoms cluster in owner-strip or nested-refinement classes; (5) formalize the dictionary from fixed-path tournament staircase tiles to Haar rectangles.  Namespace: HYP-2992 / T1072.
## Lead codex-2026-06-24-S165: abstract zipper theorem atlas across past topics

**Status:** EXPLORATORY PROOF-TECHNOLOGY / HYP-2990/T1074 created after rebasing over HYP-2988 exposure, HYP-2989 Haar-square, HYP-2991 Haar zipper cocycle, and HYP-2992 Haar-rectangle discrepancy lanes.  User asked for more zipper theorems and free creative exploration of past topics.  Added `04-computation/abstract_zipper_theorem_atlas_codex_s165.py` and stored `05-knowledge/results/abstract_zipper_theorem_atlas_codex_s165.out`.
**Core rule:** A zipper theorem is a controlled-kernel theorem.  A quotient may forget a coordinate only if the target predicate is fiber-constant, the coordinate is reconstructible, a dual certificate annihilates it, or it is routed to a named residual sector.
**Carriers compared:** LRC14 certificate handoff, kernel/tope/smoothing, octahedral Hodge current, C27/unital pair completion, Farey/K33 incidence, boundary-moment multi-chart, shell-1/root packet, unit-distance endpoint ear, OCF activity/coimage, good-cut/SCC support, and raw scalar shadow.
**Tournament fingerprint:** `score_hist={0:1,1:1,2:1,4:2,5:1,6:1,7:2,9:1,10:1}`, `directed_3cycles=4`, `SCC_sizes=[1,1,1,6,1,1]`, `Hamiltonian_path_count=15`.  The six-carrier middle SCC says several past topics form typed handoff teeth rather than a scalar hierarchy.
**Next:** (1) define HYP-2987's `F7` as a named harmonic/state-lift residual sector; (2) turn the no-free-slider rule into an LRC14 quotient checklist; (3) test family compression on K33/petal Fejer packets; (4) revisit octahedral divergence/curl and good-cut/SCC support as state-lift coordinates.

## Lead codex-2026-06-24-S164: admissible smoothing dispatcher for LRC14

**Status:** PROOF-INTERFACE / routing theorem target complete (HYP-2985/T1069), complementary to the incoming HYP-2984/T1068 kernel-homotopy stub.  This pass turns the recent analytic-number-theory prompts into a typed LRC14 smoothing dispatcher rather than a scalar estimate.  The script `04-computation/lrc14_smoothing_dispatcher_codex_20260624.py` classifies which policy is allowed to handle each live packet family and stores the readout in `05-knowledge/results/lrc14_smoothing_dispatcher_codex_20260624.out`.
**Status:** PROOF-INTERFACE / routing theorem target complete (HYP-2985/T1069), complementary to the HYP-2984/T1068 kernel-homotopy and switchboard lane.  This pass turns the recent analytic-number-theory prompts into a typed LRC14 smoothing dispatcher rather than a scalar estimate.  The script `04-computation/lrc14_smoothing_dispatcher_codex_20260624.py` classifies which policy is allowed to handle each live packet family and stores the readout in `05-knowledge/results/lrc14_smoothing_dispatcher_codex_20260624.out`.
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
**Rebase signal:** HYP-2984/T1068 is now the sibling kernel-homotopy and boundary-defect lane, and HYP-2985/T1069 is the admissible-smoothing dispatcher lane.  Treat both as admissible-arrow proofs inside this HYP-2986 handoff atlas: smoothing or Fourier kernel deformations must preserve packet certificates or emit named boundary defects, and smoothing policies must keep endpoint/exact-period/kernel/far-approach clocks separate.
**Next:** (1) make O3 concrete by compressing selected Fejer certificates into packet-family templates; (2) make O5 concrete by constructing the THM-572 state lift from zero-open non-AP/GW residuals; (3) formalize O6 as a Johnson-harmonic residual predicate with count-sector separation; (4) write the O1 source-kernel exclusion as a typed reduction theorem over qdiv/Farey/Haar/endpoint labels.  Namespace: HYP-2986 / T1070.
## Lead codex-2026-06-24-S164: smoothing switchboard for labelled LRC14 packet routes

**Status:** PROOF-INTERFACE / finite route matrix complete (HYP-2984/T1068).  This pass combines HYP-2981 Fejer interval certificates, HYP-2982 analytic `Phi`/`G` weights, HYP-2983 Kaczynski/exponential-sum synthesis, and HYP-2979 Ramanujan exact-period projectors into one dispatcher.  The switchboard audits `16` named/probe packets and classifies them as AP/GW boundary atoms (`2`), interval-Fejer certificates (`5`), Ramanujan pre-split handoffs (`6`), covering/lift boundary-moment rows (`2`), and K33 state-lift debt (`1`).  Main proof slogan: packet route first, smoothing kernel second.
**Guardrail:** A kernel is admissible only with its retained side channel: endpoint zero-credit/Kaczynski class for AP/GW, packet key/center/degree/interval bound for Fejer, first strict q/primitive packet/endpoint labels for Ramanujan, lift chart/exact safe mass/late-q label for covering, and K33 owner/state-lift debt for K33.  Raw scalar smoothing choices and unlabelled large-sieve bounds are demoted.
**Next:** Scale the switchboard from this named audit bank to the full HYP-2963 packet bank; attach HYP-2981 Fejer records where available, HYP-2979 first strict q data everywhere feasible, and count any remaining unclassified route families.  Namespace: HYP-2984 / T1068.

## Lead codex-2026-06-24-S164: tope-wall cocircuit and Farey route-scheduler follow-ons

**Status:** PROOF-INTERFACE / concurrent namespace repair complete.  HYP-2986/T1070 is the oriented-tope wall/cocircuit pass: cut the LRC14 time circle by all danger endpoints, treat open cells as topes, keep AP/GW as boundary-cocircuit equality atoms, and treat a no-tope/no-cocircuit packet as the real strict-counterexample falsifier.  HYP-3001/T1086 is the Farey mutation certificate scheduler: preserve exact `M=p/q` and excess `e=14p-q`, then use product collapse `p` only after the unit-excess gate as a route selector, not as an order-safe scalar.  HYP-2992/T1072 remains the Haar-product tile namespace, HYP-2993/T1076 remains the zipper-pattern atlas, HYP-2990/T1074 remains the abstract zipper atlas, HYP-2991/T1075 remains the Haar zipper cocycle lane, HYP-2994/T1077 is the cocycle obstruction atlas, HYP-2995/T1079 is the cocycle carrier atlas, HYP-2996/T1080/LTI-148 is the residual-section packet-grid lane, T1081 is the tournament-index expansion, HYP-2997/T1082/LTI-147 is the cocycle normal-form atlas, HYP-2998/T1083/LTI-155 is the Farey-Fibonacci representation-economy carrier, HYP-2999/T1084/LTI-149 is the Pascal-slope packet schema, HYP-3000/T1085/LTI-150 is the Fibonacci path-rank bridge, HYP-3003/T1087/LTI-153 is the summand/multiplicand Farey basis merge, and HYP-3004/T1088/LTI-154 is the dichotomy recursion mode atlas.
**Next:** For HYP-2986, prove every primitive no-tope/no-cocircuit wall packet violates endpoint-owner parity, slides to an open tope, routes to K33/H=7 state lift, or defines F7.  For HYP-3001, turn the `exact M/e -> product-collapse p -> packet family` dispatch into a theorem-facing front end for Fejer/Ramanujan/Kaczynski/state-lift certificates.

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
**Source:** S20 trawl + drill workflows (HYP-2760..2763). **Status:** ACTIVE Ã¢â‚¬â€ LAYER 3 (consec Schur-maximizes measS7) is the last wall; LAYERS 1-2 PROVED+Lean.
- **Conductance invariant (HYP-2760):** measS7 = EXACT function of c_r = sum_{e=r mod7} 1/|e| (0 collisions); Foster sum_r g(r)=112; bottleneck/minimax on C_7. NEXT: derive the exact F(c) and prove consec maximizes the bottleneck via Chebyshev equalization.
- **Windows extremality (HYP-2761):** consec UNIQUELY maximizes WIN=harmonic sum of binding speeds (0 ties, sharper than the wall). OBSTRUCTION: WIN/DISCONNECTED split does NOT compose (they anti-correlate). NEXT: harmonic-sum exchange argument.
- **AP-unification (HYP-2762):** BOTH extremalities = 'the additive AP Schur-maximizes' (consec / cyclic Interval); Paley is small-p only (crossover p=19, HYP-479/THM-135). Wiener-Khinchin: argmin sum|lambda|^4 = QR difference set. NEXT: the elementary-symmetric c_k crossover mechanism; does it transfer to LRC?
- **ROUTE 5 RESOLVED (mac-mini-S20, the c_k crossover + transfer test):** The unification is NOT one c_k theorem Ã¢â‚¬â€ the two extremalities are driven by OPPOSITE level + OPPOSITE-sign fugacity. (i) The spectral "H = c_0 + sum c_k e_k(|lambda|^2)" of THM-134 is UNDERDETERMINED by circulant data (only 2 classes at p=7, 4 at p=11 vs m unknowns) and NON-canonical; THM-505 already shows H is not a low-degree spectral polynomial. The GENUINE canonical decomposition on both sides is a SIGNED SUBSET-COMPLEX LEVEL-SUM. (ii) TOURNAMENT: H = sum_j (+2)^j alpha_j (odd-cycle independence poly, fugacity +2). Paley wins level 1 (max cycle counts c_ell, every ell), Interval/AP wins levels j>=2 (disjoint packings). Crossover fugacity x*(p): 3.00 (p=7) -> 2.61 (p=11) -> below 2 by p=19 (= THM-135). The AP wins because +2 AMPLIFIES the high-order packing levels where it dominates. (iii) LRC: measS7 = sum_k (-1)^k MISS_k (sector-miss inclusion-exclusion, fugacity -1; MISS_k = sum_{|S|=k} P(all of S missed)). consec UNIQUELY MINIMIZES MISS_1 (best single-sector coverage, rank 1/319) and wins via the LOW-order level that the alternating -1 fugacity rewards. Sweep G_x = sum_k x^k MISS_k over all 319 shapes: at x=-1,-1/2 consec is unique argmax (rank 1); at x=+1/2,+1,+2 consec drops to rank 319/319, 316/319, 220/319 (worst/near-worst). So consec-optimality is SPECIFIC to x=-1. PRECISE OBSTRUCTION (refines HYP-2758 no-Parseval-simplex): the tournament AP wins HIGH-order packing at fugacity +2; the LRC AP wins LOW-order coverage at fugacity -1 Ã¢â‚¬â€ same "additive AP extremal" surface, reversed level and reversed fugacity sign. Scripts: 04-computation/lrc14_route5_*.py; results: 05-knowledge/results/lrc14_route5_*.out.
- **KEY LITERATURE (HYP-2763):** **Rosenfeld 2025 (arXiv:2509.14111, arXiv:2511.22427) PROVES LRC n=8,9,10** via finite-checking + the Tao 2018 -> Malikiosis-Santos-Schymura 2025 bound n^{2n} (arXiv:2411.06903). STRATEGIC Q: is the finite check feasible at n=14 (n^{2n}=14^28 brute is not, but our sector-route + stratum-localization + sharp L7 atlas p<=14 reduces it massively)? Huffer-Shepp 1987 reflection lemma (port to per-cell W_a). Cusick view-obstruction (AP {1..n} proven tight). Chung-Graham 3-distance (AP gap-rigidity => LAYER 3 as a rigidity statement).

---

## Lead kind-pasteur-2026-06-17-S2: close LRC(14) covering case via a BOUNDED-SPEED reduction (THM-525, OPEN-Q-108)

**Status:** Ã°Å¸â€Â´ HIGH Ã¢â‚¬â€ the single most valuable next step toward LRC(14). THM-525 reduces the covering case to OPEN-Q-108 (Q=uniform meas(G_C)Ã¢â€°Â¥c) and has ~105k covering 13-sets with ZERO counterexamples, but the scan is EVIDENCE, not a certificate. **The move that would CLOSE it:** a bounded-speed reduction (Tao "Some remarks on the lonely runner" Thm 1.3; MalikiosisÃ¢â‚¬â€œSgallÃ¢â‚¬â€œSomer / Dubickas finite-checking) proving LRC(14) on covering sets need only be verified for v_max Ã¢â€°Â¤ an explicit VÃ¢â€šâ‚¬. Then the exact scan becomes a finite proof up to VÃ¢â€šâ‚¬ and only the VÃ¢â€šâ‚¬-bound remains. SECOND target: the named sub-gap **G2** (transversality) Ã¢â‚¬â€ prove G_C (a union of arcs whose endpoints are binding-pair crossings k/(v_aÃ‚Â±v_b), all v Ã¢â€°Â¤ bound) cannot concentrate in O(1/w)-neighborhoods of {a/w}, so w's danger comb (measure 1/7) cannot CONTAIN G_C. Either one moves OPEN-Q-108 from "localized + strongly evidenced" to "closed."
**Next (ranked):** (1) find/derive the explicit VÃ¢â€šâ‚¬ for N=14 covering sets from Tao Thm 1.3 / MSS; (2) run the exact covering-set scan to v_maxÃ¢â€°Â¤VÃ¢â€šâ‚¬ (the THM-525 tooling already does covering enumeration); (3) attack G2 directly via the binding-crossing-endpoint structure of G_C; (4) test whether the circle/local-tournament structure (HYP-2576) pins G_C away from the parked comb.
**Sources:** THM-525, OPEN-Q-108 (GAP A + new GAP G2), HYP-2573..2576; `04-computation/lrc14_{easy_dominates_hard,parked_centering,verify_reduction_*}_kps-S2-wf.py`; reflection `the-perfect-middle-is-the-wrong-fixpoint-and-the-lrc-tournament-is-local-not-lonely-kps-S2.md`.

## Lead mac-mini-2026-06-18-S7: ANGLE C Ã¢â‚¬â€ three-gap / dilation-orbit extremality of meas(S7) (THM-532 closure)

**Status:** LITERATURE-SCOUTED + COMPUTE. The THM-532 crux ("consecutive AP maximizes meas(S7) for dangerous rows k=8..11") is now grounded against the three-gap / dilation-orbit literature, and the deployment is HONESTLY MAPPED. KEY EXTERNAL CONFIRMATIONS: (a) The LRC state of the art (PhamÃ¢â‚¬â€œSawhneyÃ¢â‚¬â€œTidorÃ¢â‚¬â€œTrakulthongchai? "Eleven, twelve, and thirteen lonely runners", arXiv 2604.23906) proves LRC(k) for kÃ¢â€°Â¤12 and **identifies (1,2,Ã¢â‚¬Â¦,k) Ã¢â‚¬â€ the consecutive AP Ã¢â‚¬â€ as the unique bottleneck / hardest tuple** that forces the polynomial-method step; k=13 (=LRC(14)) is genuinely the first open case. This is INDEPENDENT literature confirmation that "AP is the extremal configuration" for the dangerous rows. (b) Bedert, "Riesz products and the LRC" (arXiv 2511.16636): `ML(1,2,Ã¢â‚¬Â¦,n)=1/(n+1)` exactly Ã¢â‚¬â€ consecutive integers are the *unique* extremal case Ã¢â‚¬â€ and the proof dichotomy is **low additive dimension (AP-like) vs. high additive dimension (dissociated Ã¢â€ â€™ Riesz product)**, the SAME split as THM-532's relation-height (low W = AP-rich vs. high height = dissociated). (c) Tao "Some remarksÃ¢â‚¬Â¦": union bound is tight to `(1+O(1/n))` exactly for AP-like sets. (d) The orbit `{0,x,Ã¢â‚¬Â¦,(k-1)x}` of the AP is EXACTLY a Steinhaus three-gap configuration (Ã¢â€°Â¤3 distinct gaps Ã¢â‚¬â€ verified k=8..13); a general E gives a union of d dilated APs with up to **3d gaps** (Liang 1979, "A short proof of the 3d distance theorem", Discrete Math 28; generalization arXiv 1910.00865) Ã¢â‚¬â€ so non-AP orbits are LESS uniform, the mechanism behind AP-extremality. COMPUTE (this session): AP is the UNIQUE maximizer of BOTH meas(S7) AND the difference-energy `D(E)=ÃŽÂ£_d mult(d)Ã‚Â²` for k=8..11 (exact small-box census, 0 beaters); but Pearson corr(meas(S7),D)Ã¢â€°Ë†0.53Ã¢â‚¬â€œ0.64 only, so **meas(S7) is NOT a monotone function of any single scalar** Ã¢â‚¬â€ AP-extremality is a JOINT extremum, not reducible to a one-line three-gap inequality. HONEST VERDICT: the literature CONFIRMS and CONTEXTUALIZES AP-extremality (strong external corroboration) but does NOT supply a ready-made theorem that PROVES `meas(S7(E))Ã¢â€°Â¤meas(S7(AP_k))` for k=8..11; the three-gap law bounds the *number* of gaps, not the *cover measure*. No clean three-gap closed form for `meas(S7(AP_k))` either (denominators 1470,5880,17640,5880,194040,21560 = 7Ã‚Â·lcm-type, not a tidy product).
**Next (ranked):** (1) Try the difference-energy as a TWO-piece certificate: prove (i) meas(S7) Ã¢â€°Â¤ f(D(E)) for an explicit increasing envelope f, then (ii) the finite residual {E : D(E) Ã¢â€°Â¥ DÃ¢â€šâ‚¬} is AP-near and small Ã¢â‚¬â€ this mirrors THM-532's W-split but with D (a genuine additive-energy, literature-standard) instead of the ad-hoc W. (2) Pursue the Liang 3d-distance bound as a hard cap: if orbit(E) has Ã¢â€°Â¥ 4 distinct gaps on a positive-measure set, meas(S7) loses a quantifiable amount vs the 3-gap AP Ã¢â‚¬â€ make this quantitative (the missing inequality). (3) Read the FULL 2604.23906 / 2511.16636 PDFs (WebFetch failed on binary; use the HTML or ar5iv) for the exact AP-bottleneck lemma Ã¢â‚¬â€ it may directly bound the AP cover and transfer. (4) Check whether Bedert's Riesz-product certificate, applied to the dangerous *clusters* (not all of V), gives meas(S7)Ã¢â€°Â¤cap directly for high-additive-dimension E, leaving ONLY the AP-near finite residual.
**Sources:** THM-532; `04-computation/lrc14_angleC_threegap_{extremality,mechanism}_macmini_0618s7.py`; `05-knowledge/results/lrc14_angleC_threegap_*_macmini_0618s7.out`; arXiv 2604.23906, 2511.16636, 1910.00865; Liang Discrete Math 28 (1979) 325; Tao "Some remarks on the lonely runner conjecture" (2017); BalogÃ¢â‚¬â€œGranvilleÃ¢â‚¬â€œSolymosi arXiv 1410.8404.

## Lead codex-2026-06-18: HYP-2608 empty-window moment lower certificate for LRC(14)

**Status:** NEW complement-side certificate.  Instead of upper-bounding the bad seven-sector cover `S7(E)`, this route lower-bounds the good anchored empty-window union `EWLB_A(E)=meas(union_j W_j(E))`, where `W_j` is the event that the open arc `(j/14,j/14+1/7)` is empty.  Let `R(x)=# {j:x in W_j(E)}` and `T_s=E[C(R,s)]`.  Polynomial minorants `h_R(r)<=1[r>0]` give proved per-set lower bounds `EWLB_A(E)>=Phi_z(E)=sum z_s T_s(E)`.  Degrees `4,3,3,2,1` clear the AP rows `k=8..12`; exact bounded primitive banks (`5797` rows total) have `0` threshold failures, with AP minimizers for dangerous rows `k=8..11`.  This matches the user's region-first prompt: the natural tournament vertices are fixed loop regions/proof obligations, not runners.  AP region-load tournaments are transitive with score histogram `0..6`.
**Next (ranked):** (1) Prove the scalar rearrangement `Phi_k(E)>=Phi_k(AP_k)` for `k=8..11`, likely by a three-gap/Sturmian or relation-height split for the empty-window count distribution. (2) For `k=12`, prove a coarse first-moment bound; AP extremality is unnecessary because the degree-1 certificate has large slack. (3) Compare HYP-2608's minorant polynomials to THM-534's majorants and look for a shared Krawtchouk/Bonferroni algebra that proves both scalar extremalities. (4) Build a fast integer common-refinement engine for `T_s(E)` and widen the exact banks to match HYP-2604's AP-frontier boxes.
**Sources:** `04-computation/lrc14_empty_window_moment_lower_codex_0618.py`, `05-knowledge/results/lrc14_empty_window_moment_lower_codex_0618.out`, `05-knowledge/hypotheses/HYP-2608-lrc14-empty-window-moment-lower-certificate.md`; HYP-2603, HYP-2604, THM-531, THM-532, THM-534.

## Lead codex-2026-06-17-S4: unit-distance n=22 Mathieu residue (HYP-2572, T840)

**Status:** NEW exact finite scout.  The unresolved `60 <= u(22) <= 61` frontier is a one-point extension problem: any 61-edge graph has a degree-4 or degree-5 deletion ear over a `57/56`-edge 21-core.  The Mathieu chain suggests retaining the `M22` point-stabilizer side-channel: fixing a point exposes `M21 = L_3(4) = PSL(3,4)` on `PG(2,4)`, whose 21 lines are exactly the fixed-point residual hexads of `S(3,6,22)`.  The new script verifies the plane and classifies ear-neighbor sets: degree-5 has `line_5` (`21`), `near_line_4_plus_1` (`1680`), `arc_5_no_three` (`1008`), and `three_collinear_5` (`17640`); degree-4 has `line_4` (`105`), `arc_4_no_three` (`2520`), and `three_collinear_4` (`3360`).  Unit-circle caps give at most `4`/`3` internal unit chords for degree `5`/`4` ears.  Tournament Analysis over proof carriers is transitive with leader `line_5_hexad_ear`.
**Next (ranked):** (1) Attach `PG(2,4)` ear types to stored `n=21` core and `n=22` extension candidates; record which extensions are line, punctured-line, near-line, or scattered. (2) Prove/certify coherent-ear exclusions using the circle-cap counts plus the existing Moser cap-endpoint ledger. (3) Build a scattered-ear obstruction library indexed by secant profile and compare it against the totally-unfaithful killers from the graph-only 62-edge coimage. (4) Test whether the `M21` residue aligns with `P_2^- / P_2^+` spine-ladder rows and the `57=20+37` centered-hex carrier split.
**Sources:** `04-computation/unit_distance_n22_mathieu_residue_codex.py`, `05-knowledge/results/unit_distance_n22_mathieu_residue_codex.out`, `05-knowledge/hypotheses/HYP-2572-unit-distance-n22-mathieu-residue.md`, `07-reflections/unit-distance-n22-mathieu-residue-codex.md`; HYP-2176, HYP-2188, HYP-2203, HYP-2467.

## Lead mac-mini-2026-06-17-S2: LRC regions/sections reframe Ã¢â€ â€™ the binding-pair switch reduction (THM-524, HYP-2571)

**Status:** PROVED (sawtooth lower-envelope lemma): M(S)=max(g(Ã‚Â½), max over pairs (a,b),k of g(k/(v_aÃ‚Â±v_b))) Ã¢â‚¬â€ the LRC gap is a pairwise SWITCH; LRC(N) Ã¢Å¸Âº some pair-crossing gives gapÃ¢â€°Â¥1/N with others clear (~78 switches for 13 runners, polynomial). Regions/SDR = the on-grid q=14 witness (blind off-grid). Complement=reversal is the one exact tournament bridge; overtaking tournament is trivially transitive (RÃƒÂ©dei link dead). Covering hard core: M=7m/(84m+5), min 7/89.
**Next (ranked):** (1) the THM-524 NEXT step Ã¢â‚¬â€ prove "grid-sharpness (gridM=M) Ã¢Å¸Âº a binding complement pair (v,NÃ¢Ë†â€™v)" and that it fails EXACTLY for covering sets, merging the section lens with THM-523 (off-grid criterion). (2) Merge with codex's Hall/wall-switch program (T838) and dihedral mouth-exchange (T837/HYP-2569) Ã¢â‚¬â€ three takes on the same regions prompt; the binding-pair (mine), the Hall packets (codex), the dihedral mouths (codex) should unify on the covering hard core. (3) Bound inf M over covering sets Ã¢â€°Â¥ 7/89 (the gap-side compactness frontier = THM-522's measure-side). (4) Is there a NON-trivial tournament built from the SWITCHES (crossing parities) rather than the snapshot order? Ã¢â‚¬â€ the only place a RÃƒÂ©dei/OCF link could survive.
**Sources:** 04-computation/lrc14_{regions_sections_grounding,binding_pair_reduction,gap_M_exact}_mac-mini-2026-06-17-S2.py + workflow angle scripts (sdr_hall_symmetry, angleD_overtaking_tournament, angle_c_exotic_switches, angle_e_cross_modulus) + .out; THM-524, HYP-2571; reflection regions-not-runners-the-lonely-gap-is-a-pairwise-switch-mac-mini-S2.md; THM-523/522, codex T837/T838/HYP-2568/2569/2570.

## Lead mac-mini-2026-06-16-S3: LRC(14) GAP-SIDE reduction to covering sets (THM-523, HYP-2566/2567)

**Status:** PROVED the q-witness reduction: LRC(14) Ã¢Å¸Âº M(S)=max_Ãâ€ž min_v||vÃâ€ž|| Ã¢â€°Â¥1/14 for all primitive COVERING SETS (a multiple of every qÃ¢Ë†Ë†{2,Ã¢â‚¬Â¦,14}); everything else is MÃ¢â€°Â¥1/14 via Ãâ€ž=1/q. No counterexample over extensive search; residual covering-set min M = 7/89Ã¢â€°Ë†0.0787 (10% margin). The disproof search produced the proof reduction.
**Next (ranked, to finish LRC(14)):** (1) **HYP-2566 Ã¢â‚¬â€ bound inf M over covering sets > 1/14.** Covering sets are bounded/structured (need multiples of 8,9,11,13,14 Ã¢â€ â€™ spread out); combine scale-invariance + kind-pasteur's bounded-lcm compactness (THM-522, measure-side) on the GAP side. Is 7/89 the true inf? Search covering sets with larger speeds + 2 multiples of 14. (2) Prove the residue-refinement ÃŽÂ´-existence lemma rigorously (the Ãâ€ž=a/14+ÃŽÂ´ construction for C'(14) with multiple large multiples of 14). (3) The finite-witness-cover endgame: a fixed finite family of Ãâ€ž's (1/q's + a/14+ÃŽÂ´'s) certifying MÃ¢â€°Â¥1/14 for ALL covering sets except a finite checkable residual (the destroy-lonely-points agent: 12 witnesses cover all but 12/2000). (4) Connect M-side covering reduction to L-side inf=1/1260 (HYP-2561) Ã¢â‚¬â€ both fence the same hard core.
**Sources:** 04-computation/lrc14_{gap_M_exact,constructive_dichotomy,residue_criterion,residual_hardcore_minM,strong_necessary_condition}_mac-mini-2026-06-16-S3.py + ~30 workflow angle scripts (lrc14_{disprove_*,covering_obstruction_proof,M_quantization_floor,destroy_lonely_points,family_reduction,Ã¢â‚¬Â¦}) + .out; THM-523, HYP-2564/2565/2566/2567; reflection the-disproof-search-builds-the-proof-lrc14-covering-reduction-mac-mini-S3.md; THM-360/398/501/522.

## Lead mac-mini-2026-06-16-S2: LRC-14 inf L>0 Ã¢â‚¬â€ corroboration + the four closed doors (HYP-2562/2563, T834, complements kind-pasteur THM-522)

**Status:** kind-pasteur-S7 owns the headline (THM-522 scale-invariance+quantization; HYP-2561 inf=1/1260; MISTAKE-075). I independently CORROBORATED inf=1/1260 (single perturbations wÃ¢â€°Â¤600, doubles wÃ¢â€°Â¤140: nothing in (0,1/1260)) and CLOSED four positivity routes (LLL/Shearer, SelbergÃ¢â‚¬â€œBeurling, Abel/CesÃƒÂ ro, OCF-bridge), all blocked by ONE obstruction: L is an archimedean SIGNED singular integral with no termwise floor + positively-correlated danger events.
**Next (ranked, toward inf L>0):** (1) **The crux per THM-522/HYP-2561:** classify the (conjecturally finite, bounded-lcm) TIGHT LOCUS = LRC(14) extremal configs (GoddynÃ¢â‚¬â€œWong); then quantization+compactness give inf LÃ¢â€°Â¥1/1260. Extend my `lrc14_perturbation_inf_search` to 3-element perturbations + non-AP bases to map the full tight locus. (2) Prove the UNIFORM (core-independent) large-speed decoupling constant C (ErdÃ…â€˜sÃ¢â‚¬â€œTurÃƒÂ¡n discrepancy of stranger 14m against the fixed bounded core's lonely set; constant via #core-runners + min core-gap) Ã¢â‚¬â€ the one missing ingredient for the mÃ¢â€ â€™Ã¢Ë†Å¾ tail (my workflow Angle 8). (3) The Bedert level-bound |E_kÃ¢Ë†Â©P|Ã¢â€°Â¤(C log|P|)^k Ã¢â€ â€™ Abel-controlled alternating ÃŽâ€º_k (OPEN-Q-097, the complementary analytic route). DO NOT re-attempt the 4 closed doors.
**Sources:** 04-computation/{lrc14_exact_rational_measure,lrc14_resonant_removal_7adic,lrc14_family_completeness_search,lrc14_perturbation_inf_search}_mac-mini-2026-06-16-S2.py + workflow angle scripts (lrc14_{7adic_archimedean_split,lll_shearer,abel_cesaro_levelmass,selberg_minorant,largespeed_decoupling,resonant_stranger_certificate,ocf_bridge_probe}) + .out; 07-reflections/why-every-lrc-positivity-method-fails-the-signed-integral-with-no-floor-mac-mini-S2.md; THM-522, HYP-2561/2562/2563, MISTAKE-075, OPEN-Q-104/097.

## Lead mac-mini-2026-06-16-S1: the octal lens on the H-spectrum + the triangular tile-count (T832, HYP-2558..2560)

**Status:** NEW probe past the kind-pasteur boat session (T830/HYP-2557, which mapped Fib/tri/square/prime + Heron Ã¢â‚¬â€ all re-verified, 3rd confirmation). The octal identity oddÃ‚Â²=8T_k+1 aimed at H(T) (always odd, RÃƒÂ©dei): at nÃ¢â€°Â¤6 the H-spectrum gaps {7,21,35,39} all AVOID residue 1 mod 8; the odd-square residue class is gap-free. Candidate necessary sieve for the forbidden-H frontier.
**Next (ranked):** (1) HYP-2558 Ã¢â‚¬â€ one Held-Karp sweep at **n=7 (and n=8 if feasible)**: is residue 1 mod 8 still gap-free? do 35,39 get realized? extend `H_mod8_octal_probe`. (2) HYP-2560 Ã¢â‚¬â€ extend the G_n/E_n invariant table to n=8 (m=21, first nontrivial Fibonacci tile-count); does anything special happen at nÃ¢Ë†Ë†{3,4,8,12}? (3) HYP-2559 Ã¢â‚¬â€ align THM-067 cÃ¢â€šÂ-vanishing Mersenne indices with the perfect-number tile-counts m=T_{2^k-1} at n=5,9,33. (4) reproduce {1,3,21,55} inside the V_4 (complementÃƒâ€”transpose) Burnside engine (`burnside_enum_v2.c`) Ã¢â‚¬â€ a positive run upgrades the Klein-four link from suggestive to structural.
**Sources:** 04-computation/{fib_tri_square_prime_heron,H_mod8_octal_probe}_mac-mini-2026-06-16-S1.py + .out; 07-reflections/the-octal-H-spectrum-and-the-triangular-tilecount-mac-mini-S1.md; T830/HYP-2557 + triangular_fibonacci_heron_boat_kps.py (prior art); THM-462/463 (H gap-freeness), THM-485/486/067/224, T232. (arXiv:1809.09936 = uniqueness Diophantine, per T830.)

## Lead mac-mini-2026-06-15-S4: LRC-14 via the lonely measure + the Riesz-product key (THM-515, T826, HYP-2540..2543, OPEN-Q-104)

**Status:** THM-515 (theta/lonely-measure form PROVED; additive-energy predictor verified; Riesz-product route set up + feasibility-probed; closed doors confirmed). inf L>0 OPEN. Adversarial check running.
**Next (ranked):** (1) OPEN-Q-104 Ã¢â‚¬â€ build the optimized dissociated Riesz product for the j=6 interior-drop extremizer cores and push Ã¢Ë†Â«MÃ‚Â·R below 1 (the certificate); port arXiv:2511.16636's construction. (2) Bonami hypercontractive level-k bound on the cores' AP relation lattice. (3) Gaussian-subordinated Selberg minorant (finitize + theta-positivity). (4) Linnik/ternary-form reduction (Pollock-style, THM-501).
**Sources:** 04-computation/{lrc14_theta_lattice,lrc14_riesz_product}_macmini_0615s4.py + .out; THM-515; reflection lrc14-is-the-lonely-measure-and-the-key-is-a-riesz-product.md; arXiv:2511.16636, Tao 1701.02048.

## Lead mac-mini-2026-06-15-S3: FKN Ã¢Å¸Â¹ Arrow Ã¢â‚¬â€ the tournament cube IS the social-choice cube (THM-512, T823, HYP-2534..2536, OPEN-Q-102)

**Status:** THM-512 proved/verified (Arrow/Condorcet bridge exact; c3=Guilbaud quadratic; MÃƒÂ¶bius sieve verified c3 nÃ¢â€°Â¤6). Builds on THM-511. Adversarial check running.
**Next:** (1) OPEN-Q-102 Ã¢â‚¬â€ is the OCF H=I(ÃŽÂ©,2) a noise-stability Stab_ÃÂ functional (mirroring Kalai ÃÂ=-1/3)? forbidden H Ã¢Å¸Â· forbidden Condorcet spectra. (2) the H multi-vertex MÃƒÂ¶bius sieve (HYP-2534, depth=band-limit D). (3) Friedgut/KKL: which invariants are juntas / which arc is decisive. (4) connect Guilbaud (arcsine/level-2) to the project's pi/Wallis constants (everything-is-the-triangle).
**Sources:** 04-computation/{fkn_tiling_cube,mobius_sieve_arrow}_macmini_0615s3.py + .out; THM-512; reflection fkn-arrow-and-the-tournament-as-a-vote.md.

## Lead mac-mini-2026-06-15-S1: baby Hodge = the det/permanent (P/#P) wall through the moment problem (THM-509, HYP-2526..2528, OPEN-Q-099, T820)

**Status:** THM-509 proved (Layer-1 det-side spectral-blindness rigorous) + all holes certified moment-interior (n=6,7); holes shown to be pure integrality gaps (flag-cut refuted). Adversarially verified.
**Next steps:**
1. **OPEN-Q-099 positivstellensatz:** prove no polynomial moment inequality (any degree, spectral or overlap) cuts a baby-Hodge hole for all n Ã¢â‚¬â€ making "hole = integer lattice point interior to the flag-feasible body, skipped by #P" a theorem.
2. **Generating function for the non-spectral defect:** does the necklace/zeta moment-cumulant pair (HYP-2526) give a clean handle on the non-spectral defect dimension A000009(n)-3 (THM-505)? The Witt/zeta ÃŽÂ (1-u^k)^{-W_k} and the restricted-partition ÃŽÂ _{k oddÃ¢â€°Â¥3}1/(1-x^k) are both Euler products Ã¢â‚¬â€ connect them.
3. **The c3-fiber hole structure:** the regular/near-regular score class is the richest hole source (c3=8 at n=6 Ã¢â€ â€™ 3 holes; c3=14 at n=7 Ã¢â€ â€™ 12). Characterize the holes as a function of the score sequence / the extremal (densest-cycle) fiber.
4. **Cross-link to dÃ¢Å Â¥H and the BSD/Hodge self-dual-middle reflection:** the det/per wall unifies them (d det-side, H per-side).
**Sources:** 04-computation/{baby_hodge_moment_region,baby_hodge_convex_certificate}_macmini_0615s1.py + .out; THM-509; reflection baby-hodge-is-the-det-permanent-wall-read-through-the-moment-problem.md.

## Lead mac-mini-2026-06-14-S1: dual-tower (self-dual) + LRC(14) singular-series structure (THM-503) + the Pascal-slope-d Pisot tower (T818/T819, HYP-2520..2525, OPEN-Q-097/098)

**Status:** session complete; THM-503 adversarially verified + adopted by the mesh (THM-504 builds on it); dual-tower answered; Pascal family identified+extended.
**Highest-value next steps:**
1. **OPEN-Q-097 (the LRC(14) prize):** prove inf L>0 over the dilated-AP cores. Now reframed (THM-503 + the concurrent THM-504) as a CROSS-LEVEL Abel-summation bound on the conditionally-convergent alternating-in-|T| series L=(6/7)^13+(6/7)^11 C_2Ã¢Ë†â€™(6/7)^10 ÃŽÂ£_3+Ã¢â‚¬Â¦ with a PÃƒÂ³lyaÃ¢â‚¬â€œVinogradov/Weil bound on each convergent signed sinc-lattice level ÃŽÂ£_k. Joint target with kind-pasteur's PZ thread.
2. **HYP-2409 indecomposability:** prove the skew tower's row code stays the INDECOMPOSABLE d+ (not e8Ã¢Å â€¢e8) for all k Ã¢â‚¬â€ the weight-4 support graph stays connected under doubling [[H,H],[Ã¢Ë†â€™HÃ¡Âµâ‚¬,HÃ¡Âµâ‚¬]]. Self-dual Type II half is provable + verified to order 64; indecomposability is the deep part (rides SO(32) forever).
3. **OPEN-Q-098 (gap-d tournament realization):** does a_d(n) count a natural gap-d tournament/staircase family for dÃ¢â€°Â¥3? Prove the realized 2Ã‚Â·Fib(m-2) circular-tournament count (d=2, S518). Define the "d-graded metagraph" whose H-level sizes are the slope-d ridge.
4. **The plastic-number bridge:** d=5 of the family (xÃ¢ÂÂµ=xÃ¢ÂÂ´+1) shares its root with Padovan (xÃ‚Â³=x+1, the monad/free-factorial THM-438 thread). Is there a slope-5 Ã¢â€ â€ monad tournament bridge, or pure algebraic coincidence?
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

## Lead kind-pasteur-2026-06-11-S3: the pentagonal product is a hub Ã¢â‚¬â€ random-sign Lyapunov ÃŽÂ³_pent, Euler-sign rigidity, the ÃŽÂ·Ã‚Â²Ã¢ÂÂ´ code-discriminant bridge, [72,36,16] localized (THM-487/488/489, HYP-2417..2423)

**Done:** THM-488 (ÃŽÂ³_pentÃ¢â€°Ë†0.206 new Lyapunov constant, validated vs Viswanath; Euler's signs the UNIQUE subexponential pattern = analytic shadow of the pentagonal #thm; IVT half proved, hard half certified on 1585 sets via argument principle). THM-489 (code discriminant PÃ¢â€šâ€šÃ¢â€šâ€ž = 16ÃŽÂ·Ã‚Â²Ã¢ÂÂ´ exactly; extremal correction cÃ¢â€šÂ(m)=Ã¢Ë†â€™42m proved; HYP-2420 MOS-mechanism CORRECTED to secular, n=3696 reproduced). THM-487 ([72,36,16] obstruction is code-combinatorial Ã¢â‚¬â€ WÃ¢â€šâ€¡Ã¢â€šâ€š all-positive, ÃŽâ€œÃ¢â€šâ€¡Ã¢â€šâ€š exists, Paley gauge stalls at d=12). Renumbered from 485/486 per claudebox-S5 first-come.
**Next:**
1. Prove the rigidity (HYP-2417): the uniform min-modulus bound on the sparse pentagonal polynomial P_S over a boundary circle Ã¢Å¸Â¹ RouchÃƒÂ© Ã¢Å¸Â¹ interior zero for all finite flip sets.
2. The ÃŽÂ·^{Ã¢Ë†â€™b} Lyapunov family: ÃŽÂ³ for general b interpolating partitions (b=1) and codes (b=24).
3. The ternary GleasonÃ¢â‚¬â€œPierce analog Ã¢â‚¬â€ extremal ternary is already negative at n=72 (Pierce); is there a ternary ÃŽÂ³-rigidity?
4. Sparse-lag extension of GoldsheidÃ¢â‚¬â€œZeitouni (arXiv:2505.00377) Ã¢Å¸Â¹ rigorous existence/positivity of ÃŽÂ³_pent.
5. HYP-2422: probe the random-extremal never-negative prefix tail at larger m (does deterministic stay the strict maximum-delay outlier?).

## Lead kind-pasteur-2026-06-11-S2: Gleason Ãƒâ€” tournaments closed (THM-481 joint w/ claudebox-S3), zigzag law opened (THM-483) Ã¢â‚¬â€ ErdÃ…â€˜sÃ¢â‚¬â€œMoser #1216 corrected

**Done:** THM-481 (merged per MSG-870): both Gleason Type II generators tournament-generated
(gÃ¢â€šâ€šÃ¢â€šâ€ž = C(I+S(PaleyÃ¢â€šâ€šÃ¢â€šÆ’)), ÃƒÂªÃ¢â€šË† = PaleyÃ¢â€šâ€¡); eQR identity proved then attributed Ã¢â‚¬â€ KimÃ¢â‚¬â€œSolÃƒÂ©, Des.
Codes Cryptogr. 49 (2008) Prop. 3 (found by the adversarial round's literature mandate);
Gleason-generation framing novel. THM-483: zigzag law trans(D(T)) = z(T); HYP-2360 +2
sandwich REFUTED by the alternating family A_l (ÃŽÂ´ = l unbounded, counterexample at n = 7,
one past the n=6 census); THM-455 addendum; ErdÃ…â€˜sÃ¢â‚¬â€œMoser bounded-increment route closed.
**Next:**
1. HYP-2413/HYP-2442: tower support gates - `trans(T127)=15` and `trans(T255)=23` are now exact, so the frontier is `T511`. Do not brute-force first; mine marked packets `X subset T127` with `trans(X)=15` and `q(X)=trans(D(D(X)))-30`. A packet with `q>=1` proves the predicted `trans(T511)>=31`; chain-only packets have `q=0`.
2. HYP-2409 claims 2Ã¢â‚¬â€œ3 (dictionary glue map; check claudebox THM-482's even-row hypothesis covers all tower levels).
3. Is RM(2,5) a tournament gauge of ANY order-32 skew-Hadamard? (finite sharp question, THM-481 remark).
4. GF(p) KimÃ¢â‚¬â€œSolÃƒÂ©: p-ary Paley tournament codes vs the GleasonÃ¢â‚¬â€œPierce self-dual rings.

## Lead mac-mini-2026-06-11-S1: REED-MULLER ON THE TILING CUBE Ã¢â‚¬â€ the blue code + the digit ladder (THM-477/478, HYP-2406..2408, T779)

**Status:** session complete; both theorems proved + adversarially verified (fresh conventions, zero corrections).
**What landed:** blue tilings = dual-containing linear code ker(1+sigma), self-dual defect = hypotenuse f, Plotkin=Mode-B recursion b(n)=b(n-2)+(n-2); digit ladder of H on the tiling cube: digit_0 = RM(0,m) (Redei), deg(digit_1) = D = 2*floor((n-1)/2) (bound PROVED via cycle-reversal pairing, equality n=4..7), digits >= 3 saturate; (D+1)-flat XOR annihilation = the literal RM(r)^perp = RM(m-r-1) statement for H mod 4. Digit-degree question for counting functions confirmed UNPUBLISHED (model case: Hamming weight digits = e_{2^i}, Canteaut-Videau).
**Next steps:**
1. **Prove HYP-2406 equality** (deg digit_1 = D for all n): exhibit a single ANF monomial of degree D with odd coefficient Ã¢â‚¬â€ candidate: the top-cycle monomial through the base path; likely a short argument on top of THM-478's pairing.
2. **HYP-2408 Ax shadow:** apply Ax divisibility (zeros of deg-D functions divisible by 2^(ceil(m/D)-1)) to joint digit level sets; target a coding-theoretic re-proof or strengthening of H Ã¢Ë†â€° {7,21}; compute achievable-residue tables of H mod 8, 16 per n.
3. **d_2 sequence 3,6,7,11** (n=4..7): extend to n=8 (2^21 tilings x H Ã¢â‚¬â€ feasible in C), find the law (HYP-2407).
4. **The extended-Hamming slot:** dual of the skew-tower code (RM(1,m)=Sylvester analog exits RM at order 16, THM-451) Ã¢â‚¬â€ what tournament object is RM(m-2,m)? Unexamined.
5. **Blue slice:** restrict digit functions to B_n; compute degrees in the sigma-invariant algebra; connect to SC/transpose-self classes (canon: transpose-self never pure black).
6. **Cite/absorb Sole-Zaslavsky 1994** (cosets of cocycle code = switching classes of signed graphs) into the THM-474 orbit Ã¢â‚¬â€ the tournament/base-path-gauge version remains the project's.
**Sources:** 04-computation/rm_digit_ladder_macmini_0611s1.py + thm477/478 adversarial rechecks; 05-knowledge/results/*0611*.out; 07-reflections/reed-muller-on-the-tiling-cube.md.

## Lead mac-mini-2026-06-10-S2: THE DETERMINANT LENS Ã¢â‚¬â€ d(T) = det(I+S)/2^(n-1) on the metagraph; tilings = switching classes (THM-468/472/473/474, HYP-2383..2389+2398..2400, T777, MISTAKE-071)

**Status:** session complete; theorems proved + adversarially verified; census exact nÃ¢â€°Â¤9; conjecture witnesses at n=13.
**What landed:** floor (d=1 Ã¢Å¸Âº local order, attribution Knuth/Babai-Cameron/Boussairi-DÃ¢â€šÂ), ceiling (DRT switching classes Ã¢Å¸Âº skew-Hadamard, verified zero-corrections), average (involutions/Hermite, attribution KMPRS LAA 707 (2025)), gauge theorem (tilings = switching classes), tournament Barba conjecture (OPEN-Q-058), mod-4 score law (HYP-2398), carousel takeover (HYP-2388/2399), d Ã¢Å Â¥ H + d-not-second-eigenvector (new metagraph coordinate), HYP-2312 universal form refuted (MISTAKE-071).
**Next steps (ranked):**
1. **Prove the tournament Barba bound** (OPEN-Q-058) Ã¢â‚¬â€ the genuinely new extremal problem; publishable with Klanderman et al. as the companion.
2. **Build the switching metagraph S_n** (A049313: 1,1,2,2,6,12,79): nodes = switching classes up to iso, per-node d + H-spectrum + SC composition; BabaiÃ¢â‚¬â€œCameron Thm 3.2 gives per-class tournament counts. New quotient layer above G_n/Z_2.
3. **Prove the mod-4 score law** (HYP-2398) via Pfaffian-minor parity expansion Ã¢â‚¬â€ clean THM candidate.
4. **OEIS submissions:** max det(I+S) 2,4,16,32,160,512,4096,8192 and d-version 1,1,2,2,5,8,32,32 (absent, checked); extend A334123 (labeled max-H counts: a(7)=240 confirmed; a(8),a(9) computable from census + |Aut| of the 6/1 max classes); extend A038375 to n=12+ via Pf=Ã‚Â±1-existential pruning (HYP-2312 amended) or circulant-first search.
5. **Klanderman arc-flip update on the metagraph:** det S' = det S(1+2SÃ¢ÂÂ»Ã‚Â¹_ij)Ã‚Â² assigns every wiggly line an exact determinant weight; flip-silent Ã¢Å¸Âº SÃ¢ÂÂ»Ã‚Â¹_ij Ã¢Ë†Ë† {0,Ã¢Ë†â€™1}. Compute the det-silent wiggly fraction per class Ã¢â‚¬â€ the d-analog of silent mutations.
6. **GrinbergÃ¢â‚¬â€œStanley mod-4 RÃƒÂ©dei (H Ã¢â€°Â¡ 1 + 2Ã‚Â·#odd-cycles mod 4, published 2023) = depth-1 OCF truncation** Ã¢â‚¬â€ formalize against THM-062/OCF; their U_D machinery vs kind-pasteur's THM-466 2-adic digit tower (flagged to kp in MSG).
7. **GLMY path homology per iso class** (literature: uncomputed for general tournaments; ÃŽÂ©Ã¢â€šâ€š spanned by transitive triangles + square differences) Ã¢â‚¬â€ new metagraph invariant; cross with THM-120's |Pf| phase separation.
8. **SchweserÃ¢â‚¬â€œStiebitzÃ¢â‚¬â€œToft arXiv:2510.10659 (RÃƒÂ©dei revisited, Feb 2026)** Ã¢â‚¬â€ read for the Dirac/Berge strengthenings; Satake 2025 (nearly-DRTs Ã¢â€ â€ almost difference sets, conditional on H-L Conjecture F) for the q Ã¢â€°Â¡ 1 mod 4 side.
**Sources:** 05-knowledge/results/{hadamard_det_census,hadamard_det_n9,circulant_odd_det,tournament_simplex,metagraph_det_gradient}_macmini_s2.out + thm472/thm473_adversarial_check.out + verify_hyp2312_and_mod4law_subagent_2026-06-11.out; reflection 07-reflections/the-determinant-lens-sgn-vs-chi-and-the-three-geometries.md.

## Lead mac-mini-2026-06-11-S2: BSD + HODGE merged into the self-dual-middle frame (T786, HYP-2420..2424, OPEN-Q-095, reflection the-self-dual-middle-where-parity-lives)

**Status:** session complete; reflection + 5 HYPs; rigorous spine adversarially verified; one duplicate (THM-490) caught & withdrawn.
**What landed:** (1) RIGOROUS SPINE Ã¢â‚¬â€ alternating pairing Ã¢Å¸Â¹ square (Pfaffian) is the shared symplectic mechanism behind |ÃÂ¨(E)|=Ã¢â€“Â¡ (Cassels), det(I+2A)=PfÃ‚Â² (THM-174), Weil/Hodge alternating forms; the square-defect is 2-ADIC in both worlds (Poonen-Stoll cÃ¢Ë†Ë†ÃÂ¨[2], the 8=2Ã‚Â³ in THM-442's 8Q) Ã¢â‚¬â€ the precise home of the user's "doubling." (2) FRAME Ã¢â‚¬â€ the self-dual-middle table (RH sÃ¢â€ â€1Ã¢Ë†â€™s / BSD root number / Hodge pÃ¢â€ â€q / Goldbach swap / project complement-ÃÆ’-SC), parity lives at each fixed middle; the one rigorous Ã‚Â½/doubling link is ÃŽÂ»(2n)=Ã¢Ë†â€™ÃŽÂ»(n). (3) MERGE Ã¢â‚¬â€ additive Goldbach (primes, segment) vs multiplicative Euler (strong-tournament H-primes); the s=2 polygonal segment bridges them; the project's PROVEN genus-2 H-semigroup R=odds\{7,21} (HYP-2271) is the tournament ÃŽÂ¶-side.
**Next steps:**
1. **OPEN-Q-095** (genuinely open): is there a finite abelian group / pairing attached to a tournament (ÃÂ¨-analog) whose order = det(I+2A) and whose alternating-vs-antisymmetric type (hence square-vs-2Ãƒâ€”square) is controlled by a tournament parity (SC / transpose-self / blue-black)? Compute the Smith normal form of I+2A and S; see whether SNF squares track det=PfÃ‚Â² and whether the 2-part carries a cÃ¢Ë†Ë†[2]-style defect. Cf. Klanderman et al. (SNF of skew D-optimal designs).
2. **HYP-2421** (polygon side = additive arity): make the sÃ¢â€ â€arity alignment a precise statement on the project's polygonalÃ¢â€ â€Goldbach ladder (HYP-1962); does s=2 segment / s=3 triangle / s=4 square map to a tournament/figurate hierarchy with a clean transition law?
3. **The root-number analog:** is there a tournament "root number" (a global Ã‚Â±1 = product of local Ã‚Â±1's over spine/ribs/sea or over strong components) whose sign predicts a parity (e.g. of some rank-like Q-grading)? The honest version must avoid the decorative local-global product; look for a genuine multiplicative Ã‚Â±1 over strong components paralleling w(E)=Ã¢Ë†ÂwÃ¡ÂµÂ¥.
4. **Cross-link** HYP-2424 to the existing band-gap/numerical-semigroup reflection (polarized-delta-fields-...-s699) Ã¢â‚¬â€ the genus-2 semigroup already has a physics (band gap) and NT (Frobenius) reading; add the Euler-product/Goldbach-duality reading.
**Sources:** 04-computation/h_realizability_goldbach_macmini_0611s2.py + .out; verify_moon_realizability_fresh_opus_0614.py/.out; reflection the-self-dual-middle-where-parity-lives.md.

## Lead kind-pasteur-2026-06-11-S1: ErdÃ…â€˜s 592 session 4 Ã¢â‚¬â€ the seam EXPLAINED (sum-free gradings), the t=7 wall, first Schipperus cutoffs, m=3 probe (THM-469/470/471, HYP-2390..2397)

**Done:** Picked up mac-mini-2026-06-10-S1's Next items 2/3/4. (3) SOLVED Ã¢â€ â€™ THM-469 (PROVED):
v_p level sets sum-free Ã¢Å¸Âº p=2 (parity closure); odd-p algebras die by 3-term-AP mono-collapse
(14-clause minimized UNSAT core at (3,4)); b=3 branching control Ã¢Å¸Â¹ the seam follows the SCHUR
ARITY, not branching (HYP-2390 CONFIRMED; THM-464 D slogan corrected); leading-digit refinement
rescues every odd p (HYP-2391 CONFIRMED: (sign,vÃ¢â€šÆ’,lead) SAT (3,4)/(3,5) verified). (2) ANSWERED
NEGATIVELY SO FAR Ã¢â€ â€™ THM-470: coarsening collapse (rungs must REFINE F2); dyadic 1-jet F2J dies
PER-SIZE at (3,7) Ã¢Å¸Â¹ F2 itself dies at t=7 (refinement) Ã¢Å¸Â¹ gap-determined algebras antitone in t:
decided alive tÃ¢â€°Â¤6/dead tÃ¢â€°Â¥7 Ã¢â‚¬â€ master run on the FULL invariant algebra at (3,7) [in flight at
write-up]; HYP-2396 candidate: R(n,2) = 2n+1 (3, 5, 7?; lower bounds Q(n,2n) SAT match at
n=1,2,3). (4) BUILT Ã¢â€ â€™ THM-471: B3 general-shape tuple grammar implemented + validated (m=1
reproduces THM-460 D); j1-march size explosion documented (BT(3,2) j1 = 3.5M leaves, vacuous
guard); first Schipperus-forced cutoffs ever computed: m=2 M=1 UNSAT at (2,2) and (2,3); m=3
probes at (2,2) in flight.
**Next:**
1. Harvest the in-flight runs (invariant wall at (3,7); F2X t=6,7; m=2 M=2-j0 sweep; m=3 M=1/M=2-j0 at (2,2); crossval).
2. If invariant dies at (3,7): prove/refute R(n,2) = 2n+1 (HYP-2396) Ã¢â‚¬â€ try the "+2 per level" induction via THM-459's L1-style row reduction; if it lives: extract the invariant witness's non-dyadic structure = the true next rung.
3. NON-invariant rungs (value-dependent features, e.g. vÃ¢â€šâ€š of coordinate values themselves) Ã¢â‚¬â€ the only route left to a t-uniform table if the wall is real.
4. m=3 vs m=2 differential at matched (s,c) with M=2-j0 and M=3 shapes; first Chang number (m=1 (3,3)) still open.
5. THM-459 hand closure; f_grid(4) optimality (inherited).

## Lead mac-mini-2026-06-10-S1: ErdÃ…â€˜s 592 session 3 Ã¢â‚¬â€ Q(3,5)/Q(3,6) SETTLED, the algebra ladder is strict, the seam is 2-adic at n=3; cubic lens placed (THM-464/465, T770, HYP-2373..2376)

**Done:** Q(3,5) = SAT and Q(3,6) = SAT, settled by explicit bi-dyadic (sign+vÃ¢â€šâ€š per coordinate) witnesses, independently verified Ã¢â‚¬â€ the 80k-clause timeout instance falls in 2.8s in the feature quotient (POKE Task 2 progressed, Task 1.3 advanced: R(3,2) > 6). Uniform-table SAT over t=4..7 is feature-UNSAT (0.3s) Ã¢Å¸Â¹ NO infinite (sign,vÃ¢â€šâ€š)-measurable strong witness: the ladder signs Ã¢Å Å  signs+vÃ¢â€šâ€š Ã¢Å Å  (open) is strict. Triadic control: sign+vÃ¢â€šÆ’ (equal-size algebra) is feature-UNSAT at (3,4) instantly Ã¢Å¸Â¹ the 2-adic seam is real at n=3; at n=2 ALL gradings share cutoff 5 (vÃ¢â€šÆ’ control). R_b(1) = R(3,b) proved+verified (classical Ramsey = height-1 row). Cubic lens (HYP-2376 census): cubes sum-free forever / non-Sidon from 1729 = first C4 of the cubic summand graph; signed from (3,4,5,6); cites kind-pasteur THM-462/463 + HYP-2370.
**Live at close:** batched (3,3) Chang (first Chang number) and ternary free sweep tÃ¢â€°Â¥7 Ã¢â‚¬â€ outputs stream to 05-knowledge/results/erdos592_{chang33_batched,ternary_seam}_macmini_s3.out.
**Next:**
1. Harvest the two live runs; R_3(2) exact; first Chang number.
2. CLIMB THE ALGEBRA LADDER: find the next rung past sign+vÃ¢â€šâ€š that admits a t-uniform table (candidates: unbounded-vÃ¢â€šâ€š with tail conditions, Larson partial-sum/scheme features, mixed row-grading Ãƒâ€” column-algebra) Ã¢â‚¬â€ a uniform rung + KÃƒÂ¶nig = constructive strong Specker at Ãâ€°Ã‚Â³.
3. Explain the vÃ¢â€šâ€š/vÃ¢â€šÆ’ asymmetry algebraically (parity closure odd+oddÃ¢â€ â€™even vs the mod-3 escape; THM-464 D open note).
4. mÃ¢â€°Â¥2 tower enumerator (THM-460 B3 grammar) Ã¢â‚¬â€ still the open-case probe (POKE Task 1.2).
5. THM-459 hand closure; f_grid(4) optimality.

## Lead mac-mini-2026-06-09-S2: ErdÃ…â€˜s 592 session 2 Ã¢â‚¬â€ Chang towers, dyadic witnesses, R(2,2)=5 lemma layer (THM-459/460, HYP-2363..2366)

**Source:** second 592 session; continues the S1 lead below.
**Done this session:** THM-459 (R(2,2)=5: lemma layer L1Ã¢â‚¬â€œL5 Ã¢â‚¬â€ doubly-dark clique, CÃ¢â€šâ€¦/triple dichotomy, trace cliques, composition-freeness Ã¢â‚¬â€ + machine closure; hand closure open); THM-460 (tower miniature for Ãâ€°^(Ãâ€°^m): full-type Ã¢Å¸Âº stacked towers PROVED, recursive binary shape grammar, KÃƒÂ¶nig bridge Ã¢Å¸Â¹ Chang/Schipperus FORCE finite tower-Ramsey cutoffs; m=1 sweep partial); invariant cutoff = free cutoff = 5 at n=2; DYADIC sufficiency (B_g through vÃ¢â€šâ€š(g)) verified through (3,4) Ã¢â‚¬â€ the 2-adic seam shows up inside 592 witnesses; f_grid = 1, 7, Ã¢â€°Â¤30; witness-tournament H flier (odd, transitive-leaning).
**HARVESTED post-close (both runs completed; results in 05-knowledge/results/):** Q(2,5) UNSAT independently RE-CONFIRMED (R(2,2)=5 double-certified by two complete-verifier implementations); Q(3,5) TIMEOUT at 3605s/80,111 clauses (free) and 1801s/62,854 (invariant) Ã¢â‚¬â€ UNDECIDED; Chang M=2: **(2,5) SAT with an 83-edge hard-fought witness (916s)** Ã¢â‚¬â€ the tower game survives at the exact ambient where the grid game dies, so Chang numbers sit strictly deeper than Specker numbers; (3,3) and all larger instances TIMEOUT Ã¢â‚¬â€ the first Chang number remains uncomputed (C2 guarantees it exists); M=3 vacuity guard: s Ã¢â€°Â¥ 3 required (height-3 grids need 3 split positions); Q(4,3) SAT (115 edges).
**Next steps:**
1. Decide Q(3,5) and the (3,3)/(2,6) Chang instances: needs symmetry breaking (row/value permutations), incremental cardinality, possibly parallel CEGAR Ã¢â‚¬â€ the instances are at the edge of naive encodings.
2. Implement the THM-460 B3 general-shape enumerator (towers of towers) Ã¢Å¸Â¹ run m=2 (Schipperus-forced cutoffs) and the **m=3 probe of the open $1000 case** Ã¢â‚¬â€ the designated centerpiece.
3. Extract a closed-form dyadic rule from the invQ(3,t)/dyadic witnesses (B_g via vÃ¢â€šâ€š(g) + within-row R) and try to prove it works for all t Ã¢Å¸Â¹ via THM-453 D1, a constructive strong form of Specker's Ãâ€°Ã‚Â³ negative.
4. Hand-close THM-459 (the two regimes are bounded case analyses; SAT core is small).
5. Witness-ensemble statistics for the H-flier before claiming anything (HYP-2364 note).

## Lead mac-mini-2026-06-09-S1: ErdÃ…â€˜s 592 Ã¢â‚¬â€ the R(n,t) tree-grid ladder and the rank-3 graded-relation frontier (THM-453, T768, HYP-2344..2346)

**Source:** ErdÃ…â€˜s 592 session (Ãâ€°^ÃŽÂ² Ã¢â€ â€™ (Ãâ€°^ÃŽÂ²,3)Ã‚Â²; open case ÃŽÂ±=Ãâ€°^(Ãâ€°Ã‚Â³) = exponent with three CNF summands).
**Status:** ACTIVE Ã¢â‚¬â€ exact finite results at nÃ¢â€°Â¤2, n=3 partially computed, Ãâ€°^(Ãâ€°^ÃŽÂ³) lift not started.
**What is established (THM-453):** witness frame (triangle-free graphs, no full-type independent set); grid characterization of full-type subsets; the compactness bridge Q(n,t)-SAT-Ã¢Ë†â‚¬t Ã¢Å¸Â¹ Ãâ€°^n Ã¢â€ â€º (Ãâ€°^n,3)Ã‚Â² (strong witness), hence positive relations force finite cutoffs R(n,2); computed **R(1,2)=3, R(2,2)=5 (exact)**, Q(3,3) SAT (35 edges), Q(3,4) SAT (346 edges, complete verifier, 692s); no pure order-pattern witness exists at Ãâ€°Ã‚Â³ (all 13 maximal triangle-free pattern sets grid-avoidable); the unique triangle-free pattern class at Ãâ€°Ã‚Â² is the SHIFT GRAPH.
**Next steps, in order of leverage:**
1. Decide Q(3,5) (was running at session close; CEGAR + complete verifier; needs either patience or a smarter verifier Ã¢â‚¬â€ bitset adjacency, subgrid search ordered by solver phase).
2. Hand-prove R(2,2)=5 (UNSAT certificate is small: 532 CEGAR clauses; extract a human Zarankiewicz/Mantel argument Ã¢â‚¬â€ would make THM-453 D+E a self-contained finite theorem with a clean constant under Specker's theorem).
3. Prove Q(3,t) SAT for all t by a UNIFORM rule (extract from the t=3,4 witnesses; the invariant quotient (R,{B_g}) format is the right ansatz) Ã¢Å¸Â¹ via THM-453 D1 a NEW constructive proof of Specker's Ãâ€°Ã‚Â³ Ã¢â€ â€º (Ãâ€°Ã‚Â³,3)Ã‚Â², SAT-discovered.
4. Lift F1 to Ãâ€°^(Ãâ€°^ÃŽÂ³): write the graded-relation algebra over the gap semigroup (CNF of ÃŽÂ³; summand count = grading rank); recover Schipperus's Ã¢â€°Â¤2-summand impossibility / Ã¢â€°Â¥4-summand constructions as the two regimes; instantiate rank 3 = the open $1000 case.
5. Import literature exacts: Darby (JCTB, "Negative partition relations for ordinals Ãâ€°^(Ãâ€°^ÃŽÂ±)"), Larson "pentagram" (JSL 2000), Schipperus APAL 2010 constructions; EHMR's presentation of Specker's witness Ã¢â‚¬â€ compare with the SAT-discovered structure.
**Files:** 04-computation/erdos592_*.py (6 scripts), 05-knowledge/results/erdos592_*.out, 03-artifacts/drafts/erdos-592-survey-and-reframes-macmini-s1.md, 07-reflections/the-additive-ladder-reaches-the-ordinals-macmini-s1.md. Mistakes: MISTAKE-066 (bridge direction), MISTAKE-067 (incomplete CEGAR verifier Ã¢â‚¬â€ read the witness).

## Lead kind-pasteur-2026-06-09-S2: ErdÃ…â€˜s problems through the doubling lens (THM-455/456/457, T769, HYP-2356..2361, MISTAKE-068/069)

**Status:** ACTIVE Ã¢â‚¬â€ all branches adversarially verified.
**ErdÃ…â€˜s 64 (#64, ErdÃ…â€˜sÃ¢â‚¬â€œGyÃƒÂ¡rfÃƒÂ¡s):** THM-456 (blowup spectrum law spec(G[KÃ¢â€šâ€š]) = gap-free [3, 2s(G)]; a single edge plants a twin CÃ¢â€šâ€ž Ã¢Å¸Â¹ blowups never counterexamples; TurÃƒÂ¡n corridor closed nÃ¢â€°Â¤9; all 71 CÃ¢â€šâ€ž-free ÃŽÂ´Ã¢â€°Â¥3 graphs n=10Ã¢â‚¬â€œ12 killed by forced CÃ¢â€šË†) + THM-457 (dyadic gate ladder: girth ladder 24 / 28 (NEW EXACT) / Ã¢â€°Â¤32 / >46 / 58; ladder principle: closing each dyadic gate inflates the next; dihedral 3-reflection Cayley family with dyadic spectrum exactly {4,32}; Exoo G78 reconstructed Ã¢â‚¬â€ contains C32, new beyond the 2014 paper). MISTAKE-069: McGee HAS 34 C8s (S710 enumeration-order artifact corrected).
**ErdÃ…â€˜sÃ¢â‚¬â€œMoser (#1216):** THM-455 Ã¢â‚¬â€ tower trans = 3,5,7,11 (extremality window 7Ã¢â‚¬â€œ31, closes by 63); sandwich trans(D(T)) Ã¢Ë†â€™ trans(T) Ã¢Ë†Ë† {1,2} over all 32768 n=6 tournaments (HYP-2360); ReidÃ¢â‚¬â€œParker FORCES the PaleyÃ¢â€šâ€¡ +2 exception (verified trans(D(PaleyÃ¢â€šâ€¡)) = 5); T31 ties published trans(QRÃ¢â€šÆ’Ã¢â€šÂ) = 7 (solver externally validated vs MomiharaÃ¢â‚¬â€œSuda Table 2).
**Next steps:** half-life ladder prediction at order 128 (HYP-2361); 2-adic character condition for the {4,32} dihedral family; classify alternating-chain +2 exceptions; hand trans data to mac-mini's 592 thread (THM-453 Part A bridge).

## Lead kind-pasteur-2026-06-09-S1: the skew-Sylvester doubling D(T) Ã¢â‚¬â€ Walsh/skew-Hadamard recursion on tournaments (THM-447/448/450..452/454, T767, HYP-2332..2337+2350..2355)

**Source:** user directive (n*2 doubling, three copies + one negated, skew-Hadamard normalization = tiling frame).
**Status:** ACTIVE Ã¢â‚¬â€ THM-447 proved + verified n=3..5; D(border(C3)) core Ã¢â€°â€¦ Paley T_7 (iso verified); branch fan-out in progress (H-formula, Mersenne tower, tiling ÃÆ’-eigenspaces, Hadamard equivalence, Clifford family classification).
**Key objects:** D(T) = [[M, M+I],[MÃ¢Ë†â€™I, Ã¢Ë†â€™M]]; spectral law M'Ã‚Â² = IÃ¢â€šâ€šÃ¢Å â€”(2MÃ‚Â²Ã¢Ë†â€™I) (Chebyshev TÃ¢â€šâ€š, ÃŽÂ»Ã‚Â²+1 doubles = Hadamard order); canonical Ham path (p, twin, reversed p); twin arcs saturate the anti-diagonal (hypotenuse) of the doubled staircase; tiling splits into ÃÆ’-symmetric copy blocks + ÃÆ’-antisymmetric cross block.
**Next steps:** (1) H(D(T)) / H(T[KÃ¢â€šâ€š]) functional in I(ÃŽÂ©(T),x) Ã¢â‚¬â€ would answer OPEN-Q-045 Q1; (2) identify tower DRT_15 and DRT_31 vs Paley; (3) ÃŽÂ©(D(T)) structure theorem; (4) engineering: skew-Walsh butterfly transform O(m log m).

---

## Lead monad-explorer-2026-06-07-S18: the two humps are one resurgent series Ã¢â‚¬â€ Watson bridge with explicit coefficients `bÃ¢Â±Â¼` (THM-438 ADD-16)

**Status:** PARTIALLY DONE (this session). The free factorial law's spectral tail `ÃÂ(x)~e^{1Ã¢Ë†â€™x}ÃŽÂ£bÃ¢Â±Â¼xÃ¢ÂÂ»ÃŠÂ²`
and moment ratio `A088368(k)/k!~eÃ‚Â·ÃŽÂ£bÃ¢Â±Â¼/(k)Ã¢Â±Â¼` are Watson-lemma images sharing the EXACT rationals
`b=1,2,10,178/3,1178/3,42494/15,Ã¢â‚¬Â¦` (`bÃ¢Â±Â¼=[tÃŠÂ²]exp(R(ÃÆ’(t))Ã¢Ë†â€™1)`, `t=ÃÆ’/Q(ÃÆ’)`). VERIFIED: tail vs
parametric density `1eÃ¢Ë†â€™16`; moments `kÃ¢â€°Â¤60`. Refines OEIS Kotesovec `a(n)~eÃ‚Â·n!`. The `bÃ¢Â±Â¼` are
Gevrey-1 (divergent) Ã¢Å¸Â¹ the `e`-overshoot on both sides is optimal truncation of one resurgent series.
**Next steps:** (1) **prove the Watson identity coefficient-by-coefficient** (leading order rigorous via
Carleman-determinacy + `e^{Ã¢Ë†â€™x}` tail; full series = standard "moment asymptotics from tail expansion"
uniform-remainder theorem, resurgent Ã¢Å¸Â¹ asymptotic only). (2) **closed form / GF for `bÃ¢Â±Â¼`** Ã¢â‚¬â€ the `bÃ¢Â±Â¼`
are built from the same Gompertz/`EÃ¢â€šÂ` `g`; is there an `EÃ¢â€šÂ`-type closed form, or an OEIS-negative
recurrence? (3) feed `bÃ¢Â±Â¼` back to the EDGE: the analogous `lnln` edge expansion (ADD-15 NEXT (3)) and
whether an analogous "edge bridge" links the `xÃ¢â€ â€™0` density to small-`k`/negative-`k` data.
**Files:** `04-computation/paley_starstar_tail_moment_watson_monad.py` (+`.out`); reflection
`07-reflections/the-two-humps-are-one-resurgent-series.md`; THM-438 ADDENDUM-16; HYP-2308/INDEX.

---

## Lead monad-explorer-2026-06-07-S17: the free factorial law's density is a closed-form parametric curve Ã¢â‚¬â€ push the crossing-`q` family `ÃŽÂ¼_q` (THM-438 ADD-15)
**Source:** THM-438 ADDENDUM-15, reflection `the-two-singularities-of-the-exponential-integral-shape-the-density.md`, HYP-2308.
**Status:** OPEN (new). ADD-15 gave the whole density as a parametric curve `x(u)=Ã¢Ë†â€™uÃ‚Â²g(u)`, `ÃÂ(u)=Ã¢Ë†â€™Im(u)/(Ãâ‚¬|u|Ã‚Â²)` on
`Im(uÃ‚Â²g(u))=0` (`u=Ã¢Ë†â€™1/G`, `g=eÃ¡ÂµËœEÃ¢â€šÂ(u)`), verified vs root-found ÃÂ to `1eÃ¢Ë†â€™12` and moments A088368 to `<0.3%`. The two
ends of the support = the two singularities of `EÃ¢â€šÂ` (log at `0` Ã¢â€ â€™ edge Ã¢Ë†Å¡log; cut `(Ã¢Ë†â€™Ã¢Ë†Å¾,0]` Ã¢â€ â€™ tail). The tail constant `e`
is DERIVED as `e^{R(0)}=e^{ÃŽÂºÃ¢â€šÂ}` via the Stokes term, with `ÃÂe^x=e^{R(G(x))}` overshooting `e` (resurgent hump = the
`m_k/k!Ã¢â€ â€™e` hump of MISTAKE-063). DONE from the S16 lead: the edge structure is now fully transparent from the
parametrization (`arg uÃ¢â€ â€™Ã¢Ë†â€™Ãâ‚¬/2`, `x~ÃŽÂµÃ‚Â²(ln(1/ÃŽÂµ)Ã¢Ë†â€™ÃŽÂ³)`); full `lnln` resummation still a clean finite task.
**Next steps:** (1) **`K_q(z)` / the parametric density of the crossing-`q` family `ÃŽÂ¼_q`** (the priority now). The
parametric machine `{x=K(w), ÃÂ=Ã¢Ë†â€™Im w/Ãâ‚¬ on Im K=0}` works for ANY `ÃŽÂ¼` with explicit `K`. If `K_q` has an `EÃ¢â€šÂ`-type
closed form, the two singularities should slide continuously across `0Ã¢â€°Â¤qÃ¢â€°Â¤1`: the `EÃ¢â€šÂ` log (free edge Ã¢Ë†Å¡log) weakening
into the classical bounded edge `Ã¢â€ â€™eÃ¢ÂÂ»Ã‚Â¹` + atom `eÃ¢ÂÂ»Ã‚Â¹ÃŽÂ´Ã¢â€šâ‚¬`. WHERE (in `q`) does the logÃ¢â€ â€™atom transition happen Ã¢â‚¬â€ `q=0Ã¢ÂÂº` or a
critical `q_c`? (2) **BelinschiÃ¢â‚¬â€œNica `B_t(ÃŽÂ¼_free)`** Ã¢â‚¬â€ `B_t` in closed `K` form; any named image. (3) Resum the full
edge `lnln` expansion (now a finite calc from the parametrization). (4) Off-diagonal `t(k,m)` THIRD deformation; `t(7,5)`;
HYP-2308 remainder (non-circulant DRT n=15).

## Lead monad-explorer-2026-06-07-S16: the free factorial law's edge is `Ã¢Ë†Å¡log`, not a constant Ã¢â‚¬â€ resum it, and find where the atom dissolves in `ÃŽÂ¼_q` (THM-438 ADD-14)
**Source:** THM-438 ADDENDUM-14, reflection `eulers-divergent-series-is-the-free-factorial-laws-r-transform.md`, HYP-2308.
**Status:** OPEN (new). ADD-14 gave the free factorial law's K-transform a CLOSED FORM `K(z)=Ã¢Ë†â€™(1/zÃ‚Â²)e^{Ã¢Ë†â€™1/z}EÃ¢â€šÂ(Ã¢Ë†â€™1/z)=Ã¢Ë†â€™(1/zÃ‚Â²)g(Ã¢Ë†â€™1/z)`
(`g`=Gompertz fn; Euler's divergent series resummed; verified `<7eÃ¢Ë†â€™17`), proved `K(Ã¢Ë†â€™1)=Ã¢Ë†â€™ÃŽÂ´` (Gompertz constant), and
showed the `xÃ¢â€ â€™0` edge has NO finite constant: `Ãâ‚¬ ÃÂÃ¢Ë†Å¡x Ã¢â€ â€™ Ã¢Ë†Å¡(ln|G|Ã¢Ë†â€™ÃŽÂ³) ~ Ã¢Ë†Å¡(Ã‚Â½ln(1/x))` (corrects ADD-12 `1/Ãâ‚¬` / ADD-13 `0.4Ã¢â‚¬â€œ0.6`).
**Next steps:** (1) **Resum the full edge asymptotic.** `Ãâ‚¬ ÃÂÃ¢Ë†Å¡x=Ã¢Ë†Å¡(ln|G|Ã¢Ë†â€™ÃŽÂ³)` is exact; solve `rÃ‚Â²x=ln rÃ¢Ë†â€™ÃŽÂ³` with `2ÃŽÂ¸Ã¢â€ â€™Ãâ‚¬` to higher
order for the complete `xÃ¢â€ â€™0` expansion of the density (the `lnln`/`ÃŽÂ³` corrections). Closed-form `K` makes this tractable.
(2) **`K_q(z)` for the crossing-`q` family.** Does `K_q` have its own exp-integral closed form? The classical end (`q=1`)
has a BOUNDED density (`Ã¢â€ â€™eÃ¢ÂÂ»Ã‚Â¹`) plus an atom `eÃ¢ÂÂ»Ã‚Â¹ÃŽÂ´Ã¢â€šâ‚¬`; the free end (`q=0`) has no atom and a `Ã¢Ë†Å¡logÃ‚Â·x^{Ã¢Ë†â€™1/2}` divergence.
Where in `0<q<1` does the atom dissolve / the divergence switch on Ã¢â‚¬â€ at `q=0Ã¢ÂÂº` or a critical `q_c`? (Ties to S15.)
(3) **BelinschiÃ¢â‚¬â€œNica `B_t(ÃŽÂ¼_free)`** is now computable in closed form (`K` explicit) Ã¢â‚¬â€ any named `B_t` image gives the
wild end a second analytic handle.

## Lead monad-explorer-2026-06-07-S15: the BelinschiÃ¢â‚¬â€œNica `B_t` semigroup vs the crossing-`q` family `ÃŽÂ¼_q` (THM-438 ADD-13)
**Source:** THM-438 ADDENDUM-13, reflection `the-two-endpoints-are-bercovici-pata-partners.md`, HYP-2308.
**Status:** OPEN (new). ADD-13 established the two named endpoints A000262/A088368 are BercoviciÃ¢â‚¬â€œPata partners
(shared cumulants `n!`, classicalÃ¢â€ â€free CP of `ÃŽÂ½=e^{-x}dx`), joined by a positive-definite measure family `ÃŽÂ¼_q`
(`m_k(q)=ÃŽÂ£_Ãâ‚¬ q^{cr(Ãâ‚¬)}Ã¢Ë†Â|B|!`, Hankel `D_n(q)` all-nonneg `q`-coeffs, `q=0` free, `q=1` classical).
**Next steps:** (1) The BelinschiÃ¢â‚¬â€œNica `B_t` transform interpolates `Ã¢Ë†â€”`/`Ã¢Å Å¾` analytically with `B_1=ÃŽâ€º` Ã¢â‚¬â€ is the
crossing-`q` family `ÃŽÂ¼_q` the SAME interpolation? If so it pins `ÃŽÂ¼_q`'s Cauchy transform to a known closed form
and may yield the free density's exact edge constant (numerics give `Ã¢â€°Ë†0.4Ã¢â‚¬â€œ0.6`, NOT `1/Ãâ‚¬`; see ADD-13 part 6).
(2) The off-diagonal `t(k,m)` columns are now ruled out as BOTH the crossing-`q` triangle AND the rate-marked
`N(k,j)` (ADD-12) Ã¢â‚¬â€ they need a THIRD deformation, transverse to both named axes. Find it.
(3) `ÃŽÂ¼_q`'s own cumulants (q-deformed) Ã¢â‚¬â€ does `ÃŽÂ¼_q` stay compound-Poisson-like for `0<q<1`?

## Lead monad-explorer-2026-06-07-S12: PROVE the diagonal `t(k,k)=A088368(k)=ÃŽÂ£_{NC(k)}Ã¢Ë†Â|B|!` (the cleanest THM-438 sub-target)
**Source:** THM-438 ADDENDUM-10, MISTAKE-063, reflection `the-cancellation-runs-between-two-named-endpoints.md`.
**Status:** VERIFIED kÃ¢â€°Â¤7 + OEIS match (A088368 = Callan "partitions of [n] into sets of noncrossing lists", arXiv:0711.4841). Closed form `t(k,k)=ÃŽÂ£_{Ãâ‚¬Ã¢Ë†Ë†NC(k)}Ã¢Ë†Â_B|B|!` VERIFIED kÃ¢â€°Â¤7. Asymptotic `~eÃ‚Â·k!` (Kotesovec) Ã¢â‚¬â€ ADD-9's "`Ã¢â€°ÂeÃ‚Â·m!`" RETRACTED (the ratio overshoots e, peaks at m=8, descends back; MISTAKE-063).
**Next steps:** (a) **Prove the bijection** diagonal even-series patterns (= max-cycle-rank = doubled plane trees, Euler tour visits v `deg(v)` times, weight `Ã¢Ë†Â(degÃ¢Ë†â€™1)!`) Ã¢â€ â€ NC(k) with block-factorial weight `Ã¢Ë†Â|B|!`. Finite, number-theory-free; would upgrade `t(k,k)=A088368(k)` to PROVED and give `h_m(m)Ã¢â€ â€™e` as a corollary. (b) The two open handoffs (#1 `(k)_m|t(k,m)`; #2 `g_m(Ã¢Ë†â€™1)=(Ã¢Ë†â€™1)^m(2^mÃ¢Ë†â€™1)`) live at the TAME end Ã¢â‚¬â€ the wild diagonal is now closed-form, so attack the tame end via the `(1+s)|Q_m` line-parity involution. (c) Off-diagonal columns/residues are all OEIS-NEGATIVE (P_m(1)=1,3,20,181; subdiag 9,72,580,4845; col3 13,72,230,560; unsigned rowsum 1,4,23,160,1262,10944; `ÃŽÂ£_{NC}Ã¢Ë†Â(|B|Ã¢Ë†â€™1)!=1,2,6,23,105,553,3311`) Ã¢â‚¬â€ do NOT re-hunt them in OEIS. (d) Still need `t(7,5)` (the k=7 background run DIED at k=6): a core-aware enumerator validated vs kÃ¢â€°Â¤6.

## Lead monad-explorer-2026-06-07-S5: prove `(Ã¢Ëœâ€¦Ã¢Ëœâ€¦)` via walk-induced ribbon genus + A215257 bijection
**Source:** THM-438 ADDENDUM-3, HYP-2308, reflection `the-drt-engine-is-S-squared-equals-J-minus-nI-the-catalan-is-genus-zero.md`.
**Status:** OPEN. `(Ã¢Ëœâ€¦Ã¢Ëœâ€¦) ÃŽÂ£_{even-series ÃÆ’}ÃŽÂ¼(0ÃŒâ€š,ÃÆ’)=(Ã¢Ë†â€™1)^kC_k` VERIFIED kÃ¢â€°Â¤5; value = free cumulant of the two-point law. Leading-order engine now `SÃ‚Â²=JÃ¢Ë†â€™nI` (DRT-universal, number-theory-free, VERIFIED pÃ¢â€°Â¤43).
**Next steps:** (a) build the rotation system the Euler walk `x_0..x_{2k}` induces on `G_ÃÆ’`, compute ribbon genus per pattern (kÃ¢â€°Â¤5 enumerated on disk), test `ÃŽÂ£_{genus 0}ÃŽÂ¼=(Ã¢Ë†â€™1)^kC_k` & `ÃŽÂ£_{genus>0}ÃŽÂ¼=0` (naive non-crossing-PARTITION version is REFUTED Ã¢â‚¬â€ do not retry). (b) write the explicit bijection even-series patterns Ã¢â€ â€ indecomposable deque-sortable permutations (OEIS A215257); its sign structure may prove `(Ã¢Ëœâ€¦Ã¢Ëœâ€¦)`. (c) HYP-2308 remainder: expander-mixing `o(n^{k+1})` bound, tested on a verified non-circulant DRT n=15. (d) extend cycle-rank triangle to k=6,7 (needs enumerator smarter than all-set-partitions) to pin its recursion.

## Priority opus-2026-05-16-S1: TRRT and All-0 Staircase

### INV-191: H=63 Unlocks at n=8 via Complete Conflict Graph
**Source:** opus-2026-05-29-S8; exact census opus-2026-05-29-S10
**Status:** EXACT at n=8; single-core mechanism identified in S11
**What:** HYP-1754 ("H=63 is universally forbidden") is refuted. A concrete n=8 tournament has H(T)=63 by both DP and direct permutation enumeration. Its odd-cycle conflict graph ÃŽÂ©(T) has 31 directed odd cycles and is complete, so OCF gives H=I(K31,2)=1+2Ã‚Â·31=63. S10 upgraded this to a finite theorem: among all 6880 n=8 isomorphism classes, exactly two have H=63; both have |Aut|=1, score sequences (1,2,2,3,3,5,6,6) and (1,1,2,4,4,5,5,6), and ÃŽÂ©(T)=K31. This explains how 63 bypasses the old disconnected K3-factor obstruction: it realizes 63 through a complete ÃŽÂ©, not through K3Ã¢Å â€2K1.
**S11 update (opus-2026-05-29-S11):** Both H=63 classes are **single-core**: every odd directed cycle contains one vertex, and deleting that core vertex leaves a transitive 7-vertex tournament. Their transitive-insertion signatures are `1001100` and `1100110`; the weighted signature count
`r(s)=ÃŽÂ£_{i<j, s_i=1, s_j=0} 2^{max(j-i-2,0)}` equals 31, matching the number of odd cycles. A complete-ÃŽÂ© census over isomorphism classes n=3..8 shows cycle counts r=3 and r=10 are absent whenever ÃŽÂ© is complete; a signature target search finds r=3 and r=10 absent for all single-core signatures up to m=16, while r=31 first appears at m=7. This gives a new focused H=21 lens: prove the r=10 gap in the single-core family, then handle non-core complete ÃŽÂ© and non-complete ÃŽÂ±-tuples separately.
**S12 update (opus-2026-05-29-S12):** Reframed the H=63 mechanism as a projection defect. The two H=63 classes are exact old-projection kills: deleting the core vertex loses 31/31 directed odd cycles and leaves `H(T-v)=1`, `alpha(T-v)=[1,0]`. A core-stratified complete-ÃŽÂ© census through n=8 confirms r=3 and r=10 are absent in every core stratum, and r=31 occurs only in core-size-1 classes. The single-core target search was extended to mÃ¢â€°Â¤40: r=3 and r=10 remain absent; r=31,42,63 appear with simple linear count laws after their first occurrence.
**S5 update (opus-2026-05-30):** The applied residue/phase contrast table classifies H=63 as the clean exact-kill endpoint of the standard anomaly set: `rho=1`, max-loss residue `(0,0,0)`, residue rank 0. This sharpens the next step: treat the single-core signature count as a finite-state target problem for absent values `r=3,10`, separate from THM-025-style near-kill failures.
**Next:**
  1. Prove the single-core signature formula and the persistent gaps r=3,10.
  2. Classify complete-ÃŽÂ© tournaments with empty cycle-family core; compare their r-support to the single-core support.
  3. Revisit H=21 by separating the complete-ÃŽÂ© case ÃŽÂ©=K10 from the non-complete ÃŽÂ±-vector cases.
  4. Locate the H=63 classes in the merged metagraph/principal-line coordinates.
  5. Prove or refute the projection-kill/near-kill hypothesis by scanning more real-root failures and complete-ÃŽÂ© classes.
**Files:** `04-computation/h63_counterexample_audit_s8.py`, `04-computation/h63_n8_isoclass_census_s10.py`, `04-computation/omega_extreme_fingerprints_s11.py`, `04-computation/projection_defect_bridge_s12.py`, `05-knowledge/results/h63_counterexample_audit_s8.out`, `05-knowledge/results/h63_n8_isoclass_census_s10.out`, `05-knowledge/results/omega_extreme_fingerprints_s11.out`, `05-knowledge/results/single_core_signature_targets_s11.out`, `05-knowledge/results/projection_defect_bridge_s12.out`, THM-344, MISTAKE-050, HYP-1757, HYP-1758, HYP-1760, HYP-1761, HYP-1762.

### INV-189: Real-Rootedness of I(ÃŽÂ©(T), x) for All Tournaments (TRRT)
**Source:** opus-2026-05-16-S1
**Status:** STALE/REFUTED AS UNIVERSAL Ã¢â‚¬â€ see THM-025
**What:** The universal conjecture "I(ÃŽÂ©(T), x) has all real, negative roots for every tournament T" is false. Canon THM-025 gives an n=9 counterexample with score sequence [1,1,3,4,4,4,6,6,7] and I(ÃŽÂ©,x)=1+94x+10xÃ‚Â²+xÃ‚Â³, violating Newton's inequality. What remains interesting is to characterize the large real-rooted subclass and understand why failures are rare.
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
**What:** H values of all-0 interleaved staircase at n=2k. 5,29,233 are Markov numbers (breaks at k=5: 2489=19Ãƒâ€”131). c3=k(k-1) confirmed through k=12.
**Full sequence (k=2..12):** 5, 29, 233, 2489, 33773, 562685, 11222321, 262755369, 7110764837, **219612027389, 7658921303353**
**Proved:** # directed 3-cycles = k(k-1) exactly.  The staircase specialization was formalized sorry-free in `eliott-monad/math-lean` commit `b5ffcde` as `Math.Tournaments.allZeroStaircase_threeCycleCount`; THM-410 records the more general interval-reversal counting theorem.
**Growth ratios:** r(k)=H(k)/H(k-1): 5.80, 8.03, 10.68, 13.57, 16.66, 19.94, 23.41, 27.06, 30.88, 34.87. Differences Ã¢â€°Ë† +3.8/step (slowly increasing). Empirical: r(k) Ã¢â€ â€™ 3k asymptotically; deficit 3kÃ¢Ë†â€™r(k) peaked at kÃ¢â€°Ë†6 then decreases toward 0.
**Ruled out:** No order-2/3 linear recurrence. Markov fails all consecutive triples. OEIS search inconclusive (oeis.org blocked all fetches; web search for specific terms found no match Ã¢â‚¬â€ likely novel sequence).
**Next:** (1) Submit to OEIS. (2) Formalize THM-410's general interval-reversal theorem in Lean. (3) Find algebraic structure (product formula, generating function). (4) Try k=13 if feasible (~5 min runtime with array.array). (5) Investigate why r(k)Ã¢â€ â€™3k Ã¢â‚¬â€ is there a combinatorial explanation?

---

## Priority S24: New Leads from Session S24

### INV-184: Stanley-Stembridge PROVED (Hikita, arXiv:2410.12758)
**Source:** opus-2026-04-05-S24 web search
**Status:** NEW Ã¢â‚¬â€ needs integration
**What:** e-positivity of chromatic symmetric functions for unit interval graphs proved by Tatsuyuki Hikita (Oct 2024, revised Dec 2025). Since Mitrovic-Stojadinovic (arXiv:2506.08841) connects Redei-Berge to chromatic functions, and Stanley-Stembridge concerns e-positivity, this potentially constrains the h-positivity of U_T.
**Next:** Read Hikita's proof. Check implications for our U_T framework (THM-062/063). Does h-positivity of U_P follow?

### INV-185: Mitrovic NC Deletion-Contraction (arXiv:2504.20968)
**Source:** opus-2026-04-05-S24 web search
**Status:** NEW
**What:** The original Redei-Berge function does NOT satisfy deletion-contraction. Mitrovic introduces a noncommutative analogue that DOES. This enables inductive proofs Ã¢â‚¬â€ exactly what we need for Claim A-type results.
**Next:** Read the NC deletion-contraction construction. Can it give a new proof of Claim A? Does it yield new tournament invariants?

### INV-186: Real-Rootedness of I(ÃŽÂ©(T), x)
**Source:** opus-2026-04-05-S24 computation
**Status:** PROVED for n Ã¢â€°Â¤ 8, DISPROVED universally at n=9 (THM-025)
**What:** All zeros of the independence polynomial of the odd-cycle conflict graph are real and negative for nÃ¢â€°Â¤8 by claw-freeness + Chudnovsky-Seymour. At n=9, THM-025 gives a counterexample. The right problem is now structural characterization, not universal proof.
**Next:** Characterize the n=9 failure mechanism and identify sufficient tournament conditions for real-rootedness.

### INV-187: E[H(T)] = n!/2^{n-1}
**Source:** opus-2026-04-05-S24 computation
**Status:** PROVED (linearity of expectation)
**What:** Clean closed form for average Hamiltonian paths. W(n) = n! Ãƒâ€” 2^{C(n-1,2)}. Probably known but not explicitly in our bibliography.
**Next:** Check literature for this result. Is it in Moon's book?
**EXTENSION (monad-explorer-2026-06-07, HYP-2307):** the Paley/circulant maximizer's RATIO over this average, `R(p)=H(T_p)/(p!/2^{pÃ¢Ë†â€™1})`, Ã¢â€ â€™ **e** Ã¢â‚¬â€ PROVEN via a character-sum cluster expansion (`R=E_ÃÆ’[Ã¢Ë†Â(1+Ãâ€¡(d_k))]Ã¢â€ â€™exp(ÃŽÂ£a_L)`, only the cherry cluster survives, `a_2=Ã¢Ë†â€™Ãâ€¡(Ã¢Ë†â€™1)=1`). Universal across circulant tournaments (cherry weight = tournament condition `g` odd). **SUB-LEMMA CLOSED Ã¢â€ â€™ THM-438 (monad-explorer-2026-06-07, same day):** `a_{2k}=0 Ã¢Ë†â‚¬kÃ¢â€°Â¥2` PROVEN UNIFORMLY (no per-k Weil): `B_L=0` Ã¢Å¸Â¹ `A_L=Ã¢Ë†â€™ÃŽÂ£`coincidence-patterns; no-leaf forces `VÃ¢â€°Â¤2k`; only the single `V=2k` pattern (one even cycle) needs Weil, all others `o(p^{2k})` trivially. **Plus the CATALAN LAW** `A_{2k}=C_k p^{k+1}+O(p^{k+1/2})` (C_k=Euler tours of plane trees=moment-method tree-walks). See `04-computation/paley_cluster_{sharp_order,catalan}_monad.py`, reflection `the-paley-cluster-integrals-are-catalan-numbers-tree-walks-and-the-moment-method.md`, OPEN-Q-013.
**RATE DIAGNOSTIC (monad-explorer, for handoff #2 / the p=31,43,47 compute node):** exact `R(p)` (p=3..23) give `(eÃ¢Ë†â€™R)Ã‚Â·p = 2.15,2.23,3.07,3.63,3.71` (climbing, decelerating Ã¢â‚¬â€ consistent with `eÃ¢Ë†â€™R~K/p`, `KÃ¢â€ â€™Ã¢â€°Ë†3.8`, i.e. `R=e(1Ã¢Ë†â€™C/p)`, `CÃ¢â€°Ë†1.4`) BUT `(eÃ¢Ë†â€™R)Ã‚Â·Ã¢Ë†Å¡p = 1.24,0.84,0.92,0.83,0.77` is actually FLATTER Ã¢â‚¬â€ so the 5 points do NOT cleanly distinguish `O(1/p)` from `O(1/Ã¢Ë†Å¡p)`. The a_4-sector cluster argument predicts a `+O(1/p)` term (`~pÃ‚Â·A_4/p^5 = 2/p+Ã¢â‚¬Â¦`), favoring `1/p`; settle + pin `C` at pÃ¢â€°Â¥31. This is HYP-2307 handoff #2 (the smooth analytic Paley signature) Ã¢â‚¬â€ now a *prediction to test*, not a blind extrapolation.
**RATE RESOLVED Ã¢â€ â€™ 1/p (monad-explorer 3rd session, MISTAKE-060 + THM-438 ADDENDUM):** the cluster-integral error term is `A_{2k}=C_k p^{k+1}+O(p^k)` Ã¢â‚¬â€ i.e. `O(p^k)`, NOT the `O(p^{k+1/2})` originally claimed. Verified: `(A_4Ã¢Ë†â€™2p^3)/p^2` is STABLE (Ã¢â€°Ë†Ã¢Ë†â€™7.1Ã¢â‚¬Â¦Ã¢Ë†â€™7.8Ã¢â€ â€™Ã¢â€°Ë†Ã¢Ë†â€™8) while `/p^{2.5}` drifts to 0. So the analytic correction to `R` is relative `O(1/p)`, confirming `R=e(1Ã¢Ë†â€™C/p+Ã¢â‚¬Â¦)` over `O(1/Ã¢Ë†Å¡p)`. (The Ã¢Ë†Å¡p-vs-1/p ambiguity above was a finite-data artifact; the Catalan-skeleton error term settles it analytically.) `C` still to be pinned at pÃ¢â€°Â¥31 Ã¢â‚¬â€ the `R=e(1Ã¢Ë†â€™C/p)` ansatz is now justified, `CÃ¢â€°Ë†1.4` from current 5 points.

### INV-188: I(ÃŽÂ©(T), x) as Tournament Invariant
**Source:** opus-2026-04-05-S24 computation
**Status:** ESTABLISHED
**What:** The full polynomial I(ÃŽÂ©(T), x) is a strictly finer tournament invariant than H(T). Distinguishes tournaments at n=6 that share the same H. Real-rooted decomposition: H = ÃŽÂ (1 + 2/r_i) where -1/r_i are the zeros.
**Next:** Classify which iso classes are distinguished. Connection to Redei-Berge function?

---

## Priority A: Key structural questions (OCF PROVED by Grinberg-Stanley)

### INV-150: Simplicial RÃƒÂ©dei Theorem (THM-220) Ã¢â‚¬â€ PROVED for all n Ã¢â€°Â¥ 4
**Source:** opus-2026-03-15-S90
**Status:** COMPLETE. sim_H Ã¢Ë†Ë† {0,1} for all tournaments on n Ã¢â€°Â¥ 4 vertices.
**Key results:** Algebraic proof (Key Lemma + case analysis), near-transitive construction (H=2^{n-2}+1), OCF proof (ÃŽÂ© = K_{2^{n-3}}), ÃŽÂ²Ã¢â€šÂ=1 for all near-transitive.
**Char poly:** p(ÃŽÂ») = ÃŽÂ»^{n-2}(ÃŽÂ»Ã‚Â²+1) - (1+ÃŽÂ»)^{n-2} (binomial coefficients!).
**Next:** Write as standalone paper. This is publishable.

### INV-151: The Cayley Monad DÃ¢â€šâ€ž Framework
**Source:** opus-2026-03-15-S90c/d/f
**Status:** ESTABLISHED. QÃ¢ÂÂ´=Id generates DÃ¢â€šâ€ž of order 8 acting on Ã¢â€žâ€šPÃ‚Â¹.
**Key results:** Cross-ratio of Q-orbit of 2 = 2 (the fugacity). Fixed points Ã‚Â±i. Bloch sphere connection: DÃ¢â€šâ€ž = Pauli group mod center. QÃ¢ÂÂ´(M)=M on transfer matrix. Spectral zeta ÃŽÂ¶_M(-3)=7, ÃŽÂ¶_M(-5)=21 (forbidden H!).
**Next:** Investigate DÃ¢â€šâ€ž action on tournament invariants. Possible paper.

### INV-152: The Ãâ€ž-Ãâ€  Clock and Quasicrystal Structure
**Source:** opus-2026-03-15-S90h
**Status:** DISCOVERED. arg(ÃŽÂ»_c)/Ãâ‚¬ Ã¢â€°Ë† ln(2) (4 sig figs, NOT exact).
**Key results:** Gear ratio Ã¢â€°Ë† ln(2), period Ã¢â€°Ë† 2/ln(2), Tr mod 8 period = 8 exactly. Bott periodicity connection. Information interpretation: 1 bit per half-turn.
**Next:** Investigate why the approximation is so good. Connection to Rauzy fractal.

### INV-153: Tournament Equidecomposability (Hilbert's Third Problem analog)
**Source:** opus-2026-03-15-S90k
**Status:** PARTIALLY PROVED. ÃŽÂ²Ã¢â€šÂ as Dehn invariant. (H,ÃŽÂ²Ã¢â€šÂ) complete at n=5.
**Key results:** 8 equidecomposability classes at n=5. Near-transitive = regular tetrahedron (non-decomposable). Within each (H,ÃŽÂ²Ã¢â€šÂ) class, I(ÃŽÂ©Ã¢â€šÆ’,x) constant.
**Next:** Verify at n=6,7. Prove or disprove completeness for general n.

### INV-154: The Golden Shadow f(n) = (n-2+Ã¢Ë†Å¡(nÃ‚Â²+4))/2
**Source:** opus-2026-03-15-S90o/p, kind-pasteur-S112
**Status:** EXPLORED. CF = [n-2; n,n,n,...]. Pell norm fÃ‚Â·f'=-n. n-tribonacci family.
**Key results:** Unifies golden ratio, Ã¢Ë†Å¡2, silver ratio in one formula. f(n) satisfies ÃŽÂ»Ã‚Â²=(n-2)ÃŽÂ»+n, transfer matrix [[n-2,n],[1,0]]. Memory correction T_n-M_{n-2} maximal at n=3.
**Next:** Investigate the n-tribonacci family more. Connection to Pisot substitution conjecture.

### INV-155: Tournament Wick Rotation and Spin-1 Ising
**Source:** opus-2026-03-15-S90g/i
**Status:** DISCOVERED. arctanh(2) = log(3)/2+iÃâ‚¬/2 Ã¢â€ â€™ complex temperature.
**Key results:** H(T) = Wick-rotated partition function. ln(2) = renormalized criticality. Arrow's theorem = arctanh pole. Discriminant ÃŽâ€(x) = 4x(xÃ‚Â²-11x-1).
**Next:** Formalize the Ising correspondence. Yang-Lee zero connections.

### INV-141: T_19 Full Omega Dims Ã¢â‚¬â€ Degrees 9-18 Pending (Engineering Priority)
**Source:** kind-pasteur-2026-03-10-S54
**Status:** PARTIALLY COMPUTED. Degrees 0-8 done: [1,9,72,540,3753,23832,136260,688266,2987622].
**What:** Complete the T_19 Omega dim sequence to degrees 0-18. Degree 9+ requires C/C++ or LinBox.
**Key constraint:** chi per eigenspace = 1 (expected). Partial chi through m=8 is 2415061. Remaining (m=9..18) must sum to -2415060 alternating.
**Blocking issue:** Python dict-based pivot storage needs ~TB RAM for degree 9 (pivot density grows: max 10Ã¢â€ â€™26Ã¢â€ â€™95Ã¢â€ â€™455 at degrees 5Ã¢â€ â€™6Ã¢â€ â€™7Ã¢â€ â€™8).
**Next steps:**
  1. Implement in C/C++ (or use SageMath/LinBox) for degree 9+
  2. Or prove palindrome conjecture: Omega_{18-m} = Omega_m (gives free half via symmetry). T_7 IS palindromic; T_11 is NOT. T_19 unknown.
  3. Use chi constraint to check partial palindrome: if Omega_9 Ã¢â€°Ë† 12.95M (extrapolating ratio 4.34Ã¢â€ â€™3.6), check if sum through m=18 = 1.
**Priority:** MEDIUM (engineering interest, pattern understanding)
**Scripts:** t19_omega_dims_sparse.py, results at 05-knowledge/results/t19_omega_dims_sparse.out

### INV-142: Engineering Product: mod_rank Library (PyPI Target)
**Source:** kind-pasteur-2026-03-10-S53 (synthesis), S54 (update)
**Status:** MOSTLY DONE. mod_rank_library.py has core functions. 24 pytest tests written (test_mod_rank_library.py, all pass). Needs: README, benchmark vs dense/scipy, create setup.py.
**What:** General-purpose small-prime modular rank library. Key functions: gauss_rank_uint8, gauss_rank_nullbasis_uint8, certified_rank, betti_number_from_boundary_ranks. Memory table: T_11 deg 9 (52550Ãƒâ€”15745) = 6.6 GB int64 Ã¢â€ â€™ 827 MB uint8.
**Engineering value:** Useful for any combinatorics/topology computation over finite fields with large sparse matrices.
**Next steps:** Write README, benchmark vs dense/scipy, create setup.py for PyPI

### INV-143: Engineering Product: circulant_homology Python Module
**Source:** kind-pasteur-2026-03-10-S54
**Status:** COMPLETE. circulant_homology.py has CirculantHomology and PaleyHomology classes. Fixed and verified.
**What:** Clean API for computing Omega dims and Betti numbers of circulant tournaments. Uses sparse column reduction (THM-125 eigenspace identity gives nÃƒâ€” speedup). Verified for T_3, T_7, T_11.
**S55 fixes:** (1) betti_numbers() had wrong formula Ã¢â‚¬â€ fixed to use correct boundary map ranks via eigenspace computation. Verified T_3=[1,1,0], T_7=[1,0,0,0,6,0,0] exactly. (2) Added 27 pytest tests (test_circulant_homology.py, all pass). (3) Added caching for omega_basis_k and face_data for performance.
**New finding (HYP-453):** T_7 eigenspace structure: k=0 = H_0 only; k=1..6 = one H_4 generator each. T_11 eigenspace: all non-trivial homology at k=0.
**monad-compute (2026-06-03) from-scratch re-verification (`04-computation/verify_t11_betti_s_monad.py`, results `..._s_monad.out`/`.._NOTES.md`):**
  - ÃŽÂ© dims `[1,5,20,70,205,460,700,690,450,180,30]` RE-CONFIRMED from scratch (use_cache=False), Ãâ€¡=1. Raw |A_m| path counts = `[1,5,25,110,430,1430,3970,8735,14395,15745,8645]`. root-field prime=23.
  - Boundary ranks rank(d_m^(k)): **k=0** `[0,0,5,15,55,150,305,390,300,150,30,0]` (388s); **k=1** `[0,1,4,16,54,151,309,390,300,150,30,0]` (382s).
  - Per-eigenspace Betti: k=0 Ã¢â€ â€™ `[1,0,0,0,0,5,5,..]`; k=1 Ã¢â€ â€™ `[..,0,1,0..]` (+1 to ÃŽÂ²_6 only).
  - **Structural finding (REFINES HYP-453):** ÃŽÂ²_5=5 lives entirely at k=0; ÃŽÂ²_6 is distributed Ã¢â‚¬â€ k=0 carries 5, each non-principal eigenspace carries +1. Predicts ÃŽÂ²_6 = 5+10Ã‚Â·1 = **15** (matches cached). HYP-453 "all T_11 homology at k=0" holds for ÃŽÂ²_5 but NOT ÃŽÂ²_6.
**monad-compute-2026-06-03-S2 Ã¢â‚¬â€ FULLY COMPLETE from scratch (all 11 eigenspaces):**
  - **Cache-bloat fixed:** clearing `_omega_basis_cache` per eigenspace held timing flat at 384Ã¢â‚¬â€œ413 s/eigenspace (no more >12 min slowdown). Script now persists each eigenspace to `verify_t11_betti_s_monad_ranks.json` (resumable) + auto-commit.
  - **All k=1..10 have IDENTICAL boundary ranks** `[0,1,4,16,54,151,309,390,300,150,30,0]`; k=0 `[0,0,5,15,55,150,305,390,300,150,30,0]`.
  - **ÃŽÂ² = `[1,0,0,0,0,5,15,0,0,0,0]` CONFIRMED.** Per-eigenspace ÃŽÂ²_6 = `[5,1,1,1,1,1,1,1,1,1,1]` Ã¢â€ â€™ 5+10Ã‚Â·1=15; ÃŽÂ²_5 = `[5,0,...]` Ã¢â€ â€™ 5 (k=0 only). The "+1 to ÃŽÂ²_6 per non-principal eigenspace" pattern holds for ALL 10.
  - **Euler check clarification:** the script's first "OVERALL MISMATCH" was a FALSE ALARM Ã¢â‚¬â€ it compared Ãâ€¡_Betti=11 against single-copy Ãâ€¡_Omega=1. By THM-125 each of the n=11 eigenspaces carries a FULL copy of ÃŽÂ©_m, so the correct identity is Ãâ€¡_Betti = nÃ‚Â·Ãâ€¡_Omega = 11Ã‚Â·1 = 11. Ã¢Å“â€œ Fixed in the script. The Betti numbers were always correct (reproduce the library's official `betti_numbers()` formula exactly).
**STATUS: ÃŽÂ²_5=5, ÃŽÂ²_6=15 mechanically re-verified from scratch across all 11 eigenspaces. INV-143 T_11 Betti re-verification CLOSED.** Remaining open engineering item: C/LinBox reimplementation for routine degree-9+ re-verification (cf. INV-141).

### INV-192: Engineering Product: Odd-Cycle Disjointness Features for Tournament TDA
**Source:** opus-2026-05-29-S11
**Status:** PARTIALLY IMPLEMENTED (opus-2026-05-29-S15)
**What:** S11's H=63 / THM-025 comparison suggests that practical tournament fingerprints should expose the disjointness geometry of odd cycles, not only H, scores, and Betti numbers. Candidate features: `alpha_vector(Omega)`, cycle-family core size, complete-ÃŽÂ© flag, disjoint-pair count, independent-triple supports, single-core signature, `r_core(s)`, and deletion profile `(H(T-v), |Omega(T-v)|, complete?)`.
**S12 update:** Extend this feature layer to projection defects: support count vs support excess, max support multiplicity, deletion loss profile `(lost, kept, loss_frac, alpha(T-v))`, exact-kill/near-kill flags, and even-graph projection weight/degree sequence for odd n.
**S15 update:** `04-computation/tournament_tda.py` now computes exact `omega_*` features for nÃ¢â€°Â¤9: `omega_alpha_1`, `omega_alpha_2`, `omega_alpha_3`, complete-ÃŽÂ© flag, disjoint-pair density, support count/excess, max support multiplicity, cycle-family core size, projection-kill vertex count, and max deletion-loss fraction. These are included in the flat ML feature vector and demo output is saved in `05-knowledge/results/tournament_tda_omega_features_s15.out`.
**Why it matters:** These features separate two phenomena that H alone compresses: (1) H=63 unlocks through a complete-core ÃŽÂ©=K31; (2) real-rootedness fails through a no-core, highly concentrated independent triple. For ranking data, this distinguishes "all inconsistency localized at one pivot" from "three disjoint inconsistency groups with lopsided coupling."
**Next steps:**
  1. Benchmark the implemented features on synthetic rankings, sports/election data, and attention-derived tournaments.
  2. Use them as prefilters before expensive full ÃŽÂ© or path-homology computations.
  3. Add odd-n even-graph projection signatures and deletion `alpha(T-v)` profiles as a second feature layer.
**Files:** `04-computation/tournament_tda.py`, `04-computation/omega_extreme_fingerprints_s11.py`, `04-computation/projection_defect_bridge_s12.py`, `05-knowledge/results/tournament_tda_omega_features_s15.out`, `05-knowledge/results/omega_extreme_fingerprints_s11.out`, `05-knowledge/results/projection_defect_bridge_s12.out`, `07-reflections/omega-extremes-as-cycle-disjointness-axis.md`, `07-reflections/projection-defect-as-common-residue.md`.

### INV-193: Projection-Defect Axis Across ÃŽÂ©, Even Graphs, and Path Homology
**Source:** opus-2026-05-29-S12
**Status:** NEW META-STRUCTURAL LEAD
**What:** Several mature threads may be instances of the same question: what residue survives a forgetful projection? S12 compared three projections:
  1. Directed odd cycles Ã¢â€ â€™ vertex supports (multiplicity defect).
  2. Vertex deletion / old-coordinate projection (cycle loss, projection kill).
  3. Tournament orientation Ã¢â€ â€™ degree-even cycle-space graph, well-defined only at odd n.
**Key examples:** H=63 is an exact old-projection kill: deleting the core removes all 31 odd cycles. THM-025 is a near-kill: one vertex supports 92/94 odd cycles, but the two surviving old cycles carry `alpha=[1,2,1]`, enough to keep the real-root failure alive. Paley T7 and interval T7 have the same 36 odd-cycle supports but different support-excess and even-graph projection weights. Path-homology HYP-408/ghost cycles asks the same structural question in chain-complex language: when do through-v-only cycles become boundaries after old projection?
**Engineering angle:** Add projection-defect features to Tournament TDA: max deletion cycle-loss fraction, old-projection kill vertices, support multiplicity defect, and odd-n even-graph projection signature. These are cheap fingerprints that may prefilter root failures, localized inconsistency, and homology anomalies.
**Next steps:**
  1. Scan known THM-025-like non-real-root examples for high max deletion loss.
  2. Compare projection-defect statistics against beta_3/beta_4 path-homology anomalies.
  3. For odd n, correlate even-graph projection fibers with ÃŽÂ© support-multiplicity defect.
  4. Extend the S15 `tournament_tda.py` `omega_*` feature block with explicit deletion-residue alpha profiles and even-graph projection signatures.
**Files:** `04-computation/projection_defect_bridge_s12.py`, `05-knowledge/results/projection_defect_bridge_s12.out`, `07-reflections/projection-defect-as-common-residue.md`, HYP-1760, HYP-1763, T282, T283.

### INV-194: Merged Tiling Bucket Constraints
**Source:** kind-pasteur-2026-05-29-S5
**Status:** THEOREM + BOOLEAN-CUBE LEAN SPECIALIZATION COMPLETE + OPEN STRUCTURAL EXTENSION
**What:** Treat quotient maps out of the tiling cube as bucket maps. THM-346 proves the general half-line balance law for any quotient `q: Q_m -> B` and any mask family `M`: `2*self_b + incident_cross_b = |q^{-1}(b)|*|M|`. THM-345 specializes this to `pi: Q_m -> G_n/Z_2`, proving that bucket size parity detects SC/NS type exactly: SC buckets are odd, NS buckets are `2 mod 4`, and `sum_M B_M=2^m`. For every Hamming layer `d`, the ordered transport matrix `W_d(M,N)` is symmetric, has row sums `B_M*C(m,d)`, has even diagonal, and therefore has forced cross-outflow parity. Lucas' theorem says the active parity layers are exactly the binary submasks of `m=C(n-1,2)`.
**S6 Lean update (kind-pasteur-2026-05-29-S6):** The local good-cut bucket gap was strengthened in `TournamentH7.GoodCuts`: nonempty good-cut support is equivalent to the existence of an upward tile, any upward tile forces at least two distinct good cuts, and every tiling satisfies `goodCutCount = 0 Ã¢Ë†Â¨ 2 Ã¢â€°Â¤ goodCutCount`. The `TournamentH7.Verify` audit confirms this core depends only on Lean foundations.
**S1 Lean update (kind-pasteur-2026-05-30-S1):** Added THM-348 / `TournamentH7.BucketBalance`, a generic axiom-free Lean proof of the oriented finite-set balance `|selfHalf|+|crossHalf|=|fiber|*|moves|`, plus the zero-cross closure criterion. This isolates the remaining Lean bridge for THM-346: prove the fixed-point-free involution pairing that turns internal half-lines into unordered self-lines counted twice.
**S2 Lean update (kind-pasteur-2026-05-30-S2):** Added THM-350 / the unordered abstract layer in `TournamentH7.BucketBalance`: the partner map `(x,u)->(step u x,u)` preserves internal half-lines for involutive moves, fixed-point-free moves forbid self-partners, and the unordered balance follows from `Even selfHalf.card`. Also added the concrete staircase-to-tournament bridge and top-good-cut/SC audit via `StaircaseConnectivity`.
**Codex 2026-05-30 orbit update:** The generic orbit-parity bridge is now Lean-proved. `BucketBalance.even_card_of_fixedPointFree_involutiveOn` proves that any finite fixed-point-free involution has even cardinality, and `BucketBalance.unordered_balance_of_involutive_fixedPointFree` removes the separate evenness assumption from the abstract unordered balance.
**Opus 2026-05-30 Boolean-cube update:** THM-351 closes the remaining Boolean-mask bridge in Lean. `BoolCube`, `xorMask`, `xorMask_involutive`, `xorMask_fixedPointFree_of_nonzero`, and `unordered_balance_boolCube_masks` specialize the abstract THM-350 layer to finite Boolean cubes with nonzero xor masks. This gives a reusable formal model for tiling-cube mask families before quotient-specific tournament maps are attached.
**Why it matters:** This turns the merged metagraph into a constrained reversible transport system, not just an unweighted graph. It also corrects the old S205 "merged tiling excess" narration: merged buckets still partition the fixed-base cube exactly.
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
- Main theorems: CÃ¢â€ â€™_5^{1,2} and CÃ¢â€ â€™_n^{1,s} with sÃ¢â€°Â 2 have ÃŽÂ²_1=1,ÃŽÂ²_2=1, otherwise ÃŽÂ²_1=1.
- Symbol matrix M_m(t) reduces rank to evaluating Laurent polynomials at roots of unity.
- Stability theorem (Thm 1.4): for large primes pÃ¢Ë†â€°Q+(S), Betti numbers stabilize.
- DOES NOT compute Betti numbers for Paley T_p (connection sets of size (p-1)/2).
- DOES NOT have Omega palindrome or top vanishing results for arbitrary connection sets.
**Impact:** The symbol-matrix framework is the right tool for T_p, but they leave Paley application open. The Stability Theorem (Thm 1.4) implies our pattern ÃŽÂ²_6(T_p)=p-1 holds for ALL large pÃ¢Ë†â€°Q+(QR_p). Computing Q+(QR_p) would be the key step.
**Next step (UPDATED S52):** COMPLETED. Tang-Yau symbol matrix applied to T_7 (deg 2-5) and T_11 (deg 2-5). Key discovery: Q+(QR_p) is EMPTY (not just "avoids p-th roots" Ã¢â‚¬â€ ALL t values work equally). Symbol matrix M_m(t) is CONSTANT (THM-125: proven algebraically). Eigenspace identity is trivial consequence. New open: prove Q+(QR_p) empty algebraically for all p.
**S71c UPDATE:** Full HTML paper read. Confirms our eigenspace approach matches their Prop 4.1. Their Conjecture 4.8 (H_m=0 for mÃ¢â€°Â¥3, no-wrap-around S) does NOT apply to Paley (QR has wrap-around). Our HYP-824 Betti concentration at d=m,m+1 is genuinely new. Key open: prove ÃŽÂ²_d^(0)=0 for 3Ã¢â€°Â¤dÃ¢â€°Â¤m-1 (the intermediate vanishing).
**Priority:** CLOSED for computational part. Proving intermediate vanishing is Priority A.

### INV-136: Schweser-Stiebitz-Toft (arXiv:2510.10659): Redei's Theorem Revisited (Oct 2025)
**Source:** Web research kind-pasteur-2026-03-10-S50
**Status:** INVESTIGATED (kind-pasteur-2026-03-12-S55).
**What:** 6-page paper unifying three equivalent stronger forms of RÃƒÂ©dei's theorem: (i) RÃƒÂ©dei's stronger theorem (adding undirected vertices creates even path count), (ii) Berge's theorem (complementary mixed graphs have same parity), (iii) Dirac's theorem (Hamiltonian paths through non-oriented edges are even). All three are interconnected.
**Relevance:** Directly relevant to H(T) parity. The mixed-graph generalization provides a structural lens complementary to OCF. THM-016/017 (our even-odd split) and DC (THM-082/083) could potentially be re-derived via Berge/Dirac. Possible route to extend OCF to mixed digraphs.
**Action:** Add to TANGENTS as "Mixed graph OCF extension via Schweser-Stiebitz-Toft stronger RÃƒÂ©dei".

### INV-137: Satake (arXiv:2502.12090): Cyclotomic Nearly-Doubly-Regular Tournaments (Feb 2025)
**Source:** Web research kind-pasteur-2026-03-10-S50
**Status:** DEFINITIVELY RESOLVED (kind-pasteur-2026-03-12-S56). Satake NDRTs do NOT maximize H.
**What:** For prime powers q Ã¢â€°Â¡ 5 (mod 8), cyclotomic tournament CT_q is NDR iff q = sÃ‚Â² + 4. When true: full adjacency spectrum computed explicitly (eigenvalues (q-1)/2 and (-1 Ã‚Â± iÃ¢Ë†Å¡(q Ã¢Ë†â€œ 2Ã¢Ë†Å¡q))/2). Under Hardy-Littlewood conjecture F: infinitely many such q. This is the n Ã¢â€°Â¡ 1 (mod 4) analog of Paley (n Ã¢â€°Â¡ 3 mod 4).
**RESOLUTION (S56 satake_analysis_ext.py, exhaustive at n=13):**
- q=5: H_sat=15 = H_max (trivially tied Ã¢â‚¬â€ ALL 4 circulants tie at n=5)
- q=13: H_sat=3,703,011, rank 40/64 Ã¢â‚¬â€ FAR from maximum (gap=8,164, HYP-456 REFUTED)
- Maximizer at n=13: cyclic interval S={7,...,12}, H=3,711,175 (unique, rank 1/64) Ã¢â‚¬â€ HYP-455 CONFIRMED
- NDR property does NOT predict H-optimality for qÃ¢â€°Â¡5 mod 8
**Pattern discovered (HYP-455):** At Paley primes pÃ¢â€°Â¡3 mod 4 Ã¢â€ â€™ Paley maximizes. At qÃ¢â€°Â¡5 mod 8 primes Ã¢â€ â€™ cyclic interval maximizes.
**Next step:** Verify cyclic interval pattern at q=29 (needs sparse/C++ Held-Karp). q=29 currently hits MemoryError in Python.

### INV-138: Ren (arXiv:2504.15126): Path Independence Complexes of Digraphs (Apr 2025)
**Source:** Web research kind-pasteur-2026-03-10-S50
**Status:** INVESTIGATED (kind-pasteur-2026-03-12-S55). HIGH PRIORITY Ã¢â‚¬â€ closely adjacent to our work.
**What:** Studies "path independence complexes" of digraphs Ã¢â‚¬â€ simplicial complexes whose faces are vertex sets with no directed path between any two. Uses GLMY infrastructure. Main results: (i) canonical embeddings from independence complex (undirected) into path independence complex (digraph), (ii) these are ÃŽÂ£_k-equivariant and isometric giving double-parametrized persistent homology, (iii) Shannon capacity consequences.
**NOT the same as Omega(T):** Omega(T) is the conflict GRAPH of directed cycles (edges = two cycles sharing a vertex). The path independence complex is a SIMPLICIAL COMPLEX on vertices with different independence condition. However, both are topological structures related to the tournament and use GLMY machinery.
**Key connection:** The embedding theorem (independence complex Ã¢â€ â€™ path independence complex) might explain the tight relationship between Omega(T) and tournament path homology. Our beta_2=0 result and seesaw mechanism (THM-095) may follow from an embedding theorem structure.
**Next step:** Read sections on Betti number bounds and explicit computations. Compare beta numbers from Ren's framework to our beta_2=0, beta_3 onset results. This is potentially a theoretical unification.

### INV-139: Tang-Yau (arXiv:2402.05682): Cellular Homology of Digraphs (Feb 2024)
**Source:** Web research kind-pasteur-2026-03-10-S50
**Status:** NEW Ã¢â‚¬â€ background.
**What:** Introduces cellular chain complexes for digraphs using admissible minimal paths. Connects to GLMY under strongly regular conditions.

### INV-140: Hepworth-Roff (arXiv:2404.06689): Magnitude-Path Spectral Sequence (Apr 2024)
**Source:** Web research kind-pasteur-2026-03-10-S50
**Status:** NEW Ã¢â‚¬â€ background.
**What:** The MPSS connects magnitude homology (page 1) to path homology (page 2). Bigraded path homology satisfies Kunneth, Mayer-Vietoris, excision. Could provide extra structure for tournament path homology via the bigraded version.

### INV-144: Awan-Bernardi B-polynomial for Digraphs (arXiv:1610.01839)
**Source:** Web research kind-pasteur-2026-03-14-S69
**Status:** NEW Ã¢â‚¬â€ high priority connection
**What:** Awan and Bernardi defined a Tutte polynomial for directed graphs satisfying deletion-contraction B(D) expressed via B(D\a) and B(D/a). Published JCTB 2020. Our F_T(x) = F_{T\e}(x) + (x-1)*F(T/e,x) is structurally the same type of relation. The B-polynomial is a 3-variable polynomial detecting acyclicity, strong connectivity, and directed paths.
**Next step:** Read paper, check if F(T,x) is a specialization of their B-polynomial.

### INV-145: Sazdanovic-Yip Categorification of Chromatic Function (arXiv:1506.03133)
**Source:** Web research kind-pasteur-2026-03-14-S69
**Status:** NEW Ã¢â‚¬â€ highest priority for categorification
**What:** Sazdanovic and Yip categorified the chromatic symmetric function using Khovanov-style techniques. Combined with Mitrovic-Stojadinovic chromatic-Redei-Berge bridge (arXiv:2506.08841), this suggests categorifying F(T,x) via long exact sequences in homology, producing "Khovanov homology for tournaments."
**Next step:** Read and check if the deletion-contraction for F(T,x) lifts to a long exact sequence in GLMY path homology.

### INV-146: Asao Magnitude-Path Spectral Sequence (arXiv:2201.08047)
**Source:** Web research kind-pasteur-2026-03-14-S69
**Status:** NEW Ã¢â‚¬â€ connects magnitude and path homology
**What:** Asao proved magnitude homology and GLMY path homology are pages of the same spectral sequence. Published Bull. London Math. Soc. 2023. Magnitude homology categorifies magnitude (Hepworth-Willerton, arXiv:1505.04125), analogous to Khovanov categorifying Jones. This creates an indirect chain: Khovanov Ã¢â€ â€ magnitude homology Ã¢â€ â€ path homology.
**Next step:** Understand what the spectral sequence looks like for tournaments.

### INV-147: Hepworth Reachability Homology (arXiv:2312.01378)
**Source:** Web research kind-pasteur-2026-03-14-S69
**Status:** NEW Ã¢â‚¬â€ unifies magnitude and path homology
**What:** Hepworth defined reachability homology of digraphs, unifying magnitude and path homology. Satisfies homotopy invariance, Kunneth, excision, Mayer-Vietoris. Published IMRN 2025. Potentially the "right" homology theory for tournaments.
**Next step:** Compute reachability homology for small tournaments, compare to GLMY.

### INV-050: Fourier Decomposition Proof of OCF Ã¢â‚¬â€ OCF PROVED AT n=5 AND n=7
**Source:** opus-2026-03-06-S11b (continued^7, ^8)
**Status:** OCF PROVED AT n=5 AND n=7 via Fourier decomposition. All identities at both n proved.
**What:** OCF decomposes into independent degree-homogeneous Fourier identities:
- **Fourier Homogeneity Theorem:** w_{n-1-2k} is a homogeneous polynomial of degree 2k in centered edge variables s_e = A_e - 1/2.
- **Degree 0:** Trivial (expected values). PROVED for all n.
- **Degree 2:** Proved via proportionality constants c_{2j+1} = C(n-3,2j-2)*(2j-2)!/2^{2j-2}. PROVED for all n.
- **Degree n-1:** Proved via path-cycle bijection: w_0 = 2*[deg-(n-1) of t_n]. PROVED for all n.
- **Degree 4 at n=7:** PROVED via counting lemmas. The degree-4 Fourier space is 2-dimensional:
  - Type P: 5-vertex spanning paths (coefficients Ã‚Â±1 in t5_d4)
  - Type Q: 6-vertex disjoint PÃ¢â€šâ€š pairs (coefficients Ã‚Â±1 in ÃŽÂ±Ã¢â€šâ€š_d4)
  - [deg-4 of tÃ¢â€šâ€¡] = (1/2)Ã‚Â·[deg-4 of tÃ¢â€šâ€¦] + [deg-4 of ÃŽÂ±Ã¢â€šâ€š] (EXACT, verified all 5985 monomials)
  - wÃ¢â€šâ€š/4 = 3Ã‚Â·[deg-4 of tÃ¢â€šâ€¦] + 6Ã‚Â·[deg-4 of ÃŽÂ±Ã¢â€šâ€š] (EXACT, verified all 5985 monomials)
  - Counting: cÃ¢â€šâ€¦=2, cÃ¢â€šâ€¡=4, paths=12 (type P); cÃ¢â€šâ€¡=8, paths=24, aÃ¢â€šâ€š=1 (type Q)
**Key insight:** Fourier supports of tÃ¢â€šâ€¦_d4 and ÃŽÂ±Ã¢â€šâ€š_d4 are DISJOINT (spanning paths vs PÃ¢â€šâ€š pairs), reducing the identity to independent counting arguments.
**Scripts:** `04-computation/fourier_homogeneity.py`, `fourier_degree2_identity.py`, `ocf_fourier_proof_framework.py`, `degree4_identity_n7.py`, `degree4_proof_n7.py`
**CRITICAL FINDING (opus-2026-03-07-S37):** At n=9, the degree-4 Fourier space has dimension **>> 200** (no saturation at 200 random tournaments with 300 probe monomials). Within the P4 type alone, coefficients are NOT proportional (117/126 vertex sets have non-constant |coeff|). The |coeff| values range from 1 to 27. This means the n=7 "two type" decomposition DOES NOT generalize. The Fourier approach is **infeasible at n=9** for middle degrees.
**Scripts:** `04-computation/degree4_n9_rank.py`, `degree4_n9_rank2.py`, `degree4_n9_saturation.py`
**Next step:** (1) The Fourier proof cannot extend to n>=9 for middle degrees. Focus on algebraic approaches (OCF already proved by Grinberg-Stanley). (2) The degree-0, degree-2, and degree-(n-1) identities still hold for all n and have clean proofs. Can they be combined differently?

### INV-123: THM-086 Universal Taylor Zeros mod 3 Ã¢â‚¬â€ PROOF SKETCH COMPLETE
**Source:** kind-pasteur-2026-03-07-S37
**Status:** PROOF SKETCH COMPLETE. Verified n=5-10, inductive structure identified.
**What:** c_j(T) = 0 mod 3 for all tournaments T on n vertices and all j < val(n), where val(n) = 2*floor((n-1)/2). This means (x-1)^{val(n)} | F(T,x) mod 3. For n odd, F(T,x) mod 3 is determined by a SINGLE parameter alpha = c_{n-1}(T) mod 3.
**Proved cases:** j=0,1,2 (THM-085, algebraic). j=3 (palindrome + THM-085). j>=4 (DC induction + palindrome, verified computationally).
**Key corollary:** Eulerian conjecture: 3|A(n,k) => 3|F_k(T) for all T. Follows from (x-1)-adic valuation of A_n(x) mod 3 being exactly val(n).
**What remains:** The "almost-tournament claim" Ã¢â‚¬â€ c_j(T\e) = 0 mod 3 for j < val(n)-1 Ã¢â‚¬â€ needs formal proof, likely via nested DC induction. Verified exhaustively at n=5, sampled at n=6-8.
**Scripts:** `04-computation/thm086_verify.py`, `dc_induction_proof.py`, `c4_induction_test.py`, `taylor_cj_mod3_analysis.py`, `eulerian_zeros_from_palindrome.py`
**Next step:** (1) Prove almost-tournament claim algebraically (N_uv formula reduces it to Taylor zeros of the "adjacent pair" polynomial). (2) Extend to mod 9. (3) Mod p for p>=5 INVESTIGATED (S38): universal zeros match Eulerian val for n >= p+2 but Eulerian conjecture FAILS for p>=5 (multiple free parameters).

### INV-124: THM-094 F_k mod 2 Tournament-Independent Ã¢â‚¬â€ PROOF SKETCH COMPLETE
**Source:** kind-pasteur-2026-03-07-S38
**Status:** PROOF SKETCH COMPLETE. Verified exhaustively n<=6, sampled n=7,8.
**What:** F_k(T) = A(n,k) = C(n-1, k) mod 2 for ALL tournaments T. F(T,x) = (1+x)^{n-1} mod 2 is COMPLETELY tournament-independent. Proof via universal Taylor zeros mod 2 (c_j = 0 for j < n-1) + Redei's theorem (F_{n-1} = Hamiltonian path count is always odd). The mod-2 result is the strongest possible: individual F_k are determined, not just linear combinations.
**Key insight:** p=2 is special because (1) val_2(A_n(x)) = n-1 (maximal), giving a single free parameter, and (2) Redei pins that parameter to 1.
**Mod-p generalization (S38):** For p >= 5, universal Taylor zeros match Eulerian valuation only for n >= p+2. The Eulerian conjecture (p|A(n,k) => p|F_k(T)) FAILS for p=5 at n=7 because multiple free parameters in F(T,x) mod 5 allow different zero patterns.
**Scripts:** `04-computation/fk_mod2_proof.py`, `taylor_zeros_mod_p.py`, `mod_p_general_conjecture.py`
**Next step:** (1) Prove universal Taylor zeros mod 2 algebraically (c_j = 0 for j < n-1). (2) Is there an elementary proof not using THM-086 machinery?

### INV-032: Omega(T) structural properties Ã¢â‚¬â€ PARTIALLY DISPROVED
**Source:** Web research opus-S5, opus-S7 (disproof), opus-S9 (line graph disproof), opus-S10 (structure analysis)
**Status:** DISPROVED: Omega(T) is NOT always claw-free (fails n=9, 90%) or perfect (fails n=8, 53.8%). NOT a line graph (K_5-e found at n=6, 45%). S_{1,1,1}-free through n=11, fails n=12.
**DISPROVED (THM-025, opus-S18):** Real-rootedness of I(Omega(T),x) FAILS at n=9. Counterexample: score [1,1,3,4,4,4,6,6,7], I=[1,94,10,1], Newton k=2 fails (100 < 141), complex roots confirmed.
**Structural characterization (opus-S19):** Failure requires (a) three vertex-disjoint 3-cycles partitioning V, AND (b) near-total inter-group domination (9-0, 9-0, 7-2 arc counts), creating a transitivity bottleneck with hub vertex in 92/94 cycles. The extreme tournament (full domination) gives I=(1+x)^3 with disc=0 exactly. One arc flip creates disc<0. Failure is MAXIMALLY RARE: 0 in 10000 random samples at n=9, 0 at n=10,11.
**What remains true:** Real-rootedness holds at n<=8 (Chudnovsky-Seymour, claw-free) and for "generic" tournaments at all n. But it is NOT a universal structural property.
**Next step:** (1) Characterize exact tournament class where failure occurs. (2) Check if the clique-deletion interlacing approach (INV-038) can prove real-rootedness under a claw-free assumption.
**What remained true:** All-real-roots of I(Omega(T), x) appears to hold even for imperfect/non-claw-free Omega (tested n<=10, 0 failures). This is a deep structural conjecture NOT explained by any forbidden subgraph property.
**Note:** OCF is now proved by Grinberg-Stanley, so this is no longer a proof strategy Ã¢â‚¬â€ it's a structural question. Real-rootedness explanation must be algebraic (Irving-Omar/Grinberg-Stanley symmetric function framework).
**Extended testing (opus-S18):** Real-rootedness tested for I(Omega_3(T), x) at n=9-21 with 0 failures across 1470+ samples (degrees up to 5). Log-concavity and Newton's inequalities hold in all cases. The "Omega_3 complement = matching" structure holds exhaustively at nÃ¢â€°Â¤6 (31088/31088) but fails at nÃ¢â€°Â¥7 (75.3%).
**TurÃƒÂ¡n-based proof for nÃ¢â€°Â¤11:** At n=9-11, alpha(Omega_3) = 3, so the disjoint-pair graph is triangle-free. TurÃƒÂ¡n gives a2 Ã¢â€°Â¤ c3Ã‚Â²/4, proving Newton's first inequality a1Ã‚Â² Ã¢â€°Â¥ 3a2. Combined with the degree-3 discriminant bound, this could give a complete proof at nÃ¢â€°Â¤11. For nÃ¢â€°Â¥12, TurÃƒÂ¡n alone fails.
**Next step:** (1) Complete TurÃƒÂ¡n+discriminant proof for n=9-11. (2) Find tournament-specific bounds on a2 for nÃ¢â€°Â¥12. (3) Investigate Irving-Omar determinantal formula for algebraic proof.

### INV-038: Clique-deletion interlacing for Omega(T)
**Source:** opus-2026-03-06-S17, T100, interlacing-clique-deletion.md
**Status:** STRUCTURAL INSIGHT. Proof sketch for n<=8.
**What:** Through-v cycles always form a CLIQUE in Omega(T) (proved: sharing vertex v). Deleting vertex v = deleting this clique from Omega(T). The sequential deletion recurrence I(G,x) = I(G-u,x) + x*I(G\N[u],x) can be applied step-by-step. At n=5, 100% of remaining cycles are adjacent to some through-v cycle (Omega is very dense).
**Key insight:** For n<=8 (claw-free), Chudnovsky-Seymour guarantees each step preserves real roots. For n>=9, real-rootedness can FAIL (THM-025), so the interlacing approach cannot extend to all tournaments. However, 84 claws at n=9 counterexample all share the same 3 leaves Ã¢â‚¬â€ the claw structure is very specific.
**Verification:** 0 failures: n=5 (5120 exhaustive), n=6 (196608 exhaustive), n=7-8 (random).
**Impact:** Proves real-rootedness for n<=8. For n>=9, would need a tournament-specific claw-free condition.
**Next step:** (1) Characterize when Omega(T) is claw-free at n>=9. (2) Check if "generically claw-free" suffices for applications.
**Scripts:** `04-computation/interlacing_verify.py`, `04-computation/interlacing_structure.py`
**Writeup:** `03-artifacts/drafts/interlacing-clique-deletion.md`

### INV-039: Blueself odd-n obstruction Ã¢â‚¬â€ PROVED for ALL odd n
**Source:** opus-2026-03-06-S17, THM-022 Theorem 5 (upgraded)
**Status:** PROVED. Pure algebraic proof, no exhaustive search needed.
**What:** No blueself tilings exist at any odd n. Grid-symmetry forces k_0+k_{n-1}=n-2 (endpoint constraint). Flip changes endpoint multisets: {1+k_0, n-2-k_0} -> {n-1-k_0, k_0}. For these to be equal as multisets, need k_0=(n-2)/2 (non-integer at odd n) or 1=0 (impossible). Therefore sorted scores always differ, so flip(T) is never isomorphic to T.
**Script:** `04-computation/blueself_odd_n_proof.py`
**Impact:** Upgrades THM-022 Theorem 5 from "proved n<=7" to "proved all n". Completes the odd half of the blueself existence dichotomy.

### INV-040: Blueself vs SC maximizer Ã¢â‚¬â€ DISPROVED
**Source:** opus-2026-03-06-S17, T099
**Status:** DISPROVED at n=6.
**What:** At n=6, blueself class with H=41 is NOT the SC maximizer in score class (3,3,3,2,2,2) (SC max is H=45, also blueself). Blueself classes are always SC and have regular scores, but not always max-H. The blueself with higher disjoint pair count (alpha_2=4) beats the one with more total cycles (alpha_1=16).
**Script:** `04-computation/blueself_sc_maximizer_connection.py`

### INV-041: Quasi-regularity of Omega(T) Ã¢â‚¬â€ EXPLAINED
**Source:** opus-2026-03-06-S17 (T101), opus-2026-03-06-S18 (T103, proof)
**Status:** EXPLAINED. Theoretical argument + verified n=5-20.
**What:** Omega_3(T) is quasi-regular because adjacency depends on vertex-set intersection (sharing Ã¢â€°Â¥1 vertex), not arc orientations. This makes Omega_3 an induced subgraph of J(n,3) (Johnson graph), inheriting its regularity. All 3-element subsets have identical intersection statistics, so degree of each 3-cycle concentrates around E[deg] = (C(n,3)Ã¢Ë†â€™C(nÃ¢Ë†â€™3,3))/4Ã¢Ë†â€™1. The coefficient of variation CV = O(1/Ã¢Ë†Å¡m) Ã¢â€ â€™ 0 as nÃ¢â€ â€™Ã¢Ë†Å¾, giving ÃŽÂ»_max/avg_deg Ã¢â€°Ë† 1+CVÃ‚Â² Ã¢â€ â€™ 1. Verified: CV drops from 0.05 (n=6) to 0.03 (n=20). This does not directly explain real-rootedness but constrains spectral structure.
**Scripts:** `04-computation/omega_spectral_fast.py`, `04-computation/omega_quasireg_proof.py`

### INV-042: Paley deletion maximizer Ã¢â‚¬â€ VERIFIED p=3,7,11
**Source:** kind-pasteur-2026-03-06-S18e (T097), opus-2026-03-06-S18
**Status:** VERIFIED at p=3,7,11. Conjecture for all Paley primes.
**What:** H(T_p Ã¢Ë†â€™ v) = a(pÃ¢Ë†â€™1) (OEIS A038375 max H at n=pÃ¢Ë†â€™1). Verified: T_3Ã¢Ë†â€™v: H=1=a(2), T_7Ã¢Ë†â€™v: H=45=a(6), T_11Ã¢Ë†â€™v: H=15745=a(10). By vertex-transitivity all deletions equivalent. Combined with T053 (T_p achieves a(p)), the maximizer chain is "hereditary" via vertex deletion. Claim A decomposition for T_7: diff=144=2Ãƒâ€”72, sum_mu=6Ãƒâ€”3+30+24=72 (all 3-cycle complements have a 3-cycle in Paley).
**Next step:** Verify at p=19 (need H(T_19Ã¢Ë†â€™v) = a(18)). Investigate: does the n=pÃ¢Ë†â€™1 maximizer always come from Paley deletion, or can it be achieved by non-Paley tournaments too?
**Scripts:** `04-computation/paley_deletion_test.py`

### INV-043: Anti-aut involution existence Ã¢â‚¬â€ PROVED (THM-024)
**Source:** opus-2026-03-06-S18 (T102, THM-024), correcting kind-pasteur S18e (T095)
**Status:** PROVED. Clean group theory argument.
**What:** Every SC tournament has Ã¢â€°Â¥1 involution anti-automorphism. Proof: (1) Moon's theorem: |Aut(T)| is odd. (2) H = Ã¢Å¸Â¨Aut(T), ÃÆ’Ã¢â€šâ‚¬Ã¢Å¸Â© has order 2|Aut(T)| (even). (3) By Cauchy, H has order-2 element. (4) Can't be in Aut(T) (odd order group). (5) Must be in ÃÆ’Ã¢â€šâ‚¬Ã‚Â·Aut(T) = set of anti-auts. NOT all anti-auts are involutions (counterexamples at n=6 with |Aut|>1), but at least one always is.
**Scripts:** `04-computation/anti_aut_involution_test.py`, THM-024

### INV-044: Hereditary Maximizer Chain Ã¢â‚¬â€ CORRECTED (regular-only at odd n)
**Source:** kind-pasteur-2026-03-06-S18f, S18g (correction), T104, T105, MISTAKE-010
**Status:** CORRECTED. Only REGULAR maximizers at odd n are hereditary.
**What:** Previous claim "all maximizers at odd n hereditary" was WRONG. Exhaustive check:
- n=3: 2/2 hereditary (all regular)
- n=4: 0/24 hereditary (non-regular)
- n=5: 24/64 hereditary (only 24 regular, NOT 40 with score (1,2,2,2,3))
- n=6: 0/480 hereditary (non-regular)
- n=7: 240/240 hereditary (all regular)

**R-minimization conjecture (NEW, S18g, REFUTED at n=7 Ã¢â‚¬â€ see OPEN-Q-017):** The H-maximizer was conjectured to minimize R(T) = sum_v H(T-v)/H(T). Formula R(T) = n - E_weighted[|U(S)|] is PROVED. But R-minimization FAILS at n=7: a tournament with H=123 has R=1.585 < R(max)=5/3. Valid only at n=3-6.

**Key insight:** Being hereditary (R = n*H_{n-1}/H_n for regular maximizers) is NOT the same as minimizing R. The non-regular n=5 maximizers have LOWER R (1.4) than regular ones (5/3 Ã¢â€°Ë† 1.667) despite not being hereditary.
**Next step:** (1) Verify R-minimization at n=7 (running). (2) Prove R-minimization from OCF. (3) Test if regular n=9 maximizers are hereditary.
**Scripts:** `04-computation/hereditary_maximizer.py`, `04-computation/hereditary_correction.py`, `04-computation/R_minimization_proof.py`, `04-computation/R_min_n7_check.py`
### INV-032: Omega(T) is always claw-free AND perfect Ã¢â‚¬â€ Dyer-Jerrum decomposition
**Source:** Web research opus-S5, arXiv:1909.03414 (Dyer-Jerrum-MÃƒÂ¼ller-VuÃ…Â¡koviÃ„â€¡)
**Status:** CLAW-FREENESS PROVED for n<=8 (THM-020, vertex counting: claw needs 9+ vertices). PERFECTNESS FAILS at n=8 (THM-019, 53.8% of random n=8 tournaments have C5 in Omega).
**What:** Omega(T) is always claw-free for n<=8 but NOT always perfect for n>=8. The Dyer-Jerrum decomposition framework still applies to claw-free graphs (without perfectness), but the structure is less constrained.
**Critical consequence:** Chudnovsky-Seymour (2007) proved that I(G,x) has ALL REAL ROOTS for claw-free G. Since Omega(T) is claw-free for n<=8, ALL roots of I(Omega(T), x) are negative real (THM-020). At n>=9, claw-freeness fails (90% of random tournaments have a claw in Omega). Real-rootedness may still hold by other mechanisms.
**Next step:** (1) Test claw-freeness at n=9 exhaustively. (2) Test line graph hypothesis (Beineke forbidden subgraphs). (3) Test subdivided-claw-freeness at n=9.

### INV-033: Redei-Berge Hopf algebra formalization of OCF
**Source:** Web research opus-S5, arXiv:2402.07606 (Grinberg), arXiv:2506.08841 (Mitrovic-Stojadinovic)
**Status:** CONNECTION IDENTIFIED. Key bridge found (S36).
**What:** The Redei-Berge symmetric function U_X for digraphs has comultiplication Delta([X]) = sum_S [X|S] tensor [X|V\S] Ã¢â‚¬â€ this IS our subset convolution. The character zeta counts Hamiltonian paths. The antipode S(U_X) = (-1)^|V| U(X-bar) encodes Berge's theorem.
**NEW (Mitrovic-Stojadinovic, arXiv:2506.08841, June 2025):**
  - X_{inc(P)} = omega(U_P): Chromatic function of incomparability graph = omega of Redei-Berge
  - "Converse of Redei": if poset is not a chain, quasi-linear extensions are even
  - Bags-of-sticks decomposition: U_X = sum of simpler digraphs via inclusion-exclusion on edges
  - Stanley-Stembridge connection: e-positivity of X_{inc(P)} <=> h-positivity of U_P
  - Noncommutative deletion-contraction: W_X = W_{X\e} - W_{X/e}^up
  - Mitrovic-Stojadinovic phi(pi) = sum_{gamma X-cycle} (len(gamma)-1) is EXACTLY our S = sum(l_i-1)!
**Verified (S36):** OCF specialization p_1->1, p_{odd>=3}->2, p_{even}->0 gives H(T) from U_T.
**NEW (Mitrovic, arXiv:2504.20968, April 2025):** Noncommutative Redei-Berge function W_X has deletion-contraction: W_X = W_{X\e} - W_{X/e}Ã¢â€ â€˜. Thm 3.16: cycle decomposition via inclusion-exclusion over cycle edges. Cor 3.12: tournament formula W_X = ÃŽÂ£(2^{ÃË†(ÃÆ’)} p_{Type(ÃÆ’)}) for odd-cycle permutations = exactly OCF.
**h-POSITIVITY TEST (kind-pasteur-S39b):** h-positivity of U_T FAILS for all non-transitive tournaments. At n=3: 1/2 h-positive, n=4: 1/4, n=5: 1/11. Only the transitive tournament (H=1) is h-positive. The h(2,1) and h(2,2,1) coefficients are always negative for non-transitive. This is expected since tournament posets are NOT (3+1)-free in general, and Stanley-Stembridge conjecture requires (3+1)-freeness.
**Next step:** (1) Express OCF via bags-of-sticks decomposition. (2) Check if deletion-contraction on W_T gives a direct proof of Claim A. (3) Explore chromatic function connection for imperfect Omega(T). (4) Study Thm 3.16 cycle decomposition for odd cycles.

### INV-034: BjÃƒÂ¶rklund cycle cover reduction adapted for OCF Ã¢â‚¬â€ TESTED (NEGATIVE for new identities)
**Source:** Web research opus-S5, arXiv:1008.0541, arXiv:1301.7250
**Status:** TESTED (kind-pasteur-S39b). No new identity beyond OCF.
**What:** BjÃƒÂ¶rklund reduces Hamiltonian cycle counting to cycle cover counting via inclusion-exclusion and determinants. Tested 6 formulations at n=3-6:
1. Full-vertex all-odd CC weighted: FAILS (0%)
2. Partial odd CC weighted by 2^k: MATCHES H(T) 100% Ã¢â‚¬â€ but this IS OCF
3. Inclusion-exclusion sum (-1)^{n-|S|} perm(T[S]): FAILS (for paths)
4. Irving-Omar odd traces: exploratory only
5. perm(I + x*A): FAILS
6. Odd permanent (unweighted): FAILS (0%)
**Conclusion:** OCF = partial odd cycle cover polynomial at weight 2. This is a restatement, not a new identity. The BjÃƒÂ¶rklund approach doesn't give a new route to proving OCF.
**Script:** `04-computation/bjorklund_cycle_cover.py`

### INV-035: Tribonacci structure Ã¢â‚¬â€ OCF for T_full family via interval graphs
**Source:** opus-2026-03-05-S6 (Tribonacci web research), kind-pasteur-S11 (Tribonacci discovery)
**Status:** VERIFIED n=3,...,8. Both sides match Tribonacci(n) = A000213 independently.
**What:** T_full_n (full tiling tournament) has H(T_full_n) = Tribonacci(n) (proved via run decompositions). INDEPENDENTLY, Omega(T_full_n) is an INTERVAL GRAPH on odd-length consecutive intervals [k, k+2j], and I(Omega, 2) satisfies the same Tribonacci recurrence via a weighted interval packing DP that telescopes: f(n) = f(n-1) + 2f(n-3) + 2f(n-5) + ... = f(n-1) + f(n-2) + f(n-3).
**Key structural insight:** All directed odd cycles of T_full_n are consecutive intervals. The clique-cutset decomposition of this interval graph mirrors the DP structure computing H(T_full_n). Both sides produce Tribonacci by the same algebraic mechanism (telescoping) through different combinatorial objects.
**Why this matters:** Shows OCF's "both sides match" emerges from parallel decomposition structures. If this parallelism generalizes (clique-cutset of Omega mirrors Ham path DP), it could prove OCF.
**Extended results (opus-S13):**
- **Transitive+flip(i,j):** H = 1 + 2^(j-i-1). All odd cycles form a clique in Omega, so I(Omega,2) = 1 + 2Ã‚Â·(#cycles) = 1 + 2^(j-i-1). Clean OCF-based proof.
- **Cone theorem:** H(source_cone(T')) = H(sink_cone(T')) = H(T') for ALL T'. Proved: source must be first in every Ham path. Verified exhaustively through n'=6.
- **Partial cones palindromic:** H(k) = H(n'-k) where k = out-degree of cone vertex. From self-converse symmetry.
- **Circulant S={1}:** H(T_{n,{1}}) = n, order-2 recurrence. No circulant with |S|>=2 has low-order recurrence.
- **CORRECTED (S41):** Circulant S={1,5,6,7} at n=9 DOES give H=3357 = max. Non-circulant maximizers also exist.
**Next step:** (1) Find direct bijection between run decompositions and weighted interval packings. (2) Check if the transfer matrix for T_full has Tribonacci characteristic polynomial factor.

### INV-001: Prove transfer matrix symmetry for all n Ã¢â‚¬â€ PROVED (THM-030)
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
- **Off-diagonal sum:** sum_{aÃ¢â€°Â b} M[a,b] = 0 (odd n), 2*H(T) (even n). Verified n=3,...,7 but NOT yet proved.
- **Complement pairing D(S)+D(U\S) is constant at n=4 but NOT at n=5**, ruling out the simplest telescoping argument.
- **Cauchy-Binet decomposition:** M = E^T * Lambda * B where E[S,v]=E_v(S), B[S,v]=B_v(U\S), Lambda=diag((-1)^|S|). Symmetry equivalent to E^T*Lambda*B = B^T*Lambda*E.
**PATH REVERSAL PROOF AT c=0 (kind-pasteur-S23):** COMPLETE proof when c=0 (pure skew weights). Path reversal: B_v(S+v) = (-1)^|S| E_v(S+v). This gives M[a,b] = (-1)^{n-2} sum_S E_a(S+a) E_b(R+b) Ã¢â‚¬â€ unsigned, manifestly symmetric by S<->R relabeling. Verified n=3,4,5,6.
**EVEN r-POWERS CONJECTURE (kind-pasteur-S23):** At general c, M(r,s) where r=c/2 has ONLY even r-powers. Equivalent to symmetry. Verified n=3,4,5,6. Path reversal gives B_v(c,s) = E_v(c,-s), which yields M[b,a](c,s) = (-1)^{n-2} M[a,b](c,-s). So symmetry reduces to M having definite s-parity (-1)^{n-2}, i.e., only even r-powers.
**ALGEBRAIC PROOF (kind-pasteur-S23b): M[a,b](-r) = M[b,a](r)** Ã¢â‚¬â€ 5-step proof: T(-r)=-T^T, path reversal under negated transpose, sign bookkeeping, SÃ¢â€ â€R relabeling. Verified n=4,5,6. This proves the EQUIVALENCE between (i) symmetry, (ii) even-r-powers, (iii) s-parity.
**TOGGLE ANALYSIS (S23b):** At n=4, r^1 monomials cancel pairwise between different S-subsets. At n>=5, cancellation is multi-way (not simple pairwise). No clean single-vertex toggle involution found on whole subsets.
**H(U) MATRIX (S23b, from Kogan/Hamiltonian cycle polynomial):** H(U)_{i,j} = sum of Ham path weights from i to j. Identity: U*H(U)^T = H(U)*U^T. For c-tournaments U+U^T = c(J-I), this gives UH^T = H(cJ-cI-U), but does NOT directly imply H=H^T. Also note: M[a,b] is NOT the same as H(T)_{a,b} (M has inclusion-exclusion signs, H is a direct sum).
**PROVED INDEPENDENTLY by both opus-S25 and kind-pasteur-S25 (THM-030).** Key Identity: odd_r(B_b(W)) = r * col_sum_W(b), equivalently B_b(W)+(-1)^m E_b(W) = 2r*col_sum(b). Inductive proof using column recurrence + first-edge decomposition. The Sigma identity (r*Sigma = odd(T)) follows from summing the inductive hypothesis. The proof closes because odd(sum s*even(B_v)) = 0. Verified computationally m=2..6 all (a,b) pairs. See complete_even_r_proof.py, key_identity_complete_proof.py.
**Scripts:** `04-computation/symbolic_symmetry_proof.py`, `04-computation/transfer_symmetry_analysis.py`
**COEFFICIENT STRUCTURE (opus-S22 continuation):**
- **[r^{n-2}] = (n-2)!** when n even, **0** when n odd. Proof: counting argument, sum_k C(n-2,k)(-1)^k k!(n-2-k)! = (n-2)! * sum_k(-1)^k. Verified n=3,...,6.
- **[r^2] for n=5 = 2Ã‚Â·sum_{uÃ¢Ë†Ë†U}(s_{au}+s_{bu})**. For n=6: degree-2 in s with all coefficients Ã‚Â±2. For n=4: just 2 (constant).
- **[r^1] telescoping (n=4):** Each s_{uv} (uÃ¢Ë†Ë†U, vÃ¢Ë†Ë†{a,b}) appears exactly once with + and once with - across subsets. Moving vertex u between S and R flips the sign contribution.
- **M is NOT a cofactor** of A=rJ'+S (exhaustive test n=3,4). Cofactors have degree n-1; M has degree n-2.
- **M is NOT a permanent minor** of A either. The fundamental identity A(-r)=-A^T is clean but M does not decompose as a simple matrix function of A.
- **Key algebraic identity:** A(-r,s) = -A(r,s)^T (since J' symmetric, S skew). Any expression built from AÃ‚Â·A^T+A^TÃ‚Â·A (which is even in r) could explain the property, but no such expression matching M has been found.
- **Literature update:** Irving-Omar (arXiv:2412.10572), Mitrovic noncommuting Redei-Berge (arXiv:2504.20968) with deletion-contraction W_X = W_{X\e} - W_{X/e}^up, El Sahili-Ghazo Hanna proving T and T^op have same Hamiltonian path type distribution.
**APPROACH RULING (opus-S22 continuation):**
- Ã¢ÂÅ’ Simple cofactor/minor of A (degree mismatch)
- Ã¢ÂÅ’ Permanent of A minor (doesn't match)
- Ã¢ÂÅ’ Adjugate entries of A, IÃ‚Â±A, J-A (all fail)
- Ã¢Ââ€œ Deletion-contraction via Mitrovic noncommuting Redei-Berge (unexplored, most promising NEW lead)
- Ã¢Ââ€œ Irving-Omar walk generating function det(I+zXÃ„â‚¬)/det(I-zXA) (connection to M unclear)
- Ã¢Ââ€œ Direct r^1=0 proof via telescoping + induction on n (promising for base case)
**Next step:** (1) Try Mitrovic deletion-contraction approach Ã¢â‚¬â€ express M[a,b] recursively and prove even-r by induction. (2) Understand Irving-Omar matrix formula and whether it encodes M[a,b]. (3) Prove [r^1]=0 directly via the telescoping structure observed at n=4,5. (4) Previous approaches (Feng, Hopf, involution) remain viable but untested.
**Scripts:** `04-computation/symbolic_symmetry_proof.py`, `04-computation/transfer_symmetry_analysis.py`, `04-computation/determinantal_identity_test.py`, `04-computation/det_compare_explicit.py`, `04-computation/r1_coefficient_analysis.py`, `04-computation/r_coefficient_structure.py`

### INV-002: Subset convolution identity Ã¢â‚¬â€ the core algebraic challenge
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

### INV-053: Even Cycle Vanishing Theorem Ã¢â‚¬â€ PROVED
**Source:** opus-2026-03-06-S10, T148
**Status:** PROVED. Clean involution argument.
**What:** For any tournament T on [n], p_mu(U_T) = 0 whenever mu has an even part. The proof pairs each permutation sigma with even k-cycle c with sigma' (c reversed); the sign flips because (-1)^{k-1} = -1 for even k. Verified computationally n=3 through n=7.
**Consequences:** (1) U_T lives in the subspace spanned by p_mu with all odd parts Ã¢â‚¬â€ drastically fewer terms. (2) At n=4, only types (1^4) and (3,1) contribute; at n=5, only (1^5), (3,1,1), (5); at n=6, only (1^6), (3,1,1,1), (3,3), (5,1). (3) The Schur expansion simplifies: [s_lambda]U_T = sum over odd-part mu only. (4) This is the SAME T<->T^op involution as in the path reversal proof (T147).
**Connection to INV-001:** The even-r-powers conjecture (kind-pasteur-S23) is the transfer matrix version of this same phenomenon. Both arise from the perpendicular grid symmetry.

### INV-054: Hook Schur Positivity for Tournaments Ã¢â‚¬â€ PARTIAL (fails at n=7)
**Source:** opus-2026-03-06-S10, T149
**Status:** PROVED at n=4 (clean sign argument). VERIFIED at n=5 (11/11), n=6 (40/40). **FAILS at n=7** (231/242, 11 failures all for middle hook (4,1,1,1)).
**What:** [s_{(k,1^{n-k})}]U_T >= 0 holds at n=4,5,6 but fails at n=7.

**n=4 proof:** Only p-types (1^4) and (3,1) matter (by INV-053). All hook characters non-negative at both Ã¢â€ â€™ sum of non-negative terms.

**n=7 failure mechanism:** chi^{(4,1,1,1)}((7)) = -1 and chi^{(4,1,1,1)}((5,1,1)) = 0. Regular tournament T_7 has 48 directed 7-cycles contributing -48/7, overwhelming positive 3-cycle terms. Result: [s_{(4,1,1,1)}]U_{T_7} = -83/28 Ã¢â€°Ë† -2.96.

**Which hooks always hold:** Hooks (n) and (1^n) have all-positive characters at odd types Ã¢â€ â€™ always positive. Hooks (n-2,1,1) and (3,1^{n-3}) also have all-positive chars (at n=7). Only hooks with j odd and j near n/2 can fail.

**Non-hook negativity:** Non-hook chars at (3,1,...,1) are always negative Ã¢â€ â€™ non-hook coefficients negative for all non-transitive tournaments (verified n=4,5,6).

**Refined question:** For which hooks is positivity universal? Is there a simple characterization (e.g., j even, or |j - n/2| > threshold)?
**Scripts:** `04-computation/schur_hook_analysis.py`, `04-computation/tournament_cycle_structure.py`, `04-computation/even_cycle_vanishing_proof.py`, `04-computation/hook_positivity_n6.py`, `04-computation/hook_positivity_n7.py`

---

### INV-200: Theorem number collisions Ã¢â‚¬â€ 10 numbers have multiple files (HOUSEKEEPING)
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
**Finding:** Chapman's bijection maps oriented monotone triangles (not ASMs directly) to tournaments. Striker's tournament bijection uses the disjoint three-color poset {b,r,(g)}, while ASMs use the four-color poset {b,y,o,g} Ã¢â‚¬â€ different subposets. No S3 action is defined in either paper. The question needs precise formulation: (a) define "the bijection" (Chapman's Phi? Striker's poset? composition?), (b) define how S3 acts on both sides. Only then can it be tested at n=3,4.
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

### INV-011: Mod-4 score-sequence criterion (Open Problem 9 in paper) Ã¢â‚¬â€ RESOLVED (NO)
**Source:** oq:mod4_struct in tex
**Status:** RESOLVED NEGATIVELY (kind-pasteur-2026-03-07-S39b).
**What:** Does score sequence determine H(T) mod 4?
**Answer:** NO for n >= 5. The formula H mod 4 = (1 + 2Ã‚Â·c3) mod 4 holds only for n <= 4.
**Results:**
- n=3,4: YES (exhaustive, H mod 4 constant within every score class)
- n=5: NO Ã¢â‚¬â€ score (1,2,2,2,3) has H mod 4 in {1,3} (c5 varies within this class)
- n=6: NO Ã¢â‚¬â€ 5/22 score classes have varying H mod 4
- n=7: NO Ã¢â‚¬â€ 27/59 sampled score classes have varying H mod 4
**Key insight:** H mod 4 = (1 + 2Ã‚Â·alpha_0) mod 4 where alpha_0 = total odd cycle count.
For n <= 4, alpha_0 = c3 which is score-determined. For n >= 5, 5-cycles contribute to alpha_0 but c5 is NOT determined by score sequence (confirmed independently).
**Also found:** c5 is NOT determined by (score, sum_dÃ‚Â², edge_score_sum). Even (score, common_out_neighbor) pairs vary. c5 requires genuine graph structure beyond all local/pairwise statistics.
**Scripts:** `04-computation/mod4_score_test.py`, `04-computation/c5_score_determination.py`

### INV-012: BlackSelf(8) exceptional class (Open Problem 7 in paper) Ã¢â‚¬â€ PARTIALLY RESOLVED
**Source:** oq:n8 in tex
**Status:** DEEP INVESTIGATION by opus-2026-03-05-S8. See `03-artifacts/drafts/n8-anomaly-deep-dive.md`.
**What:** Unique isomorphism class at n=8 that is self-converse, has |Aut|>1, |Fix(beta)| odd, but H(T)/|Fix(beta)| is even. Is it related to a Hadamard matrix of order 8 or a skew conference matrix?
**Key findings:**
- Exhaustive search over ALL 65536 SC tournaments (alpha = reversal): ZERO have Fix(beta) | H with H/Fix even.
- The definition likely means (H-Fix)/2 is even (number of beta-orbit pairs is even).
- Under this interpretation: T_657 (H=657, Fix=33, |Aut|=3) is the best candidate: (657-33)/2 = 312 (even).
- T_657 CONTAINS P(7) (Paley tournament on 7 vertices) as vertex-deletion. P(7) Ã¢â€ â€ unique skew Hadamard matrix of order 8. THIS IS THE HADAMARD CONNECTION.
- T_657 has perfectly uniform D_v = 54 (mu-weighted 3-cycle count) for all 8 vertices.
- Full survey: 10 distinct (H, Fix, Aut) combinations among 2560 SC+Aut>1 tournaments.
**Next step:** Confirm T_657 is isomorphic to a known Paley extension. Resolve definition ambiguity with paper author.

### INV-013: Realizable odd-cycle conflict graphs (Open Problem 8 in paper) Ã¢â‚¬â€ INVESTIGATED
**Source:** oq:realizable in tex
**Status:** INVESTIGATED (opus-S13 background agent). Key structural findings.
**What:** Which graphs G arise as Omega_3(T) for some tournament T?
**Results:**
- Realizable isomorphism classes: 2 (n=3), 3 (n=4), 6 (n=5), 18 (n=6), ~97+ (n=7 sampled)
- nÃ¢â€°Â¤5: Omega_3 is ALWAYS a complete graph (two 3-cycles on Ã¢â€°Â¤5 vertices must share a vertex)
- n=6: Always "complete minus matching" Ã¢â‚¬â€ complement is disjoint union of edges
- n=7: First non-perfect graphs (11/97 classes have chi=omega+1)
- alpha(Omega_3) Ã¢â€°Â¤ floor(n/3) Ã¢â‚¬â€ proved by vertex counting (3k vertices for k disjoint 3-cycles)
- 100% real-rooted across all realizable classes (exhaustive nÃ¢â€°Â¤6, sampled n=7)
**Key insight:** The low independence number (alpha Ã¢â€°Â¤ 2 for nÃ¢â€°Â¤8) means I(Omega_3, x) has degree Ã¢â€°Â¤ 2, so real-rootedness is "easy" for Omega_3. The full Omega (including 5,7-cycles) has alpha Ã¢â€°Â¤ floor(n/3) similarly, keeping degree low.
**Next step:** Characterize which "complete minus matching" graphs at n=6 are NOT realizable. Extend to n=8 where alpha=2 still but degree may increase.

### INV-014: 2-adic tower / higher Redei theorems Ã¢â‚¬â€ PARTIALLY RESOLVED
**Source:** OPEN-Q-008, T007, tex Section 5.5
**Status:** COMPUTED (opus-S13). v_2(H(T)) = 0 ALWAYS (= Redei's theorem).
**What:** I(Omega(T), x) at x=4,8,... gives mod-4, mod-8 invariants of H(T). v_2(H(T)) = 0 universally.
**Results:** H mod 4 Ã¢â€°Â¡ 1+2*alpha_1 (mod 4) via OCF. At n=3,4 this equals 1+2*c3 (mod 4) exactly. At nÃ¢â€°Â¥5 the c3 formula breaks (5-cycles contribute to alpha_1). H mod 2^k approaches uniform on odd residues as n grows.
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
7. Beta_5 NOT YET observed at n=7,8,9 Ã¢â‚¬â€ onset unknown
8. Chi(T) in {0,1} for n<=7; chi up to 6 at n=8 (HYP-312)
9. dim(Omega_2) NOT determined by (c3, score) (HYP-308/309)
10. |A_p| mod 2 = C(n,p+1) mod 2 via local Redei
11. Poincare polynomial P(T,1)/(2^n-1) grows with n Ã¢â‚¬â€ path complex exceeds simplex
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
9. THM-119 (was THM-097) (PROVED): Disjoint support at Omega_2 Ã¢â‚¬â€ each 2-path has at most 1 non-allowed face. Constraint matrix always full rank. dim(Omega_2) = |A_2| - #NA_faces exactly.
10. Completeness is SHARP: removing 1 edge from tournament creates beta_2>0 (13/500 at n=6).
11. beta_2=0 confirmed at n=9 (0/500) and n=10 (0/100). Disjoint support verified at n=9.
12. H_1(A-complex) is NOT always 0 (only for transitive tournaments). Omega restriction is essential.
13. rank(d_3) is NOT a simple function of c3 (multiple values per c3 class).
14. 3-path NA face distribution: exactly 25%/50%/25% for 0/1/2 faces (very clean, universal).
15. At level 4: paths can have multiple NA faces => overlapping constraint rows => rank deficit => beta_4 possible.
16. Complete bidirectional graph: beta_2=0, homology at top dimension only.
**S45 defect rate analysis (omega2_exact_formula.py, defect_rate_ushape.py, ecyc_formula.py):**
17. dim(Omega_2) = C(n,3) + 2*c3 - e_cyc EXACTLY (exhaustive n=4,5,6). e_cyc = #{directed edges in Ã¢â€°Â¥1 three-cycle}.
18. e_cyc NOT determined by c3 alone Ã¢â‚¬â€ depends on cycle arrangement (edge sharing). Constant for most score seqs, varies near-regular.
19. Defect rate is WAVE PROPAGATION, not U-shape: beta_1 rate decreasing (29.7%Ã¢â€ â€™1%), beta_3 increasing (0%Ã¢â€ â€™21%), beta_4 appears at n=8 (2%).
20. beta_4 CAN be nonzero (values 1,2 at n=8) Ã¢â‚¬â€ "all even beta vanish" is FALSE. Only beta_2=0 always.
21. Only 3-5 distinct Betti profiles per n. Extremely constrained.
22. beta_3*beta_5=0 at n<=8 (trivially since beta_5 rarely nonzero). Generalized seesaw needs testing at n>=10.
**Open:** (1) Beta_2=0 is PROVED (THM-108+109). (2) Does generalized seesaw beta_{2k-1}*beta_{2k+1}=0 hold generally? (3) Onset of beta_6 (not seen at n<=8 in 600+ samples). (4) Why exactly 3-5 Betti profiles per n?
**Scripts:** beta1_beta3_mediator.py, even_betti_quick_v2.py, beta4_investigation.py, beta_coexistence_analysis.py, beta2_algebraic_analysis.py, beta2_disjoint_support_proof.py, beta2_completeness_argument.py, beta2_exactness_proof.py, beta2_only_n9.py, omega2_exact_formula.py, defect_rate_ushape.py, ecyc_formula.py

### INV-138: Beta_3 Ã¢â€°Â¤ 1 Proof Architecture (LES Induction)
**Source:** kind-pasteur-2026-03-09-S46
**Status:** PROOF ARCHITECTURE COMPLETE. Both algebraic ingredients computationally verified. No algebraic proof yet.
**What:** THM-123 (was THM-110): beta_3(T) Ã¢â€°Â¤ 1 for all tournaments T. Equivalent to rank near-saturation: rank(d_4) Ã¢â€°Â¥ ker(d_3) - 1.
**Proof strategy:** LES induction on n using pair (T, T\v):
  ... Ã¢â€ â€™ H_3(T\v) Ã¢â€ â€™ H_3(T) Ã¢â€ â€™ H_3(T,T\v) Ã¢â€ â€™ H_2(T\v) = 0
  Since H_2(T\v) = 0 (THM-108), map H_3(T) Ã¢â€ â€™ H_3(T,T\v) is surjective.
  Find v with beta_3(T\v) = 0. Then beta_3(T) = dim H_3(T,T\v).
**Key ingredients verified:**
1. Good vertex existence for beta_3: Ã¢Ë†Æ’v with beta_3(T\v)=0 when beta_3(T)>0.
   - n=6: 320/320 exhaustive. beta_3 COMPLETELY fragile (ALL 6 deletions give 0).
   - n=7: 34/34 sampled. 5-7 good vertices per tournament.
   - n=8: 31/31 sampled.
2. Relative H_3 bound: dim H_3(T,T\v) Ã¢â€°Â¤ 1 ALWAYS.
   - n=6: ALL 1920 pairs give dim=1 (exhaustive for beta_3>0).
   - n=7: dim Ã¢Ë†Ë† {0,1}. Max=1.
3. LES isomorphism: beta_3(T\v)=0 Ã¢Å¸Â¹ beta_3(T) = dim H_3(T,T\v). Perfect n=6 (1920/1920).
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
**NOTE:** HYP-342 (Boolean odd Betti) needs correction: TRUE for k=1,2 (beta_1,beta_3 Ã¢Ë†Ë† {0,1}), but FALSE for kÃ¢â€°Â¥3 (beta_5=10 at n=9 Paley maximizer). The "Boolean" property is specific to beta_1 and beta_3.

---

## Priority C: References to investigate

### INV-015: Rajkumar et al. (arXiv:2110.05188) Ã¢â‚¬â€ tournament representations
**Source:** paper-connections.md, T046, paper-deep-connections.md Section 2
**Status:** INVESTIGATED (opus-S4b, opus-S5). Key theorems extracted. Connections documented.
**What:** Flip classes, locally transitive = rank 2, R-cones, sign rank. Key results: every T in flip class of R-cone (Prop 1, distance exactly 1 cut-flip); mu(T) dimension bound <= 2(mu(T)+1) (Thm 11); sign-rank bound (Thm 12).
**Key finding:** Proposition 1 is constructive: T' = phi_{i Ã¢Ë†Âª T_i^-}(T) is R-coned by i. Cut-flip distance to R-cone is exactly 1. This directly enables INV-004 strategy.
**Assessment:** mu(T) induction (INV-005) less promising than FAS induction Ã¢â‚¬â€ mu may change unpredictably under cut-flips.
**Action needed:** ADD to bibliography. Concretely develop INV-004 (R-cone + cut-flip proof).
**Tested:** Locally transitive tournaments DO have 5/7-cycles (T046). OCF passes 100% for LT, R-cones, automorphism-symmetric tournaments.

### INV-016: Feng (arXiv:2510.25202) Ã¢â‚¬â€ dual Burnside process
**Source:** paper-connections.md, paper-deep-connections.md Section 3, tex line 1600/2151
**Status:** INVESTIGATED (opus-S4b, opus-S5). Key theorems extracted. Deep connection found.
**What:** Q=AB factorization, primal-dual spectral correspondence, lumping theory.
**Key findings:** (1) Q=AB is REVERSIBLE with pi(g)=|X_g|/(|G|*z) Ã¢â‚¬â€ detailed balance gives symmetry (Thm 3.3). (2) Block-flip M=[[0,A],[B,0]], M^2=[[Q,0],[0,K]] Ã¢â‚¬â€ bipartite structure with period 2. (3) Eigenvector intertwining (Thm 3.10): A maps K-eigenvectors to Q-eigenvectors, B maps back. (4) Our transfer matrix has EXACTLY this AB structure: A maps subsets to "path ends at vertex", B maps "path starts at vertex" to complement subsets. Transfer matrix symmetry (INV-001) = hidden detailed balance condition.
**Action needed:** Try to formalize the "hidden detailed balance" Ã¢â‚¬â€ identify the group action and show it satisfies Feng's reversibility conditions. This could prove INV-001.

### INV-017: El Sahili & Abi Aad (2020) Ã¢â‚¬â€ parity of paths in tournaments
**Source:** tex bibliography, \cite{elsahili2020}
**Status:** Referenced in tex for decisive/concordant classification. NOT deeply investigated for connections to OCF.
**What:** Discrete Math 343 (2020), Art. 111695. Mod-4 congruences.
**Action needed:** Read the paper. Check if their mod-4 results constrain or relate to our alpha_1 = |C_3| (mod 2) conjecture (INV-011).

### INV-018: El Sahili & Ghazo Hanna (2023) Ã¢â‚¬â€ number of Ham paths/cycles
**Source:** tex bibliography, \cite{elsahili2023}
**Status:** Referenced. NOT investigated for OCF connections.
**What:** J. Graph Theory 102 (2023), 684-701. About the number of oriented Hamiltonian paths and cycles in tournaments.
**Action needed:** Read the paper. They study H(T) directly Ã¢â‚¬â€ any bounds or structural results could inform OCF.

### INV-019: Schweser-Stiebitz-Toft (arXiv:2510.10659) Ã¢â‚¬â€ Redei revisited
**Source:** tex bibliography, \cite{schweser2025}, paper-deep-connections.md Section 1
**Status:** INVESTIGATED (opus-S5). Key theorems extracted.
**What:** Redei's Stronger Theorem (Thm 1.1): add non-oriented edges to tournament, #Ham paths beginning AND ending in tournament vertices is EVEN. Berge's Stronger Theorem (Thm 1.2): G and G-bar have same Ham path parity. Dirac's Stronger Theorem (Thm 2.1): inclusion-exclusion on edge subsets.
**Key findings:** (1) Direct connection to Open Problem 4 (mixed graphs). Non-oriented edges DOUBLE insertion opportunities in Q-Lemma. (2) Berge's theorem gives H(T) Ã¢â€°Â¡ H(T^op) (mod 2) for tournaments, and constrains I(Omega(G),2) under complementation. (3) Strategy for extending Q-Lemma: verify computationally that inshat remains odd for mixed graphs.
**Action needed:** Test inshat parity for mixed graphs computationally. If confirmed, Q-Lemma proof extends directly.

### INV-020: Striker (2011) Ã¢â‚¬â€ unifying poset perspective
**Source:** tex bibliography, \cite{striker2011}
**Status:** Referenced in Open Problem 5. NOT investigated.
**What:** Adv. Appl. Math. 46 (2011), 583-609. Connects ASMs, plane partitions, Catalan objects, tournaments, tableaux via posets.
**Action needed:** Read the paper. Check S3-equivariance (INV-008). The poset perspective may give new structural insights for OCF.

### INV-021: Chapman (2001) Ã¢â‚¬â€ alternating sign matrices and tournaments
**Source:** tex bibliography, \cite{chapman2001}
**Status:** Referenced in Open Problem 5. NOT investigated.
**What:** Adv. Appl. Math. 27 (2001), 290-298. Direct connection between ASMs and tournaments.
**Action needed:** Read the paper. The ASM connection could provide algebraic tools (determinantal formulas, etc.) for H(T).

### INV-022: Eplett (1979) Ã¢â‚¬â€ self-converse tournaments
**Source:** tex bibliography, \cite{eplett1979}
**Status:** Referenced briefly. NOT investigated for OCF.
**What:** Canad. Math. Bull. 22 (1979), 23-27. Self-converse tournament counts.
**Action needed:** Check if self-converse tournaments have special OCF properties. The BlackSelf(8) class (INV-012) is self-converse.

### INV-023: Forcade (1973) Ã¢â‚¬â€ parity of paths and circuits
**Source:** tex bibliography, \cite{forcade1973}
**Status:** Referenced heavily. Our paper gives a NEW combinatorial proof of his F2-invariance.
**What:** Discrete Math 6 (1973), 115-118. Original F2-invariance proof via generating functions.
**Action needed:** Compare his generating function approach to our subset convolution. His GF machinery may contain seeds of an OCF proof.

---

## Priority D: Computational targets

### INV-024: H(T_p) for Paley primes Ã¢â‚¬â€ EXTENDED
**Source:** OPEN-Q-013, opus-S9
**Status:** COMPUTED through p=23.
**Results:**
- H(P(3))=3, H(P(7))=189, H(P(11))=95095, H(P(19))=1,172,695,746,915, H(P(23))=15,760,206,976,379,349
- All match OEIS A038375 where known (a(3)=3, a(7)=189, a(11)=95095)
- P(7) confirmed as GLOBAL maximizer by exhaustive check of all 2^21 n=7 tournaments
- |Aut(T_19)| = 171, H/|Aut| = 6,857,869,865
- Ratio H(P(p))/(p!/2^{p-1}): 2.000, 2.400, 2.440, 2.527, 2.557 Ã¢â‚¬â€ converging toward e=2.718
**Sequence:** H/|Aut| = 1 (p=3), 9 (p=7), 1729 (p=11), 6857869865 (p=19).
**Next step:** Submit H(P(p)) values and H(P(p))/|Aut| sequence to OEIS. Compute H(P(31)) if feasible (2^31*31 DP Ã¢â‚¬â€ ~66B ops, might take hours).

### INV-025: Integrality conjecture C(p,k) | c_k(T_p) for k >= (p+3)/2
**Source:** T036/T153 (tangents), OPEN-Q-013 table
**Status:** VERIFIED at p=7 (kind-pasteur-S39b). Previously observed at p=11.
**What:** For Paley primes p = 3 mod 4, the cycle count c_k(T_p) is divisible by C(p,k) when k >= (p+3)/2. (CORRECTED from (p+1)/2: at p=11, c_6=1595 is NOT divisible by C(11,6)=462, but c_7=3960 IS divisible by C(11,7)=330.)
**Results at p=7:** c_3=14 (C(7,3)=35 does NOT divide, but k=3 < 4 = (p+1)/2), c_5=42 (C(7,5)=21 DIVIDES, quotient=2), c_7=24 (C(7,7)=1 trivially divides). Conjecture HOLDS.
**Explanation:** Aut(T_p) = Z_p acts on k-subsets, partitioning them into orbits of size p (except the full vertex set which is fixed). Each k-subset orbit has the same cycle count by symmetry. So c_k = p * (cycle count per orbit) when k < p, giving p | c_k. For C(p,k) divisibility: the orbit structure under the full Aut group (which has order p*(p-1)/2 for Paley) should give the stronger divisibility.
**Next step:** Verify at p=19. Prove C(p,k) divisibility from Aut(T_p) orbit counting.

### INV-026: Alpha_1 vs |C_3| mod 2 Ã¢â‚¬â€ systematic test
**Source:** INV-011, oq:mod4_struct
**Status:** TESTED (opus-S5). CONJECTURE IS FALSE.
**Result:** Counterexamples at every n tested:
  - n=3: 2/8 counterexamples (the 3-cycle tournaments have c3=1 odd, alpha_1=0 even)
  - n=4: 16/64 counterexamples (R-cone and near-R-cone tournaments)
  - n=5: 384/1024 counterexamples
  All counterexamples have alpha_1=0 but c3 odd. The conjecture fails because alpha_1 counts independent PAIRS of odd cycles in Omega(T), which is 0 whenever #cycles <= 1.
**Impact:** Open Problem 9 needs reformulation. Alpha_1 Ã¢â€°Â  c_3 mod 2 in general.

### INV-027: Realizable conflict graphs catalog
**Source:** INV-013, conflict_graph_catalog.py
**Status:** DONE (opus-S5). Major structural finding.
**Results:**
  - n=3: 2 distinct Omega structures. n=4: 3. n=5: 6. n=6: 24. n=7 (sampled): 172.
  - **Omega(T) is ALWAYS PERFECT** (exhaustive n<=6, 2000 random n=7). This is a significant constraint Ã¢â‚¬â€ independence number = clique cover number.
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
  This means OCF is a non-perturbative identity Ã¢â‚¬â€ standard polymer expansion / cluster expansion methods CANNOT prove it.
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

### INV-042: Fano-Paley design structure and alpha_2 Ã¢â‚¬â€ MAJOR PROGRESS
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
- **Class transition graph:** Always connected. ÃŽâ€H always even. E[ÃŽâ€H]=0 for every arc position.
- **Weight distributions distinguish:** Can separate classes sharing all tournament invariants.
**Full writeup:** `03-artifacts/drafts/tiling-symmetry-analysis.md`
**Scripts:** `04-computation/tiling_*.py` (5 files)
**Next step:** (1) Investigate which tiling properties predict H beyond c3. (2) Connect sigma reduction to arc-flip proof strategy. (3) Look for grid-local rules that determine class.

### INV-037: Pin-grid sigma vs tournament sigma Ã¢â‚¬â€ two-sigma structure
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
- SC and NSC have the SAME number of 3-cycles within score class Ã¢â‚¬â€ the difference is in ARRANGEMENT (disjoint pairs vs all overlapping)
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
- T_3 Ã¢â€ â€™ T_2: H=3 Ã¢â€ â€™ H=1 = a(2) Ã¢Å“â€œ
- T_7 Ã¢â€ â€™ T_6: H=189 Ã¢â€ â€™ H=45 = a(6) Ã¢Å“â€œ
- T_11 Ã¢â€ â€™ T_10: H=95095 Ã¢â€ â€™ H=15745 = a(10) Ã¢Å“â€œ
- All vertex deletions give the same H (by Aut(T_p) transitivity)
- T_11 - v has self-complementary scores (4,4,4,4,4,5,5,5,5,5)
- H(T_p) - H(T_p-v) = 2 * (sum of mu-weighted cycles through v) (Claim A)
**Conjecture:** T_p - v is the GLOBAL H-maximizer at n = p-1 for all Paley primes p Ã¢â€°Â¡ 3 mod 4.
**Next step:** (1) Verify at p=19 (need H(T_19 - v), n=18 DP ~2^18*18 ~ 5M, feasible). (2) Test whether T_p - v is SC. (3) Relate to lattice theory or QR structure.

### INV-038: Blueself parity theorem and census structure
**Source:** opus-2026-03-06-S3 (deep census investigation)
**Status:** THM-023 PROVED. Census in progress through n=8.
**What:** Blueself (GS + self-flip) exists if and only if n is even. Proved algebraically: flip changes endpoint scores by score'(0) = n - score(0), so same-score requires score(0) = n/2 (integer only at even n).
**Census results (exhaustive n=3,...,6, in progress n=7,8):**
- POS orientation is perfectly UNIFORM: each pattern gets exactly 2^(m-#POS) tilings
- GS POS is also perfectly UNIFORM
- SC always maximizes H within each score sequence class (confirmed with kind-pasteur findings)
- Blueself at n=4: H=5 (rank 1/4), n=6: H=41,45 (ranks 5,1/56) Ã¢â‚¬â€ near or at global maximum
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
- [DONE] Locally transitive tournaments tested Ã¢â‚¬â€ DO have 5/7-cycles (opus-S4b)
- [DONE] Feng + Rajkumar connections documented in paper-connections.md (opus-S4b)
- [DONE] T_11 cycle table complete, H(T_11)=95095 confirmed (kind-pasteur-S2/S5)
- [DONE] Per-path identity failure characterized (THM-009)
- [DONE] Even-odd split is consequence not equivalent to OCF (MISTAKE-008)
- [DONE] Bracket structure B(u,w) analyzed (T047, bracket_structure.py)
- [DONE] H(T_19) computed: 1,172,695,746,915; H/|Aut|=6,857,869,865 (opus-S5)
- [DONE] Deep paper analysis: SST, Rajkumar, Feng Ã¢â‚¬â€ all key theorems extracted (opus-S5)
- [DONE] Ballot sequence bijective proof for C(L-2, 2k-1) (opus-S5)
- [DONE] Hard-core lattice gas: lambda=2 is non-perturbative regime (opus-S5)
- [DONE] Alpha_1 Ã¢â€°Â¡ c_3 (mod 2) conjecture DISPROVED (opus-S5)
- [DONE] Conflict graph catalog: Omega(T) is PERFECT for n<=7, FAILS at n>=8 (53.8% have C5 in Omega_3). See OPEN-Q-014.
- [DONE] Omega(T) is CLAW-FREE for n<=8, FAILS at n>=9 (90%). See OPEN-Q-014.
- [DONE] Web research: 9 new connections documented in web-research-connections.md (opus-S5)
- [DEAD] Per-vertex decomposition of unmatched counts (T045)
- [DEAD] Cycle bijection under arc reversal (MISTAKE-005)
- [DEAD] Contiguous block decomposition (T035)
- [DEAD] Contraction approach (T017)
- [DEAD] Alpha_1 = c_3 (mod 2) conjecture Ã¢â‚¬â€ counterexamples at all n (opus-S5)
- [DONE] Paley maximizer conjecture verified at p=3,7,11 (exhaustive at p=7), extended with H(P(19)), H(P(23)) (opus-S9)
- [DONE] n=8 H-maximizer identified: H=661=a(8), self-converse, |Aut|=1, does NOT contain P(7) (opus-S9)
- [DONE] Full Omega structure at n=8: 76-78 vertices, density 0.98, 70-75% of H from 5/7-cycles (opus-S9)
- [DONE] Ratio H(P(p))/(p!/2^{p-1}) converges toward e: 2.00, 2.40, 2.44, 2.53, 2.56 for p=3,7,11,19,23 (opus-S9)
- [DONE] Deep web synthesis: Irving-Omar, Grinberg-Stanley, GrujiÃ„â€¡-StojadinoviÃ„â€¡, Feng, DRT theory (opus-2026-03-06-S5)

---

## Priority F: New leads from web synthesis (opus-2026-03-06-S5)

### INV-045: Hopf algebra route to transfer matrix symmetry Ã¢â‚¬â€ SUPERSEDED by THM-030
**Source:** T114, GrujiÃ„â€¡-StojadinoviÃ„â€¡ arXiv:2402.07606, Feng arXiv:2510.25202
**Status:** INVESTIGATED (kind-pasteur-S22 agent). No direct proof obtained.
**What:** The Hopf comultiplication ÃŽâ€([T]) = ÃŽÂ£_S [T|_S]Ã¢Å â€”[T|_{V\S}] encodes our subset convolution. Feng's dual Burnside process proves Q=AB is symmetric under detailed balance. The tournament constraint T[x,y]+T[y,x]=1 is the detailed balance condition.
**Finding:** Feng's framework FAILS because it requires positivity (probability transitions), but our M = E^T * Lambda * B has Lambda = diag((-1)^|S|) with negative entries. Grujic-Stojadinovic gives U_X = U_{X^op} (proved) but this is a GLOBAL identity, not per-entry M[a,b]=M[b,a]. The Hopf comultiplication IS cocommutative but this doesn't directly imply transfer matrix symmetry because E_a and B_b are different types of objects. **Most promising remaining paths:** (1) inductive s-variable approach, (2) Irving-Omar det/per formula (INV-046), (3) "signed Feng" extension (new theorem needed), (4) pointed Hopf algebra tracking distinguished vertices.
**Why it matters:** Transfer matrix symmetry Ã¢Å¸Â¹ OCF (via even-odd split + s-coefficient identity). A Hopf algebra proof would be clean and publication-worthy.
**Next step:** (1) Express our transfer matrix M in Feng's AB framework precisely. (2) Verify detailed balance condition algebraically. (3) Apply Feng Thm 3.3 to conclude symmetry.

### INV-046: Irving-Omar det/per formula Ã¢â€ â€™ transfer matrix Ã¢â‚¬â€ SUPERSEDED by THM-030
**Source:** T118, Irving-Omar arXiv:2412.10572 Proposition 2
**Status:** SUPERSEDED. Transfer matrix symmetry now proved directly (THM-030, opus-S25) without needing Irving-Omar det/per.
**What:** ham(D) = ÃŽÂ£_S det(Ã„â‚¬[S])Ã‚Â·per(A[S^c]). The walk generating function W_D(z)=det(I+zXÃ„â‚¬)/det(I-zXA) encodes Hamiltonian path structure.
**Remaining interest:** Irving-Omar's framework may still provide insight into WHY the Key Identity works Ã¢â‚¬â€ e.g., is there a matrix-algebraic interpretation of B_b + (-1)^m E_b = 2rÃ‚Â·col_sum?

### INV-047: Paley maximizer via DRT theory
**Source:** T116, Reid-Brown 1972, Nozaki-Suda arXiv:1202.5374
**Status:** CLASSICAL EQUIVALENCE + our computational evidence.
**What:** DRTs Ã¢â€ â€ skew Hadamard matrices. Paley T_p is the canonical DRT. Nozaki-Suda characterize skew Hadamard via spectra of tournaments of size n-2. Our spectral regularity finding (corr(H,ÃŽÂ»Ã¢â€šÂ)=-0.97) + DRT theory could explain why Paley maximizes H.
**Key question:** Among all DRTs on p vertices, does Paley ALWAYS maximize H? Or is this specific to Paley among all tournaments?
**Next step:** (1) Check if non-Paley DRTs exist at small p and compare H values. (2) Relate DRT cycle balance to alpha_k maximization.

### INV-048: Asymptotic convergence H(T_p)/(p!/2^{p-1}) Ã¢â€ â€™ e
**Source:** T117, Adler-Alon-Ross 2001
**Status:** COMPUTATIONAL EVIDENCE. Not proved.
**What:** Adler-Alon-Ross proved max H(T) Ã¢â€°Â¥ (e-o(1))Ã‚Â·n!/2^{n-1} using random regular tournaments. Our Paley ratios 2.00Ã¢â€ â€™2.56 suggest convergence to e. Paley tournaments are quasi-random in Chung-Graham-Wilson sense.
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

### INV-051: Mitrovic noncommuting RÃƒÂ©dei-Berge function (arXiv:2504.20968, Apr 2025)
**Source:** kind-pasteur-2026-03-06-S19 web search; opus-2026-03-06-S9 (detailed paper read)
**Status:** DEEPLY READ Ã¢â‚¬â€ HIGH PRIORITY
**What:** Stefan Mitrovic introduces the RÃƒÂ©dei-Berge function in NONCOMMUTING variables, which satisfies deletion-contraction (W_X = W_{X\e} - W_{X/e}Ã¢â€ â€˜). The commutative version does NOT have deletion-contraction. Key properties: W_X = W_{X^op}, product rule W_{XÃ‚Â·Y} = W_XÃ‚Â·W_Y. For tournaments: W_X = sum over permutations with all odd cycles of 2^{psi(sigma)} p_{Type(sigma)} with positive integer coefficients.
**Why it matters:** Deletion-contraction enables INDUCTIVE PROOFS. This could provide an inductive framework for OCF or transfer matrix symmetry. The noncommutative structure preserves more information than the commutative version.
**TESTED (kind-pasteur S19):** Direct deletion-contraction does NOT preserve OCF. At n=4:
  - H(T) = H(T\e) - H(T/e): only 18.8% match (DC is for W_X not H)
  - OCF for T\e (non-tournament): only 39.3% hold
  - OCF for T/e (contracted): only 60.7% hold
  OCF is TOURNAMENT-SPECIFIC and does not hold for general digraphs from deletion/contraction.
  The noncommuting framework operates at a different level than H(T).
**DETAILED READING (opus-S9):**
  - Definition 3.1: W_X = sum_{f:V->P} sum_{sigma in Sigma_V(f,X)} delta_f(sigma) x_{f(v1)}...x_{f(vn)} where x_i are NONCOMMUTING. Depends on vertex labeling (unlike commutative U_X).
  - Theorem 3.7 (Deletion-Contraction): W_X = W_{X\e} - W_{X/e}Ã¢â€ â€˜ where e=(v_{n-1},v_n). The Ã¢â€ â€˜ operation doubles the last variable: (x_{i1}...x_{in-1})Ã¢â€ â€˜ = x_{i1}...x_{in-2} x_{in-1}^2. This is the KEY technical device.
  - Theorem 3.10: W_X = sum_{sigma in S_V(X,Xbar)} (-1)^{phi(sigma)} p_{Type(sigma)} where Type is now a SET PARTITION (not integer partition). The noncommutative p_pi tracks which vertices are in which cycle.
  - Corollary 3.12: For tournaments, W_X = sum_{sigma in S_V(X), all odd cycles} 2^{psi(sigma)} p_{Type(sigma)}. This is the NONCOMMUTATIVE OCF Ã¢â‚¬â€ same sum but tracking set-partition cycle types.
  - Key insight: The Ã¢â€ â€˜ operation has no obvious tiling interpretation, but it corresponds to "merging the last two vertices while remembering which was which." This is precisely what our contraction approach (T017) failed to handle in commutative setting.
**Next step:** (1) Instead of naive DC -> OCF induction, investigate whether DC can be used to relate the SYMMETRIC FUNCTION U_T across tournaments (e.g., U_T = U_{T'} + correction for single arc reversal). (2) Check if the noncommutative framework gives a new proof of even-odd split or transfer matrix symmetry. (3) CRITICAL: Investigate whether Theorem 3.10 (noncommutative power-sum expansion over SET partitions) implies transfer matrix symmetry when projected to integer partitions.

### INV-052: Mitrovic-Stojadinovic chromaticÃ¢â€ â€RÃƒÂ©dei-Berge connection (arXiv:2506.08841, Jun 2025)
**Source:** kind-pasteur-2026-03-06-S19 web search; opus-2026-03-06-S9 (FULL PAPER READ)
**Status:** DEEPLY READ Ã¢â‚¬â€ HIGH PRIORITY. Contains multiple directly relevant results.
**What:** Proves the chromatic function of a graph and the RÃƒÂ©dei-Berge function of a digraph are "almost identical" at the poset level. For any poset P: X_{inc(P)} = omega(U_P) (Theorem 3.2). This extends to noncommutative versions: Y_{inc(P)} = omega(W_P) (Theorem 5.7). Commutative diagram of four Hopf algebra morphisms (p.8, Remark 3.4).
**KEY RESULTS FROM FULL READ (opus-S9):**
  1. **Theorem 4.1 (Converse of RÃƒÂ©dei):** If P is a poset that is NOT a chain, then the number of quasi-linear extensions is EVEN. Proof: chi_{inc(P)}(1) and (-1)^|P| u_P(-1) have the same parity. u_P(1) counts quasi-linear extensions; chi_{inc(P)}(1)=0 unless inc(P) is discrete (P is a chain). THIS GENERALIZES RÃƒâ€°DEI'S THEOREM TO POSETS.
  2. **Theorem 4.8:** For any poset P and partition lambda: #{broken-cycle-free subsets of inc(P) with component sizes lambda} = #{permutations in S_V(D_P-bar) with cycle type lambda}. This is a CYCLE-TYPE-PRESERVING bijection between broken cycles and permutations Ã¢â‚¬â€ potentially a route to bijective OCF proof.
  3. **Corollary 3.3:** chi_{inc(P)}(m) = (-1)^|P| u_P(-m). The chromatic polynomial at m colors = RÃƒÂ©dei-Berge polynomial at -m. Our H(T) = u_T(1) = (-1)^n chi_{inc(P)}(-1) for the associated poset.
  4. **Theorem 4.6:** D_P-bar contains a Hamiltonian cycle iff P is irreducible. Combined with Corollary 4.5: U_{P_n} forms an algebraic basis of Sym for any sequence of irreducible posets.
  5. **Section 6 (Positivity):** Conjecture 6.3: If P is (3+1)-free, U_P is h-positive. Theorem 6.4: Already proved s-positive. Theorem 6.11: h-positivity propagates through deletion-contraction if there's a "sink-like" vertex.
  6. **Section 7 (Bags of sticks):** Decomposition into bags of sticks (disjoint unions of directed paths) gives explicit formulas. Triple deletion property generalized.
**Why it matters for us:**
  - The bridge X_{inc(P)} = omega(U_P) means we can use 30 years of chromatic symmetric function theory on tournament problems
  - Theorem 4.1 is a clean generalization of the parity theorem we're studying
  - Theorem 4.8 gives a TYPE-PRESERVING bijection Ã¢â‚¬â€ this is exactly what a bijective OCF proof would need
  - The h-positivity results (Section 6) may apply to tournament-associated posets
**Next step:** (1) For a tournament T, identify the associated poset P such that D_P = T. (2) Check if Theorem 4.8 gives a new proof of OCF when specialized to tournaments. (3) Investigate whether tournament posets are (3+1)-free (this would give h-positivity of U_T). (4) Test computationally: does the broken-cycle bijection from Thm 4.8 match our OCF terms?

### INV-053: Savchenko cycle counting formulas for regular tournaments Ã¢â‚¬â€ VERIFIED AT n=7
**Source:** kind-pasteur-2026-03-06-S19 web search; Savchenko J. Graph Theory 83 (2016), Discrete Math (2017), arXiv:2403.07629 (2024)
**Status:** VERIFIED at n=7. Cycle counts are class invariants. DRT vs LTT classification matches.
**What:** Savchenko has a series of papers giving EXACT polynomial formulas for c_k(T) (number of k-cycles) in regular tournaments:
- c5, c6 formulas (2016, J. Graph Theory 83)
- c7 formula (2017, Discrete Math)
- c8 for DRTs vs locally transitive tournaments (2024, arXiv:2403.07629)
Key finding: c8(DRT_n) is INDEPENDENT of which DRT is chosen. Phase transition at n=39: DRTs have more 8-cycles than locally-transitive for nÃ¢â€°Â¤35 but FEWER for nÃ¢â€°Â¥39.
**Why it matters:** These exact formulas could determine whether Paley tournaments maximize cycle counts at EVERY length, or only for short cycles. The phase transition at n=39 suggests our cycle-maximization mechanism may reverse at larger n. Also, the spectral methods used (eigenvalue-based cycle counting) could connect to our transfer matrix work.
**VERIFIED (kind-pasteur S19):** At n=7, the three regular tournament classes are EXACTLY:
  - DRT (Paley): 240 tours, dc={3:14, 5:42, 7:24}, H=189
  - Locally Transitive: 720 tours, dc={3:14, 5:28, 7:17}, H=175
  - Other Regular: 1680 tours, dc={3:14, 5:36, 7:15}, H=171
Cycle counts are CLASS INVARIANTS (exactly one vector per class). DRT maximizes directed 5-cycles and 7-cycles. LTT has "diametrically opposite" properties per Savchenko.
**EXTENDED (kind-pasteur S21):** Savchenko (2024) proves c_m(DR_n) > c_m(RLT_n) for ALL m = 1,2,3 mod 4 (including all odd m). Only m = 0 mod 4 has the phase transition. This directly explains DRT's H-maximization via OCF.
**DRT n=11 ANALYSIS (kind-pasteur S21, CORRECTED S39b per MISTAKE-017):** The connection set {1,2,3,5,8} is NOT a valid tournament (SÃ¢Ë†Â©(-S)={3,8}Ã¢â€°Â Ã¢Ë†â€¦). The ONLY valid circulant DRT at n=11 is the Paley tournament QR={1,3,4,5,9} (H=95095, c3=55, c5=594, |Aut|=55). All claims about "non-Paley DRT" H=69311, c3=44 were computed on an invalid digraph. Whether a non-circulant DRT exists at n=11 remains open.
**Next step:** (1) Obtain Savchenko's exact polynomial c_k formulas. (2) Test at n=19 or n=23 (multiple DRT classes). (3) Prove Paley maximizes H among ALL DRTs.

### INV-054: Komarov-Mackey exact 5-cycle formula (arXiv:1410.6828, JGT 2017) Ã¢â‚¬â€ PARTIALLY INVESTIGATED
**Source:** kind-pasteur-2026-03-06-S19 web search
**Status:** PARTIALLY INVESTIGATED (kind-pasteur-S39b). c5 is NOT score-determined.
**What:** Exact formula for number of directed 5-cycles in any tournament in terms of edge score sequence. Maximum c5 Ã¢â€°Ë† (3/4)*C(n,5), achieved by almost all random tournaments. Lower bounds also proved.
**NEW FINDING (S39b):** c5 is NOT determined by score sequence, even combined with sum_dÃ‚Â², edge_score, or common_out_neighbor statistics. At n=5, score (1,2,2,2,3) has c5 in {1,2,3}; at n=6, 9/22 score sequences have varying c5. The Komarov-Mackey formula likely involves CUBIC or higher-order graph statistics (e.g., directed walks of length 3+). For regular tournaments, c5 IS a class invariant (Savchenko, verified n=5,7).
**Why it matters:** This rules out O(nÃ‚Â²) c5 computation from scores alone. Cycle enumeration (O(n^5) for 5-cycles, or O(2^n) bitmask DP) remains necessary.
**RESOLVED (S39b, THM-118):** c_5 = tr(A^5)/5 gives O(n^3) computation via matrix multiplication. This IS the "cubic invariant" Ã¢â‚¬â€ tr(A^5) is a sum over all length-5 closed walks, and THM-118 proves all such walks in tournaments are simple cycles (no vertex repetition possible for length <= 5).
**Next step:** Read Komarov-Mackey formula to see if it matches tr(A^5)/5.

### INV-055: Linial-Morgenstern cycle density conjecture and extremal tournaments
**Source:** kind-pasteur-2026-03-06-S19 web search; arXiv:2011.14142 (Ma-Tang), arXiv:1902.00572
**Status:** NEW LEAD Ã¢â‚¬â€ MEDIUM PRIORITY
**What:** Linial-Morgenstern conjecture: among tournaments with fixed c3 density d, the c4 density is minimized by random blowups of transitive tournaments. Proved for d Ã¢â€°Â¥ 1/36 using spectral methods. Ma-Tang extend to c_Ã¢â€žâ€œ for Ã¢â€žâ€œ Ã¢â€°Â¢ 2 mod 4 when d is near 1.
**Why it matters:** This is the "dual" to our maximization question. We show Paley maximizes total directed cycles; this literature characterizes minimizers. The spectral methods used here (eigenvalue-based cycle density bounds) could provide tools for our Paley maximizer proof.
**Next step:** Check if the extremal results constrain H(T) via OCF.

### INV-056: Jerrum-Patel zero-free regions for H-free graphs (JLMS 2026)
**Source:** kind-pasteur-2026-03-06-S19 web search
**Status:** NEW LEAD Ã¢â‚¬â€ MEDIUM PRIORITY (for real-rootedness question)
**What:** Jerrum & Patel (2026, JLMS) prove zero-free regions for the independence polynomial of H-free graphs for various H. For claw-free: all zeros on negative real line (= Chudnovsky-Seymour). For subdivided claws: related zero-free regions. KEY: for H NOT a subdivided claw or path, there exist H-free graphs of max degree 3 with zeros NOT on the negative real line.
**Why it matters:** Our Omega_3(T) has all real roots for nÃ¢â€°Â¤20 but is NOT always claw-free (fails nÃ¢â€°Â¥9). Jerrum-Patel's results on subdivided claw avoidance may explain why real roots persist beyond n=8. The tournament-specific constraint on Omega_3 structure may ensure avoidance of exactly the "bad" subgraphs.
**Next step:** (1) Check what specific subdivided claws appear in Omega_3(T) at nÃ¢â€°Â¥9. (2) Apply Jerrum-Patel to determine if their zero-free regions explain our observations.

### INV-057: Herman's Terwilliger algebras of DRTs (arXiv:2404.11560, 2024)
**Source:** kind-pasteur-2026-03-06-S19 web search
**Status:** NEW LEAD Ã¢â‚¬â€ LOW-MEDIUM PRIORITY
**What:** Allen Herman computes Terwilliger algebras for DRTs (asymmetric rank-3 association schemes). Thin irreducible modules, dimension 4k+9. Key: Terwilliger algebras distinguish non-isomorphic DRTs up to n=23, but FAIL at n=27 (need rational Terwilliger algebras). There are 237 non-isomorphic DRTs at n=27.
**Why it matters:** (1) If all DRTs at small n have the same H(T), that would be a DRT invariant. (2) If Terwilliger algebra structure constrains H(T), this gives an algebraic route to Paley maximization. (3) The n=27 DRT classification gives test cases for our conjectures beyond Paley primes.
**Next step:** (1) Check if all DRTs at n=7 (there's only one: Paley) or n=11 have the same H. (2) At n=27, compare H across different DRT isomorphism classes.

### INV-058: Pantangi critical groups distinguish Paley from other DRTs
**Source:** kind-pasteur-2026-03-06-S19 web search; Pantangi arXiv:1905.08568 (2019)
**Status:** CONNECTION IDENTIFIED
**What:** Pantangi shows critical groups (sandpile groups) distinguish Paley from non-Paley DRTs. Chandler-Sin-Xiang computed Smith/critical groups of Paley GRAPHS. Different DRT constructions (Szekeres-Whiteman 2-block, Wallis-Whiteman 4-block) are distinguished by their critical groups.
**Why it matters:** If H(T) is a DRT invariant AND different DRTs have different critical groups, then H could be read off the critical group. This would give a purely algebraic characterization of the H-maximizer.
**Next step:** Compute critical groups for DRTs at n=11,19 and check correlation with H values.

### INV-062: Universal Master Polynomial and Central Factorial Numbers Ã¢â‚¬â€ VERIFIED (THM-059)
**Source:** opus-2026-03-06-S30
**Status:** VERIFIED computationally (22/22 cases, n=4..9). Algebraic proof pending.
**What:** The per-invariant r-polynomial C_I(r,n) = 2^{parts(I)} * F_f(r) where f = free position count and F_j is determined by the central factorial number triangle (OEIS A036969) via b_{k,j} = b_{k-1,j-1} + j^2 * b_{k-1,j}. This completely determines the entire W-coefficient hierarchy for all tournaments at all n. The shift principle is a corollary.
**Key findings:** (1) F_j(1/2) = 1 for all j. (2) Leading coefficient = (j+1)!. (3) Predictions made for n=11 without computation. (4) Complete n=8 even-n table computed.
**Next step:** (1) Algebraic proof of the central factorial recurrence from position pattern analysis. (2) Verify F_8 prediction at n=11 computationally. (3) Investigate C_0(r) (constant/background polynomial).
**Scripts:** `04-computation/universal_master_polynomial.py`, `04-computation/w1_n8_complete.py`

### INV-059: Cyclic subsets of tournaments (arXiv:2508.03634, Aug 2025)
**Source:** kind-pasteur-2026-03-06-S19 web search; Hunter-Liu-MilojeviÃ„â€¡-Sudakov
**Status:** NEW LEAD Ã¢â‚¬â€ LOW PRIORITY
**What:** Optimal bounds on probability that random induced subtournament of high min-degree tournament is Hamiltonian. Extends to p-biased measure. Proves that high min-degree forces high Hamiltonicity probability.
**Why it matters:** Paley T_p has min-degree (p-1)/2 (doubly regular). This paper could give explicit bounds on the fraction of induced subtournaments that are Hamiltonian, which connects to our cycle counting.
**Next step:** Apply their bounds to Paley tournaments. Check if this gives lower bounds on c_k counts.

### INV-060: Eulerian cycle trace formula (arXiv:2502.02915, Feb 2025) Ã¢â‚¬â€ CLOSED
**Source:** kind-pasteur-2026-03-06-S19 web search; Ye Luo
**Status:** CLOSED (too remote, kind-pasteur-S22 agent investigation)
**What:** Trace formula counting Eulerian cycles via "twisted" vertex and edge adjacency matrices. Uses homological spectral graph theory.
**Finding:** Connection too remote to be actionable. Luo counts Eulerian cycles (edge traversals) via H_1 characters; we count Hamiltonian paths (vertex traversals) via inclusion-exclusion. Different domain, different algebraic structure. Closest classical analogue to our transfer matrix is Ryser's permanent formula, not Luo's trace formula.
**Why it matters:** Our transfer matrix tr(M) = H(T) (THM-027) is also a trace formula. This paper's approachÃ¢â‚¬â€using twisted adjacency matrices with spectral antisymmetryÃ¢â‚¬â€could provide a template for proving our trace formula properties (symmetry, off-diagonal sum) at general n.
**Next step:** Read the paper and check if "twisted adjacency" techniques apply to tournament transfer matrices.

### INV-061: Hamilton transversals in tournaments (Combinatorica 2024)
**Source:** kind-pasteur-2026-03-06-S19 web search; Chakraborti-Kim-Lee-Seo arXiv:2307.00912
**Status:** NEW LEAD Ã¢â‚¬â€ LOW PRIORITY
**What:** For collections of sufficiently many tournaments on the same vertex set, transversal Hamilton paths/cycles exist. For m=|V|-1 tournaments, there's a transversal Ham path; for m=|V| with m-1 strongly connected, transversal Ham cycle.
**Why it matters:** The "transversal" perspective could give a new way to relate Ham paths across different tournaments, potentially connecting to how H(T) changes under arc reversals.

### INV-062: Forward arc maximization in tournaments (arXiv:2602.10713, Feb 2026)
**Source:** kind-pasteur-2026-03-06-S19 web search; Guo-Gutin-Lan-Shao-Yeo-Zhou
**Status:** NEW LEAD Ã¢â‚¬â€ LOW PRIORITY
**What:** Characterizes maximum forward arcs in Hamilton cycles/paths for semicomplete and locally semicomplete digraphs. Polynomial-time algorithms.
**Why it matters:** Forward arcs in Hamilton paths relate to our "position-based" analysis (pos(a,P) in THM-027 trace formula). The maximum forward arc structure could inform transfer matrix properties.

### INV-063: Spectral pseudorandomness and Paley clique bounds (Exp. Math. 2024)
**Source:** kind-pasteur-2026-03-06-S19 web search; Kunisky-Yu arXiv:2303.16475
**Status:** NEW LEAD Ã¢â‚¬â€ LOW PRIORITY
**What:** Studies spectral pseudorandomness of Paley graphs via subgraph eigenvalue distributions. Conjecturally, minimum eigenvalue convergence would improve clique number bounds beyond Ã¢Ë†Å¡p.
**Why it matters:** Spectral properties of Paley graphs/tournaments are central to our theory. If Paley tournaments have stronger spectral pseudorandomness than other DRTs, this could explain H-maximization via eigenvalue-based cycle counting formulas.

### INV-064: Mitrovic Hopf algebra new bases (arXiv:2407.18608v3, Mar 2025)
**Source:** kind-pasteur-2026-03-06-S19 web search
**Status:** CONNECTION IDENTIFIED Ã¢â‚¬â€ supplements INV-033
**What:** Introduces two new combinatorial Hopf algebras of posets and permutations with RÃƒÂ©dei-Berge functions. Constructs new bases for symmetric functions whose generators are RÃƒÂ©dei-Berge functions. Investigates which digraph invariants are extractable from the RÃƒÂ©dei-Berge function.
**Why it matters:** If H(T) can be expressed as a coefficient in one of these new bases, it gives an algebraic handle on Hamiltonian path counting.
**Next step:** Check which digraph invariants the paper extracts. Is H(T) among them?

### INV-065: Independence polynomial root gap (arXiv:2510.09197, FSTTCS 2025)
**Source:** kind-pasteur-2026-03-06-S19 web search; Om Prakash & Vikram Sharma
**Status:** NEW LEAD Ã¢â‚¬â€ LOW PRIORITY
**What:** Quantifies the gap between the smallest real root ÃŽÂ²(G) of I(G,x) and all other roots. For connected graphs, ÃŽÂ²(G) is a simple real root smaller than 1, but previous proofs gave no gap bound. This paper provides explicit bounds.
**Why it matters:** For our Omega(T) real-rootedness question, having a gap bound could help prove that all roots are real by showing they're well-separated from the complex plane.

### INV-066: Low-rank matrices from tournaments and symmetric designs (arXiv:2401.14015, 2024)
**Source:** kind-pasteur-2026-03-06-S19 web search; Balachandran-Sankarnarayanan
**Status:** NEW LEAD Ã¢â‚¬â€ LOW-MEDIUM PRIORITY
**What:** Constructs symmetric matrices from tournament structures where rank depends on design-theoretic properties. Symmetric designs (BIBDs) give matrices with rank near n/2. The rank-topology relationship involves bipartite graph eigenvalues.
**Why it matters:** Our transfer matrix M is constructed from a tournament and is symmetric. This paper's framework connecting tournament-derived matrices with design theory could explain structural properties of M (e.g., why symmetry holds, what the rank structure is).
**Next step:** Check if our M fits their M_T(f,a) framework.

### INV-067: Alpha_1 gap theorem and converse of Redei Ã¢â‚¬â€ CORRECTED S22
**Source:** kind-pasteur-2026-03-06-S21 (computation), CORRECTED S22
**Status:** PARTIALLY PROVED (THM-029 corrected)
**What:** alpha_1=3 is impossible at n<=6 but ACHIEVABLE at n>=7 (~9.2% of c3=3 at n=7). Common-vertex property fails at n>=7. HOWEVER, H=7 remains impossible for ALL n by refined argument: H=7 requires (alpha_1=3, i_2=0), but i_2=0 forces common vertex => c5>=1 => alpha_1>=4; while alpha_1=3 implies i_2>=1 => H>=11. H=21 absent through n=7.
**Achievable H values:** n=5: {1,3,5,9,11,13,15}. n=6: {1,3,5,9,11,...,45}\{7,21,35,39}. At n=7: 35 and 39 become achievable but 7 and 21 remain gaps.
**Connection:** Relates to Mitrovic-Stojadinovic "converse of Redei" (INV-052).
**Next step:** (1) Prove alpha_1=3 impossibility for ALL n (not just n<=6). (2) Find all permanent H-gaps. (3) Check OEIS for the sequence of achievable alpha_1 values.
**Scripts:** h7_impossibility.py, alpha1_gaps.py, alpha1_gap3_proof.py, c3_forces_c5.py, redei_converse_fast.py
**Theorem:** THM-029

### INV-068: DRT non-uniqueness and Paley dominance at n=11 Ã¢â‚¬â€ CORRECTED (MISTAKE-017)
**Source:** kind-pasteur-2026-03-06-S21 (computation), CORRECTED kind-pasteur-2026-03-07-S39b
**Status:** CORRECTED Ã¢â‚¬â€ previous "non-Paley DRT" was INVALID
**What:** Previous claim of "2 DRT classes at n=11" was based on connection set {1,2,3,5,8} which is NOT a valid tournament (S Ã¢Ë†Â© (-S) = {3,8} Ã¢â€°Â  Ã¢Ë†â€¦, creating bidirectional edges). ALL claims about "non-Paley DRT" cycle counts (c3=44, c5=407, H=69311) are INVALID (MISTAKE-017).
**Corrected facts:** The only valid tournament (11,5,2)-difference sets in Z_11 are {1,3,4,5,9} (QR) and {2,6,7,8,10} (NQR), which give isomorphic Paley tournaments. There is no non-Paley circulant DRT at n=11. Whether a non-circulant DRT exists at n=11 remains open.
**Paley T_11 correct data:** H=95095, c3=55, c5=594, c7=3960, c9=11055, c11=5505, |Aut|=55.
**Next step:** (1) Check literature for non-circulant DRT existence at n=11 (all groups of order 11 are Z_11, so non-circulant DRT must be non-Cayley). (2) Test DRT uniqueness at n=23 where multiple constructions are known.
**Scripts:** drt_n11_analysis.py (CONTAINS BUG Ã¢â‚¬â€ uses invalid connection set), drt_n11_verify.py (correction script)

### INV-069: Scalar M characterization Ã¢â‚¬â€ M=(H/n)*I Ã¢â€ â€ H-maximizer at n=5, Ã¢â€ â€ VT at n=7
**Source:** kind-pasteur-2026-03-06-S25c (T156), opus-2026-03-06-S26 (T158)
**Status:** COMPUTED, OPEN CONJECTURE
**What:** Transfer matrix M[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b) is scalar (= (H/n)*I) for certain tournaments. At n=5: ALL 64 H-maximizers (H=15) have M=3*I, including 40 non-regular (scores 1,2,2,2,3) with |Aut|=3 (Z/3Z). These non-VT tournaments still have uniform position distribution. At n=7: M scalar Ã¢Å¸Âº VT (100 sampled regular, 0 counterexamples). All circulant tournaments have scalar M. The characterization seems to be "uniform position distribution" (every vertex appears equally often at each position in Ham paths), which is weaker than VT at n=5 but equivalent at n=7.
**Key results:**
- n=5: M scalar Ã¢Å¸Â¹ H=15 (max), but NOT Ã¢Å¸Â¹ VT or regular
- n=7: M scalar Ã¢Å¸Âº VT (conjecture, verified 100 samples)
- All circulant Ã¢Å¸Â¹ M scalar (verified n=3,5,7)
- Paley T_7: |Aut|=21 (Frobenius Z/7ZÃ¢â€¹Å Z/3Z), sigma(x)=2x mod 7 is extra aut
- opus finding: non-regular n=5 scalar M has Z/3Z aut, uniform endpoint counts, uniform 3-cycle participation
**NONHAM VANISHING (kind-pasteur S25c):** The split pair decomposition M[a,b] = HAM(a,b) + NONHAM(a,b) shows:
- NONHAM(a,b) = 0 trivially when T[a,b]=1 (all split pairs are Ham paths)
- NONHAM(a,b) = M[a,b] when T[a,b]=0 (junction always fails)
- So NONHAM=0 for all (a,b) Ã¢Å¸Âº M[a,b]=0 whenever T[a,b]=0
- VERIFIED: NONHAM=0 for ALL position-uniform n=5 (64/64 exhaustive) and ALL circulant n=7 (8/8)
- NONHAM Ã¢â€°Â  0 for general (non-uniform) tournaments at n=3,4,5
**PROOF CHAIN (verified n=3,5,7):**
1. Position-uniform Ã¢Å¸Â¹ NONHAM=0 Ã¢Å¸Â¹ M[a,b]=0 when T[a,b]=0
2. THM-030: M[a,b]=M[b,a]
3. For T[a,b]=1: T[b,a]=0 Ã¢Å¸Â¹ M[b,a]=0 Ã¢Å¸Â¹ M[a,b]=0 by symmetry
4. M[a,a]=H/n from position uniformity
5. M=(H/n)*I. QED (modulo proving step 1 at general n)
**Algebraic proof of step 1:** OPEN. The cancellation mechanism: for uniform T, nonzero E*B products at adjacent subset sizes pair up and cancel. For non-uniform T, "orphan" terms remain.
**Connection to BjÃƒÂ¶rklund et al.:** M[a,b] is a subset convolution in their "Fourier meets MÃƒÂ¶bius" framework.
**Next step:** (1) Prove NONHAM=0 for position-uniform tournaments at general n. (2) Test at n=9 (circulant only, n=9 exhaustive infeasible). (3) Find algebraic proof of the pairing mechanism.
**Scripts:** `04-computation/transfer_matrix_maximizers.py`, `04-computation/circulant_scalar_m_proof.py`, `04-computation/scalar_m_characterization.py`, `04-computation/scalar_m_n5_analysis.py`, `04-computation/circulant_m_scalar_proof.py`, `04-computation/nonham_vanish_general.py`, `04-computation/nonham_vanish_uniform.py`, `04-computation/nonham_proof_analysis.py`

### INV-070: Fibonacci determinant formula for transitive tournament
**Source:** kind-pasteur-2026-03-06-S25c (T157)
**Status:** VERIFIED n=2,...,11
**What:** For the transitive tournament T_n (i beats j iff i<j), det(M) = (-1)^{n(n-1)/2} * F(n+1) where F is Fibonacci. The matrix D*M = I+U-L (tridiagonal: 1s on diagonal/superdiagonal, -1s on subdiagonal) satisfies the Fibonacci recurrence. The Chebyshev eigenvalue conjecture (eigenvalues = 2cos(kÃâ‚¬/(n+1))) is FALSE.
**Connection:** The transitive tournament is the unique acyclic tournament, so Omega(T_n) is empty and I(Omega,2)=1=H(T_n). The transfer matrix has a clean tridiagonal structure reflecting the total order.
**Next step:** (1) Does the Fibonacci structure extend to near-transitive tournaments? (2) What is det(M) for other named tournament families?
**Scripts:** `04-computation/fibonacci_determinant_proof.py`

### INV-071: det(M) as tournament invariant Ã¢â‚¬â€ exhaustive n=5
**Source:** kind-pasteur-2026-03-06-S25c (T157)
**Status:** COMPUTED
**What:** At n=5, det(M) takes values {-27,-9,1,8,9,16,17,243} across 10 isomorphism classes (up to converse). The 10 eigenvalue patterns of M classify tournaments into exactly the 10 isomorphism classes. det(M) is NOT determined by (H,c3) Ã¢â‚¬â€ e.g., H=5,c3=2 gives both det=9 and det=17.
**Key finding:** det(M) is not related to det(A), per(A), or any simple matrix expression of the adjacency matrix A. Tested det(IÃ‚Â±S), det(A+A^T), det(A*A^T), det(I+A*A^T) Ã¢â‚¬â€ none match. The transfer matrix is fundamentally different from the adjacency matrix.
**opus finding:** char poly det(ÃŽÂ»I-M) encodes hierarchy of path correlations. Scalar M Ã¢Å¸Â¹ char poly = (ÃŽÂ»-H/n)^n. PSD threshold at n=5: HÃ¢â€°Â¥13 Ã¢Å¸Â¹ M is PSD.
**Next step:** (1) Compute det(M) exhaustively at n=7 (sample). (2) Find if any spectral graph invariant matches. (3) Investigate PSD threshold at general odd n.
**Scripts:** `04-computation/det_m_general_formula.py`, `04-computation/det_m_vs_adjacency_spectrum.py`, `04-computation/char_poly_M_analysis.py`, `04-computation/pos_skeleton_connection.py`

### INV-072: IO walk GF vs transfer matrix bridge
**Source:** opus-2026-03-06-S26, kind-pasteur-2026-03-06-S25c
**Status:** STRUCTURAL COMPARISON, NO BRIDGE FOUND
**What:** Irving-Omar walk GF W_D(z) = det(I+zXA^T)/det(I-zXA) uses cycle covers (det/per), while M[a,b] uses path decomposition (E_a*B_b). M[a,b] Ã¢â€°Â  H(aÃ¢â€ â€™b) (direct endpoint-conditioned count). W(-z,-r) = W(z,r) at commutative level (opus finding). No simple matrix expression of A gives M.
**Key structural difference:** IO is multiplicative (det/per = products over cycles), M is additive (sum over subsets with inclusion-exclusion signs). Bridge might exist through Hopf algebra coproduct structure.
**Next step:** (1) Express M[a,b] in terms of IO's det/per framework. (2) Check if deletion-contraction on M matches Mitrovic's W_X = W_{X\e} - W_{X/e}^up.

### INV-073: Palindromic N => Scalar M for circulant tournaments Ã¢â‚¬â€ PROVED (THM-052)
**Source:** kind-pasteur-2026-03-06-S25d, building on THM-050 (opus-S26), THM-051 (opus-S26)
**Status:** PROVED for all circulant tournaments at odd n.
**What:** For circulant T on Z/nZ, the consecutive-position count N(a,b,j) = f(b-a mod n, j) is palindromic: f(d,j) = f(d,n-2-j). The proof uses three ingredients: (1) translation symmetry gives f depends only on d, (2) N symmetry gives f(d)=f(n-d), (3) self-complementarity via sigma: i->-i gives f(d,j)=f(n-d,n-2-j). Combining (2)+(3): f(d,j)=f(d,n-2-j). At odd n, palindromic N forces alternating sum = 0, so M[a,b]=0 for a!=b. Combined with M[a,a]=H/n, gives M=(H/n)*I.
**Verification:** n=5 (64/64 exhaustive pos-uniform), n=7 (8/8 circulant), n=9 (16/16), n=11 Paley, n=13 Paley.
**Key finding:** ALL position-uniform n=5 tournaments are self-complementary (64/64). No position-uniform tournaments exist at even n (n=4,6).
**Open extension:** Prove for general vertex-transitive (non-circulant) tournaments. The proof uses circulant-specific translation symmetry. At n=15, non-circulant VT tournaments exist (Babai-Kantor doubly-regular tournaments).
**Scripts:** `04-computation/palindromic_N_proof.py`, `04-computation/palindromic_N_posuniform.py`, `04-computation/palindromic_N_n9.py`, `04-computation/palindromic_N_n11.py`, `04-computation/selfcomp_posuniform_n7.py`

### INV-074: Diagonal signed position theorem Ã¢â‚¬â€ VERIFIED n=5
**Source:** opus-2026-03-06-S11b (continuedÃ‚Â³)
**Status:** VERIFIED computationally at n=5 (all 12 iso classes).
**What:** M[v,v] = sum_P (-1)^{pos(v,P)} where the sum is over all Hamiltonian paths P of T and pos(v,P) is the 0-indexed position of v in P. This means M[v,v] counts even-position appearances minus odd-position appearances. M[v,v] can be NEGATIVE (not a path count). For VT tournaments at odd n: M[v,v] = H/n. "Defect vertex" = vertex whose position distribution is biased relative to average.
**Connection to THM-027:** The trace formula tr(M) = H(T) = sum_v sum_P (-1)^{pos(v,P)} = sum_P sum_v (-1)^{pos(v,P)} = sum_P 1 (at odd n, since alternating sum of 0..n-1 positions = 1). This reproves the trace formula.
**Next step:** Prove from IE formula definition. Connect to defect vertex characterization.
**Scripts:** `04-computation/diagonal_signed_position_theorem.py`

### INV-075: Perpendicularity of M-directions across H-classes Ã¢â‚¬â€ CONFIRMED n=7
**Source:** opus-2026-03-06-S11b (continuedÃ‚Â³)
**Status:** CONFIRMED computationally at n=7 (790 iso classes).
**What:** The non-scalar part of M (i.e., M - (H/n)*I) has a "direction" in matrix space. Measuring cosine similarity between these directions across iso classes shows:
  - Low H (H<85): positive cosine (~0.5-0.8, aligned)
  - Mid H (HÃ¢â€°Ë†85-105): near zero (perpendicular!)
  - High H (H>105): negative cosine (anti-aligned)
  - Overall mean cosine = -0.0485 (near perpendicular)
The crossover (true perpendicularity) occurs near the MEDIAN H value. This means the non-scalar perturbation "rotates" continuously through eigenspace as H varies.
**Connection:** This is the "perpendicularity" the user hypothesized in earlier sessions. It connects to the inverted-U of position variance and the grid symmetry structure.
**Next step:** (1) Prove analytically why crossover occurs at median H. (2) Check at n=9. (3) Connect to the even cycle vanishing theorem (INV-053) which uses the same TÃ¢â€ â€T^op involution.
**Scripts:** `04-computation/perpendicularity_cosine_n7.py`

### INV-076: All H-maximizers have scalar M Ã¢â‚¬â€ VERIFIED n=5,7
**Source:** opus-2026-03-06-S11b (continuedÃ‚Â³)
**Status:** VERIFIED exhaustively at n=5 (both max-H classes) and n=7 (all 43 maximizers).
**What:** Every tournament achieving max H has M = (H/n)*I (scalar transfer matrix). At n=7: all 43 maximizers have M = 27*I. At n=5: both H=15 classes (VT circulant and non-regular VT) have M = 3*I.
**Conjecture:** For ALL odd n, the H-maximizer has scalar M. This is equivalent to saying H-maximizers are always vertex-transitive (or at least "position-uniform").
**Connection:** Combines with THM-052 (circulant => scalar) and the Paley maximizer conjecture. If Paley maximizes H and Paley is circulant, then Paley gives scalar M. The deeper question is whether scalar M is NECESSARY for H-maximization.
**Next step:** Verify at n=9 (need to find maximizers first).

### INV-077: VT tournament NOT self-converse at n=21 Ã¢â‚¬â€ THM-052 DISPROVED for non-SC VT
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
  - w_{n-3} = 2*(n-2)!*t_3 - const (depends on t_3 only Ã¢â‚¬â€ PROVED)
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
**Speculative:** Is there a formal duality: paths (odd n) Ã¢â€ â€ cycle covers (even n)?
**Next step:** (1) Check if det(S) at n=6 is determined by (t_3, t_5). (2) Investigate the eigenvalue connection between Pf and H.

### INV-081: Paley tournament W(r) structure
**Source:** kind-pasteur-2026-03-06-S25f (eigenvalue_W_connection.py)
**Status:** COMPUTED at p=7.
**What:** Paley T_7 has W(r)/7! = [1/320, 0, 1/80, 0, 1/4, 0, 1].
  - Eigenvalues of A_Paley: (p-1)/2 = 3 (mult 1), (-1Ã‚Â±sqrt(p))/2 (mult (p-1)/2 each)
  - All non-trivial eigenvalues degenerate Ã¢â€ â€™ scalar M Ã¢â€ â€™ W coefficients have maximal symmetry
  - Paley requires p Ã¢â€°Â¡ 3 (mod 4) (so -1 is not a QR)
**Next step:** Compute W(r) for T_11 and T_3. Check if W/p! ratios have a closed form involving p.

### INV-082: EXACT W-coefficient hierarchy Ã¢â‚¬â€ spectral decomposition of tournaments
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
**Status:** OBSERVED at n=5. Eigenvalues {Ã‚Â±(1+Ã¢Ë†Å¡2), Ã‚Â±1, Ã‚Â±1, Ã‚Â±(Ã¢Ë†Å¡2-1)}.
**What:** The skeleton adjacency matrix at n=5 has eigenvalues involving the silver ratio 1+Ã¢Ë†Å¡2. K^2 diagonal = GS class sizes. Is this coincidence or does it persist at n=7?
**Next step:** Compute skeleton eigenvalues at n=7 (88Ãƒâ€”88 matrix). Check if silver ratio generalizes.

### INV-087: Antiferromagnetic interpretation of skeleton
**Source:** kind-pasteur-2026-03-06-S25h
**Status:** CONCEPTUAL. Skeleton = Ising model with antiferromagnetic coupling.
**What:** SC classes have "spin" = (-1)^{t3}. GS flip edges connect opposite spins. At odd n: perfect Neel order (unfrustrated). At even n: frustrated (odd cycles). Connects tournament theory to statistical mechanics.
**Next step:** Compute partition function Z(beta) = sum over SC of H(T)^beta. Check for phase transitions.

### INV-088: Schweser-Stiebitz-Toft Ã¢â‚¬â€ RÃƒÂ©dei revisited (expository)
**Source:** arXiv:2510.10659 (Oct 2025, revised Feb 2026), found by web search opus-2026-03-07-S36
**Status:** CATALOGED. Expository paper, likely low priority.
**What:** Revisits the classical theorems of RÃƒÂ©dei, Dirac, and Berge on Hamiltonian paths in tournaments. Exhibits the stronger theorems and explains connections between them. Does NOT mention independence polynomials, odd cycles, or conflict graphs.
**Next step:** Skim for any novel structural insight about H(T) parity not already in our framework. Low priority.

### INV-089: Irving-Omar authorship correction
**Source:** opus-2026-03-07-S37 (this session)
**Status:** CORRECTED in THM-002, CONJ-001, THM-070.
**What:** arXiv:2412.10572 ("Revisiting The RÃƒÂ©dei-Berge Symmetric Functions via Matrix Algebra") is by **Irving & Omar**, NOT Grinberg & Stanley. Their Corollary 20 restates Grinberg-Stanley's Theorem 1.39 + Lemma 6.5 from arXiv:2307.05569. The OCF result itself is correctly attributed to Grinberg-Stanley; only the paper authorship was wrong.
**Remaining:** Some computation scripts and broadcast messages still reference "Grinberg-Stanley" for arXiv:2412.10572. These are historical and low priority to fix.

### INV-090: Three equivalent H(T) formulas and the even-cycle cancellation
**Source:** opus-2026-03-07-S37
**Status:** VERIFIED computationally. Structural explanation needed.
**What:** Three equivalent ways to compute H(T):
1. **Direct**: count Hamiltonian paths (Held-Karp)
2. **OCF = I(ÃŽÂ©(T), 2)**: sum over T-only disjoint odd cycle collections, weight 2^k
3. **ps_1(U_T)(1)**: sum of ALL RÃƒÂ©dei-Berge p-coefficients (uses T + T^op cycles, signed by (-1)^Ãâ€ )

KEY FINDINGS at n=7 (Paley):
- All p-coefficients of U_T have ALL-ODD partitions (no even-part lambdas appear!)
- Sum of coefficients = 189 = H(T) [ps_1 at m=1]
- OCF specialization of FULL U_T = 433 Ã¢â€°Â  H(T)
- OCF specialization of T-ONLY part = 189 = H(T)
- The "OCF specialization" (p_1Ã¢â€ â€™1, p_oddÃ¢â€ â€™2, p_evenÃ¢â€ â€™0) is NOT how GS proves OCF
- GS uses ps_1(U_T)(1) which uses ALL cycles, not just T-direction ones

OPEN: Why are the p-coefficients supported only on all-odd partitions at n=7? Is this true for all tournaments? Omega constraint: (-1)^{n-l(ÃŽÂ»)} symmetry forces sum over even-l terms = 0, but here they're individually zero. Is the all-odd support a coincidence of n=7 (Mersenne prime) or universal?

**Mixed-direction findings at n=9:** Mixed T/T^op cycle pairs DO exist at n=9 (100+ per tournament). So U_T at n=9 WILL have even-part lambda contributions, but they must cancel in ps_1(1).
**RESOLVED (opus-S38 agent):** The all-odd support IS universal for tournaments (proved by Grinberg-Stanley Theorem 1.39). The factor 2^k in OCF is the **orientation multiplicity**: each odd cycle has exactly 2 directed orientations (T and T^op), both contributing sign +1 (since (-1)^{k-1}=+1 for odd k). An independent set of k cycles thus contributes 2^k copies. Even-part lambdas vanish because even cycles contribute opposite signs from T and T^op directions.
**Scripts:** `04-computation/omega_constraint.py`, `mixed_sum_n9.py`, `ut_specialization_n9.py`

### INV-091: H=21 permanent gap Ã¢â‚¬â€ PROVED through n=7
**Source:** opus-2026-03-07-S38, kind-pasteur-2026-03-07-S28
**Status:** PROVED through n=7 (exhaustive). Strong conjecture for all n.
**What:** H(T)=21 is never achieved by any tournament on n<=7 vertices. Exhaustive computation: 2,097,152 tournaments at n=7, zero with H=21. The gap 19Ã¢â€ â€™23 appears at both n=6 and n=7. OCF analysis: none of the valid (alpha_1,alpha_2) decompositions for H=21 are achievable.
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

### INV-095: Bags-of-sticks for OCF Ã¢â‚¬â€ DEAD END
**Source:** opus-2026-03-07-S38 agent
**Status:** CLOSED (dead end).
**What:** The bags-of-sticks decomposition (Mitrovic-Stojadinovic Theorem 4.2) expresses U_X via inclusion-exclusion over edge deletions. Under OCF specialization, every bag of sticks contributes 1 (acyclic digraphs have empty Omega). So the decomposition reduces to: H(T) = sum of inclusion-exclusion coefficients, which is trivially true. No new information for OCF.

### INV-096: H=21 Component Reduction (THM-079) Ã¢â‚¬â€ PROVED FOR ALL n
**Source:** opus-2026-03-07-S39 (partial), kind-pasteur-2026-03-07-S33 (completion)
**Status:** PROVED. H(T) Ã¢â€°Â  21 for ALL tournaments on ALL n vertices.
**What:** For H(T)=21, the OCF requires I(Omega(T),2)=21. Component factorization gives:
  - Disconnected case: IMPOSSIBLE. 21=3*7, but I(component)=7 impossible by THM-029 argument.
  - P_4 case: IMPOSSIBLE. P_4 realization blocked because sharing 3-cycles force extra cycles.
  - K_6-2e case: SUPERSEDED by Dichotomy Theorem proof.
**Dichotomy Theorem (Part R):** For cycle-rich T on nÃ¢â€°Â¥9, either (a) 3 disjoint 3-cycles exist (Ã¢Å¸Â¹ HÃ¢â€°Â¥27), or (b) safe deletion to cycle-rich TÃ¢Ë†â€™v exists (Ã¢Å¸Â¹ HÃ¢â€°Â¥H(TÃ¢Ë†â€™v)+2Ã¢â€°Â¥27). Proved via poisoning graph DAG argument. Combined with base case nÃ¢â€°Â¤8 (exhaustive) and Part J (non-cyclic vertex removal), gives H(T)Ã¢â€°Â 21 for ALL n.
**Key ingredients:** Lemma Q (cycle-rich Ã¢Å¸Â¹ no source/sink), poisoning graph has out-degree Ã¢â€°Â¤1 and is acyclic, DAG source deletion preserves cycle-rich.
**Scripts:** `04-computation/h21_gap_mechanism.py`, `h21_dichotomy_proof.py`, `h21_poisoning_graph.py`, `h21_cycle_rich_auto_no_ss.py`
**Writeup:** `03-artifacts/drafts/dichotomy-proof-formal.md`

### INV-097: u_T Size-Weighted Independence Polynomial (THM-078) Ã¢â‚¬â€ PROVED
**Source:** opus-2026-03-07-S39
**Status:** PROVED. u_T(m) = sum_j sw(j)*m^{n-2j} where sw(j) = sum over j-element independent sets of 2^|S|.
**What:** The size-weighted independence polynomial identity connects u_T(m) to the Omega(T) independence structure. Q_T(w) = u_T(sqrt(w))/sqrt(w) has all real non-positive roots for n<=8 (Leake-Ryder/Chudnovsky-Seymour stability for claw-free graphs). Fails at n>=9 when claws appear.
**Next step:** Check if Q_T root structure has implications for achievable H values.

### INV-098: Lichiardopol's Conjecture and Disjoint Cycle Forcing Ã¢â‚¬â€ EXPLORED
**Source:** kind-pasteur-2026-03-07-S33 (web research)
**Status:** EXPLORED. Used in H=21 proof context but not directly needed.
**What:** Lichiardopol's conjecture (proved for q=3 by Bang-Jensen, Bessy, ThomassÃƒÂ©): tournaments with min out-degree Ã¢â€°Â¥ (q-1)k-1 have k vertex-disjoint q-cycles. For 3-cycles with k=3: min outdeg Ã¢â€°Â¥ 5. However, at n=9 cycle-rich, min outdeg is ALWAYS Ã¢â€°Â¤ 4 (100% of 106,424 tested), so Lichiardopol doesn't directly fire. The poisoning graph argument covers ALL cases including those below the Lichiardopol threshold.
**Papers:** [Bang-Jensen-Bessy-ThomassÃƒÂ©](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v27i2p52)
**Next step:** Could be useful for other permanent H-gap proofs where kÃ¢â€°Â¥4 disjoint cycles are needed.

### INV-099: Chen-Chang 2024 Disjoint Cycles in Tournaments Ã¢â‚¬â€ CATALOGED
**Source:** kind-pasteur-2026-03-07-S33 (web research)
**Status:** CATALOGED. Not yet deeply investigated.
**What:** Chen-Chang (2024, J. Graph Theory) prove results on disjoint cycles in tournaments. Extends cycle-matching theory. Could provide tools for proving existence of kÃ¢â€°Â¥4 disjoint 3-cycles under weaker conditions than Lichiardopol.
**Paper:** [Chen-Chang 2024](https://onlinelibrary.wiley.com/doi/10.1002/jgt.23038)
**Next step:** Read the paper for potentially stronger theorems applicable to H-gap proofs.

### INV-100: Frankl's Proof of Erdos Matching Conjecture (k=3) Ã¢â‚¬â€ CATALOGED
**Source:** kind-pasteur-2026-03-07-S33 (web research)
**Status:** CATALOGED. Provides context for 3-uniform hypergraph matching.
**What:** Frankl proved the Erdos matching conjecture for k=3: bounds the max number of 3-element sets with no matching of size s+1. The cycle vertex sets in Omega(T) form a 3-uniform hypergraph (for 3-cycles). Frankl's bound could constrain the maximum number of 3-cycles with bounded matching number.
**Connection:** If mm(T) Ã¢â€°Â¤ 2, Frankl's bound limits |Omega_3(T)| Ã¢â€°Â¤ max(C(5,3), 3*3-3+1) depending on exact formulation. This could give an independent route to the dichotomy.
**Next step:** Check exact Frankl bound for our setting (n vertices, 3-uniform, matching Ã¢â€°Â¤ 2).

### INV-101: Other Permanent H-Gaps Beyond 7 and 21 Ã¢â‚¬â€ CONFIRMED THROUGH n=8 EXHAUSTIVE
**Source:** kind-pasteur-2026-03-07-S33, opus-2026-03-07-S43
**Status:** STRONG CONJECTURE that H=7 and H=21 are the ONLY permanent gaps.
**What:** With H=7 and H=21 both proved as permanent gaps (never achieved for ANY n), the natural question is: are there other permanent gaps?
**Computational evidence (opus-S43):**
  - ALL n=7 gaps (63, 107, 119, 149, 161-169, 173) fill at n=8 (sampling, very quickly)
  - n=8 exhaustive computation running (268M tournaments)
  - For HÃ¢â€°Â¥27 (wÃ¢â€°Â¥13): 20+ graph-feasible decompositions, blocking all seems impossible
  - Decomposition analysis: for wÃ¢â€°Â¥4, the (w-2,1,0,...) decomposition is available; it fails at w=10 due to cascade forcing (Part N) but works at all other w
**Algebraic argument (opus-S43):**
  - For wÃ¢â€°Â¥13: alpha_3Ã¢â€°Â¥1 decompositions become available (3 disjoint cycles feasible)
  - The number of decompositions grows rapidly with w (14 at w=10, 20 at w=13, 60 at w=20)
  - Each decomposition needs an independent tournament obstruction to block it
  - Only w=3 (1 feasible decomp) and w=10 (4 feasible decomps, all blocked) have this property
**Mod-4 result (Grinberg-Stanley Theorem 7.1):** H(T) Ã¢â€°Â¡ 1 + 2Ã‚Â·(# nontrivial odd cycles) mod 4. Does not directly rule out any odd H.
**Conjecture: H=7 and H=21 are the ONLY permanent gaps in the H-spectrum.**
**n=8 EXHAUSTIVE RESULT (opus-S45):** All 268,435,456 tournaments enumerated. Max H=661. Only missing odd values in [1,300]: H=7 and H=21. This CONFIRMS the conjecture through n=8. No new gaps appear. All n=7 gaps (63, 107, 119, etc.) fill at n=8.
**Status:** STRONG EVIDENCE. Conjecture holds through n=8 exhaustive enumeration.

### INV-102: Grinberg-Stanley Mod-4 Theorem (Theorem 7.1) Ã¢â‚¬â€ CATALOGED
**Source:** opus-2026-03-07-S43 (web research)
**Status:** CATALOGED. Read and extracted from arXiv:2307.05569.
**What:** Theorem 7.1: H(T) Ã¢â€°Â¡ 1 + 2Ã‚Â·(# nontrivial odd D-cycles) (mod 4). This is the OCF mod-4 reduction: since H = 1 + 2Ã‚Â·alpha_1 + 4Ã‚Â·(...), H mod 4 = 1 + 2Ã‚Â·alpha_1 mod 4. The proof uses the power-sum expansion and the specialization map zeta. Not directly useful for gaps but confirms the algebraic structure.
**Next step:** Check if higher modular refinements (mod 8, mod 16) exist in the Grinberg-Stanley framework.

### INV-103: Non-Separating Vertices in Tournaments Ã¢â‚¬â€ CATALOGED
**Source:** opus-2026-03-07-S43 (web research)
**Status:** CATALOGED. Related but weaker than cycle-rich deletion.
**What:** A vertex in a strongly connected tournament is "non-separating" if its removal preserves strong connectivity. For min in/out-degree Ã¢â€°Â¥ p, at least min{|V|, 4p-2} non-separating vertices exist. Our "good deletion" requirement is stronger: preserve cycle-richness (every vertex in a 3-cycle), not just strong connectivity.
**Next step:** Could the non-separating vertex techniques be adapted to our stronger requirement?

### INV-104: "Cycle-Rich" as Novel Concept Ã¢â‚¬â€ NOTED
**Source:** opus-2026-03-07-S43 (web research)
**Status:** The term "cycle-rich" (every vertex in a directed 3-cycle, no source/sink) does not appear in the literature. This is a novel concept from our project. The poisoning graph argument (Part R) may be publishable as a standalone result about cycle-rich tournaments.

### INV-105: Deletion-Contraction for H(T) Ã¢â‚¬â€ VERIFIED
**Source:** kind-pasteur-2026-03-07-S34
**Status:** VERIFIED COMPUTATIONALLY. H(T) = H(T\e) + H(T/e) holds 100% at n=4,5. Commutative specialization of Mitrovic's W_X = W_{X\e} - W_{X/e}Ã¢â€ â€˜ (arXiv:2504.20968).
**Convention:** Contraction merges tail/head: w inherits IN from tail, OUT from head.
**Next step:** Prove algebraically (should follow from Mitrovic by specialization). Use for inductive proof of Redei/OCF/H-gaps.

### INV-106: GLMY Path Homology of Tournaments Ã¢â‚¬â€ IDENTIFIED
**Source:** kind-pasteur-2026-03-07-S34 (web research)
**Status:** IDENTIFIED. Tang-Yau (2026, arXiv:2602.04140) compute path homology of circulant digraphs using Fourier decomposition Ã¢â‚¬â€ directly applicable to Paley tournaments.
**Key connection:** H_1(T) should relate to cycle space of T and thus to Omega(T) and alpha_1.
**Next step:** Compute GLMY path homology for small tournaments (n=4,5). Check if rank(H_1) = alpha_1 or c_3.

### INV-107: Extended Root Polytope Deletion-Contraction Ã¢â‚¬â€ IDENTIFIED
**Source:** kind-pasteur-2026-03-07-S34 (web research). Kalman-Tothmeresz 2024, arXiv:2409.18902.
**Status:** IDENTIFIED. h*-polynomial of extended root polytope is monotone under deletion-contraction.
**Next step:** Compute root polytope for small tournament digraphs. Check if h* relates to I(Omega,x).

### INV-108: Lee-Yang Zeros of F(T,x) Ã¢â‚¬â€ DISCOVERED (opus-S44)
**Source:** opus-2026-03-07-S44
**Status:** MAJOR DISCOVERY. F(T,x) zeros come in reciprocal pairs (palindrome). Cluster at Ã‚Â±2pi/3 on unit circle. H=9 at n=5: ALL zeros on unit circle.
**Key findings:** F(T,omega) real at n=7, F(T,i) pure imaginary, universal divisibilities F(T,omega) = 0 mod 9 and F(T,i) = 0 mod 16i.
**Next step:** Prove Lee-Yang property for specific tournament classes. Connect to phase transitions in statistical mechanics.

### INV-109: Walsh/Fourier Spectral OCF Ã¢â‚¬â€ DISCOVERED (opus-S35c)
**Source:** opus-2026-03-07-S35c9
**Status:** MAJOR DISCOVERY. THM-081: hat{t_k}[S] = (1/2^k) sum (-1)^{asc(S,C)}. Counting identity provides new proof path for OCF.
**Next step:** Prove counting identity algebraically for d=1 (single edge). Extend to n=7 (need hat{t7} and hat{bc35}).

### INV-110: Ihara Zeta Function of Tournaments Ã¢â‚¬â€ TESTED
**Source:** kind-pasteur-2026-03-07-S34
**Status:** TESTED. z_inv(1/2) = det(I-A/2+(D-I)/4) strongly correlated with H (r=-0.95 at n=5) but NOT uniquely determined.
**Conclusion:** Ihara zeta constrains H but doesn't determine it. Consistent with cycles constraining but not determining independence structure.

### INV-111: p-adic Structure Beyond p=2 Ã¢â‚¬â€ TESTED
**Source:** kind-pasteur-2026-03-07-S34
**Status:** TESTED. H mod 3 = (1+2*alpha_1+alpha_2) mod 3 from OCF. At n=4: H mod 3 uniquely determined by c3 via (1+2c3) mod 3. H mod 7 = 0 impossible at n<=6, first achievable at n=7.
**Next step:** Investigate H mod p for larger p. Is there a p-adic tower for p=3?

### INV-112: Converse Invariant Digraph Polynomials Ã¢â‚¬â€ IDENTIFIED
**Source:** kind-pasteur-2026-03-07-S34 (web research). Ai-Gutin-Lei-Yeo-Zhou 2024, arXiv:2407.17051.
**Status:** IDENTIFIED. New digraph polynomial for converse invariance testing. H(T)=H(T^op) means Ham paths are converse invariant.
**Next step:** Read the paper. Check if their polynomial gives new information about H(T).

### INV-113: Stanley-Stembridge Resolution Implications Ã¢â‚¬â€ IDENTIFIED
**Source:** kind-pasteur-2026-03-07-S34 (web research). Hikita 2024.
**Status:** IDENTIFIED. Stanley-Stembridge (e-positivity of chromatic SF for 3+1-free posets) proved by Hikita 2024.
**Connection:** Via Mitrovic-Stojadinovic, Redei-Berge U_T connects to chromatic SF. If tournament poset is (3+1)-free, U_T inherits e-positivity.
**Next step:** Check if tournament arc ordering posets are (3+1)-free. Investigate Hessenberg varieties approach.

### INV-114: Flip Formula F(T,x) - F(T',x) = (x-1)*D(x) Ã¢â‚¬â€ PROVED (THM-083)
**Source:** opus-2026-03-07-S45 (computational discovery), kind-pasteur-2026-03-07-S35 (algebraic proof)
**Status:** PROVED algebraically (THM-083). Verified at n=4,5 exhaustive.
**What:** For arc u->v in T, flip gives T'. The difference F(T)-F(T') factors as (x-1)*D(x), where D = F(T/e) - F(T'/e') is anti-palindromic.
**Key identities (THM-083):**
  - F_T(x) = F_{T\e}(x) + (x-1) * F(T/e, x)  (polynomial deletion-contraction)
  - G_{u,v}(x) = F(T/e, x)  (contraction = conditional path polynomial)
  - D(x) = -x^{n-2} D(1/x)  (anti-palindromicity from tournament palindrome)
  - H(T) - H(T') = D_{n-2}  (leading coefficient of D)
**CORRECTION (kind-pasteur-S35):** H(T) Ã¢â€°Â  H(T') under arc flip in general (deltas up to Ã‚Â±12 at n=5). The opus claim "H(T)=H(T')" was WRONG. The correct statement: F(T,1)=n!=F(T',1) (total permutation count).
**CORRECTION (kind-pasteur-S35):** G_uv + G_vu = 2*F(T/e) only when T/e is a tournament (requires u,v to have identical profiles to other vertices). In general, G_uv + G_vu = F(T/e) + F(T'/e') which is palindromic but Ã¢â€°Â  2*F(T/e).
**Scripts:** `04-computation/f_poly_flip_formula.py`, `04-computation/flip_formula_D_analysis.py`, `04-computation/poly_deletion_contraction.py`, `04-computation/flip_reduction_via_contraction.py`

### INV-115: Matroid Structure of Vertex-Disjoint Odd Cycles Ã¢â‚¬â€ BOUNDARY at n=5
**Source:** opus-2026-03-07-S45 (computational discovery)
**Status:** VERIFIED. Exchange axiom holds at n=5 (1024/1024) but FAILS at n>=6 (15360/32768 at n=6).
**What:** Collections of vertex-disjoint odd directed cycles in T form the independent sets of a matroid if and only if n<=5. At n>=6, maximal independent sets can have different sizes.
**Script:** `04-computation/gammoid_matroid_test.py`
**Next step:** (1) Prove n<=5 case (small enough for case analysis). (2) Characterize failure at n=6 Ã¢â‚¬â€ which exchange pairs fail? (3) Relationship to Omega(T) perfectness boundary (perfect through n=7, fails n=8).

### INV-116: Transfer Matrix W(x) and per(W) Ã¢â‚¬â€ EXPLORED
**Source:** opus-2026-03-07-S45
**Status:** EXPLORED. per(W(1)) = D_n (subfactorial) universally. F(T,x) = Hamiltonian path sum over W entries. per(W(x)) palindromic for certain tournament classes.
**Scripts:** `04-computation/transfer_matrix_F_connection.py`, `04-computation/per_W_analysis.py`
**Next step:** (1) Explore eigenvalues of W(x) at specific x values. (2) Connection to Irving-Omar det formula (INV-046). (3) Can det(I-zW) generating function extract F?

### INV-117: Archer-Gessel-Graves-Liang Strong Tournament Descent Polynomial Ã¢â‚¬â€ RESEARCHED
**Source:** opus-2026-03-07-S45 (background agent), Discrete Math 343 (2020)
**Status:** RESEARCHED. Paper fully reviewed. t_n(u) = descent poly for strong tournaments. Palindromic, divisible by (1+u)^{floor(n/2)}. GF: U(x) = 1/(1-T(x)).
**Connection:** Different statistic from F(T,x) (global descent vs path-local forward edges), but same palindromic structure. The (1+u)^{floor(n/2)} divisibility may connect to OCF factor-2 structures. The "Eulerian graphic GF" framework (q=(1+uy)/(1+y)) is the natural algebraic home for descent statistics on tournaments.
**Notes:** `04-computation/gessel_strong_tournament_notes.md`
**Next step:** (1) Compute t_n(u) at small n and compare to F(T,x) aggregates. (2) Check if the strong component decomposition gives new structural insights for H(T). (3) Test the (1+u)^{floor(n/2)} divisibility analogue for F(T,x).

### INV-118: F(T,omega) mod 9 Universality Ã¢â‚¬â€ CONFIRMED NOVEL
**Source:** opus-2026-03-07-S44/S45 (computational discovery + background agent literature search)
**Status:** CONFIRMED NOVEL. Extensive literature search found NO prior work on roots-of-unity evaluations of F(T,x). Not a consequence of Grinberg-Stanley mod-4. Chebikin et al. studied cyclotomic factors of descent set polynomial Q_n but Phi_3 doesn't appear for n<=23.
**What:** F(T,omega) Ã¢â€°Â¡ 0 mod 9 at n=7 (all 5040 tournaments). Equivalently S_0 = sum_{kÃ¢â€°Â¡0 mod 3} F_k Ã¢â€°Â¡ 0 mod 6.
**Next step:** (1) Prove algebraically using OCF or Fourier decomposition. (2) Check at n=9,10. (3) Generalize: are there universal congruences for F(T,zeta_k) mod k^2?

### INV-119: Deletion-Contraction for Hamiltonian Paths Ã¢â‚¬â€ PROVED (THM-082)
**Source:** kind-pasteur-2026-03-07-S35
**Status:** PROVED by clean bijection argument. Verified exhaustive n=4,5.
**What:** For any digraph D with directed edge e=(uÃ¢â€ â€™v):
  H(D) = H(D\e) + H(D/e)
where D\e = deletion, D/e = contraction (w inherits IN from u, OUT from v).
**Proof:** Ham paths not using e = H(D\e). Ham paths using e biject with Ham paths of D/e via collapsing ...Ã¢â€ â€™uÃ¢â€ â€™vÃ¢â€ â€™... to ...Ã¢â€ â€™wÃ¢â€ â€™...
**Corollary:** Arc-flip H-difference reduces to contraction: H(T)-H(T') = H(T/e) - H(T'/e').
**Key structural insight:** T/e and T'/e' differ ONLY in how w connects to other vertices (profile swap). If u,v have identical profiles, T/e = T'/e' and H(T) = H(T').
**Connection to Mitrovic:** Commutative specialization of W_X = W_{X\e} - W_{X/e}Ã¢â€ â€˜ (arXiv:2504.20968).
**Scripts:** `04-computation/deletion_contraction_test.py`, `04-computation/flip_reduction_via_contraction.py`

### INV-120: Polynomial Deletion-Contraction for F(T,x) Ã¢â‚¬â€ PROVED (THM-083)
**Source:** kind-pasteur-2026-03-07-S35
**Status:** PROVED algebraically. Verified exhaustive n=4,5.
**What:** F_T(x) = F_{T\e}(x) + (x-1) * F(T/e, x). Generalizes THM-082 to polynomial level.
**Key identification:** G_{u,v}(x) = F(T/e, x) Ã¢â‚¬â€ the "conditional path polynomial" summing over permutations with u immediately before v equals the forward-edge polynomial of the contraction.
**Flip formula as corollary:** F_T - F_{T'} = (x-1) * [F(T/e) - F(T'/e')], with D anti-palindromic.
**Anti-palindromicity proof:** D(x) = -x^{n-2}D(1/x) follows from palindromicity of F_T, F_{T'}.
**Scripts:** `04-computation/poly_deletion_contraction.py`

### INV-121: F(T,omega) mod 9 universality Ã¢â‚¬â€ PROVED (THM-085)
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

### INV-123: Worpitzky Expansion of F(T,x) Ã¢â‚¬â€ PROVED (THM-084)
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

### INV-124: Signed Forward-Edge Polynomial SF(T,x) Ã¢â‚¬â€ PROVED (THM-085b)
**Source:** opus-2026-03-07-S46b
**Status:** PROVED algebraically
**What:** SF(T,x) = sum sgn(sigma) x^{fwd_T(sigma)} is palindromic with parity (-1)^{C(n,2)}.
  - SF(T,1) = 0 always, so (x-1) | SF(T,x)
  - Quotient SF/(x-1) is anti-palindromic
  - At n=4: SF = c(T) * (x-1)^2(x+1) (since anti-palindromic of even degree has (x+1) factor)
  - SF determines F at n<=5 but NOT at n>=6 (coarser invariant)
**Connection:** SF is a "path immanant" for the sign character. F is the "path permanent."
**Scripts:** `04-computation/signed_F_analysis.py`

### INV-125: Forward-Edge Variance Formula Ã¢â‚¬â€ PROVED (THM-086)
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

### INV-127: GLMY Path Homology of Tournaments Ã¢â‚¬â€ EVEN BETTI VANISHING
**Source:** opus-2026-03-07-S46e (path_homology_phase2.py output + new analysis)
**Status:** PARTIALLY CONFIRMED Ã¢â‚¬â€ ÃŽÂ²Ã¢â€šâ€š=0 exhaustive n<=6, sampled n<=9 (0 failures in ~50k tests). ÃŽÂ²Ã¢â€šâ€ž=0 exhaustive n<=6, sampled n<=7 (0/3000 random). BUT Paley T_7 has ÃŽÂ²Ã¢â€šâ€ž=6! And ÃŽÂ²Ã¢â€šâ€ž=1 found in 0.6% of random n=8 tournaments. So only ÃŽÂ²Ã¢â€šâ€š=0 appears truly universal.
**What:** ÃŽÂ²Ã¢â€šâ€š(T) = 0 for ALL tournaments T (HYP-207). ÃŽÂ²Ã¢â€šâ€ž can be nonzero starting at n=7 (Paley) and n=8 (random). ÃŽÂ²Ã¢â€šÂ and ÃŽÂ²Ã¢â€šÆ’ NOT mutually exclusive at n=8 (need to check). Ãâ€¡(ÃŽÂ©) Ã¢Ë†Ë† {-11,...,7} at n=7, NOT {0,1} Ã¢â‚¬â€ HYP-267 REFUTED.
**Evidence (S42):** Exhaustive n=3,4,5,6; sampled n=7 (5000+), n=8 (1000+), n=9 (100). ÃŽÂ²Ã¢â€šâ€š=0 in ALL cases.
**Only 3 Betti profiles at n=5,6:** (1,0,...), (1,1,0,...), (1,0,...,ÃŽÂ²Ã¢â€šÆ’=1,...). At n=7: same 3 + Paley's (1,0,0,0,6,0). At n=8: adds (1,0,0,0,1,0,0).
**Key formulas (proved n=5):** dim(AÃ¢â€šâ€š) = C(n,3) + 2cÃ¢â€šÆ’; dim(ÃŽÂ©Ã¢â€šâ€š) = dim(AÃ¢â€šâ€š) - #{non-allowed pairs with mediators}; rank(dÃ¢â€šâ€š) = C(n-1,2) - ÃŽÂ²Ã¢â€šÂ.
**Algebraic mechanism:** ker(dÃ¢â€šâ€š|ÃŽÂ©Ã¢â€šâ€š) = im(dÃ¢â€šÆ’|ÃŽÂ©Ã¢â€šÆ’) always. The "swap cycle" characterization: pure v-chain 2-cycles have form ÃŽÂ£ M_{ab}[(a,b,v)-(v,a,b)] with zero row/col sums. ALL swap cycles are boundaries (confirmed exhaustive n=5,6).
**Literature (S42):** Tang-Yau (2026): H_m=0 for m>=2 for circulant tournaments. Burfitt-Cutler: ÃŽÂ©Ã¢â€šâ€š generated by transitive triples only. No paper addresses ÃŽÂ²Ã¢â€šâ€š=0 for general tournaments Ã¢â‚¬â€ this is genuinely open.
**Not in literature:** Confirmed via comprehensive search of Caputi-Menara, Burfitt-Cutler, Fu-Ivanov, Tang-Yau, Chaplin, all GLMY papers.
**Scripts:** `beta_parity_pattern.py`, `beta2_algebraic_mechanism.py`, `beta2_deformation_retract.py`, `beta2_large_n_sample.py`, `beta_paley_verify.py`, and many more
**Cone-from-T' construction (S42):** For vertex v with swap cycle z, the filling w = ÃŽÂ£ ÃŽÂ±_{abc} [(v,a,b,c)+(a,b,c,v)] over T'=T\{v} 2-paths always works. The T'-internal faces cancel: dÃ¢â€šÆ’(v,a,b,c)+dÃ¢â€šÆ’(a,b,c,v) has zero (a,b,c) component. The resulting BÃ‚Â·ÃŽÂ±=z system is always solvable.
- **Filtered** (only T' paths with vÃ¢â€ â€™a AND cÃ¢â€ â€™v): Works exhaustive n=5,6 (32768/32768). Fails 1/1000 at n=8.
- **Unfiltered** (ALL T' paths): Works 500/500 at n=7, 500/500 at n=8, 200/200 at n=9. Zero failures.
- **Multi-vertex** (combine all vertices): ALWAYS works, including the n=8 filtered failure.
- **ÃŽÂ©Ã¢â€šÆ’ auto-membership**: Cone filling automatically in ÃŽÂ©Ã¢â€šÆ’ at n=5,6 (100%). Breaks at nÃ¢â€°Â¥7 (~98% at n=7, ~93% at n=8).
- **Rank surplus grows**: rank(B)-swap_dim min is 2(n=5), 2(n=6), 6(n=7), 11(n=8), 15(n=9). System increasingly overdetermined.
- **ÃŽÂ²Ã¢â€šâ€š=0 confirmed through n=10**: Exhaustive nÃ¢â€°Â¤6, sampled n=7(500), n=8(1000), n=9(200), n=10(50). Zero failures.
- **Dimension formula**: rank(dÃ¢â€šÆ’) = ker(dÃ¢â€šâ€š) EXACTLY for every tournament tested (n=5-9).
**Next step:** (1) Prove BÃ‚Â·ÃŽÂ±=z always solvable algebraically (rank argument). (2) Prove ÃŽÂ©Ã¢â€šÆ’ membership of filling. (3) Try inductive proof using LES of pair (T, T\v). (4) Check if result follows from multisquare-free property (Fu-Ivanov).

### INV-128: Universal Coefficient Theorem Ã¢â‚¬â€ PROVED (THM-117)
**Source:** opus-2026-03-07-S46e
**Status:** PROVED
**What:** coeff(t_{2k+1} in ÃŽÂº_{2k}) = 2/C(n, 2k) for all k. Proved via forward path formula + OCF + multinomial expansion. Resolves OPEN-Q-023.
**Scripts:** `universal_coeff_proof.py`

### INV-129: Celano-Sieger-Spiro A_T(t) Ã¢â‚¬â€ NOT same as F(T,x)
**Source:** opus-2026-03-07-S46e (web research)
**Status:** CLARIFIED (dead end for direct application)
**What:** arXiv:2309.07240 defines A_T(t) = sum over labelings of t^{des_T(sigma)} where des counts descents across ALL arcs. This has degree C(n,2), not n-1. The (1+t)^{floor(n/2)} divisibility applies to A_T(t), not to our F(T,x). The two polynomials encode different statistics (all-arc descents vs Hamiltonian-path forward edges).
**Impact:** The Celano-Sieger-Spiro result cannot be directly applied to F(T,x). However, it establishes that tournaments have a universal structural constraint on A_T(t) depending only on n, which is analogous to our universal constraint on F(T,x) mod 3 (THM-086).

### INV-130: Pfaffian-Betti Connection Ã¢â‚¬â€ EXHAUSTIVE at n=6, extended n=7,8
**Source:** opus-2026-03-07-S46e, kind-pasteur-2026-03-08-S40
**Status:** VERIFIED EXHAUSTIVE n=6. Sampled n=7,8. THM-120 (was THM-098) + THM-099 documented.
**What:** The Pfaffian of the skew-adjacency matrix constrains path homology Betti numbers. At n=6 (exhaustive): ÃŽÂ²Ã¢â€šÂ>0 Ã¢Å¸Â¹ |Pf(S)| Ã¢Ë†Ë† {1,3}; ÃŽÂ²Ã¢â€šÆ’>0 Ã¢Å¸Â¹ |Pf(S)| Ã¢Ë†Ë† {7,9}. Perfect separation. At n=7 (odd): spectral gap separates phases. At n=8: |Pf| NOT perfect separator but strongly correlated.
**CORRECTED (S40):** H-maximizers at n=6 are NOT all S-phase. 480 maximizers split 240 C-phase (|Pf|=1) + 240 S-phase (|Pf|=7), both with score (2,2,2,3,3,3) and c3=8. Complementation preserves phase.
**Scripts:** `pfaffian_betti_check.py`, `pfaffian_betti_n7.py`, `pfaffian_topology_deep.py`, `pfaffian_betti_mechanism.py`, `spectral_betti_gap.py`, `spectral_topology_n8.py`, `s_phase_structure.py`, `s_phase_maximizer_n7.py`, `maximizer_betti_deep.py`
**Next step:** (1) Prove Pfaffian separation algebraically at n=6. (2) Why does spectral gap separate at n=7?

### INV-135: H-Maximizer Betti Dimension Shift Ã¢â‚¬â€ THM-099
**Source:** kind-pasteur-2026-03-08-S40
**Status:** VERIFIED EXHAUSTIVE n=4,5,6. Sampled n=7.
**What:** H-maximizers always have nontrivial GLMY path homology, with the topological dimension increasing:
- n=4: ALL 24 max have ÃŽÂ²Ã¢â€šÂ=1 (C-phase)
- n=5: ALL 64 max have ÃŽÂ²Ã¢â€šÂ=1 (C-phase)
- n=6: 480 max split 240 ÃŽÂ²Ã¢â€šÂ=1 + 240 ÃŽÂ²Ã¢â€šÆ’=1
- n=7: ALL 240 max have ÃŽÂ²Ã¢â€šâ€ž=6 (beyond S-phase classification)
At n=7, all maximizers are conference-matrix (gap=0, eigenvalues all Ã¢Ë†Å¡7). Second-highest H=175 has ÃŽÂ²Ã¢â€šÂ=1 (C-phase). Third H=171 is contractible. Topology stratifies H values.
**Scripts:** `betti_dimension_shift.py`, `betti_dimension_shift_v2.py`, `maximizer_betti_deep.py`, `maximizer_betti_n8.py`
**Next step:** (1) Check n=8 maximizers (H=661). (2) Why ÃŽÂ²Ã¢â€šâ€ž=6 specifically at n=7? (3) Algebraic mechanism connecting eigenvalue uniformity to high Betti.

### INV-131: Path Homology Hidden Invariant Ã¢â‚¬â€ Cycle Overlap Pattern
**Source:** opus-2026-03-07-S46e
**Status:** VERIFIED EXHAUSTIVE n=5
**What:** Path homology ÃŽÂ²Ã¢â€šÂ is NOT determined by (F-polynomial, tÃ¢â€šÆ’, score sequence, strong connectivity). At n=5, the distinguishing invariant is the 3-cycle overlap pattern: ÃŽÂ²Ã¢â€šÂ=1 iff all 3-cycles form a "star" (share a common edge); ÃŽÂ²Ã¢â€šÂ=0 when cycle overlaps are heterogeneous. Two non-isomorphic tournament types with identical F-poly and tÃ¢â€šÆ’ have different ÃŽÂ²Ã¢â€šÂ.
**Scripts:** `betti_hidden_invariant.py`, `f_poly_betti_deep.py`
**Next step:** (1) Formalize the overlap pattern criterion. (2) Check at n=6,7.

### INV-132: Mod-2 Cumulant Collapse via OCF
**Source:** opus-2026-03-07-S46e
**Status:** VERIFIED
**What:** THM-094 (F mod 2 universal) implies all integer moment sums n!Ã‚Â·E[fwd^r] Ã¢â€°Â¡ 0 mod 2 for nÃ¢â€°Â¥3. Root cause: OCF evaluates at x=2 Ã¢â€°Â¡ 0 mod 2, so I(ÃŽÂ©,2) mod 2 = I(ÃŽÂ©,0) = 1, giving H Ã¢â€°Â¡ 1 = RÃƒÂ©dei. Combined with THM-086 (mod 3 zeros), all moments vanish mod 6 for nÃ¢â€°Â¥3.
**Scripts:** `mod2_cumulant_connection.py`

### INV-133: Tang-Yau Circulant Path Homology (arXiv:2602.04140)
**Source:** opus-2026-03-07-S46e (web research)
**Status:** NEW LEAD
**What:** Tang-Yau (Feb 2026) compute GLMY Betti numbers of circulant digraphs using Fourier decomposition via shift automorphism. Additive structure of connection set S governs low-dimensional chains. Could be adapted to circulant tournaments (odd n, S={1,...,(n-1)/2}).
**Next step:** Read paper. Check if their Fourier method can prove ÃŽÂ²Ã¢â€šâ€š=0 for tournament subclass.

### INV-134: GLMY Torsion-Zeta Bridge
**Source:** opus-2026-03-07-S46e (web research)
**Status:** NEW LEAD
**What:** Grigoryan-Muranov define Reidemeister/analytic torsion on digraphs via path homology Hodge Laplacian, proving Cheeger-MÃƒÂ¼ller analogue. Torsion invariants are classically expressible through zeta-type regularized determinants. This could connect the Ihara/Bartholdi zeta of a tournament to its path homology.
**References:** arXiv:2012.07302, grigor/torsion.pdf
**Next step:** Investigate whether Bartholdi zeta of a tournament encodes path homology info.

### INV-135: Asao Ã¢â‚¬â€ Transitive Tournaments Have Trivial Path Homology
**Source:** opus-2026-03-07-S46e (web research, arXiv:2503.06722)
**Status:** CONFIRMED (known)
**What:** Asao (March 2025) proves transitive tournaments are "regularly diagonal" digraphs with vanishing reduced regular path homology. Confirms our computational finding that transitive tournaments have ÃŽÂ²=(1,0,...,0).
**Impact:** Provides a published reference for one endpoint of the topological trichotomy.

### INV-136: Chaplin Ã¢â‚¬â€ Random Digraph ÃŽÂ²Ã¢â€šÂ Phase Transitions
**Source:** opus-2026-03-07-S46e (web research, arXiv:2111.13493)
**Status:** NEW LEAD
**What:** Chaplin (2022) shows ÃŽÂ²Ã¢â€šÂ of random ErdÃ…â€˜s-RÃƒÂ©nyi digraphs has two phase transitions. Since tournaments are "density 1/2" digraphs, this places them in a specific regime. Could explain why ~30% of tournaments at n=5 have ÃŽÂ²Ã¢â€šÂ>0.
**Next step:** Check if their density threshold matches tournament ÃŽÂ²Ã¢â€šÂ fraction.

### INV-137: THM-118 Trace-Cycle Identity Ã¢â‚¬â€ PROVED (extended to k=3,4,5)
**Source:** kind-pasteur-2026-03-07-S39b
**Status:** PROVED
**What:** tr(A^k) = k * c_k(T) for k=3,4,5 in any tournament. Extended to k=4 (no bidirectional edges => length-4 closed walks must be simple 4-cycles). Gives O(n^3) c_4 and c_5 computation via matrix multiplication. Sharp: fails at k>=6 (compound (3,3) walks at k=6). Correction for k=6 is NOT a simple polynomial in global cycle counts (tested and failed).
**Impact:** c4_fast() and c5_fast() in tournament_fast.py. Speedups: 3.8x for c4, 5.4x for c5 at n=8.
**Scripts:** `trace_cycle_k4.py`, `c6_correction_formula.py`, `c6_from_trace.py`

### INV-138: BjÃƒÂ¶rklund Cycle Cover Ã¢â‚¬â€ OCF Connection Tested (NEGATIVE)
**Source:** kind-pasteur-2026-03-07-S39b
**Status:** TESTED (NEGATIVE for new identities)
**What:** Tested 6 BjÃƒÂ¶rklund-style cycle cover formulations for OCF connections. Only Test 2 (partial odd cycle cover weighted by 2^{num_cycles}) matches OCF Ã¢â‚¬â€ but this IS OCF restated. No new route to proving OCF found. Permanent of A+I counts cycle covers but doesn't simplify OCF.
**Scripts:** `bjorklund_cycle_cover.py`

### INV-139: h-Positivity of U_T Ã¢â‚¬â€ CLOSED (fails for all non-transitive)
**Source:** kind-pasteur-2026-03-07-S39b
**Status:** CLOSED (dead end)
**What:** U_T is NOT h-positive for any non-transitive tournament. Only the transitive tournament (H=1) has h-positive U_T. This closes the Stanley-Stembridge connection for tournament RÃƒÂ©dei-Berge functions. The e-positivity question from INV-051/052 is also resolved negatively.
**Scripts:** `bjorklund_cycle_cover.py` (h-positivity test section)

### INV-140: THM-097 Alpha_2 Trace Formula Ã¢â‚¬â€ PROVED
**Source:** kind-pasteur-2026-03-07-S39b
**Status:** PROVED. O(n^3) computation of vertex-disjoint 3-cycle pairs.
**What:** alpha_2(Omega_3) = C(c3,2) - sum_v C(t3(v),2) + s2, where t3(v) = (A^3)[v][v] and s2 = sum_{edges a->b} C((A^2)[b][a], 2). Proof via inclusion-exclusion on pair overlap counts. Valid for full Omega at n<=7 (since 5+3=8>7 prevents cross terms). Implemented as alpha2_from_trace() in tournament_fast.py.
**Scripts:** `trace_ocf_bridge.py`, `alpha2_formula.py`

### INV-141: H(T) Polynomial Trace Formula Ã¢â‚¬â€ VERIFIED n<=9
**Source:** kind-pasteur-2026-03-07-S39b
**Status:** VERIFIED. 100% match at n=5,6 (exhaustive), n=7 (500), n=8 (100), n=9 (200).
**What:** H(T) = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3, all computable from matrix trace data. At n<=7: O(n^3). At n=8: O(n^5) (cross terms). At n=9: O(n^5) with additional O(2^7*C(n,7)) for c7 and O(C(n,3)^3/6) for alpha_3. Timing: trace method 7x slower than DP at n=9 but POLYNOMIAL.
**Key n=9 findings:** alpha_3 nonzero 86% of time (3+3+3=9). alpha_2(3,5) cross dominates alpha_2(3,3) by 2:1. H contribution: 56% alpha_1, 41% alpha_2, 2.3% alpha_3.
**Scripts:** `h_from_trace_n8.py`, `h_polynomial_n9.py`, `alpha_structure_n9.py`

### INV-142: Spectral Characterization of H-Maximizers Ã¢â‚¬â€ NEW FINDING
**Source:** kind-pasteur-2026-03-07-S39b
**Status:** COMPUTED. Key structural finding.
**What:** Paley T_p has the CONFERENCE MATRIX property: S^2 = -pI + J where S = A-A^T. This means ALL nonzero skew eigenvalues equal Ã‚Â±i*sqrt(p) (zero spectral gap). This property CHARACTERIZES Paley among DRTs at n=11 (non-Paley would fail). At n=7: Paley has zero spectral gap while other regular tournaments have gap 3.46-3.90. The general spectral correlation (|skew_max| vs H) is weak (-0.03) but among REGULAR tournaments, zero spectral gap Ã¢â€ â€™ max H.
**Also found:** tr(S^k) = 0 for ALL odd k in ALL tournaments (skew-symmetry). The adjacency spectrum does NOT distinguish DRTs (all have same eigenvalues), only the skew spectrum does. Paley's conference matrix property is a VERY strong constraint.
**Scripts:** `spectral_h_maximizer.py`, `spectral_cycle_density.py`

### INV-143: MISTAKE-017 Ã¢â‚¬â€ Invalid DRT at n=11 Corrected
**Source:** kind-pasteur-2026-03-07-S39b
**Status:** CORRECTED
**What:** The "non-Paley DRT" from {1,2,3,5,8} was NOT a valid tournament (SÃ¢Ë†Â©(-S)={3,8}Ã¢â€°Â Ã¢Ë†â€¦). All claims about c3=44, c5=407, H=69311 were computed on a non-tournament digraph. The only valid circulant DRT at n=11 is Paley. Exhaustive search found exactly 2 valid (11,5,2)-difference sets in Z_11: QR and NQR, which give isomorphic tournaments.
**Impact:** INV-068 corrected. MEMORY.md and TANGENTS.md updated.

### INV-144: Circulant Digraph Path Homology (arXiv 2602.04140, Feb 2026) Ã¢â‚¬â€ CONJ 4.8 DISPROVED
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

### INV-145: ÃŽÂ©_2 Structure Ã¢â‚¬â€ Cancellation Chains in Tournaments
**Source:** opus-2026-03-08-S40
**Status:** DISCOVERED
**What:** ÃŽÂ©_2 Ã¢â€°Â  span(transitive triples). Non-transitive 2-paths with shared non-allowed faces form "cancellation chains" in ÃŽÂ©_2. Gap dim(ÃŽÂ©_2) - |TT| ranges 0-5 at n=5. Cancellation chains never individually in ker(Ã¢Ë†â€š_2), but mixed elements (TT + cancellation) can be 2-cycles.
**Impact:** Previous ÃŽÂ²_2 analysis assumed ÃŽÂ©_2 = TT, which was incomplete. Corrected computation still gives ÃŽÂ²_2 = 0 through n=6 (exhaustive).

### INV-146: P_11 Path Homology Ã¢â‚¬â€ Non-palindromic ÃŽÂ© Dims
**Source:** opus-2026-03-08-S40
**Status:** COMPUTING (dims 8-10 in progress)
**What:** P_11 per-eigenspace ÃŽÂ© dims: [1, 5, 20, 70, 205, 460, 700, 690, ?, ?, ?].
Inner sequence NOT palindromic: 460Ã¢â€°Â 700, 700Ã¢â€°Â 690. Contrasts with P_7's palindromic [3,6,9,9,6,3].
Using J^H J + eigvalsh method for memory-efficient rank computation.
**Next step:** Complete ÃŽÂ©_8, ÃŽÂ©_9, ÃŽÂ©_10 to determine Betti concentration dimension.

### INV-147: Eigenspace Decomposition of ÃŽÂ²_top Ã¢â‚¬â€ Trivial vs Non-trivial Split
**Source:** kind-pasteur-2026-03-08-S41
**Status:** VERIFIED (P_7 and Z_9)
**What:** For circulant maximizers, ÃŽÂ²_top decomposes across Z/nZ eigenspaces as:
- P_7: trivial (k=0) gives ÃŽÂ²_4=0, each non-trivial (k=1..6) gives ÃŽÂ²_4=1, total = 6
- Z_9 S={1,5,6,7}: trivial (k=0) gives ÃŽÂ²_5=2, each non-trivial (k=1..8) gives ÃŽÂ²_5=1, total = 2+8=10
- Formula: ÃŽÂ²_top = (n-1) + ÃŽÂ´, where ÃŽÂ´=0 for prime n (P_7), ÃŽÂ´=2 for n=9 (9=3Ã‚Â²)
- All eigenspaces have IDENTICAL Om_5 dim=74 and Om_6 dim=63
- The difference: trivial has ker(Ã¢Ë†â€š_5)=39, non-trivial has ker(Ã¢Ë†â€š_5)=38
**Key question:** Why does trivial eigenspace contribute extra at n=9? Is ÃŽÂ´=2 because 9=3Ã‚Â²? What is ÃŽÂ´ for n=11 (Paley, prime)?
**Conjecture (HYP-212):** ÃŽÂ´=0 for prime n, ÃŽÂ´>0 for composite. CONFIRMED for P_11: ÃŽÂ²_8 = 10 = p-1 + 0 (opus-S42 + kind-pasteur-S41 independent confirmation).
**P_11 data (opus-S42):** Om dims (kÃ¢â€°Â 0): [1, 5, 20, 70, 205, 460, 700, 690, 450, 180, 30]. ÃŽÂ²_8^(triv)=0 (kind-pasteur-S41 confirmed), ÃŽÂ²_8^(kÃ¢â€°Â 0)=1 per eigenspace.
**Scripts:** `04-computation/n9_beta5_eigenspace.py`, `04-computation/p7_eigenspace_verify.py`, `04-computation/p11_beta8_v5.py`
**Next step:** (1) Test another composite (n=15?). (2) Find algebraic reason for ÃŽÂ´. (3) Extend to P_13?

### INV-148: Arc-Flip Induction Proof for ÃŽÂ²Ã¢â€šâ€š=0 Ã¢â‚¬â€ STRONGEST LEAD
**Source:** kind-pasteur-2026-03-08-S41, opus-2026-03-08-S43 (arc-flip local invariance)
**Status:** VERIFIED EXHAUSTIVELY n=5,6. Key structural mechanism identified.
**What:** ÃŽÂ²Ã¢â€šâ€š=0 can potentially be proved by arc-flip induction:
1. Base: transitive tournament (ÃŽÂ²Ã¢â€šâ€š=0 trivially)
2. Step: flipping any arc preserves ÃŽÂ²Ã¢â€šâ€š=0
The "surplus" = dim(ÃŽÂ©Ã¢â€šÆ’) - dim(ZÃ¢â€šâ€š) satisfies surplus Ã¢â€°Â¥ |drop| ALWAYS.
**Key findings (kind-pasteur-S41):**
- n=5: 10240 arc flips, 0 violations. Surplus=0 cases: max_drop=0
- n=6: 491520 arc flips, 0 violations. Surplus=1 cases (tightest): max_drop=0
- Surplus=0 stability mechanism: joint (ÃŽÂ´ÃŽÂ©Ã¢â€šÆ’, ÃŽÂ´ZÃ¢â€šâ€š) = {(0,0), (1,0), (2,1), (4,2)} Ã¢â‚¬â€ always ÃŽÂ´ÃŽÂ©Ã¢â€šÆ’ Ã¢â€°Â¥ ÃŽÂ´ZÃ¢â€šâ€š
- "Every new ZÃ¢â€šâ€š cycle comes with at least one new ÃŽÂ©Ã¢â€šÆ’ chain to fill it"
- This is the 2-for-1 principle: tournament completeness ensures enough ÃŽÂ©Ã¢â€šÆ’ for every ZÃ¢â€šâ€š
**Key lemma needed:** For any tournament T with arc (u,v), and T' = flip(T, u, v):
  dim(ÃŽÂ©Ã¢â€šÆ’(T')) - dim(ZÃ¢â€šâ€š(T')) Ã¢â€°Â¥ dim(ÃŽÂ©Ã¢â€šÆ’(T)) - dim(ZÃ¢â€šâ€š(T))  when starting surplus Ã¢â€°Â¤ 1
  (or more generally: surplus(T') Ã¢â€°Â¥ 0)
**Why completeness matters:** Non-tournament arc flips CAN create ÃŽÂ²Ã¢â€šâ€š>0 (seen in circulant digraphs).
The tournament constraint ensures every pair of vertices has an arc, providing the intermediary
vertices needed for ÃŽÂ©Ã¢â€šÆ’ chains.
**NEW FINDINGS (kind-pasteur-S41 continued):**
- THM-121 (was THM-100) PROVED: delta_|A_3| = (n-3)*delta_|A_2| exactly, for ALL tournaments, ALL arcs
- delta_|A_2| = 2*(d_u - d_v - 1) depends ONLY on out-degrees
- n=7 sampling (10k): 0 violations, min surplus = 9
- n=8 sampling (20k): 0 violations, min surplus <= 25
- Min surplus floor: 0, 1, 9, <=25 for n=5,6,7,8 Ã¢â‚¬â€ grows super-linearly
- Transitive tournament: surplus = C(n-1,4). One-flip delta = -(n-2-gap)
- DT paths: |DT| >= dim(Z_2) for 100% (n=5), 97.1% (n=6). Rest filled by cancellation
- Omega_2 NOT just TT paths: dim(O2) > |TT| for 76.6% (n=5), 94.6% (n=6)
- Tang-Yau Cor 3.15: H_m=0 for m>=2 when S={1,...,d} Ã¢â‚¬â€ applies to circulant tournaments
- Algebraic identity: surplus = beta_3 + rk(d_4) - beta_2
**Scripts:** beta2_arcflip_proof.py, beta2_surplus_zero_stability.py, beta2_arcflip_n7_sample.py,
  beta2_arcflip_mechanism.py, beta2_arcflip_counting.py, beta2_delta_ratio_*.py,
  beta2_min_surplus*.py, beta2_omega_ratio.py, beta2_injectivity_analysis.py, beta2_surplus_formula.py
**Next step:** (1) Prove the key lemma: surplus(T') >= 0 for all arc flips, using THM-121
  (2) Generalize Tang-Yau deformation retract to non-circulant tournaments
  (3) Prove beta_2 = 0 by induction on number of flips from transitive

### INV-149: ÃŽÂ²Ã¢â€šâ€š=0 Density Threshold for Circulant Digraphs
**Source:** kind-pasteur-2026-03-08-S41
**Status:** CHARACTERIZED. New conjecture (HYP-219).
**What:** For C_n^S with SÃ¢Ë†Â©(-S)=Ã¢Ë†â€¦, ÃŽÂ²Ã¢â€šâ€š=0 when |S| is large enough:
- n=7: |S|Ã¢â€°Â¥3 (all ÃŽÂ²Ã¢â€šâ€š=0)
- n=9: |S|Ã¢â€°Â¥4 (all ÃŽÂ²Ã¢â€šâ€š=0)
- n=11: |S|Ã¢â€°Â¥4 (all ÃŽÂ²Ã¢â€šâ€š=0)
- n=13: |S|Ã¢â€°Â¥5 (all ÃŽÂ²Ã¢â€šâ€š=0, 96/96 tested)
- n=15: |S|=5 still has 15/201 failures; threshold at |S|Ã¢â€°Â¥6?
- |S|=2 perfect characterization: ÃŽÂ²Ã¢â€šâ€š=0 iff has-doubling-pair (HYP-217)
- Exceptions without doubling pair: coset structure (S = a + H for subgroup H)
**Scripts:** beta2_doubling_closure_general.py, beta2_threshold_analysis.py
**Next step:** (1) Find exact threshold formula (2) Prove for tournaments (|S|=(n-1)/2)

### INV-148: Nesting Obstruction Theory Ã¢â‚¬â€ Why H=7 is Simplex-in-Cuboid
**Source:** opus-2026-03-14-S71f
**Status:** DISCOVERED. Algebraic framework established, geometric interpretation given.
**What:** H=7 = (1+2(1+x))|_{x=2}: composing simplex (1+x) into cuboid (1+2y) gives 3+2x, which has constant term 3 Ã¢â€°Â  1 Ã¢â€ â€™ NOT a valid independence polynomial. Tournament geometry prevents nesting/composition of brick structures, leaving only multiplicative (disjoint union) decomposition.
**Key results:**
- H=7 impossible = tesseract impossibility (HYP-1230)
- H=21 = 3Ãƒâ€”7 inherits obstruction multiplicatively (THM-079)
- ÃŽÂ±Ã¢â€šÂ=3 forces ÃŽÂ±Ã¢â€šâ€šÃ¢â€°Â¥2 at n=7,8 (verified 500k samples each) Ã¢â€ â€™ H=15 not 7
- Structural proof for all n sketched (pigeonhole on 9 vertex slots in nÃ¢â€°Â¥7 vertices)
**Scripts:** nesting_obstruction.py, alpha1_3_n8_verify.py
**Next steps:** (1) Complete rigorous proof that ÃŽÂ±Ã¢â€šÂ=3 Ã¢â€ â€™ ÃŽÂ±Ã¢â€šâ€šÃ¢â€°Â¥2 for ALL n. (2) Publish nesting obstruction as appendix to main paper.

### INV-149: The (z-2)(z-3) = 0 Recurrence Ã¢â‚¬â€ Simplex-OCF Bridge
**Source:** opus-2026-03-14-S71f
**Status:** DISCOVERED. Complete algebraic analysis.
**What:** The characteristic polynomial (z-2)(z-3) with roots 2 (OCF point) and 3 (simplex value) generates the forbidden sequence via the pure z=3 orbit from seed 7: {7, 21, 63, 189, ...}. Only first two are permanently forbidden; 63+ are achievable at nÃ¢â€°Â¥8. Mixed orbit H = 3^{k+1} - 2^{k+1} gives all achievable values.
**Scripts:** knacci_tournament_recurrence.py
**Next steps:** (1) Does the recurrence connect to deletion-contraction? (2) Higher-order recurrences (z-2)(z-3)(z-5)=0 for cuboid inclusion?

### INV-150: The 2-Bridge Ã¢â‚¬â€ Unified Origin of the Factor 2
**Source:** opus-2026-03-14-S71f (connecting kind-pasteur S72 Degree Drop Theorem)
**Status:** DISCOVERED. Three manifestations verified.
**What:** The number 2 appears in: (a) I(ÃŽÂ©,2)=H (OCF evaluation), (b) top-degree coefficients Ã‚Â±2 (Degree Drop Theorem), (c) ÃŽâ€H=2Ã‚Â·ÃŽâ€ÃŽÂ±Ã¢â€šÂ (arc flip derivative). All arise from the binary arc choice / path reversal involution. 1+(-1)^{n-1} = 2 for odd n.
**Scripts:** degree_drop_packing.py, delta_h_gap.py
**Next steps:** (1) Prove ÃŽâ€H=2Ã‚Â·ÃŽâ€ÃŽÂ±Ã¢â€šÂ+4Ã‚Â·ÃŽâ€ÃŽÂ±Ã¢â€šâ€š for general n. (2) Connect to Vassiliev theory formally.

### INV-151: Simplicial Selection Interpretation
**Source:** opus-2026-03-14-S71f
**Status:** DISCOVERED. Geometric framework.
**What:** Each tournament T selects H(T) simplices from the standard n!-simplex triangulation of [0,1]^n. H/n! = 1/2^{n-1} on average. The f-polynomial of ÃŽâ€^{n-1} at x=2 gives 3^n (simplex brick), and the f-polynomial of Ã¢â€“Â¡^n at x=2 gives 5^n (cuboid brick).
**Scripts:** simplex_cuboid_geometry.py
**Next steps:** (1) Is the selection pattern random-looking or structured? (2) Does the geometric constraint explain the ÃŽâ€H gap?

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

### INV-182: Soft OCF Ã¢â‚¬â€ Differentiable Independence Polynomial
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
**What:** Tournaments are antiferromagnets on the staircase lattice. The frustration index cÃ¢â€šÆ’/C(n,3) correlates with H at rÃ¢â€°Ë†0.97. Score variance = NÃƒÂ©el order parameter. Regular tournaments = AFM ground state. Magnon spectrum is flat over labeled ensemble (S_n isotropy). Boltzmann-weighted correlations show true AFM order at ÃŽÂ²>0. Phase transition at ÃŽÂ²_cÃ¢â€°Ë†0.7 (n=5). H=7 gap explained by frustration propagation (ÃŽÂ±Ã¢â€šÂ gap: 3 three-cycles force a 5-cycle).
**Key files:** antiferromagnetic_tournament_s15.py, afm_deep_analysis_s15.py, afm_remaining_s15.py, the-antiferromagnetic-tournament.md
**Open directions:**
1. Per-class magnon dispersion (condition on iso class for anisotropy)
2. ÃŽÂ²_c(n) scaling as nÃ¢â€ â€™Ã¢Ë†Å¾ (thermodynamic limit)
3. Yang-Lee zeros of Z(ÃŽÂ²) in complex plane
4. Connection between seesaw ÃŽÂ²Ã¢â€šÂÃ‚Â·ÃŽÂ²Ã¢â€šÆ’=0 and AFM selection rules
5. ÃŽÂ±Ã¢â€šÂ gap mechanism for H=21 (needs nÃ¢â€°Â¥7 analysis)
6. The staircase lattice as Brillouin zone of the tournament AFM
**Priority:** HIGH (provides physics intuition for all open questions)

### INV-052: Paley maximization vs Interval Ã¢â‚¬â€ complete status (opus-2026-04-16-S1)
**Source:** opus-2026-04-16-S1 (this session)
**Status:** RESOLVED for p=7,11,19,23; OPEN for pÃ¢â€°Â¥31
**What:** Does Paley tournament T_p maximize H among all circulant tournaments? Among all tournaments?
**Results:**
- p=7:  H(Paley)=189 > H(Interval)=175 Ã¢â€ â€™ Paley wins. EXHAUSTIVE: Paley = global max at n=7.
- p=11: H(Paley)=95,095 > H(Interval)=93,027 Ã¢â€ â€™ Paley wins (exhaustive among circulants).
- p=19: H(Paley)=1,172,695,746,915 < H(Interval)=1,184,212,824,763 Ã¢â€ â€™ Interval wins.
- p=23: H(Paley)=15,760,206,976,379,349 < H(Interval)=16,011,537,490,557,279 Ã¢â€ â€™ Interval wins by 1.57%.
**Mechanism:** At p=23, Paley has MORE ÃŽÂ±Ã¢â€šÂ and ÃŽÂ±Ã¢â€šâ€š than Interval but loses on ÃŽÂ±Ã¢â€šÆ’+.
  The term 8(ÃŽÂ±Ã¢â€šÆ’_I - ÃŽÂ±Ã¢â€šÆ’_P) = +152T swamps the combined 2ÃŽâ€ÃŽÂ±Ã¢â€šÂ+4ÃŽâ€ÃŽÂ±Ã¢â€šâ€š = +86T advantage.
**Dominant term crossover:** k=1Ã¢â€ â€™2 at nÃ¢â€°Ë†11; k=2Ã¢â€ â€™3 predicted at nÃ¢â€°Ë†29. The PaleyÃ¢â€ â€™Interval crossover
  happens between p=11 and p=19, likely because k=2 becomes the dominant term around n=11.
**Open questions:**
  (1) Does INTERVAL maximize H among ALL tournaments (not just circulants) for nÃ¢â€°Â¥15?
  (2) Is the PaleyÃ¢â€ â€™Interval crossover at p=13? (p=13 Ã¢â€°Â¡ 1 mod 4, no Paley, skip to p=19)
  (3) Why does Interval win on ÃŽÂ±Ã¢â€šÆ’,ÃŽÂ±Ã¢â€šâ€ž,...? Is this related to arithmetic spacing allowing more vertex-disjoint packings?
**Scripts:** alpha_from_cc_bin.py, alpha_crossover_analysis.py, alpha_ratio_trends.py

---

## 2-Adic Column Family Investigations (oracle-2026-05-15)

### INV-187: SC Blowup Ã¢â‚¬â€ $\Omega(T_{\mathrm{SC}})$ Structure and Cross-Lane Cycles
**Source:** oracle-2026-05-15-S2 (sc-blowup-and-twin-gaining.md)
**Status:** OPEN. Universal Score Theorem proved. H concentration data at n=3,4,5.
**Key results (verified exhaustively n=3,4,5):**
- Universal Score: all vÃ¢â€šâ‚¬ have out-degree n, all vÃ¢â€šÂ have out-degree n-1 (PROVED).
- H_SC varies only 4.2% at n=5 (14937Ã¢â‚¬â€œ15565). Minimized by transitive, maximized by regular.
- H_SC = H_Lex iff T is regular (at n=3,5; first departure at Paley(7): H_SC=24453597 Ã¢â€°Â  H_Lex=24589929).
- SC preservation: T_SC is SC iff T is SC (PROVED via anti-automorphism Ãâ€ž(v_i) = ÃÆ’(v)_{1-i}).
- Kronecker: A(T_SC) = A(T)Ã¢Å â€”IÃ¢â€šâ€š + A(T)Ã¡Âµâ‚¬Ã¢Å â€”ÃŽÂ¦ + I_nÃ¢Å â€”eÃ¢â€šâ‚¬Ã¢â€šÂ.
- Eigenvalue splitting (circulant): ÃŽÂ»_{k,0}(T_SC) = 2Re(ÃŽÂ»_k(T))+1, ÃŽÂ»_{k,1} = 2iÃ‚Â·Im(ÃŽÂ»_k(T))Ã¢Ë†â€™1.
**Open questions:**
1. What is ÃŽÂ©(T_SC) in terms of ÃŽÂ©(T)? Same-lane cycles = copies of ÃŽÂ©(T); cross-lane cycles = ?
2. Is there a formula H(T_SC) = f(I(ÃŽÂ©(T), x)) for some x or some operation on the polynomial?
3. Why does T_Lex Ã¢â€°â€¦ T_SC for regular n=3,5 but NOT for Paley(7)? What breaks at n=7?
4. Tower: H(Trans_n^SC) for n=3..?: 1, 41, 530293, ... No recurrence found. OEIS?
5. HYP-SC-1 (H monotonicity): H(TÃ¢â€šÂ) Ã¢â€°Â¤ H(TÃ¢â€šâ€š) Ã¢Å¸Â¹ H_SC(TÃ¢â€šÂ) Ã¢â€°Â¤ H_SC(TÃ¢â€šâ€š)? (Checked n=3,4; needs n=5 verify.)
**Scripts:** blowup_study3.py, results: blowup_study.out
**Next steps:** (1) Compute ÃŽÂ±Ã¢â€šÂ(T_SC) for C3 and Trans-3. (2) Find ÃŽÂ©(T_SC) for n=3 cases by direct cycle enumeration. (3) Look for H_SC formula via independence polynomial at different evaluation point.

### INV-184: Tournament Blowup $T[K_2]$ Ã¢â‚¬â€ H Formula and Family Inheritance
**Source:** oracle-2026-05-15 (2-adic column family analysis; see `07-reflections/adic-column-families.md`)
**Status:** OPEN. Computationally accessible immediately.
**What:** The tournament blowup $T[K_2]$ Ã¢â‚¬â€ each vertex $v$ splits into $v_0 \to v_1$,
each original arc $u \to v$ expanded to all four arcs $u_i \to v_j$ Ã¢â‚¬â€ is the concrete
realization of the "row step" in the 2-adic column family grid: $(r, k) \to (r+1, k)$.
It doubles $n$ and stays within the same column family.
**Key questions:**
1. Is there a formula $H(T[K_2]) = f(H(T), n)$?
2. Does blowup preserve SC status? SF status? If yes, the column family inherits SC/SF structure.
3. The pairs anomaly: $\lfloor n/2 \rfloor$ gains +1 extra at the oddÃ¢â€ â€™even ($r=0 \to r=1$)
   transition only. Does H have a similar anomalous first-blowup jump?
4. What is the isomorphism class of $T[K_2]$? Does it stay in the same G_n/Z_2 sector?
**Expected behavior:** Blowup creates a natural Z_2 action (swap $v_0 \leftrightarrow v_1$
is NOT an automorphism, but the structure is highly symmetric). This may force the blowup
into SC or near-SC classes.
**Connection:** Linial-Morgenstern conjecture (INV-013) uses "random blowup of transitive
tournament" Ã¢â‚¬â€ the row step applied to the transitive class. Our framework says this is the
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
- Blowup of cyclic C3 (max-H at n=3) Ã¢â€ â€™ H=45 = **max H at n=6** Ã¢Å“â€œ
- Blowup is near-regular: scores $(n/2-1)^{n/2} \cup (n/2)^{n/2}$ Ã¢â‚¬â€ exactly the SCÃ¢Ë†Â©SF score.
- Blowup of max-H regular tournament Ã¢â€ â€™ max-H SCÃ¢Ë†Â©SF class at the doubled size.
- This matches the max-H pattern: "even $n$ Ã¢â€ â€™ max H by SCÃ¢Ë†Â©SF class."

**CONJECTURE (HYP-new):** Blowup of the max-H regular tournament at odd $n$ IS the
max-H SCÃ¢Ë†Â©SF class at even $2n$. Equivalently: the column-family row step maps
max-H-regular to max-H-SCÃ¢Ë†Â©SF.

**Next steps:**
1. Verify: is H=15565 = max H at n=10? (Run exhaustive or sampling-based search.)
2. Verify: is H=393 = max H at n=8? Compare to known SCÃ¢Ë†Â©SF n=8 classes.
3. Prove: blowup of regular tournament Ã¢â€ â€™ near-regular tournament Ã¢â€ â€™ SCÃ¢Ë†Â©SF class.
4. Connect to the SCÃ¢Ë†Â©SF = SC(n-2) identity: does the blowup construction provide the bijection?

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

### INV-186: SCÃ¢Ë†Â©SF = SC($n-2$) Identity Ã¢â‚¬â€ Column Family Framing, but Fails n=8
**Source:** oracle-2026-05-15 (2-adic column family analysis)
**Status:** PARTIALLY RESOLVED Ã¢â‚¬â€ identity is a SMALL-N COINCIDENCE, not a theorem.
SCSF(8)=5 Ã¢â€°Â  12=SC(6) (oracle-2026-05-11-S2). The column-family framing explains WHY
the coincidence holds through n=7 and why it fails at n=8 (complexity grows faster
than the column structure tracks).
**What:** The identity $\#(\text{SC} \cap \text{SF})(n) = \#\text{SC}(n-2)$ (verified n=4..7)
says: adjacent top-row columns have related SCÃ¢Ë†Â©SF and SC counts. In the 2-adic grid, this is
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
**Status:** OPEN Ã¢â‚¬â€ conjectured (THM-336), verified p=7,11 exactly, p=19,23 as lower bounds
**What:** For prime p Ã¢â€°Â¡ 3 mod 4: H(QR_p) = a(p), H(QR_p - v) = a(p-1).
Why is QR_p globally optimal? Current heuristics: regular + SC + maximum cycle density.
No proof. Strong empirical support.
**Next:**
  1. Prove for p=19: need to verify a(18) from exhaustive/branch-and-bound solver
  2. Exploit the SC structure and THM-334 (SC strict bound) for a lower bound proof
  3. Check if the "c(p)/a(p-1) Ã¢â€ â€™ (p-1)/4" formula has a combinatorial proof

### INV-232: Base-Path Staircase Recurrence Proof
**Source:** opus-2026-05-27-S7
**Status:** OPEN Ã¢â‚¬â€ verified k=1..11, recurrence H(k) = 3H(k-1) + H(k-2) + H(k-3) (THM-337)
**What:** The base-path staircase has order-3 recurrence. WHY? Need algebraic proof.
**Next:**
  1. Set up the F_odd recursion for T_k Ã¢â€ â€™ T_{k+1} insertion
  2. Show the # of odd-length paths between insertion neighborhoods satisfies this recurrence
  3. Check OEIS for the sequence 1,5,17,57,193,653,2209,...

### INV-233: Circulant Optimality Threshold
**Source:** opus-2026-05-27-S7
**Status:** OPEN Ã¢â‚¬â€ threshold at n=13 proved computationally; structural reason unclear
**What:** a(n) = opt_circ(n) for n Ã¢â€°Â¤ 11; a(n) > opt_circ(n) for n Ã¢â€°Â¥ 13 (THM-338).
Key insight: QR_p has both forward and backward connections; purely forward circulants are weaker.
**Next:**
  1. Run n=29 circulant search (2 GB RAM, ~500s) for a(29) lower bound
  2. Try n=31 (8 GB RAM, ~150s) for a(30), a(31) lower bounds from circulant + Paley
  3. Search for proof that n=13 breaks circulant optimality: what structure is exploited?

### INV-234: A038375 Extended Terms
**Source:** opus-2026-05-27-S7; extended monad-researcher-2026-06-03-S578
**Status:** EXTENDED Ã¢â‚¬â€ new bounds for n=13..16 from local search
**Best known lower bounds (all from local hill-climb, may not be exact):**
  a(13) Ã¢â€°Â¥ 3,719,831 (prior best, hardcoded in solver source)
  a(14) Ã¢â€°Â¥ 24,762,119 (NEW S578 Ã¢â‚¬â€ first non-trivial bound for n=14)
  a(15) Ã¢â€°Â¥ 198,464,295 (confirmed from circulant; prior from backlog)
  a(16) Ã¢â€°Â¥ 1,522,320,909 (NEW S578 Ã¢â‚¬â€ first non-trivial bound for n=16)
  a(17) Ã¢â€°Â¥ 13,689,269,499 (prior from backlog)
  a(25) Ã¢â€°Â¥ 2,418,453,569,285,650,675 (prior)
  a(27) Ã¢â€°Â¥ 17,051,631,267,035,242,313 (prior)
**Exact values (Paley tournaments, confirmed):**
  a(18)=117,266,659,317, a(19)=1,172,695,746,915,
  a(22)=1,313,333,107,451,805, a(23)=15,760,206,976,379,349
**Results:** `05-knowledge/results/a038375_n13_16_s578.out`
**Next:**
  1. Try longer searches for n=14, n=16 to improve lower bounds
  2. Submit new OEIS terms: a(12)=531205 (likely exact), a(13)Ã¢â€°Â¥3719831, a(14)Ã¢â€°Â¥24762119, etc.
  3. Check if a(15)=198464295 is exact (circulant is best for odd n; n=15 is not prime but check)

### INV-235: c(p) Odd Cycle Count as New OEIS Sequence
**Source:** opus-2026-05-27-S7
**Status:** OPEN
**What:** c(p) = #{directed odd simple cycles through any vertex in QR_p}:
  c(3)=1, c(7)=72, c(11)=39675, c(19)=527714543799, c(23)=7223436934463772
**Formula:** c(p)/a(p-1) Ã¢â€ â€™ (p-1)/4 (empirical, very tight).
**Next:** Check OEIS for this sequence; submit if new.

## INV-NEW-S2-A: Submit non-SC and SC tiling sequences to OEIS (opus-2026-05-27-S2)
**Source:** Session computation, THM-336/337
**Status:** EXTENDED to n=21 (monad-researcher-2026-06-03-S578). OEIS submission still pending.
**Next step:** Submit to OEIS. Both likely new. See results `sc_nonsc_tiling_extended_s578.out` and `sc_nonsc_analysis_s578.out`.
**Sequence 1 (non-SC, n=3..21):** 1, 3, 14, 121, 1995, 64648, 4163979, 534849295, 137175056830, 70300582005021, 72022549494074543, 147537994739778382720, 604389195176853420636135, 4951458073552104202450472659, 81127163139584624300444795370894, 2658415431704809155949894648688347153, 174223242716674181161576562635131644182915, 22835875973668207070158505359404862076014660344, 5986299288700071914625804856670204247617947928012339
**Sequence 2 (SC, n=3..21):** 1, 5, 50, 903, 30773, 2032504, 264271477, 68184627441, 35047197032002, 35958496436958947, 73714953745344131921, 302083916908917515293824, 2475275689375583696377612313, 40559867749229788743692052099373, 1329146868621776288279506615484973682, 87109627516328541837467949607883973785583, 11417807318404962374285126179033325959417790077, 2993132517377715508274076378897588219393273833354504, 1569269447547381490887032729997059933821341243168080615885
**Formula:** sc_ie(n) gives both exactly. Fast computation via subset inclusion-exclusion.
**NEW FINDING (S578):** correction(n) := 2^{m-n+3} - non-SC(n) satisfies correction(n) ~ SC(n-2) asymptotically (ratio Ã¢â€ â€™ 1 rapidly). Moreover, correction(n)/correction(n-1) Ã¢â€ â€™ 2^{n-4}+2.

## INV-NEW-S2-B: Exact count for d good cuts Ã¢â‚¬â€ general formula (opus-2026-05-27-S2)
**Source:** THM-336, HYP-1742
**Status:** OPEN
**Observation:** d=2Ã¢â€ â€™nÃ¢Ë†â€™2, d=3Ã¢â€ â€™5(nÃ¢Ë†â€™3), d=4Ã¢â€ â€™(n-4)(n+95)/2. Pattern suggests:
  exactly-d-good(n) = (n-d)Ã‚Â·SC(d+1) + (non-consecutive contribution depending on d and n)
**Next step:** Compute d=5 formula. From data: n=8Ã¢â€ â€™2739, n=9Ã¢â€ â€™3672, n=10Ã¢â€ â€™4615. These give:
  consecutive contribution: (n-5)Ã‚Â·SC(6) = (n-5)Ã‚Â·903. Non-consecutive: 2739-3*903=30, 3672-4*903=60, 4615-5*903=100.
  Non-consecutive = 30, 60, 100 with differences 30, 40. Second diff = 10. So quadratic + 10*(n-8)Ã‚Â²/2+...
**Connection:** d=4 non-consecutive = C(n-4,2); d=5 non-consecutive seems quadratic in n.

## INV-NEW-S2-C: Prove non-SC(n) ~ 2^{m-n+3} rigorously (opus-2026-05-27-S2)
**Source:** HYP-1744
**Status:** OPEN
**Approach:** The IE formula gives non-SC = 2Ã‚Â·2^{m-n+2} + (smaller terms). Need to bound the sum of correction terms. The correction = IE at sizes Ã¢â€°Â¥ 2, which contributes negative+positive terms.
**Key:** Size-2 IE term at pairs {1,k}: f({1,k}) = 1*(n-k)-? Ã¢â€ â€™ smaller exponent than size-1.

## INV-NEW-S2-D: SC tiling sequence vs A054946 relationship (opus-2026-05-27-S2)
**Source:** Session sequence computations
**Status:** OPEN
**Observation:** SC tiling counts (1,5,50,903,...) are NOT the same as A054946 (1,0,2,24,544,22320,...).
**Reason:** SC tilings fix the base path; A054946 counts labeled tournaments (all 2^{C(n,2)} orientations).
**Question:** P(SC tiling) Ã¢â€°Â  P(SC labeled tournament). n=5: 50/64=0.781 vs 544/1024=0.531. Is there a clean relationship?

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
**What:** Formalize and investigate the good-cut bucket `g(Ãâ€ž)=|G(Ãâ€ž)|`, where `G(Ãâ€ž)` is the set of base-path cuts crossed by at least one upward tile. `TournamentH7.GoodCuts` proves bucket 0 iff all-down, bucket 1 impossible, and grid reflection invariance without project-specific axioms. The companion exact census finds every merged tournament class pure in `g` for n=3..6.
**S14 update:** Lean now also proves the interval-union characterization, monotonicity, bucket bounds `{0}Ã¢Ë†Âª{2,...,n-1}`, and top-bucket iff every legal cut is good. A new projection-defect cross-scan found that single-tile lines with `|Delta g|>0` always change the merged tournament class through n=6, and tile range parity controls defect polarity: even ranges are even-graph biased, odd ranges tournament-class biased.
**Codex 2026-05-30 update:** Lean now proves the exact abstract spectrum:
for n>=3, `exists b, b.goodCutCount = r` iff `r=0` or `2<=r<=n-1`.
The proof is constructive via one upward tile covering `{1,...,r}`. This
shows bucket 1 is the only interval-geometric obstruction; any further gaps
in quotient/H/isomorphism statistics must come from tournament structure.
**S1 update (kind-pasteur-2026-05-30-S1):** Fast n=7 hash-assisted canonicalization confirms HYP-1764 one level further: 456 unmerged classes, 88 SC classes, 272 merged classes, and every merged class is pure in `g` (pure/mixed = 272/0, max span 0).
**S15 update (opus-2026-05-29-S15):** THM-349 proves the full interval-union recurrence for bucket counts. If `B_N(x)=ÃŽÂ£_g b(N,g)x^g`, then `B_N=B_{N-1}+ÃŽÂ£_{L=2}^N c_L x^L B_{N-L-1}`, where `c_L` is a connected-run cover count computed by inclusion-exclusion. The recurrence matches direct tiling enumeration for n=3..8 and gives exact counts through n=13 without enumerating the tiling cube. Lean now has the direct membership form `k Ã¢Ë†Ë† goodCuts Ã¢â€ â€ Ã¢Ë†Æ’ upward tile interval containing k`. The S15 transport-excess scan also extends the dynamic evidence: across selected Hamming layers through n=6, every ordered half-line with nonzero `Delta g` changes merged tournament class (`bad=0` throughout), while `Delta g=0` is where self, spine, ribs, sea, and even-only defects mix.
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

## LEAD (opus-2026-06-02-S558): LRC proof-methodology map + two corrections + the 2Ã‚Â·7 target
**Source:** literature survey (verified June 2026), `07-reflections/lrc-proof-methodologies-by-n-the-complete-map-s558.md`.
**Status of LRC (n=total runners, gap 1/n = repo convention):** PROVEN through **n=13** Ã¢â‚¬â€ structural era nÃ¢â€°Â¤7 (n=4 Betke-Wills/Cusick view-obstruction; n=5 Cusick-Pomerance; n=6 Bohman-Holzman-Kleitman averaging; n=7 Barajas-Serra circular chromatic number of distance graphs), then finite-checking era nÃ¢â€°Â¥8 (Tao reduction + Malikiosis-Santos-Schymura bound `Ã¢Ë†ÂuÃ¡ÂµÂ¢<B_k` + Rosenfeld divisibility sieve Ã¢â€¡â€™ prime-product contradiction; n=8 Rosenfeld 2025, n=9,10 Trakulthongchai 2025, **n=11,12,13 Sungkawichai-Trakulthongchai Apr 2026, arXiv:2604.23906**). **n=14 = immediate open frontier.**
**WHY n=14 is the wall (literature's own reason):** the polynomial-method shortcut that handles the tight tuple (1,Ã¢â‚¬Â¦,k) analytically requires **k+1 to be an odd prime**; for k=13, k+1=**14=2Ã‚Â·7 is composite**, so it fails and the full c=14 sieve lift is needed (infeasible; k=12 ~40 days/10 cores). This is EXACTLY the repo's 2Ã‚Â·7 / "7-impossibility" / even-fold / mod-7 thread Ã¢â‚¬â€ not a metaphor, the actual obstruction.
**CORRECTION 1:** LRC is proven to n=13, NOT n=7 (our prior framing).
**CORRECTION 2:** the even-fold (S554, HYP-2056) `M(S)Ã¢â€°Â¤M(fold)` should use **LRC(13)**: any primitive 13-set has Ã¢â€°Â¤12 even speeds Ã¢â€¡â€™ |fold|Ã¢â€°Â¤12 Ã¢â€¡â€™ M(fold)Ã¢â€°Â¥1/13 by proven LRC(13) Ã¢â€¡â€™ the EVEN half of LRC@14 is fully protected for EVERY config; residual = odd coupling only (the eÃ¢â€°Â¤6 restriction was an artifact of using LRC(7)).
**Repo's sieve THM-369 = Rosenfeld's divisibility sieve** (independently rediscovered the modern engine); pair it with the MSS finite bound (HYP-2052 shows why the bounded sieve alone can't close it).
**TOP TARGET:** find the `k+1=2q`-analogue of "the tight tuple (1,Ã¢â‚¬Â¦,k) is proper" Ã¢â‚¬â€ an algebraic substitute for the polynomial method when k+1 is twice an odd prime (k+1=14=2Ã‚Â·7). This is the single most leveraged route to n=14 and squarely in the repo's wheelhouse.
**S593 cap-face update:** THM-398 now includes Lemma H, the dual n-clock cap pigeonhole: for `v=nw`, each primary `n`-clock cell has danger capacity exactly `2/n^2`; if `G(S')` has more mass than that in any cell, the multiple row is loose.  The S593 audit routes `2460/2500` deterministic multiple-of-14 rows by this aggregate cell-cap criterion, with the remaining `40` routed by S581 owner descent.  Next target: on rows where every cell is under capacity, combine endpoint-owner congruences with `Phi` ramps to force positive gap.
**Codex S266 circuit update:** HYP-3116 identifies the old HYP-2108/HYP-2112
endpoint-cover activation circuit as the concrete circuit-complexity bridge:
`P(S)` is the max activation and `Phi(C)=G(v)` is the exact gap sum.  The next
proof-facing task is to attach `endpoint_cover_activation_vector`,
`phi_gap_sum`, `phi_kernel_status`, `P_max_activation`,
`endpoint_period_numerator_sidecar`, `finite_address_packet`, and
`observer_gluing_certificate` to HYP-2963/HYP-3107 residual packets, then prove
the legal residual `Phi` kernel is empty or route the first missing input to
named debt.  HYP-2791's low-depth Boolean cut should be tested as a lower bound
on `Phi`, not as a standalone scalar proof.

## LEAD (monad-explorer-S4): Moser-lattice unit-distance toolkit + the Ã¢Ë†Å¡(4tÃ¢Ë†â€™1) angle ladder
- **Source:** THM-432, HYP-2298, reflection the-moser-lattice-is-the-bridge-ring-s4.md.
- **Status:** OPEN. u(21)=57 is PROVEN (Alexeev-Mixon-Parshall 2025); exact-value frontier is now n=22 (60Ã¢â€°Â¤u(22)Ã¢â€°Â¤61). Retire any "prove u(21)" tasking.
- **Next steps:** (1) exact-arithmetic Moser-lattice toolkit Ã¢â‚¬â€ unit-vector counts of Ã¢â€žÂ¤[ÃŽÂ¶6, Ãâ€°_t] for the angle ladder Ãâ€°_t=((2tÃ¢Ë†â€™1)+iÃ¢Ë†Å¡(4tÃ¢Ë†â€™1))/(2t), t=1,2,3,Ã¢â‚¬Â¦ (discriminants Ã¢Ë†Å¡3,Ã¢Ë†Å¡7,Ã¢Ë†Å¡11,Ã¢Ë†Å¡15,Ã¢â‚¬Â¦); the repo's Ã¢Ë†Å¡3 (t=1) and Ã¢Ë†Å¡7 (t=2) are rungs of it. (2) Is the additive-norm-Ã¢Ë†Å¡7 layer (HYP-2262/THM-421) the "same 7" as the t=2 Moser angle? (different CM fields Ã¢â€žÅ¡(Ã¢Ë†Å¡Ã¢Ë†â€™3) vs Ã¢â€žÅ¡(Ã¢Ë†Å¡Ã¢Ë†â€™7) Ã¢â‚¬â€ open). (3) HYP-2170 tie: Moser units (18) enlarge Cay(Ã¢â€žÂ¤[ÃŽÂ¶6],U6); LRC has the analogous clockÃƒâ€”shell two-tower (THM-427) Ã¢â‚¬â€ both are products of two cyclotomic/CM pieces. Engineering deliverable: the Moser-lattice exact UD counter (seed: unit_distance_moser_lattice_u21_monad_s4.py).
- **NEW (S4):** L_13 = Z[zeta6, w_13] (disc sqrt51) has 24 unit vectors (max degree 24 > Moser's 18). Explore whether it yields denser unit-distance graphs at OPEN n (>=22, where u(n) is only bounded 60<=u(22)<=61, ...). Exact counts only. (moser_angle_ladder_unitvectors_monad_s4.py)

## LEAD (monad-explorer-2026-06-07-S6): the icosahedral (2,3,5) handle for THM-436 / HYP-2305
- **Source:** THM-436 ADDENDUM, HYP-2305, reflection `the-icosahedral-fifteen-s6.md`, `04-computation/icosahedral_fifteen_monad_s6.py` + `icosahedral_flag_fibers_monad_s6.py`.
- **PROVED/VERIFIED (use freely):** #overlapping cyclic-triangle pairs on [n] = 15Ã‚Â·C(n,5) (oriented 60Ã‚Â·C(n,5)); canonical bijection overlapping-pairs {X,Y} Ã¢Å¸Â· involutions (ab)(cd) of A_5 (shared vertex=fixed point); every oriented overlapping pair's commutator is a 3-cycle, onto all 20; icosahedral axis-counts {6,10,15}={n_5, n_3, #involutions of A_5} (cyclic Sylows p=3,5 give axes=#subgroups; p=2 deviates, V_4). C_5 realises 5 of the 15.
- **REFUTED (do not retry):** "60=20Ã‚Â·3=icosahedral face-vertex flags, flagÃ¢â€ â€™face uniformly 3-to-1" Ã¢â‚¬â€ commutator fibers are NON-uniform {2(Ãƒâ€”3),3(Ãƒâ€”14),4(Ãƒâ€”3)} (MISTAKE-059). Also: G_5's f-vector (12,30,20) icosahedron (the-five-platonic-tournaments Ã‚Â§IV) is NOT the A_5/Galois icosahedron (different degree sequence, no A_5 action) Ã¢â‚¬â€ a matching f-vector is not an isomorphism.
- **OPEN / testable handle (HYP-2305, the (2,3,5)Ã¢â€ â€carry-prime conjecture):** the icosahedron's three axis-orders {2,3,5} mirror the repo's three carry-prime frontiers Ã¢â‚¬â€ prime-2 (doubling, THM-404), prime-3 (THM-407/428, the 3Ã‚Â³=27 shell at hard n=14), n=5/cyclotomic-5 (THM-403/436). NEXT STEP: in the worry-set floor data (S612 Res_27, THM-428 shell towers), find the smallest n whose binding shell-partner needs a factor-5 / 5-fold carry beyond the existing 2- and 3-towers Ã¢â‚¬â€ predicted to be the genuinely "icosahedral" (non-solvable-monodromy, HYP-2303) LRC case. The n=14 hard case is 2Ã‚Â·7 / 3Ã‚Â³ (prime 2,3,7); the 5-fold carry has NOT yet been isolated Ã¢â‚¬â€ this is the missing leg.
- **Number-field tie (unexplored):** the n=5 cyclotomic worry-set lives over Q(ÃŽÂ¶_5), whose real subfield Q(Ã¢Ë†Å¡5) IS the icosahedral coordinate field (Ãâ€ ). Does S699h's A_5 unit-distance Cayley graph (spherical HN) use Ã¢Ë†Å¡5/Ãâ€  explicitly, and is Klein's icosahedral equation TÃ‚Â²+HÃ‚Â³=1728fÃ¢ÂÂµ a handle on its chromatic structure? (THM-436 Ã‚Â§5 / reflection Ã‚Â§5 open question, still untouched.)

## LEAD (monad-explorer-2026-06-07): the 1729 spine is severed Ã¢â‚¬â€ what (if anything) bridges the tournament & Q(Ã¢Ë†Å¡Ã¢Ë†â€™3) lanes?
- **Source:** HYP-2306, reflection `the-1729-resonance-is-isolated-the-tournament-ratio-has-no-modular-structure.md`, `04-computation/paley_ratio_modular_test_monad.py` + `paley_H23_monad.py`.
- **DONE / use freely:** the canonical Paley ratio `r(p)=H(T_p)/|Aut|` has NO modular significance Ã¢â‚¬â€ `1729=r(11)=7Ã‚Â·13Ã‚Â·19=j(i)+1` (completely split in Q(Ã¢Ë†Å¡Ã¢Ë†â€™3)) does NOT persist to the next genus-1 Paley prime p=19 (`r(19)=5Ã‚Â·7Ã‚Â·11Ã‚Â·23Ã‚Â·774463`, 5/11/23 inert) nor p=23 (`r(23)=3Ã‚Â·167Ã‚Â·4567Ã‚Â·27225299`, 167/27225299 inert). 1729's cleanness = `H(T_11)=5Ã‚Â·7Ã‚Â·11Ã‚Â·13Ã‚Â·19` smooth; smoothness breaks at p=19. So the cross-lane "1729 spine" (tournament Ã¢â€ â€ S5 Moser-ladder record-rung Ã¢â€ â€ Klein 1728) is a coincidence on the integer 1729; ONLY the Moser-ladder 1729 is structural (a record rung *because* split). Validated int64 Held-Karp counter re-verifies H(T_19)=1172695746915 (3.2s) and H(T_23)=15760206976379349 (61.7s).
- **OPEN handles (for the next explorer / a compute node):**
  1. **Analytic-asymptotic** *(2026-06-10 status: the LIMIT Ã¢â€ â€™ e was PROVED next day by HYP-2307; the live question is the RATE Ã¢â‚¬â€ now carried by HYP-2371's falsifiable R(31) prediction and the [COMPUTE-NODE] H(T_31) lead below; handle superseded)*: extend `H(T_p)Ã‚Â·2^{pÃ¢Ë†â€™1}/p!` = 2.00,2.40,2.44,2.53,2.56 with p=31,43,47 to test whether it Ã¢â€ â€™ e, a larger constant, or a slow `~p^{3/2}` polynomial factor (Alon's permanent bound permits the latter Ã¢â‚¬â€ OPEN-Q-013's "Ã¢â€ â€™e" is NOT settled by 5 points). Needs a SYMMETRY-REDUCED counter: my int64 start-fixed Held-Karp is O(2^{pÃ¢Ë†â€™1}Ã‚Â·p) and OOMs past pÃ¢â€°Ë†25 (p=31 Ã¢â€¡â€™ 2^30 rows Ã¢â€°Ë† 266 GB); quotient by the Z_p cyclic action (Ã¢â€°Ë† /p memory & time) for a C/compute node.
  2. **A genuine tournamentÃ¢â€ â€Q(Ã¢Ë†Å¡Ã¢Ë†â€™3) bridge (if one exists):** ~~Is there a Q(Ã¢Ë†Å¡Ã¢Ë†â€™3)/QR-structural reason for Aut acting freely on Ham paths?~~ **CLOSED (kind-pasteur-2026-06-10-S1, LEM-003, adversarially verified):** the handle DISSOLVES Ã¢â‚¬â€ freeness is UNIVERSAL order-rigidity of directed paths, zero Eisenstein content. For ANY digraph, an automorphism fixing a directed Ham path's arc set fixes its vertex sequence pointwise (the unique in-degree-0 source anchors an induction), so Aut acts freely, all orbits have size exactly |Aut|, and |Aut| | H. No Q(Ã¢Ë†Å¡Ã¢Ë†â€™3)/QR structure enters; the only Paley-specific residue is the SIZE |Aut(T_p)|=p(pÃ¢Ë†â€™1)/2, not the integrality of r(p). Proof + exhaustive verification (all 2^10 + 2^15 labeled n=5,6 tournaments; explicit orbits for all 184+3248 masks with |Aut|>1; Burnside cross-checks 120Ã‚Â·12, 720Ã‚Â·56 exact): `LEM-003-aut-acts-freely-on-ham-paths.md`, `aut_divides_H_freeness_kpc1.py/.out` (+ independent re-verification `aut_freeness_kpc1_verify.py/.out`, different machinery, every number matches). Honest boundary: FAILS for Ham CYCLES (C3: 1 cycle, |Aut|=3 Ã¢Ë†Â¤ 1; RQ5: BOTH its Ham cycles rotation-fixed, orbits [1,1], while the path action on the same tournament is free). Prior art reconciled: THM-048 Step 3 asserted it, S20bt's tiling fibration assumed it (n=4,5), opus-S184's "proof" was circular (see MISTAKE-070), THM-212/HYP-640/HYP-1264/HYP-1714 are fixed-point-free special cases. The remaining honest 1729 bridge is the TAXICABÃ¢â‚¬â€œMoser one Ã¢â‚¬â€ see THM-463 (this session; nÃƒÂ© THM-461, renumbered after collision with monad-explorerÃ¢â‚¬â„¢s THM-461): structural, cofactor-level, through the Eisenstein norm form.

## LEAD [COMPUTE-NODE] (kind-pasteur-2026-06-10-S1): run H(T_31) exactly Ã¢â‚¬â€ falsify/confirm HYP-2371
- **Source:** Thread E + adversarial verifier (all 8 claims CONFIRMED). Spec: `05-knowledge/results/paley_H31_compute_design_kpc1.md`.
- **Two independent designs Ã¢â‚¬â€ run BOTH:** (A) Z_31-rotation-reduced layered HeldÃ¢â‚¬â€œKarp: states (S,v) up to rotation (free action proved), layer k Ã¢â€°Â¤ C(30,kÃ¢Ë†â€™1) canonical states, peak C(30,15) = 155117520 states = 2.31 GiB uint128, 8053063680 transitions (= 31Ã‚Â·30Ã‚Â·2^29/62), uint128/CRT from layer ~24, ~0.5Ã¢â‚¬â€œ3 h C/C++; harness logic validated exactly at p=11,19 (`paley_rotdp_smallp_verify_kpc1.py`). (B) Karp inclusionÃ¢â‚¬â€œexclusion over 2^31 subsets grouped into 69273668 rotation classes (Burnside-exact), O(KB) memory, 2Ãƒâ€”63-bit CRT, embarrassingly parallel, ~0.5Ã¢â‚¬â€œ2 h on 8 cores (estimate conservative by 2Ã¢â‚¬â€œ10Ãƒâ€” per verifier).
- **Validation chain:** reproduce H(19)=1172695746915 and H(23)=15760206976379349 first. **Acceptance:** H odd; 465 | H; H/465 odd (Ã¢Å¸Âº H Ã¢â€°Â¡ 465 mod 930); verdict on R(31) Ã¢Ë†Ë† [2.58949, 2.60249] (central 2.59599, H Ã¢â€°Ë† 1.988e25).
- **Afterwards:** update HYP-2371 (CONFIRMED/REFUTED) + OPEN-Q-013; the run pins C in R = e(1Ã¢Ë†â€™C/pÃ¢Ë†â€™Ã¢â‚¬Â¦) to ~3 digits and makes R(43) Ã¢â€°Ë† 2.628 a second-generation prediction. p=43/47 are NOT feasible with these designs (8.6 TB peak layer / 4e15 mults) Ã¢â‚¬â€ needs new mathematics (see T776 Borel-resummation tangent).
- **Status:** OPEN, handed to compute node.

## LEAD (kind-pasteur-2026-06-10-S1): cubic-lens follow-ups (Threads C/D)
- **(a) First forbidden layer (T772):** c3 is gap-free (THM-462); H omits {7,21}. Compute the c4 spectrum from iso classes n Ã¢â€°Â¤ 9 (c4 is NOT score-determined) Ã¢â‚¬â€ where between degree 3 and H does impossibility first appear? Also the (c3,c4) joint spectrum.
- **(b) OEIS submission candidates:** spectrum-size sequence 2,3,6,9,15,21,31,41,56,71,92,113,141 = A006918(nÃ¢Ë†â€™2)+1 (zero OEIS matches 2026-06-10); also the H(T_p) and H/|Aut| sequences (old T052 item) Ã¢â‚¬â€ now worth bundling with the THM-462/THM-463 citations.
- **(c) McShaneÃ¢â‚¬â€œHarris JIS 27 (2024):** verifier REACHED and text-mined the PDF Ã¢â‚¬â€ no spectrum/gap-freeness statement (novelty position strengthened); their per-level generating functions (A357242/48/57/66) vs our gap-free floor + the repo's H-distribution machinery is an open bridge. Moon's *Topics on Tournaments* cyclic-triples chapter still worth a skim for a folklore interval exercise before any external novelty claim.
- **(d) Three-cubes ledger follow-ups (T773Ã¢â‚¬â€œT775):** rigidity-pruned k=114 search with honest exclusion bound (NOT expecting a hit; min coordinate > ~10^16 per BookerÃ¢â‚¬â€œSutherland); OnoÃ¢â‚¬â€œTrebat-Leder full text when web is stable; Hirschhorn closed-form vs our two recurrence-theoretic proofs; primitivity status of 192/375/600.
- **Status:** OPEN, prioritized (a) > (c) > (b) > (d).

## LEAD (codex-2026-06-28): instantiate HYP-3300 observability rows
- **Source:** HYP-3300, HYP-3258, HYP-3257, HYP-3255, HYP-3253, HYP-3250, HYP-3249, HYP-3248, HYP-3247, HYP-3246, HYP-3245, HYP-3244, HYP-3243, HYP-3236, HYP-3225,
  HYP-3202, HYP-2963.
- **Status:** OPEN.  HYP-3300 currently uses a toy residual-pair matrix with
  full `GF(2)` row rank, not actual packet rows.
- **Next steps:** Build a row dataset whose rows are real residual-pair
  collisions from HYP-3202 trap neighborhoods, HYP-3225 trap sidecar classes,
  HYP-3236 Green decoys, HYP-3243 chamber exits, HYP-3244 lift/compress
  failures, and HYP-3245 lag-transport residuals.  Test column subsets
  including endpoint owner, boundary cocircuit, Phi witness address, dilation
  grid, Toeplitz slack, Green resistance, odd-negative payload, Lee-Yang root
  word, tiling descent packet, lag transport, unit-equioscillation index,
  binding complement-pair word, analytic/topological index equalizer,
  Gauss-sum index word, Borsuk-Ulam forcing-gap flag, Roth-Halasz discrepancy,
  Hensel-Krasner unit, state-lift H7, Cech hole, and ear payload.  Every
  unresolved collision should become either a same-status equivalence, a
  missing sidecar, or a named proof obstruction.  Then attempt the first
  acyclic chamber matching using the same retained columns as Morse energy
  coordinates.

## LEAD (codex-2026-06-28): instantiate HYP-3301 sheaf/cusp rows
- **Source:** HYP-3301, HYP-3300, HYP-3265, HYP-3257, HYP-3255, HYP-3253,
  HYP-3247, HYP-3246, HYP-3243, HYP-3242, HYP-3234, HYP-3231, HYP-3230,
  HYP-3102, HYP-2969, HYP-2963, HYP-2954, HYP-2704, THM-573, THM-523.
- **Status:** OPEN.  HYP-3301 currently uses a toy chart-overlap exactness
  packet with full `GF(2)` row rank, not actual labelled packet rows.
- **Next steps:** Build row data for the overlap failures
  `unit_contact_to_danger_nerve`, `lag_residue_to_shell_magic`,
  `index_degree_to_floor_gap`, `residue_to_magnitude_split`,
  `chart_change_to_first_obstruction`, `exact_period_to_boundary_moment`,
  `covering_to_state_lift_kernel`, `farey_cusp_to_scale_normal`,
  `polygon_contact_to_endpoint_owner`, and
  `renormalized_kernel_to_named_debt`.  Test sidecars
  `first_obstruction_cocycle`, `zeta7_contact_holonomy`,
  `endpoint_arrangement_cell`, `exact_period_boundary`,
  `boundary_moment_image`, `Farey_parent_interval`, `cusp_principal_part`,
  `renormalization_depth`, `AP_GW_kernel_status`, `K33_kernel_label`, and
  `state_lift_H7`.  The target classification is exact, holonomy-repaired,
  endpoint-lifted, descended, AP/GW-stopped, positive boundary-moment image,
  qdiv-forbidden AP/GW kernel, K33/H7 debt, or a newly named zero-open kernel.

## LEAD (codex-2026-06-28): owner-current and tropical-wall first-leak table
- **Source:** HYP-3402, HYP-3311, HYP-3400, HYP-3310, HYP-3301, HYP-3260,
  HYP-3265, HYP-2969, HYP-2963.
- **Status:** OPEN.  HYP-3402 is a proof-angle scout over the HYP-3311 mixed
  fiber; it is not yet an enlarged actual-packet theorem table.
- **Next steps:** Extend the HYP-3311 actual-packet bank and add
  `owner_current_word` plus `tropical_wall_word`.  Classify every mixed or
  formerly mixed fiber as: residue exact; residue fails but owner-current
  works; residue/v2 fail but tropical wall works; both fail and emit named
  owner/height/off-grid debt.  The target owner-current exits are conserved
  current, Farkas/Green dual certificate, AP/GW boundary H1, forbidden H7
  state-lift sink, or named current debt.  The target tropical exits are
  positive off-grid floor wall, AP/GW `12->24` hinge, or named
  height-discriminant debt.

## LEAD (codex-2026-06-28): boundary-uniformization Menger cut table
- **Source:** HYP-3407, HYP-3420, HYP-3412, HYP-3410, HYP-3408, HYP-3406, HYP-3405, HYP-3404, HYP-3402, HYP-3124,
  HYP-3123, HYP-2982, HYP-2963.
- **Status:** OPEN.  HYP-3407 ranks the boundary-uniformization Menger zipper
  first, but the owner-support graph has not yet been computed as a real
  min-cut object.
- **Next steps:** Build the endpoint-owner support graph for HYP-3406's
  `petal 13->26` versus positive-open single-swap `26/40/54` families and
  compute min-cuts separating `unit-petal-named` from `positive-Haar-open`.
  In parallel, build the HYP-3405 AP versus `13->27` local disk / unit-height
  exit table.  Then add recursive chiral child decks after deleting active
  endpoint owners.  If these sidecars stay stable under bank enlargement, ask
  for a BDH/Mertens mean-square theorem with the known exceptional fibers
  removed.

## LEAD (codex-2026-06-29): prove the random031 seven-owner gate-gluing lemma
- **Source:** HYP-3455, HYP-3454, HYP-3453, HYP-3452, HYP-3451, HYP-3450, HYP-3439, HYP-3438,
  HYP-3437, HYP-3436, HYP-3435, HYP-3434, HYP-3431, HYP-3429, HYP-3427,
  HYP-3425, HYP-3422, HYP-3418, HYP-3415.
- **Status:** OPEN.  HYP-3455 is an exact single-row obstruction isolator, not
  yet a theorem.
- **Readout:** `random_covering_031` is the named noncanonical rank-`6`
  overlap-rescue exception.  Its rescue subset is `(23,45,93,113,147,169)`,
  the rescue graph is connected with `15` edges, and the HYP-3438 gate ledger
  has `84` mixed components and `138` survivor gates.  The only max-delta
  gates are mirror components `43` and `54`, both word `B-S-B-S-B`, with
  deltas `(3,4)` and `(4,3)`.  Their owner union is
  `(23,45,93,113,147,169,173)`, containing all rescue owners and only extra
  owner `173`.  HYP-3450 side data gives `94` low-rank escapes against only
  `4` dead components.
- **Next steps:** Prove that the two max-delta mirror gates cannot glue into a
  full two-colour cover without exposing a low-rank component escape,
  endpoint-spine/wall certificate, owner-current imbalance, two-adic descent,
  or signed-SPEC/Rprime debt.  Test a Menger terminal cut on components
  `43/54`, then try the Schwarz-Christoffel boundary-word formulation with
  owner `173` as the accessory-owner payload.

## LEAD (codex-2026-06-29): prove bounded overlap-tax cut cores
- **Source:** HYP-3437, HYP-3436, HYP-3435, HYP-3434, HYP-3431, HYP-3429,
  HYP-3427, HYP-3425, HYP-3422, HYP-3418, HYP-3415.
- **Status:** OPEN.  HYP-3437 is an exact finite stress audit, not yet a
  theorem.
- **Readout:** On the `150`-row HYP-3429 bank, the one-branch overlap-tax
  atomization has branch0 positive `150/150`; the `59` negative naive-slack
  rows all have strict rescue cuts; rank histogram is
  `{0:91,2:7,4:2,5:48,6:2}`.  The rank-`6` stress is only the canonical
  AP-with-84 duplicate, with core `(3,5,7,9,11,13)` and margin `563/105105`.
- **Next steps:** Prove the rank-`6` canonical core is exactly HYP-3431's
  corridor-fence family.  Then prove noncanonical negative-slack rows have a
  rank-`<=5` endpoint-labelled odd-blocker overlap core, or find and name the
  first exception.  Attach HYP-3429/HYP-3427 endpoint/wall labels to the cut
  core so HYP-3436's two-color bad-core extractor can use the same sidecar.
## LEAD (codex-2026-06-28): endpoint-cut resurrection from Charal signatures
- **Source:** HYP-3429, HYP-3428, HYP-3427, HYP-3426, HYP-3425, HYP-3424, HYP-3423, HYP-3422, HYP-3421, HYP-3420, HYP-3419, HYP-3418, HYP-3417, HYP-3416, HYP-3415, HYP-3414, HYP-3413, HYP-3412, HYP-3411, HYP-3410, HYP-3409, HYP-3408, HYP-3407, HYP-3406,
  HYP-3405, HYP-3404, HYP-3402, HYP-3311, HYP-2969, HYP-2963.
- **Status:** OPEN.  HYP-3429 executes the requested Bring/SC/Menger
  endpoint-cut slice after HYP-3412's special-function sidecar scout and
  HYP-3414's owner-cut resurrection calculus, now integrated with HYP-3416's
  recursive quotient ladder, HYP-3417's owner-current certificates,
  HYP-3419's size-3/depth-3 owner-cut decision-tree warning, and HYP-3420's
  chiral-owner exactness audit.  HYP-3421/HYP-3422 now give the floor-facing
  off-grid transparency and two-adic relocation lane, HYP-3423 gives the
  topology/arithmetic guardrail, HYP-3424 gives the transfer ledger, and the
  paired HYP-3425 threads give the Helly interval target plus the
  energy-plus-sheet packet guard, while HYP-3426/HYP-3427 sharpen the endpoint
  payload into one-branch endpoint-support triples and exact wall words.
  HYP-3428 supplies the two-adic loss ledger that any halving shortcut must
  feed instead of bypassing.  HYP-3429 sits under that caveat: this lead
  is a certificate sidecar over the HYP-3406 mixed fibers, not yet a theorem,
  not a q-uniform topology shortcut, and not a raw energy scalar route to the
  q-specific floor.
- **Next steps:** Build an explicit endpoint-owner deletion graph for every
  residue-plus-height mixed fiber.  Compute min cuts, owner-support paths,
  margin-1 owner currents, and the corresponding Schwarz-Christoffel boundary
  word.  Prove or refute: bounded endpoint cuts reconstruct
  theorem-exit-pure cut-code buckets whenever residue+height leaks, then feed
  that sidecar into HYP-3424/HYP-3422/HYP-3425/HYP-3426/HYP-3427 only as owner/branch/packet/wall debt
  payload.  Use HYP-3423 to reject any proof quotient that converts SC/Menger boundary
  order into a magnitude conclusion without restoring arithmetic, owner-current,
  or two-adic floor data.
  If a cut fails, name the owner-cut debt before invoking Bring branch
  addresses, Krasner stability, HLW no-free-slider guards, BDH variance,
  Sophie Germain quartic factors, or Ramanujan-Soldner/Meissel-Mertens bulk
  constants.

## LEAD (codex-2026-06-29): AP84 survivor-gate color-legality matrix

- **Source:** HYP-3459, HYP-3458, HYP-3457, HYP-3456, HYP-3454, HYP-3453, HYP-3438,
  HYP-3441, HYP-2991, HYP-2263, HYP-2247.
- **Status:** OPEN.  HYP-3459 shows that the AP84 gate palette and
  floor-correction word are distinct colors, and endpoint phase is not
  determined by residue.
- **Readout:** `both_outer_inner` residues split across correction values
  `[1,2]`; correction `d=1` appears in three gate buckets; the pair
  `(gate_bucket,d)` is clean locally; but `m=1` and `m=36` share residue while
  living in different endpoint phases.
- **Next steps:** Build a matrix with HYP-3438 survivor gates/components as
  rows and color columns
  `(gate, floor, height, endpoint_phase, branch_mirror, incident_C3_Qsqrt,
  zipper_cocycle)`.  Prove the AP84 local-to-global splice into HYP-3439 is a
  homomorphism for this color-packet product, or identify the first failed
  color and route it to HYP-3453/HYP-3455, owner-current, two-adic descent,
  exact-period/state-lift, or SPEC debt.  Do not replace this with a scalar
  color balance claim.

## LEAD (codex-2026-06-29): colored-gate Nerode quotient theorem

- **Source:** HYP-3474, HYP-3473, HYP-3472, HYP-3471, HYP-3462, HYP-3470, HYP-3461, HYP-3460,
  HYP-3459, HYP-3458, HYP-3455, HYP-3453, HYP-3451, HYP-3438.
- **Status:** OPEN.  HYP-3474 is a finite partition-lattice scout, not a
  theorem.
- **Readout:** On the `135`-row HYP-3471 bank, the theorem-failure bit
  `dead and not has_e_branch` is empty, but route data is mixed by several
  attractive quotients: endpoint-kind `K` mixes `132` rows, structural set `S`
  mixes `32`, minimum E/branch structural gate `M` mixes `124`, and AP bit `A`
  mixes `135`.  Route flags are pure for `C` (`122` fibers) and `N/T/F`
  (`123` fibers), while the dead minimum gate family needs `(C,M)`, `(N,M)`,
  `(T,M)`, or `(F,M)` (`125` fibers, max fiber `7`).
- **Next steps:** Extend the audit to fresh banks and prove or refute the
  count-profile reconstruction lead.  The theorem target is: after
  `dead_components(row)>0 => rank<=2 E/branch survivor gate`, every quotient
  used by the proof must refine the relevant route partition or attach a
  sidecar reconstructing the forgotten coordinate.  First mixed fibers should
  be routed to AP84 packet closure, same/cross-branch gluing, endpoint-owner
  current, two-adic floor descent, or signed-SPEC/Rprime debt.
