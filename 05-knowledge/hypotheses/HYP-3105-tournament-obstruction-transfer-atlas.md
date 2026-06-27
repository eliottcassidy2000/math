---
id: HYP-3105
title: Tournament obstruction-transfer atlas
status: SYNTHESIS / technique atlas; not a proof
source: codex-2026-06-27-S259b
tangent: T1182
technique: LTI-243
tournament_technique: LTT-141
script: 04-computation/tournament_obstruction_transfer_atlas_codex_s259.py
result: 05-knowledge/results/tournament_obstruction_transfer_atlas_codex_s259.out
related:
  - HYP-3104
  - HYP-3103
  - HYP-3102
  - HYP-3101
  - HYP-3100
  - THM-002
  - THM-029
  - THM-079
  - THM-115
  - THM-264
  - THM-454
  - THM-027
  - THM-030
  - THM-255
  - THM-577
  - HYP-3099
  - THM-576
  - THM-575
  - THM-573
  - HYP-3098
  - HYP-3097
  - HYP-3096
  - HYP-3095
  - HYP-3094
  - HYP-3078
  - HYP-3076
  - HYP-3074
  - HYP-2963
  - OPEN-Q-108
---

# HYP-3105: Tournament Obstruction-Transfer Atlas

## Claim

The old H=7 and H=21 contradiction proofs should be generalized as an
obstruction-transfer method:

```text
constructed subproblem
  -> faithful tournament/OCF/conflict carrier
  -> forbidden local spectrum or forced-expansion skeleton
  -> contradiction, named residual debt, or a sharper proof obligation
```

The method is proof-facing only when the transfer functor states exactly what
LRC predicate it preserves, which coordinates it destroys, and which sidecar
repairs each destroyed coordinate.  A raw statement like "this subproblem
looks like H=7" is a guardrail failure.

The incoming tournament-contradiction grammar and Lean observer-gluing
frontier make this stricter.  HYP-3100 supplies the legality chain for
certificate pullbacks; the theorem target is `ObserverGluingCoverage`, and
`coarseWinding_degenerate_not_terminal` says coarse mod-14 winding H is not a
terminal proof carrier.  Therefore H-values can be reused only as local
forbidden spectra inside a certified observer-gluing or packet-transfer
ledger.

The latest incoming proof angles make the first targets sharper: HYP-3101 asks
for bounded normalized direct components via normal-fan, Cech, and barcode
packets; HYP-3102 asks that hidden observer-chart payloads be generated,
annihilated, descended, AP/GW-stopped, or routed to F7/THM-572; HYP-3103's
perspective-groupoid layer asks whether a quotient functor has kept the
sidecar needed by its next operation; and the incoming HYP-3104 maximizer
signal atlas separates LRC value currencies from tournament H-rigidity.  This
makes H a transfer-risk signal, not a scalar maximizer certificate.  The
obstruction-transfer atlas should be applied to those packets first, not to an
untyped scalar H shadow.

## Core Moves

1. `permanent_H_gap_transfer`: transfer a subproblem to an OCF carrier whose
   independence polynomial at 2 would equal a forbidden value, then import the
   H-gap theorem.  Required payload: complete transfer to the OCF object,
   connected-component status, and retained LRC predicate.
2. `component_value_factorization`: split a proposed obstruction by component
   values.  A 7-factor inside an H=21-looking object is illegal by THM-029 and
   THM-079 style factor blocking.
3. `forced_expansion_closure`: assume the attractive small skeleton, then
   prove the real carrier forces extra cycles, sidecars, or overlap payloads
   that jump the invariant out of the forbidden value.  This is the most
   reusable move for LRC14.
4. `deletion_invariant_induction`: show that deleting a low-load coordinate
   preserves the obstruction carrier, reducing to an already-blocked base.
5. `typed_ocf_evaluation_lift`: replace scalar H by a typed vector of
   evaluations, such as component H, cap inclusion-exclusion order, chart
   overlap type, and terminal-exit class.

The generated atlas also tested ranked proof-carrier tournaments, cycle-class
observability matrices, median-center legality tournaments, edge-flip stress
disprovers, hypertournament blocker allocation, quotient round trips, tropical
capacity gaps, valuation-residue lift audits, and semantic tournaments of
statements.

After rebasing over HYP-3099/S31ah/S65 and the KPS single-component ladder
update, the atlas imports seven sharper tournament proof-engine levers:

1. `redei_parity_certificate`: a count encoded as tournament Hamiltonian
   paths must be odd.
2. `landau_score_feasibility`: an encoded win profile must satisfy Landau's
   score inequalities.
3. `cycle_census_hole_transfer`: a missing `(c3,c5,...)` fiber is a forbidden
   spectrum distinct from forbidden H.
4. `score_exchange_nontransitivity`: an improvement tournament can prove a
   greedy exchange only when local minima collapse to the global sink; otherwise
   it exposes a finite check.
5. `winding_tie_apex_audit`: the LRC14 apex 7 is an antipodal matching count,
   not automatically an `Omega=K3` / H=7 event.
6. `certificate_engine_spectrum_generator`: enumerate tournament invariants
   and promote persistent missing values to candidate certificates instead of
   handpicking H=7/H=21.
7. `single_component_H_ladder_certificate`: identify H=7 and H=21 as exactly
   the single-component clique-Omega gaps `K3` and `K10`, with small
   Omega-realizability sparse and clique-dominated.

## Programmatic Readout

The script records the H-gap skeletons:

```text
K3:                         I(G,2)=7
P4:                         I(G,2)=21
K6 minus two disjoint edges: I(G,2)=21
K8 minus one edge:           I(G,2)=21
K10:                        I(G,2)=21
```

These are not proposed LRC carriers by themselves.  They are regression tests
for transfer discipline: if a future proof imports one of these spectra, it
must also show why the true packet cannot force the expansion mechanisms that
blocked H=7 and H=21 in THM-029, THM-079, and THM-115.

Top atlas applications after integrating the incoming certificate-engine work:

```text
single-component H ladder @ KPS spectrum discovery       score 28
certificate-engine spectrum generator @ S31ah engine     score 27
winding tie/apex audit @ apex7-forbidden-H bridge        score 26
score-exchange nontransitivity @ S65 cap optimality      score 26
cycle-census hole transfer @ baby-Hodge cycle hole       score 25
ranked proof-carrier tournament @ LRC14 observer gluing  score 24
median-center legality tournament @ route-state hull     score 23
forced-expansion closure @ THM-577 cap dip               score 20
edge-flip stress disprover @ HYP-2963 packet bank        score 20
```

## Application Notes

For LRC14 observer gluing, vertices should be chart-overlap obligations, not
runners, residues, denominators, raw arcs, or Pascal numbers.  The useful
contradiction shape is:

```text
assume a residual packet has no normalized arc exit,
no moment/Perron exit,
no nested-refinement covering exit,
and no K33/state-lift exit;
then its retained chart-overlap tournament must realize an impossible local
obstruction-transfer ledger.
```

For THM-577, the tournament-like object is the inclusion-exclusion order
vector.  The cap value alone is too coarse.  The forced-expansion question is
whether every would-be dip-free skeleton at `j=4,5` forces a higher-order
overlap term that either becomes the known finite dip or is discharged by a
moment/Perron certificate.

S65 adds a second cap-facing tournament.  The single-swap improvement
tournament is useful, but not as a clean exchange proof: it finds bounded
minimizers through `{1..20}` and greedy descent works for `j=2,3,4`, while
`j=5` gets stuck at `{1,10,11,12,13}` instead of the true minimizer
`{1,5,7,8,9}`.  The result is progress because it turns the obstruction into
finite local-minimum elimination, not because it proves optimality by
transitivity.

For HYP-3094 K33 state lift, positivity is not the right edge direction.
Nested-refinement covering rows and K33 cross-handoff rows are both
positive-open.  The tournament should compare which sidecar gives a legal
terminal exit: branch binders, endpoint-owner transitions, median centers, or
the first named state-lift debt.

For HYP-2963, the atlas suggests an automatic stale-quotient rejection tool.
Given a new scalar proposal, build an edge-flip tournament over packet rows:
if small sidecar perturbations reverse the scalar order while preserving the
LRC predicate, the scalar is only telemetry.

For q-cusp, sixth-power, and p-adic/Roth/Steiner analogies, tournament use is
mostly disproof hygiene.  The question is not whether the analogy has the
right flavor, but whether finite principal parts, support-six lane rank,
valuation fibers, discrepancy bins, and metric payloads survive the quotient.

For the S31ah certificate engine, the key refinement is generative: run the
same spectrum-gap search over H, Hamiltonian-path parity, score sequences,
cycle-count fibers, Omega shapes, Newton/real-rootedness shadows, transfer
symmetry, tournament spectra, H-maximality, and winding encodings.  This is
how the project should look for the next H=7/H=21-like obstruction.

The latest KPS spectrum-discovery pass validates the generator on the
original frontier.  It mechanically rediscovers H=7 and H=21 as persistent
gaps, rediscovers the known `Omega` non-realizability facts, and reframes
`K_m = Omega(T)` as the clique ladder `H=1+2m`: the missing clique sizes are
`m=3` and `m=10`, corresponding to H=7 and H=21.  That makes the usable
obstruction more precise: transfer to a single connected `Omega` component
whose clique-like conflict shape would require `K3` or `K10`, not merely to a
loose scalar H.

Post-rebase KPS artifacts make this a concrete regression test rather than a
metaphor.  `tournament_I21_omega_miner_kps.py` lists the five connected
`I(G,2)=21` targets `K10`, `K8-e`, `K6-M`, `K6-P3`, and `P4`, and routes their
non-realizability through the single-component H gap plus forced expansion.
`tournament_certify_applications_kps.py` upgrades the certificate battery to
the spectrum statement `#HamPaths(T) in {odd >= 1} \ {7,21}`, with Redei
parity, power-of-two collapse to transitive tournaments, strong-connectivity
lower bounds, and the `n=7` regular values as reusable guards.
The KPS close-out messages call this clean-spectrum package `HYP-3101`, but
this navigation already uses HYP-3101 for the normal-fan Cech component-bound
route; this atlas therefore cites the KPS certificate package by script/result
name and by the S31ah/KPS label.

For the apex-7 bridge audit, literal forbidden-H transfer fails.  Winding
tournaments avoid H=7 vacuously because every tournament does, and the real
apex event at denominator 14 is a perfect matching of seven antipodal diameter
pairs.  That matching is triangle-free, so it is structurally opposite to an
`Omega=K3` component.  Any future apex-to-H bridge needs a new functor.

For the baby-Hodge crosscheck, the `(c3,c5)=(8,10)` hole is not the H=7/H=21
alpha obstruction.  Neighboring fibers carry realized H-values such as `41`
and `43`; the hole is a `c5`/power-sum spectral exclusion inside the regular
score class.  This validates cycle-census holes as their own obstruction
layer.

## Ledger Schema

New obstruction-transfer attempts should record:

```text
surrogate_vertex_set
transfer_functor
preserved_lrc_predicate
target_H_or_indpoly_value
forbidden_spectrum_source
minimal_skeleton
forced_expansion_payload
component_factorization
deletion_inert_coordinate
typed_ocf_vector
certificate_invariant_family
single_component_H_gap
clique_omega_realizability
omega_sparsity
score_sequence_or_Landau_status
cycle_count_fiber
improvement_tournament_local_minima
apex_tie_matching_status
destroyed_coordinate
required_sidecar
edge_flip_stress_result
terminal_exit_or_named_debt
verdict
```

The challenged assumption is that tournament analogy is proof currency.  The
replacement principle is that tournament construction is a functorial
obstruction language: it can prove, disprove, rank, or reject only after its
payload map has been audited.
