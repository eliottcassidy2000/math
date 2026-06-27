---
id: HYP-3085
title: LRC14 covering-moment and K33 state-lift shuttle
status: EVIDENCE / two-obligation proof-interface shuttle; not a proof
source: codex-2026-06-27-S254
technique: LTI-232
tournament_technique: LTT-130
script: 04-computation/lrc14_o2_o3_shuttle_codex_s254.py
result: 05-knowledge/results/lrc14_o2_o3_shuttle_codex_s254.out
related:
  - HYP-3083
  - HYP-3084
  - HYP-3082
  - HYP-2965
  - HYP-2969
  - HYP-2963
  - HYP-3037
  - HYP-3045
  - THM-573
  - THM-572
  - OPEN-Q-108
---

# HYP-3085: LRC14 Covering-Moment and K33 State-Lift Shuttle

This hypothesis records a focused back-and-forth between the two live
obligations from HYP-3083:

```text
O2 = covering-moment / OPEN-Q-108 gamma-Node3 positive-open discharge
O3 = K33 / THM-572 state-lift construction
```

The claim is not that either obligation is closed.  The claim is that each
obligation now tells us what kind of data the other must keep.

This is rebased over the incoming `codex-2026-06-27-S253` expansion of
HYP-3083.  In that spine, O2 is the covering-family gamma/Node3 terminal
discharge for low-apex, top-balanced finite-address packets, not a scalar
safe-mass shortcut.

It is also rebased over the incoming `kind-pasteur-2026-06-27-S31af`
THM-573 level-7 lift sieve.  THM-573 closes every row with at least seven
multiples of `7`, so the O2 residual is now a smaller `<=6`-multiples-of-7
core.  Its covering-margin stress tests also refute the tempting scalar idea
that all covering rows have a uniform margin after dilation/aliasing; S254's
nested-refinement rows should therefore be read as local interface evidence,
not as a global margin theorem.

The incoming mac-mini S60 gK8 artifact
`05-knowledge/hypotheses/HYP-3084-gk8-concentration-is-low-order-moment-s2-clebsch.md`
is the complementary O2 certificate candidate: it localizes the gK8/p0
covering-moment bound to low-order factorial moments led by pairwise `S2`,
then corrects the literal Clebsch-design guess to a reflection-symmetric
Perron-mode bound.  Read together, that S60 artifact supplies a possible
moment-LP carrier for O2, while HYP-3085 records the finite-address guardrails
that prevent that carrier from collapsing into a scalar-margin claim or
swallowing K33 cross-handoff debt.

## Exact Scout

Script:

```text
04-computation/lrc14_o2_o3_shuttle_codex_s254.py
```

Stored output:

```text
05-knowledge/results/lrc14_o2_o3_shuttle_codex_s254.out
```

The scout computes exact `M`, optimizer `tau`, q-threshold, strict safe mass,
safe-component count, active binders at `tau`, endpoint-owner transitions, and
grid class on three covering-moment representatives and the three low-frontier
K33 representatives:

```text
cover tail 12->84                 O2  M=7/89   tau=37/89  safe_mu=563/105105  grid=nested_refinement
cover tail 6->98                  O2  M=9/109  tau=19/109  safe_mu=1543/294294 grid=nested_refinement
cover tail 12->168                O2  M=14/173 tau=72/173  safe_mu=263/30030   grid=nested_refinement
near/K33 12->36                   O3  M=3/41   tau=17/41   safe_mu=1/1260      grid=cross_handoff
P10+K33                           O3  M=2/27   tau=8/27    safe_mu=4/2205      grid=cross_handoff
two drop(12,13)->add(26,36)       O3  M=3/37   tau=3/37    safe_mu=79/8190     grid=cross_handoff
```

The active binder words show the split the proof has to preserve:

```text
cover tails: bind a non-14 core owner against a 14-multiple owner
K33 rows:   bind a nonunit cross-handoff owner such as 36, 20, or 26/36
```

More explicitly, the three K33 binder supports are:

```text
near/K33 12->36              (5,36)
P10+K33                      (7,20)
two drop(12,13)->add(26,36)  (1,36)
```

## Back-And-Forth Readout

### O2 progress

The selected covering rows all have exact positive strict-safe mass and
nested-refinement exits.  This supports the HYP-2965/HYP-2969 view that the
covering-moment family is locally positive; the hard part is making the
positive-open certificate uniform over the family.

### O2 failure inspiring O3

Raw positivity is not a theorem separator.  The K33 rows are also positive
open.  Therefore O2 cannot be stated as "positive safe mass exists" without
retaining grid class and owner transitions.  Covering discharge must be a
family theorem for nested-refinement / owner-transition packets, with
cross-handoff packets routed out to O3.

### O3 progress

The K33 packets are localized: in the low-frontier rows they are exactly
cross-handoff packets with small active binder supports.  This gives a more
concrete input target for THM-572:

```text
cross_handoff packet + active binder owner word + endpoint transition word
  -> TournamentStateLift
```

### O3 failure inspiring O2

The script still does not build the THM-572 lift.  If that lift cannot be
constructed from the current cross-handoff fields, then the missing field is
not a mysterious new scalar.  It should be added as a localized
covering-moment sidecar: owner transition, active-binder support,
endpoint-current page, or cycle-class row.

## Theorem-Facing Update

The two obligations should now be stated as a coupled pair:

```text
O2 theorem target:
  every low-apex, top-balanced, <=6-multiples-of-7 covering
  nested-refinement packet has a uniform
  gamma-Node3 / positive-open / moment / finite-ruler certificate,
  unless it is rerouted to cross_handoff debt.

O3 theorem target:
  every K33 cross_handoff packet with its active-binder and endpoint-owner
  address constructs the TournamentStateLift needed by THM-572, or else
  names the first missing localized sidecar.
```

This is a stricter interface than the old route split.  It forbids using safe
mass as a scalar shortcut and forbids treating K33 as solved just because it
is positive-open.

## Tournament Analysis

Vertices are proof carriers / finite-address fields, not runners:

```text
finite_address_packet
covering_nested_refinement
k33_cross_handoff_lift
endpoint_owner_transition
exact_M_and_binders
safe_mass_scalar
raw_route_label
raw_runner_set
```

Pairwise observable: retained exact scale, safe topology, endpoint owners,
grid class, terminal exit, and named debt.

Switch/gauge: prefer the carrier that retains more address data for
`M(S)>=1/14`; use the displayed finite-address path to break incomparability
between covering and K33 carriers.

Stored fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1,1]
hamiltonian_path_count=1
tie_path=finite_address_packet > covering_nested_refinement
  > k33_cross_handoff_lift > endpoint_owner_transition
  > exact_M_and_binders > safe_mass_scalar > raw_route_label
  > raw_runner_set
```

## Assumption Challenge

The scout considered runners, residues, safe components, fixed circle
sections, section boundaries, endpoint owners, cover arcs, Fourier/moment
modes, grid classes, and proof obligations.  The selected vertices are proof
obligations because they preserve the predicate `M(S)>=1/14` together with the
route needed to prove it.  The quotient destroys raw runner identity and raw
safe-mass scalar data; that destruction is legal only after the row is routed
to nested-refinement discharge, cross-handoff lift debt, or a named missing
sidecar.

## Next Pull

Turn this shuttle into a family ledger:

```text
for each low-apex/top-balanced residual packet:
  if grid_class = nested_refinement:
      attach owner-transition moment certificate or finite-ruler bridge
  if grid_class = cross_handoff:
      attempt THM-572 lift from binder/owner transition data
  if either fails:
      record the first missing sidecar rather than merging it into a scalar
```

The immediate concrete next step is to run the same shuttle fields on the full
HYP-2963 packet bank cache, not just the six representative rows.
