---
id: HYP-3479
title: LRC14 hard mirror-orbit / boundary-current join
status: EVIDENCE / finite ledger join; not an LRC14 proof
source: codex-2026-06-29 continuation of HYP-3476 exception-frontier routing, HYP-3477 hard-orbit discharge, HYP-3475 hard mirror-orbit debt, and HYP-3472 boundary-current cuts
tangent: T1439
technique: LTI-439
tournament_technique: LTT-339
script: 04-computation/lrc14_hard_orbit_current_join_codex_20260629.py
result: 05-knowledge/results/lrc14_hard_orbit_current_join_codex_20260629.out
lean_module: 04-computation/lean/TournamentH7/TournamentH7/LRCHardOrbitCurrentJoin.lean
lean_result: 05-knowledge/results/lrc14_hard_orbit_current_join_lean_codex_20260629.out
reflection: 07-reflections/lrc14-hard-orbit-current-join-codex-20260629.md
related:
  - HYP-3479
  - HYP-3477
  - HYP-3476
  - HYP-3475
  - HYP-3474
  - HYP-3473
  - HYP-3472
  - HYP-3471
  - HYP-3470
  - HYP-3462
  - HYP-3461
  - HYP-3460
  - HYP-3459
  - HYP-3458
  - HYP-3455
  - HYP-3453
  - HYP-3451
  - HYP-3438
  - HYP-3436
  - THM-523
  - OPEN-Q-108
---

# HYP-3479: LRC14 Hard Mirror-Orbit / Boundary-Current Join

## Claim

HYP-3475 isolated the hard colored mirror-orbit debt:

```text
8 mirror orbits with cover-delta >= 7
on 7 random rows
```

HYP-3472 isolated the rows where low-rank E/branch gates fail to become
projection-edge or separating boundary currents.  HYP-3476 then identified
the singleton hard/currentless overlap `random_covering_031`, and HYP-3477
audited the hard-orbit discharge split.  HYP-3479 should be read as the
Lean-backed ledger join for that shared split, not as a competing namespace.
The key finite-bank result is:

```text
hard mirror-orbit debt and current-cut failure intersect in exactly one row:
random_covering_031.
```

Thus the hard-orbit discharge target is no longer an independent seven-row
problem.  On the current bank it reduces to:

```text
hard_orbit_discharge
  <= separating_current_transfer + random031_named_clause.
```

## Exact Readout

Script:

```text
04-computation/lrc14_hard_orbit_current_join_codex_20260629.py
```

Stored output:

```text
05-knowledge/results/lrc14_hard_orbit_current_join_codex_20260629.out
```

Lean ledger:

```text
04-computation/lean/TournamentH7/TournamentH7/LRCHardOrbitCurrentJoin.lean
05-knowledge/results/lrc14_hard_orbit_current_join_lean_codex_20260629.out
```

The script parses the stored HYP-3475 hard-orbit ledger and recomputes
HYP-3472 row-local currents on the hard rows plus the named exception rows.
The recheck fields are clean:

```text
edge_exception_recheck_failures=[]
separating_exception_recheck_failures=[]
hard_edge_nonexception_recheck_failures=[]
hard_separating_nonexception_recheck_failures=[]
```

The Lean module records the same dispatch arithmetic as a sorry-free finite
ledger.  The stored build output reports:

```text
hyp3479_hard_orbit_class_partition: no axioms
hyp3479_separating_current_orbit_partition: no axioms
hyp3479_dispatch_complete: no axioms
hyp3479_dispatch_matches_join_counts: no axioms
```

Joined ledger:

```text
hard_orbits_delta_ge_7=8
hard_rows=[
  random_covering_022,
  random_covering_031,
  random_covering_049,
  random_covering_078,
  random_covering_080,
  random_covering_085,
  random_covering_113
]
hard_orbit_class_hist={cross_branch:2, same_branch:6}

hard_rows_with_projection_edge_cut=6/7
hard_rows_without_projection_edge_cut=[random_covering_031]
hard_rows_with_separating_current=6/7
hard_rows_without_separating_current=[random_covering_031]
hard_orbits_with_separating_current=7/8
hard_orbits_without_separating_current=[
  (random_covering_031, (43,54))
]
ap84_hard_rows=[]
terminal_by_hard_orbit={
  boundary_current_transfer: 7,
  random031_gluing_or_phase_branch_debt: 1
}
```

The current-exception rows that are not hard-orbit debt are:

```text
edge_exceptions_without_hard_orbit=[
  random_covering_001,
  random_covering_039,
  random_covering_062,
  random_covering_074,
  random_covering_086,
  random_covering_101
]
separating_exceptions_without_hard_orbit=[
  covering_AP_with_84,
  ap_omit_12_tail_84x01,
  random_covering_001,
  random_covering_039,
  random_covering_062,
  random_covering_074,
  random_covering_086,
  random_covering_101
]
```

## Proof Pull

HYP-3479 splits two debts that previously looked entangled:

1. Hard mirror-orbit debt is almost entirely current-transfer debt.  Seven of
   the eight hard orbits already sit on rows with separating E/branch currents.
2. The only hard orbit without a projection cut or separating current is the
   HYP-3455 row `random_covering_031`, with hard components `(43,54)`.
3. AP84 has no hard orbit debt.  Its two nonseparating base rows remain AP
   sidecar debt through HYP-3462/HYP-3470 rather than mirror-orbit debt.
4. The six remaining random current exceptions are low-delta current
   exceptions, not hard-orbit exceptions.  They should be routed through
   owner-current, two-adic, signed-SPEC, or state-lift sidecars.

Finite lemma target:

```text
If a current-bank row has a hard mirror orbit,
then either it has a separating E/branch boundary-current gate
or it is random_covering_031.
```

Global proof target:

```text
hard_orbit_discharge
  follows from
    separating_current_transfer
    + HYP-3455/HYP-3460 random031 clause.
```

## Tournament Analysis

Vertices are joined proof carriers, not runners, arcs, or raw row names.

```text
pairwise_observable =
  predicate retention + hard payload + current payload + exception localization
switch =
  higher proof-facing carrier score; ties use declared route order
score_hist={19:1,29:1,42:1,47:1,54:2,58:1,60:1}
directed_3cycles=0
hamiltonian_path=
  J00_hard_orbit_current_join
  -> J01_singleton_intersection_ledger
  -> J02_separating_current_transfer
  -> J03_random031_named_gluing_debt
  -> J04_hard_mirror_orbit_ledger
  -> J05_dead_touch_gate_universal_lemma
  -> J06_raw_exception_set_overlap
  -> J07_raw_hard_delta_count
```

## Assumption Challenge

Alternate vertices considered: runners, residues, individual survivor gates,
mirror orbits, dead-cover components, blocker labels, section boundaries,
wall-crossing events, owner-current vectors, and proof obligations.

Chosen carrier vertices: hard mirror orbits joined to row-level E/branch
boundary-current cuts.

Preserved LRC predicate: whether a hard gate debt has an immediate
current-transfer exit or must be routed to a named terminal clause.

Destroyed information: exact within-row interval geometry and component order,
except for retained typed orbit, structural orbit, and best-current sidecars.

This is not an LRC14 proof.  It is a finite dispatch lemma that makes the next
hard-orbit proof obligation smaller and more explicit.
