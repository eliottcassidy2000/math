---
id: HYP-3477
title: LRC14 hard mirror-orbit discharge audit
status: EVIDENCE / exact hard-family discharge ledger; not an LRC14 proof
source: codex-2026-06-29 continuation of HYP-3475 hard-delta mirror-orbit ledger after incoming HYP-3476 pair-current and exception-frontier packets
tangent: T1437
technique: LTI-437
tournament_technique: LTT-337
script: 04-computation/lrc14_hard_mirror_orbit_discharge_codex_20260629.py
result: 05-knowledge/results/lrc14_hard_mirror_orbit_discharge_codex_20260629.out
reflection: 07-reflections/lrc14-hard-mirror-orbit-discharge-codex-20260629.md
related:
  - HYP-3476
  - HYP-3475
  - HYP-3474
  - HYP-3473
  - HYP-3472
  - HYP-3471
  - HYP-3462
  - HYP-3461
  - HYP-3460
  - HYP-3455
  - HYP-3453
  - HYP-3451
  - HYP-3438
  - HYP-3436
  - THM-523
  - OPEN-Q-108
---

# HYP-3477: LRC14 Hard Mirror-Orbit Discharge Audit

HYP-3477 audits the eight HYP-3475 hard mirror orbits with cover-delta at
least `7`, joining exact q=`14V` phase-grid witnesses, HYP-3472
dead-projection currents, and a bounded pair-current sidecar.  It is the
orbit-level companion to the HYP-3476 pair-current audit and exception-frontier
router.  It is a finite-bank discharge ledger, not a proof of LRC14.

## Rebase Integration

The active HYP-3476 packets are compatible.  The pair-current audit closes the
two AP separating-only rows and shows the seven random edge exceptions are
zero-edge singleton projections, so larger E/branch pair cuts are not the
right invariant there.  The exception-frontier router then shows that the hard
mirror-orbit frontier and the boundary-currentless frontier intersect only at
`random_covering_031`.  HYP-3477 refines the hard side orbit by orbit: seven
hard orbits already have lower-delta projection-current exits, while
`random031` is the unique hard/currentless overlap that must route through
HYP-3455/HYP-3460 and the later gluing sidecars.

## Exact Readout

Script:

```text
04-computation/lrc14_hard_mirror_orbit_discharge_codex_20260629.py
```

Stored output:

```text
05-knowledge/results/lrc14_hard_mirror_orbit_discharge_codex_20260629.out
```

Hard family:

```text
random_covering_022, random_covering_031, random_covering_049,
random_covering_078, random_covering_080, random_covering_085,
random_covering_113
```

Aggregate result:

```text
hard_orbits=8
hard_rows=7
delta_hist={7: 7, 8: 1}
route_hist={
  single_e_branch_projection_cut: 7,
  phase_branch_bypass_with_lower_delta_component_hits: 1
}
phase_bypass_orbits=1/8
single_e_branch_projection_edge_cuts=7/8
pair_current_projection_edge_cuts=7/8
mirror_e_branch_projection_edge_cuts=7/8
```

The split is sharper than HYP-3475:

```text
7/8 hard orbits:
  hit by q=14V witnesses sometimes, but each row has a lower-delta
  E/branch gate giving a dead-cover projection-edge cut.

1/8 hard orbits:
  random_covering_031.  Its max-delta mirror pair has zero compatible
  q=14V phase-grid hits, while its hard components are hit through
  lower-delta compatible gates.  This remains the HYP-3455/HYP-3476
  gluing exception.
```

The unique delta-`8` hard orbit is on `random_covering_022`; it is not a
new obstruction because a lower-delta E/branch gate on the same row removes
`14` projection edges and a bounded pair/mirror sidecar removes `28`.

For `random_covering_031`, no single, pair, or mirror E/branch current in
the bounded HYP-3477 sidecar removes projection edges.  This confirms the
HYP-3476 exception-frontier router's singleton overlap: `random031` is not a
generic hard-orbit problem but the named finite gluing clause already isolated
by HYP-3455.

## Proof Pull

The next theorem target should split the hard-orbit debt as:

```text
non-random031 hard orbit
  -> lower-delta E/branch projection-current discharge

random031 hard orbit
  -> HYP-3455 seven-owner gluing clause
  -> HYP-3476 route sidecar
  -> owner-current / two-adic / signed-SPEC exit
```

The phase-grid test should be cited carefully.  It does not avoid every hard
orbit; it avoids exactly the `random031` max-delta pair.  For the other seven
orbits, the useful discharge is the lower-delta projection current.

## Tournament Analysis

Vertices are hard-orbit discharge proof carriers, not runners or raw gates.

```text
pairwise_observable =
  predicate retention + phase bypass + projection current
  + pair-current bridge + scalar firewall
score_hist={6:1,52:1,61:1,63:1,64:2}
directed_3cycles=0
hamiltonian_path =
  H00_phase_branch_hard_orbit_bypass
  -> H01_lower_delta_e_branch_projection_cut
  -> H02_pair_current_bridge_to_HYP3476
  -> H03_hard_mirror_orbit_gluing_clause
  -> H04_dead_projection_touch_sidecar
  -> H05_raw_hard_delta_scalar
```

Assumption challenge: runners, residues, endpoint colors, individual gates,
mirror orbits, q=`14V` phase-grid witnesses, dead-cover components, blocker
labels, pair currents, and proof obligations were considered.  The chosen
carrier preserves hard-gluing debt plus lower-delta projection-current exits.
It destroys raw runner order, scalar gate count, and untyped color mass.
