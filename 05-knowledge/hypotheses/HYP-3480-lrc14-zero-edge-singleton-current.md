---
id: HYP-3480
title: LRC14 zero-edge singleton-current audit
status: RESERVED STUB / computation in progress; not an LRC14 proof
source: codex-2026-06-29 continuation of HYP-3479 hard-orbit current join, HYP-3478 small-touch geometry, HYP-3476 pair-current/router, and HYP-3477 hard-orbit discharge
tangent: T1440
technique: LTI-440
tournament_technique: LTT-340
script: 04-computation/lrc14_zero_edge_singleton_current_codex_20260629.py
result: 05-knowledge/results/lrc14_zero_edge_singleton_current_codex_20260629.out
reflection: 07-reflections/lrc14-zero-edge-singleton-current-codex-20260629.md
related:
  - HYP-3479
  - HYP-3478
  - HYP-3477
  - HYP-3476
  - HYP-3475
  - HYP-3472
  - HYP-3471
  - HYP-3460
  - HYP-3455
  - HYP-3453
  - HYP-3451
  - HYP-3438
  - THM-523
  - OPEN-Q-108
---

# HYP-3480: LRC14 Zero-Edge Singleton-Current Audit

HYP-3476 shows that the random boundary-current exceptions are not missing a
larger E/branch pair cut: their dead-cover projections are edgeless singleton
packets.  HYP-3477 and the HYP-3476 router remove the hard-currented rows and
isolate `random_covering_031` as the unique hard/currentless overlap.  This
stub follows the incoming HYP-3478 small-touch/no-hard geometry reservation and
HYP-3479 hard-orbit/current join, reserving the executable singleton-current
packet for the remaining non-hard zero-edge rows.

## Planned Scope

Primary rows:

```text
random_covering_001
random_covering_039
random_covering_062
random_covering_074
random_covering_086
random_covering_101
```

Control row:

```text
random_covering_031
```

The planned script will join HYP-3472 dead-cover projections, HYP-3476
pair-current support labels, HYP-3476 route labels, HYP-3477 hard-orbit data,
and the S319/HYP-3472 unit-delta split when available.  It will test whether
the zero-edge rows are carried by a small owner-label involution, unit-delta
singleton-current packet, endpoint-spine certificate, two-adic residue class,
or named hard/gluing exception.

Existing HYP-3478 companions already split the six primary rows into three
clean best-touch rows (`062`, `086`, `101`), asymmetric branch-unit row `001`,
and cover-delta sidecar rows `039`, `074`.  HYP-3480 should preserve that split
while adding the `random_covering_031` hard/currentless control and route flags.
HYP-3481 now supplies that control as a mirror-punctured annulus plus bypassed
saddle atlas, so the combined audit should keep random031 topological payloads
separate from the six-row singleton-current packet.

## Tournament Analysis Plan

Vertices are terminal singleton-current proof carriers, not runners, raw row
names, or scalar gate counts.  Candidate vertices include owner-label pairs,
dead singleton components, unit-delta gate words, route labels, hard-overlap
flags, two-adic residue payloads, signed-SPEC exits, and formal proof
obligations.  The quotient must preserve whether a zero-edge row is a pure
small-touch singleton packet or the named random031 hard/gluing overlap.
