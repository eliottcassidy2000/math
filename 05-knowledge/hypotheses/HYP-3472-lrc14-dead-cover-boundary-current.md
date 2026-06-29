---
id: HYP-3472
title: LRC14 dead-cover boundary-current audit
status: RESERVED STUB / graph-current audit pending; not a proof
source: codex-2026-06-29 continuation of HYP-3451/HYP-3453/HYP-3471 after the HYP-3462 AP84 corridor splice
tangent: T1432
technique: LTI-432
tournament_technique: LTT-332
script: 04-computation/lrc14_dead_cover_boundary_current_codex_20260629.py
result: 05-knowledge/results/lrc14_dead_cover_boundary_current_codex_20260629.out
reflection: 07-reflections/lrc14-dead-cover-boundary-current-codex-20260629.md
related:
  - HYP-3472
  - HYP-3471
  - HYP-3462
  - HYP-3470
  - HYP-3461
  - HYP-3460
  - HYP-3459
  - HYP-3458
  - HYP-3455
  - HYP-3453
  - HYP-3451
  - HYP-3450
  - HYP-3438
  - HYP-3437
  - HYP-3436
  - HYP-3417
  - HYP-3129
  - THM-523
  - OPEN-Q-108
---

# HYP-3472: LRC14 Dead-Cover Boundary-Current Audit

This reserves the graph/current sibling of HYP-3471. HYP-3471 showed that
every audited row with a dead component has a rank-`<=2` E/branch survivor
gate. HYP-3451 projected dead components into a branch-coloured blocker graph.
The missing bridge is to check whether the E/branch gates touch the dead-cover
projection as small boundary-current cuts rather than merely existing in the
same row.

Known input:

```text
HYP-3451: dead components have a blocker projection graph with low-rank escape data.
HYP-3453: dead_components(row)>0 implies a rank<=2 survivor gate on the 135-row bank.
HYP-3471: dead_components(row)>0 implies a rank<=2 E/branch survivor gate on the same bank.
HYP-3462/HYP-3470: AP84 sidecars are closed enough to serve as named AP packets.
```

Planned computation:

1. Join every HYP-3453 row to its HYP-3451 dead-cover blocker projection.
2. For each low-rank E/branch gate, record the adjacent blocker labels it
   exposes on the two bad sides.
3. Remove those labels from the dead-cover incidence graph and test whether the
   gate is a projection cut, a separating cut, or only a local gate touch.
4. Record the current vector of each gate as branch-labelled blocker flow:
   `B0` labels, `B1` labels, balance, removed projection edges, change in
   largest dead-cover component, and row exceptions.
5. Run Tournament Analysis over proof carriers, not runners: projection cut,
   separating current, low-rank E/branch implication, AP84 closed packet,
   random031 seven-owner clause, typed gate word, raw gate count, and raw dead
   fraction.

Assumption challenge:

- Alternate vertices considered: runners, gaps, fixed circle sections, section
  boundaries, wall-crossing events, residues, cover arcs, Fourier modes,
  matroid circuits, survivor gates, dead components, blocker labels,
  owner-current vectors, graph cuts, and proof obligations.
- Chosen vertices for this audit: dead-cover blocker labels, E/branch survivor
  gates, projection-cut witnesses, and proof obligations.
- Preserved LRC predicate: whether a dead-cover obstruction has a graph-current
  route to a low-rank E/branch survivor gate that can feed HYP-3451/HYP-3455.
- Destroyed information: full wall geometry, exact interval lengths, scalar
  color counts, and same-row gate multiplicity beyond the cut/current sidecar.

This is not an LRC14 proof. It is a finite audit meant to determine whether
the HYP-3471 implication can be promoted into a boundary-current transfer
lemma, or whether the remaining non-AP rows must be named as owner-current,
two-adic, signed-SPEC, or state-lift debt.
