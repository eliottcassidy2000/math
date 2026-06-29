---
id: HYP-3481
title: LRC14 random031 topology atlas
status: EVIDENCE / topology reframing; not an LRC14 proof
source: codex-2026-06-29 continuation of HYP-3455, HYP-3460, HYP-3476, HYP-3477, HYP-3479, and HYP-3480
tangent: T1441
technique: LTI-441
tournament_technique: LTT-341
script: 04-computation/lrc14_random031_topology_atlas_codex_20260629.py
result: 05-knowledge/results/lrc14_random031_topology_atlas_codex_20260629.out
reflection: 07-reflections/lrc14-random031-topology-atlas-codex-20260629.md
related:
  - HYP-3480
  - HYP-3479
  - HYP-3478
  - HYP-3477
  - HYP-3476
  - HYP-3475
  - HYP-3460
  - HYP-3455
  - HYP-3453
  - HYP-3451
  - HYP-3450
  - HYP-3438
  - THM-523
  - OPEN-Q-108
---

# HYP-3481: LRC14 Random031 Topology Atlas

## Claim

`random_covering_031` should be treated as a small topological object, not as
a large scalar exception.

HYP-3479 reduces the hard mirror-orbit / boundary-current join to
`separating_current_transfer + random031_named_clause`, and HYP-3480 reserves
the zero-edge singleton-current lane with random031 as its hard control.
HYP-3481 is the local topology atlas for that named random031 clause.

Two reframings are now exact on the audited bank:

1. **Mirror-punctured annulus.** The dead-cover projection has four isolated
   singleton islands, paired by mirror.  It is not a connected graph waiting
   for a Menger edge cut.
2. **Bypassed saddle seam.** The max-delta mirror pair on components `(43,54)`
   is a short antipodal seam carrying the seven-owner boundary, but the
   q=`14V` phase flow has zero hits on that seam and twelve hits on
   lower-delta mirror-bypass gates on the same two components.

## Exact Readout

Script:

```text
04-computation/lrc14_random031_topology_atlas_codex_20260629.py
```

Stored output:

```text
05-knowledge/results/lrc14_random031_topology_atlas_codex_20260629.out
```

Dead-cover topology:

```text
dead_projection=components=4 largest=1 edges=0
dead_island_mirror_pairs=[(0,3),(1,2)]
```

The four dead islands are:

```text
island 0: [113/840,55/406], labels=('B0:45','B1:23'), rank=2
island 1: [57/406,69/490], labels=('B0:113','B1:93'), rank=2
island 2: [421/490,349/406], labels=('B0:93','B1:113'), rank=2
island 3: [351/406,727/840], labels=('B0:23','B1:45'), rank=2
```

The hard seam is the HYP-3455 max-delta mirror pair:

```text
component 43: [71/161,349/791], branch0, delta=(3,4)
component 54: [442/791,90/161], branch1, delta=(4,3)
central_mirror_gap=93/791
owners=(23,45,93,113,147,169,173)
```

The same two components carry a lower-delta mirror bypass:

```text
component 43: [344/791,286/651], branch1, delta=(0,2)
component 54: [365/651,447/791], branch0, delta=(2,0)
owners=(23,93,113)
```

Phase-flow readout:

```text
q=14V witnesses=282
hard_gate_hits=0
lower_delta_hard_component_hits=12
hard_component_hit_counter={(43,'branch1',2):6,(54,'branch0',2):6}
```

## Proof Pull

The row has two different local mechanisms living on the same mirror pair of
components:

```text
max-delta seam       = continuous seven-owner gluing debt
lower-delta bypass   = phase-compatible same-component escape
```

This suggests two lemma targets:

1. **Mirror-puncture lemma.** Four edgeless dead-cover islands with
   mirror-paired singleton blockers should route by local island current, not
   by projection-edge cuts.
2. **Bypassed-saddle lemma.** A max-delta mirror seam with zero q=`14V` phase
   flux and lower-delta same-component flux is removable unless its
   seven-owner boundary produces explicit owner-current/two-adic/SPEC debt.

This reframes random031 as a removable saddle with a named owner-boundary
sidecar.  The proof should avoid trying to cut the dead-cover projection:
there are no edges to cut.

## Tournament Analysis

Vertices are proof carriers rather than runners:

```text
mirror_punctured_annulus
bypassed_saddle_seam
seven_owner_boundary
phase_branch_flow
dead_projection_scalar
raw_rank6_scalar
```

Fingerprint:

```text
score_hist={10:1,21:1,61:1,63:1,66:1,68:1}
directed_3cycles=0
hamiltonian_path=
  mirror_punctured_annulus
  -> bypassed_saddle_seam
  -> seven_owner_boundary
  -> phase_branch_flow
  -> dead_projection_scalar
  -> raw_rank6_scalar
```

## Assumption Challenge

Considered vertices:

```text
runners, rescue owners, dead islands, mirror seams, phase witnesses,
survivor gates, blocker labels, cover arcs, Fourier modes, proof obligations
```

Chosen vertices are dead islands, mirror seams, and phase-flow bypasses.  This
preserves the terminal discharge predicate and destroys raw runner order only
after replacing it by mirror, owner-boundary, and phase-flux sidecars.
