---
id: HYP-3486
title: LRC14 random031 seam-complement fiber graph
status: EVIDENCE / exact fiber-graph scout; not an LRC14 proof
source: codex-2026-06-29
tangent: T1446
technique: LTI-446
tournament_technique: LTT-346
script: 04-computation/lrc14_random031_seam_complement_fiber_graph_codex_20260629.py
result: 05-knowledge/results/lrc14_random031_seam_complement_fiber_graph_codex_20260629.out
reflection: 07-reflections/lrc14-random031-seam-complement-fiber-graph-codex-20260629.md
related:
  - HYP-3490
  - HYP-3485
  - HYP-3484
  - HYP-3483
  - HYP-3482
  - HYP-3481
  - HYP-3480
  - HYP-3479
  - HYP-3477
  - HYP-3460
  - HYP-3455
  - HYP-3438
  - THM-523
  - OPEN-Q-108
---

# HYP-3486: LRC14 Random031 Seam-Complement Fiber Graph

## Claim

After deleting the max-delta hard seam of `random_covering_031`, the q=`14V`
phase witnesses split into a finite mirror-run graph with no hidden wall:

```text
242 gate-routed witness cells -> endpoint-rank-2 seam-complement gates
40 free-hole witness cells    -> 14 mirror-closed free-hole packets
```

The lower-delta bypass is not diffuse.  It is one pure `12`-cell mirror
component on the hard seam components `(43,54)`, with owners `(23,93,113)`.

## Exact Readout

Computed by
`04-computation/lrc14_random031_seam_complement_fiber_graph_codex_20260629.py`
and stored in
`05-knowledge/results/lrc14_random031_seam_complement_fiber_graph_codex_20260629.out`.

Use the two-adic cylinder coordinate

```text
u = 2t mod 1, node=(u_index, branch), q=2422, base=q/2=1211.
```

Fiber census:

```text
witness_cells=282
occupied_fibers=258
fiber_size_hist={1:234,2:24}
cell_class_counts={bypass:12, free_hole:40, ordinary:230}
branch_class_counts:
  branch0: bypass 6, free_hole 20, ordinary 115
  branch1: bypass 6, free_hole 20, ordinary 115
hit_count_hist={0:40,1:242}
```

Seam-deleted routing:

```text
hard_seam_hits_after_delete=0
gate_routed_cells=242
free_hole_cells=40
hit_endpoint_rank_hist={2:242}
unique_hit_components=80
unique_hit_component_class_hist={both_alive:46, branch0_only:17, branch1_only:17}
```

Legal mirror-run graph, with horizontal u-neighbor edges on each branch and
mirror edges `(u,b)->(-u,1-b)`:

```text
component_count=79
component_size_hist={2:44,4:21,6:6,8:5,10:1,12:2}
component_type_hist={ordinary:64, bypass:1, free_hole:14}
routed_component_count=65
free_hole_component_count=14
free_hole_component_size_hist={2:10,4:3,8:1}
```

Pure bypass component:

```text
bypass_cells=12
bypass_components=(43,54)
bypass_owners=(23,93,113)
bypass_u_blocks_by_branch={0:(679..684), 1:(527..532)}
bypass_phase_blocks_by_branch={0:(7..12), 1:(2..7)}
```

## Quotient Guardrail

Mirror is legal and class-preserving:

```text
mirror_missing=0
mirror_pair_class_counts={
  bypass/bypass:6,
  free_hole/free_hole:20,
  ordinary/ordinary:115
}
```

Vertical half-turn gluing is not legal without a sidecar:

```text
vertical_halfturn_present_cells=48/282
vertical_pair_class_counts={
  free_hole/free_hole:2,
  free_hole/ordinary:14,
  ordinary/ordinary:8
}
vertical_glued_component_count=69
vertical_mixed_component_count=7
```

So `n*2` is the correct address coordinate, but not a free sheet quotient.  It
mixes free-hole and ordinary cells and never carries the bypass.

## Proof Pull

The HYP-3482 target should be sharpened:

```text
delete hard seam
=> every routed witness reaches endpoint-rank-2
   and every unrouted witness is in a mirror-closed free-hole packet
   and the hard-component bypass is one pure 12-cell mirror component.
```

This suggests three local lemmas:

1. **Rank-2 seam-complement routing.**  A witness cell hitting a surviving
   branch-compatible gate after hard-seam deletion discharges through the
   endpoint-rank-2 E/branch gate packet.
2. **Mirror free-hole packets.**  A witness cell with no compatible survivor
   gate is already a phase witness and should discharge through a free-hole
   lemma, not through a forced gate route.
3. **Pure bypass component.**  The 12-cell bypass on components `(43,54)` is
   a lower-delta rank-2 channel carrying owners `(23,93,113)`; the remaining
   seam-only owners `(45,147,169,173)` are boundary debt, not phase flow.

## Compatibility With HYP-3485 And HYP-3490

Incoming HYP-3485 is the connection atlas explaining why the seam-complement
picture matches older zipper, Cech, Menger, two-adic, and PGF routes.
Incoming HYP-3490 identifies `random_covering_031` as the unique
private-label firewall row that also carries hard mirror debt.  HYP-3486 is
complementary: the dead-projection label-union carrier cannot cut the row, but
the phase fiber graph shows that the seam-complement witnesses already split
into rank-2 exits, free-hole packets, and one pure bypass component.  The
packets should meet in the final proof as:

```text
private-label firewall says projection-current cannot discharge random031;
fiber graph says phase flow does not need the forbidden seam.
```

## Tournament Analysis

Vertices are fiber-graph proof carriers, not runners or raw gates:

```text
rank2_seam_complement_discharge
legal_mirror_run_graph
pure_twelve_cell_bypass_component
free_hole_mirror_packets
vertical_halfturn_guardrail
fiber_occupancy_word
raw_282_witness_count
```

Fingerprint:

```text
score_hist={12:1,65:1,71:1,75:1,82:1,88:1,96:1}
directed_3cycles=0
hamiltonian_path=
  rank2_seam_complement_discharge
  -> legal_mirror_run_graph
  -> pure_twelve_cell_bypass_component
  -> free_hole_mirror_packets
  -> vertical_halfturn_guardrail
  -> fiber_occupancy_word
  -> raw_282_witness_count
```

## Assumption Challenge

Considered vertices: runners, raw gates, u-fibers, branch cells, horizontal
flow runs, mirror pairs, vertical half-turn fibers, component exits, endpoint
ranks, free-hole packets, and proof obligations.

Chosen vertices: witness cells and legal mirror-run components.

Preserved predicate: terminal random031 discharge after deleting the hard
seam.

Destroyed by bad quotients: branch legality, free-hole status, bypass purity,
mirror class preservation, endpoint rank, and the distinction between address
projection (`n*2`) and legal topology.
