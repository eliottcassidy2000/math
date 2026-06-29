---
id: HYP-3511
title: LRC14 random031 free-hole bracket atlas
status: EVIDENCE / exact free-hole packet atlas; not an LRC14 proof
source: monad-explorer-2026-06-29 continuation of HYP-3510, HYP-3486, HYP-3485, HYP-3484, HYP-3483, HYP-3482, HYP-3481, HYP-3490, HYP-3477, HYP-3460, and HYP-3455
tangent: T1511
technique: LTI-511
tournament_technique: LTT-411
script: 04-computation/lrc14_random031_free_hole_bracket_atlas_codex_20260629.py
result: 05-knowledge/results/lrc14_random031_free_hole_bracket_atlas_codex_20260629.out
reflection: 07-reflections/lrc14-random031-free-hole-bracket-atlas-codex-20260629.md
related:
  - HYP-3510
  - HYP-3490
  - HYP-3486
  - HYP-3485
  - HYP-3484
  - HYP-3483
  - HYP-3482
  - HYP-3481
  - HYP-3477
  - HYP-3460
  - HYP-3455
  - THM-523
  - OPEN-Q-108
---

# HYP-3511: LRC14 Random031 Free-Hole Bracket Atlas

## Claim

HYP-3486 showed that deleting the max-delta seam leaves `40` no-gate witness
cells in `14` legal mirror free-hole packets.  HYP-3511 sharpens that packet:
the free-hole cells are not a diffuse wall.

Exact split:

```text
14 legal mirror packets
= 10 individually ordinary-bracketed packets
+ 4 half-open packets in exactly 2 same-branch doublets.
```

Every exposed side of every packet or doublet is an ordinary endpoint-rank-`2`
boundary.  So every no-gate witness already lives inside a controlled
branch-order bead of the seam complement.

## Exact Readout

Computed by:

```text
04-computation/lrc14_random031_free_hole_bracket_atlas_codex_20260629.py
05-knowledge/results/lrc14_random031_free_hole_bracket_atlas_codex_20260629.out
```

Packet census:

```text
free_hole_packets=14
packet_size_hist={2:10,4:3,8:1}
individually_bracketed_packets=10
individually_bracketed_size_hist={2:7,4:3}
half_open_packets=4
half_open_size_hist={2:3,8:1}
half_open_cluster_count=2
half_open_cluster_packet_size_hist={2:2}
half_open_cluster_cell_size_hist={4:1,10:1}
```

Boundary rank:

```text
all_exposed_boundaries_ordinary=True
exposed_boundary_endpoint_rank_hist={(2,):48}
```

So the free-hole packet is locally rigid in the relevant sense: the only
nonordinary adjacency inside the branch-order carrier is the same-branch
gluing inside two small doublets.

Exact half-open doublets:

```text
cluster=(2,3)
  packet 2: branch0 u=(1149) phase=(1), branch1 u=(62) phase=(13)
            exposed ordinary components (93) and (4)
  packet 3: branch0 u=(1139) phase=(5), branch1 u=(72) phase=(9)
            exposed ordinary components (90) and (7)

cluster=(8,13)
  packet 8:  branch0 u=(447) phase=(13), branch1 u=(764) phase=(1)
             exposed ordinary components (35) and (62)
  packet 13: branch0 u=(453,454,455,456) phase=(5,6,7,8),
             branch1 u=(755,756,757,758) phase=(6,7,8,9)
             exposed ordinary components (38) and (59)
```

The ten remaining packets are individually ordinary-bracketed on both branch
sheets.

## Reframe

### 1. Free-hole packet is local, not global, debt

The free-hole carrier now has a finite lemma shape:

```text
10 ordinary-bracketed single packets
+ 2 ordinary-bracketed doublets.
```

This is much sharper than the scalar statement "`40` no-gate cells".  The
no-gate witnesses are interior seam-complement beads, not a hidden transport
wall.

### 2. Puncture shadow is suggestive, not literal

HYP-3481's continuous topology has `4` dead islands, while the q-grid free-hole
carrier has `14` mirror packets.  So the correct bridge is not equality of
counts.  The exact shared structure found here is smaller: two mirror-paired
half-open doublets carrying all same-branch free-hole adjacency.

## Compatibility With HYP-3486, HYP-3510, And HYP-3490

HYP-3486 gives the legal fiber trichotomy:

```text
242 rank-2 routed cells
40 free-hole cells in 14 packets
12 bypass cells in one pure mirror component.
```

HYP-3510 gives the opposite coarse collapse:

```text
branch order + survivor ports => one seam-complement incidence component.
```

HYP-3511 sits exactly between them.  It proves that the `40` free-hole cells
already have local bracket structure inside the seam complement.  HYP-3490 is
orthogonal: the pair-current / projection-edge carrier is blocked by a
private-label firewall, so the remaining transport object really is the
phase-side seam complement with this refined free-hole packet.

## Proof Pull

The random031 terminal packet should now split into three local lemmas:

1. **Rank-2 route lemma** from HYP-3486 for the `242` routed cells.
2. **Free-hole bracket lemma** from HYP-3511 for the `40` no-gate cells:
   `10` bracketed single packets plus `2` bracketed doublets.
3. **Pure bypass lemma** from HYP-3486 for the `12` lower-delta bypass cells.

This leaves owner-boundary puncture debt localized to the seam sidecar rather
than spread across the whole free-hole packet.

## Tournament Analysis

Vertices are proof carriers, not runners or raw phase counts:

```text
ordinary_bracketed_single_packet
half_open_doublet_packet
exposed_rank2_boundary
same_branch_free_adjacency_cluster
mirror_packet_pairing
free_hole_count_shadow
raw_40_cell_shadow
```

Fingerprint:

```text
score_hist={14:1,40:1,68:1,76:1,84:1,92:1,98:1}
directed_3cycles=0
hamiltonian_path=ordinary_bracketed_single_packet
  -> half_open_doublet_packet
  -> exposed_rank2_boundary
  -> same_branch_free_adjacency_cluster
  -> mirror_packet_pairing
  -> free_hole_count_shadow
  -> raw_40_cell_shadow
```

## Assumption Challenge

Do not treat the free-hole packet as a scalar count or as a raw list of q-grid
cells.  The useful vertices are packet types and exposed boundaries.

Chosen carrier: legal mirror packets plus branch-order boundary sidecars.
This preserves the seam-complement free-hole discharge predicate and destroys
only the low-signal cell-by-cell naming once local bracketing is retained.
