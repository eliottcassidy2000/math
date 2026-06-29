---
id: HYP-3476
title: LRC14 exception-frontier router
status: EVIDENCE / finite packet-router audit; not an LRC14 proof
source: codex-2026-06-29 continuation of HYP-3475, HYP-3474, HYP-3473, HYP-3472, and HYP-3455
tangent: T1436
technique: LTI-436
tournament_technique: LTT-336
script: 04-computation/lrc14_exception_frontier_router_codex_20260629.py
result: 05-knowledge/results/lrc14_exception_frontier_router_codex_20260629.out
reflection: 07-reflections/lrc14-exception-frontier-router-codex-20260629.md
related:
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
  - THM-523
  - OPEN-Q-108
---

# HYP-3476: LRC14 Exception-Frontier Router

## Claim

HYP-3472 left a finite boundary-current exception frontier.  HYP-3475 left a
finite hard mirror-orbit frontier.  These are not the same obstruction.

On the current `135`-row bank, the two frontiers intersect in exactly one row:

```text
random_covering_031
```

This is the HYP-3455/HYP-3460 row already known as the seven-owner
mirror-gluing and phase-branch bypass clause.  Everything else in the union
routes into one of three simpler terminal packets:

```text
hard_mirror_orbit_with_separating_current
small_touch_no_hard_current_exception
ap84_edge_cut_without_separating_current
```

Thus the next proof target should be a labelled exception-packet theorem, not
a single scalar hard-gate inequality.

## Exact Readout

Executable scout:

```text
04-computation/lrc14_exception_frontier_router_codex_20260629.py
```

Stored output:

```text
05-knowledge/results/lrc14_exception_frontier_router_codex_20260629.out
```

Frontier ledger:

```text
boundary_rows_without_projection_edge_cut=7
boundary_rows_without_separating_current=9
hard_mirror_orbit_rows_delta_ge_7=7
frontier_union_rows=15
intersection_hard_and_currentless=1 ['random_covering_031']
hard_rows_with_separating_current=6/7
currentless_rows_without_hard_orbit=8/9
```

The route partition over all `135` rows is:

```text
ordinary_separating_current: 115
nondead: 5
hard_mirror_orbit_with_separating_current: 6
small_touch_no_hard_current_exception: 6
ap84_edge_cut_without_separating_current: 2
random031_overlap_hard_and_currentless: 1
```

The six hard rows outside `random031` all already have a separating
E/branch boundary current:

```text
random_covering_022
random_covering_049
random_covering_078
random_covering_080
random_covering_085
random_covering_113
```

The six non-AP currentless rows outside `random031` have no hard mirror orbit:

```text
random_covering_001
random_covering_039
random_covering_062
random_covering_074
random_covering_086
random_covering_101
```

Their best touching gates have total adjacent-delta range:

```text
small_touch_best_delta_range=(2,3)
small_touch_best_current_hist={(1,1):5, (2,1):1}
```

The two AP84 currentless rows are edge-only rather than touch-only:

```text
covering_AP_with_84
ap_omit_12_tail_84x01
```

Both have the same projection-edge cut:

```text
gate=[8/49,97/588]
kind=B1|E
labels=('B0:13','B1:7')
removed_edges=20
current=(1,1)
delta=(1,1)
```

## Quotient Guardrail

HYP-3474's colored-gate axes were:

```text
K,N,T,S,F,C,M,A
```

HYP-3476 tests whether any subset of those axes preserves the new six-label
terminal route packet.  The answer on the current bank is negative:

```text
fewest_axis_route_pure=none among axes=('K','N','T','S','F','C','M','A')
most_compressed_route_pure=none among existing colored-gate axes
required_sidecar_R=terminal_route_packet fibers=6 max_fiber=115
```

So the route label `R` is not optional bookkeeping.  A later proof must either
carry `R` as a sidecar or prove a reconstruction theorem for it.

## Rebase Integration: S319 Gate Unit-Delta Split

After this scout was committed locally, incoming coordination added an S319
gate unit-delta split handoff.  That handoff uses a conflicting local
`HYP-3472` name, so this note refers to it as S319 until the namespace is
reconciled.

S319 splits minimum E/B gate kinds into:

```text
branch_unit_delta = 110
both_unit_delta = 1  (random_covering_088)
delta_sidecar = 19
```

Its delta-sidecar rows intersect the HYP-3476 frontier but do not contain or
replace the HYP-3476 route partition:

```text
S319_delta_sidecar cap HYP3476_frontier =
  {random_covering_022, random_covering_039, random_covering_074,
   random_covering_085, random_covering_113}

S319_delta_sidecar cap HYP3476_small_touch =
  {random_covering_039, random_covering_074}

S319_delta_sidecar cap HYP3476_hard_currented =
  {random_covering_022, random_covering_085, random_covering_113}
```

So S319 is an upstream gate-kind coordinate.  HYP-3476's route sidecar `R`
remains necessary because it also separates AP84 edge-only rows, random031's
unique hard/currentless overlap, ordinary current rows, nondead rows, and
frontier rows not in the S319 delta-sidecar packet.

## Proof Pull

The finite audited form is:

```text
currentless_frontier cap hard_orbit_frontier = {random_covering_031}
```

This splits the remaining LRC14 route into four packet lemmas:

1. Ordinary rows: use the HYP-3472 separating E/branch boundary current.
2. Hard-current rows: discharge hard mirror debt after the lower-rank
   separating current has been named.
3. Small-touch currentless rows: prove a bounded owner-current/two-adic/SPEC
   lemma for touch-only rows with adjacent delta `<=3`.
4. AP84 edge-only rows: route through the closed AP84
   endpoint/corridor/color packet.

The only row that carries both hard mirror debt and failed separating current
is `random_covering_031`, already isolated by HYP-3455 and softened by
HYP-3460.

## Tournament Analysis

Vertices are terminal exception packets, not runners, raw exception names, or
individual gates.

```text
pairwise_observable =
  frontier partition + boundary-current payload + mirror-orbit payload
  + AP84 splice + random031 gluing + quotient firewall
score_hist={12:1,45:1,46:2,49:1,54:1,58:1}
directed_3cycles=0
hamiltonian_path =
  R00_orthogonal_exception_packet_theorem
  -> R01_unique_overlap_random031_clause
  -> R02_hard_orbit_with_separating_current
  -> R04_ap84_edge_only_closed_packet
  -> R05_partition_lattice_route_label
  -> R03_small_touch_no_hard_boundary_rows
  -> R06_raw_exception_name_list
```

Assumption challenge: considered runners, gaps, fixed circle sections, section
boundaries, wall-crossing events, residues, cover arcs, individual survivor
gates, mirror orbits, dead-cover labels, quotient fibers, and proof
obligations.  The chosen carrier preserves the terminal discharge predicate
after the universal E/branch gate producer has fired.  It destroys full row
geometry unless the terminal route sidecar `R` is retained.
