---
id: HYP-3490
title: LRC14 private-label firewall audit
status: EVIDENCE / all-union currentless-row firewall; not an LRC14 proof
source: codex-2026-06-29 continuation of HYP-3484, HYP-3483, HYP-3482, HYP-3481, HYP-3480, HYP-3479, HYP-3478, HYP-3472, HYP-3476, HYP-3477, HYP-3475, and HYP-3451
tangent: T1450
technique: LTI-450
tournament_technique: LTT-350
script: 04-computation/lrc14_private_label_firewall_codex_20260629.py
result: 05-knowledge/results/lrc14_private_label_firewall_codex_20260629.out
reflection: 07-reflections/lrc14-private-label-firewall-codex-20260629.md
related:
  - HYP-3484
  - HYP-3483
  - HYP-3482
  - HYP-3481
  - HYP-3480
  - HYP-3479
  - HYP-3478
  - HYP-3477
  - HYP-3476
  - HYP-3475
  - HYP-3474
  - HYP-3473
  - HYP-3472
  - HYP-3471
  - HYP-3462
  - HYP-3460
  - HYP-3455
  - HYP-3453
  - HYP-3451
  - HYP-3450
  - HYP-3438
  - THM-523
  - OPEN-Q-108
---

# HYP-3490: LRC14 Private-Label Firewall Audit

## Claim

HYP-3476 proves that two low-rank E/branch gates do not repair the seven random
HYP-3472 projection-edge exceptions.  HYP-3490 explains why: on exactly those
seven rows, every blocker label touched by every low-rank E/branch gate is
private to a single dead component.

Therefore no union of adjacent E/branch gate labels, of any size, can delete a
dead-cover projection edge.  The HYP-3476 pair failure is not a two-gate
accident; it is an all-union private-label firewall for this carrier.

## Exact Readout

Script:

```text
04-computation/lrc14_private_label_firewall_codex_20260629.py
```

Stored output:

```text
05-knowledge/results/lrc14_private_label_firewall_codex_20260629.out
```

Aggregate ledger:

```text
rows_audited=135
dead_rows=130
dead_rows_with_h3472_projection_edge_cut=123/130
dead_rows_with_h3472_separating_current=121/130
private_firewall_rows=7/130 [
  random_covering_001,
  random_covering_031,
  random_covering_039,
  random_covering_062,
  random_covering_074,
  random_covering_086,
  random_covering_101
]
shared_touch_rows=123/130
mismatch_private_firewall_vs_no_edge_cut=[]
mismatch_shared_touch_vs_edge_cut=[]
```

Frontier split:

```text
route_hist={
  ordinary_or_hard_currented:121,
  ap84_edge_only_pair_current_terminal:2,
  random031_private_firewall_hard_overlap:1,
  small_touch_private_firewall:6
}
```

The private firewall rows have no projection edge labels at all in the
dead-cover projection touched by E/branch gates:

```text
random_covering_001: dead_components=4, touch_multiplicities={1:8}
random_covering_031: dead_components=4, touch_multiplicities={1:8}
random_covering_039: dead_components=2, touch_multiplicities={1:4}
random_covering_062: dead_components=2, touch_multiplicities={1:4}
random_covering_074: dead_components=2, touch_multiplicities={1:4}
random_covering_086: dead_components=2, touch_multiplicities={1:4}
random_covering_101: dead_components=2, touch_multiplicities={1:4}
```

Only `random_covering_031` is both private-firewalled and hard-orbit-bearing.
The other six are the small-touch/no-hard packet from the HYP-3476 frontier
router.  After rebasing over the completed HYP-3478 audit, this six-row packet
is sharper: all `14` dead components are singleton `B0`/`B1` owner-pair
components with exact mirror-swapped partners, all six dead-cover projections
are edgeless, there is no row-wise owner imbalance, and every row has complete
or touching E/branch sidecars.  HYP-3478 splits the six into cover-delta
sidecar rows `random_covering_039`, `random_covering_074` and clean
branch-unit-delta rows `random_covering_001`, `random_covering_062`,
`random_covering_086`, `random_covering_101`.

After the completed HYP-3480 singleton-current audit, the six-row terminal
packet is stronger again: all `14/14` target dead components have complete
branch-unit touches and all `7/7` mirror pairs have mirror-compatible
branch-unit gate pairs.  The executable split is now clean rows
`random_covering_062`, `random_covering_086`, `random_covering_101`,
asymmetric row `random_covering_001`, and cover-delta minimum-shadow rows
`random_covering_039`, `random_covering_074`; `random_covering_031` remains
the `0/4` complete-branch-unit hard control routed through HYP-3484/HYP-3483/HYP-3482/HYP-3481.

## Assumption Challenge

Candidate vertices considered: runners, gaps, residues, survivor gates, mirror
orbits, gate pairs, dead components, blocker labels, labelled projection
edges, quotient fibers, and proof obligations.

The chosen carrier vertices are dead-cover blocker labels with multiplicity
sidecars, especially labels touched by low-rank E/branch gates.  This preserves
the predicate:

```text
whether any union of adjacent E/branch gate labels can remove an edge from
the dead-cover projection.
```

It destroys raw interval geometry, runner order, and scalar gate counts, but it
retains exactly the shared-label multiplicity that HYP-3476 pair-current needs
and lacks on the seven random rows.

## Proof Pull

The current frontier should now split as:

```text
AP currentless rows
  -> non-private labels
  -> HYP-3476 two-gate balanced current terminal

random031
  -> private-label firewall + hard mirror orbit
  -> HYP-3484 forbidden-seam flow geometry
  -> HYP-3483 recursion-flow comparator
  -> HYP-3482 forbidden-seam atlas
  -> HYP-3481 random031 topology atlas
  -> HYP-3455 seven-owner gluing / owner-current / two-adic / SPEC

six remaining random currentless rows
  -> private-label firewall + no hard orbit
  -> HYP-3478 mirror-singleton packet
  -> HYP-3480 executable singleton-current audit
  -> mirror-unit singleton certificate on 14/14 dead components and 7/7 mirror pairs
  -> clean rows 062/086/101 + asymmetric 001 + cover-delta shadows 039/074
  -> finite singleton-current, endpoint-spine, exact-period, or state-lift exit
```

HYP-3490 is especially useful for formalization: it replaces a large pair
enumeration by a label-multiplicity lemma.  If every touched label has
multiplicity `1`, deleting any set of such labels cannot alter the unlabelled
dead-cover projection graph.

## Tournament Analysis

Vertices are label-multiplicity proof carriers, not runners or raw exception
names.

```text
pairwise_observable =
  predicate retention + private-label payload + projection-cut firewall
  + frontier split + quotient guardrail
score_hist={7:1,17:1,58:1,59:1,64:1,67:2}
directed_3cycles=0
hamiltonian_path =
  F00_private_label_firewall_theorem
  -> F01_random031_private_hard_overlap
  -> F02_small_touch_private_packet
  -> F03_ap84_nonprivate_pair_current
  -> F04_edge_support_label_axis
  -> F05_raw_pair_current_count
  -> F06_raw_exception_name
```
