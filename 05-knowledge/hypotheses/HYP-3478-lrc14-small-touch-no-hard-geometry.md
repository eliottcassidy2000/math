# HYP-3478: LRC14 Small-Touch / No-Hard Geometry

**Status:** EVIDENCE / exact six-row companion audit; not a proof.

**Claimed by:** codex-2026-06-29  
**Tangent:** T1438  
**Technique handles:** LTI-438, LTT-338  
**Open question:** OPEN-Q-108
**Atlas script:** `04-computation/lrc14_small_touch_geometry_atlas_codex_20260629.py`
**Atlas result:** `05-knowledge/results/lrc14_small_touch_geometry_atlas_codex_20260629.out`
**Atlas reflection:** `07-reflections/lrc14-small-touch-geometry-atlas-codex-20260629.md`
**Companion script:** `04-computation/lrc14_small_touch_no_hard_geometry_codex_20260629.py`
**Companion result:** `05-knowledge/results/lrc14_small_touch_no_hard_geometry_codex_20260629.out`
**Companion reflection:** `07-reflections/lrc14-small-touch-no-hard-geometry-codex-20260629.md`

## Question

HYP-3476 leaves six non-AP boundary-current exceptions that have no HYP-3475
hard mirror orbit and no E/branch pair-current edge support:

```text
random_covering_001
random_covering_039
random_covering_062
random_covering_074
random_covering_086
random_covering_101
```

This file was the no-hard reservation for that geometry audit.  The active
HYP-3478 atlas is
`05-knowledge/hypotheses/HYP-3478-lrc14-small-touch-geometry-atlas.md`; this
companion keeps the exact no-hard row/current readout from the second audit.

## Exact Readout

Both HYP-3478 scripts use the HYP-3450/HYP-3451/HYP-3453/HYP-3438 row bank and
the HYP-3472/HYP-3476 current utilities.

The atlas readout:

```text
rows=6
dead_component_count_hist={2:5,4:1}
mirror_pair_count_hist={1:5,2:1}
projection_edge_hist={0:6}
edge_support_label_hist={0:6}
all_dead_components_singleton=True
all_hard_orbit_count_hist={0:6}
s319_min_gate_kind_hist={'branch_unit_delta':4,'delta_sidecar_packet':2}
component_owner_pair_hist={(165,179):2,(81,153):2,(63,129):2,
                           (9,81):2,(15,99):2,(133,169):2,(7,175):2}
component_complete_gate_count_hist={2:2,3:6,4:2,6:2,9:2}
```

The companion row/current readout:

```text
total_dead_components=14
projection_edge_rows=[]
singleton_cover_fail_rows=[]
mirror_failure_rows=[]
owner_unbalanced_rows=[]
min_gate_kind_hist={'branch_unit_delta': 4, 'delta_sidecar_packet': 2}
touching_e_branch_counts={
  random_covering_001: 56,
  random_covering_039: 12,
  random_covering_062: 26,
  random_covering_074: 42,
  random_covering_086: 34,
  random_covering_101: 30
}
min_gate_not_touching_dead_labels=['random_covering_101']
```

Every dead component is a singleton `B0`/`B1` owner-pair component, and the
exact interval mirror maps it to the component with the owner pair swapped.
The row-wise signed owner balance is therefore zero for all six rows.  The
dead-cover projection is edgeless in every row because the branch-coloured
dead labels are singletons, not because E/branch gates miss the packet.

The finite terminal split is:

```text
cover_delta_sidecar_rows = random_covering_039, random_covering_074
clean_unit_delta_rows = random_covering_001, random_covering_062,
                        random_covering_086, random_covering_101
```

`random_covering_101` has a shortest E/branch unit-delta gate that does not
touch dead labels, but it also has a shortest touching unit-delta E/branch
gate at component `21`, interval `[31/210,181/1225]`, with labels
`('B0:175','B1:7')`.

## Assumption Challenge

Do not assume the tournament vertices are runners, arcs, or raw row names.
The audit compares alternate vertex sets: singleton dead components,
mirror-paired dead components, blocker-owner pairs, fixed circle intervals,
section boundaries, touching gate events, residues, cover arcs, Fourier/color
sidecars, and proof obligations.  The chosen quotient preserves the LRC
predicate "which terminal discharge packet remains after the low-rank
E/branch producer" while explicitly recording what it destroys: interval
order, branch orientation, endpoint wall labels, and owner-current locality.

## Current Pull

The six rows reduce to mirror-paired singleton components.  Five rows have one
mirror pair and `random_covering_001` has two.  All dead-cover projections
have zero edges, all hard-orbit counts are zero, and every pocket has complete
or touching E/branch gate sidecars.  The next proof target is a finite
singleton-current lemma for the four branch-unit rows (`001`, `062`, `086`,
`101`), followed by the cover-delta sidecar clause for `039` and `074`.

Do not keep enlarging E/branch pair-current searches unless a new invariant
creates projection edges.

## Pointers

HYP-3478, HYP-3479, HYP-3477, HYP-3476, HYP-3475, HYP-3472, HYP-3471,
HYP-3455, HYP-3453, HYP-3451, HYP-3450, HYP-3438, THM-523, T1438, LTI-438,
LTT-338, OPEN-Q-108.
