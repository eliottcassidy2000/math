# LRC14 Dead-Cover Pair-Current Reflection

**Date:** 2026-06-29
**Hypothesis:** HYP-3476
**Artifacts:**
- `04-computation/lrc14_dead_cover_pair_current_codex_20260629.py`
- `05-knowledge/results/lrc14_dead_cover_pair_current_codex_20260629.out`

## What Changed

The two-gate idea was worth testing because HYP-3472's single-gate current
left seven projection-edge exceptions and two additional separating-current
exceptions.  The result is not a universal two-gate closure, but it gives a
cleaner split:

```text
audited_dead_rows=9
audited_single_e_branch_cuts=330
audited_pair_e_branch_cuts=9113
exact_mirror_pairs=165
edge_exception_rows_closed_by_pair_edge_cut=0/7
separating_exception_rows_closed_by_pair_separating_current=2/9
```

The AP rows are genuine pair-current rows.  `covering_AP_with_84` and
`ap_omit_12_tail_84x01` each get a non-mirror two-gate separating current with
`removed_edges=36`, `largest_drop=2`, and `component_gain=2`.  The AP exact
mirror pair is still useful, removing `40` edges, but it does not separate.

The seven random rows are different.  Their dead-cover projections are already
edgeless singleton graphs:

```text
random exception projection edges = 0
edge_support_labels = ()
gate_edge_support_labels = ()
```

So no E/branch pair can create a projection-edge cut under the HYP-3472
deletion rule.  The failure is not "need a third gate"; it is "wrong current
carrier."

## Assumption Challenge

I considered runners, gaps, fixed circle sections, section boundaries,
wall-crossing events, residues, cover arcs, Fourier modes, matroid circuits,
survivor gates, mirror orbits, gate pairs, blocker labels, graph cuts, and
proof obligations as possible tournament vertices.  The script fixes vertices
as low-rank E/branch gate pairs plus finite cut/current obligations.

Preserved predicate:

```text
target dead row has a <=2-gate boundary-current carrier after adjacent-label deletion
```

Destroyed information:

```text
interval geometry, ordered branch orientation, endpoint sidecars, and row-level
owner-current data unless separately retained
```

HYP-3474 is the guardrail here: the pair quotient is legal for the AP rows
because it retains the actual projection/current predicate.  It is not legal
for random rows because the projection itself has no edge-support coordinate.

After pulling the next upstream work, the zero-edge random debt is already
partly sorted.  HYP-3477 reserves the hard mirror-orbit discharge lane and
intersects this HYP-3476 set at `random_covering_031`.  The incoming colored
gate unit-delta split intersects this HYP-3476 set at `random_covering_039`
and `random_covering_074`, both in its 19-row cover-delta sidecar packet.
That leaves `random_covering_001`, `random_covering_062`,
`random_covering_086`, and `random_covering_101` as the cleanest first test
for a unit-delta zero-edge singleton-current clause.

## Next Hook

Do not spend the next pass enumerating larger E/branch gate tuples unless the
tuple carries a new coordinate beyond adjacent labels.  The random exceptions
need a zero-edge singleton-current packet.  Good next tests:

- classify the isolated dead components by owner-current imbalance;
- compare each singleton pair against HYP-3455/HYP-3451 conductance and gluing;
- pull HYP-3460 phase-branch witnesses onto each isolated dead component;
- test whether the singleton labels form a two-adic or signed-SPEC/Rprime
  descent obstruction;
- expose a Lean-facing terminal exit:
  `pairCurrent AP | singletonCurrent random | namedDebt`.

Tournament fingerprint for HYP-3476:

```text
score_hist={7:1,11:1,50:1,56:1,58:2,60:1,66:2}
hamiltonian_path=P00_two_gate_pair_current_terminal
 -> P01_mirror_pair_current_terminal
 -> P03_colored_gate_formal_packet
 -> P02_single_gate_boundary_current
 -> P05_hard_mirror_orbit_discharge
 -> P04_partition_lattice_guardrail
 -> P06_typed_gate_pair_word
 -> P07_raw_pair_count
 -> P08_raw_label_union_size
```
