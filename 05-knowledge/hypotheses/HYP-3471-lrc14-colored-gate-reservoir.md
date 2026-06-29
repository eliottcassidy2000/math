---
id: HYP-3471
title: LRC14 colored gate-reservoir
status: EVIDENCE / exact colored gate quotient audit; not an LRC14 proof
source: codex-2026-06-29 continuation of HYP-3453 with older colored-reservoir/discrepancy routes, layered above HYP-3462 AP84 corridor-splice certificate, HYP-3470 AP84 CRT-placement bridge, HYP-3461 colored-extension gate carrier, HYP-3460 phase-branch color pullback, HYP-3459 AP84 color-packet legality, and HYP-3458 AP84 coloring-recursion
tangent: T1431
technique: LTI-431
tournament_technique: LTT-331
script: 04-computation/lrc14_colored_gate_reservoir_codex_20260629.py
result: 05-knowledge/results/lrc14_colored_gate_reservoir_codex_20260629.out
reflection: 07-reflections/lrc14-colored-gate-reservoir-codex-20260629.md
related:
  - HYP-3470
  - HYP-3462
  - HYP-3461
  - HYP-3460
  - HYP-3459
  - HYP-3458
  - HYP-3457
  - HYP-3456
  - HYP-3455
  - HYP-3454
  - HYP-3453
  - HYP-3452
  - HYP-3451
  - HYP-3450
  - HYP-3439
  - HYP-3438
  - HYP-3437
  - HYP-3436
  - HYP-3435
  - HYP-3434
  - HYP-3431
  - HYP-3429
  - HYP-3427
  - HYP-3426
  - HYP-3425
  - HYP-3422
  - HYP-3418
  - HYP-3415
  - HYP-2595
  - HYP-2594
  - HYP-2593
  - THM-523
  - OPEN-Q-108
---

# HYP-3471: LRC14 Colored Gate-Reservoir

## Claim

The older "14-colored reservoir" is not just a color-count heuristic.  In the
current component-cover/gate stack it becomes a typed gate word.

On the HYP-3453 bank, the exact sharpened finite lemma is:

```text
dead_components(row) > 0
  => row has a rank <= 2 survivor gate with one E endpoint
     and one B0/B1 endpoint.
```

Thus a dead-cover obstruction must touch a bichromatic even/branch boundary,
not merely some low-rank survivor interval.  The legal quotient is:

```text
colored gate word = typed mod-14 endpoint residues
                  + branch mask
                  + adjacency
                  + adjacent cover deltas.
```

Endpoint kind alone and numeric residue alone are useful diagnostics, but they
are not legal proof quotients for arbitrary rows.

## Exact Readout

Script:

```text
04-computation/lrc14_colored_gate_reservoir_codex_20260629.py
```

Stored output:

```text
05-knowledge/results/lrc14_colored_gate_reservoir_codex_20260629.out
```

Aggregate audit:

```text
rows_audited=135
rows_with_dead_components=130/135
low_rank_gates=8666
e_branch_low_rank_gates=7002
same_branch_low_rank_gates=1482
cross_branch_low_rank_gates=182
dead_rows_with_e_branch_low_rank_gate=130/130
dead_rows_without_e_branch_low_rank_gate=[]
dead_rows_only_e_branch_low_rank_gates=28
dead_rows_with_same_branch_low_rank_gate=102
dead_rows_with_cross_branch_low_rank_gate=54
```

Compression ladder:

```text
endpoint_kind_palette_size=8
numeric_mod14_palette_size=147
typed_mod14_palette_size=360
structural_palette_size=161
full_color_palette_size=1727
row_kind_set_count=8
row_typed_set_count=123
row_structural_set_count=107
```

Endpoint-kind histogram:

```text
B1|E: 1762
E|B0: 1762
E|B1: 1739
B0|E: 1739
B0|B0: 741
B1|B1: 741
B1|B0: 98
B0|B1: 84
```

The AP84 four-color packet is real but not universal:

```text
canonical_ap_palette=[
  B0:5|E:0,
  B1:7|E:0,
  E:0|B0:7,
  E:0|B1:5
]
rows_with_canonical_ap_palette=68/135
dead_rows_with_canonical_ap_palette=67/130
dead_rows_without_canonical_ap_palette_count=63
```

After the HYP-3462/HYP-3470 rebases, the split is sharper: HYP-3462 closes the
AP84 corridor-splice carrier, HYP-3470 carries exact `q=14V` CRT color-grid
witnesses for the canonical AP84 tail, and HYP-3471 carries the local typed
gate word once a dead-cover obstruction has been compressed to a survivor-gate
boundary.  The first is a structural AP splice, the second is a placement
sidecar, and the third is a proof-carrier sidecar.

The minimum E/branch gate per dead row is usually one of the four unit-delta
structural types:

```text
B0|E branch0 left_bad_edge  (1,1): 31 rows
B1|E branch1 left_bad_edge  (1,1): 30 rows
E|B1 branch1 right_bad_edge (1,1): 28 rows
E|B0 branch0 right_bad_edge (1,1): 21 rows
```

The remaining minimum gates carry small cover-delta sidecars; these are the
finite repair coordinates for random/non-AP rows.

## Prior Coloring Connections

HYP-2593's phase-color reservoir already made CRT witnesses into color hits in
`G_P cap C_b(E)`.  HYP-3471 identifies the local wall-level descendant of that
idea: the color is now attached to the two endpoint labels of a survivor gate.

HYP-2594/HYP-2595 warned that colored discrepancy scalars lose component and
boundary information.  HYP-3471 pinpoints the corresponding current loss:
numeric mod-14 colors collapse `360` typed residue gate words to `147` numeric
words, while endpoint kind collapses the whole bank to only `8` gate colors and
`8` row-level kind sets.  Those quotients are too coarse for the random rows.

The old "beyond 3-coloring" analogy is useful only as a compression warning.
This is not an ordinary proper-coloring theorem: there are same-branch gates
`B0|B0` and `B1|B1`, and cross-branch gates `B0|B1` / `B1|B0`.  The proof
payload is not "proper coloring exists"; it is "dead obstruction forces an
E/branch gate, while same-branch gates require owner-delta/gluing sidecars."

## Proof Pull

HYP-3453 gave:

```text
dead_components(row) > 0 => rank <= 2 survivor gate.
```

HYP-3471 sharpens it to:

```text
dead_components(row) > 0 => rank <= 2 E/branch survivor gate.
```

The next rigorous finite lemma should prove this implication without relying on
the sample bank.  A plausible route:

1. Use HYP-3451's branch-coloured blocker graph to show a dead-cover component
   has boundary current crossing from an even wall to one branch wall.
2. Use HYP-3438 gate words to convert that crossing into a survivor gate.
3. Use HYP-3471's color ladder to retain typed residue and sidecar data until
   the crossing is routed into one of:

```text
AP84 packet and exact color-grid placement:
  HYP-3462/HYP-3470/HYP-3460/HYP-3459/HYP-3458/HYP-3454/HYP-3456/HYP-3457
random031 same-branch/gluing stress: HYP-3455
canonical corridor fence: HYP-3431
general component conductance/Menger route: HYP-3451
owner-current/two-adic/signed-SPEC debt: HYP-3417/HYP-3422/HYP-3129
```

Tournament Analysis uses colored proof carriers, not runners or raw color
counts.  The Hamiltonian path begins:

```text
dead_positive_e_branch_gate
-> full_colored_gate_word
-> structural_color_sidecar
-> typed_mod14_gate_word
-> ap84_four_color_packet
-> endpoint_kind_coloring
-> numeric_14_residue
-> raw_color_count
```

Assumption challenge: runners, residue colors, endpoint-kind colors, typed
mod-14 colors, gaps, fixed circle sections, section boundaries, wall-crossing
events, cover arcs, survivor gates, dead-cover components, branch-coloured
blockers, and proof obligations were considered.  The chosen quotient preserves
the dead-cover-to-E/branch-gate predicate and the graph-composable endpoint
payload.  It destroys row-specific full wall geometry, so it is a finite lemma
target and not a global LRC14 proof.
