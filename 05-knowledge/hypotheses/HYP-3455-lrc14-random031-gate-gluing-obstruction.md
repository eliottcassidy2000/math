---
id: HYP-3455
title: LRC14 random031 gate-gluing obstruction
status: EVIDENCE / exact single-row obstruction isolator; not an LRC14 proof
source: codex-2026-06-29 continuation of HYP-3439/HYP-3438/HYP-3450/HYP-3451
tangent: T1415
technique: LTI-415
tournament_technique: LTT-315
script: 04-computation/lrc14_random031_gate_gluing_obstruction_codex_20260629.py
result: 05-knowledge/results/lrc14_random031_gate_gluing_obstruction_codex_20260629.out
reflection: 07-reflections/lrc14-random031-gate-gluing-obstruction-codex-20260629.md
related:
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
  - HYP-3425
  - HYP-3422
  - HYP-3418
  - HYP-3415
  - THM-523
  - OPEN-Q-108
---

# HYP-3455: LRC14 Random031 Gate-Gluing Obstruction

## Claim

`random_covering_031` is the named noncanonical rank-`6` obstruction left by
the broad HYP-3439 bridge search.  It is not a diffuse scalar obstruction:
HYP-3455 isolates it as a finite gate-gluing clause joining:

```text
HYP-3437 rank-6 overlap rescue graph
HYP-3438 survivor-gate words
HYP-3450/HYP-3451 component-cover escape router
```

## Exact Readout

The row is:

```text
speeds=(12,23,45,55,58,70,84,93,113,120,147,169,173)
odd=(23,45,55,93,113,147,169,173)
even_half=(6,29,35,42,60)
```

One-branch overlap-tax data:

```text
branch0_measure=18860031217607/184101593891215
naive_slack=-130807232086254/1288711157238505
rescue_rank=6
rescue_subset=(23,45,93,113,147,169)
rescue_margin=2280704077944161/92787203321172360
rescue_graph_edges=15
rescue_graph_connected=True
```

Survivor-gate data:

```text
mixed_components=84
survivor_gates=138
branch_mask_hist={both:8, branch0:65, branch1:65}
max_total_delta=7
max_delta_gate_count=2
max_delta_owner_union=(23,45,93,113,147,169,173)
rescue_subset_inside_max_delta_owner_union=True
extra_max_delta_owners=(173,)
max_delta_mirror_pairs=[(43,54)]
```

The two maximal gates are mirror partners:

```text
component 43:
  word=B-S-B-S-B
  survivor=[71/161,349/791]
  len=4/18193
  branch_mask=branch0
  labels=B0:23|B0:113
  delta=(3,4)
  covers=B0(23):B1(93)->B0(45,113):B1(147,169,173)

component 54:
  word=B-S-B-S-B
  survivor=[442/791,90/161]
  len=4/18193
  branch_mask=branch1
  labels=B1:113|B1:23
  delta=(4,3)
  covers=B0(147,169,173):B1(45,113)->B0(93):B1(23)
```

Component-cover graph data:

```text
components=98
dead_components=4
alive_components=94
low_rank_escape=94
max_dead_pair_rank=2
class_hist={both_alive:50, branch0_only:22, branch1_only:22, dead_both:4}
best_rank=2
best_survivor=[29/385,97/1211]
best_labels=L[B0:55] R[B0:173]
danger_score=0.000868
```

## Proof Pull

The noncanonical rank-`6` failure is now a seven-owner local clause:

```text
rescue core = (23,45,93,113,147,169)
gate owner union = (23,45,93,113,147,169,173)
```

Any proof route for arbitrary primitive rows should split:

```text
canonical rank-6 -> HYP-3431/HYP-3452 AP corridor and tail phase
noncanonical rank<=5 -> bounded overlap core plus HYP-3436/HYP-3438 sidecars
random031 rank-6 -> finite mirror gate-gluing clause
```

The next lemma target is: the two max-delta survivor gates cannot be glued into
a full two-colour cover without exposing one of the `94` low-rank component
escapes, an endpoint-spine/window route, an owner-current imbalance, a
two-adic descent exit, or signed-SPEC/Rprime debt.

Tournament Analysis uses proof obligations and finite gate/gluing carriers as
vertices.  The Hamiltonian path begins:

```text
max_delta_survivor_gate_pair
-> rank6_rescue_overlap_graph
-> two_color_component_escape_router
-> owner_delta_gluing_clause
-> mirror_involution_cut_word
```

Assumption challenge: HYP-3455 explicitly rejects the shortcut that a
noncanonical rank-`6` rescue core represents an uncontrolled family.  In the
audited row, it is a finite mirror-symmetric gate-gluing clause with exact
endpoint labels, owner deltas, component-cover escapes, and a named extra owner
`173`.
