---
id: HYP-3453
title: LRC14 gate-escape transversal router
status: EVIDENCE / exact component-gate join; not an LRC14 proof
source: codex-2026-06-29 continuation of HYP-3438 and HYP-3450/HYP-3451
tangent: T1413
technique: LTI-413
tournament_technique: LTT-313
script: 04-computation/lrc14_gate_escape_transversal_router_codex_20260629.py
result: 05-knowledge/results/lrc14_gate_escape_transversal_router_codex_20260629.out
reflection: 07-reflections/lrc14-gate-escape-transversal-router-codex-20260629.md
related:
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
  - HYP-3417
  - HYP-3129
  - HYP-2963
  - THM-523
  - OPEN-Q-108
---

# HYP-3453: LRC14 Gate-Escape Transversal Router

## Claim

HYP-3450 proved, on the audited bank, that every row has an endpoint-rank
`<= 2` survivor component.  HYP-3438 separately decomposed mixed `E_safe`
components into exact survivor gate words.  HYP-3453 joins these two ledgers.

Rebase integration: incoming HYP-3452 proves the AP-with-`84m` component phase
audit through `m=70`, and incoming HYP-3439 bridges one-branch rescue cores to
survivor gates and component escapes.  HYP-3453 is the broader bank-level
companion, showing when an escape is a graph-composable gate rather than a
clean-component exit.

The key compression failure is that a "low-rank survivor component" is not a
single proof object.  It splits into:

```text
clean escape: no bad-core/dead-cover obstruction is present in the component;
gate escape: a survivor gate lies on the boundary of bad-core blocks.
```

Only the gate escape composes directly with the HYP-3451 dead-cover graph,
because it retains the parent even wall, endpoint labels, branch mask, and
adjacent cover deltas.  Clean escapes are already no-obstruction exits.

## Exact Readout

Script:

```text
04-computation/lrc14_gate_escape_transversal_router_codex_20260629.py
```

Stored result:

```text
05-knowledge/results/lrc14_gate_escape_transversal_router_codex_20260629.out
```

Aggregate join over the same `135` primitive covering rows:

```text
rows_with_low_rank_component=135/135
rows_with_low_rank_gate=133/135
rows_with_dead_components=130/135
rows_with_dead_components_and_low_rank_gate=130/130
rows_dead_without_low_rank_gate=[]
dead_zero_rows=[
  random_covering_007, random_covering_014, random_covering_044,
  random_covering_053, random_covering_089
]
dead_zero_clean_only_rows=[random_covering_044, random_covering_053]
all_survivor_gates=8702
low_rank_gates=8666
gate_endpoint_rank_hist={2:8666, 3:36}
low_rank_component_total=13352
low_rank_gate_parent_components=6204
clean_low_rank_components=7148
```

The proof-facing implication is exact on the audited bank:

```text
dead_components(row) > 0  =>  row has a rank <= 2 survivor gate.
no rank <= 2 survivor gate => dead_components(row) = 0.
```

The two rows without a low-rank gate are not hidden counterexample candidates.
They have no dead components and no mixed survivor gates:

```text
random_covering_044:
  components=108, dead=0, gates=0, low_rank_components=108
random_covering_053:
  components=86, dead=0, gates=0, low_rank_components=86
```

## AP-Tail Base Case

The AP-with-`84m` danger family remains finite and structured.  For
`{1,...,11,13,84}`:

```text
components=26, dead=22, gates=4, low_rank_components=4, low_rank_gates=4
gates:
  [8/49,97/588]     labels=B1:7|E:84  delta=(1,1)
  [33/196,6/35]     labels=E:84|B1:5  delta=(1,1)
  [29/35,163/196]   labels=B0:5|E:84  delta=(1,1)
  [491/588,41/49]   labels=E:84|B0:7  delta=(1,1)
```

For `84m`, `m=2,3,4`, the same four gate pattern persists.  At `m=5`, the
row has `104` dead components, two rank-`2` boundary gates, and two additional
rank-`1` clean components bounded by `E:420|E:420`.  This is a useful stress:
the proof must allow clean exits without losing the gate transversal.

## Proof Pull

HYP-3451 asked for a proof that full two-colour blocker saturation of
`E_safe` is impossible.  HYP-3453 sharpens the local target:

```text
If a dead-cover obstruction exists at all, it has a rank <= 2 gate escape.
```

Thus a counterexample cannot hide behind the scalar statement "there is a
low-rank component."  It must either have no dead-cover obstruction, or pass
through an exact survivor gate with endpoint labels, branch mask, adjacent
cover deltas, and parent even wall retained.

The next rigorous step is to prove this gate-transversal implication for
primitive covering rows, then feed the gate into the HYP-3451 Menger,
Green-current, or algebraic-connectivity proof route.

## Tournament Analysis

Vertices are escape proof carriers, not runners or scalar survivor masses.

Pairwise observable:

```text
predicate retention + component/gate join + cut-graph composability
+ clean-case exit + endpoint payload
```

Fingerprint:

```text
score_hist={10:1, 38:1, 52:1, 55:2, 59:1}
directed_3cycles=0
hamiltonian_path =
  E00_dead_positive_rank2_gate_transversal
  -> E01_clean_obstruction_empty_exit
  -> E02_gate_endpoint_delta_spine
  -> E03_component_cover_conductance_join
  -> E04_low_rank_component_shadow
  -> E05_raw_survivor_measure
```

Assumption challenge: runners, gaps, fixed sections, section boundaries,
wall-crossing events, residues, cover arcs, endpoint walls, bad-core blocks,
survivor gates, clean `E_safe` components, dead components, branch-coloured
blockers, and proof obligations were considered.  The chosen carrier preserves
the component-cover obstruction predicate while recording whether an escape is
clean-only or graph-composable.
