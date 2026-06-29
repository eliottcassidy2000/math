---
id: HYP-3438
title: LRC14 survivor-gate word audit
status: EVIDENCE / exact survivor-gap gate-word classification; not a proof
source: codex-2026-06-29 continuation of HYP-3436 minimal bad-core cover extractor
tangent: T1399
technique: LTI-399
tournament_technique: LTT-299
script: 04-computation/lrc14_survivor_gate_word_audit_codex_20260629.py
result: 05-knowledge/results/lrc14_survivor_gate_word_audit_codex_20260629.out
reflection: 07-reflections/lrc14-survivor-gate-word-audit-codex-20260629.md
related:
  - HYP-3451
  - HYP-3450
  - HYP-3437
  - HYP-3436
  - HYP-3435
  - HYP-3434
  - HYP-3431
  - HYP-3429
  - HYP-3428
  - HYP-3427
  - HYP-3426
  - HYP-3425
  - HYP-3422
  - HYP-3417
  - HYP-3129
  - THM-523
---

# HYP-3438: LRC14 Survivor-Gate Word Audit

HYP-3438 is the immediate local-to-global follow-up to HYP-3436.

Known input: HYP-3436 classifies

```text
E_safe cap B0_odd cap B1_odd
```

component-by-component and shows that every audited primitive covering row
still has survivor gaps.  The tight row `{1,...,11,13,84}` has only four mixed
even-safe components, with survivor mass `1/105`, while most other even-safe
components are fully bad or clean.

The executed audit decomposes each mixed even-safe component into a gate word
whose letters are survivor gaps and adjacent bad-core blocks.  Each survivor
gate retains interval endpoints, endpoint wall labels, left/right bad-core
cover pairs, cover-owner deltas, branch mask, and the parent even wall.

## Exact Readout

Script:

```text
04-computation/lrc14_survivor_gate_word_audit_codex_20260629.py
```

Stored result:

```text
05-knowledge/results/lrc14_survivor_gate_word_audit_codex_20260629.out
```

Aggregate audit on the same `135` primitive covering rows as HYP-3436:

```text
mixed_E_components=6228
clean_E_components=7166
fully_bad_E_components=3770
survivor_gates_inside_mixed_components=8702
```

Gate words are alternating bad/survivor words.  The largest observed mixed
component has five survivor gates:

```text
max_survivors_in_one_mixed_component=5
  random_covering_036 component 32 word S-B-S-B-S-B-S-B-S
max_segments_in_one_mixed_component=10
  random_covering_078 component 40 word B-S-B-S-B-S-B-S-B-S
```

Branch masks are symmetric:

```text
survivor_branch_mask_hist={both:1064, branch0:3819, branch1:3819}
gate_adjacency_hist={left_bad_edge:3515, right_bad_edge:3515, two_sided:1672}
```

Endpoint kinds stay small and structured:

```text
B1|E:1769, E|B0:1769, E|B1:1746, B0|E:1746,
B0|B0:744, B1|B1:744, B1|B0:98, B0|B1:86
parent_endpoint_kind_hist={E|E:8702}
```

The dominant adjacent-cover delta is `(1,1)`, but hard random gates can require
up to `(4,4)`:

```text
cover_delta_hist has (1,1):5950
hardest displayed deltas: (4,4) on random_covering_022
```

## Tight Canonical Row

For

```text
S={1,2,3,4,5,6,7,8,9,10,11,13,84}
```

the four HYP-3436 mixed components become four one-gate edge words:

```text
component 3:  word=B-S, survivor=[8/49,97/588], len=1/588,
              branch_mask=branch1, labels=B1:7|E:84
component 4:  word=S-B, survivor=[33/196,6/35], len=3/980,
              branch_mask=branch1, labels=E:84|B1:5
component 21: word=B-S, survivor=[29/35,163/196], len=3/980,
              branch_mask=branch0, labels=B0:5|E:84
component 22: word=S-B, survivor=[491/588,41/49], len=1/588,
              branch_mask=branch0, labels=E:84|B0:7
```

This is the exact HYP-3431 corridor-fence shadow: all four survivor gates are
edge gates between an even wall and a low odd branch wall, with adjacent cover
delta `(1,1)`.

## Proof Target

The finite lemma should be stated as a gate-word obstruction:

```text
Every mixed E_safe component emits an exact survivor gate word.
Any proof using survivor mass must reconstruct or route:
  adjacent bad-core covers
  endpoint labels
  branch mask
  parent even wall.
```

Legal exits are now finite: direct branch relocation, HYP-3431 corridor-fence,
HYP-3429/HYP-3427 endpoint-spine/wall certificates, HYP-3417/HYP-3420
owner-current exceptions, HYP-3428 two-adic loss, HYP-3437 overlap-cut bridge,
or HYP-3129 signed-SPEC.  Raw survivor measure, mixed-component counts,
harmonic budgets, and topology labels are negative controls unless they
reconstruct one of those sidecars.

Rebase integration: HYP-3450/HYP-3451 are the whole-component sibling of this
gate-word quotient.  They choose obstruction cuts in the branch-alive/dead
component-cover projection; HYP-3438 then classifies the mixed components
inside those cuts by endpoint labels, adjacent cover deltas, branch masks, and
parent even walls.

## Tournament Analysis

Vertices are survivor-gate proof carriers, not runners or raw gaps.

```text
pairwise_observable =
  predicate retention + gate exactness + endpoint/cover-delta payload
  + scalar-firewall safety
score_hist={12:1, 52:1, 57:2, 58:1, 61:1, 64:1, 66:1, 67:1}
directed_3cycles=0
hamiltonian_path =
  G00_exact_survivor_gate_word
  -> G01_adjacent_cover_delta_ledger
  -> G02_endpoint_wall_alternation
  -> G03_branch_mask_relocation_witness
  -> G04_corridor_fence_recognizer
  -> G05_overlap_cut_bridge
  -> G06_owner_current_exception_router
  -> G07_signed_SPEC_route
  -> G08_raw_survivor_measure
```

Assumption challenge: runners, gaps, fixed sections, section boundaries,
wall-crossing events, residues, cover arcs, endpoint walls, branch-bad owners,
bad-core blocks, survivor gaps, cover-pair deltas, mixed even components,
owner-current labels, and proof obligations were considered.  The chosen
carrier preserves the two-adic relocation predicate and records exactly what
scalarization destroys.
