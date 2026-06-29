---
id: HYP-3472
title: LRC14 colored gate unit-delta split
status: EVIDENCE / exact HYP-3471 sharpening; not an LRC14 proof
source: monad-explorer-2026-06-29 continuation of HYP-3471, separating the H3435 structured bank from the random bank
tangent: T1432
technique: LTI-432
tournament_technique: LTT-332
script: 04-computation/lrc14_colored_gate_unit_delta_split_codex_20260629.py
result: 05-knowledge/results/lrc14_colored_gate_unit_delta_split_codex_20260629.out
reflection: 07-reflections/lrc14-colored-gate-unit-delta-split-codex-20260629.md
related:
  - HYP-3471
  - HYP-3462
  - HYP-3470
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
  - HYP-3435
  - HYP-3431
  - THM-523
  - OPEN-Q-108
---

# HYP-3472: LRC14 Colored Gate Unit-Delta Split

## Claim

HYP-3471's finite lemma

```text
dead_components(row) > 0
  => row has a rank <= 2 E/branch survivor gate
```

is much sharper on the audited `135`-row bank than the headline suggests.

The minimum such gate splits into three exact packets:

```text
structured dead rows:  branch-specific unit-delta edge gates only
random dead rows:      mostly the same
                       + exactly one both-branch unit-delta row
                       + exactly 19 single-branch cover-delta exceptions
```

So the residual is not "random rows behave arbitrarily."  It is a small
cover-delta sidecar packet living on the same `E/B` single-bad-edge skeleton.

## Exact Readout

Script:

```text
04-computation/lrc14_colored_gate_unit_delta_split_codex_20260629.py
```

Stored output:

```text
05-knowledge/results/lrc14_colored_gate_unit_delta_split_codex_20260629.out
```

Aggregate split:

```text
rows_audited=135
dead_rows=130
structured_dead_rows=15
random_dead_rows=115
min_kind_hist={
  branch_unit_delta: 110,
  both_unit_delta: 1,
  delta_sidecar_packet: 19
}
```

Structured packet:

```text
structured_unit_delta_rows=15/15
structured_min_typed_hist={
  B1:7|E:0: 10,
  E:0|B1:5: 4,
  E:0|B0:5: 1
}
structured_min_struct_hist={
  (B1|E, branch1, left_bad_edge, 1,1): 10,
  (E|B1, branch1, right_bad_edge, 1,1): 4,
  (E|B0, branch0, right_bad_edge, 1,1): 1
}
```

Random residual:

```text
random_branch_unit_delta_rows=95/115
random_both_unit_delta_rows=1/115
random_delta_sidecar_rows=19/115
```

The `19`-row sidecar packet is still very rigid:

```text
all 19 rows:
  - are single-bad-edge minima
  - are branch-specific (10 branch0, 9 branch1)
  - only need delta vectors in {(1,2),(2,1),(2,2),(1,3)}
```

Structural histogram on that packet:

```text
(B1|E, branch1, left_bad_edge, 1,2): 4
(E|B0, branch0, right_bad_edge, 1,2): 3
(B0|E, branch0, left_bad_edge, 1,2): 2
(B0|E, branch0, left_bad_edge, 2,1): 2
(E|B0, branch0, right_bad_edge, 2,1): 2
(E|B1, branch1, right_bad_edge, 2,2): 1
(B1|E, branch1, left_bad_edge, 2,1): 1
(E|B0, branch0, right_bad_edge, 2,2): 1
(E|B1, branch1, right_bad_edge, 2,1): 1
(E|B1, branch1, right_bad_edge, 1,2): 1
(E|B1, branch1, right_bad_edge, 1,3): 1
```

Exception rows are exactly

```text
random_covering_022, 036, 039, 047, 051, 054, 056, 058, 059, 063,
069, 074, 083, 085, 090, 095, 104, 107, 113.
```

## Interpretation

HYP-3471's proof carrier survives another round of controlled forgetting.
The AP84 corridor-splice packet from HYP-3462 was already rigid; HYP-3472
shows that the full structured H3435 bank is also rigid at the minimum-gate
level.  The only non-unit behavior is a finite random sidecar packet, not a
breakdown of the `E/B` boundary alphabet.

That changes the next proof target.  The theorem now does not need to explain
arbitrary colored survivor geometry.  It needs only:

1. a structural proof of the unit-delta packet for the structured branch, and
2. a conductance/Menger or gate-gluing route for the `19` random delta-sidecar
   rows and the lone both-branch unit row `random_covering_088`.

## Proof Pull

The next finite theorem target should be:

```text
dead_components(row) > 0
  => minimum E/branch gate is either
     (a) a branch-specific unit-delta edge gate,
     (b) the unique both-branch unit-delta row, or
     (c) one of the 19 named cover-delta sidecar rows.
```

This reduces the post-HYP-3471 work to a finite packet theorem.

- Route the structured unit-delta packet through HYP-3462/HYP-3431/HYP-3454/HYP-3456/HYP-3457.
- Route the `19`-row sidecar packet through HYP-3451 conductance/Menger or HYP-3455 gluing sidecars.
- Keep HYP-3470/HYP-3460/HYP-3459/HYP-3458 as color-placement / color-legality sidecars, not as substitutes for the local gate theorem.
