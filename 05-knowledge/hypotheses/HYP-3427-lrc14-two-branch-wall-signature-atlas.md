---
id: HYP-3427
title: LRC14 two-branch wall-signature atlas
status: EVIDENCE / exact wall-certificate scout; not an LRC14 proof
source: codex-2026-06-28 continuation after HYP-3426 mirror audit, HYP-3425 Helly audit, and energy-sheet sidecar warning
tangent: T1388
technique: LTI-388
tournament_technique: LTT-288
script: 04-computation/lrc14_two_branch_wall_signature_atlas_codex_20260628.py
result: 05-knowledge/results/lrc14_two_branch_wall_signature_atlas_codex_20260628.out
reflection: 07-reflections/lrc14-two-branch-wall-signature-atlas-codex-20260628.md
related:
  - HYP-3426
  - HYP-3425
  - HYP-3424
  - HYP-3423
  - HYP-3422
  - HYP-3421
  - HYP-3420
  - HYP-3419
  - HYP-3418
  - HYP-3417
  - HYP-3415
  - HYP-3140
  - HYP-3129
  - HYP-2963
  - THM-523
  - OPEN-Q-108
---

# HYP-3427: LRC14 Two-Branch Wall-Signature Atlas

## Claim

HYP-3425 gives a positive two-branch survivor set:

```text
good = E_safe \ (B0_odd cap B1_odd).
```

HYP-3426 removes the branch-choice ambiguity by the mirror involution and
endpoint-owner support audit.  HYP-3427 asks for the next proof-facing
certificate: a survivor interval should not be recorded only by its positive
measure; it should carry an exact wall signature:

```text
branch mask + left wall labels + right wall labels + midpoint binders.
```

The walls are:

```text
E:s    even-speed wall  ||(s/2)u|| = 1/14
O0:o   branch-0 odd wall ||o u/2|| = 1/14
O1:o   branch-1 odd wall ||o u/2|| = 3/7
```

The theorem target becomes:

```text
every primitive covering row has a survivor window with a bounded wall
signature, or else emits named owner-current, sheet, exact-period, or
state-lift debt.
```

## Exact Readout

Script:

```text
04-computation/lrc14_two_branch_wall_signature_atlas_codex_20260628.py
```

Stored result:

```text
05-knowledge/results/lrc14_two_branch_wall_signature_atlas_codex_20260628.out
```

The exact audit covers `67` rows: curated HYP-3422/HYP-3425 rows, the
canonical `{1..11,13,84m}` tower through `m=20`, eight two-tail probes, and
`35` deterministic primitive covering rows.

```text
rows with survivor windows:       67/67
total survivor windows:           5524
global signature types:           27
branch mask histogram:            b0=2255, b1=2255, both=1014
```

The top wall-type signatures are dominated by even walls:

```text
(E,E)       2694
(O1,E)      585
(E,O0)      585
(E,O1)      557
(O0,E)      557
(O1,O1)     250
```

Midpoint binders across all audited windows remain mixed but floor-facing:

```text
14Q=570, even_R=2425, odd_unit=2055, seven_R=474
```

## Tight Row Certificate

For the tight canonical row

```text
S = {1,2,3,4,5,6,7,8,9,10,11,13,84}
```

the four survivor windows are exactly:

```text
(8/49,97/588)     width 1/588  branch1  O1:7 -> E:84
(33/196,6/35)     width 3/980  branch1  E:84 -> O1:5
(29/35,163/196)   width 3/980  branch0  O0:5 -> E:84
(491/588,41/49)   width 1/588  branch0  E:84 -> O0:7
```

So the visible one-tail obstruction is not only positive.  It has a four-window
certificate with one even wall `E:84` and odd walls `5` and `7`.

## Proof Meaning

This is a stricter target than HYP-3425.  HYP-3425 asks for a positive
component after removing the two-color bad core.  HYP-3427 asks future agents
to prove that such a component can be named by a compact wall word.

The canonical tower changes shape: for `m<=5` all survivor windows are
branch-only; starting at `m=6`, `both` windows appear and eventually dominate
the tower.  The top canonical signature becomes `('both','E','E')`, suggesting
that the large-tail regime may be governed by even-wall-to-even-wall windows,
while the base tight row is governed by odd/even boundary alternation.

## Guardrail

Do not use this as a new scalar statistic.  The wall word is a certificate
object.  A route that forgets which exact even wall, odd branch wall, and
midpoint binders define the survivor interval has lost the proof payload.

This also reads HYP-3425's energy-sheet warning conservatively: additive energy
or SPEC data may prioritize a packet, but the interval proof still needs a
wall certificate or a named sidecar debt.

## Tournament Analysis

Vertices are proof carriers / wall certificates, not runners or raw intervals.

```text
score_hist={-28:1, 49:1, 65:1, 73:1, 76:1, 81:1, 85:1}
directed_3cycles=0
hamiltonian_path=
  W00_wall_signature_certificate
  -> W01_survivor_component_normal_form
  -> W02_branch_mask_descent_router
  -> W03_owner_current_wall_exception
  -> W04_energy_sheet_sidecar_join
  -> W05_raw_positive_measure_audit
  -> W06_named_analogy_without_walls
```

## Next Hook

Try to prove a finite certificate lemma of the form:

```text
Every primitive covering S=O union 2E has at least one survivor window whose
wall word lies in a bounded legal alphabet.
```

If the bounded alphabet fails in a larger search, record the first failure by
the missing wall type: even-cover, odd branch, owner-current, sheet, exact
period, or state-lift debt.
