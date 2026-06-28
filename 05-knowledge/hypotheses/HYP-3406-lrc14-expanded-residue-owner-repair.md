---
id: HYP-3406
title: LRC14 expanded-bank residue failure and owner-support repair
status: EVIDENCE / exact enlarged-bank audit; not an LRC14 proof
source: monad-explorer-2026-06-28
tangent: T1367
technique: LTI-367
tournament_technique: LTT-267
script: 04-computation/lrc14_expanded_residue_owner_repair_codex_20260628.py
result: 05-knowledge/results/lrc14_expanded_residue_owner_repair_codex_20260628.out
reflection: 07-reflections/lrc14-expanded-residue-owner-repair-codex-20260628.md
related:
  - HYP-3405
  - HYP-3404
  - HYP-3403
  - HYP-3402
  - HYP-3311
  - HYP-3310
  - HYP-3301
  - HYP-3265
  - HYP-3260
  - HYP-3259
  - HYP-3258
  - HYP-3257
  - HYP-3253
  - HYP-2975
  - HYP-2969
  - HYP-2963
  - THM-523
  - OPEN-Q-108
---

# HYP-3406: LRC14 Expanded-Bank Residue Failure And Owner-Support Repair

## Claim

HYP-3404's top-ranked residue-word breakpoint lead is executable, and it is
the expanded-bank companion to HYP-3405's AP-collar finite-lemma certificate.
HYP-3311's nonunit residue word is only a curated-bank separator.  On broader
HYP-2963 banks, residue-word exactness fails in two qualitatively different
ways:

1. a same-residue height leak;
2. a stronger endpoint-owner leak that survives both `v2` and exact nonunit
   height data.

On the scanned enlarged banks, the exact next repair is not
`residue + v2` or `residue + height`.  It is

```text
coarse sheaf base + residue word + endpoint-owner support word.
```

So HYP-3402's owner-current route is the next exact theorem-facing sidecar on
the enlarged packet bank, while tropical height data explains only the first
subcase.

## Exact Readout

The executable audit enlarges the bank from the curated HYP-2969 sample to
four HYP-2963 banks:

```text
(single_limit, two_swap_limit) =
  (20,4), (30,8), (48,12), (60,16).
```

Controlled-forgetting summary:

```text
(20,4):  residue 1 mixed, residue+v2 0, residue+height 0,
         residue+unit_qsqrt7 1, residue+owner_support 0

(30,8):  residue 2 mixed, residue+v2 1, residue+height 1,
         residue+unit_qsqrt7 1, residue+owner_support 0

(48,12): residue 2 mixed, residue+v2 1, residue+height 1,
         residue+unit_qsqrt7 1, residue+owner_support 0

(60,16): residue 2 mixed, residue+v2 1, residue+height 1,
         residue+unit_qsqrt7 1, residue+owner_support 0
```

So the owner-support repair stays exact through the largest scanned bank
(`872` rows).

## First Leak: Height

The first residue-word failure already appears at `(20,4)`:

```text
P10+GW
GW-shell alias 12->132
```

Both lie in the same coarse-plus-residue fiber, but have different theorem
exits:

```text
P10+GW                 -> unit-petal-named
GW-shell alias 12->132 -> positive-Haar-open
```

This collision is repaired by `v2` / exact nonunit height data.  It is the
height-flex leak predicted by HYP-3310/HYP-3402.

## Stronger Leak: Endpoint Owner

By `(30,8)`, a stronger collision appears:

```text
petal 13->26
single swap 1->26
single swap 3->26
single swap 5->26
```

These rows share:

```text
same coarse sheaf base
same nonunit residue word
same nonunit v2 word
same exact nonunit height word
```

but split theorem exits:

```text
petal 13->26           -> unit-petal-named
single swap *->26      -> positive-Haar-open
```

So the next missing coordinate is not nonunit height.  It is the endpoint
owner pattern carried by the boundary layer.

This family grows on larger banks:

```text
(48,12): add single swaps *->40
(60,16): add single swaps *->54 and 9->54
```

The same phenomenon persists: residue, `v2`, and exact nonunit height remain
blind, while endpoint-owner support separates the exits.

## Why Owner Support Works

The owner-support word is extracted from HYP-2975's taut-bridge boundary layer:
the set of endpoint-owner labels that support positive bridges or taut vertices.

Examples:

```text
petal 13->26      owner_support = (12:g2, 1:g1)
single swap 1->26 owner_support = (12:g2, 13:g1)
single swap 3->26 owner_support = (11:g1, 12:g2, 13:g1, 6:g2)
P10+GW            owner_support = (6:g2, 7:g7)
GW-shell 12->132  owner_support = (11:g1, 13:g1, 5:g1, 6:g2, 7:g7)
```

This is exactly the coordinate that the residue-only packet forgot: not merely
which covering residues exist, but which endpoint-owner channels are active at
the positive boundary.

## Interpretation

The enlarged-bank continuation of HYP-3311 is now concrete:

```text
curated bank:
  residue word is exact

expanded bank:
  first failure = height-driven
  next persistent failure = endpoint-owner driven
  residue + owner-support is exact on the scanned banks
```

This sharpens HYP-3402.  The true order of next-sidecar repairs is:

```text
1. nonunit residue word is only bank-local;
2. some failures are repaired by height/v2;
3. the stronger surviving failures require endpoint-owner support.
```

## Caveat

This is exact enlarged-bank evidence, not a proof that owner support is
globally terminal.  The next task is to enlarge the bank further and test
whether:

```text
residue + owner_support
```

eventually fails, and whether the next leak after that is tropical/off-grid,
unit-contact holonomy, or a new named debt.

## Tournament Analysis

Use sidecar repairs as vertices:

```text
residue_plus_owner_support
residue_plus_unit_qsqrt7
residue_plus_v2
residue_plus_height
residue_plus_unit_slot
residue_only
```

Pairwise observable:

```text
how many mixed theorem-exit fibers survive on the enlarged packet banks
```

Exact fingerprint:

```text
directed_3cycles = 0
hamiltonian_path_count = 1
priority_path =
  residue_plus_owner_support
  -> residue_plus_unit_qsqrt7
  -> residue_plus_v2
  -> residue_plus_height
  -> residue_plus_unit_slot
  -> residue_only
```

The point is not that owner support is "smallest".  It is that it is the first
repair that stays exact once the bank leaves the curated sample.
