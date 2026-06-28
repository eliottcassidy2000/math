---
id: HYP-3416
title: LRC14 recursive quotient ladder
status: SYNTHESIS / executable recursion-pattern atlas; not an LRC14 proof
source: codex-2026-06-28
tangent: T1375
technique: LTI-375
tournament_technique: LTT-275
script: 04-computation/lrc14_recursive_quotient_ladder_codex_20260628.py
result: 05-knowledge/results/lrc14_recursive_quotient_ladder_codex_20260628.out
reflection: 07-reflections/lrc14-recursive-quotient-ladder-codex-20260628.md
related:
  - HYP-3410
  - HYP-3409
  - HYP-3408
  - HYP-3407
  - HYP-3406
  - HYP-3405
  - HYP-3404
  - HYP-3403
  - HYP-3402
  - HYP-3401
  - HYP-3311
  - HYP-3301
  - HYP-3234
  - HYP-3231
  - HYP-3230
  - HYP-3244
  - HYP-3124
  - HYP-3123
  - HYP-3039
  - HYP-3023
  - HYP-2963
  - OPEN-Q-108
---

# HYP-3416: LRC14 Recursive Quotient Ladder

## Claim

The recent LRC14 proof stack is repeatedly discovering one recursive law:

```text
compress a labelled packet
-> declare the LRC predicate preserved by the compression
-> audit theorem-exit purity on compressed fibers
-> if a fiber mixes exits, name the first forgotten coordinate
-> add that sidecar, cut, child deck, local disk, chart change, or debt
-> recurse on the smaller/purer packet.
```

This should be treated as a proof target, not just as an organizing slogan.

Concurrent-mainline bridge: HYP-3408 is the exact guardrail-sidecar router,
HYP-3409 is the tighter recursive sidecar pattern atlas, and HYP-3410 is the
charal/Menger mixed-fiber slice.  HYP-3416 is the wider quotient-ladder view
that folds those into older zippers, tail/tip decks, scale-normal descent,
tiling lifts, signed chart changes, analytic ledgers, and branch-sheet
payloads.  Its lesson is that after a quotient fails the theorem-exit purity
test, the first illegal forgetting event becomes the next packet coordinate.

> Recursive Packet Ladder Lemma target:
> every primitive LRC14 quotient used in the final proof is either
> theorem-exit pure, or its first impurity is discharged by owner support,
> a unit-height disk exit, tail/tip child recursion, a scale-normal sidecar,
> a chart-change sidecar, an analytic exceptional ledger, a branch-sheet
> payload, or a named finite base case.

The strongest new abstraction is that the recursive vertex is often not a
runner, an arc, a residue, or even a row.  The vertex is the first legal
proof-compression state.

Local side-signal preserved during close-out: the untracked
`lrc_c6_is_bu_times_index_macmini_S85.py` / `.out` pair verifies the
`C6 = C2 x C3` residue reading of the binding layer and explicitly warns that
the index-theorem/C6 quotient still forgets the 2-adic height.  That is exactly
the HYP-3416 guardrail in miniature: a beautiful quotient is useful only after
the magnitude layer it destroys is named as a sidecar.

## Executable Atlas

The script

```text
04-computation/lrc14_recursive_quotient_ladder_codex_20260628.py
```

ranks thirteen recursion carriers:

```text
R01 119  Controlled-forgetting purity ladder
R02 115  Residue-height-owner ladder
R03 113  Endpoint-owner Menger descent
R04 105  AP-collar local disk recursion
R05 104  Tail-tip child-deck recursion
R06  97  Scale-normal primitive descent
R07  94  Fiber zipper magnitude-cocycle ladder
R08  88  Tiling / half-tiling interlock
R12  81  Chiral mirror recursion
R09  81  Signed chart-change recursion
R10  77  Three-gap / Stern-Brocot cap recursion
R11  69  Analytic exceptional-fiber recursion
R13  57  Branch-sheet / resolvent recursion
```

Tournament fingerprint:

```text
vertices=13
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1,11:1,12:1}
directed_3cycles=0
hamiltonian_path_count=1
priority_path=R01 -> R02 -> R03 -> R04 -> R05 -> R06 -> R07 -> R08 -> R12 -> R09 -> R10 -> R11 -> R13
```

The pairwise observable is:

```text
score, theorem purity, first-failure power, HYP-3406 bridge,
recursion depth, computability, low payload cost, low external risk.
```

## Patterns Abstractly Noticed

### 1. The Quotient Itself Recurred

HYP-3039, HYP-3301, HYP-3404, HYP-3407, HYP-3408, HYP-3409, HYP-3410, and now
HYP-3416 keep saying the same thing in different dialects: a quotient is
proof-legal only if `exit_status` is a function of the quotient.  Any mixed
fiber is not noise; it is the next coordinate.

### 2. The HYP-3406 Ladder Is The Current Data Spine

The newest exact data spine is:

```text
coarse sheaf -> residue -> v2/height -> owner_support.
```

At `(single_limit,two_swap_limit)=(72,20)`, residue and height still leave
owner-driven mixed fibers, while `residue + owner_support` kills them.  This
makes endpoint owner the first live recursive child after residue/height.

### 3. Menger Cuts Are The Owner Recursion Made Geometric

Endpoint-owner support should not remain a word forever.  It should become a
cut object:

```text
boundary channels -> min cut -> child packet on each shore.
```

The immediate test is the HYP-3406 fiber containing petal `13->26`, the
positive-open swaps into `26`, later swaps into `40/54/68`, and the
`petal 10->20` two-drop/add-20 frontier.

### 4. HYP-3405 Is The Local Base Case

The AP-collar leak is not just another counterexample to residue-only data.
It is the local disk base case of the ladder:

```text
AP/GW boundary disk -> strict-open collar -> first unit-height exit.
```

The missing coordinate is the unit-height lift `(13,0)->(13,1)`.  This is the
template for what a "first exit" should look like in the larger owner graph.

### 5. Tail/Tip Deletion Is Endpoint Owner In Tournament Form

HYP-3124's edge recursion says an edge is a two-ended obligation:

```text
edge -> tail child + tip child -> recursive child edge decks.
```

HYP-3416 proposes attaching this directly to owner-support words.  If an
owner-support quotient is exact, its tail/tip deletion children should either
remain theorem-exit pure or name which child carries the next obstruction.

### 6. Scale, Zippers, Tiling, And Charts Are The Same Move

HYP-3231 scale-normal recursion, HYP-3023 automatic-fiber zippers, HYP-3244
tiling/half-tiling interlocks, and HYP-3234 signed chart changes all implement
the same pattern:

```text
lift to a richer cover,
compress to a quotient,
record what the quotient destroys,
recurse only after sidecars restore theorem-exit purity.
```

The apparent variety is real, but the proof protocol is shared.

### 7. Analytic Estimates Should Recurse Only After Exceptions Are Named

BDH/Mertens-style mean-square language belongs after the fiber structure is
known.  It should average over pure packets or over packets whose exceptional
fibers have names.  Otherwise the smoothing function may hide exactly the
owner, height, period, or branch coordinate that changes theorem exit.

## Assumption Challenge

Alternate tournament vertices considered:

```text
runners, arcs, gaps, residue classes, cover arcs, endpoint owners,
fixed circle sections, wall-crossing events, Fourier modes, child decks,
Menger cuts, local disks, quotient fibers, branch sheets, scale orbits,
chart changes, analytic exceptional fibers, proof obligations.
```

Chosen vertices:

```text
recursive proof-compression carriers.
```

Preserved LRC predicate:

```text
theorem-exit purity for boundary-tight / strict-open / positive-Haar-open /
unit-petal / K33-H7 / named residual exits, hence eventual M(S) >= 1/14.
```

Destroyed information if quotients are used naively:

```text
unit height, exact magnitude, endpoint owner, tail/tip role, scale depth,
path presentation, automorphism fiber, signed chart basis, exact period,
continued-fraction address, smoothing exception identity, branch sheet,
and chiral mirror side.
```

## Next Tests

1. Attach HYP-3124 tail/tip child decks to HYP-3406 `owner_support` fibers.
2. Build the HYP-3406 owner-support Menger graph and compute cut shores.
3. Treat HYP-3405 AP-vs-`13->27` as the local-disk base case.
4. Re-run fiber-purity counts after adding chiral mirror and cut-shore labels.
5. Only then attempt BDH/Mertens averaging over pure or named-exception fibers.
