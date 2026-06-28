---
id: HYP-3431
title: LRC14 canonical corridor-fence certificate
status: PARTIAL PROOF / all-m canonical two-branch relocation certificate; not an LRC14 proof
source: codex-2026-06-28; renumbered after mainline HYP-3428/HYP-3429/HYP-3430 claims
tangent: T1392
technique: LTI-392
tournament_technique: LTT-292
script: 04-computation/lrc14_canonical_corridor_fence_codex_20260628.py
result: 05-knowledge/results/lrc14_canonical_corridor_fence_codex_20260628.out
reflection: 07-reflections/lrc14-canonical-corridor-fence-codex-20260628.md
related:
  - HYP-3430
  - HYP-3429
  - HYP-3428
  - HYP-3427
  - HYP-3426
  - HYP-3425
  - HYP-3424
  - HYP-3423
  - HYP-3422
  - HYP-3421
  - HYP-3418
  - HYP-3415
  - HYP-3140
  - HYP-3129
  - HYP-2272
  - THM-523
  - OPEN-Q-108
---

# HYP-3431: LRC14 Canonical Corridor-Fence Certificate

## Claim

For the canonical resonant one-tail tower

```text
S_m = {1,2,3,4,5,6,7,8,9,10,11,13,84m},
```

HYP-3425's two-branch relocation target, sharpened by HYP-3426's one-branch
mirror audit, HYP-3427's wall-signature atlas, HYP-3428's two-adic loss ledger,
HYP-3429's endpoint-spine certificate, and HYP-3430's harmonic-intercept
firewall, is provable for every `m >= 1` by a fixed-corridor / moving-grid
argument on this tower.

Write `S_m = O union 2E` and `u=2t`.  The low core

```text
O={1,3,5,7,9,11,13},    E_low={1,2,3,4,5}
```

leaves exactly two fixed branch-good corridors:

```text
[8/49, 6/35]      length 2/245
[29/35, 41/49]    length 2/245
```

The only moving obstruction in the tower is the high even half-speed
`N=42m`.  Its bad intervals are disjoint grid intervals of width
`1/(7N)=1/(294m)`, with gaps `1/(49m)`.  Since

```text
2/245 > 1/(294m)       for every m >= 1,
```

no connected fixed corridor can be covered by the disjoint high-grid bad
intervals.  Therefore every `S_m` has a positive two-branch relocation window.

Relative to HYP-3428, this is the clean canonical ledger case: the halved even
child has one moving high-grid packet, the odd blockers are exactly the fixed
endpoint walls `{5,7}`, and no owner-current/even-hinge debt remains after the
component-width inequality.

Relative to HYP-3429, this proves the canonical endpoint-spine pattern for all
`m`: the survivor components are bounded by the fixed odd walls `5,7` together
with the moving even wall `E:84m`, and after the high-grid fence opens a
component the rank is at most the endpoint-spine rank seen in the finite audit.

Relative to HYP-3430, this is the exact sidecar the harmonic firewall demands:
the proof uses wall labels and component widths, not the scalar
`H_N - log N` tail intercept that only calibrates scale.

This proves the all-`m` canonical tower certificate, but it does not prove the
full LRC14 covering floor.

## Evidence

Script:

```text
04-computation/lrc14_canonical_corridor_fence_codex_20260628.py
```

Stored output:

```text
05-knowledge/results/lrc14_canonical_corridor_fence_codex_20260628.out
```

The executable scout verifies the low-core corridor identity exactly, then
audits the closed-form tower through `m<=40` and checks the symbolic margin
through `m<=400`.

Key readout:

```text
fixed corridors: [8/49,6/35], [29/35,41/49]
fixed corridor length: 2/245
high bad width: 1/(294m)
symbolic margin: (12m-5)/(1470m) > 0
proof-margin failures for m<=400: []
worst good measure for m<=40: 1/105 at m=1
```

Endpoint owners:

```text
left corridor:  B1 odd 7 wall -> B1 odd 5 wall
right corridor: B0 odd 5 wall -> B0 odd 7 wall
```

The endpoint ownership matters.  A scalar measure table hides why the
canonical tower is easy: the low odd/even core creates fixed corridors, and
the high resonant speed only punches a separated grid through them.

## Proof Route

The canonical tower no longer needs a finite Helly search:

1. Prove the low-core corridor identity once.
2. Prove the high-grid bad intervals are pairwise disjoint with component
   width `1/(294m)`.
3. Use the connected-interval fence lemma: a connected interval of length
   greater than every component of a disjoint union cannot be covered by that
   union.
4. Lift the surviving `u` point to `t=u/2` or `t=(u+1)/2`, obtaining
   `M(S_m) >= 1/14`.

This is a proof-facing base lemma for
HYP-3430/HYP-3429/HYP-3428/HYP-3427/HYP-3426/HYP-3425/HYP-3422.  The next
generalization is to
search for non-canonical rows whose low core leaves a fixed or slowly moving
corridor longer than the widest remaining moving bad component.  Rows that fail
that test should return to the HYP-3430 scalar firewall, the HYP-3429
endpoint-spine lemma, the HYP-3428 loss ledger, the HYP-3425 component-gap
Helly certificate, or owner-current exception routing.

After the HYP-3425 additive-energy sidecar and HYP-3428 descent ledger, the
non-canonical search should also keep sheet and loss-class data.  A row that
fails the corridor-fence test should not be summarized by raw `fullE`, `RE`, or
`oddE`; it should carry an
energy-plus-sheet packet such as `(RE,q_zero_mass)` or `(fullE,q_range_hi)`
before feeding SPEC, phase-cover debt, or a named terminal exception.

## Tournament Analysis

Vertices are proof carriers and wall certificates, not runners or raw
intervals.

```text
score_hist={28:1,54:1,56:1,58:1,59:1,61:1,62:1}
directed_3cycles=0
hamiltonian_path=
  high_grid_fence_lemma
  -> canonical_84m_all_m_certificate
  -> two_branch_helly_generalization
  -> fixed_low_corridor_identity
  -> endpoint_wall_ownership_dictionary
  -> owner_current_exception_router
  -> raw_measure_table
```

Assumption challenge: runners, fixed corridors, corridor walls, high-grid wall
events, odd branch walls, even half-speeds, owner labels, Fourier modes, and
proof obligations were considered.  The chosen quotient preserves the
two-branch relocation predicate for the canonical tower and destroys endpoint
ownership if collapsed to a scalar measure.
