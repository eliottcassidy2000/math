---
id: HYP-3422
title: LRC14 two-adic off-grid relocation lemma
status: EVIDENCE / exact interval scout; not an LRC14 proof
source: codex-2026-06-28 response to resonant-transparency prompt after HYP-3418
tangent: T1383
technique: LTI-383
tournament_technique: LTT-283
script: 04-computation/lrc14_two_adic_offgrid_relocation_codex_20260628.py
result: 05-knowledge/results/lrc14_two_adic_offgrid_relocation_codex_20260628.out
reflection: 07-reflections/lrc14-two-adic-offgrid-relocation-codex-20260628.md
related:
  - HYP-3421
  - HYP-3420
  - HYP-3419
  - HYP-3418
  - HYP-3417
  - HYP-3416
  - HYP-3415
  - HYP-3410
  - HYP-3409
  - HYP-3408
  - HYP-3407
  - HYP-3406
  - HYP-3129
  - HYP-2963
  - THM-523
  - OPEN-Q-108
---

# HYP-3422: LRC14 Two-Adic Off-Grid Relocation Lemma

## Claim

The resonant-transparency observation should be recast as a two-adic
relocation lemma, not as a coprime-to-14 reduction.

Split a covering row as

```text
S = O union 2E
```

with `O` odd and `2E` even.  Put `u = 2t`.  The even speeds are safe at `t`
exactly when the halved set `E` is safe at `u`.  The two lifts of `u` give two
odd filters:

```text
t = u/2       requires ||o*u/2|| >= 1/14 for every odd o
t = (u+1)/2   requires ||o*u/2|| <= 3/7 for every odd o
```

So the exact finite target is:

```text
E_safe(1/14) intersects odd_branch_0_good union odd_branch_1_good.
```

This keeps the useful part of the user's prompt: the full optimizer moves
off-grid so resonant speeds do not obstruct at the grid.  It also keeps
HYP-3418's correction: the nonresonant/coprime witness itself is usually
`t=1/2`, where every even speed dies.

Incoming HYP-3421 is the companion full-optimum transparency ledger and
Rprime-closure route.  HYP-3422 is narrower: it names the exact two-adic
half-lift interval statement that should make HYP-3418's even-speed descent
rigorous.  Incoming HYP-3420 remains the owner-cut/chiral sidecar lane for
finite exceptional packets.

## Executable Readout

Script:

```text
04-computation/lrc14_two_adic_offgrid_relocation_codex_20260628.py
```

Stored output:

```text
05-knowledge/results/lrc14_two_adic_offgrid_relocation_codex_20260628.out
```

Exact audit on `24` rows (`4` curated covering rows plus `20` deterministic
primitive covering rows with speed cap `120`):

```text
full exact M(S) >= 1/14:              24/24
full optimizer selected off 14-grid:  24/24
naive nonresonant witness fails:      24/24
branch-0 relocation certificates:     24/24
branch-1 relocation certificates:     24/24
either-branch certificates verified:  24/24
```

The tightest displayed relocation score is still above the LRC14 threshold:

```text
covering_AP_with_12_and_84:
  selected branch = 1
  t = 136/147
  score = 11/147 > 1/14
```

The smallest branch-union row in the scout is the canonical covering row
`(1,2,3,4,5,6,7,8,9,10,11,13,84)`, with exact branch measures
`563/105105` for each lift branch.  This is small but positive, which is the
kind of finite-ruler lower bound the proof needs.

## Proof Meaning

The old slogan:

```text
resonance puts danger on the grid; Q-lonely lives off-grid
```

is true only after choosing the correct lift.  In the `u=2t` formulation,
branch `1` is the off-grid lift near `t=1/2`; it keeps odd speeds near
half-integers while the even speeds solve the halved problem.

This suggests a proof route:

1. Use LRC<=13, or a sharper covering floor, on the halved even packet `E`.
2. Prove that the resulting safe set in `u` is not fully covered by the odd
   branch-bad intervals.
3. Iterate when `E` still has even speeds, producing a two-adic descent.
4. Stop at odd-only rows, where the branch near `t=1/2` is naturally safe.

This is not yet a proof because step `2` is a real interval-overlap theorem.
But it is more precise than "nonresonant decorrelation": it names the exact
bad intervals and the exact branch lift.

## Tournament Analysis

Vertices are proof obligations, not runners or constants.

```text
score_hist={-30:1, 36:1, 52:1, 61:1, 65:1, 72:1, 74:1, 85:1}
directed_3cycles=0
hamiltonian_path=
  R00_two_adic_lift_identity
  -> R01_branch_one_odd_tolerance
  -> R02_even_half_descent
  -> R03_owner_current_cut_sidecar
  -> R04_nonresonant_decorrelation_floor
  -> R05_resonant_grid_transparency
  -> R06_apex7_offpath_guardrail
  -> R07_raw_named_analogy
```

The ordering is intentional: the exact two-adic lift identity and branch
tolerance should be proved before using apex-7/Galois/census structure, which
HYP-3418 now places off the covering-floor critical path.

## Next Hook

Replace the sample audit by a theorem over interval families:

```text
For every primitive covering 13-set S = O union 2E,
measure(E_safe cap branch_good) > 0.
```

The most promising first bound is not a global scalar average.  It is a
finite-ruler / Helly-style interval statement on `E_safe` after removing the
odd branch-bad intervals, with HYP-3417 owner-current labels used only to name
finite exceptional packets.
