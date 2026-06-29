---
id: HYP-3456
title: LRC14 AP84 mod-35 floor-count closure
status: EVIDENCE / floor-count derivation of AP-tail escape clock; not an LRC14 proof
source: codex-2026-06-29 continuation of HYP-3431, HYP-3452, HYP-3454, and HYP-3453
tangent: T1416
technique: LTI-416
tournament_technique: LTT-316
script: 04-computation/lrc14_ap84_mod35_floor_count_codex_20260629.py
result: 05-knowledge/results/lrc14_ap84_mod35_floor_count_codex_20260629.out
reflection: 07-reflections/lrc14-ap84-mod35-floor-count-codex-20260629.md
related:
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
  - HYP-3433
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

# HYP-3456: LRC14 AP84 Mod-35 Floor-Count Closure

## Claim

HYP-3456 closes the sampled mod-`35` clock left by HYP-3454 into an exact
floor-count formula over the HYP-3431 low corridors.

For

```text
S_m = {1,2,...,11,13,84m},
```

HYP-3431 gives the two fixed low branch-union corridors

```text
C1=[8/49,6/35]
C0=[29/35,41/49].
```

The moving even half-speed `42m` has safe gaps

```text
G_k(m)=[(14k+1)/(588m),(14k+13)/(588m)].
```

An AP-tail escape component is counted by a high safe gap intersecting one of
the fixed corridors.  For `C1`, the intersection inequalities are

```text
(14k+13)/(588m) > 8/49
(14k+1)/(588m)  < 6/35,
```

so the count per corridor is

```text
N(m)=floor((504m-6)/70)-floor((96m-13)/14).
```

The second corridor is its mirror, hence

```text
escapes(m)=2*N(m).
```

This is the period-`35` clock from HYP-3452/HYP-3454, now derived from fixed
corridor endpoints and the moving `E:84m` grid.

## Exact Readout

Script:

```text
04-computation/lrc14_ap84_mod35_floor_count_codex_20260629.py
```

Stored output:

```text
05-knowledge/results/lrc14_ap84_mod35_floor_count_codex_20260629.out
```

Validation against the exact HYP-3452 component audit through `m=70`:

```text
mirror_failures=[]
formula_failures=[]
component_audit_failures=[]
```

The correction vector is

```text
[2,2,1,1,1,1,1,2,1,1,2,1,1,1,1,2,2,1,1,1,1,2,2,1,1,2,1,1,2,1,2,2,1,1,0].
```

and

```text
N(m)=floor(12m/35)+d[(m-1) mod 35].
```

Because

```text
floor((504(m+35)-6)/70)=floor((504m-6)/70)+252
floor((96(m+35)-13)/14)=floor((96m-13)/14)+240,
```

we get

```text
N(m+35)-N(m)=12
escapes(m+35)-escapes(m)=24.
```

The script checks this shift for `m=1..210` with no failures.

## Proof Pull

HYP-3454 reduced the AP-tail endpoint phase to:

```text
finite transients m=1..4
closed-form endpoint interval for m>=5
mod-35 boundary-count lemma
```

HYP-3456 supplies the boundary-count lemma, conditional on the HYP-3431 fixed
low-corridor identity.  The AP-tail bridge is now down to:

```text
1. use HYP-3431 to keep C1 and C0 as the complete low branch-union carrier;
2. use HYP-3454 for the rank-one endpoint interval inside those corridors;
3. use HYP-3456 to count every moving high-grid escape;
4. handle m=1..4 as finite mixed E/B1 transients.
```

This does not prove arbitrary primitive rows.  It closes the AP84 tail clock
needed by HYP-3439/HYP-3454 and hands HYP-3453/HYP-3451 a cleaner local gate
input.

Rebase integration: incoming HYP-3455 is the noncanonical
`random_covering_031` gate-gluing clause.  HYP-3456 is the canonical AP-tail
floor-count sibling: use HYP-3455 for the named rank-`6` noncanonical gate
case and HYP-3456 for the AP84 residue clock.

## Tournament Analysis

Vertices are floor-count proof obligations, not runners or raw components.

```text
pairwise_observable =
  fixed-corridor retention + high-gap count exactness + residue payload
  + audit validation
score_hist={13:1,28:1,47:1,52:1,55:1,57:1,59:1}
directed_3cycles=0
hamiltonian_path =
  gap_intersection_floor_count
  -> fixed_low_corridor_identity
  -> period35_residue_vector
  -> mirror_two_corridor_doubling
  -> component_audit_validation
  -> raw_beatty_fit
  -> raw_dead_fraction_peak
```

Assumption challenge: runners, `m`-values, residues mod `35`, high-grid gaps,
fixed low corridor sections, endpoint walls, component-cover graph nodes, and
proof obligations were considered.  The chosen quotient preserves the exact
AP-tail escape count and endpoint-wall clock, but destroys non-AP wall geometry
and dead-cover adjacency.  It is an AP-tail floor-count lemma, not a global
primitive-row theorem.
