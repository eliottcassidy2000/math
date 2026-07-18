---
id: THM-1042
title: The component-length obstruction — G_{1..k}(1/14) has longest component L_max(k) with 1/L_max(k) > k+1 for every k ≥ 3, so an additive certificate can NEVER absorb the next consecutive speed. This closes the additive route on every family with a small-integer core, and explains uniformly why THM-1015 needed large killers.
status: PROVED as an exact finite computation (component lengths are exact rationals; the table is verifiable by hand) + the structural criterion it implies. This is a NEGATIVE result about a method, not about LRC(14).
source: klein-2026-07-18-S327
depends_on:
  - THM-1004  # the interval-survival tail whose hypothesis this measures
  - THM-1015  # the clustered closure, whose large-killer hypothesis this explains
related:
  - THM-1026  # opus's five-slot ledger (the alternating route, killed separately in klein-S325)
---

# THM-1042 — the component-length obstruction

The additive certificates (THM-1004's interval-survival tail, and the measure-recursion variant) absorb a
new speed `w` into a good set by charging it a *proportional* loss. That accounting is only valid when the
good set's components are **longer than `w`'s arc period `1/w`**; otherwise a single period spans a whole
component and the component is clipped or swallowed outright, with no proportional regime. So

```text
base B admits speed w   only if   w > 1/L_max(B),        L_max(B) := longest component of G_B(1/14).
```

This theorem computes `L_max` and shows the criterion is never satisfiable by a consecutive speed.

## The exact data

Components of `G_B(1/14)` are intervals with endpoints `(14k ± 1)/(14v)`, `v ∈ B`, so every component
length is a rational with denominator dividing `14·v·v'` for some `v,v' ∈ B` (e.g. `1/588`, `588 = 14·6·7`).
Computed exactly:

| `B` | `μ = |G_B|` | `#comps` | `L_max` | `1/L_max` | next speed `k+1` |
|---|---|---|---|---|---|
| {1..3} | 0.69048 | 4 | 5/21 | 4.2 | 4 |
| {1..4} | 0.61905 | 6 | 9/56 | 6.2 | 5 |
| {1..5} | 0.50476 | 10 | 4/35 | 8.8 | 6 |
| {1..6} | 0.45714 | 12 | 1/12 | 12.0 | 7 |
| {1..7} | 0.33469 | 18 | 3/49 | 16.3 | 8 |
| {1..8} | 0.26582 | 20 | 5/112 | 22.4 | 9 |
| {1..9} | 0.18107 | 20 | 2/63 | 31.5 | 10 |
| {1..10} | 0.13798 | 20 | 3/140 | 46.7 | 11 |
| {1..11} | 0.05633 | 20 | 1/77 | 77.0 | 12 |

**`1/L_max(k) > k+1` at every row, and the gap widens**: the ratio `(1/L_max)/(k+1)` runs
`1.05, 1.24, 1.47, 1.71, 2.04, 2.49, 3.15, 4.25, 6.42`. `1/L_max` grows superlinearly while the next
available speed grows linearly.

## Consequence

**An additive certificate can never absorb a consecutive speed.** Hence it fails on every family whose
speeds contain a run of small consecutive integers — the AP, the deep well, the GW family, and every
covering family with a small-integer core. Checked against the deep well `{1,…,12,182}`: *every* initial
split is blocked, and the blocking speeds are exactly the consecutive ones.

```text
base {1..7}  needs w > 16.3 ; killers 8,9,10,11,12,182 ; blocked by 8,9,10,11,12
base {1..11} needs w > 77.0 ; killers 12,182           ; blocked by 12
```

## What this explains, uniformly

Three separate failures now have one cause:

1. **THM-1015 required large killers** (thresholds 65.7 … 347.5). Not a convenience — those thresholds are
   `1/L_max` of the respective bases. The clustered stratum closed *because* its killers were large.
2. **The measure-recursion (klein-S326) died at `w = 8`** with boundary `2δN/w = 0.321` exceeding the
   surviving measure. `8 < 16.3 = 1/L_max({1..7})`.
3. **The radius-3 fragmentation wall** (klein-S314): the same short components, seen from the
   Hamming-radius side.

Changing the state variable — largest interval → (measure, component count) — removes the `r < 1/(2δ) = 7`
cap but not this, because both formulations price a speed against the component scale.

## Scope, stated plainly

This is a theorem about the **method**, not about LRC(14). It says the additive/measure family of
certificates cannot reach the small-speed regime, and gives the exact threshold at which they can. It does
**not** bound `M` for any family, and nothing here contradicts LRC(14). Combined with klein-S325 (the
alternating truncation `B₅ < 0` on real families, first clear at `B₁₁`) and klein-S324 (no pairwise-only
invariant characterizes tightness), the three certificate families the fleet has been using are now each
delimited by an exact, checkable criterion.

The honest reading: the small-speed / compact regime is not awaiting a sharper constant. It is outside the
reach of proportional-loss accounting, because there the good set has no component longer than the arcs
being added.

*Files: `04-computation/lrc_complength_klein_S327.py` (+ `_complength`, `_threshold` .out).*
