---
id: THM-1132
title: Complete search-free proof that the sharp measure horn certifies the r=6 worst shape — R_sharp(b) < 1 for the core {1,2,4,7,9,11,12} with consecutive killers {b,…,b+8}, all b ≥ 157
status: PROVED (rigorous). Finite part = exact rational check 157≤b≤399; tail part b≥400 = an explicit rational witness with exact margins. No search, no floating point in the certificate. This closes the "alignment" and "error-control" gaps THM-1123 left open, FOR THE WORST CORE + consecutive-killer shape (the empirical maximizer of R_sharp at r=6). It does NOT yet cover other cores/shapes.
source: death-star-2026-07-18-S58 (continued)
depends_on:
  - THM-1123   # the sharp horn 1/(7L); this discharges its worst-shape uniform bound
related:
  - THM-1102   # the conservative-constant horn whose "wall" this dissolves
script: 04-computation/r6_worstcore_finite_exact_deathstar_S58.py
output: 05-knowledge/results/r6_worstcore_finite_exact_deathstar_S58.out
---

# THM-1132 — the r=6 worst shape is closed with no search

## Statement

Let `P₀ = {1,2,4,7,9,11,12}` (the empirical R_sharp-maximizing 7-core) and, for an integer
`b ≥ 157 = 13·max(P₀)+1`, let `K_b = {b, b+2, b+4, b+6, b+8}` (consecutive step-2 killers —
the empirical worst shape). Let `L(b)` be the length of the largest component of
`G_b = S(P₀) ∖ ⋃_{k∈K_b} danger(k)`, where `danger(k)={t:‖kt‖<1/14}` and
`S(P₀)={t:‖pt‖≥1/14 ∀p∈P₀}`.

> **THM-1132.** `L(b) > 1/(7(b+8))` for every integer `b ≥ 157`.

Consequently, for any 6th (dodged) killer `k₆ > b+8`, `L(b) > 1/(7k₆)`, so the sharp measure
horn (THM-1123) exhibits a real `1/14`-lonely time: **R_sharp(b) < 1 for all b ≥ 157**, with
no search and no floating-point in the certificate.

## Proof

**Part A — finite range `157 ≤ b ≤ 399` (exact).** `L(b)` is computed in exact rational
arithmetic (`r6_worstcore_finite_exact`, cross-checked by an independent full-circle exact
routine). Result: `R_sharp(b) = 1/(7 L(b)(b+8)) < 1` for all 243 values, with maximum
`R_sharp = 10325/12888 = 0.801133 at b = 171` (`L(171) = 72/72275`). ∎(A)

**Part B — tail `b ≥ 400` (explicit rational witness).** Use the core-safe arc

  `A* = (71/154, 13/28)`, a component of `S(P₀)` of width `13/28 − 71/154 = 1/308`,

whose left endpoint is the speed-11 danger boundary and right endpoint the speed-2 danger
boundary (so **no** core speed has a danger arc in the interior of `A*`). Note
`2·A* = (71/77, 13/14) ⊂ (6/7, 1)`, a band on which the AP-arc gap function `G(σ) > 1/7`
(THM-1123's lemma).

Fix `b ≥ 400` and set `ρ = 9/(77b)`. Define
`A_good(b) = {t∈A* : ‖pt‖ − 1/14 ≥ p·ρ  ∀ p∈P₀}`. Because the only core constraints on `A*`
are the two endpoints (slopes 11 and 2), moving in by `ρ` from each endpoint clears them:

  `A_good(b) ⊇ (71/154 + ρ, 13/28 − ρ)`,  so  `|A_good(b)| ≥ 1/308 − 18/(77b)`.

For `b ≥ 400`, `(b+8)(1/308 − 18/(77b)) ≥ 1` (value `1.086` at `b=400`, increasing in `b`),
so the scaled interval `(b+8)·A_good(b)` has length `≥ 1` and therefore contains a
half-integer `n + 1/2`. Put

  `t* = (n + 1/2)/(b+8) ∈ A_good(b)`.

*Killers.* `‖(b+8)t*‖ = ‖n+1/2‖ = 1/2`. For `m=0,…,3`,
`(b+2m)t* = (n+1/2) − (8−2m)t*`, so `‖(b+2m)t*‖ = ‖(8−2m)t* − 1/2‖`. Since `t*∈A*`,
`8t* mod 1 ∈ (53/77, 5/7)`, giving `‖bt*‖ ∈ [29/154, 3/14]`; the other three are strictly
larger. Hence the smallest killer margin is
`(‖bt*‖ − 1/14)/b ≥ (29/154 − 11/154)/b = 9/(77b) = ρ`.

*Core.* `t*∈A_good(b)` gives `(‖pt*‖ − 1/14)/p ≥ ρ` for every `p∈P₀`.

Therefore every one of the 13 speeds `v ∈ P₀ ∪ K_b` satisfies `‖vt*‖ ≥ 1/14 + vρ`, so the
interval `[t*−ρ, t*+ρ]` avoids every danger arc: it is `1/14`-safe for all 13 speeds. Thus

  `L(b) ≥ 2ρ = 18/(77b)`,  and  `18/(77b) > 1/(7(b+8)) ⟺ 126b + 1008 > 77b`, always true. ∎(B)

Parts A and B together prove `L(b) > 1/(7(b+8))` for all `b ≥ 157`. ∎

## What this closes (and what it does not)

- **Closed:** THM-1123's two flagged gaps, *for this core and shape*. The "alignment" gap is
  discharged by the explicit `A* = (71/154,13/28)` with `2·A* ⊂ (6/7,1)`; the "error-control"
  gap is discharged by replacing the `G(σ)` approximation with an **exact rational witness**
  `t*` carrying explicit margins `9/(77b)` — no `σ`-drift estimate is needed. The construction
  is verified failure-free for `400 ≤ b ≤ 4000` and holds for all `b ≥ 400` by the length bound.
- **Mechanism confirmed:** `L(b)·b → 4/7` as `b→∞` (large killers: the 4 non-origin arcs cluster
  at `σ≈0.928`, leaving a 4/7-period gap), so `R_sharp(b) → 1/4`; the maximum `0.8011` is a
  small-`b` (`b=171`) phenomenon, which is why a finite check plus a generous tail suffices.
- **Not yet done:** the other 791 cores and the non-consecutive killer shapes. Each is the same
  two-part template (exact finite head + explicit-witness tail in a core-safe arc whose double
  lies in a good `G`-band); the tail witness generalizes, but the good-band arc and the finite
  bound `B₀` are core-specific. Automating "find a core-safe arc `A` with `2A` in a good band and
  run the two-part proof" over all 792 cores is the path to the full r=6 uniform theorem.

This is the first **search-free** closure of any part of the r=6 covering-killer stratum.
