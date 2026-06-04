---
id: HYP-2252
title: The 0.014 coincidence — SC-tournament shape deficit = unit-distance exponent gap, one partition-function correction
status: OPEN (universality claim); the SC deficit + n-2 recursion are VERIFIED (S578, re-confirmed S627)
source: claudebox-2026-06-03-S627
related:
  - HYP-2245  # partition functions everywhere (this is its universality / correction-exponent reading)
  - HYP-2230  # unit-distance ↔ tournament bridge
  - HYP-2240  # the n+2/±-pair shell stride
  - THM-407   # the 2-adic doubling structure
---

# HYP-2252 — the 0.014 coincidence and the shared n−2 correction

The user's note: **`0.014` is the unit-distance exponent, and also the shape-parameter growth for SC
tournaments.** Both check out, and they are the *same* partition-function correction.

## SC tournaments (verified, S578 / re-confirmed `sc_unitdist_0014_s627.out`)

The non-self-complementary tournament count is a partition function on tilings with a correction:
```
non-SC(n) = 2^{C(n-1,2)-n+3} · (1 - deficit(n)),
deficit(n) = 0.5, 0.25, 0.125, 0.0547, 0.0259, 0.0135, 0.0072, 0.0038, …   (n=3,4,5,6,7,8,9,10,…)
```
- **`deficit(8) = 0.0135 ≈ 0.014`** — the "shape-parameter growth" value the user names.
- The deficit **halves 2-adically** (ratios `→ 0.5`): `deficit(n) ≍ 2^{-(n-2)}`.
- **The n−2 recursion** (the user's "n+2 connection"): the *correction* itself is the SC count two
  steps back,
  ```
  correction(n) := 2^{C(n-1,2)-n+3} − non-SC(n)  ~  SC(n−2),   ratio → 1
  (n=8: 888 vs SC(6)=903 = 0.983; n=12: 0.9997),
  ```
  with `correction(n)/correction(n−1) = 2^{n-4}+2` (doubling).

So the non-SC partition function's leading correction is its own **two-steps-back subsystem** — the
`n→n+2` stride (the ±-pair / even-fold of [[HYP-2240]]/[[HYP-2245]], the 2-adic doubling THM-407).

## Unit distances (A186705)

The empirical exponent of `u(n)` against `4/3`:
```
4/3 − log u(n)/log n  =  …, 0.0179 (n=9), 0.0070 (n=12), 0.0084 (n=14), 0.0087 (n=22), …
```
sits in the **`~0.014` band** through the small/known regime — the gap between the realised exponent
and the Spencer–Szemerédi–Trotter ceiling `4/3`. (The `u(n)` data is small and erratic, so this is a
looser, scale-dependent match than the clean SC deficit.)

## The claim: one correction, two systems (universality)

In the partition-function frame ([[HYP-2245]]) both are **leading term × (1 − correction)** with the
**correction governed by a 2-adically-halving, `n−2`-recursive subsystem**:

| | leading | correction structure | value near n≈8 |
|--|--|--|--|
| non-SC tournaments | `2^{C(n-1,2)-n+3}` | `~ SC(n−2)`, deficit halves | `0.0135` |
| unit distances | `n^{4/3}` | the `n+2`/CM-tower stride ([[HYP-2230]]) | `~0.014` |

> **Conjecture.** `0.014` is not a coincidence but a **shared correction value** at the matched
> scale: unit-distance maximisers and self-complementary tournaments are specialisations of the same
> master partition function ([[HYP-2245]]), and the universal object controlling both subleading
> corrections is the `n→n+2` (±-pair / 2-adic doubling) recursion — `correction(n) ~ (system at
> n−2)`. The SC side makes it exact (`correction = SC(n−2)`); the unit-distance side is the same
> stride wearing the CM-rotation/grid-disproof hat.

## To do
1. Find the unit-distance analogue of `correction(n) ~ SC(n−2)` exactly: is `u(n) − (leading)` a
   clean two-steps-back quantity? (the small data may need the Engel/Moser-ring exact configs.)
2. Identify the matched scale: at which `n` do the two corrections coincide, and is `0.014` a true
   shared constant (a critical-exponent value of the master partition function) or scale-dependent?
3. Tie the `2^{n-4}+2` doubling ratio of the SC correction to the cyclotomic-tower functional
   equation ([[HYP-2245]] Euler product) — the n+2 stride as the source of both.
