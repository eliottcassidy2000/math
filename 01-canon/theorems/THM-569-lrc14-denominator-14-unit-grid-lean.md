---
id: THM-569
title: LRC(14) denominator-14 unit-grid criterion — in Lean, every unit apex point `a/14` is lonely exactly when no speed is divisible by 14
status: PROVED / FORMALIZED in Lean (`TournamentH7.LRCUnitGrid14`)
source: codex-2026-06-22 (formalization session)
depends_on:
  - THM-523   # q-witness reduction / covering-set necessary condition
  - THM-568   # apex-denominator lemma
related:
  - HYP-2910  # tightness-star exact atlas companion
  - OPEN-Q-108
---

# THM-569 — LRC(14) denominator-14 unit-grid criterion

## Statement

Let `v : ι → ℤ` be any integer speed family and let `a` be a unit modulo `14`.
Then

```lean
Lonely 14 v ((a : ℝ) / 14) ↔ ∀ i, ¬ ((14 : ℤ) ∣ v i).
```

The six concrete unit residues `1, 3, 5, 9, 11, 13` are exported as named Lean
theorems, and the general q-sieve corollary is specialized:

```lean
(∀ t : ℝ, ¬ Lonely 14 v t) → ∃ i, (14 : ℤ) ∣ v i.
```

## Proof

The forward direction is the apex-floor observation in the concrete `Lonely`
predicate: if `14 ∣ v_i`, then at time `a/14`,

```text
v_i · a/14 = (v_i/14) · a ∈ ℤ,
```

so the runner is exactly on the observer and cannot satisfy
`1/14 ≤ ‖v_i a/14‖`.

The reverse direction is `LonelyRunner.sieve_frac` with `n = q = 14`: if `a` is
coprime to `14` and no speed is divisible by `14`, the numerator
`v_i a - m·14` is a nonzero integer for every runner and every integer `m`,
hence the distance is at least `1/14`.

## Role in LRC(14)

This theorem formalizes the exact apex-grid split used in the current proof
template:

- `14`-free rows are settled at the denominator-14 unit grid, exactly.
- `14`-covering rows block every denominator-14 unit point and must be handled
  off-grid by the remaining analytic covering-strictness argument.

Assumption challenge: the useful carrier here is not a tournament on runners.
It is the quotient of the six unit residues modulo `14` by the divisibility
predicate `14 ∣ v_i`. This quotient preserves the exact `Lonely 14` truth value
on apex-grid times and destroys the off-grid Weyl/equidistribution information,
which is why it cleanly belongs to the finite side of the proof template only.
