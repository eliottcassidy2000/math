---
id: THM-386
name: lrc-lonely-central-box-grounding
status: PROVED (Lean-checked)
date: 2026-06-01
session: oracle-2026-06-01-S513
depends_on:
  - THM-369
lean: 04-computation/lean/TournamentH7/TournamentH7/LonelyRunner.lean
---

# THM-386: `Lonely` is the standard central-box / nearest-integer predicate

This grounds the proof-friendly `Lonely` predicate used by the repo's Lean LRC
sieve (THM-369) against the **standard Lonely Runner statement**: every
`v_i · t` has fractional part in the central box `[1/n, 1-1/n]` (Cusick's
view-obstruction box). Machine-checked, axiom-clean.

## Statements (all in `LonelyRunner.lean`)

Recall `Lonely n v t := ∀ i, ∀ m : ℤ, 1/n ≤ |v i · t − m|` ("every `v i · t` is
at distance `≥ 1/n` from every integer").

- **`far_iff_fract`** (the core real-analysis lemma). For all `c, x : ℝ`,
  ```
  (∀ m : ℤ, c ≤ |x − m|)  ↔  (c ≤ Int.fract x ∧ Int.fract x ≤ 1 − c).
  ```
  "`x` is `≥ c` from every integer" iff "its fractional part lies in the band
  `[c, 1−c]`." Proof: forward, evaluate at `⌊x⌋` (gives `c ≤ fract x`) and at
  `⌊x⌋+1` (gives `c ≤ 1 − fract x`); backward, for any `m` write
  `x − m = fract x − (m − ⌊x⌋)` and bound `|fract x − k|` by cases `k ≤ 0` /
  `k ≥ 1` using `0 ≤ fract x < 1`.

- **`lonely_iff_fract_mem`** (the view-obstruction form).
  ```
  Lonely n v t  ↔  ∀ i, 1/n ≤ Int.fract (v i · t) ∧ Int.fract (v i · t) ≤ 1 − 1/n.
  ```
  Immediate from `far_iff_fract` via `forall_congr'`. This is exactly the
  statement "the orbit `t ↦ (v i · t)` lies in the central cube `[1/n, 1−1/n]^k`."

- **`lonely_add_one`** (periodicity). For integer speeds, loneliness is
  `1`-periodic: `Lonely n v (t + 1) ↔ Lonely n v t` (reindex `m ↦ m ± v i`).

## Significance

THM-369 introduced `Lonely` in the "far from every integer" form because it makes
lower bounds on the nearest-integer norm into clean universally-quantified
inequalities (no `Int.fract`/measure machinery in the sieve proofs). THM-386
proves that this form is **literally** the textbook predicate — fractional part in
the central box — so the Lean development is now formally anchored to the standard
Lonely Runner Conjecture, not merely a convenient surrogate. Periodicity makes the
time variable a genuine `ℝ/ℤ` quantity.

All three carry a `#print axioms` audit showing only the foundational
`propext / Classical.choice / Quot.sound`.

## Related

- THM-369: the Lean denominator-sieve master lemma (uses `Lonely`).
- HYP-1981 / HYP-1987: the marked-A000568 source-reachability reformulation.
- `07-reflections/lrc-denominator-sieve-lean-and-two-regimes-s18.md`.
