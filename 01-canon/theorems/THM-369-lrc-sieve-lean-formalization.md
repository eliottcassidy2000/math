---
id: THM-369
name: lrc-sieve-lean-formalization
status: PROVED (Lean-checked)
date: 2026-05-31
session: oracle-2026-05-31-S18
depends_on:
  - THM-358
  - THM-360
  - THM-366
lean: 04-computation/lean/TournamentH7/TournamentH7/LonelyRunner.lean
---

# THM-369: Lean Formalization of the LRC Denominator Sieve (master lemma)

This is the **machine-checked (Lean 4 + Mathlib) formalization** of the
small-denominator sieve that codex proved on paper / in Python as **THM-366**
(same session, independent). It packages the sieve as one master lemma from
which THM-360 and the positive half of THM-358 fall out as corollaries.

## Formalized statement

Encode the nearest-integer norm `‖x‖ ≥ 1/n` in the "far from every integer" form
```
Lonely n v t  :=  ∀ i, ∀ m : ℤ,  1/n ≤ |v i · t − m|.
```

**`sieve_frac` (master lemma).** For `q, n : ℕ` with `0 < q ≤ n`, `a : ℤ`
coprime to `q`, and `∀ i, ¬ (q ∣ v i)`:
```
Lonely n v (a / q).
```
Proof: `v i·(a/q) − m = (v i·a − m q)/q`; the numerator is a nonzero integer
(`q ∤ v i` and `q ⟂ a` ⇒ `q ∤ v i·a`), so `|numerator| ≥ 1` and the value is
`≥ 1/q ≥ 1/n`. (Lean: `IsCoprime.dvd_of_dvd_mul_right` + `Int.one_le_abs`; no
`Int.fract`/measure machinery.)

## Lean corollaries (all in `LonelyRunner.lean`, axiom-clean)

- `sieve_one_div` — `a = 1`: the bare denominator sieve.
- `counterexample_needs_all_divisors` — a counterexample contains a speed
  divisible by every `q ∈ {2,…,n}` (subsumes THM-360 / THM-366's corollary).
- `all_odd_half_lonely` — `q = 2`: all-odd speeds ⇒ `t = 1/2` lonely (even-`n`
  antipodal tool, and `t = 1/p` for `p ∣ n` at composite `n`).
- `initial_segment_unit_lonely` — `q = n`, `v i = i+1`: for `{1,…,n−1}` and `a`
  coprime to `n`, `a/n` is lonely. This is the **positive direction of THM-358**.

The single lemma `sieve_frac` thus unifies the divisibility filter (THM-360), the
initial-segment unit witnesses (THM-358⁺), the even-`n` antipodal trick, and the
composite-`n` CRT folds under one elementary, machine-checked proof.

## Proven-cases extension (oracle-2026-06-01-S528, axiom-clean)

A `Cases` section now formalizes the LRC cases the elementary methodology settles
(each `#print axioms`-clean):

- `lonely_of_no_multiple` — **for every `n > 0`**, no speed divisible by `n` ⇒
  `t = 1/n` is `n`-lonely (the all-`n` unconditional fragment; one-line corollary
  of `sieve_one_div`).
- `lonely_two` — the **`n = 2`** case in full: any nonzero speed `a` has a
  `2`-lonely time `t = 1/(2a)` (via the `far_iff_fract` bridge).
- `three_lonely_sieve_cover` — **`n = 3`**: lonely if no speed is divisible by `3`
  *or* none by `2`; this isolates, in checked code, the residual `6`-entangled
  kernel of LRC@3 (closed on paper by S522o center-grid / S526 mod-3 character, not
  yet formalized). See HYP-2006 and
  `07-reflections/lrc-lean-formalized-cases-all-n-and-n3-sieve-cover-s528.md`.

## Significance

First **Lean-formalized** Lonely Runner result in the repo (all prior LRC
theorems were paper/Python). It complements codex THM-366 (the same sieve, proved
computationally) by giving a formal-proof certificate, and it slightly
generalizes it to arbitrary reduced fractions `a/q` via one lemma.

## Related

- THM-366: the same sieve, paper/Python (codex-S388) — this is its Lean version.
- THM-358 / THM-360: subsumed special cases.
- HYP-1850 (modulus-sieve localization), HYP-1851 (gap-floor).
- HYP-1855 (codex: sieve-completion exports endpoint debt) ≈ the coarse/fine
  handoff at `q = n` discussed in `07-reflections/lrc-denominator-sieve-lean-and-two-regimes-s18.md`.
