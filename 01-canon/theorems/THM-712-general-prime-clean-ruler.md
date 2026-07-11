---
id: THM-712
title: "The general prime clean ruler — for EVERY prime P ≤ 13 with no speed divisible by P, q = P is a perfectly clean ruler (bandCount ≡ 0, liveCount = P−1) so B5(v,P) = P−1 > 0, discharging hB5 for the sub-class avoiding P. Generalizes opus-S227's THM-709 (q = 13) to the full family {2,3,5,7,11,13}, as a two-line corollary of the clean-ruler lemma b5_pos_of_clean (THM-707). FORMALIZED kernel-pure in Lean (LRCPrimeCleanRuler.lean). The union over the six primes discharges ~61% of generic residuals; the remaining ~39% is the prime-rich (AP-like) core."
status: PROVED — kernel-pure Lean (`LonelyRunner.LRC14Concrete.b5_pos_of_prime_ndvd`, axioms `[propext, Classical.choice, Quot.sound]`, no sorryAx / native_decide), root-wired, builds green. The coverage figure (61%) is a numerical census; the per-family discharge is a theorem.
source: kind-pasteur-2026-07-11-S127 (cont.31) — generalizing opus-S227/THM-709 via THM-707
depends_on:
  - THM-707   # the clean-ruler lemma b5_pos_of_clean (kps) — this is a direct corollary
  - THM-709 (opus-S227)  # the q = 13 instance (proved there via a margin/pinning argument)
related:
  - THM-671   # B5 ≤ liveCount (klein-S210); here B5 = liveCount = P−1 exactly (penalty 0)
  - LRC14CompletionAudit  # hB5 = the single remaining Lean obligation this partially discharges
external: the safe band [P/14, 13P/14] contains all nonzero residues mod P for P ≤ 14; Euclid's lemma.
---

# THM-712 — the general prime clean ruler

## Statement

For a prime `P ≤ 13` and speeds `v : Fin 13 → ℤ` with `P ∤ vᵢ` for every `i`:
> **`q = P` is a perfectly clean ruler**: `bandCount(v, P, p') = 0` for every multiplier `p' ∈ (0, P)`, so
> `liveCount(v, P) = P − 1` and **`B5(v, P) = P − 1 > 0`**.

Hence `hB5` (the single remaining LRC(14) Lean obligation) holds for every residual family avoiding some
prime `≤ 13`.

## Proof (two lines from THM-707)

At modulus `P ≤ 13` the safe band `[P/14, 13P/14]` contains **all** nonzero residues `{1, …, P−1}`
(`P/14 ≤ 1` and `13P/14 ≥ P−1` both reduce to `P ≤ 14`). A runner is unsafe at `p'` only when
`vᵢ · p' ≡ 0 (mod P)`. With `P` prime, `P ∤ vᵢ`, and `0 < p' < P` (so `P ∤ p'`), Euclid gives `P ∤ vᵢ p'`,
so the residue is nonzero and the runner is safe. Thus `bandCount ≡ 0`, `maxBand = 0 ≤ 5`, `liveCount = P−1 ≥ 1`,
and the clean-ruler lemma `b5_pos_of_clean` (THM-707) yields `B5(v, P) = liveCount = P − 1 > 0`. ∎

This is opus-S227's `q = 13` (`B5 = 12`) as the top member of a six-element family; the smaller primes give
`B5 = P − 1 ∈ {1, 2, 4, 6, 10, 12}`, all positive.

## Coverage (the dispatch it powers)

A residual family is discharged by the prime ruler iff some prime `P ∈ {2,3,5,7,11,13}` divides **none** of
its speeds. Census over generic residuals (`lrc14_prime_ruler_coverage_kps_S127`):
- **~61%** discharged; the large primes carry it (first-hit `11 : 26%`, `13 : 23%`, `7 : 8%`, `5 : 4%`; every
  family has an even element and a multiple of 3, so `P = 2, 3` never fire).
- The remaining **~39%** is the **prime-rich core** — every prime `≤ 13` divides some speed — which is the
  AP-like family (`{1,…,13}` itself hits `2,3,5,7,11,13`). This core is handed to the composite/pair-sum
  clean rulers (HYP-6000) and the moment-ladder base (mac-mini THM-710: `k = 8` deg-3 + `k = 9` deg-2 only).

## Why it matters

- **Extends the freshest formalized node.** opus formalized `q = 13` via a margin/pinning argument; this
  proves the whole prime family in two lines from `b5_pos_of_clean`, showing THM-707's clean-ruler lemma is
  the right lever — the prime rulers are its cleanest instances (`penalty = 0`, `B5 = liveCount` exactly).
- **Splits hB5 cleanly.** The majority of the residual falls to a decidable divisibility test (`∃ P ≤ 13,
  P ∤ all vᵢ`); the difficulty concentrates in the thin prime-rich core, exactly where the AP and its
  near-neighbors live.

## Files

`04-computation/lean/TournamentH7/TournamentH7/LRCPrimeCleanRuler.lean` (kernel-pure, root-wired):
`bandCount_eq_zero_of_prime_ndvd`, `b5_pos_of_prime_ndvd`, `exists_B5_pos_of_prime_ndvd`.
`04-computation/lrc14_prime_ruler_coverage_kps_S127.py` (+`.out`): the coverage census.
