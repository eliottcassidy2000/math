---
id: THM-709
title: "The prime clean ruler — cleanRuler_of_not_dvd_13: any 13-speed family with NO speed divisible by 13 has q=13 as a clean ruler (bandCount ≡ 0), because the safe band [⌈13/14⌉,⌊13·13/14⌋]=[1,12] at q=13 is exactly the nonzero residues mod 13. This discharges the single LRC(14) Lean obligation hB5 for the entire sub-class avoiding 13∣v_i (via kps THM-707's exists_B5_pos_of_cleanRuler); the hard residuals are exactly those WITH a speed ≡0 mod 13 — the '13 is prime' phenomenon."
status: PROVED and FORMALIZED, kernel-pure. LRCPrimeRuler.lean (cleanRuler_of_not_dvd_13), root-wired, builds green, axioms [propext, Classical.choice, Quot.sound] — no sorryAx, no native_decide. Elementary residue arithmetic (Prime.dvd_mul + Int.dvd_of_emod_eq_zero + omega). Verified numerically: q=13 clean on 2000/2000 random 13-speed families avoiding 13∣·.
source: opus-2026-07-11-S227
depends_on:
  - THM-707   # kps: the clean-ruler reduction of hB5 (exists_B5_pos_of_cleanRuler)
related:
  - THM-671   # klein: B5 ≤ liveCount
  - HYP-4382  # mac-mini: the tight locus is the AP because 13 is prime (residue_pinning_13)
  - LRCPairSumDispatch  # the pair-sum ruler for the hard case (speed ≡ 0 mod 13)
external: the prime residue field 𝔽₁₃ (units = nonzero residues); the LRC danger arc 1/14.
---

# THM-709 — the prime clean ruler `q = 13`

## Statement (Lean: `cleanRuler_of_not_dvd_13`)

> For `v : Fin 13 → ℤ`, if `13 ∤ v i` for every runner `i`, then `CleanRuler v` — i.e. `q = 13` has a live
> multiplier and `bandCount v 13 p ≤ 5` for every `p ∈ (0,13)` (in fact `bandCount v 13 p = 0`).

By kps's `exists_B5_pos_of_cleanRuler` (THM-707) this yields the per-family obligation `∃ q, 0 < B5 v q`,
so **`hB5` — the single remaining analytic obligation of the LRC(14) grand assembly — is discharged for the
entire residual sub-class with no speed divisible by 13.**

## Why (one paragraph)

At `q = 13` the safe band `[⌈13/14⌉, ⌊13·13/14⌋] = [1, 12]` is *exactly* the nonzero residues mod 13.  So a
runner is *out* of band at `p` iff its residue `(v i · p) mod 13 = 0` iff `13 ∣ v i · p`.  Since `13` is
prime and `13 ∤ p` for `p ∈ (0,13)`, having `13 ∤ v i` gives `13 ∤ v i · p`, so every residue is nonzero,
every runner is in band, and `bandCount ≡ 0` across all 12 multipliers: `liveCount = 12`, `maxBand = 0`.

## Scope — and the hard complement

The families this does **not** cover are exactly those with a speed `≡ 0 mod 13` (a multiple of 13): the AP
`{1,…,13}` (speed 13), `{1,…,12, 26}` (`26 = 2·13`), etc.  These are the tight / near-tight residuals, where
`q = 13` fails (the multiple-of-13 runner sits at residue 0 = danger) and the **pair-sum ruler** (kps THM-707,
binding `{1..12,26}` at `q = 27`) is required.  This is the "13 is prime" phenomenon (mac-mini HYP-4382): the
extremal tight locus is pinned at the multiples of 13.

The same residue argument gives a clean ruler at **any prime `q ∈ {7, 11, 13}`** (each has safe band
`[1, q−1]` = all nonzero residues, and for prime `q ≤ 13` and `p ∈ (0,q)`, `bandCount = #{i : q ∣ v_i}`).  So
a family is clean-ruled unless it has a speed divisible by **each** of 7, 11, 13 — shrinking the genuinely
hard class further.

## Generalized + composite (opus-2026-07-11-S228, FORMALIZED)

`LRCPrimeRuler.lean` now proves the general lemma and the composite, all kernel-pure
(`[propext, Classical.choice, Quot.sound]`), root-wired, build green:

- **`cleanRuler_of_prime_not_dvd (q) (hq : q.Prime) (hq14 : q ≤ 14) (h : ∀ i, ¬ (q:ℤ) ∣ v i) : CleanRuler v`**
  — the general prime ruler (band `[1,q−1]` for any prime `q ≤ 14`; the inequalities close by `omega` from
  `q ≤ 14`, `Euclid` from `Nat.prime_iff_prime_int`).
- **`cleanRuler_of_not_dvd_7 / _11 / _13`** — the three instances (`by norm_num` on primality and `q ≤ 14`).
- **`cleanRuler_of_avoids_some_prime`** — if the family avoids `0 mod q` for *some* `q ∈ {7,11,13}`, it has a
  clean ruler.  ⟹ **the only residuals NOT discharged by a prime clean ruler are those with a speed divisible
  by each of 7, 11, 13** — a structured hard core (measured ≈ 40% of random 13-speed families, far smaller
  than "all near-AP") handled by the pair-sum ruler (kps THM-707 / cont.30).

This is a clean two-branch split of `hB5`: `[avoids 0 mod some prime in {7,11,13} → prime ruler, maxBand 0]`
∪ `[multiple of each of 7,11,13 → pair-sum ruler]`.

## Files
`04-computation/lean/TournamentH7/TournamentH7/LRCPrimeRuler.lean` (root-wired);
`04-computation/lrc14_clean_ruler_nearAP_opus_S226.py` (the S226 verification that motivated it).
