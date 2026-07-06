/-
  TournamentH7.LRCGapCandidate — GAP CANDIDATES ARE DIVISIBILITY-RICH
  (kind-pasteur-2026-07-06-S24, HYP-4417 follow-on).

  A "gap candidate" is a 12-family that is NOT loose at `2/25` — i.e. for every
  `t` some runner is within `2/25` of an integer (the danger arcs cover).  By
  opus's divisor protection (`LRCDivisorProtection.int_far_of_not_dvd_k`), a
  family with NO multiple of `k` is lonely at `t = 1/k` with margin `1/k`; and
  `1/k > 2/25` for every `k ≤ 12`.  Contrapositive, collected over all `k`:

  > a gap candidate contains a multiple of EVERY `k ∈ {2,…,12}`.

  This is the covering-system structure of gap candidates (the parity-covering
  lens, S708): the runners must "cover" divisibility by each `k ≤ 12`.  The AP
  `{1,…,12}` is the MINIMAL such family (runner `k` is the multiple of `k`).  In
  particular a gap candidate has a runner divisible by `11`, one by `7`, one by
  `8 = 2³`, one by `9 = 3²`, one by `5` — the prime-power obstructions force
  height and additive spread.  It welds opus's `coverer_height` (HYP-4406, a
  single-`k` CRT dichotomy) into the full divisibility profile, and it sharpens
  the residue-collision (additive) side of the S23 split: the collision families
  are exactly the divisibility-covering ones.

  Kernel-pure; no `sorry`, no `native_decide`.
-/
import Mathlib
import TournamentH7.LRCDivisorProtection

namespace LonelyRunner
namespace GapCandidate

open DivisorProtection
open scoped Classical

/-- **GAP CANDIDATES ARE DIVISIBILITY-RICH**: a 12-family whose danger arcs cover
(NOT loose at `2/25`) contains a multiple of every `k ∈ {2,…,12}`.  For if no
runner were divisible by `k`, the family would be lonely at `t = 1/k` with margin
`1/k ≥ 1/12 > 2/25`, contradicting coverage at `t = 1/k`. -/
theorem gap_candidate_has_multiple (v : Fin 12 → ℤ)
    (hcov : ∀ t : ℝ, ∃ i : Fin 12, ∃ m : ℤ, |(v i : ℝ) * t - m| < 2 / 25)
    (k : ℕ) (hk2 : 2 ≤ k) (hk12 : k ≤ 12) :
    ∃ i : Fin 12, (k : ℤ) ∣ v i := by
  by_contra hcon
  push_neg at hcon
  have hk0 : 0 < k := by omega
  have hkR : (0 : ℝ) < (k : ℝ) := by exact_mod_cast hk0
  -- coverage at t = 1/k gives a runner within 2/25 of an integer there
  obtain ⟨i, m, hm⟩ := hcov (1 / (k : ℝ))
  -- but that runner is NOT divisible by k, so it is ≥ 1/k away — contradiction
  have hfar := int_far_of_not_dvd_k k hk0 (v i) (hcon i) m
  have hval : (v i : ℝ) * (1 / (k : ℝ)) - m = (v i : ℝ) / (k : ℝ) - m := by ring
  rw [hval] at hm
  -- 2/25 ≤ 1/12 ≤ 1/k
  have h12k : (1 : ℝ) / 12 ≤ 1 / (k : ℝ) := by
    apply one_div_le_one_div_of_le hkR
    exact_mod_cast hk12
  have hlt : (2 : ℝ) / 25 < 1 / (k : ℝ) := by
    have : (2 : ℝ) / 25 < 1 / 12 := by norm_num
    linarith
  linarith

/-- **The covering profile in one statement**: a gap candidate has a runner
divisible by `k` for each `k` in the danger range `{2,…,12}` simultaneously —
the "danger-covering system" the runners must realize. -/
theorem gap_candidate_covers_all (v : Fin 12 → ℤ)
    (hcov : ∀ t : ℝ, ∃ i : Fin 12, ∃ m : ℤ, |(v i : ℝ) * t - m| < 2 / 25) :
    ∀ k : ℕ, 2 ≤ k → k ≤ 12 → ∃ i : Fin 12, (k : ℤ) ∣ v i :=
  fun k hk2 hk12 => gap_candidate_has_multiple v hcov k hk2 hk12

/-- **The prime-power obstructions** — the maximal divisibility demands a gap
candidate must meet (each implies its divisors' demands): multiples of
`5, 7, 8, 9, 11` and `12`.  These are the "hard" coverers (`8 = 2³`, `9 = 3²`,
and the large primes `5, 7, 11`), forcing five spread runners; the AP meets each
minimally at `5, 7, 8, 9, 11, 12`. -/
theorem gap_candidate_prime_powers (v : Fin 12 → ℤ)
    (hcov : ∀ t : ℝ, ∃ i : Fin 12, ∃ m : ℤ, |(v i : ℝ) * t - m| < 2 / 25) :
    (∃ i, (5 : ℤ) ∣ v i) ∧ (∃ i, (7 : ℤ) ∣ v i) ∧ (∃ i, (8 : ℤ) ∣ v i) ∧
    (∃ i, (9 : ℤ) ∣ v i) ∧ (∃ i, (11 : ℤ) ∣ v i) ∧ (∃ i, (12 : ℤ) ∣ v i) :=
  ⟨gap_candidate_has_multiple v hcov 5 (by norm_num) (by norm_num),
   gap_candidate_has_multiple v hcov 7 (by norm_num) (by norm_num),
   gap_candidate_has_multiple v hcov 8 (by norm_num) (by norm_num),
   gap_candidate_has_multiple v hcov 9 (by norm_num) (by norm_num),
   gap_candidate_has_multiple v hcov 11 (by norm_num) (by norm_num),
   gap_candidate_has_multiple v hcov 12 (by norm_num) (by norm_num)⟩

#print axioms gap_candidate_has_multiple
#print axioms gap_candidate_covers_all
#print axioms gap_candidate_prime_powers

end GapCandidate
end LonelyRunner
