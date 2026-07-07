/-
  TournamentH7.LRCBandFloor — THE GENERAL AVOID-BAND 2/25-FLOOR CERTIFICATE
  (kind-pasteur-2026-07-06-S51, HYP-4677).

  The single tool underneath every layer of the (C) covering system.  Clearing a
  family at modulus `q` means finding a rotation `c ∈ (ℤ/q)*` that puts every speed's
  residue into the *safe band* `[μ, q−μ]` (distance `≥ μ` from `0`), where `μ` is the
  band half-width `⌈2q/25⌉`.  Whenever `μ/q ≥ 2/25` (i.e. `2q ≤ 25μ`) this forces
  `M ≥ 2/25` — the family is LOOSE.

      **`∀i, μ ≤ (vᵢ·c) % q ≤ q−μ`  with  `2q ≤ 25μ`   ⟹   M(v) ≥ 2/25`  (at `t = c/q`).**

  This subsumes the whole covering:
    · `q = 25, μ = 2`  — the mod-25 non-transversal layer (`LRCMod25Floor`);
    · `q ≤ 12, μ = 1`  — the small-modulus / Fan-Sun gcd layer (`LRCSmallModFloor`);
    · `13 ≤ q ≤ 32, μ = ⌈2q/25⌉` — the near-AP moat layer (kps-S50), e.g. `q=17, μ=2`
      (`2·17 = 34 ≤ 50 = 25·2`), `q=32, μ=3` (`64 ≤ 75`).
  So every covering certificate is one instance of `loose_of_band`.  Direct instance of
  `rational_point_margin`; kernel-pure.
-/
import Mathlib
import TournamentH7.LRCHarmonicGate

namespace LonelyRunner
namespace BandFloor

open HarmonicGate

/-- **The general avoid-band loose certificate.**  If a rotation `c` puts every
`vᵢ·c` into the band `[μ, q−μ]` mod `q`, and the band clears `2/25` (`2q ≤ 25μ`), then
at `t = c/q` every runner is `≥ 2/25` from every integer — LOOSE. -/
theorem loose_of_band {ι : Type*} (v : ι → ℤ) (q μ c : ℤ) (hq : 0 < q) (hband : 2 * q ≤ 25 * μ)
    (hcov : ∀ i, μ ≤ (v i * c) % q ∧ (v i * c) % q ≤ q - μ) :
    ∃ t : ℝ, ∀ i, ∀ m : ℤ, (2 : ℝ) / 25 ≤ |(v i : ℝ) * t - m| := by
  refine ⟨(c : ℝ) / q, fun i m => ?_⟩
  have h := rational_point_margin v c q μ hq hcov i m
  have hqR : (0 : ℝ) < (q : ℝ) := by exact_mod_cast hq
  have h2 : (2 : ℝ) / 25 ≤ (μ : ℝ) / q := by
    rw [div_le_div_iff₀ (by norm_num) hqR]
    have : (2 : ℝ) * q ≤ 25 * μ := by exact_mod_cast hband
    linarith
  exact le_trans h2 (by simpa using h)

/-- **Moat instance, `q = 17, μ = 2`** (`2·17 = 34 ≤ 50`) — the dominant near-AP moat
modulus (kps-S50 histogram peak).  A rotation off the band `[2,15]` mod 17 clears. -/
theorem loose_at_17 {ι : Type*} (v : ι → ℤ) (c : ℤ)
    (hcov : ∀ i, 2 ≤ (v i * c) % 17 ∧ (v i * c) % 17 ≤ 15) :
    ∃ t : ℝ, ∀ i, ∀ m : ℤ, (2 : ℝ) / 25 ≤ |(v i : ℝ) * t - m| :=
  loose_of_band v 17 2 c (by norm_num) (by norm_num) hcov

/-- **Moat endpoint, `q = 32, μ = 3`** (`2·32 = 64 ≤ 75`) — the top of the moat range
`[14,32]` (klein-S144).  A rotation off `[3,29]` mod 32 clears. -/
theorem loose_at_32 {ι : Type*} (v : ι → ℤ) (c : ℤ)
    (hcov : ∀ i, 3 ≤ (v i * c) % 32 ∧ (v i * c) % 32 ≤ 29) :
    ∃ t : ℝ, ∀ i, ∀ m : ℤ, (2 : ℝ) / 25 ≤ |(v i : ℝ) * t - m| :=
  loose_of_band v 32 3 c (by norm_num) (by norm_num) hcov

#print axioms loose_of_band
#print axioms loose_at_17
#print axioms loose_at_32

end BandFloor
end LonelyRunner
