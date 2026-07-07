/-
  TournamentH7.LRCTranslateSpectrum — THE CONSECUTIVE-BLOCK (TRANSLATE) SPECTRUM
  (kind-pasteur-2026-07-06-S48, HYP-4637).

  Closes the *uniform-k* half of opus-S127's "covering escape": the families that
  evade the whole finite covering system are the `lcm`-lifts of the AP, `{i + L·k_i}`.
  When all `k_i` are equal these are **translates** — the consecutive block
  `{m, m+1, …, m+11}`.  Their loneliness spectrum is exact and clean:

      M({m, …, m+11}) = m / (2m + 11),   at the witness  t = 1/(2m+11).

  Proof (this file, the lower bound, which is all that is needed): at `t = 1/(2m+11)`
  every speed `m+i` (`0 ≤ i ≤ 11`) has residue `m+i` (no wraparound, since
  `m+11 < 2m+11`), lying in `[m, m+11] = [μ, s−μ]`, so its distance to every integer
  is `≥ m/(2m+11)`.  Hence `M ≥ m/(2m+11)`, and since `m/(2m+11) ≥ 2/15 ⟺ 11m ≥ 22`,

      **for `m ≥ 2`,  `M({m,…,m+11}) ≥ 2/15 > 2/25`  — LOOSE.**

  Only `m = 1` (the AP `{1,…,12}`) is tight (`M = 1/13`).  So the uniform-k escape is
  loose except the AP — the translate-rigidity leg of opus-S127.  Direct instance of
  `rational_point_margin` (`s = 2m+11`, `k = 1`, `μ = m`); kernel-pure.
-/
import Mathlib
import TournamentH7.LRCHarmonicGate

namespace LonelyRunner
namespace TranslateSpectrum

open HarmonicGate

/-- **The translate margin.**  The consecutive block `v i = m + i` (`i : Fin 12`),
`m ≥ 1`, has margin `≥ m/(2m+11)` at `t = 1/(2m+11)`. -/
theorem translate_margin (m : ℤ) (hm : 1 ≤ m) (v : Fin 12 → ℤ)
    (hv : ∀ i : Fin 12, v i = m + (i : ℤ)) :
    ∀ i : Fin 12, ∀ k : ℤ,
      (m : ℝ) / (2 * m + 11) ≤ |(v i : ℝ) * ((1 : ℝ) / (2 * m + 11)) - k| := by
  have hs : 0 < 2 * m + 11 := by linarith
  have hcond : ∀ i : Fin 12, m ≤ (v i * 1) % (2 * m + 11) ∧
      (v i * 1) % (2 * m + 11) ≤ (2 * m + 11) - m := by
    intro i
    have hi : (i : ℤ) ≤ 11 := by
      have : (i : ℕ) < 12 := i.isLt
      have : ((i : ℕ) : ℤ) ≤ 11 := by exact_mod_cast Nat.lt_succ_iff.mp this
      simpa using this
    have hi0 : 0 ≤ (i : ℤ) := by positivity
    rw [hv i, mul_one]
    have hnn : 0 ≤ m + (i : ℤ) := by linarith
    have hlt : m + (i : ℤ) < 2 * m + 11 := by linarith
    rw [Int.emod_eq_of_lt hnn hlt]
    exact ⟨by linarith, by linarith⟩
  have h := rational_point_margin v 1 (2 * m + 11) m hs hcond
  intro i k
  have := h i k
  simpa using this

/-- **The translate is loose for `m ≥ 2`** (`m/(2m+11) ≥ 2/15 > 2/25`).  Only the AP
(`m = 1`) is tight; every other consecutive block clears — the uniform-k escape leg. -/
theorem translate_loose (m : ℤ) (hm : 2 ≤ m) (v : Fin 12 → ℤ)
    (hv : ∀ i : Fin 12, v i = m + (i : ℤ)) :
    ∃ t : ℝ, ∀ i : Fin 12, ∀ k : ℤ, (2 : ℝ) / 25 ≤ |(v i : ℝ) * t - k| := by
  refine ⟨(1 : ℝ) / (2 * m + 11), fun i k => ?_⟩
  have hmargin := translate_margin m (by linarith) v hv i k
  have hsR : (0 : ℝ) < 2 * (m : ℝ) + 11 := by
    have : (2 : ℝ) ≤ (m : ℝ) := by exact_mod_cast hm
    linarith
  have hstep : (2 : ℝ) / 25 ≤ (m : ℝ) / (2 * m + 11) := by
    rw [div_le_div_iff₀ (by norm_num) hsR]
    have : (2 : ℝ) ≤ (m : ℝ) := by exact_mod_cast hm
    nlinarith
  exact le_trans hstep hmargin

#print axioms translate_margin
#print axioms translate_loose

end TranslateSpectrum
end LonelyRunner
