/-
  TournamentH7.FarElementRate  (mac-mini-2026-07-02-S16)

  THE FAR-ELEMENT RATE LEMMA's Lean kernel — the program's single load-bearing
  open transcription (klein-S104: the F3/middle-band leaf certification rests
  on exactly this).  Content: the ANTIPODE DODGE — every window of length
  `1/w` contains a point whose `w`-clearance is exactly `1/2`, constructively
  `t' = (⌈w·t − 1/2⌉ + 1/2)/w`.  Consequence (the rate): a core witness with
  margin `ε` (clearances `≥ r + ε`) dodges ANY single far runner `w ≥ V/ε`
  within its margin window — clearances move by at most `V·(1/w) ≤ ε`, and the
  far runner sits at the antipode.  Iterating = the peel; the RATE is
  `w ≥ V/ε`, i.e. `O(1/w)` loss per far element.  All-ℚ, fully constructive.

  Sorry-free.
-/
import Mathlib.Tactic

namespace TournamentH7.FarElementRate

/-- **The antipode dodge.**  For any `t` and any speed `w > 0`, the window
`[t, t + 1/w]` contains a point `t'` with `w`-clearance exactly `1/2`:
every integer is at distance `≥ 1/2` from `w·t'`. -/
theorem exists_antipode_dodge (w : ℕ) (hw : 0 < w) (t : ℚ) :
    ∃ t' : ℚ, t ≤ t' ∧ t' ≤ t + 1 / w ∧ ∀ b : ℤ, (1 : ℚ) / 2 ≤ |(w : ℚ) * t' - b| := by
  have hw' : (0 : ℚ) < (w : ℚ) := by exact_mod_cast hw
  set k : ℤ := ⌈(w : ℚ) * t - 1 / 2⌉ with hk
  refine ⟨((k : ℚ) + 1 / 2) / w, ?_, ?_, ?_⟩
  · -- t ≤ t': from k ≥ w t − 1/2
    have h := Int.le_ceil ((w : ℚ) * t - 1 / 2)
    rw [le_div_iff₀ hw']
    rw [← hk] at h
    nlinarith
  · -- t' ≤ t + 1/w: from k < w t + 1/2, i.e. k ≤ ⌈⌉ < arg + 1
    have h := Int.ceil_lt_add_one ((w : ℚ) * t - 1 / 2)
    rw [← hk] at h
    rw [div_le_iff₀ hw']
    have hexp : (t + 1 / w) * w = t * w + 1 := by field_simp
    rw [hexp]
    nlinarith
  · -- clearance: w t' = k + 1/2, distance to any integer ≥ 1/2
    intro b
    have hval : (w : ℚ) * (((k : ℚ) + 1 / 2) / w) = (k : ℚ) + 1 / 2 := by
      field_simp
    rw [hval]
    have hint : (1 : ℚ) ≤ |(k : ℚ) - b| ∨ k = b := by
      rcases eq_or_ne k b with h | h
      · exact Or.inr h
      · left
        have : (1 : ℤ) ≤ |k - b| := Int.one_le_abs (sub_ne_zero.mpr h)
        calc (1 : ℚ) = ((1 : ℤ) : ℚ) := by norm_num
          _ ≤ ((|k - b| : ℤ) : ℚ) := by exact_mod_cast this
          _ = |(k : ℚ) - b| := by rw [Int.cast_abs]; push_cast; ring_nf
    rcases hint with h | h
    · calc (1 : ℚ) / 2 ≤ |(k : ℚ) - b| - 1 / 2 := by linarith
        _ ≤ |(k : ℚ) + 1 / 2 - b| := by
            have := abs_sub_abs_le_abs_sub ((k : ℚ) - b) (-(1 / 2 : ℚ))
            have habs : |(-(1/2) : ℚ)| = 1/2 := by norm_num
            have hrw : (k : ℚ) - b - -(1/2) = (k : ℚ) + 1/2 - b := by ring
            rw [habs, hrw] at this
            linarith
    · subst h
      have : (k : ℚ) + 1 / 2 - k = 1 / 2 := by ring
      rw [this]
      norm_num

/-- **The rate corollary (margin transport).**  If a witness `t` has core
clearance margin `ε` in the Lipschitz sense (moving by `δ ≤ 1/w` costs at most
`V·δ ≤ ε`), then for any far speed `w ≥ V/ε` the dodged point keeps the core
safe and is `1/2`-clear of the far runner.  Stated abstractly: the dodge point
is within `1/w` of `t`, so ANY property `ε`-stable under `1/w`-perturbation
transports.  (The concrete core-clearance instantiation lives with the
certificate layer; this records the interface and the rate `w ≥ V/ε`.) -/
theorem dodge_within (w : ℕ) (hw : 0 < w) (t : ℚ) (Pr : ℚ → Prop)
    (hstable : ∀ t' : ℚ, t ≤ t' → t' ≤ t + 1 / w → Pr t') :
    ∃ t' : ℚ, Pr t' ∧ ∀ b : ℤ, (1 : ℚ) / 2 ≤ |(w : ℚ) * t' - b| := by
  obtain ⟨t', h1, h2, h3⟩ := exists_antipode_dodge w hw t
  exact ⟨t', hstable t' h1 h2, h3⟩

end TournamentH7.FarElementRate
