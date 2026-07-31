import Mathlib

namespace LRC14

/-- The algebraic heart of the exact non-strict integer-apex gate.

If `p < 7`, `h, γ > 0`, and an integer label lies strictly beyond
`T = 7pγ / ((7-p)h)`, then the component majorant attached to that label is
strictly below the average `h/p` required of a `p`-cover.
-/
theorem componentMajorant_lt_coverAverage
    (p w : ℕ) (h γ : ℝ)
    (hp : 0 < p) (hp7 : p < 7) (hw : 0 < w)
    (hh : 0 < h) (hγ : 0 < γ)
    (hfar :
      7 * (p : ℝ) * γ / (((7 : ℕ) - p : ℕ) * h) < (w : ℝ)) :
    h / 7 + γ / w < h / p := by
  have hpR : (0 : ℝ) < p := by exact_mod_cast hp
  have hwR : (0 : ℝ) < w := by exact_mod_cast hw
  have hgapN : 0 < 7 - p := by omega
  have hgapR : (0 : ℝ) < (7 - p : ℕ) := by exact_mod_cast hgapN
  have hcross :
      7 * (p : ℝ) * γ <
        (w : ℝ) * ((7 - p : ℕ) * h) := by
    exact (div_lt_iff₀ (mul_pos hgapR hh)).mp hfar
  have htail :
      γ / w < ((7 - p : ℕ) * h) / (7 * p) := by
    apply (div_lt_div_iff₀ hwR (mul_pos (by norm_num) hpR)).2
    convert hcross using 1 <;> ring
  calc
    h / 7 + γ / w <
        h / 7 + ((7 - p : ℕ) * h) / (7 * p) :=
      add_lt_add_left htail _
    _ = h / p := by
      field_simp [hpR.ne']
      push_cast
      ring

/-- A finite family of `p` labels all beyond the exact apex cannot cover
mass `h` when every restricted coverage obeys the component majorant.
-/
theorem sum_componentMajorant_lt_mass
    (p : ℕ) (h γ : ℝ) (S : Finset ℕ) (c : ℕ → ℝ)
    (hp : 0 < p) (hp7 : p < 7) (hh : 0 < h) (hγ : 0 < γ)
    (hcard : S.card = p)
    (hw : ∀ w ∈ S, 0 < w)
    (hfar :
      ∀ w ∈ S,
        7 * (p : ℝ) * γ / (((7 : ℕ) - p : ℕ) * h) < (w : ℝ))
    (hc : ∀ w ∈ S, c w ≤ h / 7 + γ / w) :
    ∑ w ∈ S, c w < h := by
  have havg :
      ∀ w ∈ S, c w < h / p := by
    intro w hwS
    exact lt_of_le_of_lt (hc w hwS)
        (componentMajorant_lt_coverAverage p w h γ hp hp7
          (hw w hwS) hh hγ (hfar w hwS))
  have hSne : S.Nonempty := by
    rw [← Finset.card_pos, hcard]
    exact hp
  calc
    ∑ w ∈ S, c w < ∑ _w ∈ S, h / p := by
      exact Finset.sum_lt_sum_of_nonempty hSne havg
    _ = h := by
      simp [hcard, hp.ne']

end LRC14
