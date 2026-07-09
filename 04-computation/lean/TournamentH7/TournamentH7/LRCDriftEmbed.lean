/-
  TournamentH7.LRCDriftEmbed — the DRIFT-ABSORBED ruler embedding (klein-2026-07-09-S205).

  This discharges `hembed` (THM-527 Part A's realization) on the regime where the drift is absorbable,
  WITH the `e_i = Vmax − v_i` binding that kps-S105 correctly flagged as missing from the abstract
  `hembed` of klein-S203/S204 (as an isolated implication that was false: `E` unrelated to `v`).

  THE MECHANISM.  Write the witness time on the `Vmax`-ruler with a FREE fast phase:

      τ = (j + φ) / Vmax,     so    Vmax·τ = j + φ    and    frac(Vmax·τ) = φ   (for φ ∈ [0,1)).

  With `v_i = Vmax − e_i` the runner phase is exactly (mod 1, since `j : ℤ`)

      v_i·τ = (φ − t_i − d_i) + j,     t_i := e_i·j/Vmax  (the tooth),  d_i := e_i·φ/Vmax  (the DRIFT).

  So the fast phase `φ` must clear each tooth `t_i` *after* it has drifted by `d_i`, and
  `|d_i| = |e_i|·φ/Vmax ≤ spread/Vmax` since `φ < 1` (kps-S105 / opus-S176's "tooth wobble").  Placing `φ`
  at the MIDPOINT of a tooth-free gap of width `g` gives a margin `g/2`; the drift eats at most
  `spread/Vmax`.  Hence:

  > **`g > 1/7 + 2·spread/Vmax`  ⟹  every runner is `> 1/14` from the origin at `τ = (j + a + g/2)/Vmax`.**

  No Lipschitz estimate is needed: the `∀ n : ℤ, r ≤ |y − n|` characterisation of `nearInt` (kps-S31
  `GapReach`) turns the drift absorption into a plain triangle inequality.

  This is the `Vmax ≫ spread` half of Part A (drift `→ 0`); the complementary `Vmax ≤ V*` corner is the
  bounded finite check (`V* ≤ 234 / 1106 / 3^12`, kps-S105).  Self-contained on `LRCCriterionC` + `LRCGapReach`.
-/
import Mathlib
import TournamentH7.LRCCriterionC
import TournamentH7.LRCGapReach

namespace LonelyRunner
namespace LRC14Concrete

open Set

/-- Any single time bounds `Mreach` from below (`Mreach` is the global sup of `minReach`). -/
theorem Mreach_ge_of_minReach (v : Fin 13 → ℤ) (τ : ℝ) (h : (1 : ℝ) / 14 ≤ minReach v τ) :
    (1 : ℝ) / 14 ≤ Mreach v := by
  rw [Mreach_eq_global_sSup]
  refine le_trans h (le_csSup ⟨1 / 2, ?_⟩ ⟨τ, rfl⟩)
  rintro y ⟨s, rfl⟩
  exact minReach_le_half v s

/-- **The drift-absorbed ruler embedding (hembed, large-`Vmax` regime).**

Hypotheses: the co-offset BINDING `v_i = Vmax − e_i`; a spread bound `|e_i| ≤ spread`; a tooth-free
open gap `(a, a+g) ⊆ [0,1]` (no integer translate of any tooth `e_i·j/Vmax` enters it); and the DRIFT
MARGIN `1/7 + 2·spread/Vmax < g`.

Conclusion: at the ruler time `τ = (j + a + g/2)/Vmax` — the gap midpoint used as the fast phase —
every runner is at least `1/14` from the origin.  This is THM-527 Part A's realization, with the
finite-`Vmax` coupling (the tooth wobble) explicitly absorbed. -/
theorem minReach_ge_of_driftGap
    (v e : Fin 13 → ℤ) (Vmax : ℤ) (j : ℤ) (a g spread : ℝ)
    (hV : 0 < Vmax)
    (hbind : ∀ i, (v i : ℝ) = (Vmax : ℝ) - (e i : ℝ))
    (hsp : ∀ i, |(e i : ℝ)| ≤ spread)
    (ha : 0 ≤ a) (hag : a + g ≤ 1) (hgpos : 0 < g)
    (hfree : ∀ i : Fin 13, ∀ n : ℤ,
      ((e i : ℝ) * (j : ℝ) / (Vmax : ℝ) + (n : ℝ)) ∉ Ioo a (a + g))
    (hgap : 1 / 7 + 2 * (spread * (a + g / 2) / (Vmax : ℝ)) < g) :
    (1 : ℝ) / 14 ≤ minReach v (((j : ℝ) + (a + g / 2)) / (Vmax : ℝ)) := by
  have hVR : (0 : ℝ) < (Vmax : ℝ) := by exact_mod_cast hV
  have hsp0 : 0 ≤ spread := le_trans (abs_nonneg _) (hsp 0)
  set φ : ℝ := a + g / 2 with hφ
  have hφ0 : 0 ≤ φ := by rw [hφ]; linarith
  have hφ1 : φ ≤ 1 := by rw [hφ]; linarith
  have hdrnn : (0 : ℝ) ≤ spread * φ / (Vmax : ℝ) := by positivity
  -- margin: the midpoint of the free gap clears every tooth-translate by g/2
  have hmargin : ∀ i : Fin 13, ∀ n : ℤ,
      g / 2 ≤ |φ - ((e i : ℝ) * (j : ℝ) / (Vmax : ℝ) + (n : ℝ))| := by
    intro i n
    by_contra hcon
    push_neg at hcon
    rw [abs_lt] at hcon
    obtain ⟨h1, h2⟩ := hcon
    exact hfree i n ⟨by rw [hφ] at h2; linarith, by rw [hφ] at h1; linarith⟩
  -- **the SHARPENED drift bound**: `|d_i| = |e_i|·φ/Vmax ≤ spread·φ/Vmax` — proportional to the fast
  -- phase `φ`, not merely `≤ spread/Vmax`.  Low fast phases drift less.
  have hdrift : ∀ i : Fin 13, |(e i : ℝ) * φ / (Vmax : ℝ)| ≤ spread * φ / (Vmax : ℝ) := by
    intro i
    rw [abs_div, abs_mul, abs_of_pos hVR, abs_of_nonneg hφ0]
    gcongr
    exact hsp i
  -- every runner clears the origin by > 1/14 at `τ = (j + φ)/Vmax`
  have hkey : ∀ i : Fin 13,
      (1 : ℝ) / 14 < nearInt ((v i : ℝ) * (((j : ℝ) + φ) / (Vmax : ℝ))) := by
    intro i
    have halg : (v i : ℝ) * (((j : ℝ) + φ) / (Vmax : ℝ))
        = (φ - (e i : ℝ) * (j : ℝ) / (Vmax : ℝ) - (e i : ℝ) * φ / (Vmax : ℝ)) + (j : ℝ) := by
      rw [hbind i]; field_simp; ring
    rw [halg, nearInt_add_int]
    apply GapReach.nearInt_gt_of_forall_int
    intro n
    have h1 := hmargin i n
    have h2 := hdrift i
    have h3 : |φ - ((e i : ℝ) * (j : ℝ) / (Vmax : ℝ) + (n : ℝ))|
          - |(e i : ℝ) * φ / (Vmax : ℝ)|
        ≤ |φ - (e i : ℝ) * (j : ℝ) / (Vmax : ℝ) - (e i : ℝ) * φ / (Vmax : ℝ) - (n : ℝ)| := by
      have hh := abs_sub_abs_le_abs_sub
        (φ - ((e i : ℝ) * (j : ℝ) / (Vmax : ℝ) + (n : ℝ))) ((e i : ℝ) * φ / (Vmax : ℝ))
      rwa [show (φ - ((e i : ℝ) * (j : ℝ) / (Vmax : ℝ) + (n : ℝ)))
              - (e i : ℝ) * φ / (Vmax : ℝ)
            = φ - (e i : ℝ) * (j : ℝ) / (Vmax : ℝ) - (e i : ℝ) * φ / (Vmax : ℝ) - (n : ℝ)
          by ring] at hh
    linarith
  unfold minReach
  exact le_ciInf fun i => le_of_lt (hkey i)

/-- **`Mreach ≥ 1/14` from a drift-absorbed good period** — THM-527 Part A's realization on the
`Vmax ≫ spread` regime, with the `e = Vmax − v` binding. -/
theorem Mreach_ge_of_driftGap
    (v e : Fin 13 → ℤ) (Vmax : ℤ) (j : ℤ) (a g spread : ℝ)
    (hV : 0 < Vmax)
    (hbind : ∀ i, (v i : ℝ) = (Vmax : ℝ) - (e i : ℝ))
    (hsp : ∀ i, |(e i : ℝ)| ≤ spread)
    (ha : 0 ≤ a) (hag : a + g ≤ 1) (hgpos : 0 < g)
    (hfree : ∀ i : Fin 13, ∀ n : ℤ,
      ((e i : ℝ) * (j : ℝ) / (Vmax : ℝ) + (n : ℝ)) ∉ Ioo a (a + g))
    (hgap : 1 / 7 + 2 * (spread * (a + g / 2) / (Vmax : ℝ)) < g) :
    (1 : ℝ) / 14 ≤ Mreach v :=
  Mreach_ge_of_minReach v _
    (minReach_ge_of_driftGap v e Vmax j a g spread hV hbind hsp ha hag hgpos hfree hgap)

/-! ## Axiom audit -/
#print axioms minReach_ge_of_driftGap
#print axioms Mreach_ge_of_driftGap

end LRC14Concrete
end LonelyRunner
