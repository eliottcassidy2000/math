/-
  TournamentH7.LRCMissMonotone -- the structural backbone of the far-element L_y
  drop (kind-pasteur S31d): the factorial-moment miss measures `J(A,E)` are
  ANTITONE in the speed list.

  In the THM-534 L_y route, `L_y(E) = Σ_r y_r S_r(E)` with factorial moments
  `S_r(E) = Σ_{|A|=r} J(A,E)` and `J(A,E) = meas{x : every sector in A is missed by
  EVERY speed of E}`.  Adding a runner can only SHRINK each miss set (one more arc
  to avoid), so

      E ⊆ E'  ⟹  J(A,E') ≤ J(A,E).

  This is the structural half of "consec maximizes L_y / far elements drop L_y"
  (mac-mini S28/S33): far elements monotonically shrink the miss measures; the
  quantitative DROP factor `(1 - |A|/7)` per decorrelated far element is the
  Fejér/spectral-concentration residual (HYP-2873, the additive energy = `∫|Ê|⁴`).
  Sorry-free.
-/

import TournamentH7.LRCDenseCovers

namespace LonelyRunner
namespace MissMonotone

open MeasureTheory

/-- The **miss set** `{x : every inner sector in `A` is missed by every speed in
`E`}` — the support of the factorial-moment term `J(A,E)`. -/
def missSet (A : Finset ℕ) (E : List ℤ) : Set ℝ :=
  {x | ∀ e ∈ E, ∀ j ∈ A, ¬ ((j : ℝ) / 7 ≤ Int.fract ((e : ℝ) * x) ∧
      Int.fract ((e : ℝ) * x) < ((j : ℝ) + 1) / 7)}

/-- Adding speeds shrinks the miss set: `E ⊆ E' ⟹ missSet A E' ⊆ missSet A E`. -/
theorem missSet_antitone {A : Finset ℕ} {E E' : List ℤ} (h : ∀ e ∈ E, e ∈ E') :
    missSet A E' ⊆ missSet A E := by
  intro x hx e heE j hjA
  exact hx e (h e heE) j hjA

/-- The factorial-moment miss measure `J(A,E) = meas(missSet A E)`. -/
noncomputable def J (A : Finset ℕ) (E : List ℤ) : ℝ :=
  (DenseCovers.slowμ (missSet A E)).toReal

theorem J_nonneg (A : Finset ℕ) (E : List ℤ) : 0 ≤ J A E := ENNReal.toReal_nonneg

/-- **`J` is antitone in the speed list** — the structural backbone of the
far-element `L_y` drop: more speeds ⟹ smaller miss measure. -/
theorem J_antitone {A : Finset ℕ} {E E' : List ℤ} (h : ∀ e ∈ E, e ∈ E') :
    J A E' ≤ J A E :=
  ENNReal.toReal_mono (measure_ne_top _ _) (measure_mono (missSet_antitone h))

/-- Cons form: adding one runner `w` to the front never increases `J`. -/
theorem J_cons_le (A : Finset ℕ) (w : ℤ) (E : List ℤ) : J A (w :: E) ≤ J A E :=
  J_antitone (fun e he => List.mem_cons_of_mem w he)

/-! ## Axiom audit -/

#print axioms missSet_antitone
#print axioms J_antitone
#print axioms J_cons_le

end MissMonotone
end LonelyRunner
