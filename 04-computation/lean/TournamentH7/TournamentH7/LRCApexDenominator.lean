/-
  TournamentH7.LRCApexDenominator -- THM-568, the apex-denominator lemma
  (kind-pasteur S31aa).  The arithmetic CORE of "a tight LRC(14) optimum has
  denominator a multiple of 14".

  At a tight optimum `t* = a/D` (lowest terms, `M(S) = 1/14`) the origin sits in an
  empty arc of length exactly `1/7`, so a binding runner `v` satisfies
  `‖v · (a/D)‖ = 1/14`, i.e.

        v · a / D = m ± 1/14   for some integer m,
        ⟺   14 · v · a = D · (14 m ± 1).

  Since `14 m ± 1` is coprime to `14`, this forces **`14 ∣ D`** — the optimum lives
  at the apex denominator `14 = 2·7`.  This is the LRC analogue of THM-079's
  cycle-forcing, made purely arithmetic.  Sorry-free; the only mathematical content
  is the coprimality + `IsCoprime.dvd_of_dvd_mul_right`.

  The companion `LRCApex7Floor.D14_never_certifies` is the *floor* (a multiple of 14
  is blocked at `D=14`); together they say: tight ⟹ `14∣D`, and 14-covering ⟹ the
  binding cannot be at `D=14`, so a tight set must be 14-free.
-/

import Mathlib.RingTheory.Coprime.Basic
import Mathlib.Tactic

namespace LonelyRunner
namespace ApexDenominator

/-- `14 m ± 1` is coprime to `14` (it is `≡ ±1 (mod 14)`). -/
theorem isCoprime_fourteen_step (m e : ℤ) (he : e = 1 ∨ e = -1) :
    IsCoprime (14 : ℤ) (14 * m + e) := by
  have h1 : (14 : ℤ) * m + e = e + 14 * m := by ring
  rw [h1]
  rcases he with rfl | rfl
  · exact (isCoprime_one_right).add_mul_left_right m
  · exact (isCoprime_one_right.neg_right).add_mul_left_right m

/-- **THM-568, arithmetic core.**  If a (binding) runner satisfies the rational
binding equation `14 · v · a = D · (14 m + e)` with `e = ±1` — i.e. `‖v·(a/D)‖ = 1/14`
— then **`14 ∣ D`**: the tight optimum denominator is a multiple of the apex `14`. -/
theorem apex_dvd_of_binding (v a D m e : ℤ) (he : e = 1 ∨ e = -1)
    (h : 14 * v * a = D * (14 * m + e)) : (14 : ℤ) ∣ D := by
  have hcop : IsCoprime (14 : ℤ) (14 * m + e) := isCoprime_fourteen_step m e he
  have hdvd : (14 : ℤ) ∣ D * (14 * m + e) := ⟨v * a, by rw [← h]; ring⟩
  exact hcop.dvd_of_dvd_mul_right hdvd

/-- Both binders of a tight optimum sum to a multiple of `D`: if `v_a·a/D = m_a + 1/14`
and `v_b·a/D = m_b − 1/14` then `(v_a+v_b)·a = D·(m_a+m_b)`, so `D ∣ (v_a+v_b)·a`. -/
theorem binders_sum_dvd (va vb a D ma mb : ℤ)
    (ha : 14 * va * a = D * (14 * ma + 1)) (hb : 14 * vb * a = D * (14 * mb - 1)) :
    (D : ℤ) ∣ (va + vb) * a := by
  refine ⟨ma + mb, ?_⟩
  have h14 : (14 : ℤ) * ((va + vb) * a) = 14 * (D * (ma + mb)) := by
    have : (14 : ℤ) * ((va + vb) * a) = 14 * va * a + 14 * vb * a := by ring
    rw [this, ha, hb]; ring
  exact mul_left_cancel₀ (by norm_num : (14 : ℤ) ≠ 0) h14

/-! ## Axiom audit -/

#print axioms apex_dvd_of_binding
#print axioms binders_sum_dvd

end ApexDenominator
end LonelyRunner
