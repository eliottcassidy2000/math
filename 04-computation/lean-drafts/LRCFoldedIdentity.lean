/- LRCFoldedIdentity.lean -- opus-2026-07-17-S346 (HYP-7290 / THM-965), v3.
   THE TWO-VARIABLE FOLDED IDENTITY over boxeph's muNum:

     14 * muNum a b = 4ab + fold((a+b) % 14) - fold((b-a) % 14),
     fold r = r * (14 - r),   for 1 ≤ a ≤ b.

   Division-free route: peel j = 0, split the shifted range at Qm = (b-a)/14
   (plateau 4a; descending tail), evaluate the tail by Nat.le_induction, and
   close with ONE linear_combination whose coefficients come from the exact
   factorizations
     S-bracket = -(S - 14Qp - r+)(S - 14 - 14Qp + r+),
     D-bracket = +(D - 14Qm - r-)(D - 14 - 14Qm + r-).
   Paper: THM-965 (verified 400/400 exact, every gcd). -/
import Mathlib
import TournamentH7.LRCTreeHunter

namespace LonelyRunner.LRC14.Hunter

/-- The descending tail, division-free (ℤ), by `Nat.le_induction`. -/
theorem tail_sum (S : ℤ) (u v : ℕ) (huv : u ≤ v) :
    (∑ j ∈ Finset.Ico u v, (2 * S - 28 * ((j : ℤ) + 1)))
      = 2 * S * ((v : ℤ) - u) - 14 * ((v : ℤ) * (v + 1) - u * (u + 1)) := by
  induction v, huv using Nat.le_induction with
  | base => simp
  | succ v huv ih =>
      rw [Finset.sum_Ico_succ_top huv, ih]
      push_cast
      ring

/-- **THM-965, muNum form.** -/
theorem muNum_folded (a b : ℕ) (ha : 1 ≤ a) (hab : a ≤ b) :
    (14 * muNum a b : ℤ)
      = 4 * (a : ℤ) * b
        + (((a + b) % 14 : ℕ) : ℤ) * (14 - (((a + b) % 14 : ℕ) : ℤ))
        - (((b - a) % 14 : ℕ) : ℤ) * (14 - (((b - a) % 14 : ℕ) : ℤ)) := by
  set S := a + b with hSdef
  set D := b - a with hDdef
  set Qp := S / 14 with hQpdef
  set Qm := D / 14 with hQmdef
  have hDS : D ≤ S := by omega
  have hQmp : Qm ≤ Qp := Nat.div_le_div_right hDS
  have h14Qm : 14 * Qm ≤ D := by
    rw [hQmdef, mul_comm]
    exact Nat.div_mul_le_self D 14
  have h14Qp : 14 * Qp ≤ S := by
    rw [hQpdef, mul_comm]
    exact Nat.div_mul_le_self S 14
  have hmin : 2 * min a b = 2 * a := by rw [Nat.min_eq_left hab]
  -- term shapes after the j = 0 peel
  have hshape : ∀ j ∈ Finset.range Qp,
      (if j + 1 = 0 then min S (2 * min a b)
       else 2 * min (S - 14 * (j + 1)) (2 * min a b))
      = if j < Qm then 4 * a else 2 * (S - 14 * (j + 1)) := by
    intro j hj
    rw [Finset.mem_range] at hj
    rw [if_neg (Nat.succ_ne_zero j), hmin]
    rcases lt_or_ge j Qm with hlt | hge
    · rw [if_pos hlt]
      have h14 : 14 * (j + 1) ≤ D := by omega
      have hge2a : 2 * a ≤ S - 14 * (j + 1) := by omega
      rw [Nat.min_eq_right hge2a]
      ring
    · rw [if_neg (not_lt.mpr hge)]
      have hDlt : D < 14 * (j + 1) := by
        have h1 : D / 14 < j + 1 := by omega
        have := (Nat.div_lt_iff_lt_mul (by norm_num : 0 < 14)).mp h1
        omega
      have hle2a : S - 14 * (j + 1) ≤ 2 * a := by omega
      rw [Nat.min_eq_left hle2a]
  have hj0 : min S (2 * min a b) = 2 * a := by
    rw [hmin]
    omega
  -- the ℕ assembly
  have hnat : muNum a b
      = 2 * a + 4 * a * Qm
        + ∑ j ∈ Finset.Ico Qm Qp, 2 * (S - 14 * (j + 1)) := by
    unfold muNum
    rw [← hSdef, ← hQpdef, Finset.sum_range_succ',
      Finset.sum_congr rfl hshape, if_pos rfl, hj0, Finset.range_eq_Ico,
      ← Finset.sum_Ico_consecutive _ (Nat.zero_le Qm) hQmp]
    have hplateau : ∀ j ∈ Finset.Ico 0 Qm,
        (if j < Qm then 4 * a else 2 * (S - 14 * (j + 1))) = 4 * a :=
      fun j hj => if_pos (Finset.mem_Ico.mp hj).2
    have htail : ∀ j ∈ Finset.Ico Qm Qp,
        (if j < Qm then 4 * a else 2 * (S - 14 * (j + 1)))
          = 2 * (S - 14 * (j + 1)) :=
      fun j hj => if_neg (not_lt.mpr (Finset.mem_Ico.mp hj).1)
    rw [Finset.sum_congr rfl hplateau, Finset.sum_congr rfl htail,
      Finset.sum_const, Nat.card_Ico, Nat.sub_zero, smul_eq_mul]
    ring
  -- cast the tail sum to ℤ (in the push_cast-distributed shape)
  have hcastsum : (∑ j ∈ Finset.Ico Qm Qp, 2 * ((S - 14 * (j + 1) : ℕ) : ℤ))
      = ∑ j ∈ Finset.Ico Qm Qp, (2 * (S : ℤ) - 28 * ((j : ℤ) + 1)) := by
    refine Finset.sum_congr rfl fun j hj => ?_
    rw [Finset.mem_Ico] at hj
    have h14 : 14 * (j + 1) ≤ S := by omega
    have hcast : ((S - 14 * (j + 1) : ℕ) : ℤ) = (S : ℤ) - 14 * ((j : ℤ) + 1) := by
      push_cast [h14]
      ring
    rw [hcast]
    ring
  have hgoal_lhs : (14 * muNum a b : ℤ)
      = 14 * (2 * (a : ℤ) + 4 * a * Qm)
        + 14 * (2 * (S : ℤ) * ((Qp : ℤ) - Qm)
            - 14 * ((Qp : ℤ) * (Qp + 1) - (Qm : ℤ) * (Qm + 1))) := by
    rw [hnat]
    push_cast
    rw [hcastsum, tail_sum (S : ℤ) Qm Qp hQmp]
    ring
  rw [hgoal_lhs]
  -- the fold algebra: one linear_combination with factored coefficients
  have h2a : (2 * a : ℤ) = (S : ℤ) - (D : ℤ) := by
    rw [hSdef, hDdef]
    push_cast
    omega
  have hrS : (S : ℤ) = 14 * (Qp : ℤ) + ((S % 14 : ℕ) : ℤ) := by
    rw [hQpdef]
    push_cast
    omega
  have hrD : (D : ℤ) = 14 * (Qm : ℤ) + ((D % 14 : ℕ) : ℤ) := by
    rw [hQmdef]
    push_cast
    omega
  have h4ab : 4 * (a : ℤ) * b = (S : ℤ) ^ 2 - (D : ℤ) ^ 2 := by
    have hS' : (S : ℤ) = (a : ℤ) + b := by rw [hSdef]; push_cast; ring
    have hD' : (D : ℤ) = (b : ℤ) - a := by rw [hDdef]; push_cast; omega
    rw [hS', hD']
    ring
  have hSb : (((a + b) % 14 : ℕ) : ℤ) = ((S % 14 : ℕ) : ℤ) := by rw [hSdef]
  have hDb : (((b - a) % 14 : ℕ) : ℤ) = ((D % 14 : ℕ) : ℤ) := by rw [hDdef]
  rw [hSb, hDb, h4ab]
  linear_combination (14 + 28 * (Qm : ℤ)) * h2a
    - ((S : ℤ) - 14 - 14 * (Qp : ℤ) + ((S % 14 : ℕ) : ℤ)) * hrS
    + ((D : ℤ) - 14 - 14 * (Qm : ℤ) + ((D % 14 : ℕ) : ℤ)) * hrD

end LonelyRunner.LRC14.Hunter
