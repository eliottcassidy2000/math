/- DRAFT (opus S345): the structural steps are COMPLETE (hterm: the cap-boundary
   term shapes; hsplit: the three-piece Ico decomposition); the remaining sorries
   are Gauss-sum plumbing only (Sum_{Ico} j closed form + the final
   linear_combination against hr/hrm).  LANDING PLAN (one sitting): replace the
   explicit Gauss route with Nat.le_induction on Qp fixing Qm (each increment
   adds 2(S - 14 Qp); the target difference closes by ring with hr), OR reuse
   boxeph's sum_shifted pattern from consecutive_closed_form.  The paper identity
   is THM-965 (400/400 exact).  Keep OUT of the build tree until sorry-free. -/
/- LRCFoldedIdentity.lean -- opus-2026-07-17-S345 (HYP-7300 / THM-965).
   THE TWO-VARIABLE FOLDED IDENTITY over boxeph's muNum (coprime-normalized
   pairs; the gcd case rescales):

     14 * muNum a b = 4*a*b + (a+b)%14 * (14 - (a+b)%14)
                            - (b-a)%14 * (14 - (b-a)%14)      (a ≤ b)

   stated in ℤ to keep the subtraction honest.  Proof: split the muNum sum
   at the cap boundary Q₋ = (b-a)/14 — terms with j ≤ Q₋ are capped at 4a
   (the j = 0 term at 2a), terms above descend as 2(a+b) - 28j — then the
   arithmetic telescope: both partial sums are quadratics in Q± and the
   residues fold.  Paper: THM-965 (verified 400/400 exact, every gcd).
   Kernel-pure target: no sorry, no native_decide. -/
import Mathlib
import TournamentH7.LRCTreeHunter

namespace LonelyRunner.LRC14.Hunter

/-- **THM-965, muNum form** (coprime-normalized; `a ≤ b`, `1 ≤ a`):
`14·muNum a b = 4ab + fold((a+b) % 14) − fold((b−a) % 14)`,
`fold r = r(14 − r)`. -/
theorem muNum_folded (a b : ℕ) (ha : 1 ≤ a) (hab : a ≤ b) :
    (14 * muNum a b : ℤ)
      = 4 * a * b + ((a + b) % 14 : ℕ) * (14 - ((a + b) % 14 : ℕ))
        - ((b - a) % 14 : ℕ) * (14 - ((b - a) % 14 : ℕ)) := by
  -- notation: S = a+b, D = b−a (ℕ), quotients/residues
  set S := a + b with hS
  set D := b - a with hD
  have hDS : D ≤ S := by omega
  have hmin : min a b = a := Nat.min_eq_left hab
  have hcap : 2 * min a b = 2 * a := by rw [hmin]
  -- the sum splits at Q₋ = D / 14; total range is Q₊ = S / 14
  set Qp := S / 14 with hQp
  set Qm := D / 14 with hQm
  have hQmp : Qm ≤ Qp := Nat.div_le_div_right hDS
  -- term shapes
  have hterm : ∀ j ∈ Finset.range (Qp + 1),
      (if j = 0 then min S (2 * min a b)
       else 2 * min (S - 14 * j) (2 * min a b))
      = if j = 0 then 2 * a
        else if j ≤ Qm then 4 * a
        else 2 * (S - 14 * j) := by
    intro j hj
    rw [Finset.mem_range] at hj
    rcases Nat.eq_zero_or_pos j with rfl | hjpos
    · simp [hcap]
      omega
    · rw [if_neg (Nat.pos_iff_ne_zero.mp hjpos)]
      rcases le_or_lt j Qm with hle | hgt
      · rw [if_neg (Nat.pos_iff_ne_zero.mp hjpos), if_pos hle, hcap]
        have : 14 * j ≤ D := le_trans (Nat.mul_le_mul_left 14 hle)
          (Nat.mul_div_le D 14 ▸ by
            have := Nat.div_mul_le_self D 14
            omega)
        omega
      · rw [if_neg (Nat.pos_iff_ne_zero.mp hjpos), if_neg (not_le.mpr hgt),
          hcap]
        have h1 : D < 14 * j := by
          have := Nat.lt_div_iff_mul_lt (by norm_num : 0 < 14) |>.mp
          omega
        omega
  unfold muNum
  rw [← hS, Finset.sum_congr rfl hterm]
  -- evaluate the three-piece sum
  have hsplit : ∑ j ∈ Finset.range (Qp + 1),
      (if j = 0 then 2 * a
       else if j ≤ Qm then 4 * a else 2 * (S - 14 * j))
      = 2 * a + 4 * a * Qm
        + ∑ j ∈ Finset.Ico (Qm + 1) (Qp + 1), 2 * (S - 14 * j) := by
    rw [Finset.range_eq_Ico,
      show Finset.Ico 0 (Qp + 1)
        = Finset.Ico 0 1 ∪ Finset.Ico 1 (Qm + 1) ∪ Finset.Ico (Qm + 1) (Qp + 1)
        from by
          rw [Finset.Ico_union_Ico_eq_Ico (by omega) (by omega),
            Finset.Ico_union_Ico_eq_Ico (by omega) (by omega)]]
    rw [Finset.sum_union, Finset.sum_union]
    · congr 1
      · congr 1
        · simp
        · rw [Finset.sum_congr rfl
            (fun j hj => by
              rw [Finset.mem_Ico] at hj
              rw [if_neg (by omega), if_pos (by omega)]),
            Finset.sum_const, Nat.card_Ico]
          ring_nf
          omega
      · exact Finset.sum_congr rfl fun j hj => by
          rw [Finset.mem_Ico] at hj
          rw [if_neg (by omega), if_neg (by omega)]
    · exact Finset.Ico_disjoint_Ico_consecutive 0 1 (Qm + 1) |>.mono_right
        (by exact Finset.Ico_subset_Ico_left (by omega))
    · exact Finset.Ico_disjoint_Ico_consecutive 0 (Qm + 1) (Qp + 1)
  rw [hsplit]
  -- the descending tail: Σ_{Qm+1}^{Qp} (2S − 28j) in ℤ, then residues fold
  have htail : (∑ j ∈ Finset.Ico (Qm + 1) (Qp + 1), 2 * (S - 14 * j) : ℤ)
      = ∑ j ∈ Finset.Ico (Qm + 1) (Qp + 1), (2 * S - 28 * j : ℤ) := by
    push_cast
    refine Finset.sum_congr rfl fun j hj => ?_
    rw [Finset.mem_Ico] at hj
    have : 14 * j ≤ S := by
      have := Nat.div_mul_le_self S 14
      have hj2 : j ≤ Qp := by omega
      calc 14 * j ≤ 14 * Qp := Nat.mul_le_mul_left 14 hj2
        _ ≤ S := by rw [hQp]; omega
    omega
  have hr : S % 14 + 14 * Qp = S := by rw [hQp]; omega
  have hrm : D % 14 + 14 * Qm = D := by rw [hQm]; omega
  have hsum_lin : (∑ j ∈ Finset.Ico (Qm + 1) (Qp + 1), (2 * S - 28 * j : ℤ))
      = 2 * S * (Qp - Qm) - 28 * ((Qp * (Qp + 1) - Qm * (Qm + 1)) / 2 : ℤ)
        := by
    rw [Finset.sum_sub_distrib, Finset.sum_const, Nat.card_Ico,
      ← Finset.sum_mul, ← Finset.mul_sum]
    have hgauss : (∑ j ∈ Finset.Ico (Qm + 1) (Qp + 1), (j : ℤ))
        = ((Qp * (Qp + 1) - Qm * (Qm + 1)) / 2 : ℤ) := by
      have := Finset.sum_range_id_mul_two (Qp + 1)
      have := Finset.sum_range_id_mul_two (Qm + 1)
      sorry
    sorry
  sorry

end LonelyRunner.LRC14.Hunter
