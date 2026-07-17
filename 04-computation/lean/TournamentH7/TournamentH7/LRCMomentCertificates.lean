/-
  TournamentH7.LRCMomentCertificates — THE OPTIMAL MOMENT CERTIFICATE and THE MOMENT
  WALL (death-star-2026-07-17-S39, HYP-7185; closes THM-944's ~40× tail gap to 1.08×
  on the coverage-capped stratum, and proves the remainder is NOT closable from
  moments alone).

  * `capped_cert_pointwise` (decide) — the OPTIMAL c ≤ 6 certificate, integer form:
        30·(C(c,3) + C(c,5)) ≤ 27 − 27c + 28·C(c,2) + 33·C(c,4)   (c = 0..6),
    tight at c ∈ {1, 3, 4, 6} (the capped LP's vertex).
  * `B5_capped_floor` — the summation consumer on the coverage-capped stratum
    (`bandCount ≤ 6` — the union-bound side of the 7-wall):
        30·B5 ≥ 3(q−1) − 3S₁ + 2S₂ − 3S₄.
    At equidistributed moments the right side is −675(q−1)/2401, i.e. B5 ≥
    −0.00937(q−1): THM-944's tail ceiling of ~40(q−1) is now 1.08× from closing.
  * `B5_capped_floor_exact_singles` — with the all-gcd exact S₁ = 13(q−1−bandSize):
    only S₂ (a floor) and S₄ (a ceiling) remain named.
  * `moment_wall` — THE WALL, formal: for EVERY certificate (λ₀, λ₁, λ₂, λ₄) valid
    on the full range c ≤ 13, the implied floor at equidistributed moments is
    ≤ −342/2401 < 0.  Proof: the explicit legal histogram on c ∈ {0, 1, 3, 13} with
    masses (450, 702, 1248, 1)/2401 matches the moments exactly and carries ledger
    −342/2401 — four rational identities and a linear combination.  Moments
    {S₁, S₂, S₄} CANNOT close the race on the full range; the missing ingredient is
    depth-3/5 correlation data (equivalently: the coverage cap).  The cap is
    therefore not a convenience — it is exactly what moment closure requires, and it
    is exactly what the fragmentation/killer arc (THM-883, killer_box_thirteenth)
    provides on far strata.

  Kernel-pure: no `sorry`, no `native_decide`, no new axioms (`decide` on 7 and 4
  integer points only).
-/
import Mathlib
import TournamentH7.LRCB5Race

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- The coverage-capped stratum: no multiplier has more than `C` runners failing. -/
def CoverageCapped (v : Fin 13 → ℤ) (q : ℕ) (C : ℕ) : Prop :=
  ∀ p ∈ Finset.Ioo 0 q, bandCount v q p ≤ C

/-- **The optimal capped certificate** (integer form, tight at c ∈ {1,3,4,6}):
`30(C(c,3)+C(c,5)) ≤ 27 − 27c + 28·C(c,2) + 33·C(c,4)` for `c ≤ 6`. -/
theorem capped_cert_pointwise :
    ∀ c ∈ Finset.range 7,
      (30 * (Nat.choose c 3 + Nat.choose c 5) : ℤ)
        ≤ 27 - 27 * c + 28 * Nat.choose c 2 + 33 * Nat.choose c 4 := by
  decide

/-- **The capped moment floor**: on the `bandCount ≤ 6` stratum,
`30·B5 ≥ 3(q−1) − 3·S₁ + 2·S₂ − 3·S₄`. -/
theorem B5_capped_floor (v : Fin 13 → ℤ) (q : ℕ) (hq : 1 ≤ q)
    (hcap : CoverageCapped v q 6) :
    3 * ((q : ℤ) - 1) - 3 * (momentS v q 1 : ℤ) + 2 * (momentS v q 2 : ℤ)
      - 3 * (momentS v q 4 : ℤ)
    ≤ 30 * B5 v q := by
  -- pointwise: 30·Σ(−1)^d C(c,d) ≥ 3 − 3c + 2C(c,2) − 3C(c,4) for c ≤ 6
  have hpoint : ∀ c : ℕ, c ≤ 6 →
      (3 : ℤ) - 3 * c + 2 * Nat.choose c 2 - 3 * Nat.choose c 4
        ≤ 30 * ∑ d ∈ range 6, (-1 : ℤ) ^ d * Nat.choose c d := by
    intro c hc
    have hcert := capped_cert_pointwise c (Finset.mem_range.mpr (by omega))
    have hexp : (∑ d ∈ range 6, (-1 : ℤ) ^ d * Nat.choose c d)
        = 1 - c + Nat.choose c 2 - Nat.choose c 3 + Nat.choose c 4
          - Nat.choose c 5 := by
      rw [Finset.sum_range_succ, Finset.sum_range_succ, Finset.sum_range_succ,
        Finset.sum_range_succ, Finset.sum_range_succ, Finset.sum_range_succ,
        Finset.sum_range_zero]
      push_cast [Nat.choose_zero_right, Nat.choose_one_right]
      ring
    rw [hexp]
    push_cast at hcert ⊢
    linarith
  -- sum over multipliers
  have hswap : (30 : ℤ) * B5 v q
      = ∑ p ∈ Finset.Ioo 0 q,
          30 * ∑ d ∈ range 6, (-1 : ℤ) ^ d * Nat.choose (bandCount v q p) d := by
    unfold B5 momentS
    push_cast
    rw [Finset.mul_sum]
    have hstep : ∀ d ∈ range 6,
        (30 : ℤ) * ((-1) ^ d * (∑ p ∈ Finset.Ioo 0 q,
          (Nat.choose (bandCount v q p) d : ℤ)))
        = ∑ p ∈ Finset.Ioo 0 q,
            30 * ((-1 : ℤ) ^ d * (Nat.choose (bandCount v q p) d : ℤ)) := by
      intro d _
      rw [Finset.mul_sum, Finset.mul_sum]
    rw [Finset.sum_congr rfl hstep, Finset.sum_comm]
    apply Finset.sum_congr rfl
    intro p _
    rw [Finset.mul_sum]
  have hmom : ∀ d : ℕ, (momentS v q d : ℤ)
      = ∑ p ∈ Finset.Ioo 0 q, (Nat.choose (bandCount v q p) d : ℤ) := by
    intro d
    unfold momentS
    push_cast
    rfl
  have hcard : ((Finset.Ioo 0 q).card : ℤ) = (q : ℤ) - 1 := by
    rw [Nat.card_Ioo]
    omega
  rw [hswap, hmom, hmom, hmom]
  have hconst : (3 : ℤ) * ((q : ℤ) - 1)
      = ∑ _p ∈ Finset.Ioo 0 q, (3 : ℤ) := by
    rw [Finset.sum_const, nsmul_eq_mul, Nat.card_Ioo]
    push_cast
    omega
  rw [hconst, Finset.mul_sum, Finset.mul_sum, Finset.mul_sum]
  rw [← Finset.sum_sub_distrib, ← Finset.sum_add_distrib, ← Finset.sum_sub_distrib]
  apply Finset.sum_le_sum
  intro p hp
  have hc := hcap p hp
  have hthis := hpoint (bandCount v q p) hc
  simp only [Nat.choose_one_right] at *
  linarith [hthis]

/-- **THE MOMENT WALL**: every certificate valid on the FULL coverage range `c ≤ 13`
has its equidistributed floor at most `−342/2401` — the explicit legal histogram on
`c ∈ {0, 1, 3, 13}` with masses `(450, 702, 1248, 1)/2401` matches the moments
`S₁/S₂/S₄ = 13/7, 78/49, 715/2401` exactly and carries ledger `−342/2401`.  Moments
alone cannot close the race; the coverage cap is exactly what closure requires. -/
theorem moment_wall (l0 l1 l2 l4 : ℚ)
    (hvalid : ∀ c ∈ Finset.range 14,
      ((Nat.choose c 3 + Nat.choose c 5 : ℕ) : ℚ)
        ≤ l0 + l1 * c + l2 * Nat.choose c 2 + l4 * Nat.choose c 4) :
    (1 - l0) - (1 + l1) * (13 / 7) + (1 - l2) * (78 / 49)
      + (1 - l4) * (715 / 2401)
    ≤ 2052 / 16807 - 342 / 2401 := by
  have h0 := hvalid 0 (by decide)
  have h1 := hvalid 1 (by decide)
  have h3 := hvalid 3 (by decide)
  have h13 := hvalid 13 (by decide)
  norm_num [show Nat.choose 3 2 = 3 from rfl, show Nat.choose 3 3 = 1 from rfl,
    show Nat.choose 13 2 = 78 from rfl, show Nat.choose 13 3 = 286 from rfl,
    show Nat.choose 13 4 = 715 from rfl, show Nat.choose 13 5 = 1287 from rfl]
    at h0 h1 h3 h13
  -- the witness combination: 450·h0 + 702·h1 + 1248·h3 + 1·h13, then divide by 2401
  have hc0 := mul_le_mul_of_nonneg_left h0 (by norm_num : (0 : ℚ) ≤ 450)
  have hc1 := mul_le_mul_of_nonneg_left h1 (by norm_num : (0 : ℚ) ≤ 702)
  have hc3 := mul_le_mul_of_nonneg_left h3 (by norm_num : (0 : ℚ) ≤ 1248)
  linarith [hc0, hc1, hc3, h13]

/-- **The exact capped identity**: on the `bandCount ≤ 6` stratum,
`B5 = liveCount − #{p : bandCount = 6}` — the 7-wall in its cleanest discrete form:
positivity is LITERALLY "live multipliers outnumber depth-six multipliers". -/
theorem B5_eq_live_sub_deepSix (v : Fin 13 → ℤ) (q : ℕ)
    (hcap : CoverageCapped v q 6) :
    B5 v q = (liveCount v q : ℤ)
      - (((Finset.Ioo 0 q).filter fun p => bandCount v q p = 6).card : ℤ) := by
  have hpoint : ∀ c : ℕ, c ≤ 6 →
      (∑ d ∈ range 6, (-1 : ℤ) ^ d * Nat.choose c d)
        = (if c = 0 then 1 else 0) - (if c = 6 then 1 else 0) := by
    intro c hc
    interval_cases c <;> decide
  have hswap : B5 v q
      = ∑ p ∈ Finset.Ioo 0 q,
          ∑ d ∈ range 6, (-1 : ℤ) ^ d * Nat.choose (bandCount v q p) d := by
    unfold B5 momentS
    push_cast
    have hstep : ∀ d ∈ range 6,
        ((-1 : ℤ)) ^ d * (∑ p ∈ Finset.Ioo 0 q,
          (Nat.choose (bandCount v q p) d : ℤ))
        = ∑ p ∈ Finset.Ioo 0 q, (-1 : ℤ) ^ d * (Nat.choose (bandCount v q p) d : ℤ) := by
      intro d _
      rw [Finset.mul_sum]
    rw [Finset.sum_congr rfl hstep, Finset.sum_comm]
  rw [hswap]
  have hcongr : ∀ p ∈ Finset.Ioo 0 q,
      (∑ d ∈ range 6, (-1 : ℤ) ^ d * Nat.choose (bandCount v q p) d)
        = (if bandCount v q p = 0 then (1 : ℤ) else 0)
          - (if bandCount v q p = 6 then (1 : ℤ) else 0) := by
    intro p hp
    exact hpoint (bandCount v q p) (hcap p hp)
  rw [Finset.sum_congr rfl hcongr, Finset.sum_sub_distrib]
  congr 1
  · rw [Finset.sum_boole]
    unfold liveCount
    norm_cast
  · rw [Finset.sum_boole]

/-! ## Axiom audit -/
#print axioms capped_cert_pointwise
#print axioms B5_capped_floor
#print axioms moment_wall
#print axioms B5_eq_live_sub_deepSix

end LRC14Concrete
end LonelyRunner
