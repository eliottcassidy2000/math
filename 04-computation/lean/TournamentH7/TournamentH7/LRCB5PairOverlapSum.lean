import TournamentH7.LRCClosedBudget
import TournamentH7.LRCRationalOpenComb

/-!
Exact Raabe/Bernoulli evaluation of the cyclic tooth-pair overlap sum.  This
discharges the analytic arithmetic after a single explicit finite reindexing:
the normalized merged-comb length must equal the gcd-scaled overlap ledger.
-/

namespace LonelyRunner
namespace LRCB5PairOverlapSum

open LRC14.ClosedBudget
open LRCRationalOpenComb RatIntervals

noncomputable section

/-- Cyclic grid sum of two interval-overlap tents. -/
def gridPairOverlap (period : ℕ) (length₁ length₂ shift : ℝ) : ℝ :=
  ∑ index ∈ Finset.range period,
    pairOverlap length₁ length₂
      (Int.fract (shift + (index : ℝ) / period))

/-- Raabe evaluation of a complete cyclic grid of overlap tents. -/
theorem gridPairOverlap_eq_B2
    (period : ℕ) (hperiod : 1 ≤ period)
    (length₁ length₂ shift : ℝ)
    (hlength₁0 : 0 ≤ length₁) (hlength₁1 : length₁ ≤ 1)
    (hlength₂0 : 0 ≤ length₂) (hlength₂1 : length₂ ≤ 1) :
    gridPairOverlap period length₁ length₂ shift =
      period * (length₁ * length₂) +
        (1 / (2 * period)) *
          (B2R (period * shift) -
            B2R (period * (shift - length₂)) -
            B2R (period * (shift + length₁)) +
            B2R (period * (shift + length₁ - length₂))) := by
  have hterm : ∀ index ∈ Finset.range period,
      pairOverlap length₁ length₂
          (Int.fract (shift + (index : ℝ) / period)) =
        length₁ * length₂ + (1 / 2) *
          (B2R (shift + (index : ℝ) / period) -
            B2R (shift - length₂ + (index : ℝ) / period) -
            B2R (shift + length₁ + (index : ℝ) / period) +
            B2R (shift + length₁ - length₂ + (index : ℝ) / period)) := by
    intro index _hindex
    have hoverlap := pair_overlap_B2 length₁ length₂
      (Int.fract (shift + (index : ℝ) / period))
      hlength₁0 hlength₁1 hlength₂0 hlength₂1
      (Int.fract_nonneg _) (Int.fract_lt_one _)
    unfold pairOverlap
    rw [hoverlap]
    have hone :
        B2R (Int.fract (shift + (index : ℝ) / period)) =
          B2R (shift + (index : ℝ) / period) := by
      simpa using B2R_fract_add (shift + (index : ℝ) / period) 0
    have htwo :
        B2R (Int.fract (shift + (index : ℝ) / period) - length₂) =
          B2R (shift - length₂ + (index : ℝ) / period) := by
      rw [sub_eq_add_neg, B2R_fract_add]
      congr 1
      ring
    have hthree :
        B2R (Int.fract (shift + (index : ℝ) / period) + length₁) =
          B2R (shift + length₁ + (index : ℝ) / period) := by
      rw [B2R_fract_add]
      congr 1
      ring
    have hfour :
        B2R (Int.fract (shift + (index : ℝ) / period) + length₁ - length₂) =
          B2R (shift + length₁ - length₂ + (index : ℝ) / period) := by
      rw [show Int.fract (shift + (index : ℝ) / period) + length₁ - length₂ =
          Int.fract (shift + (index : ℝ) / period) + (length₁ - length₂) by ring,
        B2R_fract_add]
      congr 1
      ring
    rw [hone, htwo, hthree, hfour]
  have hsum₁ := raabe_B2 period hperiod shift
  have hsum₂ := raabe_B2 period hperiod (shift - length₂)
  have hsum₃ := raabe_B2 period hperiod (shift + length₁)
  have hsum₄ := raabe_B2 period hperiod (shift + length₁ - length₂)
  unfold gridPairOverlap
  rw [Finset.sum_congr rfl hterm, Finset.sum_add_distrib,
    Finset.sum_const, Finset.card_range, nsmul_eq_mul]
  have hcorrection :
      (∑ index ∈ Finset.range period,
        (B2R (shift + (index : ℝ) / period) -
          B2R (shift - length₂ + (index : ℝ) / period) -
          B2R (shift + length₁ + (index : ℝ) / period) +
          B2R (shift + length₁ - length₂ + (index : ℝ) / period))) =
        (1 / period) *
          (B2R (period * shift) -
            B2R (period * (shift - length₂)) -
            B2R (period * (shift + length₁)) +
            B2R (period * (shift + length₁ - length₂))) := by
    rw [Finset.sum_add_distrib, Finset.sum_sub_distrib,
      Finset.sum_sub_distrib, hsum₁, hsum₂, hsum₃, hsum₄]
    ring
  rw [← Finset.mul_sum, hcorrection]
  have hperiod0 : (period : ℝ) ≠ 0 := by positivity
  field_simp

/-- The cyclic overlap ledger for speeds `g*p` and `g*q`; `g` is the
multiplicity of each difference class and `g*p*q` is their lcm when `p,q`
are coprime. -/
def scaledPairOverlapLedger (g p q : ℕ) : ℝ :=
  (g : ℝ) *
    gridPairOverlap (g * p * q)
      (1 / (7 * (g * p : ℕ) : ℝ))
      (1 / (7 * (g * q : ℕ) : ℝ))
      (1 / (14 * (g * q : ℕ) : ℝ) -
        1 / (14 * (g * p : ℕ) : ℝ))

/-- Exact B₂ evaluation of the gcd-scaled cyclic tooth-pair ledger. -/
theorem scaledPairOverlapLedger_sub_eq_B2
    (g p q : ℕ) (hg : 0 < g) (hp : 0 < p) (hq : 0 < q) :
    scaledPairOverlapLedger g p q - 1 / 49 =
      (B2R (((p : ℝ) - q) / 14) -
        B2R (((p : ℝ) + q) / 14)) /
          ((p : ℝ) * q) := by
  have hgpNat : 0 < g * p := Nat.mul_pos hg hp
  have hgqNat : 0 < g * q := Nat.mul_pos hg hq
  have hperiod : 1 ≤ g * p * q := Nat.mul_pos hgpNat hq
  have hgp : (0 : ℝ) < (g * p : ℕ) := by exact_mod_cast hgpNat
  have hgq : (0 : ℝ) < (g * q : ℕ) := by exact_mod_cast hgqNat
  have hlength₁0 :
      0 ≤ (1 / (7 * (g * p : ℕ) : ℝ)) := by positivity
  have hlength₂0 :
      0 ≤ (1 / (7 * (g * q : ℕ) : ℝ)) := by positivity
  have hlength₁1 :
      (1 / (7 * (g * p : ℕ) : ℝ)) ≤ 1 := by
    have hdenNat : 1 ≤ 7 * (g * p) := by omega
    have hden : (1 : ℝ) ≤ 7 * (g * p : ℕ) := by exact_mod_cast hdenNat
    exact (div_le_one (by positivity)).2 hden
  have hlength₂1 :
      (1 / (7 * (g * q : ℕ) : ℝ)) ≤ 1 := by
    have hdenNat : 1 ≤ 7 * (g * q) := by omega
    have hden : (1 : ℝ) ≤ 7 * (g * q : ℕ) := by exact_mod_cast hdenNat
    exact (div_le_one (by positivity)).2 hden
  rw [scaledPairOverlapLedger,
    gridPairOverlap_eq_B2 (g * p * q) hperiod
      (1 / (7 * (g * p : ℕ) : ℝ))
      (1 / (7 * (g * q : ℕ) : ℝ))
      (1 / (14 * (g * q : ℕ) : ℝ) -
        1 / (14 * (g * p : ℕ) : ℝ))
      hlength₁0 hlength₁1 hlength₂0 hlength₂1]
  have harg₁ :
      ((g * p * q : ℕ) : ℝ) *
          (1 / (14 * (g * q : ℕ) : ℝ) -
            1 / (14 * (g * p : ℕ) : ℝ)) =
        ((p : ℝ) - q) / 14 := by
    push_cast
    field_simp
  have harg₂ :
      ((g * p * q : ℕ) : ℝ) *
          ((1 / (14 * (g * q : ℕ) : ℝ) -
            1 / (14 * (g * p : ℕ) : ℝ)) -
              1 / (7 * (g * q : ℕ) : ℝ)) =
        -(((p : ℝ) + q) / 14) := by
    push_cast
    field_simp
    ring
  have harg₃ :
      ((g * p * q : ℕ) : ℝ) *
          ((1 / (14 * (g * q : ℕ) : ℝ) -
            1 / (14 * (g * p : ℕ) : ℝ)) +
              1 / (7 * (g * p : ℕ) : ℝ)) =
        ((p : ℝ) + q) / 14 := by
    push_cast
    field_simp
    ring
  have harg₄ :
      ((g * p * q : ℕ) : ℝ) *
          ((1 / (14 * (g * q : ℕ) : ℝ) -
            1 / (14 * (g * p : ℕ) : ℝ)) +
              1 / (7 * (g * p : ℕ) : ℝ) -
              1 / (7 * (g * q : ℕ) : ℝ)) =
        -(((p : ℝ) - q) / 14) := by
    push_cast
    field_simp
    ring
  rw [harg₁, harg₂, harg₃, harg₄,
    B2R_neg (((p : ℝ) + q) / 14),
    B2R_neg (((p : ℝ) - q) / 14)]
  push_cast
  field_simp
  ring

/-- The full rational-list identity is reduced to one exact finite
reindexing statement: the merged pair-region length equals the cyclic
tooth-pair overlap ledger. Everything after that reindexing is proved. -/
theorem length_ratOpenPairRegion_sub_eq_B2_of_scaledLedger
    (first second : ℕ) (hfirst : 0 < first) (hsecond : 0 < second)
    (hledger :
      ((length (ratOpenPairRegion first second) : ℚ) : ℝ) =
        scaledPairOverlapLedger (Nat.gcd first second)
          (first / Nat.gcd first second)
          (second / Nat.gcd first second)) :
    ((length (ratOpenPairRegion first second) : ℚ) : ℝ) - 1 / 49 =
      (B2R ((((first / Nat.gcd first second : ℕ) : ℝ) -
          (second / Nat.gcd first second : ℕ)) / 14) -
        B2R ((((first / Nat.gcd first second : ℕ) : ℝ) +
          (second / Nat.gcd first second : ℕ)) / 14)) /
      ((first / Nat.gcd first second : ℕ) *
        (second / Nat.gcd first second : ℕ)) := by
  have hg : 0 < Nat.gcd first second :=
    Nat.gcd_pos_of_pos_left second hfirst
  have hp : 0 < first / Nat.gcd first second :=
    Nat.div_pos (Nat.gcd_le_left second hfirst) hg
  have hq : 0 < second / Nat.gcd first second :=
    Nat.div_pos (Nat.gcd_le_right first hsecond) hg
  rw [hledger]
  exact scaledPairOverlapLedger_sub_eq_B2
    (Nat.gcd first second)
    (first / Nat.gcd first second)
    (second / Nat.gcd first second) hg hp hq

#print axioms gridPairOverlap_eq_B2
#print axioms scaledPairOverlapLedger_sub_eq_B2
#print axioms length_ratOpenPairRegion_sub_eq_B2_of_scaledLedger

end
end LRCB5PairOverlapSum
end LonelyRunner
