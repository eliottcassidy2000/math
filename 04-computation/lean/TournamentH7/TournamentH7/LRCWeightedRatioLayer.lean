import TournamentH7.LRCPairRatioLayerArithmetic

namespace LonelyRunner
namespace LRCWeightedRatioLayer

open Finset
open LRCPairRatioLayerArithmetic

def tau2 : ℚ := 6 / 637
def tau3 : ℚ := 2 / 441
def tau4 : ℚ := 5 / 2646
def tau5 : ℚ := 1 / 1764
def tau6 : ℚ := 1 / 3136
def tau7 : ℚ := 5 / 37632
def tau8 : ℚ := 1 / 9996
def tau9 : ℚ := 1 / 14112
def ratio11Level : ℚ := 4 / 539
def ratio12Level : ℚ := 5 / 588

def countAbove {α : Type*} [DecidableEq α] (edges : Finset α)
    (weight : α → ℚ) (threshold : ℚ) : ℕ :=
  (edges.filter fun edge => threshold < weight edge).card

def stepAbove (threshold weight : ℚ) : ℚ :=
  if threshold < weight then 1 else 0

def indicatorAbove {α : Type*} (weight : α → ℚ)
    (threshold : ℚ) (edge : α) : ℚ :=
  stepAbove threshold (weight edge)

def layerCeiling (weight : ℚ) : ℚ :=
  tau9
    + (tau8 - tau9) * stepAbove tau9 weight
    + (tau7 - tau8) * stepAbove tau8 weight
    + (tau6 - tau7) * stepAbove tau7 weight
    + (tau5 - tau6) * stepAbove tau6 weight
    + (tau4 - tau5) * stepAbove tau5 weight
    + (tau3 - tau4) * stepAbove tau4 weight
    + (ratio11Level - tau3) * stepAbove tau3 weight
    + (ratio12Level - ratio11Level) *
        stepAbove ratio11Level weight
    + (tau2 - ratio12Level) *
        stepAbove ratio12Level weight

theorem weight_le_layerCeiling (weight : ℚ) (hmax : weight ≤ tau2) :
    weight ≤ layerCeiling weight := by
  have h98 : tau9 < tau8 := by norm_num [tau9, tau8]
  have h87 : tau8 < tau7 := by norm_num [tau8, tau7]
  have h76 : tau7 < tau6 := by norm_num [tau7, tau6]
  have h65 : tau6 < tau5 := by norm_num [tau6, tau5]
  have h54 : tau5 < tau4 := by norm_num [tau5, tau4]
  have h43 : tau4 < tau3 := by norm_num [tau4, tau3]
  have h3r11 : tau3 < ratio11Level := by
    norm_num [tau3, ratio11Level]
  have hr1112 : ratio11Level < ratio12Level := by
    norm_num [ratio11Level, ratio12Level]
  by_cases h12 : ratio12Level < weight
  · have h11 : ratio11Level < weight := hr1112.trans h12
    have h3 : tau3 < weight := h3r11.trans h11
    have h4 : tau4 < weight := h43.trans h3
    have h5 : tau5 < weight := h54.trans h4
    have h6 : tau6 < weight := h65.trans h5
    have h7 : tau7 < weight := h76.trans h6
    have h8 : tau8 < weight := h87.trans h7
    have h9 : tau9 < weight := h98.trans h8
    calc
      weight ≤ tau2 := hmax
      _ = layerCeiling weight := by
        simp only [layerCeiling, stepAbove, if_pos h12, if_pos h11,
          if_pos h3, if_pos h4, if_pos h5, if_pos h6, if_pos h7,
          if_pos h8, if_pos h9]
        ring
  · by_cases h11 : ratio11Level < weight
    · have h3 : tau3 < weight := h3r11.trans h11
      have h4 : tau4 < weight := h43.trans h3
      have h5 : tau5 < weight := h54.trans h4
      have h6 : tau6 < weight := h65.trans h5
      have h7 : tau7 < weight := h76.trans h6
      have h8 : tau8 < weight := h87.trans h7
      have h9 : tau9 < weight := h98.trans h8
      calc
        weight ≤ ratio12Level := le_of_not_gt h12
        _ = layerCeiling weight := by
          simp only [layerCeiling, stepAbove, if_neg h12, if_pos h11,
            if_pos h3, if_pos h4, if_pos h5, if_pos h6, if_pos h7,
            if_pos h8, if_pos h9]
          ring
    · by_cases h3 : tau3 < weight
      · have h4 : tau4 < weight := h43.trans h3
        have h5 : tau5 < weight := h54.trans h4
        have h6 : tau6 < weight := h65.trans h5
        have h7 : tau7 < weight := h76.trans h6
        have h8 : tau8 < weight := h87.trans h7
        have h9 : tau9 < weight := h98.trans h8
        calc
          weight ≤ ratio11Level := le_of_not_gt h11
          _ = layerCeiling weight := by
            simp only [layerCeiling, stepAbove, if_neg h12, if_neg h11,
              if_pos h3, if_pos h4, if_pos h5, if_pos h6, if_pos h7,
              if_pos h8, if_pos h9]
            ring
      · by_cases h4 : tau4 < weight
        · have h5 : tau5 < weight := h54.trans h4
          have h6 : tau6 < weight := h65.trans h5
          have h7 : tau7 < weight := h76.trans h6
          have h8 : tau8 < weight := h87.trans h7
          have h9 : tau9 < weight := h98.trans h8
          calc
            weight ≤ tau3 := le_of_not_gt h3
            _ = layerCeiling weight := by
              simp only [layerCeiling, stepAbove, if_neg h12, if_neg h11,
                if_neg h3, if_pos h4, if_pos h5, if_pos h6, if_pos h7,
                if_pos h8, if_pos h9]
              ring
        · by_cases h5 : tau5 < weight
          · have h6 : tau6 < weight := h65.trans h5
            have h7 : tau7 < weight := h76.trans h6
            have h8 : tau8 < weight := h87.trans h7
            have h9 : tau9 < weight := h98.trans h8
            calc
              weight ≤ tau4 := le_of_not_gt h4
              _ = layerCeiling weight := by
                simp only [layerCeiling, stepAbove, if_neg h12, if_neg h11,
                  if_neg h3, if_neg h4, if_pos h5, if_pos h6, if_pos h7,
                  if_pos h8, if_pos h9]
                ring
          · by_cases h6 : tau6 < weight
            · have h7 : tau7 < weight := h76.trans h6
              have h8 : tau8 < weight := h87.trans h7
              have h9 : tau9 < weight := h98.trans h8
              calc
                weight ≤ tau5 := le_of_not_gt h5
                _ = layerCeiling weight := by
                  simp only [layerCeiling, stepAbove, if_neg h12,
                    if_neg h11, if_neg h3, if_neg h4, if_neg h5,
                    if_pos h6, if_pos h7, if_pos h8, if_pos h9]
                  ring
            · by_cases h7 : tau7 < weight
              · have h8 : tau8 < weight := h87.trans h7
                have h9 : tau9 < weight := h98.trans h8
                calc
                  weight ≤ tau6 := le_of_not_gt h6
                  _ = layerCeiling weight := by
                    simp only [layerCeiling, stepAbove, if_neg h12,
                      if_neg h11, if_neg h3, if_neg h4, if_neg h5,
                      if_neg h6, if_pos h7, if_pos h8, if_pos h9]
                    ring
              · by_cases h8 : tau8 < weight
                · have h9 : tau9 < weight := h98.trans h8
                  calc
                    weight ≤ tau7 := le_of_not_gt h7
                    _ = layerCeiling weight := by
                      simp only [layerCeiling, stepAbove, if_neg h12,
                        if_neg h11, if_neg h3, if_neg h4, if_neg h5,
                        if_neg h6, if_neg h7, if_pos h8, if_pos h9]
                      ring
                · by_cases h9 : tau9 < weight
                  · calc
                      weight ≤ tau8 := le_of_not_gt h8
                      _ = layerCeiling weight := by
                        simp only [layerCeiling, stepAbove, if_neg h12,
                          if_neg h11, if_neg h3, if_neg h4, if_neg h5,
                          if_neg h6, if_neg h7, if_neg h8, if_pos h9]
                        ring
                  · calc
                      weight ≤ tau9 := le_of_not_gt h9
                      _ = layerCeiling weight := by
                        simp only [layerCeiling, stepAbove, if_neg h12,
                          if_neg h11, if_neg h3, if_neg h4, if_neg h5,
                          if_neg h6, if_neg h7, if_neg h8, if_neg h9]
                        ring

theorem sum_indicatorAbove_eq_countAbove {α : Type*} [DecidableEq α]
    (edges : Finset α) (weight : α → ℚ) (threshold : ℚ) :
    ∑ edge ∈ edges, indicatorAbove weight threshold edge =
      (countAbove edges weight threshold : ℚ) := by
  classical
  simpa only [indicatorAbove, stepAbove, countAbove] using
    (Finset.sum_boole (R := ℚ)
      (fun edge => threshold < weight edge) edges)

theorem sum_layerCeiling_eq {α : Type*} [DecidableEq α]
    (edges : Finset α) (weight : α → ℚ) :
    ∑ edge ∈ edges, layerCeiling (weight edge) =
      (edges.card : ℚ) * tau9
        + (tau8 - tau9) * countAbove edges weight tau9
        + (tau7 - tau8) * countAbove edges weight tau8
        + (tau6 - tau7) * countAbove edges weight tau7
        + (tau5 - tau6) * countAbove edges weight tau6
        + (tau4 - tau5) * countAbove edges weight tau5
        + (tau3 - tau4) * countAbove edges weight tau4
        + (ratio11Level - tau3) * countAbove edges weight tau3
        + (ratio12Level - ratio11Level) *
            countAbove edges weight ratio11Level
        + (tau2 - ratio12Level) * countAbove edges weight ratio12Level := by
  classical
  simp only [layerCeiling]
  simp_rw [Finset.sum_add_distrib, ← Finset.mul_sum]
  simp only [stepAbove, Finset.sum_boole, countAbove]
  simp only [Finset.sum_const, nsmul_eq_mul]

/-- Abstract finite weighted-layer adapter.  The seven middle hypotheses are
the Turan edge caps coming from the strict ratio-clique certificates.  The top
two are only the union-of-two-paths and one-path bounds for ratios `12,13`. -/
theorem weighted_pair_sum_le_path_bound {α : Type*} [DecidableEq α]
    (edges : Finset α) (weight : α → ℚ)
    (hmax : ∀ edge ∈ edges, weight edge ≤ tau2)
    (hcard : edges.card ≤ 78)
    (h9 : countAbove edges weight tau9 ≤ 73)
    (h8 : countAbove edges weight tau8 ≤ 72)
    (h7 : countAbove edges weight tau7 ≤ 70)
    (h6 : countAbove edges weight tau6 ≤ 67)
    (h5 : countAbove edges weight tau5 ≤ 63)
    (h4 : countAbove edges weight tau4 ≤ 56)
    (h3 : countAbove edges weight tau3 ≤ 42)
    (h11 : countAbove edges weight ratio11Level ≤ 24)
    (h12 : countAbove edges weight ratio12Level ≤ 12) :
    ∑ edge ∈ edges, weight edge ≤ negativePairTierBoundPathOnly := by
  have hpoint : ∀ edge ∈ edges, weight edge ≤ layerCeiling (weight edge) :=
    fun edge hedge => weight_le_layerCeiling (weight edge) (hmax edge hedge)
  have hsum :
      ∑ edge ∈ edges, weight edge ≤
        ∑ edge ∈ edges, layerCeiling (weight edge) :=
    Finset.sum_le_sum hpoint
  rw [sum_layerCeiling_eq] at hsum
  have hcardQ : (edges.card : ℚ) ≤ 78 := by exact_mod_cast hcard
  have h9Q : (countAbove edges weight tau9 : ℚ) ≤ 73 := by exact_mod_cast h9
  have h8Q : (countAbove edges weight tau8 : ℚ) ≤ 72 := by exact_mod_cast h8
  have h7Q : (countAbove edges weight tau7 : ℚ) ≤ 70 := by exact_mod_cast h7
  have h6Q : (countAbove edges weight tau6 : ℚ) ≤ 67 := by exact_mod_cast h6
  have h5Q : (countAbove edges weight tau5 : ℚ) ≤ 63 := by exact_mod_cast h5
  have h4Q : (countAbove edges weight tau4 : ℚ) ≤ 56 := by exact_mod_cast h4
  have h3Q : (countAbove edges weight tau3 : ℚ) ≤ 42 := by exact_mod_cast h3
  have h11Q : (countAbove edges weight ratio11Level : ℚ) ≤ 24 := by
    exact_mod_cast h11
  have h12Q : (countAbove edges weight ratio12Level : ℚ) ≤ 12 := by
    exact_mod_cast h12
  apply hsum.trans
  calc
    _ ≤ 78 * tau9
        + (tau8 - tau9) * 73
        + (tau7 - tau8) * 72
        + (tau6 - tau7) * 70
        + (tau5 - tau6) * 67
        + (tau4 - tau5) * 63
        + (tau3 - tau4) * 56
        + (ratio11Level - tau3) * 42
        + (ratio12Level - ratio11Level) * 24
        + (tau2 - ratio12Level) * 12 := by
          gcongr <;>
            norm_num [tau2, tau3, tau4, tau5, tau6, tau7, tau8, tau9,
              ratio11Level, ratio12Level]
    _ = negativePairTierBoundPathOnly := by
      norm_num [tau2, tau3, tau4, tau5, tau6, tau7, tau8, tau9,
        ratio11Level, ratio12Level, negativePairTierBoundPathOnly]

theorem weighted_pair_sum_lt_target {α : Type*} [DecidableEq α]
    (edges : Finset α) (weight : α → ℚ)
    (hmax : ∀ edge ∈ edges, weight edge ≤ tau2)
    (hcard : edges.card ≤ 78)
    (h9 : countAbove edges weight tau9 ≤ 73)
    (h8 : countAbove edges weight tau8 ≤ 72)
    (h7 : countAbove edges weight tau7 ≤ 70)
    (h6 : countAbove edges weight tau6 ≤ 67)
    (h5 : countAbove edges weight tau5 ≤ 63)
    (h4 : countAbove edges weight tau4 ≤ 56)
    (h3 : countAbove edges weight tau3 ≤ 42)
    (h11 : countAbove edges weight ratio11Level ≤ 24)
    (h12 : countAbove edges weight ratio12Level ≤ 12) :
    ∑ edge ∈ edges, weight edge < 13 / 30 := by
  have hbound := weighted_pair_sum_le_path_bound edges weight hmax hcard
    h9 h8 h7 h6 h5 h4 h3 h11 h12
  exact hbound.trans_lt negativePairTierBoundPathOnly_lt_target

#print axioms weight_le_layerCeiling
#print axioms sum_layerCeiling_eq
#print axioms weighted_pair_sum_le_path_bound
#print axioms weighted_pair_sum_lt_target

end LRCWeightedRatioLayer
end LonelyRunner
