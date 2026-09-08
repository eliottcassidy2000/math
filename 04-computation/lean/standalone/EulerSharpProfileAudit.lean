import Mathlib.Tactic

/-! Independent scalar certificate for a refinement of the explicit periodic
profile in Euler/EulerProof.lean. This file does not import the Euler project
and does not certify a PDE theorem. -/

noncomputable section

namespace SharpProfileAudit

def denominator (δ c : ℝ) : ℝ := (1 + δ)^2 - 2*(1 + δ)*c + 1

def slope (δ c : ℝ) : ℝ := ((1 + δ)*c - 1) / denominator δ c

theorem denominator_pos {δ c : ℝ} (hδ : 0 < δ) (hc : c ≤ 1) :
    0 < denominator δ c := by
  have h := mul_nonneg (by positivity : 0 ≤ 2*(1+δ)) (sub_nonneg.mpr hc)
  unfold denominator
  nlinarith [sq_pos_of_pos hδ]

theorem slope_lower {δ c : ℝ} (hδ : 0 < δ) (hcl : -1 ≤ c) (hcu : c ≤ 1) :
    (-1 : ℝ)/(δ+2) ≤ slope δ c := by
  unfold slope
  apply (div_le_div_iff₀ (by linarith : 0 < δ+2) (denominator_pos hδ hcu)).2
  have h := mul_nonneg (by positivity : 0 ≤ δ*(1+δ)) (by linarith : 0 ≤ c+1)
  unfold denominator
  nlinarith

theorem slope_upper {δ c : ℝ} (hδ : 0 < δ) (hc : c ≤ 1) :
    slope δ c ≤ 1/δ := by
  unfold slope
  apply (div_le_div_iff₀ (denominator_pos hδ hc) hδ).2
  have h := mul_nonneg (by positivity : 0 ≤ (1+δ)*(2+δ)) (sub_nonneg.mpr hc)
  unfold denominator
  nlinarith

theorem slope_at_neg_one {δ : ℝ} (hδ : 0 < δ) :
    slope δ (-1) = (-1 : ℝ)/(δ+2) := by
  have h1 : δ+2 ≠ 0 := by linarith
  have hd : denominator δ (-1) = (δ+2)^2 := by unfold denominator; ring
  unfold slope
  rw [hd]
  field_simp
  ring

theorem slope_at_one {δ : ℝ} (hδ : 0 < δ) : slope δ 1 = 1/δ := by
  have h1 : δ ≠ 0 := ne_of_gt hδ
  have hd : denominator δ 1 = δ^2 := by unfold denominator; ring
  unfold slope
  rw [hd]
  field_simp
  ring

theorem sharp_scalar_pressure_cost {a flux δ c : ℝ}
    (ha : 0 ≤ a) (hf : 0 ≤ flux) (hδ : 0 < δ) (hcl : -1 ≤ c) (hcu : c ≤ 1) :
    -2*a*flux*slope δ c ≤ 2*a*flux/(δ+2) := by
  have h := mul_le_mul_of_nonneg_left (slope_lower hδ hcl hcu)
    (by positivity : 0 ≤ 2*a*flux)
  calc
    -2*a*flux*slope δ c = -(2*a*flux*slope δ c) := by ring
    _ ≤ -(2*a*flux*((-1 : ℝ)/(δ+2))) := neg_le_neg h
    _ = 2*a*flux/(δ+2) := by ring

end SharpProfileAudit

#print axioms SharpProfileAudit.slope_lower
#print axioms SharpProfileAudit.slope_upper
#print axioms SharpProfileAudit.sharp_scalar_pressure_cost
