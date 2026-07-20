/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Codex (LRC multi-agent project, 2026-07-19)
-/
import TournamentH7.LRCCenteredSurvivorProtrusion

/-!
# Terminal endpoint transfer and gcd tax (THM-1283)

This module checks the arithmetic consumer behind the paper endpoint topology:
the signed endpoint residual, the outward-width bound, the endpoint-suffix
quantile, normalization of the exterior seam, the strengthened rational tax,
its gcd weakening, and integer rounding.  It also checks the complementary
endpoint/flood split: the inward and outward pieces of one endpoint-owner tooth
cancel the residual numerator `Q` exactly when their two independent budgets
are combined.  Finally it records the small-branch inverse used by the nested
multi-owner endpoint forest.

Selection of an endpoint owner, proper crossing of the two actual teeth,
containment of the five-comb survivor beyond the endpoint-owner wall, and
extraction of the terminal chronological occurrence remain explicit paper
providers.  No covering statement is hidden in an axiom here.
-/

namespace LRC14.TerminalEndpointTransferGcdTax

/-- A strict signed endpoint residual `-c < z < c` gives a positive endpoint
numerator `Q=c+z` strictly below `2c`. -/
theorem endpoint_residual_bounds (c z : ℚ)
    (hzLower : -c < z) (hzUpper : z < c) :
    0 < c + z ∧ c + z < 2 * c := by
  constructor <;> linarith

/-- The exact endpoint numerator formula makes the outward part of the owner
tooth shorter than one full owner-tooth width. -/
theorem endpoint_width_lt_owner_tooth (c x Q q : ℚ)
    (hc : 0 < c) (hx : 0 < x) (hQUpper : Q < 2 * c)
    (hq : q = Q / (14 * c * x)) :
    q < 1 / (7 * x) := by
  rw [hq]
  have hden : 0 < 14 * c * x := by positivity
  calc
    Q / (14 * c * x) < (2 * c) / (14 * c * x) :=
      div_lt_div_of_pos_right hQUpper hden
    _ = 1 / (7 * x) := by
      field_simp
      ring

/-- Since every endpoint owner is faster than the carrier, its outward wall
is reached before the far wall of the adjacent carrier tooth. -/
theorem owner_tooth_inside_carrier_tooth (c x : ℚ)
    (hc : 0 < c) (hcx : c < x) :
    1 / (7 * x) < 1 / (7 * c) := by
  have hx : 0 < x := lt_trans hc hcx
  apply (div_lt_div_iff₀ (by positivity : 0 < 7 * x)
    (by positivity : 0 < 7 * c)).2
  nlinarith

/-- A protected margin of size `1/(2c)` contains every endpoint seam, whose
supremum is below `1/(7x)`. -/
theorem endpoint_seam_inside_half_carrier_margin (c x q : ℚ)
    (hc : 0 < c) (hcx : c < x) (hq : q < 1 / (7 * x)) :
    q < 1 / (2 * c) := by
  have hOwner := owner_tooth_inside_carrier_tooth c x hc hcx
  have hCarrier : 1 / (7 * c) < 1 / (2 * c) := by
    have hc7 : 0 < 7 * c := by positivity
    have hc2 : 0 < 2 * c := by positivity
    apply (div_lt_div_iff₀ hc7 hc2).2
    nlinarith
  exact lt_trans hq (lt_trans hOwner hCarrier)

/-- The endpoint density inverse remains valid without a prior small-tail
case split.  If the outer suffix is longer than `1/6` the conclusion is
automatic; otherwise its mass is at most `(3/4)y`. -/
theorem outer_suffix_quantile (y mass : ℚ)
    (hmass : 11 / 360 < mass)
    (hSmall : y ≤ 1 / 6 → mass ≤ (3 / 4) * y) :
    11 / 270 < y := by
  by_cases hy : y ≤ 1 / 6
  · have hUpper := hSmall hy
    linarith
  · have hyLarge : 1 / 6 < y := lt_of_not_ge hy
    linarith

/-- Normalizing a physical endpoint seam by the complete `a`-safe component
length `6/(7a)` gives the exact tax `aQ/(12cx)`. -/
theorem normalized_endpoint_seam (c a x Q : ℚ)
    (hc : c ≠ 0) (hx : x ≠ 0) :
    (7 * a / 6) * (Q / (14 * c * x)) = a * Q / (12 * c * x) := by
  field_simp
  ring

/-- Algebraic identity converting the strict normalized tail margin into the
integer-coefficient endpoint tax. -/
theorem endpoint_tax_margin_identity (c a x Q : ℚ)
    (hc : c ≠ 0) (hx : x ≠ 0) :
    540 * c * x *
        ((13 / 12 - a / (2 * c)) -
          (11 / 270 + a * Q / (12 * c * x))) =
      563 * c * x - 270 * a * x - 45 * a * Q := by
  field_simp
  ring

/-- The survivor suffix inequality `ell-eta>11/270`, the centered protrusion
formula, and `rho≤1/2` imply the exact residue tax. -/
theorem endpoint_tax_consumer
    (c a x Q ell eta rho : ℚ)
    (hc : 0 < c) (hx : 0 < x)
    (hEll : ell = 1 / 2 + 7 * rho / 6 - a / (2 * c))
    (hEta : eta = a * Q / (12 * c * x))
    (hSuffix : 11 / 270 < ell - eta)
    (hRho : rho ≤ 1 / 2) :
    270 * a * x + 45 * a * Q < 563 * c * x := by
  have hEllUpper : ell ≤ 13 / 12 - a / (2 * c) := by
    rw [hEll]
    linarith
  have hTail : 11 / 270 + a * Q / (12 * c * x) < ell := by
    rw [hEta] at hSuffix
    linarith
  have hMargin :
      0 < (13 / 12 - a / (2 * c)) -
        (11 / 270 + a * Q / (12 * c * x)) := by
    linarith
  have hScale : 0 < 540 * c * x := by positivity
  have hPositive := mul_pos hScale hMargin
  rw [endpoint_tax_margin_identity c a x Q (ne_of_gt hc) (ne_of_gt hx)] at hPositive
  linarith

/-- Retaining the exact centered error numerator `e=2c*rho` sharpens the
coarse `563c` tax to `(248c+315e)`. -/
theorem exact_centered_error_tax_consumer
    (c a x Q e ell eta : ℚ)
    (hc : 0 < c) (hx : 0 < x)
    (hEll : ell = (6 * c + 7 * e - 6 * a) / (12 * c))
    (hEta : eta = a * Q / (12 * c * x))
    (hSuffix : 11 / 270 < ell - eta) :
    270 * a * x + 45 * a * Q < (248 * c + 315 * e) * x := by
  rw [hEll, hEta] at hSuffix
  have hScale : 0 < 540 * c * x := by positivity
  have h := mul_lt_mul_of_pos_left hSuffix hScale
  field_simp at h
  nlinarith

/-- Replacing the exact endpoint numerator by any nonnegative lower bound
(in particular `gcd(c,x)`) weakens the exact tax safely. -/
theorem endpoint_gcd_tax_consumer
    (c a x Q g : ℚ)
    (ha : 0 ≤ a) (hg : g ≤ Q)
    (hExact : 270 * a * x + 45 * a * Q < 563 * c * x) :
    270 * a * x + 45 * a * g < 563 * c * x := by
  have hMonotone : 45 * a * g ≤ 45 * a * Q := by
    exact mul_le_mul_of_nonneg_left hg (by positivity)
  linarith

/-- Strict rational endpoint tax rounds to the displayed closed integer cut. -/
theorem integer_endpoint_tax_rounding (c a x Q : ℤ)
    (hTax : 270 * a * x + 45 * a * Q < 563 * c * x) :
    270 * a * x + 45 * a * Q ≤ 563 * c * x - 1 := by
  omega

/-! ## Complementary endpoint/flood cancellation -/

/-- In the functional normalization, the inward part of an endpoint-owner
tooth and the appropriately rescaled outward endpoint seam sum to one complete
`c/(8h)` fastest-tooth quantum.  The endpoint numerator `Q` disappears. -/
theorem functional_endpoint_partition_identity (c a h Q : ℚ)
    (hc : c ≠ 0) (ha : a ≠ 0) (hh : h ≠ 0) :
    (2 * c - Q) / (16 * h) +
        (3 * c / (4 * a)) * (a * Q / (12 * c * h)) =
      c / (8 * h) := by
  field_simp
  ring

/-- In the harmonic normalization, the same inward/outward partition is one
complete `7/(6h)` fastest-tooth quantum, again independently of `Q`. -/
theorem harmonic_endpoint_partition_identity (c a h Q : ℚ)
    (hc : c ≠ 0) (ha : a ≠ 0) (hh : h ≠ 0) :
    7 * (2 * c - Q) / (12 * c * h) +
        (7 / a) * (a * Q / (12 * c * h)) =
      7 / (6 * h) := by
  field_simp
  ring

/-- Abstract functional two-budget consumer.  The first hypothesis is the
internal flood invoice (chronological seams plus the inward tooth piece); the
second is the strict exterior survivor invoice.  Their supports need not be
identified: once the paper proves the two inequalities independently, the
complementary partition upgrades them to a full fastest-tooth quantum. -/
theorem functional_endpoint_flood_consumer
    (c a h Q seam internal tail : ℚ)
    (hc : 0 < c) (ha : 0 < a) (hh : 0 < h)
    (hInternal : (2 * c - Q) / (16 * h) + seam ≤ internal)
    (hExterior : a * Q / (12 * c * h) < tail) :
    seam + c / (8 * h) <
      internal + (3 * c / (4 * a)) * tail := by
  have hScale : 0 < 3 * c / (4 * a) := by positivity
  have hExteriorScaled := mul_lt_mul_of_pos_left hExterior hScale
  have hPartition := functional_endpoint_partition_identity c a h Q
    (ne_of_gt hc) (ne_of_gt ha) (ne_of_gt hh)
  linarith

/-- Harmonic counterpart of `functional_endpoint_flood_consumer`. -/
theorem harmonic_endpoint_flood_consumer
    (c a h Q seam internal tail : ℚ)
    (hc : 0 < c) (ha : 0 < a) (hh : 0 < h)
    (hInternal : 7 * (2 * c - Q) / (12 * c * h) + seam ≤ internal)
    (hExterior : a * Q / (12 * c * h) < tail) :
    seam + 7 / (6 * h) < internal + (7 / a) * tail := by
  have hScale : 0 < 7 / a := by positivity
  have hExteriorScaled := mul_lt_mul_of_pos_left hExterior hScale
  have hPartition := harmonic_endpoint_partition_identity c a h Q
    (ne_of_gt hc) (ne_of_gt ha) (ne_of_gt hh)
  linarith

/-! ## Nested endpoint-owner forest: small density branch -/

/-- If the star-forest overlap credit contributes `(3/4)S` to the lower mass
bound while the surviving suffix has density at most `3/4`, inversion gives
the exact additive endpoint tax `11/270 + S`. -/
theorem multiowner_small_branch_inverse (S y mass : ℚ)
    (hLower : 11 / 360 + (3 / 4) * S < mass)
    (hUpper : mass ≤ (3 / 4) * y) :
    11 / 270 + S < y := by
  linarith

/-- Geometric form of `multiowner_small_branch_inverse`: if `etaMax` is the
longest nested endpoint seam and `S` is the sum of the remaining seams, the
small-density branch charges their full sum against the carrier tail. -/
theorem multiowner_small_branch_tail (ell etaMax S mass : ℚ)
    (hLower : 11 / 360 + (3 / 4) * S < mass)
    (hUpper : mass ≤ (3 / 4) * (ell - etaMax)) :
    11 / 270 + S + etaMax < ell := by
  have hInverse := multiowner_small_branch_inverse S (ell - etaMax) mass
    hLower hUpper
  linarith

#print axioms endpoint_residual_bounds
#print axioms endpoint_width_lt_owner_tooth
#print axioms owner_tooth_inside_carrier_tooth
#print axioms endpoint_seam_inside_half_carrier_margin
#print axioms outer_suffix_quantile
#print axioms normalized_endpoint_seam
#print axioms endpoint_tax_margin_identity
#print axioms endpoint_tax_consumer
#print axioms exact_centered_error_tax_consumer
#print axioms endpoint_gcd_tax_consumer
#print axioms integer_endpoint_tax_rounding
#print axioms functional_endpoint_partition_identity
#print axioms harmonic_endpoint_partition_identity
#print axioms functional_endpoint_flood_consumer
#print axioms harmonic_endpoint_flood_consumer
#print axioms multiowner_small_branch_inverse
#print axioms multiowner_small_branch_tail

end LRC14.TerminalEndpointTransferGcdTax
