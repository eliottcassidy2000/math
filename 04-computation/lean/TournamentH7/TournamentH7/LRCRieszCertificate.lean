/-
  TournamentH7.LRCRieszCertificate — soundness of the Riesz-product loneliness certificate
  (opus-2026-07-09-S173).

  The LRC singular-series / lonely-measure route (THM-515, HYP-2540): a speed set `S` is LOOSE
  (its lonely set has positive measure) iff the covering multiplicity `M(τ) = #{v : ‖vτ‖ ≤ 1/14}`
  is NOT `≥ 1` almost everywhere.  The Riesz-product method certifies looseness by pairing `M`
  against a NONNEGATIVE test density `R` (a Riesz product `∏(1 + aₘcos2πmτ) ≥ 0`):

      if  `R ≥ 0`  and  `∫ M·R < ∫ R`,  then  `M < 1` on a positive-measure set  ⟹  `S` is loose.

  This file proves that soundness in full generality (kernel-pure).  It is the LOGICAL CORE of the
  singular-series lower-bound program — the POSITIVE-DEFINITE route (`ĥ = 1_safe ≥ 0`) that sidesteps
  the signed/Mertens cancellation blocking the covering-`W` side (opus-S172).  The analytic content —
  constructing an `R` with `∫M·R < ∫R` for every loose `S` (the tuned dissociated Bedert-2025
  construction, `inf L > 0`) — is exactly what stays open; here the certificate's VALIDITY is checked.
-/
import Mathlib

namespace LonelyRunner.Riesz

open MeasureTheory

/-- **Riesz-certificate soundness.**  For a nonnegative test density `R` with `M·R` and `R`
integrable, the strict pairing inequality `∫ M·R < ∫ R` rules out `M ≥ 1` almost everywhere —
i.e. `M < 1` on a positive-measure set.  (In the LRC application `M ≥ 0` is the covering
multiplicity and `{M < 1} = {M = 0}` is the lonely set, so this certifies positive lonely
measure ⟹ loose.)  Pure monotonicity of the integral: if `1 ≤ M` a.e. then `R ≤ M·R` a.e.
(as `R ≥ 0`), forcing `∫R ≤ ∫M·R`, contradicting the certificate. -/
theorem riesz_certificate {α : Type*} [MeasurableSpace α] {μ : Measure α} {M R : α → ℝ}
    (hR : 0 ≤ᵐ[μ] R)
    (hMR : Integrable (fun x => M x * R x) μ) (hRi : Integrable R μ)
    (hcert : ∫ x, M x * R x ∂μ < ∫ x, R x ∂μ) :
    ¬ (∀ᵐ x ∂μ, (1 : ℝ) ≤ M x) := by
  intro hM
  have hmono : ∫ x, R x ∂μ ≤ ∫ x, M x * R x ∂μ := by
    apply integral_mono_ae hRi hMR
    filter_upwards [hR, hM] with x hRx hMx
    simp only [Pi.zero_apply] at hRx
    calc R x = 1 * R x := (one_mul _).symm
      _ ≤ M x * R x := mul_le_mul_of_nonneg_right hMx hRx
  linarith

/-- Contrapositive packaging: if the covering multiplicity `M` IS `≥ 1` a.e. (`S` tight), then no
nonnegative `R` can witness `∫ M·R < ∫ R` — the VALIDITY guarantee (no false positive), matching the
numeric check that the tight `{1..12}∪{182}` gives ratio `≥ 1` (opus-S173). -/
theorem no_certificate_of_ae_covered {α : Type*} [MeasurableSpace α] {μ : Measure α} {M R : α → ℝ}
    (hR : 0 ≤ᵐ[μ] R)
    (hMR : Integrable (fun x => M x * R x) μ) (hRi : Integrable R μ)
    (hcov : ∀ᵐ x ∂μ, (1 : ℝ) ≤ M x) :
    ∫ x, R x ∂μ ≤ ∫ x, M x * R x ∂μ := by
  apply integral_mono_ae hRi hMR
  filter_upwards [hR, hcov] with x hRx hMx
  simp only [Pi.zero_apply] at hRx
  calc R x = 1 * R x := (one_mul _).symm
    _ ≤ M x * R x := mul_le_mul_of_nonneg_right hMx hRx

end LonelyRunner.Riesz
