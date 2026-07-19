/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Codex (LRC multi-agent project, 2026-07-19)
-/
import Mathlib.Tactic

/-!
# Coherent blocker selection on one chronological tooth word (THM-1254)

The paper layer chooses every centered-spoke blocker from one deletion-minimal
tooth cover.  A maximum-speed vertex on a blocker cycle has a speed-descent
outgoing edge, whose THM-1248 relative digit is binary.  The order of that
edge's two consecutive marked teeth selects either the original word or its
reflection, making one full chronological-drift invoice unconditional.  This
module checks that same-edge dichotomy, the address-product collapse, the
mixed-circuit invoice, reflection, and the secondary absolute carrier cutoff.
-/

namespace LRC14
namespace CoherentBlockerChronology

/-- When the centered closing blocker tooth is the initial chronological
tooth, the generic address-product gap collapses to one relative digit. -/
theorem coherent_address_product_collapse
    (P n₀ nᵣ k delta : ℤ)
    (hP : P = k + nᵣ + delta) :
    P * n₀ - n₀ * nᵣ = n₀ * (k + delta) := by
  rw [hP]
  ring

/-- Substituting the coherent product gap into the mixed chronological /
centered identity removes the free affine sign term. -/
theorem coherent_mixed_residual_identity
    (P N s₀ sᵣ n₀ nᵣ k delta residual Delta : ℚ)
    (hn₀ : n₀ ≠ 0) (hnᵣ : nᵣ ≠ 0)
    (hN : N = n₀)
    (hP : P = k + nᵣ + delta)
    (hDelta : Delta = s₀ / n₀ - sᵣ / nᵣ)
    (hresidual : residual = P * s₀ - N * sᵣ) :
    residual = n₀ * nᵣ * Delta + s₀ * (k + delta) := by
  rw [hresidual, hN, hP, hDelta]
  field_simp [hn₀, hnᵣ]
  ring

/-- A nonnegative coherent digit factor makes the unshifted residual pay the
entire positive chronological drift. -/
theorem coherent_nonnegative_digit_invoice
    {residual n₀ nᵣ Delta s₀ k delta : ℚ}
    (hidentity : residual = n₀ * nᵣ * Delta + s₀ * (k + delta))
    (hn₀ : 0 < n₀) (hnᵣ : 0 < nᵣ) (hDelta : 0 < Delta)
    (hs₀ : 0 ≤ s₀) (hgap : 0 ≤ k + delta) :
    0 < n₀ * nᵣ * Delta ∧ n₀ * nᵣ * Delta ≤ residual := by
  constructor
  · positivity
  · rw [hidentity]
    exact le_add_of_nonneg_right (mul_nonneg hs₀ hgap)

/-- Reflection sends the lower carrier-gap address `k` to `c-k-1` and a
target tooth address `N` to `d-N`. -/
theorem reflected_middle_address (c d k N : ℤ) :
    (c - k - 1) + (d - N) = (c + d) - 1 - (k + N) := by
  ring

/-- The THM-1248 relative digit transforms by `delta ↦ 1-delta` under circle
reflection. -/
theorem reflected_relative_digit
    (Pj Qj M delta Pj' M' delta' : ℤ)
    (hdelta : delta = Pj - M)
    (hPj' : Pj' = Qj - Pj)
    (hM' : M' = Qj - 1 - M)
    (hdelta' : delta' = Pj' - M') :
    delta' = 1 - delta := by
  rw [hdelta', hPj', hM', hdelta]
  ring

/-- The reflected coherent sign factor is `c-k-delta`. -/
theorem reflected_gap_digit (c k delta : ℤ) :
    (c - k - 1) + (1 - delta) = c - k - delta := by
  ring

/-- On a THM-1248 binary speed-descent edge, both possible same-edge invoice
factors are nonnegative.  The original word uses `k+delta`; the reflected word
uses `c-k-delta`. -/
theorem binary_speed_descent_factors
    (c k delta : ℤ)
    (hk : 0 ≤ k) (hkc : k < c)
    (hbinary : delta = 0 ∨ delta = 1) :
    0 ≤ k + delta ∧ 0 ≤ c - k - delta := by
  rcases hbinary with rfl | rfl <;> omega

/-- The marked teeth are distinct.  If the successor tooth is earlier, the
original binary factor pays; if it is later, reflection reverses the order and
the reflected binary factor pays.  The same cycle edge is used in both cases.
-/
theorem binary_speed_descent_same_edge_factor
    (c k delta a aNext : ℤ)
    (hk : 0 ≤ k) (hkc : k < c)
    (hbinary : delta = 0 ∨ delta = 1)
    (hne : aNext ≠ a) :
    (aNext < a ∧ 0 ≤ k + delta) ∨
      (a < aNext ∧ 0 ≤ c - k - delta) := by
  have hfactors := binary_speed_descent_factors c k delta hk hkc hbinary
  omega

/-- Full arithmetic consumer of the unconditional same-edge dichotomy.  The
paper topology supplies a binary speed-descent edge and distinct marked
positions; according to their order, either the original or reflected mixed
residual pays its entire strictly positive chronological drift. -/
theorem binary_speed_descent_same_edge_full_invoice
    {c k delta residual residualR n₀ nᵣ Delta s₀
      n₀R nᵣR DeltaR s₀R : ℚ}
    {a aNext : ℤ}
    (hk : 0 ≤ k) (hkc : k < c)
    (hbinary : delta = 0 ∨ delta = 1)
    (hne : aNext ≠ a)
    (hidentity : residual = n₀ * nᵣ * Delta + s₀ * (k + delta))
    (hidentityR :
      residualR = n₀R * nᵣR * DeltaR + s₀R * (c - k - delta))
    (hn₀ : 0 < n₀) (hnᵣ : 0 < nᵣ) (hDelta : 0 < Delta)
    (hs₀ : 0 ≤ s₀)
    (hn₀R : 0 < n₀R) (hnᵣR : 0 < nᵣR) (hDeltaR : 0 < DeltaR)
    (hs₀R : 0 ≤ s₀R) :
    (aNext < a ∧ 0 < n₀ * nᵣ * Delta ∧
      n₀ * nᵣ * Delta ≤ residual) ∨
    (a < aNext ∧ 0 < n₀R * nᵣR * DeltaR ∧
      n₀R * nᵣR * DeltaR ≤ residualR) := by
  have hpos : aNext < a ∨ a < aNext := by omega
  have hfactors : 0 ≤ k + delta ∧ 0 ≤ c - k - delta := by
    rcases hbinary with rfl | rfl <;> constructor <;> linarith
  rcases hpos with hdown | hup
  · left
    refine ⟨hdown, ?_⟩
    exact coherent_nonnegative_digit_invoice hidentity hn₀ hnᵣ hDelta hs₀
      hfactors.1
  · right
    refine ⟨hup, ?_⟩
    constructor
    · positivity
    · rw [hidentityR]
      exact le_add_of_nonneg_right (mul_nonneg hs₀R hfactors.2)

/-- Secondary general-digit guardrail: failure of the original drift invoice
can occur only in one of the first 586 carrier gaps. -/
theorem lower_gap_failure_bounded
    (k delta : ℤ)
    (hk : 0 ≤ k) (hdelta : -586 ≤ delta)
    (hfail : k + delta < 0) :
    k ≤ 585 := by
  omega

/-- Secondary general-digit guardrail: if both an arbitrary original cyclic
descent and an arbitrary reflected cyclic descent fail their drift invoices,
the carrier itself is absolutely bounded. -/
theorem two_sided_invoice_failure_bounds_carrier
    (c k deltaDown deltaUp : ℤ)
    (hk : 0 ≤ k)
    (hdeltaDown : -586 ≤ deltaDown)
    (hdeltaUp : deltaUp ≤ 587)
    (hdown : k + deltaDown < 0)
    (hup : c - k - deltaUp < 0) :
    c ≤ 1171 := by
  omega

/-- Secondary sign-only scale split for arbitrary (not necessarily binary)
edge choices: one circuit has a nonnegative factor, or `c≤1171`. -/
theorem coherent_invoice_or_bounded_carrier
    (c k deltaDown deltaUp : ℤ)
    (hk : 0 ≤ k)
    (hdeltaDown : -586 ≤ deltaDown)
    (hdeltaUp : deltaUp ≤ 587) :
    0 ≤ k + deltaDown ∨ 0 ≤ c - k - deltaUp ∨ c ≤ 1171 := by
  omega

#print axioms coherent_address_product_collapse
#print axioms coherent_mixed_residual_identity
#print axioms coherent_nonnegative_digit_invoice
#print axioms reflected_middle_address
#print axioms reflected_relative_digit
#print axioms reflected_gap_digit
#print axioms binary_speed_descent_factors
#print axioms binary_speed_descent_same_edge_factor
#print axioms binary_speed_descent_same_edge_full_invoice
#print axioms lower_gap_failure_bounded
#print axioms two_sided_invoice_failure_bounds_carrier
#print axioms coherent_invoice_or_bounded_carrier

end CoherentBlockerChronology
end LRC14
