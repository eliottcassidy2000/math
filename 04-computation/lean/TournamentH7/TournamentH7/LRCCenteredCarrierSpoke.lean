/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Codex (LRC multi-agent project, 2026-07-19)
-/
import Mathlib.Tactic

/-!
# Arithmetic consumers for THM-1240 centered carrier-spoke blockers

The nearest-integer theorem supplies a beat within half a numerator of the
scaled slow-gap centre.  This module checks the resulting interior/depth and
curvature inequalities, the nonzero-clock divisibility core, orbit collision,
and the cut-clock separation identity.
-/

namespace LRC14
namespace CenteredCarrierSpoke

/-- Since `q=c+d>2c`, the half-numerator approximation lies strictly inside
the slow gap's half-width `3/(7c)`. -/
theorem centered_beat_is_interior_scale
    {c d q : ℝ} (hc : 0 < c) (hcd : c < d) (hq : q = c + d) :
    1 / (2 * q) < 3 / (7 * c) := by
  have hqpos : 0 < q := by rw [hq]; linarith
  apply (div_lt_div_iff₀ (by positivity : (0 : ℝ) < 2 * q)
    (by positivity : (0 : ℝ) < 7 * c)).2
  rw [hq]
  nlinarith

/-- Exact deep-distance invoice from the centred nearest-integer error. -/
theorem centered_depth_invoice
    {c d q delta : ℝ} (hc : 0 < c) (hcd : c < d)
    (hq : q = c + d) (hdelta : |delta| ≤ c / (2 * q)) :
    d / (2 * q) ≤ 1 / 2 - |delta| ∧
      1 / 4 < d / (2 * q) := by
  have hqpos : 0 < q := by rw [hq]; linarith
  have hsum : c / (2 * q) + d / (2 * q) = 1 / 2 := by
    rw [hq]
    field_simp
    ring
  constructor
  · linarith
  · apply (div_lt_div_iff₀ (by norm_num : (0 : ℝ) < 4)
      (by positivity : (0 : ℝ) < 2 * q)).2
    rw [hq]
    nlinarith

/-- A determinant at least `d/2` gives the explicit positive curvature
slack `14D-q ≥ 6d-c`. -/
theorem deep_determinant_slack
    {c d q D : ℝ} (hcd : c < d) (hq : q = c + d)
    (hD : d / 2 ≤ D) :
    6 * d - c ≤ 14 * D - q ∧ 0 < 14 * D - q := by
  constructor <;> rw [hq] <;> nlinarith

/-- If the master residue were zero, every speed divisible by the common
scale would be integral at that beat. -/
theorem zero_master_residue_forces_integral
    {q d₀ L c c' p r : ℕ}
    (hq : q = d₀ * L) (hc : c = d₀ * c') (hp : p = L * r) :
    q ∣ c * p := by
  refine ⟨c' * r, ?_⟩
  omega

/-- Seven iterates of a self-map on six labels have an orbit collision, the
finite core behind the blocker-cycle argument. -/
theorem six_label_orbit_collision (f : Fin 6 → Fin 6) (start : Fin 6) :
    ∃ i j : Fin 7, i ≠ j ∧
      (f^[i.val]) start = (f^[j.val]) start := by
  exact Fintype.exists_ne_map_eq_of_card_lt
    (fun i : Fin 7 ↦ (f^[i.val]) start) (by norm_num)

/-- A loopless blocker selection excludes one-cycles. -/
theorem loopless_excludes_fixed_point (f : Fin 6 → Fin 6)
    (hloopless : ∀ i, f i ≠ i) :
    ∀ i, ¬f i = i := by
  intro i hi
  exact hloopless i hi

/-- Adjacent carrier-spoke clocks inherit exactly the macroscopic speed cut. -/
theorem spoke_clock_cut_identity (c x y : ℝ) :
    (c + y) - (c + x) = y - x := by
  ring

/-- Exact lower endpoint for the spoke slack. -/
theorem spoke_slack_is_positive {c d : ℝ} (hcd : c < d) :
    0 < 6 * d - c := by
  nlinarith

#print axioms centered_beat_is_interior_scale
#print axioms centered_depth_invoice
#print axioms deep_determinant_slack
#print axioms zero_master_residue_forces_integral
#print axioms six_label_orbit_collision
#print axioms loopless_excludes_fixed_point
#print axioms spoke_clock_cut_identity
#print axioms spoke_slack_is_positive

end CenteredCarrierSpoke
end LRC14
