/-
  TournamentH7.LRCWitnessPartA -- finite-Vmax error-budget glue for the LRC(14)
  witness route.

  KPS S30 observed that THM-527 Part A's finite-Vmax correction has the form

      rho_K = witnessG2 + O(#arcs(GOOD(E)) / Vmax).

  Once the Bonferroni/p0 route supplies a positive margin `delta <= witnessG2`,
  Part A only needs the finite-ruler error to be smaller than that margin.  This
  module records that arithmetic as sorry-free Lean glue.  It deliberately does
  not define `rho_K`, `GOOD(E)`, or the arc count; those are the remaining
  analytic/event-level objects.
-/

import Mathlib.Tactic
import TournamentH7.LRCWitnessBonferroni

namespace LonelyRunner
namespace LRC14
namespace PartA

/-- If a finite witness density `rhoK` is within `eps` of an asymptotic density
`rho`, and `rho` is at least a margin `delta > eps`, then `rhoK` is positive. -/
theorem finite_witness_pos_from_abs_error
    (rhoK rho delta eps : ℝ)
    (hfloor : delta ≤ rho)
    (herror : |rhoK - rho| ≤ eps)
    (hbudget : eps < delta) :
    0 < rhoK := by
  have hlow : -eps ≤ rhoK - rho := (abs_le.mp herror).1
  linarith

/-- The same positivity criterion with the THM-527 Part-A error budget
`eps = #arcs / Vmax`. -/
theorem finite_witness_pos_from_arc_error
    (rhoK rho delta : ℝ) (arcCount vmax : ℕ)
    (hfloor : delta ≤ rho)
    (herror : |rhoK - rho| ≤ (arcCount : ℝ) / (vmax : ℝ))
    (hbudget : (arcCount : ℝ) / (vmax : ℝ) < delta) :
    0 < rhoK :=
  finite_witness_pos_from_abs_error rhoK rho delta
    ((arcCount : ℝ) / (vmax : ℝ)) hfloor herror hbudget

/-- Multiplicative form of the arc budget: `#arcs < delta * Vmax` implies
`#arcs / Vmax < delta`. -/
theorem arc_div_lt_delta_of_lt_mul
    (arcCount vmax : ℕ) (delta : ℝ) (hvmax : 0 < vmax)
    (hbudget : (arcCount : ℝ) < delta * (vmax : ℝ)) :
    (arcCount : ℝ) / (vmax : ℝ) < delta := by
  have hvmaxR : (0 : ℝ) < (vmax : ℝ) := by exact_mod_cast hvmax
  have hvmax_ne : (vmax : ℝ) ≠ 0 := ne_of_gt hvmaxR
  have hcancel :
      ((arcCount : ℝ) / (vmax : ℝ)) * (vmax : ℝ) = (arcCount : ℝ) := by
    field_simp [hvmax_ne]
  nlinarith

/-- Multiplicative version of `finite_witness_pos_from_arc_error`. -/
theorem finite_witness_pos_from_arc_error_mul
    (rhoK rho delta : ℝ) (arcCount vmax : ℕ)
    (hvmax : 0 < vmax)
    (hfloor : delta ≤ rho)
    (herror : |rhoK - rho| ≤ (arcCount : ℝ) / (vmax : ℝ))
    (hbudget : (arcCount : ℝ) < delta * (vmax : ℝ)) :
    0 < rhoK :=
  finite_witness_pos_from_arc_error rhoK rho delta arcCount vmax hfloor herror
    (arc_div_lt_delta_of_lt_mul arcCount vmax delta hvmax hbudget)

/-- Uniform-bound form of the arc budget.  If `arcCount ≤ arcBound`, the margin
`delta` is at least a positive floor `deltaFloor`, and the finite ruler satisfies
`arcBound < deltaFloor * Vmax`, then the actual error term is below `delta`. -/
theorem arc_div_lt_delta_of_le_bound
    (arcCount arcBound vmax : ℕ) (deltaFloor delta : ℝ) (hvmax : 0 < vmax)
    (harc : arcCount ≤ arcBound)
    (hdelta : deltaFloor ≤ delta)
    (hbudget : (arcBound : ℝ) < deltaFloor * (vmax : ℝ)) :
    (arcCount : ℝ) / (vmax : ℝ) < delta := by
  have harcR : (arcCount : ℝ) ≤ (arcBound : ℝ) := by exact_mod_cast harc
  have hbudget' : (arcCount : ℝ) < deltaFloor * (vmax : ℝ) :=
    lt_of_le_of_lt harcR hbudget
  have hdiv : (arcCount : ℝ) / (vmax : ℝ) < deltaFloor :=
    arc_div_lt_delta_of_lt_mul arcCount vmax deltaFloor hvmax hbudget'
  exact lt_of_lt_of_le hdiv hdelta

/-- Finite witness positivity from a uniform arc bound and a uniform margin
floor.  This is the algebraic shape suggested by the period-bounded
`#arcs(GOOD)` signal: once `Vmax` clears the uniform threshold, the finite
correction cannot erase the asymptotic witness density. -/
theorem finite_witness_pos_from_arc_bound
    (rhoK rho delta deltaFloor : ℝ) (arcCount arcBound vmax : ℕ)
    (hvmax : 0 < vmax)
    (hfloor : delta ≤ rho)
    (hdelta : deltaFloor ≤ delta)
    (harc : arcCount ≤ arcBound)
    (herror : |rhoK - rho| ≤ (arcCount : ℝ) / (vmax : ℝ))
    (hbudget : (arcBound : ℝ) < deltaFloor * (vmax : ℝ)) :
    0 < rhoK :=
  finite_witness_pos_from_arc_error rhoK rho delta arcCount vmax hfloor herror
    (arc_div_lt_delta_of_le_bound arcCount arcBound vmax deltaFloor delta
      hvmax harc hdelta hbudget)

/-- The p0-wide-bound margin itself lower-bounds `witnessG2`, shape by shape.
This is the margin version of the HYP-2832 unification, before comparing that
margin with the conservative global `m_P` floor. -/
theorem p0_margin_le_witnessG2_shapes
    (nuShape measGP p0Shape cap delta : Shape → ℝ)
    (hbonf : ∀ s, nuShape s + measGP s - 1 ≤ witnessG2 s)
    (hDp0  : ∀ s, (1 - nuShape s) ≤ p0Shape s)
    (hp0cap : ∀ s, p0Shape s ≤ cap s - delta s)
    (hmeasGP : ∀ s, cap s ≤ measGP s)
    (s : Shape) :
    delta s ≤ witnessG2 s := by
  have h1 : nuShape s + measGP s - 1 ≤ witnessG2 s := hbonf s
  have h2 : (1 - nuShape s) ≤ p0Shape s := hDp0 s
  have h3 : p0Shape s ≤ cap s - delta s := hp0cap s
  have h4 : cap s ≤ measGP s := hmeasGP s
  linarith

/-- Shape-indexed finite witness positivity from an already supplied
`delta <= witnessG2` floor and an arc-count error budget. -/
theorem finite_witness_pos_from_floor_shapes
    (finiteRho delta : Shape → ℝ) (arcCount vmax : Shape → ℕ)
    (hfloor : ∀ s, delta s ≤ witnessG2 s)
    (herror : ∀ s, |finiteRho s - witnessG2 s| ≤
      (arcCount s : ℝ) / (vmax s : ℝ))
    (hbudget : ∀ s, (arcCount s : ℝ) / (vmax s : ℝ) < delta s)
    (s : Shape) :
    0 < finiteRho s :=
  finite_witness_pos_from_arc_error (finiteRho s) (witnessG2 s) (delta s)
    (arcCount s) (vmax s) (hfloor s) (herror s) (hbudget s)

/-- Shape-indexed finite witness positivity from a uniform arc bound and a
uniform positive margin floor. -/
theorem finite_witness_pos_from_uniform_arc_bound_shapes
    (finiteRho delta : Shape → ℝ) (arcCount vmax : Shape → ℕ)
    (arcBound : ℕ) (deltaFloor : ℝ)
    (hfloor : ∀ s, delta s ≤ witnessG2 s)
    (hdelta : ∀ s, deltaFloor ≤ delta s)
    (hvmax : ∀ s, 0 < vmax s)
    (harc : ∀ s, arcCount s ≤ arcBound)
    (herror : ∀ s, |finiteRho s - witnessG2 s| ≤
      (arcCount s : ℝ) / (vmax s : ℝ))
    (hbudget : ∀ s, (arcBound : ℝ) < deltaFloor * (vmax s : ℝ))
    (s : Shape) :
    0 < finiteRho s :=
  finite_witness_pos_from_arc_bound (finiteRho s) (witnessG2 s) (delta s)
    deltaFloor (arcCount s) arcBound (vmax s) (hvmax s) (hfloor s)
    (hdelta s) (harc s) (herror s) (hbudget s)

/-- Shape-indexed finite witness positivity directly from the p0-wide-bound
margin, `D <= p0`, and the arc-count error budget. -/
theorem finite_witness_pos_from_p0_margin_shapes
    (nuShape measGP p0Shape cap delta finiteRho : Shape → ℝ)
    (arcCount vmax : Shape → ℕ)
    (hbonf : ∀ s, nuShape s + measGP s - 1 ≤ witnessG2 s)
    (hDp0  : ∀ s, (1 - nuShape s) ≤ p0Shape s)
    (hp0cap : ∀ s, p0Shape s ≤ cap s - delta s)
    (hmeasGP : ∀ s, cap s ≤ measGP s)
    (herror : ∀ s, |finiteRho s - witnessG2 s| ≤
      (arcCount s : ℝ) / (vmax s : ℝ))
    (hbudget : ∀ s, (arcCount s : ℝ) / (vmax s : ℝ) < delta s)
    (s : Shape) :
    0 < finiteRho s :=
  finite_witness_pos_from_floor_shapes finiteRho delta arcCount vmax
    (p0_margin_le_witnessG2_shapes nuShape measGP p0Shape cap delta
      hbonf hDp0 hp0cap hmeasGP)
    herror hbudget s

/-- Multiplicative-budget version of
`finite_witness_pos_from_p0_margin_shapes`. -/
theorem finite_witness_pos_from_p0_margin_shapes_mul
    (nuShape measGP p0Shape cap delta finiteRho : Shape → ℝ)
    (arcCount vmax : Shape → ℕ)
    (hbonf : ∀ s, nuShape s + measGP s - 1 ≤ witnessG2 s)
    (hDp0  : ∀ s, (1 - nuShape s) ≤ p0Shape s)
    (hp0cap : ∀ s, p0Shape s ≤ cap s - delta s)
    (hmeasGP : ∀ s, cap s ≤ measGP s)
    (hvmax : ∀ s, 0 < vmax s)
    (herror : ∀ s, |finiteRho s - witnessG2 s| ≤
      (arcCount s : ℝ) / (vmax s : ℝ))
    (hbudget : ∀ s, (arcCount s : ℝ) < delta s * (vmax s : ℝ))
    (s : Shape) :
    0 < finiteRho s := by
  have hfloor : delta s ≤ witnessG2 s :=
    p0_margin_le_witnessG2_shapes nuShape measGP p0Shape cap delta
      hbonf hDp0 hp0cap hmeasGP s
  exact finite_witness_pos_from_arc_error_mul
    (finiteRho s) (witnessG2 s) (delta s) (arcCount s) (vmax s)
    (hvmax s) hfloor (herror s) (hbudget s)

/-- Direct p0-margin version with a uniform arc bound and a uniform margin
floor.  This packages the incoming `#arcs(GOOD)` period-bounded signal as a
Lean-facing threshold obligation: prove `arcCount s ≤ arcBound`, prove
`deltaFloor ≤ delta s`, then require only
`arcBound < deltaFloor * vmax s`. -/
theorem finite_witness_pos_from_p0_margin_uniform_arc_bound_shapes
    (nuShape measGP p0Shape cap delta finiteRho : Shape → ℝ)
    (arcCount vmax : Shape → ℕ) (arcBound : ℕ) (deltaFloor : ℝ)
    (hbonf : ∀ s, nuShape s + measGP s - 1 ≤ witnessG2 s)
    (hDp0  : ∀ s, (1 - nuShape s) ≤ p0Shape s)
    (hp0cap : ∀ s, p0Shape s ≤ cap s - delta s)
    (hmeasGP : ∀ s, cap s ≤ measGP s)
    (hdelta : ∀ s, deltaFloor ≤ delta s)
    (hvmax : ∀ s, 0 < vmax s)
    (harc : ∀ s, arcCount s ≤ arcBound)
    (herror : ∀ s, |finiteRho s - witnessG2 s| ≤
      (arcCount s : ℝ) / (vmax s : ℝ))
    (hbudget : ∀ s, (arcBound : ℝ) < deltaFloor * (vmax s : ℝ))
    (s : Shape) :
    0 < finiteRho s := by
  have hfloor : ∀ s, delta s ≤ witnessG2 s :=
    p0_margin_le_witnessG2_shapes nuShape measGP p0Shape cap delta
      hbonf hDp0 hp0cap hmeasGP
  exact finite_witness_pos_from_uniform_arc_bound_shapes
    finiteRho delta arcCount vmax arcBound deltaFloor hfloor hdelta hvmax
    harc herror hbudget s

/-- Conditional LRC14 assembly through a finite-ruler Part-A node.  This is the
formal shape suggested by the S30 arc-count signal: p0 margin plus a
`#arcs/Vmax` error budget gives positive finite witness density; a finite
positive-witness criterion then gives `Mreach >= 1/14`. -/
theorem lrc14_from_finite_partA_p0_margin_shapes
    (nuShape measGP p0Shape cap delta finiteRho : Shape → ℝ)
    (arcCount vmax : Shape → ℕ)
    (hbonf : ∀ s, nuShape s + measGP s - 1 ≤ witnessG2 s)
    (hDp0  : ∀ s, (1 - nuShape s) ≤ p0Shape s)
    (hp0cap : ∀ s, p0Shape s ≤ cap s - delta s)
    (hmeasGP : ∀ s, cap s ≤ measGP s)
    (herror : ∀ s, |finiteRho s - witnessG2 s| ≤
      (arcCount s : ℝ) / (vmax s : ℝ))
    (hbudget : ∀ s, (arcCount s : ℝ) / (vmax s : ℝ) < delta s)
    (hfinitePartA : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      0 < finiteRho (shapeOf v) → (1 : ℝ) / 14 ≤ Mreach v) :
    LRC14Statement := by
  intro v hv
  have hpos : 0 < finiteRho (shapeOf v) :=
    finite_witness_pos_from_p0_margin_shapes
      nuShape measGP p0Shape cap delta finiteRho arcCount vmax
      hbonf hDp0 hp0cap hmeasGP herror hbudget (shapeOf v)
  exact lonely_of_Mreach_ge v hv (hfinitePartA v hv hpos)

/-- Conditional LRC14 assembly through finite-ruler Part A with the current
small/large branch split.  Small clusters use the already isolated `m_P` witness
floor; large clusters use the p0-wide-bound margin.  Thus the finite-ruler
error budget also splits into the two natural thresholds:
`#arcs/Vmax < m_P` on the small branch and `#arcs/Vmax < delta` on the large
branch. -/
theorem lrc14_from_finite_partA_p0_margin_split_shapes
    (nuShape measGP p0Shape cap delta finiteRho : Shape → ℝ)
    (arcCount vmax : Shape → ℕ)
    (hsmall : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      clusterSize (shapeOf v) ≤ 7 →
      (witnessMP : ℝ) ≤ witnessG2 (shapeOf v))
    (hbonf : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      nuShape s + measGP s - 1 ≤ witnessG2 s)
    (hDp0  : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (1 - nuShape s) ≤ p0Shape s)
    (hp0cap : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      p0Shape s ≤ cap s - delta s)
    (hmeasGP : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      cap s ≤ measGP s)
    (herror : ∀ s, |finiteRho s - witnessG2 s| ≤
      (arcCount s : ℝ) / (vmax s : ℝ))
    (hsmallBudget : ∀ s, clusterSize s ≤ 7 →
      (arcCount s : ℝ) / (vmax s : ℝ) < (witnessMP : ℝ))
    (hlargeBudget : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (arcCount s : ℝ) / (vmax s : ℝ) < delta s)
    (hsize : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      clusterSize (shapeOf v) ≤ 13)
    (hfinitePartA : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      0 < finiteRho (shapeOf v) → (1 : ℝ) / 14 ≤ Mreach v) :
    LRC14Statement := by
  intro v hv
  have hpos : 0 < finiteRho (shapeOf v) := by
    by_cases h7 : clusterSize (shapeOf v) ≤ 7
    · exact finite_witness_pos_from_arc_error
        (finiteRho (shapeOf v)) (witnessG2 (shapeOf v)) (witnessMP : ℝ)
        (arcCount (shapeOf v)) (vmax (shapeOf v))
        (hsmall v hv h7) (herror (shapeOf v))
        (hsmallBudget (shapeOf v) h7)
    · have h8 : 8 ≤ clusterSize (shapeOf v) := by omega
      have h13 : clusterSize (shapeOf v) ≤ 13 := hsize v hv
      have hfloor : delta (shapeOf v) ≤ witnessG2 (shapeOf v) := by
        have h1 : nuShape (shapeOf v) + measGP (shapeOf v) - 1 ≤
            witnessG2 (shapeOf v) := hbonf (shapeOf v) h8 h13
        have h2 : (1 - nuShape (shapeOf v)) ≤ p0Shape (shapeOf v) :=
          hDp0 (shapeOf v) h8 h13
        have h3 : p0Shape (shapeOf v) ≤
            cap (shapeOf v) - delta (shapeOf v) :=
          hp0cap (shapeOf v) h8 h13
        have h4 : cap (shapeOf v) ≤ measGP (shapeOf v) :=
          hmeasGP (shapeOf v) h8 h13
        linarith
      exact finite_witness_pos_from_arc_error
        (finiteRho (shapeOf v)) (witnessG2 (shapeOf v)) (delta (shapeOf v))
        (arcCount (shapeOf v)) (vmax (shapeOf v))
        hfloor (herror (shapeOf v)) (hlargeBudget (shapeOf v) h8 h13)
  exact lonely_of_Mreach_ge v hv (hfinitePartA v hv hpos)

/-! ## Axiom audit -/

#print axioms finite_witness_pos_from_abs_error
#print axioms finite_witness_pos_from_arc_error
#print axioms arc_div_lt_delta_of_lt_mul
#print axioms finite_witness_pos_from_arc_error_mul
#print axioms arc_div_lt_delta_of_le_bound
#print axioms finite_witness_pos_from_arc_bound
#print axioms p0_margin_le_witnessG2_shapes
#print axioms finite_witness_pos_from_floor_shapes
#print axioms finite_witness_pos_from_uniform_arc_bound_shapes
#print axioms finite_witness_pos_from_p0_margin_shapes
#print axioms finite_witness_pos_from_p0_margin_shapes_mul
#print axioms finite_witness_pos_from_p0_margin_uniform_arc_bound_shapes
#print axioms lrc14_from_finite_partA_p0_margin_shapes
#print axioms lrc14_from_finite_partA_p0_margin_split_shapes

end PartA
end LRC14
end LonelyRunner
