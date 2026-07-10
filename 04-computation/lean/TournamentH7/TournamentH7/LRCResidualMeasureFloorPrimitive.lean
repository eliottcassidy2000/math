/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-10-S207)
-/
import Mathlib
import TournamentH7.LRCResidualMeasureFloor
import TournamentH7.LRCPrimitiveAssembly

/-!
# The measure floor, restated on the PRIMITIVE residual (the sharpest form)

kps-S127's `SafeMeasureFloor` (`LRCResidualMeasureFloor`) reduces LRC(14) to "every residual family has a
safe set of positive measure". opus-S206's primitivity peel (`lrc14_grand_assembly_primitive`) reduces
LRC(14) to PRIMITIVE families. Composing them gives the **weakest floor hypothesis in the corpus**:

  `SafeMeasureFloorPrimitive` — every residual family that is ALSO primitive (`tupleGcd v = 1`) has
  `0 < volume (safePeriod v)` — and `lrc14_of_measureFloor_primitive : LRCUpTo13 → SafeMeasureFloorPrimitive
  → LRC14Statement`.

**Why primitivity is the right restriction (opus-S205, machine-checked).** The unrestricted residual
domain admits non-primitive dilates: `2·[1,2,3,4,5,6,7,8,9,11,12,13,20]` satisfies every residual clause
with `gcd = 2`, and since `α ↦ c·α` is measure-preserving, `μ = μ(core) = 1/980` with the core already
window-censused (`Vmax = 20`). Those dilates drive `inf μ → 0` — so a UNIFORM floor `μ ≥ μ₀` is FALSE on
the unrestricted residual but the per-family floor `μ > 0` still holds. On the PRIMITIVE residual the
dilates are gone: `inf μ` is bounded away from `0` (adversarial search: `inf μ ≈ 0.0085`, vs iid
`(6/7)^13 ≈ 0.1348`), with the minimizers concentrated at small `Vmax` and `μ → (6/7)^13` as `Vmax → ∞`
(decorrelation; `Vmax > 30 ⟹ μ ≥ 0.044`). So on the primitive residual BOTH the per-family floor `μ > 0`
(needed here, via kps's direct bridge) AND a uniform floor `μ ≥ μ₀` (which would additionally power klein's
THM-685 liveness route) are true and well-posed.

This is a *reduction*, not a proof of the floor: the floor on the primitive residual is the open analytic
core. But it is the tightest target — strictly fewer families than `SafeMeasureFloor`.

Kernel-pure: no `sorry`, no `native_decide`. Axioms: `[propext, Classical.choice, Quot.sound]`.
-/

namespace LonelyRunner
namespace LRC14Grand

open MeasureTheory

/-- **The sole remaining analytic obligation, on the PRIMITIVE residual.** Every residual family that is
also primitive (`tupleGcd v = 1`) has a positive-measure safe set. Strictly weaker than kps-S127's
`SafeMeasureFloor` (the non-primitive dilates need not be floored — they reduce to their cores). -/
def SafeMeasureFloorPrimitive : Prop :=
  ∀ v : Fin 13 → ℤ, LRC14.tupleGcd v = 1 → IsResidual v → 0 < volume (safePeriod v)

/-- The primitive residual obligation follows from the primitive measure floor (kps's direct bridge
`lonely_of_safePeriod_measure_pos`, now handed `tupleGcd v = 1`). -/
theorem residualObligationPrimitive_of_measureFloorPrimitive
    (hfloor : SafeMeasureFloorPrimitive) : ResidualObligationPrimitive := by
  intro v hv hgcd hcov hgap hcomp hdist hlarge hdiv hcoarse hres
  exact lonely_of_safePeriod_measure_pos
    (hfloor v hgcd ⟨hv, hcov, hgap, hcomp, hdist, hlarge, hdiv, hcoarse, hres⟩)

/-- **LRC(14) from the citation and the PRIMITIVE measure floor** — the sharpest reduction: the safe set
need be shown positive-measure only for the primitive residual (the dilates, which forced `inf μ = 0`, are
peeled by `lrc14_grand_assembly_primitive`). -/
theorem lrc14_of_measureFloor_primitive (cite : LRCUpTo13)
    (hfloor : SafeMeasureFloorPrimitive) : LRC14.LRC14Statement :=
  lrc14_grand_assembly_primitive cite
    (residualObligationPrimitive_of_measureFloorPrimitive hfloor)

/-- The primitive floor is genuinely weaker: kps's `SafeMeasureFloor` implies it (forget primitivity). So
switching to the primitive target loses nothing and strictly shrinks the family set to be floored. -/
theorem safeMeasureFloorPrimitive_of_safeMeasureFloor (h : SafeMeasureFloor) :
    SafeMeasureFloorPrimitive :=
  fun v _ hres => h v hres

end LRC14Grand
end LonelyRunner
