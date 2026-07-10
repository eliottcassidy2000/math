/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: kind-pasteur (LRC multi-agent project, 2026-07-10-S127)
-/
import Mathlib
import TournamentH7.LRCTightRigidity
import TournamentH7.LRCResidualMeasureFloorPrimitive
import TournamentH7.LRCDissociatedAssembly

/-!
# Wiring the difference-primitive collapse into opus's dissociated peel (kps-S127)

opus-S209 (`LRCDissociatedAssembly`) peels the `d = 2,3` near-dilate minimizers to THM-678
(`MultiDetunedDispatch`) and leaves the analytic floor obligation on the DISSOCIATED residual
(`ResidualObligationDissoc`).  Both that surface and opus-S207's weaker `SafeMeasureFloorPrimitive` carry
the clause `tupleGcd v = 1` — which is **exactly** my `Primitive v` (`primitive_of_tupleGcd_one`).  So the
S127 collapse `dilated_primitive_eq_range` applies verbatim on opus's surface, and this file composes them.

**The composition.** On the primitive residual the collapse says the tight locus is *literally* `{1,…,13}`
(no dilate); but every residual family is `GapFamily` (ratio `> 13`), and `{1,…,13}` has ratio exactly `13`
(`gapFamily_image_ne_range`).  So **no primitive residual family is tight** — the S127
`PrimitiveTightRigidity` discharges opus's weakest floor:

  `PrimitiveTightRigidity → SafeMeasureFloorPrimitive → ResidualObligationPrimitive`
                                                       → `ResidualObligationDissoc` (a fortiori)
                                                       → `LRC14Statement` (opus's dissoc assembly).

**Honest boundary — and what the wiring reveals.**  `PrimitiveTightRigidity` is the S127 open conjecture
(`≥ LRC(14)`); this file is a *reduction*, not a proof of the floor.  What it makes precise is the *relative
strength* of the two open routes: the tight rigidity is strictly **stronger** than the floor opus needs —
it discharges the WHOLE primitive residual, so the THM-678 dispatch (`hMD`) is **redundant on this route**
(`lrc14_of_primitiveTightRigidity`, `hMD`-free, proves the same `LRC14Statement`).  The moral for the
fleet: `SafeMeasureFloorPrimitive` (opus) is the minimal analytic target; the AP extremal rigidity is
sufficient overkill.  Prove the floor, not the rigidity.

Kernel-pure: no `sorry`, no `native_decide`.  Axioms: `[propext, Classical.choice, Quot.sound]`.
-/

namespace LonelyRunner
namespace LRC14Grand

open MeasureTheory

/-- **The bridge.**  opus's tuple-primitivity `tupleGcd v = 1` is exactly the S127 `Primitive v` (no
`d ≥ 2` divides every speed).  If some `d ≥ 2` divided all `vᵢ`, then `d.natAbs` would divide every
`(vᵢ).natAbs`, hence their `gcd = tupleGcd v = 1` — impossible for `d ≥ 2`. -/
theorem primitive_of_tupleGcd_one {v : Fin 13 → ℤ} (h : LRC14.tupleGcd v = 1) : Primitive v := by
  intro d hd hall
  have hnat : ∀ i, d.natAbs ∣ (v i).natAbs := fun i => Int.natAbs_dvd_natAbs.mpr (hall i)
  have hdvd : d.natAbs ∣ LRC14.tupleGcd v := Finset.dvd_gcd (fun i _ => hnat i)
  rw [h] at hdvd
  have h1 : d.natAbs = 1 := Nat.dvd_one.mp hdvd
  omega

/-- **The AP is excluded by the scale gap.**  The collapse target `{1,…,13}` is a dilated interval
(`c = 1`), so a `GapFamily` (ratio `> 13`) can never have `image |v| = {1,…,13}`.  This is the concrete
content of `not_dilated_of_gapFamily` at the AP. -/
theorem gapFamily_image_ne_range {v : Fin 13 → ℤ} (hgap : GapFamily v) :
    (Finset.univ.image fun i => |v i|) ≠ Finset.Icc 1 13 := by
  intro heq
  apply not_dilated_of_gapFamily hgap
  refine ⟨1, one_pos, ?_⟩
  rw [heq]
  ext x
  simp

/-- **The S127 primitive rigidity discharges opus's weakest floor.**  For a primitive residual family, if
the safe set had measure zero the collapse would force `image |v| = {1,…,13}` — contradicting `GapFamily`.
Hence every primitive residual family has a positive-measure safe set. -/
theorem safeMeasureFloorPrimitive_of_primitiveTightRigidity
    (h : PrimitiveTightRigidity) : SafeMeasureFloorPrimitive := by
  intro v hgcd hres
  obtain ⟨hv, _hcov, hgap, _hcomp, _hdist, _hlarge, _hdiv, _hcoarse, _hcres⟩ := hres
  by_cases hz : volume (safePeriod v) = 0
  · exact absurd (h v hv (primitive_of_tupleGcd_one hgcd) hz) (gapFamily_image_ne_range hgap)
  · exact pos_iff_ne_zero.mpr hz

/-- The primitive residual obligation, from the S127 rigidity (via opus's floor bridge). -/
theorem residualObligationPrimitive_of_primitiveTightRigidity
    (h : PrimitiveTightRigidity) : ResidualObligationPrimitive :=
  residualObligationPrimitive_of_measureFloorPrimitive
    (safeMeasureFloorPrimitive_of_primitiveTightRigidity h)

/-- **The requested wire: the S127 collapse discharges opus's DISSOCIATED residual obligation.**  A
fortiori — the rigidity floors the whole primitive residual, so the extra dissociation clause is free. -/
theorem residualObligationDissoc_of_primitiveTightRigidity
    (h : PrimitiveTightRigidity) : ResidualObligationDissoc := by
  intro v hv hgcd hcov hgap hcomp hdist hlarge hdiv hcoarse hcres _hdissoc
  exact residualObligationPrimitive_of_primitiveTightRigidity h
    v hv hgcd hcov hgap hcomp hdist hlarge hdiv hcoarse hcres

/-- **LRC(14) through opus's dissociated assembly, from the S127 primitive rigidity.**  The explicit
composition requested: cite + THM-678 dispatch + the primitive collapse ⟹ `LRC14Statement`. -/
theorem lrc14_of_primitiveTightRigidity_dissoc (cite : LRCUpTo13)
    (hMD : MultiDetunedDispatch) (h : PrimitiveTightRigidity) : LRC14.LRC14Statement :=
  lrc14_grand_assembly_dissoc cite hMD (residualObligationDissoc_of_primitiveTightRigidity h)

/-- **The `hMD`-free route — the honest headline.**  The S127 rigidity discharges the entire primitive
residual, so the THM-678 near-dilate dispatch is *not needed* to close LRC(14) on this route.  This locates
`PrimitiveTightRigidity` strictly above opus's dissociated floor in strength: it is sufficient overkill,
and `SafeMeasureFloorPrimitive` (opus) is the minimal analytic target. -/
theorem lrc14_of_primitiveTightRigidity (cite : LRCUpTo13)
    (h : PrimitiveTightRigidity) : LRC14.LRC14Statement :=
  lrc14_grand_assembly_primitive cite (residualObligationPrimitive_of_primitiveTightRigidity h)

end LRC14Grand
end LonelyRunner
