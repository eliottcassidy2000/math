/-
  death-star-2026-07-07-S3 -- THE Good ⊆ goodSet BRIDGE.

  There are two "good set" formalisms for "max-gap of {frac(e·x)} > 1/7":
   · `LRCTailDiameter.Good (1/7) E`  (real-anchored: ∃ real a, ∀ e, frac(e·x−a) > 1/7) — what the
     muGood diameter floors (death-star AP20/muGood_diam_floor, boxeph AP30/AP44/AP76) bound;
   · `GoodSet.goodSet E`             (phase-anchored: ∃ c∈E, ∀ b, frac((b−c)x) ∉ (0,1/7]) — the
     CONCRETE measurable event that the skeleton's `witnessG2` measures (via LRCEventMeasureBridge).

  This file proves `Good (1/7) E.toFinset ⊆ goodSet E`, so `volume(goodSet ∩ [0,1]) ≥ muGood ≥ m_P`
  transfers EVERY diameter floor to the concrete event.  The bridge: take `c = argmax_e frac(e·x−a)`;
  then for every b, `u_b − u_c ≤ 0` with `u_b, u_c ∈ (1/7, 1)`, so `frac((b−c)x) = frac(u_b−u_c)`
  is `0` (if `u_b = u_c`) or `u_b − u_c + 1 ∈ (1/7, 1)` — never in `(0, 1/7]`.

  Kernel-pure: no sorry, no native_decide.
-/
import Mathlib
import TournamentH7.LRCTailDiameter
import TournamentH7.LRCGoodSet
import TournamentH7.LRCAP20Certificate

namespace LonelyRunner
namespace GoodSetBridge

open TailDiameter
open TournamentH7.GoodSet
open MeasureTheory
open scoped ENNReal

/-- **THE BRIDGE.**  The real-anchored `Good` set embeds in the concrete phase-anchored `goodSet`. -/
theorem Good_subset_goodSet (E : List ℤ) (hE : E.toFinset.Nonempty) :
    Good (1/7) E.toFinset ⊆ goodSet E := by
  intro x hx
  obtain ⟨a, ha⟩ := hx
  -- ha : ∀ e ∈ E.toFinset, 1/7 < Int.fract ((e:ℝ)*x - a)
  obtain ⟨c, hcmem, hcmax⟩ :=
    Finset.exists_max_image E.toFinset (fun e => Int.fract ((e:ℝ) * x - a)) hE
  -- unfold goodSet membership
  simp only [goodSet, Set.mem_iUnion, Set.mem_iInter, Set.mem_setOf_eq, exists_prop]
  refine ⟨c, hcmem, ?_⟩
  intro b hbmem
  -- reduce frac((b-c)x) to frac(u_b - u_c)
  have hcast : (((b - c : ℤ)) : ℝ) * x = ((b:ℝ) * x - a) - ((c:ℝ) * x - a) := by push_cast; ring
  have hkey : Int.fract (((b - c : ℤ) : ℝ) * x)
      = Int.fract (Int.fract ((b:ℝ) * x - a) - Int.fract ((c:ℝ) * x - a)) := by
    rw [hcast, fract_sub_eq_fract_fract_sub_fract]
  rw [hkey]
  have hub : (1:ℝ)/7 < Int.fract ((b:ℝ) * x - a) := ha b hbmem
  have huc : (1:ℝ)/7 < Int.fract ((c:ℝ) * x - a) := ha c hcmem
  have hub1 : Int.fract ((b:ℝ) * x - a) < 1 := Int.fract_lt_one _
  have huc1 : Int.fract ((c:ℝ) * x - a) < 1 := Int.fract_lt_one _
  have hle : Int.fract ((b:ℝ) * x - a) ≤ Int.fract ((c:ℝ) * x - a) := hcmax b hbmem
  rcases eq_or_lt_of_le hle with heq | hlt
  · -- u_b = u_c ⟹ frac(0) = 0 ∉ (0, 1/7]
    have : Int.fract ((b:ℝ) * x - a) - Int.fract ((c:ℝ) * x - a) = 0 := by rw [heq]; ring
    rw [this, Int.fract_zero]
    simp only [Set.mem_Ioc]
    rintro ⟨h, _⟩; exact absurd h (lt_irrefl 0)
  · -- u_b < u_c ⟹ diff ∈ (-6/7, 0), frac = diff + 1 ∈ (1/7, 1) > 1/7
    set d : ℝ := Int.fract ((b:ℝ) * x - a) - Int.fract ((c:ℝ) * x - a) with hd
    have hdlo : (-1:ℝ) ≤ d := by rw [hd]; linarith
    have hdhi : d < 0 := by rw [hd]; linarith
    have hfloor : ⌊d⌋ = (-1 : ℤ) := by
      refine (Int.floor_eq_iff).mpr ⟨?_, ?_⟩
      · push_cast; linarith
      · push_cast; linarith
    have hfr : Int.fract d = d + 1 := by unfold Int.fract; rw [hfloor]; push_cast; ring
    rw [hfr]
    simp only [Set.mem_Ioc, not_and, not_le]
    intro _; rw [hd]; linarith

/-- Measure form: the concrete good event has `volume ≥ muGood`. -/
theorem volume_goodSet_ge_muGood (E : List ℤ) (hE : E.toFinset.Nonempty) :
    muGood (1/7) E.toFinset ≤ volume (goodSet E ∩ Set.Icc (0:ℝ) 1) := by
  unfold muGood
  exact measure_mono (Set.inter_subset_inter_left _ (Good_subset_goodSet E hE))

/-- **THE CONCRETE-EVENT FLOOR (bounded diameter).**  For a nonempty family `E` whose translate
lies in `{0,…,D}` with `2 ≤ D ≤ 30`, the concrete `goodSet` event has measure `≥ m_P` — the
lower bound the skeleton's `witnessG2` (= `μ(goodSet ∩ G_P)`) needs, in the pure-cluster case. -/
theorem goodSet_measure_ge_mP (E : List ℤ) (hE : E.toFinset.Nonempty) (m D : ℤ)
    (hD2 : 2 ≤ D) (hD30 : D ≤ 30) (hEd : ∀ e ∈ E.toFinset, e - m ∈ Finset.Icc (0 : ℤ) D) :
    ENNReal.ofReal ((14249:ℝ)/252252) ≤ volume (goodSet E ∩ Set.Icc (0:ℝ) 1) :=
  le_trans (AP20Cert.hlarge_floor_of_diam_le E.toFinset m D hD2 hD30 hEd)
    (volume_goodSet_ge_muGood E hE)

/-! ### The `slowμ` form — the exact measure `witnessG2` uses (`= volume.restrict (Ico 0 1)`). -/

/-- The concrete good event's `slowμ`-measure is `≥ m_P` (bounded diameter).  `slowμ` is the
skeleton's period probability measure `volume.restrict (Ico 0 1)`; this is EXACTLY the measure the
`witnessG2` identification (`LRCEventMeasureBridge.hwitness`) uses. -/
theorem slowμ_goodSet_ge_mP (E : List ℤ) (hE : E.toFinset.Nonempty) (m D : ℤ)
    (hD2 : 2 ≤ D) (hD30 : D ≤ 30) (hEd : ∀ e ∈ E.toFinset, e - m ∈ Finset.Icc (0 : ℤ) D) :
    ENNReal.ofReal ((14249:ℝ)/252252) ≤ DenseCovers.slowμ (goodSet E) := by
  have hmeas : MeasurableSet (goodSet E) := measurableSet_goodSet E
  have hrestrict : DenseCovers.slowμ (goodSet E) = volume (goodSet E ∩ Set.Ico (0:ℝ) 1) := by
    rw [DenseCovers.slowμ, Measure.restrict_apply hmeas]
  rw [hrestrict]
  have haeq : volume (Good (1/7) E.toFinset ∩ Set.Icc (0:ℝ) 1)
      = volume (Good (1/7) E.toFinset ∩ Set.Ico (0:ℝ) 1) := by
    apply measure_congr
    exact (Filter.EventuallyEq.refl _ (Good (1/7) E.toFinset)).inter Ico_ae_eq_Icc.symm
  have hfloor := AP20Cert.hlarge_floor_of_diam_le E.toFinset m D hD2 hD30 hEd
  rw [show muGood (1/7) E.toFinset = volume (Good (1/7) E.toFinset ∩ Set.Icc (0:ℝ) 1) from rfl,
    haeq] at hfloor
  exact le_trans hfloor
    (measure_mono (Set.inter_subset_inter_left _ (Good_subset_goodSet E hE)))

/-- The witnessG2-ready real form: `(slowμ(goodSet E)).toReal ≥ m_P`.  In the pure-cluster case
(`G_P = univ`), `witnessG2 = (slowμ(goodSet cluster)).toReal`, so this IS the census floor
`witnessMP ≤ witnessG2` for a nonempty cluster of diameter ≤ 30. -/
theorem slowμ_goodSet_toReal_ge_mP (E : List ℤ) (hE : E.toFinset.Nonempty) (m D : ℤ)
    (hD2 : 2 ≤ D) (hD30 : D ≤ 30) (hEd : ∀ e ∈ E.toFinset, e - m ∈ Finset.Icc (0 : ℤ) D) :
    (14249:ℝ)/252252 ≤ (DenseCovers.slowμ (goodSet E)).toReal := by
  have h := slowμ_goodSet_ge_mP E hE m D hD2 hD30 hEd
  have hne : DenseCovers.slowμ (goodSet E) ≠ ⊤ := measure_ne_top _ _
  have := ENNReal.toReal_mono hne h
  rwa [ENNReal.toReal_ofReal (by norm_num)] at this

end GoodSetBridge
end LonelyRunner
