/-
  TournamentH7.LRCSlowFast — the slow-fast identity, algebraic core of the ruler embedding hembed
  (kind-pasteur-2026-07-09-S105).

  THM-527 Part A's "slow-fast change of variables" views time from the fastest runner `Vmax`.
  Writing `Vmax·τ = j + φ` (`j = ⌊Vmax·τ⌋` the ruler period, `φ = fract(Vmax·τ)` the fast phase),
  every runner `u = Vmax − e` (co-offset `e`) has

      ‖u·τ‖  =  ‖φ − e·τ‖                                     (SLOW-FAST IDENTITY)

  and, expanding `e·τ = e·j/Vmax + e·φ/Vmax`, `‖u·τ‖ = ‖(φ − c) − e·φ/Vmax‖` where `c = frac(e·j/Vmax)`
  is the tooth (klein's `teeth`).  So the runner distance is the fast-phase clearance `φ − c` MINUS the
  **slow-fast drift** `e·φ/Vmax`.  This file proves the exact identity and makes the drift explicit —
  the foundational reduction of `hembed` (`minReach v τ` ↔ teeth clearance), and the precise place the
  finite-`Vmax` obstruction lives: the drift is `≤ spread/Vmax`, negligible for small spread but `O(1)`
  in the good-period window (`Vmax ≈ 7·spread/6`), which is exactly why the embedding needs the
  density/arc-count argument there, not a single fast phase.  Self-contained (imports `LRCMreachConcrete`).
-/
import Mathlib
import TournamentH7.LRCMreachConcrete

namespace LonelyRunner
namespace LRC14Concrete

/-- `nearInt` is invariant under adding an integer (nearest-integer distance is `1`-periodic). -/
theorem nearInt_int_add (m : ℤ) (x : ℝ) : nearInt ((m : ℝ) + x) = nearInt x := by
  simp [nearInt, Int.fract_intCast_add]

/-- **The slow-fast identity.**  For a runner speed `u = Vmax − e` (integers) and any time `τ`, with
fast phase `φ = fract(Vmax·τ)`, the runner's nearest-integer distance equals the fast-phase distance
to the drifted co-offset: `‖u·τ‖ = ‖φ − e·τ‖`.  (Because `u·τ = Vmax·τ − e·τ = ⌊Vmax·τ⌋ + (φ − e·τ)`
and `nearInt` erases the integer `⌊Vmax·τ⌋`.) -/
theorem nearInt_speed_eq_phase_sub (Vmax e : ℤ) (τ : ℝ) :
    nearInt (((Vmax - e : ℤ) : ℝ) * τ) =
      nearInt (Int.fract ((Vmax : ℝ) * τ) - (e : ℝ) * τ) := by
  have hsplit : (((Vmax - e : ℤ) : ℝ) * τ)
      = ((⌊(Vmax : ℝ) * τ⌋ : ℤ) : ℝ) + (Int.fract ((Vmax : ℝ) * τ) - (e : ℝ) * τ) := by
    have hfr : (Vmax : ℝ) * τ = ((⌊(Vmax : ℝ) * τ⌋ : ℤ) : ℝ) + Int.fract ((Vmax : ℝ) * τ) := by
      rw [Int.floor_add_fract]
    push_cast
    push_cast at hfr
    linarith
  rw [hsplit, nearInt_int_add]

/-- **The drift is exactly `e·φ/Vmax`.**  Writing `Vmax·τ = j + φ`, the co-offset phase `e·τ` exceeds
the integer-grid tooth position `e·j/Vmax` by exactly the slow-fast drift `e·φ/Vmax`.  Combined with
the identity, `‖u·τ‖ = ‖(φ − c) − e·φ/Vmax‖` — the teeth clearance minus this drift.  The drift is
`≤ e/Vmax ≤ spread/Vmax`: negligible for small spread, `O(1)` in the good-period window
(`Vmax ≈ 7·spread/6`), which is the precise finite-`Vmax` obstruction. -/
theorem drift_eq (Vmax e j : ℤ) (φ : ℝ) (hV : (Vmax : ℝ) ≠ 0) :
    (e : ℝ) * ((j : ℝ) + φ) / Vmax - (e : ℝ) * j / Vmax = (e : ℝ) * φ / Vmax := by
  field_simp
  ring

end LRC14Concrete
end LonelyRunner
