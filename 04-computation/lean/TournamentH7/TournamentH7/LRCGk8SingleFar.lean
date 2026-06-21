/-
  TournamentH7.LRCGk8SingleFar -- arithmetic kernel for the gK8 single-far reduction.

  CONTEXT (HYP-2829, mac-mini-2026-06-22-S24).  The gK8 concentration extremality
  `max_E L_yK8 ≤ 10·cap_k` (the open node `gK8_concentration_extremality` of the
  LRC(14) skeleton) has its BINDING wide case at far-count `r = #{e>14} = 1`
  (single far runner).  Stratified by `r`, the maximum of `L_yK8 = 10q0+q3+10q6`
  is dominated by the bounded `r=0` configurations (consec), and the binding wide
  value is the single-far supremum -- which clears `10·cap` with a COMFORTABLE
  margin `~1` (vs THM-563's razor-thin `0.13` for `p0`).  The analytic content
  (the THM-563 periodicity: `far·(deviation of each qₜ)` is periodic in `far`, so
  the supremum is a finite window-max plus a `1/far` tail) lives in
  `04-computation/lrc_gk8_singlefar_bound_macmini_0622s24.py`.

  This module records the EXACT rational margins as a `native_decide` checksum.
  It does NOT formalize the periodicity scan; it certifies the arithmetic that the
  bounded binding `L_yK8(consec_k)` and the single-far window supremum both sit
  strictly below `10·cap_k`, for the binding rows `k = 8, 9, 10`.

  Rationals as `(num, den)` with the comparison `a/b < c/d  ↔  a·d < c·b`
  (all denominators positive), checked over `ℕ`.

  BUILD NOTE: arithmetic verified exact in Python (all 9 ℕ-inequalities true); follows the
  proven `LRCQ6Contraction.lean` native_decide pattern.  Author lacked a local Lean toolchain --
  please `lake build TournamentH7.LRCGk8SingleFar` to confirm and add to the root import.
-/

namespace LonelyRunner
namespace Gk8SingleFar

/-! ### Bounded binding `L_yK8(consec_k)` vs `10·cap_k`. -/

/-- k=8: `L_yK8(consec_8) = 2633/735 < 10·cap_8 = 2243/588`. -/
theorem bounded_k8_below_cap : (2633 : Nat) * 588 < 2243 * 735 := by native_decide

/-- k=9: `L_yK8(consec_9) = 3259/735 < 10·cap_9 = 9895/2002`. -/
theorem bounded_k9_below_cap : (3259 : Nat) * 2002 < 9895 * 735 := by native_decide

/-- k=10: `L_yK8(consec_10) = 37/7 < 10·cap_10 = 40/7`. -/
theorem bounded_k10_below_cap : (37 : Nat) * 7 < 40 * 7 := by native_decide

/-! ### Single-far window supremum `sup_{15≤f≤120} L_yK8(consec_{k-1} ∪ {f})` vs `10·cap_k`.
The window supremum is the binding wide value (attained at `f=21,21,22`). -/

/-- k=8: single-far sup `2323/980 < 10·cap_8 = 2243/588` (margin `≈ 1.44`). -/
theorem singlefar_k8_below_cap : (2323 : Nat) * 588 < 2243 * 980 := by native_decide

/-- k=9: single-far sup `2876/735 < 10·cap_9 = 9895/2002` (margin `≈ 1.03`). -/
theorem singlefar_k9_below_cap : (2876 : Nat) * 2002 < 9895 * 735 := by native_decide

/-- k=10: single-far sup `62267/12936 < 10·cap_10 = 40/7` (margin `≈ 0.90`). -/
theorem singlefar_k10_below_cap : (62267 : Nat) * 7 < 40 * 12936 := by native_decide

/-! ### The single-far value sits strictly BELOW the bounded binding (far-count
monotonicity `r=1 < r=0`), confirming the binding is the bounded check. -/

/-- k=8: single-far `2323/980 < L_yK8(consec_8) = 2633/735`. -/
theorem singlefar_k8_below_bounded : (2323 : Nat) * 735 < 2633 * 980 := by native_decide

/-- k=9: single-far `2876/735 < L_yK8(consec_9) = 3259/735`. -/
theorem singlefar_k9_below_bounded : (2876 : Nat) < 3259 := by native_decide

/-- k=10: single-far `62267/12936 < L_yK8(consec_10) = 37/7`. -/
theorem singlefar_k10_below_bounded : (62267 : Nat) * 7 < 37 * 12936 := by native_decide

end Gk8SingleFar
end LonelyRunner
