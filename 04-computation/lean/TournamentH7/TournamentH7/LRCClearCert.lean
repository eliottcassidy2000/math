/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-05-S93)
-/
import TournamentH7.LRCDescentSurface

/-!
# Clear-component certificates: the kernel gate for the descent surface (HYP-4216)

`strict_lonely_of_clear_component` consumes a pointwise-strict base hypothesis on an
interval. This file provides its KERNEL-CHECKABLE certificate: the interval
`(A/Q, (A+LL)/Q)` is strictly `1/13`-clear for integer speed `w` iff every tooth index
in the finite window misses it — a pure-INTEGER disjunction per index (the S47 lesson:
never make the kernel normalize rationals). `clear_of_cert` turns the decidable
certificate into the pointwise hypothesis; a base of ≤ 5 runners needs one certificate
per runner on the chosen component, then the descent surface fires.
-/

namespace LonelyRunner
namespace ClearCert

/-- The integer certificate: tooth `m` of speed `w` misses `(A/Q, (A+LL)/Q)` —
entirely right of it or entirely left of it, in `13Q`-scaled integers. -/
def toothMiss (w A LL Q : ℤ) (m : ℤ) : Prop :=
  13 * w * (A + LL) < Q * (13 * m - 1) ∨ Q * (13 * m + 1) < 13 * w * A

instance (w A LL Q m : ℤ) : Decidable (toothMiss w A LL Q m) := by
  unfold toothMiss; infer_instance

/-- **The certificate gate**: if every tooth index in the covering window misses, the
open interval is strictly `1/13`-clear for `w`, pointwise. -/
theorem clear_of_cert (w A LL Q : ℤ) (hw : 0 < w) (hQ : 0 < Q) (_hLL : 0 < LL)
    (hcert : ∀ m : ℤ, w * A - Q ≤ Q * m → Q * m ≤ w * (A + LL) + Q →
      toothMiss w A LL Q m) :
    ∀ t : ℝ, (A : ℝ) / Q < t → t < ((A : ℝ) + LL) / Q →
      ∀ m : ℤ, (1 : ℝ) / 13 < |(w : ℝ) * t - m| := by
  intro t ht1 ht2 m
  have hQR : (0 : ℝ) < (Q : ℝ) := by exact_mod_cast hQ
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  rw [div_lt_iff₀ hQR] at ht1
  rw [lt_div_iff₀ hQR] at ht2
  -- w*t ∈ (w*A/Q, w*(A+LL)/Q)
  have hwt1 : (w : ℝ) * A < (w : ℝ) * (t * Q) := by
    apply mul_lt_mul_of_pos_left ht1 hwR
  have hwt2 : (w : ℝ) * (t * Q) < (w : ℝ) * ((A : ℝ) + LL) := by
    apply mul_lt_mul_of_pos_left ht2 hwR
  by_cases hm : w * A - Q ≤ Q * m ∧ Q * m ≤ w * (A + LL) + Q
  · have hmiss := hcert m hm.1 hm.2
    rw [lt_abs]
    rcases hmiss with h | h
    · -- interval entirely LEFT of the tooth: m - w t > 1/13
      right
      have hR : 13 * (w : ℝ) * ((A : ℝ) + LL) < (Q : ℝ) * (13 * m - 1) := by
        exact_mod_cast h
      -- w t Q < w (A + LL): 13 w t Q < 13 w (A+LL) < Q(13m − 1): 13 w t < 13 m − 1
      nlinarith [hwt2, hQR]
    · -- entirely RIGHT: w t - m > 1/13
      left
      have hR : (Q : ℝ) * (13 * m + 1) < 13 * (w : ℝ) * A := by
        exact_mod_cast h
      nlinarith [hwt1, hQR]
  · -- m outside the covering window: distance > 1 > 1/13
    push Not at hm
    rw [lt_abs]
    rcases lt_or_ge (Q * m) (w * A - Q) with hout | hge
    · -- Q m ≤ w A − Q − 1 ⟹ m + 1 < w A / Q < w t: distance > 1 > 1/13
      left
      have hZ : Q * m + 1 ≤ w * A - Q := hout
      have hZR : (Q : ℝ) * m + 1 ≤ (w : ℝ) * A - Q := by exact_mod_cast hZ
      nlinarith [hwt1, hQR, hZR]
    · have hout2 := hm hge
      right
      have hZ : w * (A + LL) + Q + 1 ≤ Q * m := hout2
      have hZR : (w : ℝ) * ((A : ℝ) + LL) + Q + 1 ≤ (Q : ℝ) * m := by exact_mod_cast hZ
      nlinarith [hwt2, hQR]

#print axioms clear_of_cert

end ClearCert
end LonelyRunner
