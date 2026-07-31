/-
TournamentH7.ArtanhSandwich — the certified artanh log-ratio sandwich.

The certified-inequality engine flagged by HYP-9023 ("the missing
certified-inequality layer"), whose external home is now CONFIRMED as AMM
Problem 12592's proposer-side certificate (27) (HYP-9061; results note
amm12592-snippet-context-confirmation-deathstar-coinC2.md).

Kernel-checked contents, no new axioms:

* `log_ratio_lower` : 0 ≤ t < 1 → 2(t + t³/3 + t⁵/5) ≤ log((1+t)/(1-t));
* `log_ratio_upper` : 0 ≤ t < 1 → log((1+t)/(1-t)) ≤
                        2(t + t³/3 + t⁵/(5(1-t²)));
* `certificate_27`  : 1/25 < (2457/6592)·log(8847357/2974400)
                        − log(1285/896)  — the decoded certificate (27);
* `log_two_lower`   : 842/1215 ≤ log 2  (the sandwich at t = 1/3);
* `mass_ordering_M62_M43` : 22/1215 ≤ 26·log 2 − 18  (THM-2000's
                        M(6,2) > M(4,3), float-free).

Mechanism: the lower/upper differences vanish at 0 and have exact
derivatives 2t⁶/(1−t²) and (4/5)t⁶/(1−t²)² on (−1,1); FTC-2 plus
nonnegativity of the integrand closes both directions.
-/
import Mathlib

open Real Set intervalIntegral

namespace ArtanhSandwich

/-- `(log((1+s)/(1-s)))' = 2/(1-t²)` on `(-1, 1)`. -/
lemma hasDerivAt_logRatio (t : ℝ) (h0 : -1 < t) (h1 : t < 1) :
    HasDerivAt (fun s : ℝ => Real.log ((1 + s) / (1 - s)))
      (2 / (1 - t ^ 2)) t := by
  have h1p : (0 : ℝ) < 1 + t := by linarith
  have h1m : (0 : ℝ) < 1 - t := by linarith
  have hnum : HasDerivAt (fun s : ℝ => 1 + s) 1 t := by
    simpa using (hasDerivAt_id t).const_add (1 : ℝ)
  have hden : HasDerivAt (fun s : ℝ => 1 - s) (-1) t := by
    simpa using (hasDerivAt_id t).const_sub (1 : ℝ)
  have hq : HasDerivAt (fun s : ℝ => (1 + s) / (1 - s))
      ((1 * (1 - t) - (1 + t) * (-1)) / (1 - t) ^ 2) t :=
    hnum.div hden h1m.ne'
  have hpos : (0 : ℝ) < (1 + t) / (1 - t) := div_pos h1p h1m
  have h := hq.log hpos.ne'
  have hne : (1 : ℝ) - t ^ 2 ≠ 0 := by nlinarith
  convert h using 1
  field_simp
  ring

/-- Lower-difference derivative:
`(log((1+s)/(1-s)) - 2(s + s³/3 + s⁵/5))' = 2t⁶/(1-t²)`. -/
lemma hasDerivAt_lowerDiff (t : ℝ) (h0 : -1 < t) (h1 : t < 1) :
    HasDerivAt
      (fun s : ℝ => Real.log ((1 + s) / (1 - s))
        - 2 * (s + s ^ 3 / 3 + s ^ 5 / 5))
      (2 * t ^ 6 / (1 - t ^ 2)) t := by
  have hpoly : HasDerivAt (fun s : ℝ => 2 * (s + s ^ 3 / 3 + s ^ 5 / 5))
      (2 * (1 + t ^ 2 + t ^ 4)) t := by
    have h1' : HasDerivAt (fun s : ℝ => s) 1 t := hasDerivAt_id t
    have h3 : HasDerivAt (fun s : ℝ => s ^ 3 / 3) (t ^ 2) t := by
      have h := (hasDerivAt_pow 3 t).div_const 3
      convert h using 1
      push_cast
      ring
    have h5 : HasDerivAt (fun s : ℝ => s ^ 5 / 5) (t ^ 4) t := by
      have h := (hasDerivAt_pow 5 t).div_const 5
      convert h using 1
      push_cast
      ring
    have hsum := (h1'.add h3).add h5
    exact hsum.const_mul (2 : ℝ)
  have h := (hasDerivAt_logRatio t h0 h1).sub hpoly
  convert h using 1
  have hne : (1 : ℝ) - t ^ 2 ≠ 0 := by nlinarith
  field_simp
  ring

/-- Upper-difference derivative:
`(2(s + s³/3 + s⁵/(5(1-s²))) - log((1+s)/(1-s)))' = (4/5)t⁶/(1-t²)²`. -/
lemma hasDerivAt_upperDiff (t : ℝ) (h0 : -1 < t) (h1 : t < 1) :
    HasDerivAt
      (fun s : ℝ => 2 * (s + s ^ 3 / 3 + s ^ 5 / (5 * (1 - s ^ 2)))
        - Real.log ((1 + s) / (1 - s)))
      (4 / 5 * t ^ 6 / (1 - t ^ 2) ^ 2) t := by
  have hne : (1 : ℝ) - t ^ 2 ≠ 0 := by nlinarith
  have hden5 : (5 : ℝ) * (1 - t ^ 2) ≠ 0 := by
    intro hc
    apply hne
    linarith [hc]
  have hd : HasDerivAt (fun s : ℝ => 5 * (1 - s ^ 2)) (5 * (-(2 * t))) t := by
    have hsq : HasDerivAt (fun s : ℝ => s ^ 2) (2 * t) t := by
      simpa using hasDerivAt_pow 2 t
    exact (hsq.neg.const_add (1 : ℝ)).const_mul (5 : ℝ)
  have hq : HasDerivAt (fun s : ℝ => s ^ 5 / (5 * (1 - s ^ 2)))
      ((5 * t ^ 4 * (5 * (1 - t ^ 2)) - t ^ 5 * (5 * (-(2 * t))))
        / (5 * (1 - t ^ 2)) ^ 2) t := by
    have hp5 : HasDerivAt (fun s : ℝ => s ^ 5) (5 * t ^ 4) t := by
      simpa using hasDerivAt_pow 5 t
    exact hp5.div hd hden5
  have hpoly : HasDerivAt
      (fun s : ℝ => 2 * (s + s ^ 3 / 3 + s ^ 5 / (5 * (1 - s ^ 2))))
      (2 * (1 + t ^ 2
        + (5 * t ^ 4 * (5 * (1 - t ^ 2)) - t ^ 5 * (5 * (-(2 * t))))
          / (5 * (1 - t ^ 2)) ^ 2)) t := by
    have h1' : HasDerivAt (fun s : ℝ => s) 1 t := hasDerivAt_id t
    have h3 : HasDerivAt (fun s : ℝ => s ^ 3 / 3) (t ^ 2) t := by
      have h := (hasDerivAt_pow 3 t).div_const 3
      convert h using 1
      push_cast
      ring
    have hsum := (h1'.add h3).add hq
    exact hsum.const_mul (2 : ℝ)
  have h := hpoly.sub (hasDerivAt_logRatio t h0 h1)
  convert h using 1
  field_simp
  ring

/-- **Lower bound of the sandwich.** -/
theorem log_ratio_lower {t : ℝ} (h0 : 0 ≤ t) (h1 : t < 1) :
    2 * (t + t ^ 3 / 3 + t ^ 5 / 5) ≤ Real.log ((1 + t) / (1 - t)) := by
  have key : (∫ s in (0:ℝ)..t, 2 * s ^ 6 / (1 - s ^ 2)) =
      (Real.log ((1 + t) / (1 - t)) - 2 * (t + t ^ 3 / 3 + t ^ 5 / 5))
      - (Real.log ((1 + 0) / (1 - 0))
          - 2 * ((0:ℝ) + 0 ^ 3 / 3 + 0 ^ 5 / 5)) := by
    apply intervalIntegral.integral_eq_sub_of_hasDerivAt
    · intro x hx
      rw [uIcc_of_le h0] at hx
      exact hasDerivAt_lowerDiff x (by linarith [hx.1])
        (lt_of_le_of_lt hx.2 h1)
    · apply ContinuousOn.intervalIntegrable
      apply ContinuousOn.div
      · fun_prop
      · fun_prop
      · intro x hx
        rw [uIcc_of_le h0] at hx
        have hx1 : x < 1 := lt_of_le_of_lt hx.2 h1
        have : 0 < 1 - x ^ 2 := by nlinarith [hx.1]
        exact this.ne'
  have hnn : 0 ≤ ∫ s in (0:ℝ)..t, 2 * s ^ 6 / (1 - s ^ 2) := by
    apply intervalIntegral.integral_nonneg h0
    intro u hu
    have hu1 : u < 1 := lt_of_le_of_lt hu.2 h1
    have hpos : 0 < 1 - u ^ 2 := by nlinarith [hu.1]
    positivity
  rw [key] at hnn
  norm_num [Real.log_one] at hnn
  linarith

/-- **Upper bound of the sandwich.** -/
theorem log_ratio_upper {t : ℝ} (h0 : 0 ≤ t) (h1 : t < 1) :
    Real.log ((1 + t) / (1 - t))
      ≤ 2 * (t + t ^ 3 / 3 + t ^ 5 / (5 * (1 - t ^ 2))) := by
  have key : (∫ s in (0:ℝ)..t, 4 / 5 * s ^ 6 / (1 - s ^ 2) ^ 2) =
      (2 * (t + t ^ 3 / 3 + t ^ 5 / (5 * (1 - t ^ 2)))
        - Real.log ((1 + t) / (1 - t)))
      - (2 * ((0:ℝ) + 0 ^ 3 / 3 + 0 ^ 5 / (5 * (1 - 0 ^ 2)))
        - Real.log ((1 + 0) / (1 - 0))) := by
    apply intervalIntegral.integral_eq_sub_of_hasDerivAt
    · intro x hx
      rw [uIcc_of_le h0] at hx
      exact hasDerivAt_upperDiff x (by linarith [hx.1])
        (lt_of_le_of_lt hx.2 h1)
    · apply ContinuousOn.intervalIntegrable
      apply ContinuousOn.div
      · fun_prop
      · fun_prop
      · intro x hx
        rw [uIcc_of_le h0] at hx
        have hx1 : x < 1 := lt_of_le_of_lt hx.2 h1
        have h2 : 0 < 1 - x ^ 2 := by nlinarith [hx.1]
        positivity
  have hnn : 0 ≤ ∫ s in (0:ℝ)..t, 4 / 5 * s ^ 6 / (1 - s ^ 2) ^ 2 := by
    apply intervalIntegral.integral_nonneg h0
    intro u hu
    have hu1 : u < 1 := lt_of_le_of_lt hu.2 h1
    have hpos : 0 < 1 - u ^ 2 := by nlinarith [hu.1]
    positivity
  rw [key] at hnn
  norm_num [Real.log_one] at hnn
  linarith

/-- **Certificate (27)** of AMM Problem 12592 (proposer-side notes,
decoded and context-confirmed): exactly fair, kernel-checked. -/
theorem certificate_27 :
    (1 / 25 : ℝ) <
      (2457 / 6592) * Real.log (8847357 / 2974400)
        - Real.log (1285 / 896) := by
  have htB0 : (0 : ℝ) ≤ 5872957 / 11821757 := by norm_num
  have htB1 : (5872957 / 11821757 : ℝ) < 1 := by norm_num
  have htA0 : (0 : ℝ) ≤ 389 / 2181 := by norm_num
  have htA1 : (389 / 2181 : ℝ) < 1 := by norm_num
  have hlow := log_ratio_lower htB0 htB1
  have hup := log_ratio_upper htA0 htA1
  have hargB : ((1 : ℝ) + 5872957 / 11821757) / (1 - 5872957 / 11821757)
      = 8847357 / 2974400 := by norm_num
  have hargA : ((1 : ℝ) + 389 / 2181) / (1 - 389 / 2181)
      = 1285 / 896 := by norm_num
  rw [hargB] at hlow
  rw [hargA] at hup
  have hmul : (2457 / 6592 : ℝ)
        * (2 * (5872957 / 11821757 + (5872957 / 11821757 : ℝ) ^ 3 / 3
          + (5872957 / 11821757 : ℝ) ^ 5 / 5))
      ≤ (2457 / 6592) * Real.log (8847357 / 2974400) := by
    apply mul_le_mul_of_nonneg_left hlow
    norm_num
  have hrat : (1 / 25 : ℝ) <
      (2457 / 6592 : ℝ)
        * (2 * (5872957 / 11821757 + (5872957 / 11821757 : ℝ) ^ 3 / 3
          + (5872957 / 11821757 : ℝ) ^ 5 / 5))
      - 2 * (389 / 2181 + (389 / 2181 : ℝ) ^ 3 / 3
          + (389 / 2181 : ℝ) ^ 5 / (5 * (1 - (389 / 2181 : ℝ) ^ 2))) := by
    norm_num
  linarith

/-- The sandwich at `t = 1/3`: a certified rational lower bound for
`log 2`. -/
theorem log_two_lower : (842 / 1215 : ℝ) ≤ Real.log 2 := by
  have h := log_ratio_lower (t := (1/3 : ℝ)) (by norm_num) (by norm_num)
  have harg : ((1 : ℝ) + 1/3) / (1 - 1/3) = 2 := by norm_num
  rw [harg] at h
  have : (842 / 1215 : ℝ)
      = 2 * ((1/3 : ℝ) + (1/3 : ℝ) ^ 3 / 3 + (1/3 : ℝ) ^ 5 / 5) := by
    norm_num
  linarith [h]

/-- THM-2000's figurate mass ordering `M(6,2) > M(4,3)`, float-free:
`M(6,2) - M(4,3) = 26 log 2 - 18 ≥ 22/1215 > 0`. -/
theorem mass_ordering_M62_M43 : (22 / 1215 : ℝ) ≤ 26 * Real.log 2 - 18 := by
  have h := log_two_lower
  nlinarith [h]

end ArtanhSandwich
