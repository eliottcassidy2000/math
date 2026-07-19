import TournamentH7.LRCEssentialRegion

/-!
# Nested carrier-window arithmetic (THM-1212)

This module kernel-checks the complete rational inequality ladder used to
close the five- and six-carrier `r = 6` strata.  The measurable cover input is
`LRC14.multi_speed_density_bound` from `LRCEssentialRegion`; the lemmas below
are its exact dimensionless consumers.

The remaining formal composition boundary is geometric: instantiate the
density theorem on each nested closed interval and extract the carrier-safe
boundary consumed by the residue-owner chart.  No numerical certificate is
hidden in this file.
-/

namespace LRC14.NestedCarrierWindow

/-- Length of `[1/(14x), 13/(14q)]`, before the common factor `1/x` is
restored: if `r=q/x`, then `x * |J(x,q)| = phaseLength r`. -/
noncomputable def phaseLength (r : ℝ) : ℝ := (13 - r) / (14 * r)

/-- The elementary geometry behind every nested window: if
`x ≤ o ≤ q` and `τ ∈ [1/(14x), 13/(14q)]`, then the protected carrier has
phase in the closed safe band `[1/14,13/14]`.  This is the no-wrap prefix
window lemma used before each density handoff. -/
theorem prefix_window_phase_bounds {x q o τ : ℝ}
    (hx : 0 < x) (hq : 0 < q) (hxo : x ≤ o) (hoq : o ≤ q)
    (hτLeft : 1 / (14 * x) ≤ τ) (hτRight : τ ≤ 13 / (14 * q)) :
    1 / 14 ≤ o * τ ∧ o * τ ≤ 13 / 14 := by
  have ho : 0 < o := lt_of_lt_of_le hx hxo
  have h14x : 0 < (14 : ℝ) * x := by positivity
  have h14q : 0 < (14 : ℝ) * q := by positivity
  constructor
  · have hbase : (1 : ℝ) / 14 ≤ o / (14 * x) := by
      rw [div_le_div_iff₀ (by norm_num : (0 : ℝ) < 14) h14x]
      nlinarith
    have hmul := mul_le_mul_of_nonneg_left hτLeft (le_of_lt ho)
    have : o / (14 * x) ≤ o * τ := by
      simpa [div_eq_mul_inv] using hmul
    exact hbase.trans this
  · have hmul := mul_le_mul_of_nonneg_left hτRight (le_of_lt ho)
    calc
      o * τ ≤ o * (13 / (14 * q)) := hmul
      _ = (13 * o) / (14 * q) := by ring
      _ ≤ 13 / 14 := by
        rw [div_le_div_iff₀ h14q (by norm_num : (0 : ℝ) < 14)]
        nlinarith

lemma phaseLength_pos {r : ℝ} (hr : 0 < r) (hr13 : r < 13) :
    0 < phaseLength r := by
  unfold phaseLength
  exact div_pos (by linarith) (by positivity)

/-- The five-carrier residual is killed in one step:
`r < 14/9` makes the protected two-carrier window too long for three later
combs to cover under the density budget. -/
theorem rho5_density_contradiction {r : ℝ} (hr : 0 < r)
    (hres : r < 14 / 9) :
    3 / r < 4 * phaseLength r := by
  have hden : 0 < 7 * r := by positivity
  have hid : 4 * phaseLength r - 3 / r = (5 - 2 * r) / (7 * r) := by
    unfold phaseLength
    field_simp
    <;> ring
  have hnum : 0 < 5 - 2 * r := by nlinarith
  have hpos : 0 < (5 - 2 * r) / (7 * r) := div_pos hnum hden
  rw [← hid] at hpos
  linarith

private lemma numerator_pos_of_positive_quotient {a b : ℝ}
    (hb : 0 < b) (h : 0 < a / b) : 0 < a := by
  rcases (div_pos_iff.mp h) with hab | hab
  · exact hab.1
  · exact False.elim ((not_lt_of_ge (le_of_lt hb)) hab.2)

/-- First six-carrier handoff.  A four-comb cover of the two-carrier window,
together with `r=y/x<35/12`, would force `s=z/x<1960/363`. -/
theorem rho6_first_handoff {r s : ℝ} (hr : 0 < r) (hs : 0 < s)
    (hres : r < 35 / 12)
    (hcover : 3 * phaseLength r < 4 / s) :
    s < 1960 / 363 := by
  have hden : 0 < 14 * r * s := by positivity
  have hid : 4 / s - 3 * phaseLength r =
      (56 * r - 3 * s * (13 - r)) / (14 * r * s) := by
    unfold phaseLength
    field_simp
    <;> ring
  have hquot : 0 <
      (56 * r - 3 * s * (13 - r)) / (14 * r * s) := by
    rw [← hid]
    linarith
  have hnum : 0 < 56 * r - 3 * s * (13 - r) :=
    numerator_pos_of_positive_quotient hden hquot
  have hrden : (121 / 12 : ℝ) < 13 - r := by nlinarith
  have hs3 : 0 < 3 * s := by positivity
  have hprod : 3 * s * (121 / 12) < 3 * s * (13 - r) := by
    exact mul_lt_mul_of_pos_left hrden hs3
  have hrnum : 56 * r < 56 * (35 / 12 : ℝ) := by nlinarith
  nlinarith

lemma first_handoff_rounding : (1960 / 363 : ℝ) < 27 / 5 := by norm_num

/-- Second handoff.  The rounded input `s<27/5` turns a three-comb cover into
the exact bound `t=w/x<567/76`. -/
theorem rho6_second_handoff {s t : ℝ} (hs : 0 < s) (ht : 0 < t)
    (hsBound : s < 27 / 5)
    (hcover : 4 * phaseLength s < 3 / t) :
    t < 567 / 76 := by
  have hden : 0 < 14 * s * t := by positivity
  have hid : 3 / t - 4 * phaseLength s =
      (42 * s - 4 * t * (13 - s)) / (14 * s * t) := by
    unfold phaseLength
    field_simp
    <;> ring
  have hquot : 0 <
      (42 * s - 4 * t * (13 - s)) / (14 * s * t) := by
    rw [← hid]
    linarith
  have hnum : 0 < 42 * s - 4 * t * (13 - s) :=
    numerator_pos_of_positive_quotient hden hquot
  have hsden : (38 / 5 : ℝ) < 13 - s := by nlinarith
  have ht4 : 0 < 4 * t := by positivity
  have hprod : 4 * t * (38 / 5) < 4 * t * (13 - s) := by
    exact mul_lt_mul_of_pos_left hsden ht4
  have hsnum : 42 * s < 42 * (27 / 5 : ℝ) := by nlinarith
  nlinarith

lemma second_handoff_inside : (567 / 76 : ℝ) < 13 := by norm_num

/-- Third handoff.  A two-comb cover of the four-carrier window forces
`v/x<15876/2105`. -/
theorem rho6_third_handoff {t v : ℝ} (ht : 0 < t) (hv : 0 < v)
    (htBound : t < 567 / 76)
    (hcover : 5 * phaseLength t < 2 / v) :
    v < 15876 / 2105 := by
  have hden : 0 < 14 * t * v := by positivity
  have hid : 2 / v - 5 * phaseLength t =
      (28 * t - 5 * v * (13 - t)) / (14 * t * v) := by
    unfold phaseLength
    field_simp
    <;> ring
  have hquot : 0 <
      (28 * t - 5 * v * (13 - t)) / (14 * t * v) := by
    rw [← hid]
    linarith
  have hnum : 0 < 28 * t - 5 * v * (13 - t) :=
    numerator_pos_of_positive_quotient hden hquot
  have htden : (421 / 76 : ℝ) < 13 - t := by nlinarith
  have hv5 : 0 < 5 * v := by positivity
  have hprod : 5 * v * (421 / 76) < 5 * v * (13 - t) := by
    exact mul_lt_mul_of_pos_left htden hv5
  have htnum : 28 * t < 28 * (567 / 76 : ℝ) := by nlinarith
  nlinarith

lemma third_handoff_below_eleven : (15876 / 2105 : ℝ) < 11 := by norm_num

/-- Once the fifth normalized carrier has ratio below `11`, the final
protected closed interval is longer than a tooth of that carrier. -/
theorem final_window_longer_than_tooth {v : ℝ} (hv : 0 < v)
    (hv11 : v < 11) :
    1 / (7 * v) < phaseLength v := by
  have hden : 0 < 14 * v := by positivity
  have hid : phaseLength v - 1 / (7 * v) = (11 - v) / (14 * v) := by
    unfold phaseLength
    field_simp
    <;> ring
  have hpos : 0 < (11 - v) / (14 * v) := by
    exact div_pos (by linarith) hden
  rw [← hid] at hpos
  linarith

/-- Kernel-checked composition of the complete dimensionless ladder. -/
theorem rho6_nested_ratio_ladder {r s t v : ℝ}
    (hr : 0 < r) (hs : 0 < s) (ht : 0 < t) (hv : 0 < v)
    (hres : r < 35 / 12)
    (h2 : 3 * phaseLength r < 4 / s)
    (h3 : 4 * phaseLength s < 3 / t)
    (h4 : 5 * phaseLength t < 2 / v) :
    v < 15876 / 2105 ∧ 1 / (7 * v) < phaseLength v := by
  have hsExact := rho6_first_handoff hr hs hres h2
  have hsRound : s < 27 / 5 := lt_trans hsExact first_handoff_rounding
  have htExact := rho6_second_handoff hs ht hsRound h3
  have hvExact := rho6_third_handoff ht hv htExact h4
  have hv11 : v < 11 := lt_trans hvExact third_handoff_below_eleven
  exact ⟨hvExact, final_window_longer_than_tooth hv hv11⟩

/-! ### Set-level consumers

The preceding lemmas isolate the dimensionless arithmetic.  The next theorem
crosses one of the geometric composition boundaries explicitly: it applies
`multi_speed_density_bound` to the actual closed two-carrier window and proves
that the three later carrier combs in the `rho = 5` residual cannot cover it.
-/

/-- In the five-carrier residual, the actual closed prefix window
`[1/(14x), 13/(14y)]` is not covered by the three later danger combs.

This is the set-level form of the argument: `multi_speed_density_bound` is
instantiated on the concrete `Finset {z,w,v}` and the concrete interval, rather
than supplied as a pre-normalized arithmetic hypothesis. -/
theorem rho5_actual_window_not_covered {x y z w v : ℕ}
    (hxy : x < y) (hyz : y < z) (hzw : z < w) (hwv : w < v)
    (hres : 9 * y < 14 * x) :
    ¬ Set.Icc (1 / (14 * (x : ℝ))) (13 / (14 * (y : ℝ))) ⊆
        ⋃ s ∈ ({z, w, v} : Finset ℕ), badArcs s (1 / 14) := by
  have hxNat : 0 < x := by
    by_contra hx
    have hx0 : x = 0 := Nat.eq_zero_of_not_pos hx
    simp [hx0] at hres
  have hyNat : 0 < y := lt_trans hxNat hxy
  have hx : (0 : ℝ) < x := by exact_mod_cast hxNat
  have hy : (0 : ℝ) < y := by exact_mod_cast hyNat
  have hresReal : 9 * (y : ℝ) < 14 * (x : ℝ) := by exact_mod_cast hres

  let left : ℝ := 1 / (14 * (x : ℝ))
  let right : ℝ := 13 / (14 * (y : ℝ))
  let L : ℝ := right - left
  have hLformula :
      L = (13 * (x : ℝ) - (y : ℝ)) / (14 * (x : ℝ) * (y : ℝ)) := by
    dsimp [L, left, right]
    field_simp
    <;> ring
  have hnum : 0 < 13 * (x : ℝ) - (y : ℝ) := by nlinarith
  have hL : 0 < L := by
    rw [hLformula]
    exact div_pos hnum (by positivity)
  have hright : left + L = right := by simp [L]

  let S : Finset ℕ := {z, w, v}
  have hzv : z < v := lt_trans hzw hwv
  have hzwNe : z ≠ w := Nat.ne_of_lt hzw
  have hzvNe : z ≠ v := Nat.ne_of_lt hzv
  have hwvNe : w ≠ v := Nat.ne_of_lt hwv
  have hcard : S.card = 3 := by simp [S, hzwNe, hzvNe, hwvNe]
  have hSpos : ∀ s ∈ S, 1 ≤ s := by
    intro s hs
    simp only [S, Finset.mem_insert, Finset.mem_singleton] at hs
    rcases hs with rfl | rfl | rfl <;> omega
  have hSy : ∀ s ∈ S, y ≤ s := by
    intro s hs
    simp only [S, Finset.mem_insert, Finset.mem_singleton] at hs
    rcases hs with rfl | rfl | rfl <;> omega

  intro hcover
  have hcover' : Set.Icc left (left + L) ⊆ ⋃ s ∈ S, badArcs s (1 / 14) := by
    rw [hright]
    simpa [left, right, S] using hcover
  have hdensity := multi_speed_density_bound S hSpos (1 / 14) L left
    (by norm_num) (le_of_lt hL) hcover'

  have hsum : ∑ s ∈ S, (1 : ℝ) / s ≤ (S.card : ℝ) / (y : ℝ) := by
    calc
      ∑ s ∈ S, (1 : ℝ) / s ≤ ∑ s ∈ S, (1 : ℝ) / y := by
        refine Finset.sum_le_sum (fun s hs => ?_)
        have hys : (y : ℝ) ≤ (s : ℝ) := by exact_mod_cast hSy s hs
        exact one_div_le_one_div_of_le hy hys
      _ = (S.card : ℝ) / (y : ℝ) := by
        rw [Finset.sum_const, nsmul_eq_mul]
        ring
  have hscale : 2 * (1 / 14 : ℝ) * (∑ s ∈ S, (1 : ℝ) / s)
      ≤ 2 * (1 / 14 : ℝ) * ((S.card : ℝ) / (y : ℝ)) :=
    mul_le_mul_of_nonneg_left hsum (by norm_num)
  have hbudget : 4 * L ≤ 3 / (y : ℝ) := by
    rw [hcard] at hdensity hscale
    norm_num at hdensity hscale ⊢
    linarith

  have hratio : (y : ℝ) / (x : ℝ) < 14 / 9 := by
    apply (div_lt_iff₀ hx).2
    nlinarith
  have hrho5 := rho5_density_contradiction (div_pos hy hx) hratio
  have hphase : phaseLength ((y : ℝ) / (x : ℝ)) = (x : ℝ) * L := by
    rw [hLformula]
    unfold phaseLength
    field_simp
    <;> ring
  have hreciprocal :
      3 / ((y : ℝ) / (x : ℝ)) = (x : ℝ) * (3 / (y : ℝ)) := by
    field_simp
    <;> ring
  have hlongScaled :
      (x : ℝ) * (3 / (y : ℝ)) < (x : ℝ) * (4 * L) := by
    rw [hreciprocal, hphase] at hrho5
    nlinarith
  have hlong : 3 / (y : ℝ) < 4 * L := by
    by_contra hnot
    have hreverse : 4 * L ≤ 3 / (y : ℝ) := le_of_not_gt hnot
    have hscaled : (x : ℝ) * (4 * L) ≤ (x : ℝ) * (3 / (y : ℝ)) :=
      mul_le_mul_of_nonneg_left hreverse (le_of_lt hx)
    exact (not_lt_of_ge hscaled) hlongScaled
  exact (not_lt_of_ge hbudget) hlong

private lemma actual_window_length_formula {x q : ℕ}
    (hx : (0 : ℝ) < x) (hq : (0 : ℝ) < q) :
    13 / (14 * (q : ℝ)) - 1 / (14 * (x : ℝ)) =
      (13 * (x : ℝ) - (q : ℝ)) / (14 * (x : ℝ) * (q : ℝ)) := by
  field_simp
  <;> ring

private lemma actual_window_length_pos {x q : ℕ}
    (hx : (0 : ℝ) < x) (hq : (0 : ℝ) < q)
    (hq13 : (q : ℝ) < 13 * (x : ℝ)) :
    0 < 13 / (14 * (q : ℝ)) - 1 / (14 * (x : ℝ)) := by
  rw [actual_window_length_formula hx hq]
  exact div_pos (by linarith) (by positivity)

private lemma phaseLength_nat_ratio_eq_scaled_window {x q : ℕ}
    (hx : (0 : ℝ) < x) (hq : (0 : ℝ) < q) :
    phaseLength ((q : ℝ) / (x : ℝ)) =
      (x : ℝ) * (13 / (14 * (q : ℝ)) - 1 / (14 * (x : ℝ))) := by
  unfold phaseLength
  field_simp
  <;> ring

private lemma ratio_reciprocal_eq_scaled {x q : ℕ} (c : ℝ)
    (hx : (0 : ℝ) < x) (hq : (0 : ℝ) < q) :
    c / ((q : ℝ) / (x : ℝ)) = (x : ℝ) * (c / (q : ℝ)) := by
  field_simp
  <;> ring

/-- Set-level first handoff for six carriers.  If the four later combs cover
`J(x,y)`, the concrete density theorem and strict speed order force the first
ratio in the nested ladder. -/
theorem rho6_first_actual_cover_forces_ratio {x y z w v u : ℕ}
    (hxy : x < y) (hyz : y < z) (hzw : z < w) (hwv : w < v) (hvu : v < u)
    (hres : 12 * y < 35 * x)
    (hcover : Set.Icc (1 / (14 * (x : ℝ))) (13 / (14 * (y : ℝ))) ⊆
      ⋃ s ∈ ({z, w, v, u} : Finset ℕ), badArcs s (1 / 14)) :
    (z : ℝ) / (x : ℝ) < 1960 / 363 := by
  have hxNat : 0 < x := by omega
  have hyNat : 0 < y := lt_trans hxNat hxy
  have hzNat : 0 < z := lt_trans hyNat hyz
  have hx : (0 : ℝ) < x := by exact_mod_cast hxNat
  have hy : (0 : ℝ) < y := by exact_mod_cast hyNat
  have hz : (0 : ℝ) < z := by exact_mod_cast hzNat
  have hresReal : 12 * (y : ℝ) < 35 * (x : ℝ) := by exact_mod_cast hres
  have hy13 : (y : ℝ) < 13 * (x : ℝ) := by nlinarith
  let left : ℝ := 1 / (14 * (x : ℝ))
  let L : ℝ := 13 / (14 * (y : ℝ)) - left
  have hL : 0 < L := by
    dsimp [L, left]
    exact actual_window_length_pos hx hy hy13
  let S : Finset ℕ := {z, w, v, u}
  have hzwNe : z ≠ w := Nat.ne_of_lt hzw
  have hzvNe : z ≠ v := Nat.ne_of_lt (lt_trans hzw hwv)
  have hzuNe : z ≠ u := Nat.ne_of_lt (lt_trans (lt_trans hzw hwv) hvu)
  have hwvNe : w ≠ v := Nat.ne_of_lt hwv
  have hwuNe : w ≠ u := Nat.ne_of_lt (lt_trans hwv hvu)
  have hvuNe : v ≠ u := Nat.ne_of_lt hvu
  have hcard : S.card = 4 := by
    simp [S, hzwNe, hzvNe, hzuNe, hwvNe, hwuNe, hvuNe]
  have hSpos : ∀ s ∈ S, 1 ≤ s := by
    intro s hs
    simp only [S, Finset.mem_insert, Finset.mem_singleton] at hs
    rcases hs with rfl | rfl | rfl | rfl <;> omega
  have hcover' : Set.Icc left (left + L) ⊆ ⋃ s ∈ S, badArcs s (1 / 14) := by
    simpa [left, L, S] using hcover
  have hdensity := multi_speed_density_bound S hSpos (1 / 14) L left
    (by norm_num) (le_of_lt hL) hcover'
  have hsum : ∑ s ∈ S, (1 : ℝ) / s < 4 / (z : ℝ) := by
    have hzwReal : (z : ℝ) < w := by exact_mod_cast hzw
    have hzvReal : (z : ℝ) < v := by exact_mod_cast (lt_trans hzw hwv)
    have hzuReal : (z : ℝ) < u := by exact_mod_cast (lt_trans (lt_trans hzw hwv) hvu)
    have hwRec : (1 : ℝ) / w < 1 / z := one_div_lt_one_div_of_lt hz hzwReal
    have hvRec : (1 : ℝ) / v < 1 / z := one_div_lt_one_div_of_lt hz hzvReal
    have huRec : (1 : ℝ) / u < 1 / z := one_div_lt_one_div_of_lt hz hzuReal
    have hsumExpand : ∑ s ∈ S, (1 : ℝ) / s =
        1 / (z : ℝ) + (1 / (w : ℝ) + (1 / (v : ℝ) + 1 / (u : ℝ))) := by
      simp [S, hzwNe, hzvNe, hzuNe, hwvNe, hwuNe, hvuNe]
    rw [hsumExpand]
    calc
      1 / (z : ℝ) + (1 / (w : ℝ) + (1 / (v : ℝ) + 1 / (u : ℝ)))
          < 1 / (z : ℝ) + (1 / (z : ℝ) + (1 / (z : ℝ) + 1 / (z : ℝ))) := by
            linarith
      _ = 4 / (z : ℝ) := by ring
  have hactual : 3 * L < 4 / (z : ℝ) := by
    have hscaled := mul_le_mul_of_nonneg_left hdensity (show (0 : ℝ) ≤ 7 by norm_num)
    rw [hcard] at hscaled
    have hbase : 3 * L ≤ ∑ s ∈ S, (1 : ℝ) / s := by
      convert hscaled using 1 <;> ring
    exact lt_of_le_of_lt hbase hsum
  have hphase := phaseLength_nat_ratio_eq_scaled_window hx hy
  have hrec := ratio_reciprocal_eq_scaled (x := x) (q := z) 4 hx hz
  have hnormalized :
      3 * phaseLength ((y : ℝ) / (x : ℝ)) < 4 / ((z : ℝ) / (x : ℝ)) := by
    rw [hphase, hrec]
    convert mul_lt_mul_of_pos_left hactual hx using 1 <;> ring
  have hratio : (y : ℝ) / (x : ℝ) < 35 / 12 := by
    apply (div_lt_iff₀ hx).2
    nlinarith
  exact rho6_first_handoff (div_pos hy hx) (div_pos hz hx) hratio hnormalized

/-- Set-level second handoff for six carriers.  A cover of `J(x,z)` by the
three remaining combs, together with the rounded first ratio, forces the
second ratio in the ladder. -/
theorem rho6_second_actual_cover_forces_ratio {x z w v u : ℕ}
    (hxz : x < z) (hzw : z < w) (hwv : w < v) (hvu : v < u)
    (hzBound : 5 * z < 27 * x)
    (hcover : Set.Icc (1 / (14 * (x : ℝ))) (13 / (14 * (z : ℝ))) ⊆
      ⋃ s ∈ ({w, v, u} : Finset ℕ), badArcs s (1 / 14)) :
    (w : ℝ) / (x : ℝ) < 567 / 76 := by
  have hxNat : 0 < x := by omega
  have hzNat : 0 < z := lt_trans hxNat hxz
  have hwNat : 0 < w := lt_trans hzNat hzw
  have hx : (0 : ℝ) < x := by exact_mod_cast hxNat
  have hz : (0 : ℝ) < z := by exact_mod_cast hzNat
  have hw : (0 : ℝ) < w := by exact_mod_cast hwNat
  have hzBoundReal : 5 * (z : ℝ) < 27 * (x : ℝ) := by exact_mod_cast hzBound
  have hz13 : (z : ℝ) < 13 * (x : ℝ) := by nlinarith
  let left : ℝ := 1 / (14 * (x : ℝ))
  let L : ℝ := 13 / (14 * (z : ℝ)) - left
  have hL : 0 < L := by
    dsimp [L, left]
    exact actual_window_length_pos hx hz hz13
  let S : Finset ℕ := {w, v, u}
  have hwvNe : w ≠ v := Nat.ne_of_lt hwv
  have hwuNe : w ≠ u := Nat.ne_of_lt (lt_trans hwv hvu)
  have hvuNe : v ≠ u := Nat.ne_of_lt hvu
  have hcard : S.card = 3 := by simp [S, hwvNe, hwuNe, hvuNe]
  have hSpos : ∀ s ∈ S, 1 ≤ s := by
    intro s hs
    simp only [S, Finset.mem_insert, Finset.mem_singleton] at hs
    rcases hs with rfl | rfl | rfl <;> omega
  have hcover' : Set.Icc left (left + L) ⊆ ⋃ s ∈ S, badArcs s (1 / 14) := by
    simpa [left, L, S] using hcover
  have hdensity := multi_speed_density_bound S hSpos (1 / 14) L left
    (by norm_num) (le_of_lt hL) hcover'
  have hsum : ∑ s ∈ S, (1 : ℝ) / s < 3 / (w : ℝ) := by
    have hwvReal : (w : ℝ) < v := by exact_mod_cast hwv
    have hwuReal : (w : ℝ) < u := by exact_mod_cast (lt_trans hwv hvu)
    have hvRec : (1 : ℝ) / v < 1 / w := one_div_lt_one_div_of_lt hw hwvReal
    have huRec : (1 : ℝ) / u < 1 / w := one_div_lt_one_div_of_lt hw hwuReal
    have hsumExpand : ∑ s ∈ S, (1 : ℝ) / s =
        1 / (w : ℝ) + (1 / (v : ℝ) + 1 / (u : ℝ)) := by
      simp [S, hwvNe, hwuNe, hvuNe]
    rw [hsumExpand]
    calc
      1 / (w : ℝ) + (1 / (v : ℝ) + 1 / (u : ℝ))
          < 1 / (w : ℝ) + (1 / (w : ℝ) + 1 / (w : ℝ)) := by linarith
      _ = 3 / (w : ℝ) := by ring
  have hactual : 4 * L < 3 / (w : ℝ) := by
    have hscaled := mul_le_mul_of_nonneg_left hdensity (show (0 : ℝ) ≤ 7 by norm_num)
    rw [hcard] at hscaled
    have hbase : 4 * L ≤ ∑ s ∈ S, (1 : ℝ) / s := by
      convert hscaled using 1 <;> ring
    exact lt_of_le_of_lt hbase hsum
  have hphase := phaseLength_nat_ratio_eq_scaled_window hx hz
  have hrec := ratio_reciprocal_eq_scaled (x := x) (q := w) 3 hx hw
  have hnormalized :
      4 * phaseLength ((z : ℝ) / (x : ℝ)) < 3 / ((w : ℝ) / (x : ℝ)) := by
    rw [hphase, hrec]
    convert mul_lt_mul_of_pos_left hactual hx using 1 <;> ring
  have hratio : (z : ℝ) / (x : ℝ) < 27 / 5 := by
    apply (div_lt_iff₀ hx).2
    nlinarith
  exact rho6_second_handoff (div_pos hz hx) (div_pos hw hx) hratio hnormalized

/-- Set-level third handoff for six carriers.  A cover of `J(x,w)` by the
last two combs, together with the second ratio, forces the final protected
prefix ratio. -/
theorem rho6_third_actual_cover_forces_ratio {x w v u : ℕ}
    (hxw : x < w) (hwv : w < v) (hvu : v < u)
    (hwBound : 76 * w < 567 * x)
    (hcover : Set.Icc (1 / (14 * (x : ℝ))) (13 / (14 * (w : ℝ))) ⊆
      ⋃ s ∈ ({v, u} : Finset ℕ), badArcs s (1 / 14)) :
    (v : ℝ) / (x : ℝ) < 15876 / 2105 := by
  have hxNat : 0 < x := by omega
  have hwNat : 0 < w := lt_trans hxNat hxw
  have hvNat : 0 < v := lt_trans hwNat hwv
  have hx : (0 : ℝ) < x := by exact_mod_cast hxNat
  have hw : (0 : ℝ) < w := by exact_mod_cast hwNat
  have hv : (0 : ℝ) < v := by exact_mod_cast hvNat
  have hwBoundReal : 76 * (w : ℝ) < 567 * (x : ℝ) := by exact_mod_cast hwBound
  have hw13 : (w : ℝ) < 13 * (x : ℝ) := by nlinarith
  let left : ℝ := 1 / (14 * (x : ℝ))
  let L : ℝ := 13 / (14 * (w : ℝ)) - left
  have hL : 0 < L := by
    dsimp [L, left]
    exact actual_window_length_pos hx hw hw13
  let S : Finset ℕ := {v, u}
  have hvuNe : v ≠ u := Nat.ne_of_lt hvu
  have hcard : S.card = 2 := by simp [S, hvuNe]
  have hSpos : ∀ s ∈ S, 1 ≤ s := by
    intro s hs
    simp only [S, Finset.mem_insert, Finset.mem_singleton] at hs
    rcases hs with rfl | rfl <;> omega
  have hcover' : Set.Icc left (left + L) ⊆ ⋃ s ∈ S, badArcs s (1 / 14) := by
    simpa [left, L, S] using hcover
  have hdensity := multi_speed_density_bound S hSpos (1 / 14) L left
    (by norm_num) (le_of_lt hL) hcover'
  have hsum : ∑ s ∈ S, (1 : ℝ) / s < 2 / (v : ℝ) := by
    have hvuReal : (v : ℝ) < u := by exact_mod_cast hvu
    have huRec : (1 : ℝ) / u < 1 / v := one_div_lt_one_div_of_lt hv hvuReal
    have hsumExpand : ∑ s ∈ S, (1 : ℝ) / s =
        1 / (v : ℝ) + 1 / (u : ℝ) := by simp [S, hvuNe]
    rw [hsumExpand]
    calc
      1 / (v : ℝ) + 1 / (u : ℝ) < 1 / (v : ℝ) + 1 / (v : ℝ) := by linarith
      _ = 2 / (v : ℝ) := by ring
  have hactual : 5 * L < 2 / (v : ℝ) := by
    have hscaled := mul_le_mul_of_nonneg_left hdensity (show (0 : ℝ) ≤ 7 by norm_num)
    rw [hcard] at hscaled
    have hbase : 5 * L ≤ ∑ s ∈ S, (1 : ℝ) / s := by
      convert hscaled using 1 <;> ring
    exact lt_of_le_of_lt hbase hsum
  have hphase := phaseLength_nat_ratio_eq_scaled_window hx hw
  have hrec := ratio_reciprocal_eq_scaled (x := x) (q := v) 2 hx hv
  have hnormalized :
      5 * phaseLength ((w : ℝ) / (x : ℝ)) < 2 / ((v : ℝ) / (x : ℝ)) := by
    rw [hphase, hrec]
    convert mul_lt_mul_of_pos_left hactual hx using 1 <;> ring
  have hratio : (w : ℝ) / (x : ℝ) < 567 / 76 := by
    apply (div_lt_iff₀ hx).2
    nlinarith
  exact rho6_third_handoff (div_pos hw hx) (div_pos hv hx) hratio hnormalized

/-- A closed interval strictly longer than one open danger tooth cannot be
contained in a single speed's danger comb.  This packages the endpoint
topology used at the last step of the six-carrier ladder. -/
theorem closed_interval_not_subset_badArcs_of_longer_than_tooth
    {speed : ℕ} (hspeed : 1 ≤ speed) {a b : ℝ} (hab : a ≤ b)
    (hlong : 1 / (7 * (speed : ℝ)) < b - a) :
    ¬ Set.Icc a b ⊆ badArcs speed (1 / 14) := by
  intro hcover
  have hle := Icc_length_le_of_subset_badArcs speed hspeed (1 / 14)
    (by norm_num) a b hab hcover
  have htooth : 2 * (1 / 14 : ℝ) / (speed : ℝ) = 1 / (7 * (speed : ℝ)) := by
    field_simp
    <;> ring
  rw [htooth] at hle
  exact (not_lt_of_ge hle) hlong

#print axioms rho5_density_contradiction
#print axioms prefix_window_phase_bounds
#print axioms rho6_first_handoff
#print axioms rho6_second_handoff
#print axioms rho6_third_handoff
#print axioms final_window_longer_than_tooth
#print axioms rho6_nested_ratio_ladder
#print axioms rho5_actual_window_not_covered
#print axioms rho6_first_actual_cover_forces_ratio
#print axioms rho6_second_actual_cover_forces_ratio
#print axioms rho6_third_actual_cover_forces_ratio
#print axioms closed_interval_not_subset_badArcs_of_longer_than_tooth

end LRC14.NestedCarrierWindow
