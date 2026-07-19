import Mathlib

/-!
# Centered blocker address compression (THM-1248)

This module checks the algebraic and rational-arithmetic core of the finite
relative-address quotient.  THM-1240 supplies the centered spokes and blocker
cycle; THM-1233 supplies the compact ratio bound.  The first-tooth-boundary
continuity argument and the interpretation as strict danger/safety sets remain
the explicit paper provider.
-/

namespace LRC14.CenteredBlockerAddressCompression

/-- Adding the carrier-gap address to the blocker tooth produces the exact
central-remainder cocycle. -/
theorem lower_gap_blocker_cocycle
    (P Qi Qj c dj k N X r M ell : ℤ)
    (hQj : Qj = c + dj)
    (hX : X = c*P - k*Qi)
    (hr : r = dj*P - N*Qi)
    (hM : M = k + N)
    (hell : ell = X + r) :
    P*Qj - M*Qi = ell := by
  rw [hQj, hX, hr, hM, hell]
  ring

/-- The lower-gap central remainder is invariant under an integral lift of
the circle time. -/
theorem lower_gap_cocycle_translation_invariant
    (P Qi Qj c dj k N n : ℤ) (hQj : Qj = c + dj) :
    (P+n*Qi)*Qj - ((k+n*c)+(N+n*dj))*Qi =
      P*Qj - (k+N)*Qi := by
  rw [hQj]
  ring

/-- The phase determinant differs from its positive central remainder by the
relative address digit times the source clock. -/
theorem relative_address_determinant
    (P Pj Qi Qj M ell delta : ℤ)
    (hell : ell = P*Qj - M*Qi)
    (hdelta : delta = Pj - M) :
    P*Qj - Pj*Qi = ell - delta*Qi := by
  rw [hell, hdelta]
  ring

/-- Any common divisor of the two spoke clocks divides the central
remainder.  This is the exact connection to the full gcd/torsion sheet. -/
theorem central_remainder_common_clock
    (P M Qi Qj ell g : ℤ)
    (hQi : g ∣ Qi) (hQj : g ∣ Qj)
    (hell : ell = P*Qj - M*Qi) :
    g ∣ ell := by
  rw [hell]
  exact dvd_sub (dvd_mul_of_dvd_right hQj P) (dvd_mul_of_dvd_right hQi M)

/-- The blocker residue expressed in centered rounding errors. -/
theorem recentered_blocker_identity
    (c dj Qi Qj t0 P Pj N k eps epsj beta delta : ℚ)
    (hQi0 : Qi ≠ 0)
    (hQj : Qj = c + dj)
    (ht0 : c*t0 = k + 1/2)
    (hP : P = Qi*t0 + eps)
    (hPj : Pj = Qj*t0 + epsj)
    (hbeta : beta = dj*(P/Qi) - N)
    (hdelta : delta = Pj - (k+N)) :
    beta = delta - 1/2 + (dj/Qi)*eps - epsj := by
  rw [hQj] at hPj
  rw [hP, hPj, hbeta, hdelta]
  field_simp [hQi0]
  linear_combination Qi * ht0

/-- The carrier central band plus strict blocker danger gives the sharp
positive remainder band `(5/28,23/28)`. -/
theorem central_remainder_band
    (Q X r ell : ℝ)
    (hQ : 0 < Q)
    (hXlo : Q/4 < X) (hXhi : X < 3*Q/4)
    (hrlo : -Q/14 < r) (hrhi : r < Q/14)
    (hell : ell = X+r) :
    5*Q/28 < ell ∧ ell < 23*Q/28 := by
  constructor <;> nlinarith

/-- Ratio-sensitive pre-integer bounds for the recentered digit. -/
theorem relative_digit_ratio_bounds
    (a eps epsj beta delta : ℝ)
    (ha : 0 ≤ a)
    (hepslo : -1/2 ≤ eps) (hepshi : eps ≤ 1/2)
    (hepsjlo : -1/2 ≤ epsj) (hepsjhi : epsj ≤ 1/2)
    (hbetalo : -1/14 < beta) (hbetahi : beta < 1/14)
    (hdelta : delta = 1/2-a*eps+epsj+beta) :
    -a/2-1/14 < delta ∧ delta < 1+a/2+1/14 := by
  have hlo := mul_le_mul_of_nonneg_left hepslo ha
  have hhi := mul_le_mul_of_nonneg_left hepshi ha
  constructor <;> nlinarith

/-- THM-1233's `a<2345/2` ratio consequence puts every real digit strictly
between `-587` and `588`. -/
theorem relative_digit_global_real_bounds
    (a eps epsj beta delta : ℝ)
    (ha : 0 ≤ a) (hacap : 2*a < 2345)
    (hepslo : -1/2 ≤ eps) (hepshi : eps ≤ 1/2)
    (hepsjlo : -1/2 ≤ epsj) (hepsjhi : epsj ≤ 1/2)
    (hbetalo : -1/14 < beta) (hbetahi : beta < 1/14)
    (hdelta : delta = 1/2-a*eps+epsj+beta) :
    -587 < delta ∧ delta < 588 := by
  have hb := relative_digit_ratio_bounds a eps epsj beta delta ha
    hepslo hepshi hepsjlo hepsjhi hbetalo hbetahi hdelta
  constructor <;> nlinarith [hb.1, hb.2]

theorem relative_digit_integer_bank
    (delta : ℤ) (hlo : (-587 : ℤ) < delta) (hhi : delta < 588) :
    -586 ≤ delta ∧ delta ≤ 587 := by
  omega

/-- Every edge with `d_j/Q_i≤1` has a binary relative digit. -/
theorem relative_digit_binary_real_bounds
    (a eps epsj beta delta : ℝ)
    (ha : 0 ≤ a) (ha1 : a ≤ 1)
    (hepslo : -1/2 ≤ eps) (hepshi : eps ≤ 1/2)
    (hepsjlo : -1/2 ≤ epsj) (hepsjhi : epsj ≤ 1/2)
    (hbetalo : -1/14 < beta) (hbetahi : beta < 1/14)
    (hdelta : delta = 1/2-a*eps+epsj+beta) :
    -4/7 < delta ∧ delta < 11/7 := by
  have hb := relative_digit_ratio_bounds a eps epsj beta delta ha
    hepslo hepshi hepsjlo hepsjhi hbetalo hbetahi hdelta
  constructor <;> nlinarith [hb.1, hb.2]

theorem relative_digit_binary_integer
    (delta : ℤ) (hlo : (-1 : ℤ) < delta) (hhi : delta < 2) :
    delta = 0 ∨ delta = 1 := by
  omega

/-- A binary digit is exactly the chronological side bit between the two
centered spoke phases. -/
theorem binary_digit_phase_separation
    (theta diff : ℝ) (delta : ℤ)
    (htheta_lo : 5/28 < theta) (htheta_hi : theta < 23/28)
    (hdiff : diff = theta-delta)
    (hdelta : delta = 0 ∨ delta = 1) :
    5/28 < |diff| ∧ |diff| < 23/28 := by
  rcases hdelta with rfl | rfl
  · simp only [Int.cast_zero, sub_zero] at hdiff
    rw [hdiff, abs_of_pos (by nlinarith)]
    exact ⟨htheta_lo, htheta_hi⟩
  · simp only [Int.cast_one] at hdiff
    rw [hdiff, abs_of_neg (by nlinarith)]
    constructor <;> nlinarith

/-- Three-edge instance of the contracting affine rounding-error transport. -/
theorem affine_triangle_transport
    (e1 e2 e3 a1 a2 a3 b1 b2 b3 : ℝ)
    (h1 : e2 = a1*e1+b1)
    (h2 : e3 = a2*e2+b2)
    (h3 : e1 = a3*e3+b3) :
    (1-a1*a2*a3)*e1 = a3*a2*b1+a3*b2+b3 := by
  linear_combination h3 + a3*h2 + a3*a2*h1

theorem contraction_gap_constant :
    1-((2345 : ℚ)/2346)^2 = 4691/5503716 := by
  norm_num

/-- The central positive defects force positive product holonomy on a
three-cycle.  The general cycle is the same telescoping identity. -/
theorem positive_central_triangle_holonomy
    (Q1 Q2 Q3 P1 P2 P3 M1 M2 M3 ell1 ell2 ell3 : ℝ)
    (hQ1 : 0 < Q1)
    (hP1 : 0 < P1) (hP2 : 0 < P2)
    (hM2 : 0 < M2) (hM3 : 0 < M3)
    (hell1 : 0 < ell1) (hell2 : 0 < ell2) (hell3 : 0 < ell3)
    (hr1 : P1*Q2-M1*Q1=ell1)
    (hr2 : P2*Q3-M2*Q2=ell2)
    (hr3 : P3*Q1-M3*Q3=ell3) :
    M1*M2*M3 < P1*P2*P3 := by
  have hid : (P1*P2*P3-M1*M2*M3)*Q1 =
      M2*M3*ell1+P1*M3*ell2+P1*P2*ell3 := by
    linear_combination M2*M3*hr1 + P1*M3*hr2 + P1*P2*hr3
  have hmul : 0 < (P1*P2*P3-M1*M2*M3)*Q1 := by
    rw [hid]
    positivity
  have hdiff : 0 < P1*P2*P3-M1*M2*M3 :=
    pos_of_mul_pos_right hmul (le_of_lt hQ1)
  exact sub_pos.mp hdiff

/-- Arithmetic consumer for the first target-tooth boundary: an ascending
source loses less than `1/7` of depth and remains above `1/14`. -/
theorem wall_event_source_margin :
    (1 : ℚ)/4-1/7 = 3/28 ∧ (1 : ℚ)/14 < 3/28 := by
  norm_num

#print axioms lower_gap_blocker_cocycle
#print axioms lower_gap_cocycle_translation_invariant
#print axioms relative_address_determinant
#print axioms central_remainder_common_clock
#print axioms recentered_blocker_identity
#print axioms central_remainder_band
#print axioms relative_digit_ratio_bounds
#print axioms relative_digit_global_real_bounds
#print axioms relative_digit_integer_bank
#print axioms relative_digit_binary_real_bounds
#print axioms relative_digit_binary_integer
#print axioms binary_digit_phase_separation
#print axioms affine_triangle_transport
#print axioms contraction_gap_constant
#print axioms positive_central_triangle_holonomy
#print axioms wall_event_source_margin

end LRC14.CenteredBlockerAddressCompression
