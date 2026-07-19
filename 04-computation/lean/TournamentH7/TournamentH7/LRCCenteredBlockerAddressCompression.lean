import Mathlib.Tactic

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
  calc
    P*Qj-M*Qi = (c*P-k*Qi)+(dj*P-N*Qi) := by
      rw [hQj, hM]
      ring
    _ = X+r := by rw [hX, hr]
    _ = ell := hell.symm

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
  calc
    beta = dj*(P/Qi)-N := hbeta
    _ = dj*t0+(dj/Qi)*eps-N := by
      rw [hP]
      field_simp [hQi0]
    _ = delta-1/2+(dj/Qi)*eps-epsj := by
      rw [hdelta, hPj, hQj]
      linear_combination -ht0

/-- The carrier central band plus strict blocker danger gives the sharp
positive remainder band `(5/28,23/28)`. -/
theorem central_remainder_band
    (Q X r ell : ℝ)
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

/-- The finite relative digit chooses the oriented tooth germ: nonpositive
digits put the target phase to the left, while digits at least one put it to
the right.  In the application `diff=Q_j(t_i-t_j)`. -/
theorem relative_digit_selects_germ
    (theta diff delta : ℝ)
    (htheta_lo : 0 < theta) (htheta_hi : theta < 1)
    (hdiff : diff = theta-delta) :
    (delta ≤ 0 → 0 < diff) ∧ (1 ≤ delta → diff < 0) := by
  constructor <;> intro h <;> nlinarith

/-- The signed left-minus-right tooth-germ imbalance is already encoded by
the translation-invariant central remainder and carrier coordinate. -/
theorem oriented_germ_imbalance
    (d Q X r ell : ℚ) (hell : ell = X+r) :
    (1/(14*d)+r/(d*Q))-(1/(14*d)-r/(d*Q)) =
      2*(ell-X)/(d*Q) := by
  rw [hell]
  ring

/-- The central carrier coordinate is part of the common centered residue
orbit: with `E=2c*eps`, the identity reads `2X=Q+E`. -/
theorem central_coordinate_from_rounding_error
    (c Q k P eps X : ℚ)
    (hX : X = c*P-k*Q)
    (heps : c*eps = c*P-Q*(k+1/2)) :
    2*X = Q+2*c*eps := by
  calc
    2*X = 2*(c*P-k*Q) := by rw [hX]
    _ = Q+2*(c*P-Q*(k+1/2)) := by ring
    _ = Q+2*(c*eps) := by rw [heps]
    _ = Q+2*c*eps := by ring

/-- The centered residue orbit has an exact cleared-denominator edge step. -/
theorem centered_residue_edge_step
    (c dj Qi Qj eps epsj beta delta ell X E Ej : ℚ)
    (hQj : Qj = c+dj)
    (hE : E = 2*c*eps) (hEj : Ej = 2*c*epsj)
    (htransport : Qi*epsj = dj*eps+Qi*(delta-1/2-beta))
    (hbeta : Qi*beta = ell-X)
    (hX : X = Qi/2+c*eps) :
    Qi*Ej-Qj*E = 2*c*(Qi*delta-ell) := by
  calc
    Qi*Ej-Qj*E = 2*c*(Qi*epsj-dj*eps-c*eps) := by
      rw [hE, hEj, hQj]
      ring
    _ = 2*c*(Qi*delta-Qi/2-Qi*beta-c*eps) := by
      rw [htransport]
      ring
    _ = 2*c*(Qi*delta-Qi/2-(ell-X)-c*eps) := by rw [hbeta]
    _ = 2*c*(Qi*delta-ell) := by
      rw [hX]
      ring

/-- Dividing rounding error by its runner speed makes every functional-tail
edge contract by the source factor `d_i/(c+d_i)`. -/
theorem normalized_tail_transport
    (di dj Qi eps epsj z zj b : ℚ)
    (hdi0 : di ≠ 0) (hdj0 : dj ≠ 0) (hQi0 : Qi ≠ 0)
    (hz : z = eps/di) (hzj : zj = epsj/dj)
    (htransport : epsj = (dj/Qi)*eps+b) :
    zj = (di/Qi)*z+b/dj := by
  rw [hz, hzj, htransport]
  field_simp [hdi0, hdj0, hQi0]

/-- The chosen wall is an exact rational function of the target centered
residue, relative digit, and side sign. -/
theorem wall_coordinate_reconstruction
    (c dj Qj Pj k N delta sigma E : ℚ)
    (hdj0 : dj ≠ 0) (hQj0 : Qj ≠ 0)
    (hQj : Qj = c+dj)
    (hPj : Pj = k+N+delta)
    (hE : E = 2*c*Pj-(2*k+1)*Qj) :
    14*dj*Qj*((14*N+sigma)/(14*dj)-Pj/Qj) =
      7*E+(7-14*delta+sigma)*Qj := by
  rw [hE, hPj]
  field_simp [hdj0, hQj0]
  rw [hQj]
  ring

/-- Arithmetic core of the sharpened target-wall clearance.  In the
application `y=σE_j`, `a=14|δ-1/2|-1≥6`, and
`C=14 d_j Q_j σ(t_j-w)`. -/
theorem target_wall_clearance_lower_bound
    (c dj Q y a C : ℝ)
    (hc : 0 < c) (hd : c < dj) (hy : y ≤ c) (ha : 6 ≤ a)
    (hQ : Q = c+dj) (hC : C = -7*y+a*Q) :
    5*dj < C := by
  rw [hC, hQ]
  nlinarith

/-- Once the target-wall leg is joined to the positive source-to-wall leg,
the cleared `> 5 d_j` estimate gives the sharp phase distance `> 5/14`.
This is the arithmetic bridge consumed by THM-1256. -/
theorem target_wall_clearance_forces_phase_distance
    {dj Q C targetLeg sourceLeg phaseDist : ℝ}
    (hdj : 0 < dj) (hQ : 0 < Q)
    (hC : 5*dj < C)
    (hCleg : C = 14*dj*Q*targetLeg)
    (hsource : 0 < sourceLeg)
    (hdist : |phaseDist| = Q*(sourceLeg+targetLeg)) :
    5/14 < |phaseDist| := by
  rw [hCleg] at hC
  have hscaled : dj*5 < dj*(14*Q*targetLeg) := by
    simpa [mul_assoc, mul_comm, mul_left_comm] using hC
  have hleg : 5 < 14*Q*targetLeg :=
    lt_of_mul_lt_mul_left hscaled hdj.le
  have hsource' : 0 < Q*sourceLeg := mul_pos hQ hsource
  rw [hdist]
  nlinarith

/-- The source runner's residue at the target spoke is reconstructed from
the determinant and the target carrier coordinate. -/
theorem reciprocal_blocker_residue
    (Pi Pj Qi Qj c di k D Xj : ℤ)
    (hD : D = Pi*Qj-Pj*Qi)
    (hXj : Xj = c*Pj-k*Qj)
    (hQi : Qi = c+di) :
    D+Xj = Qj*(Pi-k)-di*Pj := by
  rw [hD, hXj, hQi]
  ring

/-- A loopless functional lasso on six vertices, whose cycle has at least
two vertices, has at most four proper tail edges. -/
theorem functional_tail_length_at_most_four
    (tail cycle : ℕ) (hcycle : 2 ≤ cycle) (htotal : tail+cycle ≤ 6) :
    tail ≤ 4 := by
  omega

/-- THM-1198's `d₁/c<13/6` makes every nonempty normalized tail contract by
strictly more than the universal factor `13/19`. -/
theorem slowest_tail_multiplier_cap
    (c d R : ℝ)
    (hc : 0 < c) (hd : 0 < d)
    (hR : R ≤ d/(c+d)) (hratio : 6*d < 13*c) :
    R < 13/19 := by
  have hden : 0 < c+d := by positivity
  have hfrac : d/(c+d) < (13 : ℝ)/19 := by
    rw [div_lt_iff₀ hden]
    nlinarith
  exact lt_of_le_of_lt hR hfrac

/-- Two distinct located seam pairs, hence a two-edge forest, contribute the
exact `49/6` positioned Hunter coefficient on the full carrier gap. -/
theorem two_seam_hunter_rearrangement
    (G H q c : ℝ)
    (hc : 0 < c)
    (hG : G = 6/(7*c))
    (hhunter : G/7+q ≤ 6*H/49) :
    1/c+(49/6)*q ≤ H := by
  have hc0 : c ≠ 0 := ne_of_gt hc
  have hscale : (49/6)*(G/7) = 1/c := by
    rw [hG]
    field_simp [hc0]
    ring
  calc
    1/c+(49/6)*q = (49/6)*(G/7+q) := by
      rw [mul_add, hscale]
    _ ≤ (49/6)*(6*H/49) :=
      mul_le_mul_of_nonneg_left hhunter (by norm_num)
    _ = H := by ring

/-- Numeric shadow of the owner-reuse pigeonhole: two ascent walls cannot
inject into the outside labels once the lasso already uses five labels. -/
theorem two_ascent_owner_pigeonhole
    (ascents orbitVertices : ℕ)
    (ha : 2 ≤ ascents) (hv : 5 ≤ orbitVertices)
    (hinj : ascents ≤ 6-orbitVertices) : False := by
  omega

/-- If both walls of every one of `v` blocker teeth choose an owner, absence
of lasso incidence and absence of repeated outside ownership would inject
`2v` wall occurrences into only `6-v` outside labels.  This is impossible as
soon as the lasso has at least three vertices. -/
theorem two_wall_owner_pigeonhole
    (v : ℕ) (hv : 3 ≤ v) (hinj : 2*v ≤ 6-v) : False := by
  omega

/-- Every positive difference of two fast-tooth endpoints has cleared
numerator at least `gcd(j,h)`.  The paper layer divides by `14jh`; unlike the
stronger abutment quantum, this universal interior-overlap quantum needs no
mod-14 compatibility hypothesis. -/
theorem positive_endpoint_numerator_gcd_quantum
    (j h A B N : ℤ) (hj : 0 < j)
    (hN : N = h*A-j*B) (hpos : 0 < N) :
    (Int.gcd j h : ℤ) ≤ N := by
  have hgj : (Int.gcd j h : ℤ) ∣ j := Int.gcd_dvd_left j h
  have hgh : (Int.gcd j h : ℤ) ∣ h := Int.gcd_dvd_right j h
  have hdiv : (Int.gcd j h : ℤ) ∣ N := by
    rw [hN]
    exact dvd_sub (dvd_mul_of_dvd_left hgh A) (dvd_mul_of_dvd_left hgj B)
  have hgpos : 0 < (Int.gcd j h : ℤ) := by
    exact_mod_cast Int.gcd_pos_of_ne_zero_left h hj.ne'
  exact Int.le_of_dvd hpos hdiv

/-- The selected target-tooth wall has a quotient-visible source residue.
After the paper layer applies circle distance, this decides whether the
source remains safe at that wall, including on admissible descent edges. -/
theorem selected_wall_source_residue
    (di dj Q E r sigma A : ℚ)
    (hdj0 : dj ≠ 0) (hQ0 : Q ≠ 0)
    (hA : A = 7*dj*(Q-E)+di*(sigma*Q-14*r)) :
    (1/2-E/(2*Q))+(di/dj)*(sigma/14-r/Q) = A/(14*dj*Q) := by
  rw [hA]
  field_simp [hdj0, hQ0]
  ring

/-- Once all six centered residues are retained, a candidate owner's exact
residue at a selected wall is gauge-invariant and needs no absolute tooth
address.  Strict ownership is the paper predicate that the circle norm of
this rational number is below `1/14`. -/
theorem sampled_owner_residue_reconstruction
    (c dj Q Eh h Ei r sigma B L : ℚ)
    (hc0 : c ≠ 0) (hdj0 : dj ≠ 0) (hQ0 : Q ≠ 0)
    (hB : B = 7*dj*Q*(c-Eh)+7*h*Ei*dj-14*c*h*r+c*h*sigma*Q)
    (hL : L = 14*c*dj*Q) :
    B/L = 1/2-Eh/(2*c)+h*Ei/(2*c*Q)-h*r/(dj*Q)+h*sigma/(14*dj) := by
  rw [hB, hL]
  field_simp [hc0, hdj0, hQ0]
  ring

/-- The forward and backward signed margins of an owning tooth add to its
full danger-tooth width. -/
theorem oriented_owner_germ_width
    (h alpha sigma forward backward : ℚ)
    (hh0 : h ≠ 0)
    (hf : forward = (1/14-sigma*alpha)/h)
    (hb : backward = (1/14+sigma*alpha)/h) :
    forward+backward = 1/(7*h) := by
  rw [hf, hb]
  field_simp [hh0]
  ring

/-- Combining the protected rank-two debt with the fully located tree debt
charges the three remaining tree edges only after the protected harmonic
surplus has been exhausted. -/
theorem protected_tree_tariff_combination
    (S B wF wR : ℝ)
    (hprotected : B+(7/12)*wF ≤ S)
    (htree : 1+(7/12)*(wF+wR) ≤ S) :
    B+(7/12)*wF+max 0 ((7/12)*wR-(B-1)) ≤ S := by
  by_cases h : (7/12)*wR-(B-1) ≤ 0
  · rw [max_eq_left h]
    linarith
  · have hp : 0 ≤ (7/12)*wR-(B-1) := le_of_lt (lt_of_not_ge h)
    rw [max_eq_right hp]
    linarith

/-- The protected-component base at the THM-1198 slowest-ratio cap. -/
theorem protected_base_at_ratio_cap :
    (9*((13 : ℚ)/6)+2)/(3*((13 : ℚ)/6)*(1+(13 : ℚ)/6)) = 258/247 := by
  norm_num

/-- Division-free identity behind the strict `258/247` protected base. -/
theorem protected_base_surplus_identity
    (x : ℝ) (hx0 : x ≠ 0) (hx1 : x+1 ≠ 0) :
    (9*x+2)/(3*x*(1+x))-258/247 =
      ((13-6*x)*(129*x+38))/(741*x*(x+1)) := by
  have hx1' : 1+x ≠ 0 := by simpa [add_comm] using hx1
  field_simp [hx0, hx1, hx1']
  ring

theorem protected_base_strict_surplus
    (x : ℝ) (hx : 0 < x) (hcap : 6*x < 13) :
    (258 : ℝ)/247 < (9*x+2)/(3*x*(1+x)) := by
  have hx0 : x ≠ 0 := ne_of_gt hx
  have hx1 : x+1 ≠ 0 := by positivity
  have hid := protected_base_surplus_identity x hx0 hx1
  have hrhs : 0 < ((13-6*x)*(129*x+38))/(741*x*(x+1)) := by positivity
  linarith

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
  have hmul' : 0 < Q1*(P1*P2*P3-M1*M2*M3) := by
    simpa [mul_comm] using hmul
  have hdiff : 0 < P1*P2*P3-M1*M2*M3 :=
    pos_of_mul_pos_right hmul' (le_of_lt hQ1)
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
#print axioms relative_digit_selects_germ
#print axioms oriented_germ_imbalance
#print axioms central_coordinate_from_rounding_error
#print axioms centered_residue_edge_step
#print axioms normalized_tail_transport
#print axioms wall_coordinate_reconstruction
#print axioms target_wall_clearance_lower_bound
#print axioms target_wall_clearance_forces_phase_distance
#print axioms reciprocal_blocker_residue
#print axioms functional_tail_length_at_most_four
#print axioms slowest_tail_multiplier_cap
#print axioms two_seam_hunter_rearrangement
#print axioms two_ascent_owner_pigeonhole
#print axioms two_wall_owner_pigeonhole
#print axioms positive_endpoint_numerator_gcd_quantum
#print axioms selected_wall_source_residue
#print axioms sampled_owner_residue_reconstruction
#print axioms oriented_owner_germ_width
#print axioms protected_tree_tariff_combination
#print axioms protected_base_at_ratio_cap
#print axioms protected_base_surplus_identity
#print axioms protected_base_strict_surplus
#print axioms affine_triangle_transport
#print axioms contraction_gap_constant
#print axioms positive_central_triangle_holonomy
#print axioms wall_event_source_margin

end LRC14.CenteredBlockerAddressCompression
