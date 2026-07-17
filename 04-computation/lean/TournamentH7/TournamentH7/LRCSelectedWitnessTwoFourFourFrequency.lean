import TournamentH7.LRCEndgameTwoEight

namespace LonelyRunner
namespace LRC14Grand

open scoped Classical

/-- The integer scalar frequency carried by the complete binary prefix code
`{1/2,1/4,1/4}` after the three reduced numerators are normalized to residue
one modulo four. -/
def twoFourFourPhaseFrequency (a₂ a₄a a₄b : ℤ) : ℤ :=
  (2 * a₂ + a₄a + a₄b) / 4

/-- The zero-frequency q244 residual is genuine at the fixed-phase interface.
The normalized numerators `(-7,5,9)` are all one modulo four and have frequency
zero, while the corresponding speed rows `(-14,5,9)` cover all four branches
at `u=5/12`.  Thus residue arithmetic alone cannot discard the exact signed
support-three relation.

The natural vertices here are the four branch residue classes: the q-two row
owns `{1,3}` and the q-four rows own `{2}` and `{0}`.  Orienting this exact
partition as a runner tournament would erase the cover predicate. -/
theorem twoFourFour_zero_frequency_fixedPhase_counterexample :
    twoFourFourPhaseFrequency (-7) 5 9 = 0 ∧
      ¬ HasThreeDetunedGoodBranch (-14) 5 9 4 ((5 : ℝ) / 12) := by
  constructor
  · norm_num [twoFourFourPhaseFrequency]
  · rintro ⟨c, hcIco, hcTwo, hcFourA, hcFourB⟩
    have hcBounds := Finset.mem_Ico.mp hcIco
    have hc0 : 0 ≤ c := hcBounds.1
    have hc4 : c < 4 := hcBounds.2
    interval_cases c
    · apply hcFourB
      rw [detunedBadBranches, Finset.mem_filter]
      refine ⟨by norm_num, 1, ?_⟩
      norm_num [abs_of_nonneg, abs_of_nonpos]
    · apply hcTwo
      rw [detunedBadBranches, Finset.mem_filter]
      refine ⟨by norm_num, -5, ?_⟩
      norm_num [abs_of_nonneg, abs_of_nonpos]
    · apply hcFourA
      rw [detunedBadBranches, Finset.mem_filter]
      refine ⟨by norm_num, 3, ?_⟩
      norm_num [abs_of_nonneg, abs_of_nonpos]
    · apply hcTwo
      rw [detunedBadBranches, Finset.mem_filter]
      refine ⟨by norm_num, -12, ?_⟩
      norm_num [abs_of_nonneg, abs_of_nonpos]

/-- The exact q244 residue cover has zero cyclic first moment. -/
theorem twoFourFour_code_offset_dvd
    (c₂ c₄a c₄b : ℤ)
    (hparityA : ¬ c₂ ≡ c₄a [ZMOD 2])
    (hparityB : ¬ c₂ ≡ c₄b [ZMOD 2])
    (hdistinct : ¬ c₄a ≡ c₄b [ZMOD 4]) :
    (4 : ℤ) ∣ 2 * c₂ + c₄a + c₄b := by
  change ¬ c₂ % 2 = c₄a % 2 at hparityA
  change ¬ c₂ % 2 = c₄b % 2 at hparityB
  change ¬ c₄a % 4 = c₄b % 4 at hdistinct
  rw [Int.dvd_iff_emod_eq_zero]
  have hc₂0 := Int.emod_nonneg c₂ (by norm_num : (2 : ℤ) ≠ 0)
  have hc₂2 := Int.emod_lt_of_pos c₂ (by norm_num : (0 : ℤ) < 2)
  have hca0 := Int.emod_nonneg c₄a (by norm_num : (4 : ℤ) ≠ 0)
  have hca4 := Int.emod_lt_of_pos c₄a (by norm_num : (0 : ℤ) < 4)
  have hcb0 := Int.emod_nonneg c₄b (by norm_num : (4 : ℤ) ≠ 0)
  have hcb4 := Int.emod_lt_of_pos c₄b (by norm_num : (0 : ℤ) < 4)
  have hcaParity : c₄a % 4 % 2 = c₄a % 2 :=
    Int.emod_emod_of_dvd c₄a (by norm_num : (2 : ℤ) ∣ 4)
  have hcbParity : c₄b % 4 % 2 = c₄b % 2 :=
    Int.emod_emod_of_dvd c₄b (by norm_num : (2 : ℤ) ∣ 4)
  omega

/-- Remove the common divisor from one bad row and absorb an independent sign
into its integer witness. -/
theorem normalizedPhaseBad_of_detunedBadBranch
    (δ g d q ε a c : ℤ) (u : ℝ)
    (hd : d ≠ 0) (hq : q ≠ 0)
    (hδ : δ = d * (ε * a)) (hg : g = d * q)
    (hε : ε = 1 ∨ ε = -1)
    (hbad : c ∈ detunedBadBranches δ g u) :
    ∃ n : ℤ,
      |(a : ℝ) * ((u + (c : ℝ)) / (q : ℝ)) - n| < 1 / 14 := by
  obtain ⟨n, hn⟩ := (Finset.mem_filter.mp hbad).2
  have hdR : (d : ℝ) ≠ 0 := by exact_mod_cast hd
  have hqR : (q : ℝ) ≠ 0 := by exact_mod_cast hq
  rcases hε with rfl | rfl
  · refine ⟨n, ?_⟩
    have heq :
        (δ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n =
          (a : ℝ) * ((u + (c : ℝ)) / (q : ℝ)) - n := by
      rw [hδ, hg]
      push_cast
      field_simp [hdR, hqR]
    rwa [heq] at hn
  · refine ⟨-n, ?_⟩
    have heq :
        (δ : ℝ) * ((u + (c : ℝ)) / (g : ℝ)) - n =
          -((a : ℝ) * ((u + (c : ℝ)) / (q : ℝ)) - (-n : ℤ)) := by
      rw [hδ, hg]
      push_cast
      field_simp [hdR, hqR]
      ring
    rw [heq, abs_neg] at hn
    exact hn

/-- A bad normalized denominator-`q` row is contained in one integer-frequency
wall of radius `q/14`.  Multiplication by `q` absorbs the branch offset into the
integer witness. -/
theorem integerFrequencyBad_of_normalizedBad
    (numerator branch denominator : ℤ) (u : ℝ)
    (hdenominator : 0 < denominator)
    (hbad : ∃ integer : ℤ,
      |(numerator : ℝ) *
          ((u + (branch : ℝ)) / (denominator : ℝ)) - integer| < 1 / 14) :
    ∃ integer : ℤ,
      |(numerator : ℝ) * u - integer| < (denominator : ℝ) / 14 := by
  obtain ⟨integer, hinteger⟩ := hbad
  refine ⟨denominator * integer - numerator * branch, ?_⟩
  have hdenominatorReal : (0 : ℝ) < denominator := by
    exact_mod_cast hdenominator
  have heq :
      (numerator : ℝ) * u -
          ((denominator * integer - numerator * branch : ℤ) : ℝ) =
        (denominator : ℝ) *
          ((numerator : ℝ) *
            ((u + (branch : ℝ)) / (denominator : ℝ)) - integer) := by
    push_cast
    field_simp [ne_of_gt hdenominatorReal]
    ring
  rw [heq, abs_mul, abs_of_pos hdenominatorReal]
  calc
    (denominator : ℝ) *
        |(numerator : ℝ) *
          ((u + (branch : ℝ)) / (denominator : ℝ)) - integer| <
      (denominator : ℝ) * (1 / 14) :=
        mul_lt_mul_of_pos_left hinteger hdenominatorReal
    _ = (denominator : ℝ) / 14 := by ring

/-- The q-two specialization has radius `1/7`.  Unlike the combined q244 wall,
this scalar wall does not suffer cancellation between selected speeds. -/
theorem integerFrequencyBad_of_normalizedQTwoBad
    (numerator branch : ℤ) (u : ℝ)
    (hbad : ∃ integer : ℤ,
      |(numerator : ℝ) * ((u + (branch : ℝ)) / 2) - integer| < 1 / 14) :
    ∃ integer : ℤ,
      |(numerator : ℝ) * u - integer| < 1 / 7 := by
  convert integerFrequencyBad_of_normalizedBad
    numerator branch 2 u (by norm_num) hbad using 1 <;> norm_num

/-- A q2 or q4 reduced numerator can be independently signed to residue one
modulo four. -/
theorem exists_signed_reducedNumerator_modFour
    (δ g q : ℤ) (hg : 0 < g)
    (hqValue : q = 2 ∨ q = 4)
    (hreduced : g / (Int.gcd δ g : ℤ) = q) :
    ∃ d ε a : ℤ,
      d ≠ 0 ∧ δ = d * (ε * a) ∧ g = d * q ∧
      (ε = 1 ∨ ε = -1) ∧ (4 : ℤ) ∣ a - 1 := by
  let d : ℤ := (Int.gcd δ g : ℤ)
  let p : ℤ := δ / d
  have hdNat : 0 < Int.gcd δ g := by
    rw [Int.gcd_pos_iff]
    right
    omega
  have hd : d ≠ 0 := by
    dsimp [d]
    exact_mod_cast (ne_of_gt hdNat)
  have hδ : δ = d * p := (Int.mul_ediv_cancel' (Int.gcd_dvd_left δ g)).symm
  have hgFactor : g = d * q := by
    have hfactor := Int.mul_ediv_cancel' (Int.gcd_dvd_right δ g)
    change d * (g / d) = g at hfactor
    rw [hreduced] at hfactor
    exact hfactor.symm
  have hcoprime : Int.gcd p q = 1 := by
    have h := Int.gcd_div_gcd_div_gcd hdNat
    change Int.gcd p (g / d) = 1 at h
    rwa [hreduced] at h
  have hpOdd : ¬ (2 : ℤ) ∣ p := by
    intro hpEven
    have hqEven : (2 : ℤ) ∣ q := by
      rcases hqValue with rfl | rfl <;> norm_num
    have hbez : (1 : ℤ) =
        p * Int.gcdA p q + q * Int.gcdB p q := by
      rw [← Int.gcd_eq_gcd_ab, hcoprime]
      norm_num
    have honeEven : (2 : ℤ) ∣ 1 := by
      rw [hbez]
      exact dvd_add (dvd_mul_of_dvd_left hpEven _) (dvd_mul_of_dvd_left hqEven _)
    norm_num at honeEven
  have hpRem : p % 4 = 1 ∨ p % 4 = 3 := by
    have hpRem0 := Int.emod_nonneg p (by norm_num : (4 : ℤ) ≠ 0)
    have hpRem4 := Int.emod_lt_of_pos p (by norm_num : (0 : ℤ) < 4)
    have hpParity : p % 4 % 2 = p % 2 :=
      Int.emod_emod_of_dvd p (by norm_num : (2 : ℤ) ∣ 4)
    have hpModTwo : p % 2 ≠ 0 := by
      intro hzero
      exact hpOdd (Int.dvd_iff_emod_eq_zero.mpr hzero)
    omega
  rcases hpRem with hpOne | hpThree
  · refine ⟨d, 1, p, hd, ?_, hgFactor, Or.inl rfl, ?_⟩
    · rw [hδ]
      ring
    · rw [Int.dvd_iff_emod_eq_zero]
      omega
  · refine ⟨d, -1, -p, hd, ?_, hgFactor, Or.inr rfl, ?_⟩
    · rw [hδ]
      ring
    · have hpDecomp := Int.emod_add_mul_ediv p 4
      rw [hpThree] at hpDecomp
      refine ⟨-(p / 4) - 1, ?_⟩
      omega

/-- Absorbing a sign into a reduced numerator gives an exact scaling identity.
This is the algebraic bridge from the normalized q244 frequency back to the
three original selected speeds. -/
theorem signed_reducedNumerator_scale
    (delta modulus divisor denominator sign numerator : ℤ)
    (hdelta : delta = divisor * (sign * numerator))
    (hmodulus : modulus = divisor * denominator)
    (hsign : sign = 1 ∨ sign = -1) :
    denominator * (sign * delta) = modulus * numerator := by
  rcases hsign with rfl | rfl
  · rw [hdelta, hmodulus]
    ring
  · rw [hdelta, hmodulus]
    ring

/-- Absolute-value form of a normalized denominator scale. -/
theorem normalized_abs_scale
    (delta modulus denominator sign numerator : ℤ)
    (hmodulus : 0 ≤ modulus)
    (hdenominator : 0 ≤ denominator)
    (hsign : sign = 1 ∨ sign = -1)
    (hscale : denominator * (sign * delta) = modulus * numerator) :
    (denominator : ℝ) * |(delta : ℝ)| =
      (modulus : ℝ) * |(numerator : ℝ)| := by
  have hscaleReal :
      (denominator : ℝ) * ((sign : ℝ) * (delta : ℝ)) =
        (modulus : ℝ) * (numerator : ℝ) := by
    exact_mod_cast hscale
  have hsignAbs : |(sign : ℝ)| = 1 := by
    rcases hsign with rfl | rfl <;> norm_num
  have hmodulusReal : (0 : ℝ) ≤ modulus := by exact_mod_cast hmodulus
  have hdenominatorReal : (0 : ℝ) ≤ denominator := by
    exact_mod_cast hdenominator
  have habs := congrArg abs hscaleReal
  rw [abs_mul, abs_mul, abs_mul, hsignAbs,
    abs_of_nonneg hmodulusReal, abs_of_nonneg hdenominatorReal] at habs
  norm_num only [one_mul] at habs
  exact habs

/-- Absolute-value form of the normalized q-two scale. -/
theorem normalizedQTwo_abs_scale
    (delta modulus sign numerator : ℤ)
    (hmodulus : 0 ≤ modulus)
    (hsign : sign = 1 ∨ sign = -1)
    (hscale : 2 * (sign * delta) = modulus * numerator) :
    (2 : ℝ) * |(delta : ℝ)| =
      (modulus : ℝ) * |(numerator : ℝ)| :=
  normalized_abs_scale delta modulus 2 sign numerator
    hmodulus (by norm_num) hsign hscale

/-- Exact q244 frequency dictionary.  Once the three primitive numerators are
signed to residue one modulo four, the scalar wall frequency is not an
auxiliary parameter: multiplying it by the common modulus gives the signed
sum of the three selected speeds. -/
theorem signed_sum_eq_modulus_mul_twoFourFourPhaseFrequency
    (deltaTwo deltaFourA deltaFourB modulus : ℤ)
    (signTwo signFourA signFourB
      numeratorTwo numeratorFourA numeratorFourB : ℤ)
    (hscaleTwo :
      2 * (signTwo * deltaTwo) = modulus * numeratorTwo)
    (hscaleFourA :
      4 * (signFourA * deltaFourA) = modulus * numeratorFourA)
    (hscaleFourB :
      4 * (signFourB * deltaFourB) = modulus * numeratorFourB)
    (hresTwo : (4 : ℤ) ∣ numeratorTwo - 1)
    (hresFourA : (4 : ℤ) ∣ numeratorFourA - 1)
    (hresFourB : (4 : ℤ) ∣ numeratorFourB - 1) :
    signTwo * deltaTwo + signFourA * deltaFourA +
        signFourB * deltaFourB =
      modulus * twoFourFourPhaseFrequency
        numeratorTwo numeratorFourA numeratorFourB := by
  obtain ⟨kTwo, hkTwo⟩ := hresTwo
  obtain ⟨kFourA, hkFourA⟩ := hresFourA
  obtain ⟨kFourB, hkFourB⟩ := hresFourB
  have hfrequencyDvd :
      (4 : ℤ) ∣ 2 * numeratorTwo + numeratorFourA + numeratorFourB := by
    refine ⟨2 * kTwo + kFourA + kFourB + 1, ?_⟩
    omega
  have hfrequency :
      twoFourFourPhaseFrequency numeratorTwo numeratorFourA numeratorFourB * 4 =
        2 * numeratorTwo + numeratorFourA + numeratorFourB :=
    Int.ediv_mul_cancel hfrequencyDvd
  have hfour :
      4 * (signTwo * deltaTwo + signFourA * deltaFourA +
          signFourB * deltaFourB) =
        4 * (modulus * twoFourFourPhaseFrequency
          numeratorTwo numeratorFourA numeratorFourB) := by
    calc
      4 * (signTwo * deltaTwo + signFourA * deltaFourA +
          signFourB * deltaFourB) =
        2 * (2 * (signTwo * deltaTwo)) +
          4 * (signFourA * deltaFourA) +
          4 * (signFourB * deltaFourB) := by ring
      _ = 2 * (modulus * numeratorTwo) +
          modulus * numeratorFourA + modulus * numeratorFourB := by
            rw [hscaleTwo, hscaleFourA, hscaleFourB]
      _ = modulus *
          (2 * numeratorTwo + numeratorFourA + numeratorFourB) := by ring
      _ = modulus *
          (twoFourFourPhaseFrequency
            numeratorTwo numeratorFourA numeratorFourB * 4) := by
              rw [hfrequency]
      _ = 4 * (modulus * twoFourFourPhaseFrequency
          numeratorTwo numeratorFourA numeratorFourB) := by ring
  omega

/-- The q244 parallel partition forces its normalized scalar frequency into
the open `3/14` neighborhood of an integer.  `hcode` is exactly the offset
congruence of the complete residue code: `2*c₂+c₄a+c₄b = 0 (mod 4)`. -/
theorem frequency_bad_of_twoFourFour_matching
    (a₂ a₄a a₄b c₂ c₄a c₄b : ℤ) (u : ℝ)
    (hres₂ : (4 : ℤ) ∣ a₂ - 1)
    (hres₄a : (4 : ℤ) ∣ a₄a - 1)
    (hres₄b : (4 : ℤ) ∣ a₄b - 1)
    (hcode : (4 : ℤ) ∣ 2 * c₂ + c₄a + c₄b)
    (hbad₂ : ∃ n : ℤ,
      |(a₂ : ℝ) * ((u + (c₂ : ℝ)) / 2) - n| < 1 / 14)
    (hbad₄a : ∃ n : ℤ,
      |(a₄a : ℝ) * ((u + (c₄a : ℝ)) / 4) - n| < 1 / 14)
    (hbad₄b : ∃ n : ℤ,
      |(a₄b : ℝ) * ((u + (c₄b : ℝ)) / 4) - n| < 1 / 14) :
    ∃ n : ℤ,
      |(twoFourFourPhaseFrequency a₂ a₄a a₄b : ℝ) * u - n| < 3 / 14 := by
  obtain ⟨n₂, hn₂⟩ := hbad₂
  obtain ⟨n₄a, hn₄a⟩ := hbad₄a
  obtain ⟨n₄b, hn₄b⟩ := hbad₄b
  obtain ⟨k₂, hk₂⟩ := hres₂
  obtain ⟨k₄a, hk₄a⟩ := hres₄a
  obtain ⟨k₄b, hk₄b⟩ := hres₄b
  obtain ⟨kcode, hkcode⟩ := hcode
  have hfrequencyDvd : (4 : ℤ) ∣ 2 * a₂ + a₄a + a₄b := by
    refine ⟨2 * k₂ + k₄a + k₄b + 1, ?_⟩
    omega
  have hoffsetDvd :
      (4 : ℤ) ∣ 2 * a₂ * c₂ + a₄a * c₄a + a₄b * c₄b := by
    refine ⟨kcode + 2 * k₂ * c₂ + k₄a * c₄a + k₄b * c₄b, ?_⟩
    have ha₂ : a₂ = 4 * k₂ + 1 := by omega
    have ha₄a : a₄a = 4 * k₄a + 1 := by omega
    have ha₄b : a₄b = 4 * k₄b + 1 := by omega
    have hcode' : 2 * c₂ + c₄a + c₄b = 4 * kcode := hkcode
    rw [ha₂, ha₄a, ha₄b]
    calc
      2 * (4 * k₂ + 1) * c₂ + (4 * k₄a + 1) * c₄a +
          (4 * k₄b + 1) * c₄b =
        4 * (2 * k₂ * c₂ + k₄a * c₄a + k₄b * c₄b) +
          (2 * c₂ + c₄a + c₄b) := by ring
      _ = 4 * (kcode + 2 * k₂ * c₂ + k₄a * c₄a + k₄b * c₄b) := by
        rw [hcode']
        ring
  let F : ℤ := (2 * a₂ + a₄a + a₄b) / 4
  let K : ℤ :=
    (2 * a₂ * c₂ + a₄a * c₄a + a₄b * c₄b) / 4
  have hF : F * 4 = 2 * a₂ + a₄a + a₄b :=
    Int.ediv_mul_cancel hfrequencyDvd
  have hK : K * 4 =
      2 * a₂ * c₂ + a₄a * c₄a + a₄b * c₄b :=
    Int.ediv_mul_cancel hoffsetDvd
  have hfrequency : twoFourFourPhaseFrequency a₂ a₄a a₄b = F := by
    rfl
  refine ⟨n₂ + n₄a + n₄b - K, ?_⟩
  have heq :
      (F : ℝ) * u - ((n₂ + n₄a + n₄b - K : ℤ) : ℝ) =
        ((a₂ : ℝ) * ((u + (c₂ : ℝ)) / 2) - n₂) +
        ((a₄a : ℝ) * ((u + (c₄a : ℝ)) / 4) - n₄a) +
        ((a₄b : ℝ) * ((u + (c₄b : ℝ)) / 4) - n₄b) := by
    push_cast
    have hFR : (F : ℝ) * 4 =
        2 * (a₂ : ℝ) + a₄a + a₄b := by exact_mod_cast hF
    have hKR : (K : ℝ) * 4 =
        2 * (a₂ : ℝ) * c₂ + a₄a * c₄a + a₄b * c₄b := by
      exact_mod_cast hK
    have hFdiv : (F : ℝ) =
        (2 * (a₂ : ℝ) + a₄a + a₄b) / 4 := by linarith
    have hKdiv : (K : ℝ) =
        (2 * (a₂ : ℝ) * c₂ + a₄a * c₄a + a₄b * c₄b) / 4 := by
      linarith
    rw [hFdiv, hKdiv]
    ring
  rw [hfrequency, heq]
  calc
    |((a₂ : ℝ) * ((u + (c₂ : ℝ)) / 2) - n₂) +
        ((a₄a : ℝ) * ((u + (c₄a : ℝ)) / 4) - n₄a) +
        ((a₄b : ℝ) * ((u + (c₄b : ℝ)) / 4) - n₄b)|
        ≤ |(a₂ : ℝ) * ((u + (c₂ : ℝ)) / 2) - n₂| +
          |(a₄a : ℝ) * ((u + (c₄a : ℝ)) / 4) - n₄a| +
          |(a₄b : ℝ) * ((u + (c₄b : ℝ)) / 4) - n₄b| := by
            nlinarith [
              abs_add_le
                ((a₂ : ℝ) * ((u + (c₂ : ℝ)) / 2) - n₂)
                ((a₄a : ℝ) * ((u + (c₄a : ℝ)) / 4) - n₄a),
              abs_add_le
                (((a₂ : ℝ) * ((u + (c₂ : ℝ)) / 2) - n₂) +
                  ((a₄a : ℝ) * ((u + (c₄a : ℝ)) / 4) - n₄a))
                ((a₄b : ℝ) * ((u + (c₄b : ℝ)) / 4) - n₄b)]
    _ < 3 / 14 := by linarith

/-- Concrete q244 failure wall.  The only extra data beyond the existing
denominator hypotheses are signed reduced normal forms whose numerators are
all `1 mod 4`; those signs always exist for the primitive q2/q4 numerators. -/
theorem frequency_bad_of_twoFourFour_failure
    (δ₂ δ₄a δ₄b g : ℤ) (u : ℝ) (hg : 2 ≤ g)
    (hq2 : g / (Int.gcd δ₂ g : ℤ) = 2)
    (hq4a : g / (Int.gcd δ₄a g : ℤ) = 4)
    (hq4b : g / (Int.gcd δ₄b g : ℤ) = 4)
    (d₂ d₄a d₄b ε₂ ε₄a ε₄b a₂ a₄a a₄b : ℤ)
    (hd₂ : d₂ ≠ 0) (hd₄a : d₄a ≠ 0) (hd₄b : d₄b ≠ 0)
    (hδ₂ : δ₂ = d₂ * (ε₂ * a₂))
    (hδ₄a : δ₄a = d₄a * (ε₄a * a₄a))
    (hδ₄b : δ₄b = d₄b * (ε₄b * a₄b))
    (hg₂ : g = d₂ * 2) (hg₄a : g = d₄a * 4)
    (hg₄b : g = d₄b * 4)
    (hε₂ : ε₂ = 1 ∨ ε₂ = -1)
    (hε₄a : ε₄a = 1 ∨ ε₄a = -1)
    (hε₄b : ε₄b = 1 ∨ ε₄b = -1)
    (hres₂ : (4 : ℤ) ∣ a₂ - 1)
    (hres₄a : (4 : ℤ) ∣ a₄a - 1)
    (hres₄b : (4 : ℤ) ∣ a₄b - 1)
    (hfail : ¬ HasThreeDetunedGoodBranch δ₂ δ₄a δ₄b g u) :
    ∃ n : ℤ,
      |(twoFourFourPhaseFrequency a₂ a₄a a₄b : ℝ) * u - n| < 3 / 14 := by
  obtain ⟨c₂, c₄a, c₄b, hc₂, hc₄a, hc₄b, -, -, -,
      hparityA, hparityB, hdistinct⟩ :=
    qTwo_four_four_failure_normal_form
      δ₂ δ₄a δ₄b g u hg hq2 hq4a hq4b hfail
  have hbad₂ := normalizedPhaseBad_of_detunedBadBranch
    δ₂ g d₂ 2 ε₂ a₂ c₂ u hd₂ (by norm_num) hδ₂ hg₂ hε₂ hc₂
  have hbad₄a := normalizedPhaseBad_of_detunedBadBranch
    δ₄a g d₄a 4 ε₄a a₄a c₄a u hd₄a (by norm_num)
      hδ₄a hg₄a hε₄a hc₄a
  have hbad₄b := normalizedPhaseBad_of_detunedBadBranch
    δ₄b g d₄b 4 ε₄b a₄b c₄b u hd₄b (by norm_num)
      hδ₄b hg₄b hε₄b hc₄b
  exact frequency_bad_of_twoFourFour_matching
    a₂ a₄a a₄b c₂ c₄a c₄b u hres₂ hres₄a hres₄b
      (twoFourFour_code_offset_dvd c₂ c₄a c₄b
        hparityA hparityB hdistinct)
      hbad₂ hbad₄a hbad₄b

/-- Canonical q244 wall with its arithmetic content retained.  In addition to
the failure wall, this version exposes the three independent signs and the
exact signed-speed identity that the earlier existential interface discarded.
Consequently a zero wall frequency is literally a support-three relation with
coefficients in `{+1,-1}`. -/
theorem exists_twoFourFour_normalized_signed_failureWall
    (deltaTwo deltaFourA deltaFourB modulus : ℤ) (hmodulus : 2 ≤ modulus)
    (hqTwo : modulus / (Int.gcd deltaTwo modulus : ℤ) = 2)
    (hqFourA : modulus / (Int.gcd deltaFourA modulus : ℤ) = 4)
    (hqFourB : modulus / (Int.gcd deltaFourB modulus : ℤ) = 4) :
    ∃ signTwo signFourA signFourB
        numeratorTwo numeratorFourA numeratorFourB : ℤ,
      (signTwo = 1 ∨ signTwo = -1) ∧
      (signFourA = 1 ∨ signFourA = -1) ∧
      (signFourB = 1 ∨ signFourB = -1) ∧
      (4 : ℤ) ∣ numeratorTwo - 1 ∧
      (4 : ℤ) ∣ numeratorFourA - 1 ∧
      (4 : ℤ) ∣ numeratorFourB - 1 ∧
      2 * (signTwo * deltaTwo) = modulus * numeratorTwo ∧
      4 * (signFourA * deltaFourA) = modulus * numeratorFourA ∧
      4 * (signFourB * deltaFourB) = modulus * numeratorFourB ∧
      signTwo * deltaTwo + signFourA * deltaFourA +
          signFourB * deltaFourB =
        modulus * twoFourFourPhaseFrequency
          numeratorTwo numeratorFourA numeratorFourB ∧
      (∀ u : ℝ,
        ¬ HasThreeDetunedGoodBranch
          deltaTwo deltaFourA deltaFourB modulus u →
        ∃ integer : ℤ, |(numeratorTwo : ℝ) * u - integer| < 1 / 7) ∧
      (∀ u : ℝ,
        ¬ HasThreeDetunedGoodBranch
          deltaTwo deltaFourA deltaFourB modulus u →
        ∃ integer : ℤ, |(numeratorFourA : ℝ) * u - integer| < 2 / 7) ∧
      (∀ u : ℝ,
        ¬ HasThreeDetunedGoodBranch
          deltaTwo deltaFourA deltaFourB modulus u →
        ∃ integer : ℤ, |(numeratorFourB : ℝ) * u - integer| < 2 / 7) ∧
      ∀ u : ℝ,
        ¬ HasThreeDetunedGoodBranch
          deltaTwo deltaFourA deltaFourB modulus u →
        ∃ integer : ℤ,
          |(twoFourFourPhaseFrequency
              numeratorTwo numeratorFourA numeratorFourB : ℝ) * u - integer| <
            3 / 14 := by
  obtain ⟨divisorTwo, signTwo, numeratorTwo, hdivisorTwo, hdeltaTwo,
      hmodulusTwo, hsignTwo, hresTwo⟩ :=
    exists_signed_reducedNumerator_modFour deltaTwo modulus 2
      (by omega) (Or.inl rfl) hqTwo
  obtain ⟨divisorFourA, signFourA, numeratorFourA, hdivisorFourA,
      hdeltaFourA, hmodulusFourA, hsignFourA, hresFourA⟩ :=
    exists_signed_reducedNumerator_modFour deltaFourA modulus 4
      (by omega) (Or.inr rfl) hqFourA
  obtain ⟨divisorFourB, signFourB, numeratorFourB, hdivisorFourB,
      hdeltaFourB, hmodulusFourB, hsignFourB, hresFourB⟩ :=
    exists_signed_reducedNumerator_modFour deltaFourB modulus 4
      (by omega) (Or.inr rfl) hqFourB
  have hscaleTwo :
      2 * (signTwo * deltaTwo) = modulus * numeratorTwo :=
    signed_reducedNumerator_scale deltaTwo modulus divisorTwo 2
      signTwo numeratorTwo hdeltaTwo hmodulusTwo hsignTwo
  have hscaleFourA :
      4 * (signFourA * deltaFourA) = modulus * numeratorFourA :=
    signed_reducedNumerator_scale deltaFourA modulus divisorFourA 4
      signFourA numeratorFourA hdeltaFourA hmodulusFourA hsignFourA
  have hscaleFourB :
      4 * (signFourB * deltaFourB) = modulus * numeratorFourB :=
    signed_reducedNumerator_scale deltaFourB modulus divisorFourB 4
      signFourB numeratorFourB hdeltaFourB hmodulusFourB hsignFourB
  have hidentity :=
    signed_sum_eq_modulus_mul_twoFourFourPhaseFrequency
      deltaTwo deltaFourA deltaFourB modulus
      signTwo signFourA signFourB numeratorTwo numeratorFourA numeratorFourB
      hscaleTwo hscaleFourA hscaleFourB hresTwo hresFourA hresFourB
  have hrowWalls : ∀ u : ℝ,
      ¬ HasThreeDetunedGoodBranch
          deltaTwo deltaFourA deltaFourB modulus u →
      (∃ integer : ℤ, |(numeratorTwo : ℝ) * u - integer| < 1 / 7) ∧
      (∃ integer : ℤ, |(numeratorFourA : ℝ) * u - integer| < 2 / 7) ∧
      (∃ integer : ℤ, |(numeratorFourB : ℝ) * u - integer| < 2 / 7) := by
    intro u hfail
    obtain ⟨branchTwo, branchFourA, branchFourB,
        hbranchTwo, hbranchFourA, hbranchFourB, -, -, -, -, -, -⟩ :=
      qTwo_four_four_failure_normal_form
        deltaTwo deltaFourA deltaFourB modulus u hmodulus
          hqTwo hqFourA hqFourB hfail
    have hbadTwo := normalizedPhaseBad_of_detunedBadBranch
      deltaTwo modulus divisorTwo 2 signTwo numeratorTwo branchTwo u
        hdivisorTwo (by norm_num) hdeltaTwo hmodulusTwo hsignTwo hbranchTwo
    have hbadFourA := normalizedPhaseBad_of_detunedBadBranch
      deltaFourA modulus divisorFourA 4 signFourA numeratorFourA branchFourA u
        hdivisorFourA (by norm_num) hdeltaFourA hmodulusFourA
          hsignFourA hbranchFourA
    have hbadFourB := normalizedPhaseBad_of_detunedBadBranch
      deltaFourB modulus divisorFourB 4 signFourB numeratorFourB branchFourB u
        hdivisorFourB (by norm_num) hdeltaFourB hmodulusFourB
          hsignFourB hbranchFourB
    refine ⟨integerFrequencyBad_of_normalizedQTwoBad
      numeratorTwo branchTwo u hbadTwo, ?_, ?_⟩
    · convert integerFrequencyBad_of_normalizedBad
        numeratorFourA branchFourA 4 u (by norm_num) hbadFourA using 1 <;>
          norm_num
    · convert integerFrequencyBad_of_normalizedBad
        numeratorFourB branchFourB 4 u (by norm_num) hbadFourB using 1 <;>
          norm_num
  refine ⟨signTwo, signFourA, signFourB,
    numeratorTwo, numeratorFourA, numeratorFourB,
    hsignTwo, hsignFourA, hsignFourB,
    hresTwo, hresFourA, hresFourB,
    hscaleTwo, hscaleFourA, hscaleFourB, hidentity,
    fun u hfail => (hrowWalls u hfail).1,
    fun u hfail => (hrowWalls u hfail).2.1,
    fun u hfail => (hrowWalls u hfail).2.2, ?_⟩
  intro u hfail
  exact frequency_bad_of_twoFourFour_failure
    deltaTwo deltaFourA deltaFourB modulus u hmodulus hqTwo hqFourA hqFourB
      divisorTwo divisorFourA divisorFourB signTwo signFourA signFourB
      numeratorTwo numeratorFourA numeratorFourB
      hdivisorTwo hdivisorFourA hdivisorFourB
      hdeltaTwo hdeltaFourA hdeltaFourB
      hmodulusTwo hmodulusFourA hmodulusFourB
      hsignTwo hsignFourA hsignFourB hresTwo hresFourA hresFourB hfail

/-- Denominator hypotheses alone canonically supply some normalized q244
scalar wall.  The frequency may vanish; that exact support-three relation is
the residual not handled by interval escape. -/
theorem exists_twoFourFour_normalized_failureWall
    (δ₂ δ₄a δ₄b g : ℤ) (hg : 2 ≤ g)
    (hq2 : g / (Int.gcd δ₂ g : ℤ) = 2)
    (hq4a : g / (Int.gcd δ₄a g : ℤ) = 4)
    (hq4b : g / (Int.gcd δ₄b g : ℤ) = 4) :
    ∃ a₂ a₄a a₄b : ℤ,
      (4 : ℤ) ∣ a₂ - 1 ∧ (4 : ℤ) ∣ a₄a - 1 ∧
      (4 : ℤ) ∣ a₄b - 1 ∧
      ∀ u : ℝ,
        ¬ HasThreeDetunedGoodBranch δ₂ δ₄a δ₄b g u →
        ∃ n : ℤ,
          |(twoFourFourPhaseFrequency a₂ a₄a a₄b : ℝ) * u - n| < 3 / 14 := by
  obtain ⟨d₂, ε₂, a₂, hd₂, hδ₂, hg₂, hε₂, hres₂⟩ :=
    exists_signed_reducedNumerator_modFour δ₂ g 2 (by omega)
      (Or.inl rfl) hq2
  obtain ⟨d₄a, ε₄a, a₄a, hd₄a, hδ₄a, hg₄a, hε₄a, hres₄a⟩ :=
    exists_signed_reducedNumerator_modFour δ₄a g 4 (by omega)
      (Or.inr rfl) hq4a
  obtain ⟨d₄b, ε₄b, a₄b, hd₄b, hδ₄b, hg₄b, hε₄b, hres₄b⟩ :=
    exists_signed_reducedNumerator_modFour δ₄b g 4 (by omega)
      (Or.inr rfl) hq4b
  refine ⟨a₂, a₄a, a₄b, hres₂, hres₄a, hres₄b, ?_⟩
  intro u hfail
  exact frequency_bad_of_twoFourFour_failure
    δ₂ δ₄a δ₄b g u hg hq2 hq4a hq4b
      d₂ d₄a d₄b ε₂ ε₄a ε₄b a₂ a₄a a₄b
      hd₂ hd₄a hd₄b hδ₂ hδ₄a hδ₄b hg₂ hg₄a hg₄b
      hε₂ hε₄a hε₄b hres₂ hres₄a hres₄b hfail

end LRC14Grand
end LonelyRunner
