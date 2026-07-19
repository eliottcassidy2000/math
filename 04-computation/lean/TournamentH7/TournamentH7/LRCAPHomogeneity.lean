/-
# Twelve-speed arithmetic-progression homogeneity (THM-1171)

For the reduced progression

`g * (A + k * D)`,  `0 <= k < 12`,  `Nat.Coprime A D`,

multiplication by an inverse of `A` modulo `D` puts every phase on the same
half-residue.  When `D >= 2` this gives clearance at least `1/3`.  When
`D = 1`, the centered time gives clearance `A / (2*A+11)`, which is strictly
larger than `1/13` exactly when `A > 1`.

The terminal theorem below is deliberately phrased only for arithmetic
progressions.  It does not assert that a general tight twelve-speed family is
an arithmetic progression, and therefore does not close the open sporadic
branch.

No `sorry`; no `native_decide`.
-/

import TournamentH7.LonelyRunner

namespace LonelyRunner
namespace AP12Homogeneity

/-- The `k`th term of the gcd-reduced arithmetic progression. -/
def apTerm (g A D k : ℕ) : ℤ :=
  (g : ℤ) * ((A : ℤ) + (k : ℤ) * (D : ℤ))

/-- The twelve-speed arithmetic progression `g * (A + kD)`, `0 <= k < 12`. -/
def ap12 (g A D : ℕ) : Fin 12 → ℤ :=
  fun k => apTerm g A D k.1

/-- Coprimality supplies a Bezout representative of the half-residue.

The equation is oriented for the later phase identity:
`floor(D/2) = A*u + D*z`. -/
theorem exists_halfResidue_bezout (A D : ℕ) (hcop : Nat.Coprime A D) :
    ∃ u z : ℤ,
      ((D / 2 : ℕ) : ℤ) = (A : ℤ) * u + (D : ℤ) * z := by
  obtain ⟨x, hx⟩ := Int.mod_coprime hcop
  have hscaled := hx.mul_right (((D / 2 : ℕ) : ℤ))
  rw [Int.modEq_iff_add_fac] at hscaled
  obtain ⟨z, hz⟩ := hscaled
  refine ⟨x * ((D / 2 : ℕ) : ℤ), z, ?_⟩
  calc
    ((D / 2 : ℕ) : ℤ) =
        ((A : ℤ) * x) * ((D / 2 : ℕ) : ℤ) + (D : ℤ) * z := by
          simpa using hz
    _ = (A : ℤ) * (x * ((D / 2 : ℕ) : ℤ)) + (D : ℤ) * z := by ring

/-- Exact rational phase collapse from explicit Bezout data.

At time `u/(gD)`, every AP term differs from the common half-residue by the
integer `k*u-z`.  This is an equality in `ℚ`, before passing to the circle. -/
theorem phaseCollapseQ
    (g A D k : ℕ) (u z : ℤ)
    (hg : 0 < g) (hD : 0 < D)
    (hbez : ((D / 2 : ℕ) : ℤ) = (A : ℤ) * u + (D : ℤ) * z) :
    (apTerm g A D k : ℚ) *
        ((u : ℚ) / ((g : ℚ) * (D : ℚ))) =
      ((D / 2 : ℕ) : ℚ) / (D : ℚ) +
        (((k : ℤ) * u - z : ℤ) : ℚ) := by
  have hgQ : (g : ℚ) ≠ 0 := by exact_mod_cast (Nat.ne_of_gt hg)
  have hDQ : (D : ℚ) ≠ 0 := by exact_mod_cast (Nat.ne_of_gt hD)
  have hbezQ :
      ((D / 2 : ℕ) : ℚ) =
        (A : ℚ) * (u : ℚ) + (D : ℚ) * (z : ℚ) := by
    calc
      ((D / 2 : ℕ) : ℚ) = ((((D / 2 : ℕ) : ℤ)) : ℚ) := by
        rw [Int.cast_natCast]
      _ = (((A : ℤ) * u + (D : ℤ) * z : ℤ) : ℚ) := by rw [hbez]
      _ = (A : ℚ) * (u : ℚ) + (D : ℚ) * (z : ℚ) := by push_cast; rfl
  unfold apTerm
  push_cast
  field_simp [hgQ, hDQ]
  linarith [hbezQ]

/-- Real-valued form of `phaseCollapseQ`, suitable for `Lonely`. -/
theorem phaseCollapseR
    (g A D k : ℕ) (u z : ℤ)
    (hg : 0 < g) (hD : 0 < D)
    (hbez : ((D / 2 : ℕ) : ℤ) = (A : ℤ) * u + (D : ℤ) * z) :
    (apTerm g A D k : ℝ) *
        ((u : ℝ) / ((g : ℝ) * (D : ℝ))) =
      ((D / 2 : ℕ) : ℝ) / (D : ℝ) +
        (((k : ℤ) * u - z : ℤ) : ℝ) := by
  have hgR : (g : ℝ) ≠ 0 := by exact_mod_cast (Nat.ne_of_gt hg)
  have hDR : (D : ℝ) ≠ 0 := by exact_mod_cast (Nat.ne_of_gt hD)
  have hbezR :
      ((D / 2 : ℕ) : ℝ) =
        (A : ℝ) * (u : ℝ) + (D : ℝ) * (z : ℝ) := by
    calc
      ((D / 2 : ℕ) : ℝ) = ((((D / 2 : ℕ) : ℤ)) : ℝ) := by
        rw [Int.cast_natCast]
      _ = (((A : ℤ) * u + (D : ℤ) * z : ℤ) : ℝ) := by rw [hbez]
      _ = (A : ℝ) * (u : ℝ) + (D : ℝ) * (z : ℝ) := by push_cast; rfl
  unfold apTerm
  push_cast
  field_simp [hgR, hDR]
  linarith [hbezR]

/-- The half-residue has the sharp uniform lower bound `1/3` for `D >= 2`. -/
theorem one_third_le_halfResidue (D : ℕ) (hD : 2 ≤ D) :
    (1 : ℝ) / 3 ≤ ((D / 2 : ℕ) : ℝ) / (D : ℝ) := by
  have hfloor : D ≤ 3 * (D / 2) := by omega
  have hfloorR : (D : ℝ) ≤ 3 * ((D / 2 : ℕ) : ℝ) := by
    exact_mod_cast hfloor
  have hDR : (0 : ℝ) < (D : ℝ) := by exact_mod_cast (lt_of_lt_of_le (by omega) hD)
  rw [le_div_iff₀ hDR]
  nlinarith

/-- The chosen half-residue lies in the first half of the unit circle. -/
theorem halfResidue_le_one_half (D : ℕ) (hD : 0 < D) :
    ((D / 2 : ℕ) : ℝ) / (D : ℝ) ≤ (1 : ℝ) / 2 := by
  have hfloor : 2 * (D / 2) ≤ D := by omega
  have hfloorR : 2 * ((D / 2 : ℕ) : ℝ) ≤ (D : ℝ) := by
    exact_mod_cast hfloor
  have hDR : (0 : ℝ) < (D : ℝ) := by exact_mod_cast hD
  rw [div_le_iff₀ hDR]
  nlinarith

/-- Exact centered phase for the consecutive (`D=1`) branch. -/
theorem centeredPhaseR (g A k : ℕ) (hg : 0 < g) :
    (apTerm g A 1 k : ℝ) *
        ((1 : ℝ) / ((g : ℝ) * (2 * (A : ℝ) + 11))) =
      ((A : ℝ) + (k : ℝ)) / (2 * (A : ℝ) + 11) := by
  have hgR : (g : ℝ) ≠ 0 := by exact_mod_cast (Nat.ne_of_gt hg)
  unfold apTerm
  push_cast
  field_simp [hgR]

/-- The centered clearance comparison is exact: it beats `1/13` iff `A>1`. -/
theorem centeredClearance_gt_one_thirteenth_iff (A : ℕ) (hA : 0 < A) :
    (1 : ℝ) / 13 < (A : ℝ) / (2 * (A : ℝ) + 11) ↔ 1 < A := by
  have hden : (0 : ℝ) < 2 * (A : ℝ) + 11 := by positivity
  rw [lt_div_iff₀ hden]
  constructor
  · intro h
    have hAR : (1 : ℝ) < (A : ℝ) := by nlinarith
    exact_mod_cast hAR
  · intro h
    have hAR : (1 : ℝ) < (A : ℝ) := by exact_mod_cast h
    nlinarith

/-- Every runner in the `D >= 2` branch clears every integer by at least
`1/3` at the modular-inverse time.  This is stronger than merely proving
`Lonely 13`. -/
theorem spreadAP_clearance
    (g A D : ℕ) (hg : 0 < g) (hD : 2 ≤ D) (hcop : Nat.Coprime A D) :
    ∃ t : ℝ, ∀ k : Fin 12, ∀ m : ℤ,
      (1 : ℝ) / 3 ≤ |(ap12 g A D k : ℝ) * t - (m : ℝ)| := by
  obtain ⟨u, z, hbez⟩ := exists_halfResidue_bezout A D hcop
  let x : ℝ := ((D / 2 : ℕ) : ℝ) / (D : ℝ)
  refine ⟨(u : ℝ) / ((g : ℝ) * (D : ℝ)), ?_⟩
  intro k m
  have hDpos : 0 < D := lt_of_lt_of_le (by omega) hD
  have hxlo : (1 : ℝ) / 3 ≤ x := by
    simpa [x] using one_third_le_halfResidue D hD
  have hxhi : x ≤ 1 - (1 : ℝ) / 3 := by
    have hhalf := halfResidue_le_one_half D hDpos
    dsimp [x]
    nlinarith
  have hx0 : (0 : ℝ) ≤ x := le_trans (by norm_num) hxlo
  have hx1 : x < 1 :=
    lt_of_le_of_lt hxhi (by norm_num : 1 - (1 : ℝ) / 3 < 1)
  have hfract : Int.fract x = x := Int.fract_eq_self.mpr ⟨hx0, hx1⟩
  have hfar : ∀ q : ℤ, (1 : ℝ) / 3 ≤ |x - (q : ℝ)| := by
    apply (far_iff_fract ((1 : ℝ) / 3) x).2
    rw [hfract]
    exact ⟨hxlo, hxhi⟩
  let j : ℤ := (k.1 : ℤ) * u - z
  have hphase :
      (ap12 g A D k : ℝ) *
          ((u : ℝ) / ((g : ℝ) * (D : ℝ))) = x + (j : ℝ) := by
    simpa [ap12, x, j] using
      phaseCollapseR g A D k.1 u z hg hDpos hbez
  rw [hphase]
  have heq : x + (j : ℝ) - (m : ℝ) = x - ((m - j : ℤ) : ℝ) := by
    push_cast
    ring
  rw [heq]
  exact hfar (m - j)

/-- Existing-predicate consumer for the spread branch. -/
theorem spreadAP_lonely
    (g A D : ℕ) (hg : 0 < g) (hD : 2 ≤ D) (hcop : Nat.Coprime A D) :
    ∃ t : ℝ, Lonely 13 (ap12 g A D) t := by
  obtain ⟨t, ht⟩ := spreadAP_clearance g A D hg hD hcop
  refine ⟨t, fun k m => ?_⟩
  exact le_trans (by norm_num : (1 : ℝ) / 13 ≤ 1 / 3) (ht k m)

/-- The spread branch actually beats the `1/13` threshold strictly. -/
theorem spreadAP_strictWitness
    (g A D : ℕ) (hg : 0 < g) (hD : 2 ≤ D) (hcop : Nat.Coprime A D) :
    ∃ t : ℝ, ∀ k : Fin 12, ∀ m : ℤ,
      (1 : ℝ) / 13 < |(ap12 g A D k : ℝ) * t - (m : ℝ)| := by
  obtain ⟨t, ht⟩ := spreadAP_clearance g A D hg hD hcop
  refine ⟨t, fun k m => ?_⟩
  exact lt_of_lt_of_le (by norm_num : (1 : ℝ) / 13 < 1 / 3) (ht k m)

/-- At the centered time, every consecutive-branch runner clears every
integer by at least the exact endpoint value `A/(2A+11)`. -/
theorem centeredAP_clearance
    (g A : ℕ) (hg : 0 < g) (hA : 0 < A) :
    ∃ t : ℝ, ∀ k : Fin 12, ∀ m : ℤ,
      (A : ℝ) / (2 * (A : ℝ) + 11) ≤
        |(ap12 g A 1 k : ℝ) * t - (m : ℝ)| := by
  let den : ℝ := 2 * (A : ℝ) + 11
  let c : ℝ := (A : ℝ) / den
  refine ⟨(1 : ℝ) / ((g : ℝ) * den), ?_⟩
  intro k m
  have hden : 0 < den := by dsimp [den]; positivity
  have hk : k.1 ≤ 11 := by omega
  have hkR : (k.1 : ℝ) ≤ 11 := by exact_mod_cast hk
  let x : ℝ := ((A : ℝ) + (k.1 : ℝ)) / den
  have hxlo : c ≤ x := by
    dsimp [c, x]
    rw [div_le_div_iff_of_pos_right hden]
    have hk0 : (0 : ℝ) ≤ (k.1 : ℝ) := by positivity
    linarith
  have hxhi : x ≤ 1 - c := by
    have hsym : 1 - c = ((A : ℝ) + 11) / den := by
      dsimp [c, den]
      field_simp
      ring
    rw [hsym]
    dsimp [x]
    rw [div_le_div_iff_of_pos_right hden]
    linarith
  have hc0 : (0 : ℝ) ≤ c := by dsimp [c]; positivity
  have hcpos : (0 : ℝ) < c := by dsimp [c]; positivity
  have hx0 : (0 : ℝ) ≤ x := le_trans hc0 hxlo
  have hx1 : x < 1 := lt_of_le_of_lt hxhi (sub_lt_self 1 hcpos)
  have hfract : Int.fract x = x := Int.fract_eq_self.mpr ⟨hx0, hx1⟩
  have hfar : ∀ q : ℤ, c ≤ |x - (q : ℝ)| := by
    apply (far_iff_fract c x).2
    rw [hfract]
    exact ⟨hxlo, hxhi⟩
  have hphase :
      (ap12 g A 1 k : ℝ) * ((1 : ℝ) / ((g : ℝ) * den)) = x := by
    simpa [ap12, den, x] using centeredPhaseR g A k.1 hg
  rw [hphase]
  simpa [c, den] using hfar m

/-- Existing-predicate consumer for the reduced-consecutive branch with
`A > 1`. -/
theorem centeredAP_lonely
    (g A : ℕ) (hg : 0 < g) (hA : 1 < A) :
    ∃ t : ℝ, Lonely 13 (ap12 g A 1) t := by
  have hApos : 0 < A := by omega
  obtain ⟨t, ht⟩ := centeredAP_clearance g A hg hApos
  have hstrict := (centeredClearance_gt_one_thirteenth_iff A hApos).2 hA
  refine ⟨t, fun k m => ?_⟩
  exact le_trans hstrict.le (ht k m)

/-- The nonhomogeneous consecutive branch beats `1/13` strictly. -/
theorem centeredAP_strictWitness
    (g A : ℕ) (hg : 0 < g) (hA : 1 < A) :
    ∃ t : ℝ, ∀ k : Fin 12, ∀ m : ℤ,
      (1 : ℝ) / 13 < |(ap12 g A 1 k : ℝ) * t - (m : ℝ)| := by
  have hApos : 0 < A := by omega
  obtain ⟨t, ht⟩ := centeredAP_clearance g A hg hApos
  have hstrict := (centeredClearance_gt_one_thirteenth_iff A hApos).2 hA
  exact ⟨t, fun k m => lt_of_lt_of_le hstrict (ht k m)⟩

/-- **Terminal AP-only classification.**  If a positive primitive reduced AP
admits no time at which all twelve runners beat `1/13` strictly, then its
reduced parameters are `A=D=1`.

This is the formal homogeneity conclusion of THM-1171.  The hypothesis is a
tightness/no-strict-improvement condition on an AP already in hand; there is
no claim that an arbitrary twelve-speed tight set has AP shape. -/
theorem reducedAP_eq_one_of_no_strictWitness
    (g A D : ℕ) (hg : 0 < g) (hA : 0 < A) (hD : 0 < D)
    (hcop : Nat.Coprime A D)
    (htight : ¬ ∃ t : ℝ, ∀ k : Fin 12, ∀ m : ℤ,
      (1 : ℝ) / 13 < |(ap12 g A D k : ℝ) * t - (m : ℝ)|) :
    A = 1 ∧ D = 1 := by
  rcases (show D = 1 ∨ 2 ≤ D by omega) with rfl | hspread
  · refine ⟨?_, rfl⟩
    by_contra hne
    have hAgt : 1 < A := by omega
    exact htight (centeredAP_strictWitness g A hg hAgt)
  · exact (htight (spreadAP_strictWitness g A D hg hspread hcop)).elim

/-- Pointwise homogeneous form of the terminal classification. -/
theorem ap12_eq_homogeneous_of_no_strictWitness
    (g A D : ℕ) (hg : 0 < g) (hA : 0 < A) (hD : 0 < D)
    (hcop : Nat.Coprime A D)
    (htight : ¬ ∃ t : ℝ, ∀ k : Fin 12, ∀ m : ℤ,
      (1 : ℝ) / 13 < |(ap12 g A D k : ℝ) * t - (m : ℝ)|) :
    ap12 g A D = fun k : Fin 12 => (g : ℤ) * ((k.1 : ℤ) + 1) := by
  obtain ⟨rfl, rfl⟩ :=
    reducedAP_eq_one_of_no_strictWitness g A D hg hA hD hcop htight
  funext k
  change
    (g : ℤ) * ((1 : ℤ) + (k.1 : ℤ) * (1 : ℤ)) =
      (g : ℤ) * ((k.1 : ℤ) + 1)
  ring

#print axioms exists_halfResidue_bezout
#print axioms phaseCollapseQ
#print axioms spreadAP_clearance
#print axioms centeredAP_clearance
#print axioms reducedAP_eq_one_of_no_strictWitness
#print axioms ap12_eq_homogeneous_of_no_strictWitness

end AP12Homogeneity
end LonelyRunner
