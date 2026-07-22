import TournamentH7.GMC2Formalization

/-!
# WIP: the `DvdK1 → NC2` capstone composition (death-star)

Assembles codex's spine into the conditional NC2 theorem: under the one-variable DvdK premise,
a moment-null polynomial is charge-one-sided.  Structure: the descent gives a finite-field torus
point `w` where every integral zero relation is preserved and the lowest-face seed is nonzero;
the normalized moment relation is `0` at `w` (null, preserved) yet `≠ 0` at `w`
(`three_case_sum_ne_zero`, nonzero seed) — contradiction.

Remaining `sorry`s mark the transport obligations still being discharged.
-/

open GMC2 GMC2NormalizedMoment GMC2ResidueAssembly Finset

namespace GMC2NC2Capstone

/-- **Conditional NC2 under one-variable DvdK.** WIP skeleton. -/
theorem nc2_of_dvdk1 (hDvdK : GMC2DvdKInterface.DvdK1) :
    ∀ P : MvPolynomial (Fin 2) ℂ, GMC2.NC2At P := by
  intro P hnull
  by_contra hnotone
  obtain ⟨lambda, delta, F, m0, hsubset, hm0, hlower, hFdef,
      A, D, w, hpprime, hchar, hfin, hfield, hunit, hwnz,
      hpreserve, hmomzero, hseednz⟩ :=
    GMC2IntegralFaceSeedDescent.exists_finite_field_moment_point_preserving_integral_lowest_face_seed
      hDvdK P hnull hnotone
  -- The residue field, its char-p field structure:
  letI : Field D.ResidueField := D.fieldStructure
  -- exponent = the support inclusion:
  set exponent : ↥P.support → Fin 2 →₀ ℕ := fun s => (s : Fin 2 →₀ ℕ) with hexp
  -- A reference channel + face height A0 (from the nonzero seed) — TODO extract.
  -- The min height floor `p*A0` at mass `p*m0` — TODO from balanced_natural_height_floor_of_reference.
  -- Then `aeval w (normalizedMomentRelationInt exponent (p*m0) (p*A0))` is:
  --   (a) = 0   (null ⟹ integral zero relation, preserved by hpreserve)
  --   (b) ≠ 0   (three_case_sum_ne_zero: non-dilated→0, off-face→0, face=w^p sum = seed ≠ 0)
  sorry

end GMC2NC2Capstone
