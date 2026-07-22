import TournamentH7.GMC2FaceSeed
import TournamentH7.GMC2DvdKUniqueChannel

/-!
# Discharging the DvdK1 input for the unique-channel class

The GMC(2) spine consumes the one-variable DvdK theorem in exactly one place: `GMC2FaceSeed.
exists_nonzero_lowest_face_seed` calls `DvdK1` to produce a nonzero constant-term seed on the
rational lowest face.  Everything upstream — the slope `λ`, level `δ`, the exact face `F`, charge
injectivity, straddling — comes from *pure Newton-polygon geometry* (`exists_rational_lowest_face_
finset`), with no analytic input.

So on any polynomial whose lowest face carries a **unique balanced channel** (death-star-S101's 84%
of supports), the seed can be supplied by the elementary `ct_ne_zero_of_unique_balanced`
(`GMC2DvdKUniqueChannel`) instead of `DvdK1`.  `exists_nonzero_lowest_face_seed_of_uniqueChannel`
below has the *same conclusion* as codex's `exists_nonzero_lowest_face_seed` but takes a
unique-channel hypothesis in place of the `DvdK1` premise — a drop-in, DvdK-axiom-free replacement
for the whole unique-channel class.  The residual is exactly the coincident-channel stratum (`card ≥
2`, the symmetric/resonant 16%), which is the irreducible DvdK content handled on paper by codex
THM-2067.
-/

open MvPolynomial Finset

namespace GMC2DvdKUniqueChannelBypass

/-- A support (charge vector `q`) has a **unique balanced channel** if some positive size `m0`
carries exactly one balanced composition.  This is the property `ct_ne_zero_of_unique_balanced`
needs, and by death-star-S101 it holds for 84% of straddling supports. -/
def HasUniqueBalancedChannel {ι : Type*} [Fintype ι] [DecidableEq ι] (q : ι → ℤ) : Prop :=
  ∃ m0 : ℕ, 1 ≤ m0 ∧ ∃ r0 : ι → ℕ,
    r0 ∈ Finset.piAntidiag (Finset.univ : Finset ι) m0 ∧
    GMC2ConstantTermRelations.totalCharge q r0 = 0 ∧
    ∀ r ∈ Finset.piAntidiag (Finset.univ : Finset ι) m0,
      GMC2ConstantTermRelations.totalCharge q r = 0 → r = r0

/-- The lowest-face hypothesis replacing `DvdK1` for one polynomial: every face cut out by a slope
`λ`/level `δ` on `P`'s support has a unique balanced channel in its charges. -/
def LowestFaceUniqueChannel (P : MvPolynomial (Fin 2) ℂ) : Prop :=
  ∀ (lambda delta : ℚ) (F : Finset (Fin 2 →₀ ℕ)),
    (∀ s, s ∈ F ↔ s ∈ P.support ∧
      GMC2.radialExponentQ s - lambda * GMC2.chargeQ s = delta) →
    HasUniqueBalancedChannel (fun s : ↥F => GMC2.charge s)

/-- **DvdK-free lowest-face seed for the unique-channel class.**  Identical conclusion to
`GMC2FaceSeed.exists_nonzero_lowest_face_seed`, but with the `DvdK1` premise replaced by the
unique-channel hypothesis `LowestFaceUniqueChannel P`.  The seed is produced by the elementary
`ct_ne_zero_of_unique_balanced`; no external analytic axiom is used. -/
theorem exists_nonzero_lowest_face_seed_of_uniqueChannel
    (P : MvPolynomial (Fin 2) ℂ) (hP : ¬GMC2.ChargeOneSided P)
    (hUC : LowestFaceUniqueChannel P) :
    ∃ lambda delta : ℚ, ∃ F : Finset (Fin 2 →₀ ℕ), ∃ m0 : ℕ,
      1 ≤ m0 ∧
      F ⊆ P.support ∧
      (∀ s ∈ P.support,
        delta ≤ GMC2.radialExponentQ s - lambda * GMC2.chargeQ s) ∧
      (∀ s, s ∈ F ↔ s ∈ P.support ∧
        GMC2.radialExponentQ s - lambda * GMC2.chargeQ s = delta) ∧
      MvPolynomial.aeval (fun s : ↥F => P.coeff s)
        (GMC2ConstantTermRelations.constantTermRelation
          (fun s : ↥F => GMC2.charge s) m0) ≠ 0 := by
  obtain ⟨lambda, delta, F, hsubset, hlower, hface, _hstraddle⟩ :=
    GMC2.exists_rational_lowest_face_finset P hP
  obtain ⟨m0, hm0, r0, hr0mem, hr0bal, huniq⟩ := hUC lambda delta F hface
  have hcoeff : ∀ s : ↥F, P.coeff s ≠ 0 := by
    intro s
    have hsSupport : (s : Fin 2 →₀ ℕ) ∈ P.support := hsubset s.property
    simpa only [MvPolynomial.mem_support_iff] using hsSupport
  have hseed := GMC2DvdKUniqueChannel.ct_ne_zero_of_unique_balanced
    (fun s : ↥F => GMC2.charge s) (fun s : ↥F => P.coeff s) hcoeff m0 r0 hr0mem hr0bal huniq
  exact ⟨lambda, delta, F, m0, hm0, hsubset, hlower, hface, hseed⟩

end GMC2DvdKUniqueChannelBypass

#print axioms GMC2DvdKUniqueChannelBypass.exists_nonzero_lowest_face_seed_of_uniqueChannel
