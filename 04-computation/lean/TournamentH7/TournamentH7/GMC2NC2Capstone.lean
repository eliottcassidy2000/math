import TournamentH7.GMC2NC2

/-!
# Compatibility surface for the NC2/GMC(2) capstone

The former work-in-progress theorem in this module contained one `sorry` and
overstated the current formal boundary as `DvdK1 -> NC2`.  The checked
capstone lives in `GMC2NC2`: after a compact normalized-height witness is
supplied, the direct finite-field specialization is simultaneously zero by
moment nullity and nonzero by the Frobenius face residue.

This module retains the historical name as a small, sorry-free compatibility
surface.  It deliberately exposes both remaining inputs in its theorem names.
-/

namespace GMC2NC2Capstone

/-- The exact remaining internal composition interface. -/
abbrev HeightWitnessSupplier := GMC2NC2.HeightWitnessSupplier

/-- Checked NC2 endpoint from the published DvdK input and the explicit
height-witness supplier. -/
theorem nc2_of_dvdk1_of_heightWitnessSupplier
    (hDvdK : GMC2DvdKInterface.DvdK1)
    (hHeight : HeightWitnessSupplier) : GMC2.NC2 :=
  GMC2NC2.nc2_of_dvdK1_of_heightWitnessSupplier hDvdK hHeight

/-- Checked GMC(2) endpoint through the same two visible inputs. -/
theorem gmc2_of_dvdk1_of_heightWitnessSupplier
    (hDvdK : GMC2DvdKInterface.DvdK1)
    (hHeight : HeightWitnessSupplier)
    (P Q : MvPolynomial (Fin 2) ℂ)
    (hnull : ∀ m : ℕ, 1 ≤ m → GMC2.E (P ^ m) = 0) :
    ∃ N : ℕ, ∀ m ≥ N, GMC2.E (Q * P ^ m) = 0 :=
  GMC2NC2.gmc2_of_dvdK1_of_heightWitnessSupplier
    hDvdK hHeight P Q hnull

end GMC2NC2Capstone

#print axioms GMC2NC2Capstone.nc2_of_dvdk1_of_heightWitnessSupplier
#print axioms GMC2NC2Capstone.gmc2_of_dvdk1_of_heightWitnessSupplier
