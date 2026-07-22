import TournamentH7.GMC2LowestFacePackage
import TournamentH7.GMC2FrobeniusFace

/-!
# Dictionary between polynomial bidegrees and Frobenius-face coordinates

The existence layer uses `GMC2.charge` on natural Finsupp bidegrees, whereas
the face arithmetic layer accepts arbitrary integer exponent functions.  The
two representations are definitionally the same after the natural casts; the
lemmas here make that cross-module transport explicit.
-/

namespace GMC2FaceDictionary

/-- Integer `Z`- and `W`-exponents of an exact bidegree. -/
def exponentA (s : Fin 2 →₀ ℕ) : ℤ := s 0

def exponentB (s : Fin 2 →₀ ℕ) : ℤ := s 1

@[simp] theorem charge_eq (s : Fin 2 →₀ ℕ) :
    GMC2FrobeniusFace.charge exponentA exponentB s = GMC2.charge s := by
  simp [GMC2FrobeniusFace.charge, exponentA, exponentB, GMC2.charge]

@[simp] theorem tiltedHeight_eq (lambda : ℚ) (s : Fin 2 →₀ ℕ) :
    GMC2FrobeniusFace.tiltedHeight exponentA exponentB lambda s =
      GMC2.radialExponentQ s - lambda * GMC2.chargeQ s := by
  simp [GMC2FrobeniusFace.tiltedHeight, GMC2.radialExponentQ, GMC2.chargeQ]
  norm_cast

end GMC2FaceDictionary

#print axioms GMC2FaceDictionary.charge_eq
#print axioms GMC2FaceDictionary.tiltedHeight_eq
