/-
  TournamentH7.LRCDilationSymmetry -- LAYER 2 of the LRC(14) consec-max wall is
  FORCED by multiplicative (dilation) symmetry.  (HYP-2756/2757, mac-mini-2026-06-21-S19.)

  Setup.  After STRATUM-LOCALIZATION (HYP-2749) the LRC(14) consec-extremality reduces
  to the full-Z/7-residue stratum.  A full-residue k=8 offset shape has residue multiset
  {0,1,2,3,4,5,6} together with exactly one REPEATED residue `r*`.  HYP-2753 localizes the
  remaining wall to three layers; LAYER 2 is the empirical fact "consec doubles the IDENTITY
  residue 0".  This file proves the structural reason: that multiset is invariant under the
  multiplicative (dilation) group Z/7* IFF `r* = 0` -- the unique dilation-FIXED residue.
  Hence LAYER 2 is forced by symmetry, mirroring Paley=QR as the unique Z/7*-symmetric
  tournament (HYP-2755).  Everything here is finite and decidable (mathlib-free core Lean).
-/
namespace TournamentH7.LRCDilationSymmetry

/-- Multiplicity of residue `r` (mod 7) in the multiset `{0,1,...,6}` with `rstar` repeated. -/
def cnt (rstar r : Nat) : Nat := 1 + (if r % 7 == rstar % 7 then 1 else 0)

/-- The residue multiset (with `rstar` doubled) is fixed by dilation by `c` for every unit
    `c ∈ {1,...,6}` of Z/7 iff `cnt (c*r) = cnt r` for all residues `r`.  Quantifier bounds are
    written `< 7` first so the proposition is a decidable bounded quantifier in core Lean. -/
def dilFixed (rstar : Nat) : Prop :=
  ∀ c, c < 7 → 0 < c → ∀ r, r < 7 → cnt rstar ((c * r) % 7) = cnt rstar r

instance (rstar : Nat) : Decidable (dilFixed rstar) := by unfold dilFixed; infer_instance

/-- **LAYER 2, forced by multiplicative symmetry.**  Among the doubled-residue choices,
    only the identity residue `0` yields a Z/7*-(dilation)-invariant residue multiset.
    So "double the identity residue" is not a coincidence of the cover measure: it is the
    unique dilation-symmetric configuration, and consec is its minimal-magnitude shape. -/
theorem layer2_only_identity :
    dilFixed 0 ∧ (∀ rstar, rstar < 7 → 0 < rstar → ¬ dilFixed rstar) := by decide

/-- Restatement: the set of dilation-fixed doubled-residues in {0,...,6} is exactly `{0}`. -/
theorem dilFixed_iff_zero : ∀ rstar, rstar < 7 → (dilFixed rstar ↔ rstar = 0) := by decide

end TournamentH7.LRCDilationSymmetry
