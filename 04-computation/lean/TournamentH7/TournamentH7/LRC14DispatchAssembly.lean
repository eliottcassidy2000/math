/-
  TournamentH7.LRC14DispatchAssembly -- honest assembly of the LRC(14)
  sieve dispatch and AP-core/dominance bridge.

  `LRCSieveDispatch` proves that only the covering class remains.
  `LRCAPCoreBridge` proves the hard/dominant class once an inverse theorem
  supplies a 13-fold dominant speed.  These facts do not, by themselves,
  prove that every covering family belongs to that hard class: covering does
  not imply dominance.  This module records the genuinely needed easy/hard
  split as an explicit hypothesis and checks the resulting composition.

  Kernel-pure: no `sorry`, no `native_decide`, and no custom axioms.
-/
import TournamentH7.LRCAPCoreBridge
import TournamentH7.LRCSieveDispatch

namespace LonelyRunner

/-- A predicate-level partition of every positive covering family into an
easy class or the compact class to which an inverse theorem applies.  In the
analytic LRC application, `Easy` is the class with an immediate positive
margin and `Compact` is the genuinely small-margin residual. -/
def CoveringSplit (Easy Compact : (Fin 13 → ℤ) → Prop) : Prop :=
  ∀ (v : Fin 13 → ℤ), (∀ i, 0 < v i) → Covering v → Easy v ∨ Compact v

/-- The easy side of a covering split already has a `1/14`-lonely time. -/
def EasyCase (Easy : (Fin 13 → ℤ) → Prop) : Prop :=
  ∀ (v : Fin 13 → ℤ), (∀ i, 0 < v i) → Easy v → ∃ t : ℝ, Lonely 14 v t

/-- Honest covering-case assembly.  Besides the LRC(≤13) citation and the
inverse theorem on `Compact`, one must supply both the covering easy/hard
split and the easy-side witness.  In particular, `INV Compact` alone is not a
proof of `CoveringCase`. -/
theorem coveringCase_of_split {Easy Compact : (Fin 13 → ℤ) → Prop}
    (cite : LRCUpTo13) (inv : INV Compact)
    (split : CoveringSplit Easy Compact) (easy : EasyCase Easy) :
    CoveringCase := by
  intro v hpos hcover
  rcases split v hpos hcover with heasy | hcompact
  · exact easy v hpos heasy
  · exact lonely14_of_INV cite inv v hpos hcompact

/-- Full dispatch assembly: non-covering families are discharged by the
sieve, easy covering families by `EasyCase`, and compact covering families by
`INV` followed by the AP-core descent bridge. -/
theorem lrc14_of_covering_split {Easy Compact : (Fin 13 → ℤ) → Prop}
    (cite : LRCUpTo13) (inv : INV Compact)
    (split : CoveringSplit Easy Compact) (easy : EasyCase Easy)
    (v : Fin 13 → ℤ) (hpos : ∀ i, 0 < v i) :
    ∃ t : ℝ, Lonely 14 v t :=
  lrc14_of_covering (coveringCase_of_split cite inv split easy) v hpos

end LonelyRunner

#print axioms LonelyRunner.coveringCase_of_split
#print axioms LonelyRunner.lrc14_of_covering_split
