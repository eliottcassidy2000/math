/-
  TournamentH7.LRCCrystallographicThetaFrontier -- S263 / HYP-3110.

  This module records the De Moivre/Jacobi/crystallographic sidecar lane for
  the LRC(14) proof frontier.  It does not prove LRC(14).  It adds:

    * exact finite catalog counts for the 17 wallpaper groups, 230 3D space
      groups, 14 3D Bravais lattices, and four Jacobi theta channels;
    * a checked De Moivre quintic normal-form cancellation over `Rat`;
    * open predicates saying how theta/orbifold sidecars must feed the existing
      HYP-3107 proof-frontier interface.
-/

import TournamentH7.LRCProofFrontier

namespace LonelyRunner
namespace LRC14

/-! ## Finite crystallographic and theta catalogs -/

/-- The classical 17 wallpaper groups, as a finite sidecar catalog. -/
inductive WallpaperGroup where
  | p1
  | p2
  | pm
  | pg
  | cm
  | pmm
  | pmg
  | pgg
  | cmm
  | p4
  | p4m
  | p4g
  | p3
  | p3m1
  | p31m
  | p6
  | p6m
deriving DecidableEq, Repr, Fintype

/-- There are 17 wallpaper groups. -/
theorem wallpaperGroup_count : Fintype.card WallpaperGroup = 17 := by
  native_decide

/-- The three-dimensional space-group catalog is represented by its 230
numbered entries.  HYP-3110 uses the count as a quotient-audit budget, not as a
geometric proof by itself. -/
abbrev SpaceGroup3D := Fin 230

/-- There are 230 three-dimensional crystallographic space groups. -/
theorem spaceGroup3D_count : Fintype.card SpaceGroup3D = 230 := by
  simp [SpaceGroup3D]

/-- The 14 three-dimensional Bravais lattice types, represented by index. -/
abbrev Bravais3D := Fin 14

/-- There are 14 three-dimensional Bravais lattice types. -/
theorem bravais3D_count : Fintype.card Bravais3D = 14 := by
  simp [Bravais3D]

/-- The Jacobi theta channels relevant to sidecar bookkeeping. -/
inductive JacobiThetaChannel where
  | theta1_odd
  | theta2_shifted
  | theta3_lattice
  | theta4_alternating
deriving DecidableEq, Repr, Fintype

/-- Four Jacobi theta channels are tracked. -/
theorem jacobiThetaChannel_count : Fintype.card JacobiThetaChannel = 4 := by
  native_decide

/-- HYP-3110 sidecars, treated as proof carriers rather than runners. -/
inductive CrystallographicThetaCarrier where
  | finite_address_exit
  | observer_gluing_certificate
  | jacobi_theta_tail
  | lee_yang_root_curve
  | de_moivre_quintic_fold
  | space_group_230_orbifold
  | wallpaper_17_orbifold
  | raw_scalar_shadow
deriving DecidableEq, Repr, Fintype

/-- Eight carrier vertices are in the HYP-3110 scout tournament. -/
theorem crystallographicThetaCarrier_count :
    Fintype.card CrystallographicThetaCarrier = 8 := by
  native_decide

/-! ## De Moivre quintic fold -/

/-- De Moivre's solvable quintic normal-form cancellation.  If
`x = u - a/u`, then the quintic part collapses to a two-term expression:

`x^5 + 5a*x^3 + 5a^2*x = u^5 - a^5/u^5`.

HYP-3110 uses this as a finite-depth cancellation detector.  It is not the
general Abel-Ruffini/A5 obstruction lane. -/
theorem deMoivreQuintic_fold_rat (a u : Rat) (hu : u ≠ 0) :
    (u - a / u) ^ 5 + 5 * a * (u - a / u) ^ 3 +
        5 * a ^ 2 * (u - a / u) = u ^ 5 - a ^ 5 / u ^ 5 := by
  field_simp [hu]
  ring

/-! ## Open sidecar obligations -/

/-- A Jacobi-theta tail certificate for the signed residue-cusp/support-six
tail after low-height wall deletion. -/
def JacobiThetaTailCertificate
    (isTailRow : (Fin 13 -> Int) -> Prop)
    (thetaBudget : (Fin 13 -> Int) -> Real) : Prop :=
  forall v : Fin 13 -> Int, (forall i, v i ≠ 0) ->
    isTailRow v -> thetaBudget v <= (1 : Real) / 14

/-- A crystallographic quotient audit: every quotient must state its wallpaper
or space-group sidecar, preserve the LRC predicate, and identify the destroyed
coordinate. -/
structure CrystallographicQuotientAudit where
  usesWallpaper : Bool
  usesSpaceGroup : Bool
  hasTranslationLattice : Bool
  hasStabilizerWord : Bool
  hasGlideOrScrewSidecar : Bool
  preservedLRCPredicate : Bool
  destroyedCoordinateNamed : Bool

/-- The quotient is legally sidecar-complete only when all required flags are
present. -/
def CrystallographicQuotientAudit.complete
    (q : CrystallographicQuotientAudit) : Prop :=
  q.hasTranslationLattice = true ∧
  q.hasStabilizerWord = true ∧
  q.hasGlideOrScrewSidecar = true ∧
  q.preservedLRCPredicate = true ∧
  q.destroyedCoordinateNamed = true

/-- A residual row producer that may use theta/crystallographic sidecars, but
must end in an existing proof-frontier packet interface. -/
def ThetaCrystallographicResidualProducer
    (isResidual : (Fin 13 -> Int) -> Prop)
    (audit : (Fin 13 -> Int) -> CrystallographicQuotientAudit) : Prop :=
  forall v : Fin 13 -> Int, (forall i, v i ≠ 0) ->
    isResidual v ->
      (audit v).complete ∧
        (Nonempty (ObserverGluingCertificate v) ∨
         Nonempty (FiniteAddressBranchPacket v))

/-- HYP-3110 closes against the existing frontier only after the residual
producer emits observer-gluing or finite-address packets. -/
theorem lrc14_from_theta_crystallographic_residuals
    (isResidual : (Fin 13 -> Int) -> Prop)
    (audit : (Fin 13 -> Int) -> CrystallographicQuotientAudit)
    (classifier : FrontierResidualClassifier isResidual)
    (producer : ThetaCrystallographicResidualProducer isResidual audit) :
    LRC14Statement := by
  refine lrc14_from_residual_observer_gluing_frontier classifier ?_
  intro v hv hres
  rcases producer v hv hres with ⟨_, hpacket⟩
  rcases hpacket with hobs | hfa
  · exact hobs
  · rcases hfa with ⟨packet⟩
    exact ⟨ObserverGluingCertificate.ofFiniteAddressBranchPacket packet⟩

#print axioms wallpaperGroup_count
#print axioms spaceGroup3D_count
#print axioms bravais3D_count
#print axioms jacobiThetaChannel_count
#print axioms crystallographicThetaCarrier_count
#print axioms deMoivreQuintic_fold_rat
#print axioms lrc14_from_theta_crystallographic_residuals

end LRC14
end LonelyRunner
