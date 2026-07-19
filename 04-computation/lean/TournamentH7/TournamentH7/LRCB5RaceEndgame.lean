/-
  TournamentH7.LRCB5RaceEndgame

  An endgame socket for the B5 race scoreboard.  The certificate below does
  not assert odd-tail equidistribution: it records the two required tail
  estimates explicitly, together with the strict numerical margin that makes
  the scoreboard positive.  Thus this module isolates the remaining analytic
  obligation without pretending to prove it.

  No `sorry`; no `native_decide`.
-/

import TournamentH7.LRCB5Race
import TournamentH7.LRCSelectedWitnessResidual

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner
open scoped Classical

/-- A proof-producing certificate for the odd-tail side of the B5 race.

The `tail3_upper` and `tail5_upper` fields are the actual finite deviation-sum
bounds owed by equidistribution.  The last field says that, after paying those
two bounds and the already-proved even-layer floors, the explicit scoreboard
still has positive margin. -/
structure B5RaceTailCertificate (v : Fin 13 → ℤ) where
  q : ℕ
  fourteen_le_q : 14 ≤ q
  primitive_mod_q : ∀ i, Int.gcd (v i) (q : ℤ) = 1
  tail3 : ℚ
  tail5 : ℚ
  tail3_upper :
    (∑ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard 3,
      LRC14Concrete.deviation v q T) ≤ tail3
  tail5_upper :
    (∑ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard 5,
      LRC14Concrete.deviation v q T) ≤ tail5
  scoreboard_rhs_pos :
    0 < ((q : ℚ) - 1) * (2052 / 16807)
      - 78 * (((q : ℚ) - 1) / 49)
      - 715 * (((q : ℚ) - 1) / 2401)
      - tail3 - tail5

/-- A race-tail certificate gives the positive integer B5 witness consumed by
the discrete Bonferroni endgame. -/
theorem B5RaceTailCertificate.b5_pos {v : Fin 13 → ℤ}
    (certificate : B5RaceTailCertificate v) :
    0 < LRC14Concrete.B5 v certificate.q := by
  have hscore := LRC14Concrete.B5_race_scoreboard v certificate.q
    certificate.fourteen_le_q certificate.primitive_mod_q
    certificate.tail3 certificate.tail5 certificate.tail3_upper
    certificate.tail5_upper
  have hposQ : (0 : ℚ) < (LRC14Concrete.B5 v certificate.q : ℚ) :=
    lt_of_lt_of_le certificate.scoreboard_rhs_pos hscore
  exact_mod_cast hposQ

/-- Odd-tail race certificates only on the primitive, dissociated,
chain-dense core.  Its tuple hypotheses are exactly those of
`DenseCoreDissociatedB5Supply`; the only change is that the output carries the
two explicit odd-tail estimates rather than bare B5 positivity. -/
def DenseCoreB5RaceTailSupply : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → LRC14.tupleGcd v = 1 →
    LRC14.CoveringFamily v → GapFamily v →
    (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) →
    (∀ i j, i ≠ j → |v i| ≠ |v j|) →
    (∃ i, 23 ≤ |v i|) →
    (∀ g : ℤ, 2 ≤ g → ∀ i₀ : Fin 13,
      (∀ j, j ≠ i₀ → g ∣ v j) → g ∣ v i₀) →
    (¬ ∃ (L : ℤ) (k a : Fin 13 → ℤ) (A : ℝ),
      (∀ i, v i = a i + L * k i) ∧ 0 < (L : ℝ) ∧
      (∀ i, |(a i : ℝ)| ≤ A) ∧ A / (L : ℝ) ≤ 1 / 13 - 1 / 14 ∧
      (∀ i, k i ≠ 0) ∧ (Finset.univ.image k).card ≤ 12) →
    (∀ g : ℤ, 2 ≤ g →
      nonMultCard v g ≠ 2 ∧ nonMultCard v g ≠ 3) →
    (∃ σ : Equiv.Perm (Fin 13), Monotone (fun i ↦ |v (σ i)|) ∧
      ChainDenseCore (fun i ↦ |v (σ i)|)) →
    Nonempty (B5RaceTailCertificate v)

/-- The explicit race-tail supplier discharges the bare positive-B5 supplier
used by the chain-dense endgame. -/
theorem denseCoreDissociatedB5Supply_of_B5RaceTail
    (hsupply : DenseCoreB5RaceTailSupply) :
    DenseCoreDissociatedB5Supply := by
  intro v hv hgcd hcov hgap hcomp hdist hlarge hdiv hcoarse hdissoc hcore
  obtain ⟨certificate⟩ :=
    hsupply v hv hgcd hcov hgap hcomp hdist hlarge hdiv hcoarse hdissoc hcore
  exact ⟨certificate.q, by omega, certificate.b5_pos⟩

/-- Residual-selected capstone with the B5 race socket.  The selected-phase
and pair-tower hypotheses are already confined to their nongeneric residue;
the dense-core hypothesis is precisely the still-open odd-tail certificate
supply. -/
theorem lrc14_from_residualSelectedWitnessSupplies_and_B5RaceTail
    (cite : LRCUpTo13)
    (hpairs : NonterminatingPairTowerSupply)
    (h22 : ResidualTwoTwoSelectedWitnessSupply)
    (h244 : ResidualTwoFourFourSelectedWitnessSupply)
    (h333 : ResidualUniformThreeSelectedWitnessSupply)
    (hsupply : DenseCoreB5RaceTailSupply) :
    LRC14.LRC14Statement :=
  lrc14_from_twoThree_detuned_and_denseCore_dissociated_B5 cite
    (deepExceptionalDetunedDispatchTwoThree_of_finalResidues cite
      (deepExceptionalDetunedDispatchFinalResidues_of_residualSelectedWitnessSupplies
        hpairs h22 h244 h333))
    (denseCoreDissociatedB5Supply_of_B5RaceTail hsupply)

/-! ## Axiom audit -/

#print axioms B5RaceTailCertificate.b5_pos
#print axioms denseCoreDissociatedB5Supply_of_B5RaceTail
#print axioms lrc14_from_residualSelectedWitnessSupplies_and_B5RaceTail

end LRC14Grand
end LonelyRunner
