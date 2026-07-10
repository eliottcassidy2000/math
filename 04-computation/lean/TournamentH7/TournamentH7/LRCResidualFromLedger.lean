/-
  TournamentH7.LRCResidualFromLedger — wiring the pair-sum ledger into the grand assembly's residual
  (kind-pasteur-2026-07-09-S118).

  monad-explorer-S6's `lrc14_grand_assembly` (THM-671) derives the full `LRC14Statement` from the
  LRC(≤13) citation and a SINGLE `ResidualObligation`: every covering, scale-gapped, compressed,
  distinct-speed family with some `|v_i| ≥ 23` has a lonely instant.  That residual is exactly the
  covering ratio>13 case the THM-668 pair-sum machinery targets.

  This file reduces that analytic obligation (`∃ t, Lonely 14 v t`) to a **discrete, number-theoretic**
  one — the pair-sum ledger's own conclusion, `HasLiveRuler`: some pair-sum modulus `q` has fewer than
  `q−1` blocked multipliers.  The bridge is my consumer chain: `mreach_ge_of_blocked_lt` (kps-S117,
  a live ruler ⟹ `Mreach ≥ 1/14`) composed with `lonely_of_Mreach_ge` (a lonely instant exists).  So

      `lrc14_from_ledger : LRC(≤13) → [every residual family has a live pair-sum ruler] → LRC14`.

  The remaining obligation is now the ledger's exact form — mac-mini's C1 gcd-ledger / THM-675
  descent-burden / klein's signed box all prove `HasLiveRuler`.  This is the concrete bridge from the
  top-level assembly down to the pair-sum liveness census.  (`fires v q p = fires |v| q p` since the
  band `[q/14,13q/14]` is symmetric under `r ↦ q−r`, so the ledger's abs-speed statement supplies this
  verbatim.)  Builds on `LRC14GrandAssembly` and `LRCLedgerConsumer`.
-/
import Mathlib
import TournamentH7.LRC14GrandAssembly
import TournamentH7.LRCLedgerConsumer

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner LRC14Concrete

/-- **The pair-sum ledger's conclusion for a family.**  Some pair-sum modulus `q > 0` is *live*: among
the nonzero multipliers `{1,…,q−1}`, fewer than `q−1` fail to fire the ruler (leave the band).  This is
exactly what mac-mini's gcd-exact ledger / THM-675 / klein's signed box certify. -/
def HasLiveRuler (v : Fin 13 → ℤ) : Prop :=
  ∃ q : ℤ, 0 < q ∧ ∃ N : ℕ, (N : ℤ) = q - 1 ∧
    ((Finset.range N).filter (fun p : ℕ => ¬ fires v q ((p : ℤ) + 1))).card < N

/-- A live pair-sum ruler yields a lonely instant (consumer chain: `mreach_ge_of_blocked_lt` then
`lonely_of_Mreach_ge`). -/
theorem lonely_of_hasLiveRuler (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) (h : HasLiveRuler v) :
    ∃ t : ℝ, Lonely 14 v t := by
  obtain ⟨q, hq, N, hN, hcount⟩ := h
  have hM := mreach_ge_of_blocked_lt v q hq N hN hcount
  exact LRC14Concrete.lonely_of_Mreach_ge v hv hM

/-- **The ledger discharges the residual obligation.**  If every residual-class family has a live
pair-sum ruler, the grand assembly's `ResidualObligation` holds. -/
theorem residualObligation_of_ledger
    (hledger : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → LRC14.CoveringFamily v → GapFamily v →
      (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) → (∀ i j, i ≠ j → |v i| ≠ |v j|) →
      (∃ i, 23 ≤ |v i|) → HasLiveRuler v) :
    ResidualObligation := by
  intro v hv hcov hgap hcomp hdist hlarge
  exact lonely_of_hasLiveRuler v hv (hledger v hv hcov hgap hcomp hdist hlarge)

/-- **LRC(14) from LRC(≤13) and the pair-sum ledger.**  Composing with `lrc14_grand_assembly`: the full
14-runner Lonely-Runner statement follows from the ≤13 citation and the single **discrete** obligation
that every residual-class family has a live pair-sum ruler.  This converts the assembly's last analytic
surface into the number-theoretic ledger that mac-mini/klein are proving. -/
theorem lrc14_from_ledger (cite : LRCUpTo13)
    (hledger : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → LRC14.CoveringFamily v → GapFamily v →
      (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) → (∀ i j, i ≠ j → |v i| ≠ |v j|) →
      (∃ i, 23 ≤ |v i|) → HasLiveRuler v) :
    LRC14.LRC14Statement :=
  lrc14_grand_assembly cite (residualObligation_of_ledger hledger)

end LRC14Grand
end LonelyRunner
