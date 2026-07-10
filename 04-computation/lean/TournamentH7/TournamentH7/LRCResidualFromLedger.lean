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
import TournamentH7.LRCDiscreteBonferroni

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
      (∃ i, 23 ≤ |v i|) →
      (∀ g : ℤ, 2 ≤ g → ∀ i₀ : Fin 13, (∀ j, j ≠ i₀ → g ∣ v j) → g ∣ v i₀) →
      (¬ ∃ (L : ℤ) (k a : Fin 13 → ℤ) (A : ℝ), (∀ i, v i = a i + L * k i) ∧ 0 < (L : ℝ) ∧
        (∀ i, |(a i : ℝ)| ≤ A) ∧ A / (L : ℝ) ≤ 1/13 - 1/14 ∧ (∀ i, k i ≠ 0) ∧
        (Finset.univ.image k).card ≤ 12) →
      HasLiveRuler v) :
    ResidualObligation := by
  intro v hv hcov hgap hcomp hdist hlarge hdiv hcoarse
  exact lonely_of_hasLiveRuler v hv (hledger v hv hcov hgap hcomp hdist hlarge hdiv hcoarse)

/-- **LRC(14) from LRC(≤13) and the pair-sum ledger.**  Composing with `lrc14_grand_assembly`: the full
14-runner Lonely-Runner statement follows from the ≤13 citation and the single **discrete** obligation
that every residual-class family has a live pair-sum ruler.  This converts the assembly's last analytic
surface into the number-theoretic ledger that mac-mini/klein are proving. -/
theorem lrc14_from_ledger (cite : LRCUpTo13)
    (hledger : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → LRC14.CoveringFamily v → GapFamily v →
      (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) → (∀ i j, i ≠ j → |v i| ≠ |v j|) →
      (∃ i, 23 ≤ |v i|) →
      (∀ g : ℤ, 2 ≤ g → ∀ i₀ : Fin 13, (∀ j, j ≠ i₀ → g ∣ v j) → g ∣ v i₀) →
      (¬ ∃ (L : ℤ) (k a : Fin 13 → ℤ) (A : ℝ), (∀ i, v i = a i + L * k i) ∧ 0 < (L : ℝ) ∧
        (∀ i, |(a i : ℝ)| ≤ A) ∧ A / (L : ℝ) ≤ 1/13 - 1/14 ∧ (∀ i, k i ≠ 0) ∧
        (Finset.univ.image k).card ≤ 12) →
      HasLiveRuler v) :
    LRC14.LRC14Statement :=
  lrc14_grand_assembly cite (residualObligation_of_ledger hledger)

/-- **monad's THM-680 conclusion → `HasLiveRuler`.**  THM-680 defines a multiplier `p` on a pair-sum
ruler `q` as *live* when `v_l·p mod q ∈ B` for every `l` — this is exactly `fires v q p` — and its
per-ruler liveness floor `LM > 0` supplies such a live `p`.  A single live `p ∈ {1,…,q−1}` is unblocked,
so the blocked-multiplier count is `< q−1`: `HasLiveRuler v`.  This is the bridge that turns THM-680's
Fourier/Parseval liveness floor into the discrete obligation the assembly consumes. -/
theorem hasLiveRuler_of_fires (v : Fin 13 → ℤ) (q : ℤ) (hq : 0 < q) (m : ℤ)
    (hm1 : 1 ≤ m) (hmq : m ≤ q - 1) (hfire : fires v q m) : HasLiveRuler v := by
  obtain ⟨N, hN⟩ : ∃ N : ℕ, (N : ℤ) = q - 1 := ⟨(q - 1).toNat, Int.toNat_of_nonneg (by omega)⟩
  obtain ⟨p0, hp0⟩ : ∃ p0 : ℕ, (p0 : ℤ) = m - 1 := ⟨(m - 1).toNat, Int.toNat_of_nonneg (by omega)⟩
  refine ⟨q, hq, N, hN, ?_⟩
  have hp0lt : p0 < N := by
    have h : (p0 : ℤ) < (N : ℤ) := by rw [hp0, hN]; omega
    exact_mod_cast h
  have hp0mem : p0 ∈ Finset.range N := Finset.mem_range.mpr hp0lt
  have hp0fires : fires v q ((p0 : ℤ) + 1) := by
    rw [hp0, show m - 1 + 1 = m from by ring]; exact hfire
  have hne : (Finset.range N).filter (fun p : ℕ => ¬ fires v q ((p : ℤ) + 1)) ≠ Finset.range N := by
    intro heq
    have hmem : p0 ∈ (Finset.range N).filter (fun p : ℕ => ¬ fires v q ((p : ℤ) + 1)) := heq.symm ▸ hp0mem
    exact (Finset.mem_filter.mp hmem).2 hp0fires
  have hlt := Finset.card_lt_card
    (Finset.ssubset_iff_subset_ne.mpr ⟨Finset.filter_subset _ _, hne⟩)
  rwa [Finset.card_range] at hlt

/-- **`HasLiveRuler` from THM-680's existence conclusion** (some pair-sum ruler carries a live
multiplier). -/
theorem hasLiveRuler_of_exists_live (v : Fin 13 → ℤ)
    (h : ∃ q : ℤ, 0 < q ∧ ∃ m : ℤ, 1 ≤ m ∧ m ≤ q - 1 ∧ fires v q m) : HasLiveRuler v := by
  obtain ⟨q, hq, m, hm1, hmq, hfire⟩ := h
  exact hasLiveRuler_of_fires v q hq m hm1 hmq hfire

/-- **LRC(14) from LRC(≤13) and monad's THM-680 (per-ruler liveness) on the residual class.**  When
THM-680 supplies a live multiplier for every residual family, this closes LRC(14) through the grand
assembly — the Fourier/Parseval liveness route to the finish line. -/
theorem lrc14_from_liveness (cite : LRCUpTo13)
    (hlive : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → LRC14.CoveringFamily v → GapFamily v →
      (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) → (∀ i j, i ≠ j → |v i| ≠ |v j|) →
      (∃ i, 23 ≤ |v i|) →
      (∀ g : ℤ, 2 ≤ g → ∀ i₀ : Fin 13, (∀ j, j ≠ i₀ → g ∣ v j) → g ∣ v i₀) →
      (¬ ∃ (L : ℤ) (k a : Fin 13 → ℤ) (A : ℝ), (∀ i, v i = a i + L * k i) ∧ 0 < (L : ℝ) ∧
        (∀ i, |(a i : ℝ)| ≤ A) ∧ A / (L : ℝ) ≤ 1/13 - 1/14 ∧ (∀ i, k i ≠ 0) ∧
        (Finset.univ.image k).card ≤ 12) →
      ∃ q : ℤ, 0 < q ∧ ∃ m : ℤ, 1 ≤ m ∧ m ≤ q - 1 ∧ fires v q m) :
    LRC14.LRC14Statement :=
  lrc14_from_ledger cite (fun v hv hcov hgap hcomp hdist hlarge hdiv hcoarse =>
    hasLiveRuler_of_exists_live v (hlive v hv hcov hgap hcomp hdist hlarge hdiv hcoarse))

/-- **The B5 discrete-Bonferroni certifier discharges the residual** (folding in death-star's
`LRCDiscreteBonferroni`).  For a residual family, a single pair-sum ruler `q` with positive depth-5
Bonferroni bound `B5 v q > 0` forces a live multiplier (`B5 ≤ liveCount`), hence `Mreach ≥ 1/14`
(`mreach_ge_of_B5_pos`, which routes through my `mreach_ge_of_pairsum_band`), hence `∃ t, Lonely`.  This
is the fold death-star requested; `B5 > 0` is the exact-integer, `native_decide`-checkable certificate
whose a-priori supply for the residual is THM-671 part 6 (monad's THM-680 Fourier floor / klein's box). -/
theorem residualObligation_of_B5
    (hB5 : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → LRC14.CoveringFamily v → GapFamily v →
      (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) → (∀ i j, i ≠ j → |v i| ≠ |v j|) →
      (∃ i, 23 ≤ |v i|) →
      (∀ g : ℤ, 2 ≤ g → ∀ i₀ : Fin 13, (∀ j, j ≠ i₀ → g ∣ v j) → g ∣ v i₀) →
      (¬ ∃ (L : ℤ) (k a : Fin 13 → ℤ) (A : ℝ), (∀ i, v i = a i + L * k i) ∧ 0 < (L : ℝ) ∧
        (∀ i, |(a i : ℝ)| ≤ A) ∧ A / (L : ℝ) ≤ 1/13 - 1/14 ∧ (∀ i, k i ≠ 0) ∧
        (Finset.univ.image k).card ≤ 12) →
      ∃ q : ℕ, 0 < q ∧ 0 < LRC14Concrete.B5 v q) :
    ResidualObligation := by
  intro v hv hcov hgap hcomp hdist hlarge hdiv hcoarse
  obtain ⟨q, hq, hB5pos⟩ := hB5 v hv hcov hgap hcomp hdist hlarge hdiv hcoarse
  exact LRC14Concrete.lonely_of_Mreach_ge v hv (LRC14Concrete.mreach_ge_of_B5_pos v q hq hB5pos)

/-- **LRC(14) from LRC(≤13) and the B5 certifier's a-priori supply on the residual class.**  The
discrete-Bonferroni finish: LRC(14) follows once every residual family carries a positive B5 bound at
some pair-sum ruler — the sole remaining a-priori obligation (THM-671 part 6). -/
theorem lrc14_from_B5 (cite : LRCUpTo13)
    (hB5 : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → LRC14.CoveringFamily v → GapFamily v →
      (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) → (∀ i j, i ≠ j → |v i| ≠ |v j|) →
      (∃ i, 23 ≤ |v i|) →
      (∀ g : ℤ, 2 ≤ g → ∀ i₀ : Fin 13, (∀ j, j ≠ i₀ → g ∣ v j) → g ∣ v i₀) →
      (¬ ∃ (L : ℤ) (k a : Fin 13 → ℤ) (A : ℝ), (∀ i, v i = a i + L * k i) ∧ 0 < (L : ℝ) ∧
        (∀ i, |(a i : ℝ)| ≤ A) ∧ A / (L : ℝ) ≤ 1/13 - 1/14 ∧ (∀ i, k i ≠ 0) ∧
        (Finset.univ.image k).card ≤ 12) →
      ∃ q : ℕ, 0 < q ∧ 0 < LRC14Concrete.B5 v q) :
    LRC14.LRC14Statement :=
  lrc14_grand_assembly cite (residualObligation_of_B5 hB5)

end LRC14Grand
end LonelyRunner
