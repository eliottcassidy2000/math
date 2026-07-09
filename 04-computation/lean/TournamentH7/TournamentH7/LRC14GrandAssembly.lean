/-
  TournamentH7.LRC14GrandAssembly — THE TOP-LEVEL SURFACE (monad-explorer-2026-07-09-S6, HYP-5757).

  ONE theorem deriving `LRC14Statement` from the LRC(≤13) citation plus ONE residual
  hypothesis, with every other branch discharged by the existing sorry-free corpus:

    (1) NON-COVERING (some q ∈ {2,…,14} divides no speed):
        `sieve_one_div` — explicit lonely time t = 1/q.            [unconditional]
    (2) ALL-COMPARABLE (every pair of |speeds| within factor 13):
        `spread13_lonely` — explicit t = 1/(min+max).              [unconditional]
    (3) DOMINANT (one runner exceeds 13× all others):
        `LRC14.hdom_discharged` — the sharp dominant peel.         [cite : LRCUpTo13]
    (4) BOUNDED WINDOW (all |speeds| ≤ 22):
        `WindowData.hwindow22_closed` — machine-checked window.    [cite : LRCUpTo13]
    (5) REPEATED |SPEED| (two runners share an absolute speed):
        `lonely_of_abs` (here) + `lonely14_of_repeat`.             [cite : LRCUpTo13]

  THE RESIDUAL (the single remaining obligation, `ResidualObligation` below):
    covering ∧ scale-gapped (some ratio > 13) ∧ compressed (no dominant runner)
    ∧ all |speeds| distinct ∧ max |speed| ≥ 23.
  This is strictly sharper than every prior surface in the corpus
  (`lrc14_of_compressed`: covering ∧ compressed; `lrc14_endgame`: opaque-witnessG2).

  Also here: `lonely_comp_perm` (permutation invariance) and `covering18_complete`
  (the kps-S115 966-list is COMPLETE over [1,18] — one native_decide over the 8568
  13-subsets), pinning the [1,18] base case end-to-end.

  Kernel-pure except where noted (`covering18_complete` uses native_decide).
-/
import Mathlib
import TournamentH7.LonelyRunner
import TournamentH7.LRCFourteenSkeleton
import TournamentH7.LRCSpread13
import TournamentH7.LRC14CertRoute
import TournamentH7.LRCEndgameAssembly
import TournamentH7.LRCWindowData22
import TournamentH7.LRC13Citation
import TournamentH7.LRCCovering966

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner

/-- A 13-speed family has a SCALE GAP if some pair of speeds differs by more than
a factor 13 in absolute value. -/
def GapFamily (v : Fin 13 → ℤ) : Prop :=
  ¬ ∀ i j : Fin 13, |v i| ≤ 13 * |v j|

/-- **THE RESIDUAL SURFACE.**  The single remaining obligation after branches (1)-(5):
covering, scale-gapped, compressed (no dominant runner), all absolute speeds distinct,
and reaching beyond the machine-checked window. -/
def ResidualObligation : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
    LRC14.CoveringFamily v →
    GapFamily v →
    (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) →
    (∀ i j, i ≠ j → |v i| ≠ |v j|) →
    (∃ i, 23 ≤ |v i|) →
    ∃ t : ℝ, Lonely 14 v t

/-- Loneliness transfers from the absolute-value family back to the signed family
(the `∀ m : ℤ` quantifier absorbs the sign). -/
theorem lonely_of_abs {n : ℕ} (v : Fin 13 → ℤ) (t : ℝ)
    (h : Lonely n (fun i => |v i|) t) : Lonely n v t := by
  intro i m
  rcases abs_choice (v i) with hcase | hcase
  · have := h i m
    simpa [hcase] using this
  · have hthis := h i (-m)
    simp only [hcase] at hthis
    have heq : |((-(v i) : ℤ) : ℝ) * t - ((-m : ℤ) : ℝ)| = |(v i : ℝ) * t - m| := by
      push_cast
      rw [show -(v i : ℝ) * t - -(m : ℝ) = -((v i : ℝ) * t - m) by ring, abs_neg]
    calc (1 : ℝ) / n ≤ |((-(v i) : ℤ) : ℝ) * t - ((-m : ℤ) : ℝ)| := by exact_mod_cast hthis
    _ = |(v i : ℝ) * t - m| := heq

/-- **THE GRAND ASSEMBLY.**  `LRC14Statement` from the LRC(≤13) citation and the
single residual obligation. -/
theorem lrc14_grand_assembly (cite : LRCUpTo13) (hresidual : ResidualObligation) :
    LRC14.LRC14Statement := by
  intro v hv
  -- (1) non-covering: the denominator sieve at the missing modulus
  by_cases hcov : LRC14.CoveringFamily v
  swap
  · rw [LRC14.CoveringFamily] at hcov
    push_neg at hcov
    obtain ⟨q, h2, h14, hdiv⟩ := hcov
    exact ⟨(1 : ℝ) / q, sieve_one_div 14 q v h14 (by omega) hdiv⟩
  -- (4) bounded window: all |v i| ≤ 22
  by_cases hwin : ∀ i, |v i| ≤ 22
  · exact WindowData.hwindow22_closed cite v hv hwin
  -- (2) all-comparable: ratio ≤ 13 throughout
  by_cases hgap : GapFamily v
  swap
  · rw [GapFamily, not_not] at hgap
    obtain ⟨imax, -, hmax⟩ :=
      Finset.exists_max_image (Finset.univ) (fun i => |v i|) ⟨0, Finset.mem_univ 0⟩
    obtain ⟨jmin, -, hmin⟩ :=
      Finset.exists_min_image (Finset.univ) (fun i => |v i|) ⟨0, Finset.mem_univ 0⟩
    exact ⟨_, LRC14.spread13_lonely v (|v jmin|) (|v imax|)
      (abs_pos.mpr (hv jmin))
      (fun i => hmin i (Finset.mem_univ i))
      (fun i => hmax i (Finset.mem_univ i))
      (hgap imax jmin)⟩
  -- (3) dominant runner: the sharp peel
  by_cases hdom : ∃ i, ∀ j, j ≠ i → 13 * |v j| < |v i|
  · exact LRC14.hdom_discharged cite v hv hcov hdom
  -- (5) repeated absolute speed: reduce to ≤ 12 distinct values, cite
  by_cases hrep : ∃ i j, i ≠ j ∧ |v i| = |v j|
  · obtain ⟨i, j, hij, heq⟩ := hrep
    obtain ⟨t, ht⟩ := lonely14_of_repeat cite (fun i => |v i|)
      (fun i => by simpa using hv i) hij heq
    exact ⟨t, lonely_of_abs v t ht⟩
  -- the residual class
  · push_neg at hdom hrep hwin
    obtain ⟨iw, hiw⟩ := hwin
    refine hresidual v hv hcov hgap ?_ ?_ ⟨iw, by omega⟩
    · intro i
      obtain ⟨j, hji, hj⟩ := hdom i
      exact ⟨j, hji, hj⟩
    · intro i j hij
      exact hrep i j hij

/-- **The kernel-pure variant** — identical surface WITHOUT the window-22 branch (which
carries two `native_decide` certificate axioms).  The residual class here additionally
contains the `max |speed| ≤ 22` families; axioms: propext, Classical.choice, Quot.sound
only. -/
def ResidualObligationPure : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
    LRC14.CoveringFamily v →
    GapFamily v →
    (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) →
    (∀ i j, i ≠ j → |v i| ≠ |v j|) →
    ∃ t : ℝ, Lonely 14 v t

theorem lrc14_grand_assembly_pure (cite : LRCUpTo13)
    (hresidual : ResidualObligationPure) : LRC14.LRC14Statement := by
  intro v hv
  by_cases hcov : LRC14.CoveringFamily v
  swap
  · rw [LRC14.CoveringFamily] at hcov
    push_neg at hcov
    obtain ⟨q, h2, h14, hdiv⟩ := hcov
    exact ⟨(1 : ℝ) / q, sieve_one_div 14 q v h14 (by omega) hdiv⟩
  by_cases hgap : GapFamily v
  swap
  · rw [GapFamily, not_not] at hgap
    obtain ⟨imax, -, hmax⟩ :=
      Finset.exists_max_image (Finset.univ) (fun i => |v i|) ⟨0, Finset.mem_univ 0⟩
    obtain ⟨jmin, -, hmin⟩ :=
      Finset.exists_min_image (Finset.univ) (fun i => |v i|) ⟨0, Finset.mem_univ 0⟩
    exact ⟨_, LRC14.spread13_lonely v (|v jmin|) (|v imax|)
      (abs_pos.mpr (hv jmin))
      (fun i => hmin i (Finset.mem_univ i))
      (fun i => hmax i (Finset.mem_univ i))
      (hgap imax jmin)⟩
  by_cases hdom : ∃ i, ∀ j, j ≠ i → 13 * |v j| < |v i|
  · exact LRC14.hdom_discharged cite v hv hcov hdom
  by_cases hrep : ∃ i j, i ≠ j ∧ |v i| = |v j|
  · obtain ⟨i, j, hij, heq⟩ := hrep
    obtain ⟨t, ht⟩ := lonely14_of_repeat cite (fun i => |v i|)
      (fun i => by simpa using hv i) hij heq
    exact ⟨t, lonely_of_abs v t ht⟩
  · push_neg at hdom hrep
    refine hresidual v hv hcov hgap ?_ ?_
    · intro i
      obtain ⟨j, hji, hj⟩ := hdom i
      exact ⟨j, hji, hj⟩
    · intro i j hij
      exact hrep i j hij

/-- `Lonely` is invariant under permuting the family. -/
theorem lonely_comp_perm {n : ℕ} (v : Fin 13 → ℤ) (t : ℝ) (σ : Equiv.Perm (Fin 13))
    (h : Lonely n v t) : Lonely n (v ∘ σ) t :=
  fun i m => h (σ i) m

/-- **The 966-list is COMPLETE over [1,18]** — every primitive covering 13-subset of
{1,…,18} appears among kps-S115's certified witnesses (one decidable sweep over the
C(18,13) = 8568 subsets).  With `coveringWitnesses_lonely`, the [1,18] slice is
machine-closed end-to-end. -/
theorem covering18_complete :
    ∀ l ∈ (List.range 18).sublistsLen 13,
      LRC14Concrete.coveringPrim (l.map (fun x => (x : ℤ) + 1)) →
        (l.map (fun x => (x : ℤ) + 1)) ∈
          LRC14Concrete.coveringWitnesses.map Prod.fst := by
  native_decide

/-! ## Axiom audit -/
#print axioms lrc14_grand_assembly
#print axioms lrc14_grand_assembly_pure
#print axioms lonely_of_abs
#print axioms lonely_comp_perm
#print axioms covering18_complete

end LRC14Grand
end LonelyRunner
