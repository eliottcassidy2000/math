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
    ∧ all |speeds| distinct ∧ max |speed| ≥ 23 ∧ no detuned harmonic ∧ no coarse
    decomposition ∧ NO NONTRIVIAL COMMON RESIDUE (branch 8, THM-682(a):
    `LRCCommonResidue.lonely_of_common_residue` — monad-S12).
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
import TournamentH7.LRCDetunedDispatch
import TournamentH7.LRCCoarseReduction
import TournamentH7.LRCCommonResidue
import TournamentH7.LRCWindow22Census

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
    (∀ g : ℤ, 2 ≤ g → ∀ i₀ : Fin 13, (∀ j, j ≠ i₀ → g ∣ v j) → g ∣ v i₀) →
    (¬ ∃ (L : ℤ) (k a : Fin 13 → ℤ) (A : ℝ), (∀ i, v i = a i + L * k i) ∧ 0 < (L : ℝ) ∧
      (∀ i, |(a i : ℝ)| ≤ A) ∧ A / (L : ℝ) ≤ 1/13 - 1/14 ∧ (∀ i, k i ≠ 0) ∧
      (Finset.univ.image k).card ≤ 12) →
    (∀ d : ℤ, 2 ≤ d → ∀ a : ℤ, (∀ i, d ∣ (v i - a)) → d ∣ a) →
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
  -- (4) bounded window: all |v i| ≤ 22 — kernel-pure 6-witness pigeonhole (LEM-024, opus-S202),
  --     replacing the `winData22` native_decide (kps-S127: the 6 witnesses; opus-S202: the pigeonhole)
  by_cases hwin : ∀ i, |v i| ≤ 22
  · exact WindowData.hwindowW_closed 22 cite Window22Census.hdistinct22_kernel v hv hwin
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
  -- (6) detuned harmonic: THM-668 (LRCDetunedDispatch)
  by_cases hdet : ∃ g : ℤ, 2 ≤ g ∧ ∃ i₀ : Fin 13, (∀ j, j ≠ i₀ → g ∣ v j) ∧ ¬ g ∣ v i₀
  · obtain ⟨g, hg2, i₀, hH, hnd⟩ := hdet
    exact DetunedDispatch.lonely14_of_detuned cite v hv g hg2 i₀ hH hnd
  -- (7) multi-scale: the coarse reduction to LRC(≤13)
  by_cases hms : ∃ (L : ℤ) (k a : Fin 13 → ℤ) (A : ℝ), (∀ i, v i = a i + L * k i) ∧
      0 < (L : ℝ) ∧ (∀ i, |(a i : ℝ)| ≤ A) ∧ A / (L : ℝ) ≤ 1/13 - 1/14 ∧
      (∀ i, k i ≠ 0) ∧ (Finset.univ.image k).card ≤ 12
  · obtain ⟨L, k, a, A, hdecomp, hL, ha, hbudget, hkne, hcard⟩ := hms
    exact CoarseReduction.lonely14_of_coarse_le12 cite v k a L A hdecomp hL ha hbudget hkne hcard
  -- (8) nontrivial common residue: THM-682(a) (LRCCommonResidue)
  by_cases hcr : ∃ d : ℤ, 2 ≤ d ∧ ∃ a : ℤ, ¬ d ∣ a ∧ ∀ i, d ∣ (v i - a)
  · obtain ⟨d, hd2, a, hna, hres⟩ := hcr
    exact CommonResidue.lonely_of_common_residue v d a hd2 hna hres
  -- the residual class
  · push_neg at hdom hrep hwin hdet
    obtain ⟨iw, hiw⟩ := hwin
    refine hresidual v hv hcov hgap ?_ ?_ ⟨iw, by omega⟩ ?_ hms ?_
    · intro i
      obtain ⟨j, hji, hj⟩ := hdom i
      exact ⟨j, hji, hj⟩
    · intro i j hij
      exact hrep i j hij
    · intro g hg2 i₀ hH
      exact hdet g hg2 i₀ hH
    · intro d hd2 a hall
      by_contra hnda
      exact hcr ⟨d, hd2, a, hnda, hall⟩

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
    (∀ g : ℤ, 2 ≤ g → ∀ i₀ : Fin 13, (∀ j, j ≠ i₀ → g ∣ v j) → g ∣ v i₀) →
    (¬ ∃ (L : ℤ) (k a : Fin 13 → ℤ) (A : ℝ), (∀ i, v i = a i + L * k i) ∧ 0 < (L : ℝ) ∧
      (∀ i, |(a i : ℝ)| ≤ A) ∧ A / (L : ℝ) ≤ 1/13 - 1/14 ∧ (∀ i, k i ≠ 0) ∧
      (Finset.univ.image k).card ≤ 12) →
    (∀ d : ℤ, 2 ≤ d → ∀ a : ℤ, (∀ i, d ∣ (v i - a)) → d ∣ a) →
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
  -- (6) detuned harmonic: THM-668 (LRCDetunedDispatch)
  by_cases hdet : ∃ g : ℤ, 2 ≤ g ∧ ∃ i₀ : Fin 13, (∀ j, j ≠ i₀ → g ∣ v j) ∧ ¬ g ∣ v i₀
  · obtain ⟨g, hg2, i₀, hH, hnd⟩ := hdet
    exact DetunedDispatch.lonely14_of_detuned cite v hv g hg2 i₀ hH hnd
  -- (7) multi-scale: the coarse reduction to LRC(≤13)
  by_cases hms : ∃ (L : ℤ) (k a : Fin 13 → ℤ) (A : ℝ), (∀ i, v i = a i + L * k i) ∧
      0 < (L : ℝ) ∧ (∀ i, |(a i : ℝ)| ≤ A) ∧ A / (L : ℝ) ≤ 1/13 - 1/14 ∧
      (∀ i, k i ≠ 0) ∧ (Finset.univ.image k).card ≤ 12
  · obtain ⟨L, k, a, A, hdecomp, hL, ha, hbudget, hkne, hcard⟩ := hms
    exact CoarseReduction.lonely14_of_coarse_le12 cite v k a L A hdecomp hL ha hbudget hkne hcard
  -- (8) nontrivial common residue: THM-682(a) (LRCCommonResidue)
  by_cases hcr : ∃ d : ℤ, 2 ≤ d ∧ ∃ a : ℤ, ¬ d ∣ a ∧ ∀ i, d ∣ (v i - a)
  · obtain ⟨d, hd2, a, hna, hres⟩ := hcr
    exact CommonResidue.lonely_of_common_residue v d a hd2 hna hres
  · push_neg at hdom hrep hdet
    refine hresidual v hv hcov hgap ?_ ?_ ?_ hms ?_
    · intro i
      obtain ⟨j, hji, hj⟩ := hdom i
      exact ⟨j, hji, hj⟩
    · intro i j hij
      exact hrep i j hij
    · intro g hg2 i₀ hH
      exact hdet g hg2 i₀ hH
    · intro d hd2 a hall
      by_contra hnda
      exact hcr ⟨d, hd2, a, hnda, hall⟩

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
