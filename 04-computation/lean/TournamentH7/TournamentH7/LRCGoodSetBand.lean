/-
  TournamentH7.LRCGoodSetBand — the goodSet band lemma (mac-mini-2026-07-09-S65 cont.21).

  CLAIM STUB (honest): statements below are the claim; proofs land this session.
  NOT root-imported until sorry-free.

  The last plumbing lemma of the witness-floor discharge pipeline: a checkable
  per-difference band condition forces `Ico α β ⊆ goodSet E` — the `hgood` input of
  cont.18's `witnessG2_pos_of_anchor`, and (composed with brick (iii)) the full
  m_P floors of `hsmall3`/`hlarge` become engine-lists-in, proof-out.

  Band shapes (e := the positive difference; endpoints rational in practice):
  * up (a < b, e = b − a):  ∃ j, j + 1/7 < e·α  ∧  e·β ≤ j + 1
      — on [α,β): e·x ∈ (j + 1/7, j + 1), so ⌊e·x⌋ = j and fract > 1/7.
  * down (b < a, e = a − b): ∃ j, j ≤ e·α  ∧  e·β ≤ j + 6/7
      — e·x ∈ [j, j + 6/7), fract(e·x) < 6/7, and the reflection
        fract(−y) ∈ {0} ∪ {1 − fract y} lands outside (0, 1/7].
  * b = a: the zero difference, free (`fract 0 = 0`).
-/
import Mathlib
import TournamentH7.LRCFourteenSkeleton
import TournamentH7.LRCGoodSetSmall
import TournamentH7.LRCUnionVolume

namespace LonelyRunner
namespace LRC14
namespace GoodSetBand

open DenseCovers MeasureTheory TournamentH7.GoodSet

/-- The per-element band certificate for tooth `a`, element `b`, interval `[α, β)`. -/
def bandCert (a b : ℤ) (α β : ℝ) : Prop :=
  b = a ∨
    (a < b ∧ ∃ j : ℤ, (j : ℝ) + 1 / 7 < ((b - a : ℤ) : ℝ) * α ∧
      ((b - a : ℤ) : ℝ) * β ≤ (j : ℝ) + 1) ∨
    (b < a ∧ ∃ j : ℤ, (j : ℝ) ≤ ((a - b : ℤ) : ℝ) * α ∧
      ((a - b : ℤ) : ℝ) * β ≤ (j : ℝ) + 6 / 7)

/-- **The goodSet band lemma.** A tooth `a ∈ E` whose every difference passes the band
check puts the whole interval inside `goodSet E`. -/
theorem Ico_subset_goodSet_of_bounds {E : List ℤ} {a : ℤ} (ha : a ∈ E) {α β : ℝ}
    (hband : ∀ b ∈ E, bandCert a b α β) :
    Set.Ico α β ⊆ goodSet E := by
  sorry

/-- Both band checks at once: the interval sits inside the intersection `witnessG2` measures. -/
theorem Ico_subset_good_inter_safe {E P : List ℤ} {a : ℤ} (ha : a ∈ E) {α β : ℝ}
    (hgood : ∀ b ∈ E, bandCert a b α β)
    (hposP : ∀ p ∈ P, (0 : ℤ) < p)
    (hsafe : ∀ p ∈ P, ∃ j : ℤ, (j : ℝ) + 1 / 14 ≤ (p : ℝ) * α ∧
      (p : ℝ) * β ≤ (j : ℝ) + 13 / 14) :
    Set.Ico α β ⊆ goodSet E ∩ safeSet P := by
  sorry

/-- **The leg-shaped floor**: a sorted-disjoint engine list, each interval carrying a tooth
and both band certificates, forces `witnessG2 s ≥ Σ lengths`. -/
theorem witnessG2_ge_of_sorted_bands (s : Shape) (l : List (ℝ × ℝ))
    (hposP : ∀ p ∈ s.1, (0 : ℤ) < p)
    (hcert : ∀ q ∈ l, ∃ a ∈ s.2, (∀ b ∈ s.2, bandCert a b q.1 q.2) ∧
      (∀ p ∈ s.1, ∃ j : ℤ, (j : ℝ) + 1 / 14 ≤ (p : ℝ) * q.1 ∧
        (p : ℝ) * q.2 ≤ (j : ℝ) + 13 / 14))
    (hin : ∀ q ∈ l, 0 ≤ q.1 ∧ q.1 ≤ q.2 ∧ q.2 ≤ 1)
    (hsorted : l.Pairwise (fun q r => q.2 ≤ r.1)) :
    (l.map (fun q => q.2 - q.1)).sum ≤ witnessG2 s := by
  sorry

end GoodSetBand
end LRC14
end LonelyRunner
