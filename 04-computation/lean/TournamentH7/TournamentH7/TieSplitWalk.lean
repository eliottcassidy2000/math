/- TieSplitWalk.lean — mac-mini-2026-07-16-S127.
   THM-866's kernel arithmetic, formalized (rung two of the LRC(14)/metagraph ladder).
   (1) The F3 tie-split identity: flipping the arc between two TIED vertices raises
       x = Σ d_v² by exactly 8 (pure ring arithmetic).
   (2) The general per-flip law (THM-855 F3): Δx = 4*(d_l - d_w) + 8.
   (3) Distinct scores on n vertices with total C(n,2) are exactly {0, …, n-1}
       (the pigeonhole half of "no tie ⟹ transitive"). -/
import Mathlib.Data.Finset.Sort
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Data.Int.Interval
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith

namespace LRC14

/-- **(2) The per-flip drop law (THM-855 F3).** Flipping arc w→l sends
    `d_w ↦ d_w - 2, d_l ↦ d_l + 2`; the change in `d_w² + d_l²` is `4*(d_l - d_w) + 8`. -/
theorem flip_delta (dw dl : ℤ) :
    (dw - 2)^2 + (dl + 2)^2 - (dw^2 + dl^2) = 4*(dl - dw) + 8 := by ring

/-- **(1) The tie-split +8 (THM-866's step).** If `d_w = d_l`, the flip raises x by 8. -/
theorem tie_split_plus_eight (d : ℤ) :
    (d - 2)^2 + (d + 2)^2 - (d^2 + d^2) = 8 := by ring

/-- **(3) Pigeonhole.** An injective score map `s : Fin n → ℤ` with `0 ≤ s v ≤ n-1`
    hits every value in `{0, …, n-1}` (so distinct scores force the transitive
    score vector). -/
theorem scores_are_range (n : ℕ) (s : Fin n → ℤ)
    (hinj : Function.Injective s)
    (hlo : ∀ v, 0 ≤ s v) (hhi : ∀ v, s v ≤ (n : ℤ) - 1) :
    ∀ k : ℤ, 0 ≤ k → k ≤ (n : ℤ) - 1 → ∃ v, s v = k := by
  classical
  -- the image is a subset of Icc 0 (n-1) of full cardinality, hence equal
  have himage : (Finset.univ.image s) ⊆ Finset.Icc 0 ((n : ℤ) - 1) := by
    intro k hk
    rcases Finset.mem_image.mp hk with ⟨v, _, rfl⟩
    exact Finset.mem_Icc.mpr ⟨hlo v, hhi v⟩
  have hcard : (Finset.univ.image s).card = n := by
    rw [Finset.card_image_of_injective _ hinj, Finset.card_univ, Fintype.card_fin]
  have hIcc : (Finset.Icc (0 : ℤ) ((n : ℤ) - 1)).card = n := by
    rw [Int.card_Icc]
    omega
  have heq : Finset.univ.image s = Finset.Icc 0 ((n : ℤ) - 1) :=
    Finset.eq_of_subset_of_card_le himage (by rw [hcard, hIcc])
  intro k hk0 hk1
  have : k ∈ Finset.univ.image s := heq ▸ Finset.mem_Icc.mpr ⟨hk0, hk1⟩
  rcases Finset.mem_image.mp this with ⟨v, _, hv⟩
  exact ⟨v, hv⟩

/-- **The walk's step count.** From level `x` to the ceiling `c` in exact +8 steps:
    if `c - x = 8 * k` then the walk length is `k` (the graded-chain arithmetic). -/
theorem walk_length (x c : ℤ) (k : ℕ) (h : c - x = 8 * k) :
    x + 8 * k = c := by omega

end LRC14
