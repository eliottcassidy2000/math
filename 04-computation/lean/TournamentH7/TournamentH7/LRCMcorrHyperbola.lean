/-
# LRCMcorrHyperbola — wiring the hyperbola box count into the per-cell correlation bound

Connects death-star's `LRCHyperbolaBox.hyperbola_box_count` (the character-free combinatorial per-cell
count) to the diagonal-suppression anchor `MultCorrelation.offdiag_mcorr_sq_le`.

The ratio-`w` correlation `zcorr A w = #{(a,b) ∈ A² : a = w·b}` (`= mcorr A w` on units, since
`a·b⁻¹ = w ⟺ a = w·b`) injects — via `(a,b) ↦ b` — into death-star's box `{k ≠ 0 : cdist k ≤ K,
cdist(w·k) ≤ K}` whenever `A` avoids `0` and is `K`-close to `0`.  Hence `hyperbola_box_count` gives the
**per-cell bound** `(zcorr A w − 1)·P ≤ 4K²` for any hyperbola-minimum `P` on the ratio `w` — exactly the
`mcorr A w ≤ M` input (`M = 1 + 4K²/P`) that `offdiag_mcorr_sq_le` aggregates over `w ≠ 1` into the total
`t₂` off-diagonal energy bound.  The three files close the character route's combinatorial core:
`hyperbola_box_count` (per-cell count, death-star) → `zcorr_percell` (this file) → `offdiag_mcorr_sq_le`
(L² aggregation, kps).

kind-pasteur-2026-07-10-S127.
-/
import Mathlib
import TournamentH7.LRCHyperbolaBox
import TournamentH7.LRCMultCorrelation

namespace LonelyRunner
namespace McorrHyperbola

open HyperbolaBox Finset

variable {q : ℕ} [NeZero q]

/-- The ratio-`w` correlation on `ZMod q`: `#{(a,b) ∈ A² : a = w·b}`.  On a set of units this is
`MultCorrelation.mcorr A w` (`a·b⁻¹ = w ⟺ a = w·b`). -/
def zcorr (A : Finset (ZMod q)) (w : ZMod q) : ℕ :=
  ((A ×ˢ A).filter (fun p => p.1 = w * p.2)).card

/-- **The correlation is bounded by the hyperbola box.**  For `A` avoiding `0` with every element within
circle-distance `K` of `0`, the map `(a,b) ↦ b` injects the ratio-`w` pairs into
`{k ≠ 0 : cdist k ≤ K, cdist(w·k) ≤ K}` (indeed `b ∈ A` and `w·b = a ∈ A` are both nonzero and `K`-close). -/
theorem zcorr_le_box (A : Finset (ZMod q)) (w : ZMod q) (K : ℕ)
    (h0 : ∀ a ∈ A, a ≠ 0) (hK : ∀ a ∈ A, cdist a ≤ K) :
    zcorr A w ≤ ((Finset.univ : Finset (ZMod q)).filter
        (fun k => k ≠ 0 ∧ cdist k ≤ K ∧ cdist (w * k) ≤ K)).card := by
  rw [zcorr]
  apply Finset.card_le_card_of_injOn (fun p => p.2)
  · rintro ⟨a, b⟩ hp
    simp only [Finset.mem_coe, Finset.mem_filter, Finset.mem_product] at hp
    obtain ⟨⟨ha, hb⟩, hab⟩ := hp
    simp only [Finset.mem_coe, Finset.mem_filter, Finset.mem_univ, true_and]
    exact ⟨h0 b hb, hK b hb, by rw [← hab]; exact hK a ha⟩
  · rintro ⟨a, b⟩ hp ⟨c, d⟩ hq hbd
    simp only [Finset.mem_coe, Finset.mem_filter, Finset.mem_product] at hp hq
    obtain ⟨_, hab⟩ := hp; obtain ⟨_, hcd⟩ := hq
    have hbd' : b = d := hbd
    exact Prod.ext_iff.mpr ⟨by rw [hab, hcd, hbd'], hbd'⟩

/-- **Wiring: the per-cell correlation bound from `hyperbola_box_count`.**  Given a hyperbola-minimum `P`
on the ratio `w` (`∀ h ≠ 0, P ≤ cdist h · cdist(w·h)` — the generic/non-resonant condition), a `K`-bounded
`A` obeys `(zcorr A w − 1)·P ≤ 4K²`.  This is the per-cell input `mcorr A w ≤ M` (`M = 1 + 4K²/P`) that
`MultCorrelation.offdiag_mcorr_sq_le` consumes over all off-diagonal ratios `w ≠ 1`. -/
theorem zcorr_percell (A : Finset (ZMod q)) (w : ZMod q) (K P : ℕ)
    (h0 : ∀ a ∈ A, a ≠ 0) (hK : ∀ a ∈ A, cdist a ≤ K)
    (hP : ∀ h : ZMod q, h ≠ 0 → P ≤ cdist h * cdist (w * h)) :
    ((zcorr A w : ℤ) - 1) * (P : ℤ) ≤ 4 * K * K := by
  have hbox := hyperbola_box_count (q := q) w K K P hP
  have hle := zcorr_le_box A w K h0 hK
  have hcast : (zcorr A w : ℤ) ≤ (((Finset.univ : Finset (ZMod q)).filter
      (fun k => k ≠ 0 ∧ cdist k ≤ K ∧ cdist (w * k) ≤ K)).card : ℤ) := by exact_mod_cast hle
  calc ((zcorr A w : ℤ) - 1) * (P : ℤ)
      ≤ ((((Finset.univ : Finset (ZMod q)).filter
          (fun k => k ≠ 0 ∧ cdist k ≤ K ∧ cdist (w * k) ≤ K)).card : ℤ) - 1) * (P : ℤ) :=
        mul_le_mul_of_nonneg_right (by omega) (by positivity)
    _ ≤ 4 * K * K := hbox

end McorrHyperbola
end LonelyRunner
