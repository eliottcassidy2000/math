/-
# LRCMultCorrelation — the multiplicative pair-correlation anchor (diagonal suppression, L² assembly)

The clean L² anchor of the multiplicative character-frame route (klein S225–226 / death-star
`LRCHyperbolaBox`).  The *multiplicative pair-correlation* of a set `A` at ratio `w` is

  `mcorr A w = #{(a,b) ∈ A² : a·b⁻¹ = w}`.

Two identities pin its shape: the **total** `∑_w mcorr A w = |A|²` (every pair has a ratio) and the
**diagonal** `mcorr A 1 = |A|` (the pairs `(a,a)`).  The *diagonal-suppression law* (klein) says the
off-diagonal mass `∑_{w≠1} mcorr² ` — the `t₂` object — is small for generic (covering, non-geometric)
sets.  This file supplies the **assembly anchor**: a uniform per-cell bound `mcorr A w ≤ M` (`w ≠ 1`) — the
equidistribution input that `LRCHyperbolaBox.hyperbola_box_count` provides — converts into the total
off-diagonal energy bound

  `∑_{w≠1} (mcorr A w)² ≤ M · (|A|² − |A|)`.

This is the safe, structural half of the character bound (an *identity*/L²-conversion, not the refuted
*absolute* Cauchy–Schwarz — death-star-S11): the per-cell smallness is the analytic content, this is its
aggregation.  It is the multiplicative twin of the additive resonance-sum aggregation, and its `w = 1`
diagonal is where the multiplicative geometric rigidity (LEM-023) concentrates.

kind-pasteur-2026-07-09-S127.
-/
import Mathlib

namespace LonelyRunner
namespace MultCorrelation

open Finset

variable {G : Type*} [Group G] [Fintype G] [DecidableEq G]

/-- The **multiplicative pair-correlation** of `A` at ratio `w`: the number of ordered pairs
`(a,b) ∈ A²` with `a · b⁻¹ = w`. -/
def mcorr (A : Finset G) (w : G) : ℕ := ((A ×ˢ A).filter (fun p => p.1 * p.2⁻¹ = w)).card

/-- **The total correlation is `|A|²`** — every ordered pair has exactly one ratio. -/
theorem sum_mcorr (A : Finset G) : ∑ w : G, mcorr A w = A.card ^ 2 := by
  have h : (A ×ˢ A).card = ∑ w : G, mcorr A w :=
    Finset.card_eq_sum_card_fiberwise (t := (Finset.univ : Finset G))
      (f := fun p : G × G => p.1 * p.2⁻¹) (fun _ _ => Finset.mem_univ _)
  rw [← h, Finset.card_product, sq]

/-- **The diagonal correlation is `|A|`** — `w = 1` counts the pairs `(a,a)`. -/
theorem mcorr_one (A : Finset G) : mcorr A 1 = A.card := by
  rw [mcorr]
  apply Finset.card_bij (fun p _ => p.1)
  · rintro ⟨a, b⟩ hp
    simp only [Finset.mem_filter, Finset.mem_product, mul_inv_eq_one] at hp
    exact hp.1.1
  · rintro ⟨a, b⟩ hp ⟨c, d⟩ hq hac
    simp only [Finset.mem_filter, Finset.mem_product, mul_inv_eq_one] at hp hq
    obtain ⟨_, hab⟩ := hp; obtain ⟨_, hcd⟩ := hq
    have hac' : a = c := hac
    exact Prod.ext hac' ((hab.symm.trans hac').trans hcd)
  · intro a ha
    refine ⟨(a, a), ?_, rfl⟩
    simp only [Finset.mem_filter, Finset.mem_product, mul_inv_eq_one, and_true]
    exact ⟨ha, ha⟩

/-- **Each correlation is at most `|A|`** — the trivial per-cell bound (`(a,b) ↦ a` injects the fiber
into `A`, since `b = w⁻¹·a` is determined). -/
theorem mcorr_le_card (A : Finset G) (w : G) : mcorr A w ≤ A.card := by
  rw [mcorr]
  apply Finset.card_le_card_of_injOn (fun p => p.1)
  · rintro ⟨a, b⟩ hp
    simp only [Finset.mem_coe, Finset.mem_filter, Finset.mem_product] at hp
    exact hp.1.1
  · rintro ⟨a, b⟩ hp ⟨c, d⟩ hq hac
    simp only [Finset.mem_coe, Finset.mem_filter, Finset.mem_product] at hp hq
    obtain ⟨_, hab⟩ := hp; obtain ⟨_, hcd⟩ := hq
    have hac' : a = c := hac
    -- a·b⁻¹ = w = c·d⁻¹ and a = c ⟹ b = d
    have hbb : a * b⁻¹ = c * d⁻¹ := by rw [hab, hcd]
    rw [hac'] at hbb
    have hbd : b = d := inv_inj.mp (mul_left_cancel hbb)
    exact Prod.ext_iff.mpr ⟨hac', hbd⟩

/-- **The diagonal-suppression assembly anchor.**  Given a uniform per-cell correlation bound
`mcorr A w ≤ M` for every off-diagonal ratio `w ≠ 1` (the equidistribution input, e.g. from
`LRCHyperbolaBox.hyperbola_box_count`), the total off-diagonal multiplicative energy is bounded:

  `∑_{w ≠ 1} (mcorr A w)² ≤ M · (|A|² − |A|)`.

The `t₂` off-line mass is thereby controlled by the per-cell bound, the diagonal `w = 1` (mass `|A|²`)
being isolated. -/
theorem offdiag_mcorr_sq_le (A : Finset G) (M : ℕ) (hM : ∀ w, w ≠ 1 → mcorr A w ≤ M) :
    ∑ w ∈ Finset.univ.erase 1, (mcorr A w) ^ 2 ≤ M * (A.card ^ 2 - A.card) := by
  have hbound : ∑ w ∈ Finset.univ.erase 1, (mcorr A w) ^ 2
      ≤ ∑ w ∈ Finset.univ.erase 1, M * mcorr A w := by
    apply Finset.sum_le_sum
    intro w hw
    rw [Finset.mem_erase] at hw
    rw [sq]
    exact Nat.mul_le_mul_right _ (hM w hw.1)
  calc ∑ w ∈ Finset.univ.erase 1, (mcorr A w) ^ 2
      ≤ ∑ w ∈ Finset.univ.erase 1, M * mcorr A w := hbound
    _ = M * ∑ w ∈ Finset.univ.erase 1, mcorr A w := by rw [Finset.mul_sum]
    _ = M * (A.card ^ 2 - A.card) := by
        have hadd := Finset.add_sum_erase Finset.univ (mcorr A) (Finset.mem_univ 1)
        rw [sum_mcorr, mcorr_one] at hadd
        congr 1
        omega

end MultCorrelation
end LonelyRunner
