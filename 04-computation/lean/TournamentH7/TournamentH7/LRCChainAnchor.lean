/-
  TournamentH7.LRCChainAnchor  (boxeph-2026-07-09-S12, HYP-5853 Lean layer 1)

  The 1/3 CHAIN ANCHOR — the constructive heart of chain-Bonferroni
  (HYP-5853): a doubling chain is never the obstruction, because
  `t = 1/(3a)` clears every element `a·2^e` with clearance exactly `1/3`
  (`3 ∤ 2^e`, so `|2^e/3 − m| ≥ 1/3` for every integer `m`).

  Main results:
    * `two_pow_third_clears` — `1/3 ≤ |2^e/3 − m|` for all `e : ℕ, m : ℤ`.
    * `lonely_of_pow_two_ratios` — any 13-tuple whose entries all have the
      form `a·2^e` (a single geometric chain, in any order, with repeats
      and any exponents) is `Lonely 14` at `t = 1/(3a)`, with the FULL
      clearance `1/3` (`≥ 1/14` a fortiori).  This is HYP-5853's
      `W₀ = 12` outright-lonely case, constructive — no measure theory.

  The two-chain case (`W₀ = 11`, the chain-Bonferroni frontier) needs the
  measure layer: `μ_L` as volumes of explicit rational interval unions
  (`μ₂ = 11/14`: pair-safe set `= [1/14,13/28] ∪ [15/28,13/14]`),
  dilation invariance (`LRCGoodDilation`), two-set Bonferroni, and
  positive-volume ⟹ nonempty (`LRCDensityFloorCert`).  Handoff documented
  in HYP-5853; this file supplies the constructive layer it rests on.

  Kernel-pure target: no `sorry`, no `native_decide` (audited below).
-/
import Mathlib
import TournamentH7.LonelyRunner

namespace LonelyRunner
namespace ChainAnchor

/-- `3 ∤ 2^e`: powers of two are never divisible by three. -/
theorem three_not_dvd_two_pow (e : ℕ) : ¬ (3 : ℤ) ∣ 2 ^ e := by
  intro h
  have h2 : Prime (3 : ℤ) := by norm_num
  have h3 : (3 : ℤ) ∣ 2 := h2.dvd_of_dvd_pow h
  norm_num at h3

/-- **The 1/3 clearance of the dyadic orbit:** for every `e : ℕ` and every
integer `m`, `|2^e/3 − m| ≥ 1/3` — the point `2^e/3` sits at distance
exactly `1/3` from the nearest integers, uniformly in `e`. -/
theorem two_pow_third_clears (e : ℕ) (m : ℤ) :
    (1 : ℝ) / 3 ≤ |(2 : ℝ) ^ e / 3 - m| := by
  have hne : (2 : ℤ) ^ e - 3 * m ≠ 0 := by
    intro h
    apply three_not_dvd_two_pow e
    exact ⟨m, by linarith [h]⟩
  have habs : (1 : ℤ) ≤ |(2 : ℤ) ^ e - 3 * m| := Int.one_le_abs hne
  have habsR : (1 : ℝ) ≤ |(2 : ℝ) ^ e - 3 * m| := by exact_mod_cast habs
  have hkey : |(2 : ℝ) ^ e / 3 - m| = |(2 : ℝ) ^ e - 3 * m| / 3 := by
    rw [show (2 : ℝ) ^ e / 3 - m = ((2 : ℝ) ^ e - 3 * m) / 3 by ring,
      abs_div]
    norm_num
  rw [hkey]
  linarith

/-- **The chain anchor:** any 13 speeds of the form `a·2^(e i)` (a single
geometric doubling chain, arbitrary exponents/order/multiplicity) are
`Lonely 14` at the explicit time `t = 1/(3a)` — with clearance `1/3`,
far above the `1/14` band.  HYP-5853's `W₀ = 12` case, constructive. -/
theorem lonely_of_pow_two_ratios (a : ℤ) (ha : 0 < a) (e : Fin 13 → ℕ) :
    Lonely 14 (fun i : Fin 13 => a * 2 ^ (e i)) (1 / (3 * (a : ℝ))) := by
  intro i m
  have haR : (0 : ℝ) < (a : ℝ) := by exact_mod_cast ha
  have hval : ((a * 2 ^ (e i) : ℤ) : ℝ) * (1 / (3 * (a : ℝ)))
      = (2 : ℝ) ^ (e i) / 3 := by
    push_cast
    first
    | (field_simp; ring)
    | field_simp
  rw [hval]
  have h13 := two_pow_third_clears (e i) m
  have h14 : (1 : ℝ) / 14 ≤ 1 / 3 := by norm_num
  linarith

/-- The same, packaged for subset/embedding consumers: every speed tuple
`v` with `v i = a·2^(e i)` admits a lonely time — the `∃ t` form the
sieve/dispatch layers consume. -/
theorem exists_lonely_of_chain (a : ℤ) (ha : 0 < a) (v : Fin 13 → ℤ)
    (hv : ∀ i, ∃ e : ℕ, v i = a * 2 ^ e) :
    ∃ t : ℝ, Lonely 14 v t := by
  choose e he using hv
  refine ⟨1 / (3 * (a : ℝ)), ?_⟩
  intro i m
  have h := lonely_of_pow_two_ratios a ha e i m
  rwa [show v i = a * 2 ^ (e i) from he i]

end ChainAnchor
end LonelyRunner

-- kernel-purity audit
#print axioms LonelyRunner.ChainAnchor.two_pow_third_clears
#print axioms LonelyRunner.ChainAnchor.lonely_of_pow_two_ratios
#print axioms LonelyRunner.ChainAnchor.exists_lonely_of_chain
