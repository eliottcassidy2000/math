/-
  TournamentH7.LRCDeepCount — THE PINNED-P COUNTING and THE CENSUS CRITERION
  (death-star-2026-07-17-S42, HYP-7196; the composition THM-949 named).

  Recon truth (4000 planted ladder-7 families, q ≤ 14·v_bot): 99.5% carry ZERO
  deep multipliers; the maximum seen is 3, only at resonant q.  This module makes
  the mechanism formal:

  * `window_unique_p` — for `q ≤ 7·v`, each witness value admits AT MOST ONE bad
    multiplier: two bad p's sharing a witness force `v·(p₂−p₁) < q/7 ≤ v`, hence
    `p₂ = p₁`.  The bad set of a runner injects into its witness range `[1, v]`
    (THM-949's window) — and through THM-949's tower determinism, the whole
    7-ladder deep set injects into the BOTTOM witness `n₁ ∈ [1, v_bot]`:
    CoverageCapped-style statements on explicit strata are FINITE CHECKS.
  * `deepSet_inj_bottom` — the injection: on the window, distinct deep multipliers
    have distinct bottom witnesses; hence `#deep ≤ v_bot`.
  * `bled_ge_neg792` (decide) — the uncapped pointwise floor
    `Σ_{d≤5} (−1)^d C(c,d) ≥ −792` for every coverage `c ≤ 13` (the worst case
    `−C(12,5)` at `c = 13`), and `= if c = 0 then 1 else …` refined:
    the ledger summand is `1` at `c = 0`, `0` for `1 ≤ c ≤ 5`, and `≥ −792`
    beyond.
  * `B5_ge_live_sub_deep` — THE CENSUS CRITERION: unconditionally,
        `B5 ≥ liveCount − 792·#{p : bandCount ≥ 6}`;
    with the deep count pinned (this module + THM-949) the race closes by census:
    `B5 > 0` whenever live multipliers outnumber `792×` the deep count — and on
    the recon-typical stratum the deep count is ZERO.

  Kernel-pure: no `sorry`, no `native_decide` (decide on 14 points only).
-/
import Mathlib
import TournamentH7.LRCQWindow
import TournamentH7.LRCMomentCertificates

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- **Window uniqueness**: for `q ≤ 7·v`, two bad multipliers sharing a witness
coincide. -/
theorem window_unique_p (v : Fin 13 → ℤ) (q : ℕ) (i : Fin 13) (p₁ p₂ : ℕ)
    (hq : 0 < q) (hvpos : 0 < v i) (hwin : (q : ℤ) ≤ 7 * v i)
    (hbad₁ : ¬ inBand v q p₁ i) (hbad₂ : ¬ inBand v q p₂ i)
    (hEq : failWitness v q p₁ i = failWitness v q p₂ i) :
    p₁ = p₂ := by
  have hw₁ := bad_at_witness v q p₁ i hq hbad₁
  have hw₂ := bad_at_witness v q p₂ i hq hbad₂
  rw [hEq] at hw₁
  set n : ℤ := failWitness v q p₂ i with hn
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  -- |v·(p₁ − p₂)| ≤ |v p₁ − n q| + |n q − v p₂| < q/7 ≤ v
  have hd : 14 * |v i * (p₁ : ℤ) - v i * (p₂ : ℤ)| < 2 * (q : ℤ) := by
    have htri : |v i * (p₁ : ℤ) - v i * (p₂ : ℤ)|
        ≤ |v i * (p₁ : ℤ) - n * (q : ℤ)| + |v i * (p₂ : ℤ) - n * (q : ℤ)| := by
      have h1 : v i * (p₁ : ℤ) - v i * (p₂ : ℤ)
          = (v i * (p₁ : ℤ) - n * (q : ℤ)) - (v i * (p₂ : ℤ) - n * (q : ℤ)) := by ring
      rw [h1]
      exact abs_sub _ _
    nlinarith [hw₁, hw₂, htri]
  have habs : |v i * (p₁ : ℤ) - v i * (p₂ : ℤ)| = v i * |(p₁ : ℤ) - (p₂ : ℤ)| := by
    rw [show v i * (p₁ : ℤ) - v i * (p₂ : ℤ) = v i * ((p₁ : ℤ) - (p₂ : ℤ)) from by ring,
      abs_mul, abs_of_pos hvpos]
  rw [habs] at hd
  -- 14·v·|Δp| < 2q ≤ 14·v ⟹ |Δp| < 1 ⟹ Δp = 0
  have hlt : v i * |(p₁ : ℤ) - (p₂ : ℤ)| < v i := by nlinarith [hd, hwin]
  have hdp : |(p₁ : ℤ) - (p₂ : ℤ)| < 1 := by
    by_contra hge
    push Not at hge
    have := mul_le_mul_of_nonneg_left hge hvpos.le
    nlinarith [this]
  have : (p₁ : ℤ) = (p₂ : ℤ) := by
    have := abs_lt.mp hdp
    omega
  exact_mod_cast this

/-- **The bottom injection**: on the window `q ≤ 7·v i`, the deep set of any family
containing runner `i` among its bad runners injects into `i`'s witness range —
`#{p ∈ (0,q) : S ⊆ bad(p)} ≤ #(witness values) ≤ v i` whenever `i ∈ S`. -/
theorem deepSet_card_le (v : Fin 13 → ℤ) (q : ℕ) (i : Fin 13)
    (S : Finset (Fin 13)) (hiS : i ∈ S)
    (hq : 0 < q) (hvpos : 0 < v i) (hwin : (q : ℤ) ≤ 7 * v i) :
    (((Finset.Ioo 0 q).filter fun p => ∀ k ∈ S, ¬ inBand v q p k)).card
      ≤ (v i).toNat := by
  -- inject p ↦ witness, landing in [1, v i] ⊆ Icc 1 (v i).toNat via toNat
  have hmap : ∀ p ∈ (Finset.Ioo 0 q).filter (fun p => ∀ k ∈ S, ¬ inBand v q p k),
      (failWitness v q p i).toNat ∈ Finset.Icc 1 (v i).toNat := by
    intro p hp
    rw [Finset.mem_filter, Finset.mem_Ioo] at hp
    obtain ⟨⟨hp0, hpq⟩, hbadAll⟩ := hp
    have hbad := hbadAll i hiS
    have hwin14 : (q : ℤ) ≤ 14 * v i := by linarith
    have h1 := failWitness_pos_of_window v q p i hq hp0 hvpos hwin14 hbad
    have h2 := failWitness_le v q p i hq hpq hvpos hbad
    rw [Finset.mem_Icc]
    constructor
    · omega
    · omega
  have hinj : ∀ p₁ ∈ (Finset.Ioo 0 q).filter (fun p => ∀ k ∈ S, ¬ inBand v q p k),
      ∀ p₂ ∈ (Finset.Ioo 0 q).filter (fun p => ∀ k ∈ S, ¬ inBand v q p k),
      (failWitness v q p₁ i).toNat = (failWitness v q p₂ i).toNat → p₁ = p₂ := by
    intro p₁ hp₁ p₂ hp₂ hEq
    rw [Finset.mem_filter, Finset.mem_Ioo] at hp₁ hp₂
    have hbad₁ := hp₁.2 i hiS
    have hbad₂ := hp₂.2 i hiS
    have hwin14 : (q : ℤ) ≤ 14 * v i := by linarith
    have h1 := failWitness_pos_of_window v q p₁ i hq hp₁.1.1 hvpos hwin14 hbad₁
    have h2 := failWitness_pos_of_window v q p₂ i hq hp₂.1.1 hvpos hwin14 hbad₂
    have hEqZ : failWitness v q p₁ i = failWitness v q p₂ i := by omega
    exact window_unique_p v q i p₁ p₂ hq hvpos hwin hbad₁ hbad₂ hEqZ
  calc (((Finset.Ioo 0 q).filter fun p => ∀ k ∈ S, ¬ inBand v q p k)).card
      ≤ (Finset.Icc 1 (v i).toNat).card :=
        Finset.card_le_card_of_injOn (fun p => (failWitness v q p i).toNat) hmap hinj
    _ = (v i).toNat := by
        rw [Nat.card_Icc]
        omega

/-- The uncapped pointwise ledger floor: for every coverage `c ≤ 13`,
`Σ_{d≤5} (−1)^d C(c,d) ≥ −792`, with value `1` at `c = 0` and `0` for
`1 ≤ c ≤ 5`. -/
theorem bled_ge_neg792 :
    ∀ c ∈ Finset.range 14,
      (-792 : ℤ) ≤ ∑ d ∈ range 6, (-1 : ℤ) ^ d * Nat.choose c d := by
  decide

/-- The pointwise ledger vanishes for `1 ≤ c ≤ 5` and is `1` at `c = 0`. -/
theorem bled_eq_low :
    ∀ c ∈ Finset.range 6,
      (∑ d ∈ range 6, (-1 : ℤ) ^ d * Nat.choose c d)
        = if c = 0 then 1 else 0 := by
  decide

/-- **THE CENSUS CRITERION**: unconditionally,
`B5 ≥ liveCount − 792·#{p : 6 ≤ bandCount}`. -/
theorem B5_ge_live_sub_deep (v : Fin 13 → ℤ) (q : ℕ) :
    (liveCount v q : ℤ)
      - 792 * (((Finset.Ioo 0 q).filter fun p => 6 ≤ bandCount v q p).card : ℤ)
    ≤ B5 v q := by
  have hswap : B5 v q
      = ∑ p ∈ Finset.Ioo 0 q,
          ∑ d ∈ range 6, (-1 : ℤ) ^ d * Nat.choose (bandCount v q p) d := by
    unfold B5 momentS
    push_cast
    have hstep : ∀ d ∈ range 6,
        ((-1 : ℤ)) ^ d * (∑ p ∈ Finset.Ioo 0 q,
          (Nat.choose (bandCount v q p) d : ℤ))
        = ∑ p ∈ Finset.Ioo 0 q, (-1 : ℤ) ^ d * (Nat.choose (bandCount v q p) d : ℤ) := by
      intro d _
      rw [Finset.mul_sum]
    rw [Finset.sum_congr rfl hstep, Finset.sum_comm]
  rw [hswap]
  -- pointwise: the summand ≥ (if c = 0 then 1 else 0) − 792·(if 6 ≤ c then 1 else 0)
  have hpoint : ∀ p ∈ Finset.Ioo 0 q,
      (if bandCount v q p = 0 then (1 : ℤ) else 0)
        - 792 * (if 6 ≤ bandCount v q p then (1 : ℤ) else 0)
      ≤ ∑ d ∈ range 6, (-1 : ℤ) ^ d * Nat.choose (bandCount v q p) d := by
    intro p _
    set c := bandCount v q p with hc
    have hc13 : c ≤ 13 := bandCount_le_thirteen v q p
    by_cases h6 : 6 ≤ c
    · have hfloor := bled_ge_neg792 c (Finset.mem_range.mpr (by omega))
      have h0 : ¬ (c = 0) := by omega
      rw [if_neg h0, if_pos h6]
      linarith
    · have hlow := bled_eq_low c (Finset.mem_range.mpr (by omega))
      rw [if_neg h6, hlow]
      split <;> simp
  calc (liveCount v q : ℤ)
      - 792 * (((Finset.Ioo 0 q).filter fun p => 6 ≤ bandCount v q p).card : ℤ)
      = ∑ p ∈ Finset.Ioo 0 q,
          ((if bandCount v q p = 0 then (1 : ℤ) else 0)
            - 792 * (if 6 ≤ bandCount v q p then (1 : ℤ) else 0)) := by
        rw [Finset.sum_sub_distrib]
        congr 1
        · rw [Finset.sum_boole]
          unfold liveCount
          norm_cast
        · rw [← Finset.mul_sum, Finset.sum_boole]
    _ ≤ _ := Finset.sum_le_sum hpoint

/-- **Positivity by census**: live multipliers outnumbering `792×` the deep count
force `B5 > 0` — on recon-typical strata the deep count is ZERO and the criterion
is simply `liveCount > 0`. -/
theorem B5_pos_of_live_beats_deep (v : Fin 13 → ℤ) (q : ℕ)
    (h : 792 * (((Finset.Ioo 0 q).filter fun p => 6 ≤ bandCount v q p).card : ℤ)
      < (liveCount v q : ℤ)) :
    0 < B5 v q := by
  have := B5_ge_live_sub_deep v q
  linarith

/-! ## Axiom audit -/
#print axioms window_unique_p
#print axioms deepSet_card_le
#print axioms B5_ge_live_sub_deep
#print axioms B5_pos_of_live_beats_deep

end LRC14Concrete
end LonelyRunner
