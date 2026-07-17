/-
  TournamentH7.LRCB5Race — THE BELOW-PAIR DILATE COUNT and THE B5 RACE SCOREBOARD
  (death-star-2026-07-17-S38, HYP-7182; assembles THM-940/942/943 into the formal
  race for B5-positivity on the dense core).

  * `dilatePairs_card_le` — THE COUNTING LEMMA: an injective positive 13-family
    carries at most 12 doubling pairs (a pair is determined by its top, and the
    minimum is never a top).
  * `dilate_top_below_pair` — THE TRAP CONFINEMENT: on a chain-dense family the top
    of every doubling pair sits at element index ≤ j+1 — dilate pairs live entirely
    in the bottom block (THM-939's A1 with the explicit `(1, −2)` coefficient
    vector).
  * `deviation_lower` — the proved floor `−(q−1)/7^|T| ≤ D_T` (from `N_T ≥ 0`).
  * `deviation_upper_of_mem` — the proved ceiling through THM-942A's exact
    singleton: `D_T ≤ (q−1)/7 − (q−1)/7^|T|` whenever some `i ∈ T` has
    `gcd(v i, q) = 1` and `q ≥ 14`.
  * `B5_race_scoreboard` — THE RACE, formal: the THM-940 identity with the proved
    floors (k = 1 HELPS: `−Σ D_i ≥ 0`; k = 2 and k = 4 floored trivially) and NAMED
    odd-tail hypotheses — every constant proved, the odd tails exactly what
    equidistribution owes:
        B5 ≥ (q−1)·2052/16807 − 78(q−1)/49 − 715(q−1)/2401 − tail3 − tail5.
    Sharper per-stratum floors (e.g. THM-943A's exact dilate mass at 14 ∣ q) slot
    into the k = 2 term without touching the frame.

  Kernel-pure: no `sorry`, no `native_decide`, no new axioms.
-/
import Mathlib
import TournamentH7.LRCDeviationPairs
import TournamentH7.LRCDenseCoreRelationTrap

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- The doubling pairs of a family. -/
def dilatePairs (w : Fin 13 → ℤ) : Finset (Fin 13 × Fin 13) :=
  Finset.univ.filter fun p => w p.2 = 2 * w p.1

/-- **The counting lemma**: an injective positive family carries at most 12 doubling
pairs. -/
theorem dilatePairs_card_le (w : Fin 13 → ℤ) (hpos : ∀ i, 0 < w i)
    (hinj : Function.Injective w) : (dilatePairs w).card ≤ 12 := by
  obtain ⟨m₀, -, hm₀⟩ :=
    Finset.exists_min_image (Finset.univ : Finset (Fin 13)) w ⟨0, Finset.mem_univ 0⟩
  have hmap : ∀ p ∈ dilatePairs w, p.2 ∈ Finset.univ.erase m₀ := by
    intro p hp
    rw [dilatePairs, Finset.mem_filter] at hp
    rw [Finset.mem_erase]
    refine ⟨?_, Finset.mem_univ _⟩
    intro hEq
    have h1 := hm₀ p.1 (Finset.mem_univ p.1)
    have h2 := hpos p.1
    have h3 : w p.2 = 2 * w p.1 := hp.2
    rw [hEq] at h3
    omega
  have hinjOn : ∀ p₁ ∈ dilatePairs w, ∀ p₂ ∈ dilatePairs w,
      p₁.2 = p₂.2 → p₁ = p₂ := by
    intro p₁ hp₁ p₂ hp₂ hEq
    rw [dilatePairs, Finset.mem_filter] at hp₁ hp₂
    have h1 : 2 * w p₁.1 = 2 * w p₂.1 := by
      rw [← hp₁.2, ← hp₂.2, hEq]
    have h2 : w p₁.1 = w p₂.1 := by omega
    exact Prod.ext (hinj h2) hEq
  calc (dilatePairs w).card
      ≤ (Finset.univ.erase m₀).card :=
        Finset.card_le_card_of_injOn (fun p => p.2) hmap hinjOn
    _ = 12 := by
        rw [Finset.card_erase_of_mem (Finset.mem_univ m₀)]
        simp

/-- **The trap confinement**: on a chain-dense family the top of every doubling pair
sits at element index ≤ j+1. -/
theorem dilate_top_below_pair (w : Fin 13 → ℤ)
    (hpos : ∀ i, 0 < w i) (hmono : Monotone w) (j : Fin 12)
    (hladder : ∀ k : Fin 12, j < k → 3 * w k.castSucc ≤ w k.succ)
    (s t : Fin 13) (hdil : w t = 2 * w s) :
    (t : ℕ) ≤ (j : ℕ) + 1 := by
  by_contra hcon
  push Not at hcon
  have ht2 : (j : ℕ) + 2 ≤ (t : ℕ) := by omega
  set coeff : Fin 13 → ℤ := fun x => if x = t then 1 else if x = s then -2 else 0
    with hcoeff
  have hst : s ≠ t := by
    intro hEq
    rw [hEq] at hdil
    have := hpos t
    omega
  have htop : coeff t ≠ 0 := by
    rw [hcoeff]
    simp
  have hhigh : ∀ x : Fin 13, (t : ℕ) < (x : ℕ) → coeff x = 0 := by
    intro x hx
    have hxt : x ≠ t := by
      intro hEq
      rw [hEq] at hx
      omega
    have hs_lt : (s : ℕ) < (t : ℕ) := by
      by_contra hge
      push Not at hge
      have hle : w t ≤ w s := hmono hge
      have := hpos s
      omega
    have hxs : x ≠ s := by
      intro hEq
      rw [hEq] at hx
      omega
    rw [hcoeff]
    simp [hxt, hxs]
  have hmass : (∑ x ∈ Finset.univ.filter (fun x : Fin 13 => (x : ℕ) < (t : ℕ)),
      |coeff x|) ≤ 2 := by
    have hpoint : ∀ x ∈ Finset.univ.filter (fun x : Fin 13 => (x : ℕ) < (t : ℕ)),
        |coeff x| = if x = s then 2 else 0 := by
      intro x hx
      have hxt : x ≠ t := by
        intro hEq
        have := (Finset.mem_filter.mp hx).2
        rw [hEq] at this
        omega
      by_cases hxs : x = s
      · rw [if_pos hxs, hcoeff]
        simp only [if_neg hxt, if_pos hxs]
        norm_num
      · rw [if_neg hxs, hcoeff]
        simp only [if_neg hxt, if_neg hxs]
        norm_num
    rw [Finset.sum_congr rfl hpoint, Finset.sum_ite_eq'
      (Finset.univ.filter (fun x : Fin 13 => (x : ℕ) < (t : ℕ))) s (fun _ => (2 : ℤ))]
    split <;> norm_num
  have hzero : (∑ x, coeff x * w x) = 0 := by
    have h1 : (∑ x, coeff x * w x)
        = ∑ x ∈ ({t, s} : Finset (Fin 13)), coeff x * w x := by
      symm
      apply Finset.sum_subset (Finset.subset_univ _)
      intro x _ hx
      have hxt : x ≠ t := by
        intro hEq
        exact hx (by rw [hEq]; exact Finset.mem_insert_self ..)
      have hxs : x ≠ s := by
        intro hEq
        exact hx (by
          rw [hEq]
          exact Finset.mem_insert_of_mem (Finset.mem_singleton_self _))
      rw [hcoeff]
      simp [hxt, hxs]
    rw [h1, Finset.sum_insert (by
      rw [Finset.mem_singleton]
      exact fun h => hst (h.symm))]
    rw [Finset.sum_singleton, hcoeff]
    have hts : ¬ (s = t) := hst
    simp only [if_pos rfl, hts, if_neg, if_pos rfl]
    push_cast
    linarith [hdil]
  exact LRC14Grand.no_low_mass_relation_above_pair w hpos hmono j hladder coeff t
    ht2 htop hhigh hmass hzero

/-- The proved deviation floor: `−(q−1)/7^|T| ≤ D_T` (joint counts are nonnegative). -/
theorem deviation_lower (v : Fin 13 → ℤ) (q : ℕ) (T : Finset (Fin 13)) :
    -(((q : ℚ) - 1) / 7 ^ T.card) ≤ deviation v q T := by
  unfold deviation
  have h0 : (0 : ℚ) ≤ (jointFail v q T : ℚ) := by positivity
  linarith

/-- The proved deviation ceiling through THM-942A's exact singleton: whenever some
`i ∈ T` has `gcd(v i, q) = 1` and `q ≥ 14`,
`D_T ≤ (q−1)/7 − (q−1)/7^|T|`. -/
theorem deviation_upper_of_mem (v : Fin 13 → ℤ) (q : ℕ) (T : Finset (Fin 13))
    (i : Fin 13) (hiT : i ∈ T) (hq : 14 ≤ q)
    (hgcd : Int.gcd (v i) (q : ℤ) = 1) :
    deviation v q T ≤ ((q : ℚ) - 1) / 7 - ((q : ℚ) - 1) / 7 ^ T.card := by
  have hq0 : 0 < q := by omega
  have hanti : jointFail v q T ≤ jointFail v q {i} :=
    jointFail_anti v q (Finset.singleton_subset_iff.mpr hiT)
  have hexact := jointFail_singleton_eq v q i hq0 hgcd
  have hbound := deviation_singleton_bounds v q i hq hgcd
  unfold deviation at hbound ⊢
  rw [Finset.card_singleton] at hbound
  have hle : (jointFail v q T : ℚ) ≤ (jointFail v q {i} : ℚ) := by exact_mod_cast hanti
  have h2 := hbound.2
  rw [pow_one] at h2
  linarith

/-- **THE RACE SCOREBOARD**: the THM-940 identity with the proved floors — k = 1
helps (`−Σ D_i ≥ 0`), k = 2 and k = 4 pay at worst the trivial floor — and the
NAMED odd-tail hypotheses.  Every constant proved; the odd tails are exactly what
equidistribution owes on the stratum. -/
theorem B5_race_scoreboard (v : Fin 13 → ℤ) (q : ℕ) (hq : 14 ≤ q)
    (hgcd : ∀ i, Int.gcd (v i) (q : ℤ) = 1)
    (tail3 tail5 : ℚ)
    (h3 : (∑ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard 3, deviation v q T)
      ≤ tail3)
    (h5 : (∑ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard 5, deviation v q T)
      ≤ tail5) :
    ((q : ℚ) - 1) * (2052 / 16807)
      - 78 * (((q : ℚ) - 1) / 49) - 715 * (((q : ℚ) - 1) / 2401)
      - tail3 - tail5
    ≤ (B5 v q : ℚ) := by
  have hid := B5_eq_equilibrium_add_deviation v q
  -- expand the six ledger terms
  rw [Finset.sum_range_succ, Finset.sum_range_succ, Finset.sum_range_succ,
    Finset.sum_range_succ, Finset.sum_range_succ, Finset.sum_range_succ,
    Finset.sum_range_zero] at hid
  -- k = 0: the empty subset contributes 0
  have hk0 : (∑ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard 0, deviation v q T)
      = 0 := by
    rw [Finset.powersetCard_zero, Finset.sum_singleton]
    exact deviation_empty v q (by omega)
  -- k = 1: every singleton deviation is ≤ 0, so −Σ ≥ 0
  have hk1 : (∑ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard 1, deviation v q T)
      ≤ 0 := by
    apply Finset.sum_nonpos
    intro T hT
    obtain ⟨hsub, hcard⟩ := Finset.mem_powersetCard.mp hT
    obtain ⟨i, hi⟩ := Finset.card_eq_one.mp hcard
    rw [hi]
    exact (deviation_singleton_bounds v q i hq (hgcd i)).2
  -- k = 2: the trivial floor per pair
  have hk2 : -(78 * (((q : ℚ) - 1) / 49))
      ≤ ∑ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard 2, deviation v q T := by
    have hcount : ((Finset.univ : Finset (Fin 13)).powersetCard 2).card = 78 := by
      rw [Finset.card_powersetCard]
      simp only [Finset.card_univ, Fintype.card_fin]
      rfl
    calc -(78 * (((q : ℚ) - 1) / 49))
        = ∑ _T ∈ (Finset.univ : Finset (Fin 13)).powersetCard 2,
            -(((q : ℚ) - 1) / 49) := by
          rw [Finset.sum_const, hcount, nsmul_eq_mul]
          ring
      _ ≤ _ := by
          apply Finset.sum_le_sum
          intro T hT
          have hcard := (Finset.mem_powersetCard.mp hT).2
          have := deviation_lower v q T
          rw [hcard] at this
          calc -(((q : ℚ) - 1) / 49) = -(((q : ℚ) - 1) / 7 ^ 2) := by norm_num
            _ ≤ deviation v q T := this
  -- k = 4: the trivial floor per quadruple
  have hk4 : -(715 * (((q : ℚ) - 1) / 2401))
      ≤ ∑ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard 4, deviation v q T := by
    have hcount : ((Finset.univ : Finset (Fin 13)).powersetCard 4).card = 715 := by
      rw [Finset.card_powersetCard]
      simp only [Finset.card_univ, Fintype.card_fin]
      rfl
    calc -(715 * (((q : ℚ) - 1) / 2401))
        = ∑ _T ∈ (Finset.univ : Finset (Fin 13)).powersetCard 4,
            -(((q : ℚ) - 1) / 2401) := by
          rw [Finset.sum_const, hcount, nsmul_eq_mul]
          ring
      _ ≤ _ := by
          apply Finset.sum_le_sum
          intro T hT
          have hcard := (Finset.mem_powersetCard.mp hT).2
          have := deviation_lower v q T
          rw [hcard] at this
          calc -(((q : ℚ) - 1) / 2401) = -(((q : ℚ) - 1) / 7 ^ 4) := by norm_num
            _ ≤ deviation v q T := this
  -- assemble
  rw [hid, hk0]
  push_cast
  nlinarith [hk1, hk2, hk4, h3, h5]

/-! ## Axiom audit -/
#print axioms dilatePairs_card_le
#print axioms dilate_top_below_pair
#print axioms deviation_lower
#print axioms deviation_upper_of_mem
#print axioms B5_race_scoreboard

end LRC14Concrete
end LonelyRunner
