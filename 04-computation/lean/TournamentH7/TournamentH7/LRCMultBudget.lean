/-
# LRCMultBudget — the multiplicative supply-chain rung (M₃ deficit off the geometric extremum)

The multiplicative twin of `LRCE3Budget`.  Where the additive rung supplies `E₃ < C(k,2)` for the
residual (covering ⟹ not a dilated interval), this supplies the multiplicative deficit

  `M₃ S < C(k,2)`   for any set that is NOT a geometric power-chain,

and — the endgame wiring — shows every **covering, compressed** residual family is not a geometric chain,
hence carries the strict multiplicative deficit.  This is the extremal anchor for klein's multiplicative
character-frame route (HYP-5835/S225–226): the character sum resonates maximally on geometric progressions
(LEM-023 `multCount_eq_choose_iff_geometric`), so a residual — being off that extremum by `M₃ < C(k,2)` —
carries the multiplicative off-line deficit the character bound consumes, exactly as `E₃ < C(k,2)` feeds
the additive resonance bound.

The residual-exclusion is the clean arithmetic dual of `dilated_max_eq_card_mul_min`: a covering geometric
chain `{a,…,aᵏ}` needs `11 ∣ a` and `13 ∣ a` (covering primes ⟹ `a ≥ 143`), but a compressed one needs
`a ≤ 13` (the top element within `13×` of a smaller one) — impossible.

kind-pasteur-2026-07-09-S127.
-/
import Mathlib
import TournamentH7.LRCMultRigidity

namespace LonelyRunner
namespace LRC14Ledger

open LRCMultRigidity Finset

/-- **`M₃ ≤ C(k,2)`** — the multiplicative Schur count never exceeds the number of 2-subsets (the
injection `(a,b) ↦ {a, a·b}` lands in `powersetCard 2`). -/
theorem multCount_le_choose {S : Finset ℕ} (hmin2 : ∀ x ∈ S, 2 ≤ x) : M3 S ≤ S.card.choose 2 := by
  rw [M3, ← Finset.card_powersetCard 2 S]
  apply Finset.card_le_card_of_injOn (fun p => ({p.1, p.1 * p.2} : Finset ℕ))
  · rintro ⟨a, b⟩ hp
    simp only [Finset.mem_coe, Finset.mem_filter, Finset.mem_product] at hp
    obtain ⟨⟨ha, hb⟩, h3⟩ := hp
    have hab : a < a * b := lt_mul_of_one_lt_right (by have := hmin2 a ha; omega) (by have := hmin2 b hb; omega)
    rw [Finset.mem_coe, Finset.mem_powersetCard]
    refine ⟨fun z hz => ?_, Finset.card_pair (Nat.ne_of_lt hab)⟩
    simp only [Finset.mem_insert, Finset.mem_singleton] at hz
    rcases hz with rfl | rfl
    · exact ha
    · exact h3
  · rintro ⟨a, b⟩ hp ⟨c, d⟩ hq h
    simp only [Finset.mem_coe, Finset.mem_filter, Finset.mem_product] at hp hq
    obtain ⟨⟨ha, hb⟩, _⟩ := hp; obtain ⟨⟨hc, hd⟩, _⟩ := hq
    have hab : a < a * b := lt_mul_of_one_lt_right (by have := hmin2 a ha; omega) (by have := hmin2 b hb; omega)
    have hcd : c < c * d := lt_mul_of_one_lt_right (by have := hmin2 c hc; omega) (by have := hmin2 d hd; omega)
    simp only [Finset.ext_iff, Finset.mem_insert, Finset.mem_singleton] at h
    have hac : a = c := by
      rcases (h a).mp (Or.inl rfl) with h1 | h1
      · exact h1
      · rcases (h c).mpr (Or.inl rfl) with h2 | h2
        · exact h2.symm
        · exfalso
          have hA : c < a := h1 ▸ hcd
          have hB : a < c := h2 ▸ hab
          omega
    have hbd : b = d := by
      rcases (h (a * b)).mp (Or.inr rfl) with h1 | h1
      · exfalso; rw [← hac] at h1; omega
      · rw [hac] at h1; exact Nat.eq_of_mul_eq_mul_left (by have := hmin2 c hc; omega) h1
    exact Prod.ext hac hbd

/-- **The multiplicative deficit off the extremum (the supply-chain anchor).**  A set of integers `≥ 2`
that is NOT a geometric power-chain has strictly fewer than the maximal number of multiplicative Schur
triples: `M₃ S < C(k,2)` — combining `multCount_le_choose` with the LEM-023 equality characterisation. -/
theorem multCount_lt_choose_of_not_geometric {S : Finset ℕ} (hmin2 : ∀ x ∈ S, 2 ≤ x)
    (hng : ¬ GeometricChain S) : M3 S < S.card.choose 2 := by
  have hne : S.Nonempty := by
    rcases S.eq_empty_or_nonempty with rfl | h
    · exact absurd (⟨2, le_refl 2, by simp⟩ : GeometricChain ∅) hng
    · exact h
  have hle := multCount_le_choose hmin2
  have hne_eq : M3 S ≠ S.card.choose 2 := fun h =>
    hng ((multCount_eq_choose_iff_geometric hne hmin2).mp h)
  omega

/-- **A covering, compressed set is not a geometric chain.**  A geometric chain `{a,…,aᵏ}` that is
*covering* (`11 ∣` and `13 ∣` some elements, so `11 ∣ a` and `13 ∣ a`, forcing `a ≥ 143`) cannot be
*compressed* (the top element `aᵏ` within `13×` of a strictly smaller `aʲ` forces `a ≤ 13`).  This is the
residual-class exclusion, dual to `dilated_max_eq_card_mul_min`. -/
theorem not_geometric_of_covering_compressed {S : Finset ℕ} (hne : S.Nonempty)
    (h11 : ∃ x ∈ S, (11 : ℕ) ∣ x) (h13 : ∃ x ∈ S, (13 : ℕ) ∣ x)
    (hcomp : ∀ x ∈ S, ∃ y ∈ S, y ≠ x ∧ x ≤ 13 * y) : ¬ GeometricChain S := by
  rintro ⟨a, ha2, hSeq⟩
  set k := S.card with hk
  have hkpos : 0 < k := hne.card_pos
  have hmemS : ∀ x, x ∈ S ↔ ∃ i, 1 ≤ i ∧ i ≤ k ∧ x = a ^ i := by
    intro x
    rw [hSeq, Finset.mem_image]
    constructor
    · rintro ⟨i, hi, rfl⟩; rw [Finset.mem_Icc] at hi; exact ⟨i, hi.1, hi.2, rfl⟩
    · rintro ⟨i, hi1, hik, rfl⟩; exact ⟨i, by rw [Finset.mem_Icc]; exact ⟨hi1, hik⟩, rfl⟩
  have h11a : (11 : ℕ) ∣ a := by
    obtain ⟨x, hx, hdvd⟩ := h11
    obtain ⟨i, _, _, rfl⟩ := (hmemS x).mp hx
    exact (by norm_num : Nat.Prime 11).dvd_of_dvd_pow hdvd
  have h13a : (13 : ℕ) ∣ a := by
    obtain ⟨x, hx, hdvd⟩ := h13
    obtain ⟨i, _, _, rfl⟩ := (hmemS x).mp hx
    exact (by norm_num : Nat.Prime 13).dvd_of_dvd_pow hdvd
  have h143 : (11 * 13 : ℕ) ∣ a := (by norm_num : Nat.Coprime 11 13).mul_dvd_of_dvd_of_dvd h11a h13a
  have hage : 143 ≤ a := by have := Nat.le_of_dvd (by omega) h143; omega
  have hakmem : a ^ k ∈ S := (hmemS _).mpr ⟨k, hkpos, le_refl k, rfl⟩
  obtain ⟨y, hy, hyne, hyle⟩ := hcomp (a ^ k) hakmem
  obtain ⟨j, hj1, hjk, rfl⟩ := (hmemS y).mp hy
  have hjk' : j < k := by
    rcases lt_or_eq_of_le hjk with h | h
    · exact h
    · exact absurd (by rw [h] : a ^ j = a ^ k) hyne
  have hsplit : a ^ k = a ^ (k - j) * a ^ j := by rw [← pow_add]; congr 1; omega
  rw [hsplit] at hyle
  have hajpos : 0 < a ^ j := pow_pos (by omega) j
  have hle13 : a ^ (k - j) ≤ 13 := Nat.le_of_mul_le_mul_right hyle hajpos
  have hale : a ≤ a ^ (k - j) := Nat.le_self_pow (by omega) a
  omega

/-- **The multiplicative deficit for the residual class (the endgame wiring).**  Every covering,
compressed set of integers `≥ 2` carries the strict multiplicative Schur deficit `M₃ S < C(k,2)` — the
exact input the multiplicative character-frame off-line bound (klein's route) consumes, dual to the
additive `E3_lt_choose_of_gap`. -/
theorem multCount_lt_choose_of_covering_compressed {S : Finset ℕ} (hne : S.Nonempty)
    (hmin2 : ∀ x ∈ S, 2 ≤ x)
    (h11 : ∃ x ∈ S, (11 : ℕ) ∣ x) (h13 : ∃ x ∈ S, (13 : ℕ) ∣ x)
    (hcomp : ∀ x ∈ S, ∃ y ∈ S, y ≠ x ∧ x ≤ 13 * y) : M3 S < S.card.choose 2 :=
  multCount_lt_choose_of_not_geometric hmin2
    (not_geometric_of_covering_compressed hne h11 h13 hcomp)

end LRC14Ledger
end LonelyRunner
