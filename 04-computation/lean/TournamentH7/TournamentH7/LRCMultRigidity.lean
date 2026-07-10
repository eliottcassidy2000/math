/-
# LRCMultRigidity — multiplicative rank-1 geometric-chain rigidity (the dual of the E₃ rigidity)

The multiplicative twin of `LRCSchurRigidity`.  With the **multiplicative Schur count**

  `M₃ S = #{(a,b) ∈ S² : a·b ∈ S}`,

a finite set of integers `≥ 2` has `M₃ S = C(|S|,2)` — the maximum — **iff** `S` is a **geometric
power-chain** `{a, a², …, aᵏ}` (`a ≥ 2`).  This is `E₃ S = C(k,2) ⟺ dilated interval` transported through
`x ↦ log x`: the multiplicative structure of `S` is the additive structure of its exponents, so the
extremal is the exp-image of the additive dilated interval `{1,2,…,k}`.  The **doubling** case `a = 2`,
`{2,4,8,…,2ᵏ}`, is the pure 2-adic chain that klein's multiplicative character frame sees at prime rulers
and that the LRC(14) final rung's dyadic carrier lives on (see the S127 diagonal split).

The proof reuses the additive rigidity (`LRCSchurRigidity.dilated_of_closedUnderDiff`) on the exponent set
`S.image (Nat.log a)`.

kind-pasteur-2026-07-09-S127.
-/
import Mathlib
import TournamentH7.LRCSchurRigidity

namespace LRCMultRigidity

open Finset

/-- `S` is **closed under exact quotients**: `x < y` in `S` implies `x ∣ y` and `y / x ∈ S`.  The
multiplicative analogue of `ClosedUnderDiff`. -/
def MultClosedUnderQuot (S : Finset ℕ) : Prop :=
  ∀ x ∈ S, ∀ y ∈ S, x < y → x ∣ y ∧ y / x ∈ S

/-- `S` is a **geometric power-chain** `{a, a², …, aᵏ}` (`k = |S|`, `a ≥ 2`). -/
def GeometricChain (S : Finset ℕ) : Prop :=
  ∃ a : ℕ, 2 ≤ a ∧ S = (Finset.Icc 1 S.card).image (fun i => a ^ i)

/-- The **multiplicative Schur count** `#{(a,b) ∈ S² : a·b ∈ S}`. -/
def M3 (S : Finset ℕ) : ℕ := ((S ×ˢ S).filter (fun p => p.1 * p.2 ∈ S)).card

/-- **Every element of a quotient-closed set (of integers ≥ 2) is a power of its minimum.** -/
theorem pow_of_multClosed {S : Finset ℕ} (hne : S.Nonempty) (hmin2 : ∀ x ∈ S, 2 ≤ x)
    (hmc : MultClosedUnderQuot S) : ∀ y ∈ S, ∃ m, 1 ≤ m ∧ y = (S.min' hne) ^ m := by
  set a := S.min' hne with ha
  have ha2 : 2 ≤ a := hmin2 a (S.min'_mem hne)
  intro y
  induction y using Nat.strong_induction_on with
  | _ y ih =>
    intro hy
    rcases eq_or_lt_of_le (S.min'_le y hy) with heq | hlt
    · exact ⟨1, le_refl 1, by rw [← heq, pow_one]⟩
    · obtain ⟨hdvd, hquot⟩ := hmc a (S.min'_mem hne) y hy hlt
      have hya : y / a < y := Nat.div_lt_self (by have := hmin2 y hy; omega) (by omega)
      obtain ⟨m, hm1, hmeq⟩ := ih (y / a) hya hquot
      exact ⟨m + 1, by omega, by rw [pow_succ, ← hmeq, Nat.div_mul_cancel hdvd]⟩

/-- **Multiplicative rigidity: a quotient-closed set of integers ≥ 2 is a geometric power-chain.**
Proved by transporting the additive rigidity to the exponent set `E = S.image (log a)`. -/
theorem geometric_of_multClosed {S : Finset ℕ} (hne : S.Nonempty) (hmin2 : ∀ x ∈ S, 2 ≤ x)
    (hmc : MultClosedUnderQuot S) : GeometricChain S := by
  set a := S.min' hne with ha
  have ha1 : 1 < a := hmin2 a (S.min'_mem hne)
  -- every element is `a ^ (log_a of it)`
  have hpow : ∀ y ∈ S, y = a ^ (Nat.log a y) := by
    intro y hy
    obtain ⟨m, _, hm⟩ := pow_of_multClosed hne hmin2 hmc y hy
    rw [hm, Nat.log_pow ha1]
  -- the exponent set
  set E := S.image (Nat.log a) with hE
  have hlogInj : Set.InjOn (Nat.log a) S := by
    intro x hx y hy hxy
    rw [hpow x hx, hpow y hy, hxy]
  have hcardE : E.card = S.card := Finset.card_image_of_injOn hlogInj
  have hEne : E.Nonempty := hne.image _
  have h0E : 0 ∉ E := by
    rw [hE, Finset.mem_image]
    rintro ⟨y, hy, hlogy⟩
    obtain ⟨m, hm1, hm⟩ := pow_of_multClosed hne hmin2 hmc y hy
    rw [hm, Nat.log_pow ha1] at hlogy
    omega
  -- `E` is closed under differences
  have hclE : LRCSchurRigidity.ClosedUnderDiff E := by
    intro e he f hf hef
    rw [hE, Finset.mem_image] at he hf ⊢
    obtain ⟨y, hy, rfl⟩ := he
    obtain ⟨z, hz, rfl⟩ := hf
    have hyz : y < z := by
      rw [hpow y hy, hpow z hz]; exact Nat.pow_lt_pow_right ha1 hef
    obtain ⟨hdvd, hquot⟩ := hmc y hy z hz hyz
    refine ⟨z / y, hquot, ?_⟩
    have hzy : z / y = a ^ (Nat.log a z - Nat.log a y) := by
      conv_lhs => rw [hpow z hz, hpow y hy]
      exact Nat.pow_div (le_of_lt hef) (by omega)
    rw [hzy, Nat.log_pow ha1]
  -- transport the additive rigidity
  obtain ⟨d, hd0, hEeq⟩ := LRCSchurRigidity.dilated_of_closedUnderDiff hEne h0E hclE
  -- `1 ∈ E` (the minimum `a = a¹`), so `d = 1`
  have h1E : 1 ∈ E := by
    rw [hE, Finset.mem_image]
    refine ⟨a, S.min'_mem hne, ?_⟩
    have hmem := hpow a (S.min'_mem hne)
    have h1 : a ^ 1 = a ^ Nat.log a a := by rw [pow_one]; exact hmem
    exact (Nat.pow_right_injective ha1 h1).symm
  have hd1 : d = 1 := by
    rw [hEeq, Finset.mem_image] at h1E
    obtain ⟨i, _, hid⟩ := h1E
    have hdvd1 : d ∣ 1 := ⟨i, by rw [mul_comm]; exact hid.symm⟩
    have := Nat.le_of_dvd Nat.one_pos hdvd1
    omega
  have hEeq' : E = Finset.Icc 1 S.card := by
    rw [hEeq, hd1, hcardE]; simp
  -- assemble `S = (Icc 1 k).image (a ^ ·)`
  refine ⟨a, ha1, ?_⟩
  have himg : S.image (fun y => a ^ (Nat.log a y)) = S := by
    have hcong : S.image (fun y => a ^ (Nat.log a y)) = S.image id :=
      Finset.image_congr (fun y hy => by rw [id_eq]; exact (hpow y hy).symm)
    rw [hcong, Finset.image_id]
  calc S = S.image (fun y => a ^ (Nat.log a y)) := himg.symm
    _ = (S.image (Nat.log a)).image (fun i => a ^ i) := by rw [Finset.image_image]; rfl
    _ = (Finset.Icc 1 S.card).image (fun i => a ^ i) := by rw [← hE, hEeq']

/-- **Converse: a geometric power-chain is quotient-closed.** -/
theorem multClosed_of_geometric {S : Finset ℕ} (h : GeometricChain S) : MultClosedUnderQuot S := by
  obtain ⟨a, ha2, hSeq⟩ := h
  have ha1 : 1 < a := ha2
  intro x hx y hy hxy
  rw [hSeq] at hx hy ⊢
  simp only [Finset.mem_image, Finset.mem_Icc] at hx hy ⊢
  obtain ⟨i, ⟨hi1, hik⟩, rfl⟩ := hx
  obtain ⟨j, ⟨hj1, hjk⟩, rfl⟩ := hy
  have hij : i < j := by
    by_contra hcon; push_neg at hcon
    exact absurd hxy (not_lt.mpr (Nat.pow_le_pow_right (le_of_lt ha1) hcon))
  refine ⟨pow_dvd_pow a (le_of_lt hij), j - i, ⟨by omega, by omega⟩, ?_⟩
  rw [Nat.pow_div (le_of_lt hij) (by omega)]

/-- **The multiplicative Schur count is maximal iff the set is quotient-closed.**  `M₃ S = C(k,2)` — the
maximum — iff every 2-subset `{x,y}` (`x<y`) is realised as `{x, x·(y/x)}` with `x ∣ y`, `y/x ∈ S`, i.e.
`MultClosedUnderQuot`.  The multiplicative twin of `schurCount_eq_choose_iff_closedUnderDiff`, via the
injection `(a,b) ↦ {a, a·b}`. -/
theorem multCount_eq_choose_iff_multClosed {S : Finset ℕ} (hmin2 : ∀ x ∈ S, 2 ≤ x) :
    M3 S = S.card.choose 2 ↔ MultClosedUnderQuot S := by
  set P : Finset (ℕ × ℕ) := (S ×ˢ S).filter (fun p => p.1 * p.2 ∈ S) with hPdef
  have hmem : ∀ p ∈ P, p.1 ∈ S ∧ p.2 ∈ S ∧ p.1 * p.2 ∈ S := by
    intro p hp
    simp only [hPdef, Finset.mem_filter, Finset.mem_product] at hp
    exact ⟨hp.1.1, hp.1.2, hp.2⟩
  have hlt : ∀ p ∈ P, p.1 < p.1 * p.2 := by
    intro p hp
    obtain ⟨h1, h2, _⟩ := hmem p hp
    exact lt_mul_of_one_lt_right (by have := hmin2 _ h1; omega) (by have := hmin2 _ h2; omega)
  set φ : ℕ × ℕ → Finset ℕ := fun p => ({p.1, p.1 * p.2} : Finset ℕ) with hφ
  have hsub : P.image φ ⊆ S.powersetCard 2 := by
    intro T hT
    rw [Finset.mem_image] at hT
    obtain ⟨p, hp, rfl⟩ := hT
    obtain ⟨h1, _, h3⟩ := hmem p hp
    rw [Finset.mem_powersetCard]
    refine ⟨?_, Finset.card_pair (Nat.ne_of_lt (hlt p hp))⟩
    intro z hz
    simp only [hφ, Finset.mem_insert, Finset.mem_singleton] at hz
    rcases hz with rfl | rfl
    · exact h1
    · exact h3
  have hinj : Set.InjOn φ P := by
    intro p hp q hq h
    have hpl := hlt p hp; have hql := hlt q hq
    simp only [hφ, Finset.ext_iff, Finset.mem_insert, Finset.mem_singleton] at h
    have hp1 : p.1 = q.1 := by
      rcases (h p.1).mp (Or.inl rfl) with h1 | h1
      · exact h1
      · rcases (h q.1).mpr (Or.inl rfl) with h2 | h2
        · exact h2.symm
        · exfalso
          have hA : q.1 < p.1 := h1 ▸ hql
          have hB : p.1 < q.1 := h2 ▸ hpl
          omega
    have hp2 : p.2 = q.2 := by
      rcases (h (p.1 * p.2)).mp (Or.inr rfl) with h1 | h1
      · exfalso; rw [← hp1] at h1; omega
      · rw [hp1] at h1; exact Nat.eq_of_mul_eq_mul_left (by have := hmin2 q.1 (hmem q hq).1; omega) h1
    exact Prod.ext hp1 hp2
  have hcardI : (P.image φ).card = M3 S := by
    rw [M3, ← hPdef]; exact Finset.card_image_of_injOn hinj
  have hcardPC : (S.powersetCard 2).card = S.card.choose 2 := Finset.card_powersetCard 2 S
  constructor
  · intro hE
    have hIeq : P.image φ = S.powersetCard 2 :=
      Finset.eq_of_subset_of_card_le hsub (by rw [hcardI, hcardPC, hE])
    intro x hx y hy hxy
    have hTpc : ({x, y} : Finset ℕ) ∈ S.powersetCard 2 := by
      rw [Finset.mem_powersetCard]
      refine ⟨?_, Finset.card_pair (Nat.ne_of_lt hxy)⟩
      intro z hz
      simp only [Finset.mem_insert, Finset.mem_singleton] at hz
      rcases hz with rfl | rfl
      · exact hx
      · exact hy
    rw [← hIeq, Finset.mem_image] at hTpc
    obtain ⟨p, hp, hpφ⟩ := hTpc
    have hpl := hlt p hp
    obtain ⟨_, hp2S, _⟩ := hmem p hp
    simp only [hφ, Finset.ext_iff, Finset.mem_insert, Finset.mem_singleton] at hpφ
    have hp1x : p.1 = x := by
      rcases (hpφ x).mpr (Or.inl rfl) with h1 | h1
      · exact h1.symm
      · exfalso
        rcases (hpφ y).mpr (Or.inr rfl) with h2 | h2
        · rw [h1, h2] at hxy; omega
        · rw [h1, h2] at hxy; omega
    have hp2y : p.1 * p.2 = y := by
      rcases (hpφ y).mpr (Or.inr rfl) with h1 | h1
      · exfalso; rw [hp1x] at h1; omega
      · exact h1.symm
    have hx0 : 0 < x := by have := hmin2 x hx; omega
    refine ⟨⟨p.2, by rw [← hp2y, hp1x]⟩, ?_⟩
    have hyx : y / x = p.2 := by rw [← hp2y, hp1x, Nat.mul_div_cancel_left _ hx0]
    rw [hyx]; exact hp2S
  · intro hcl
    have hsup : S.powersetCard 2 ⊆ P.image φ := by
      intro T hT
      rw [Finset.mem_powersetCard] at hT
      obtain ⟨hTS, hT2⟩ := hT
      obtain ⟨a, b, hab, rfl⟩ := Finset.card_eq_two.mp hT2
      have haS : a ∈ S := hTS (Finset.mem_insert_self a {b})
      have hbS : b ∈ S := hTS (Finset.mem_insert_of_mem (Finset.mem_singleton_self b))
      rw [Finset.mem_image]
      rcases Nat.lt_or_ge a b with hlt' | hge
      · obtain ⟨hdvd, hquot⟩ := hcl a haS b hbS hlt'
        refine ⟨(a, b / a), ?_, ?_⟩
        · simp only [hPdef, Finset.mem_filter, Finset.mem_product]
          exact ⟨⟨haS, hquot⟩, by rw [Nat.mul_div_cancel' hdvd]; exact hbS⟩
        · simp only [hφ]; rw [Nat.mul_div_cancel' hdvd]
      · have hlt' : b < a := lt_of_le_of_ne hge (Ne.symm hab)
        obtain ⟨hdvd, hquot⟩ := hcl b hbS a haS hlt'
        refine ⟨(b, a / b), ?_, ?_⟩
        · simp only [hPdef, Finset.mem_filter, Finset.mem_product]
          exact ⟨⟨hbS, hquot⟩, by rw [Nat.mul_div_cancel' hdvd]; exact haS⟩
        · simp only [hφ]; rw [Nat.mul_div_cancel' hdvd, Finset.pair_comm]
    have hIeq : P.image φ = S.powersetCard 2 := Finset.Subset.antisymm hsub hsup
    rw [← hcardI, hIeq, hcardPC]

/-- **Multiplicative rank-1 rigidity (the dual of `schurCount_eq_choose_iff_dilated`).**  For a finite
set of integers `≥ 2`, the multiplicative Schur count is maximal — `M₃ S = C(k,2)` — **iff** `S` is a
geometric power-chain `{a, a², …, aᵏ}`.  The `a = 2` instance `{2,4,…,2ᵏ}` is klein's ratio-2 doubling
chain — the 2-adic carrier of the LRC(14) final rung. -/
theorem multCount_eq_choose_iff_geometric {S : Finset ℕ} (hne : S.Nonempty) (hmin2 : ∀ x ∈ S, 2 ≤ x) :
    M3 S = S.card.choose 2 ↔ GeometricChain S := by
  rw [multCount_eq_choose_iff_multClosed hmin2]
  exact ⟨geometric_of_multClosed hne hmin2, multClosed_of_geometric⟩

end LRCMultRigidity
