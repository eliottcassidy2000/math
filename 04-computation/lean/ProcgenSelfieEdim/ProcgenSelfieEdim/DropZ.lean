import ProcgenSelfieEdim.CollatzDrop
import ProcgenSelfieEdim.ListLemmas
import ProcgenSelfieEdim.PolyNorm

set_option autoImplicit false
set_option linter.unusedSimpArgs false

/-!
# THM-4530 (a): Syracuse drops over all odd integers — "two copies over Z" (Prop. 2.6)

Extend the Syracuse map and the drop to all odd integers `A` (both signs), keeping the sign
in the odd part: `3A + 1 = 2^v · S` with `S = syrZ A` odd (of the sign of `3A + 1`), and
`K(A) = dropZ A = (A - S) / 2`.

* `dropZ_identity`: `6 K(A) + 1 = (2^v - 3) S(A)` for every odd `A ∈ Z`.
* `dropZ_unit_pos`, `dropZ_unit_neg`: every `d ∈ Z` is the drop of `8d + 1` (branch `v = 2`,
  `S = 6d + 1`) and of `-4d - 1` (branch `v = 1`, `S = -6d - 1`); the two unit copies have
  opposite signs (`unit_copies_opposite_signs`).
* `dropZ_fibre` (Prop. 2.6): the odd `A ∈ Z` with drop `d` are exactly
  `A = (1 + 2^(v+1) d) / (2^v - 3)` over the `v ≥ 1` with `(2^v - 3) ∣ 6d + 1`; distinct `v`
  give distinct `A` (`dropZ_fibre_inj`).
* `fibreZ_spec`, `fibreZ_length`: a duplicate-free list of exactly these `A` has length
  `2 + N`, where `N = #{v ≥ 3 : (2^v - 3) ∣ 6d + 1}` (`extraZ`, independent of the cutoff by
  `extraZ_stable`). So every `d ∈ Z` is a drop exactly `2 + N(6d+1)` times.
* `extra_copy_sign`: the extra copies (`v ≥ 3`) have the sign of `6d + 1`.
* `labels_two_copies` (owner's labels `M = (A + 1)/2`): among the preimages of `d`, the
  even label is exactly `-2d` (`A = -4d - 1`), the label `≡ 1 (mod 4)` is exactly `4d + 1`
  (`A = 8d + 1`), and all extra copies have labels `≡ 3 (mod 4)`.
-/

namespace ProcgenSelfieEdim

/-! ## Signed odd part -/

/-- The 2-adic valuation of a nonzero integer. -/
def v2Z (n : Int) : Nat := v2 n.natAbs

/-- The signed odd part: `n = 2^(v2Z n) · oddPartZ n`, with `oddPartZ n` odd and of the
sign of `n`. -/
def oddPartZ (n : Int) : Int :=
  if n < 0 then -((oddPart n.natAbs : Nat) : Int) else ((oddPart n.natAbs : Nat) : Int)

theorem decompZ (n : Int) (hn : n ≠ 0) : n = 2 ^ v2Z n * oddPartZ n ∧ oddPartZ n % 2 = 1 := by
  have hpos : 0 < n.natAbs := Int.natAbs_pos.2 hn
  obtain ⟨h1, h2⟩ := decomp n.natAbs hpos
  have hc : (n.natAbs : Int) = 2 ^ v2 n.natAbs * ((oddPart n.natAbs : Nat) : Int) := by
    have e := congrArg (fun x : Nat => (x : Int)) h1
    simp only [Int.natCast_mul, Int.natCast_pow] at e
    exact e
  unfold v2Z oddPartZ
  by_cases hneg : n < 0
  · rw [if_pos hneg]
    have habs : (n.natAbs : Int) = -n := Int.ofNat_natAbs_of_nonpos (by omega)
    constructor
    · rw [Int.mul_neg, ← hc]
      omega
    · omega
  · rw [if_neg hneg]
    have habs : (n.natAbs : Int) = n := Int.natAbs_of_nonneg (by omega)
    constructor
    · rw [← hc]
      omega
    · omega

/-- Uniqueness of the signed decomposition `n = 2^v · s` with `s` odd. -/
theorem v2Z_oddPartZ_of_eq (n : Int) (v : Nat) (s : Int) (h : n = 2 ^ v * s) (hs : s % 2 = 1) :
    v2Z n = v ∧ oddPartZ n = s := by
  have habs : n.natAbs = 2 ^ v * s.natAbs := by
    rw [h, Int.natAbs_mul, Int.natAbs_pow]
    rfl
  have hsodd : s.natAbs % 2 = 1 := by omega
  obtain ⟨hv, ho⟩ := v2_oddPart_of_eq _ _ _ habs hsodd
  refine ⟨hv, ?_⟩
  unfold oddPartZ
  rw [ho]
  have hp : (0 : Int) < 2 ^ v := Int.pow_pos (by decide)
  by_cases hs0 : s < 0
  · have hn : n < 0 := by
      rw [h]
      exact Int.mul_neg_of_pos_of_neg hp hs0
    rw [if_pos hn]
    omega
  · have hn : ¬ n < 0 := by
      rw [h]
      exact Int.not_lt.2 (Int.mul_nonneg (Int.le_of_lt hp) (by omega))
    rw [if_neg hn]
    omega

/-! ## Syracuse map and drop on all odd integers -/

/-- The Syracuse map on odd integers: `S(A) = oddpart(3A + 1)`, sign kept. -/
def syrZ (A : Int) : Int := oddPartZ (3 * A + 1)

/-- The drop `K(A) = (A - S(A)) / 2`. -/
def dropZ (A : Int) : Int := (A - syrZ A) / 2

theorem syrZ_decomp (A : Int) (hA : A % 2 = 1) :
    3 * A + 1 = 2 ^ v2Z (3 * A + 1) * syrZ A ∧ syrZ A % 2 = 1 :=
  decompZ (3 * A + 1) (by omega)

theorem v2Z_pos (A : Int) (hA : A % 2 = 1) : 1 ≤ v2Z (3 * A + 1) := by
  obtain ⟨h1, h2⟩ := syrZ_decomp A hA
  cases hv : v2Z (3 * A + 1) with
  | zero =>
    rw [hv, Int.pow_zero, Int.one_mul] at h1
    omega
  | succ w => omega

/-- **The drop identity over Z.** `6 K(A) + 1 = (2^v - 3) S(A)`, `v = v_2(3A + 1)`. -/
theorem dropZ_identity (A : Int) (hA : A % 2 = 1) :
    6 * dropZ A + 1 = ((2 : Int) ^ v2Z (3 * A + 1) - 3) * syrZ A := by
  obtain ⟨h1, h2⟩ := syrZ_decomp A hA
  unfold dropZ
  rw [Int.sub_mul]
  omega

/-- Every admissible branch `v ≥ 1` gives an odd preimage, `A = 2d + (6d + 1)/(2^v - 3)`,
with `(2^v - 3) A = 1 + 2^(v+1) d`. -/
theorem dropZ_backward (d : Int) (v : Nat) (hv : 1 ≤ v) (hdvd : ((2 : Int) ^ v - 3) ∣ (6 * d + 1)) :
    ∃ A : Int, A % 2 = 1 ∧ v2Z (3 * A + 1) = v ∧ syrZ A = (6 * d + 1) / ((2 : Int) ^ v - 3) ∧
      dropZ A = d ∧ ((2 : Int) ^ v - 3) * A = 1 + 2 ^ (v + 1) * d := by
  have hq := Int.mul_ediv_cancel' hdvd
  generalize hs : (6 * d + 1) / ((2 : Int) ^ v - 3) = s at hq
  obtain ⟨w, rfl⟩ : ∃ w, v = w + 1 := ⟨v - 1, by omega⟩
  have hpow : (2 : Int) ^ (w + 1) = 2 * (2 : Int) ^ w := by rw [Int.pow_succ, Int.mul_comm]
  have hpow2 : (2 : Int) ^ (w + 1 + 1) = 2 * (2 * (2 : Int) ^ w) := by
    rw [Int.pow_succ, hpow, Int.mul_comm]
  rw [hpow] at hq
  have hq' : 2 * ((2 : Int) ^ w * s) - 3 * s = 6 * d + 1 := by
    rw [← hq, Int.sub_mul, Int.mul_assoc]
  have hsodd : s % 2 = 1 := by omega
  refine ⟨2 * d + s, by omega, ?_⟩
  have key : 3 * (2 * d + s) + 1 = 2 ^ (w + 1) * s := by
    rw [hpow, Int.mul_assoc]
    omega
  obtain ⟨hv2, hodd⟩ := v2Z_oddPartZ_of_eq _ _ _ key hsodd
  refine ⟨hv2, ?_, ?_, ?_⟩
  · unfold syrZ
    exact hodd
  · unfold dropZ syrZ
    rw [hodd]
    omega
  · rw [hpow, hpow2]
    generalize (2 : Int) ^ w = X at hq ⊢
    pnorm at hq ⊢
    omega

/-- Forward direction: the branch `v = v_2(3A + 1)` of an odd `A` is admissible for its drop. -/
theorem dropZ_forward (A : Int) (hA : A % 2 = 1) :
    1 ≤ v2Z (3 * A + 1) ∧ ((2 : Int) ^ v2Z (3 * A + 1) - 3) ∣ (6 * dropZ A + 1) ∧
      ((2 : Int) ^ v2Z (3 * A + 1) - 3) * A = 1 + 2 ^ (v2Z (3 * A + 1) + 1) * dropZ A := by
  obtain ⟨h1, h2⟩ := syrZ_decomp A hA
  have hid := dropZ_identity A hA
  refine ⟨v2Z_pos A hA, ⟨syrZ A, hid⟩, ?_⟩
  have hAd : A = 2 * dropZ A + syrZ A := by
    unfold dropZ
    omega
  have hpow2 : (2 : Int) ^ (v2Z (3 * A + 1) + 1) = 2 * (2 : Int) ^ v2Z (3 * A + 1) := by
    rw [Int.pow_succ, Int.mul_comm]
  generalize hX : (2 : Int) ^ v2Z (3 * A + 1) = X at hid h1 hpow2 ⊢
  rw [hpow2]
  conv => lhs; rw [hAd]
  pnorm at hid ⊢
  omega

/-- A nonzero `M` divides `6d + 1` only for `|M| ≤ |6d + 1|`. -/
theorem two_pow_lt_of_dvd (d : Int) (v V : Nat) (hd1 : 6 * d + 1 < (2 : Int) ^ V - 3)
    (hd2 : -(6 * d + 1) < (2 : Int) ^ V - 3) (hdvd : ((2 : Int) ^ v - 3) ∣ (6 * d + 1)) : v < V := by
  by_cases hvV : V ≤ v
  · exfalso
    have hmono := two_pow_mono V v hvV
    by_cases hpos : 0 < 6 * d + 1
    · have := Int.le_of_dvd hpos hdvd
      omega
    · have hneg : 0 < -(6 * d + 1) := by omega
      have := Int.le_of_dvd hneg (Int.dvd_neg.2 hdvd)
      omega
  · omega

/-! ## The fibre of a drop value -/

/-- **Prop. 2.6.** The odd `A ∈ Z` with drop `d` are exactly
`A = (1 + 2^(v+1) d) / (2^v - 3)` over the `v ≥ 1` with `(2^v - 3) ∣ 6d + 1`. -/
theorem dropZ_fibre (d A : Int) : (A % 2 = 1 ∧ dropZ A = d) ↔
    ∃ v : Nat, 1 ≤ v ∧ ((2 : Int) ^ v - 3) ∣ (6 * d + 1) ∧
      A = (1 + 2 ^ (v + 1) * d) / ((2 : Int) ^ v - 3) := by
  constructor
  · rintro ⟨hA, rfl⟩
    obtain ⟨f1, f2, f3⟩ := dropZ_forward A hA
    refine ⟨v2Z (3 * A + 1), f1, f2, ?_⟩
    rw [← f3, Int.mul_ediv_cancel_left _ (two_pow_sub_three_ne_zero _)]
  · rintro ⟨v, hv, hdvd, rfl⟩
    obtain ⟨A, hA1, _, _, hA4, hA5⟩ := dropZ_backward d v hv hdvd
    rw [← hA5, Int.mul_ediv_cancel_left _ (two_pow_sub_three_ne_zero _)]
    exact ⟨hA1, hA4⟩

/-- Distinct branches give distinct preimages: the branch of a preimage is `v_2(3A + 1)`. -/
theorem dropZ_fibre_inj (d : Int) (v : Nat) (hv : 1 ≤ v) (hdvd : ((2 : Int) ^ v - 3) ∣ (6 * d + 1)) :
    v2Z (3 * ((1 + 2 ^ (v + 1) * d) / ((2 : Int) ^ v - 3)) + 1) = v := by
  obtain ⟨A, _, hA2, _, _, hA5⟩ := dropZ_backward d v hv hdvd
  rw [← hA5, Int.mul_ediv_cancel_left _ (two_pow_sub_three_ne_zero _)]
  exact hA2

/-- **The unit copy `8d + 1`** (branch `v = 2`, `2^2 - 3 = 1`). -/
theorem dropZ_unit_pos (d : Int) :
    (8 * d + 1) % 2 = 1 ∧ v2Z (3 * (8 * d + 1) + 1) = 2 ∧ syrZ (8 * d + 1) = 6 * d + 1 ∧
      dropZ (8 * d + 1) = d := by
  have key : 3 * (8 * d + 1) + 1 = 2 ^ 2 * (6 * d + 1) := by
    show _ = 4 * (6 * d + 1)
    omega
  obtain ⟨h1, h2⟩ := v2Z_oddPartZ_of_eq _ _ _ key (by omega)
  refine ⟨by omega, h1, h2, ?_⟩
  unfold dropZ syrZ
  rw [h2]
  omega

/-- **The unit copy `-4d - 1`** (branch `v = 1`, `2^1 - 3 = -1`). -/
theorem dropZ_unit_neg (d : Int) :
    (-4 * d - 1) % 2 = 1 ∧ v2Z (3 * (-4 * d - 1) + 1) = 1 ∧ syrZ (-4 * d - 1) = -6 * d - 1 ∧
      dropZ (-4 * d - 1) = d := by
  have key : 3 * (-4 * d - 1) + 1 = 2 ^ 1 * (-6 * d - 1) := by
    show _ = 2 * (-6 * d - 1)
    omega
  obtain ⟨h1, h2⟩ := v2Z_oddPartZ_of_eq _ _ _ key (by omega)
  refine ⟨by omega, h1, h2, ?_⟩
  unfold dropZ syrZ
  rw [h2]
  omega

/-- The two unit copies have opposite signs. -/
theorem unit_copies_opposite_signs (d : Int) :
    (0 < 8 * d + 1 ∧ -4 * d - 1 < 0) ∨ (8 * d + 1 < 0 ∧ 0 < -4 * d - 1) := by
  by_cases h : 0 ≤ d
  · exact Or.inl ⟨by omega, by omega⟩
  · exact Or.inr ⟨by omega, by omega⟩

/-- The extra copies (`v ≥ 3`) have the sign of `6d + 1`. -/
theorem extra_copy_sign (d : Int) (v : Nat) (hv : 3 ≤ v) (hdvd : ((2 : Int) ^ v - 3) ∣ (6 * d + 1)) :
    (0 < 6 * d + 1 → 0 < (1 + 2 ^ (v + 1) * d) / ((2 : Int) ^ v - 3)) ∧
      (6 * d + 1 < 0 → (1 + 2 ^ (v + 1) * d) / ((2 : Int) ^ v - 3) < 0) := by
  obtain ⟨A, _, _, _, _, hA5⟩ := dropZ_backward d v (by omega) hdvd
  rw [← hA5, Int.mul_ediv_cancel_left _ (two_pow_sub_three_ne_zero _)]
  have h8 := two_pow_mono 3 v hv
  have h8' : (2 : Int) ^ 3 = 8 := by decide
  have hpow2 : (2 : Int) ^ (v + 1) = 2 * (2 : Int) ^ v := by rw [Int.pow_succ, Int.mul_comm]
  rw [hpow2] at hA5
  generalize hX : (2 : Int) ^ v = X at hA5 h8
  constructor
  · intro hpos
    by_cases hA : 0 < A
    · exact hA
    · exfalso
      have hd : 0 ≤ d := by omega
      have h1 : (X - 3) * A ≤ 0 := Int.mul_nonpos_of_nonneg_of_nonpos (by omega) (by omega)
      have h2 : 0 ≤ X * d := Int.mul_nonneg (by omega) hd
      rw [Int.mul_assoc] at hA5
      omega
  · intro hneg
    by_cases hA : A < 0
    · exact hA
    · exfalso
      have hd : d ≤ -1 := by omega
      have h1 : 0 ≤ (X - 3) * A := Int.mul_nonneg (by omega) (by omega)
      have h2 : X * d ≤ X * (-1) := Int.mul_le_mul_of_nonneg_left hd (by omega)
      rw [Int.mul_assoc] at hA5
      omega

/-! ## Counting: exactly `2 + N(6d + 1)` copies -/

/-- The branch test `v ≥ 1 ∧ (2^v - 3) ∣ 6d + 1` as a Boolean. -/
def branchB (d : Int) (v : Nat) : Bool :=
  decide (1 ≤ v) && decide ((6 * d + 1) % ((2 : Int) ^ v - 3) = 0)

/-- The extra-branch test `v ≥ 3 ∧ (2^v - 3) ∣ 6d + 1`. -/
def extraB (d : Int) (v : Nat) : Bool :=
  decide (3 ≤ v) && decide ((6 * d + 1) % ((2 : Int) ^ v - 3) = 0)

/-- The preimages of `d` from the branches `v < V`. -/
def fibreZ (d : Int) (V : Nat) : List Int :=
  ((List.range V).filter (branchB d)).map (fun v => (1 + 2 ^ (v + 1) * d) / ((2 : Int) ^ v - 3))

/-- `N(6d + 1)` truncated at `V`: the number of `3 ≤ v < V` with `(2^v - 3) ∣ 6d + 1`. -/
def extraZ (d : Int) (V : Nat) : Nat := ((List.range V).filter (extraB d)).length

theorem branchB_iff (d : Int) (v : Nat) :
    branchB d v = true ↔ 1 ≤ v ∧ ((2 : Int) ^ v - 3) ∣ (6 * d + 1) := by
  unfold branchB
  rw [Bool.and_eq_true, decide_eq_true_iff, decide_eq_true_iff, Int.dvd_iff_emod_eq_zero]

/-- **The fibre, as a list.** If `|6d + 1| < 2^V - 3`, then `fibreZ d V` is duplicate-free and
contains exactly the odd `A ∈ Z` with `K(A) = d`. -/
theorem fibreZ_spec (d : Int) (V : Nat) (hd1 : 6 * d + 1 < (2 : Int) ^ V - 3)
    (hd2 : -(6 * d + 1) < (2 : Int) ^ V - 3) :
    (fibreZ d V).Nodup ∧ ∀ A, A ∈ fibreZ d V ↔ (A % 2 = 1 ∧ dropZ A = d) := by
  constructor
  · unfold fibreZ
    apply nodup_map_of_inj_on _ _ (List.Nodup.sublist List.filter_sublist (nodup_range V))
    intro a ha b hb hab
    have ha' := (branchB_iff d a).1 (List.mem_filter.1 ha).2
    have hb' := (branchB_iff d b).1 (List.mem_filter.1 hb).2
    have e1 := dropZ_fibre_inj d a ha'.1 ha'.2
    have e2 := dropZ_fibre_inj d b hb'.1 hb'.2
    simp only at hab
    rw [hab] at e1
    rw [← e1, e2]
  · intro A
    rw [dropZ_fibre]
    unfold fibreZ
    rw [List.mem_map]
    constructor
    · rintro ⟨v, hv, rfl⟩
      rw [List.mem_filter, branchB_iff] at hv
      exact ⟨v, hv.2.1, hv.2.2, rfl⟩
    · rintro ⟨v, hv, hdvd, rfl⟩
      refine ⟨v, ?_, rfl⟩
      rw [List.mem_filter, branchB_iff, List.mem_range]
      exact ⟨two_pow_lt_of_dvd d v V hd1 hd2 hdvd, hv, hdvd⟩

theorem branch_one_two (d : Int) : branchB d 1 = true ∧ branchB d 2 = true := by
  constructor
  · rw [branchB_iff]
    refine ⟨Nat.le_refl 1, ?_⟩
    have : (2 : Int) ^ 1 - 3 = -1 := by decide
    rw [this]
    exact ⟨-(6 * d + 1), by omega⟩
  · rw [branchB_iff]
    refine ⟨by decide, ?_⟩
    have : (2 : Int) ^ 2 - 3 = 1 := by decide
    rw [this]
    exact Int.one_dvd _

theorem length_filter_range_succ (p : Nat → Bool) (n : Nat) :
    ((List.range (n + 1)).filter p).length = ((List.range n).filter p).length + (if p n then 1 else 0) := by
  rw [List.range_succ, List.filter_append, List.length_append]
  cases h : p n
  · rw [List.filter_cons_of_neg (by rw [h]; decide), List.filter_nil]
    rfl
  · rw [List.filter_cons_of_pos h, List.filter_nil]
    rfl

/-- **Exactly `2 + N(6d + 1)` copies.** For `V ≥ 3`, the fibre list has length `2 + extraZ d V`. -/
theorem fibreZ_length (d : Int) : ∀ (V : Nat), 3 ≤ V → (fibreZ d V).length = 2 + extraZ d V
  | 0, h => absurd h (by decide)
  | 1, h => absurd h (by decide)
  | 2, h => absurd h (by decide)
  | 3, _ => by
    unfold fibreZ extraZ
    rw [List.length_map]
    have h0 : branchB d 0 = false := rfl
    obtain ⟨h1, h2⟩ := branch_one_two d
    have e0 : extraB d 0 = false := rfl
    have e1 : extraB d 1 = false := rfl
    have e2 : extraB d 2 = false := rfl
    rw [length_filter_range_succ, length_filter_range_succ, length_filter_range_succ,
      length_filter_range_succ, length_filter_range_succ, length_filter_range_succ, h0, h1, h2,
      e0, e1, e2]
    rfl
  | V + 4, _ => by
    have ih := fibreZ_length d (V + 3) (by omega)
    unfold fibreZ extraZ at ih ⊢
    rw [List.length_map] at ih ⊢
    rw [length_filter_range_succ, length_filter_range_succ (extraB d)]
    have hb : branchB d (V + 3) = extraB d (V + 3) := by
      unfold branchB extraB
      have h1 : decide (1 ≤ V + 3) = true := decide_eq_true (by omega)
      have h3 : decide (3 ≤ V + 3) = true := decide_eq_true (by omega)
      rw [h1, h3]
    rw [hb, ih]
    omega

/-- The count does not depend on the cutoff once `2^V - 3 > |6d + 1|`. -/
theorem extraZ_stable (d : Int) (V : Nat) (hd1 : 6 * d + 1 < (2 : Int) ^ V - 3)
    (hd2 : -(6 * d + 1) < (2 : Int) ^ V - 3) : ∀ k, extraZ d (V + k) = extraZ d V
  | 0 => rfl
  | k + 1 => by
    have ih := extraZ_stable d V hd1 hd2 k
    unfold extraZ at ih ⊢
    rw [← Nat.add_assoc, length_filter_range_succ, ih]
    have : extraB d (V + k) = false := by
      cases h : extraB d (V + k)
      · rfl
      · exfalso
        unfold extraB at h
        rw [Bool.and_eq_true, decide_eq_true_iff, decide_eq_true_iff] at h
        have := two_pow_lt_of_dvd d (V + k) V hd1 hd2 (Int.dvd_iff_emod_eq_zero.2 h.2)
        omega
    rw [this]
    rfl

/-! ## The owner's labels `M = (A + 1)/2` -/

/-- The branch of an odd `A` is read off `A mod 8`. -/
theorem v2Z_mod_eight (A : Int) (hA : A % 2 = 1) :
    (v2Z (3 * A + 1) = 1 ↔ A % 4 = 3) ∧ (v2Z (3 * A + 1) = 2 ↔ A % 8 = 1) ∧
      (3 ≤ v2Z (3 * A + 1) ↔ A % 8 = 5) := by
  obtain ⟨h1, h2⟩ := syrZ_decomp A hA
  have hv := v2Z_pos A hA
  generalize hS : syrZ A = S at h1 h2
  generalize hw : v2Z (3 * A + 1) = w at h1 hv
  have hcase : (w = 1 ∧ A % 4 = 3) ∨ (w = 2 ∧ A % 8 = 1) ∨ (3 ≤ w ∧ A % 8 = 5) := by
    rcases Nat.lt_or_ge w 3 with hw3 | hw3
    · rcases (show w = 1 ∨ w = 2 by omega) with rfl | rfl
      · rw [show (2 : Int) ^ 1 = 2 from rfl] at h1
        exact Or.inl ⟨rfl, by omega⟩
      · rw [show (2 : Int) ^ 2 = 4 from rfl] at h1
        exact Or.inr (Or.inl ⟨rfl, by omega⟩)
    · obtain ⟨k, rfl⟩ : ∃ k, w = k + 3 := ⟨w - 3, by omega⟩
      rw [Int.pow_add, show (2 : Int) ^ 3 = 8 from rfl, Int.mul_comm _ 8, Int.mul_assoc] at h1
      exact Or.inr (Or.inr ⟨by omega, by omega⟩)
  refine ⟨⟨fun h => ?_, fun h => ?_⟩, ⟨fun h => ?_, fun h => ?_⟩, ⟨fun h => ?_, fun h => ?_⟩⟩
  all_goals rcases hcase with ⟨c1, c2⟩ | ⟨c1, c2⟩ | ⟨c1, c2⟩
  all_goals omega

/-- **The owner's two copies, in labels.** For odd `A` with `K(A) = d` and label
`M = (A + 1)/2`: `M` is even iff `A = -4d - 1` (`M = -2d`), `M ≡ 1 (mod 4)` iff `A = 8d + 1`
(`M = 4d + 1`), and every other preimage has `M ≡ 3 (mod 4)`. -/
theorem labels_two_copies (d A : Int) (hA : A % 2 = 1) (hd : dropZ A = d) :
    ((A + 1) / 2 % 2 = 0 ↔ A = -4 * d - 1) ∧ ((A + 1) / 2 % 4 = 1 ↔ A = 8 * d + 1) ∧
      (A ≠ -4 * d - 1 → A ≠ 8 * d + 1 → (A + 1) / 2 % 4 = 3) := by
  obtain ⟨m1, m2, m3⟩ := v2Z_mod_eight A hA
  obtain ⟨f1, _, f3⟩ := dropZ_forward A hA
  rw [hd] at f3
  have hpow2 : (2 : Int) ^ (v2Z (3 * A + 1) + 1) = 2 * (2 : Int) ^ v2Z (3 * A + 1) := by
    rw [Int.pow_succ, Int.mul_comm]
  rw [hpow2] at f3
  -- the branch determines `A` from `d`
  have b1 : v2Z (3 * A + 1) = 1 → A = -4 * d - 1 := by
    intro h
    rw [h, show (2 : Int) ^ 1 = 2 from rfl] at f3
    omega
  have b2 : v2Z (3 * A + 1) = 2 → A = 8 * d + 1 := by
    intro h
    rw [h, show (2 : Int) ^ 2 = 4 from rfl] at f3
    omega
  -- the unit preimages have branches 1 and 2
  have u1 : A = -4 * d - 1 → v2Z (3 * A + 1) = 1 := by
    rintro rfl
    exact (dropZ_unit_neg d).2.1
  have u2 : A = 8 * d + 1 → v2Z (3 * A + 1) = 2 := by
    rintro rfl
    exact (dropZ_unit_pos d).2.1
  refine ⟨⟨fun h => b1 (m1.2 (by omega)), fun h => by have := m1.1 (u1 h); omega⟩,
    ⟨fun h => b2 (m2.2 (by omega)), fun h => by have := m2.1 (u2 h); omega⟩, fun h1 h2 => ?_⟩
  have hv3 : 3 ≤ v2Z (3 * A + 1) := by
    rcases Nat.lt_or_ge (v2Z (3 * A + 1)) 3 with h | h
    · rcases (show v2Z (3 * A + 1) = 1 ∨ v2Z (3 * A + 1) = 2 by omega) with h' | h'
      · exact absurd (b1 h') h1
      · exact absurd (b2 h') h2
    · exact h
  have := m3.1 hv3
  omega

/-! ## Worked examples -/

/-- `d = 2` (`6d + 1 = 13 = 2^4 - 3`): the preimages are `-9` (`v = 1`), `17` (`v = 2`) and
`5` (`v = 4`), so `2 + N(13) = 3`. `d = 24` (`6d + 1 = 145 = 5 · 29`): `-97`, `193`, `77`
(`v = 3`) and `53` (`v = 5`), so `2 + N(145) = 4`. -/
theorem fibreZ_examples :
    fibreZ 2 6 = [-9, 17, 5] ∧ fibreZ 24 9 = [-97, 193, 77, 53] ∧
      (dropZ (-9) = 2 ∧ dropZ 17 = 2 ∧ dropZ 5 = 2) ∧
      (dropZ (-97) = 24 ∧ dropZ 193 = 24 ∧ dropZ 77 = 24 ∧ dropZ 53 = 24) := by
  decide

end ProcgenSelfieEdim
