import ProcgenSelfieEdim.Hypercube
import ProcgenSelfieEdim.ArcParity

set_option autoImplicit false

/-!
# THM-4525, Lemma L3: the alternating sum of an edge histogram

For the edge `e = {u, u + 2^i}` of `Q_d` (bit `i` of `u` is `0`) and `U` = all coordinates
except `i`, write `χ(w) = (-1)^{Σ_{j ≠ i} w_j}`. Then

  `Σ_r (-1)^r H_e(r) = Σ_{s ∈ S} (-1)^{d(e,s)} = χ(u) · Σ_{s ∈ S} χ(s)`

(`alt_sum`, `alt_sum_hist`). The key step is the projection lemma
`d(e, s) = #{j ≠ i : u_j ≠ s_j}` (`edgeDist_projection`). Within one direction `i` the
alternating sum therefore takes only the two values `±Σ_s χ(s)`.
-/

namespace ProcgenSelfieEdim

/-- `sgn n = (-1)^n`. -/
def sgn (n : Nat) : Int := if n % 2 = 0 then 1 else -1

theorem sgn_add (a b : Nat) : sgn (a + b) = sgn a * sgn b := by
  unfold sgn
  by_cases ha : a % 2 = 0 <;> by_cases hb : b % 2 = 0
  · rw [if_pos (by omega), if_pos ha, if_pos hb]
    rfl
  · rw [if_neg (by omega), if_pos ha, if_neg hb]
    rfl
  · rw [if_neg (by omega), if_neg ha, if_pos hb]
    rfl
  · rw [if_pos (by omega), if_neg ha, if_neg hb]
    rfl

theorem sgn_congr (a b : Nat) (h : a % 2 = b % 2) : sgn a = sgn b := by
  unfold sgn
  by_cases ha : a % 2 = 0
  · rw [if_pos ha, if_pos (by omega)]
  · rw [if_neg ha, if_neg (by omega)]

/-- Hamming distance as a sum over coordinates. -/
theorem hamming_eq_rsum : ∀ (d u v : Nat),
    hamming d u v = rsum d (fun j => if bit u j = bit v j then 0 else 1)
  | 0, _, _ => rfl
  | d + 1, u, v => by
    rw [hamming_succ, rsum_shift, hamming_eq_rsum d (u / 2) (v / 2), bit_zero, bit_zero]
    congr 1
    apply rsum_congr
    intro j _
    rw [bit_succ, bit_succ]

theorem bit_lt_two (u j : Nat) : bit u j < 2 := by
  unfold bit
  omega

/-- Flipping bit `i` (from `0` to `1`) changes no other bit. -/
theorem bit_flip : ∀ (i u j : Nat), bit u i = 0 → bit (u + 2 ^ i) j = if j = i then 1 else bit u j
  | 0, u, j, hb => by
    rw [bit_zero] at hb
    cases j with
    | zero =>
      rw [if_pos rfl, bit_zero, Nat.pow_zero]
      omega
    | succ j =>
      rw [if_neg (by omega), bit_succ, bit_succ, Nat.pow_zero, show (u + 1) / 2 = u / 2 by omega]
  | i + 1, u, j, hb => by
    rw [bit_succ] at hb
    have hp := two_pow_succ i
    cases j with
    | zero =>
      rw [if_neg (by omega), bit_zero, bit_zero]
      omega
    | succ j =>
      rw [bit_succ, bit_succ, show (u + 2 ^ (i + 1)) / 2 = u / 2 + 2 ^ i by omega,
        bit_flip i (u / 2) j hb]
      by_cases hj : j = i
      · rw [if_pos hj, if_pos (by omega)]
      · rw [if_neg hj, if_neg (by omega)]

/-- **Projection lemma.** For the edge `{u, u + 2^i}`, `d(e, s) = #{j ≠ i : u_j ≠ s_j}`. -/
theorem edgeDist_projection (d i u s : Nat) (hi : i < d) (hb : bit u i = 0) :
    edgeDist d u (u + 2 ^ i) s =
      rsum d (fun j => if j = i then 0 else (if bit u j = bit s j then 0 else 1)) := by
  unfold edgeDist
  rw [hamming_eq_rsum, hamming_eq_rsum]
  -- the two sums differ only at `j = i`
  have split : ∀ (f : Nat → Nat), rsum d f =
      rsum d (fun j => if j = i then 0 else f j) + f i := by
    intro f
    have := rsum_single d i f hi
    rw [← this, ← rsum_add]
    apply rsum_congr
    intro j _
    by_cases hj : j = i
    · subst hj
      rw [if_pos rfl, if_pos rfl, Nat.zero_add]
    · rw [if_neg hj, if_neg (Ne.symm hj), Nat.add_zero]
  rw [split (fun j => if bit u j = bit s j then 0 else 1),
    split (fun j => if bit (u + 2 ^ i) j = bit s j then 0 else 1)]
  have hrest : rsum d (fun j => if j = i then 0 else
      (if bit (u + 2 ^ i) j = bit s j then 0 else 1)) =
      rsum d (fun j => if j = i then 0 else (if bit u j = bit s j then 0 else 1)) := by
    apply rsum_congr
    intro j _
    by_cases hj : j = i
    · rw [if_pos hj, if_pos hj]
    · rw [if_neg hj, if_neg hj, bit_flip i u j hb, if_neg hj]
  rw [hrest, bit_flip i u i hb, if_pos rfl, hb]
  have := bit_lt_two s i
  by_cases hs : bit s i = 0
  · rw [hs, if_pos rfl, if_neg (by decide)]
    omega
  · rw [if_neg (Ne.symm hs), if_pos (by omega)]
    omega

/-- `χ(w) = (-1)^{Σ_{j ≠ i} w_j}` (bits `0, …, d-1`). -/
def chi (d i w : Nat) : Int := sgn (rsum d (fun j => if j = i then 0 else bit w j))

/-- **L3.** `Σ_{s ∈ S} (-1)^{d(e,s)} = χ(u) · Σ_{s ∈ S} χ(s)` for `e = {u, u + 2^i}`. -/
theorem alt_sum (d i u : Nat) (hi : i < d) (hb : bit u i = 0) :
    ∀ (S : List Nat), (S.map (fun s => sgn (edgeDist d u (u + 2 ^ i) s))).sum =
      chi d i u * (S.map (chi d i)).sum
  | [] => by
    show (0 : Int) = chi d i u * 0
    rw [Int.mul_zero]
  | s :: t => by
    rw [List.map_cons, List.map_cons, List.sum_cons, List.sum_cons, alt_sum d i u hi hb t, Int.mul_add]
    congr 1
    rw [edgeDist_projection d i u s hi hb]
    unfold chi
    rw [← sgn_add]
    apply sgn_congr
    rw [← rsum_add]
    apply rsum_mod_two_congr
    intro j _
    by_cases hj : j = i
    · rw [if_pos hj, if_pos hj, if_pos hj]
    · rw [if_neg hj, if_neg hj, if_neg hj]
      have h1 := bit_lt_two u j
      have h2 := bit_lt_two s j
      by_cases he : bit u j = bit s j
      · rw [if_pos he]
        omega
      · rw [if_neg he]
        omega

/-- **L3, histogram form.** `Σ_{r even} H_e(r) - Σ_{r odd} H_e(r) = χ(u) · Σ_{s ∈ S} χ(s)`. -/
theorem alt_sum_hist (d i u : Nat) (hi : i < d) (hb : bit u i = 0) (S : List Nat) :
    (rsum d (fun r => if r % 2 = 0 then hist d S u (u + 2 ^ i) r else 0) : Int) -
      (rsum d (fun r => if r % 2 = 0 then 0 else hist d S u (u + 2 ^ i) r) : Int) =
      chi d i u * (S.map (chi d i)).sum := by
  have he : hamming d u (u + 2 ^ i) = 1 := hamming_flip d i u hi hb
  rw [← alt_sum d i u hi hb S]
  have e1 := weighted_sum_eq d u (u + 2 ^ i) he (fun r => if r % 2 = 0 then 1 else 0) S
  have e2 := weighted_sum_eq d u (u + 2 ^ i) he (fun r => if r % 2 = 0 then 0 else 1) S
  have c1 : rsum d (fun r => (if r % 2 = 0 then 1 else 0) * hist d S u (u + 2 ^ i) r) =
      rsum d (fun r => if r % 2 = 0 then hist d S u (u + 2 ^ i) r else 0) :=
    rsum_congr d _ _ (fun r _ => by by_cases h : r % 2 = 0 <;> simp only [h, if_true, if_false] <;> omega)
  have c2 : rsum d (fun r => (if r % 2 = 0 then 0 else 1) * hist d S u (u + 2 ^ i) r) =
      rsum d (fun r => if r % 2 = 0 then 0 else hist d S u (u + 2 ^ i) r) :=
    rsum_congr d _ _ (fun r _ => by by_cases h : r % 2 = 0 <;> simp only [h, if_true, if_false] <;> omega)
  rw [← c1, ← c2, ← e1, ← e2]
  -- per landmark: `sgn k = [k even] - [k odd]`
  clear e1 e2 c1 c2
  induction S with
  | nil => rfl
  | cons s t ih =>
    rw [List.map_cons, List.map_cons, List.map_cons, List.sum_cons, List.sum_cons, List.sum_cons, ← ih]
    unfold sgn
    by_cases h : edgeDist d u (u + 2 ^ i) s % 2 = 0
    · rw [if_pos h, if_pos h, if_pos h]
      push_cast
      omega
    · rw [if_neg h, if_neg h, if_neg h]
      push_cast
      omega

theorem sgn_cases (n : Nat) : sgn n = 1 ∨ sgn n = -1 := by
  unfold sgn
  by_cases h : n % 2 = 0
  · rw [if_pos h]
    exact Or.inl rfl
  · rw [if_neg h]
    exact Or.inr rfl

/-- **L3, corollary.** Within one direction `i`, the alternating sum of an edge histogram takes
only the two values `± Σ_{s ∈ S} χ(s)`. -/
theorem alt_sum_two_values (d i u : Nat) (hi : i < d) (hb : bit u i = 0) (S : List Nat) :
    (S.map (fun s => sgn (edgeDist d u (u + 2 ^ i) s))).sum = (S.map (chi d i)).sum ∨
      (S.map (fun s => sgn (edgeDist d u (u + 2 ^ i) s))).sum = -(S.map (chi d i)).sum := by
  rw [alt_sum d i u hi hb S]
  unfold chi
  rcases sgn_cases (rsum d (fun j => if j = i then 0 else bit u j)) with h | h
  · rw [h, Int.one_mul]
    exact Or.inl rfl
  · rw [h, Int.neg_one_mul]
    exact Or.inr rfl

end ProcgenSelfieEdim
