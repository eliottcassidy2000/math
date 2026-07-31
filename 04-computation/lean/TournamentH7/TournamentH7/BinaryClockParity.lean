/-
THM-2976 (repo canon): binary-clock parity for critical-run deficit ledgers
(AMM 12592 / THM-2966 spine normal form).

Lean core, over F_2 = ZMod 2, with `beta M d` the forced-parity polynomial
  beta = (1+X)^(M+d+1) + (1+X^(M+1)) * (1+X)^d :

* `checkpoint_vanishing` (T1): if M+1 = 2^r then beta = 0 — for EVERY d.
* `clock_coeff_lt`  (T2, part i): if M+1 = 2^v * c with c odd, 3 ≤ c, then
  beta.coeff o = 0 for all o < 2^v.
* `clock_coeff_top` (T2, part ii): beta.coeff (2^v) = 1.

So the minimal forced-odd position is exactly 2^(v_2(M+1)) whenever M+1 is
not a power of two, and there is none at all when it is: the parity ledger
runs on a binary clock. Mechanism: Frobenius (1+X)^(2^v) = 1+X^(2^v) and
the factorization beta = expand_(2^v)[(1+Y)^c + 1 + Y^c] * (1+X)^d, whose
bracket has zero constant term and coefficient c ≡ 1 (mod 2) at Y^1.
-/
import Mathlib

open Polynomial

namespace BinaryClockParity

abbrev F2 := ZMod 2

/-- The forced-parity polynomial of THM-2976 Lemma 0 (lane D closed form). -/
noncomputable def beta (M d : ℕ) : Polynomial F2 :=
  (1 + X) ^ (M + d + 1) + (1 + X ^ (M + 1)) * (1 + X) ^ d

lemma coeff_oneX (n k : ℕ) :
    ((1 + X : Polynomial F2) ^ n).coeff k = (n.choose k : F2) := by
  rw [show (1 + X : Polynomial F2) = X + 1 from add_comm 1 X]
  exact coeff_X_add_one_pow F2 n k

/-- Frobenius over F_2[X]: `(1+X)^(2^r) = 1 + X^(2^r)`. -/
theorem one_add_X_pow_two_pow (r : ℕ) :
    (1 + X : Polynomial F2) ^ 2 ^ r = 1 + X ^ 2 ^ r := by
  induction r with
  | zero => simp
  | succ n ih =>
    rw [pow_succ 2 n, pow_mul, ih, CharTwo.add_sq, one_pow, ← pow_mul,
      ← pow_succ 2 n]

/-- **T1 (checkpoint vanishing).** If `M + 1 = 2^r` then `beta M d = 0`,
for every depth value `d`. -/
theorem checkpoint_vanishing (M d r : ℕ) (hM : M + 1 = 2 ^ r) :
    beta M d = 0 := by
  unfold beta
  have hA : M + d + 1 = 2 ^ r + d := by omega
  rw [hA, pow_add, one_add_X_pow_two_pow, ← hM]
  exact CharTwo.add_self_eq_zero _

/-- The clock bracket `(1+Y)^c + 1 + Y^c`. -/
noncomputable def bracket (c : ℕ) : Polynomial F2 :=
  (1 + X) ^ c + 1 + X ^ c

theorem bracket_coeff_zero (c : ℕ) (hc : 1 ≤ c) :
    (bracket c).coeff 0 = 0 := by
  unfold bracket
  rw [coeff_add, coeff_add, coeff_oneX, coeff_X_pow,
    if_neg (by omega : (0 : ℕ) ≠ c)]
  simp only [Nat.choose_zero_right, Nat.cast_one, coeff_one, add_zero]
  decide

theorem bracket_coeff_one (c : ℕ) (hodd : Odd c) (hc : 3 ≤ c) :
    (bracket c).coeff 1 = 1 := by
  unfold bracket
  rw [coeff_add, coeff_add, coeff_oneX, coeff_X_pow,
    if_neg (by omega : (1 : ℕ) ≠ c), Nat.choose_one_right]
  have hcast : (c : F2) = 1 := by
    rcases hodd with ⟨k, hk⟩
    subst hk
    push_cast
    have h2 : (2 : F2) = 0 := by decide
    rw [h2]
    ring
  have h1 : ((1 : Polynomial F2)).coeff 1 = 0 := by
    rw [coeff_one]
    norm_num
  rw [hcast, h1]
  ring

/-- Structural identity: for `M + 1 = 2^v * c`,
`beta M d = expand_(2^v)(bracket c) * (1+X)^d`. -/
theorem beta_eq_expand_bracket (M d v c : ℕ) (hM : M + 1 = 2 ^ v * c) :
    beta M d = (expand F2 (2 ^ v) (bracket c)) * (1 + X) ^ d := by
  unfold beta bracket
  have hA : M + d + 1 = 2 ^ v * c + d := by omega
  rw [hA, pow_add, pow_mul, one_add_X_pow_two_pow]
  rw [map_add, map_add, map_pow, map_pow, map_one, map_add, map_one,
    expand_X, ← pow_mul, ← hM]
  ring

lemma expand_bracket_coeff_small (v c : ℕ) (hc1 : 1 ≤ c) (i : ℕ)
    (hi : i < 2 ^ v) :
    (expand F2 (2 ^ v) (bracket c)).coeff i = 0 := by
  have hvpos : 0 < 2 ^ v := pow_pos (by norm_num : (0:ℕ) < 2) v
  rw [coeff_expand hvpos]
  split_ifs with hdvd
  · rcases hdvd with ⟨t, ht⟩
    rcases Nat.eq_zero_or_pos t with h0 | hpos
    · subst h0
      simp only [Nat.mul_zero] at ht
      subst ht
      simpa using bracket_coeff_zero c hc1
    · exfalso
      have h1 : 2 ^ v * 1 ≤ 2 ^ v * t := Nat.mul_le_mul_left _ hpos
      have h2 : 2 ^ v ≤ i := by rw [ht]; simpa using h1
      omega
  · rfl

/-- **T2(i) (silence below the clock).** If `M + 1 = 2^v * c` with `c` odd
and `3 ≤ c`, then `(beta M d).coeff o = 0` for every `o < 2^v`. -/
theorem clock_coeff_lt (M d v c : ℕ) (hM : M + 1 = 2 ^ v * c)
    (_hodd : Odd c) (hc : 3 ≤ c) (o : ℕ) (ho : o < 2 ^ v) :
    (beta M d).coeff o = 0 := by
  rw [beta_eq_expand_bracket M d v c hM, coeff_mul]
  apply Finset.sum_eq_zero
  rintro ⟨i, j⟩ hij
  rw [Finset.mem_antidiagonal] at hij
  have hi : i < 2 ^ v := by omega
  rw [expand_bracket_coeff_small v c (by omega) i hi, zero_mul]

/-- **T2(ii) (the clock strikes).** Under the same hypotheses,
`(beta M d).coeff (2^v) = 1`: the minimal forced-odd position is exactly
`2^(v_2(M+1))`. -/
theorem clock_coeff_top (M d v c : ℕ) (hM : M + 1 = 2 ^ v * c)
    (hodd : Odd c) (hc : 3 ≤ c) :
    (beta M d).coeff (2 ^ v) = 1 := by
  have hvpos : 0 < 2 ^ v := pow_pos (by norm_num : (0:ℕ) < 2) v
  rw [beta_eq_expand_bracket M d v c hM, coeff_mul]
  rw [Finset.sum_eq_single ((2 : ℕ) ^ v, (0 : ℕ))]
  · rw [coeff_expand hvpos, if_pos dvd_rfl, Nat.div_self hvpos,
      bracket_coeff_one c hodd hc, coeff_oneX, Nat.choose_zero_right]
    simp
  · rintro ⟨i, j⟩ hij hne
    rw [Finset.mem_antidiagonal] at hij
    rcases Nat.lt_or_ge i (2 ^ v) with hi | hi
    · rw [expand_bracket_coeff_small v c (by omega) i hi, zero_mul]
    · have hieq : i = 2 ^ v := by omega
      have hj : j = 0 := by omega
      exact absurd (Prod.ext hieq hj) hne
  · intro hmem
    exact absurd (Finset.mem_antidiagonal.mpr (by omega)) hmem

end BinaryClockParity
