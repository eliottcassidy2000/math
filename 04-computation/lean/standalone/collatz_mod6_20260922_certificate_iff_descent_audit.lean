/- Audit scratch (NOT part of any package): the free-parameter descent certificate of
   collatz_mod6_20260922_descent_certificate.lean is EQUIVALENT to plain strict descent
   `∃ t > 0, iterate collatz t n < n` for every positive start. Core Lean 4.30.0, no Mathlib. -/
def iterate {α : Type} (f : α → α) : Nat → α → α
  | 0, x => x
  | k + 1, x => iterate f k (f x)

def collatz (n : Nat) : Nat := if n % 2 = 0 then n / 2 else 3 * n + 1

def DescentCertificate (n : Nat) : Prop :=
  ∃ t K L B : Nat, 0 < t ∧
    2 ^ K * iterate collatz t n = 3 ^ L * n + B ∧
    3 ^ L * n + B < 2 ^ K * n

theorem collatz_pos (n : Nat) (hn : 0 < n) : 0 < collatz n := by
  unfold collatz
  split <;> omega

theorem iterate_pos (t : Nat) : ∀ n, 0 < n → 0 < iterate collatz t n := by
  induction t with
  | zero => intro n hn; simpa [iterate] using hn
  | succ k ih => intro n hn; simp only [iterate]; exact ih _ (collatz_pos n hn)

theorem lt_two_pow (n : Nat) : n < 2 ^ n := by
  induction n with
  | zero => decide
  | succ k ih => rw [Nat.pow_succ]; omega

theorem certificate_iff_descent (n : Nat) (hn : 0 < n) :
    DescentCertificate n ↔ ∃ t, 0 < t ∧ iterate collatz t n < n := by
  constructor
  · rintro ⟨t, K, L, B, ht, hid, hlt⟩
    refine ⟨t, ht, ?_⟩
    have h : 2 ^ K * iterate collatz t n < 2 ^ K * n := by rw [hid]; exact hlt
    exact (Nat.mul_lt_mul_left (Nat.two_pow_pos K)).1 h
  · rintro ⟨t, ht, hlt⟩
    have hy : 0 < iterate collatz t n := iterate_pos t n hn
    have h3 : n ≤ 2 ^ n * iterate collatz t n := by
      have h1 := lt_two_pow n
      have h2 := Nat.mul_le_mul_left (2 ^ n) hy
      rw [Nat.mul_one] at h2
      omega
    have h4 : 2 ^ n * iterate collatz t n < 2 ^ n * n :=
      (Nat.mul_lt_mul_left (Nat.two_pow_pos n)).2 hlt
    refine ⟨t, n, 0, 2 ^ n * iterate collatz t n - n, ht, ?_, ?_⟩
    · simp only [Nat.pow_zero, Nat.one_mul]; omega
    · simp only [Nat.pow_zero, Nat.one_mul]; omega
