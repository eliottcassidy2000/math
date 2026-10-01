import ProcgenSelfieEdim.Redei
import ProcgenSelfieEdim.CollatzDrop

set_option autoImplicit false

/-!
# THM-4524 A3(e): switching classes with constant `H` need `N = 2^k`

If `H` is constant on the switching class of a tournament `T` on `N ≥ 1` vertices, the
loop-set sum `Σ_L H(switch_L T) = 2 · N!` (A2, `switch_sum`) gives `2^(N-1) · H(T) = N!`.
Since `H(T)` is odd (Rédei, `redei`), `v_2(N!) = N - 1`; Legendre's formula
`v_2(N!) + s_2(N) = N` (`legendre_two`, `s_2` = binary digit sum) then forces `s_2(N) = 1`,
i.e. `N` is a power of two (`constant_switching_class`).
-/

namespace ProcgenSelfieEdim

/-- The binary digit sum `s_2(n)`. -/
def s2 (n : Nat) : Nat := if n = 0 then 0 else n % 2 + s2 (n / 2)
termination_by n
decreasing_by omega

theorem s2_eq (n : Nat) : s2 n = if n = 0 then 0 else n % 2 + s2 (n / 2) := by
  rw [s2]

theorem s2_zero : s2 0 = 0 := by
  rw [s2_eq, if_pos rfl]

theorem s2_double (q : Nat) : s2 (2 * q) = s2 q := by
  by_cases hq : q = 0
  · subst hq
    rfl
  · rw [s2_eq (2 * q), if_neg (by omega), show 2 * q / 2 = q by omega, show 2 * q % 2 = 0 by omega,
      Nat.zero_add]

theorem s2_double_add_one (q : Nat) : s2 (2 * q + 1) = s2 q + 1 := by
  rw [s2_eq (2 * q + 1), if_neg (by omega), show (2 * q + 1) / 2 = q by omega,
    show (2 * q + 1) % 2 = 1 by omega, Nat.add_comm]

theorem s2_eq_zero : ∀ (k n : Nat), n < k → s2 n = 0 → n = 0
  | 0, _, h, _ => absurd h (Nat.not_lt_zero _)
  | k + 1, n, hk, h => by
    by_cases hn : n = 0
    · exact hn
    · exfalso
      rw [s2_eq, if_neg hn] at h
      have h2 : s2 (n / 2) = 0 := by omega
      have := s2_eq_zero k (n / 2) (by omega) h2
      omega

/-- `s_2(n) = 1` only for powers of two. -/
theorem pow_two_of_s2_eq_one : ∀ (k n : Nat), n < k → s2 n = 1 → ∃ e, n = 2 ^ e
  | 0, _, h, _ => absurd h (Nat.not_lt_zero _)
  | k + 1, n, hk, h => by
    by_cases hn : n = 0
    · subst hn
      rw [s2_zero] at h
      exact absurd h (by decide)
    · rw [s2_eq, if_neg hn] at h
      by_cases hodd : n % 2 = 1
      · have h0 : s2 (n / 2) = 0 := by omega
        have := s2_eq_zero n (n / 2) (by omega) h0
        exact ⟨0, by rw [Nat.pow_zero]; omega⟩
      · have h1 : s2 (n / 2) = 1 := by omega
        obtain ⟨e, he⟩ := pow_two_of_s2_eq_one k (n / 2) (by omega) h1
        exact ⟨e + 1, by rw [Nat.pow_succ, ← he]; omega⟩

/-- `v_2` is additive on positive numbers. -/
theorem v2_mul (a b : Nat) (ha : 0 < a) (hb : 0 < b) : v2 (a * b) = v2 a + v2 b := by
  obtain ⟨h1, h2⟩ := decomp a ha
  obtain ⟨h3, h4⟩ := decomp b hb
  have e : a * b = 2 ^ (v2 a + v2 b) * (oddPart a * oddPart b) := by
    rw [Nat.pow_add]
    conv => lhs; rw [h1, h3]
    rw [Nat.mul_assoc, Nat.mul_left_comm (oddPart a), ← Nat.mul_assoc]
  have hodd : oddPart a * oddPart b % 2 = 1 := by
    rw [Nat.mul_mod, h2, h4]
  exact (v2_oddPart_of_eq _ _ _ e hodd).1

theorem v2_odd (n : Nat) (h : n % 2 = 1) : v2 n = 0 :=
  (v2_oddPart_of_eq n 0 n (by rw [Nat.pow_zero, Nat.one_mul]) h).1

theorem v2_double (q : Nat) (hq : 0 < q) : v2 (2 * q) = v2 q + 1 := by
  rw [v2_mul 2 q (by decide) hq, Nat.add_comm]
  have : v2 2 = 1 := (v2_oddPart_of_eq 2 1 1 rfl rfl).1
  rw [this]

/-- Adding one turns the trailing ones of `n` into zeros: `s_2(n+1) + v_2(n+1) = s_2(n) + 1`. -/
theorem s2_succ : ∀ (k n : Nat), n < k → s2 (n + 1) + v2 (n + 1) = s2 n + 1
  | 0, _, h => absurd h (Nat.not_lt_zero _)
  | k + 1, n, hk => by
    by_cases hev : n % 2 = 0
    · -- `n = 2q`: `n + 1` is odd
      obtain ⟨q, rfl⟩ : ∃ q, n = 2 * q := ⟨n / 2, by omega⟩
      rw [s2_double_add_one, s2_double, v2_odd (2 * q + 1) (by omega)]
    · -- `n = 2q + 1`: `n + 1 = 2(q + 1)`
      obtain ⟨q, rfl⟩ : ∃ q, n = 2 * q + 1 := ⟨n / 2, by omega⟩
      have ih := s2_succ k q (by omega)
      rw [show 2 * q + 1 + 1 = 2 * (q + 1) by omega, s2_double, v2_double (q + 1) (by omega),
        s2_double_add_one]
      omega

theorem fact_pos : ∀ (n : Nat), 0 < fact n
  | 0 => Nat.one_pos
  | n + 1 => Nat.mul_pos (Nat.succ_pos n) (fact_pos n)

/-- **Legendre's formula for `p = 2`.** `v_2(n!) + s_2(n) = n`. -/
theorem legendre_two : ∀ (n : Nat), v2 (fact n) + s2 n = n
  | 0 => by
    rw [s2_zero, show fact 0 = 1 from rfl, v2_odd 1 rfl]
  | n + 1 => by
    have ih := legendre_two n
    have hs := s2_succ (n + 1) n (Nat.lt_succ_self n)
    show v2 ((n + 1) * fact n) + s2 (n + 1) = n + 1
    rw [v2_mul (n + 1) (fact n) (Nat.succ_pos n) (fact_pos n)]
    omega

/-- **THM-4524 A3(e).** If `H` is constant on the switching class of a tournament on
`N ≥ 1` vertices, then `2^(N-1) · H(T) = N!` and `N` is a power of two. -/
theorem constant_switching_class (N : Nat) (hN : 1 ≤ N) (T : Nat → Nat → Bool)
    (hT : IsTournament N T) (hconst : ∀ L : Nat → Bool, hpCount N (switch L T) = hpCount N T) :
    2 ^ (N - 1) * hpCount N T = fact N ∧ ∃ k, N = 2 ^ k := by
  have hsum := switch_sum N hN T hT
  rw [List.map_congr_left (g := fun _ => hpCount N T) (fun bs _ => hconst (nthB bs)), sum_map_const,
    length_allBools] at hsum
  have hpow : 2 ^ N = 2 * 2 ^ (N - 1) := by
    rw [← Nat.pow_succ']
    congr 1
    omega
  have hfact : 2 ^ (N - 1) * hpCount N T = fact N := by
    rw [hpow] at hsum
    have h2 : 2 * (2 ^ (N - 1) * hpCount N T) = 2 * fact N := by
      rw [← hsum, Nat.mul_comm (hpCount N T), Nat.mul_assoc]
    exact Nat.eq_of_mul_eq_mul_left (by decide) h2
  refine ⟨hfact, ?_⟩
  have hodd := redei N T hT
  have hv : v2 (fact N) = N - 1 := (v2_oddPart_of_eq _ _ _ hfact.symm hodd).1
  have hl := legendre_two N
  exact pow_two_of_s2_eq_one (N + 1) N (Nat.lt_succ_self N) (by omega)

/-- The bound is attained at `N = 4`: every switch of `tourC 16` has exactly `3 = 4!/2^3`
Hamiltonian paths (kernel check over the 16 loop sets). -/
theorem constant_class_check :
    allR (fun bs => Nat.beq (lenR (hamPathsR 4 (switch (nthB bs) (tourC 16)))) 3) (allBools 4) = true := by
  decide +kernel

theorem constant_class_example (L : Nat → Bool) : hpCount 4 (switch L (tourC 16)) = 3 := by
  have hmem : ascList L 4 ∈ allBools 4 := (mem_allBools 4 _).2 (ascList_length L 4)
  have h := allR_eq_true _ _ constant_class_check _ hmem
  rw [Nat.beq_eq, lenR_eq] at h
  have hsame : SameOn 4 (switch L (tourC 16)) (switch (nthB (ascList L 4)) (tourC 16)) :=
    switch_loops_congr 4 L _ (tourC 16) (fun x hx => (nthB_ascList L 4 x hx).symm)
  rw [hpCount_congr 4 _ _ hsame]
  exact h

end ProcgenSelfieEdim
