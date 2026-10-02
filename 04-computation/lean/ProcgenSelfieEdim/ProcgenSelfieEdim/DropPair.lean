import ProcgenSelfieEdim.CollatzDrop
import ProcgenSelfieEdim.PolyNorm

set_option autoImplicit false
set_option linter.unusedSimpArgs false

/-!
# THM-4530 (b): two consecutive drops determine the point (Theorem 4.5)

On the odd `A ≥ 1`, the map `A ↦ (K(A), K(S(A)))` is injective, for the Syracuse map
`S(A) = oddpart(3A + 1)` (`drop_pair_injective`) and for the `3x - 1` map
`T(A) = oddpart(3A - 1)` (`dropM_pair_injective`).

The proof is the note's elimination. Write `M_v = 2^v - 3`, `S = S(A)`, `S' = S(A')` and
`S(S) = S₂`, `S(S') = S₂'`, with valuations `v, v'` (first step) and `u, u'` (second step).
Equal first drops give `M_v S = M_v' S'`, equal second drops give `S - S₂ = S' - S₂'`.

* If `v = v'` the points coincide. Otherwise both `v, v' ≥ 2` (a sign argument: `M_1 = -1`),
  likewise `u ≠ u'` and `u, u' ≥ 2`; by symmetry `u < u'`, `δ = u' - u`, `P = 2^δ`.
* `elim_identities`: with `E = P M_u M_v' - M_u' M_v`, `S E = ε (P - 1) M_v'` and
  `S' E = ε (P - 1) M_v` (`ε = ±1` the sheet). These are polynomial identities, checked by
  the `pnorm` normaliser and `omega`.
* `3x + 1` (`ε = 1`, so `E > 0`): `v' > v` (`plus_gt`); then `M_u (2^e - 1) < 4` forces
  `u = 2` and `e = v' - v ∈ {1, 2}`; `e = 2` is too big, and `e = 1` leaves `δ = 1`, `v = 3`,
  `S = 13`, where `3A + 1 = 8 · 13` is impossible (`plus_lt`; the note's `n₁ = 65 ≢ 1 mod 6`).
* `3x - 1` (`ε = -1`, so `E < 0`): `v' < v` is impossible (`minus_gt`); `v' > v` forces
  `u = 2`, `e = 1`, and `2^v ∈ (9(P - 1)/(2P - 3), 10(P - 1)/(2P - 3)]`, which contains no
  power of two (`minus_lt`).

The divisibility step `E ∣ (P - 1)` uses that `E` is prime to 3 (`three_dvd_of_mul`).
-/

namespace ProcgenSelfieEdim

/-! ## Small arithmetic facts -/

theorem two_pow_mod_three (n : Nat) : (2 : Int) ^ n % 3 = 1 ∨ (2 : Int) ^ n % 3 = 2 := by
  induction n with
  | zero => exact Or.inl (by decide)
  | succ k ih =>
    rw [Int.pow_succ]
    omega

theorem two_pow_cases (v : Nat) (hv : 2 ≤ v) :
    (2 : Int) ^ v = 4 ∨ (2 : Int) ^ v = 8 ∨ 16 ≤ (2 : Int) ^ v := by
  rcases (show v = 2 ∨ v = 3 ∨ 4 ≤ v by omega) with rfl | rfl | h
  · exact Or.inl (by decide)
  · exact Or.inr (Or.inl (by decide))
  · have h1 := two_pow_mono 4 v h
    have h16 : (2 : Int) ^ 4 = 16 := by decide
    exact Or.inr (Or.inr (by omega))

theorem two_pow_cases' (e : Nat) (he : 1 ≤ e) :
    (2 : Int) ^ e = 2 ∨ (2 : Int) ^ e = 4 ∨ 8 ≤ (2 : Int) ^ e := by
  rcases (show e = 1 ∨ e = 2 ∨ 3 ≤ e by omega) with rfl | rfl | h
  · exact Or.inl (by decide)
  · exact Or.inr (Or.inl (by decide))
  · have h1 := two_pow_mono 3 e h
    have h8 : (2 : Int) ^ 3 = 8 := by decide
    exact Or.inr (Or.inr (by omega))

theorem two_pow_lt (v w : Nat) (h : v < w) : (2 : Int) ^ v < (2 : Int) ^ w := by
  have h1 := two_pow_mono (v + 1) w h
  have h2 : (2 : Int) ^ (v + 1) = (2 : Int) ^ v * 2 := Int.pow_succ _ _
  have h3 : (0 : Int) < (2 : Int) ^ v := Int.pow_pos (by decide)
  omega

theorem le_of_mul_eq_pos (a b c : Int) (ha : 1 ≤ a) (hb : 0 ≤ b) (h : a * b = c) : b ≤ c := by
  have h1 : b * 1 ≤ b * a := Int.mul_le_mul_of_nonneg_left ha hb
  rw [Int.mul_one, Int.mul_comm] at h1
  omega

theorem pos_of_mul_pos (a b : Int) (hb : 0 < b) (h : 0 < a * b) : 0 < a := by
  by_cases ha : 0 < a
  · exact ha
  · have : a * b ≤ 0 := Int.mul_nonpos_of_nonpos_of_nonneg (by omega) (Int.le_of_lt hb)
    omega

theorem neg_of_mul_pos_neg (a b : Int) (ha : 0 < a) (h : a * b < 0) : b < 0 := by
  by_cases hb : b < 0
  · exact hb
  · have : 0 ≤ a * b := Int.mul_nonneg (Int.le_of_lt ha) (by omega)
    omega

theorem pos_of_mul_pos_left (a b : Int) (ha : 0 < a) (h : 0 < a * b) : 0 < b := by
  by_cases hb : 0 < b
  · exact hb
  · have : a * b ≤ 0 := Int.mul_nonpos_of_nonneg_of_nonpos (Int.le_of_lt ha) (by omega)
    omega

/-- If `E` is prime to 3 and `c E = 3k`, then `3 ∣ c`. -/
theorem three_dvd_of_mul (c E k : Int) (hE : E % 3 = 1 ∨ E % 3 = 2) (h : c * E = 3 * k) :
    c % 3 = 0 := by
  have h1 := Int.mul_emod c E 3
  rw [h] at h1
  rcases (show c % 3 = 0 ∨ c % 3 = 1 ∨ c % 3 = 2 by omega) with hc | hc | hc
  · exact hc
  · rw [hc] at h1
    rcases hE with hE | hE <;> rw [hE] at h1 <;> omega
  · rw [hc] at h1
    rcases hE with hE | hE <;> rw [hE] at h1 <;> omega

/-- `X (2P - 3)` is prime to 3 when `X` and `P` are. -/
theorem mod_three_prod (X P : Int) (hX : X % 3 = 1 ∨ X % 3 = 2) (hP : P % 3 = 1 ∨ P % 3 = 2) :
    X * (2 * P - 3) % 3 = 1 ∨ X * (2 * P - 3) % 3 = 2 := by
  have h1 := Int.mul_emod X (2 * P - 3) 3
  rcases hP with hP | hP
  · have h2 : (2 * P - 3) % 3 = 2 := by omega
    rw [h2] at h1
    rcases hX with hX | hX <;> rw [hX] at h1 <;> omega
  · have h2 : (2 * P - 3) % 3 = 1 := by omega
    rw [h2] at h1
    rcases hX with hX | hX <;> rw [hX] at h1 <;> omega

/-- Equal values `M_a S = M_b S'` with `S, S' > 0` and `a ≠ b` force `a, b ≥ 2`
(`M_1 = -1` is the only negative modulus with exponent `≥ 1`). -/
theorem both_ge_two (a b : Nat) (S Sp : Int) (ha : 1 ≤ a) (hb : 1 ≤ b) (hS : 0 < S) (hSp : 0 < Sp)
    (hne : a ≠ b) (h : ((2 : Int) ^ a - 3) * S = ((2 : Int) ^ b - 3) * Sp) : 2 ≤ a ∧ 2 ≤ b := by
  have key : ∀ (x y : Nat) (Sx Sy : Int), 1 ≤ y → x = 1 → 0 < Sx → 0 < Sy →
      ((2 : Int) ^ x - 3) * Sx = ((2 : Int) ^ y - 3) * Sy → y = 1 := by
    intro x y Sx Sy hy hx hSx hSy hxy
    subst hx
    rw [show (2 : Int) ^ 1 - 3 = -1 by decide] at hxy
    by_cases hy2 : 2 ≤ y
    · have h4 := pow_two_ge_four y hy2
      have : 0 < ((2 : Int) ^ y - 3) * Sy := Int.mul_pos (by omega) hSy
      omega
    · omega
  constructor
  · by_cases ha1 : a = 1
    · exact absurd ((key a b S Sp hb ha1 hS hSp h).trans ha1.symm) (Ne.symm hne)
    · omega
  · by_cases hb1 : b = 1
    · exact absurd ((key b a Sp S ha hb1 hSp hS h.symm).trans hb1.symm) hne
    · omega

/-! ## The elimination identities -/

/-- With `n₂ = M_u S₂ = M_u' S₂'` eliminated: `S E = ε (P - 1)(Y - 3)` and
`S' E = ε (P - 1)(X - 3)`, where `X = 2^v`, `Y = 2^v'`, `U = 2^u`, `U P = 2^u'` and
`E = P (U - 3)(Y - 3) - (U P - 3)(X - 3)`. -/
theorem elim_identities (ε X Y U P S Sp S2 S2p : Int)
    (H1 : (X - 3) * S = (Y - 3) * Sp)
    (H3 : 3 * S + ε = U * S2) (H3p : 3 * Sp + ε = U * P * S2p)
    (H4 : S - S2 = Sp - S2p) :
    S * (P * (U - 3) * (Y - 3) - (U * P - 3) * (X - 3)) = ε * (P - 1) * (Y - 3) ∧
      Sp * (P * (U - 3) * (Y - 3) - (U * P - 3) * (X - 3)) = ε * (P - 1) * (X - 3) := by
  have n2 : (U - 3) * S2 = (U * P - 3) * S2p := by
    pnorm at H3 H3p ⊢
    omega
  have K1 : P * (U - 3) * (3 * S + ε) = (U * P - 3) * (3 * Sp + ε) := by
    have a := congrArg (fun z => P * (U - 3) * z) H3
    have b := congrArg (fun z => (U * P - 3) * z) H3p
    have c := congrArg (fun z => P * U * z) n2
    simp only at a b c
    pnorm at a b c ⊢
    omega
  constructor
  · have a := congrArg (fun z => (Y - 3) * z) K1
    have b := congrArg (fun z => 3 * (U * P - 3) * z) H1
    simp only at a b
    pnorm at a b ⊢
    omega
  · have a := congrArg (fun z => (X - 3) * z) K1
    have b := congrArg (fun z => 3 * P * (U - 3) * z) H1
    simp only at a b
    pnorm at a b ⊢
    omega

/-! ## The `3x + 1` sheet -/

/-- `3x + 1`, case `v' < v`: then `E < 0`, contradicting `S E = (P - 1)(Y - 3) > 0`. -/
theorem plus_gt (X Y U P S : Int) (hY : 4 ≤ Y) (hYX : Y < X) (hU : 4 ≤ U) (hP : 2 ≤ P)
    (hS : 0 < S)
    (K2 : S * (P * (U - 3) * (Y - 3) - (U * P - 3) * (X - 3)) = 1 * (P - 1) * (Y - 3)) : False := by
  generalize hE : P * (U - 3) * (Y - 3) - (U * P - 3) * (X - 3) = E at K2
  have hpos : 0 < 1 * (P - 1) * (Y - 3) := by
    rw [Int.one_mul]
    exact Int.mul_pos (by omega) (by omega)
  have hEpos : 0 < E := pos_of_mul_pos_left S E hS (by rw [K2]; exact hpos)
  have m1 : P * (U - 3) * (Y - X) < 0 :=
    Int.mul_neg_of_pos_of_neg (Int.mul_pos (by omega) (by omega)) (by omega)
  have m2 : 0 < (P - 1) * (X - 3) := Int.mul_pos (by omega) (by omega)
  pnorm at hE m1 m2
  omega

/-- `3x + 1`, case `v < v'` (`Y = X Q`): only `u = 2`, `Q = 2`, `P = 2`, `X = 8`, `S = 13`
survives the bounds, and then `3A + 1 = 8 · 13` has no solution. -/
theorem plus_lt (X Q U P S Sp : Int)
    (hX : X = 4 ∨ X = 8 ∨ 16 ≤ X) (hX3 : X % 3 = 1 ∨ X % 3 = 2)
    (hU : U = 4 ∨ 8 ≤ U) (hP : 2 ≤ P) (hP3 : P % 3 = 1 ∨ P % 3 = 2)
    (hQ : Q = 2 ∨ Q = 4 ∨ 8 ≤ Q) (hS : 0 < S) (hSp : 0 < Sp)
    (K2 : S * (P * (U - 3) * (X * Q - 3) - (U * P - 3) * (X - 3)) = 1 * (P - 1) * (X * Q - 3))
    (K3 : Sp * (P * (U - 3) * (X * Q - 3) - (U * P - 3) * (X - 3)) = 1 * (P - 1) * (X - 3))
    (hA : (X * S - 1) % 3 = 0) : False := by
  have hX4 : 4 ≤ X := by omega
  have hQ2 : 2 ≤ Q := by omega
  have hU4 : 4 ≤ U := by omega
  generalize hE : P * (U - 3) * (X * Q - 3) - (U * P - 3) * (X - 3) = E at K2 K3
  have hXQ : 8 ≤ X * Q := by
    have := Int.mul_le_mul hX4 hQ2 (by decide) (by omega)
    omega
  have hEpos : 0 < E := by
    apply pos_of_mul_pos_left S E hS
    rw [K2, Int.one_mul]
    exact Int.mul_pos (by omega) (by omega)
  have hEle : E ≤ (P - 1) * (X - 3) := by
    rw [← Int.one_mul ((P - 1) * (X - 3)), ← Int.mul_assoc]
    exact le_of_mul_eq_pos Sp E _ (by omega) (Int.le_of_lt hEpos) K3
  -- `(U - 3)(Q - 1) ≤ 3`
  have hW : (U - 3) * (Q - 1) ≤ 3 := by
    by_cases hW : 4 ≤ (U - 3) * (Q - 1)
    · exfalso
      have m := Int.mul_le_mul_of_nonneg_left hW (show 0 ≤ P * X from Int.mul_nonneg (by omega) (by omega))
      have hE' := hE
      pnorm at m hE' hEle
      omega
    · omega
  have hU3 : U - 3 ≤ 3 := by
    have m := Int.mul_le_mul_of_nonneg_left (show 1 ≤ Q - 1 by omega) (show 0 ≤ U - 3 by omega)
    rw [Int.mul_one] at m
    omega
  have hQ3 : Q - 1 ≤ 3 := by
    have m := Int.mul_le_mul_of_nonneg_right (show 1 ≤ U - 3 by omega) (show 0 ≤ Q - 1 by omega)
    rw [Int.one_mul] at m
    omega
  have hU' : U = 4 := by omega
  subst hU'
  rcases (show Q = 2 ∨ Q = 4 by omega) with rfl | rfl
  · -- `e = 1`
    have hc : (S - 2 * Sp) * E = 3 * (P - 1) := by
      pnorm at K2 K3 ⊢
      omega
    have hEX : E = 9 * (P - 1) - X * (2 * P - 3) := by
      rw [← hE]
      pnorm
      omega
    have hZ := mod_three_prod X P hX3 hP3
    have hE3 : E % 3 = 1 ∨ E % 3 = 2 := by omega
    have hcpos : 0 < S - 2 * Sp := pos_of_mul_pos _ E hEpos (by rw [hc]; omega)
    have hc3 := three_dvd_of_mul _ E _ hE3 hc
    obtain ⟨c', hc'⟩ : ∃ c', S - 2 * Sp = 3 * c' := ⟨(S - 2 * Sp) / 3, by omega⟩
    rw [hc', Int.mul_assoc] at hc
    have hc'E : c' * E = P - 1 := by omega
    have hEle' : E ≤ P - 1 := le_of_mul_eq_pos c' E _ (by omega) (Int.le_of_lt hEpos) hc'E
    rcases hX with rfl | rfl | hX16
    · omega
    · -- `X = 8`: `P = 2`, `E = 1`, `S = 13`
      have hP2 : P = 2 := by omega
      subst hP2
      have hE1 : E = 1 := by omega
      subst hE1
      rw [Int.mul_one] at K2
      omega
    · have m := Int.mul_le_mul_of_nonneg_right hX16 (show 0 ≤ 2 * P - 3 by omega)
      omega
  · -- `e = 2`: `E = 3X + 9(P - 1)` but `E ≤ 9(P - 1)`
    have hc : (S - 4 * Sp) * E = 9 * (P - 1) := by
      pnorm at K2 K3 ⊢
      omega
    have hcpos : 0 < S - 4 * Sp := pos_of_mul_pos _ E hEpos (by rw [hc]; omega)
    have hle := le_of_mul_eq_pos _ E _ (by omega) (Int.le_of_lt hEpos) hc
    pnorm at hE
    omega

/-- **Core of Theorem 4.5 (`3x + 1`).** -/
theorem pair_core_plus (v v' u u' : Nat) (S Sp S2 S2p : Int)
    (hv : 1 ≤ v) (hv' : 1 ≤ v') (hu : 1 ≤ u) (hu' : 1 ≤ u')
    (hS : 0 < S) (hSp : 0 < Sp) (hS2 : 0 < S2) (hS2p : 0 < S2p)
    (H1 : ((2 : Int) ^ v - 3) * S = ((2 : Int) ^ v' - 3) * Sp)
    (H3 : 3 * S + 1 = (2 : Int) ^ u * S2) (H3p : 3 * Sp + 1 = (2 : Int) ^ u' * S2p)
    (H4 : S - S2 = Sp - S2p)
    (hA : ((2 : Int) ^ v * S - 1) % 3 = 0) (hA' : ((2 : Int) ^ v' * Sp - 1) % 3 = 0) :
    v = v' ∧ S = Sp := by
  by_cases hvv : v = v'
  · subst hvv
    exact ⟨rfl, Int.eq_of_mul_eq_mul_left (two_pow_sub_three_ne_zero v) H1⟩
  exfalso
  obtain ⟨hv2, hv2'⟩ := both_ge_two v v' S Sp hv hv' hS hSp hvv H1
  have n2 : ((2 : Int) ^ u - 3) * S2 = ((2 : Int) ^ u' - 3) * S2p := by
    rw [Int.sub_mul, Int.sub_mul]
    omega
  have huu : u ≠ u' := by
    intro huu
    subst huu
    have h22 := Int.eq_of_mul_eq_mul_left (two_pow_sub_three_ne_zero u) n2
    have hSS : S = Sp := by omega
    subst hSS
    have hXY : (2 : Int) ^ v = (2 : Int) ^ v' := by
      have := Int.eq_of_mul_eq_mul_left (show S ≠ 0 by omega)
        ((Int.mul_comm _ _).trans (H1.trans (Int.mul_comm _ _)))
      omega
    rcases Nat.lt_or_gt_of_ne hvv with h | h
    · have := two_pow_lt v v' h
      omega
    · have := two_pow_lt v' v h
      omega
  obtain ⟨hu2, hu2'⟩ := both_ge_two u u' S2 S2p hu hu' hS2 hS2p huu n2
  -- the general step, with `u < u'`
  have step : ∀ (v v' u δ : Nat) (S Sp S2 S2p : Int), 2 ≤ v → 2 ≤ v' → 2 ≤ u → 1 ≤ δ →
      0 < S → 0 < Sp → v ≠ v' →
      ((2 : Int) ^ v - 3) * S = ((2 : Int) ^ v' - 3) * Sp →
      3 * S + 1 = (2 : Int) ^ u * S2 → 3 * Sp + 1 = (2 : Int) ^ u * (2 : Int) ^ δ * S2p →
      S - S2 = Sp - S2p → ((2 : Int) ^ v * S - 1) % 3 = 0 → False := by
    intro v v' u δ S Sp S2 S2p hv hv' hu hδ hS hSp hne H1 H3 H3p H4 hA
    obtain ⟨K2, K3⟩ := elim_identities 1 _ _ _ _ S Sp S2 S2p H1 H3 H3p H4
    have hU := two_pow_cases u hu
    have hP : 2 ≤ (2 : Int) ^ δ := by
      have := two_pow_mono 1 δ hδ
      exact this
    rcases Nat.lt_or_gt_of_ne hne with h | h
    · obtain ⟨e, rfl⟩ : ∃ e, v' = v + e := ⟨v' - v, by omega⟩
      rw [Int.pow_add] at K2 K3
      exact plus_lt _ _ _ _ S Sp (two_pow_cases v hv) (two_pow_mod_three v) (by omega) hP
        (two_pow_mod_three δ) (two_pow_cases' e (by omega)) hS hSp K2 K3 hA
    · exact plus_gt _ _ _ _ S (pow_two_ge_four v' hv') (two_pow_lt v' v h)
        (pow_two_ge_four u hu) hP hS K2
  rcases Nat.lt_or_gt_of_ne huu with h | h
  · obtain ⟨δ, rfl⟩ : ∃ δ, u' = u + δ := ⟨u' - u, by omega⟩
    rw [Int.pow_add] at H3p
    exact step v v' u δ S Sp S2 S2p hv2 hv2' hu2 (by omega) hS hSp hvv H1 H3 H3p H4 hA
  · obtain ⟨δ, rfl⟩ : ∃ δ, u = u' + δ := ⟨u - u', by omega⟩
    rw [Int.pow_add] at H3
    exact step v' v u' δ Sp S S2p S2 hv2' hv2 hu2' (by omega) hSp hS (Ne.symm hvv) H1.symm H3p H3
      H4.symm hA'

/-! ## The `3x - 1` sheet -/

/-- `3x - 1`, case `v' < v`: `-E ≤ (P - 1)(X - 3)` but `-E > 3 (P - 1)(X - 3)`. -/
theorem minus_gt (X Y U P S Sp : Int) (hY : 4 ≤ Y) (hYX : Y < X) (hU : 4 ≤ U) (hP : 2 ≤ P)
    (hS : 0 < S) (hSp : 0 < Sp)
    (K2 : S * (P * (U - 3) * (Y - 3) - (U * P - 3) * (X - 3)) = -1 * (P - 1) * (Y - 3))
    (K3 : Sp * (P * (U - 3) * (Y - 3) - (U * P - 3) * (X - 3)) = -1 * (P - 1) * (X - 3)) :
    False := by
  generalize hE : P * (U - 3) * (Y - 3) - (U * P - 3) * (X - 3) = E at K2 K3
  have hneg : -1 * (P - 1) * (Y - 3) < 0 := by
    have : 0 < (P - 1) * (Y - 3) := Int.mul_pos (by omega) (by omega)
    rw [Int.mul_assoc]
    omega
  have hEneg : E < 0 := neg_of_mul_pos_neg S E hS (by rw [K2]; exact hneg)
  have hF : Sp * (-E) = (P - 1) * (X - 3) := by
    rw [Int.mul_neg, K3, Int.mul_assoc]
    omega
  have hFle := le_of_mul_eq_pos Sp (-E) _ (by omega) (by omega) hF
  have m1 : 0 < P * (U - 3) * (X - Y) :=
    Int.mul_pos (Int.mul_pos (by omega) (by omega)) (by omega)
  have m2 : 0 < (P - 1) * (X - 3) := Int.mul_pos (by omega) (by omega)
  pnorm at hE m1 m2 hFle
  omega

/-- `3x - 1`, case `v < v'` (`Y = X Q`): forced to `u = 2`, `Q = 2`, and then
`9(P - 1) < 2^v (2P - 3) ≤ 10(P - 1)` has no solution with `2^v` a power of two. -/
theorem minus_lt (X Q U P S Sp : Int)
    (hX : X = 4 ∨ X = 8 ∨ 16 ≤ X) (hX3 : X % 3 = 1 ∨ X % 3 = 2)
    (hU : U = 4 ∨ 8 ≤ U) (hP : 2 ≤ P) (hP3 : P % 3 = 1 ∨ P % 3 = 2)
    (hQ : Q = 2 ∨ Q = 4 ∨ 8 ≤ Q) (hS : 0 < S)
    (K2 : S * (P * (U - 3) * (X * Q - 3) - (U * P - 3) * (X - 3)) = -1 * (P - 1) * (X * Q - 3))
    (K3 : Sp * (P * (U - 3) * (X * Q - 3) - (U * P - 3) * (X - 3)) = -1 * (P - 1) * (X - 3)) :
    False := by
  have hX4 : 4 ≤ X := by omega
  have hQ2 : 2 ≤ Q := by omega
  have hU4 : 4 ≤ U := by omega
  generalize hE : P * (U - 3) * (X * Q - 3) - (U * P - 3) * (X - 3) = E at K2 K3
  have hXQ : 8 ≤ X * Q := by
    have := Int.mul_le_mul hX4 hQ2 (by decide) (by omega)
    omega
  have hEneg : E < 0 := by
    apply neg_of_mul_pos_neg S E hS
    rw [K2]
    have : 0 < (P - 1) * (X * Q - 3) := Int.mul_pos (by omega) (by omega)
    rw [Int.mul_assoc]
    omega
  -- `(U - 3)(Q - 1) ≤ 2`
  have hW : (U - 3) * (Q - 1) ≤ 2 := by
    by_cases hW : 3 ≤ (U - 3) * (Q - 1)
    · exfalso
      have m := Int.mul_le_mul_of_nonneg_left hW (show 0 ≤ P * X from Int.mul_nonneg (by omega) (by omega))
      have hE' := hE
      pnorm at m hE'
      omega
    · omega
  have hU3 : U - 3 ≤ 2 := by
    have m := Int.mul_le_mul_of_nonneg_left (show 1 ≤ Q - 1 by omega) (show 0 ≤ U - 3 by omega)
    rw [Int.mul_one] at m
    omega
  have hQ3 : Q - 1 ≤ 2 := by
    have m := Int.mul_le_mul_of_nonneg_right (show 1 ≤ U - 3 by omega) (show 0 ≤ Q - 1 by omega)
    rw [Int.one_mul] at m
    omega
  have hU' : U = 4 := by omega
  have hQ' : Q = 2 := by omega
  subst hU' hQ'
  have hc : (S - 2 * Sp) * (-E) = 3 * (P - 1) := by
    pnorm at K2 K3 ⊢
    omega
  have hEX : -E = X * (2 * P - 3) - 9 * (P - 1) := by
    rw [← hE]
    pnorm
    omega
  have hZ := mod_three_prod X P hX3 hP3
  have hE3 : (-E) % 3 = 1 ∨ (-E) % 3 = 2 := by omega
  have hcpos : 0 < S - 2 * Sp := pos_of_mul_pos _ (-E) (by omega) (by rw [hc]; omega)
  have hc3 := three_dvd_of_mul _ (-E) _ hE3 hc
  obtain ⟨c', hc'⟩ : ∃ c', S - 2 * Sp = 3 * c' := ⟨(S - 2 * Sp) / 3, by omega⟩
  rw [hc', Int.mul_assoc] at hc
  have hc'E : c' * (-E) = P - 1 := by omega
  have hFle : -E ≤ P - 1 := le_of_mul_eq_pos c' (-E) _ (by omega) (by omega) hc'E
  rcases hX with rfl | rfl | hX16
  · omega
  · omega
  · have m := Int.mul_le_mul_of_nonneg_right hX16 (show 0 ≤ 2 * P - 3 by omega)
    omega

/-- **Core of Theorem 4.5 (`3x - 1`).** -/
theorem pair_core_minus (v v' u u' : Nat) (S Sp S2 S2p : Int)
    (hv : 1 ≤ v) (hv' : 1 ≤ v') (hu : 1 ≤ u) (hu' : 1 ≤ u')
    (hS : 0 < S) (hSp : 0 < Sp) (hS2 : 0 < S2) (hS2p : 0 < S2p)
    (H1 : ((2 : Int) ^ v - 3) * S = ((2 : Int) ^ v' - 3) * Sp)
    (H3 : 3 * S - 1 = (2 : Int) ^ u * S2) (H3p : 3 * Sp - 1 = (2 : Int) ^ u' * S2p)
    (H4 : S - S2 = Sp - S2p) :
    v = v' ∧ S = Sp := by
  by_cases hvv : v = v'
  · subst hvv
    exact ⟨rfl, Int.eq_of_mul_eq_mul_left (two_pow_sub_three_ne_zero v) H1⟩
  exfalso
  obtain ⟨hv2, hv2'⟩ := both_ge_two v v' S Sp hv hv' hS hSp hvv H1
  have n2 : ((2 : Int) ^ u - 3) * S2 = ((2 : Int) ^ u' - 3) * S2p := by
    rw [Int.sub_mul, Int.sub_mul]
    omega
  have huu : u ≠ u' := by
    intro huu
    subst huu
    have h22 := Int.eq_of_mul_eq_mul_left (two_pow_sub_three_ne_zero u) n2
    have hSS : S = Sp := by omega
    subst hSS
    have hXY : (2 : Int) ^ v = (2 : Int) ^ v' := by
      have := Int.eq_of_mul_eq_mul_left (show S ≠ 0 by omega)
        ((Int.mul_comm _ _).trans (H1.trans (Int.mul_comm _ _)))
      omega
    rcases Nat.lt_or_gt_of_ne hvv with h | h
    · have := two_pow_lt v v' h
      omega
    · have := two_pow_lt v' v h
      omega
  obtain ⟨hu2, hu2'⟩ := both_ge_two u u' S2 S2p hu hu' hS2 hS2p huu n2
  have step : ∀ (v v' u δ : Nat) (S Sp S2 S2p : Int), 2 ≤ v → 2 ≤ v' → 2 ≤ u → 1 ≤ δ →
      0 < S → 0 < Sp → v ≠ v' →
      ((2 : Int) ^ v - 3) * S = ((2 : Int) ^ v' - 3) * Sp →
      3 * S + -1 = (2 : Int) ^ u * S2 → 3 * Sp + -1 = (2 : Int) ^ u * (2 : Int) ^ δ * S2p →
      S - S2 = Sp - S2p → False := by
    intro v v' u δ S Sp S2 S2p hv hv' hu hδ hS hSp hne H1 H3 H3p H4
    obtain ⟨K2, K3⟩ := elim_identities (-1) _ _ _ _ S Sp S2 S2p H1 H3 H3p H4
    have hU := two_pow_cases u hu
    have hP : 2 ≤ (2 : Int) ^ δ := two_pow_mono 1 δ hδ
    rcases Nat.lt_or_gt_of_ne hne with h | h
    · obtain ⟨e, rfl⟩ : ∃ e, v' = v + e := ⟨v' - v, by omega⟩
      rw [Int.pow_add] at K2 K3
      exact minus_lt _ _ _ _ S Sp (two_pow_cases v hv) (two_pow_mod_three v) (by omega) hP
        (two_pow_mod_three δ) (two_pow_cases' e (by omega)) hS K2 K3
    · exact minus_gt _ _ _ _ S Sp (pow_two_ge_four v' hv') (two_pow_lt v' v h)
        (pow_two_ge_four u hu) hP hS hSp K2 K3
  rcases Nat.lt_or_gt_of_ne huu with h | h
  · obtain ⟨δ, rfl⟩ : ∃ δ, u' = u + δ := ⟨u' - u, by omega⟩
    rw [Int.pow_add] at H3p
    exact step v v' u δ S Sp S2 S2p hv2 hv2' hu2 (by omega) hS hSp hvv H1 (by omega) (by omega) H4
  · obtain ⟨δ, rfl⟩ : ∃ δ, u = u' + δ := ⟨u - u', by omega⟩
    rw [Int.pow_add] at H3
    exact step v' v u' δ Sp S S2p S2 hv2' hv2 hu2' (by omega) hSp hS (Ne.symm hvv) H1.symm
      (by omega) (by omega) H4.symm

/-! ## The theorems on the odd naturals -/

/-- Casting a decomposition `n = 2^v · s` to `Int`. -/
theorem cast_decomp (n v s : Nat) (h : n = 2 ^ v * s) : (n : Int) = (2 : Int) ^ v * (s : Int) := by
  have e := congrArg (fun x : Nat => (x : Int)) h
  simp only [Int.natCast_mul, Int.natCast_pow] at e
  exact e

theorem v2_pos_of_odd (n : Nat) (hn : n % 2 = 0) (hpos : 0 < n) : 1 ≤ v2 n := by
  obtain ⟨h1, h2⟩ := decomp n hpos
  cases hv : v2 n with
  | zero =>
    rw [hv, Nat.pow_zero, Nat.one_mul] at h1
    omega
  | succ w => omega

/-- **Theorem 4.5 (`3x + 1`).** Two consecutive drops determine the odd `A ≥ 1`. -/
theorem drop_pair_injective (A A' : Nat) (hA : A % 2 = 1) (hA' : A' % 2 = 1)
    (h1 : drop A = drop A') (h2 : drop (syr A) = drop (syr A')) : A = A' := by
  obtain ⟨d1, o1⟩ := syr_decomp A
  obtain ⟨d1', o1'⟩ := syr_decomp A'
  obtain ⟨d2, o2⟩ := syr_decomp (syr A)
  obtain ⟨d2', o2'⟩ := syr_decomp (syr A')
  have c1 := cast_decomp _ _ _ d1
  have c1' := cast_decomp _ _ _ d1'
  have c2 := cast_decomp _ _ _ d2
  have c2' := cast_decomp _ _ _ d2'
  simp only [Int.natCast_add, Int.natCast_mul] at c1 c1' c2 c2'
  unfold drop at h1 h2
  have e1 : (A : Int) - syr A = A' - syr A' := by omega
  have e2 : (syr A : Int) - syr (syr A) = syr A' - syr (syr A') := by omega
  have hv := v2_pos_of_odd (3 * A + 1) (by omega) (by omega)
  have hv' := v2_pos_of_odd (3 * A' + 1) (by omega) (by omega)
  have hu := v2_pos_of_odd (3 * syr A + 1) (by omega) (by omega)
  have hu' := v2_pos_of_odd (3 * syr A' + 1) (by omega) (by omega)
  have H1 : ((2 : Int) ^ v2 (3 * A + 1) - 3) * (syr A : Int) =
      ((2 : Int) ^ v2 (3 * A' + 1) - 3) * (syr A' : Int) := by
    rw [Int.sub_mul, Int.sub_mul]
    omega
  obtain ⟨hvv, hSS⟩ := pair_core_plus _ _ _ _ (syr A : Int) (syr A' : Int) (syr (syr A) : Int)
    (syr (syr A') : Int) hv hv' hu hu' (by omega) (by omega) (by omega) (by omega) H1
    (by omega) (by omega) e2 (by omega) (by omega)
  omega

/-- The `3x - 1` map `T(A) = oddpart(3A - 1)` on odd `A ≥ 1`. -/
def syrM (A : Nat) : Nat := oddPart (3 * A - 1)

/-- Its drop `(A - T(A)) / 2`. -/
def dropM (A : Nat) : Int := ((A : Int) - (syrM A : Int)) / 2

theorem syrM_decomp (A : Nat) (hA : A % 2 = 1) :
    3 * A - 1 = 2 ^ v2 (3 * A - 1) * syrM A ∧ syrM A % 2 = 1 :=
  decomp (3 * A - 1) (by omega)

/-- **Theorem 4.5 (`3x - 1`).** Two consecutive drops determine the odd `A ≥ 1`. -/
theorem dropM_pair_injective (A A' : Nat) (hA : A % 2 = 1) (hA' : A' % 2 = 1)
    (h1 : dropM A = dropM A') (h2 : dropM (syrM A) = dropM (syrM A')) : A = A' := by
  obtain ⟨d1, o1⟩ := syrM_decomp A hA
  obtain ⟨d1', o1'⟩ := syrM_decomp A' hA'
  obtain ⟨d2, o2⟩ := syrM_decomp (syrM A) o1
  obtain ⟨d2', o2'⟩ := syrM_decomp (syrM A') o1'
  have c1 := cast_decomp _ _ _ d1
  have c1' := cast_decomp _ _ _ d1'
  have c2 := cast_decomp _ _ _ d2
  have c2' := cast_decomp _ _ _ d2'
  have s1 : ((3 * A - 1 : Nat) : Int) = 3 * (A : Int) - 1 := by omega
  have s1' : ((3 * A' - 1 : Nat) : Int) = 3 * (A' : Int) - 1 := by omega
  have s2 : ((3 * syrM A - 1 : Nat) : Int) = 3 * (syrM A : Int) - 1 := by omega
  have s2' : ((3 * syrM A' - 1 : Nat) : Int) = 3 * (syrM A' : Int) - 1 := by omega
  rw [s1] at c1
  rw [s1'] at c1'
  rw [s2] at c2
  rw [s2'] at c2'
  unfold dropM at h1 h2
  have e1 : (A : Int) - syrM A = A' - syrM A' := by omega
  have e2 : (syrM A : Int) - syrM (syrM A) = syrM A' - syrM (syrM A') := by omega
  have hv := v2_pos_of_odd (3 * A - 1) (by omega) (by omega)
  have hv' := v2_pos_of_odd (3 * A' - 1) (by omega) (by omega)
  have hu := v2_pos_of_odd (3 * syrM A - 1) (by omega) (by omega)
  have hu' := v2_pos_of_odd (3 * syrM A' - 1) (by omega) (by omega)
  have H1 : ((2 : Int) ^ v2 (3 * A - 1) - 3) * (syrM A : Int) =
      ((2 : Int) ^ v2 (3 * A' - 1) - 3) * (syrM A' : Int) := by
    rw [Int.sub_mul, Int.sub_mul]
    omega
  obtain ⟨hvv, hSS⟩ := pair_core_minus _ _ _ _ (syrM A : Int) (syrM A' : Int)
    (syrM (syrM A) : Int) (syrM (syrM A') : Int) hv hv' hu hu' (by omega) (by omega) (by omega)
    (by omega) H1 (by omega) (by omega) e2
  omega

end ProcgenSelfieEdim
