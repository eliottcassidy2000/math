import ProcgenSelfieEdim.DropPair
import ProcgenSelfieEdim.SwitchSum

set_option autoImplicit false
set_option linter.unusedSimpArgs false

/-!
# THM-4530 (c): unit branches of the `qx + 1` drop map (Theorem 5.2)

For odd `q`, the `qx + 1` map on odd `A ≥ 1` is `T(A) = oddpart(qA + 1)`, with drop
`K(A) = (A - T(A)) / 2` and `2qK + 1 = (2^v - q) T` (`dropQ_identity`, Theorem 5.1, one step).

* `dropQ_unit_desc`: if `q = 2^a - 1`, every `d ≥ 0` is the drop of `2^(a+1) d + 1`.
* `dropQ_unit_asc`: if `q = 2^a + 1` (`a ≥ 1`), every `d = -n < 0` is the drop of `2^(a+1) n - 1`.
* `unit_of_zero_drop`: if `0` is a drop, then `q = 2^a - 1` (a drop `0` is a fixed point,
  `(2^v - q) A = 1`).
* `unit_of_neg_drop`: if `-q!` is a drop, then `q = 2^a + 1` (the ascending modulus
  `q - 2^v ≤ q` divides `q!`, hence divides `2q·q! - (q - 2^v) T = 1`).
* `units_nonneg_iff`, `units_neg_iff`, `units_all_iff` (**Theorem 5.2**): every `d ≥ 0` is a
  drop iff `q = 2^a - 1`; every `d < 0` iff `q = 2^a + 1`; every `d ∈ Z` iff `q = 3`.

The note's density statement (a side without a unit misses a set of density `> 0.38`) is not
formalized; the "only if" directions here use single explicit witnesses (`d = 0`, `d = -q!`).
-/

namespace ProcgenSelfieEdim

/-- The `qx + 1` map: `T(A) = oddpart(qA + 1)`. -/
def syrQ (q A : Nat) : Nat := oddPart (q * A + 1)

/-- Its drop `(A - T(A)) / 2`. -/
def dropQ (q A : Nat) : Int := ((A : Int) - (syrQ q A : Int)) / 2

theorem syrQ_decomp (q A : Nat) :
    q * A + 1 = 2 ^ v2 (q * A + 1) * syrQ q A ∧ syrQ q A % 2 = 1 :=
  decomp (q * A + 1) (by omega)

theorem odd_mul_odd (q A : Nat) (hq : q % 2 = 1) (hA : A % 2 = 1) : q * A % 2 = 1 := by
  rw [Nat.mul_mod, hq, hA]

/-- **Theorem 5.1 (one step).** `2 q K(A) + 1 = (2^v - q) T(A)`, `v = v_2(qA + 1)`. -/
theorem dropQ_identity (q A : Nat) (hA : A % 2 = 1) :
    2 * (q : Int) * dropQ q A + 1 = ((2 : Int) ^ v2 (q * A + 1) - q) * (syrQ q A : Int) := by
  obtain ⟨h1, h2⟩ := syrQ_decomp q A
  have c := cast_decomp _ _ _ h1
  pcast at c
  have hK : 2 * dropQ q A = (A : Int) - syrQ q A := by
    unfold dropQ
    omega
  have e := congrArg (fun z => (q : Int) * z) hK
  simp only at e
  pnorm at e c ⊢
  omega

/-- **Descending unit.** If `q + 1 = 2^a`, then `d ≥ 0` is the drop of `2^(a+1) d + 1`
(branch `v = a`, `T = 2qd + 1`). -/
theorem dropQ_unit_desc (q a d : Nat) (hqa : q + 1 = 2 ^ a) :
    (2 ^ (a + 1) * d + 1) % 2 = 1 ∧ dropQ q (2 ^ (a + 1) * d + 1) = d := by
  have e : 2 ^ (a + 1) * d = 2 * ((q + 1) * d) := by
    rw [Nat.pow_succ, ← hqa, Nat.mul_comm (q + 1) 2, Nat.mul_assoc]
  rw [e]
  refine ⟨by omega, ?_⟩
  have hq' : (q : Int) + 1 = (2 : Int) ^ a := by
    have := congrArg (fun x : Nat => (x : Int)) hqa
    pcast at this
    exact this
  have key : q * (2 * ((q + 1) * d) + 1) + 1 = 2 ^ a * (2 * (q * d) + 1) := by
    apply Int.ofNat_inj.1
    pcast
    rw [← hq']
    pnorm
    omega
  have hs := (v2_oddPart_of_eq _ _ _ key (by omega)).2
  unfold dropQ syrQ
  rw [hs]
  pcast
  pnorm
  omega

/-- **Ascending unit.** If `q = 2^a + 1` with `a ≥ 1`, then `d = -n < 0` is the drop of
`2^(a+1) n - 1` (branch `v = a`, `T = 2qn - 1`). -/
theorem dropQ_unit_asc (q a n : Nat) (ha : 1 ≤ a) (hqa : q = 2 ^ a + 1) (hn : 1 ≤ n) :
    (2 ^ (a + 1) * n - 1) % 2 = 1 ∧ dropQ q (2 ^ (a + 1) * n - 1) = -(n : Int) := by
  have hR : 2 ≤ 2 ^ a := by
    have := Nat.pow_le_pow_right (show 0 < 2 by decide) ha
    exact this
  have e : 2 ^ (a + 1) * n = 2 * 2 ^ a * n := by rw [Nat.pow_succ, Nat.mul_comm (2 ^ a) 2]
  have hpos : 1 ≤ 2 ^ a * n := Nat.mul_pos (by omega) (by omega)
  obtain ⟨A, hA⟩ : ∃ A, A + 1 = 2 * 2 ^ a * n := ⟨2 * 2 ^ a * n - 1, by rw [Nat.mul_assoc]; omega⟩
  have hAe : 2 ^ (a + 1) * n - 1 = A := by
    rw [e]
    omega
  rw [hAe]
  have hAodd : A % 2 = 1 := by
    rw [Nat.mul_assoc] at hA
    omega
  refine ⟨hAodd, ?_⟩
  have key : q * A + 1 = 2 ^ a * (A + 2 * n) := by
    apply Int.ofNat_inj.1
    have hA' := congrArg (fun x : Nat => (x : Int)) hA
    pcast at hA'
    subst hqa
    pcast
    pnorm at hA' ⊢
    omega
  have hs := (v2_oddPart_of_eq _ _ _ key (by omega)).2
  unfold dropQ syrQ
  rw [hs]
  pcast
  omega

/-- **Only if (descending).** A drop `0` is a fixed point, so `q = 2^v - 1`. -/
theorem unit_of_zero_drop (q A : Nat) (hA : A % 2 = 1) (h : dropQ q A = 0) :
    q + 1 = 2 ^ v2 (q * A + 1) := by
  obtain ⟨h1, h2⟩ := syrQ_decomp q A
  have hT : syrQ q A = A := by
    unfold dropQ at h
    omega
  rw [hT] at h1
  generalize 2 ^ v2 (q * A + 1) = X at h1 ⊢
  -- `X A = q A + 1`, so `(X - q) A = 1` and `A = 1`
  have hlt : q < X := by
    by_cases hqX : X ≤ q
    · have := Nat.mul_le_mul_right A hqX
      omega
    · omega
  have hm : (X - q) * A = 1 := by
    rw [Nat.sub_mul]
    omega
  have hA1 : A = 1 := Nat.eq_one_of_mul_eq_one_right (by rw [Nat.mul_comm]; exact hm)
  subst hA1
  omega

theorem dvd_fact (m : Nat) (hm : 1 ≤ m) : ∀ q, m ≤ q → m ∣ fact q
  | 0, h => absurd h (by omega)
  | q + 1, h => by
    show m ∣ (q + 1) * fact q
    by_cases hmq : m = q + 1
    · subst hmq
      exact ⟨fact q, rfl⟩
    · exact Nat.dvd_mul_left_of_dvd (dvd_fact m hm q (by omega)) (q + 1)

/-- `Nat` divisibility casts to `Int` (core's `Int.natCast_dvd_natCast` uses
`Classical.choice`). -/
theorem natCast_dvd (m n : Nat) (h : m ∣ n) : (m : Int) ∣ (n : Int) := by
  obtain ⟨c, hc⟩ := h
  exact ⟨c, by rw [hc, Int.natCast_mul]⟩

theorem fact_pos' : ∀ q, 1 ≤ fact q
  | 0 => Nat.le_refl 1
  | q + 1 => by
    show 1 ≤ (q + 1) * fact q
    exact Nat.mul_pos (by omega) (fact_pos' q)

/-- **Only if (ascending).** If `-q!` is a drop, then `q = 2^v + 1` with `v ≥ 1`. -/
theorem unit_of_neg_drop (q A : Nat) (hq : q % 2 = 1) (hA : A % 2 = 1)
    (h : dropQ q A = -((fact q : Nat) : Int)) : ∃ a, 1 ≤ a ∧ q = 2 ^ a + 1 := by
  have hid := dropQ_identity q A hA
  obtain ⟨h1, h2⟩ := syrQ_decomp q A
  have hv : 1 ≤ v2 (q * A + 1) := v2_pos_of_odd _ (by have := odd_mul_odd q A hq hA; omega) (by omega)
  rw [h] at hid
  generalize hT : syrQ q A = T at hid h2
  generalize hV : v2 (q * A + 1) = v at hid hv
  have hpv : ((2 ^ v : Nat) : Int) = (2 : Int) ^ v := by pcast
  rw [← hpv] at hid
  generalize hW : 2 ^ v = W at hid
  have hF := fact_pos' q
  generalize hG : fact q = G at hid hF
  -- the left side is negative, so `W = 2^v < q`
  have hprod : 2 * (q : Int) * -((G : Nat) : Int) + 1 < 0 := by
    have : (1 : Int) ≤ (q : Int) * ((G : Nat) : Int) :=
      Int.mul_le_mul (show (1 : Int) ≤ q by omega) (show (1 : Int) ≤ ((G : Nat) : Int) by omega)
        (by decide) (by omega)
    rw [Int.mul_neg, Int.mul_assoc]
    omega
  have hTpos : (0 : Int) < T := by omega
  have hlt : W < q := by
    by_cases hc : q ≤ W
    · have : 0 ≤ ((W : Int) - q) * T := Int.mul_nonneg (by omega) (by omega)
      omega
    · omega
  -- `m = q - W` divides `q!` and `2 q q! - 1`, so `m = 1`
  have hmd : ((q - W : Nat) : Int) ∣ ((G : Nat) : Int) := by
    rw [← hG]
    exact natCast_dvd _ _ (dvd_fact (q - W) (by omega) q (by omega))
  have hm1 : ((q - W : Nat) : Int) ∣ 2 * (q : Int) * ((G : Nat) : Int) - 1 := by
    refine ⟨T, ?_⟩
    have hc : ((q - W : Nat) : Int) = (q : Int) - (W : Int) := by omega
    rw [hc]
    have e := hid
    rw [Int.mul_neg] at e
    pnorm at e ⊢
    omega
  have hm2 : ((q - W : Nat) : Int) ∣ 2 * (q : Int) * ((G : Nat) : Int) :=
    Int.dvd_trans hmd (Int.dvd_mul_left _ _)
  have hone := Int.dvd_sub hm2 hm1
  rw [show 2 * (q : Int) * ((G : Nat) : Int) - (2 * (q : Int) * ((G : Nat) : Int) - 1) = 1
    by omega] at hone
  have := Int.eq_one_of_dvd_one (by omega) hone
  exact ⟨v, hv, by omega⟩

/-! ## Theorem 5.2 -/

/-- **Theorem 5.2 (descending side).** For odd `q`, every `d ≥ 0` is a drop of the `qx + 1`
map on the odd naturals iff `q = 2^a - 1`. -/
theorem units_nonneg_iff (q : Nat) :
    (∀ d : Nat, ∃ A : Nat, A % 2 = 1 ∧ dropQ q A = d) ↔ ∃ a : Nat, q + 1 = 2 ^ a := by
  constructor
  · intro h
    obtain ⟨A, hA, h0⟩ := h 0
    exact ⟨_, unit_of_zero_drop q A hA h0⟩
  · rintro ⟨a, ha⟩ d
    obtain ⟨h1, h2⟩ := dropQ_unit_desc q a d ha
    exact ⟨_, h1, h2⟩

/-- **Theorem 5.2 (ascending side).** For odd `q`, every `d < 0` is a drop iff `q = 2^a + 1`
with `a ≥ 1`. -/
theorem units_neg_iff (q : Nat) (hq : q % 2 = 1) :
    (∀ n : Nat, 1 ≤ n → ∃ A : Nat, A % 2 = 1 ∧ dropQ q A = -(n : Int)) ↔
      ∃ a : Nat, 1 ≤ a ∧ q = 2 ^ a + 1 := by
  constructor
  · intro h
    obtain ⟨A, hA, h0⟩ := h (fact q) (fact_pos' q)
    exact unit_of_neg_drop q A hq hA h0
  · rintro ⟨a, ha, hqa⟩ n hn
    obtain ⟨h1, h2⟩ := dropQ_unit_asc q a n ha hqa hn
    exact ⟨_, h1, h2⟩

/-- **Theorem 5.2.** For odd `q`, the drop map of `qx + 1` hits every integer iff `q = 3`. -/
theorem units_all_iff (q : Nat) (hq : q % 2 = 1) :
    (∀ d : Int, ∃ A : Nat, A % 2 = 1 ∧ dropQ q A = d) ↔ q = 3 := by
  constructor
  · intro h
    obtain ⟨a, ha⟩ := (units_nonneg_iff q).1 (fun d => h d)
    obtain ⟨b, hb, hqb⟩ := (units_neg_iff q hq).1 (fun n _ => h (-(n : Int)))
    -- `2^a = 2^b + 2` forces `b = 1`
    rcases (show b = 1 ∨ 2 ≤ b by omega) with rfl | hb2
    · rw [Nat.pow_one] at hqb
      exact hqb
    · exfalso
      obtain ⟨c, rfl⟩ : ∃ c, b = c + 2 := ⟨b - 2, by omega⟩
      rw [Nat.pow_add] at hqb
      have h4 : (2 : Nat) ^ 2 = 4 := rfl
      rw [h4] at hqb
      rcases (show a ≤ 1 ∨ 2 ≤ a by omega) with ha1 | ha2
      · have : 2 ^ a ≤ 2 ^ 1 := Nat.pow_le_pow_right (by decide) ha1
        have : 1 ≤ 2 ^ c := Nat.one_le_two_pow
        omega
      · obtain ⟨e, rfl⟩ : ∃ e, a = e + 2 := ⟨a - 2, by omega⟩
        rw [Nat.pow_add, h4] at ha
        omega
  · rintro rfl d
    by_cases hd : 0 ≤ d
    · obtain ⟨n, rfl⟩ : ∃ n : Nat, d = n := ⟨d.toNat, by omega⟩
      obtain ⟨h1, h2⟩ := dropQ_unit_desc 3 2 n rfl
      exact ⟨_, h1, h2⟩
    · obtain ⟨n, hn⟩ : ∃ n : Nat, d = -(n : Int) := ⟨(-d).toNat, by omega⟩
      subst hn
      obtain ⟨h1, h2⟩ := dropQ_unit_asc 3 1 n (Nat.le_refl 1) rfl (by omega)
      exact ⟨_, h1, h2⟩

/-- Examples: `5x + 1` misses some `d ≥ 0` (namely `d = 0`; `5 + 1` is not a power of two),
and `7x + 1` misses some `d < 0` (`7 - 1` is not a power of two). -/
theorem units_examples :
    (¬ ∀ d : Nat, ∃ A : Nat, A % 2 = 1 ∧ dropQ 5 A = d) ∧
      (¬ ∀ n : Nat, 1 ≤ n → ∃ A : Nat, A % 2 = 1 ∧ dropQ 7 A = -(n : Int)) := by
  constructor
  · intro h
    obtain ⟨a, ha⟩ := (units_nonneg_iff 5).1 h
    rcases (show a ≤ 2 ∨ 3 ≤ a by omega) with h' | h'
    · rcases (show a = 0 ∨ a = 1 ∨ a = 2 by omega) with rfl | rfl | rfl <;>
        exact absurd ha (by decide)
    · have : 2 ^ 3 ≤ 2 ^ a := Nat.pow_le_pow_right (by decide) h'
      have h8 : (2 : Nat) ^ 3 = 8 := rfl
      omega
  · intro h
    obtain ⟨a, _, ha⟩ := (units_neg_iff 7 rfl).1 h
    rcases (show a ≤ 2 ∨ 3 ≤ a by omega) with h' | h'
    · rcases (show a = 0 ∨ a = 1 ∨ a = 2 by omega) with rfl | rfl | rfl <;>
        exact absurd ha (by decide)
    · have : 2 ^ 3 ≤ 2 ^ a := Nat.pow_le_pow_right (by decide) h'
      have h8 : (2 : Nat) ^ 3 = 8 := rfl
      omega

end ProcgenSelfieEdim
