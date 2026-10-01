import ProcgenSelfieEdim.ListLemmas
import ProcgenSelfieEdim.Tournament

set_option autoImplicit false

/-!
# Coding labeled tournaments by natural numbers

For `a < b`, bit number `b (b - 1) / 2 + a` (`= pairs b + a`) of a code `c` says
whether `a → b` in the tournament `tourC c`. `tourC_code` shows that every labeled
tournament on `{0, …, n-1}` is `tourC c` (on all ordered pairs of distinct vertices)
for some `c < 2 ^ C(n, 2)`. Finite statements about all tournaments on `n` vertices
therefore reduce to checks over the codes `c < 2 ^ C(n, 2)`.
-/

namespace ProcgenSelfieEdim

/-- Bit `i` of `c`. -/
def bitB (c i : Nat) : Bool := Nat.beq (c / 2 ^ i % 2) 1

/-- The tournament coded by `c`: for `a < b`, `a → b` iff bit `b (b-1)/2 + a` of `c`
is set; no loops. -/
def tourC (c a b : Nat) : Bool :=
  cond (Nat.blt a b) (bitB c (b * (b - 1) / 2 + a))
    (cond (Nat.blt b a) (!bitB c (a * (a - 1) / 2 + b)) false)

theorem tourC_lt (c a b : Nat) (h : a < b) : tourC c a b = bitB c (b * (b - 1) / 2 + a) := by
  unfold tourC
  rw [Nat.blt_eq.mpr h, cond_true]

theorem tourC_gt (c a b : Nat) (h : b < a) : tourC c a b = !bitB c (a * (a - 1) / 2 + b) := by
  unfold tourC
  have h1 : Nat.blt a b = false := by
    cases hb : Nat.blt a b
    · rfl
    · exact absurd (Nat.blt_eq.mp hb) (by omega)
  rw [h1, cond_false, Nat.blt_eq.mpr h, cond_true]

theorem tourC_self (c a : Nat) : tourC c a a = false := by
  unfold tourC
  have h1 : Nat.blt a a = false := by
    cases hb : Nat.blt a a
    · rfl
    · exact absurd (Nat.blt_eq.mp hb) (Nat.lt_irrefl a)
  rw [h1, cond_false, cond_false]

theorem tourC_isTournament (n c : Nat) : IsTournament n (tourC c) := by
  intro a b _ _ hab
  by_cases h : a < b
  · rw [tourC_lt c a b h, tourC_gt c b a h]
  · rw [tourC_gt c a b (by omega), tourC_lt c b a (by omega), Bool.not_not]

/-! ## Bit lists -/

/-- Indexing into a bit list (out-of-range bits are `false`). -/
def nthB : List Bool → Nat → Bool
  | [], _ => false
  | x :: _, 0 => x
  | _ :: l, i + 1 => nthB l i

theorem nthB_append_left : ∀ (l1 l2 : List Bool) (i : Nat), i < l1.length →
    nthB (l1 ++ l2) i = nthB l1 i
  | [], _, _, h => absurd h (Nat.not_lt_zero _)
  | _ :: _, _, 0, _ => rfl
  | _ :: l1, l2, i + 1, h => by
    show nthB (l1 ++ l2) i = nthB l1 i
    exact nthB_append_left l1 l2 i (by simp only [List.length_cons] at h; omega)

theorem nthB_append_right : ∀ (l1 l2 : List Bool) (i : Nat), nthB (l1 ++ l2) (l1.length + i) = nthB l2 i
  | [], _, i => by rw [List.nil_append, List.length_nil, Nat.zero_add]
  | _ :: l1, l2, i => by
    rw [List.cons_append, List.length_cons, Nat.add_right_comm]
    exact nthB_append_right l1 l2 i

/-- Horner encoding of a bit list (first bit = least significant). -/
def bitsToNat : List Bool → Nat
  | [] => 0
  | x :: l => (if x = true then 1 else 0) + 2 * bitsToNat l

theorem bitB_bitsToNat : ∀ (l : List Bool) (i : Nat), bitB (bitsToNat l) i = nthB l i
  | [], i => by
    show Nat.beq (0 / 2 ^ i % 2) 1 = false
    rw [Nat.zero_div]
    rfl
  | x :: l, 0 => by
    show Nat.beq (((if x = true then 1 else 0) + 2 * bitsToNat l) / 2 ^ 0 % 2) 1 = x
    rw [Nat.pow_zero, Nat.div_one]
    cases x
    · rw [if_neg (by decide)]
      have : (0 + 2 * bitsToNat l) % 2 = 0 := by omega
      rw [this]
      rfl
    · rw [if_pos rfl]
      have : (1 + 2 * bitsToNat l) % 2 = 1 := by omega
      rw [this]
      rfl
  | x :: l, i + 1 => by
    show Nat.beq (((if x = true then 1 else 0) + 2 * bitsToNat l) / 2 ^ (i + 1) % 2) 1 = nthB l i
    rw [← bitB_bitsToNat l i]
    unfold bitB
    rw [show 2 ^ (i + 1) = 2 * 2 ^ i from by rw [Nat.pow_succ, Nat.mul_comm],
      ← Nat.div_div_eq_div_mul]
    have : ((if x = true then 1 else 0) + 2 * bitsToNat l) / 2 = bitsToNat l := by
      split <;> omega
    rw [this]

theorem bitsToNat_lt : ∀ (l : List Bool), bitsToNat l < 2 ^ l.length
  | [] => by decide
  | x :: l => by
    have ih := bitsToNat_lt l
    show (if x = true then 1 else 0) + 2 * bitsToNat l < 2 ^ (l.length + 1)
    rw [Nat.pow_succ]
    split <;> omega

/-! ## The code of a tournament -/

/-- `[f 0, f 1, …, f (k-1)]`. -/
def ascList (f : Nat → Bool) : Nat → List Bool
  | 0 => []
  | k + 1 => ascList f k ++ [f k]

theorem ascList_length (f : Nat → Bool) : ∀ (k : Nat), (ascList f k).length = k
  | 0 => rfl
  | k + 1 => by rw [ascList, List.length_append, ascList_length f k]; rfl

theorem nthB_ascList (f : Nat → Bool) : ∀ (k a : Nat), a < k → nthB (ascList f k) a = f a
  | 0, _, h => absurd h (Nat.not_lt_zero _)
  | k + 1, a, h => by
    rw [ascList]
    by_cases ha : a < k
    · rw [nthB_append_left _ _ a (by rw [ascList_length]; exact ha), nthB_ascList f k a ha]
    · have hak : a = k := by omega
      subst hak
      have := nthB_append_right (ascList f a) [f a] 0
      rw [ascList_length, Nat.add_zero] at this
      rw [this]
      rfl

/-- The bits `T a b` (`a < b < n`) in the order of `pairs b + a`. -/
def codeList (T : Nat → Nat → Bool) : Nat → List Bool
  | 0 => []
  | b + 1 => codeList T b ++ ascList (fun a => T a b) b

theorem codeList_length (T : Nat → Nat → Bool) : ∀ (n : Nat), (codeList T n).length = pairs n
  | 0 => rfl
  | n + 1 => by
    rw [codeList, List.length_append, codeList_length T n, ascList_length]
    rfl

theorem pairs_mono : ∀ (a b : Nat), a ≤ b → pairs a ≤ pairs b
  | a, 0, h => by
    have : a = 0 := by omega
    subst this
    exact Nat.le_refl _
  | a, b + 1, h => by
    by_cases hab : a = b + 1
    · subst hab
      exact Nat.le_refl _
    · have := pairs_mono a b (by omega)
      show pairs a ≤ pairs b + b
      omega

theorem codeList_nthB (T : Nat → Nat → Bool) :
    ∀ (n a b : Nat), a < b → b < n → nthB (codeList T n) (pairs b + a) = T a b
  | 0, _, _, _, hb => absurd hb (Nat.not_lt_zero _)
  | n + 1, a, b, hab, hb => by
    rw [codeList]
    by_cases hbn : b < n
    · have hlt : pairs b + a < pairs n := by
        have := pairs_mono (b + 1) n (by omega)
        show pairs b + a < pairs n
        have h2 : pairs (b + 1) = pairs b + b := rfl
        omega
      rw [nthB_append_left _ _ _ (by rw [codeList_length]; exact hlt), codeList_nthB T n a b hab hbn]
    · have hbn' : b = n := by omega
      subst hbn'
      have := nthB_append_right (codeList T b) (ascList (fun a => T a b) b) a
      rw [codeList_length] at this
      rw [this, nthB_ascList _ b a hab]

theorem pairs_formula (b : Nat) : b * (b - 1) / 2 = pairs b := by
  have := pairs_eq b
  omega

/-- **Every labeled tournament has a code.** -/
theorem tourC_code (n : Nat) (T : Nat → Nat → Bool) (hT : IsTournament n T) :
    bitsToNat (codeList T n) < 2 ^ pairs n ∧ SameOn n T (tourC (bitsToNat (codeList T n))) := by
  refine ⟨by have := bitsToNat_lt (codeList T n); rwa [codeList_length] at this, ?_⟩
  intro a b ha hb hab
  by_cases h : a < b
  · rw [tourC_lt _ a b h, pairs_formula, bitB_bitsToNat, codeList_nthB T n a b h hb]
  · rw [tourC_gt _ a b (by omega), pairs_formula, bitB_bitsToNat,
      codeList_nthB T n b a (by omega) ha]
    exact hT b a hb ha (Ne.symm hab)

end ProcgenSelfieEdim
