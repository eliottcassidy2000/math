/-
Rectangle-only recursion over the pinned Rudich--Rousseau engine.
The rectangle identity `sum_Zι` is inherited, not a new theorem here.
This module proves a new adapter and all-depth semantic recurrence;
it does not replace the external TM2 compiler or prove a new time bound.
-/
module

public import FactoringAux.EngineGroup

@[expose] public section

namespace Factoring.Aux.RectangleOnly

open Finset

variable {m G lam : ℕ}

/-- Extend a rectangle choice by zero on all unused quadruples. -/
def extendChoice (ι : Choice G lam) : Quad G → Fin (tB lam + 1) :=
  fun r => if h : r ∈ rects G then ι ⟨r, h⟩ else 0

@[simp] theorem restrict_extend (ι : Choice G lam) :
    restrictC (extendChoice ι) = ι := by
  funext r
  simp [restrictC, extendChoice, r.2]

def choiceCode (ι : Choice G lam) :
    Fin ((tB lam + 1) ^ (G * G * (G * G))) :=
  phiEquiv G (tB lam) (extendChoice ι)

@[simp] theorem decode_choiceCode (ι : Choice G lam) :
    restrictC ((phiEquiv G (tB lam)).symm (choiceCode ι)) = ι := by
  simp [choiceCode]

theorem choiceCode_injective : Function.Injective (@choiceCode G lam) := by
  intro a b h
  have := congrArg (fun c => restrictC ((phiEquiv G (tB lam)).symm c)) h
  simpa using this

/-- The entire child table, not only its scalar value, ignores unused choices. -/
theorem nextT_same_on_fibres (K : ℕ) (hK : 0 < K)
    (T : Clique (ZMod m) (G * G) lam)
    (c d : Fin ((tB lam + 1) ^ (G * G * (G * G))))
    (h : restrictC ((phiEquiv G (tB lam)).symm c) =
         restrictC ((phiEquiv G (tB lam)).symm d)) :
    nextT K hK T c = nextT K hK T d := by
  unfold nextT
  rw [h]

def next (K : ℕ) (hK : 0 < K)
    (T : Clique (ZMod m) (G * G) lam) (ι : Choice G lam) :=
  nextT K hK T (choiceCode ι)

theorem next_value (K : ℕ) (hK : 2 ≤ K) (hG : G = 3 * K ^ 8)
    (T : Clique (ZMod m) (G * G) lam) (ι : Choice G lam) :
    (next K (by omega) T ι).value = Zι (toGrid T) ι := by
  unfold next
  rw [value_nextT K hK hG]
  simp

theorem sum_next (K : ℕ) (hK : 2 ≤ K) (hG : G = 3 * K ^ 8)
    (T : Clique (ZMod m) (G * G) lam) :
    ∑ ι : Choice G lam, (next K (by omega) T ι).value =
      (2 : ZMod m) ^ (tB lam * (rects G).card) * T.value := by
  simp_rw [next_value K hK hG T]
  rw [← sum_Zι, ← value_toGrid]

theorem choice_card (G lam : ℕ) :
    Fintype.card (Choice G lam) = (tB lam + 1) ^ (rects G).card := by
  simp [Choice]

def exponent (K0 : ℕ) : (d L : ℕ) → ℕ
  | 0, _ => 0
  | d + 1, L => tB L * (rects (Gs K0 d)).card +
      exponent K0 d (Lnext (Kb K0 d) (Gs K0 d) L)

def padding (K0 : ℕ) : (d L : ℕ) → ℕ
  | 0, _ => 1
  | d + 1, L => (tB L + 1) ^ nonRects (Gs K0 d) *
      padding K0 d (Lnext (Kb K0 d) (Gs K0 d) L)

def value (K0 : ℕ) (hK : 2 ≤ K0) :
    (d L : ℕ) → Clique (ZMod m) (Kb K0 d) L → ZMod m
  | 0, _, T => T.value
  | d + 1, L, T => ∑ ι : Choice (Gs K0 d) L,
      value K0 hK d (Lnext (Kb K0 d) (Gs K0 d) L)
        (next (Kb K0 d) (by have := Kb_ge K0 hK d; omega) T ι)

theorem value_eq (K0 : ℕ) (hK : 2 ≤ K0) :
    ∀ (d L : ℕ) (T : Clique (ZMod m) (Kb K0 d) L),
      value K0 hK d L T = (2 : ZMod m) ^ exponent K0 d L * T.value
  | 0, _, T => by simp [value, exponent]
  | d + 1, L, T => by
      simp only [value, exponent]
      simp only [value_eq K0 hK d]
      rw [← Finset.mul_sum]
      rw [sum_next (m := m) (G := Gs K0 d) (Kb K0 d) (Kb_ge K0 hK d) rfl T]
      rw [pow_add]
      ring_nf
      rfl

theorem mult_factorization (K0 : ℕ) : ∀ d L,
    deepMult K0 d L = padding K0 d L * 2 ^ exponent K0 d L
  | 0, _ => by simp [deepMult, padding, exponent]
  | d + 1, L => by
      simp only [deepMult, padding, exponent, mult, mult_factorization K0 d, pow_add]
      ring

theorem old_value_eq_padding (K0 : ℕ) (hK : 2 ≤ K0)
    (d L : ℕ) (T : Clique (ZMod m) (Kb K0 d) L) :
    deepVal K0 hK d L T = (padding K0 d L : ZMod m) * value K0 hK d L T := by
  rw [deepVal_eq, value_eq, mult_factorization]
  push_cast
  ring

/-- Boundary assignments remain named; this is an equality of whole functions. -/
theorem boundary_values {B : Type} (K0 : ℕ) (hK : 2 ≤ K0)
    (d L : ℕ) (T : B → Clique (ZMod m) (Kb K0 d) L) :
    (fun b => value K0 hK d L (T b)) =
      (fun b => (2 : ZMod m) ^ exponent K0 d L * (T b).value) := by
  funext b
  exact value_eq K0 hK d L (T b)

/-- Number of leaves in the unpruned nested sums; not a TM2 time bound. -/
def fullLeaves (K0 : ℕ) : (d L : ℕ) → ℕ
  | 0, _ => 1
  | d + 1, L => (tB L + 1) ^ (Gs K0 d * Gs K0 d * (Gs K0 d * Gs K0 d)) *
      fullLeaves K0 d (Lnext (Kb K0 d) (Gs K0 d) L)

def reducedLeaves (K0 : ℕ) : (d L : ℕ) → ℕ
  | 0, _ => 1
  | d + 1, L => (tB L + 1) ^ (rects (Gs K0 d)).card *
      reducedLeaves K0 d (Lnext (Kb K0 d) (Gs K0 d) L)

theorem quad_partition (G : ℕ) :
    nonRects G + (rects G).card = G * G * (G * G) := by
  have h := Finset.card_filter_add_card_filter_not
    (s := (Finset.univ : Finset (Quad G))) (fun r => r ∈ rects G)
  have hf : (univ.filter (fun r : Quad G => r ∈ rects G)) = rects G := by
    ext r
    simp
  rw [hf] at h
  simpa [nonRects, Quad, add_comm] using h

theorem leaves_eq_padding (K0 : ℕ) : ∀ d L,
    fullLeaves K0 d L = padding K0 d L * reducedLeaves K0 d L
  | 0, _ => by simp [fullLeaves, padding, reducedLeaves]
  | d + 1, L => by
      simp only [fullLeaves, padding, reducedLeaves, leaves_eq_padding K0 d]
      rw [← quad_partition, pow_add]
      ring

end Factoring.Aux.RectangleOnly

#print axioms Factoring.Aux.RectangleOnly.nextT_same_on_fibres
#print axioms Factoring.Aux.RectangleOnly.choiceCode_injective
#print axioms Factoring.Aux.RectangleOnly.sum_next
#print axioms Factoring.Aux.RectangleOnly.choice_card
#print axioms Factoring.Aux.RectangleOnly.value_eq
#print axioms Factoring.Aux.RectangleOnly.old_value_eq_padding
#print axioms Factoring.Aux.RectangleOnly.boundary_values
#print axioms Factoring.Aux.RectangleOnly.leaves_eq_padding
