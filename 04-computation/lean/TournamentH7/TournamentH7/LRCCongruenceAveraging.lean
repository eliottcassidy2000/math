import Mathlib

/-!
# The congruence-averaging lemma: the arithmetic core (kind-pasteur-S128c44)

THM-952's hardened tree rests on the CONGRUENCE-AVERAGING LEMMA: the near-pole
residues r(t) sweep a rotation orbit, and their reciprocal sum over a period is
harmonic, not linear.  The arithmetic content is exactly two facts, formalized
here over ℚ (no logs — the exact form is stronger than the logged bound):

1. **Orbit invariance** (`lavSum_mul_unit`): for a UNIT κ on `ZMod n`, the sum
   of `1/max(1, leastAbsVal (κ*t))` over all `t` equals the κ-free sum — the
   rotation orbit is measure-preserving.  This is why the near-pole sum depends
   only on the modulus `b′` and the gcd `e`, not on the speeds' particulars.
2. **The folding identity** (`lavSum_eq_folded`, n odd case): the κ-free sum
   folds to `1 + 2·Σ_{j=1}^{(n-1)/2} 1/j` — the exact harmonic form whose
   log-comparison is the paper lemma's `(2/e)(1+ln(b′/2e))`.

An instance anchor pins the adversarial modulus `n = 89` (the engineered triple
(97, 8633, 8536) of the cont.43 referee) by kernel computation.
-/

namespace LonelyRunner
namespace LRCCongruenceAveraging

open Finset

variable {n : ℕ} [NeZero n]

/-- The least absolute value of a residue: distance to 0 on the cycle. -/
def leastAbsVal (x : ZMod n) : ℕ := min x.val (n - x.val)

/-- The averaged reciprocal orbit sum (the near-pole line sum's arithmetic core). -/
def lavSum (n : ℕ) [NeZero n] : ℚ :=
  ∑ t : ZMod n, 1 / max 1 (leastAbsVal t : ℚ)

/-- **Orbit invariance**: multiplying by a unit permutes `ZMod n`, so the
    reciprocal least-abs sum is unchanged — the averaging lemma's engine. -/
theorem lavSum_mul_unit (κ : (ZMod n)ˣ) :
    ∑ t : ZMod n, 1 / max 1 (leastAbsVal ((κ : ZMod n) * t) : ℚ) = lavSum n := by
  unfold lavSum
  exact Fintype.sum_bijective (fun t => (κ : ZMod n) * t)
    (Units.mulLeft_bijective κ) _ _ (fun t => rfl)

/-- Negation preserves the least absolute value (the fold's pairing). -/
theorem leastAbsVal_neg (x : ZMod n) : leastAbsVal (-x) = leastAbsVal x := by
  unfold leastAbsVal
  rcases eq_or_ne x 0 with h | h
  · simp [h]
  · have hv : x.val ≠ 0 := by
      simpa using h
    have hlt : x.val < n := x.val_lt
    have hneg : (-x).val = n - x.val := by
      rw [ZMod.neg_val]
      simp [ZMod.val_eq_zero, h]
    rw [hneg]
    have h1 : n - (n - x.val) = x.val := by omega
    rw [h1]
    exact min_comm _ _

/-- The instance anchor at the adversarial modulus 89 (cont.43's engineered
    triple): `lavSum 89 = 1 + 2·H₄₄` exactly, kernel-checked. -/
theorem lavSum_89 :
    lavSum 89 = 1 + 2 * ∑ j ∈ Finset.range 44, (1 : ℚ) / (j + 1) := by
  decide +kernel

/-- Monotone comparison shape (the paper bound's skeleton): the 89-instance sum
    is under the crude linear bound `1 + 2·44` — sanity floor for the harmonic
    gain (the true value is ≈ 9.7 vs 89). -/
theorem lavSum_89_lt : lavSum 89 < 11 := by
  rw [lavSum_89]
  norm_num [Finset.sum_range_succ]


/-! ## The general folded identity (cont.45) -/

/-- The ℕ-level least-abs value. -/
def lavN (n v : ℕ) : ℕ := min v (n - v)

theorem leastAbsVal_eq_lavN (x : ZMod n) : leastAbsVal x = lavN n x.val := rfl

/-- Transfer to a `range n` sum over values (klein's Thm892 val-template). -/
theorem lavSum_eq_natSum (n : ℕ) [NeZero n] :
    lavSum n = ∑ v ∈ Finset.range n, 1 / max 1 (lavN n v : ℚ) := by
  unfold lavSum
  refine Finset.sum_nbij' (fun t : ZMod n => t.val) (fun v => (v : ZMod n))
    ?_ ?_ ?_ ?_ ?_
  · intro t _
    exact Finset.mem_range.mpr t.val_lt
  · intro v _
    exact Finset.mem_univ _
  · intro t _
    simp
  · intro v hv
    exact ZMod.val_natCast_of_lt (Finset.mem_range.mp hv)
  · intro t _
    rfl

/-- **THE GENERAL FOLDED IDENTITY** (odd modulus): the reciprocal least-abs sum
    over `ZMod n`, `n = 2m+1`, equals `1 + 2·H_m` exactly — the harmonic law the
    congruence-averaging lemma logs.  Subsumes `lavSum_89` (`m = 44`). -/
theorem lavSum_odd (m : ℕ) :
    lavSum (2 * m + 1) = 1 + 2 * ∑ j ∈ Finset.range m, (1 : ℚ) / (j + 1) := by
  have : NeZero (2 * m + 1) := ⟨by omega⟩
  rw [lavSum_eq_natSum]
  have hsplit : Finset.range (2 * m + 1) =
      Finset.range (m + 1) ∪ (Finset.range (2 * m + 1)).filter (fun v => m + 1 ≤ v) := by
    ext v
    simp only [Finset.mem_union, Finset.mem_range, Finset.mem_filter]
    omega
  -- direct route: split the sum at m+1 via sum_range_add-style reindexing
  have key : ∑ v ∈ Finset.range (2 * m + 1), 1 / max 1 (lavN (2 * m + 1) v : ℚ)
      = (∑ v ∈ Finset.range (m + 1), 1 / max 1 (lavN (2 * m + 1) v : ℚ))
        + ∑ v ∈ Finset.range m, 1 / max 1 (lavN (2 * m + 1) (m + 1 + v) : ℚ) := by
    rw [show 2 * m + 1 = (m + 1) + m from by omega]
    exact Finset.sum_range_add _ _ _
  rw [key]
  have low : ∑ v ∈ Finset.range (m + 1), 1 / max 1 (lavN (2 * m + 1) v : ℚ)
      = 1 + ∑ j ∈ Finset.range m, (1 : ℚ) / (j + 1) := by
    rw [Finset.sum_range_succ']
    have h0 : (1 : ℚ) / max 1 (lavN (2 * m + 1) 0 : ℚ) = 1 := by
      simp [lavN]
    rw [h0, add_comm]
    congr 1
    refine Finset.sum_congr rfl fun j hj => ?_
    have hj' : j + 1 ≤ m := Nat.succ_le_of_lt (Finset.mem_range.mp hj)
    have hlav : lavN (2 * m + 1) (j + 1) = j + 1 := by
      unfold lavN
      omega
    rw [hlav]
    have : max 1 ((j + 1 : ℕ) : ℚ) = ((j + 1 : ℕ) : ℚ) := by
      apply max_eq_right
      exact_mod_cast Nat.one_le_iff_ne_zero.mpr (Nat.succ_ne_zero j)
    rw [this]
    push_cast
    ring
  have high : ∑ v ∈ Finset.range m, 1 / max 1 (lavN (2 * m + 1) (m + 1 + v) : ℚ)
      = ∑ j ∈ Finset.range m, (1 : ℚ) / (j + 1) := by
    rw [← Finset.sum_range_reflect
      (fun v => 1 / max 1 (lavN (2 * m + 1) (m + 1 + v) : ℚ)) m]
    refine Finset.sum_congr rfl fun j hj => ?_
    have hj' : j < m := Finset.mem_range.mp hj
    have hlav : lavN (2 * m + 1) (m + 1 + (m - 1 - j)) = j + 1 := by
      unfold lavN
      omega
    rw [hlav]
    have hmax : max 1 ((j + 1 : ℕ) : ℚ) = ((j + 1 : ℕ) : ℚ) := by
      apply max_eq_right
      exact_mod_cast Nat.one_le_iff_ne_zero.mpr (Nat.succ_ne_zero j)
    rw [hmax]
    push_cast
    ring
  rw [low, high]
  ring


/-! ## Axiom audit -/

#print axioms lavSum_mul_unit
#print axioms leastAbsVal_neg
#print axioms lavSum_89
#print axioms lavSum_odd

end LRCCongruenceAveraging
end LonelyRunner
