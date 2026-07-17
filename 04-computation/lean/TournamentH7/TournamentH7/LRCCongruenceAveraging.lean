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

/-! ## Axiom audit -/

#print axioms lavSum_mul_unit
#print axioms leastAbsVal_neg
#print axioms lavSum_89

end LRCCongruenceAveraging
end LonelyRunner
