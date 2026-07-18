import Mathlib

/-!
# The grid-count core (kind-pasteur-S128c48)

The weight-1 sampling bridge's kernel piece (THM-984): integers strictly
between two rationals number at least the gap minus one.  Applied at modulus
`q` to a good interval `(α, β)`, the integers `p ∈ (qα, qβ)` give good grid
points `p/q` with count ≥ `q(β−α) − 1` — the per-interval live floor.  The
union assembly over a packet's disjoint good intervals sums these floors
(disjoint intervals give disjoint integer ranges), yielding
`liveCount ≥ q·μ₀ − (#intervals)`, the explicit-modulus engine of the
live-route certificate.
-/

namespace LonelyRunner
namespace LRCGridCount

/-- Every integer in `Ioo ⌊x⌋ ⌈y⌉` lies strictly between `x` and `y`. -/
theorem mem_Ioo_strict (x y : ℚ) (p : ℤ) (hp : p ∈ Finset.Ioo ⌊x⌋ ⌈y⌉) :
    x < (p : ℚ) ∧ (p : ℚ) < y := by
  obtain ⟨h1, h2⟩ := Finset.mem_Ioo.mp hp
  constructor
  · calc x < (⌊x⌋ : ℚ) + 1 := Int.lt_floor_add_one x
      _ ≤ (p : ℚ) := by exact_mod_cast Int.add_one_le_iff.mpr h1
  · have hpInt : p ≤ (⌈y⌉ : ℤ) - 1 := by omega
    have hpRat : (p : ℚ) ≤ (⌈y⌉ : ℚ) - 1 := by exact_mod_cast hpInt
    have hceil : (⌈y⌉ : ℚ) < y + 1 := Int.ceil_lt_add_one y
    linarith

/-- **The strict-count floor**: integers strictly between `x` and `y` number
    at least `y − x − 1`. -/
theorem card_Ioo_floor_ceil (x y : ℚ) :
    y - x - 1 ≤ ((Finset.Ioo ⌊x⌋ ⌈y⌉).card : ℚ) := by
  by_cases h : (⌈y⌉ : ℤ) ≤ ⌊x⌋
  · have hcard : (Finset.Ioo ⌊x⌋ ⌈y⌉).card = 0 := by
      rw [Finset.card_eq_zero]
      exact Finset.Ioo_eq_empty (not_lt.mpr h)
    rw [hcard]
    have hy : (⌈y⌉ : ℚ) ≤ (⌊x⌋ : ℚ) := by exact_mod_cast h
    have h1 : y ≤ (⌈y⌉ : ℚ) := Int.le_ceil y
    have h2 : (⌊x⌋ : ℚ) ≤ x := Int.floor_le x
    push_cast
    linarith
  · have hlt : (⌊x⌋ : ℤ) < ⌈y⌉ := by omega
    have hcard : ((Finset.Ioo ⌊x⌋ ⌈y⌉).card : ℤ) = ⌈y⌉ - ⌊x⌋ - 1 := by
      rw [Int.card_Ioo]
      exact Int.toNat_of_nonneg (by omega)
    have hq : ((Finset.Ioo ⌊x⌋ ⌈y⌉).card : ℚ) = (⌈y⌉ : ℚ) - (⌊x⌋ : ℚ) - 1 := by
      exact_mod_cast hcard
    rw [hq]
    have h1 : y ≤ (⌈y⌉ : ℚ) := Int.le_ceil y
    have h2 : (⌊x⌋ : ℚ) ≤ x := Int.floor_le x
    linarith

/-- **The per-interval live floor** (the sampling bridge's application shape):
    at modulus `q > 0`, the good interval `(α, β)` contains at least
    `q(β−α) − 1` grid points `p/q` — each certified strictly interior. -/
theorem interval_live_floor (α β : ℚ) (q : ℕ) (hq : 0 < q) :
    (q : ℚ) * (β - α) - 1 ≤ ((Finset.Ioo ⌊(q : ℚ) * α⌋ ⌈(q : ℚ) * β⌉).card : ℚ) ∧
    ∀ p ∈ Finset.Ioo ⌊(q : ℚ) * α⌋ ⌈(q : ℚ) * β⌉,
      α < (p : ℚ) / q ∧ (p : ℚ) / q < β := by
  constructor
  · have := card_Ioo_floor_ceil ((q : ℚ) * α) ((q : ℚ) * β)
    linarith [this, mul_sub (q : ℚ) β α]
  · intro p hp
    obtain ⟨h1, h2⟩ := mem_Ioo_strict ((q : ℚ) * α) ((q : ℚ) * β) p hp
    have hq' : (0 : ℚ) < q := by exact_mod_cast hq
    constructor
    · rw [lt_div_iff₀ hq']
      linarith [h1]
    · rw [div_lt_iff₀ hq']
      linarith [h2]

/-! ## Axiom audit -/

#print axioms card_Ioo_floor_ceil
#print axioms interval_live_floor

end LRCGridCount
end LonelyRunner
