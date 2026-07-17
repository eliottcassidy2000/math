/- Thm892Jordan.lean — klein-2026-07-16-S318.
   THM-892 (C*), the Möbius/Jordan layer (boxeph manifest item 6, second lemma):
   with J₂(q) = Σ_{d∣q} μ(q/d)·d² (the second Jordan totient), the divisor sums
   close:  Σ_{d∣q} J₂(d) = q²  for every q ≠ 0.
   This is the arithmetic driving THM-892's Σ_{u∈(Z/q)*} csc²(πu/q) = J₂(q)/3:
   Möbius inversion of MI0 over the divisor lattice is exactly this identity.
   Proof: pure Dirichlet algebra — J₂ = μ * pow2, so
   ζ * J₂ = (ζ * μ) * pow2 = 1 * pow2 = pow2.  No sorries, no native_decide. -/
import Mathlib.NumberTheory.ArithmeticFunction
import Mathlib.Tactic

namespace LonelyRunner
namespace LRC14
namespace Thm892

/-- The second Jordan totient as a Möbius-weighted divisor sum. -/
def J2 (q : ℕ) : ℤ :=
  ∑ d ∈ q.divisors, (ArithmeticFunction.moebius (q / d)) * (d : ℤ) ^ 2

/-- `J₂` is the Dirichlet convolution `μ * pow 2`, pointwise. -/
theorem J2_eq_conv (m : ℕ) :
    (ArithmeticFunction.moebius *
      ((ArithmeticFunction.pow 2 : ArithmeticFunction ℕ) : ArithmeticFunction ℤ)) m
      = J2 m := by
  rw [ArithmeticFunction.mul_apply]
  unfold J2
  rw [← Nat.sum_divisorsAntidiagonal'
      (fun a b => (ArithmeticFunction.moebius a : ℤ) * (b : ℤ) ^ 2)]
  apply Finset.sum_congr rfl
  intro x hx
  rw [Nat.mem_divisorsAntidiagonal] at hx
  have hx2 : x.2 ≠ 0 := by
    intro h
    exact hx.2 (by rw [← hx.1, h, mul_zero])
  rw [ArithmeticFunction.natCoe_apply, ArithmeticFunction.pow_apply,
      if_neg (by exact fun hc => hx2 hc.2)]
  push_cast
  ring

/-- **THM-892 (C*), the Möbius/Jordan layer**: `Σ_{d∣q} J₂(d) = q²` for `q ≠ 0`. -/
theorem J2_divisor_sum (q : ℕ) (hq : q ≠ 0) :
    ∑ d ∈ q.divisors, J2 d = (q : ℤ) ^ 2 := by
  have hpow : ((ArithmeticFunction.pow 2 : ArithmeticFunction ℕ) :
      ArithmeticFunction ℤ) q = (q : ℤ) ^ 2 := by
    rw [ArithmeticFunction.natCoe_apply, ArithmeticFunction.pow_apply,
        if_neg (by exact fun hc => hq hc.2)]
    push_cast
    ring
  calc ∑ d ∈ q.divisors, J2 d
      = ∑ d ∈ q.divisors, (ArithmeticFunction.moebius *
          ((ArithmeticFunction.pow 2 : ArithmeticFunction ℕ) :
            ArithmeticFunction ℤ)) d :=
        Finset.sum_congr rfl (fun d _ => (J2_eq_conv d).symm)
    _ = (((ArithmeticFunction.zeta : ArithmeticFunction ℕ) : ArithmeticFunction ℤ) *
          (ArithmeticFunction.moebius *
            ((ArithmeticFunction.pow 2 : ArithmeticFunction ℕ) :
              ArithmeticFunction ℤ))) q := by
        rw [ArithmeticFunction.coe_zeta_mul_apply]
    _ = ((((ArithmeticFunction.zeta : ArithmeticFunction ℕ) : ArithmeticFunction ℤ) *
          ArithmeticFunction.moebius) *
            ((ArithmeticFunction.pow 2 : ArithmeticFunction ℕ) :
              ArithmeticFunction ℤ)) q := by
        rw [mul_assoc]
    _ = ((1 : ArithmeticFunction ℤ) *
          ((ArithmeticFunction.pow 2 : ArithmeticFunction ℕ) :
            ArithmeticFunction ℤ)) q := by
        rw [ArithmeticFunction.coe_zeta_mul_moebius]
    _ = ((ArithmeticFunction.pow 2 : ArithmeticFunction ℕ) :
          ArithmeticFunction ℤ) q := by
        rw [one_mul]
    _ = (q : ℤ) ^ 2 := hpow

/- Concrete anchors (J₂(7) = 48, J₂(13) = 168, J₂(14) = 144) are left to the
   paper side: kernel `decide` cannot reduce `moebius` (its `minFac` recursion
   is well-founded, the known decide-blocker), and the identity above is the
   manifest item. -/

end Thm892
end LRC14
end LonelyRunner
