import Mathlib

/-
# The algebraic kernel of the `H ≥ disc` SCC reduction (THM-1950)

klein-2026-07-21-S400.  death-star's WOWII candidate `H(T) ≥ disc(T)` (Rédei Hamiltonian-path
count ≥ the poly-time skew-determinant `disc = |det(I+K)|/2^{n-1}`, HYP-8636) reduces to the
strongly-connected case by the strong-component decomposition, using the total inverse-response
`s(T) = 𝟙ᵀ(I+K)⁻¹𝟙`.  The proved SCC-composition law is
  `disc(C₁ ⇒ C₂) = disc(C₁)·disc(C₂)·(1 + s₁s₂)/2`,  `s(C₁ ⇒ C₂) = (s₁+s₂)/(1+s₁s₂)`,
and the induction on the invariant `P(T) = max(1, s(T))·disc(T)` (peel the top strong component,
using `H(T)=H(C₁)H(T')`, THM-1860) closes iff the following **two-variable kernel inequality**
holds.  It is the algebraic heart of the reduction — the analogue of THM-1860's `sum ≤ product`
kernel (`SumLeProd.lean`), and the piece the WOWII "formalize" step asks for.  Proved here
`sorry`-free.

For `x = s(C₁) ≥ 0`, `y = s(T') ≥ 0`:  `max(1,x)·max(1,y) ≥ max(1+xy, x+y)/2`, and since
`s(C₁ ⇒ T')·(1+xy) = x+y` and `disc(C₁ ⇒ T') = disc(C₁)disc(T')(1+xy)/2`, this is exactly the
step `H(C₁)H(T') ≥ P(C₁)P(T') ≥ P(C₁ ⇒ T') ≥ disc(C₁ ⇒ T')`.
-/

namespace HgeDiscKernel

/-- **The kernel inequality.**  For `x, y ≥ 0`,
`max (1 + x*y) (x + y) / 2 ≤ max 1 x * max 1 y`.  This is the algebraic core of the `H ≥ disc`
SCC reduction (THM-1950): with `x = s(C₁)`, `y = s(T')` it powers the peel of the top strong
component. -/
theorem kernel_ineq {x y : ℝ} (hx : 0 ≤ x) (hy : 0 ≤ y) :
    max (1 + x * y) (x + y) / 2 ≤ max 1 x * max 1 y := by
  have h2 : max (1 + x * y) (x + y) ≤ 2 * (max 1 x * max 1 y) := by
    rcases le_total 1 x with hx1 | hx1 <;> rcases le_total 1 y with hy1 | hy1
    · -- 1 ≤ x, 1 ≤ y
      rw [max_eq_right hx1, max_eq_right hy1]
      apply max_le
      · nlinarith [mul_nonneg (sub_nonneg.2 hx1) (sub_nonneg.2 hy1)]
      · nlinarith [mul_nonneg hx (sub_nonneg.2 hy1), mul_nonneg hy (sub_nonneg.2 hx1)]
    · -- 1 ≤ x, y ≤ 1
      rw [max_eq_right hx1, max_eq_left hy1]
      apply max_le
      · nlinarith [sub_nonneg.2 hx1, mul_nonneg hx (sub_nonneg.2 hy1)]
      · nlinarith [sub_nonneg.2 hx1, sub_nonneg.2 hy1]
    · -- x ≤ 1, 1 ≤ y
      rw [max_eq_left hx1, max_eq_right hy1]
      apply max_le
      · nlinarith [sub_nonneg.2 hy1, mul_nonneg hy (sub_nonneg.2 hx1)]
      · nlinarith [sub_nonneg.2 hy1, sub_nonneg.2 hx1]
    · -- x ≤ 1, y ≤ 1
      rw [max_eq_left hx1, max_eq_left hy1]
      apply max_le
      · nlinarith [sub_nonneg.2 hx1, mul_nonneg hx (sub_nonneg.2 hy1)]
      · nlinarith [sub_nonneg.2 hx1, sub_nonneg.2 hy1]
  linarith

/-- The peel step in the exact form used by the reduction.  With `dc, dt ≥ 0` the component
discs, `x = s(C₁) ≥ 0`, `y = s(T') ≥ 0`, and the strong/inductive bases
`H₁ ≥ max(1,x)·dc`, `Ht ≥ max(1,y)·dt`, the product `H₁·Ht` dominates the composite
`P = max(1, (x+y)/(1+xy))·(dc·dt·(1+xy)/2)`.  Stated with the composite already simplified
via `max(1, s)·(1+xy) = max(1+xy, x+y)` (since `s(1+xy)=x+y`). -/
theorem peel_step {dc dt x y H1 Ht : ℝ}
    (hdc : 0 ≤ dc) (hdt : 0 ≤ dt) (hx : 0 ≤ x) (hy : 0 ≤ y)
    (hH1 : max 1 x * dc ≤ H1) (hHt : max 1 y * dt ≤ Ht) :
    max (1 + x * y) (x + y) / 2 * (dc * dt) ≤ H1 * Ht := by
  have hk := kernel_ineq hx hy
  have hbase : (0:ℝ) ≤ max 1 x * dc := mul_nonneg (le_trans zero_le_one (le_max_left _ _)) hdc
  have hbaset : (0:ℝ) ≤ max 1 y * dt := mul_nonneg (le_trans zero_le_one (le_max_left _ _)) hdt
  have step : max (1 + x * y) (x + y) / 2 * (dc * dt)
      ≤ (max 1 x * dc) * (max 1 y * dt) := by
    have : max (1 + x * y) (x + y) / 2 * (dc * dt)
        ≤ (max 1 x * max 1 y) * (dc * dt) := by
      apply mul_le_mul_of_nonneg_right hk (mul_nonneg hdc hdt)
    calc max (1 + x * y) (x + y) / 2 * (dc * dt)
        ≤ (max 1 x * max 1 y) * (dc * dt) := this
      _ = (max 1 x * dc) * (max 1 y * dt) := by ring
  calc max (1 + x * y) (x + y) / 2 * (dc * dt)
      ≤ (max 1 x * dc) * (max 1 y * dt) := step
    _ ≤ H1 * Ht := by
        apply mul_le_mul hH1 hHt hbaset (le_trans hbase hH1)

#print axioms kernel_ineq
#print axioms peel_step

end HgeDiscKernel
