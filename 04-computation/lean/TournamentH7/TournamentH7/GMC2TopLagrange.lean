import Mathlib

/-!
# The top Lagrange identity `∑ α^{|s|-1} / nodal'(α) = 1`

codex's `GMC2LaurentShiftCheckA.GMC2FullRootLagrange.sum_pow_div_derivative_nodal_eq_zero` proves the
*vanishing* half of the Lagrange/residue family: for a finite set `s` of distinct points and `k ≤ |s|-2`,
`∑_{α∈s} α^k / nodal'(α) = 0` (the residue-at-infinity identity behind `hfull`).  This module supplies
the **complementary top case `k = |s|-1`, whose value is `1`** — so the full family reads

  `∑_{α∈s} α^k / nodal'(α) = δ_{k, |s|-1}`   for `0 ≤ k ≤ |s|-1`.

Both halves are one appeal to `Lagrange.coeff_eq_sum`: the sum equals the `(|s|-1)`-coefficient of the
interpolant of `X^k` through `s`, which is `X^k`'s own `(|s|-1)`-coefficient (degree `< |s|`), i.e.
`0` for `k < |s|-1` and `1` for `k = |s|-1`.

**Why this is the right elementary building block.**  For `Φ = X^M − t·R` over `F(t)`, the small-root
factor is the monic degree-`M` Weierstrass distinguished polynomial `P` (mac-mini's
`GMC2DvdKWeierstrass.smallRootFactor`).  The DvdK "b = 1 wrapper" is the residue identity
`∑_{β root of P} β^{M-1}/Φ'(β) = F(t)` (`= 1` under constant-term vanishing).  Since `Φ'(β) = P'(β)·h(β)`
at a root `β` of `P`, the `h ≡ 1` (leading / `D_0`) part of that sum is exactly
`∑_β β^{M-1}/P'(β) = ∑_β β^{|roots|-1}/nodal'(β) = 1` — this lemma.  It is *not* the deep
`[x^0]`/annulus identity (the shared valuation crux, owned by the Weierstrass lane); it is the
uncontested classical value that identity sits on top of.  General field, kernel-pure.
-/

open Polynomial
open scoped BigOperators

namespace GMC2TopLagrange

/-- **The top Lagrange identity.**  For a finite set `s` of distinct field elements,
`∑_{α∈s} α^{|s|-1} / (nodal s)'(α) = 1`.  The `k = |s|-1` companion of codex's
`sum_pow_div_derivative_nodal_eq_zero` (`k ≤ |s|-2 ⟹ 0`); together they give the Lagrange family
`∑ α^k/nodal'(α) = δ_{k,|s|-1}`.  Proof: `Lagrange.coeff_eq_sum` equates the sum with the top
(`|s|-1`) coefficient of `X^{|s|-1}`, which is `1`. -/
theorem sum_pow_pred_div_derivative_nodal_eq_one
    {K : Type*} [Field K] (s : Finset K) (hs : 0 < s.card) :
    ∑ α ∈ s, α ^ (s.card - 1) /
        (Polynomial.derivative (Lagrange.nodal s id)).eval α = 1 := by
  classical
  have hdeg : ((Polynomial.X : K[X]) ^ (s.card - 1)).degree < s.card := by
    rw [Polynomial.degree_X_pow]
    exact_mod_cast (by omega : s.card - 1 < s.card)
  have hlag := Lagrange.coeff_eq_sum (s := s) (v := id)
    (P := (Polynomial.X : K[X]) ^ (s.card - 1)) Function.injective_id.injOn hdeg
  have hcoeff : ((Polynomial.X : K[X]) ^ (s.card - 1)).coeff (s.card - 1) = 1 := by
    rw [Polynomial.coeff_X_pow]; simp
  rw [hcoeff] at hlag
  calc
    ∑ α ∈ s, α ^ (s.card - 1) /
        (Polynomial.derivative (Lagrange.nodal s id)).eval α =
        ∑ α ∈ s, ((Polynomial.X : K[X]) ^ (s.card - 1)).eval α /
          ∏ β ∈ s.erase α, (α - β) := by
      apply Finset.sum_congr rfl
      intro α hα
      have hderiv := Lagrange.eval_nodal_derivative_eval_node_eq
        (s := s) (v := id) hα
      simp only [id_eq] at hderiv
      rw [hderiv, Lagrange.eval_nodal]
      simp
    _ = 1 := hlag.symm

end GMC2TopLagrange

#print axioms GMC2TopLagrange.sum_pow_pred_div_derivative_nodal_eq_one
