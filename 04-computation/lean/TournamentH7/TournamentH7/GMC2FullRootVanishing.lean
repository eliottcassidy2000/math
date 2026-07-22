import TournamentH7.GMC2LaurentShiftCheckA

/-!
# The concrete full-root Lagrange vanishing — the `hfull` input of the additive route

`GMC2RootPacketConcrete.root_packet_eq_zero` (boxeph — the additive analog of the orbit-product core)
takes as a hypothesis `hfull`, that the **full-root** derivative-weighted sum vanishes:
`∑_{α ∈ rootSet Φ} αᵏ/Φ'(α) = 0` for `k ≤ deg Φ − 2`.  This is the residue-at-infinity / Lagrange
identity, and it is the additive analog of Vieta.

This module discharges `hfull` concretely.  `GMC2FullRootLagrange.sum_pow_div_derivative_nodal_eq_zero`
(codex) proves the identity for a `Finset` of distinct roots against the **monic** nodal polynomial.
The bridge to `Φ`'s actual (non-monic) derivative over its splitting field is:

* the product over the `rootSet` subtype becomes a `Finset` sum over the root set (separable ⇒ nodup,
  `Finset.sum_attach`) — the same reindexing used for the multiplicative-route Vieta
  (`GMC2PhiVieta.prod_rootSet_Phi`);
* `Φ.map = C(lc)·nodal` (`Splits.eq_prod_roots` + the nodup multiset→Finset product collapse), so
  `Φ'(α) = lc · nodal'(α)` and the leading coefficient `lc` factors out of the whole sum.

So the additive route's `hfull` is a theorem — leaving only `hA` (the small-cluster packet sum equals
the base-field value) as the remaining input, the additive analog of THM-1550's `hS`.
-/

open Polynomial
open scoped BigOperators

namespace GMC2FullRootVanishing

variable {F L : Type*} [Field F] [Field L] [Algebra F L]

/-- **The full-root Lagrange sum vanishes** (`hfull` of `root_packet_eq_zero`).  For `Φ` splitting
separably in `L` and `k + 1 < deg Φ`, the derivative-weighted sum over the *distinct* roots is `0`.
The leading coefficient of `Φ.map` factors out, reducing to codex's monic nodal identity. -/
theorem full_root_lagrange_sum_eq_zero
    (Φ : F[X]) (k : ℕ)
    (hsplit : (Φ.map (algebraMap F L)).Splits)
    (hsep : Separable (Φ.map (algebraMap F L)))
    (hk : k + 1 < Φ.natDegree) :
    ∑ α : Φ.rootSet L, ((α : L) ^ k / aeval (α : L) (derivative Φ)) = 0 := by
  classical
  set Ψ := Φ.map (algebraMap F L) with hΨ
  have hnd : Ψ.roots.Nodup := nodup_roots hsep
  set s := Ψ.roots.toFinset with hs
  -- `Φ.map = C(lc) · nodal(roots)` (separable ⇒ nodup collapses the multiset product to the nodal one)
  have hfact : Ψ = C Ψ.leadingCoeff * Lagrange.nodal s id := by
    have key : Ψ = C Ψ.leadingCoeff * (Ψ.roots.map (fun a => X - C a)).prod :=
      hsplit.eq_prod_roots
    have hnodal : Lagrange.nodal s id = (Ψ.roots.map (fun a => X - C a)).prod := by
      rw [Lagrange.nodal, Finset.prod_eq_multiset_prod, hs, Multiset.toFinset_val, hnd.dedup]
      simp only [id_eq]
    rw [hnodal]; exact key
  have hcard : s.card = Φ.natDegree := by
    have h1 : Ψ.roots.card = Ψ.natDegree := Polynomial.splits_iff_card_roots.mp hsplit
    have h2 : Ψ.natDegree = Φ.natDegree := by
      simp [hΨ, Polynomial.natDegree_map_eq_of_injective (algebraMap F L).injective]
    rw [hs, Multiset.toFinset_card_of_nodup hnd, h1, h2]
  -- termwise: `Φ'(x) = lc · nodal'(x)`
  have hterm : ∀ x ∈ s, aeval x (derivative Φ)
      = Ψ.leadingCoeff * (derivative (Lagrange.nodal s id)).eval x := by
    intro x hx
    have h1 : aeval x (derivative Φ) = eval x (derivative Ψ) := by
      rw [Polynomial.aeval_def, Polynomial.eval₂_eq_eval_map, hΨ, Polynomial.derivative_map]
    rw [h1]
    conv_lhs => rw [hfact, Polynomial.derivative_C_mul, Polynomial.eval_mul, Polynomial.eval_C]
  -- reindex the subtype sum to the root `Finset`
  have hsum1 : (∑ α : Φ.rootSet L, ((α : L) ^ k / aeval (α : L) (derivative Φ)))
      = ∑ x ∈ s, (x ^ k / aeval x (derivative Φ)) := by
    have hbridge : (∑ α : Φ.rootSet L, ((α : L) ^ k / aeval (α : L) (derivative Φ)))
        = ∑ x ∈ (Φ.aroots L).toFinset, (x ^ k / aeval x (derivative Φ)) := by
      first
      | exact Finset.sum_attach _ _
      | exact Finset.sum_attach _ (fun x => x ^ k / aeval x (derivative Φ))
    rw [hbridge]
  rw [hsum1]
  -- pull out `lc⁻¹` and apply codex's monic Lagrange identity
  have hstep : ∀ x ∈ s, x ^ k / aeval x (derivative Φ)
      = (Ψ.leadingCoeff)⁻¹ * (x ^ k / (derivative (Lagrange.nodal s id)).eval x) := by
    intro x hx
    rw [hterm x hx, div_eq_mul_inv, div_eq_mul_inv, mul_inv]; ring
  have hz : (∑ x ∈ s, x ^ k / (derivative (Lagrange.nodal s id)).eval x) = 0 :=
    GMC2FullRootLagrange.sum_pow_div_derivative_nodal_eq_zero s k (by rw [hcard]; exact hk)
  rw [Finset.sum_congr rfl hstep, ← Finset.mul_sum, hz, mul_zero]

end GMC2FullRootVanishing

#print axioms GMC2FullRootVanishing.full_root_lagrange_sum_eq_zero
