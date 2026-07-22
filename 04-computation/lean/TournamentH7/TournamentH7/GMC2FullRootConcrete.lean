import Mathlib
import TournamentH7.GMC2LaurentShiftCheckA
import TournamentH7.GMC2RootPacketConcrete

/-!
# The concrete full-root Lagrange sum vanishes (`hfull`)

Discharges the `hfull` premise of `GMC2RootPacketConcrete.root_packet_eq_zero`: for `Φ` irreducible
over a characteristic-zero field `F` and splitting in `L`, with `0 < M < deg Φ`,
`∑_{α ∈ Φ.rootSet L} α^{M-1} / Φ'(α) = 0`.

Proof: `Φ` is separable (irreducible over a perfect field), so its roots in `L` are distinct;
`Φ.map = C(lc)·nodal(roots)` and hence `Φ'(α) = lc·nodal'(α)`, reducing the sum to
`lc⁻¹ · (codex's `sum_pow_pred_div_derivative_nodal_eq_zero`) = 0`.
-/

open Polynomial

namespace GMC2FullRootConcrete

variable {F L : Type*} [Field F] [CharZero F] [Field L] [Algebra F L]

theorem fullRootSum_eq_zero (Φ : F[X]) (hΦ : Irreducible Φ)
    (hsplit : Polynomial.Splits (Φ.map (algebraMap F L))) (M : ℕ) (hM0 : 0 < M)
    (hMdeg : M < Φ.natDegree) :
    ∑ α : Φ.rootSet L, ((α : L) ^ (M - 1) / aeval (α : L) (derivative Φ)) = 0 := by
  classical
  set φ := algebraMap F L
  set Ψ : L[X] := Φ.map φ with hΨ
  have hΦ0 : Φ ≠ 0 := hΦ.ne_zero
  have hΨ0 : Ψ ≠ 0 := by rw [hΨ]; exact (Polynomial.map_ne_zero_iff (algebraMap F L).injective).mpr hΦ0
  -- separability ⇒ distinct roots
  have hsep : Φ.Separable := (PerfectField.ofCharZero (K := F)).separable_of_irreducible hΦ
  have hΨsep : Ψ.Separable := hsep.map
  have hnodup : Ψ.roots.Nodup := Polynomial.nodup_roots hΨsep
  have hΨsplit : Polynomial.Splits Ψ := hsplit
  set s : Finset L := Ψ.roots.toFinset with hs
  -- s.card = deg Φ
  have hcard : s.card = Φ.natDegree := by
    rw [hs, Multiset.toFinset_card_of_nodup hnodup,
      Polynomial.splits_iff_card_roots.mp hΨsplit, hΨ, Polynomial.natDegree_map]
  have hMcard : M < s.card := by rw [hcard]; exact hMdeg
  -- rootSet = ↑s
  have hrootSet : Φ.rootSet L = ↑s := by
    rw [Polynomial.rootSet_def, hs, hΨ, Polynomial.aroots_def]
  -- Ψ = C(lc) * nodal s id
  have hsval : Ψ.roots = s.val := by rw [hs, Multiset.toFinset_val, Multiset.dedup_eq_self.mpr hnodup]
  have hprod : (Ψ.roots.map (fun a => X - C a)).prod = Lagrange.nodal s id := by
    rw [hsval]; simp [Lagrange.nodal]
  have hlc0 : Ψ.leadingCoeff ≠ 0 := leadingCoeff_ne_zero.mpr hΨ0
  have hΨeq : Ψ = C Ψ.leadingCoeff * Lagrange.nodal s id := by
    rw [← hprod]; exact hΨsplit.eq_prod_roots
  -- Φ'(α) mapped = lc * nodal'(α)
  have hderiv : ∀ α : L, aeval α (derivative Φ) = Ψ.leadingCoeff *
      (derivative (Lagrange.nodal s id)).eval α := by
    intro α
    have hΦΨ : aeval α (derivative Φ) = (derivative Ψ).eval α := by
      rw [aeval_def, eval₂_eq_eval_map, ← derivative_map, ← hΨ]
    rw [hΦΨ]
    conv_lhs => rw [hΨeq]
    rw [derivative_mul, derivative_C, zero_mul, zero_add, eval_mul, eval_C]
  -- assemble: convert the subtype sum to a `Finset` sum over `s`, rewrite the denominator, factor `lc⁻¹`
  have hconv : (∑ α : Φ.rootSet L, ((α : L) ^ (M - 1) / aeval (α : L) (derivative Φ)))
      = ∑ α ∈ s, (α ^ (M - 1) /
          (Ψ.leadingCoeff * (derivative (Lagrange.nodal s id)).eval α)) := by
    rw [← Finset.sum_coe_sort s]
    refine Fintype.sum_equiv (Equiv.setCongr hrootSet) _ _ (fun α => ?_)
    simp only [Equiv.setCongr_apply]
    rw [hderiv]
  rw [hconv]
  have hfactor : ∀ α : L, α ^ (M - 1) /
      (Ψ.leadingCoeff * (derivative (Lagrange.nodal s id)).eval α)
      = (α ^ (M - 1) / (derivative (Lagrange.nodal s id)).eval α) / Ψ.leadingCoeff :=
    fun α => by rw [mul_comm, div_mul_eq_div_div]
  simp_rw [hfactor, ← Finset.sum_div]
  rw [GMC2FullRootLagrange.sum_pow_pred_div_derivative_nodal_eq_zero s M hM0 hMcard, zero_div]

/-- **The additive route reduces to the b=1 wrapper.**  Combining the now-proved full-root Lagrange
sum (`fullRootSum_eq_zero`, discharging `hfull`) with the root-packet lemma
(`GMC2RootPacketConcrete.root_packet_eq_zero`): for `Φ` irreducible over a char-0 field `F` splitting
in `L`, if the small-cluster residue sum equals `1` (`hb`, the b=1 wrapper — the one remaining input),
then `False`.  So the additive one-variable DvdK contradiction is complete **modulo** the b=1 wrapper,
with **no** THM-1550 / small-root product / Hensel / Wiener–Hopf. -/
theorem additive_dvdk_reduces_to_smallSum [FiniteDimensional F L] [Normal F L]
    (Φ : F[X]) (hΦ : Irreducible Φ) (hsplit : Polynomial.Splits (Φ.map (algebraMap F L)))
    (M : ℕ) (hM0 : 0 < M) (hMdeg : M < Φ.natDegree)
    (x0 : Φ.rootSet L) (S : Finset (Φ.rootSet L))
    (hb : algebraMap F L 1
        = ∑ β ∈ S, ((β : L) ^ (M - 1) / aeval (β : L) (derivative Φ))) :
    False := by
  have h1 : algebraMap F L (1 : F) = 0 :=
    GMC2RootPacketConcrete.root_packet_eq_zero Φ hΦ (M - 1) x0 S
      (fullRootSum_eq_zero Φ hΦ hsplit M hM0 hMdeg) 1 hb
  rw [map_one] at h1
  exact one_ne_zero h1

end GMC2FullRootConcrete

#print axioms GMC2FullRootConcrete.fullRootSum_eq_zero
#print axioms GMC2FullRootConcrete.additive_dvdk_reduces_to_smallSum
