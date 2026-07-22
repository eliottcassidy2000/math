import TournamentH7.GMC2LaurentShiftCheckA
import TournamentH7.GMC2RootPacketConcrete

/-!
# The full-root Lagrange sum for `Φ` (non-monic) — discharging `hfull` via the Weierstrass product

`GMC2RootPacketConcrete.root_packet_eq_zero` takes the full-root sum `∑_{α∈roots} αᵏ/Φ'(α) = 0` as a
hypothesis (`hfull`).  codex's `GMC2FullRootLagrange.sum_pow_div_derivative_nodal_eq_zero` proves it
for the *monic nodal* polynomial.  This module bridges to the actual (non-monic) `Φ` via the
**Weierstrass product form** `Φ.map = C(leadingCoeff)·∏(X − root)`: differentiating gives
`Φ'(α) = leadingCoeff · nodal'(α)` at a root, so the `Φ`-sum is `(1/leadingCoeff)` times the nodal sum
`= 0`.
-/

open Polynomial

namespace GMC2FullRootPhi

variable {F L : Type*} [Field F] [Field L] [Algebra F L] [DecidableEq L]

/-- **Weierstrass product form.**  A polynomial that splits with distinct roots is its leading
coefficient times the nodal polynomial of its roots. -/
theorem phi_map_eq (Φ : F[X]) (hsplit : Splits (Φ.map (algebraMap F L)))
    (hnd : (Φ.aroots L).Nodup) :
    Φ.map (algebraMap F L)
      = C (Φ.map (algebraMap F L)).leadingCoeff * Lagrange.nodal (Φ.aroots L).toFinset id := by
  conv_lhs => rw [hsplit.eq_prod_roots]
  congr 1
  rw [Lagrange.nodal, Finset.prod_eq_multiset_prod, Multiset.toFinset_val,
    Multiset.dedup_eq_self.mpr hnd, aroots_def]
  simp only [id_eq]

/-- **The derivative-at-root proportionality (Weierstrass).**  `Φ'(α) = leadingCoeff · nodal'(α)`. -/
theorem aeval_deriv_eq (Φ : F[X]) (hsplit : Splits (Φ.map (algebraMap F L)))
    (hnd : (Φ.aroots L).Nodup) (α : L) :
    aeval α (derivative Φ)
      = (Φ.map (algebraMap F L)).leadingCoeff *
          (derivative (Lagrange.nodal (Φ.aroots L).toFinset id)).eval α := by
  rw [Polynomial.aeval_def, Polynomial.eval₂_eq_eval_map, ← Polynomial.derivative_map]
  conv_lhs => rw [phi_map_eq Φ hsplit hnd]
  rw [Polynomial.derivative_C_mul, Polynomial.eval_mul, Polynomial.eval_C]

/-- **The full-root Lagrange sum for the non-monic `Φ`.**  `∑_{α∈roots} αᵏ/Φ'(α) = 0` for
`k + 1 < deg Φ`.  This discharges `hfull` in `root_packet_eq_zero`. -/
theorem full_root_sum_eq_zero (Φ : F[X]) (hsplit : Splits (Φ.map (algebraMap F L)))
    (hnd : (Φ.aroots L).Nodup) (hsep : Φ.Separable) (k : ℕ) (hk : k + 1 < Φ.natDegree) :
    ∑ α : Φ.rootSet L, ((α : L) ^ k / aeval (α : L) (derivative Φ)) = 0 := by
  set lc := (Φ.map (algebraMap F L)).leadingCoeff with hlc
  have hΦ0 : Φ ≠ 0 := by
    intro h; rw [h] at hk; simp at hk
  -- reindex the subtype sum to the roots Finset
  have hmem : ∀ x : L, x ∈ (Φ.aroots L).toFinset ↔ x ∈ Φ.rootSet L := by
    intro x
    simp only [Multiset.mem_toFinset, Polynomial.mem_aroots, Polynomial.mem_rootSet,
      Polynomial.map_ne_zero_iff (FaithfulSMul.algebraMap_injective F L)]
  rw [← Finset.sum_subtype (Φ.aroots L).toFinset hmem
    (fun α => α ^ k / aeval α (derivative Φ))]
  -- replace `Φ'(α)` by `lc · nodal'(α)` and factor out `1/lc`
  have hstep : ∀ α ∈ (Φ.aroots L).toFinset,
      (α ^ k / aeval α (derivative Φ))
        = (α ^ k / (derivative (Lagrange.nodal (Φ.aroots L).toFinset id)).eval α) / lc := by
    intro α _
    rw [aeval_deriv_eq Φ hsplit hnd α, hlc]
    rw [div_div, mul_comm]
  rw [Finset.sum_congr rfl hstep, ← Finset.sum_div]
  -- the nodal sum vanishes (codex's Lagrange), so the whole thing is `0 / lc = 0`
  have hcard : (Φ.aroots L).toFinset.card = Φ.natDegree := by
    rw [Multiset.toFinset_card_of_nodup hnd, aroots_def,
      Polynomial.splits_iff_card_roots.mp hsplit,
      Polynomial.natDegree_map_eq_of_injective (FaithfulSMul.algebraMap_injective F L)]
  rw [GMC2FullRootLagrange.sum_pow_div_derivative_nodal_eq_zero (Φ.aroots L).toFinset k
    (by rw [hcard]; exact hk), zero_div]

/-- **The self-contained additive root-packet lemma** (`hfull` discharged).  For `Φ` irreducible over
a char-0 field, splitting with distinct roots in `L`, the barycentric packet sum `b_k(S) ∈ F ⟹ = 0` —
no `hfull` hypothesis, no THM-1550, no small-root product/Hensel/Vieta.  This is the additive-route
core made self-contained, ready for the `b = 1` wrapper. -/
theorem root_packet_eq_zero_selfcontained [CharZero F] [FiniteDimensional F L] [Normal F L]
    (Φ : F[X]) (hΦ : Irreducible Φ) (hsplit : Splits (Φ.map (algebraMap F L)))
    (hnd : (Φ.aroots L).Nodup) (hsep : Φ.Separable) (k : ℕ) (hk : k + 1 < Φ.natDegree)
    (x0 : Φ.rootSet L) (S : Finset (Φ.rootSet L))
    (b0 : F) (hb : algebraMap F L b0 = ∑ β ∈ S, ((β : L) ^ k / aeval (β : L) (derivative Φ))) :
    algebraMap F L b0 = 0 :=
  GMC2RootPacketConcrete.root_packet_eq_zero Φ hΦ k x0 S
    (full_root_sum_eq_zero Φ hsplit hnd hsep k hk) b0 hb

end GMC2FullRootPhi

#print axioms GMC2FullRootPhi.phi_map_eq
#print axioms GMC2FullRootPhi.aeval_deriv_eq
#print axioms GMC2FullRootPhi.full_root_sum_eq_zero
#print axioms GMC2FullRootPhi.root_packet_eq_zero_selfcontained
