import TournamentH7.GMC2GalRootAction
import TournamentH7.GMC2LaurentShiftCheckA

/-!
# The concrete root-packet lemma — the additive-route core (bypasses THM-1550)

THM-2101's **root-packet lemma**: for `Φ` irreducible over a char-0 field `F` splitting in `L`, and a
subset `S` of its roots, the barycentric sum `b_k(S) = ∑_{α∈S} αᵏ/Φ'(α)` — *if it lies in `F`* — is
zero (`k ≤ deg Φ − 2`).  This is the additive analog of the orbit-product/Vieta core, and it removes
the small-root product, Hensel factorization, and Wiener–Hopf (THM-1550) from the DvdK proof.

The abstract additive machinery is codex's (`GMC2LaurentShiftCheckA`): `card_nsmul_translateSum_eq`
(orbit sum) and `sum_pow_pred_div_derivative_nodal_eq_zero` (Lagrange full-root sum = 0).  This module
instantiates it concretely at the **Galois action on the roots** (`GMC2GalRootAction`, my transitive
direct action), with the equivariant derivative weight — the additive analog of my S236 concrete
orbit-product instantiation.
-/

open Polynomial

namespace GMC2RootPacketConcrete

variable {F L : Type*} [Field F] [Field L] [Algebra F L]

/-- The barycentric derivative weight `w(α) = αᵏ / Φ'(α)` is Galois-equivariant: `w(σ•α) = σ(w α)`
(because `Φ' ∈ F[X]`, so `σ` commutes with evaluating it). -/
theorem weight_equivariant (Φ : F[X]) (k : ℕ) (σ : L ≃ₐ[F] L) (β : Φ.rootSet L) :
    ((σ • β : Φ.rootSet L) : L) ^ k / (aeval ((σ • β : Φ.rootSet L) : L) (derivative Φ))
      = σ ((β : L) ^ k / aeval (β : L) (derivative Φ)) := by
  rw [GMC2GalRootAction.coe_smul, map_div₀, map_pow]
  congr 1
  exact Polynomial.aeval_algHom_apply σ.toAlgHom (β : L) (derivative Φ)

/-- **The root-packet lemma.**  If the barycentric packet sum `b = ∑_{β∈S} βᵏ/Φ'(β)` lies in the base
field `F` (here presented as the image of `b0 : F`) and the *full*-root sum vanishes, then `b = 0`.

Assembles: transitivity (`isPretransitive_rootAction`, from irreducibility) + equivariance
(`weight_equivariant`) ⟹ every Galois translate of the packet sum equals `b`
(`card_nsmul_translateSum_eq`) ⟹ `|G|·b = |S|·|Stab|·(full sum) = 0` ⟹ `b = 0` in char 0. -/
theorem root_packet_eq_zero [CharZero F] [FiniteDimensional F L] [Normal F L]
    (Φ : F[X]) (hΦ : Irreducible Φ) (k : ℕ)
    (x0 : Φ.rootSet L) (S : Finset (Φ.rootSet L))
    (hfull : ∑ α : Φ.rootSet L, ((α : L) ^ k / aeval (α : L) (derivative Φ)) = 0)
    (b0 : F)
    (hb : algebraMap F L b0 = ∑ β ∈ S, ((β : L) ^ k / aeval (β : L) (derivative Φ))) :
    algebraMap F L b0 = 0 := by
  classical
  haveI : MulAction.IsPretransitive (L ≃ₐ[F] L) (Φ.rootSet L) :=
    GMC2GalRootAction.isPretransitive_rootAction Φ hΦ
  set w : Φ.rootSet L → L := fun α => (α : L) ^ k / aeval (α : L) (derivative Φ) with hw
  have htranslate : ∀ σ : L ≃ₐ[F] L, ∑ β ∈ S, w (σ • β) = algebraMap F L b0 := by
    intro σ
    calc ∑ β ∈ S, w (σ • β)
        = ∑ β ∈ S, σ (w β) :=
          Finset.sum_congr rfl (fun β _ => weight_equivariant Φ k σ β)
      _ = σ (∑ β ∈ S, w β) := by rw [← map_sum]
      _ = σ (algebraMap F L b0) := by rw [← hb]
      _ = algebraMap F L b0 := AlgEquiv.commutes σ b0
  have hkey := GMC2AdditiveOrbitSum.card_nsmul_translateSum_eq w S x0 (algebraMap F L b0) htranslate
  rw [hfull, nsmul_zero, nsmul_eq_mul] at hkey
  -- hkey : (|G| : L) * algebraMap F L b0 = 0, with (|G| : L) ≠ 0
  have hcardL : (Fintype.card (L ≃ₐ[F] L) : L) ≠ 0 := by
    rw [← map_natCast (algebraMap F L)]
    exact (map_ne_zero_iff _ (FaithfulSMul.algebraMap_injective F L)).mpr
      (Nat.cast_ne_zero.mpr Fintype.card_ne_zero)
  exact (mul_eq_zero.mp hkey).resolve_left hcardL

end GMC2RootPacketConcrete

#print axioms GMC2RootPacketConcrete.weight_equivariant
#print axioms GMC2RootPacketConcrete.root_packet_eq_zero
