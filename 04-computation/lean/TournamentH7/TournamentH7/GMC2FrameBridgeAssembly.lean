import TournamentH7.GMC2FrameBridge
import TournamentH7.GMC2FrameBridgePacket
import TournamentH7.GMC2Thm2067HSonly

/-!
# Frame-bridge assembly: `hS` from divisibility + the Weierstrass factor's constant coefficient

Combining the packet construction (`GMC2FrameBridgePacket.exists_packet_prod_eq`), Vieta
(`Splits.coeff_zero_eq_prod_roots_of_monic`), and the embedding pullback
(`GMC2FrameBridge.prod_eq_algebraMap_of_embedding`), the splitting-field `hS` identity follows from two
**concrete** inputs on the Weierstrass distinguished factor `Pω` (the image of `smallRootFactor` in the
Laurent frame `Ω`):

* `Pω ∣ Φ.map` over `Ω` — the divisibility supplied by the transpose of the Weierstrass factorization
  `Φ = P·h` (death-star's shared transpose hom), and
* `Pω.coeff 0 = algebraMap v` — the constant coefficient descends to the base field, with value `v`
  determined by the Weierstrass computation (`= ±c·t` under `hderiv`).

Then `∏_{β∈S} β = algebraMap ((-1)^{deg Pω}·v)`, exactly the `hS` fed to
`GMC2Thm2067HSonly.thm2067_reduced_to_hS`.  This is the whole bridge modulo those two frame-side facts —
kernel-pure, no valuation.
-/

open Polynomial

namespace GMC2FrameBridgeAssembly

variable {F : Type*} [Field F]

/-- **The frame bridge, reduced to divisibility + coefficient value.**  For `Φ ≠ 0` over `RatFunc F`,
an embedding `ψ` into `Ω`, and a monic, split, separable (nodup roots) polynomial `Pω` over `Ω`
dividing `Φ` there with `Pω.coeff 0 = algebraMap v`, the product of the packet of splitting-field roots
landing on `Pω` equals `algebraMap ((-1)^{deg Pω}·v)`. -/
theorem hS_of_dvd_value (Φ : (RatFunc F)[X]) (hΦ0 : Φ ≠ 0)
    {Ω : Type*} [Field Ω] [Algebra (RatFunc F) Ω]
    (ψ : Φ.SplittingField →ₐ[RatFunc F] Ω)
    (Pω : Polynomial Ω) (hmonic : Pω.Monic) (hPωsplit : Pω.Splits)
    (hPωnd : Pω.roots.Nodup)
    (hdvd : Pω ∣ Φ.map (algebraMap (RatFunc F) Ω))
    (v : RatFunc F) (hval : Pω.coeff 0 = algebraMap (RatFunc F) Ω v) :
    ∃ S : Finset (Φ.rootSet Φ.SplittingField),
      (∏ β ∈ S, (β : Φ.SplittingField))
        = algebraMap (RatFunc F) Φ.SplittingField ((-1) ^ Pω.natDegree * v) := by
  obtain ⟨S, hS⟩ := GMC2FrameBridgePacket.exists_packet_prod_eq Φ hΦ0 ψ Pω hPωnd hdvd
  refine ⟨S, ?_⟩
  refine GMC2FrameBridge.prod_eq_algebraMap_of_embedding Φ ψ S ((-1) ^ Pω.natDegree * v) ?_
  rw [hS]
  have hvieta : Pω.roots.prod = (-1) ^ Pω.natDegree * Pω.coeff 0 := by
    rw [hPωsplit.coeff_zero_eq_prod_roots_of_monic hmonic, ← mul_assoc, ← pow_add, ← two_mul,
      pow_mul, neg_one_sq, one_pow, one_mul]
  rw [hvieta, hval, map_mul, map_pow, map_neg, map_one]

/-- **The frame bridge, completed into the DvdK contradiction.**  For the concrete `Φ = Phi R M`
(`1 ≤ M < deg R`, `R(0) ≠ 0`) over a characteristic-zero field, given the two frame-side facts —
`Pω ∣ Φ` over `Ω` (the transpose of the Weierstrass factorization) and `Pω.coeff 0 = algebraMap
((-1)^{deg Pω}·(c·t))` (the Weierstrass value under `hderiv`) — the orbit-product contradiction
`GMC2Thm2067HSonly.thm2067_reduced_to_hS` closes: `False`.  This is `SinglePolyCrux`'s contradiction
reduced to **exactly** the transpose-provided divisibility and the `hderiv`-provided value; everything
else is kernel-pure and valuation-free.  (The value is packaged so that `(-1)^{deg Pω}` cancels, giving
the `hS` shape `∏_{β∈S} β = algebraMap (c·t)` directly.) -/
theorem false_of_frame_data [CharZero F] (R : Polynomial F) (M : ℕ)
    (hM : 1 ≤ M) (hMd : M < R.natDegree) (hR0 : R.coeff 0 ≠ 0)
    {Ω : Type*} [Field Ω] [Algebra (RatFunc F) Ω]
    (ψ : (GMC2PhiVieta.Phi R M).SplittingField →ₐ[RatFunc F] Ω)
    (Pω : Polynomial Ω) (hmonic : Pω.Monic) (hPωsplit : Pω.Splits) (hPωnd : Pω.roots.Nodup)
    (hdvd : Pω ∣ (GMC2PhiVieta.Phi R M).map (algebraMap (RatFunc F) Ω))
    (c : F) (hc : c ≠ 0)
    (hval : Pω.coeff 0
        = algebraMap (RatFunc F) Ω ((-1) ^ Pω.natDegree * (RatFunc.C c * RatFunc.X)))
    (x0 : (GMC2PhiVieta.Phi R M).rootSet (GMC2PhiVieta.Phi R M).SplittingField) :
    False := by
  have hΦ0 : GMC2PhiVieta.Phi R M ≠ 0 := (GMC2DvdKAssembly.irreducible_Phi R M hM hR0).ne_zero
  obtain ⟨S, hS⟩ := hS_of_dvd_value (GMC2PhiVieta.Phi R M) hΦ0 ψ Pω hmonic hPωsplit hPωnd hdvd
    ((-1) ^ Pω.natDegree * (RatFunc.C c * RatFunc.X)) hval
  rw [show (-1 : RatFunc F) ^ Pω.natDegree * ((-1) ^ Pω.natDegree * (RatFunc.C c * RatFunc.X))
        = RatFunc.C c * RatFunc.X by
      rw [← mul_assoc, ← pow_add, ← two_mul, pow_mul, neg_one_sq, one_pow, one_mul]] at hS
  exact GMC2Thm2067HSonly.thm2067_reduced_to_hS R M hM hMd hR0 S x0 c hc hS

end GMC2FrameBridgeAssembly

#print axioms GMC2FrameBridgeAssembly.hS_of_dvd_value
#print axioms GMC2FrameBridgeAssembly.false_of_frame_data
