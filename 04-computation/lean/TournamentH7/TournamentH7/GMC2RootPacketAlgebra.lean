import TournamentH7.GMC2LaurentShiftCheckA

/-!
# Algebraic rigidity of an additive DvdK root packet

This module isolates the algebraic half of THM-2101's additive proof.  A
transitive finite group action cannot carry a proper packet whose equivariant
weight has a nonzero base-field sum while the full-orbit sum is zero.  The
specialization below uses the Galois action on the roots of an irreducible
polynomial.

The intended downstream use is deliberately narrower than
`RootPacketContourPremise`: the remaining analytic or t-adic input should only
select a finite root packet and identify its barycentric sum with `1`.  It
should not assume the desired nonvanishing conclusion.
-/

open Polynomial
open scoped BigOperators
open MulAction

namespace GMC2RootPacketAlgebra

noncomputable section

variable {G Ω K : Type*} [Group G] [Fintype G] [MulAction G Ω]
  [Fintype Ω] [DecidableEq Ω] [Field K] [CharZero K]
  [MulSemiringAction G K]

/-- A transitive equivariant packet whose sum is fixed by the group and whose
full-orbit sum is zero must itself have sum zero.  This is the additive
orbit-incidence engine with the contradiction removed from its statement. -/
theorem packetSum_eq_zero_of_fixed [IsPretransitive G Ω]
    (f : Ω → K) (S : Finset Ω) (x : Ω)
    (hf : ∀ (g : G) (a : Ω), f (g • a) = g • f a)
    (hfixed : ∀ g : G, g • (∑ a ∈ S, f a) = ∑ a ∈ S, f a)
    (hfull : ∑ a : Ω, f a = 0) :
    ∑ a ∈ S, f a = 0 := by
  let b : K := ∑ a ∈ S, f a
  have htranslate : ∀ g : G, ∑ a ∈ S, f (g • a) = b := by
    intro g
    calc
      ∑ a ∈ S, f (g • a) = ∑ a ∈ S, g • f a := by
        apply Finset.sum_congr rfl
        intro a ha
        exact hf g a
      _ = g • (∑ a ∈ S, f a) := by
        change (∑ a ∈ S, g • f a) = g • (∑ a ∈ S, f a)
        rw [Finset.smul_sum]
      _ = b := by simpa [b] using hfixed g
  have hincidence := GMC2AdditiveOrbitSum.card_nsmul_translateSum_eq
    f S x b htranslate
  rw [hfull, nsmul_zero] at hincidence
  have hcard : (Fintype.card G : K) ≠ 0 :=
    Nat.cast_ne_zero.mpr Fintype.card_ne_zero
  rw [nsmul_eq_mul] at hincidence
  simpa [b] using (mul_eq_zero.mp hincidence).resolve_left hcard

variable {F : Type*} [Field F] [CharZero F]

/-- Galois specialization of `packetSum_eq_zero_of_fixed`.  If an equivariant
weight on the roots of an irreducible polynomial has total sum zero, then any
root packet whose sum descends to the base field has base value zero. -/
theorem galoisPacket_baseValue_eq_zero
    (p : F[X]) (hp : Irreducible p)
    (f : p.rootSet p.SplittingField → p.SplittingField)
    (hf : ∀ (g : p.Gal) (a : p.rootSet p.SplittingField),
      f (g • a) = g • f a)
    (hfull : ∑ a : p.rootSet p.SplittingField, f a = 0)
    (S : Finset (p.rootSet p.SplittingField)) (b : F)
    (hpacket : ∑ a ∈ S, f a = algebraMap F p.SplittingField b) :
    b = 0 := by
  classical
  letI : IsPretransitive p.Gal (p.rootSet p.SplittingField) := ⟨by
    intro x y
    have hx := minpoly.eq_of_irreducible hp (mem_rootSet.mp x.2).2
    have hy := minpoly.eq_of_irreducible hp (mem_rootSet.mp y.2).2
    obtain ⟨g, hg⟩ :=
      (Normal.minpoly_eq_iff_mem_orbit p.SplittingField).mp (hy.symm.trans hx)
    exact ⟨g, Subtype.ext hg⟩
  ⟩
  have hrootCard : 0 < Fintype.card (p.rootSet p.SplittingField) := by
    rw [Polynomial.card_rootSet_eq_natDegree hp.separable (SplittingField.splits p)]
    exact hp.natDegree_pos
  letI : Nonempty (p.rootSet p.SplittingField) :=
    Fintype.card_pos_iff.mp hrootCard
  have hfixed : ∀ g : p.Gal,
      g • (∑ a ∈ S, f a) = ∑ a ∈ S, f a := by
    intro g
    rw [hpacket]
    exact g.commutes b
  have hzero : ∑ a ∈ S, f a = 0 :=
    packetSum_eq_zero_of_fixed f S (Classical.arbitrary _) hf hfixed hfull
  rw [hpacket] at hzero
  apply (algebraMap F p.SplittingField).injective
  simpa using hzero

end

end GMC2RootPacketAlgebra

#print axioms GMC2RootPacketAlgebra.packetSum_eq_zero_of_fixed
#print axioms GMC2RootPacketAlgebra.galoisPacket_baseValue_eq_zero
