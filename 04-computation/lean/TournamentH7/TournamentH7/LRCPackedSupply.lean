/-
  TournamentH7.LRCPackedSupply -- THE PACKED-SUPPLY DICHOTOMY
  (klein-2026-07-11-S250, HYP-5975, THM-698's Lean half): the certificate
  supply's conclusion proved for EVERY packed two-scale shape, by combining
  the two complementary witness routes through the P-DICHOTOMY:

    for |P| <= 5 with |P| + |E| = 13:  EITHER some q* in [8,13] avoids P
    with |E| < q* (the TEST-POINT route, THM-693/696), OR the room is denied
    -- which forces P = [|E|+1, 13] EXACTLY (interval containment + card
    equality), hence pmin >= 9 and the FIRST-WINDOW route (THM-697) applies.

  packed_supply: P subset [1,13], |P| <= 5, |P| + |E| = 13, co-offsets
  0 <= e <= 90 with 168e < V, V > 2366  ==>  the THM527ACertificateSupply
  conclusion (an interval with thirteen SafeIvStrict certificates).

  With THM-698's audit (the supply's domain IS uniformly two-scale;
  witnessG2 = meas(G_P minus D(E)); its positivity = the dead-zone theorem),
  the supply's remaining content is exactly the spread-cluster realization
  at moderate V.

  Kernel-pure: no native_decide, no sorry.
-/
import Mathlib
import TournamentH7.LRCTestDataSupply
import TournamentH7.LRCFirstWindowWitness

namespace LonelyRunner
namespace PackedSupply

open LRC14Grand MeasureTransfer

/-- **The window-route supply** for top-block-type shapes (`pmin ≥ 9`):
the leftmost 6/7-phase multiplier serves every small speed in `[9,13]` and
every co-offset `≤ 90`, then the fattening bridge produces the data. -/
theorem firstWindow_supply (v : Fin 13 → ℤ) (V : ℤ) (P E : Finset ℤ)
    (hV : 2366 < V)
    (hPsub : ∀ p ∈ P, 9 ≤ p ∧ p ≤ 13)
    (hE : ∀ e ∈ E, 0 ≤ e ∧ e ≤ 90)
    (hv : ∀ i, v i ∈ P ∨ ∃ e ∈ E, v i = V - e) :
    ∃ D x y q : ℤ, 0 < D ∧ 0 < q ∧ 0 ≤ x ∧ y < D ∧ D < (y - x) * q ∧
      ∀ i, SafeIvStrict (v i) D x y := by
  obtain ⟨j, hj0, hjlo, hjhi⟩ :=
    FirstWindowWitness.firstWindow_j_exists V 9 (by norm_num) (by omega)
  have hw0 : (0 : ℤ) < 7 * j + 6 := by omega
  have hwq : 7 * j + 6 < 7 * V := by nlinarith [hjhi]
  have hlive : StrictlyLive v (7 * V) (7 * j + 6) := by
    refine FirstWindowWitness.firstWindow_strictlyLive v V j (by omega) hj0
      hwq ?_
    intro i
    rcases hv i with hp | ⟨e, he, hve⟩
    · obtain ⟨hp9, hp13⟩ := hPsub _ hp
      left
      constructor
      · nlinarith [hjlo]
      · nlinarith [hjhi]
    · obtain ⟨he0, he90⟩ := hE e he
      exact Or.inr ⟨e, he0, hve, by nlinarith [hjhi]⟩
  refine TestDataSupply.supply_of_strictlyLive hlive (fun i => ?_)
    (by omega : (0 : ℤ) < V)
  rcases hv i with hp | ⟨e, he, hve⟩
  · obtain ⟨h9, h13⟩ := hPsub _ hp
    exact ⟨by omega, by omega⟩
  · obtain ⟨he0, he90⟩ := hE e he
    rw [hve]
    exact ⟨by omega, by omega⟩

/-- **The P-dichotomy**: either the test-point route has room, or `P` is
exactly the top block `[|E|+1, 13]` — forcing every small speed `≥ 9`. -/
theorem packed_dichotomy (P E : Finset ℤ) (hcard : P.card ≤ 5)
    (hfull : (P.card : ℤ) + (E.card : ℤ) = 13) :
    (∃ qs : ℤ, 8 ≤ qs ∧ qs ≤ 13 ∧ qs ∉ P ∧ (E.card : ℤ) < qs) ∨
    (∀ p ∈ P, 9 ≤ p) := by
  by_cases h : ∃ qs : ℤ, 8 ≤ qs ∧ qs ≤ 13 ∧ qs ∉ P ∧ (E.card : ℤ) < qs
  · exact Or.inl h
  · right
    push_neg at h
    have hE8 : (8 : ℤ) ≤ (E.card : ℤ) := by omega
    have hsub : Finset.Icc ((E.card : ℤ) + 1) 13 ⊆ P := by
      intro x hx
      rw [Finset.mem_Icc] at hx
      by_contra hxP
      have := h x (by omega) (by omega) hxP
      omega
    have hcardIcc : (Finset.Icc ((E.card : ℤ) + 1) 13).card = P.card := by
      rw [Int.card_Icc]
      have h1 : (13 + 1 - ((E.card : ℤ) + 1)) = (P.card : ℤ) := by omega
      rw [h1]
      exact Int.toNat_natCast _
    have hPeq : P = Finset.Icc ((E.card : ℤ) + 1) 13 :=
      (Finset.eq_of_subset_of_card_le hsub hcardIcc.ge).symm
    intro p hp
    rw [hPeq, Finset.mem_Icc] at hp
    omega

/-- **THE PACKED SUPPLY (THM-698's Lean half)**: the certificate-supply
conclusion for EVERY packed two-scale shape — the two witness routes,
dispatched by the dichotomy. -/
theorem packed_supply (v : Fin 13 → ℤ) (V : ℤ) (P E : Finset ℤ)
    (hV : 2366 < V)
    (hPsub : ∀ p ∈ P, 1 ≤ p ∧ p ≤ 13) (hPcard : P.card ≤ 5)
    (hfull : (P.card : ℤ) + (E.card : ℤ) = 13)
    (hE : ∀ e ∈ E, 0 ≤ e ∧ e ≤ 90 ∧ 168 * e < V)
    (hv : ∀ i, v i ∈ P ∨ ∃ e ∈ E, v i = V - e) :
    ∃ D x y q : ℤ, 0 < D ∧ 0 < q ∧ 0 ≤ x ∧ y < D ∧ D < (y - x) * q ∧
      ∀ i, SafeIvStrict (v i) D x y := by
  rcases packed_dichotomy P E hPcard hfull with hqs | h9
  · exact TestDataSupply.twoScale_supply v V P E (by omega) hPsub hqs
      (fun e he => ⟨(hE e he).1, (hE e he).2.2⟩) hv
  · exact firstWindow_supply v V P E hV
      (fun p hp => ⟨h9 p hp, (hPsub p hp).2⟩)
      (fun e he => ⟨(hE e he).1, (hE e he).2.1⟩) hv

end PackedSupply
end LonelyRunner

-- kernel-purity audit (fleet convention)
#print axioms LonelyRunner.PackedSupply.firstWindow_supply
#print axioms LonelyRunner.PackedSupply.packed_dichotomy
#print axioms LonelyRunner.PackedSupply.packed_supply
