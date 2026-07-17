/-
  TournamentH7.LRCDeepReflectionParity

  Reflection pairs every exact-depth multiplier except the possible midpoint
  of an even modulus.  Thus the parity of an exact-depth stratum is exactly
  its midpoint membership.  At q = 2m the midpoint depth is, more concretely,
  the number of even speeds.  In particular, at odd modulus the two
  exceptional weighted-census strata at depths six and seven occur in
  two-event quanta.  This does not pay those strata, but removes singleton
  tails from any subsequent transport argument.

  Tournament-analysis audit: vertices are multiplier events, the observable
  is exact band depth, and `p ↦ q-p` is the gauge involution.  The resulting
  orbit graph is a perfect matching at odd modulus; a runner tournament would
  forget both multiplier activity and the reflection orbit.

  No `sorry`; no `native_decide`.
-/

import TournamentH7.BucketBalance
import TournamentH7.LRCResonanceNucleus
import TournamentH7.LRCWeightedDeepCensus

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- At every modulus, reflection pairs all members of an exact-depth stratum
except its possible fixed point.  The right side is therefore either zero or
one, but retaining it as a cardinality makes the fixed-point statement valid
without separate odd/even hypotheses. -/
theorem exactDepthCount_mod_two_eq_midpointCount
    (v : Fin 13 → ℤ) (q depth : ℕ) :
    exactDepthCount v q depth % 2 =
      ((Finset.Ioo 0 q).filter fun p =>
        bandCount v q p = depth ∧ 2 * p = q).card % 2 := by
  classical
  let s := (Finset.Ioo 0 q).filter fun p => bandCount v q p = depth
  let paired := s.filter fun p => 2 * p ≠ q
  have hpaired : Even paired.card := by
    let reflect : ℕ → ℕ := fun p => q - p
    apply Tournament.BucketBalance.even_card_of_fixedPointFree_involutiveOn
      paired reflect
    · intro p hp
      unfold paired s at hp ⊢
      simp only [Finset.mem_filter, Finset.mem_Ioo] at hp ⊢
      obtain ⟨⟨⟨hp0, hpq⟩, hdepth⟩, hnotmid⟩ := hp
      have hq : 0 < q := hp0.trans hpq
      refine ⟨⟨⟨Nat.sub_pos_of_lt hpq, Nat.sub_lt hq hp0⟩, ?_⟩, ?_⟩
      · rw [bandCount_reflect v q p hq hp0 hpq]
        exact hdepth
      · unfold reflect
        omega
    · intro p hp
      unfold reflect
      unfold paired s at hp
      simp only [Finset.mem_filter, Finset.mem_Ioo] at hp
      omega
    · intro p hp hfixed
      unfold reflect at hfixed
      unfold paired s at hp
      simp only [Finset.mem_filter, Finset.mem_Ioo] at hp
      omega
  have hsplit :
      (s.filter fun p => 2 * p = q).card + paired.card = s.card := by
    simpa [paired] using
      (Finset.card_filter_add_card_filter_not
        (s := s) (p := fun p => 2 * p = q))
  have hmid :
      s.filter (fun p => 2 * p = q) =
        (Finset.Ioo 0 q).filter fun p =>
          bandCount v q p = depth ∧ 2 * p = q := by
    ext p
    simp [s, and_assoc]
  unfold exactDepthCount
  rw [← hmid]
  change s.card % 2 = (s.filter fun p => 2 * p = q).card % 2
  rcases hpaired with ⟨k, hk⟩
  omega

/-- A runner is bad at the half multiplier precisely when its speed is even.
For an odd speed the residue is `m`; for an even speed it is zero. -/
theorem midpoint_not_inBand_iff_evenSpeed
    (v : Fin 13 → ℤ) (m : ℕ) (hm : 0 < m) (i : Fin 13) :
    ¬ inBand v (m + m) m i ↔ (2 : ℤ) ∣ v i := by
  have hq : ((m + m : ℕ) : ℤ) = 2 * (m : ℤ) := by
    push_cast
    ring
  constructor
  · intro hbad
    rcases Int.even_or_odd (v i) with heven | hodd
    · exact heven.two_dvd
    · exfalso
      apply hbad
      obtain ⟨c, hc⟩ : ∃ c, v i = 2 * c + 1 := hodd
      have hres :
          (v i * (m : ℤ)) % ((m + m : ℕ) : ℤ) = (m : ℤ) := by
        rw [hq, hc]
        have hsplit :
            (2 * c + 1) * (m : ℤ) =
              (m : ℤ) + (2 * (m : ℤ)) * c := by
          ring
        rw [hsplit, Int.add_mul_emod_self_left]
        exact Int.emod_eq_of_lt (by positivity) (by omega)
      unfold inBand
      rw [hres]
      constructor <;> push_cast <;> omega
  · rintro ⟨c, hc⟩ hband
    have hres :
        (v i * (m : ℤ)) % ((m + m : ℕ) : ℤ) = 0 := by
      rw [hq, hc]
      have hmultiple : 2 * c * (m : ℤ) = c * (2 * (m : ℤ)) := by
        ring
      rw [hmultiple]
      exact Int.mul_emod_left _ _
    obtain ⟨hlower, _⟩ := hband
    rw [hres, hq] at hlower
    omega

/-- The exact midpoint depth at `q = 2m` is the number of even speeds. -/
theorem bandCount_midpoint_eq_evenSpeedCount
    (v : Fin 13 → ℤ) (m : ℕ) (hm : 0 < m) :
    bandCount v (m + m) m =
      (Finset.univ.filter fun i : Fin 13 => (2 : ℤ) ∣ v i).card := by
  unfold bandCount
  apply congrArg Finset.card
  ext i
  simp only [Finset.mem_filter, Finset.mem_univ, true_and]
  exact midpoint_not_inBand_iff_evenSpeed v m hm i

/-- At a positive even modulus, an exact-depth stratum is odd exactly when
the midpoint itself has that depth. -/
theorem exactDepthCount_mod_two_even
    (v : Fin 13 → ℤ) (m depth : ℕ) (hm : 0 < m) :
    exactDepthCount v (m + m) depth % 2 =
      if bandCount v (m + m) m = depth then 1 else 0 := by
  rw [exactDepthCount_mod_two_eq_midpointCount]
  by_cases hdepth : bandCount v (m + m) m = depth
  · rw [if_pos hdepth]
    have hset :
        ((Finset.Ioo 0 (m + m)).filter fun p =>
          bandCount v (m + m) p = depth ∧ 2 * p = m + m) = {m} := by
      ext p
      simp only [Finset.mem_filter, Finset.mem_Ioo, Finset.mem_singleton]
      constructor
      · rintro ⟨⟨hp0, hpq⟩, _, hfixed⟩
        omega
      · rintro rfl
        exact ⟨⟨hm, by omega⟩, hdepth, by omega⟩
    rw [hset]
    simp
  · rw [if_neg hdepth]
    have hset :
        ((Finset.Ioo 0 (m + m)).filter fun p =>
          bandCount v (m + m) p = depth ∧ 2 * p = m + m) = ∅ := by
      apply Finset.eq_empty_iff_forall_notMem.mpr
      intro p hp
      simp only [Finset.mem_filter, Finset.mem_Ioo] at hp
      obtain ⟨_, hpdepth, hfixed⟩ := hp
      have hp : p = m := by omega
      subst p
      exact hdepth hpdepth
    rw [hset]
    simp

/-- Substituting the midpoint formula: the only unpaired exact-depth event
is present precisely when `depth` is the number of even speeds. -/
theorem exactDepthCount_mod_two_even_eq_evenSpeedCount
    (v : Fin 13 → ℤ) (m depth : ℕ) (hm : 0 < m) :
    exactDepthCount v (m + m) depth % 2 =
      if (Finset.univ.filter fun i : Fin 13 => (2 : ℤ) ∣ v i).card = depth
      then 1 else 0 := by
  rw [exactDepthCount_mod_two_even v m depth hm,
    bandCount_midpoint_eq_evenSpeedCount v m hm]

/-- Every exact-depth multiplier stratum has even cardinality at odd modulus.
The reflection action has no midpoint because an odd integer is not twice a
natural number. -/
theorem exactDepthCount_even_of_odd
    (v : Fin 13 → ℤ) (q depth : ℕ) (hqodd : Odd q) :
    Even (exactDepthCount v q depth) := by
  unfold exactDepthCount
  let s := (Finset.Ioo 0 q).filter fun p => bandCount v q p = depth
  let reflect : ℕ → ℕ := fun p => q - p
  apply Tournament.BucketBalance.even_card_of_fixedPointFree_involutiveOn
    s reflect
  · intro p hp
    unfold s at hp ⊢
    obtain ⟨hpWindow, hdepth⟩ := Finset.mem_filter.mp hp
    have hp0 : 0 < p := (Finset.mem_Ioo.mp hpWindow).1
    have hpq : p < q := (Finset.mem_Ioo.mp hpWindow).2
    have hq : 0 < q := hp0.trans hpq
    refine Finset.mem_filter.mpr ⟨?_, ?_⟩
    · rw [Finset.mem_Ioo]
      exact ⟨Nat.sub_pos_of_lt hpq, Nat.sub_lt hq hp0⟩
    · rw [bandCount_reflect v q p hq hp0 hpq]
      exact hdepth
  · intro p hp
    unfold reflect
    unfold s at hp
    have hpWindow := (Finset.mem_filter.mp hp).1
    rw [Finset.mem_Ioo] at hpWindow
    omega
  · intro p hp hfixed
    unfold reflect at hfixed
    unfold s at hp
    have hpWindow := (Finset.mem_filter.mp hp).1
    rw [Finset.mem_Ioo] at hpWindow
    obtain ⟨k, hk⟩ := hqodd
    omega

/-- The depth-six residue in the weighted B5 census is even at odd modulus. -/
theorem exactDepthCount_six_even_of_odd
    (v : Fin 13 → ℤ) (q : ℕ) (hqodd : Odd q) :
    Even (exactDepthCount v q 6) :=
  exactDepthCount_even_of_odd v q 6 hqodd

/-- The depth-seven residue in the weighted B5 census is even at odd modulus. -/
theorem exactDepthCount_seven_even_of_odd
    (v : Fin 13 → ℤ) (q : ℕ) (hqodd : Odd q) :
    Even (exactDepthCount v q 7) :=
  exactDepthCount_even_of_odd v q 7 hqodd

/-! ## Axiom audit -/

#print axioms exactDepthCount_mod_two_eq_midpointCount
#print axioms midpoint_not_inBand_iff_evenSpeed
#print axioms bandCount_midpoint_eq_evenSpeedCount
#print axioms exactDepthCount_mod_two_even
#print axioms exactDepthCount_mod_two_even_eq_evenSpeedCount
#print axioms exactDepthCount_even_of_odd
#print axioms exactDepthCount_six_even_of_odd
#print axioms exactDepthCount_seven_even_of_odd

end LRC14Concrete
end LonelyRunner
