/-
  TournamentH7.LRCParityPairing -- the parity pairing law (klein-2026-07-09-S227,
  formalizing klein-S222 / HYP-5825; the Redei-involution transplant).

  THE LAW: p ↦ q − p is a live-preserving involution on multipliers (the safe
  band is mirror-symmetric), whose unique possible fixed point is q/2 -- and
  q/2 is live iff EVERY speed is odd.  Hence liveCount ≡ [q even ∧ all speeds
  odd] (mod 2), so for covering sets (which contain an even speed) the live
  count is EVEN at every modulus: live multipliers come in exact ± pairs
  (every certificate exhibiting p yields q − p free), and an odd live count
  at a covering set is a computation-bug detector.

  Kernel-pure: no native_decide, no sorry.
-/

import Mathlib
import TournamentH7.LRCDiscreteBonferroni

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- The live set as a Finset (bandCount = 0 on the open window (0, q)). -/
def liveSet (v : Fin 13 → ℤ) (q : ℕ) : Finset ℕ :=
  (Finset.Ioo 0 q).filter fun p => bandCount v q p = 0

theorem liveCount_eq_card (v : Fin 13 → ℤ) (q : ℕ) :
    liveCount v q = (liveSet v q).card := rfl

/-- Band mirror symmetry: if `v·p` has residue `x` in the band, then `v·(q−p)`
has residue `q − x`, also in the band.  Route: the two residues sum to `0`
mod `q`, both lie in `[0, q)`, and the band forces `x ≠ 0` -- omega closes. -/
theorem inBand_mirror (v : Fin 13 → ℤ) (q p : ℕ) (hq : 0 < q) (hp : p < q)
    (i : Fin 13) (h : inBand v q p i) : inBand v q (q - p) i := by
  obtain ⟨h1, h2⟩ := h
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  set x : ℤ := (v i * (p : ℤ)) % (q : ℤ) with hxdef
  set y : ℤ := (v i * ((q - p : ℕ) : ℤ)) % (q : ℤ) with hydef
  have hx0 : 0 ≤ x := Int.emod_nonneg _ (by omega)
  have hxq : x < (q : ℤ) := Int.emod_lt_of_pos _ hqZ
  have hy0 : 0 ≤ y := Int.emod_nonneg _ (by omega)
  have hyq : y < (q : ℤ) := Int.emod_lt_of_pos _ hqZ
  have hxpos : 0 < x := by
    rcases lt_or_eq_of_le hx0 with h | h
    · exact h
    · exfalso; rw [← h] at h1; omega
  -- the residues sum to 0 mod q
  have hsum : ((x + y) % (q : ℤ)) = 0 := by
    have hcast : ((q - p : ℕ) : ℤ) = (q : ℤ) - (p : ℤ) := by
      exact_mod_cast Nat.cast_sub (le_of_lt hp)
    have htot : v i * (p : ℤ) + v i * ((q - p : ℕ) : ℤ) = v i * (q : ℤ) := by
      rw [hcast]; ring
    calc (x + y) % (q : ℤ)
        = (v i * (p : ℤ) + v i * ((q - p : ℕ) : ℤ)) % (q : ℤ) := by
          rw [hxdef, hydef, Int.add_emod, Int.emod_emod_of_dvd _ dvd_rfl,
              Int.emod_emod_of_dvd _ dvd_rfl, ← Int.add_emod]
      _ = (v i * (q : ℤ)) % (q : ℤ) := by rw [htot]
      _ = 0 := Int.mul_emod_left _ _
  obtain ⟨k, hk⟩ : (q : ℤ) ∣ (x + y) := Int.dvd_of_emod_eq_zero hsum
  have hy : y = (q : ℤ) - x := by
    have hk01 : k = 0 ∨ k = 1 := by
      rcases lt_trichotomy k 0 with h | h | h
      · exfalso; nlinarith
      · exact Or.inl h
      · have : k = 1 := by nlinarith
        exact Or.inr this
    rcases hk01 with rfl | rfl
    · exfalso; simp at hk; omega
    · omega
  rw [hydef] at hy
  simp only [inBand]
  constructor <;> omega

/-- Liveness is mirror-symmetric on the open window. -/
theorem live_mirror (v : Fin 13 → ℤ) (q p : ℕ) (hp0 : 0 < p) (hpq : p < q)
    (h : bandCount v q p = 0) : bandCount v q (q - p) = 0 := by
  have hq : 0 < q := lt_trans hp0 hpq
  unfold bandCount at h ⊢
  rw [Finset.card_eq_zero, Finset.filter_eq_empty_iff] at h ⊢
  intro i hi
  have hband : inBand v q p i := by
    have := h hi
    simpa using this
  simpa using inBand_mirror v q p hq hpq i hband

/-- The half-point law: for `q = m + m` (`m > 0`), the multiplier `m` is live
iff every speed is odd. -/
theorem half_live_iff (v : Fin 13 → ℤ) (m : ℕ) (hm : 0 < m) :
    bandCount v (m + m) m = 0 ↔ ∀ i, ¬ (2 : ℤ) ∣ v i := by
  have hq : ((m + m : ℕ) : ℤ) = 2 * (m : ℤ) := by push_cast; ring
  constructor
  · intro h i h2
    unfold bandCount at h
    rw [Finset.card_eq_zero, Finset.filter_eq_empty_iff] at h
    have hband : inBand v (m + m) m i := by
      have := h (Finset.mem_univ i)
      simpa using this
    obtain ⟨c, hc⟩ := h2
    obtain ⟨h1, _⟩ := hband
    have hres : (v i * (m : ℤ)) % ((m + m : ℕ) : ℤ) = 0 := by
      rw [hq, hc]
      have : 2 * c * (m : ℤ) = c * (2 * (m : ℤ)) := by ring
      rw [this]
      exact Int.mul_emod_left _ _
    rw [hres] at h1
    have : ((m + m : ℕ) : ℤ) ≤ 0 := by simpa using h1
    rw [hq] at this
    omega
  · intro hodd
    unfold bandCount
    rw [Finset.card_eq_zero, Finset.filter_eq_empty_iff]
    intro i _
    simp only [not_not, inBand]
    obtain ⟨c, hc⟩ : ∃ c, v i = 2 * c + 1 := by
      rcases Int.even_or_odd (v i) with he | ho
      · exact absurd he.two_dvd (hodd i)
      · exact ho
    have hres : (v i * (m : ℤ)) % ((m + m : ℕ) : ℤ) = (m : ℤ) := by
      rw [hq, hc]
      have hsplit : (2 * c + 1) * (m : ℤ) = (m : ℤ) + (2 * (m : ℤ)) * c := by ring
      rw [hsplit, Int.add_mul_emod_self_left]
      exact Int.emod_eq_of_lt (by positivity) (by push_cast; omega)
    constructor <;> rw [hres] <;> push_cast <;> omega

/-- **THE PARITY PAIRING LAW**: the live count's parity equals the half-point's
liveness (the ± involution pairs everything else). -/
theorem liveCount_parity (v : Fin 13 → ℤ) (q : ℕ) :
    liveCount v q % 2
      = ((Finset.Ioo 0 q).filter fun p => bandCount v q p = 0 ∧ 2 * p = q).card % 2 := by
  classical
  rw [liveCount_eq_card]
  set L := liveSet v q with hL
  have hsplit1 : (L.filter fun p => 2 * p < q).card
      + (L.filter fun p => ¬ 2 * p < q).card = L.card :=
    Finset.filter_card_add_filter_neg_card_eq_card _
  have hsplit2 : ((L.filter fun p => ¬ 2 * p < q).filter fun p => 2 * p = q).card
      + ((L.filter fun p => ¬ 2 * p < q).filter fun p => ¬ 2 * p = q).card
      = (L.filter fun p => ¬ 2 * p < q).card :=
    Finset.filter_card_add_filter_neg_card_eq_card _
  -- the strict-upper half bijects with the strict-lower half via p ↦ q − p
  have hbij : (L.filter fun p => 2 * p < q).card
      = ((L.filter fun p => ¬ 2 * p < q).filter fun p => ¬ 2 * p = q).card := by
    apply Finset.card_bij' (i := fun p _ => q - p) (j := fun p _ => q - p)
    · intro p hp
      simp only [hL, liveSet, Finset.mem_filter, Finset.mem_Ioo] at hp ⊢
      obtain ⟨⟨⟨hp0, hpq⟩, hlive⟩, hlt⟩ := hp
      exact ⟨⟨⟨⟨by omega, by omega⟩, live_mirror v q p hp0 hpq hlive⟩, by omega⟩, by omega⟩
    · intro p hp
      simp only [hL, liveSet, Finset.mem_filter, Finset.mem_Ioo] at hp ⊢
      obtain ⟨⟨⟨⟨hp0, hpq⟩, hlive⟩, hge⟩, hne⟩ := hp
      exact ⟨⟨⟨by omega, by omega⟩, live_mirror v q p hp0 hpq hlive⟩, by omega⟩
    · intro p hp
      simp only [hL, liveSet, Finset.mem_filter, Finset.mem_Ioo] at hp
      omega
    · intro p hp
      simp only [hL, liveSet, Finset.mem_filter, Finset.mem_Ioo] at hp
      omega
  have hmid : (L.filter fun p => ¬ 2 * p < q).filter (fun p => 2 * p = q)
      = (Finset.Ioo 0 q).filter fun p => bandCount v q p = 0 ∧ 2 * p = q := by
    rw [hL]
    unfold liveSet
    rw [Finset.filter_filter, Finset.filter_filter]
    apply Finset.filter_congr
    intro p _
    constructor
    · rintro ⟨hlive, _, heq⟩; exact ⟨hlive, heq⟩
    · rintro ⟨hlive, heq⟩; exact ⟨hlive, by omega, heq⟩
  have hmidcard : ((L.filter fun p => ¬ 2 * p < q).filter fun p => 2 * p = q).card
      = ((Finset.Ioo 0 q).filter fun p => bandCount v q p = 0 ∧ 2 * p = q).card := by
    rw [hmid]
  omega

/-- **Covering ⟹ even live count**: a family with an even speed has a dead
half-point, so live multipliers come in exact ± pairs at EVERY modulus. -/
theorem liveCount_even_of_even_speed (v : Fin 13 → ℤ) (q : ℕ)
    (heven : ∃ i, (2 : ℤ) ∣ v i) : liveCount v q % 2 = 0 := by
  rw [liveCount_parity]
  have hzero : ((Finset.Ioo 0 q).filter
      fun p => bandCount v q p = 0 ∧ 2 * p = q).card = 0 := by
    rw [Finset.card_eq_zero, Finset.filter_eq_empty_iff]
    intro p hp
    rw [Finset.mem_Ioo] at hp
    rintro ⟨hlive, heq⟩
    rcases Nat.even_or_odd q with hev | hodd
    · obtain ⟨m, hm⟩ := hev
      have hm0 : 0 < m := by omega
      have hpm : p = m := by omega
      rw [hpm, hm] at hlive
      rw [half_live_iff v m hm0] at hlive
      obtain ⟨i, hi⟩ := heven
      exact hlive i hi
    · obtain ⟨k, hk⟩ := hodd
      omega
  rw [hzero]

end LRC14Concrete
end LonelyRunner

-- kernel-purity audit (fleet convention)
#print axioms LonelyRunner.LRC14Concrete.inBand_mirror
#print axioms LonelyRunner.LRC14Concrete.live_mirror
#print axioms LonelyRunner.LRC14Concrete.half_live_iff
#print axioms LonelyRunner.LRC14Concrete.liveCount_parity
#print axioms LonelyRunner.LRC14Concrete.liveCount_even_of_even_speed
