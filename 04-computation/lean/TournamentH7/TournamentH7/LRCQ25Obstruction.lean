/-
  TournamentH7.LRCQ25Obstruction — exact formal obstruction to a universal
  denominator-q <= 25 LRC(14) witness (codex-S58, THM-764).

  The explicit covering family

      S0 = {43,55,61,70,73,79,83,99,103,104,109,113,156}

  has no weak safe-band multiplier at any q in [15,25].  Here "weak" means
  the actual LRC threshold, with inclusive inequalities

      q <= 14 * ((v_i p) mod q) <= 13q.

  Thus this is stronger than merely ruling out the strict-live predicate.
  It does NOT refute LRC(14): p/q = 2/27 is an explicit witness.  The file
  also formalizes the infinite common-translation orbit by
  L = lcm(2,...,25) = 26771144400.  Translation by kL preserves every
  residue through q=25, so every translated family is still covering and
  still has no q-witness in [15,25].

  All finite certificates use plain `decide`, not `native_decide`; there are
  no `sorry`s and no new axioms.
-/
import TournamentH7.LRCB5Periodic
import TournamentH7.LRC14WindowWiring

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- The codex-S58 13-speed obstruction row. -/
def q25ObstructionSpeeds : Fin 13 → ℤ :=
  ![43, 55, 61, 70, 73, 79, 83, 99, 103, 104, 109, 113, 156]

/-- A denominator-`q` rational witness in the exact inclusive LRC(14) band. -/
def HasQWitness (v : Fin 13 → ℤ) (q : ℕ) : Prop :=
  ∃ p ∈ Finset.Ioo 0 q, ∀ i, inBand v q p i

/-- Vanishing `bandCount` is exactly the pointwise safe-band condition. -/
theorem bandCount_eq_zero_iff_all_inBand (v : Fin 13 → ℤ) (q p : ℕ) :
    bandCount v q p = 0 ↔ ∀ i, inBand v q p i := by
  unfold bandCount
  rw [Finset.card_eq_zero, Finset.filter_eq_empty_iff]
  simp

/-- `HasQWitness` is equivalently positivity of the exact live-multiplier count. -/
theorem hasQWitness_iff_liveCount_pos (v : Fin 13 → ℤ) (q : ℕ) :
    HasQWitness v q ↔ 0 < liveCount v q := by
  constructor
  · rintro ⟨p, hp, hband⟩
    unfold liveCount
    apply Finset.card_pos.mpr
    exact ⟨p, Finset.mem_filter.mpr
      ⟨hp, (bandCount_eq_zero_iff_all_inBand v q p).2 hband⟩⟩
  · intro hlive
    unfold liveCount at hlive
    obtain ⟨p, hp⟩ := Finset.card_pos.mp hlive
    exact ⟨p, (Finset.mem_filter.mp hp).1,
      (bandCount_eq_zero_iff_all_inBand v q p).1 (Finset.mem_filter.mp hp).2⟩

/-- A zero owner blocks every multiplier at a positive modulus, even for the
inclusive (weak) safe band. -/
theorem not_hasQWitness_of_dvd (v : Fin 13 → ℤ) (q : ℕ) (hq : 0 < q)
    (i : Fin 13) (hi : (q : ℤ) ∣ v i) : ¬ HasQWitness v q := by
  rintro ⟨p, _hp, hband⟩
  have hz : (v i * (p : ℤ)) % (q : ℤ) = 0 :=
    Int.emod_eq_zero_of_dvd (hi.mul_right (p : ℤ))
  have hlo := (hband i).1
  rw [hz] at hlo
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  omega

/-- The row is a genuine covering family: every q in [2,14] divides a speed. -/
theorem q25Obstruction_covering : LRC14.CoveringFamily q25ObstructionSpeeds := by
  intro q hq2 hq14
  interval_cases q
  · exact ⟨3, by decide⟩
  · exact ⟨7, by decide⟩
  · exact ⟨9, by decide⟩
  · exact ⟨1, by decide⟩
  · exact ⟨12, by decide⟩
  · exact ⟨3, by decide⟩
  · exact ⟨9, by decide⟩
  · exact ⟨7, by decide⟩
  · exact ⟨3, by decide⟩
  · exact ⟨1, by decide⟩
  · exact ⟨12, by decide⟩
  · exact ⟨9, by decide⟩
  · exact ⟨3, by decide⟩

/-- No modulus q in [15,25] has a zero residue in the obstruction row. -/
theorem q25Obstruction_zero_free (q : ℕ) (hq15 : 15 ≤ q) (hq25 : q ≤ 25) :
    ∀ i, q25ObstructionSpeeds i % (q : ℤ) ≠ 0 := by
  interval_cases q <;> decide

/-- Exact finite checksum: the live-multiplier count is zero for 15 <= q <= 25. -/
theorem q25Obstruction_liveCount_zero (q : ℕ) (hq15 : 15 ≤ q) (hq25 : q ≤ 25) :
    liveCount q25ObstructionSpeeds q = 0 := by
  interval_cases q <;> decide

/-- **The q <= 25 obstruction.**  No rational p/q at the actual inclusive
LRC threshold witnesses this covering, zero-free family for 15 <= q <= 25. -/
theorem q25Obstruction_no_qWitness (q : ℕ) (hq15 : 15 ≤ q) (hq25 : q ≤ 25) :
    ¬ HasQWitness q25ObstructionSpeeds q := by
  rw [hasQWitness_iff_liveCount_pos,
    q25Obstruction_liveCount_zero q hq15 hq25]
  simp

/-- The same row in fact has no denominator-q witness anywhere in the full
range 2 <= q <= 25.  Covering supplies a zero owner through q=14, and the
exact S58 checksum handles q=15,...,25. -/
theorem q25Obstruction_no_qWitness_through25
    (q : ℕ) (hq2 : 2 ≤ q) (hq25 : q ≤ 25) :
    ¬ HasQWitness q25ObstructionSpeeds q := by
  by_cases hq14 : q ≤ 14
  · obtain ⟨i, hi⟩ := q25Obstruction_covering q hq2 hq14
    exact not_hasQWitness_of_dvd q25ObstructionSpeeds q (by omega) i hi
  · exact q25Obstruction_no_qWitness q (by omega) hq25

/-! ## Honest upper edge: q = 27 works -/

/-- The multiplier p=2 at q=27 clears every runner in the inclusive band. -/
theorem q25Obstruction_q27_multiplier_two :
    ∀ i, inBand q25ObstructionSpeeds 27 2 i := by
  decide

/-- Hence the obstruction is bounded-denominator only: it has a q=27 witness. -/
theorem q25Obstruction_hasQWitness_27 : HasQWitness q25ObstructionSpeeds 27 :=
  ⟨2, by decide, q25Obstruction_q27_multiplier_two⟩

/-- The q=27 certificate reaches the genuine LRC(14) floor. -/
theorem q25Obstruction_mreach_ge :
    (1 : ℝ) / 14 ≤ Mreach q25ObstructionSpeeds :=
  mreach_ge_of_pairsum_band q25ObstructionSpeeds 2 27 (by norm_num)
    q25Obstruction_q27_multiplier_two

/-! ## The infinite common-translation obstruction orbit -/

/-- `lcm(2,...,25)`, recorded as the exact integer used by the S58 orbit. -/
def q25TranslationPeriod : ℕ := 26771144400

/-- Add a common multiple of every q <= 25 to all thirteen speeds. -/
def q25TranslatedSpeeds (k : ℕ) : Fin 13 → ℤ := fun i =>
  q25ObstructionSpeeds i + (k : ℤ) * (q25TranslationPeriod : ℤ)

/-- Every modulus in [2,25] divides the recorded translation period. -/
theorem q_dvd_q25TranslationPeriod (q : ℕ) (hq2 : 2 ≤ q) (hq25 : q ≤ 25) :
    q ∣ q25TranslationPeriod := by
  interval_cases q <;> norm_num [q25TranslationPeriod]

/-- A common translation by `kL` preserves each speed modulo q <= 25. -/
theorem q25Translated_modEq (k q : ℕ) (hq2 : 2 ≤ q) (hq25 : q ≤ 25) (i : Fin 13) :
    q25TranslatedSpeeds k i ≡ q25ObstructionSpeeds i [ZMOD (q : ℤ)] := by
  rw [Int.modEq_iff_dvd]
  have hqN : q ∣ q25TranslationPeriod :=
    q_dvd_q25TranslationPeriod q hq2 hq25
  have hqZ : (q : ℤ) ∣ (q25TranslationPeriod : ℤ) :=
    Int.natCast_dvd_natCast.mpr hqN
  have hmul : (q : ℤ) ∣ (k : ℤ) * (q25TranslationPeriod : ℤ) :=
    hqZ.mul_left (k : ℤ)
  simpa [q25TranslatedSpeeds] using hmul.neg_right

/-- Consequently every band count is identical through q=25. -/
theorem q25Translated_bandCount (k q p : ℕ) (hq2 : 2 ≤ q) (hq25 : q ≤ 25) :
    bandCount (q25TranslatedSpeeds k) q p = bandCount q25ObstructionSpeeds q p := by
  apply bandCount_congr_mod
  intro i
  exact (q25Translated_modEq k q hq2 hq25 i).mul_right (p : ℤ)

/-- The entire live-multiplier count is translation-invariant through q=25. -/
theorem q25Translated_liveCount (k q : ℕ) (hq2 : 2 ≤ q) (hq25 : q ≤ 25) :
    liveCount (q25TranslatedSpeeds k) q = liveCount q25ObstructionSpeeds q := by
  unfold liveCount
  apply congrArg Finset.card
  apply Finset.filter_congr
  intro p _
  rw [q25Translated_bandCount k q p hq2 hq25]

/-- The orbit is uniformly zero-free at every q in [15,25]. -/
theorem q25Translated_zero_free (k q : ℕ) (hq15 : 15 ≤ q) (hq25 : q ≤ 25) :
    ∀ i, q25TranslatedSpeeds k i % (q : ℤ) ≠ 0 := by
  intro i
  have hmod := q25Translated_modEq k q (by omega) hq25 i
  rw [hmod.eq]
  exact q25Obstruction_zero_free q hq15 hq25 i

/-- Every translated row remains covering. -/
theorem q25Translated_covering (k : ℕ) :
    LRC14.CoveringFamily (q25TranslatedSpeeds k) := by
  intro q hq2 hq14
  obtain ⟨i, hi⟩ := q25Obstruction_covering q hq2 hq14
  refine ⟨i, ?_⟩
  have hqN : q ∣ q25TranslationPeriod :=
    q_dvd_q25TranslationPeriod q hq2 (by omega)
  have hqZ : (q : ℤ) ∣ (q25TranslationPeriod : ℤ) :=
    Int.natCast_dvd_natCast.mpr hqN
  have hshift : (q : ℤ) ∣ (k : ℤ) * (q25TranslationPeriod : ℤ) :=
    hqZ.mul_left (k : ℤ)
  simpa [q25TranslatedSpeeds] using hi.add hshift

/-- Every translated row is primitive.  The common shift cancels from the
differences `70-43=27` and `83-43=40`, and `3*27-2*40=1`. -/
theorem q25Translated_tupleGcd_one (k : ℕ) :
    LRC14.tupleGcd (q25TranslatedSpeeds k) = 1 := by
  have h27 : (LRC14.tupleGcd (q25TranslatedSpeeds k) : ℤ) ∣ 27 := by
    have h := (LRC14.tupleGcd_dvd (q25TranslatedSpeeds k) (3 : Fin 13)).sub
      (LRC14.tupleGcd_dvd (q25TranslatedSpeeds k) (0 : Fin 13))
    simpa [q25TranslatedSpeeds, q25ObstructionSpeeds] using h
  have h40 : (LRC14.tupleGcd (q25TranslatedSpeeds k) : ℤ) ∣ 40 := by
    have h := (LRC14.tupleGcd_dvd (q25TranslatedSpeeds k) (6 : Fin 13)).sub
      (LRC14.tupleGcd_dvd (q25TranslatedSpeeds k) (0 : Fin 13))
    simpa [q25TranslatedSpeeds, q25ObstructionSpeeds] using h
  have h1 : (LRC14.tupleGcd (q25TranslatedSpeeds k) : ℤ) ∣ 1 := by
    have h := (h27.mul_left (3 : ℤ)).sub (h40.mul_left (2 : ℤ))
    norm_num at h ⊢
    exact h
  exact Nat.dvd_one.mp (Int.natCast_dvd_natCast.mp h1)

/-- In particular the base obstruction row itself is primitive. -/
theorem q25Obstruction_tupleGcd_one :
    LRC14.tupleGcd q25ObstructionSpeeds = 1 := by
  simpa [q25TranslatedSpeeds] using q25Translated_tupleGcd_one 0

/-- **Infinite uniform obstruction.**  For every k, the translated covering
family has no denominator-q witness for any 15 <= q <= 25. -/
theorem q25Translated_no_qWitness (k q : ℕ) (hq15 : 15 ≤ q) (hq25 : q ≤ 25) :
    ¬ HasQWitness (q25TranslatedSpeeds k) q := by
  rw [hasQWitness_iff_liveCount_pos,
    q25Translated_liveCount k q (by omega) hq25,
    q25Obstruction_liveCount_zero q hq15 hq25]
  simp

/-- Uniform full-window form for the infinite orbit: no translated primitive
covering row has a denominator-q witness for any 2 <= q <= 25. -/
theorem q25Translated_no_qWitness_through25
    (k q : ℕ) (hq2 : 2 ≤ q) (hq25 : q ≤ 25) :
    ¬ HasQWitness (q25TranslatedSpeeds k) q := by
  by_cases hq14 : q ≤ 14
  · obtain ⟨i, hi⟩ := q25Translated_covering k q hq2 hq14
    exact not_hasQWitness_of_dvd (q25TranslatedSpeeds k) q (by omega) i hi
  · exact q25Translated_no_qWitness k q (by omega) hq25

/-! ## Axiom audit -/
#print axioms bandCount_eq_zero_iff_all_inBand
#print axioms q25Obstruction_covering
#print axioms q25Obstruction_no_qWitness
#print axioms q25Obstruction_no_qWitness_through25
#print axioms q25Obstruction_hasQWitness_27
#print axioms q25Obstruction_mreach_ge
#print axioms q25Translated_covering
#print axioms q25Translated_tupleGcd_one
#print axioms q25Translated_zero_free
#print axioms q25Translated_no_qWitness
#print axioms q25Translated_no_qWitness_through25

end LRC14Concrete
end LonelyRunner
