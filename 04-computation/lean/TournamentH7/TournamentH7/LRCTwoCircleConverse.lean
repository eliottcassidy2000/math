/-
  TournamentH7.LRCTwoCircleConverse — THE TWO-CIRCLE DEEP CERTIFICATE,
  CONVERSE PART II, THE SPEED-ONE AND SPEED-TWO CASES.

  `LRCTwoCircle` proves that both resonance circles lie in the canonical
  depth-six set.  This module starts the converse case split and closes its
  first two genuine branches.  If the canonical phase is deep and speed one
  fails, then the phase lies on the integer circle.  If speed one is safe and
  speed two fails, then the phase lies on the half circle.

  The proof keeps the witnesses.  Among six distinct failing speeds one has
  speed at least six.  The coefficient-weight-14 relation lock forces that
  runner's witness to be the corresponding multiple of the speed-one
  witness.  The latter is either zero or one, and the large runner's strict
  band inequality sharpens the speed-one `1/14` neighborhood to the required
  `1/84` resonance circle.

  On the speed-two branch the witness is forced to one; odd speeds 3,...,11
  are excluded by relation locking, while a failing speed thirteen is
  incompatible with every even speed 4,...,12 by a one-unit branch jump.
  Depth six therefore forces speed twelve to fail with witness six.

  These are exactly the `k₀ = 1` and `k₀ = 2` branches of the structured Part
  II plan.  The finite endpoint `k₀ >= 9` is also excluded.  This module does
  not claim the remaining `k₀ = 3,...,8` cases or the full converse.

  Kernel-pure: no `sorry`, no `native_decide`, no new axioms.
-/
import Mathlib
import TournamentH7.LRCTwoCircle

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- If speed one and one speed `k >= 6` both fail, relation locking places
their witnesses on the same integer ray and the phase lies on Circle I. -/
theorem circleI_of_one_and_large_witness
    (q p k : ℕ) (wOne wLarge : ℤ)
    (hq : 0 < q) (hp : 0 < p) (hpq : p < q)
    (hk6 : 6 ≤ k) (hk13 : k ≤ 13)
    (hOne : 14 * |(p : ℤ) - wOne * q| < q)
    (hLarge : 14 * |(k : ℤ) * (p : ℤ) - wLarge * q| < q) :
    84 * (p : ℤ) < q ∨ 84 * ((q : ℤ) - p) < q := by
  have hk1 : 1 ≤ k := by omega
  have hlock : (k : ℤ) * wOne = wLarge := by
    have h := rational_lock_weight14
      1 (k : ℤ) wOne wLarge 1 k q p
      (by omega) hk1 (by omega) (by norm_num) (by simpa using hOne) hLarge
    simpa using h
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  have hpZ : (0 : ℤ) < (p : ℤ) := by exact_mod_cast hp
  have hpqZ : (p : ℤ) < q := by exact_mod_cast hpq
  have hwLower : 0 ≤ wOne := by
    by_contra hneg
    have hw : wOne ≤ -1 := by omega
    have hdist : (q : ℤ) ≤ |(p : ℤ) - wOne * q| := by
      have hpos : (0 : ℤ) < (p : ℤ) - wOne * q := by nlinarith
      rw [abs_of_pos hpos]
      nlinarith
    nlinarith [hOne]
  have hwUpper : wOne ≤ 1 := by
    by_contra hhigh
    have hw : 2 ≤ wOne := by omega
    have hdist : (q : ℤ) ≤ |(p : ℤ) - wOne * q| := by
      have hneg : (p : ℤ) - wOne * q < 0 := by nlinarith
      rw [abs_of_neg hneg]
      nlinarith
    nlinarith [hOne]
  rcases (show wOne = 0 ∨ wOne = 1 by omega) with rfl | rfl
  · left
    have hwLarge : wLarge = 0 := by simpa using hlock.symm
    rw [hwLarge, zero_mul, sub_zero, abs_of_pos (mul_pos (by positivity) hpZ)] at hLarge
    have hkp : 84 * (p : ℤ) ≤ 14 * (k : ℤ) * (p : ℤ) := by nlinarith
    exact lt_of_le_of_lt hkp hLarge
  · right
    have hwLarge : wLarge = (k : ℤ) := by simpa using hlock.symm
    rw [hwLarge] at hLarge
    have hrewrite : (k : ℤ) * (p : ℤ) - (k : ℤ) * (q : ℤ)
        = -((k : ℤ) * ((q : ℤ) - p)) := by ring
    rw [hrewrite, abs_neg, abs_of_pos (mul_pos (by positivity) (by linarith))] at hLarge
    have hkqp : 84 * ((q : ℤ) - p)
        ≤ 14 * (k : ℤ) * ((q : ℤ) - p) := by nlinarith
    exact lt_of_le_of_lt hkqp hLarge

/-- **Part II, `k₀ = 1`.**  A canonical depth-six phase at which speed one
fails lies on the integer resonance circle. -/
theorem circleI_of_deep_and_speedOne_bad
    (q p : ℕ) (hq : 0 < q) (hp : 0 < p) (hpq : p < q)
    (hdeep : 6 ≤ bandCount canonical q p)
    (hOneBad : ¬ inBand canonical q p (0 : Fin 13)) :
    84 * (p : ℤ) < q ∨ 84 * ((q : ℤ) - p) < q := by
  let bad : Finset (Fin 13) :=
    Finset.univ.filter fun i : Fin 13 => ¬ inBand canonical q p i
  have hbadCard : 6 ≤ bad.card := by
    simpa [bad, bandCount] using hdeep
  have hlarge : ∃ i : Fin 13, i ∈ bad ∧ 5 ≤ (i : ℕ) := by
    by_contra hnone
    push_neg at hnone
    have hsub : bad ⊆ Finset.univ.filter fun i : Fin 13 => (i : ℕ) < 5 := by
      intro i hi
      rw [Finset.mem_filter]
      exact ⟨Finset.mem_univ i, by omega⟩
    have hsmallCard :
        (Finset.univ.filter fun i : Fin 13 => (i : ℕ) < 5).card = 5 := by
      decide
    have hcardLe := Finset.card_le_card hsub
    rw [hsmallCard] at hcardLe
    omega
  obtain ⟨i, hiBad, hi5⟩ := hlarge
  have hiBad' : ¬ inBand canonical q p i := by
    exact (Finset.mem_filter.mp hiBad).2
  let k : ℕ := (i : ℕ) + 1
  have hk6 : 6 ≤ k := by omega
  have hk13 : k ≤ 13 := by
    exact Nat.succ_le_iff.mpr i.isLt
  have hOne := bad_at_witness canonical q p (0 : Fin 13) hq hOneBad
  have hLarge := bad_at_witness canonical q p i hq hiBad'
  have hOne' :
      14 * |(p : ℤ) - failWitness canonical q p (0 : Fin 13) * q| < q := by
    simpa [canonical] using hOne
  have hLarge' :
      14 * |(k : ℤ) * (p : ℤ) - failWitness canonical q p i * q| < q := by
    simpa [canonical, k] using hLarge
  exact circleI_of_one_and_large_witness q p k
    (failWitness canonical q p (0 : Fin 13))
    (failWitness canonical q p i) hq hp hpq hk6 hk13 hOne' hLarge'

/-- A depth-six canonical phase has a failing speed among `1,...,8`.  This is
the finite endpoint of the Part II smallest-failing-speed split: `k₀ >= 9`
cannot occur because only five canonical speeds remain. -/
theorem exists_bad_speed_le_eight_of_deep
    (q p : ℕ) (hdeep : 6 ≤ bandCount canonical q p) :
    ∃ i : Fin 13, ¬ inBand canonical q p i ∧ (i : ℕ) < 8 := by
  let bad : Finset (Fin 13) :=
    Finset.univ.filter fun i : Fin 13 => ¬ inBand canonical q p i
  have hbadCard : 6 ≤ bad.card := by
    simpa [bad, bandCount] using hdeep
  by_contra hnone
  push_neg at hnone
  have hsub : bad ⊆ Finset.univ.filter fun i : Fin 13 => 8 ≤ (i : ℕ) := by
    intro i hi
    rw [Finset.mem_filter]
    exact ⟨Finset.mem_univ i, hnone i (Finset.mem_filter.mp hi).2⟩
  have hhighCard :
      (Finset.univ.filter fun i : Fin 13 => 8 ≤ (i : ℕ)).card = 5 := by
    decide
  have hcardLe := Finset.card_le_card hsub
  rw [hhighCard] at hcardLe
  omega

/-- If speed one is safe while speed two fails, the speed-two witness is the
central integer `1`; the endpoint witnesses `0` and `2` would already make
speed one fail. -/
theorem speedTwo_witness_eq_one_of_speedOne_safe
    (q p : ℕ) (wTwo : ℤ) (hq : 0 < q) (hp : 0 < p) (hpq : p < q)
    (hOneSafe : inBand canonical q p (0 : Fin 13))
    (hTwo : 14 * |2 * (p : ℤ) - wTwo * q| < q) :
    wTwo = 1 := by
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  have hpZ : (0 : ℤ) < (p : ℤ) := by exact_mod_cast hp
  have hpqZ : (p : ℤ) < q := by exact_mod_cast hpq
  have hwLower : 0 ≤ wTwo := by
    by_contra hneg
    have hw : wTwo ≤ -1 := by omega
    have hdist : (q : ℤ) ≤ |2 * (p : ℤ) - wTwo * q| := by
      have hpos : (0 : ℤ) < 2 * (p : ℤ) - wTwo * q := by nlinarith
      rw [abs_of_pos hpos]
      nlinarith
    nlinarith [hTwo]
  have hwUpper : wTwo ≤ 2 := by
    by_contra hhigh
    have hw : 3 ≤ wTwo := by omega
    have hdist : (q : ℤ) ≤ |2 * (p : ℤ) - wTwo * q| := by
      have hneg : 2 * (p : ℤ) - wTwo * q < 0 := by nlinarith
      rw [abs_of_neg hneg]
      nlinarith
    nlinarith [hTwo]
  rcases (show wTwo = 0 ∨ wTwo = 1 ∨ wTwo = 2 by omega) with rfl | rfl | rfl
  · exfalso
    apply (not_inBand_of_witness canonical q p (0 : Fin 13) hq 0) hOneSafe
    have hsmall : 14 * (p : ℤ) < q := by
      rw [zero_mul, sub_zero, abs_of_pos (by positivity : (0 : ℤ) < 2 * (p : ℤ))] at hTwo
      nlinarith
    simpa [canonical, abs_of_pos hpZ] using hsmall
  · rfl
  · exfalso
    apply (not_inBand_of_witness canonical q p (0 : Fin 13) hq 1) hOneSafe
    have hsmall : 14 * ((q : ℤ) - p) < q := by
      have hrewrite : 2 * (p : ℤ) - 2 * (q : ℤ) = -2 * ((q : ℤ) - p) := by ring
      rw [hrewrite, abs_mul, abs_neg, abs_of_pos (by linarith : (0 : ℤ) < (q : ℤ) - p)] at hTwo
      norm_num at hTwo
      nlinarith
    have hrewrite : (p : ℤ) - 1 * q = -((q : ℤ) - p) := by ring
    rw [canonical, hrewrite, abs_neg, abs_of_pos (by linarith : (0 : ℤ) < (q : ℤ) - p)]
    exact hsmall

/-- With central witness `w₂=1`, no odd speed from three through eleven can
fail: the light ratio relation would make an odd integer twice an integer. -/
theorem odd_speed_not_bad_of_speedTwo_witness_one
    (q p k : ℕ) (wOdd : ℤ)
    (hk3 : 3 ≤ k) (hk11 : k ≤ 11) (hkOdd : k % 2 = 1)
    (hTwo : 14 * |2 * (p : ℤ) - q| < q)
    (hOdd : 14 * |(k : ℤ) * (p : ℤ) - wOdd * q| < q) : False := by
  have hlock := rational_lock_weight14
    2 (k : ℤ) 1 wOdd 2 k q p
    (by omega) (by omega) (by omega) (by norm_num) hTwo hOdd
  have hparity : k % 2 = 0 := by
    have hkEq : (k : ℤ) = 2 * wOdd := by simpa using hlock
    omega
  omega

/-- A failing speed thirteen on the central speed-two branch excludes every
even speed at least four.  The pair cross-bound forces the `13` witness onto
one of the two adjacent branches `2w₁₃-13=+-1`; the resulting exact `q` jump
is too large to be paid by the two strict `1/14` errors. -/
theorem speedThirteen_incompatible_with_even_ge_four
    (q p r : ℕ) (wThirteen wEven : ℤ)
    (hq : 0 < q) (hr2 : 2 ≤ r) (hr6 : r ≤ 6)
    (hTwo : 14 * |2 * (p : ℤ) - q| < q)
    (hThirteen : 14 * |13 * (p : ℤ) - wThirteen * q| < q)
    (hEven : 14 * |(2 * r : ℕ) * (p : ℤ) - wEven * q| < q) : False := by
  have hcross := witness_cross_bound
    2 13 1 wThirteen q p hq (by norm_num) (by norm_num) hTwo hThirteen
  have habsBranch : |2 * wThirteen - 13| = 1 := by
    have hbound : |2 * wThirteen - 13| ≤ 1 := by
      rw [abs_sub_comm]
      simpa using (show |1 * 13 - wThirteen * 2| ≤ 1 by omega)
    have hnonzero : 2 * wThirteen - 13 ≠ 0 := by omega
    have hone := Int.one_le_abs hnonzero
    omega
  have hEvenLock := rational_lock_weight14
    2 (2 * (r : ℤ)) 1 wEven 1 r q p
    (by omega) (by omega) (by omega) (by ring) hTwo
    (by simpa only [Nat.cast_mul, Nat.cast_ofNat] using hEven)
  have hwEven : wEven = (r : ℤ) := by
    have : (r : ℤ) = wEven := by simpa using hEvenLock
    exact this.symm
  have hEvenError :
      14 * (r : ℤ) * |2 * (p : ℤ) - q| < q := by
    rw [hwEven] at hEven
    have hfactor : ((2 * r : ℕ) : ℤ) * (p : ℤ) - (r : ℤ) * q
        = (r : ℤ) * (2 * (p : ℤ) - q) := by
      push_cast
      ring
    rw [hfactor, abs_mul, abs_of_pos (by exact_mod_cast (by omega : 0 < r))] at hEven
    nlinarith
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  have hkey :
      13 * (2 * (p : ℤ) - q) - 2 * (13 * (p : ℤ) - wThirteen * q)
        = (2 * wThirteen - 13) * q := by ring
  have htriangle :
      (q : ℤ) ≤ 13 * |2 * (p : ℤ) - q|
        + 2 * |13 * (p : ℤ) - wThirteen * q| := by
    have h := abs_sub (13 * (2 * (p : ℤ) - q))
      (2 * (13 * (p : ℤ) - wThirteen * q))
    rw [abs_mul, abs_of_pos (by norm_num : (0 : ℤ) < 13),
      abs_mul, abs_of_pos (by norm_num : (0 : ℤ) < 2), hkey,
      abs_mul, habsBranch, abs_of_pos hqZ, one_mul] at h
    exact h
  have hTwoSharp : 28 * |2 * (p : ℤ) - q| < q := by
    have habs : (0 : ℤ) ≤ |2 * (p : ℤ) - q| := abs_nonneg _
    have hrZ : (2 : ℤ) ≤ r := by exact_mod_cast hr2
    nlinarith [hEvenError]
  nlinarith [htriangle, hTwoSharp, hThirteen,
    abs_nonneg (2 * (p : ℤ) - q),
    abs_nonneg (13 * (p : ℤ) - wThirteen * q)]

/-- **Part II, `k₀ = 2`.**  If the canonical phase is deep, speed one is
safe, and speed two fails, then all six even speeds must fail coherently and
the phase lies on the half-resonance circle. -/
theorem circleII_of_deep_speedOne_safe_speedTwo_bad
    (q p : ℕ) (hq : 0 < q) (hp : 0 < p) (hpq : p < q)
    (hdeep : 6 ≤ bandCount canonical q p)
    (hOneSafe : inBand canonical q p (0 : Fin 13))
    (hTwoBad : ¬ inBand canonical q p (1 : Fin 13)) :
    84 * |2 * (p : ℤ) - q| < q := by
  let bad : Finset (Fin 13) :=
    Finset.univ.filter fun i : Fin 13 => ¬ inBand canonical q p i
  have hbadCard : 6 ≤ bad.card := by
    simpa [bad, bandCount] using hdeep
  let wTwo := failWitness canonical q p (1 : Fin 13)
  have hTwoRaw := bad_at_witness canonical q p (1 : Fin 13) hq hTwoBad
  have hTwo : 14 * |2 * (p : ℤ) - wTwo * q| < q := by
    simpa [canonical, wTwo] using hTwoRaw
  have hwTwo : wTwo = 1 :=
    speedTwo_witness_eq_one_of_speedOne_safe q p wTwo hq hp hpq hOneSafe hTwo
  have hTwoOne : 14 * |2 * (p : ℤ) - q| < q := by simpa [hwTwo] using hTwo
  have hOddSafe : ∀ i : Fin 13, 2 ≤ (i : ℕ) → (i : ℕ) ≤ 10 →
      (i : ℕ) % 2 = 0 → inBand canonical q p i := by
    intro i hi2 hi10 hiEven
    by_contra hiBad
    have hiRaw := bad_at_witness canonical q p i hq hiBad
    have hiBound :
        14 * |(((i : ℕ) + 1 : ℕ) : ℤ) * (p : ℤ)
          - failWitness canonical q p i * q| < q := by
      simpa [canonical] using hiRaw
    exact odd_speed_not_bad_of_speedTwo_witness_one
      q p ((i : ℕ) + 1) (failWitness canonical q p i)
      (by omega) (by omega) (by omega) hTwoOne hiBound
  have hThirteenSafe : inBand canonical q p (12 : Fin 13) := by
    by_contra hThirteenBad
    have hThirteenRaw :=
      bad_at_witness canonical q p (12 : Fin 13) hq hThirteenBad
    have hThirteen :
        14 * |13 * (p : ℤ)
          - failWitness canonical q p (12 : Fin 13) * q| < q := by
      simpa [canonical] using hThirteenRaw
    have hEvenSafe : ∀ i : Fin 13, 3 ≤ (i : ℕ) → (i : ℕ) ≤ 11 →
        (i : ℕ) % 2 = 1 → inBand canonical q p i := by
      intro i hi3 hi11 hiOdd
      by_contra hiBad
      have hiRaw := bad_at_witness canonical q p i hq hiBad
      let r : ℕ := ((i : ℕ) + 1) / 2
      have hr2 : 2 ≤ r := by omega
      have hr6 : r ≤ 6 := by omega
      have hspeed : (i : ℕ) + 1 = 2 * r := by omega
      have hiBound :
          14 * |(2 * r : ℕ) * (p : ℤ)
            - failWitness canonical q p i * q| < q := by
        simpa [canonical, hspeed] using hiRaw
      exact speedThirteen_incompatible_with_even_ge_four
        q p r (failWitness canonical q p (12 : Fin 13))
        (failWitness canonical q p i) hq hr2 hr6 hTwoOne hThirteen hiBound
    have hsub : bad ⊆
        Finset.univ.filter fun i : Fin 13 => (i : ℕ) = 1 ∨ (i : ℕ) = 12 := by
      intro i hi
      have hiBad : ¬ inBand canonical q p i := (Finset.mem_filter.mp hi).2
      rw [Finset.mem_filter]
      refine ⟨Finset.mem_univ i, ?_⟩
      by_cases hallowed : (i : ℕ) = 1 ∨ (i : ℕ) = 12
      · exact hallowed
      · exfalso
        by_cases hiParity : (i : ℕ) % 2 = 0
        · by_cases hiZero : (i : ℕ) = 0
          · have : i = (0 : Fin 13) := Fin.ext hiZero
            exact hiBad (this ▸ hOneSafe)
          · exact hiBad (hOddSafe i (by omega) (by omega) hiParity)
        · have hiOdd : (i : ℕ) % 2 = 1 := by omega
          exact hiBad (hEvenSafe i (by omega) (by omega) hiOdd)
    have htwoCard :
        (Finset.univ.filter fun i : Fin 13 =>
          (i : ℕ) = 1 ∨ (i : ℕ) = 12).card = 2 := by
      decide
    have hcardLe := Finset.card_le_card hsub
    rw [htwoCard] at hcardLe
    omega
  have hTwelveBad : ¬ inBand canonical q p (11 : Fin 13) := by
    intro hTwelveSafe
    have hsub : bad ⊆ Finset.univ.filter fun i : Fin 13 =>
        (i : ℕ) % 2 = 1 ∧ (i : ℕ) < 10 := by
      intro i hi
      have hiBad : ¬ inBand canonical q p i := (Finset.mem_filter.mp hi).2
      rw [Finset.mem_filter]
      refine ⟨Finset.mem_univ i, ?_⟩
      by_cases hiZero : (i : ℕ) = 0
      · have : i = (0 : Fin 13) := Fin.ext hiZero
        exact False.elim (hiBad (this ▸ hOneSafe))
      by_cases hiTwelve : (i : ℕ) = 11
      · have : i = (11 : Fin 13) := Fin.ext hiTwelve
        exact False.elim (hiBad (this ▸ hTwelveSafe))
      by_cases hiThirteen : (i : ℕ) = 12
      · have : i = (12 : Fin 13) := Fin.ext hiThirteen
        exact False.elim (hiBad (this ▸ hThirteenSafe))
      have hiOdd : (i : ℕ) % 2 = 1 := by
        by_contra hiNotOdd
        have hiEven : (i : ℕ) % 2 = 0 := by omega
        exact hiBad (hOddSafe i (by omega) (by omega) hiEven)
      exact ⟨hiOdd, by omega⟩
    have hfiveCard :
        (Finset.univ.filter fun i : Fin 13 =>
          (i : ℕ) % 2 = 1 ∧ (i : ℕ) < 10).card = 5 := by
      decide
    have hcardLe := Finset.card_le_card hsub
    rw [hfiveCard] at hcardLe
    omega
  have hTwelveRaw :=
    bad_at_witness canonical q p (11 : Fin 13) hq hTwelveBad
  have hTwelve :
      14 * |12 * (p : ℤ)
        - failWitness canonical q p (11 : Fin 13) * q| < q := by
    simpa [canonical] using hTwelveRaw
  have hTwelveLock := rational_lock_weight14
    2 12 1 (failWitness canonical q p (11 : Fin 13)) 1 6 q p
    (by omega) (by omega) (by omega) (by norm_num) hTwoOne hTwelve
  have hwTwelve : failWitness canonical q p (11 : Fin 13) = 6 := by
    have : (6 : ℤ) = failWitness canonical q p (11 : Fin 13) := by
      simpa using hTwelveLock
    exact this.symm
  rw [hwTwelve] at hTwelve
  have hfactor : 12 * (p : ℤ) - 6 * (q : ℤ)
      = 6 * (2 * (p : ℤ) - q) := by ring
  rw [hfactor, abs_mul, abs_of_pos (by norm_num : (0 : ℤ) < 6)] at hTwelve
  nlinarith

/-! ## Axiom audit -/
#print axioms circleI_of_one_and_large_witness
#print axioms circleI_of_deep_and_speedOne_bad
#print axioms exists_bad_speed_le_eight_of_deep
#print axioms speedTwo_witness_eq_one_of_speedOne_safe
#print axioms odd_speed_not_bad_of_speedTwo_witness_one
#print axioms speedThirteen_incompatible_with_even_ge_four
#print axioms circleII_of_deep_speedOne_safe_speedTwo_bad

end LRC14Concrete
end LonelyRunner
