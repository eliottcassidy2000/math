/-
  TournamentH7.LRCPrimitiveCarrierINV — the primitive q=14 carrier inverse
  and its global LRC(14) consumer (codex-2026-07-18-S67).

  Literal `INVcov` is false because Covering is not invariant under division by
  the common gcd: the doubled AP covers q=14, while its primitive core does not.
  The correct normalized object separates the sole divisor not already forced by
  absence of a 1/13-lonely time.

  If a family has no `Lonely 13` time, the denominator sieve forces a carrier for
  every q=2,...,13.  Thus full Covering is equivalent to the single remaining
  predicate `HasFourteenCarrier`.  This yields the primitive trichotomy

      Lonely13  OR  no q=14 carrier  OR  thirteen-fold dominance.

  Each branch has a proved LRC(14) exit: band monotonicity, the exact time 1/14,
  or `ap_core_bridge`.  Dividing by `tupleGcd` before applying the trichotomy and
  transporting the resulting time back by dilation proves the global consumer.

  The structural supplier `PrimitiveCarrierINV` remains an explicit hypothesis.
  Everything connecting it to LRC(14) below is kernel-checked and sorry-free.
-/
import Mathlib
import TournamentH7.LRCDilation
import TournamentH7.LRC14WindowWiring
import TournamentH7.LRCAPCoreBridge
import TournamentH7.LRCSieveDispatch

namespace LonelyRunner
namespace PrimitiveCarrierInverse

/-- The one divisor obligation not forced by absence of a `Lonely 13` time. -/
def HasFourteenCarrier (v : Fin 13 → ℤ) : Prop :=
  ∃ i, (14 : ℤ) ∣ v i

/-- Some speed dominates all the other speeds by a factor of thirteen. -/
def ThirteenDominance (v : Fin 13 → ℤ) : Prop :=
  ∃ vstar : Fin 13, ∀ i, i ≠ vstar → 13 * v i ≤ v vstar

/-- The honest primitive version of the historical Covering inverse. -/
def PrimitiveINVcov : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, 0 < v i) → LRC14.tupleGcd v = 1 →
    Covering v → (¬ ∃ t : ℝ, Lonely 13 v t) → ThirteenDominance v

/-- The q=14-carrier form of the primitive inverse supplier.  Under
`¬ ∃ t, Lonely 13 v t`, carriers for q=2,...,13 are automatic, so this is
equivalent to `PrimitiveINVcov`. -/
def PrimitiveCarrierINV : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, 0 < v i) → LRC14.tupleGcd v = 1 →
    HasFourteenCarrier v → (¬ ∃ t : ℝ, Lonely 13 v t) → ThirteenDominance v

/-- The operational three-exit statement on a primitive positive family. -/
def PrimitiveFourteenTrichotomy : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, 0 < v i) → LRC14.tupleGcd v = 1 →
    (∃ t : ℝ, Lonely 13 v t) ∨
      ¬ HasFourteenCarrier v ∨ ThirteenDominance v

/-- Under absence of a `Lonely 13` time, Covering through 14 is exactly the
single q=14-carrier bit.  The denominator sieve already forces q=2,...,13. -/
theorem covering_iff_hasFourteenCarrier_of_no_lonely13
    (v : Fin 13 → ℤ) (hno : ¬ ∃ t : ℝ, Lonely 13 v t) :
    Covering v ↔ HasFourteenCarrier v := by
  have hcex : ∀ t : ℝ, ¬ Lonely 13 v t := by
    intro t ht
    exact hno ⟨t, ht⟩
  constructor
  · intro hcover
    exact hcover 14 (by norm_num) (by norm_num)
  · intro h14 q hq2 hq14
    by_cases hq : q = 14
    · subst q
      exact h14
    · exact counterexample_needs_all_divisors 13 v hcex q hq2 (by omega)

/-- The carrier supplier and the primitive Covering supplier are the same
statement once the proved q=2,...,13 sieve information is retained. -/
theorem primitiveCarrierINV_iff_primitiveINVcov :
    PrimitiveCarrierINV ↔ PrimitiveINVcov := by
  constructor
  · intro hcarrier v hpos hprim hcover hno
    apply hcarrier v hpos hprim
    · exact (covering_iff_hasFourteenCarrier_of_no_lonely13 v hno).mp hcover
    · exact hno
  · intro hcover v hpos hprim h14 hno
    apply hcover v hpos hprim
    · exact (covering_iff_hasFourteenCarrier_of_no_lonely13 v hno).mpr h14
    · exact hno

/-- The conditional carrier inverse is exactly the three-exit trichotomy. -/
theorem primitiveCarrierINV_iff_trichotomy :
    PrimitiveCarrierINV ↔ PrimitiveFourteenTrichotomy := by
  constructor
  · intro hinv v hpos hprim
    by_cases h13 : ∃ t : ℝ, Lonely 13 v t
    · exact Or.inl h13
    · by_cases h14 : HasFourteenCarrier v
      · exact Or.inr (Or.inr (hinv v hpos hprim h14 h13))
      · exact Or.inr (Or.inl h14)
  · intro htri v hpos hprim h14 h13
    rcases htri v hpos hprim with hlonely | hmiss | hdom
    · exact False.elim (h13 hlonely)
    · exact False.elim (hmiss h14)
    · exact hdom

/-- Existence of a lonely time is unchanged by division by the tuple gcd. -/
theorem exists_lonely_primPart_iff (n : ℕ) (v : Fin 13 → ℤ)
    (hv : ∀ i, v i ≠ 0) :
    (∃ t : ℝ, Lonely n v t) ↔
      ∃ t : ℝ, Lonely n (LRC14.primPart v) t := by
  have hgpos : 0 < LRC14.tupleGcd v := LRC14.tupleGcd_pos v hv
  have hg : (LRC14.tupleGcd v : ℤ) ≠ 0 := by
    exact_mod_cast hgpos.ne'
  simpa [LRC14.primPart] using
    (exists_lonely_of_dvd n v (LRC14.tupleGcd v : ℤ) hg
      (LRC14.tupleGcd_dvd v))

/-- In the actual no-`Lonely14` branch, full Covering can safely be rederived
after gcd normalization.  This deliberately uses `14`, not merely `13`. -/
theorem covering_primPart_of_no_lonely14 (v : Fin 13 → ℤ)
    (hv : ∀ i, v i ≠ 0) (hno : ¬ ∃ t : ℝ, Lonely 14 v t) :
    Covering (LRC14.primPart v) := by
  by_contra hcover
  apply hno
  exact (exists_lonely_primPart_iff 14 v hv).mpr
    (sieve_dispatch (LRC14.primPart v) hcover)

/-- A primitive carrier inverse supplies a `Lonely 14` time for every positive
family.  Normalize first, take one of the three primitive exits, then transport
the witness back through the common gcd. -/
theorem lonely14_of_primitiveCarrierINV (cite : LRCUpTo13)
    (hinv : PrimitiveCarrierINV) (v : Fin 13 → ℤ) (hpos : ∀ i, 0 < v i) :
    ∃ t : ℝ, Lonely 14 v t := by
  have hv : ∀ i, v i ≠ 0 := fun i => (hpos i).ne'
  let u : Fin 13 → ℤ := LRC14.primPart v
  have hupos : ∀ i, 0 < u i := by
    intro i
    exact LRC14.primPart_pos v hpos i
  have huprim : LRC14.tupleGcd u = 1 := by
    exact LRC14.primPart_gcd_eq_one v hv
  have htri : PrimitiveFourteenTrichotomy :=
    primitiveCarrierINV_iff_trichotomy.mp hinv
  have hu14 : ∃ t : ℝ, Lonely 14 u t := by
    rcases htri u hupos huprim with h13 | hmiss | hdom
    · obtain ⟨t, ht⟩ := h13
      exact ⟨t, lonely14_of_lonely_le (by norm_num) (by norm_num) ht⟩
    · refine ⟨1 / 14, lonely_of_no_multiple 14 (by norm_num) u ?_⟩
      intro i hi
      exact hmiss ⟨i, hi⟩
    · obtain ⟨vstar, hstar⟩ := hdom
      exact ap_core_bridge cite u hupos vstar hstar
  exact (exists_lonely_primPart_iff 14 v hv).mpr hu14

/-- The positive working LRC(14) statement from the primitive carrier inverse. -/
theorem LRC14_of_primitiveCarrierINV (cite : LRCUpTo13)
    (hinv : PrimitiveCarrierINV) :
    ∀ v : Fin 13 → ℤ, (∀ i, 0 < v i) → ∃ t : ℝ, Lonely 14 v t :=
  fun v hpos => lonely14_of_primitiveCarrierINV cite hinv v hpos

/-- The signed/nonzero global `LRC14Statement` follows after absolute-value
normalization. -/
theorem lrc14Statement_of_primitiveCarrierINV (cite : LRCUpTo13)
    (hinv : PrimitiveCarrierINV) : LRC14.LRC14Statement := by
  intro v hv
  let va : Fin 13 → ℤ := fun i => |v i|
  have hvapos : ∀ i, 0 < va i := fun i => abs_pos.mpr (hv i)
  obtain ⟨t, ht⟩ := lonely14_of_primitiveCarrierINV cite hinv va hvapos
  exact ⟨t, (LRC14.lonely_abs_iff 14 v t).mp ht⟩

end PrimitiveCarrierInverse
end LonelyRunner

#print axioms LonelyRunner.PrimitiveCarrierInverse.covering_iff_hasFourteenCarrier_of_no_lonely13
#print axioms LonelyRunner.PrimitiveCarrierInverse.primitiveCarrierINV_iff_primitiveINVcov
#print axioms LonelyRunner.PrimitiveCarrierInverse.primitiveCarrierINV_iff_trichotomy
#print axioms LonelyRunner.PrimitiveCarrierInverse.exists_lonely_primPart_iff
#print axioms LonelyRunner.PrimitiveCarrierInverse.covering_primPart_of_no_lonely14
#print axioms LonelyRunner.PrimitiveCarrierInverse.lonely14_of_primitiveCarrierINV
#print axioms LonelyRunner.PrimitiveCarrierInverse.LRC14_of_primitiveCarrierINV
#print axioms LonelyRunner.PrimitiveCarrierInverse.lrc14Statement_of_primitiveCarrierINV
