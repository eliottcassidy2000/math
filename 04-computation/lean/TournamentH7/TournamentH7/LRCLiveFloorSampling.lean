/-
  TournamentH7.LRCLiveFloorSampling -- THM-984's live-floor composition

  `LRCGridSampling` proves the interface-free weight-one estimate for a
  separated family of open intervals.  This file supplies the LRC-specific
  composition:

    * a lonely rational point `p/q` lies in every inclusive residue band;
    * therefore every sampled point of an interval table inside `safePeriod`
      contributes to `liveCount`;
    * hence `q * mu0 - error <= liveCount`, and the explicit ruler
      `ceil ((error + 1) / mu0)` has positive live count;
    * positive live count feeds both the direct loneliness adapter and the
      capped-five census consumer.

  The one missing producer is stated honestly by the premises: the current
  measure files do not expose a theorem decomposing `safePeriod v` into a
  separated half-open/open interval table with a component bound.  Once that
  producer supplies at most `E` components and total length at least `mu0`,
  the theorems below give the headline `q * mu0 - E` bridge.

  Assumption/carrier audit: the useful vertices here are interval components
  and ruler points, not runners or arcs.  This quotient preserves safe-set
  membership and its count, but deliberately forgets endpoint owners and all
  pairwise/tournament data; Tournament Analysis therefore adds no clean
  observable to this purely one-dimensional counting step.

  Kernel-pure: no `sorry`, no `native_decide`.
-/
import Mathlib
import TournamentH7.LRCGridSampling
import TournamentH7.LRCLiveCountLonely
import TournamentH7.LRCDeepCertificate
import TournamentH7.LRCResidualMeasureFloor

namespace LonelyRunner
namespace LRC14Grand

open Finset

noncomputable section

/-- The exact converse needed by the sampler: a lonely rational point has all
residues in the inclusive safe band. -/
theorem inBand_of_lonely_rational (v : Fin 13 → ℤ) (q p : ℕ) (hq : 0 < q)
    (hlonely : Lonely 14 v ((p : ℝ) / (q : ℝ))) :
    ∀ i, LRC14Concrete.inBand v q p i := by
  intro i
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  have hqR : (0 : ℝ) < (q : ℝ) := by exact_mod_cast hq
  set x : ℤ := v i * (p : ℤ) with hx
  set r : ℤ := x % (q : ℤ) with hr
  set d : ℤ := x / (q : ℤ) with hd
  have hr0 : 0 ≤ r := by
    rw [hr]
    exact Int.emod_nonneg _ hqZ.ne'
  have hrq : r < (q : ℤ) := by
    rw [hr]
    exact Int.emod_lt_of_pos _ hqZ
  have hx0 : x - d * (q : ℤ) = r := by
    rw [hd, hr, Int.emod_def]
    ring
  have hx1 : x - (d + 1) * (q : ℤ) = r - (q : ℤ) := by
    rw [hd, hr, Int.emod_def]
    ring
  have hmargin : ∀ m : ℤ,
      ((q : ℝ) / 14) ≤ ((|x - m * (q : ℤ)| : ℤ) : ℝ) := by
    intro m
    have hm := hlonely i m
    have heq :
        (v i : ℝ) * ((p : ℝ) / (q : ℝ)) - m =
          ((x - m * (q : ℤ) : ℤ) : ℝ) / (q : ℝ) := by
      rw [hx]
      push_cast
      field_simp
    rw [heq, abs_div, abs_of_pos hqR, le_div_iff₀ hqR, ← Int.cast_abs] at hm
    norm_num at hm ⊢
    linarith
  have hloR : (q : ℝ) ≤ 14 * (r : ℝ) := by
    have h := hmargin d
    rw [hx0, abs_of_nonneg hr0] at h
    linarith
  have hhiR : 14 * (r : ℝ) ≤ 13 * (q : ℝ) := by
    have h := hmargin (d + 1)
    rw [hx1, abs_of_nonpos (by exact_mod_cast (sub_nonpos.mpr hrq.le)), neg_sub] at h
    push_cast at h
    linarith
  change (q : ℤ) ≤ 14 * ((v i * (p : ℤ)) % (q : ℤ)) ∧
    14 * ((v i * (p : ℤ)) % (q : ℤ)) ≤ 13 * (q : ℤ)
  rw [← hx, ← hr]
  constructor
  · exact_mod_cast hloR
  · exact_mod_cast hhiR

/-- A lonely rational point is counted as a live multiplier. -/
theorem bandCount_eq_zero_of_lonely_rational (v : Fin 13 → ℤ) (q p : ℕ) (hq : 0 < q)
    (hlonely : Lonely 14 v ((p : ℝ) / (q : ℝ))) :
    LRC14Concrete.bandCount v q p = 0 := by
  apply Finset.card_eq_zero.mpr
  rw [Finset.filter_eq_empty_iff]
  intro i _
  exact not_not.mpr (inBand_of_lonely_rational v q p hq hlonely i)

/-- **The LRC-specific weight-one sampler.**  A separated interval table
inside `safePeriod v` contributes its full grid sample to `liveCount`, losing
at most one point per interval. -/
theorem liveCount_ge_interval_table (v : Fin 13 → ℤ) (q : ℕ) (hq : 0 < q)
    {n : ℕ} (a b : Fin n → ℝ)
    (h0 : ∀ i, 0 ≤ a i) (h1 : ∀ i, b i ≤ 1)
    (hsep : ∀ i j, i ≠ j → b i ≤ a j ∨ b j ≤ a i)
    (hsub : ∀ i, Set.Ioo (a i) (b i) ⊆ safePeriod v) :
    (q : ℝ) * (∑ i, (b i - a i)) - n ≤ (LRC14Concrete.liveCount v q : ℝ) := by
  let sampled : Finset ℕ :=
    (Finset.Ioo 0 q).filter fun p : ℕ =>
      ∃ i, a i < (p : ℝ) / q ∧ (p : ℝ) / q < b i
  have hsample :
      (q : ℝ) * (∑ i, (b i - a i)) - n ≤ (sampled.card : ℝ) := by
    exact LRC14.Sampling.card_grid_family q hq n a b h0 h1 hsep
  have hsampledSub : sampled ⊆
      (Finset.Ioo 0 q).filter (fun p => LRC14Concrete.bandCount v q p = 0) := by
    intro p hp
    change p ∈ (Finset.Ioo 0 q).filter (fun p : ℕ =>
      ∃ i, a i < (p : ℝ) / q ∧ (p : ℝ) / q < b i) at hp
    rw [Finset.mem_filter] at hp
    obtain ⟨hpwindow, i, hpi⟩ := hp
    have hsafe := hsub i hpi
    exact Finset.mem_filter.mpr ⟨hpwindow,
      bandCount_eq_zero_of_lonely_rational v q p hq hsafe.2⟩
  calc
    (q : ℝ) * (∑ i, (b i - a i)) - n ≤ (sampled.card : ℝ) := hsample
    _ ≤ (LRC14Concrete.liveCount v q : ℝ) := by
      unfold LRC14Concrete.liveCount
      exact_mod_cast Finset.card_le_card hsampledSub

/-- A certified interval-length floor gives `q * mu0 - n <= liveCount`. -/
theorem liveCount_ge_of_interval_floor (v : Fin 13 → ℤ) (q : ℕ) (hq : 0 < q)
    {n : ℕ} (a b : Fin n → ℝ) (mu0 : ℝ)
    (h0 : ∀ i, 0 ≤ a i) (h1 : ∀ i, b i ≤ 1)
    (hsep : ∀ i j, i ≠ j → b i ≤ a j ∨ b j ≤ a i)
    (hsub : ∀ i, Set.Ioo (a i) (b i) ⊆ safePeriod v)
    (hfloor : mu0 ≤ ∑ i, (b i - a i)) :
    (q : ℝ) * mu0 - n ≤ (LRC14Concrete.liveCount v q : ℝ) := by
  calc
    (q : ℝ) * mu0 - n ≤ (q : ℝ) * (∑ i, (b i - a i)) - n := by
      gcongr
    _ ≤ (LRC14Concrete.liveCount v q : ℝ) :=
      liveCount_ge_interval_table v q hq a b h0 h1 hsep hsub

/-- The headline form with an externally supplied component/endpoint budget. -/
theorem liveCount_ge_of_endpoint_budget (v : Fin 13 → ℤ) (q : ℕ) (hq : 0 < q)
    {n : ℕ} (a b : Fin n → ℝ) (mu0 : ℝ) (error : ℕ)
    (h0 : ∀ i, 0 ≤ a i) (h1 : ∀ i, b i ≤ 1)
    (hsep : ∀ i j, i ≠ j → b i ≤ a j ∨ b j ≤ a i)
    (hsub : ∀ i, Set.Ioo (a i) (b i) ⊆ safePeriod v)
    (hfloor : mu0 ≤ ∑ i, (b i - a i)) (herror : n ≤ error) :
    (q : ℝ) * mu0 - error ≤ (LRC14Concrete.liveCount v q : ℝ) := by
  have hbase := liveCount_ge_of_interval_floor v q hq a b mu0
    h0 h1 hsep hsub hfloor
  have herrorR : (n : ℝ) ≤ (error : ℝ) := by exact_mod_cast herror
  linarith

/-- Positivity once the ruler beats the component error. -/
theorem liveCount_pos_of_interval_floor (v : Fin 13 → ℤ) (q : ℕ) (hq : 0 < q)
    {n : ℕ} (a b : Fin n → ℝ) (mu0 : ℝ)
    (h0 : ∀ i, 0 ≤ a i) (h1 : ∀ i, b i ≤ 1)
    (hsep : ∀ i j, i ≠ j → b i ≤ a j ∨ b j ≤ a i)
    (hsub : ∀ i, Set.Ioo (a i) (b i) ⊆ safePeriod v)
    (hfloor : mu0 ≤ ∑ i, (b i - a i))
    (hlarge : (n : ℝ) < (q : ℝ) * mu0) :
    0 < LRC14Concrete.liveCount v q := by
  have hlower := liveCount_ge_of_interval_floor v q hq a b mu0
    h0 h1 hsep hsub hfloor
  have hposR : (0 : ℝ) < (LRC14Concrete.liveCount v q : ℝ) := by linarith
  exact_mod_cast hposR

/-- A slightly sharper explicit modulus than `ceil (2E/mu0)`: adding one in
the numerator makes the endpoint-error inequality strict uniformly. -/
def liveFloorModulus (error : ℕ) (mu0 : ℝ) : ℕ :=
  ⌈((error : ℝ) + 1) / mu0⌉₊

theorem liveFloorModulus_pos (error : ℕ) {mu0 : ℝ} (hmu0 : 0 < mu0) :
    0 < liveFloorModulus error mu0 := by
  have h : 1 ≤ ⌈((error : ℝ) + 1) / mu0⌉₊ :=
    (Nat.one_le_ceil_iff).2 (by positivity)
  simpa [liveFloorModulus] using h

theorem error_lt_liveFloorModulus_mul (error : ℕ) {mu0 : ℝ} (hmu0 : 0 < mu0) :
    (error : ℝ) < (liveFloorModulus error mu0 : ℝ) * mu0 := by
  have hceil : ((error : ℝ) + 1) / mu0 ≤ (liveFloorModulus error mu0 : ℝ) :=
    Nat.le_ceil _
  have hmul := mul_le_mul_of_nonneg_right hceil hmu0.le
  have hcancel : (((error : ℝ) + 1) / mu0) * mu0 = (error : ℝ) + 1 := by
    field_simp
  rw [hcancel] at hmul
  linarith

/-- The explicit-`q0` positive-live-count form with an endpoint budget. -/
theorem liveCount_pos_at_explicit_modulus (v : Fin 13 → ℤ)
    {n : ℕ} (a b : Fin n → ℝ) (mu0 : ℝ) (error : ℕ)
    (hmu0 : 0 < mu0)
    (h0 : ∀ i, 0 ≤ a i) (h1 : ∀ i, b i ≤ 1)
    (hsep : ∀ i j, i ≠ j → b i ≤ a j ∨ b j ≤ a i)
    (hsub : ∀ i, Set.Ioo (a i) (b i) ⊆ safePeriod v)
    (hfloor : mu0 ≤ ∑ i, (b i - a i)) (herror : n ≤ error) :
    0 < LRC14Concrete.liveCount v (liveFloorModulus error mu0) := by
  have hbound := liveCount_ge_of_endpoint_budget v (liveFloorModulus error mu0)
    (liveFloorModulus_pos error hmu0) a b mu0 error
    h0 h1 hsep hsub hfloor herror
  have hgap := error_lt_liveFloorModulus_mul error hmu0
  have hposR :
      (0 : ℝ) < (LRC14Concrete.liveCount v (liveFloorModulus error mu0) : ℝ) := by
    linarith
  exact_mod_cast hposR

/-- End-to-end loneliness at the explicit modulus, using the direct adapter
from `LRCLiveCountLonely`. -/
theorem lonely_at_explicit_modulus_of_interval_floor (v : Fin 13 → ℤ)
    {n : ℕ} (a b : Fin n → ℝ) (mu0 : ℝ) (error : ℕ)
    (hmu0 : 0 < mu0)
    (h0 : ∀ i, 0 ≤ a i) (h1 : ∀ i, b i ≤ 1)
    (hsep : ∀ i j, i ≠ j → b i ≤ a j ∨ b j ≤ a i)
    (hsub : ∀ i, Set.Ioo (a i) (b i) ⊆ safePeriod v)
    (hfloor : mu0 ≤ ∑ i, (b i - a i)) (herror : n ≤ error) :
    ∃ t : ℝ, Lonely 14 v t := by
  exact LRC14Concrete.lonely_of_liveCount_pos
    (liveCount_pos_at_explicit_modulus v a b mu0 error hmu0
      h0 h1 hsep hsub hfloor herror)

/-- The same positive live count wired through the existing cap-five census
consumer, for callers already carrying that decidable gate. -/
theorem lonely_of_interval_floor_capped5 (v : Fin 13 → ℤ) (q : ℕ) (hq : 0 < q)
    (hv : ∀ i, v i ≠ 0) {n : ℕ} (a b : Fin n → ℝ) (mu0 : ℝ)
    (h0 : ∀ i, 0 ≤ a i) (h1 : ∀ i, b i ≤ 1)
    (hsep : ∀ i j, i ≠ j → b i ≤ a j ∨ b j ≤ a i)
    (hsub : ∀ i, Set.Ioo (a i) (b i) ⊆ safePeriod v)
    (hfloor : mu0 ≤ ∑ i, (b i - a i))
    (hlarge : (n : ℝ) < (q : ℝ) * mu0)
    (hcap : LRC14Concrete.CoverageCapped v q 5) :
    ∃ t : ℝ, Lonely 14 v t := by
  exact LRC14Concrete.lonely_of_census_capped5 v q hq hv hcap
    (liveCount_pos_of_interval_floor v q hq a b mu0
      h0 h1 hsep hsub hfloor hlarge)

end


/-! ## Axiom audit -/

#print axioms inBand_of_lonely_rational
#print axioms bandCount_eq_zero_of_lonely_rational
#print axioms liveCount_ge_interval_table
#print axioms liveCount_ge_of_interval_floor
#print axioms liveCount_ge_of_endpoint_budget
#print axioms liveCount_pos_of_interval_floor
#print axioms liveCount_pos_at_explicit_modulus
#print axioms lonely_at_explicit_modulus_of_interval_floor
#print axioms lonely_of_interval_floor_capped5

end LRC14Grand
end LonelyRunner
