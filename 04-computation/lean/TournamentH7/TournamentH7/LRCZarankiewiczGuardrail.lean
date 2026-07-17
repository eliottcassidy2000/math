import Mathlib

/-!
# A Zarankiewicz guardrail for small relation supports

The parallel-class/crossing quotient records which unordered pairs of runner
indices are jointly owned by a small relation support.  When distinct supports
never own the same pair, their pair sets are disjoint, so the total pair load is
at most `choose 13 2 = 78`.  This is the exact useful counting content of the
Zarankiewicz viewpoint for the depth-five relation strata.

This quotient deliberately does not encode relation coefficients, Fourier
weights, or signs.  Accordingly the theorem below is only a multiplicity
guardrail, not a positivity statement for `B5`.
-/

namespace LonelyRunner
namespace LRCZarankiewiczGuardrail

open scoped Classical

/-- The unordered index pairs owned by a support. -/
def supportPairs (support : Finset (Fin 13)) : Finset (Finset (Fin 13)) :=
  support.powersetCard 2

/-- Distinct supports are pair-unique when they never share an unordered pair. -/
def PairUnique (supports : Finset (Finset (Fin 13))) : Prop :=
  (supports : Set (Finset (Fin 13))).PairwiseDisjoint supportPairs

/-- The supports owning a fixed unordered index pair. -/
def pairOwners (supports : Finset (Finset (Fin 13))) (pair : Finset (Fin 13)) :
    Finset (Finset (Fin 13)) :=
  supports.filter fun support => pair ∈ supportPairs support

theorem supportPairs_subset_univPairs (support : Finset (Fin 13)) :
    supportPairs support ⊆ (Finset.univ : Finset (Fin 13)).powersetCard 2 := by
  intro pair hpair
  have hcard : pair.card = 2 := (Finset.mem_powersetCard.mp
    (show pair ∈ support.powersetCard 2 from hpair)).2
  exact Finset.mem_powersetCard.mpr ⟨Finset.subset_univ pair, hcard⟩

/-- Exact double counting: total pair load across supports equals total owner
multiplicity across the `78` unordered runner pairs. -/
theorem pair_load_eq_sum_owners (supports : Finset (Finset (Fin 13))) :
    ∑ support ∈ supports, (support.card.choose 2) =
      ∑ pair ∈ (Finset.univ : Finset (Fin 13)).powersetCard 2,
        (pairOwners supports pair).card := by
  calc
    ∑ support ∈ supports, (support.card.choose 2) =
        ∑ support ∈ supports,
          ∑ pair ∈ (Finset.univ : Finset (Fin 13)).powersetCard 2,
            if pair ∈ supportPairs support then 1 else 0 := by
      apply Finset.sum_congr rfl
      intro support _hsupport
      calc
        support.card.choose 2 = (supportPairs support).card := by
          simp [supportPairs]
        _ = ∑ pair ∈ (Finset.univ : Finset (Fin 13)).powersetCard 2,
              if pair ∈ supportPairs support then 1 else 0 :=
          Finset.card_eq_sum_ite (supportPairs_subset_univPairs support)
    _ = ∑ pair ∈ (Finset.univ : Finset (Fin 13)).powersetCard 2,
          ∑ support ∈ supports, if pair ∈ supportPairs support then 1 else 0 := by
      rw [Finset.sum_comm]
    _ = ∑ pair ∈ (Finset.univ : Finset (Fin 13)).powersetCard 2,
          (pairOwners supports pair).card := by
      apply Finset.sum_congr rfl
      intro pair _hpair
      simp [pairOwners]

/-- If every unordered pair has at most `multiplicity` owners, total relation
pair load is at most `78 * multiplicity`. -/
theorem pair_load_le_78_mul (supports : Finset (Finset (Fin 13)))
    (multiplicity : ℕ)
    (howners : ∀ pair ∈ (Finset.univ : Finset (Fin 13)).powersetCard 2,
      (pairOwners supports pair).card ≤ multiplicity) :
    ∑ support ∈ supports, (support.card.choose 2) ≤ 78 * multiplicity := by
  rw [pair_load_eq_sum_owners]
  calc
    ∑ pair ∈ (Finset.univ : Finset (Fin 13)).powersetCard 2,
        (pairOwners supports pair).card ≤
        ∑ _pair ∈ (Finset.univ : Finset (Fin 13)).powersetCard 2,
          multiplicity := by
      exact Finset.sum_le_sum fun pair hpair => howners pair hpair
    _ = 78 * multiplicity := by
      norm_num [Finset.card_powersetCard, Nat.choose]

theorem pair_load_le_78 (supports : Finset (Finset (Fin 13)))
    (hunique : PairUnique supports) :
    ∑ support ∈ supports, (support.card.choose 2) ≤ 78 := by
  have hunion : supports.biUnion supportPairs ⊆
      (Finset.univ : Finset (Fin 13)).powersetCard 2 := by
    intro pair hpair
    rw [Finset.mem_biUnion] at hpair
    obtain ⟨support, _hsupport, hpairSupport⟩ := hpair
    exact supportPairs_subset_univPairs support hpairSupport
  calc
    ∑ support ∈ supports, (support.card.choose 2)
        = ∑ support ∈ supports, (supportPairs support).card := by
            apply Finset.sum_congr rfl
            intro support _hsupport
            simp [supportPairs]
    _ = (supports.biUnion supportPairs).card := by
          rw [Finset.card_biUnion hunique]
    _ ≤ ((Finset.univ : Finset (Fin 13)).powersetCard 2).card :=
          Finset.card_le_card hunion
    _ = 78 := by decide

/-- A uniform lower bound on support size converts directly to a lower bound
on pair load. -/
theorem choose_mul_card_le_pair_load (supports : Finset (Finset (Fin 13))) (size : ℕ)
    (hsize : ∀ support ∈ supports, size ≤ support.card) :
    (size.choose 2) * supports.card ≤
      ∑ support ∈ supports, (support.card.choose 2) := by
  calc
    (size.choose 2) * supports.card =
        ∑ _support ∈ supports, (size.choose 2) := by simp [Nat.mul_comm]
    _ ≤ ∑ support ∈ supports, (support.card.choose 2) := by
      exact Finset.sum_le_sum fun support hsupport =>
        Nat.choose_le_choose 2 (hsize support hsupport)

/-- Pair-unique supports of size at least three number at most `26`. -/
theorem card_le_26_of_three_le (supports : Finset (Finset (Fin 13)))
    (hunique : PairUnique supports) (hsize : ∀ support ∈ supports, 3 ≤ support.card) :
    supports.card ≤ 26 := by
  have hthree := choose_mul_card_le_pair_load supports 3 hsize
  have hload := pair_load_le_78 supports hunique
  norm_num at hthree
  omega

/-- Pair-unique supports of size at least four number at most `13`. -/
theorem card_le_13_of_four_le (supports : Finset (Finset (Fin 13)))
    (hunique : PairUnique supports) (hsize : ∀ support ∈ supports, 4 ≤ support.card) :
    supports.card ≤ 13 := by
  have hfour := choose_mul_card_le_pair_load supports 4 hsize
  have hload := pair_load_le_78 supports hunique
  norm_num [Nat.choose] at hfour
  omega

/-- Pair-unique supports of size at least five number at most `7`. -/
theorem card_le_7_of_five_le (supports : Finset (Finset (Fin 13)))
    (hunique : PairUnique supports) (hsize : ∀ support ∈ supports, 5 ≤ support.card) :
    supports.card ≤ 7 := by
  have hfive := choose_mul_card_le_pair_load supports 5 hsize
  have hload := pair_load_le_78 supports hunique
  norm_num [Nat.choose] at hfive
  omega

/-- More than `26` support-`≥3` relations force a repeated owner pair, hence
a collision of parallel classes in the incidence quotient. -/
theorem not_pairUnique_of_26_lt_card (supports : Finset (Finset (Fin 13)))
    (hsize : ∀ support ∈ supports, 3 ≤ support.card) (hmany : 26 < supports.card) :
    ¬ PairUnique supports := by
  intro hunique
  exact (Nat.not_lt_of_ge (card_le_26_of_three_le supports hunique hsize)) hmany

/-- Failure of pair uniqueness is an explicit parallel-class collision: two
distinct supports own the same unordered index pair. -/
theorem exists_shared_pair_of_not_pairUnique (supports : Finset (Finset (Fin 13)))
    (hcollision : ¬ PairUnique supports) :
    ∃ first ∈ supports, ∃ second ∈ supports, first ≠ second ∧
      ∃ pair, pair ∈ supportPairs first ∧ pair ∈ supportPairs second := by
  simp only [PairUnique, Set.PairwiseDisjoint, Set.Pairwise, Function.onFun,
    not_forall] at hcollision
  obtain ⟨first, hfirst, second, hsecond, hne, hnotDisjoint⟩ := hcollision
  obtain ⟨pair, hpairFirst, hpairSecond⟩ :=
    Finset.not_disjoint_iff.mp hnotDisjoint
  exact ⟨first, hfirst, second, hsecond, hne, pair, hpairFirst, hpairSecond⟩

/-- The support-`≥3` collision alternative in witness form. -/
theorem exists_shared_pair_of_26_lt_card (supports : Finset (Finset (Fin 13)))
    (hsize : ∀ support ∈ supports, 3 ≤ support.card) (hmany : 26 < supports.card) :
    ∃ first ∈ supports, ∃ second ∈ supports, first ≠ second ∧
      ∃ pair, pair ∈ supportPairs first ∧ pair ∈ supportPairs second :=
  exists_shared_pair_of_not_pairUnique supports
    (not_pairUnique_of_26_lt_card supports hsize hmany)

/-! ## The tiny-scale relation floor -/

/-- Thirteen distinct positive speeds below `40` cannot be Sidon: among the
`78` unordered index pairs, two distinct pairs have the same sum.  The two
pairs may share one endpoint (a support-three relation) or be disjoint (a
support-four relation), exactly matching THM-935's small-scale floor. -/
theorem exists_small_support_relation_of_lt_40 (v : Fin 13 → ℕ)
    (hpositive : ∀ i, 0 < v i) (hsmall : ∀ i, v i < 40)
    (hinjective : Function.Injective v) :
    ∃ first second : Finset (Fin 13),
      first ∈ (Finset.univ : Finset (Fin 13)).powersetCard 2 ∧
      second ∈ (Finset.univ : Finset (Fin 13)).powersetCard 2 ∧
      first ≠ second ∧
      (∑ i ∈ first, v i) = ∑ i ∈ second, v i := by
  let pairs : Finset (Finset (Fin 13)) :=
    (Finset.univ : Finset (Fin 13)).powersetCard 2
  let possibleSums : Finset ℕ := Finset.Icc 2 77
  have hcard : possibleSums.card < pairs.card := by
    dsimp [possibleSums, pairs]
    norm_num [Finset.card_powersetCard, Nat.choose]
  have hmaps : Set.MapsTo (fun pair : Finset (Fin 13) => ∑ i ∈ pair, v i)
      (pairs : Set (Finset (Fin 13))) (possibleSums : Set ℕ) := by
    intro pair hpair
    have hpairCard : pair.card = 2 :=
      (Finset.mem_powersetCard.mp (show pair ∈ pairs from hpair)).2
    obtain ⟨first, second, hne, rfl⟩ := Finset.card_eq_two.mp hpairCard
    have hvaluesNe : v first ≠ v second := hinjective.ne hne
    simp only [Finset.sum_insert, Finset.sum_singleton, Finset.mem_Icc]
    constructor <;> omega
  obtain ⟨first, hfirst, second, hsecond, hne, heq⟩ :=
    Finset.exists_ne_map_eq_of_card_lt_of_maps_to hcard hmaps
  exact ⟨first, hfirst, second, hsecond, hne, heq⟩

/-! ## Axiom audit -/

#print axioms pair_load_le_78
#print axioms pair_load_eq_sum_owners
#print axioms pair_load_le_78_mul
#print axioms card_le_26_of_three_le
#print axioms card_le_13_of_four_le
#print axioms card_le_7_of_five_le
#print axioms not_pairUnique_of_26_lt_card
#print axioms exists_shared_pair_of_26_lt_card
#print axioms exists_small_support_relation_of_lt_40

end LRCZarankiewiczGuardrail
end LonelyRunner
