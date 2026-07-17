/-
  Exact 4-by-13 Zarankiewicz incidence bound used by the depth-five relation
  budget. Pair-codegree at most one forces at most nineteen incidences, so four
  supports of size at least five force a repeated runner pair. No \`sorry\`;
  no \`native_decide\`.
-/

import TournamentH7.LRCZarankiewiczGuardrail

namespace LonelyRunner
namespace LRCZarankiewiczFourThirteen

open Finset
open LRCZarankiewiczGuardrail

/-- Right degree in a `4 × 13` incidence graph. -/
def rightDegree (neighbor : Fin 4 → Finset (Fin 13)) (runner : Fin 13) : ℕ :=
  ((Finset.univ : Finset (Fin 4)).filter fun support => runner ∈ neighbor support).card

/-- Runners common to every left vertex in an unordered left pair. -/
def commonRunners (neighbor : Fin 4 → Finset (Fin 13))
    (pair : Finset (Fin 4)) : Finset (Fin 13) :=
  (Finset.univ : Finset (Fin 13)).filter fun runner =>
    ∀ support ∈ pair, runner ∈ neighbor support

theorem total_incidence_eq_sum_rightDegree (neighbor : Fin 4 → Finset (Fin 13)) :
    ∑ support, (neighbor support).card = ∑ runner, rightDegree neighbor runner := by
  simp only [rightDegree, Finset.card_filter]
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro support _
  simp

/-- Double-count length-two paths between left vertices. -/
theorem sum_choose_rightDegree_eq_sum_commonRunners
    (neighbor : Fin 4 → Finset (Fin 13)) :
    ∑ runner, (rightDegree neighbor runner).choose 2 =
      ∑ pair ∈ (Finset.univ : Finset (Fin 4)).powersetCard 2,
        (commonRunners neighbor pair).card := by
  have hpoint : ∀ runner : Fin 13,
      (rightDegree neighbor runner).choose 2 =
        (((Finset.univ : Finset (Fin 4)).powersetCard 2).filter fun pair =>
          ∀ support ∈ pair, runner ∈ neighbor support).card := by
    intro runner
    rw [rightDegree, ← Finset.card_powersetCard]
    congr 1
    ext pair
    simp only [Finset.mem_powersetCard, Finset.mem_filter]
    constructor
    · rintro ⟨hsub, hcard⟩
      exact ⟨⟨Finset.subset_univ pair, hcard⟩,
        fun support hs => (Finset.mem_filter.mp (hsub hs)).2⟩
    · rintro ⟨⟨_, hcard⟩, hall⟩
      exact ⟨fun support hs => Finset.mem_filter.mpr
        ⟨Finset.mem_univ support, hall support hs⟩, hcard⟩
  rw [Finset.sum_congr rfl fun runner _ => hpoint runner]
  simp only [Finset.card_filter, commonRunners]
  rw [Finset.sum_comm]

theorem rightDegree_le_four (neighbor : Fin 4 → Finset (Fin 13)) (runner : Fin 13) :
    rightDegree neighbor runner ≤ 4 := by
  unfold rightDegree
  exact le_trans (Finset.card_filter_le _ _) (by simp)

theorem degree_le_one_add_choose_two (degree : ℕ) (hdegree : degree ≤ 4) :
    degree ≤ 1 + degree.choose 2 := by
  interval_cases degree <;> decide

/-- Exact Zarankiewicz upper value `z(4,13;2,2) ≤ 19`. -/
theorem z_four_thirteen_le_nineteen (neighbor : Fin 4 → Finset (Fin 13))
    (hK22 : ∀ pair ∈ (Finset.univ : Finset (Fin 4)).powersetCard 2,
      (commonRunners neighbor pair).card ≤ 1) :
    ∑ support, (neighbor support).card ≤ 19 := by
  rw [total_incidence_eq_sum_rightDegree]
  have hcollision : ∑ runner, (rightDegree neighbor runner).choose 2 ≤ 6 := by
    rw [sum_choose_rightDegree_eq_sum_commonRunners]
    calc
      ∑ pair ∈ (Finset.univ : Finset (Fin 4)).powersetCard 2,
          (commonRunners neighbor pair).card ≤
        ∑ _pair ∈ (Finset.univ : Finset (Fin 4)).powersetCard 2, 1 := by
          exact Finset.sum_le_sum fun pair hpair => hK22 pair hpair
      _ = 6 := by decide
  calc
    ∑ runner, rightDegree neighbor runner ≤
        ∑ runner, (1 + (rightDegree neighbor runner).choose 2) := by
      exact Finset.sum_le_sum fun runner _ =>
        degree_le_one_add_choose_two _ (rightDegree_le_four neighbor runner)
    _ = 13 + ∑ runner, (rightDegree neighbor runner).choose 2 := by
      rw [Finset.sum_add_distrib]
      simp
    _ ≤ 13 + 6 := Nat.add_le_add_left hcollision 13
    _ = 19 := by decide

/-- Four supports of size at least five force a `K_{2,2}` collision. -/
theorem four_five_supports_force_shared_pair
    (neighbor : Fin 4 → Finset (Fin 13))
    (hsize : ∀ support, 5 ≤ (neighbor support).card) :
    ¬ (∀ pair ∈ (Finset.univ : Finset (Fin 4)).powersetCard 2,
      (commonRunners neighbor pair).card ≤ 1) := by
  intro hK22
  have hupper := z_four_thirteen_le_nineteen neighbor hK22
  have hlower : 20 ≤ ∑ support, (neighbor support).card := by
    calc
      20 = ∑ _support : Fin 4, 5 := by decide
      _ ≤ ∑ support, (neighbor support).card :=
        Finset.sum_le_sum fun support _ => hsize support
  omega

/-- Sharp support-five cap: a pair-unique family on thirteen runners contains
at most three supports of size at least five. -/
theorem card_le_three_of_five_le
    (supports : Finset (Finset (Fin 13)))
    (hunique : PairUnique supports)
    (hsize : ∀ support ∈ supports, 5 ≤ support.card) :
    supports.card ≤ 3 := by
  by_contra hcard
  have hfour : 4 ≤ supports.card := by omega
  obtain ⟨chosen, hchosen, hchosenCard⟩ :=
    Finset.exists_subset_card_eq hfour
  let enumerate : Fin 4 ≃ ↑chosen :=
    (finCongr hchosenCard.symm).trans chosen.equivFin.symm
  let neighbor : Fin 4 → Finset (Fin 13) := fun index => (enumerate index).1
  have hneighborMem (index : Fin 4) : neighbor index ∈ supports := by
    exact hchosen (enumerate index).2
  have hneighborSize (index : Fin 4) : 5 ≤ (neighbor index).card :=
    hsize (neighbor index) (hneighborMem index)
  have hneighborInj : Function.Injective neighbor := by
    intro first second heq
    apply enumerate.injective
    apply Subtype.ext
    exact heq
  have hK22 : ∀ pair ∈ (Finset.univ : Finset (Fin 4)).powersetCard 2,
      (commonRunners neighbor pair).card ≤ 1 := by
    intro pair hpair
    rw [Finset.card_le_one_iff]
    intro first second hfirst hsecond
    by_contra hne
    have hpairCard : pair.card = 2 :=
      (Finset.mem_powersetCard.mp hpair).2
    obtain ⟨left, right, hleftRight, rfl⟩ :=
      Finset.card_eq_two.mp hpairCard
    have hleft : first ∈ neighbor left ∧ first ∈ neighbor right := by
      simpa [commonRunners, hleftRight] using hfirst
    have hright : second ∈ neighbor left ∧ second ∈ neighbor right := by
      simpa [commonRunners, hleftRight] using hsecond
    have hsupportNe : neighbor left ≠ neighbor right :=
      hneighborInj.ne hleftRight
    have hdisjoint : Disjoint (supportPairs (neighbor left))
        (supportPairs (neighbor right)) :=
      hunique (hneighborMem left) (hneighborMem right) hsupportNe
    have hownedLeft : {first, second} ∈ supportPairs (neighbor left) := by
      rw [supportPairs, Finset.mem_powersetCard]
      refine ⟨?_, by simpa [hne]⟩
      intro runner hrunner
      simp only [Finset.mem_insert, Finset.mem_singleton] at hrunner
      rcases hrunner with rfl | rfl
      · exact hleft.1
      · exact hright.1
    have hownedRight : {first, second} ∈ supportPairs (neighbor right) := by
      rw [supportPairs, Finset.mem_powersetCard]
      refine ⟨?_, by simpa [hne]⟩
      intro runner hrunner
      simp only [Finset.mem_insert, Finset.mem_singleton] at hrunner
      rcases hrunner with rfl | rfl
      · exact hleft.2
      · exact hright.2
    exact (Finset.disjoint_left.mp hdisjoint hownedLeft) hownedRight
  exact four_five_supports_force_shared_pair neighbor hneighborSize hK22


/-! ## Axiom audit -/

#print axioms z_four_thirteen_le_nineteen
#print axioms four_five_supports_force_shared_pair
#print axioms card_le_three_of_five_le

end LRCZarankiewiczFourThirteen
end LonelyRunner
