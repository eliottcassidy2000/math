/-
  Exact K_{2,2}-free incidence values for six through nine left vertices
  and thirteen runners. The degree-collision upper bound is paired with
  explicit sharp witnesses. No `sorry`; no `native_decide`.
-/

import TournamentH7.LRCZarankiewiczGuardrail

namespace LonelyRunner
namespace LRCZarankiewiczSixNineThirteen

open Finset

def rightDegreeM {m : ℕ} (neighbor : Fin m → Finset (Fin 13))
    (runner : Fin 13) : ℕ :=
  ((Finset.univ : Finset (Fin m)).filter fun support =>
    runner ∈ neighbor support).card

def commonRunnersM {m : ℕ} (neighbor : Fin m → Finset (Fin 13))
    (pair : Finset (Fin m)) : Finset (Fin 13) :=
  (Finset.univ : Finset (Fin 13)).filter fun runner =>
    ∀ support ∈ pair, runner ∈ neighbor support

theorem total_incidence_eq_sum_rightDegreeM {m : ℕ}
    (neighbor : Fin m → Finset (Fin 13)) :
    ∑ support, (neighbor support).card =
      ∑ runner, rightDegreeM neighbor runner := by
  simp only [rightDegreeM, Finset.card_filter]
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro support _
  simp

theorem sum_choose_rightDegreeM_eq_sum_commonRunnersM {m : ℕ}
    (neighbor : Fin m → Finset (Fin 13)) :
    ∑ runner, (rightDegreeM neighbor runner).choose 2 =
      ∑ pair ∈ (Finset.univ : Finset (Fin m)).powersetCard 2,
        (commonRunnersM neighbor pair).card := by
  have hpoint : ∀ runner : Fin 13,
      (rightDegreeM neighbor runner).choose 2 =
        (((Finset.univ : Finset (Fin m)).powersetCard 2).filter fun pair =>
          ∀ support ∈ pair, runner ∈ neighbor support).card := by
    intro runner
    rw [rightDegreeM, ← Finset.card_powersetCard]
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
  simp only [Finset.card_filter, commonRunnersM]
  rw [Finset.sum_comm]

theorem rightDegreeM_le {m : ℕ} (neighbor : Fin m → Finset (Fin 13))
    (runner : Fin 13) : rightDegreeM neighbor runner ≤ m := by
  unfold rightDegreeM
  calc
    ((Finset.univ : Finset (Fin m)).filter fun support =>
        runner ∈ neighbor support).card ≤
        (Finset.univ : Finset (Fin m)).card := Finset.card_filter_le _ _
    _ = m := by simp

theorem two_mul_degree_le_three_add_choose_two (degree : ℕ)
    (hdegree : degree ≤ 9) :
    2 * degree ≤ 3 + degree.choose 2 := by
  interval_cases degree <;> decide

theorem two_mul_incidence_le_39_add_choose {m : ℕ} (hm : m ≤ 9)
    (neighbor : Fin m → Finset (Fin 13))
    (hK22 : ∀ pair ∈ (Finset.univ : Finset (Fin m)).powersetCard 2,
      (commonRunnersM neighbor pair).card ≤ 1) :
    2 * (∑ support, (neighbor support).card) ≤ 39 + m.choose 2 := by
  have hcollision :
      ∑ runner, (rightDegreeM neighbor runner).choose 2 ≤ m.choose 2 := by
    rw [sum_choose_rightDegreeM_eq_sum_commonRunnersM]
    calc
      ∑ pair ∈ (Finset.univ : Finset (Fin m)).powersetCard 2,
          (commonRunnersM neighbor pair).card ≤
        ∑ _pair ∈ (Finset.univ : Finset (Fin m)).powersetCard 2, 1 := by
          exact Finset.sum_le_sum fun pair hpair => hK22 pair hpair
      _ = m.choose 2 := by simp [Finset.card_powersetCard]
  rw [total_incidence_eq_sum_rightDegreeM]
  calc
    2 * ∑ runner, rightDegreeM neighbor runner =
        ∑ runner, 2 * rightDegreeM neighbor runner := by
          rw [Finset.mul_sum]
    _ ≤ ∑ runner, (3 + (rightDegreeM neighbor runner).choose 2) := by
      exact Finset.sum_le_sum fun runner _ =>
        two_mul_degree_le_three_add_choose_two _
          (le_trans (rightDegreeM_le neighbor runner) hm)
    _ = 39 + ∑ runner, (rightDegreeM neighbor runner).choose 2 := by
      rw [Finset.sum_add_distrib]
      simp
    _ ≤ 39 + m.choose 2 := Nat.add_le_add_left hcollision 39

theorem z_six_thirteen_le_27 (neighbor : Fin 6 → Finset (Fin 13))
    (hK22 : ∀ pair ∈ (Finset.univ : Finset (Fin 6)).powersetCard 2,
      (commonRunnersM neighbor pair).card ≤ 1) :
    ∑ support, (neighbor support).card ≤ 27 := by
  have h := two_mul_incidence_le_39_add_choose (by decide : 6 ≤ 9) neighbor hK22
  norm_num [Nat.choose] at h ⊢
  omega

theorem z_seven_thirteen_le_30 (neighbor : Fin 7 → Finset (Fin 13))
    (hK22 : ∀ pair ∈ (Finset.univ : Finset (Fin 7)).powersetCard 2,
      (commonRunnersM neighbor pair).card ≤ 1) :
    ∑ support, (neighbor support).card ≤ 30 := by
  have h := two_mul_incidence_le_39_add_choose (by decide : 7 ≤ 9) neighbor hK22
  norm_num [Nat.choose] at h ⊢
  omega

theorem z_eight_thirteen_le_33 (neighbor : Fin 8 → Finset (Fin 13))
    (hK22 : ∀ pair ∈ (Finset.univ : Finset (Fin 8)).powersetCard 2,
      (commonRunnersM neighbor pair).card ≤ 1) :
    ∑ support, (neighbor support).card ≤ 33 := by
  have h := two_mul_incidence_le_39_add_choose (by decide : 8 ≤ 9) neighbor hK22
  norm_num [Nat.choose] at h ⊢
  omega

theorem z_nine_thirteen_le_37 (neighbor : Fin 9 → Finset (Fin 13))
    (hK22 : ∀ pair ∈ (Finset.univ : Finset (Fin 9)).powersetCard 2,
      (commonRunnersM neighbor pair).card ≤ 1) :
    ∑ support, (neighbor support).card ≤ 37 := by
  have h := two_mul_incidence_le_39_add_choose (by decide : 9 ≤ 9) neighbor hK22
  norm_num [Nat.choose] at h ⊢
  omega

def witness6 : Fin 6 → Finset (Fin 13) := ![
  {0, 1, 2, 3},
  {0, 4, 5, 6},
  {0, 7, 8, 9},
  {1, 4, 7, 10, 11},
  {2, 5, 8, 10, 12},
  {3, 6, 9, 11, 12}
]

def witness7 : Fin 7 → Finset (Fin 13) := ![
  {0, 1, 2},
  {0, 3, 4, 5},
  {0, 6, 7, 8, 9},
  {1, 3, 6, 10},
  {1, 4, 7, 11, 12},
  {2, 3, 8, 11},
  {2, 5, 9, 10, 12}
]

def witness8 : Fin 8 → Finset (Fin 13) := ![
  {0, 1, 2, 7},
  {0, 3, 4, 8},
  {0, 5, 6, 9},
  {1, 3, 5, 10},
  {1, 4, 6, 11},
  {2, 3, 6, 12},
  {2, 4, 5},
  {7, 8, 9, 10, 11, 12}
]

def witness9 : Fin 9 → Finset (Fin 13) := ![
  {0, 3, 6, 9},
  {0, 4, 7, 10},
  {0, 5, 8, 11, 12},
  {1, 3, 8, 10},
  {1, 4, 6, 11},
  {1, 5, 7, 9},
  {2, 3, 7, 12},
  {2, 4, 8, 9},
  {2, 5, 6, 10}
]

theorem witness6_K22 :
    ∀ pair ∈ (Finset.univ : Finset (Fin 6)).powersetCard 2,
      (commonRunnersM witness6 pair).card ≤ 1 := by decide

theorem witness7_K22 :
    ∀ pair ∈ (Finset.univ : Finset (Fin 7)).powersetCard 2,
      (commonRunnersM witness7 pair).card ≤ 1 := by
  set_option maxRecDepth 10000 in decide

theorem witness8_K22 :
    ∀ pair ∈ (Finset.univ : Finset (Fin 8)).powersetCard 2,
      (commonRunnersM witness8 pair).card ≤ 1 := by
  set_option maxRecDepth 10000 in decide

theorem witness9_K22 :
    ∀ pair ∈ (Finset.univ : Finset (Fin 9)).powersetCard 2,
      (commonRunnersM witness9 pair).card ≤ 1 := by
  set_option maxRecDepth 10000 in decide

theorem witness6_edges : ∑ support, (witness6 support).card = 27 := by decide
theorem witness7_edges : ∑ support, (witness7 support).card = 30 := by decide
theorem witness8_edges : ∑ support, (witness8 support).card = 33 := by decide
theorem witness9_edges : ∑ support, (witness9 support).card = 37 := by decide

#print axioms z_six_thirteen_le_27
#print axioms z_seven_thirteen_le_30
#print axioms z_eight_thirteen_le_33
#print axioms z_nine_thirteen_le_37
#print axioms witness6_K22
#print axioms witness9_edges

end LRCZarankiewiczSixNineThirteen
end LonelyRunner
