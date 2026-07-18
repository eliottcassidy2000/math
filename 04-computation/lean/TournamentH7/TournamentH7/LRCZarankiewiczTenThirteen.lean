/-
  Sharp K_{2,2}-free incidence bounds with thirteen right vertices and
  ten through thirteen left vertices.  The ten-by-thirteen upper bound has
  one non-convex step: a hypothetical 41-edge family would have thirteen
  blocks of size three or four, exactly two of size four, and would cover
  every left pair exactly once.  Counting companions at each left vertex
  forces that vertex into exactly one of the two four-blocks, which is
  impossible because those blocks have only eight incidences.

  This is an incidence guardrail only.  It does not by itself exclude any
  lonely-runner threshold graph whose edge count lies below the bound.

  Tournament-analysis audit: use left supports, rather than runners, as the
  vertices and pairwise intersection cardinality as the observable.  The
  collision switch is `2 <= |A ∩ B|`; index order supplies the tie Hamiltonian
  path.  On a K22-free family the switch is identically off, so the resulting
  tie tournament is transitive: its scores are `0, ..., m - 1`, it has no
  directed cycle, every SCC is a singleton, and its Hamiltonian path is unique.
  This quotient preserves exactly the repeated-right-pair obstruction and
  destroys coefficients, signs, circle phase, and divisibility data.  The
  challenged assumption is that taking more supports merely strengthens a
  density argument: the Singer family instead saturates a parallel-class
  geometry, while the ten-support 41-edge equality case fails locally.
-/

import TournamentH7.LRCZarankiewiczSixNineThirteen

namespace LonelyRunner
namespace LRCZarankiewiczTenThirteen

open Finset
open LRCZarankiewiczSixNineThirteen
open LRCZarankiewiczGuardrail

theorem three_mul_degree_le_six_add_choose_two (degree : ℕ)
    (hdegree : degree ≤ 13) :
    3 * degree ≤ 6 + degree.choose 2 := by
  interval_cases degree <;> decide

theorem collision_sum_le_choose {m : ℕ}
    (neighbor : Fin m → Finset (Fin 13))
    (hK22 : ∀ pair ∈ (Finset.univ : Finset (Fin m)).powersetCard 2,
      (commonRunnersM neighbor pair).card ≤ 1) :
    ∑ runner, (rightDegreeM neighbor runner).choose 2 ≤ m.choose 2 := by
  rw [sum_choose_rightDegreeM_eq_sum_commonRunnersM]
  calc
    ∑ pair ∈ (Finset.univ : Finset (Fin m)).powersetCard 2,
        (commonRunnersM neighbor pair).card ≤
      ∑ _pair ∈ (Finset.univ : Finset (Fin m)).powersetCard 2, 1 := by
        exact Finset.sum_le_sum fun pair hpair => hK22 pair hpair
    _ = m.choose 2 := by simp [Finset.card_powersetCard]

theorem three_mul_incidence_le_78_add_collision {m : ℕ} (hm : m ≤ 13)
    (neighbor : Fin m → Finset (Fin 13)) :
    3 * (∑ support, (neighbor support).card) ≤
      78 + ∑ runner, (rightDegreeM neighbor runner).choose 2 := by
  rw [total_incidence_eq_sum_rightDegreeM]
  calc
    3 * ∑ runner, rightDegreeM neighbor runner =
        ∑ runner, 3 * rightDegreeM neighbor runner := by
          rw [Finset.mul_sum]
    _ ≤ ∑ runner, (6 + (rightDegreeM neighbor runner).choose 2) := by
      exact Finset.sum_le_sum fun runner _ =>
        three_mul_degree_le_six_add_choose_two _
          (le_trans (rightDegreeM_le neighbor runner) hm)
    _ = 78 + ∑ runner, (rightDegreeM neighbor runner).choose 2 := by
      rw [Finset.sum_add_distrib]
      simp

theorem three_mul_incidence_le_78_add_choose {m : ℕ} (hm : m ≤ 13)
    (neighbor : Fin m → Finset (Fin 13))
    (hK22 : ∀ pair ∈ (Finset.univ : Finset (Fin m)).powersetCard 2,
      (commonRunnersM neighbor pair).card ≤ 1) :
    3 * (∑ support, (neighbor support).card) ≤ 78 + m.choose 2 := by
  exact (three_mul_incidence_le_78_add_collision hm neighbor).trans
    (Nat.add_le_add_left (collision_sum_le_choose neighbor hK22) 78)

theorem z_ten_thirteen_le_41 (neighbor : Fin 10 → Finset (Fin 13))
    (hK22 : ∀ pair ∈ (Finset.univ : Finset (Fin 10)).powersetCard 2,
      (commonRunnersM neighbor pair).card ≤ 1) :
    ∑ support, (neighbor support).card ≤ 41 := by
  have h := three_mul_incidence_le_78_add_choose
    (by decide : 10 ≤ 13) neighbor hK22
  norm_num [Nat.choose] at h ⊢
  omega

theorem degree_eq_three_or_four_of_tight (degree : ℕ)
    (hdegree : degree ≤ 10)
    (htight : 3 * degree = 6 + degree.choose 2) :
    degree = 3 ∨ degree = 4 := by
  interval_cases degree <;> norm_num [Nat.choose] at htight
  all_goals simp

theorem eq_one_of_nine_eq_two_mul_add (base remainder : ℕ)
    (hsum : 9 = 2 * base + remainder) (hremainder : remainder ≤ 2) :
    remainder = 1 := by
  omega

theorem rightDegree_sub_one_eq_erase_card
    (neighbor : Fin 10 → Finset (Fin 13))
    (support : Fin 10) (runner : Fin 13)
    (hrunner : runner ∈ neighbor support) :
    rightDegreeM neighbor runner - 1 =
      ((Finset.univ.erase support).filter fun other =>
        runner ∈ neighbor other).card := by
  unfold rightDegreeM
  have hsupport : support ∈
      ((Finset.univ : Finset (Fin 10)).filter fun other =>
        runner ∈ neighbor other) := by
    simp [hrunner]
  rw [← Finset.card_erase_of_mem hsupport]
  congr 1
  ext other
  simp

theorem incident_companion_sum_eq_pair_intersections
    (neighbor : Fin 10 → Finset (Fin 13)) (support : Fin 10) :
    ∑ runner ∈ neighbor support, (rightDegreeM neighbor runner - 1) =
      ∑ other ∈ (Finset.univ : Finset (Fin 10)).erase support,
        (neighbor support ∩ neighbor other).card := by
  calc
    ∑ runner ∈ neighbor support, (rightDegreeM neighbor runner - 1) =
        ∑ runner ∈ neighbor support,
          (((Finset.univ : Finset (Fin 10)).erase support).filter fun other =>
            runner ∈ neighbor other).card := by
      apply Finset.sum_congr rfl
      intro runner hrunner
      exact rightDegree_sub_one_eq_erase_card neighbor support runner hrunner
    _ = ∑ runner ∈ neighbor support,
          ∑ other ∈ (Finset.univ : Finset (Fin 10)).erase support,
            if runner ∈ neighbor other then 1 else 0 := by
      apply Finset.sum_congr rfl
      intro runner _
      simp
    _ = ∑ other ∈ (Finset.univ : Finset (Fin 10)).erase support,
          ∑ runner ∈ neighbor support,
            if runner ∈ neighbor other then 1 else 0 := by
      rw [Finset.sum_comm]
    _ = ∑ other ∈ (Finset.univ : Finset (Fin 10)).erase support,
        (neighbor support ∩ neighbor other).card := by
      apply Finset.sum_congr rfl
      intro other _
      simp

theorem commonRunnersM_pair_eq_inter {m : ℕ}
    (neighbor : Fin m → Finset (Fin 13))
    (first second : Fin m) :
    commonRunnersM neighbor {first, second} =
      neighbor first ∩ neighbor second := by
  ext runner
  simp [commonRunnersM]

theorem k22_of_pairwise_inter_card_le_one {m : ℕ}
    (neighbor : Fin m → Finset (Fin 13))
    (hpairwise : ∀ first second, first ≠ second →
      (neighbor first ∩ neighbor second).card ≤ 1) :
    ∀ pair ∈ (Finset.univ : Finset (Fin m)).powersetCard 2,
      (commonRunnersM neighbor pair).card ≤ 1 := by
  intro pair hpair
  have hpairCard : pair.card = 2 :=
    (Finset.mem_powersetCard.mp hpair).2
  obtain ⟨first, second, hne, rfl⟩ := Finset.card_eq_two.mp hpairCard
  rw [commonRunnersM_pair_eq_inter]
  exact hpairwise first second hne

theorem z_ten_thirteen_ne_41 (neighbor : Fin 10 → Finset (Fin 13))
    (hK22 : ∀ pair ∈ (Finset.univ : Finset (Fin 10)).powersetCard 2,
      (commonRunnersM neighbor pair).card ≤ 1) :
    ∑ support, (neighbor support).card ≠ 41 := by
  intro hedges
  have hsumDegree : ∑ runner, rightDegreeM neighbor runner = 41 := by
    rw [← total_incidence_eq_sum_rightDegreeM]
    exact hedges
  have hcollisionLe :
      ∑ runner, (rightDegreeM neighbor runner).choose 2 ≤ 45 := by
    have h := collision_sum_le_choose neighbor hK22
    norm_num [Nat.choose] at h ⊢
    exact h
  have hincidence := three_mul_incidence_le_78_add_collision
    (by decide : 10 ≤ 13) neighbor
  rw [hedges] at hincidence
  have hcollisionEq :
      ∑ runner, (rightDegreeM neighbor runner).choose 2 = 45 := by
    omega
  have hsumTight :
      ∑ runner, 3 * rightDegreeM neighbor runner =
        ∑ runner, (6 + (rightDegreeM neighbor runner).choose 2) := by
    rw [← Finset.mul_sum, hsumDegree, Finset.sum_add_distrib]
    simp [hcollisionEq]
  have hdegreeTight (runner : Fin 13) :
      3 * rightDegreeM neighbor runner =
        6 + (rightDegreeM neighbor runner).choose 2 := by
    exact (Finset.sum_eq_sum_iff_of_le fun runner _ =>
      three_mul_degree_le_six_add_choose_two _
        (le_trans (rightDegreeM_le neighbor runner) (by decide : 10 ≤ 13))).mp
          hsumTight runner (Finset.mem_univ runner)
  have hdegreeCases (runner : Fin 13) :
      rightDegreeM neighbor runner = 3 ∨
        rightDegreeM neighbor runner = 4 :=
    degree_eq_three_or_four_of_tight _
      (rightDegreeM_le neighbor runner) (hdegreeTight runner)
  have hcommonSum :
      ∑ pair ∈ (Finset.univ : Finset (Fin 10)).powersetCard 2,
          (commonRunnersM neighbor pair).card = 45 := by
    rw [← sum_choose_rightDegreeM_eq_sum_commonRunnersM]
    exact hcollisionEq
  have hcommonSumOnes :
      ∑ pair ∈ (Finset.univ : Finset (Fin 10)).powersetCard 2,
          (commonRunnersM neighbor pair).card =
        ∑ _pair ∈ (Finset.univ : Finset (Fin 10)).powersetCard 2, 1 := by
    rw [hcommonSum]
    norm_num [Finset.card_powersetCard, Nat.choose]
  have hcommonOne :
      ∀ pair ∈ (Finset.univ : Finset (Fin 10)).powersetCard 2,
        (commonRunnersM neighbor pair).card = 1 := by
    exact (Finset.sum_eq_sum_iff_of_le fun pair hpair =>
      hK22 pair hpair).mp hcommonSumOnes
  let high : Finset (Fin 13) :=
    (Finset.univ : Finset (Fin 13)).filter fun runner =>
      rightDegreeM neighbor runner = 4
  have hdegreeExpand (runner : Fin 13) :
      rightDegreeM neighbor runner =
        3 + if rightDegreeM neighbor runner = 4 then 1 else 0 := by
    rcases hdegreeCases runner with hthree | hfour
    · simp [hthree]
    · simp [hfour]
  have hhighCard : high.card = 2 := by
    have hdecompose :
        ∑ runner, rightDegreeM neighbor runner = 39 + high.card := by
      calc
        ∑ runner, rightDegreeM neighbor runner =
            ∑ runner,
              (3 + if rightDegreeM neighbor runner = 4 then 1 else 0) := by
          exact Finset.sum_congr rfl fun runner _ => hdegreeExpand runner
        _ = 39 + high.card := by
          simp [high, Finset.sum_add_distrib]
    omega
  have hcompanion (support : Fin 10) :
      ∑ runner ∈ neighbor support,
        (rightDegreeM neighbor runner - 1) = 9 := by
    rw [incident_companion_sum_eq_pair_intersections]
    calc
      ∑ other ∈ (Finset.univ : Finset (Fin 10)).erase support,
          (neighbor support ∩ neighbor other).card =
          ∑ _other ∈ (Finset.univ : Finset (Fin 10)).erase support, 1 := by
        apply Finset.sum_congr rfl
        intro other hother
        have hne : support ≠ other := by
          exact fun heq => (Finset.mem_erase.mp hother).1 heq.symm
        have hpair : {support, other} ∈
            (Finset.univ : Finset (Fin 10)).powersetCard 2 := by
          rw [Finset.mem_powersetCard]
          exact ⟨Finset.subset_univ _, by simp [hne]⟩
        rw [← commonRunnersM_pair_eq_inter neighbor support other]
        exact hcommonOne _ hpair
      _ = 9 := by simp
  have hsubExpand (runner : Fin 13) :
      rightDegreeM neighbor runner - 1 =
        2 + if rightDegreeM neighbor runner = 4 then 1 else 0 := by
    rcases hdegreeCases runner with hthree | hfour
    · simp [hthree]
    · simp [hfour]
  have hhighIncident (support : Fin 10) :
      (neighbor support ∩ high).card = 1 := by
    have hfilter :
        neighbor support ∩ high =
          (neighbor support).filter fun runner =>
            rightDegreeM neighbor runner = 4 := by
      ext runner
      simp [high]
    have hlocal :
        9 = 2 * (neighbor support).card + (neighbor support ∩ high).card := by
      rw [← hcompanion support]
      calc
        ∑ runner ∈ neighbor support,
            (rightDegreeM neighbor runner - 1) =
            ∑ runner ∈ neighbor support,
              (2 + if rightDegreeM neighbor runner = 4 then 1 else 0) := by
          exact Finset.sum_congr rfl fun runner _ => hsubExpand runner
        _ = 2 * (neighbor support).card + (neighbor support ∩ high).card := by
          rw [hfilter]
          simp [Finset.sum_add_distrib, Nat.mul_comm]
    have hle : (neighbor support ∩ high).card ≤ 2 := by
      exact (Finset.card_le_card (Finset.inter_subset_right)).trans_eq hhighCard
    exact eq_one_of_nine_eq_two_mul_add _ _ hlocal hle
  let highNeighbor : Fin 10 → Finset (Fin 13) :=
    fun support => neighbor support ∩ high
  have hleft : ∑ support, (highNeighbor support).card = 10 := by
    calc
      ∑ support, (highNeighbor support).card = ∑ _support : Fin 10, 1 := by
        exact Finset.sum_congr rfl fun support _ => hhighIncident support
      _ = 10 := by simp
  have hrestrictedDegree (runner : Fin 13) :
      rightDegreeM highNeighbor runner =
        if runner ∈ high then rightDegreeM neighbor runner else 0 := by
    by_cases hrunner : runner ∈ high <;>
      simp [rightDegreeM, highNeighbor, hrunner]
  have hright : ∑ runner, rightDegreeM highNeighbor runner = 8 := by
    calc
      ∑ runner, rightDegreeM highNeighbor runner =
          ∑ runner,
            if runner ∈ high then rightDegreeM neighbor runner else 0 := by
        exact Finset.sum_congr rfl fun runner _ => hrestrictedDegree runner
      _ = ∑ runner ∈ high, rightDegreeM neighbor runner := by
        simp
      _ = ∑ _runner ∈ high, 4 := by
        apply Finset.sum_congr rfl
        intro runner hrunner
        exact (Finset.mem_filter.mp hrunner).2
      _ = 8 := by simp [hhighCard]
  have hdouble := total_incidence_eq_sum_rightDegreeM highNeighbor
  omega

theorem z_ten_thirteen_le_40 (neighbor : Fin 10 → Finset (Fin 13))
    (hK22 : ∀ pair ∈ (Finset.univ : Finset (Fin 10)).powersetCard 2,
      (commonRunnersM neighbor pair).card ≤ 1) :
    ∑ support, (neighbor support).card ≤ 40 := by
  have hle := z_ten_thirteen_le_41 neighbor hK22
  have hne := z_ten_thirteen_ne_41 neighbor hK22
  omega

theorem z_eleven_thirteen_le_44 (neighbor : Fin 11 → Finset (Fin 13))
    (hK22 : ∀ pair ∈ (Finset.univ : Finset (Fin 11)).powersetCard 2,
      (commonRunnersM neighbor pair).card ≤ 1) :
    ∑ support, (neighbor support).card ≤ 44 := by
  have h := three_mul_incidence_le_78_add_choose
    (by decide : 11 ≤ 13) neighbor hK22
  norm_num [Nat.choose] at h ⊢
  omega

theorem z_twelve_thirteen_le_48 (neighbor : Fin 12 → Finset (Fin 13))
    (hK22 : ∀ pair ∈ (Finset.univ : Finset (Fin 12)).powersetCard 2,
      (commonRunnersM neighbor pair).card ≤ 1) :
    ∑ support, (neighbor support).card ≤ 48 := by
  have h := three_mul_incidence_le_78_add_choose
    (by decide : 12 ≤ 13) neighbor hK22
  norm_num [Nat.choose] at h ⊢
  omega

theorem z_thirteen_thirteen_le_52 (neighbor : Fin 13 → Finset (Fin 13))
    (hK22 : ∀ pair ∈ (Finset.univ : Finset (Fin 13)).powersetCard 2,
      (commonRunnersM neighbor pair).card ≤ 1) :
    ∑ support, (neighbor support).card ≤ 52 := by
  have h := three_mul_incidence_le_78_add_choose
    (by decide : 13 ≤ 13) neighbor hK22
  norm_num [Nat.choose] at h ⊢
  omega

/-! ## Sharp collision contrapositives -/

/-- Forty-one incidences on ten left supports force one unordered support
pair to share at least two right vertices.  This is purely an incidence
statement; it retains no coefficient, sign, or phase data. -/
theorem exists_repeated_pair_of_41_le_ten_thirteen
    (neighbor : Fin 10 → Finset (Fin 13))
    (hlower : 41 ≤ ∑ support, (neighbor support).card) :
    ∃ pair ∈ (Finset.univ : Finset (Fin 10)).powersetCard 2,
      2 ≤ (commonRunnersM neighbor pair).card := by
  by_contra hcollision
  push Not at hcollision
  have hK22 : ∀ pair ∈ (Finset.univ : Finset (Fin 10)).powersetCard 2,
      (commonRunnersM neighbor pair).card ≤ 1 := by
    intro pair hpair
    have hlt := hcollision pair hpair
    omega
  have hupper := z_ten_thirteen_le_40 neighbor hK22
  omega

/-- Forty-five incidences on eleven left supports force one unordered support
pair to share at least two right vertices. -/
theorem exists_repeated_pair_of_45_le_eleven_thirteen
    (neighbor : Fin 11 → Finset (Fin 13))
    (hlower : 45 ≤ ∑ support, (neighbor support).card) :
    ∃ pair ∈ (Finset.univ : Finset (Fin 11)).powersetCard 2,
      2 ≤ (commonRunnersM neighbor pair).card := by
  by_contra hcollision
  push Not at hcollision
  have hK22 : ∀ pair ∈ (Finset.univ : Finset (Fin 11)).powersetCard 2,
      (commonRunnersM neighbor pair).card ≤ 1 := by
    intro pair hpair
    have hlt := hcollision pair hpair
    omega
  have hupper := z_eleven_thirteen_le_44 neighbor hK22
  omega

/-- Forty-nine incidences on twelve left supports force one unordered support
pair to share at least two right vertices. -/
theorem exists_repeated_pair_of_49_le_twelve_thirteen
    (neighbor : Fin 12 → Finset (Fin 13))
    (hlower : 49 ≤ ∑ support, (neighbor support).card) :
    ∃ pair ∈ (Finset.univ : Finset (Fin 12)).powersetCard 2,
      2 ≤ (commonRunnersM neighbor pair).card := by
  by_contra hcollision
  push Not at hcollision
  have hK22 : ∀ pair ∈ (Finset.univ : Finset (Fin 12)).powersetCard 2,
      (commonRunnersM neighbor pair).card ≤ 1 := by
    intro pair hpair
    have hlt := hcollision pair hpair
    omega
  have hupper := z_twelve_thirteen_le_48 neighbor hK22
  omega

/-- Fifty-three incidences on thirteen left supports force one unordered
support pair to share at least two right vertices. -/
theorem exists_repeated_pair_of_53_le_thirteen_thirteen
    (neighbor : Fin 13 → Finset (Fin 13))
    (hlower : 53 ≤ ∑ support, (neighbor support).card) :
    ∃ pair ∈ (Finset.univ : Finset (Fin 13)).powersetCard 2,
      2 ≤ (commonRunnersM neighbor pair).card := by
  by_contra hcollision
  push Not at hcollision
  have hK22 : ∀ pair ∈ (Finset.univ : Finset (Fin 13)).powersetCard 2,
      (commonRunnersM neighbor pair).card ≤ 1 := by
    intro pair hpair
    have hlt := hcollision pair hpair
    omega
  have hupper := z_thirteen_thirteen_le_52 neighbor hK22
  omega

/-! ## Mixed-stratum zero-slack consumers -/

/-- Among ten supports of size at least four, one support of size at least five
forces a repeated right-vertex pair. -/
theorem exists_repeated_pair_of_four_le_of_exists_five_ten_thirteen
    (neighbor : Fin 10 → Finset (Fin 13))
    (hfour : ∀ support, 4 ≤ (neighbor support).card)
    (hfive : ∃ support, 5 ≤ (neighbor support).card) :
    ∃ pair ∈ (Finset.univ : Finset (Fin 10)).powersetCard 2,
      2 ≤ (commonRunnersM neighbor pair).card := by
  obtain ⟨special, hspecial⟩ := hfive
  have hstrict :
      (∑ _support : Fin 10, 4) <
        ∑ support, (neighbor support).card := by
    exact Finset.sum_lt_sum (fun support _ => hfour support)
      ⟨special, Finset.mem_univ special, by omega⟩
  have hbaseline : (∑ _support : Fin 10, 4) = 40 := by decide
  rw [hbaseline] at hstrict
  exact exists_repeated_pair_of_41_le_ten_thirteen neighbor (by omega)

/-- Among eleven supports of size at least four, one support of size at least
five forces a repeated right-vertex pair. -/
theorem exists_repeated_pair_of_four_le_of_exists_five_eleven_thirteen
    (neighbor : Fin 11 → Finset (Fin 13))
    (hfour : ∀ support, 4 ≤ (neighbor support).card)
    (hfive : ∃ support, 5 ≤ (neighbor support).card) :
    ∃ pair ∈ (Finset.univ : Finset (Fin 11)).powersetCard 2,
      2 ≤ (commonRunnersM neighbor pair).card := by
  obtain ⟨special, hspecial⟩ := hfive
  have hstrict :
      (∑ _support : Fin 11, 4) <
        ∑ support, (neighbor support).card := by
    exact Finset.sum_lt_sum (fun support _ => hfour support)
      ⟨special, Finset.mem_univ special, by omega⟩
  have hbaseline : (∑ _support : Fin 11, 4) = 44 := by decide
  rw [hbaseline] at hstrict
  exact exists_repeated_pair_of_45_le_eleven_thirteen neighbor (by omega)

/-- Among twelve supports of size at least four, one support of size at least
five forces a repeated right-vertex pair. -/
theorem exists_repeated_pair_of_four_le_of_exists_five_twelve_thirteen
    (neighbor : Fin 12 → Finset (Fin 13))
    (hfour : ∀ support, 4 ≤ (neighbor support).card)
    (hfive : ∃ support, 5 ≤ (neighbor support).card) :
    ∃ pair ∈ (Finset.univ : Finset (Fin 12)).powersetCard 2,
      2 ≤ (commonRunnersM neighbor pair).card := by
  obtain ⟨special, hspecial⟩ := hfive
  have hstrict :
      (∑ _support : Fin 12, 4) <
        ∑ support, (neighbor support).card := by
    exact Finset.sum_lt_sum (fun support _ => hfour support)
      ⟨special, Finset.mem_univ special, by omega⟩
  have hbaseline : (∑ _support : Fin 12, 4) = 48 := by decide
  rw [hbaseline] at hstrict
  exact exists_repeated_pair_of_49_le_twelve_thirteen neighbor (by omega)

/-- Among thirteen supports of size at least four, one support of size at
least five forces a repeated right-vertex pair. -/
theorem exists_repeated_pair_of_four_le_of_exists_five_thirteen_thirteen
    (neighbor : Fin 13 → Finset (Fin 13))
    (hfour : ∀ support, 4 ≤ (neighbor support).card)
    (hfive : ∃ support, 5 ≤ (neighbor support).card) :
    ∃ pair ∈ (Finset.univ : Finset (Fin 13)).powersetCard 2,
      2 ≤ (commonRunnersM neighbor pair).card := by
  obtain ⟨special, hspecial⟩ := hfive
  have hstrict :
      (∑ _support : Fin 13, 4) <
        ∑ support, (neighbor support).card := by
    exact Finset.sum_lt_sum (fun support _ => hfour support)
      ⟨special, Finset.mem_univ special, by omega⟩
  have hbaseline : (∑ _support : Fin 13, 4) = 52 := by decide
  rw [hbaseline] at hstrict
  exact exists_repeated_pair_of_53_le_thirteen_thirteen neighbor (by omega)

/-! ## Sharp finite-family mixed-stratum bound -/

/-- Reindexing a pair-unique finite support family preserves the statement
that every unordered pair of indices has at most one common runner. -/
theorem commonRunnersM_enumeration_card_le_one {m : ℕ}
    (supports : Finset (Finset (Fin 13)))
    (enumeration : Fin m ≃ ↥supports)
    (hunique : PairUnique supports)
    (pair : Finset (Fin m))
    (hpair : pair ∈ (Finset.univ : Finset (Fin m)).powersetCard 2) :
    (commonRunnersM
      (fun index => (enumeration index).1) pair).card ≤ 1 := by
  let neighbor : Fin m → Finset (Fin 13) :=
    fun index => (enumeration index).1
  have hneighborMem (index : Fin m) : neighbor index ∈ supports := by
    exact (enumeration index).2
  have hneighborInj : Function.Injective neighbor := by
    exact Subtype.val_injective.comp enumeration.injective
  change (commonRunnersM neighbor pair).card ≤ 1
  rw [Finset.card_le_one_iff]
  intro firstRunner secondRunner hfirst hsecond
  by_contra hrunnerNe
  have hpairCard : pair.card = 2 :=
    (Finset.mem_powersetCard.mp hpair).2
  obtain ⟨firstIndex, secondIndex, hindexNe, rfl⟩ :=
    Finset.card_eq_two.mp hpairCard
  have hfirstMembership :
      firstRunner ∈ neighbor firstIndex ∧
        firstRunner ∈ neighbor secondIndex := by
    simpa [commonRunnersM, hindexNe] using hfirst
  have hsecondMembership :
      secondRunner ∈ neighbor firstIndex ∧
        secondRunner ∈ neighbor secondIndex := by
    simpa [commonRunnersM, hindexNe] using hsecond
  have hsupportNe : neighbor firstIndex ≠ neighbor secondIndex :=
    hneighborInj.ne hindexNe
  have hdisjoint :
      Disjoint (supportPairs (neighbor firstIndex))
        (supportPairs (neighbor secondIndex)) :=
    hunique (hneighborMem firstIndex) (hneighborMem secondIndex) hsupportNe
  have hownedFirst :
      {firstRunner, secondRunner} ∈ supportPairs (neighbor firstIndex) := by
    rw [supportPairs, Finset.mem_powersetCard]
    refine ⟨?_, by simp [hrunnerNe]⟩
    intro runner hrunner
    simp only [Finset.mem_insert, Finset.mem_singleton] at hrunner
    rcases hrunner with rfl | rfl
    · exact hfirstMembership.1
    · exact hsecondMembership.1
  have hownedSecond :
      {firstRunner, secondRunner} ∈ supportPairs (neighbor secondIndex) := by
    rw [supportPairs, Finset.mem_powersetCard]
    refine ⟨?_, by simp [hrunnerNe]⟩
    intro runner hrunner
    simp only [Finset.mem_insert, Finset.mem_singleton] at hrunner
    rcases hrunner with rfl | rfl
    · exact hfirstMembership.2
    · exact hsecondMembership.2
  exact (Finset.disjoint_left.mp hdisjoint hownedFirst) hownedSecond

/-- A pair-unique family of supports on thirteen runners, with every support
of size at least four and at least one support of size at least five, has at
most nine members.  Existing `witness9` shows that the numerical bound is
sharp at the indexed-incidence level. -/
theorem card_le_nine_of_four_le_of_exists_five
    (supports : Finset (Finset (Fin 13)))
    (hunique : PairUnique supports)
    (hfour : ∀ support ∈ supports, 4 ≤ support.card)
    (hfive : ∃ support ∈ supports, 5 ≤ support.card) :
    supports.card ≤ 9 := by
  by_contra hnot
  have hten : 10 ≤ supports.card := by omega
  have hthirteen : supports.card ≤ 13 :=
    card_le_13_of_four_le supports hunique hfour
  have hcases :
      supports.card = 10 ∨ supports.card = 11 ∨
        supports.card = 12 ∨ supports.card = 13 := by
    omega
  rcases hcases with hcard | hcard | hcard | hcard
  · let enumeration : Fin 10 ≃ ↥supports :=
      (finCongr hcard.symm).trans supports.equivFin.symm
    let neighbor : Fin 10 → Finset (Fin 13) :=
      fun index => (enumeration index).1
    have hneighborFour (index : Fin 10) :
        4 ≤ (neighbor index).card :=
      hfour _ (enumeration index).2
    have hneighborFive : ∃ index, 5 ≤ (neighbor index).card := by
      obtain ⟨special, hspecialMem, hspecialSize⟩ := hfive
      obtain ⟨index, hindex⟩ :=
        enumeration.surjective ⟨special, hspecialMem⟩
      have hvalue : neighbor index = special := congrArg Subtype.val hindex
      exact ⟨index, by simpa [hvalue] using hspecialSize⟩
    obtain ⟨pair, hpair, hcommon⟩ :=
      exists_repeated_pair_of_four_le_of_exists_five_ten_thirteen
        neighbor hneighborFour hneighborFive
    have hle := commonRunnersM_enumeration_card_le_one
      supports enumeration hunique pair hpair
    change (commonRunnersM neighbor pair).card ≤ 1 at hle
    omega
  · let enumeration : Fin 11 ≃ ↥supports :=
      (finCongr hcard.symm).trans supports.equivFin.symm
    let neighbor : Fin 11 → Finset (Fin 13) :=
      fun index => (enumeration index).1
    have hneighborFour (index : Fin 11) :
        4 ≤ (neighbor index).card :=
      hfour _ (enumeration index).2
    have hneighborFive : ∃ index, 5 ≤ (neighbor index).card := by
      obtain ⟨special, hspecialMem, hspecialSize⟩ := hfive
      obtain ⟨index, hindex⟩ :=
        enumeration.surjective ⟨special, hspecialMem⟩
      have hvalue : neighbor index = special := congrArg Subtype.val hindex
      exact ⟨index, by simpa [hvalue] using hspecialSize⟩
    obtain ⟨pair, hpair, hcommon⟩ :=
      exists_repeated_pair_of_four_le_of_exists_five_eleven_thirteen
        neighbor hneighborFour hneighborFive
    have hle := commonRunnersM_enumeration_card_le_one
      supports enumeration hunique pair hpair
    change (commonRunnersM neighbor pair).card ≤ 1 at hle
    omega
  · let enumeration : Fin 12 ≃ ↥supports :=
      (finCongr hcard.symm).trans supports.equivFin.symm
    let neighbor : Fin 12 → Finset (Fin 13) :=
      fun index => (enumeration index).1
    have hneighborFour (index : Fin 12) :
        4 ≤ (neighbor index).card :=
      hfour _ (enumeration index).2
    have hneighborFive : ∃ index, 5 ≤ (neighbor index).card := by
      obtain ⟨special, hspecialMem, hspecialSize⟩ := hfive
      obtain ⟨index, hindex⟩ :=
        enumeration.surjective ⟨special, hspecialMem⟩
      have hvalue : neighbor index = special := congrArg Subtype.val hindex
      exact ⟨index, by simpa [hvalue] using hspecialSize⟩
    obtain ⟨pair, hpair, hcommon⟩ :=
      exists_repeated_pair_of_four_le_of_exists_five_twelve_thirteen
        neighbor hneighborFour hneighborFive
    have hle := commonRunnersM_enumeration_card_le_one
      supports enumeration hunique pair hpair
    change (commonRunnersM neighbor pair).card ≤ 1 at hle
    omega
  · let enumeration : Fin 13 ≃ ↥supports :=
      (finCongr hcard.symm).trans supports.equivFin.symm
    let neighbor : Fin 13 → Finset (Fin 13) :=
      fun index => (enumeration index).1
    have hneighborFour (index : Fin 13) :
        4 ≤ (neighbor index).card :=
      hfour _ (enumeration index).2
    have hneighborFive : ∃ index, 5 ≤ (neighbor index).card := by
      obtain ⟨special, hspecialMem, hspecialSize⟩ := hfive
      obtain ⟨index, hindex⟩ :=
        enumeration.surjective ⟨special, hspecialMem⟩
      have hvalue : neighbor index = special := congrArg Subtype.val hindex
      exact ⟨index, by simpa [hvalue] using hspecialSize⟩
    obtain ⟨pair, hpair, hcommon⟩ :=
      exists_repeated_pair_of_four_le_of_exists_five_thirteen_thirteen
        neighbor hneighborFour hneighborFive
    have hle := commonRunnersM_enumeration_card_le_one
      supports enumeration hunique pair hpair
    change (commonRunnersM neighbor pair).card ≤ 1 at hle
    omega

/-! ## Cyclic Singer witnesses -/

/-- A cyclic line of the order-three projective plane, using the difference
set `{0,1,3,9}` in `Z/13Z`. -/
def singerBlock (index : Fin 13) : Finset (Fin 13) :=
  {index, index + 1, index + 3, index + 9}

/-- The thirteen cyclic translates of `{0,1,3,9}`, written extensionally so
the finite K22 checks reduce without normalizing modular arithmetic. -/
def witness13 : Fin 13 → Finset (Fin 13) := ![
  {0, 1, 3, 9},
  {1, 2, 4, 10},
  {2, 3, 5, 11},
  {3, 4, 6, 12},
  {0, 4, 5, 7},
  {1, 5, 6, 8},
  {2, 6, 7, 9},
  {3, 7, 8, 10},
  {4, 8, 9, 11},
  {5, 9, 10, 12},
  {0, 6, 10, 11},
  {1, 7, 11, 12},
  {0, 2, 8, 12}
]

def singerWitness {m : ℕ} (hm : m ≤ 13) : Fin m → Finset (Fin 13) :=
  fun support => witness13 (Fin.castLE hm support)

def witness10 : Fin 10 → Finset (Fin 13) := singerWitness (by decide)
def witness11 : Fin 11 → Finset (Fin 13) := singerWitness (by decide)
def witness12 : Fin 12 → Finset (Fin 13) := singerWitness (by decide)

theorem witness13_pairwise :
    ∀ first second : Fin 13, first ≠ second →
      (witness13 first ∩ witness13 second).card ≤ 1 := by
  set_option maxRecDepth 10000 in decide

/-- The four Singer lines through a point.  The point is a right vertex in the
incidence graph and the returned indices are the four incident left blocks. -/
def singerSupportsThrough (runner : Fin 13) : Finset (Fin 13) :=
  Finset.univ.filter fun support => runner ∈ witness13 support

/-- Delete the common point from one of its incident Singer lines. -/
def singerTail (runner support : Fin 13) : Finset (Fin 13) :=
  (witness13 support).erase runner

/-- Local parallel-class form of the Singer plane: at every point, the four
incident lines leave four disjoint three-point tails partitioning the other
twelve points.  This is the saturated counterpart of the impossible local
partition forced by a hypothetical forty-one-edge ten-point design. -/
theorem singer_parallel_class_partition :
    ∀ runner : Fin 13,
      (singerSupportsThrough runner).card = 4 ∧
      (∀ support ∈ singerSupportsThrough runner,
        (singerTail runner support).card = 3) ∧
      (∀ first ∈ singerSupportsThrough runner,
        ∀ second ∈ singerSupportsThrough runner, first ≠ second →
          Disjoint (singerTail runner first) (singerTail runner second)) ∧
      (singerSupportsThrough runner).biUnion (singerTail runner) =
        (Finset.univ : Finset (Fin 13)).erase runner := by
  set_option maxRecDepth 10000 in decide

theorem singerWitness_K22 {m : ℕ} (hm : m ≤ 13) :
    ∀ pair ∈ (Finset.univ : Finset (Fin m)).powersetCard 2,
      (commonRunnersM (singerWitness hm) pair).card ≤ 1 := by
  apply k22_of_pairwise_inter_card_le_one
  intro first second hne
  exact witness13_pairwise _ _ ((Fin.castLE_injective hm).ne hne)

theorem witness10_K22 :
    ∀ pair ∈ (Finset.univ : Finset (Fin 10)).powersetCard 2,
      (commonRunnersM witness10 pair).card ≤ 1 :=
  singerWitness_K22 (by decide)

theorem witness11_K22 :
    ∀ pair ∈ (Finset.univ : Finset (Fin 11)).powersetCard 2,
      (commonRunnersM witness11 pair).card ≤ 1 :=
  singerWitness_K22 (by decide)

theorem witness12_K22 :
    ∀ pair ∈ (Finset.univ : Finset (Fin 12)).powersetCard 2,
      (commonRunnersM witness12 pair).card ≤ 1 :=
  singerWitness_K22 (by decide)

theorem witness13_K22 :
    ∀ pair ∈ (Finset.univ : Finset (Fin 13)).powersetCard 2,
      (commonRunnersM witness13 pair).card ≤ 1 :=
  k22_of_pairwise_inter_card_le_one witness13 witness13_pairwise

theorem witness10_edges : ∑ support, (witness10 support).card = 40 := by decide
theorem witness11_edges : ∑ support, (witness11 support).card = 44 := by decide
theorem witness12_edges : ∑ support, (witness12 support).card = 48 := by decide
theorem witness13_edges : ∑ support, (witness13 support).card = 52 := by decide

/-! ## Representation-level sharpness packages -/

theorem ten_thirteen_upper_and_witness :
    (∀ neighbor : Fin 10 → Finset (Fin 13),
      (∀ pair ∈ (Finset.univ : Finset (Fin 10)).powersetCard 2,
        (commonRunnersM neighbor pair).card ≤ 1) →
      ∑ support, (neighbor support).card ≤ 40) ∧
    (∀ pair ∈ (Finset.univ : Finset (Fin 10)).powersetCard 2,
      (commonRunnersM witness10 pair).card ≤ 1) ∧
    ∑ support, (witness10 support).card = 40 := by
  exact ⟨z_ten_thirteen_le_40, witness10_K22, witness10_edges⟩

theorem eleven_thirteen_upper_and_witness :
    (∀ neighbor : Fin 11 → Finset (Fin 13),
      (∀ pair ∈ (Finset.univ : Finset (Fin 11)).powersetCard 2,
        (commonRunnersM neighbor pair).card ≤ 1) →
      ∑ support, (neighbor support).card ≤ 44) ∧
    (∀ pair ∈ (Finset.univ : Finset (Fin 11)).powersetCard 2,
      (commonRunnersM witness11 pair).card ≤ 1) ∧
    ∑ support, (witness11 support).card = 44 := by
  exact ⟨z_eleven_thirteen_le_44, witness11_K22, witness11_edges⟩

theorem twelve_thirteen_upper_and_witness :
    (∀ neighbor : Fin 12 → Finset (Fin 13),
      (∀ pair ∈ (Finset.univ : Finset (Fin 12)).powersetCard 2,
        (commonRunnersM neighbor pair).card ≤ 1) →
      ∑ support, (neighbor support).card ≤ 48) ∧
    (∀ pair ∈ (Finset.univ : Finset (Fin 12)).powersetCard 2,
      (commonRunnersM witness12 pair).card ≤ 1) ∧
    ∑ support, (witness12 support).card = 48 := by
  exact ⟨z_twelve_thirteen_le_48, witness12_K22, witness12_edges⟩

theorem thirteen_thirteen_upper_and_witness :
    (∀ neighbor : Fin 13 → Finset (Fin 13),
      (∀ pair ∈ (Finset.univ : Finset (Fin 13)).powersetCard 2,
        (commonRunnersM neighbor pair).card ≤ 1) →
      ∑ support, (neighbor support).card ≤ 52) ∧
    (∀ pair ∈ (Finset.univ : Finset (Fin 13)).powersetCard 2,
      (commonRunnersM witness13 pair).card ≤ 1) ∧
    ∑ support, (witness13 support).card = 52 := by
  exact ⟨z_thirteen_thirteen_le_52, witness13_K22, witness13_edges⟩

/-! ## Axiom audit -/

#print axioms z_ten_thirteen_le_40
#print axioms z_eleven_thirteen_le_44
#print axioms z_twelve_thirteen_le_48
#print axioms z_thirteen_thirteen_le_52
#print axioms witness10_K22
#print axioms witness13_K22
#print axioms singer_parallel_class_partition
#print axioms ten_thirteen_upper_and_witness
#print axioms exists_repeated_pair_of_41_le_ten_thirteen
#print axioms exists_repeated_pair_of_53_le_thirteen_thirteen
#print axioms exists_repeated_pair_of_four_le_of_exists_five_ten_thirteen
#print axioms exists_repeated_pair_of_four_le_of_exists_five_thirteen_thirteen
#print axioms card_le_nine_of_four_le_of_exists_five

end LRCZarankiewiczTenThirteen
end LonelyRunner
