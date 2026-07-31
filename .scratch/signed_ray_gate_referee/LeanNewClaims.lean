import TournamentH7.LRCGridAttainment

open Filter Topology
open TournamentH7.LRCWitness

namespace LonelyRunner.LRC14.SignedRayReferee

/-- Circular distance is subadditive. -/
theorem distZ_sub_le_add (x y : ℝ) :
    distZ (x - y) ≤ distZ x + distZ y := by
  have hnear :=
    LonelyRunner.GridAttainment.distZ_le_abs_sub
      (x - y) (round x - round y)
  have hrearrange :
      (x - y) - (((round x - round y : ℤ) : ℝ)) =
        (x - (round x : ℝ)) - (y - (round y : ℝ)) := by
    push_cast
    ring
  rw [hrearrange] at hnear
  calc
    distZ (x - y)
        ≤ |(x - (round x : ℝ)) - (y - (round y : ℝ))| := hnear
    _ ≤ |x - (round x : ℝ)| + |(round y : ℝ) - y| := by
      simpa [sub_zero, zero_sub, abs_neg] using
        (abs_sub_le (x - (round x : ℝ)) 0 (y - (round y : ℝ)))
    _ = |x - (round x : ℝ)| + |y - (round y : ℝ)| := by
      rw [abs_sub_comm (round y : ℝ) y]
    _ = distZ x + distZ y := by
      rw [LonelyRunner.GridAttainment.distZ_eq_round,
        LonelyRunner.GridAttainment.distZ_eq_round]

theorem not_both_strict_bad_of_separation (x y : ℝ)
    (hsep : (1 : ℝ) / 7 ≤ distZ (x - y)) :
    ¬ (distZ x < (1 : ℝ) / 14 ∧ distZ y < (1 : ℝ) / 14) := by
  intro hbad
  have htri := distZ_sub_le_add x y
  linarith

theorem escape_of_equiv_bad_card_lt {α : Type*} [DecidableEq α]
    (S B : Finset α) (e : α ≃ α) (hcard : B.card < S.card) :
    ∃ x ∈ S, e x ∉ B := by
  by_contra hno
  have hall : ∀ x ∈ S, e x ∈ B := by
    intro x hx
    by_contra hbad
    exact hno ⟨x, hx, hbad⟩
  have hsub : S.map e.toEmbedding ⊆ B := by
    intro y hy
    rw [Finset.mem_map] at hy
    obtain ⟨x, hx, rfl⟩ := hy
    exact hall x hx
  have hle := Finset.card_le_card hsub
  rw [Finset.card_map] at hle
  omega

theorem escape_of_invariant_stratum_card_lt {α : Type*} [DecidableEq α]
    (S B O : Finset α) (e : α ≃ α)
    (hinv : ∀ x, x ∈ O ↔ e x ∈ O)
    (hcard : (B ∩ O).card < (S ∩ O).card) :
    ∃ x ∈ S, e x ∉ B := by
  by_contra hno
  have hall : ∀ x ∈ S, e x ∈ B := by
    intro x hx
    by_contra hbad
    exact hno ⟨x, hx, hbad⟩
  have hsub : (S ∩ O).map e.toEmbedding ⊆ B ∩ O := by
    intro y hy
    rw [Finset.mem_map] at hy
    obtain ⟨x, hx, rfl⟩ := hy
    have hxS : x ∈ S := (Finset.mem_inter.mp hx).1
    have hxO : x ∈ O := (Finset.mem_inter.mp hx).2
    exact Finset.mem_inter.mpr ⟨hall x hxS, (hinv x).mp hxO⟩
  have hle := Finset.card_le_card hsub
  rw [Finset.card_map] at hle
  omega

def beta (d : ℕ) : ℕ := 2 * ((d - 1) / 14) + 1

/-- The unit-count threshold never exceeds the order-seven density scale. -/
theorem beta_le_ceil_seventh (d : ℕ) (hd : 2 ≤ d) :
    beta d ≤ (d + 6) / 7 := by
  unfold beta
  omega

/-- Arithmetic core of the translated open-interval count.  `span` is the
integer distance between the first and last occupied lattice sites;
`hspacing` is unit spacing and `hopen` is strict containment in an interval
of length `d/7`. -/
theorem strict_lattice_interval_capacity (d k span : ℕ)
    (hspacing : k ≤ span + 1) (hopen : 7 * span < d) :
    k ≤ (d + 6) / 7 := by
  omega

/-- The arithmetic bound is sharp: the target number of lattice sites has a
unit-spaced span still strictly shorter than `d/7`. -/
theorem strict_lattice_interval_capacity_sharp (d : ℕ) (hd : 0 < d) :
    ∃ span : ℕ, (d + 6) / 7 ≤ span + 1 ∧ 7 * span < d := by
  refine ⟨(d + 6) / 7 - 1, ?_, ?_⟩ <;> omega

def badBand (d : ℕ) : Finset ℕ :=
  (Finset.range d).filter (fun r => 14 * min r (d - r) < d)

lemma badBand_eq_interval_union (d : ℕ) (hd : 2 ≤ d) :
    badBand d =
      Finset.range ((d - 1) / 14 + 1) ∪
        Finset.Ico (d - (d - 1) / 14) d := by
  ext r
  simp only [badBand, Finset.mem_filter, Finset.mem_range,
    Finset.mem_union, Finset.mem_Ico]
  constructor
  · rintro ⟨hr, hbad⟩
    have hnear : min r (d - r) ≤ (d - 1) / 14 := by
      omega
    rcases (min_le_iff).mp hnear with hleft | hright
    · exact Or.inl (by omega)
    · exact Or.inr ⟨by omega, hr⟩
  · rintro (hleft | hright)
    · exact ⟨by omega, by
        have hnear : min r (d - r) ≤ (d - 1) / 14 :=
          le_trans (min_le_left _ _) (by omega)
        omega⟩
    · exact ⟨hright.2, by
        have hnear : min r (d - r) ≤ (d - 1) / 14 :=
          le_trans (min_le_right _ _) (by omega)
        omega⟩

theorem badBand_card (d : ℕ) (hd : 2 ≤ d) :
    (badBand d).card = beta d := by
  let b := (d - 1) / 14
  have hb : 2 * b < d := by
    dsimp [b]
    omega
  have hdisj : Disjoint (Finset.range (b + 1)) (Finset.Ico (d - b) d) := by
    rw [Finset.disjoint_left]
    intro r hrange hico
    simp only [Finset.mem_range] at hrange
    simp only [Finset.mem_Ico] at hico
    omega
  rw [badBand_eq_interval_union d hd]
  change (Finset.range (b + 1) ∪ Finset.Ico (d - b) d).card = beta d
  rw [Finset.card_union_of_disjoint hdisj, Finset.card_range, Nat.card_Ico]
  dsimp [beta, b]
  omega

def badBandFin (d : ℕ) : Finset (Fin d) :=
  Finset.univ.filter (fun r => 14 * min r.val (d - r.val) < d)

lemma badBandFin_map (d : ℕ) :
    (badBandFin d).map Fin.valEmbedding = badBand d := by
  ext r
  constructor
  · intro hr
    rw [Finset.mem_map] at hr
    obtain ⟨x, hx, rfl⟩ := hr
    simp only [badBandFin, Finset.mem_filter, Finset.mem_univ, true_and] at hx
    exact Finset.mem_filter.mpr ⟨Finset.mem_range.mpr x.isLt, hx⟩
  · intro hr
    obtain ⟨hrd, hbad⟩ := Finset.mem_filter.mp hr
    let x : Fin d := ⟨r, Finset.mem_range.mp hrd⟩
    rw [Finset.mem_map]
    exact ⟨x, by simpa [badBandFin] using hbad, rfl⟩

theorem badBandFin_card (d : ℕ) (hd : 2 ≤ d) :
    (badBandFin d).card = beta d := by
  have hmap := congrArg Finset.card (badBandFin_map d)
  simpa [badBand_card d hd] using hmap

theorem escape_badBandFin_of_beta_lt (d : ℕ) (hd : 2 ≤ d)
    (S : Finset (Fin d)) (e : Equiv.Perm (Fin d))
    (hcard : beta d < S.card) :
    ∃ x ∈ S, e x ∉ badBandFin d := by
  apply escape_of_equiv_bad_card_lt S (badBandFin d) e
  rw [badBandFin_card d hd]
  exact hcard

theorem equivalence_bad_preimage_sharp {α : Type*} [DecidableEq α]
    (B : Finset α) (e : Equiv.Perm α) :
    ∃ S : Finset α, S.card = B.card ∧ ∀ x ∈ S, e x ∈ B := by
  refine ⟨B.map e.symm.toEmbedding, by simp, ?_⟩
  intro x hx
  rw [Finset.mem_map] at hx
  obtain ⟨y, hy, hxy⟩ := hx
  rw [← hxy]
  change e (e.symm y) ∈ B
  simpa using hy

/-! ## Arithmetic cores for the complete affine flag classification -/

/-- A clique with more vertices than the common facet capacity cannot itself
be a face contained in one such facet. -/
theorem not_subset_of_card_gt_capacity {α : Type*} [DecidableEq α]
    (S C : Finset α) (kappa : ℕ) (hC : C.card ≤ kappa)
    (hS : kappa < S.card) : ¬ S ⊆ C := by
  intro hsub
  have hle := Finset.card_le_card hsub
  omega

/-- The equality cases for `R=1` in
`q = ceil(Rq/7)` are only `q=0,1`. -/
theorem ceil_seventh_R1_eq (q : ℕ) :
    q = (q + 6) / 7 ↔ q ≤ 1 := by
  omega

/-- The equality cases for `R=2` in `q = ceil(Rq/7)`. -/
theorem ceil_seventh_R2_eq (q : ℕ) :
    q = (2 * q + 6) / 7 ↔ q ≤ 1 := by
  omega

/-- The equality cases for `R=3` in `q = ceil(Rq/7)`. -/
theorem ceil_seventh_R3_eq (q : ℕ) :
    q = (3 * q + 6) / 7 ↔ q ≤ 1 := by
  omega

/-- The equality cases for `R=4` in `q = ceil(Rq/7)`. -/
theorem ceil_seventh_R4_eq (q : ℕ) :
    q = (4 * q + 6) / 7 ↔ q ≤ 2 := by
  omega

/-- The equality cases for `R=5` in `q = ceil(Rq/7)`. -/
theorem ceil_seventh_R5_eq (q : ℕ) :
    q = (5 * q + 6) / 7 ↔ q ≤ 3 := by
  omega

/-- The equality cases for `R=6` in `q = ceil(Rq/7)`. -/
theorem ceil_seventh_R6_eq (q : ℕ) :
    q = (6 * q + 6) / 7 ↔ q ≤ 6 := by
  omega

/-- Multiples of seven are always on the cardinality-equality boundary. -/
theorem ceil_seventh_R7_eq (q : ℕ) :
    q = (7 * q + 6) / 7 := by
  omega

/-- The pair counts `N_1=m-3`, `N_2=m-2` force index separations three and
two in any length-`m` arithmetic-progression facet. -/
theorem pair_counts_force_indices (m h1 h2 : ℕ)
    (hm : 3 ≤ m) (hh1 : h1 ≤ m) (hh2 : h2 ≤ m)
    (hN1 : m - h1 = m - 3) (hN2 : m - h2 = m - 2) :
    h1 = 3 ∧ h2 = 2 := by
  omega

/-- A natural modulus at least nine cannot divide a signed four or eight. -/
theorem not_dvd_signed_four_or_eight (d : ℕ) (hd : 9 ≤ d) (w : ℤ)
    (hw : w = 4 ∨ w = -4 ∨ w = 8 ∨ w = -8) :
    ¬ (d : ℤ) ∣ w := by
  intro hdvd
  rcases hw with rfl | rfl | rfl | rfl
  all_goals
    have hle := Int.natAbs_le_of_dvd_ne_zero hdvd (by norm_num)
    norm_num at hle
    omega

/-- Elimination core of the `N_1/N_2` obstruction.  If pair counts forced
`3v=+-1` and `2v=+-2` modulo `d`, then `d` would divide a signed four or
eight, impossible once `d>=9`. -/
theorem pair_count_sign_contradiction (d : ℕ) (hd : 9 ≤ d)
    (v eps1 eps2 : ℤ)
    (heps1 : eps1 = 1 ∨ eps1 = -1)
    (heps2 : eps2 = 1 ∨ eps2 = -1)
    (h1 : (d : ℤ) ∣ 3 * v - eps1)
    (h2 : (d : ℤ) ∣ 2 * v - 2 * eps2) : False := by
  have hsmall : (d : ℤ) ∣ 2 * eps1 - 6 * eps2 := by
    obtain ⟨k1, hk1⟩ := h1
    obtain ⟨k2, hk2⟩ := h2
    refine ⟨3 * k2 - 2 * k1, ?_⟩
    calc
      2 * eps1 - 6 * eps2
          = 3 * (2 * v - 2 * eps2) - 2 * (3 * v - eps1) := by ring
      _ = 3 * ((d : ℤ) * k2) - 2 * ((d : ℤ) * k1) := by rw [hk2, hk1]
      _ = (d : ℤ) * (3 * k2 - 2 * k1) := by ring
  apply not_dvd_signed_four_or_eight d hd (2 * eps1 - 6 * eps2) ?_ hsmall
  rcases heps1 with rfl | rfl <;> rcases heps2 with rfl | rfl <;> norm_num

end LonelyRunner.LRC14.SignedRayReferee
