/-
  TournamentH7.LRCRatWitness — the RATIONAL-WITNESS GATE (kind-pasteur-2026-07-02-S12,
  HYP-3969).  Single-writer satellite.

  THE CENSUS-LEG ENGINE: the bounded branch of the census+peel surface ultimately rests on
  concrete families with concrete lonely times.  For an INTEGER family `v` and a RATIONAL
  time `t`, loneliness is a FINITE rational computation: each speed needs one comparison
  `1/14 ≤ |v i · t − round (v i · t)|` (round-point minimality covers every other integer,
  since distinct integers are ≥ 1 apart and the round distance is ≤ 1/2 ≤ 1 − 1/14).

  `ratWitnessOK v t` bundles the 13 comparisons into ONE decidable predicate;
  `lonely_of_ratWitness` is its soundness (kernel-pure).  Any witness table
  (family, rational time) — e.g. produced by the Python census engines — becomes a
  `Lonely` theorem by a single `decide`/`native_decide` per row.
-/
import TournamentH7.LRC14CertRoute

namespace LonelyRunner
namespace LRC14

/-- One-speed rational loneliness check at threshold `1/n`: a single exact rational
comparison against the round point. -/
def ratLonelyAt (n : ℕ) (s : ℤ) (t : ℚ) : Prop :=
  (1 : ℚ) / n ≤ |(s : ℚ) * t - round ((s : ℚ) * t)|

instance (n : ℕ) (s : ℤ) (t : ℚ) : Decidable (ratLonelyAt n s t) := by
  unfold ratLonelyAt
  infer_instance

/-- Round-point minimality over ℚ: the check at the round point covers every integer
(`n ≥ 2`; distinct integers are ≥ 1 apart, the round distance is ≤ 1/2). -/
theorem ratLonelyAt_forall_int {n : ℕ} (hn : 2 ≤ n) {s : ℤ} {t : ℚ}
    (h : ratLonelyAt n s t) : ∀ m : ℤ, (1 : ℚ) / n ≤ |(s : ℚ) * t - m| := by
  intro m
  unfold ratLonelyAt at h
  set x : ℚ := (s : ℚ) * t with hx
  rcases eq_or_ne m (round x) with rfl | hne
  · exact h
  · have hr : |x - (round x : ℚ)| ≤ 1 / 2 := abs_sub_round x
    have hmr : (1 : ℚ) ≤ |(m : ℚ) - (round x : ℚ)| := by
      have hne' : m - round x ≠ 0 := sub_ne_zero.mpr hne
      have h1 : (1 : ℤ) ≤ |m - round x| := Int.one_le_abs hne'
      have h1' : ((1 : ℤ) : ℚ) ≤ ((|m - round x| : ℤ) : ℚ) := by exact_mod_cast h1
      rw [Int.cast_abs] at h1'
      push_cast at h1'
      linarith
    have htri : |(m : ℚ) - (round x : ℚ)| ≤ |(m : ℚ) - x| + |x - (round x : ℚ)| :=
      abs_sub_le _ _ _
    have habs : |x - (m : ℚ)| = |(m : ℚ) - x| := abs_sub_comm _ _
    have hn2 : (1 : ℚ) / n ≤ 1 / 2 := by
      apply div_le_div_of_nonneg_left (by norm_num) (by norm_num)
      exact_mod_cast hn
    rw [habs]
    linarith

/-- **The decidable family gate**: all 13 speeds clear the round-point check at `t`. -/
def ratWitnessOK (v : Fin 13 → ℤ) (t : ℚ) : Prop :=
  ∀ i, ratLonelyAt 14 (v i) t

instance (v : Fin 13 → ℤ) (t : ℚ) : Decidable (ratWitnessOK v t) := by
  unfold ratWitnessOK
  infer_instance

/-- The ℚ→ℝ loneliness lift: an all-integers rational bound at a rational time casts to
`Lonely` over ℝ.  Shared by the witness gate and the good-set pipeline. -/
theorem lonely_of_rat_forall {v : Fin 13 → ℤ} {t : ℚ}
    (h : ∀ i, ∀ m : ℤ, (1 : ℚ) / 14 ≤ |(v i : ℚ) * t - m|) :
    Lonely 14 v ((t : ℚ) : ℝ) := by
  intro i m
  have hq := h i m
  have hcast : ((|((v i : ℚ)) * t - (m : ℚ)| : ℚ) : ℝ) = |(v i : ℝ) * ((t : ℚ) : ℝ) - (m : ℝ)| := by
    push_cast
    ring_nf
  calc (1 : ℝ) / (14 : ℕ) = (((1 : ℚ) / 14 : ℚ) : ℝ) := by norm_num
    _ ≤ ((|((v i : ℚ)) * t - (m : ℚ)| : ℚ) : ℝ) := by exact_mod_cast hq
    _ = |(v i : ℝ) * ((t : ℚ) : ℝ) - (m : ℝ)| := hcast

/-- **Gate soundness (kernel-pure)**: a passing rational witness IS a lonely time. -/
theorem lonely_of_ratWitness {v : Fin 13 → ℤ} {t : ℚ} (h : ratWitnessOK v t) :
    Lonely 14 v ((t : ℚ) : ℝ) :=
  lonely_of_rat_forall fun i => ratLonelyAt_forall_int (by norm_num) (h i)

/-- Existence form: the census-table consumer. -/
theorem exists_lonely_of_ratWitness {v : Fin 13 → ℤ} (t : ℚ) (h : ratWitnessOK v t) :
    ∃ s : ℝ, Lonely 14 v s :=
  ⟨((t : ℚ) : ℝ), lonely_of_ratWitness h⟩

/-- Sanity row: the initial segment `{1,…,13}` at its classical witness `t = 1/14`. -/
example : ∃ s : ℝ, Lonely 14 (fun i : Fin 13 => (i : ℤ) + 1) s :=
  exists_lonely_of_ratWitness (1 / 14) (by native_decide)

end LRC14
end LonelyRunner
