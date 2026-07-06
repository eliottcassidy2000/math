/-
  TournamentH7.LRCLadderLoose — THE m/(12m+1) FAREY LADDER, PARAMETRIC LOOSE
  CERTIFICATE (kind-pasteur-2026-07-06-S21, HYP-4357).

  The broad spectral-gap census (lrc_spectral_gap_census_kps_S21,
  lrc_gap_boundary_sweep_kps_S21b) found the open gap `(1/13, 2/25)` EMPTY over
  all primitive 12-families of bounded height, and the STRUCTURAL reason: the
  family `{1,…,11, 12m}` has `M = m/(12m+1)` EXACTLY, a Farey/Stern-Brocot ladder
  converging to `1/12` whose first two rungs are the gap edges:
    · `m = 1`  → `1/13`  (the AP `{1,…,12}` — the gap's LOWER edge, TIGHT)
    · `m = 2`  → `2/25`  (the gap's UPPER edge)
    · `m ≥ 3`  → `3/37, 4/49, …` (all `> 2/25`).

  This file formalizes the LOOSE half: for `m ≥ 2`, `{1,…,11, 12m}` carries a
  `2/25`-margin point at `t = m/(12m+1)` — so it is loose (clears), NOT a gap
  member.  Via the `rational_point_margin` atom: at `k/s = m/(12m+1)`, every
  runner's residue lands in `[m, 11m+1] = [μ, s−μ]` (the small runners `i·m`
  for `i ≤ 11`, the big runner `12m·m ≡ 11m+1`), giving margin `≥ m/(12m+1) ≥
  2/25`.

  This COMPLEMENTS opus HYP-4356 (`{1,…,11, v}` loose for `12 ∤ v`): together the
  `{1,…,11, v}` slice is fully classified — loose for every `v` except `v = 12`
  (the unique tight AP at `1/13`).  The resonant `12 ∣ v` case, which opus's
  protection lemma does not reach, sits on this ladder and still avoids the open
  gap.

  Kernel-pure; no `sorry`, no `native_decide`.
-/
import Mathlib
import TournamentH7.LRCHarmonicGate

namespace LonelyRunner
namespace LadderLoose

open HarmonicGate
open scoped Classical

/-- The ladder family `{1,…,11, 12m}` as a map `Fin 12 → ℤ`: index `< 11` is the
speed `i+1 ∈ {1,…,11}`, index `11` is the big speed `12m`. -/
def ladderSpeed (m : ℤ) (i : Fin 12) : ℤ :=
  if (i : ℕ) < 11 then ((i : ℕ) + 1 : ℤ) else 12 * m

/-- **THE LADDER LOOSE CERTIFICATE**: for `m ≥ 2`, the family `{1,…,11, 12m}`
has a `2/25`-margin point at `t = m/(12m+1)` — it is loose (clears), hence NOT a
member of the open gap `(1/13, 2/25)`.  The resonant `12 ∣ v` case of the
`{1,…,11, v}` classification (complement to opus HYP-4356's `12 ∤ v`). -/
theorem ladder_family_loose (m : ℤ) (hm : 2 ≤ m) :
    ∃ tstar : ℝ, ∀ i : Fin 12, ∀ p : ℤ,
      (2 : ℝ) / 25 ≤ |(ladderSpeed m i : ℝ) * tstar - p| := by
  have hm0 : 0 < m := by omega
  set s : ℤ := 12 * m + 1 with hs_def
  have hs : 0 < s := by rw [hs_def]; omega
  -- the residue condition at k = m, μ = m
  have hcond : ∀ i : Fin 12,
      m ≤ (ladderSpeed m i * m) % s ∧ (ladderSpeed m i * m) % s ≤ s - m := by
    intro i
    rw [ladderSpeed]
    by_cases hi : (i : ℕ) < 11
    · -- small runner: speed = i+1 ≤ 11, residue = (i+1)*m (no wrap)
      rw [if_pos hi]
      set w : ℤ := ((i : ℕ) + 1 : ℤ) with hw
      have hw1 : 1 ≤ w := by rw [hw]; omega
      have hw11 : w ≤ 11 := by rw [hw]; omega
      have hlt : w * m < s := by rw [hs_def]; nlinarith [hw11, hm0]
      have hnonneg : 0 ≤ w * m := by positivity
      have hmod : (w * m) % s = w * m := Int.emod_eq_of_lt hnonneg hlt
      rw [hmod]
      constructor
      · nlinarith [hw1, hm0]
      · rw [hs_def]; nlinarith [hw11, hm0]
    · -- big runner: speed = 12m, residue = 11m+1
      rw [if_neg hi]
      have hbig : (12 * m * m) % s = 11 * m + 1 := by
        have hrw : 12 * m * m = (11 * m + 1) + (m - 1) * s := by rw [hs_def]; ring
        rw [hrw, Int.add_mul_emod_self_right]
        apply Int.emod_eq_of_lt
        · omega
        · rw [hs_def]; omega
      rw [hbig]
      constructor
      · omega
      · rw [hs_def]; omega
  -- apply the atom at t = m/s
  have hatom := rational_point_margin (ladderSpeed m) m s m hs hcond
  refine ⟨(m : ℝ) / s, ?_⟩
  intro i p
  have h := hatom i p
  -- margin ≥ m/s ≥ 2/25
  have hsR : (0 : ℝ) < (s : ℝ) := by exact_mod_cast hs
  have hmR : (0 : ℝ) < (m : ℝ) := by exact_mod_cast hm0
  have hlb : (2 : ℝ) / 25 ≤ (m : ℝ) / s := by
    rw [le_div_iff₀ hsR]
    have hsval : (s : ℝ) = 12 * m + 1 := by rw [hs_def]; push_cast; ring
    rw [hsval]
    have hmR2 : (2 : ℝ) ≤ (m : ℝ) := by exact_mod_cast hm
    nlinarith [hmR2]
  exact le_trans hlb h

#print axioms ladder_family_loose

end LadderLoose
end LonelyRunner
