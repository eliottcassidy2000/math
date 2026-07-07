/-
  TournamentH7.LRCCoarseReduction — THE COARSE / SCALE REDUCTION as a Route-1 tool.
  (kind-pasteur-2026-07-06-S53, HYP-4707; formalizes the kps-S52 bound `M(v) ≥ M(K) − A/L`.)

  Post opus-S130: Route 2 is retired (J-K bounds accumulation points, not the supremum);
  the correctly-aimed route is **Route 1** — bound the loneliness threshold `1/14` DIRECTLY.
  The coarse reduction is a general fact about the loneliness margin, so it lives on Route 1
  and **grounds the MULTI-SCALE families in the settled LRC(≤13)**.

  If `v i = a i + L·k i` with `|a i| ≤ A` and the coarse family `k` is lonely at `s₀` with
  margin `μ`, then `v` is lonely at `s₀/L` with margin `≥ μ − A/L`: the fine offsets `a i`,
  seen at the fast time `s₀/L`, perturb each phase by only `|a i|·|s₀|/L ≤ A/L`, while the
  coarse phase `k i·s₀` carries the full margin `μ`.  With `μ = 1/13` (LRC(≤13) on a
  `≤ 12`-distinct coarse family) and `A/L ≤ 1/13 − 1/14 = 1/182`, the transferred margin is
  `≥ 1/14` — LONELY.

  So every 13-speed family that clusters into `≤ 12` groups at some scale `L` (equivalently,
  is a bounded perturbation of a dilation of a `≤ 12`-speed family) is lonely by the settled
  cases — no new analysis.  The density node (Route 1) is thereby reduced to the SINGLE-SCALE
  (compressed) families.  Kernel-pure, IVT-free: one 1-Lipschitz triangle inequality.
-/
import Mathlib
import TournamentH7.LonelyRunner
import TournamentH7.LRC13Citation

namespace LonelyRunner
namespace CoarseReduction

/-- `Lonely` is invariant under replacing the time `t` by its fractional part
(`‖v i · t‖` only depends on `t mod 1`).  Used to normalize the LRC(≤13) witness into
`[0,1)`, so `|s₀| ≤ 1` for the transfer. -/
theorem lonely_fract {ι : Type*} (n : ℕ) (v : ι → ℤ) (t : ℝ) (h : Lonely n v t) :
    Lonely n v (Int.fract t) := by
  intro i m
  have hkey : (v i : ℝ) * Int.fract t - m
      = (v i : ℝ) * t - ((m + v i * ⌊t⌋ : ℤ) : ℝ) := by
    have hf : Int.fract t = t - (⌊t⌋ : ℝ) := rfl
    rw [hf]; push_cast; ring
  rw [hkey]; exact h i _

/-- **The coarse margin transfer (core, elementary).**  If every coarse phase `k i · s₀`
is `≥ μ` from every integer, `|a i| ≤ A`, `|s₀| ≤ 1`, and `0 < L`, then every fine phase
`(a i + L·k i)·(s₀/L)` is `≥ μ − A/L` from every integer.  This is `M(v) ≥ M(K) − A/L`
at the level of a fixed witness `s₀`. -/
theorem reach_transfer_coarse {ι : Type*} (k a : ι → ℤ) (L : ℤ) (A μ s₀ : ℝ)
    (hL : 0 < (L : ℝ)) (hs0 : |s₀| ≤ 1) (ha : ∀ i, |(a i : ℝ)| ≤ A)
    (hk : ∀ i, ∀ m : ℤ, μ ≤ |(k i : ℝ) * s₀ - m|) :
    ∀ i, ∀ m : ℤ, μ - A / L ≤ |((a i + L * k i : ℤ) : ℝ) * (s₀ / L) - m| := by
  intro i m
  set X : ℝ := (k i : ℝ) * s₀ - m with hX
  set Y : ℝ := (a i : ℝ) * s₀ / L with hY
  have hLne : (L : ℝ) ≠ 0 := ne_of_gt hL
  have hsplit : ((a i + L * k i : ℤ) : ℝ) * (s₀ / L) - m = X + Y := by
    rw [hX, hY]; push_cast; field_simp; ring
  have hXbound : μ ≤ |X| := hk i m
  have hYbound : |Y| ≤ A / L := by
    have hAnn : 0 ≤ A := le_trans (abs_nonneg _) (ha i)
    have hYabs : |Y| = |(a i : ℝ)| * |s₀| / L := by
      rw [hY, abs_div, abs_mul, abs_of_pos hL]
    rw [hYabs]
    have hnum : |(a i : ℝ)| * |s₀| ≤ A := by
      calc |(a i : ℝ)| * |s₀| ≤ A * 1 :=
            mul_le_mul (ha i) hs0 (abs_nonneg _) hAnn
        _ = A := mul_one A
    gcongr
  have hrev : |X| - |Y| ≤ |X + Y| := by
    have h := abs_sub_abs_le_abs_sub X (-Y)
    rw [abs_neg, sub_neg_eq_add] at h
    exact h
  rw [hsplit]; linarith

/-- **Coarse reduction ⟹ Lonely 14.**  If the coarse family `k` is `1/13`-lonely at `s₀`
(the LRC(≤13) margin), the offsets satisfy `|a i| ≤ A`, and the scale clears the budget
`A/L ≤ 1/13 − 1/14`, then `v i = a i + L·k i` is `1/14`-lonely at `s₀/L`. -/
theorem lonely14_of_coarse_of_lonely13 {ι : Type*} (k a : ι → ℤ) (L : ℤ) (A s₀ : ℝ)
    (hL : 0 < (L : ℝ)) (hs0 : |s₀| ≤ 1) (ha : ∀ i, |(a i : ℝ)| ≤ A)
    (hbudget : A / (L : ℝ) ≤ 1 / 13 - 1 / 14) (hk : Lonely 13 k s₀) :
    Lonely 14 (fun i => a i + L * k i) (s₀ / L) := by
  intro i m
  have hstep := reach_transfer_coarse k a L A ((1 : ℝ) / ((13 : ℕ) : ℝ)) s₀ hL hs0 ha hk i m
  have e13 : (1 : ℝ) / ((13 : ℕ) : ℝ) = 1 / 13 := by norm_num
  have e14 : (1 : ℝ) / ((14 : ℕ) : ℝ) = 1 / 14 := by norm_num
  rw [e13] at hstep
  rw [e14]
  have h14 : (1 : ℝ) / 14 ≤ 1 / 13 - A / (L : ℝ) := by linarith
  simpa using le_trans h14 hstep

/-- A 13-tuple whose distinct nonzero values number `≤ 12` is `1/13`-lonely at some time,
by the LRC(≤13) citation (the coarse-family input to the transfer). -/
theorem lonely13_of_card_le12 (cite : LRCUpTo13) (k : Fin 13 → ℤ)
    (hne : ∀ i, k i ≠ 0) (hcard : (Finset.univ.image k).card ≤ 12) :
    ∃ t : ℝ, Lonely 13 k t := by
  set S : Finset ℤ := Finset.univ.image k with hS
  let e : Fin S.card ≃ S := S.equivFin.symm
  set w : Fin S.card → ℤ := fun mm => (e mm : ℤ) with hw
  have hwne : ∀ mm, w mm ≠ 0 := by
    intro mm
    have hmem : (e mm : ℤ) ∈ S := (e mm).2
    simp only [hS, Finset.mem_image] at hmem
    obtain ⟨i0, -, hi0⟩ := hmem
    rw [hw]; simp only; rw [← hi0]; exact hne i0
  obtain ⟨t, hL⟩ := cite S.card hcard w hwne
  refine ⟨t, fun i0 m => ?_⟩
  have hmem : k i0 ∈ S := by simp [hS, Finset.mem_image]
  set m0 : Fin S.card := e.symm ⟨k i0, hmem⟩ with hm0
  have hwm0 : w m0 = k i0 := by
    rw [hw, hm0]; simp only; rw [Equiv.apply_symm_apply]
  have hstep : (1 : ℝ) / ((S.card : ℝ) + 1) ≤ |(k i0 : ℝ) * t - m| := by
    have hcite := hL m0 m
    rw [hwm0] at hcite
    have hcast : ((S.card + 1 : ℕ) : ℝ) = (S.card : ℝ) + 1 := by push_cast; ring
    rwa [hcast] at hcite
  have h13 : (1 : ℝ) / 13 ≤ 1 / ((S.card : ℝ) + 1) := by
    have hpos : (0 : ℝ) < (S.card : ℝ) + 1 := by positivity
    have hle : ((S.card : ℝ) + 1) ≤ 13 := by
      have : (S.card : ℝ) ≤ 12 := by exact_mod_cast hcard
      linarith
    exact one_div_le_one_div_of_le hpos hle
  have ecast : (1 : ℝ) / ((13 : ℕ) : ℝ) = 1 / 13 := by norm_num
  rw [ecast]; linarith

/-- **THE HEADLINE (grounds the multi-scale branch of Route 1 in LRC(≤13)).**  A 13-speed
family `v` with a scale decomposition `v i = a i + L·k i` whose coarse part `k` has `≤ 12`
distinct nonzero values (`v` clusters into `≤ 12` groups at scale `L`), with bounded offsets
`|a i| ≤ A` and the scale clearing `A/L ≤ 1/13 − 1/14`, is LONELY at `1/14` — by the settled
LRC(≤13).  No new analysis; the density node is left only the single-scale families. -/
theorem lonely14_of_coarse_le12 (cite : LRCUpTo13) (v k a : Fin 13 → ℤ) (L : ℤ) (A : ℝ)
    (hdecomp : ∀ i, v i = a i + L * k i)
    (hL : 0 < (L : ℝ)) (ha : ∀ i, |(a i : ℝ)| ≤ A)
    (hbudget : A / (L : ℝ) ≤ 1 / 13 - 1 / 14)
    (hkne : ∀ i, k i ≠ 0) (hcard : (Finset.univ.image k).card ≤ 12) :
    ∃ t : ℝ, Lonely 14 v t := by
  obtain ⟨s₀, hs₀L⟩ := lonely13_of_card_le12 cite k hkne hcard
  have hfrac : Lonely 13 k (Int.fract s₀) := lonely_fract 13 k s₀ hs₀L
  have hs0abs : |Int.fract s₀| ≤ 1 := by
    rw [abs_of_nonneg (Int.fract_nonneg s₀)]; linarith [Int.fract_lt_one s₀]
  have hlon : Lonely 14 (fun i => a i + L * k i) (Int.fract s₀ / L) :=
    lonely14_of_coarse_of_lonely13 k a L A (Int.fract s₀) hL hs0abs ha hbudget hfrac
  refine ⟨Int.fract s₀ / L, fun i m => ?_⟩
  have hthis := hlon i m
  rw [hdecomp i]
  simpa using hthis

/-! ## Axiom audit -/

#print axioms lonely_fract
#print axioms reach_transfer_coarse
#print axioms lonely14_of_coarse_of_lonely13
#print axioms lonely13_of_card_le12
#print axioms lonely14_of_coarse_le12

end CoarseReduction
end LonelyRunner
