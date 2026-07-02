/-
klein-2026-07-02-S112 (HYP-4017) — INTERVAL TRANSPORT: the epsilon-schedule primitives.

The band margin analysis (03-artifacts/drafts/interval-census-base-klein-S112.md) shows
the covering band is uniformly bounded away from tightness (optimum ≥ 1/12 everywhere,
ZERO tight rows at W ≤ 22 — the classical tight family (1,…,13) misses q = 14 and lives
on the sieve leg), so the length-invariant peel wrapper (kps-S17 design note) has an
unobstructed base. This file lands the transport primitives that upgrade point
certificates to interval certificates:

* `lonely_interval_of_margin` — loneliness at radius ρ at a point, plus a speed bound B,
  transports to 1/14-loneliness on the interval of half-length (ρ − 1/14)/B.
* `repeat13_lonely` — the citation node yields radius 1/13 (not just 1/14) for tuples
  with a repeated entry; margin 1/13 − 1/14 = 1/182 is the repeat leg's interval budget.
* `repeat_lonely_interval` — the two composed: bounded repeated-entry tuples have an
  explicit lonely interval of half-length 1/(182·B).
-/
import Mathlib
import TournamentH7.LRC13Citation

namespace LonelyRunner

/-- **Interval transport**: point-loneliness at radius ρ transports to 1/14-loneliness
on an interval, at the rate of one margin unit per unit speed. -/
theorem lonely_interval_of_margin {v : Fin 13 → ℤ} {t₀ ρ B δ : ℝ}
    (hpt : ∀ i : Fin 13, ∀ m : ℤ, ρ ≤ |(v i : ℝ) * t₀ - m|)
    (hB : ∀ i, |(v i : ℝ)| ≤ B)
    (hδ : 0 ≤ δ) (hmargin : 1 / 14 + B * δ ≤ ρ) :
    ∀ t : ℝ, |t - t₀| ≤ δ → Lonely 14 v t := by
  intro t ht i m
  have key : |(v i : ℝ) * t₀ - m| ≤ |(v i : ℝ) * t - m| + |(v i : ℝ)| * |t₀ - t| := by
    calc |(v i : ℝ) * t₀ - m|
        = |((v i : ℝ) * t - m) + (v i : ℝ) * (t₀ - t)| := by congr 1; ring
      _ ≤ |(v i : ℝ) * t - m| + |(v i : ℝ) * (t₀ - t)| := abs_add_le _ _
      _ = |(v i : ℝ) * t - m| + |(v i : ℝ)| * |t₀ - t| := by rw [abs_mul]
  have h1 : |(v i : ℝ)| * |t₀ - t| ≤ B * δ := by
    apply mul_le_mul (hB i) ?_ (abs_nonneg _) ((abs_nonneg _).trans (hB i))
    rwa [abs_sub_comm]
  have h2 := hpt i m
  push_cast
  linarith

/-- **The repeat leg at radius 1/13**: a 13-tuple with a repeated entry is lonely at the
citation's full gap 1/(k+1) with k ≤ 12 distinct speeds — not merely at the band radius
1/14. The margin 1/182 funds the interval form below. (Body adapted from
`lonely14_of_repeat`, keeping the un-relaxed bound.) -/
theorem repeat13_lonely (cite : LRCUpTo13) (v : Fin 13 → ℤ)
    (hv : ∀ i, v i ≠ 0) {i j : Fin 13} (hij : i ≠ j) (heq : v i = v j) :
    ∃ t : ℝ, ∀ i0 : Fin 13, ∀ m : ℤ, (1 : ℝ) / 13 ≤ |(v i0 : ℝ) * t - m| := by
  set S : Finset ℤ := Finset.univ.image v with hS
  have hcard : S.card ≤ 12 := by
    have hle : S.card ≤ 13 := by
      simpa [hS] using Finset.card_image_le (s := (Finset.univ : Finset (Fin 13))) (f := v)
    rcases Nat.lt_or_ge S.card 13 with h | h
    · omega
    · exfalso
      have heq13 : S.card = 13 := le_antisymm hle h
      have hinj : Set.InjOn v (Finset.univ : Finset (Fin 13)) := by
        apply Finset.injOn_of_card_image_eq
        simpa [hS] using heq13
      exact hij (hinj (Finset.mem_coe.mpr (Finset.mem_univ i))
        (Finset.mem_coe.mpr (Finset.mem_univ j)) heq)
  let e : Fin S.card ≃ S := S.equivFin.symm
  set w : Fin S.card → ℤ := fun m => (e m : ℤ) with hw
  have hwne : ∀ m, w m ≠ 0 := by
    intro m
    have hmem : (e m : ℤ) ∈ S := (e m).2
    simp only [hS, Finset.mem_image] at hmem
    obtain ⟨i0, -, hi0⟩ := hmem
    rw [hw]
    simp only
    rw [← hi0]
    exact hv i0
  obtain ⟨t, hL⟩ := cite S.card hcard w hwne
  refine ⟨t, fun i0 m => ?_⟩
  have hmem : v i0 ∈ S := by
    simp [hS, Finset.mem_image]
  set m0 : Fin S.card := e.symm ⟨v i0, hmem⟩ with hm0
  have hwm0 : w m0 = v i0 := by
    rw [hw, hm0]
    simp only
    rw [Equiv.apply_symm_apply]
  have hstep : (1 : ℝ) / ((S.card : ℝ) + 1) ≤ |(v i0 : ℝ) * t - m| := by
    have := hL m0 m
    rw [hwm0] at this
    have hcast : ((S.card + 1 : ℕ) : ℝ) = (S.card : ℝ) + 1 := by push_cast; ring
    rwa [hcast] at this
  have h13 : (1 : ℝ) / 13 ≤ 1 / ((S.card : ℝ) + 1) := by
    have hpos : (0 : ℝ) < (S.card : ℝ) + 1 := by positivity
    have hle13 : ((S.card : ℝ) + 1) ≤ 13 := by
      have h12 : (S.card : ℝ) ≤ 12 := by exact_mod_cast hcard
      linarith
    exact one_div_le_one_div_of_le hpos hle13
  linarith

/-- **The repeat leg's interval form**: a bounded 13-tuple with a repeated entry has an
explicit lonely INTERVAL of half-length 1/(182·B) — the length-positivity the peel
induction consumes, delivered by the citation node's 1/13-vs-1/14 margin. -/
theorem repeat_lonely_interval (cite : LRCUpTo13) (v : Fin 13 → ℤ)
    (hv : ∀ i, v i ≠ 0) {i j : Fin 13} (hij : i ≠ j) (heq : v i = v j)
    {B : ℝ} (hB : ∀ i, |(v i : ℝ)| ≤ B) (hBpos : 0 < B) :
    ∃ t₀ : ℝ, ∀ t : ℝ, |t - t₀| ≤ 1 / (182 * B) → Lonely 14 v t := by
  obtain ⟨t₀, h13⟩ := repeat13_lonely cite v hv hij heq
  refine ⟨t₀, lonely_interval_of_margin h13 hB (by positivity) ?_⟩
  have hBδ : B * (1 / (182 * B)) = 1 / 182 := by
    field_simp
  rw [hBδ]
  norm_num

end LonelyRunner
