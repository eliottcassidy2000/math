/-
  TournamentH7.LRCGoodPipeline — THE GOOD-SET PIPELINE (kind-pasteur-2026-07-02-S13,
  HYP-3970).  Single-writer satellite.

  The last named wiring piece of the census leg: the DANGER-COMB COMPLETENESS BRIDGE.
  A point of `[0,1)` avoiding the wrapped danger comb of speed `s` is `h`-safe for `s`
  against EVERY integer — the half-open comb teeth cover the whole danger set, with the
  wrap layer supplying the tooth that straddles `1`.

  Composed with the S12 engines this closes the pipeline:

      0 < length (goodRegion speeds h)          [ONE computable rational per class]
        → exists_mem_of_length_pos              [explicit rational point]
        → mem_diff                              [avoids every danger tooth]
        → not_mem_wrap_comb_forall              [safe for every speed, every integer]
        → lonely_of_rat_forall                  [Lonely over ℝ]

  `exists_lonely_of_goodRegion_pos` is the composite gate: a 13-family whose good region
  has positive computed length is lonely — no witness search needed, ONE length
  computation per family/class.  This is the second census-leg gate (dual to the
  `ratWitnessOK` witness-table gate).
-/
import TournamentH7.LRCRegionDiff
import TournamentH7.LRCRatWitness

namespace LonelyRunner
namespace RatIntervals

/-! ## The danger-comb completeness bridge -/

/-- Tooth widths of the phase-0 comb are `2h/s ≤ 1`. -/
private theorem comb_widths {s : ℕ} (hs : 0 < s) {h : ℚ} (hh2 : 2 * h ≤ 1) :
    ∀ p ∈ comb s h 0, p.2 - p.1 ≤ 1 := by
  intro p hp
  unfold comb at hp
  simp only [List.map_map, List.mem_map, List.mem_range, Function.comp_apply] at hp
  obtain ⟨k, _, rfl⟩ := hp
  simp only
  have hsQ : (0 : ℚ) < (s : ℚ) := by exact_mod_cast hs
  have hdiff : ((k : ℚ) + 0 + h) / s - ((k : ℚ) + 0 - h) / s = 2 * h / s := by
    field_simp
    try ring
  rw [hdiff, div_le_one hsQ]
  have hs1 : (1 : ℚ) ≤ (s : ℚ) := by exact_mod_cast hs
  linarith

/-- **The completeness bridge**: a point of `[0,1)` outside the wrapped phase-0 comb of
speed `s` is `h`-far from every integer multiple position — the comb teeth (with the
wrap tooth at the seam) cover the entire danger set. -/
theorem not_mem_wrap_comb_forall {s : ℕ} (hs : 0 < s) {h : ℚ} (hh2 : 2 * h ≤ 1)
    {x : ℚ} (hx0 : 0 ≤ x) (hx1 : x < 1)
    (hnot : ¬ mem x (wrap (comb s h 0))) :
    ∀ m : ℤ, h ≤ |(s : ℚ) * x - m| := by
  intro m
  by_contra hcon
  apply hnot
  have habs : |(s : ℚ) * x - m| < h := lt_of_not_ge hcon
  rw [abs_lt] at habs
  have hsZ : (0 : ℤ) < (s : ℤ) := by exact_mod_cast hs
  have hsZne : ((s : ℕ) : ℤ) ≠ 0 := by exact_mod_cast hs.ne'
  rw [mem_wrap hx0 hx1 (comb_widths hs hh2)]
  refine ⟨-(m / (s : ℤ)), ?_⟩
  rw [TournamentH7.CombPatterns.mem_comb hs]
  have hem0 : 0 ≤ m % (s : ℤ) := Int.emod_nonneg m hsZne
  have hemS : m % (s : ℤ) < (s : ℤ) := Int.emod_lt_of_pos m hsZ
  have hme : m % (s : ℤ) + (s : ℤ) * (m / (s : ℤ)) = m := Int.emod_add_ediv m (s : ℤ)
  refine ⟨(m % (s : ℤ)).toNat, ?_, ?_, ?_⟩
  · omega
  · -- (k : ℚ) + 0 − h ≤ s·(x + n)
    have htn : (((m % (s : ℤ)).toNat : ℕ) : ℚ) = ((m % (s : ℤ) : ℤ) : ℚ) := by
      exact_mod_cast Int.toNat_of_nonneg hem0
    have hmeQ : ((m % (s : ℤ) : ℤ) : ℚ) + (s : ℚ) * ((m / (s : ℤ) : ℤ) : ℚ) = (m : ℚ) := by
      exact_mod_cast hme
    rw [htn]
    push_cast
    nlinarith [habs.1, hmeQ]
  · have htn : (((m % (s : ℤ)).toNat : ℕ) : ℚ) = ((m % (s : ℤ) : ℤ) : ℚ) := by
      exact_mod_cast Int.toNat_of_nonneg hem0
    have hmeQ : ((m % (s : ℤ) : ℤ) : ℚ) + (s : ℚ) * ((m / (s : ℤ) : ℤ) : ℚ) = (m : ℚ) := by
      exact_mod_cast hme
    rw [htn]
    push_cast
    nlinarith [habs.2, hmeQ]

/-! ## The good region and the composite gate -/

/-- The good region of a speed list at band `h`: the unit window minus every wrapped
danger comb.  Built with the COMPUTATIONAL difference (`diffF`, live pieces only) so
per-class positivity is one tractable `native_decide`. -/
def goodRegion (speeds : List ℤ) (h : ℚ) : Region :=
  diffF [((0 : ℚ), 1)] (speeds.flatMap fun s => wrap (comb s.toNat h 0))

/-- A member of the good region is inside `[0,1)` and `h`-safe for every speed against
every integer. -/
theorem good_mem_safe {speeds : List ℤ} {h : ℚ} (hh2 : 2 * h ≤ 1)
    (hpos : ∀ s ∈ speeds, 0 < s) {x : ℚ} (hx : mem x (goodRegion speeds h)) :
    (0 ≤ x ∧ x < 1) ∧ ∀ s ∈ speeds, ∀ m : ℤ, h ≤ |(s : ℚ) * x - m| := by
  obtain ⟨hwin, havoid⟩ := mem_diffF.mp hx
  obtain ⟨p, hp, hx1, hx2⟩ := hwin
  rcases List.mem_singleton.mp hp with rfl
  simp only at hx1 hx2
  refine ⟨⟨hx1, hx2⟩, ?_⟩
  intro s hs m
  have hspos := hpos s hs
  have hsnat : 0 < s.toNat := by omega
  have hbridge := not_mem_wrap_comb_forall hsnat hh2 hx1 hx2 (fun hmem => by
    obtain ⟨q, hq, hxq⟩ := hmem
    exact havoid q (List.mem_flatMap.mpr ⟨s, hs, hq⟩) hxq) m
  have hcast : ((s.toNat : ℕ) : ℚ) = ((s : ℤ) : ℚ) := by
    exact_mod_cast Int.toNat_of_nonneg hspos.le
  rwa [hcast] at hbridge

/-- **THE COMPOSITE GATE**: a 13-family of positive speeds whose good region has positive
length is lonely.  One computable rational length per family/class — no witness search. -/
theorem exists_lonely_of_goodRegion_pos (v : Fin 13 → ℤ) (hv : ∀ i, 0 < v i)
    (hlen : 0 < length (goodRegion (List.ofFn v) (1 / 14))) :
    ∃ t : ℝ, Lonely 14 v t := by
  obtain ⟨x, hx⟩ := exists_mem_of_length_pos hlen
  have hpos : ∀ s ∈ List.ofFn v, 0 < s := by
    intro s hs
    obtain ⟨i, rfl⟩ := List.mem_ofFn.mp hs
    exact hv i
  have hsafe := good_mem_safe (by norm_num) hpos hx
  exact ⟨((x : ℚ) : ℝ), LRC14.lonely_of_rat_forall fun i m =>
    hsafe.2 (v i) (List.mem_ofFn.mpr ⟨i, rfl⟩) m⟩

/-- End-to-end pipeline test: the consecutive block `{3,…,15}` is lonely because its
good region has positive computed length — ONE `native_decide`, no witness search. -/
example : ∃ t : ℝ, Lonely 14 (LRC14.blockFamily 15) t := by
  apply exists_lonely_of_goodRegion_pos
  · intro i
    unfold LRC14.blockFamily
    have := i.isLt
    omega
  · native_decide

end RatIntervals
end LonelyRunner
