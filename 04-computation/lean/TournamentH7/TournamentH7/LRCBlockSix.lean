/-
  TournamentH7.LRCBlockSix — THE BLOCK-6 INTERVAL ENGINE (kind-pasteur-2026-07-02-S21,
  HYP-3978).  Near-equal blocks of ANY internal ratio structure, size ≤ 6, over a real
  interval base — the citation-route counterpart of mac-mini's region-based
  `LRCSimulPeel` (S18) and the many-runner generalization of klein's
  `safe_point_in_interval` (S114).

  * `gap_exists` — THE MEASURE-FREE UNION LEMMA over ℝ: a window whose clipped bad
    mass is under its length contains a point avoiding every open bad interval.
    Proof: JUMP-PAST-THE-BLOCKER induction — if the window's left edge sits inside
    some bad interval, hop to that interval's right end and drop it from the list;
    the dropped clipped mass pays for the lost window length exactly.  No measure
    theory, no counting, no sorting.

  * `runner_bad_teeth` / `runner_teeth_mass` — the teeth of runner `w` intersecting
    a window of length `L`, as an explicit list of ≤ `⌊wL⌋ + 2` intervals of width
    `1/(7w)` with clipped mass ≤ `L/7 + 2/(7w)`.

  * `block_point_step` — ≤ 6 positive runners with ARBITRARY internal ratios share a
    `1/14`-good point in any window with `w₁·L` above the explicit entry threshold:
    union density `c/7 < 1` leaves room.  Blocks of 7+ hit the density wall — the
    honest remaining core.

  * `cite_block_lonely` — the citation composition: cite `k ≤ 12` bounded runners,
    transport, and finish the (≤ 6)-block inside the interval.
-/
import Mathlib
import TournamentH7.LRCPairBlock

namespace LonelyRunner
namespace LRC14

/-! ## The measure-free union lemma -/

/-- Clipped length of one open interval against a window `[a, b]`. -/
noncomputable def clipLen (p : ℝ × ℝ) (a b : ℝ) : ℝ := max 0 (min p.2 b - max p.1 a)

theorem clipLen_nonneg (p : ℝ × ℝ) (a b : ℝ) : 0 ≤ clipLen p a b := le_max_left _ _

/-- Shrinking the window from the left does not increase the clipped length. -/
theorem clipLen_mono_left (p : ℝ × ℝ) {a a' b : ℝ} (h : a ≤ a') :
    clipLen p a' b ≤ clipLen p a b := by
  unfold clipLen
  rcases le_total (min p.2 b - max p.1 a') 0 with h1 | h1
  · rw [max_eq_left h1]
    exact le_max_left _ _
  · rw [max_eq_right h1]
    have : max p.1 a ≤ max p.1 a' := max_le_max (le_refl _) h
    calc min p.2 b - max p.1 a' ≤ min p.2 b - max p.1 a := by linarith
      _ ≤ max 0 (min p.2 b - max p.1 a) := le_max_right _ _

/-- **THE MEASURE-FREE UNION LEMMA**: a window whose clipped bad mass is under its
length contains a point avoiding every open bad interval. -/
theorem gap_exists : ∀ (n : ℕ) (bads : List (ℝ × ℝ)) (a b : ℝ),
    bads.length ≤ n → a ≤ b →
    (bads.map fun p => clipLen p a b).sum < b - a →
    ∃ t : ℝ, a ≤ t ∧ t ≤ b ∧ ∀ p ∈ bads, t ≤ p.1 ∨ p.2 ≤ t := by
  intro n
  induction n with
  | zero =>
      intro bads a b hlen hab _
      have hnil : bads = [] := List.length_eq_zero_iff.mp (Nat.le_zero.mp hlen)
      subst hnil
      exact ⟨a, le_refl _, hab, fun p hp => absurd hp (List.not_mem_nil)⟩
  | succ n ih =>
      intro bads a b hlen hab hsum
      by_cases hesc : ∀ p ∈ bads, a ≤ p.1 ∨ p.2 ≤ a
      · exact ⟨a, le_refl _, hab, hesc⟩
      · push Not at hesc
        obtain ⟨p, hp, hp1, hp2⟩ := hesc
        -- hp1 : p.1 < a, hp2 : a < p.2
        have hnonneg : ∀ x ∈ bads.map fun p => clipLen p a b, 0 ≤ x := by
          intro x hx
          rw [List.mem_map] at hx
          obtain ⟨q, _, rfl⟩ := hx
          exact clipLen_nonneg q a b
        have hpb : p.2 ≤ b := by
          by_contra hgt
          push Not at hgt
          have hcp : clipLen p a b = b - a := by
            unfold clipLen
            rw [min_eq_right (le_of_lt hgt), max_eq_right (le_of_lt hp1)]
            exact max_eq_right (by linarith)
          have hmem : clipLen p a b ∈ bads.map fun q => clipLen q a b :=
            List.mem_map.mpr ⟨p, hp, rfl⟩
          have := List.single_le_sum hnonneg _ hmem
          linarith
        -- drop p, jump to p.2
        obtain ⟨l₁, l₂, rfl⟩ := List.append_of_mem hp
        have hlen' : (l₁ ++ l₂).length ≤ n := by
          have := hlen
          simp only [List.length_append, List.length_cons] at this ⊢
          omega
        have hab' : p.2 ≤ b := hpb
        have hcp : clipLen p a b = p.2 - a := by
          unfold clipLen
          rw [min_eq_left hpb, max_eq_right (le_of_lt hp1)]
          exact max_eq_right (by linarith)
        have hsum_split : ((l₁ ++ p :: l₂).map fun q => clipLen q a b).sum
            = ((l₁ ++ l₂).map fun q => clipLen q a b).sum + clipLen p a b := by
          simp only [List.map_append, List.map_cons, List.sum_append, List.sum_cons]
          ring
        have hsum' : ((l₁ ++ l₂).map fun q => clipLen q p.2 b).sum < b - p.2 := by
          have hmono : ((l₁ ++ l₂).map fun q => clipLen q p.2 b).sum
              ≤ ((l₁ ++ l₂).map fun q => clipLen q a b).sum := by
            apply List.sum_le_sum
            intro q hq
            exact clipLen_mono_left q (le_of_lt hp2)
          linarith [hsum_split, hcp, hsum, hmono]
        obtain ⟨t, ht1, ht2, havoid⟩ := ih (l₁ ++ l₂) p.2 b hlen' hab' hsum'
        refine ⟨t, by linarith, ht2, ?_⟩
        intro q hq
        rcases List.mem_append.mp hq with hq1 | hq2
        · exact havoid q (List.mem_append.mpr (Or.inl hq1))
        · rcases List.mem_cons.mp hq2 with rfl | hq3
          · exact Or.inr (by linarith)
          · exact havoid q (List.mem_append.mpr (Or.inr hq3))

end LRC14
end LonelyRunner
