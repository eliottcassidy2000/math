/-
  TournamentH7.LRCUEscape  (death-star-2026-07-12-S14)

  THE U-ESCAPE FLOOR (THM-721, Parts 2+3) — the strict-margin repair of the descent at
  near-dilates, in witness form. For a 13-speed family `V` with `vᵢ = bᵢ + L·kᵢ`
  (`|bᵢ| ≤ B`, `0 < L`), split PURE (`bᵢ = 0`) vs IMPURE (`bᵢ ≠ 0`, count `j`).

  On the closed line `t = (m + s)/L` the fast coordinate is EXACT: `distZ(vᵢ t) =
  distZ(kᵢ s)` for pure runners (no loss!), and `= distZ(kᵢ s + bᵢ t)` for impure ones.
  Fixing `s = s*` (a `1/13`-safe time of the pure lift set, supplied by the LRC(≤13)
  citation: ≤ 12 distinct pure lifts), each impure runner forbids a `u`-set which on the
  grid `{g/N}` occupies ≤ `(2/13)N + 2|bᵢ|` points (fibering `g ↦ g mod (N/|bᵢ|)` +
  at most two clipped windows per period). With `j ≤ 6`: total bad ≤ `(12/13)N + 12B < N`
  for `N = 157·B·B!` — a survivor `u* = g*/N` exists (finite pigeonhole, NO measure
  theory). Rounding `m = round(L·u* − s*)` gives `|t − u*| ≤ 1/(2L)`, so

        margin V t ≥ 1/13 − B/(2L)  >  1/14   once   91·B < L.

  This is the leg the 1D atom (LRCDecorrelation13) cannot reach: at a near-dilate the
  lift family is the tight AP (`reach(K) = 1/14` exactly) and `reach(K) − B/L` is below
  threshold; here the pure sub-family has ≤ 12 distinct speeds, so the LRC(≤13) floor
  `1/13` applies with NO loss on the fast coordinate. Sharp: `{L, 2L, …, 12L, 13L+1}`
  has exact `M = 1/13` at every diameter (verified, THM-721 Part 4).

  The 2D atom underneath is HYP-4342 (mac-mini-S10, LRCTorusRate.lean); this file is its
  new consumer: the u-escape union bound + the j ≤ 6 corollary, kernel-pure.

  Stage plan (staged like LRCFourierCompletion): STAGE A (this file, target green):
  distZ_neg, the window count, the per-period grid count `count_grid_near_int`.
  STAGE B: the fibered count for general |b| (b ∣ N), the pigeonhole survivor
  `exists_good_u`, and the witness assembly `lonely14_of_compressed_j6`.

  Toolchain note: `set`-introduced let-bindings poison `linarith` in this Mathlib
  (minimal repro: an ldecl + a cast-division atom in a hypothesis) — this file
  deliberately uses NO `set`, only literals and `obtain`-bound constants.
-/
import TournamentH7.LRCDecorrelation13

open TournamentH7.LRCWitness
open scoped Classical

namespace TournamentH7.LRCUEscape

/-- `distZ` is even. -/
lemma distZ_neg (y : ℝ) : distZ (-y) = distZ y := by
  apply le_antisymm
  · rw [le_distZ_iff (distZ (-y)) y]
    intro m
    have h := (le_distZ_iff (distZ (-y)) (-y)).1 le_rfl (-m)
    rwa [show -y - ((-m : ℤ) : ℝ) = -(y - (m : ℝ)) by push_cast; ring, abs_neg] at h
  · rw [le_distZ_iff (distZ y) (-y)]
    intro m
    have h := (le_distZ_iff (distZ y) y).1 le_rfl (-m)
    rwa [show y - ((-m : ℤ) : ℝ) = -(-y - (m : ℝ)) by push_cast; ring, abs_neg] at h

/-- Naturals `r < N'` satisfying `x ≤ r < x + ℓ` (a real half-open window): at most
`ℓ + 1` of them. -/
lemma card_range_filter_window (N' : ℕ) (x ℓ : ℝ) (hℓ : 0 ≤ ℓ)
    (P : ℕ → Prop) [DecidablePred P]
    (hP : ∀ r : ℕ, r < N' → P r → x ≤ (r : ℝ) ∧ (r : ℝ) < x + ℓ) :
    (((Finset.range N').filter P).card : ℝ) ≤ ℓ + 1 := by
  have hsub : ((Finset.range N').filter P).card ≤ (Finset.Icc (⌈x⌉) (⌈x⌉ + ⌊ℓ⌋)).card := by
    refine Finset.card_le_card_of_injOn (fun r : ℕ => (r : ℤ)) ?_ ?_
    · intro r hr
      have hrr := Finset.mem_filter.1 hr
      have hPr := hP r (Finset.mem_range.1 hrr.1) hrr.2
      have h1 : ⌈x⌉ ≤ (r : ℤ) := Int.ceil_le.2 hPr.1
      have h2 : (r : ℤ) ≤ ⌈x⌉ + ⌊ℓ⌋ := by
        have hcx : x ≤ ((⌈x⌉ : ℤ) : ℝ) := Int.le_ceil x
        have h3 : ((r : ℤ) - ⌈x⌉ : ℤ) ≤ ⌊ℓ⌋ := by
          apply Int.le_floor.2
          push_cast
          linarith [hPr.2]
        omega
      exact Finset.mem_Icc.2 ⟨h1, h2⟩
    · intro p _ q _ hpq
      exact Nat.cast_injective (by simpa using hpq)
  have h2 : (0 : ℤ) ≤ ⌊ℓ⌋ + 1 := by
    have := Int.floor_nonneg.2 hℓ
    omega
  have hcardZ : ((Finset.Icc (⌈x⌉) (⌈x⌉ + ⌊ℓ⌋)).card : ℤ) = ⌊ℓ⌋ + 1 := by
    rw [Int.card_Icc]
    rw [show (⌈x⌉ + ⌊ℓ⌋ + 1 - ⌈x⌉ : ℤ) = ⌊ℓ⌋ + 1 from by ring]
    exact Int.toNat_of_nonneg h2
  have hcardR : ((Finset.Icc (⌈x⌉) (⌈x⌉ + ⌊ℓ⌋)).card : ℝ) = ((⌊ℓ⌋ : ℝ) + 1) := by
    exact_mod_cast congrArg (fun z : ℤ => (z : ℝ)) hcardZ
  have hfl : (⌊ℓ⌋ : ℝ) ≤ ℓ := Int.floor_le ℓ
  calc (((Finset.range N').filter P).card : ℝ)
      ≤ ((Finset.Icc (⌈x⌉) (⌈x⌉ + ⌊ℓ⌋)).card : ℝ) := by exact_mod_cast hsub
    _ = (⌊ℓ⌋ : ℝ) + 1 := hcardR
    _ ≤ ℓ + 1 := by linarith

/-- **The per-period grid count.** Among the `N'` grid points `θ + r/N'`, `r < N'`, at most
`(2/13)·N' + 2` are within `1/13` of an integer. (The bad set is one circular arc of length
`2/13`; unrolled it is at most two windows whose lengths SUM to `(2/13)·N'`, each window
costing `+1`.) -/
lemma count_grid_near_int (θ : ℝ) (N' : ℕ) :
    (((Finset.range N').filter
        (fun r : ℕ => distZ (θ + (r : ℝ) / (N' : ℝ)) < 1/13)).card : ℝ)
      ≤ 2/13 * (N' : ℝ) + 2 := by
  rcases Nat.eq_zero_or_pos N' with hN0 | hN'
  · subst hN0
    simp only [Finset.range_zero, Finset.filter_empty, Finset.card_empty, Nat.cast_zero]
    norm_num
  have hNR : (0 : ℝ) < (N' : ℝ) := by exact_mod_cast hN'
  obtain ⟨m₁, hm₁⟩ : ∃ z : ℤ, z = ⌊θ - 1/13⌋ + 1 := ⟨_, rfl⟩
  have hm₁gt : θ - 1/13 < (m₁ : ℝ) := by
    have h := Int.lt_floor_add_one (θ - 1/13)
    rw [hm₁]; push_cast; linarith
  have hm₁le : (m₁ : ℝ) ≤ θ - 1/13 + 1 := by
    have h := Int.floor_le (θ - 1/13)
    rw [hm₁]; push_cast; linarith
  -- every bad point is (1/13)-close to m₁ or to m₁ + 1
  have hsplit : (Finset.range N').filter
      (fun r : ℕ => distZ (θ + (r : ℝ) / (N' : ℝ)) < 1/13) ⊆
      ((Finset.range N').filter
        (fun r : ℕ => |θ + (r : ℝ) / (N' : ℝ) - (m₁ : ℝ)| < 1/13)) ∪
      ((Finset.range N').filter
        (fun r : ℕ => |θ + (r : ℝ) / (N' : ℝ) - ((m₁ : ℝ) + 1)| < 1/13)) := by
    intro r hr
    have hr' := Finset.mem_filter.1 hr
    have hrange := hr'.1
    obtain ⟨m, hm⟩ : ∃ m : ℤ, |θ + (r : ℝ) / (N' : ℝ) - (m : ℝ)| < 1/13 := by
      by_contra hcon
      push_neg at hcon
      exact absurd ((le_distZ_iff (1/13) _).2 hcon) (not_le.2 hr'.2)
    have hrN : (r : ℝ) / (N' : ℝ) < 1 := by
      rw [div_lt_one hNR]
      exact_mod_cast Finset.mem_range.1 hrange
    have hr0 : (0 : ℝ) ≤ (r : ℝ) / (N' : ℝ) := by positivity
    obtain ⟨habs1, habs2⟩ := abs_lt.1 hm
    have hlow : m₁ ≤ m := by
      have h1 : θ - 1/13 < (m : ℝ) := by linarith
      have h2 : ⌊θ - 1/13⌋ < m := by
        by_contra hcon2
        push_neg at hcon2
        have h2a : (m : ℝ) ≤ θ - 1/13 :=
          le_trans (by exact_mod_cast hcon2) (Int.floor_le _)
        linarith
      omega
    have hhigh : m ≤ m₁ + 1 := by
      have h3 : (m : ℝ) < θ + 1 + 1/13 := by linarith
      have h4 : (m : ℝ) < (m₁ : ℝ) + 2 := by linarith
      have h6 : m < m₁ + 2 := by exact_mod_cast h4
      omega
    rcases (by omega : m = m₁ ∨ m = m₁ + 1) with hmm | hmm
    · subst hmm
      exact Finset.mem_union_left _ (Finset.mem_filter.2 ⟨hrange, hm⟩)
    · subst hmm
      refine Finset.mem_union_right _ (Finset.mem_filter.2 ⟨hrange, ?_⟩)
      have hcast : ((m₁ + 1 : ℤ) : ℝ) = (m₁ : ℝ) + 1 := by push_cast; ring
      rwa [hcast] at hm
  -- window bounds; split on whether the arc wraps (−1/13 ≤ a) or not, a := θ − m₁
  have htotal : ((((Finset.range N').filter
        (fun r : ℕ => |θ + (r : ℝ) / (N' : ℝ) - (m₁ : ℝ)| < 1/13)).card : ℝ) +
      (((Finset.range N').filter
        (fun r : ℕ => |θ + (r : ℝ) / (N' : ℝ) - ((m₁ : ℝ) + 1)| < 1/13)).card : ℝ))
      ≤ 2/13 * (N' : ℝ) + 2 := by
    rcases le_or_gt (-(1/13) : ℝ) (θ - (m₁ : ℝ)) with hwrap | hnowrap
    · -- wrapping: A₁ ⊆ [0, (1/13−a)N'), A₂ ⊆ [(1−1/13−a)N', N'); lengths sum to (2/13)N'
      have hc₁ : ((((Finset.range N').filter
          (fun r : ℕ => |θ + (r : ℝ) / (N' : ℝ) - (m₁ : ℝ)| < 1/13)).card : ℝ))
          ≤ (1/13 - (θ - (m₁ : ℝ))) * (N' : ℝ) + 1 := by
        apply card_range_filter_window N' 0 ((1/13 - (θ - (m₁ : ℝ))) * (N' : ℝ))
          (by nlinarith [hNR.le, hm₁gt])
        intro r _ hrP
        obtain ⟨habs1, habs2⟩ := abs_lt.1 hrP
        refine ⟨Nat.cast_nonneg r, ?_⟩
        have h5 : (r : ℝ) / (N' : ℝ) < 1/13 - (θ - (m₁ : ℝ)) := by linarith
        have h7 : (r : ℝ) < (1/13 - (θ - (m₁ : ℝ))) * (N' : ℝ) := by
          rw [div_lt_iff₀ hNR] at h5
          exact h5
        linarith
      have hc₂ : ((((Finset.range N').filter
          (fun r : ℕ => |θ + (r : ℝ) / (N' : ℝ) - ((m₁ : ℝ) + 1)| < 1/13)).card : ℝ))
          ≤ (1/13 + (θ - (m₁ : ℝ))) * (N' : ℝ) + 1 := by
        apply card_range_filter_window N'
          ((1 - 1/13 - (θ - (m₁ : ℝ))) * (N' : ℝ)) ((1/13 + (θ - (m₁ : ℝ))) * (N' : ℝ))
          (by nlinarith [hNR.le, hwrap])
        intro r hrN' hrP
        obtain ⟨habs1, habs2⟩ := abs_lt.1 hrP
        constructor
        · have h8 : 1 - 1/13 - (θ - (m₁ : ℝ)) < (r : ℝ) / (N' : ℝ) := by linarith
          have h9 : (1 - 1/13 - (θ - (m₁ : ℝ))) * (N' : ℝ) < (r : ℝ) := by
            rw [lt_div_iff₀ hNR] at h8
            exact h8
          linarith
        · have h10 : (1 - 1/13 - (θ - (m₁ : ℝ))) * (N' : ℝ)
              + (1/13 + (θ - (m₁ : ℝ))) * (N' : ℝ) = (N' : ℝ) := by ring
          rw [h10]
          exact_mod_cast hrN'
      linarith
    · -- non-wrapping: A₂ empty, A₁ one window of length (2/13)·N'
      have hA₂empty : ((Finset.range N').filter
          (fun r : ℕ => |θ + (r : ℝ) / (N' : ℝ) - ((m₁ : ℝ) + 1)| < 1/13)) = ∅ := by
        rw [Finset.filter_eq_empty_iff]
        intro r hrange
        rw [not_lt]
        have hrN : (r : ℝ) / (N' : ℝ) < 1 := by
          rw [div_lt_one hNR]
          exact_mod_cast Finset.mem_range.1 hrange
        have hval : θ + (r : ℝ) / (N' : ℝ) - ((m₁ : ℝ) + 1) < -(1/13) := by linarith
        calc (1/13 : ℝ) ≤ -(θ + (r : ℝ) / (N' : ℝ) - ((m₁ : ℝ) + 1)) := by linarith
          _ ≤ |θ + (r : ℝ) / (N' : ℝ) - ((m₁ : ℝ) + 1)| := neg_le_abs _
      have hc₁ : ((((Finset.range N').filter
          (fun r : ℕ => |θ + (r : ℝ) / (N' : ℝ) - (m₁ : ℝ)| < 1/13)).card : ℝ))
          ≤ 2/13 * (N' : ℝ) + 1 := by
        apply card_range_filter_window N'
          ((-(1/13) - (θ - (m₁ : ℝ))) * (N' : ℝ)) (2/13 * (N' : ℝ))
          (by nlinarith [hNR.le])
        intro r _ hrP
        obtain ⟨habs1, habs2⟩ := abs_lt.1 hrP
        constructor
        · have h13 : -(1/13) - (θ - (m₁ : ℝ)) < (r : ℝ) / (N' : ℝ) := by linarith
          have h14 : (-(1/13) - (θ - (m₁ : ℝ))) * (N' : ℝ) < (r : ℝ) := by
            rw [lt_div_iff₀ hNR] at h13
            exact h13
          linarith
        · have h16 : (r : ℝ) / (N' : ℝ) < 1/13 - (θ - (m₁ : ℝ)) := by linarith
          have h18 : (r : ℝ) < (1/13 - (θ - (m₁ : ℝ))) * (N' : ℝ) := by
            rw [div_lt_iff₀ hNR] at h16
            exact h16
          have h19 : (-(1/13) - (θ - (m₁ : ℝ))) * (N' : ℝ) + 2/13 * (N' : ℝ)
              = (1/13 - (θ - (m₁ : ℝ))) * (N' : ℝ) := by ring
          linarith
      rw [hA₂empty]
      simp only [Finset.card_empty, Nat.cast_zero]
      linarith
  have hfinal : (((Finset.range N').filter
      (fun r : ℕ => distZ (θ + (r : ℝ) / (N' : ℝ)) < 1/13)).card : ℝ)
      ≤ ((((Finset.range N').filter
        (fun r : ℕ => |θ + (r : ℝ) / (N' : ℝ) - (m₁ : ℝ)| < 1/13)).card : ℝ) +
      (((Finset.range N').filter
        (fun r : ℕ => |θ + (r : ℝ) / (N' : ℝ) - ((m₁ : ℝ) + 1)| < 1/13)).card : ℝ)) := by
    have h20 := Finset.card_le_card hsplit
    have h21 := Finset.card_union_le
      ((Finset.range N').filter
        (fun r : ℕ => |θ + (r : ℝ) / (N' : ℝ) - (m₁ : ℝ)| < 1/13))
      ((Finset.range N').filter
        (fun r : ℕ => |θ + (r : ℝ) / (N' : ℝ) - ((m₁ : ℝ) + 1)| < 1/13))
    exact_mod_cast le_trans h20 h21
  linarith

#print axioms distZ_neg
#print axioms card_range_filter_window
#print axioms count_grid_near_int

end TournamentH7.LRCUEscape
