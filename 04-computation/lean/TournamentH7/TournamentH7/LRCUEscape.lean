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

open TournamentH7.LRCWitness TournamentH7.LRCDecorr
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

/-! ### Stage B — the fibered count, the pigeonhole survivor, and the witness assembly -/

/-- **The fibered grid count.** For `β ≠ 0` and `N = |β|·N'`, among the `N` grid points
`θ + β·g/N`, `g < N`, at most `(2/13)·N + 2|β|` are within `1/13` of an integer
(each residue class `g ≡ r (mod N')` reduces to the unit-speed count at phase `±θ`). -/
lemma count_grid_scaled (θ : ℝ) (β : ℤ) (hβ : β ≠ 0) (N N' : ℕ)
    (hN : N = β.natAbs * N') (hN' : 0 < N') :
    (((Finset.range N).filter
        (fun g : ℕ => distZ (θ + (β : ℝ) * (g : ℝ) / (N : ℝ)) < 1/13)).card : ℝ)
      ≤ 2/13 * (N : ℝ) + 2 * (β.natAbs : ℝ) := by
  have hA : 0 < β.natAbs := Int.natAbs_pos.2 hβ
  have hNpos : 0 < N := by rw [hN]; positivity
  have hNR : (0 : ℝ) < (N : ℝ) := by exact_mod_cast hNpos
  have hN'R : (0 : ℝ) < (N' : ℝ) := by exact_mod_cast hN'
  have hAR : (0 : ℝ) < (β.natAbs : ℝ) := by exact_mod_cast hA
  have hNfac : (N : ℝ) = (β.natAbs : ℝ) * (N' : ℝ) := by exact_mod_cast hN
  obtain ⟨θ', hθ'⟩ : ∃ x : ℝ, x = (if 0 < β then θ else -θ) := ⟨_, rfl⟩
  -- key reduction: the value at g depends only on g % N', at unit speed and phase θ'
  have key : ∀ g : ℕ, g < N →
      distZ (θ + (β : ℝ) * (g : ℝ) / (N : ℝ))
        = distZ (θ' + ((g % N' : ℕ) : ℝ) / (N' : ℝ)) := by
    intro g _
    have hgR : (g : ℝ) = (N' : ℝ) * ((g / N' : ℕ) : ℝ) + ((g % N' : ℕ) : ℝ) := by
      exact_mod_cast congrArg (fun n : ℕ => (n : ℝ)) (Nat.div_add_mod g N').symm
    rcases lt_or_gt_of_ne hβ with hneg | hpos
    · have hβR : (β : ℝ) = -((β.natAbs : ℝ)) := by
        exact_mod_cast congrArg (fun z : ℤ => (z : ℝ)) (by omega : β = -(β.natAbs : ℤ))
      have hθ'v : θ' = -θ := by rw [hθ', if_neg (by omega : ¬ (0 < β))]
      have hval : θ + (β : ℝ) * (g : ℝ) / (N : ℝ)
          = (θ - ((g % N' : ℕ) : ℝ) / (N' : ℝ)) + ((-(g / N' : ℕ) : ℤ) : ℝ) := by
        rw [Int.cast_neg, Int.cast_natCast, hgR, hNfac, hβR]
        field_simp
        ring
      rw [hval, distZ_add_int, ← distZ_neg (θ - ((g % N' : ℕ) : ℝ) / (N' : ℝ)), hθ'v]
      congr 1
      ring
    · have hβR : (β : ℝ) = ((β.natAbs : ℝ)) := by
        have h1 : ((β.natAbs : ℕ) : ℤ) = β := Int.natAbs_of_nonneg hpos.le
        have h2 : (((β.natAbs : ℕ) : ℤ) : ℝ) = ((β.natAbs : ℕ) : ℝ) := Int.cast_natCast _
        rw [← h2, h1]
      have hθ'v : θ' = θ := by rw [hθ', if_pos hpos]
      have hval : θ + (β : ℝ) * (g : ℝ) / (N : ℝ)
          = (θ + ((g % N' : ℕ) : ℝ) / (N' : ℝ)) + (((g / N' : ℕ) : ℤ) : ℝ) := by
        rw [Int.cast_natCast, hgR, hNfac, hβR]
        field_simp
        ring
      rw [hval, distZ_add_int, hθ'v]
  -- fiber the count over g % N'
  have hmaps : ∀ g ∈ (Finset.range N).filter
      (fun g : ℕ => distZ (θ + (β : ℝ) * (g : ℝ) / (N : ℝ)) < 1/13),
      g % N' ∈ (Finset.range N').filter
        (fun r : ℕ => distZ (θ' + (r : ℝ) / (N' : ℝ)) < 1/13) := by
    intro g hg
    have hg' := Finset.mem_filter.1 hg
    refine Finset.mem_filter.2 ⟨Finset.mem_range.2 (Nat.mod_lt g hN'), ?_⟩
    have := hg'.2
    rwa [key g (Finset.mem_range.1 hg'.1)] at this
  have hcount := Finset.card_eq_sum_card_fiberwise hmaps
  have hfiber : ∀ r ∈ (Finset.range N').filter
      (fun r : ℕ => distZ (θ' + (r : ℝ) / (N' : ℝ)) < 1/13),
      (((Finset.range N).filter
        (fun g : ℕ => distZ (θ + (β : ℝ) * (g : ℝ) / (N : ℝ)) < 1/13)).filter
          (fun g => g % N' = r)).card ≤ β.natAbs := by
    intro r _
    have hsub : (((Finset.range N).filter
        (fun g : ℕ => distZ (θ + (β : ℝ) * (g : ℝ) / (N : ℝ)) < 1/13)).filter
          (fun g => g % N' = r)).card ≤ (Finset.range β.natAbs).card := by
      refine Finset.card_le_card_of_injOn (fun g => g / N') ?_ ?_
      · intro g hg
        have hg1 := (Finset.mem_filter.1 (Finset.mem_filter.1 hg).1).1
        refine Finset.mem_range.2 ?_
        have hgN : g < N := Finset.mem_range.1 hg1
        have : g < β.natAbs * N' := by rw [← hN]; exact hgN
        exact Nat.div_lt_of_lt_mul (by rw [Nat.mul_comm] at this; exact this)
      · intro g₁ hg₁ g₂ hg₂ hdiv
        have hdiv' : g₁ / N' = g₂ / N' := hdiv
        have hm₁ : g₁ % N' = r := (Finset.mem_filter.1 hg₁).2
        have hm₂ : g₂ % N' = r := (Finset.mem_filter.1 hg₂).2
        calc g₁ = N' * (g₁ / N') + g₁ % N' := (Nat.div_add_mod g₁ N').symm
          _ = N' * (g₂ / N') + g₂ % N' := by rw [hdiv', hm₁, hm₂]
          _ = g₂ := Nat.div_add_mod g₂ N'
    simpa using hsub
  have hsum : ((Finset.range N).filter
      (fun g : ℕ => distZ (θ + (β : ℝ) * (g : ℝ) / (N : ℝ)) < 1/13)).card
      ≤ ((Finset.range N').filter
        (fun r : ℕ => distZ (θ' + (r : ℝ) / (N' : ℝ)) < 1/13)).card * β.natAbs := by
    rw [hcount]
    calc ∑ r ∈ (Finset.range N').filter
          (fun r : ℕ => distZ (θ' + (r : ℝ) / (N' : ℝ)) < 1/13),
        (((Finset.range N).filter
          (fun g : ℕ => distZ (θ + (β : ℝ) * (g : ℝ) / (N : ℝ)) < 1/13)).filter
            (fun g => g % N' = r)).card
        ≤ ∑ _r ∈ (Finset.range N').filter
          (fun r : ℕ => distZ (θ' + (r : ℝ) / (N' : ℝ)) < 1/13), β.natAbs :=
          Finset.sum_le_sum hfiber
      _ = ((Finset.range N').filter
          (fun r : ℕ => distZ (θ' + (r : ℝ) / (N' : ℝ)) < 1/13)).card * β.natAbs := by
          rw [Finset.sum_const, smul_eq_mul]
  have hunit := count_grid_near_int θ' N'
  have hcast : (((Finset.range N).filter
      (fun g : ℕ => distZ (θ + (β : ℝ) * (g : ℝ) / (N : ℝ)) < 1/13)).card : ℝ)
      ≤ (((Finset.range N').filter
        (fun r : ℕ => distZ (θ' + (r : ℝ) / (N' : ℝ)) < 1/13)).card : ℝ) * (β.natAbs : ℝ) := by
    exact_mod_cast hsum
  calc (((Finset.range N).filter
      (fun g : ℕ => distZ (θ + (β : ℝ) * (g : ℝ) / (N : ℝ)) < 1/13)).card : ℝ)
      ≤ (((Finset.range N').filter
        (fun r : ℕ => distZ (θ' + (r : ℝ) / (N' : ℝ)) < 1/13)).card : ℝ) * (β.natAbs : ℝ) := hcast
    _ ≤ (2/13 * (N' : ℝ) + 2) * (β.natAbs : ℝ) := by nlinarith [hAR.le]
    _ = 2/13 * ((β.natAbs : ℝ) * (N' : ℝ)) + 2 * (β.natAbs : ℝ) := by ring
    _ = 2/13 * (N : ℝ) + 2 * (β.natAbs : ℝ) := by rw [← hNfac]

/-- **The u-escape survivor (finite pigeonhole).** For 13 phases `θᵢ` and offsets `bᵢ`
with `|bᵢ| ≤ B`, at most `j ≤ 6` of them nonzero, some `u ∈ [0,1)` keeps EVERY nonzero-offset
phase `θᵢ + bᵢ·u` at distance `≥ 1/13` from the integers. (Grid `{g/N}`, `N = 157·B·B!`:
total bad ≤ `(12/13)·N + 12·B < N`.) -/
lemma exists_good_u (θ : Fin 13 → ℝ) (b : Fin 13 → ℤ) (B : ℤ) (hB : 0 < B)
    (hb : ∀ i, |b i| ≤ B)
    (hj : ((Finset.univ : Finset (Fin 13)).filter (fun i => b i ≠ 0)).card ≤ 6) :
    ∃ u : ℝ, 0 ≤ u ∧ u < 1 ∧ ∀ i, b i ≠ 0 → 1/13 ≤ distZ (θ i + (b i : ℝ) * u) := by
  obtain ⟨N, hNdef⟩ : ∃ n : ℕ, n = 157 * B.toNat * Nat.factorial B.toNat := ⟨_, rfl⟩
  have hBt : 0 < B.toNat := by omega
  have hNpos : 0 < N := by
    rw [hNdef]
    exact Nat.mul_pos (Nat.mul_pos (by norm_num) hBt) (Nat.factorial_pos _)
  have hNR : (0 : ℝ) < (N : ℝ) := by exact_mod_cast hNpos
  have hBR : (0 : ℝ) < (B : ℝ) := by exact_mod_cast hB
  -- per-impure-runner bad sets and their counts
  have hbadcard : ∀ i : Fin 13, b i ≠ 0 →
      (((Finset.range N).filter
        (fun g : ℕ => distZ (θ i + (b i : ℝ) * (g : ℝ) / (N : ℝ)) < 1/13)).card : ℝ)
      ≤ 2/13 * (N : ℝ) + 2 * (B : ℝ) := by
    intro i hi
    have hAi : 0 < (b i).natAbs := Int.natAbs_pos.2 hi
    have hdvd : (b i).natAbs ∣ N := by
      have hle : (b i).natAbs ≤ B.toNat := by
        have h := hb i
        rw [abs_le] at h
        omega
      rw [hNdef]
      exact Dvd.dvd.mul_left (Nat.dvd_factorial hAi hle) _
    obtain ⟨N', hN'⟩ := hdvd
    have hN'pos : 0 < N' := by
      rcases Nat.eq_zero_or_pos N' with h0 | h
      · exfalso; rw [h0, Nat.mul_zero] at hN'; omega
      · exact h
    have hcount := count_grid_scaled (θ i) (b i) hi N N' hN' hN'pos
    have hcast : ((b i).natAbs : ℝ) ≤ (B : ℝ) := by
      have h := hb i
      rw [abs_le] at h
      have h1 : (((b i).natAbs : ℕ) : ℤ) ≤ B := by omega
      have h2 : ((((b i).natAbs : ℕ) : ℤ) : ℝ) ≤ ((B : ℤ) : ℝ) := Int.cast_le.2 h1
      rwa [Int.cast_natCast] at h2
    linarith
  -- the union of bad sets misses some grid point
  obtain ⟨g, hgrange, hggood⟩ : ∃ g ∈ Finset.range N,
      g ∉ ((Finset.univ : Finset (Fin 13)).filter (fun i => b i ≠ 0)).biUnion
        (fun i => (Finset.range N).filter
          (fun g : ℕ => distZ (θ i + (b i : ℝ) * (g : ℝ) / (N : ℝ)) < 1/13)) := by
    by_contra hcon
    push_neg at hcon
    have hsub : Finset.range N ⊆ ((Finset.univ : Finset (Fin 13)).filter
        (fun i => b i ≠ 0)).biUnion
        (fun i => (Finset.range N).filter
          (fun g : ℕ => distZ (θ i + (b i : ℝ) * (g : ℝ) / (N : ℝ)) < 1/13)) := hcon
    have hcard := Finset.card_le_card hsub
    rw [Finset.card_range] at hcard
    have hbiUnion := Finset.card_biUnion_le
      (s := (Finset.univ : Finset (Fin 13)).filter (fun i => b i ≠ 0))
      (t := fun i => (Finset.range N).filter
        (fun g : ℕ => distZ (θ i + (b i : ℝ) * (g : ℝ) / (N : ℝ)) < 1/13))
    -- real-valued chain: N ≤ card(biUnion) ≤ Σ ≤ 6·(2/13·N + 2B) = 12/13·N + 12B < N
    have hsumR : ((∑ i ∈ (Finset.univ : Finset (Fin 13)).filter (fun i => b i ≠ 0),
        ((Finset.range N).filter
          (fun g : ℕ => distZ (θ i + (b i : ℝ) * (g : ℝ) / (N : ℝ)) < 1/13)).card : ℕ) : ℝ)
        ≤ 6 * (2/13 * (N : ℝ) + 2 * (B : ℝ)) := by
      push_cast
      calc ∑ i ∈ (Finset.univ : Finset (Fin 13)).filter (fun i => b i ≠ 0),
            (((Finset.range N).filter
              (fun g : ℕ => distZ (θ i + (b i : ℝ) * (g : ℝ) / (N : ℝ)) < 1/13)).card : ℝ)
          ≤ ∑ i ∈ (Finset.univ : Finset (Fin 13)).filter (fun i => b i ≠ 0),
            (2/13 * (N : ℝ) + 2 * (B : ℝ)) := by
            refine Finset.sum_le_sum ?_
            intro i hi
            exact hbadcard i (Finset.mem_filter.1 hi).2
        _ = (((Finset.univ : Finset (Fin 13)).filter (fun i => b i ≠ 0)).card : ℝ)
            * (2/13 * (N : ℝ) + 2 * (B : ℝ)) := by
            rw [Finset.sum_const, nsmul_eq_mul]
        _ ≤ 6 * (2/13 * (N : ℝ) + 2 * (B : ℝ)) := by
            have hjR : (((Finset.univ : Finset (Fin 13)).filter
                (fun i => b i ≠ 0)).card : ℝ) ≤ 6 := by exact_mod_cast hj
            nlinarith [hNR.le, hBR.le]
    have hNbig : 156 * (B : ℝ) < (N : ℝ) := by
      have h1 : 157 * B.toNat ≤ N := by
        rw [hNdef]
        calc 157 * B.toNat = 157 * B.toNat * 1 := by ring
          _ ≤ 157 * B.toNat * Nat.factorial B.toNat := by
              exact Nat.mul_le_mul_left _ (Nat.factorial_pos _)
      have h1Z : (157 : ℤ) * B ≤ (N : ℤ) := by omega
      have h2 : (157 : ℝ) * (B : ℝ) ≤ (N : ℝ) := by exact_mod_cast h1Z
      linarith
    have hchain : (N : ℝ) ≤ 12/13 * (N : ℝ) + 12 * (B : ℝ) := by
      have hc1 : ((N : ℕ) : ℝ) ≤ ((((Finset.univ : Finset (Fin 13)).filter
          (fun i => b i ≠ 0)).biUnion
          (fun i => (Finset.range N).filter
            (fun g : ℕ => distZ (θ i + (b i : ℝ) * (g : ℝ) / (N : ℝ)) < 1/13))).card : ℝ) := by
        exact_mod_cast hcard
      have hc2 : ((((Finset.univ : Finset (Fin 13)).filter (fun i => b i ≠ 0)).biUnion
          (fun i => (Finset.range N).filter
            (fun g : ℕ => distZ (θ i + (b i : ℝ) * (g : ℝ) / (N : ℝ)) < 1/13))).card : ℝ)
          ≤ ((∑ i ∈ (Finset.univ : Finset (Fin 13)).filter (fun i => b i ≠ 0),
            ((Finset.range N).filter
              (fun g : ℕ => distZ (θ i + (b i : ℝ) * (g : ℝ) / (N : ℝ)) < 1/13)).card : ℕ) : ℝ) := by
        exact_mod_cast hbiUnion
      linarith
    linarith
  -- the survivor
  refine ⟨(g : ℝ) / (N : ℝ), by positivity, ?_, ?_⟩
  · rw [div_lt_one hNR]
    exact_mod_cast Finset.mem_range.1 hgrange
  · intro i hi
    have hnotbad : g ∉ (Finset.range N).filter
        (fun g : ℕ => distZ (θ i + (b i : ℝ) * (g : ℝ) / (N : ℝ)) < 1/13) := by
      intro hmem
      exact hggood (Finset.mem_biUnion.2
        ⟨i, Finset.mem_filter.2 ⟨Finset.mem_univ i, hi⟩, hmem⟩)
    have : ¬ (distZ (θ i + (b i : ℝ) * (g : ℝ) / (N : ℝ)) < 1/13) := by
      intro hlt
      exact hnotbad (Finset.mem_filter.2 ⟨hgrange, hlt⟩)
    rw [not_lt] at this
    rwa [mul_div_assoc] at this

/-- **THM-721 Part 3 (witness form) — the compressed `j ≤ 6` stratum is loose.**
`Vᵢ = bᵢ + L·Kᵢ`, `|bᵢ| ≤ B`, at most 6 nonzero `bᵢ`, `L > 91·B`, and the PURE lift
sub-family admits a common `1/13`-safe time (the LRC(≤13) citation: ≤ 12 distinct pure
lifts). Then some real time has margin `> 1/14`: the family is LOOSE at floor
`1/13 − B/(2L)`. NO LRC(14) input. -/
theorem margin_uescape_j6 (V K b : Fin 13 → ℤ) (L B : ℤ) (hL : 0 < L) (hB : 0 < B)
    (hV : ∀ i, V i = b i + L * K i) (hb : ∀ i, |b i| ≤ B)
    (hj : ((Finset.univ : Finset (Fin 13)).filter (fun i => b i ≠ 0)).card ≤ 6)
    (hLB : 91 * B < L)
    (hpure : ∃ s : ℝ, ∀ i, b i = 0 → 1/13 ≤ distZ ((K i : ℝ) * s)) :
    ∃ t : ℝ, (1 : ℝ) / 14 < margin V t := by
  obtain ⟨s, hs⟩ := hpure
  obtain ⟨u, hu0, hu1, hugood⟩ := exists_good_u (fun i => (K i : ℝ) * s) b B hB hb hj
  have hLR : (0 : ℝ) < (L : ℝ) := by exact_mod_cast hL
  have hBR : (0 : ℝ) < (B : ℝ) := by exact_mod_cast hB
  obtain ⟨m, hm⟩ : ∃ z : ℤ, z = round ((L : ℝ) * u - s) := ⟨_, rfl⟩
  have hround : |(L : ℝ) * u - s - (m : ℝ)| ≤ 1/2 := by
    rw [hm]
    exact abs_sub_round _
  refine ⟨((m : ℝ) + s) / (L : ℝ), ?_⟩
  have hclose : |((m : ℝ) + s) / (L : ℝ) - u| ≤ 1 / (2 * (L : ℝ)) := by
    have hval : ((m : ℝ) + s) / (L : ℝ) - u = -(((L : ℝ) * u - s - (m : ℝ)) / (L : ℝ)) := by
      field_simp
      ring
    rw [hval, abs_neg, abs_div, abs_of_pos hLR]
    rw [div_le_div_iff₀ hLR (by positivity)]
    calc |(L : ℝ) * u - s - (m : ℝ)| * (2 * (L : ℝ))
        ≤ (1/2) * (2 * (L : ℝ)) := by nlinarith [hLR.le]
      _ = 1 * (L : ℝ) := by ring
  -- per-runner floor
  have hfloor : ∀ i : Fin 13,
      1/13 - (B : ℝ) / (2 * (L : ℝ)) ≤ distZ ((V i : ℝ) * (((m : ℝ) + s) / (L : ℝ))) := by
    intro i
    have hexp : (V i : ℝ) * (((m : ℝ) + s) / (L : ℝ))
        = ((K i : ℝ) * s + (b i : ℝ) * (((m : ℝ) + s) / (L : ℝ))) + ((K i * m : ℤ) : ℝ) := by
      have hVi : ((V i : ℤ) : ℝ) = (b i : ℝ) + (L : ℝ) * (K i : ℝ) := by
        rw [hV i]; push_cast; ring
      rw [hVi]
      push_cast
      field_simp
      ring
    rw [hexp, distZ_add_int]
    rcases eq_or_ne (b i) 0 with hbi | hbi
    · -- pure runner: the fast coordinate is exact, no loss
      have h1 := hs i hbi
      rw [hbi]
      push_cast
      simp only [zero_mul, add_zero]
      have hpos : (0 : ℝ) < (B : ℝ) / (2 * (L : ℝ)) := by positivity
      linarith
    · -- impure runner: pay the Lipschitz loss on the slow coordinate
      have h2 := hugood i hbi
      have hlip := distZ_lipschitz ((K i : ℝ) * s + (b i : ℝ) * (((m : ℝ) + s) / (L : ℝ)))
        ((K i : ℝ) * s + (b i : ℝ) * u)
      have hdiff : |((K i : ℝ) * s + (b i : ℝ) * (((m : ℝ) + s) / (L : ℝ)))
          - ((K i : ℝ) * s + (b i : ℝ) * u)| ≤ (B : ℝ) / (2 * (L : ℝ)) := by
        rw [show ((K i : ℝ) * s + (b i : ℝ) * (((m : ℝ) + s) / (L : ℝ)))
            - ((K i : ℝ) * s + (b i : ℝ) * u)
            = (b i : ℝ) * ((((m : ℝ) + s) / (L : ℝ)) - u) from by ring]
        rw [abs_mul]
        have hbiB : |(b i : ℝ)| ≤ (B : ℝ) := by
          have := hb i
          exact_mod_cast this
        calc |(b i : ℝ)| * |(((m : ℝ) + s) / (L : ℝ)) - u|
            ≤ (B : ℝ) * (1 / (2 * (L : ℝ))) := by
              apply mul_le_mul hbiB hclose (abs_nonneg _) hBR.le
          _ = (B : ℝ) / (2 * (L : ℝ)) := by ring
      linarith
  -- assemble the margin
  have hmargin : 1/13 - (B : ℝ) / (2 * (L : ℝ)) ≤ margin V (((m : ℝ) + s) / (L : ℝ)) := by
    rw [margin, Finset.le_inf'_iff]
    intro i _
    exact hfloor i
  have hgap : (1 : ℝ) / 14 < 1/13 - (B : ℝ) / (2 * (L : ℝ)) := by
    have hLBR : (91 : ℝ) * (B : ℝ) < (L : ℝ) := by exact_mod_cast hLB
    have hBL : (B : ℝ) / (2 * (L : ℝ)) < 1/182 := by
      rw [div_lt_iff₀ (by positivity : (0:ℝ) < 2 * (L : ℝ))]
      linarith
    linarith
  linarith [hmargin, hgap]

/-- Corollary in kps's reach form: the compressed `j ≤ 6` family has
`sSup (margin V '' [0,1]) > 1/14` (wrap the witness into `[0,1]` by periodicity). -/
theorem reach_uescape_j6 (V K b : Fin 13 → ℤ) (L B : ℤ) (hL : 0 < L) (hB : 0 < B)
    (hV : ∀ i, V i = b i + L * K i) (hb : ∀ i, |b i| ≤ B)
    (hj : ((Finset.univ : Finset (Fin 13)).filter (fun i => b i ≠ 0)).card ≤ 6)
    (hLB : 91 * B < L)
    (hpure : ∃ s : ℝ, ∀ i, b i = 0 → 1/13 ≤ distZ ((K i : ℝ) * s)) :
    (1 : ℝ) / 14 < sSup (margin V '' Set.Icc (0 : ℝ) 1) := by
  obtain ⟨t, ht⟩ := margin_uescape_j6 V K b L B hL hB hV hb hj hLB hpure
  -- wrap t into [0,1] by integer shifts: margin V (t - ⌊t⌋) = margin V t
  obtain ⟨t', ht'mem, ht'val⟩ : ∃ t' ∈ Set.Icc (0:ℝ) 1, margin V t' = margin V t := by
    refine ⟨t - (⌊t⌋ : ℝ), ?_, ?_⟩
    · constructor
      · linarith [Int.floor_le t]
      · linarith [Int.lt_floor_add_one t]
    · -- margin is 1-periodic; induct the integer shift via margin_periodic
      have hshift : ∀ n : ℕ, ∀ x : ℝ, margin V (x + (n : ℝ)) = margin V x := by
        intro n
        induction n with
        | zero => intro x; norm_num
        | succ k ih =>
            intro x
            have h1 : x + ((k+1 : ℕ) : ℝ) = (x + (k : ℕ)) + 1 := by push_cast; ring
            rw [h1, margin_periodic, ih]
      rcases le_or_gt 0 ⌊t⌋ with hf | hf
      · obtain ⟨n, hn⟩ : ∃ n : ℕ, (n : ℤ) = ⌊t⌋ := ⟨⌊t⌋.toNat, by omega⟩
        have : t = (t - (⌊t⌋ : ℝ)) + (n : ℝ) := by
          have : ((n : ℤ) : ℝ) = (⌊t⌋ : ℝ) := by exact_mod_cast hn
          push_cast at this ⊢
          linarith
        rw [show margin V t = margin V ((t - (⌊t⌋ : ℝ)) + (n : ℝ)) from by rw [← this]]
        rw [hshift n]
      · obtain ⟨n, hn⟩ : ∃ n : ℕ, (n : ℤ) = -⌊t⌋ := ⟨(-⌊t⌋).toNat, by omega⟩
        have hteq : t - (⌊t⌋ : ℝ) = t + (n : ℝ) := by
          have : ((n : ℤ) : ℝ) = -(⌊t⌋ : ℝ) := by exact_mod_cast congrArg (fun z : ℤ => (z:ℝ)) hn
          push_cast at this ⊢
          linarith
        rw [hteq, hshift n]
  have hbdd : BddAbove (margin V '' Set.Icc (0 : ℝ) 1) :=
    ⟨1 / 2, by rintro _ ⟨x, _, rfl⟩; exact LRCDecorr13.margin_le_half13 V x⟩
  have hle : margin V t' ≤ sSup (margin V '' Set.Icc (0 : ℝ) 1) :=
    le_csSup hbdd ⟨t', ht'mem, rfl⟩
  rw [ht'val] at hle
  linarith

#print axioms distZ_neg
#print axioms card_range_filter_window
#print axioms count_grid_near_int
#print axioms count_grid_scaled
#print axioms exists_good_u
#print axioms margin_uescape_j6
#print axioms reach_uescape_j6

end TournamentH7.LRCUEscape
