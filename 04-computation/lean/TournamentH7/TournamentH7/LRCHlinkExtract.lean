/-
  TournamentH7.LRCHlinkExtract — the mergeSort argmax + wrapping-gap extraction for klein-S203's
  `hlink` (opus-2026-07-09-S175).

  kps-S101 (`LRCTeethGap`) discharged the NON-WRAPPING freeness core of `hlink`.  This file supplies
  the remaining two pieces the owner named:
    (1) the `foldl max` ARGMAX extraction — `maxCircGap = G > 0` is achieved by an ADJACENT pair of the
        sorted-and-closed residue cycle `cyc = ps ++ [p0+Vmax]`;
    (2) the WRAPPING-gap case — when the max gap is the wraparound `(ps.last, p0+Vmax)`, the arc past `1`
        is still free of every integer tooth-translate (uses that `ps.last` is the max and `p0` the min).

  Built bottom-up: the `foldl max` membership helpers first (self-contained, reusable), then the
  adjacent-pair extraction, then the two freeness branches, then `hlink`.
-/
import Mathlib
import TournamentH7.LRCTeethGap

namespace LRC14

open LonelyRunner

/-! ### 1. `foldl max` argmax helpers -/

/-- The `foldl max` of a `ℕ`-list is either the initial accumulator or a member of the list. -/
theorem foldl_max_eq_or_mem : ∀ (L : List ℕ) (a : ℕ), L.foldl max a = a ∨ L.foldl max a ∈ L := by
  intro L
  induction L with
  | nil => intro a; left; rfl
  | cons x xs ih =>
    intro a
    simp only [List.foldl_cons]
    rcases ih (max a x) with h | h
    · rw [h]
      rcases max_choice a x with hx | hx
      · left; exact hx
      · right; rw [hx]; exact List.mem_cons_self
    · right; exact List.mem_cons_of_mem x h

/-- A POSITIVE `foldl max 0` is achieved by a member of the list. -/
theorem foldl_max_mem (L : List ℕ) (h : 0 < L.foldl max 0) : L.foldl max 0 ∈ L := by
  rcases foldl_max_eq_or_mem L 0 with h0 | hm
  · rw [h0] at h; exact absurd h (lt_irrefl 0)
  · exact hm

/-- Every element `d` of `List.zipWith (fun a b => b - a) c c.tail` is `y − x` for an ADJACENT pair
`x, y` of `c`: `c = l₁ ++ x :: y :: l₂` with `y − x = d`. -/
theorem mem_zipWith_sub_tail : ∀ (c : List ℕ) (d : ℕ),
    d ∈ List.zipWith (fun a b => b - a) c c.tail →
    ∃ (l₁ l₂ : List ℕ) (x y : ℕ), c = l₁ ++ x :: y :: l₂ ∧ y - x = d := by
  intro c
  induction c with
  | nil => intro d hd; simp at hd
  | cons x rest ihx =>
    cases rest with
    | nil => intro d hd; simp at hd
    | cons y rest2 =>
      intro d hd
      simp only [List.tail_cons, List.zipWith_cons_cons, List.mem_cons] at hd
      rcases hd with h | h
      · exact ⟨[], rest2, x, y, by simp, h.symm⟩
      · have h' : d ∈ List.zipWith (fun a b => b - a) (y :: rest2) (y :: rest2).tail := by
          simpa using h
        obtain ⟨l₁, l₂, x', y', heq, hd'⟩ := ihx d h'
        exact ⟨x :: l₁, l₂, x', y', by rw [heq]; simp, hd'⟩

/-- From a positive `maxCircGap` (nonempty residue list `ps = p0 :: rest`), the gap is realized by an
ADJACENT pair `x, y` of the closed cycle `cyc = ps ++ [p0 + Vmax]`. -/
theorem exists_gap_decomp (E : List ℕ) (Vmax j p0 : ℕ) (rest : List ℕ)
    (hps : (E.map (fun e => e * j % Vmax)).mergeSort (· ≤ ·) = p0 :: rest)
    (hpos : 0 < maxCircGap E Vmax j) :
    ∃ (l₁ l₂ : List ℕ) (x y : ℕ),
      (p0 :: rest ++ [p0 + Vmax]) = l₁ ++ x :: y :: l₂ ∧ y - x = maxCircGap E Vmax j := by
  have hval : maxCircGap E Vmax j
      = List.foldl max 0 (List.zipWith (fun a b => b - a)
          (p0 :: rest ++ [p0 + Vmax]) (p0 :: rest ++ [p0 + Vmax]).tail) := by
    simp only [maxCircGap, hps]
  rw [hval] at hpos ⊢
  exact mem_zipWith_sub_tail _ _ (foldl_max_mem _ hpos)

/-! ### 2. Sortedness / residue facts of the closed cycle -/

/-- Every residue lies in `[0, Vmax)`. -/
theorem residue_lt (E : List ℕ) (Vmax j r : ℕ) (hV : 0 < Vmax)
    (hr : r ∈ (E.map (fun e => e * j % Vmax)).mergeSort (· ≤ ·)) : r < Vmax := by
  have hperm := (List.mergeSort_perm (E.map (fun e => e * j % Vmax)) (· ≤ ·))
  have hr' : r ∈ E.map (fun e => e * j % Vmax) := hperm.mem_iff.mp hr
  simp only [List.mem_map] at hr'
  obtain ⟨e, -, rfl⟩ := hr'
  exact Nat.mod_lt _ hV

/-- The sorted residues are `Pairwise (≤)`. -/
theorem sorted_ps (E : List ℕ) (Vmax j : ℕ) :
    List.Pairwise (· ≤ ·) ((E.map (fun e => e * j % Vmax)).mergeSort (· ≤ ·)) := by
  exact List.pairwise_mergeSort' (r := (· ≤ ·)) _

/-! ### 3. `hlink` — the free residue gap from a good period -/

/-- **The free-gap extraction (`hlink`).**  A good period `j` (`Vmax < 7·maxCircGap`) yields a real
interval `(a, a+g)` of length `g > 1/7` containing no integer translate of any tooth — the missing link
of `LRCGoodPeriodReach.mreach_ge_of_goodPeriod` (klein-S203).  Handles both the interior gap (kps-S101
`free_translate_of_free_subInterval`) and the WRAPPING gap `(ps.last, p0+Vmax)` (directly). -/
theorem hlink_extract (E : List ℕ) (Vmax : ℕ) (hV : 0 < Vmax)
    (hgp : HasGoodPeriod E Vmax) :
    ∃ (j : ℕ) (a g : ℝ), 1 / 7 < g ∧
      ∀ c ∈ teeth E Vmax j, ∀ n : ℤ, (c + (n : ℝ)) ∉ Set.Ioo a (a + g) := by
  obtain ⟨j, -, hgood⟩ := hgp
  rw [LRC14.IsGoodPeriod] at hgood            -- hgood : Vmax < 7 * maxCircGap E Vmax j
  have hVR : (0 : ℝ) < (Vmax : ℝ) := by exact_mod_cast hV
  -- residue list empty ⟹ teeth empty ⟹ vacuous
  rcases hps_cases : (E.map (fun e => e * j % Vmax)).mergeSort (· ≤ ·) with _ | ⟨p0, rest⟩
  · refine ⟨j, 0, 1, by norm_num, ?_⟩
    intro c hc
    exfalso
    have hlen : E.length = 0 := by
      have h := (List.mergeSort_perm (E.map (fun e => e * j % Vmax)) (· ≤ ·)).length_eq
      rw [hps_cases] at h
      simp only [List.length_nil, List.length_map] at h
      omega
    rw [List.length_eq_zero_iff] at hlen
    subst hlen
    simp [teeth] at hc
  · -- nonempty: extract the adjacent gap pair
    have hpos : 0 < maxCircGap E Vmax j := by omega
    obtain ⟨l₁, l₂, x, y, hdecomp, hxy⟩ :=
      exists_gap_decomp E Vmax j p0 rest hps_cases hpos
    -- sortedness of cyc = ps ++ [p0+Vmax]
    have hsps : List.Pairwise (· ≤ ·) (p0 :: rest) := by
      rw [← hps_cases]; exact sorted_ps E Vmax j
    have hall_lt : ∀ r ∈ (p0 :: rest), r < Vmax := by
      intro r hr; exact residue_lt E Vmax j r hV (by rw [hps_cases]; exact hr)
    have hsorted_cyc : List.Pairwise (· ≤ ·) (p0 :: rest ++ [p0 + Vmax]) := by
      rw [List.pairwise_append]
      refine ⟨hsps, List.pairwise_singleton _ _, ?_⟩
      intro a ha b hb
      simp only [List.mem_singleton] at hb; subst hb
      exact le_of_lt (lt_of_lt_of_le (hall_lt a ha) (by omega))
    rw [hdecomp] at hsorted_cyc
    -- unpack pairwise-≤ of l₁ ++ x :: y :: l₂
    rw [List.pairwise_append] at hsorted_cyc
    obtain ⟨-, hsm, hl₁x⟩ := hsorted_cyc
    have hxley : x ≤ y := by
      rw [List.pairwise_cons] at hsm
      exact hsm.1 y (by simp)
    have hyl₂ : ∀ b ∈ l₂, y ≤ b := by
      rw [List.pairwise_cons, List.pairwise_cons] at hsm
      exact hsm.2.1
    have hl₁x' : ∀ a ∈ l₁, a ≤ x := by
      intro a ha; exact hl₁x a ha x (by simp)
    -- residue freeness: no r ∈ ps with x < r < y
    have hfree_res : ∀ r ∈ (p0 :: rest), ¬ (x < r ∧ r < y) := by
      intro r hr ⟨hxr, hry⟩
      have hrcyc : r ∈ l₁ ++ x :: y :: l₂ := by rw [← hdecomp]; exact List.mem_append_left _ hr
      rw [List.mem_append, List.mem_cons, List.mem_cons] at hrcyc
      rcases hrcyc with h | h | h | h
      · exact absurd (hl₁x' r h) (by omega)
      · omega
      · omega
      · exact absurd (hyl₂ r h) (by omega)
    -- g = (y - x)/Vmax > 1/7
    have hgap : y - x = maxCircGap E Vmax j := hxy
    have hgpos : 1 / 7 < ((y - x : ℕ) : ℝ) / (Vmax : ℝ) := by
      rw [lt_div_iff₀ hVR, hgap]
      have : (Vmax : ℝ) < 7 * ((maxCircGap E Vmax j : ℕ) : ℝ) := by exact_mod_cast hgood
      push_cast at this ⊢; linarith
    refine ⟨j, (x : ℝ) / Vmax, ((y - x : ℕ) : ℝ) / Vmax, hgpos, ?_⟩
    -- teeth membership: c = r/Vmax with r ∈ ps
    have htooth : ∀ c ∈ teeth E Vmax j, ∃ r ∈ (p0 :: rest), c = (r : ℝ) / Vmax := by
      intro c hc
      simp only [teeth, List.mem_toFinset, List.mem_map] at hc
      obtain ⟨e, he, rfl⟩ := hc
      refine ⟨e * j % Vmax, ?_, rfl⟩
      have hperm := List.mergeSort_perm (E.map (fun e => e * j % Vmax)) (· ≤ ·)
      rw [hps_cases] at hperm
      exact hperm.mem_iff.mpr (List.mem_map_of_mem he)
    have haeq : (x : ℝ) / Vmax + ((y - x : ℕ) : ℝ) / Vmax = (y : ℝ) / Vmax := by
      rw [← add_div]; congr 1
      have : ((y - x : ℕ) : ℝ) = (y : ℝ) - (x : ℝ) := by
        rw [Nat.cast_sub hxley]
      rw [this]; ring
    -- CASE: y is a residue (interior) or y = p0 + Vmax (wraparound)
    have hy_mem : y ∈ (p0 :: rest) ++ [p0 + Vmax] := by rw [hdecomp]; simp
    rw [List.mem_append, List.mem_singleton] at hy_mem
    rcases hy_mem with hy_res | hy_wrap
    · -- INTERIOR: y < Vmax, non-wrapping, use kps-S101
      have hyV : y < Vmax := hall_lt y hy_res
      apply free_translate_of_free_subInterval (teeth E Vmax j) _ _
        (teeth_subset_Ico E Vmax j hV) (by positivity)
      · rw [haeq, div_le_one hVR]; exact_mod_cast le_of_lt hyV
      · intro c hc hmem
        obtain ⟨r, hr, rfl⟩ := htooth c hc
        rw [haeq] at hmem
        obtain ⟨hlo, hhi⟩ := hmem
        have hxr : x < r := by
          by_contra hcon; push_neg at hcon
          have hcast : (r : ℝ) ≤ (x : ℝ) := by exact_mod_cast hcon
          have : (r : ℝ) / Vmax ≤ (x : ℝ) / Vmax := by gcongr
          linarith
        have hry : r < y := by
          by_contra hcon; push_neg at hcon
          have hcast : (y : ℝ) ≤ (r : ℝ) := by exact_mod_cast hcon
          have : (y : ℝ) / Vmax ≤ (r : ℝ) / Vmax := by gcongr
          linarith
        exact hfree_res r hr ⟨hxr, hry⟩
    · -- WRAPAROUND: y = p0 + Vmax
      subst hy_wrap
      intro c hc n hmem
      obtain ⟨r, hr, rfl⟩ := htooth c hc
      obtain ⟨hlo, hhi⟩ := hmem
      -- bounds: p0 ≤ r ≤ x
      have hp0r : p0 ≤ r := by
        rw [List.pairwise_cons] at hsps
        have hr' := hr; rw [List.mem_cons] at hr'
        rcases hr' with h | h
        · omega
        · exact hsps.1 r h
      have hrx : r ≤ x := by
        have hrcyc : r ∈ l₁ ++ x :: (p0 + Vmax) :: l₂ := by
          rw [← hdecomp]; exact List.mem_append_left _ hr
        rw [List.mem_append, List.mem_cons, List.mem_cons] at hrcyc
        rcases hrcyc with h | h | h | h
        · exact hl₁x' r h
        · omega
        · have := hall_lt r hr; omega
        · have := hyl₂ r h; have := hall_lt r hr; omega
      -- a = x/Vmax, a+g = (p0+Vmax)/Vmax = p0/Vmax + 1
      rw [haeq] at hhi
      -- hlo : x/Vmax < r/Vmax + n ; hhi : r/Vmax + n < (p0+Vmax)/Vmax = p0/Vmax + 1
      have hsplit : ((p0 + Vmax : ℕ) : ℝ) / (Vmax : ℝ) = (p0 : ℝ) / Vmax + 1 := by
        push_cast; field_simp
      have hrx_div : (r : ℝ) / Vmax ≤ (x : ℝ) / Vmax := by
        have : (r : ℝ) ≤ (x : ℝ) := by exact_mod_cast hrx
        gcongr
      have hp0_div : (p0 : ℝ) / Vmax ≤ (r : ℝ) / Vmax := by
        have : (p0 : ℝ) ≤ (r : ℝ) := by exact_mod_cast hp0r
        gcongr
      rw [hsplit] at hhi
      have hn_pos : (0 : ℝ) < (n : ℝ) := by linarith [hlo, hrx_div]
      have hn_lt : (n : ℝ) < 1 := by linarith [hhi, hp0_div]
      have h1 : (0 : ℤ) < n := by exact_mod_cast hn_pos
      have h2 : n < 1 := by exact_mod_cast hn_lt
      omega

/-- **`hlink` discharged.**  With the free-gap extraction proved (`hlink_extract`), klein-S203's
good-period → reach chain (`mreach_ge_of_goodPeriod`) needs ONLY the ruler embedding `hembed`
(THM-527 Part A, the shared blocker).  This closes the `hlink` link of the good-period leg. -/
theorem mreach_ge_of_goodPeriod_of_hembed (E : List ℕ) (Vmax : ℕ) (v : Fin 13 → ℤ) (hV : 0 < Vmax)
    (hgp : HasGoodPeriod E Vmax)
    (hembed : (∃ (j : ℕ) (φ : ℝ), ∀ c ∈ teeth E Vmax j,
        1 / 14 < LonelyRunner.LRC14Concrete.nearInt (φ - c)) →
      ∃ τ : ℝ, τ ∈ Set.Icc (0 : ℝ) 1 ∧ (1 : ℝ) / 14 ≤ LonelyRunner.LRC14Concrete.minReach v τ) :
    (1 : ℝ) / 14 ≤ LonelyRunner.LRC14Concrete.Mreach v :=
  mreach_ge_of_goodPeriod E Vmax v hgp (hlink_extract E Vmax hV) hembed

end LRC14
