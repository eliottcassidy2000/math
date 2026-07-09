/-
  TournamentH7.LRCGoodPeriodFreeGap — discharging `hlink` (klein-2026-07-09-S204).

  klein-S203 `LRCGoodPeriodReach.mreach_ge_of_goodPeriod` left two links; kps-S101 `LRCTeethGap`
  discharged the geometric core of the FIRST (`hlink`) on the non-wrapping branch.  This file finishes
  `hlink`: the **mergeSort argmax extraction** (a good period's maximal circular gap IS a tooth-free
  consecutive-residue interval) and the **wrapping case**.

  KEY SIMPLIFICATION of the wrap (the LRC co-offset convention `0 ∈ E`): the residue of `0` is `0`, the
  minimum, so `p0 = ps.head = 0` and `cyc.last = p0 + Vmax = Vmax`.  Hence EVERY consecutive gap
  interval `(cyc[i]/Vmax, cyc[i+1]/Vmax)` has right endpoint `≤ Vmax/Vmax = 1` — a subinterval of
  `[0,1]`.  So kps-S101's non-wrapping `free_translate_of_free_subInterval` covers the wrap too; no
  separate cyclic argument is needed.  Self-contained on `LRCTeethGap`.
-/
import Mathlib
import TournamentH7.LRCTeethGap

namespace LRC14

open LonelyRunner

/-- `l.foldl max a` is either the seed `a` or a member of `l`. -/
theorem foldl_max_mem : ∀ (l : List ℕ) (a : ℕ), l.foldl max a = a ∨ l.foldl max a ∈ l := by
  intro l
  induction l with
  | nil => intro a; left; rfl
  | cons x xs ih =>
    intro a
    simp only [List.foldl_cons]
    rcases ih (max a x) with h | h
    · rw [h]; rcases max_choice a x with hc | hc
      · left; exact hc
      · right; rw [hc]; simp
    · right; exact List.mem_cons_of_mem x h

/-- **The good period gives a free integer gap.**  For a good period `j` (with `0 ∈ E`), the maximal
circular gap yields naturals `lo ≤ hi ≤ Vmax` with `Vmax < 7·(hi−lo)` such that NO residue
`e·j mod Vmax` lies strictly between `lo` and `hi`.  The interval `(lo/Vmax, hi/Vmax) ⊆ [0,1]` is the
tooth-free window feeding `GapReach`. -/
theorem goodPeriod_intFreeGap (E : List ℕ) (Vmax j : ℕ) (hV : 0 < Vmax) (h0 : 0 ∈ E)
    (hgp : IsGoodPeriod E Vmax j) :
    ∃ lo hi : ℕ, lo ≤ hi ∧ hi ≤ Vmax ∧ Vmax < 7 * (hi - lo) ∧
      ∀ e ∈ E, ¬ (lo < e * j % Vmax ∧ e * j % Vmax < hi) := by
  classical
  set res : List ℕ := E.map (fun e => e * j % Vmax) with hres
  set ps : List ℕ := res.mergeSort (· ≤ ·) with hps
  -- residues are sorted
  have hsorted : ps.Pairwise (· ≤ ·) := by
    have := List.pairwise_mergeSort (le := (· ≤ ·))
      (by intro a b c; simp; omega) (by intro a b; simp; omega) res
    simpa [hps] using this
  -- 0 is a residue, hence in ps
  have h0res : (0 : ℕ) ∈ res := by
    rw [hres]; exact List.mem_map.mpr ⟨0, h0, by simp⟩
  have h0ps : (0 : ℕ) ∈ ps := by rw [hps]; exact (List.mem_mergeSort).mpr h0res
  -- ps ≠ []
  have hpsne : ps ≠ [] := by intro h; rw [h] at h0ps; exact List.not_mem_nil h0ps
  -- residues are < Vmax
  have hreslt : ∀ r ∈ res, r < Vmax := by
    intro r hr; rw [hres] at hr; obtain ⟨e, -, rfl⟩ := List.mem_map.mp hr; exact Nat.mod_lt _ hV
  have hpslt : ∀ r ∈ ps, r < Vmax := by
    intro r hr; exact hreslt r ((List.mem_mergeSort).mp (by rwa [hps] at hr))
  -- head of ps is 0 (minimum, since 0 ∈ ps and ps sorted)
  obtain ⟨p0, rest, hpeq⟩ : ∃ p0 rest, ps = p0 :: rest := by
    cases hp : ps with
    | nil => exact absurd hp hpsne
    | cons a l => exact ⟨a, l, rfl⟩
  have hp0 : p0 = 0 := by
    -- p0 ≤ 0 : head ≤ every element, and 0 ∈ ps
    have hle : ∀ r ∈ ps, p0 ≤ r := by
      intro r hr
      rw [hpeq] at hr hsorted
      rcases List.mem_cons.mp hr with h | h
      · omega
      · exact (List.pairwise_cons.mp hsorted).1 r h
    have : p0 ≤ 0 := hle 0 h0ps
    omega
  -- the cyclic list
  set cyc : List ℕ := ps ++ [p0 + Vmax] with hcyc
  -- cyc is sorted
  have hcycsorted : cyc.Pairwise (· ≤ ·) := by
    rw [hcyc]
    rw [List.pairwise_append]
    refine ⟨hsorted, by simp, ?_⟩
    intro a ha b hb
    simp only [List.mem_singleton] at hb
    subst hb
    have := hpslt a ha
    omega
  -- length facts
  have hcyclen : cyc.length = ps.length + 1 := by rw [hcyc]; simp
  -- get-monotone on cyc
  have hmono : ∀ (a b : ℕ) (ha : a < cyc.length) (hb : b < cyc.length), a < b →
      cyc[a] ≤ cyc[b] := by
    intro a b ha hb hab
    exact (List.pairwise_iff_getElem.mp hcycsorted) a b ha hb hab
  -- every cyc element ≤ Vmax (residues < Vmax, and the wrap term p0+Vmax = Vmax)
  have hcyc_mem_le : ∀ x ∈ cyc, x ≤ Vmax := by
    intro x hx
    rw [hcyc] at hx
    rcases List.mem_append.mp hx with h | h
    · exact le_of_lt (hpslt x h)
    · simp only [List.mem_singleton] at h; omega
  have hcycle : ∀ (a : ℕ) (ha : a < cyc.length), cyc[a] ≤ Vmax := by
    intro a ha; exact hcyc_mem_le _ (List.getElem_mem _)
  -- the gap list and the max
  set gaps : List ℕ := List.zipWith (fun a b => b - a) cyc cyc.tail with hgaps
  have hGdef : maxCircGap E Vmax j = gaps.foldl max 0 := by
    rw [hgaps, hcyc, maxCircGap, ← hres, ← hps, hpeq]
  have hG : Vmax < 7 * (gaps.foldl max 0) := by rw [← hGdef]; exact hgp
  -- the max gap is a member of gaps (positive, so not the seed)
  have hGpos : 0 < gaps.foldl max 0 := by omega
  have hGmem : gaps.foldl max 0 ∈ gaps := by
    rcases foldl_max_mem gaps 0 with h | h
    · omega
    · exact h
  -- extract the index
  obtain ⟨i, hi, hval⟩ := List.mem_iff_getElem.mp hGmem
  have hgapslen : gaps.length = cyc.tail.length := by rw [hgaps]; simp
  have htaillen : cyc.tail.length = cyc.length - 1 := by simp
  have hilt : i < cyc.tail.length := by rw [← hgapslen]; exact hi
  have hi1 : i + 1 < cyc.length := by rw [htaillen] at hilt; omega
  have hi0 : i < cyc.length := by omega
  -- gaps[i] = cyc[i+1] - cyc[i]
  have hgapval : cyc[i+1] - cyc[i] = maxCircGap E Vmax j := by
    have h1 : gaps[i]'hi = cyc[i+1] - cyc[i] := by
      simp only [hgaps, List.getElem_zipWith, List.getElem_tail]
    rw [hGdef, ← hval, h1]
  refine ⟨cyc[i], cyc[i+1], hmono i (i+1) hi0 hi1 (by omega), hcycle (i+1) hi1, ?_, ?_⟩
  · rw [hgapval]; exact hgp
  · -- freeness: no residue strictly between cyc[i] and cyc[i+1]
    intro e he ⟨hlo, hhi⟩
    set r : ℕ := e * j % Vmax with hr
    have hrres : r ∈ res := by rw [hres, hr]; exact List.mem_map.mpr ⟨e, he, rfl⟩
    have hrps : r ∈ ps := by rw [hps]; exact (List.mem_mergeSort).mpr hrres
    have hrcyc : r ∈ cyc := by rw [hcyc]; exact List.mem_append_left _ hrps
    obtain ⟨k, hk, hkval⟩ := List.mem_iff_getElem.mp hrcyc
    rcases Nat.lt_or_ge k (i+1) with hki | hki
    · -- k ≤ i ⟹ cyc[k] ≤ cyc[i] = lo, contradicting lo < r
      rcases Nat.lt_or_eq_of_le (Nat.le_of_lt_succ hki) with hki2 | hki2
      · have := hmono k i hk hi0 hki2; rw [hkval] at this; omega
      · have : cyc[k] = cyc[i] := by congr 1
        rw [hkval] at this; omega
    · -- k ≥ i+1 ⟹ cyc[k] ≥ cyc[i+1] = hi, contradicting r < hi
      rcases Nat.lt_or_eq_of_le hki with hki2 | hki2
      · have := hmono (i+1) k hi1 hk hki2; rw [hkval] at this; omega
      · have : cyc[i+1] = cyc[k] := by congr 1
        rw [hkval] at this; omega

/-- **`hlink` DISCHARGED (klein-S204 + kps-S101).**  A good period yields a free real interval
`(a, a+g)` of length `g > 1/7`, contained in `[0,1]`, that no integer translate of any tooth enters —
exactly the hypothesis `mreach_ge_of_goodPeriod` abstracts as `hlink`.  Combines the integer free-gap
extraction (this file) with kps-S101's `teeth_subset_Ico` + `free_translate_of_free_subInterval`. -/
theorem hlink_of_goodPeriod (E : List ℕ) (Vmax : ℕ) (h0 : 0 ∈ E) :
    HasGoodPeriod E Vmax →
      ∃ (j : ℕ) (a g : ℝ), 1 / 7 < g ∧
        ∀ c ∈ teeth E Vmax j, ∀ n : ℤ, (c + (n : ℝ)) ∉ Set.Ioo a (a + g) := by
  rintro ⟨j, hjmem, hgp⟩
  have hV : 0 < Vmax := by
    rw [Finset.mem_Ioo] at hjmem; omega
  have hVR : (0 : ℝ) < (Vmax : ℝ) := by exact_mod_cast hV
  obtain ⟨lo, hi, hlohi, hhiV, hgapbig, hfree⟩ := goodPeriod_intFreeGap E Vmax j hV h0 hgp
  have hcast : ((hi - lo : ℕ) : ℝ) = (hi : ℝ) - (lo : ℝ) := by
    rw [Nat.cast_sub hlohi]
  refine ⟨j, (lo : ℝ) / (Vmax : ℝ), ((hi : ℝ) - (lo : ℝ)) / (Vmax : ℝ), ?_, ?_⟩
  · -- 1/7 < g  ⟺  Vmax < 7·(hi−lo)
    have h1 : (Vmax : ℝ) < 7 * ((hi : ℝ) - (lo : ℝ)) := by
      rw [← hcast]; exact_mod_cast hgapbig
    rw [lt_div_iff₀ hVR]; linarith
  · -- the free-translate conclusion, via kps-S101's non-wrapping reduction
    have hag : (lo : ℝ) / (Vmax : ℝ) + ((hi : ℝ) - (lo : ℝ)) / (Vmax : ℝ) ≤ 1 := by
      rw [← add_div, show ((lo : ℝ) + ((hi : ℝ) - (lo : ℝ))) = (hi : ℝ) from by ring,
        div_le_one hVR]; exact_mod_cast hhiV
    apply free_translate_of_free_subInterval (teeth E Vmax j)
      ((lo : ℝ) / (Vmax : ℝ)) (((hi : ℝ) - (lo : ℝ)) / (Vmax : ℝ))
      (teeth_subset_Ico E Vmax j hV) (by positivity) hag
    -- no tooth in the open interval
    intro c hc hmem
    simp only [teeth, List.mem_toFinset, List.mem_map] at hc
    obtain ⟨e, he, rfl⟩ := hc
    obtain ⟨hlt1, hlt2⟩ := hmem
    rw [← add_div, show ((lo : ℝ) + ((hi : ℝ) - (lo : ℝ))) = (hi : ℝ) from by ring] at hlt2
    -- lo < e·j%Vmax < hi as naturals
    have hloR : (lo : ℝ) < ((e * j % Vmax : ℕ) : ℝ) := (div_lt_div_iff_of_pos_right hVR).mp hlt1
    have hhiR : ((e * j % Vmax : ℕ) : ℝ) < (hi : ℝ) := (div_lt_div_iff_of_pos_right hVR).mp hlt2
    exact hfree e he ⟨by exact_mod_cast hloR, by exact_mod_cast hhiR⟩

/-- **The good-period leg, with `hlink` discharged (klein-S204).**  `HasGoodPeriod E Vmax` (for the LRC
co-offset convention `0 ∈ E`) gives `1/14 ≤ Mreach v` modulo ONLY the ruler embedding `hembed`
(= THM-527 Part A).  The finite free-gap link is now proven; the sole remaining hypothesis is the
shared Part-A node. -/
theorem mreach_ge_of_goodPeriod_of_embed
    (E : List ℕ) (Vmax : ℕ) (v : Fin 13 → ℤ) (h0 : 0 ∈ E)
    (hgp : HasGoodPeriod E Vmax)
    (hembed : (∃ (j : ℕ) (φ : ℝ), ∀ c ∈ teeth E Vmax j,
        1 / 14 < LRC14Concrete.nearInt (φ - c)) →
      ∃ τ : ℝ, τ ∈ Set.Icc (0 : ℝ) 1 ∧ (1 : ℝ) / 14 ≤ LRC14Concrete.minReach v τ) :
    (1 : ℝ) / 14 ≤ LRC14Concrete.Mreach v :=
  mreach_ge_of_goodPeriod E Vmax v hgp (hlink_of_goodPeriod E Vmax h0) hembed

/-! ## Axiom audit -/
#print axioms goodPeriod_intFreeGap
#print axioms hlink_of_goodPeriod
#print axioms mreach_ge_of_goodPeriod_of_embed

end LRC14
