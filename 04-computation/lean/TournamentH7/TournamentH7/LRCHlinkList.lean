/-
  TournamentH7.LRCHlinkList — the two list lemmas that finish klein-S203's `hlink`
  (kind-pasteur-2026-07-09-S103).

  `maxCircGap = (zipWith (·-·) cyc cyc.tail).foldl max 0`.  After `foldl_max_pos_mem`
  (kps-S102) the widest gap is a member of the zipWith list.  These two lemmas turn that into a
  free residue interval:

    * `mem_zipWith_sub_adjacency` — a member of `zipWith (b−a) L L.tail` is a difference `q − p`
      of ADJACENT list entries `L = l₁ ++ p :: q :: l₂` (positional zipWith → adjacency);
    * `sorted_adjacency_sep` — in a SORTED list `l₁ ++ p :: q :: l₂`, every entry `r` satisfies
      `r ≤ p ∨ q ≤ r` (nothing lands strictly between adjacent sorted entries).

  Together: the widest gap `q − p` has NO residue in the open interval `(p, q)`, which normalized by
  `Vmax` is the free arc `> 1/7` that `hlink` needs (dispatched to `LRCTeethGap`/`LRCMaxgapArgmax`'s
  translate reductions).  Self-contained (Mathlib only).
-/
import Mathlib

namespace LRC14HlinkList

/-- **zipWith consecutive-difference → adjacency.**  If `x` is a member of
`zipWith (fun a b => b − a) L L.tail`, then `x = q − p` for a pair `p, q` ADJACENT in `L`
(`L = l₁ ++ p :: q :: l₂`). -/
theorem mem_zipWith_sub_adjacency : ∀ (L : List ℕ) (x : ℕ),
    x ∈ List.zipWith (fun a b => b - a) L L.tail →
    ∃ (l₁ l₂ : List ℕ) (p q : ℕ), L = l₁ ++ p :: q :: l₂ ∧ q - p = x
  | [], x, hx => by simp at hx
  | [_], x, hx => by simp at hx
  | a :: b :: t, x, hx => by
    simp only [List.tail_cons, List.zipWith_cons_cons, List.mem_cons] at hx
    rcases hx with h | h
    · exact ⟨[], t, a, b, rfl, h.symm⟩
    · have h' : x ∈ List.zipWith (fun a b => b - a) (b :: t) (b :: t).tail := by
        simpa using h
      obtain ⟨l₁, l₂, p, q, heq, hpq⟩ := mem_zipWith_sub_adjacency (b :: t) x h'
      refine ⟨a :: l₁, l₂, p, q, ?_, hpq⟩
      rw [List.cons_append, heq]

/-- **Sorted adjacency separation.**  In a sorted list `l₁ ++ p :: q :: l₂`, every element `r`
satisfies `r ≤ p ∨ q ≤ r`: nothing lies strictly between the adjacent entries `p, q`. -/
theorem sorted_adjacency_sep {l₁ l₂ : List ℕ} {p q : ℕ}
    (hL : List.Pairwise (· ≤ ·) (l₁ ++ p :: q :: l₂)) :
    ∀ r ∈ l₁ ++ p :: q :: l₂, r ≤ p ∨ q ≤ r := by
  intro r hr
  obtain ⟨-, hp2, hcross⟩ := List.pairwise_append.mp hL
  rw [List.mem_append] at hr
  rcases hr with hr | hr
  · exact Or.inl (hcross r hr p (List.mem_cons_self ..))
  · rw [List.mem_cons] at hr
    rcases hr with rfl | hr
    · exact Or.inl (le_refl _)
    · rw [List.mem_cons] at hr
      rcases hr with rfl | hr
      · exact Or.inr (le_refl _)
      · rw [List.pairwise_cons] at hp2
        rw [List.pairwise_cons] at hp2
        exact Or.inr (hp2.2.1 r hr)

/-- **The cyclic list is sorted.**  Appending a point `x` no smaller than every entry to a sorted
list keeps it sorted.  Applied with `l = ps` (sorted residues) and `x = p0 + Vmax` (all residues
`< Vmax ≤ p0 + Vmax`), this is the sortedness of `cyc = ps ++ [p0 + Vmax]` needed for
`sorted_adjacency_sep`. -/
theorem pairwise_append_singleton_of_le (l : List ℕ) (x : ℕ)
    (hl : List.Pairwise (· ≤ ·) l) (hx : ∀ a ∈ l, a ≤ x) :
    List.Pairwise (· ≤ ·) (l ++ [x]) := by
  rw [List.pairwise_append]
  refine ⟨hl, List.pairwise_singleton _ _, ?_⟩
  intro a ha b hb
  rw [List.mem_singleton] at hb
  subst hb
  exact hx a ha

/-- **The unified translate reduction (no internal/wrap dispatch).**  A tooth `r/Vmax` (residue
`r < Vmax`), with the sorted-separation `r ≤ p ∨ q ≤ r` and the circular bound `q ≤ r + Vmax`,
has NO integer translate in the open gap `(p/Vmax, q/Vmax)`.  This covers BOTH the internal gap
(`q < Vmax`) and the wraparound gap (`q = p0 + Vmax`) in one shot: the bound `q ≤ r + Vmax` (from
`q ≤ p0 + Vmax` and `p0 ≤ r`, the minimum residue) is what makes the `n = ±1` translates miss.
This is the single reduction that finishes `hlink`, replacing the S101/S102 case split. -/
theorem tooth_not_in_gap (Vmax p q r : ℕ) (n : ℤ)
    (hV : 0 < Vmax) (hsep : r ≤ p ∨ q ≤ r) (hqrV : q ≤ r + Vmax) (hrV : r < Vmax) :
    (r : ℝ) / Vmax + (n : ℝ) ∉ Set.Ioo ((p : ℝ) / Vmax) ((q : ℝ) / Vmax) := by
  rintro ⟨hlo, hhi⟩
  have hVR : (0 : ℝ) < (Vmax : ℝ) := by exact_mod_cast hV
  have hVne : (Vmax : ℝ) ≠ 0 := ne_of_gt hVR
  -- clear the denominator: p < r + n·Vmax < q  (over ℝ)
  have h1 : (p : ℝ) < (r : ℝ) + (n : ℝ) * Vmax := by
    have h := mul_lt_mul_of_pos_right hlo hVR
    rwa [div_mul_cancel₀ _ hVne, add_mul, div_mul_cancel₀ _ hVne] at h
  have h2 : (r : ℝ) + (n : ℝ) * Vmax < (q : ℝ) := by
    have h := mul_lt_mul_of_pos_right hhi hVR
    rwa [add_mul, div_mul_cancel₀ _ hVne, div_mul_cancel₀ _ hVne] at h
  have hrR : (r : ℝ) < (Vmax : ℝ) := by exact_mod_cast hrV
  have hqR : (q : ℝ) ≤ (r : ℝ) + (Vmax : ℝ) := by exact_mod_cast hqrV
  rcases hsep with hs | hs
  · have hsR : (r : ℝ) ≤ (p : ℝ) := by exact_mod_cast hs
    by_cases hn : (n : ℝ) ≤ 0
    · have : (n : ℝ) * Vmax ≤ 0 := mul_nonpos_of_nonpos_of_nonneg hn (le_of_lt hVR)
      linarith
    · have hn' : 0 < (n : ℝ) := not_le.mp hn
      have hn1 : (1 : ℝ) ≤ (n : ℝ) := by
        have hz : 0 < n := by exact_mod_cast hn'
        have : (1 : ℤ) ≤ n := hz
        exact_mod_cast this
      have : (Vmax : ℝ) ≤ (n : ℝ) * Vmax := le_mul_of_one_le_left (le_of_lt hVR) hn1
      linarith
  · have hsR : (q : ℝ) ≤ (r : ℝ) := by exact_mod_cast hs
    by_cases hn : 0 ≤ (n : ℝ)
    · have : (0 : ℝ) ≤ (n : ℝ) * Vmax := mul_nonneg hn (le_of_lt hVR)
      linarith
    · have hn' : (n : ℝ) < 0 := not_le.mp hn
      have hn1 : (n : ℝ) ≤ -1 := by
        have hz : n < 0 := by exact_mod_cast hn'
        have : n ≤ -1 := by omega
        exact_mod_cast this
      have : (n : ℝ) * Vmax ≤ (-1 : ℝ) * Vmax := mul_le_mul_of_nonneg_right hn1 (le_of_lt hVR)
      linarith

end LRC14HlinkList
