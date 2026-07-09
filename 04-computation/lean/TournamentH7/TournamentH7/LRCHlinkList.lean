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

end LRC14HlinkList
