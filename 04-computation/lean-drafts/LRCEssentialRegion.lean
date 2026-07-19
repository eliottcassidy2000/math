/- LRCEssentialRegion.lean -- opus-2026-07-17-S378.

   THE ESSENTIAL-REGION CRITERION (THM-1120), proved.

   S377 found the criterion empirically: replacing speed `i` by `r` in a family
   preserves tightness exactly when the ESSENTIAL REGION of `i` -- the part of the
   interval covered by no OTHER speed -- is contained in the new speed's arcs.
   The criterion predicted the swappable speed of {1,...,13} to be exactly {12},
   matching an independent substitution search with no exceptions.

   The proof is a set identity: removing `i` and adjoining `r` leaves uncovered
   exactly  (essential region of i) \ (arcs of r).  De Morgan, i.e. Set.diff_diff.
   The measure and emptiness forms follow at once.

   The final lemma justifies the INTERVAL-CONTAINMENT test used by the S377 search:
   consecutive arcs of a single modulus are strictly separated (gap (1-2*lam)/w),
   so a connected piece of the essential region cannot straddle two arcs -- it must
   sit inside a single one, which is exactly what the script checked. -/
import Mathlib.MeasureTheory.Measure.Lebesgue.Basic
import TournamentH7.FragmentationLemma

open MeasureTheory

namespace LRC14

/-- The **essential region** of speed `i` in the family `S`, relative to the
    ambient interval `I`: the part of `I` that NO other speed of `S` covers.
    Only `i` (or its replacement) can cover this set. -/
def essentialRegion (I : Set ℝ) (S : Finset ℕ) (i : ℕ) (lam : ℝ) : Set ℝ :=
  I \ ⋃ j ∈ S.erase i, badArcs j lam

/-- **THE ESSENTIAL-REGION CRITERION, set form.**  After deleting `i` from the
    family and adjoining `r`, the uncovered part of `I` is exactly the essential
    region of `i` minus the arcs of `r`.  Everything else follows from this. -/
theorem uncovered_swap_eq (I : Set ℝ) (S : Finset ℕ) (i r : ℕ) (lam : ℝ) :
    I \ ((⋃ j ∈ S.erase i, badArcs j lam) ∪ badArcs r lam)
      = essentialRegion I S i lam \ badArcs r lam := by
  rw [essentialRegion, Set.diff_diff]

/-- **Tightness form.**  The swapped family is tight (uncovered measure zero)
    exactly when the essential region of `i` is covered by `r` up to a null set. -/
theorem tight_swap_iff (I : Set ℝ) (S : Finset ℕ) (i r : ℕ) (lam : ℝ) :
    volume (I \ ((⋃ j ∈ S.erase i, badArcs j lam) ∪ badArcs r lam)) = 0
      ↔ volume (essentialRegion I S i lam \ badArcs r lam) = 0 := by
  rw [uncovered_swap_eq]

/-- **Containment form** -- the direction the S377 search actually tests.
    If the essential region of `i` sits inside the arcs of `r`, the swapped
    family covers `I` completely. -/
theorem swap_covers_of_subset (I : Set ℝ) (S : Finset ℕ) (i r : ℕ) (lam : ℝ)
    (h : essentialRegion I S i lam ⊆ badArcs r lam) :
    I \ ((⋃ j ∈ S.erase i, badArcs j lam) ∪ badArcs r lam) = ∅ := by
  rw [uncovered_swap_eq]
  exact Set.diff_eq_empty.mpr h

/-- Containment of the essential region is also NECESSARY for full covering:
    the criterion is an iff at the level of sets, not merely sufficient. -/
theorem subset_of_swap_covers (I : Set ℝ) (S : Finset ℕ) (i r : ℕ) (lam : ℝ)
    (h : I \ ((⋃ j ∈ S.erase i, badArcs j lam) ∪ badArcs r lam) = ∅) :
    essentialRegion I S i lam ⊆ badArcs r lam := by
  rw [uncovered_swap_eq] at h
  exact Set.diff_eq_empty.mp h

/-- **Separation of consecutive arcs.**  For `0 < w` and `2*lam < 1`, the arc
    around `a/w` ends strictly before the arc around `(a+1)/w` begins; the gap
    has length `(1 - 2*lam)/w > 0`.

    This is what licenses the interval-containment test: since distinct arcs of a
    single modulus are strictly separated, a connected component of the essential
    region that meets no gap must lie inside ONE arc.  So checking "each component
    of E_i sits inside some component of D_r" is not merely sufficient for the
    criterion -- for connected pieces it is the only way containment can hold. -/
theorem badArcs_consecutive_separated (w : ℕ) (hw : 0 < w) (lam : ℝ)
    (h2 : 2 * lam < 1) (a : ℤ) :
    (a : ℝ) / (w : ℝ) + lam / (w : ℝ) < ((a : ℝ) + 1) / (w : ℝ) - lam / (w : ℝ) := by
  have hw' : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  rw [← add_div, ← sub_div]
  apply div_lt_div_of_pos_right _ hw'
  linarith

/-- The gap between consecutive arcs has length exactly `(1 - 2*lam)/w`. -/
theorem badArcs_gap_length (w : ℕ) (lam : ℝ) (a : ℤ) :
    (((a : ℝ) + 1) / (w : ℝ) - lam / (w : ℝ)) - ((a : ℝ) / (w : ℝ) + lam / (w : ℝ))
      = (1 - 2 * lam) / (w : ℝ) := by
  field_simp
  ring

/-! ### The swap bound

    The separation lemma has a quantitative consequence that turns the S377
    substitution search into a proof.  A connected piece of the essential region
    contained in `badArcs r lam` must lie inside a SINGLE arc, and arcs have
    length `2*lam/r`; so the longest component of `E_i` bounds the admissible
    replacement:  `r ≤ 2*lam / ellMax(E_i)`.  Beyond that bound no `r` can work,
    so checking up to it is exhaustive. -/

/-- The right endpoint of arc `a` lies in no arc at all: arcs are strictly
    separated, so the gap point escapes the whole union. -/
theorem endpoint_not_mem_badArcs (w : ℕ) (hw : 0 < w) (lam : ℝ) (h2 : 2 * lam < 1)
    (a : ℤ) :
    (a : ℝ) / (w : ℝ) + lam / (w : ℝ) ∉ badArcs w lam := by
  have hw' : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  simp only [badArcs, Set.mem_iUnion, Set.mem_Ioo, not_exists, not_and, not_lt]
  intro b hb
  by_cases hba : b ≤ a
  · -- arcs at or to the left: their right end is at or before ours
    have hcast : ((b : ℤ) : ℝ) ≤ ((a : ℤ) : ℝ) := by exact_mod_cast hba
    have hd : ((b : ℤ) : ℝ) / (w : ℝ) ≤ ((a : ℤ) : ℝ) / (w : ℝ) := by gcongr
    linarith
  · -- arcs to the right: separation puts our point before their left end
    exfalso
    have hba' : a + 1 ≤ b := by omega
    have hab : ((a : ℤ) : ℝ) + 1 ≤ ((b : ℤ) : ℝ) := by exact_mod_cast hba'
    have hstep : (((a : ℤ) : ℝ) + 1) / (w : ℝ) ≤ ((b : ℤ) : ℝ) / (w : ℝ) := by gcongr
    have hsep := badArcs_consecutive_separated w hw lam h2 a
    linarith

/-- **The swap bound.**  If a closed interval `[x,y]` is contained in
    `badArcs r lam`, then its length is at most the arc length `2*lam/r`.
    Consequently the longest component of an essential region bounds the
    admissible replacement speed. -/
theorem Icc_length_le_of_subset_badArcs (r : ℕ) (hr : 0 < r) (lam : ℝ)
    (h2 : 2 * lam < 1) (x y : ℝ) (hxy : x ≤ y)
    (h : Set.Icc x y ⊆ badArcs r lam) :
    y - x ≤ 2 * lam / (r : ℝ) := by
  by_contra hcon
  push_neg at hcon
  -- `2*lam/r` and `lam/r` must be related as atoms before linarith can use both
  have hsplit : 2 * lam / (r : ℝ) = lam / (r : ℝ) + lam / (r : ℝ) := by ring
  rw [hsplit] at hcon
  -- x lies in some arc; that arc's right endpoint lies in [x,y] yet not in badArcs
  have hx : x ∈ badArcs r lam := h (Set.left_mem_Icc.mpr hxy)
  simp only [badArcs, Set.mem_iUnion, Set.mem_Ioo] at hx
  obtain ⟨a, hxl, hxr⟩ := hx
  have hpx : x ≤ ((a : ℤ) : ℝ) / (r : ℝ) + lam / (r : ℝ) := le_of_lt hxr
  have hpy : ((a : ℤ) : ℝ) / (r : ℝ) + lam / (r : ℝ) ≤ y := by linarith
  exact endpoint_not_mem_badArcs r hr lam h2 a (h ⟨hpx, hpy⟩)

end LRC14
