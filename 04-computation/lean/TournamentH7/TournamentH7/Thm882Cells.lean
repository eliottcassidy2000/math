/- Thm882Cells.lean — klein-2026-07-16-S318 (boxeph manifest item 1).
   THM-882 (flat = twice good), the per-cell certificates: the 12 Farey-12
   adjacencies (L = a/i, R = b/j) with i + j = 13 carry the entire flat set;
   on each cell the good window is [L + 1/(14i), R − 1/(14j)] of length
   1/(14ij), the flat interval [flo, fhi] has length 1/(14·min(i,j)·|i−j|),
   and cl(G) ⊆ F ⊆ [L, R].  Data: lrc14_flat_2x_law_boxeph_S23.py Part 2
   (all rows re-extracted and re-verified exactly, klein-S318).
   The global masses: m(F) = 6617/97020 = 2 · m(G).  All pure ℚ.
   No sorries, no native_decide. -/
import Mathlib.Tactic

namespace LonelyRunner
namespace LRC14
namespace Thm882

/-- One Farey cell's certificate: adjacency gap, good-window endpoints and
    length, flat-interval length, and the containment chain
    `L ≤ flo ≤ gl ≤ gh ≤ fhi ≤ R`. -/
def cellOK (L R gl gh flo fhi flen : ℚ) (i j : ℚ) : Prop :=
  R - L = 1 / (i * j) ∧
  gl = L + 1 / (14 * i) ∧
  gh = R - 1 / (14 * j) ∧
  gh - gl = 1 / (14 * i * j) ∧
  fhi - flo = flen ∧
  L ≤ flo ∧ flo ≤ gl ∧ gl ≤ gh ∧ gh ≤ fhi ∧ fhi ≤ R

theorem cell01 : cellOK 0 (1/12) (1/14) (13/168) (1/14) (6/77) (1/154) 1 12 := by
  unfold cellOK; norm_num

theorem cell02 : cellOK (1/7) (1/6) (15/98) (13/84) (1/7) (13/84) (1/84) 7 6 := by
  unfold cellOK; norm_num

theorem cell03 : cellOK (2/9) (1/4) (29/126) (13/56) (8/35) (13/56) (1/280) 9 4 := by
  unfold cellOK; norm_num

theorem cell04 : cellOK (3/10) (1/3) (43/140) (13/42) (15/49) (13/42) (1/294) 10 3 := by
  unfold cellOK; norm_num

theorem cell05 : cellOK (3/8) (2/5) (43/112) (27/70) (8/21) (27/70) (1/210) 8 5 := by
  unfold cellOK; norm_num

theorem cell06 : cellOK (5/11) (1/2) (71/154) (13/28) (29/63) (13/28) (1/252) 11 2 := by
  unfold cellOK; norm_num

theorem cell07 : cellOK (1/2) (6/11) (15/28) (83/154) (15/28) (34/63) (1/252) 2 11 := by
  unfold cellOK; norm_num

theorem cell08 : cellOK (3/5) (5/8) (43/70) (69/112) (43/70) (13/21) (1/210) 5 8 := by
  unfold cellOK; norm_num

theorem cell09 : cellOK (2/3) (7/10) (29/42) (97/140) (29/42) (34/49) (1/294) 3 10 := by
  unfold cellOK; norm_num

theorem cell10 : cellOK (3/4) (7/9) (43/56) (97/126) (43/56) (27/35) (1/280) 4 9 := by
  unfold cellOK; norm_num

theorem cell11 : cellOK (5/6) (6/7) (71/84) (83/98) (71/84) (6/7) (1/84) 6 7 := by
  unfold cellOK; norm_num

theorem cell12 : cellOK (11/12) 1 (155/168) (13/14) (71/77) (13/14) (1/154) 12 1 := by
  unfold cellOK; norm_num

/-- **THM-882, the twelve cells together** (the i + j = 13 Farey-12 sites carry
    the flat set; the mirror pairing 1↔12, …, 6↔7 is visible in the lengths). -/
theorem thm882_cells :
    cellOK 0 (1/12) (1/14) (13/168) (1/14) (6/77) (1/154) 1 12 ∧
    cellOK (1/7) (1/6) (15/98) (13/84) (1/7) (13/84) (1/84) 7 6 ∧
    cellOK (2/9) (1/4) (29/126) (13/56) (8/35) (13/56) (1/280) 9 4 ∧
    cellOK (3/10) (1/3) (43/140) (13/42) (15/49) (13/42) (1/294) 10 3 ∧
    cellOK (3/8) (2/5) (43/112) (27/70) (8/21) (27/70) (1/210) 8 5 ∧
    cellOK (5/11) (1/2) (71/154) (13/28) (29/63) (13/28) (1/252) 11 2 ∧
    cellOK (1/2) (6/11) (15/28) (83/154) (15/28) (34/63) (1/252) 2 11 ∧
    cellOK (3/5) (5/8) (43/70) (69/112) (43/70) (13/21) (1/210) 5 8 ∧
    cellOK (2/3) (7/10) (29/42) (97/140) (29/42) (34/49) (1/294) 3 10 ∧
    cellOK (3/4) (7/9) (43/56) (97/126) (43/56) (27/35) (1/280) 4 9 ∧
    cellOK (5/6) (6/7) (71/84) (83/98) (71/84) (6/7) (1/84) 6 7 ∧
    cellOK (11/12) 1 (155/168) (13/14) (71/77) (13/14) (1/154) 12 1 :=
  ⟨cell01, cell02, cell03, cell04, cell05, cell06, cell07, cell08, cell09,
   cell10, cell11, cell12⟩

/-- **THM-882's global masses**: the twelve flat lengths sum to
    m(F) = 6617/97020, the twelve good lengths to m(G) = 6617/194040,
    and m(F) = 2·m(G) — the flat-equals-twice-good law's arithmetic face. -/
theorem thm882_masses :
    (1/154 + 1/84 + 1/280 + 1/294 + 1/210 + 1/252
      + 1/252 + 1/210 + 1/294 + 1/280 + 1/84 + 1/154 : ℚ) = 6617/97020 ∧
    (1/(14*1*12) + 1/(14*7*6) + 1/(14*9*4) + 1/(14*10*3) + 1/(14*8*5)
      + 1/(14*11*2) + 1/(14*2*11) + 1/(14*5*8) + 1/(14*3*10) + 1/(14*4*9)
      + 1/(14*6*7) + 1/(14*12*1) : ℚ) = 6617/194040 ∧
    (6617/97020 : ℚ) = 2 * (6617/194040) := by
  refine ⟨?_, ?_, ?_⟩ <;> norm_num

end Thm882
end LRC14
end LonelyRunner
