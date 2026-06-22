/-
  TournamentH7.LRCLowerThresholdNeighborhood -- exact arithmetic checksum for
  the HYP-2847 lower-threshold Farey-neighborhood scout.

  The exact engine is:

    04-computation/lrc_lower_threshold_nbhd_width_codex_s88.py

  with stored output:

    05-knowledge/results/lrc_lower_threshold_nbhd_width_codex_s88.out

  It tests the HYP-2842 redirect for even thresholds N=2q<=14 in the dense
  bounded branch.  Proper Farey neighborhoods alone have zero certified
  residual in the deepest bounded branch for every checked N.  Adding the
  b=1 origin-collapse neighborhood repairs the certificate:

      N= 6 : 1/15
      N= 8 : 5/168
      N=10 : 1/60
      N=12 : 5/264
      N=14 : 2/273

  This file records only those rational readouts and comparisons.  It is not
  the interval-containment theorem.  The next Lean atom should prove that the
  Farey/origin neighborhoods produced by the script are subsets of
  `GOOD(E) ∩ G_P`, and that the exact interval subtraction ledger lower-bounds
  `slowμ (GOOD(E) ∩ G_P)`.
-/

namespace LonelyRunner
namespace LowerThresholdNeighborhood

/-! ## Proper-neighborhood failure readouts -/

def properBoundedN6MinNum : Nat := 0
def properBoundedN8MinNum : Nat := 0
def properBoundedN10MinNum : Nat := 0
def properBoundedN12MinNum : Nat := 0
def properBoundedN14MinNum : Nat := 0

/-! ## Origin+proper bounded dense residual readouts -/

def originBoundedN6MinNum : Nat := 1
def originBoundedN6MinDen : Nat := 15

def originBoundedN8MinNum : Nat := 5
def originBoundedN8MinDen : Nat := 168

def originBoundedN10MinNum : Nat := 1
def originBoundedN10MinDen : Nat := 60

def originBoundedN12MinNum : Nat := 5
def originBoundedN12MinDen : Nat := 264

def originBoundedN14MinNum : Nat := 2
def originBoundedN14MinDen : Nat := 273

/-- Proper Farey neighborhoods alone vanish in the deepest bounded dense branch
for every checked lower threshold N=6,8,10,12,14. -/
theorem proper_bounded_deep_branch_zero_all :
    properBoundedN6MinNum = 0 ∧
      properBoundedN8MinNum = 0 ∧
      properBoundedN10MinNum = 0 ∧
      properBoundedN12MinNum = 0 ∧
      properBoundedN14MinNum = 0 := by
  native_decide

/-- Adding the b=1 origin-collapse neighborhood gives a positive rational
residual in every checked bounded dense branch N=6,8,10,12,14. -/
theorem origin_bounded_dense_positive_all :
    0 < originBoundedN6MinNum ∧
      0 < originBoundedN8MinNum ∧
      0 < originBoundedN10MinNum ∧
      0 < originBoundedN12MinNum ∧
      0 < originBoundedN14MinNum := by
  native_decide

/-- The N=14 residual `2/273` is the smallest of the reported
origin+proper bounded dense floors. -/
theorem origin_n14_is_smallest_reported_floor :
    originBoundedN14MinNum * originBoundedN6MinDen ≤
        originBoundedN6MinNum * originBoundedN14MinDen ∧
      originBoundedN14MinNum * originBoundedN8MinDen ≤
        originBoundedN8MinNum * originBoundedN14MinDen ∧
      originBoundedN14MinNum * originBoundedN10MinDen ≤
        originBoundedN10MinNum * originBoundedN14MinDen ∧
      originBoundedN14MinNum * originBoundedN12MinDen ≤
        originBoundedN12MinNum * originBoundedN14MinDen := by
  native_decide

/-- Packaged checksum for the HYP-2847 exact scout. -/
theorem audited_lower_threshold_neighborhood_readouts :
    properBoundedN6MinNum = 0 ∧
      properBoundedN8MinNum = 0 ∧
      properBoundedN10MinNum = 0 ∧
      properBoundedN12MinNum = 0 ∧
      properBoundedN14MinNum = 0 ∧
      0 < originBoundedN6MinNum ∧
      0 < originBoundedN8MinNum ∧
      0 < originBoundedN10MinNum ∧
      0 < originBoundedN12MinNum ∧
      0 < originBoundedN14MinNum ∧
      originBoundedN14MinNum * originBoundedN6MinDen ≤
        originBoundedN6MinNum * originBoundedN14MinDen ∧
      originBoundedN14MinNum * originBoundedN8MinDen ≤
        originBoundedN8MinNum * originBoundedN14MinDen ∧
      originBoundedN14MinNum * originBoundedN10MinDen ≤
        originBoundedN10MinNum * originBoundedN14MinDen ∧
      originBoundedN14MinNum * originBoundedN12MinDen ≤
        originBoundedN12MinNum * originBoundedN14MinDen := by
  native_decide

/-! ### Axiom audit -/

#print axioms proper_bounded_deep_branch_zero_all
#print axioms origin_bounded_dense_positive_all
#print axioms origin_n14_is_smallest_reported_floor
#print axioms audited_lower_threshold_neighborhood_readouts

end LowerThresholdNeighborhood
end LonelyRunner
