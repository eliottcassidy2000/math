/-
  TournamentH7.LRCPeelThreeCombHybrid -- arithmetic kernel for THM-1213.

  The geometric inputs are external: the eight-core component floor from
  THM-1148, the exact one-period transfer from THM-1137, and the sharp
  three-comb mass/component lemma from THM-1097.  This file checks that the
  first transfer is in its saturated regime, that the resulting three-comb
  gate is scale-free, that the complete early-candidate table lies below the
  five-piece exact envelope on `25/4 <= d/a < 49/6`, and that its
  positive-dispersion decomposition makes the remaining cone
  `49*a <= 6*d` strictly positive for `a<b<c<d`.

  Kernel-pure: no `sorry`, no `native_decide`, and no new axioms.
-/
import TournamentH7.LRCSharpCombArithmetic

namespace LonelyRunner
namespace PeelThreeCombHybrid

/-- THM-1148's normalized core floor lies strictly beyond the saturation
threshold in the exact one-period transfer. -/
theorem coreEntry_saturates : (13 : ℝ) / 7 < 72 / 35 := by
  norm_num

/-- The normalized floor remains saturated after passing to any larger legal
first-killer product. -/
theorem legalCoreEntry_saturates {x : ℝ} (hx : 72 / 35 ≤ x) :
    (13 : ℝ) / 7 < x := by
  linarith [coreEntry_saturates]

/-- Substitution of the complete first safe gap `ell=6/(7*m*a)` into the
three-comb analytic quantity from THM-1097 cancels the scale `m`. -/
theorem threeComb_after_firstGap_scale_free {m a b c d : ℝ}
    (hm : m ≠ 0) (ha : a ≠ 0) (hb : b ≠ 0) (hc : c ≠ 0) :
    21 * (6 / (7 * m * a)) * (m * d)
          - 7 * (6 / (7 * m * a)) * (m * b + m * c)
          - 6 * ((m * d) / (m * b) + (m * d) / (m * c)) - 37 =
      6 * (3 * d - b - c) / a - 6 * d / b - 6 * d / c - 37 := by
  field_simp
  ring

/-- Exact positive-dispersion identity for the hybrid gate.  The two extra
terms measure how far the interior shape coordinates are from the endpoints.
-/
theorem hybridQ_dispersion_identity {a b c d : ℝ}
    (ha : a ≠ 0) (hb : b ≠ 0) (hc : c ≠ 0) :
    6 * (3 * d - b - c) / a - 6 * d / b - 6 * d / c - 37 =
      6 * d / a - 49
        + 6 * (b - a) * (d - b) / (a * b)
        + 6 * (c - a) * (d - c) / (a * c) := by
  field_simp
  ring

/-- In the ordered positive shape cone, both dispersion energies are
strictly positive.  Thus the weak boundary `49*a <= 6*d` implies the strict
analytic gate needed by the three-comb theorem. -/
theorem hybridQ_positive_of_cone {a b c d : ℝ}
    (ha : 0 < a) (hab : a < b) (hbc : b < c) (hcd : c < d)
    (hcone : 49 * a ≤ 6 * d) :
    0 < 6 * (3 * d - b - c) / a - 6 * d / b - 6 * d / c - 37 := by
  have hb : 0 < b := lt_trans ha hab
  have hc : 0 < c := lt_trans hb hbc
  have hbd : b < d := lt_trans hbc hcd
  have hac : a < c := lt_trans hab hbc
  have hbase : 0 ≤ 6 * d / a - 49 := by
    apply sub_nonneg.mpr
    exact (le_div_iff₀ ha).2 (by nlinarith)
  have henergy_b : 0 < 6 * (b - a) * (d - b) / (a * b) := by
    positivity
  have henergy_c : 0 < 6 * (c - a) * (d - c) / (a * c) := by
    positivity
  rw [hybridQ_dispersion_identity ha.ne' hb.ne' hc.ne']
  linarith

/-- The first-transfer floor and the clean shape cone packaged as the two
arithmetic obligations used by THM-1213's geometric proof. -/
theorem peelThreeComb_arithmetic_package {x a b c d : ℝ}
    (hx : 72 / 35 ≤ x) (ha : 0 < a) (hab : a < b)
    (hbc : b < c) (hcd : c < d) (hcone : 49 * a ≤ 6 * d) :
    (13 : ℝ) / 7 < x ∧
      0 < 6 * (3 * d - b - c) / a - 6 * d / b - 6 * d / c - 37 := by
  exact ⟨legalCoreEntry_saturates hx,
    hybridQ_positive_of_cone ha hab hbc hcd hcone⟩

/-! ## Exact fractional-discrepancy strengthening

The external geometric proof defines three fractional discrepancy fees
`eb,ec,ed` and integer incidence costs.  The next lemma checks the exact
scale cancellation in the mass score; it is independent of how those fees
are produced. -/

/-- With `L=6/(7*m*a)`, the inequality `7*m*d*A>C` is scale-free even when
the three exact fractional discrepancy fees are retained separately. -/
theorem exactMassScore_scale_free {m a b c d eb ec ed : ℝ}
    (hm : m ≠ 0) (ha : a ≠ 0) (hb : b ≠ 0) (hc : c ≠ 0) (hd : d ≠ 0) :
    7 * (m * d) *
        (4 * (6 / (7 * m * a)) / 7
          - (1 / m) * (eb / b + ec / c + ed / d)) =
      24 * d / (7 * a) - 7 * d * (eb / b + ec / c + ed / d) := by
  field_simp
  ring

/-- Abstract assembly of the sawtooth envelope: if the two interior phase
costs are bounded by a champion `K`, a strict champion margin implies the
exact one-peel gate. -/
theorem exactGate_of_champion {r kx ky kr K : ℝ}
    (hx : kx ≤ K) (hy : ky ≤ K)
    (hmargin : 1 + 2 * K + kr < 24 * r / 7) :
    1 + kx + ky + kr < 24 * r / 7 := by
  linarith

/-- At a zero envelope margin, distinct interior phases still give strictness
when the champion is unique and at least one phase has smaller cost. -/
theorem twoCosts_lt_twice_champion {kx ky K : ℝ}
    (hx : kx ≤ K) (hy : ky ≤ K) (hstrict : kx < K ∨ ky < K) :
    kx + ky < 2 * K := by
  rcases hstrict with h | h <;> linarith

/-! ## Complete early-candidate audit

The sawtooth proof of the sharp `25/4` cone has nine candidate values below
the first old champion `z=6`.  The next predicate records exactly the table
in THM-1213.  Its five lemmas kernel-check that every entry is strictly below
the appropriate affine envelope piece.  Together with the later champion
comparisons, this removes any need to restrict the two interior coordinates
by the asymptotic `Q4` slope gate.
-/

def EarlyCandidatesBelow (r K : ℝ) : Prop :=
  2 + r / 7 < K ∧
  2 + 9 * r / 14 < K ∧
  3 + 6 * r / 91 < K ∧
  3 + 12 * r / 35 < K ∧
  4 + 3 * r / 70 < K ∧
  4 + 18 * r / 77 < K ∧
  5 + 2 * r / 63 < K ∧
  5 + 36 * r / 203 < K ∧
  6 + 3 * r / 119 < K

/-- Every early endpoint lies below `K_6=6+r/7` on the first envelope
piece. -/
theorem earlyCandidates_below_K6 {r : ℝ}
    (hlo : 25 / 4 ≤ r) (hi : r ≤ 41 / 6) :
    EarlyCandidatesBelow r (6 + r / 7) := by
  unfold EarlyCandidatesBelow
  repeat' apply And.intro
  all_goals linarith

/-- Every early endpoint lies below `K_41=7+6r/287` on the second envelope
piece. -/
theorem earlyCandidates_below_K41 {r : ℝ}
    (_hlo : 41 / 6 ≤ r) (hi : r ≤ 246 / 35) :
    EarlyCandidatesBelow r (7 + 6 * r / 287) := by
  unfold EarlyCandidatesBelow
  repeat' apply And.intro
  all_goals linarith

/-- Every early endpoint lies below the moving-endpoint value `36r/7-29`
on the middle envelope piece. -/
theorem earlyCandidates_below_moving {r : ℝ}
    (hlo : 246 / 35 ≤ r) (_hi : r ≤ 43 / 6) :
    EarlyCandidatesBelow r (36 * r / 7 - 29) := by
  unfold EarlyCandidatesBelow
  repeat' apply And.intro
  all_goals linarith

/-- Every early endpoint lies below `K_43=7+36r/301` on the fourth envelope
piece. -/
theorem earlyCandidates_below_K43 {r : ℝ}
    (hlo : 43 / 6 ≤ r) (hi : r ≤ 8) :
    EarlyCandidatesBelow r (7 + 36 * r / 301) := by
  unfold EarlyCandidatesBelow
  repeat' apply And.intro
  all_goals linarith

/-- Every early endpoint lies below `K_8=8+r/56` on the final exact-envelope
piece. -/
theorem earlyCandidates_below_K8 {r : ℝ}
    (_hlo : 8 ≤ r) (hi : r ≤ 49 / 6) :
    EarlyCandidatesBelow r (8 + r / 56) := by
  unfold EarlyCandidatesBelow
  repeat' apply And.intro
  all_goals linarith

/-- First affine margin in the `25/4` envelope. -/
theorem margin_25_4_to_41_6 {r : ℝ} (hr : 25 / 4 ≤ r) :
    0 ≤ 4 * r - 25 := by
  linarith

/-- Second affine margin; it is already at least one at `r=41/6`. -/
theorem margin_41_6_to_7 {r : ℝ} (hr : 41 / 6 ≤ r) :
    1 ≤ (1218 / 287) * r - 28 := by
  nlinarith

/-- Third affine margin, decreasing to `58/35` at `246/35`. -/
theorem margin_7_to_246_35 {r : ℝ} (hr : r ≤ 246 / 35) :
    58 / 35 ≤ 14 - (504 / 287) * r := by
  nlinarith

/-- Fourth affine margin, with its sole zero at the upper endpoint `43/6`. -/
theorem margin_246_35_to_43_6 {r : ℝ} (hr : r ≤ 43 / 6) :
    0 ≤ 86 - 12 * r := by
  linarith

/-- Fifth affine margin, with its sole zero at the lower endpoint `43/6`. -/
theorem margin_43_6_to_8 {r : ℝ} (hr : 43 / 6 ≤ r) :
    0 ≤ (1218 / 301) * r - 29 := by
  nlinarith

/-- Sixth affine margin; the incidence jump at `r=8` still leaves margin two.
-/
theorem margin_8_to_49_6 {r : ℝ} (hr : 8 ≤ r) :
    2 ≤ (119 / 28) * r - 32 := by
  nlinarith

/-! ## First-gap address repair

For the phase-blind sharp family, `b = 6*a`.  A complete `a`-safe gap has
endpoints `(14*j+1)/(14*a)` and `(14*j+13)/(14*a)`; the next lemma checks
their exact images in the `b`-coordinate.  The interval-occupancy statement
(exactly five danger teeth) remains geometric, while its resulting strict
score is kernel-checked below.
-/

theorem sixfold_firstGap_endpoints {a j : ℝ} (ha : a ≠ 0) :
    (6 * a) * ((14 * j + 1) / (14 * a)) = 6 * j + 3 / 7 ∧
    (6 * a) * ((14 * j + 13) / (14 * a)) = 6 * j + 39 / 7 := by
  constructor <;> field_simp <;> ring

/-- Strict positivity of the repaired score for
`(a,b,c,d)=(4n,24n,24n+1,25n-1)`, `n≥3`. -/
theorem addressRepair_margin_pos {n : ℝ} (hn : 3 ≤ n) :
    0 < ((15 * n + 1) * (40 * n - 31)) / (24 * n * (24 * n + 1)) := by
  have hn0 : 0 < n := by linarith
  have hleft : 0 < 15 * n + 1 := by nlinarith
  have hright : 0 < 40 * n - 31 := by nlinarith
  have hdenLeft : 0 < 24 * n := by positivity
  have hdenRight : 0 < 24 * n + 1 := by nlinarith
  exact div_pos (mul_pos hleft hright) (mul_pos hdenLeft hdenRight)

end PeelThreeCombHybrid
end LonelyRunner

#print axioms LonelyRunner.PeelThreeCombHybrid.threeComb_after_firstGap_scale_free
#print axioms LonelyRunner.PeelThreeCombHybrid.hybridQ_dispersion_identity
#print axioms LonelyRunner.PeelThreeCombHybrid.hybridQ_positive_of_cone
#print axioms LonelyRunner.PeelThreeCombHybrid.peelThreeComb_arithmetic_package
#print axioms LonelyRunner.PeelThreeCombHybrid.exactMassScore_scale_free
#print axioms LonelyRunner.PeelThreeCombHybrid.exactGate_of_champion
#print axioms LonelyRunner.PeelThreeCombHybrid.twoCosts_lt_twice_champion
#print axioms LonelyRunner.PeelThreeCombHybrid.earlyCandidates_below_K6
#print axioms LonelyRunner.PeelThreeCombHybrid.earlyCandidates_below_K41
#print axioms LonelyRunner.PeelThreeCombHybrid.earlyCandidates_below_moving
#print axioms LonelyRunner.PeelThreeCombHybrid.earlyCandidates_below_K43
#print axioms LonelyRunner.PeelThreeCombHybrid.earlyCandidates_below_K8
#print axioms LonelyRunner.PeelThreeCombHybrid.margin_25_4_to_41_6
#print axioms LonelyRunner.PeelThreeCombHybrid.margin_41_6_to_7
#print axioms LonelyRunner.PeelThreeCombHybrid.margin_7_to_246_35
#print axioms LonelyRunner.PeelThreeCombHybrid.margin_246_35_to_43_6
#print axioms LonelyRunner.PeelThreeCombHybrid.margin_43_6_to_8
#print axioms LonelyRunner.PeelThreeCombHybrid.margin_8_to_49_6
#print axioms LonelyRunner.PeelThreeCombHybrid.sixfold_firstGap_endpoints
#print axioms LonelyRunner.PeelThreeCombHybrid.addressRepair_margin_pos
