import TournamentH7.LRCNestedCarrierWindow

/-!
# Five-killer carrier-window arithmetic (THM-1214)

This module isolates the reusable arithmetic boundary of the clustered
five-killer proof.  It checks:

* the twelve-chart owner capacity for every carrier cardinality `rho ≤ 5`;
* the endpoint contradiction which eliminates the residual two-carrier
  one-tooth cover;
* the full- and truncated-window density margins, including the five exact
  first-good thresholds in the three-carrier dispatch;
* the uniform four-carrier density margin and its concrete set-level cover
  exclusion; and
* the first ratio handoff and nested-window contradiction for five carriers.

The geometric inputs remaining outside this file are extraction of a
carrier-safe boundary from a compact uncovered set, the residue-owner chart
lemma itself, and the finite open wall-cell reduction.  The four-carrier row
below now instantiates the multi-comb measure bound on the actual closed first
window, so that this important arithmetic-to-set bridge is kernel checked.
Thus no finite search is hidden behind the arithmetic API.
-/

namespace LRC14.FiveKillerCarrierWindow

/-! ## Owner capacity -/

/-- With five killers and `rho ≤ 5` carriers, the chart ledger
`d + 2(5-rho) < 12` is exactly the owner bound `d < 2rho+2`. -/
theorem owner_budget_iff {rho d : ℕ} (hrho : rho ≤ 5) :
    d + 2 * (5 - rho) < 12 ↔ d < 2 * rho + 2 := by
  omega

/-- The convenient integral cap `d ≤ 2rho+1` always leaves a chart. -/
theorem owner_cap_survives {rho d : ℕ} (hrho : rho ≤ 5)
    (hd : d ≤ 2 * rho + 1) :
    d + 2 * (5 - rho) < 12 := by
  apply (owner_budget_iff hrho).2
  omega

/-- The capped owner count used by THM-1214 also covers the `rho=4,5` rows,
where all eight core owners may be crossed. -/
theorem capped_owner_survives {rho d : ℕ} (hrho : rho ≤ 5)
    (hd : d ≤ min 8 (2 * rho + 1)) :
    d + 2 * (5 - rho) < 12 := by
  apply owner_cap_survives hrho
  exact le_trans hd (min_le_right _ _)

/-! ## Two carriers: endpoint arithmetic -/

/-- Eliminating the putative tooth speed `y` from the two strict endpoint
inequalities gives THM-1214 equation (20). -/
theorem rho2_endpoint_elimination {p x y n : ℝ}
    (hleft : x * (14 * n - 1) < y)
    (hright : y < p * (14 * n + 1)) :
    14 * n * (x - p) < x + p := by
  nlinarith

/-- Already from `p ≤ 7` and `x ≥ p+6`, no positive-centre tooth can satisfy
both strict containment inequalities.  This is slightly stronger than the
`x < 3p` residual formulation used in the paper proof.  The proof is valid
over the reals and hence applies in particular to the integer tooth centre. -/
theorem rho2_endpoint_contradiction {p x y n : ℝ}
    (hp : p ≤ 7) (hgap : p + 6 ≤ x) (hn : 1 ≤ n)
    (hleft : x * (14 * n - 1) < y)
    (hright : y < p * (14 * n + 1)) : False := by
  have helim := rho2_endpoint_elimination hleft hright
  have hn' : 0 ≤ n - 1 := by linarith
  have hgap' : 0 ≤ x - p - 6 := by linarith
  have hproduct : 0 ≤ (n - 1) * (x - p - 6) :=
    mul_nonneg hn' hgap'
  nlinarith

/-! ## Three carriers: density margins -/

/-- Exact full-window identity behind THM-1214 equation (24). -/
theorem rho3_full_window_identity {x : ℝ} (hx : 0 < x) :
    6 / (7 * x) - 2 / (5 * (x + 1)) =
      (16 * x + 30) / (35 * x * (x + 1)) := by
  have hx0 : x ≠ 0 := ne_of_gt hx
  have hx10 : x + 1 ≠ 0 := by linarith
  field_simp
  ring

/-- The full least-carrier window is too long for the two-comb cover budget. -/
theorem rho3_full_window_margin {x : ℝ} (hx : 0 < x) :
    2 / (5 * (x + 1)) < 6 / (7 * x) := by
  have hid := rho3_full_window_identity hx
  have hden : 0 < 35 * x * (x + 1) := by positivity
  have hnum : 0 < 16 * x + 30 := by positivity
  have hpos : 0 < (16 * x + 30) / (35 * x * (x + 1)) :=
    div_pos hnum hden
  linarith

/-- Exact numerator of the truncated owner-window density margin. -/
theorem rho3_truncated_window_identity {p x : ℝ}
    (hp : 0 < p) (hx : 0 < x) :
    5 * ((x - p) / (14 * x * p)) - 2 / (x + 1) =
      (5 * (x - p) * (x + 1) - 28 * x * p) /
        (14 * x * p * (x + 1)) := by
  have hp0 : p ≠ 0 := ne_of_gt hp
  have hx0 : x ≠ 0 := ne_of_gt hx
  have hx10 : x + 1 ≠ 0 := by linarith
  field_simp
  ring

/-- Positivity of the cleared polynomial is precisely the usable truncated
window margin. -/
theorem rho3_truncated_window_margin {p x : ℝ}
    (hp : 0 < p) (hx : 0 < x)
    (hpoly : 28 * x * p < 5 * (x - p) * (x + 1)) :
    2 / (x + 1) < 5 * ((x - p) / (14 * x * p)) := by
  have hid := rho3_truncated_window_identity hp hx
  have hden : 0 < 14 * x * p * (x + 1) := by positivity
  have hnum : 0 < 5 * (x - p) * (x + 1) - 28 * x * p := by
    linarith
  have hpos : 0 <
      (5 * (x - p) * (x + 1) - 28 * x * p) /
        (14 * x * p * (x + 1)) := div_pos hnum hden
  linarith

/-- The exact five first-good thresholds in the three-carrier dispatch.
For `p=1`, `x=9` is already the least legal row (`x ≥ p+8`). -/
theorem rho3_density_thresholds {p x : ℝ}
    (h : (p = 1 ∧ 9 ≤ x) ∨
      (p = 2 ∧ 13 ≤ x) ∨
      (p = 3 ∧ 19 ≤ x) ∨
      (p = 4 ∧ 26 ≤ x) ∨
      (p = 5 ∧ 33 ≤ x)) :
    28 * x * p < 5 * (x - p) * (x + 1) := by
  rcases h with h | h | h | h | h
  · rcases h with ⟨rfl, hx⟩
    nlinarith [sq_nonneg (x - 9)]
  · rcases h with ⟨rfl, hx⟩
    nlinarith [sq_nonneg (x - 13)]
  · rcases h with ⟨rfl, hx⟩
    nlinarith [sq_nonneg (x - 19)]
  · rcases h with ⟨rfl, hx⟩
    nlinarith [sq_nonneg (x - 26)]
  · rcases h with ⟨rfl, hx⟩
    nlinarith [sq_nonneg (x - 33)]

/-- The predecessor rows show that the displayed thresholds for `p=2,...,5`
cannot be lowered within the integer table. -/
theorem rho3_threshold_predecessors :
    ¬ 28 * (12 : ℤ) * 2 < 5 * (12 - 2) * (12 + 1) ∧
    ¬ 28 * (18 : ℤ) * 3 < 5 * (18 - 3) * (18 + 1) ∧
    ¬ 28 * (25 : ℤ) * 4 < 5 * (25 - 4) * (25 + 1) ∧
    ¬ 28 * (32 : ℤ) * 5 < 5 * (32 - 5) * (32 + 1) := by
  norm_num

/-- The three terminal rho-three candidates are geometrically trivial after
the exact finite reduction: every speed in `[14,16]` is safe throughout the
whole residual window `[1/182,1/70]`.  This removes the terminal wall-cell
bank from the formal boundary. -/
theorem rho3_terminal_candidate_safe_band {q u : ℝ}
    (hqLower : 14 ≤ q) (hqUpper : q ≤ 16)
    (huLower : (1 : ℝ) / 182 ≤ u) (huUpper : u ≤ 1 / 70) :
    (1 : ℝ) / 14 < q * u ∧ q * u < 13 / 14 := by
  have hq0 : 0 ≤ q := by linarith
  have hu0 : 0 ≤ u := by
    have : (0 : ℝ) < 1 / 182 := by norm_num
    linarith
  have hlow : (14 : ℝ) * (1 / 182) ≤ q * u :=
    mul_le_mul hqLower huLower (by norm_num) hq0
  have hupp : q * u ≤ (16 : ℝ) * (1 / 70) :=
    mul_le_mul hqUpper huUpper hu0 (by norm_num)
  constructor <;> norm_num at hlow hupp ⊢ <;> linarith

/-! ## Four and five carriers -/

/-- Exact four-carrier dimensionless margin `4*(6/7)-3=3/7`. -/
theorem rho4_dimensionless_margin :
    4 * ((6 : ℝ) / 7) - 3 = 3 / 7 := by
  norm_num

/-- Restoring the least-carrier scale preserves the strict rho-four margin. -/
theorem rho4_density_margin {x : ℝ} (hx : 0 < x) :
    3 / x < 4 * (6 / (7 * x)) := by
  have hden : 0 < 7 * x := by positivity
  have hid : 4 * (6 / (7 * x)) - 3 / x = 3 / (7 * x) := by
    field_simp
    ring
  have hpos : 0 < 3 / (7 * x) := div_pos (by norm_num) hden
  linarith

/-- The concrete first safe window of the least carrier cannot be covered by
the other three carrier danger combs.  This is the set-level `rho = 4` step
from THM-1214: it instantiates `multi_speed_density_bound` on the actual
Finset `{y,z,w}` and retains the strict reciprocal saving coming from the
ordered carrier speeds. -/
theorem rho4_actual_first_window_not_covered {x y z w : ℕ}
    (hxNat : 0 < x) (hxy : x < y) (hyz : y < z) (hzw : z < w) :
    ¬ Set.Icc (1 / (14 * (x : ℝ))) (13 / (14 * (x : ℝ))) ⊆
        ⋃ s ∈ ({y, z, w} : Finset ℕ), badArcs s (1 / 14) := by
  have hx : (0 : ℝ) < x := by exact_mod_cast hxNat
  let left : ℝ := 1 / (14 * (x : ℝ))
  let L : ℝ := 6 / (7 * (x : ℝ))
  have hL : 0 < L := by
    dsimp [L]
    positivity
  have hright : left + L = 13 / (14 * (x : ℝ)) := by
    dsimp [left, L]
    field_simp
    ring

  let S : Finset ℕ := {y, z, w}
  have hyzNe : y ≠ z := Nat.ne_of_lt hyz
  have hywNe : y ≠ w := Nat.ne_of_lt (lt_trans hyz hzw)
  have hzwNe : z ≠ w := Nat.ne_of_lt hzw
  have hcard : S.card = 3 := by simp [S, hyzNe, hywNe, hzwNe]
  have hSpos : ∀ s ∈ S, 1 ≤ s := by
    intro s hs
    simp only [S, Finset.mem_insert, Finset.mem_singleton] at hs
    rcases hs with rfl | rfl | rfl <;> omega

  intro hcover
  have hcover' : Set.Icc left (left + L) ⊆
      ⋃ s ∈ S, badArcs s (1 / 14) := by
    rw [hright]
    simpa [left, S] using hcover
  have hdensity := multi_speed_density_bound S hSpos (1 / 14) L left
    (by norm_num) (le_of_lt hL) hcover'

  have hsum : ∑ s ∈ S, (1 : ℝ) / s < 3 / (x : ℝ) := by
    have hxyReal : (x : ℝ) < y := by exact_mod_cast hxy
    have hxzReal : (x : ℝ) < z := by exact_mod_cast (lt_trans hxy hyz)
    have hxwReal : (x : ℝ) < w := by
      exact_mod_cast (lt_trans (lt_trans hxy hyz) hzw)
    have hyRec : (1 : ℝ) / y < 1 / x := one_div_lt_one_div_of_lt hx hxyReal
    have hzRec : (1 : ℝ) / z < 1 / x := one_div_lt_one_div_of_lt hx hxzReal
    have hwRec : (1 : ℝ) / w < 1 / x := one_div_lt_one_div_of_lt hx hxwReal
    have hsumExpand : ∑ s ∈ S, (1 : ℝ) / s =
        1 / (y : ℝ) + (1 / (z : ℝ) + 1 / (w : ℝ)) := by
      simp [S, hyzNe, hywNe, hzwNe]
    rw [hsumExpand]
    calc
      1 / (y : ℝ) + (1 / (z : ℝ) + 1 / (w : ℝ))
          < 1 / (x : ℝ) + (1 / (x : ℝ) + 1 / (x : ℝ)) := by
            linarith
      _ = 3 / (x : ℝ) := by ring

  have hbudget : 4 * L ≤ ∑ s ∈ S, (1 : ℝ) / s := by
    have hscaled :=
      mul_le_mul_of_nonneg_left hdensity (show (0 : ℝ) ≤ 7 by norm_num)
    rw [hcard] at hscaled
    convert hscaled using 1 <;> ring
  have hmargin : 3 / (x : ℝ) < 4 * L := by
    dsimp [L]
    exact rho4_density_margin hx
  exact (not_lt_of_ge hbudget) (lt_trans hsum hmargin)

/-- A first-window four-comb cover budget forces `y/x < 14/9`.  The variable
`budget` is the normalized reciprocal sum supplied by the density theorem. -/
theorem rho5_first_ratio_of_cover_budget {r budget : ℝ} (hr : 0 < r)
    (hcover : 3 * ((6 : ℝ) / 7) ≤ budget)
    (hbudget : budget < 4 / r) :
    r < 14 / 9 := by
  have hraw : 3 * ((6 : ℝ) / 7) < 4 / r := lt_of_le_of_lt hcover hbudget
  have hmul : 3 * ((6 : ℝ) / 7) * r < 4 :=
    (lt_div_iff₀ hr).mp hraw
  nlinarith

/-- Re-export of THM-1212's exact nested-window inequality at the stronger
five-killer first-ratio bound. -/
theorem rho5_nested_window_margin {r : ℝ} (hr : 0 < r)
    (hres : r < 14 / 9) :
    3 / r < 4 * LRC14.NestedCarrierWindow.phaseLength r :=
  LRC14.NestedCarrierWindow.rho5_density_contradiction hr hres

/-- Consequently the strict density inequality required by a cover of the
protected two-carrier prefix window is impossible. -/
theorem rho5_nested_cover_impossible {r : ℝ} (hr : 0 < r)
    (hres : r < 14 / 9) :
    ¬ 4 * LRC14.NestedCarrierWindow.phaseLength r < 3 / r := by
  intro hcover
  have hmargin := rho5_nested_window_margin hr hres
  linarith

/-- Arithmetic composition of both five-carrier cover handoffs. -/
theorem rho5_cover_budget_contradiction {r budget : ℝ} (hr : 0 < r)
    (hfirstCover : 3 * ((6 : ℝ) / 7) ≤ budget)
    (hfirstBudget : budget < 4 / r)
    (hnestedCover : 4 * LRC14.NestedCarrierWindow.phaseLength r < 3 / r) :
    False := by
  have hratio := rho5_first_ratio_of_cover_budget hr hfirstCover hfirstBudget
  exact rho5_nested_cover_impossible hr hratio hnestedCover

#print axioms owner_budget_iff
#print axioms rho2_endpoint_contradiction
#print axioms rho3_full_window_margin
#print axioms rho3_truncated_window_margin
#print axioms rho3_density_thresholds
#print axioms rho3_terminal_candidate_safe_band
#print axioms rho4_density_margin
#print axioms rho4_actual_first_window_not_covered
#print axioms rho5_first_ratio_of_cover_budget
#print axioms rho5_cover_budget_contradiction

end LRC14.FiveKillerCarrierWindow
