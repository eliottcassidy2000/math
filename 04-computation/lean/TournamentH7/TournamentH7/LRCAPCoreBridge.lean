/-
  TournamentH7.LRCAPCoreBridge — the AP-core bridge, elementary half of LRC(14).
  boxeph-2026-07-18-S105.

  The Lean payoff of the COMPLETED ELEMENTARY HALF of LRC(14) (THM-1017(II)):
  the compact / AP-core covering case reduces to LRC(≤13) via the descent floor.

  `ap_core_bridge` (PROVED, kernel-pure): if the largest speed `vstar` dominates
  every other speed 13-fold (`13·v_i ≤ v_max`, i.e. ρ = v_max/v_2nd ≥ 13), then the
  13-family is 1/14-lonely.  Proof: the LRC(≤13) citation makes the 12 non-max speeds
  1/13-lonely at some `t₀`; `descent_dominant` (THM-1008) lifts that to a 1/14-lonely
  time for the whole family.

  The open inverse theorem is intended to deliver this dominance on a compact
  small-margin class.  `lonely14_of_INV` records exactly that scoped
  composition: given `INV Compact` and the LRC(≤13) citation, every family in
  `Compact` is 1/14-lonely.  It does not prove that every covering family is in
  `Compact`; the easy/compact split and its easy-side witness are explicit in
  `LRC14DispatchAssembly`.
-/
import Mathlib
import TournamentH7.LonelyRunner
import TournamentH7.LRCDescentFloor
import TournamentH7.LRC13Citation

namespace LonelyRunner

/-- **The AP-core bridge (THM-1017(II), elementary half — PROVED).**
If the largest speed `vstar` dominates (`13·v_i ≤ v_max` for every other `i`, i.e.
ρ = v_max/v_2nd ≥ 13), the 13-family is `1/14`-lonely.  LRC(≤13) makes the 12 non-max
speeds `1/13`-lonely at some `t₀`; the descent floor `descent_dominant` lifts it. -/
theorem ap_core_bridge (cite : LRCUpTo13) (v : Fin 13 → ℤ)
    (hpos : ∀ i, 0 < v i) (vstar : Fin 13)
    (hdom : ∀ i, i ≠ vstar → 13 * v i ≤ v vstar) :
    ∃ t : ℝ, Lonely 14 v t := by
  -- the 12 non-max speeds, reindexed through `vstar.succAbove : Fin 12 ↪ Fin 13`
  set w : Fin 12 → ℤ := fun j => v (vstar.succAbove j) with hw
  have hwne : ∀ j, w j ≠ 0 := fun j => (hpos _).ne'
  -- LRC(≤13) citation: 12 speeds ⟹ a 1/13-lonely time
  obtain ⟨t0, ht0⟩ := cite 12 (le_refl 12) w hwne
  -- translate the citation's loneliness to the descent-floor's `hkept`
  have hkept : ∀ i, i ≠ vstar → ∀ m : ℤ, (1 : ℝ) / 13 ≤ |(v i : ℝ) * t0 - m| := by
    intro i hi m
    obtain ⟨j, hj⟩ := Fin.exists_succAbove_eq hi
    have h := ht0 j m
    simp only [hw, hj] at h
    have hcast : (1 : ℝ) / ((12 + 1 : ℕ) : ℝ) = (1 : ℝ) / 13 := by norm_num
    rwa [hcast] at h
  exact descent_dominant v vstar t0 hkept hpos hdom

/-- **The AP-core bridge in explicit deep-well shape (THM-1017(II) mechanism).**
If the 12 non-max speeds are bounded by `12·d` (a dilated-AP core `d·{1,…,12}` has
`v_i = d·i ≤ 12·d`) and the far element clears `156·d = 13·(12·d)`, then the family is
`1/14`-lonely.  The inverse theorem forces the far element to be a multiple of
`lcm(13,14)·d = 182·d ≥ 156·d`, so this is precisely the compact-case discharge. -/
theorem ap_core_bridge_of_shape (cite : LRCUpTo13) (v : Fin 13 → ℤ)
    (hpos : ∀ i, 0 < v i) (vstar : Fin 13) (d : ℤ)
    (hcore : ∀ i, i ≠ vstar → v i ≤ 12 * d)
    (hfar : 156 * d ≤ v vstar) :
    ∃ t : ℝ, Lonely 14 v t := by
  refine ap_core_bridge cite v hpos vstar ?_
  intro i hi
  have h1 : 13 * v i ≤ 13 * (12 * d) := by
    have := hcore i hi; nlinarith [this]
  have h2 : (13 : ℤ) * (12 * d) = 156 * d := by ring
  linarith [h1, h2, hfar]

/-- **The inverse theorem `INV` (OPEN, dominance form).**  Stated as the conclusion
that a compact small-margin class needs: every 13-family in the class encoded
by the abstract predicate `Compact` has its non-max speeds
dominated 13-fold by `v_max`.  This is the intended rho-at-least-13 form of
the compact inverse theorem.  It enters as a hypothesis, never a `sorry`;
this definition itself makes no
claim that `Compact` is the whole covering class. -/
def INV (Compact : (Fin 13 → ℤ) → Prop) : Prop :=
  ∀ (v : Fin 13 → ℤ), (∀ i, 0 < v i) → Compact v →
    ∃ vstar : Fin 13, ∀ i, i ≠ vstar → 13 * v i ≤ v vstar

/-- **The compact-case composition.**  Given the LRC(≤13) citation and the inverse
theorem `INV`, every compact 13-family is `1/14`-lonely.  This is the elementary half
assembled: `INV` supplies the dominance, `ap_core_bridge` supplies the loneliness. -/
theorem lonely14_of_INV {Compact : (Fin 13 → ℤ) → Prop}
    (cite : LRCUpTo13) (inv : INV Compact)
    (v : Fin 13 → ℤ) (hpos : ∀ i, 0 < v i) (hc : Compact v) :
    ∃ t : ℝ, Lonely 14 v t := by
  obtain ⟨vstar, hdom⟩ := inv v hpos hc
  exact ap_core_bridge cite v hpos vstar hdom

end LonelyRunner

#print axioms LonelyRunner.ap_core_bridge
#print axioms LonelyRunner.ap_core_bridge_of_shape
#print axioms LonelyRunner.lonely14_of_INV
