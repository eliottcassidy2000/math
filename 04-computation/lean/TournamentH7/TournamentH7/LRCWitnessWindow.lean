/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-02-S35)
-/
import TournamentH7.LRCWitnessCert
import TournamentH7.LonelyRunnerMathlib
import Mathlib.Tactic.IntervalCases

/-!
# Below-V* finite windows for the witness-certificate shapes (DAG surface item 1 of 3)

`LRCWitnessCert` certifies each shape `(P, offs)` for the V-TAIL (`V ≥ V*`).  This file closes
the finite windows BELOW `V*`: for every admissible `V` in `[max offs + 1, V*)`, an explicit
rational witness time `a/b`, checked by kernel `decide` through the closed Nat criterion
`b ≤ 14 · min((v·a) % b, b − (v·a) % b)` (the exact form of `‖v·(a/b)‖ ≥ 1/14` via
`LonelyRunner.norm_natCast_div`).  Combined `_all` theorems (window + tail) remove the V-gap
entirely: each certified 13-runner family now holds for EVERY `V` with positive speeds.

Witness tables generated and independently verified by
`04-computation/lrc14_window_witness_generator_opus_20260702_S35.py` (exact max-min over the
complete breakpoint grid; 53 rows, all pass).
-/

namespace LonelyRunner
namespace WitnessCert

/-- The per-runner window check: speed `w ≥ 0` is `1/14`-safe at time `a/b`.  This is the
closed Nat form of `‖w·(a/b)‖ ≥ 1/14`. -/
theorem norm_ge_of_natCheck {w : ℤ} (hw : 0 ≤ w) {a b : ℕ} (hb : 0 < b)
    (hchk : b ≤ 14 * min ((w.toNat * a) % b) (b - (w.toNat * a) % b)) :
    ((h14 : ℚ) : ℝ) ≤ ‖(((w : ℝ) * ((a : ℝ) / (b : ℝ)) : ℝ) : UnitAddCircle)‖ := by
  have hw' : ((w.toNat : ℕ) : ℝ) = (w : ℝ) := by exact_mod_cast Int.toNat_of_nonneg hw
  rw [← hw']
  have hcast : ((w.toNat : ℕ) : ℝ) * ((a : ℝ) / (b : ℝ)) = ((w.toNat * a : ℕ) : ℝ) / (b : ℝ) := by
    push_cast; ring
  rw [hcast, norm_natCast_div]
  have hb' : (0 : ℝ) < (b : ℝ) := by exact_mod_cast hb
  have hchk' : ((b : ℕ) : ℝ) ≤
      14 * ((min ((w.toNat * a) % b) (b - (w.toNat * a) % b) : ℕ) : ℝ) := by
    exact_mod_cast hchk
  have h14v : ((h14 : ℚ) : ℝ) = 1 / 14 := by norm_num [h14]
  rw [h14v, le_div_iff₀ hb']
  linarith

/-- A window row: witness `a/b` certifies the full shape `(P, offs)` at reference speed `V`.
Decidable by kernel `decide` (small numerals). -/
def windowRowOK (P offs : List ℤ) (V : ℤ) (a b : ℕ) : Prop :=
  0 < b ∧
  (∀ s ∈ P, 0 ≤ s ∧ b ≤ 14 * min ((s.toNat * a) % b) (b - (s.toNat * a) % b)) ∧
  (∀ o ∈ offs, 0 ≤ V - o ∧
    b ≤ 14 * min (((V - o).toNat * a) % b) (b - ((V - o).toNat * a) % b))

instance (P offs : List ℤ) (V : ℤ) (a b : ℕ) : Decidable (windowRowOK P offs V a b) := by
  unfold windowRowOK; infer_instance

/-- A checked window row yields the lonely time for its shape instance. -/
theorem lonely_of_windowRowOK {P offs : List ℤ} {V : ℤ} {a b : ℕ} (h : windowRowOK P offs V a b) :
    ∃ τ : ℝ,
      (∀ s ∈ P, ((h14 : ℚ) : ℝ) ≤ ‖(((s : ℝ) * τ : ℝ) : UnitAddCircle)‖) ∧
      (∀ o ∈ offs, ((h14 : ℚ) : ℝ) ≤ ‖((((V - o : ℤ) : ℝ) * τ : ℝ) : UnitAddCircle)‖) := by
  obtain ⟨hb, hP, hO⟩ := h
  refine ⟨(a : ℝ) / (b : ℝ), fun s hs => ?_, fun o ho => ?_⟩
  · obtain ⟨h0, hchk⟩ := hP s hs
    exact norm_ge_of_natCheck h0 hb hchk
  · obtain ⟨h0, hchk⟩ := hO o ho
    exact norm_ge_of_natCheck h0 hb hchk

/-! ## The AP shape: window `V ∈ [13, 15)`, then `certAP_tail` -/

theorem certAP_window (V : ℤ) (h1 : 13 ≤ V) (h2 : V < 15) :
    ∃ τ : ℝ, ∀ o ∈ offsAP, ((h14 : ℚ) : ℝ) ≤
      ‖((((V - o : ℤ) : ℝ) * τ : ℝ) : UnitAddCircle)‖ := by
  interval_cases V
  · obtain ⟨τ, -, hO⟩ := lonely_of_windowRowOK (P := ([] : List ℤ)) (offs := offsAP)
      (V := 13) (a := 1) (b := 14) (by decide)
    exact ⟨τ, hO⟩
  · obtain ⟨τ, -, hO⟩ := lonely_of_windowRowOK (P := ([] : List ℤ)) (offs := offsAP)
      (V := 14) (a := 1) (b := 16) (by decide)
    exact ⟨τ, hO⟩

/-- **The dilated-AP family for every `V ≥ 13`** — window and tail combined; no V-gap. -/
theorem certAP_all (V : ℤ) (hV : 13 ≤ V) :
    ∃ τ : ℝ, ∀ o ∈ offsAP, ((h14 : ℚ) : ℝ) ≤
      ‖((((V - o : ℤ) : ℝ) * τ : ℝ) : UnitAddCircle)‖ := by
  rcases lt_or_ge V 15 with h | h
  · exact certAP_window V hV h
  · exact certAP_tail V h

/-! ## Shape 3: window `V ∈ [21, 47)`, then `cert3_tail` -/

theorem cert3_window (V : ℤ) (h1 : 21 ≤ V) (h2 : V < 47) :
    ∃ τ : ℝ,
      (∀ s ∈ P3, ((h14 : ℚ) : ℝ) ≤ ‖(((s : ℝ) * τ : ℝ) : UnitAddCircle)‖) ∧
      (∀ o ∈ offs3, ((h14 : ℚ) : ℝ) ≤ ‖((((V - o : ℤ) : ℝ) * τ : ℝ) : UnitAddCircle)‖) := by
  interval_cases V
  · exact lonely_of_windowRowOK (P := P3) (offs := offs3) (V := 21) (a := 3) (b := 17) (by decide)
  · exact lonely_of_windowRowOK (P := P3) (offs := offs3) (V := 22) (a := 1) (b := 6) (by decide)
  · exact lonely_of_windowRowOK (P := P3) (offs := offs3) (V := 23) (a := 1) (b := 8) (by decide)
  · exact lonely_of_windowRowOK (P := P3) (offs := offs3) (V := 24) (a := 2) (b := 17) (by decide)
  · exact lonely_of_windowRowOK (P := P3) (offs := offs3) (V := 25) (a := 1) (b := 9) (by decide)
  · exact lonely_of_windowRowOK (P := P3) (offs := offs3) (V := 26) (a := 4) (b := 27) (by decide)
  · exact lonely_of_windowRowOK (P := P3) (offs := offs3) (V := 27) (a := 4) (b := 23) (by decide)
  · exact lonely_of_windowRowOK (P := P3) (offs := offs3) (V := 28) (a := 1) (b := 6) (by decide)
  · exact lonely_of_windowRowOK (P := P3) (offs := offs3) (V := 29) (a := 3) (b := 19) (by decide)
  · exact lonely_of_windowRowOK (P := P3) (offs := offs3) (V := 30) (a := 3) (b := 26) (by decide)
  · exact lonely_of_windowRowOK (P := P3) (offs := offs3) (V := 31) (a := 1) (b := 8) (by decide)
  · exact lonely_of_windowRowOK (P := P3) (offs := offs3) (V := 32) (a := 7) (b := 39) (by decide)
  · exact lonely_of_windowRowOK (P := P3) (offs := offs3) (V := 33) (a := 5) (b := 29) (by decide)
  · exact lonely_of_windowRowOK (P := P3) (offs := offs3) (V := 34) (a := 1) (b := 6) (by decide)
  · exact lonely_of_windowRowOK (P := P3) (offs := offs3) (V := 35) (a := 4) (b := 25) (by decide)
  · exact lonely_of_windowRowOK (P := P3) (offs := offs3) (V := 36) (a := 7) (b := 32) (by decide)
  · exact lonely_of_windowRowOK (P := P3) (offs := offs3) (V := 37) (a := 11) (b := 52) (by decide)
  · exact lonely_of_windowRowOK (P := P3) (offs := offs3) (V := 38) (a := 3) (b := 17) (by decide)
  · exact lonely_of_windowRowOK (P := P3) (offs := offs3) (V := 39) (a := 6) (b := 35) (by decide)
  · exact lonely_of_windowRowOK (P := P3) (offs := offs3) (V := 40) (a := 1) (b := 6) (by decide)
  · exact lonely_of_windowRowOK (P := P3) (offs := offs3) (V := 41) (a := 8) (b := 37) (by decide)
  · exact lonely_of_windowRowOK (P := P3) (offs := offs3) (V := 42) (a := 13) (b := 62) (by decide)
  · exact lonely_of_windowRowOK (P := P3) (offs := offs3) (V := 43) (a := 11) (b := 61) (by decide)
  · exact lonely_of_windowRowOK (P := P3) (offs := offs3) (V := 44) (a := 7) (b := 40) (by decide)
  · exact lonely_of_windowRowOK (P := P3) (offs := offs3) (V := 45) (a := 7) (b := 41) (by decide)
  · exact lonely_of_windowRowOK (P := P3) (offs := offs3) (V := 46) (a := 1) (b := 6) (by decide)

/-- **The `{1,2,3} ∪ {V−o}` family for every `V ≥ 21`** — window + tail; no V-gap. -/
theorem cert3_all (V : ℤ) (hV : 21 ≤ V) :
    ∃ τ : ℝ,
      (∀ s ∈ P3, ((h14 : ℚ) : ℝ) ≤ ‖(((s : ℝ) * τ : ℝ) : UnitAddCircle)‖) ∧
      (∀ o ∈ offs3, ((h14 : ℚ) : ℝ) ≤ ‖((((V - o : ℤ) : ℝ) * τ : ℝ) : UnitAddCircle)‖) := by
  rcases lt_or_ge V 47 with h | h
  · exact cert3_window V hV h
  · exact cert3_tail V h

/-! ## Shape 4: window `V ∈ [15, 40)`, then `cert4_tail` -/

theorem cert4_window (V : ℤ) (h1 : 15 ≤ V) (h2 : V < 40) :
    ∃ τ : ℝ,
      (∀ s ∈ P4, ((h14 : ℚ) : ℝ) ≤ ‖(((s : ℝ) * τ : ℝ) : UnitAddCircle)‖) ∧
      (∀ o ∈ offs4, ((h14 : ℚ) : ℝ) ≤ ‖((((V - o : ℤ) : ℝ) * τ : ℝ) : UnitAddCircle)‖) := by
  interval_cases V
  · exact lonely_of_windowRowOK (P := P4) (offs := offs4) (V := 15) (a := 1) (b := 8) (by decide)
  · exact lonely_of_windowRowOK (P := P4) (offs := offs4) (V := 16) (a := 1) (b := 9) (by decide)
  · exact lonely_of_windowRowOK (P := P4) (offs := offs4) (V := 17) (a := 1) (b := 8) (by decide)
  · exact lonely_of_windowRowOK (P := P4) (offs := offs4) (V := 18) (a := 1) (b := 20) (by decide)
  · exact lonely_of_windowRowOK (P := P4) (offs := offs4) (V := 19) (a := 1) (b := 8) (by decide)
  · exact lonely_of_windowRowOK (P := P4) (offs := offs4) (V := 20) (a := 1) (b := 11) (by decide)
  · exact lonely_of_windowRowOK (P := P4) (offs := offs4) (V := 21) (a := 1) (b := 8) (by decide)
  · exact lonely_of_windowRowOK (P := P4) (offs := offs4) (V := 22) (a := 1) (b := 24) (by decide)
  · exact lonely_of_windowRowOK (P := P4) (offs := offs4) (V := 23) (a := 1) (b := 8) (by decide)
  · exact lonely_of_windowRowOK (P := P4) (offs := offs4) (V := 24) (a := 2) (b := 17) (by decide)
  · exact lonely_of_windowRowOK (P := P4) (offs := offs4) (V := 25) (a := 1) (b := 8) (by decide)
  · exact lonely_of_windowRowOK (P := P4) (offs := offs4) (V := 26) (a := 2) (b := 19) (by decide)
  · exact lonely_of_windowRowOK (P := P4) (offs := offs4) (V := 27) (a := 1) (b := 8) (by decide)
  · exact lonely_of_windowRowOK (P := P4) (offs := offs4) (V := 28) (a := 7) (b := 19) (by decide)
  · exact lonely_of_windowRowOK (P := P4) (offs := offs4) (V := 29) (a := 1) (b := 8) (by decide)
  · exact lonely_of_windowRowOK (P := P4) (offs := offs4) (V := 30) (a := 2) (b := 21) (by decide)
  · exact lonely_of_windowRowOK (P := P4) (offs := offs4) (V := 31) (a := 1) (b := 8) (by decide)
  · exact lonely_of_windowRowOK (P := P4) (offs := offs4) (V := 32) (a := 3) (b := 25) (by decide)
  · exact lonely_of_windowRowOK (P := P4) (offs := offs4) (V := 33) (a := 1) (b := 8) (by decide)
  · exact lonely_of_windowRowOK (P := P4) (offs := offs4) (V := 34) (a := 1) (b := 9) (by decide)
  · exact lonely_of_windowRowOK (P := P4) (offs := offs4) (V := 35) (a := 1) (b := 8) (by decide)
  · exact lonely_of_windowRowOK (P := P4) (offs := offs4) (V := 36) (a := 10) (b := 27) (by decide)
  · exact lonely_of_windowRowOK (P := P4) (offs := offs4) (V := 37) (a := 1) (b := 8) (by decide)
  · exact lonely_of_windowRowOK (P := P4) (offs := offs4) (V := 38) (a := 11) (b := 29) (by decide)
  · exact lonely_of_windowRowOK (P := P4) (offs := offs4) (V := 39) (a := 1) (b := 8) (by decide)

/-- **The `{2,…,6} ∪ {V−even}` family for every `V ≥ 15`** — window + tail; no V-gap. -/
theorem cert4_all (V : ℤ) (hV : 15 ≤ V) :
    ∃ τ : ℝ,
      (∀ s ∈ P4, ((h14 : ℚ) : ℝ) ≤ ‖(((s : ℝ) * τ : ℝ) : UnitAddCircle)‖) ∧
      (∀ o ∈ offs4, ((h14 : ℚ) : ℝ) ≤ ‖((((V - o : ℤ) : ℝ) * τ : ℝ) : UnitAddCircle)‖) := by
  rcases lt_or_ge V 40 with h | h
  · exact cert4_window V hV h
  · exact cert4_tail V h

end WitnessCert
end LonelyRunner
