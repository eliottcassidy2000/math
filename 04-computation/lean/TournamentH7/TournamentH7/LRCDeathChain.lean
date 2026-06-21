/-
  TournamentH7.LRCDeathChain -- finite death-chain algebra for LRC14.

  This module formalizes the proof-order quotient used in HYP-2708: before
  doing analytic estimates, record which missed-sector depths can change under
  r far hits for a chosen row currency.
-/

import Mathlib.Tactic

namespace LonelyRunner

/-- A row currency depends only on the number of missed inner sectors. -/
abbrev Currency := Nat -> Int

/-- The direct cover currency for `p0`: one unit exactly at depth zero. -/
def coverCurrency : Currency :=
  fun t => if t = 0 then 1 else 0

/-- The survival currency from HYP-2701:
`C = p1 + p2 + p3 + p4 - 4*p6`, indexed by missed depth. -/
def survivalCurrency : Currency :=
  fun t => if 1 <= t /\ t <= 4 then 1 else if t = 6 then -4 else 0

/-- A hit count reachable from depth `t` after `r` far hits. -/
def ReachableHit (r t : Nat) := Fin (Nat.min r t + 1)

/-- Boolean live-depth witness search. -/
def liveDepthBool (C : Currency) (r t : Nat) : Bool :=
  (List.range (Nat.min r t + 1)).any (fun h => decide (C (t - h) ≠ C t))

/-- Boolean silent-depth witness search. -/
def silentDepthBool (C : Currency) (r t : Nat) : Bool :=
  (List.range (Nat.min r t + 1)).all (fun h => decide (C (t - h) = C t))

/-- A depth is live when some reachable hit count changes the currency. -/
abbrev LiveDepth (C : Currency) (r t : Nat) : Prop :=
  liveDepthBool C r t = true

/-- A depth is silent when every reachable hit count preserves the currency. -/
abbrev SilentDepth (C : Currency) (r t : Nat) : Prop :=
  silentDepthBool C r t = true

/-- Live depths are monotone in the number of far hits. -/
theorem liveDepth_mono {C : Currency} {r s t : Nat} (hrs : r <= s) :
    LiveDepth C r t -> LiveDepth C s t := by
  unfold LiveDepth liveDepthBool
  intro h
  rw [List.any_eq_true] at h ⊢
  rcases h with ⟨a, ha_mem, ha_change⟩
  refine ⟨a, ?_, ha_change⟩
  have hmin : Nat.min r t <= Nat.min s t :=
    le_min (le_trans (min_le_left r t) hrs) (min_le_right r t)
  exact List.mem_range.mpr (lt_of_lt_of_le (List.mem_range.mp ha_mem)
    (Nat.succ_le_succ hmin))

/-- Direct `p0`, one far hit: only one missed sector can change the row. -/
theorem cover_oneFar_live_iff (t : Nat) (ht : t <= 6) :
    LiveDepth coverCurrency 1 t <-> t = 1 := by
  interval_cases t <;> native_decide

/-- Survival currency, one far hit: shallow closure plus the two high-tail depths. -/
theorem survival_oneFar_live_iff (t : Nat) (ht : t <= 6) :
    LiveDepth survivalCurrency 1 t <-> t = 1 \/ t = 5 \/ t = 6 := by
  interval_cases t <;> native_decide

/-- Survival currency, two far hits: HYP-2708 has exactly four live depths. -/
theorem survival_twoFar_live_iff (t : Nat) (ht : t <= 6) :
    LiveDepth survivalCurrency 2 t <-> t = 1 \/ t = 2 \/ t = 5 \/ t = 6 := by
  interval_cases t <;> native_decide

/-- Survival currency, three far hits: only depth four remains silent above zero. -/
theorem survival_threeFar_live_iff (t : Nat) (ht : t <= 6) :
    LiveDepth survivalCurrency 3 t <-> t = 1 \/ t = 2 \/ t = 3 \/ t = 5 \/ t = 6 := by
  interval_cases t <;> native_decide

/-- Survival currency, four far hits: every positive depth through six is live. -/
theorem survival_fourFar_live_iff (t : Nat) (ht : t <= 6) :
    LiveDepth survivalCurrency 4 t <-> 1 <= t := by
  interval_cases t <;> native_decide

/-- The two-far middle plateau is exactly silent. -/
theorem survival_twoFar_middle_silent :
    SilentDepth survivalCurrency 2 3 /\ SilentDepth survivalCurrency 2 4 := by
  native_decide

/-! ### Axiom audit

These finite computations should use no project-specific axioms.  In a working
Lean environment, `#print axioms` should report only Lean's ordinary
foundational axioms, and typically none for the purely decidable theorems.
-/

#print axioms liveDepth_mono
#print axioms cover_oneFar_live_iff
#print axioms survival_oneFar_live_iff
#print axioms survival_twoFar_live_iff
#print axioms survival_threeFar_live_iff
#print axioms survival_fourFar_live_iff
#print axioms survival_twoFar_middle_silent

end LonelyRunner
