import Mathlib.Data.Nat.Basic
import Mathlib.Tactic

/-
# TNC detection depth: D consecutive zeros propagate forward

kind-pasteur-2026-07-20-S128c124.  The Lean core of THM-1710(i).

The toral moment sequence `a_m = [u^{Mm}] R^m` obeys an order-`D` linear recurrence
`Σ_{i=0}^D P_i(m) a_{m+i} = 0` (THM-1670).  When the leading coefficient `P_D(m)` is
nonzero, the recurrence can be solved forward for `a_{m+D}`, and in particular:

  if the `D` terms `a_m, …, a_{m+D-1}` all vanish, then `a_{m+D}` vanishes.

The lemma below is the pure consequence of that single step: `D` consecutive initial zeros
`a_1 = ⋯ = a_D = 0` force `a_m = 0` for every `m ≥ 1`.  Hence the toral nullcone is cut out
by the FIRST `D` moments — the "detection depth `≤ D`" cap — with no analysis, no DvdK.

`step` is stated abstractly (it is exactly what a nonvanishing leading coefficient plus the
recurrence provide), so this file is `sorry`-free and self-contained; the identification of
`step` with the recurrence is in the prose of THM-1710.  It is the forward-propagation dual
of `GMC2Hermite.ThreeTerm.no_common_root`.
-/

namespace TNCDepth

/-- **Forward zero-propagation.**  If, for every `m ≥ 1`, the vanishing of the `D`
consecutive terms `a m, …, a (m+D-1)` forces `a (m+D) = 0` (`step` — supplied by the
order-`D` recurrence with nonvanishing leading coefficient), and the first `D` terms
`a 1, …, a D` all vanish (`init`), then `a m = 0` for all `m ≥ 1`.

This is THM-1710(i): the toral nullcone is determined by the first `D` moments. -/
theorem zeros_propagate {α : Type*} [Zero α] {a : ℕ → α} {D : ℕ} (_hD : 1 ≤ D)
    (step : ∀ m, 1 ≤ m → (∀ i, i < D → a (m + i) = 0) → a (m + D) = 0)
    (init : ∀ j, 1 ≤ j → j ≤ D → a j = 0) :
    ∀ m, 1 ≤ m → a m = 0 := by
  intro m
  induction m using Nat.strong_induction_on with
  | _ m ih =>
    intro hm
    by_cases hmD : m ≤ D
    · -- m ≤ D : covered by the initial data
      exact init m hm hmD
    · -- m > D : write m = (m - D) + D and apply `step` at m - D
      have hmD' : D ≤ m := by omega
      have hm' : 1 ≤ m - D := by omega
      have hstep := step (m - D) hm' ?_
      · -- a ((m - D) + D) = 0, and (m - D) + D = m
        have hcancel : (m - D) + D = m := Nat.sub_add_cancel hmD'
        rwa [hcancel] at hstep
      · -- the D predecessors a (m-D+i), i < D, all vanish by strong induction
        intro i hi
        have hlt : (m - D) + i < m := by omega
        have hge : 1 ≤ (m - D) + i := by omega
        exact ih ((m - D) + i) hlt hge

/-- The contrapositive packaging used in THM-1710(iv): if some `a m` with `m ≥ 1` is
nonzero, then one of the first `D` moments is already nonzero.  So a two-sided `R` escapes
the toral nullcone already at depth `≤ D`. -/
theorem nonzero_within_depth {α : Type*} [Zero α] {a : ℕ → α} {D : ℕ} (_hD : 1 ≤ D)
    (step : ∀ m, 1 ≤ m → (∀ i, i < D → a (m + i) = 0) → a (m + D) = 0)
    (hne : ∃ m, 1 ≤ m ∧ a m ≠ 0) :
    ∃ j, 1 ≤ j ∧ j ≤ D ∧ a j ≠ 0 := by
  by_contra hcon
  obtain ⟨m, hm, ham⟩ := hne
  refine ham (zeros_propagate _hD step (fun j hj hjD => ?_) m hm)
  by_contra haj
  exact hcon ⟨j, hj, hjD, haj⟩

end TNCDepth
