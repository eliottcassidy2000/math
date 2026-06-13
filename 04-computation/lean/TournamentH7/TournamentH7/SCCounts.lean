/-
  TournamentH7.SCCounts — SC tiling counts and THM-342 diagonal formulas

  ─── What this module records ──────────────────────────────────────────
  The project canon defines, for each n ≥ 2, the integer
        SC(n)  :=  #{path-fixed strongly connected tilings on n vertices}.
  Computed values (project canon, opus-2026-05-27-S1 / opus-2026-05-28-S3b):
        SC(2)  = 1,
        SC(3)  = 1,
        SC(4)  = 5,
        SC(5)  = 50,
        SC(6)  = 903,
        SC(7)  = 30773.

  ─── THM-340 (Composition Formula, project-novel) ─────────────────────
  Define `Q(d, k)` to be the number of labelled arrangements of `k`
  non-overlapping SC-tiling blocks of total width `d`.  Then
        Q(d, k)  =  [xᵈ] B(x)ᵏ,    where   B(x) = ∑_{a ≥ 2} SC(a+1) · xᵃ.

  ─── THM-342 (Diagonal Formulas, project-novel) ───────────────────────
  Letting `C(x) = B(x)/x²`:
        Q(2k,   k) = 1,                                   (j = 0 diagonal)
        Q(2k+1, k) = 5k,                                   (j = 1)
        Q(2k+2, k) = 25·k·(k+3)/2,                         (j = 2)
        Q(2k+3, k) = 903k + 250·k·(k−1) + 125·C(k, 3),     (j = 3)
        Q(2k+4, k) = 30773k + 4515·k·(k−1) + 2500·C(k, 2)
                     + 3750·C(k, 3) + 625·C(k, 4)          (j = 4)

  ─── Lean strategy ────────────────────────────────────────────────────
  Defining `Q(d, k)` rigorously requires formal power series in Mathlib
  (`Mathlib.RingTheory.PowerSeries`) plus a concrete tile-counting
  setup.  Here we *axiomatise* the diagonal formulas above and prove
  small derivations (e.g. the j = 0 case is just the definition of
  Q(2k, k); the j = 1 case is `5·k`, an arithmetic statement).
-/

import TournamentH7.Basic
import Mathlib.Data.Nat.Choose.Basic

namespace Tournament

/-! ### SC tiling counts (axiomatic constants) -/

/-- `SC n` = number of path-fixed strongly connected tilings on n vertices.
    Project canon constant; treated abstractly here. -/
opaque SCcount : ℕ → ℕ

@[simp] axiom SCcount_2 : SCcount 2 = 1
@[simp] axiom SCcount_3 : SCcount 3 = 1
@[simp] axiom SCcount_4 : SCcount 4 = 5
@[simp] axiom SCcount_5 : SCcount 5 = 50
@[simp] axiom SCcount_6 : SCcount 6 = 903
@[simp] axiom SCcount_7 : SCcount 7 = 30773

/-! ### Q-triangle: composition counts -/

/-- `Q d k` = number of labelled arrangements of k non-overlapping SC
    tiling blocks of total width d.  Project canon definition;
    axiomatised here. -/
opaque Qcount : ℕ → ℕ → ℕ

/-! ### THM-342 diagonal formulas (axioms) -/

/-- **Axiom (THM-342, j = 0 diagonal).** -/
axiom thm342_diag0 (k : ℕ) (hk : 1 ≤ k) : Qcount (2 * k) k = 1

/-- **Axiom (THM-342, j = 1 diagonal).** -/
axiom thm342_diag1 (k : ℕ) (hk : 1 ≤ k) : Qcount (2 * k + 1) k = 5 * k

/-- **Axiom (THM-342, j = 2 diagonal).** -/
axiom thm342_diag2 (k : ℕ) (hk : 1 ≤ k) :
    2 * Qcount (2 * k + 2) k = 25 * k * (k + 3)

/-- **Axiom (THM-342, j = 3 diagonal).** -/
axiom thm342_diag3 (k : ℕ) (hk : 1 ≤ k) :
    Qcount (2 * k + 3) k = 903 * k + 250 * k * (k - 1) + 125 * Nat.choose k 3

/-! ### Derived corollaries -/

/-- Q(2, 1) = 1.  (k = 1, j = 0.) -/
example : Qcount 2 1 = 1 := by
  have := thm342_diag0 1 (by omega)
  simp at this; exact this

/-- Q(3, 1) = 5.  (k = 1, j = 1.) -/
example : Qcount 3 1 = 5 := by
  have := thm342_diag1 1 (by omega)
  simpa using this

/-- Q(5, 2) = 10.  (k = 2, j = 1.) -/
example : Qcount 5 2 = 10 := by
  have := thm342_diag1 2 (by omega)
  -- 5 * 2 = 10
  simpa using this

/-- Q(4, 1) = 50, matching SC(5).  (k = 1, j = 2 reduces to 25·1·4/2 = 50.) -/
example : Qcount 4 1 = 50 := by
  have h := thm342_diag2 1 (by omega)
  -- h : 2 * Q(2·1 + 2) 1 = 25 * 1 * (1 + 3) = 100, hence Q(4, 1) = 50.
  show Qcount 4 1 = 50
  have h' : 2 * Qcount 4 1 = 100 := by
    have : 2 * Qcount (2 * 1 + 2) 1 = 25 * 1 * (1 + 3) := h
    simpa using this
  omega

end Tournament
