/-
klein-2026-07-02-S99 (HYP-4006b) — the Lean WINDOW RE-ENUMERATION pilot (N1 [LEAN] flip).

Fully computable re-enumeration of the THM-534 L_y windows inside Lean: exact ℚ interval
arithmetic (bad-interval list → sorted merge → measure), the sector moments S₁, S₂, the
row-11/12/13 dual L_y = 1 − ½S₁ + ⅙S₂, and `native_decide` over every shape in the
bounded-spread windows (k = 11: 286 shapes, k = 12: 78, k = 13: 13). Mirrors the Python
censuses (ly_windows_k10_13_klein.py) — same maxima, now kernel^native-checked.
-/
import Mathlib

namespace LyWindowEnum

/-- Sorted-merge measure of a union of intervals. -/
def measUnion (l : List (ℚ × ℚ)) : ℚ :=
  let s := l.mergeSort (fun a b => a.1 ≤ b.1)
  let f : (ℚ × ℚ) → (ℚ × ℚ) → (ℚ × ℚ) := fun acc iv =>
    -- acc = (total so far EXCLUDING current block, current block right end packed via fst)
    acc  -- placeholder, replaced by explicit recursion below
  (go s 0 0 0).1
where
  go : List (ℚ × ℚ) → ℚ → ℚ → ℚ → ℚ × Unit
    | [], tot, cl, ch => (tot + (ch - cl), ())
    | (lo, hi) :: t, tot, cl, ch =>
      if lo ≤ ch then go t tot cl (max ch hi)
      else go t (tot + (ch - cl)) lo hi

/-- Bad intervals for speeds `E`, avoided sectors `A` (sector j = [j/7,(j+1)/7)). -/
def bad (E A : List ℕ) : List (ℚ × ℚ) :=
  E.flatMap fun e =>
    if e = 0 then [] else
    A.flatMap fun j =>
      (List.range e).map fun a =>
        (((a : ℚ) + j / 7) / e, ((a : ℚ) + (j + 1) / 7) / e)

/-- The avoidance measure J(A, E). -/
def J (E A : List ℕ) : ℚ := 1 - measUnion (bad E A)

def sectors (r : ℕ) : List (List ℕ) := ((List.range 6).map (· + 1)).sublistsLen r

def S (E : List ℕ) (r : ℕ) : ℚ := ((sectors r).map (J E)).sum

/-- The rows-11/12/13 dual: L_y = 1 − ½S₁ + ⅙S₂. -/
def Ly (E : List ℕ) : ℚ := 1 - S E 1 / 2 + S E 2 / 6

def shapes (k spread : ℕ) : List (List ℕ) :=
  (((List.range spread).map (· + 1)).sublistsLen (k - 1)).map (0 :: ·)

/-- k = 13 window (13 shapes): every L_y ≤ cap₁₃ = 1. -/
theorem window13 : (shapes 13 13).all (fun E => Ly E ≤ 1) = true := by native_decide
/-- k = 12 window (78 shapes): every L_y ≤ cap₁₂ = 78/91. -/
theorem window12 : (shapes 12 13).all (fun E => Ly E ≤ 78/91) = true := by native_decide
/-- k = 11 window (286 shapes): every L_y ≤ cap₁₁ = 66/91. -/
theorem window11 : (shapes 11 13).all (fun E => Ly E ≤ 66/91) = true := by native_decide

/-- Third sector moment (for the rows-9/10 degree-3 dual). -/
def S3 (E : List ℕ) : ℚ := ((sectors 3).map (J E)).sum

/-- The rows-9/10 dual: L_y = 1 − (13/18)S₁ + (4/9)S₂ − (1/6)S₃. -/
def Ly3 (E : List ℕ) : ℚ := 1 - 13 * S E 1 / 18 + 4 * S E 2 / 9 - S3 E / 6

/-- k = 10 window (2002 shapes, degree-3 dual): every L_y ≤ cap₁₀ = 55/91. -/
theorem window10 : (shapes 10 14).all (fun E => Ly3 E ≤ 55/91) = true := by native_decide

/-- k = 9 window (6435 shapes, same degree-3 dual): every L_y ≤ cap₉ = 1979/4004. -/
theorem window9 : (shapes 9 15).all (fun E => Ly3 E ≤ 1979/4004) = true := by native_decide

/-- Fourth sector moment (for the row-8 degree-4 dual). -/
def S4 (E : List ℕ) : ℚ := ((sectors 4).map (J E)).sum

/-- The row-8 dual: L_y = 1 − S₁ + S₂ − (9/10)S₃ + (3/5)S₄. -/
def Ly4 (E : List ℕ) : ℚ := 1 - S E 1 + S E 2 - 9 * S3 E / 10 + 3 * S4 E / 5

/-- k = 8 window, chunked by second element: shapes with tail starting at m. -/
def shapes8 (m : ℕ) : List (List ℕ) :=
  ((((List.range 16).map (· + 1)).filter (· > m)).sublistsLen 6).map (fun t => 0 :: m :: t)

/-- k = 8 window (11440 shapes total, degree-4 dual), verified in two chunks. -/
theorem window8a : ((List.range 4).map (· + 1)).all
    (fun m => (shapes8 m).all (fun E => Ly4 E ≤ 2243/5880)) = true := by native_decide
theorem window8b : (((List.range 12).map (· + 5))).all
    (fun m => (shapes8 m).all (fun E => Ly4 E ≤ 2243/5880)) = true := by native_decide

end LyWindowEnum


