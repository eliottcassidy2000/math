/-
  TournamentH7.LRCGapSurplusLedger — the section-formula ledger, machine-checked.
  kind-pasteur-2026-07-01-S34 (HYP-3957; THM-599 + HYP-3956).

  THM-599 (torus-band) restated Lean-first (the SECTION FORMULA): the c-averaged lonely measure is
      A(U) = ∫₀¹ F_U(x) dx,   F_U(x) = Σ over circular gaps g of {frac(u·x)} of (g − w)⁺,  w = 1/7,
  and F_U is piecewise-affine in x with breakpoints only at { j/m } ∪ { (7j±1)/(7m) } over the
  pairwise differences m of U.  Hence A(U) equals the EXACT midpoint sum over breakpoint intervals
  — a finite rational computation.  This file:

   * implements `Fgap` (the gap surplus) and `Aledger` (the breakpoint midpoint sum) as computable
     functions on ℚ — no reals, no measure theory;
   * discharges the LEDGER ROWS by `native_decide`: the exact rationals of the k = 2..13 argmin
     patterns (36/49, 61/98, 11/21, …, 144699147/1267426160) and the floor comparisons against
     witnessMP = 14249/252252 — the finite table of the (⋆)-census (HYP-3953/3955);
   * states the ANALYTIC BRIDGE (Aledger = the true ∫F; THM-599 (SF)+(GT), proved on paper) as an
     explicit hypothesis-parameter Prop, and proves the DAG-node glue: bridge + table ⟹ the
     ledger floor for the true measure.  No axioms, no sorries; the bridge is the one named
     remaining analytic node on this branch.

  Division of labor: mac-mini HYP-3857 owns DangerousPatterns (the ℓ¹-simplex enumeration) and
  BonferroniTruncation (the binomial-sign lemma); klein HYP-4000 owns NestDecidable (the
  homogeneous side).  This file is the c-averaged column.

  Cross-validation: `Aledger` agrees with TWO independent Python enumerations (c-side interval
  engine, x-side gap engine) on every row — lrc14_{symbolic_ledger,section_formula_xengine}_kps.out.
-/
import Mathlib.Data.Rat.Floor

namespace LonelyRunner
namespace GapSurplusLedger

/-- The band width `w = 2h` at LRC(14): `w = 1/7`. -/
def w : ℚ := 1/7

/-- Fractional part of a rational. -/
def fracQ (x : ℚ) : ℚ := x - ⌊x⌋

/-- Insertion sort (self-contained; avoids sort-API drift). -/
def insertQ (a : ℚ) : List ℚ → List ℚ
  | [] => [a]
  | b :: l => if a ≤ b then a :: b :: l else b :: insertQ a l

def sortQ : List ℚ → List ℚ
  | [] => []
  | a :: l => insertQ a (sortQ l)

/-- Sorted phase list of the pattern `U` at time `x`: `{frac (u x)}`. -/
def phases (U : List ℤ) (x : ℚ) : List ℚ :=
  sortQ (U.map fun u => fracQ ((u : ℚ) * x))

/-- Circular gaps of a sorted list of points of `[0,1)`. -/
def gapsOf (l : List ℚ) : List ℚ :=
  match l with
  | [] => []
  | a :: _ => List.zipWith (fun q p => q - p) (l.tail ++ [a + 1]) l

/-- The gap surplus `F_U(x) = Σ (g − w)⁺` — the exact `c`-section measure of the free set
(THM-599 section formula). -/
def Fgap (U : List ℤ) (x : ℚ) : ℚ :=
  ((gapsOf (phases U x)).map fun g => max 0 (g - w)).sum

/-- Positive pairwise differences of the pattern. -/
def diffs (U : List ℤ) : List ℤ :=
  (U.flatMap fun a => U.map fun b => b - a).filter (0 < ·)

/-- The breakpoint list of `F_U` on `[0,1]`: for every difference `m`, the collision points `j/m`
and the gap-through-`w` points `(7j ± 1)/(7m)`, plus the endpoints. -/
def breakpoints (U : List ℤ) : List ℚ :=
  let core : List ℚ := (diffs U).flatMap fun (m : ℤ) =>
    (List.range m.toNat).flatMap fun (j : ℕ) =>
      [ ((j : ℤ) : ℚ) / ((m : ℤ) : ℚ),
        ((7*(j : ℤ) + 1 : ℤ) : ℚ) / ((7*m : ℤ) : ℚ),
        (((7*(j : ℤ) - 1).emod (7*m) : ℤ) : ℚ) / ((7*m : ℤ) : ℚ) ]
  sortQ (((0 : ℚ) :: (1 : ℚ) :: core).dedup)

/-- The EXACT ledger value: the midpoint sum over breakpoint intervals.  By THM-599 (SF)+(GT)
(piecewise-affineness of `F_U` between breakpoints), this equals `∫₀¹ F_U`. -/
def Aledger (U : List ℤ) : ℚ :=
  let bps := breakpoints U
  (List.zip bps bps.tail).foldl
    (fun acc p => acc + (p.2 - p.1) * Fgap U ((p.1 + p.2)/2)) 0

/-- The admissible witness floor `m_P = 14249/252252` (THM-530). -/
def witnessMP : ℚ := 14249/252252

/-! ## The ledger rows (the finite table, discharged) -/

theorem ledger_pair   : Aledger [1, 2] = 36/49 := by native_decide
theorem ledger_ap3    : Aledger [1, 2, 3] = 61/98 := by native_decide
theorem ledger_ap3'   : Aledger [10, 24, 38] = 61/98 := by native_decide
theorem ledger_k4     : Aledger [1, 2, 3, 4] = 11/21 := by native_decide
theorem ledger_k5     : Aledger [3, 5, 6, 7, 9] = 127/294 := by native_decide
theorem ledger_k6     : Aledger [2, 4, 6, 7, 8, 10] = 141/392 := by native_decide
theorem ledger_k7     : Aledger [3, 5, 6, 7, 8, 9, 11] = 173/588 := by native_decide
theorem ledger_k8     : Aledger [1, 3, 5, 6, 7, 8, 9, 11] = 169/686 := by native_decide
theorem ledger_k9     : Aledger [1, 3, 4, 5, 6, 7, 8, 9, 11] = 4267/20580 := by native_decide
theorem ledger_k10    : Aledger [1, 3, 4, 5, 6, 7, 8, 9, 11, 13] = 4279/24696 := by native_decide
theorem ledger_k11    : Aledger [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14] = 9451/61740 := by
  native_decide
theorem ledger_k12    :
    Aledger [4, 8, 9, 10, 12, 13, 14, 15, 16, 17, 22, 24] = 4615489/35315280 := by native_decide
theorem ledger_k13    :
    Aledger [3, 7, 8, 9, 11, 13, 15, 16, 17, 19, 23, 24, 27] = 144699147/1267426160 := by
  native_decide

/-- The argmin table (the S30/S32 census minimizer patterns, k = 2..13). -/
def argminTable : List (List ℤ) :=
  [ [1, 2], [1, 2, 3], [1, 2, 3, 4], [3, 5, 6, 7, 9], [2, 4, 6, 7, 8, 10],
    [3, 5, 6, 7, 8, 9, 11], [1, 3, 5, 6, 7, 8, 9, 11], [1, 3, 4, 5, 6, 7, 8, 9, 11],
    [1, 3, 4, 5, 6, 7, 8, 9, 11, 13], [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14],
    [4, 8, 9, 10, 12, 13, 14, 15, 16, 17, 22, 24],
    [3, 7, 8, 9, 11, 13, 15, 16, 17, 19, 23, 24, 27] ]

/-- **The finite-table node, discharged**: every argmin pattern clears the admissible witness
floor.  (The margins run ×13.01 at k = 2 down to ×2.02 at k = 13.) -/
theorem ledger_argmins_clear_MP : ∀ U ∈ argminTable, witnessMP ≤ Aledger U := by native_decide

/-! ## The DAG wiring (no axioms; the bridge is a named hypothesis parameter) -/

/-- **The analytic bridge** (THM-599 (SF) + (GT), proved on paper, not yet mechanized): a
function `A_true` (the genuine c-averaged lonely measure `∫₀¹ F_U`) agrees with the computable
breakpoint sum.  This is the single remaining analytic node of the ledger branch; its paper
proof is one generic lemma ("the integral of a piecewise-affine function is its midpoint sum on
a refining grid") plus the breakpoint-completeness of `breakpoints` ("no collision ⟹ phase
order constant ⟹ affine"). -/
def LedgerBridge (A_true : List ℤ → ℚ) : Prop :=
  ∀ U ∈ argminTable, A_true U = Aledger U

/-- **DAG-node glue**: the bridge plus the discharged finite table give the ledger floor for the
true measure — the (⋆)-census interface consumed by the witness route (HYP-3953 §5). -/
theorem ledger_floor_of_bridge (A_true : List ℤ → ℚ) (hb : LedgerBridge A_true) :
    ∀ U ∈ argminTable, witnessMP ≤ A_true U := by
  intro U hU
  rw [hb U hU]
  exact ledger_argmins_clear_MP U hU

end GapSurplusLedger
end LonelyRunner
