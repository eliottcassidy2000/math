/-
  TournamentH7.LRCL7Discrepancy -- the finite integer core of the LRC(14) L7
  torus-line cell-discrepancy (kind-pasteur HYP-2730/2733/2736).

  The balanced two-far cover correction R(p/q) is bounded by the L1 cell-discrepancy
  D_{p,q} of the (q,p) torus geodesic against the 7x7 sector grid.  On the common
  7*p*q breakpoint grid the occupancy is an integer 7x7 matrix `c` whose row sums and
  column sums all equal `S = p*q` (the uniform 1D marginals, PROVED elementarily).
  Writing the integer discrepancy `Ddef(c,S) = sum_{i,j} |7*c i j - S|` (a normalization
  of `D_{p,q}` vanishing on the same locus), the apex-prime law HYP-2733 says
  `Ddef = 0 <=> 7 | p*q`.

  MATHLIB-FREE (`Int.natAbs` + `omega` + `native_decide`).  We formalize:
   * `cell_apex_necessity` : a balanced cell with zero deviation forces `7 | S`
       (the scalar core of HYP-2733's easy half; cellwise => the matrix law);
   * `cell_uniform_zero`   : the uniform value has zero cell deviation;
   * concrete `native_decide` instances on real atlas ratios: the marginal balance
       (row/col sums = p*q) and the apex law `Ddef = 0 <=> 7 | p*q` for (q,p) = (7,8)
       [apex-aligned, 7|q, Ddef=0] and (2,3),(1,2) [non-apex, Ddef>0].
  The 49-cell universally-quantified omega is intentionally avoided (too heavy);
  the general matrix law is the scalar lemma applied cellwise.
-/

namespace LonelyRunner
namespace L7Discrepancy

/-- A 7x7 integer matrix as a flat row-major lookup. -/
abbrev Mat := Fin 7 → Fin 7 → Int

/-- Build a matrix from a flat row-major list of 49 entries. -/
def ofList (xs : List Int) : Mat := fun i j => xs.getD (i.val * 7 + j.val) 0

/-- Explicit sum over the 49 cells. -/
def sum49 (f : Fin 7 → Fin 7 → Int) : Int :=
  (f 0 0 + f 0 1 + f 0 2 + f 0 3 + f 0 4 + f 0 5 + f 0 6) +
  (f 1 0 + f 1 1 + f 1 2 + f 1 3 + f 1 4 + f 1 5 + f 1 6) +
  (f 2 0 + f 2 1 + f 2 2 + f 2 3 + f 2 4 + f 2 5 + f 2 6) +
  (f 3 0 + f 3 1 + f 3 2 + f 3 3 + f 3 4 + f 3 5 + f 3 6) +
  (f 4 0 + f 4 1 + f 4 2 + f 4 3 + f 4 4 + f 4 5 + f 4 6) +
  (f 5 0 + f 5 1 + f 5 2 + f 5 3 + f 5 4 + f 5 5 + f 5 6) +
  (f 6 0 + f 6 1 + f 6 2 + f 6 3 + f 6 4 + f 6 5 + f 6 6)

/-- Absolute value as an `Int`, via core `Int.natAbs`. -/
def iabs (x : Int) : Int := (x.natAbs : Int)

/-- The integer discrepancy `Ddef(c,S) = sum_{i,j} |7*c i j - S|`. -/
def Ddef (c : Mat) (S : Int) : Int := sum49 (fun i j => iabs (7 * c i j - S))

/-- The `i`-th row sum. -/
def rowSum (c : Mat) (i : Fin 7) : Int :=
  c i 0 + c i 1 + c i 2 + c i 3 + c i 4 + c i 5 + c i 6

/-- The `j`-th column sum. -/
def colSum (c : Mat) (j : Fin 7) : Int :=
  c 0 j + c 1 j + c 2 j + c 3 j + c 4 j + c 5 j + c 6 j

/-! ### The scalar apex-prime core (general, omega) -/

/-- **Apex necessity (scalar core of HYP-2733).** A balanced cell with zero deviation
    `7*c - S = 0` forces `7 | S`.  Applied cellwise this gives the matrix law
    `Ddef = 0 -> 7 | S` (a zero sum of nonnegative `|7*c i j - S|` forces each to vanish). -/
theorem cell_apex_necessity (c S : Int) (h : iabs (7 * c - S) = 0) : (7 : Int) ∣ S := by
  simp only [iabs] at h
  omega

/-- **Apex sufficiency (uniform value).** The uniform occupancy value has zero cell
    deviation. -/
theorem cell_uniform_zero (m : Int) : iabs (7 * m - 7 * m) = 0 := by
  simp only [iabs]; omega

/-! ### Concrete atlas instances (native_decide) -/

/-- Occupancy of the `(q,p)=(2,3)` geodesic (ratio 3/2, the `sup D*q` ratio); `S=pq=6`. -/
def c32 : Mat := ofList
  [2,1,0,1,2,0,0, 0,1,2,0,0,2,1, 2,0,0,2,1,0,1, 0,2,1,0,1,2,0,
   1,0,1,2,0,0,2, 1,2,0,0,2,1,0, 0,0,2,1,0,1,2]

/-- Occupancy of the apex-aligned `(q,p)=(7,8)` geodesic (`7 | q`); uniform, `S=pq=56`. -/
def c87 : Mat := ofList (List.replicate 49 8)

/-- Occupancy of `(q,p)=(1,2)` (ratio 2/1, the `sup D*p` ratio); `S=pq=2`. -/
def c21 : Mat := ofList
  [1,1,0,0,0,0,0, 0,0,1,1,0,0,0, 0,0,0,0,1,1,0, 1,0,0,0,0,0,1,
   0,1,1,0,0,0,0, 0,0,0,1,1,0,0, 0,0,0,0,0,1,1]

/-- Marginal balance: every row/column sum of `c32` equals `S = p*q = 6`. -/
theorem c32_rowSums : ∀ i, rowSum c32 i = 6 := by decide
theorem c32_colSums : ∀ j, colSum c32 j = 6 := by decide

/-- **Apex law, non-apex instance:** `7 ∤ p*q = 6`, discrepancy positive. -/
theorem c32_disc : Ddef c32 6 = 252 := by native_decide

/-- **Apex law, apex-aligned instance:** `7 | p*q = 56` (since `7|q`), discrepancy `= 0`. -/
theorem c87_disc_zero : Ddef c87 56 = 0 := by native_decide
theorem c87_rowSums : ∀ i, rowSum c87 i = 56 := by decide

/-- The `2/1` ratio: balanced (`S=2`), non-apex, discrepancy `= 140`. -/
theorem c21_disc : Ddef c21 2 = 140 := by native_decide
theorem c21_rowSums : ∀ i, rowSum c21 i = 2 := by decide

end L7Discrepancy
end LonelyRunner
