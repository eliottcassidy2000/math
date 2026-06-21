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

/-- `iabs` is nonnegative. -/
theorem iabs_nonneg (x : Int) : 0 ≤ iabs x := by simp only [iabs]; omega

/-- The `(0,0)` summand is at most the whole sum, when all summands are nonnegative.
    Proved by linear arithmetic treating each `g i j` as an opaque nonnegative integer
    (no `natAbs` here, so `omega` stays tractable over the 49 terms). -/
theorem cell_le_sum49 (g : Fin 7 → Fin 7 → Int) (hg : ∀ i j, 0 ≤ g i j) :
    g 0 0 ≤ sum49 g := by
  have h00 := hg 0 0; have h01 := hg 0 1; have h02 := hg 0 2; have h03 := hg 0 3
  have h04 := hg 0 4; have h05 := hg 0 5; have h06 := hg 0 6
  have h10 := hg 1 0; have h11 := hg 1 1; have h12 := hg 1 2; have h13 := hg 1 3
  have h14 := hg 1 4; have h15 := hg 1 5; have h16 := hg 1 6
  have h20 := hg 2 0; have h21 := hg 2 1; have h22 := hg 2 2; have h23 := hg 2 3
  have h24 := hg 2 4; have h25 := hg 2 5; have h26 := hg 2 6
  have h30 := hg 3 0; have h31 := hg 3 1; have h32 := hg 3 2; have h33 := hg 3 3
  have h34 := hg 3 4; have h35 := hg 3 5; have h36 := hg 3 6
  have h40 := hg 4 0; have h41 := hg 4 1; have h42 := hg 4 2; have h43 := hg 4 3
  have h44 := hg 4 4; have h45 := hg 4 5; have h46 := hg 4 6
  have h50 := hg 5 0; have h51 := hg 5 1; have h52 := hg 5 2; have h53 := hg 5 3
  have h54 := hg 5 4; have h55 := hg 5 5; have h56 := hg 5 6
  have h60 := hg 6 0; have h61 := hg 6 1; have h62 := hg 6 2; have h63 := hg 6 3
  have h64 := hg 6 4; have h65 := hg 6 5; have h66 := hg 6 6
  simp only [sum49]; omega

/-- **The general matrix apex law (HYP-2733, easy half), now PROVED.** If the integer
    discrepancy of any occupancy matrix vanishes then `7 | S`.  No 49-term `natAbs`
    blowup: bound the `(0,0)` summand by the (zero) total, force it to vanish, then
    apply the scalar `cell_apex_necessity`. -/
theorem matrix_apex_necessity (c : Mat) (S : Int) (h : Ddef c S = 0) : (7 : Int) ∣ S := by
  have hle : iabs (7 * c 0 0 - S) ≤ Ddef c S :=
    cell_le_sum49 (fun i j => iabs (7 * c i j - S)) (fun i j => iabs_nonneg _)
  have hz : iabs (7 * c 0 0 - S) = 0 := by
    have := iabs_nonneg (7 * c 0 0 - S); omega
  exact cell_apex_necessity (c 0 0) S hz

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

/-! ### The sharp residue closed form (HYP-2739)

The thread-C result: the L1 cell-discrepancy is RESIDUE-ONLY,
`D_{p,q} = 4 * f(‖p‖₇, ‖q‖₇) / (7*p*q)`, where `‖x‖₇ = min(x%7, 7-x%7) ∈ {0,1,2,3}` and
`f a b = 0` if `a*b=0`, `= a*b+3` if `a≠b`, `= a*b+4-2|a-2|` if `a=b`.  In the integer
normalization here (`Ddef = Σ|7c-S| = 7 * S`, since `49c-7pq = 7(7c-S)`), this reads
`Ddef = 7 * Sres p q`.  We machine-check it on the two sharp faces and the apex case. -/

/-- `‖x‖₇ = min(x mod 7, 7 - x mod 7) ∈ {0,1,2,3}`. -/
def norm7 (x : Int) : Int := min (x % 7) (7 - x % 7)

/-- The residue function `f`. -/
def fres (a b : Int) : Int :=
  if a * b = 0 then 0 else if a = b then a * b + 4 - 2 * iabs (a - 2) else a * b + 3

/-- The residue invariant `Sres p q = 4 * f(‖p‖₇, ‖q‖₇)`; `D_{p,q} = Sres/(7pq)` and
    `Ddef = 7 * Sres`. -/
def Sres (p q : Int) : Int := 4 * fres (norm7 p) (norm7 q)

/-- **Sharp residue closed form (HYP-2739), apex face `(q,p)=(7,8)`:** `7 | q`, so
    `Sres = 0` and the discrepancy vanishes. -/
theorem c87_residue : Ddef c87 56 = 7 * Sres 8 7 := by native_decide

/-- **Sharp residue closed form, the `sup D*q = 12/7` face `(q,p)=(2,3)`:** `Sres 3 2 = 36`,
    `Ddef = 7*36 = 252`. -/
theorem c32_residue : Ddef c32 6 = 7 * Sres 3 2 := by native_decide

/-- **Sharp residue closed form, the `sup D*p = 20/7` face `(q,p)=(1,2)`:** `Sres 2 1 = 20`,
    `Ddef = 7*20 = 140`. -/
theorem c21_residue : Ddef c21 2 = 7 * Sres 2 1 := by native_decide

/-- The residue invariant is bounded by `44` on the four residue classes `{0,1,2,3}^2`
    (max `f(3,3)=11`), so `D_{p,q} ≤ 44/(7pq)` (HYP-2739).  Checked on the residue table
    by `decide`; the reduction `norm7 x ∈ {0,1,2,3}` for general `x` is left to prose. -/
theorem Sres_table_le_44 : ∀ a : Fin 4, ∀ b : Fin 4,
    4 * fres (a.val : Int) (b.val : Int) ≤ 44 := by decide

end L7Discrepancy
end LonelyRunner
