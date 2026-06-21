/-
  TournamentH7.LRCFactorialAtom -- finite factorial atom identities for LRC14.

  This module isolates the algebra behind HYP-2720..HYP-2728.  It does not
  formalize analytic LRC estimates; it proves the seven-coordinate binomial
  identities that make the `Q_0` origin atom an alternating factorial boundary.
-/

namespace LonelyRunner
namespace FactorialAtom

/-- LRC14 has six inner sectors, so the missed-count atom has coordinates 0..6. -/
def atomCount : Nat := 7

/-- A signed missed-count atom vector `q_t`, `t=0..6`. -/
abbrev Atom := Fin atomCount -> Int

/-- Tiny local binomial coefficient, avoiding a Mathlib dependency for this
finite audit module. -/
def choose : Nat -> Nat -> Nat
  | _, 0 => 1
  | 0, _ + 1 => 0
  | n + 1, k + 1 => choose n k + choose n (k + 1)

/-- Natural binomial coefficients, cast to `Int`. -/
def chooseInt (n k : Nat) : Int := (choose n k : Int)

/-- Sum over the seven atom coordinates. -/
def sum7 (f : Nat -> Int) : Int :=
  f 0 + f 1 + f 2 + f 3 + f 4 + f 5 + f 6

/-- The alternating sign `(-1)^k` as an integer. -/
def altSign (k : Nat) : Int := if k % 2 = 0 then 1 else -1

/-- The finite-difference atom packet created by a unit factorial moment `W_j`. -/
def basis (j : Fin atomCount) : Atom :=
  fun t =>
    if t.val <= j.val then altSign (j.val - t.val) * chooseInt j.val t.val else 0

/-- Factorial moment readout `W_j(q)=sum_{t>=j} binom(t,j) q_t`. -/
def moment (q : Atom) (j : Fin atomCount) : Int :=
  sum7 fun t =>
    if h : t < atomCount then
      if j.val <= t then chooseInt t j.val * q ⟨t, h⟩ else 0
    else
      0

/-- The origin atom coordinate. -/
def q0 (q : Atom) : Int := q ⟨0, by decide⟩

/-- Coefficient of atom coordinate `q_t` in the alternating boundary of moments. -/
def originCoeff (t : Fin atomCount) : Int :=
  sum7 fun j =>
    if j <= t.val then altSign j * chooseInt t.val j else 0

/-- The high-tail readout `q5+5q6`. -/
def tail45 (q : Atom) : Int :=
  q ⟨5, by decide⟩ + 5 * q ⟨6, by decide⟩

/-- The Bonferroni4 high-tail readout `U4=q0+q5+5q6`. -/
def U4 (q : Atom) : Int :=
  q ⟨0, by decide⟩ + tail45 q

/-- Low leakage readout `W1+W2` used by HYP-2722/HYP-2726. -/
def low12 (q : Atom) : Int :=
  moment q ⟨1, by decide⟩ + moment q ⟨2, by decide⟩

/-- HYP-2721 cheap directions, scaled to integral atom vectors.

The index `0..4` represents cheap directions `r=1..5`; denominators are cleared
by `cheapScale`. -/
def cheapScaled (r : Fin 5) : Atom :=
  fun t =>
    match r.val, t.val with
    | 0, 0 => 5
    | 0, 1 => -6
    | 0, 6 => 1
    | 1, 0 => 5
    | 1, 1 => -9
    | 1, 3 => 5
    | 1, 6 => -1
    | 2, 0 => 2
    | 2, 1 => -6
    | 2, 2 => 5
    | 2, 5 => -2
    | 2, 6 => 1
    | 3, 0 => 2
    | 3, 1 => -9
    | 3, 2 => 15
    | 3, 3 => -10
    | 3, 5 => 3
    | 3, 6 => -1
    | 4, 0 => 1
    | 4, 1 => -6
    | 4, 2 => 15
    | 4, 3 => -20
    | 4, 4 => 15
    | 4, 5 => -6
    | 4, 6 => 1
    | _, _ => 0

/-- Denominator-clearing scale for `cheapScaled`. -/
def cheapScale (r : Fin 5) : Int :=
  match r.val with
  | 0 => 5
  | 1 => 5
  | 2 => 2
  | 3 => 2
  | _ => 1

/-- Integral tail45 target for the scaled cheap directions. -/
def cheapTailTarget (r : Fin 5) : Int :=
  match r.val with
  | 0 => 5
  | 1 => -5
  | 2 => 3
  | 3 => -2
  | _ => -1

/-- Lower numerator of the generated tail45 strip from HYP-2731. -/
def tailFloorNum : Int := 182

/-- Lower denominator of the generated tail45 strip from HYP-2731. -/
def tailFloorDen : Int := 2005

/-- Upper numerator of the generated tail45 strip from HYP-2731. -/
def tailCeilNum : Int := 10910

/-- Upper denominator of the generated tail45 strip from HYP-2731. -/
def tailCeilDen : Int := 21539

/-- Cross-multiplied statement that a scaled packet is below the tail strip. -/
def belowTailFloor (q : Atom) (scale : Int) : Prop :=
  tail45 q * tailFloorDen < tailFloorNum * scale

/-- Cross-multiplied statement that a scaled packet is above the tail strip. -/
def aboveTailCeil (q : Atom) (scale : Int) : Prop :=
  tailCeilNum * scale < tail45 q * tailCeilDen

/-- Boolean form of the cheap-side tail-strip exclusion check. -/
def outsideTailStripBool (q : Atom) (scale : Int) : Bool :=
  (tail45 q * tailFloorDen < tailFloorNum * scale) ||
    (tailCeilNum * scale < tail45 q * tailCeilDen)

/-- The finite-difference packets are exactly dual to the factorial moments. -/
theorem basis_moment_delta :
    ∀ i j : Fin atomCount, moment (basis j) i = if i = j then 1 else 0 := by
  native_decide

/-- The `Q_0` coordinate is the alternating boundary of factorial moments. -/
theorem originCoeff_delta :
    ∀ t : Fin atomCount, originCoeff t = if t.val = 0 then 1 else 0 := by
  native_decide

/-- On a unit packet, `Q_0(B_j)=(-1)^j`. -/
theorem basis_q0_sign :
    ∀ j : Fin atomCount, q0 (basis j) = altSign j.val := by
  native_decide

/-- `U4` sees packets through degree four and is blind to `B_5,B_6`. -/
theorem U4_basis :
    ∀ j : Fin atomCount, U4 (basis j) = if j.val <= 4 then altSign j.val else 0 := by
  native_decide

/-- The low `W1+W2` leakage sees exactly the `B_1` and `B_2` packet axes. -/
theorem low12_basis :
    ∀ j : Fin atomCount, low12 (basis j) =
      if j.val = 1 ∨ j.val = 2 then 1 else 0 := by
  native_decide

/-- The scaled cheap directions have the declared positive origin scale. -/
theorem cheapScaled_q0 :
    ∀ r : Fin 5, q0 (cheapScaled r) = cheapScale r := by
  native_decide

/-- The tail45 covector separates the scaled cheap directions into
`5,-5,3,-2,-1`; after normalization these are `1,-1,3/2,-1,-1`. -/
theorem cheapScaled_tail45 :
    ∀ r : Fin 5, tail45 (cheapScaled r) = cheapTailTarget r := by
  native_decide

/-! ### The Delsarte / moment-LP bound (THM-534, HYP-2726)

The LRC(14) cover `p0 = q0` is bounded per-shape by a **Delsarte linear program**:
`q0 <= L_y` whenever the dual `g(t) = sum_r y_r * C(t,r)` dominates `[t=0]` and the
missed-count atom `q` is a genuine (nonnegative) distribution.  Here we formalize the
`k=8` case, whose THM-534 dual `y = (1,-1,1,-9/10,3/5,0,0)` cleared by `*10` gives the
integer covector `yK8 = (10,-10,10,-9,6,0,0)`.  Its dual readout collapses to the clean
vector `g = (10,0,0,1,0,0,10)`, i.e. `L_y(q) = 10*q0 + q3 + 10*q6`. -/

/-- Scaled (`*10`) integer moment-LP dual coefficients for `k=8` (THM-534). -/
def yK8 : Fin atomCount -> Int := fun r =>
  match r.val with
  | 0 => 10
  | 1 => -10
  | 2 => 10
  | 3 => -9
  | 4 => 6
  | _ => 0

/-- The scaled Delsarte functional `L_y(q) = sum_r y_r * W_r(q)` (`= 10 *` THM-534's `L_y`). -/
def LyK8 (q : Atom) : Int :=
  (yK8 ⟨0, by decide⟩) * moment q ⟨0, by decide⟩
  + (yK8 ⟨1, by decide⟩) * moment q ⟨1, by decide⟩
  + (yK8 ⟨2, by decide⟩) * moment q ⟨2, by decide⟩
  + (yK8 ⟨3, by decide⟩) * moment q ⟨3, by decide⟩
  + (yK8 ⟨4, by decide⟩) * moment q ⟨4, by decide⟩
  + (yK8 ⟨5, by decide⟩) * moment q ⟨5, by decide⟩
  + (yK8 ⟨6, by decide⟩) * moment q ⟨6, by decide⟩

/-- The dual readout `g(t) = sum_{r<=t} y_r * C(t,r)`; for the `k=8` dual it is the
Krawtchouk-nonnegative vector `(10,0,0,1,0,0,10)` (HYP-2726a). -/
def gK8 (t : Fin atomCount) : Int :=
  sum7 fun r =>
    if h : r < atomCount then
      (if r <= t.val then yK8 ⟨r, h⟩ * chooseInt t.val r else 0)
    else 0

/-- The dual `g` is exactly `(10,0,0,1,0,0,10)` -- Delsarte dual feasibility,
Krawtchouk-nonnegative, dominating `10*[t=0]`. -/
theorem gK8_values :
    ∀ t : Fin atomCount,
      gK8 t = (if t.val = 0 then 10 else if t.val = 3 then 1 else if t.val = 6 then 10 else 0) := by
  native_decide

/-- Delsarte dual feasibility: `g(t) >= 10*[t=0]` for every depth `t`. -/
theorem gK8_dominates :
    ∀ t : Fin atomCount, (if t.val = 0 then (10 : Int) else 0) <= gK8 t := by
  native_decide

/-- **The Delsarte readout identity.**  The scaled moment-LP functional collapses to
the dual covector: `L_y(q) = 10*q0 + q3 + 10*q6`. -/
theorem LyK8_readout (q : Atom) :
    LyK8 q = 10 * q ⟨0, by decide⟩ + q ⟨3, by decide⟩ + 10 * q ⟨6, by decide⟩ := by
  simp only [LyK8, moment, sum7, atomCount, yK8, chooseInt, choose,
    Nat.reduceLT, Nat.reduceLeDiff, reduceDIte, reduceIte, Nat.reduceAdd]
  omega

/-- **The per-shape Delsarte / moment-LP bound (THM-534, HYP-2726).**  For any
nonnegative missed-count atom `q` (a genuine row distribution), the origin/cover atom
satisfies `10 * q0 <= L_y(q)`. -/
theorem delsarte_bound_k8 (q : Atom) (hq : ∀ t : Fin atomCount, 0 <= q t) :
    10 * q0 q <= LyK8 q := by
  rw [LyK8_readout]
  have h3 : 0 <= q ⟨3, by decide⟩ := hq _
  have h6 : 0 <= q ⟨6, by decide⟩ := hq _
  unfold q0
  omega

/-! #### The other binding rows: `k=9,10` and `k=11,12,13`

Same Delsarte structure at every binding row, with THM-534 duals
`y(9,10)=(1,-13/18,4/9,-1/6,0,0,0)` (×18 → `(18,-13,8,-3,0,0,0)`) and
`y(11,12,13)=(1,-1/2,1/6,0,0,0,0)` (×6 → `(6,-3,1,0,0,0,0)`).  Readouts:
`L_y(9)=18*q0+5*q1+2*q4+3*q5`, `L_y(11)=6*q0+3*q1+q2+q5+3*q6`. -/

/-- ×18 integer moment-LP dual for `k=9,10` (THM-534). -/
def yK9 : Fin atomCount -> Int := fun r =>
  match r.val with
  | 0 => 18 | 1 => -13 | 2 => 8 | 3 => -3 | _ => 0

/-- Scaled Delsarte functional for `k=9,10`. -/
def LyK9 (q : Atom) : Int :=
  (yK9 ⟨0, by decide⟩) * moment q ⟨0, by decide⟩
  + (yK9 ⟨1, by decide⟩) * moment q ⟨1, by decide⟩
  + (yK9 ⟨2, by decide⟩) * moment q ⟨2, by decide⟩
  + (yK9 ⟨3, by decide⟩) * moment q ⟨3, by decide⟩
  + (yK9 ⟨4, by decide⟩) * moment q ⟨4, by decide⟩
  + (yK9 ⟨5, by decide⟩) * moment q ⟨5, by decide⟩
  + (yK9 ⟨6, by decide⟩) * moment q ⟨6, by decide⟩

theorem LyK9_readout (q : Atom) :
    LyK9 q = 18 * q ⟨0, by decide⟩ + 5 * q ⟨1, by decide⟩
      + 2 * q ⟨4, by decide⟩ + 3 * q ⟨5, by decide⟩ := by
  simp only [LyK9, moment, sum7, atomCount, yK9, chooseInt, choose,
    Nat.reduceLT, Nat.reduceLeDiff, reduceDIte, reduceIte, Nat.reduceAdd, Nat.reduceMul]
  omega

/-- Delsarte / moment-LP bound at `k=9,10`: `18 * q0 <= L_y(q)` for nonnegative `q`. -/
theorem delsarte_bound_k9 (q : Atom) (hq : ∀ t : Fin atomCount, 0 <= q t) :
    18 * q0 q <= LyK9 q := by
  rw [LyK9_readout]
  have h1 : 0 <= q ⟨1, by decide⟩ := hq _
  have h4 : 0 <= q ⟨4, by decide⟩ := hq _
  have h5 : 0 <= q ⟨5, by decide⟩ := hq _
  unfold q0
  omega

/-- ×6 integer moment-LP dual for `k=11,12,13` (THM-534). -/
def yK11 : Fin atomCount -> Int := fun r =>
  match r.val with
  | 0 => 6 | 1 => -3 | 2 => 1 | _ => 0

/-- Scaled Delsarte functional for `k=11,12,13`. -/
def LyK11 (q : Atom) : Int :=
  (yK11 ⟨0, by decide⟩) * moment q ⟨0, by decide⟩
  + (yK11 ⟨1, by decide⟩) * moment q ⟨1, by decide⟩
  + (yK11 ⟨2, by decide⟩) * moment q ⟨2, by decide⟩
  + (yK11 ⟨3, by decide⟩) * moment q ⟨3, by decide⟩
  + (yK11 ⟨4, by decide⟩) * moment q ⟨4, by decide⟩
  + (yK11 ⟨5, by decide⟩) * moment q ⟨5, by decide⟩
  + (yK11 ⟨6, by decide⟩) * moment q ⟨6, by decide⟩

theorem LyK11_readout (q : Atom) :
    LyK11 q = 6 * q ⟨0, by decide⟩ + 3 * q ⟨1, by decide⟩
      + q ⟨2, by decide⟩ + q ⟨5, by decide⟩ + 3 * q ⟨6, by decide⟩ := by
  simp only [LyK11, moment, sum7, atomCount, yK11, chooseInt, choose,
    Nat.reduceLT, Nat.reduceLeDiff, reduceDIte, reduceIte, Nat.reduceAdd, Nat.reduceMul]
  omega

/-- Delsarte / moment-LP bound at `k=11,12,13`: `6 * q0 <= L_y(q)` for nonnegative `q`. -/
theorem delsarte_bound_k11 (q : Atom) (hq : ∀ t : Fin atomCount, 0 <= q t) :
    6 * q0 q <= LyK11 q := by
  rw [LyK11_readout]
  have h1 : 0 <= q ⟨1, by decide⟩ := hq _
  have h2 : 0 <= q ⟨2, by decide⟩ := hq _
  have h5 : 0 <= q ⟨5, by decide⟩ := hq _
  have h6 : 0 <= q ⟨6, by decide⟩ := hq _
  unfold q0
  omega

/-- Every scaled cheap direction is outside the generated tail45 strip.

This is a Boolean audit theorem so the module remains self-contained without
extra order/typeclass imports. -/
theorem cheapScaled_outside_tailStrip_bool :
    ∀ r : Fin 5, outsideTailStripBool (cheapScaled r) (cheapScale r) = true := by
  native_decide

/-! ### Axiom audit

These are closed finite identities over seven coordinates.  In a working Lean
environment, the axiom output should be empty or contain only Lean foundations.
-/

#print axioms basis_moment_delta
#print axioms originCoeff_delta
#print axioms basis_q0_sign
#print axioms U4_basis
#print axioms low12_basis
#print axioms cheapScaled_q0
#print axioms cheapScaled_tail45
#print axioms gK8_values
#print axioms gK8_dominates
#print axioms LyK8_readout
#print axioms delsarte_bound_k8
#print axioms cheapScaled_outside_tailStrip_bool

end FactorialAtom
end LonelyRunner
