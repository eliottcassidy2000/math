/-
  TournamentH7.LRCBooleanTypeCut -- finite arithmetic for the HYP-2751
  Boolean/type signed aggregate cut.

  The companion Python script validates the cut against the full k=8 bounded
  bank.  This Lean module records the exact rational arithmetic of the three
  active low-depth type atoms:

    T1    = type (1,(1))
    T2sep = type (2,(1,1))
    T2adj = type (2,(2)).

  Positive coefficients make AP a strict minimizer of the low-depth miss mass;
  negative coefficients give the signed aggregate on which AP is maximal.
-/

import Mathlib.Tactic

namespace LonelyRunner
namespace BooleanTypeCut

/-- A row of low-depth type differences, with common denominator `den`.
The represented rational triple is `(t1/den, t2sep/den, t2adj/den)`. -/
structure TypeDiff where
  den : Int
  t1 : Int
  t2sep : Int
  t2adj : Int
deriving DecidableEq, Repr

/-- Integer coefficient triple for `(T1,T2sep,T2adj)`. -/
structure Coeff where
  c1 : Int
  c2sep : Int
  c2adj : Int
deriving DecidableEq, Repr

/-- Numerator of the coefficient/difference pairing. -/
def evalNum (c : Coeff) (d : TypeDiff) : Int :=
  c.c1 * d.t1 + c.c2sep * d.t2sep + c.c2adj * d.t2adj

/-- Sum of the three integer coefficients. -/
def coeffSum (c : Coeff) : Int := c.c1 + c.c2sep + c.c2adj

/-- The sharp L1-normalized three-atom cut from the active LP. -/
def optimalCoeff : Coeff := { c1 := 60601, c2sep := 164633, c2adj := 5815 }

/-- A compact integer cut with nearly the same normalized margin. -/
def compactCoeff : Coeff := { c1 := 21, c2sep := 57, c2adj := 2 }

/-- Active row `(0,1,2,4,5,6,7,10)`, differences
`(27/490, 11/980, 44/735) = (162,33,176)/2940`. -/
def active0 : TypeDiff := { den := 2940, t1 := 162, t2sep := 33, t2adj := 176 }

/-- Active row `(0,2,3,4,5,6,7,8)`, differences
`(25/1176, 149/5880, 19/1470) = (125,149,76)/5880`. -/
def active1 : TypeDiff := { den := 5880, t1 := 125, t2sep := 149, t2adj := 76 }

/-- Active row `(0,4,5,6,8,9,10,14)`, differences
`(157/8820, 23/980, 899/8820) = (157,207,899)/8820`. -/
def active2 : TypeDiff := { den := 8820, t1 := 157, t2sep := 207, t2adj := 899 }

/-- The three LP-active rows. -/
def active (r : Fin 3) : TypeDiff :=
  match r.val with
  | 0 => active0
  | 1 => active1
  | _ => active2

/-- Exact sharp integer margin for `optimalCoeff`, before dividing by coefficient sum:
`60601*T1+164633*T2sep+5815*T2adj >= 2324813/420`. -/
def optimalMarginNum : Int := 2324813
def optimalMarginDen : Int := 420

/-- Exact compact-cut margin:
`21*T1+57*T2sep+2*T2adj >= 8447/4410`. -/
def compactMarginNum : Int := 8447
def compactMarginDen : Int := 4410

/-- The sharp coefficient denominator from the L1-normalized LP. -/
theorem optimalCoeff_sum :
    coeffSum optimalCoeff = 231049 := by
  native_decide

/-- The compact coefficient sum used for normalization. -/
theorem compactCoeff_sum :
    coeffSum compactCoeff = 80 := by
  native_decide

/-- The three active rows exactly saturate the sharp LP cut. -/
theorem optimal_active_rows_equal :
    ∀ r : Fin 3,
      optimalMarginDen * evalNum optimalCoeff (active r) =
        optimalMarginNum * (active r).den := by
  native_decide

/-- The compact cut is valid on the sharp active rows.  Only `active2` is tight;
the full 3430-row validation is carried by the companion exact script. -/
theorem compact_active_rows_above :
    ∀ r : Fin 3,
      compactMarginDen * evalNum compactCoeff (active r) >=
        compactMarginNum * (active r).den := by
  native_decide

/-- The compact cut's minimum witness among the sharp active rows is `active2`. -/
theorem compact_active2_equal :
    compactMarginDen * evalNum compactCoeff active2 =
      compactMarginNum * active2.den := by
  native_decide

/-- The sharp and compact margins are genuinely positive. -/
theorem margins_positive :
    0 < optimalMarginNum ∧ 0 < optimalMarginDen ∧
      0 < compactMarginNum ∧ 0 < compactMarginDen := by
  native_decide

/-! ### Axiom audit -/

#print axioms optimalCoeff_sum
#print axioms compactCoeff_sum
#print axioms optimal_active_rows_equal
#print axioms compact_active_rows_above
#print axioms compact_active2_equal
#print axioms margins_positive

end BooleanTypeCut
end LonelyRunner
