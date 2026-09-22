import Std

/- Closed arithmetic witnesses only. The rational-to-cleared-polynomial
   interpretation is described in README; no global curve or cycle theorem
   is asserted by this package. -/
namespace CatalanEllipticAudit

def aa : Nat :=
  154476802108746166441951315019919837485664325669565431700026634898253202035277999

def bb : Nat :=
  368751317941299998271978115652254748254929799689719709962831374716372246340555790

def cc : Nat :=
  43736126779286972578612526023713901528165375581616136186214379933784234677720360

def fruitNumerator (a b c : Nat) : Nat :=
  a * (a + b) * (a + c) + b * (a + b) * (b + c) + c * (a + c) * (b + c)

def fruitDenominator (a b c : Nat) : Nat :=
  (a + b) * (a + c) * (b + c)

theorem pasted_not_four :
    fruitNumerator aa bb cc ≠ 4 * fruitDenominator aa bb cc := by decide

theorem pasted_denominators_positive :
    0 < aa + bb ∧ 0 < aa + cc ∧ 0 < bb + cc := by decide

theorem repair_divisions_exact :
    10 * (bb / 10) = bb ∧ 10 * (cc / 10) = cc := by decide

theorem repaired_equals_four :
    fruitNumerator aa (bb / 10) (cc / 10) =
      4 * fruitDenominator aa (bb / 10) (cc / 10) := by decide

theorem repaired_denominators_positive :
    0 < aa + bb / 10 ∧ 0 < aa + cc / 10 ∧ 0 < bb / 10 + cc / 10 := by decide

theorem polynomial_point_ten_fortynine :
    (49 : Nat)^2 = 2 * 10^3 + 4 * 10^2 + 1 := by decide

theorem collatz_gap :
    (2 : Int)^11 - 3^7 = -139 := by decide

theorem collatz_ordered_carry :
    (3 : Nat)^6 + 3^5 * 2 + 3^4 * 2^2 + 3^3 * 2^3 +
      3^2 * 2^5 + 3 * 2^6 + 2^7 = 2363 := by decide

theorem collatz_carry_cancellation : (2363 : Nat) = 17 * 139 := by decide

theorem collatz_fixed_point_equation : (-139 : Int) * (-17) = 2363 := by decide

theorem collatz_seven_cleared_steps :
    (3 : Int) * (-17) + 1 = 2 * (-25) ∧
    (3 : Int) * (-25) + 1 = 2 * (-37) ∧
    (3 : Int) * (-37) + 1 = 2 * (-55) ∧
    (3 : Int) * (-55) + 1 = 4 * (-41) ∧
    (3 : Int) * (-41) + 1 = 2 * (-61) ∧
    (3 : Int) * (-61) + 1 = 2 * (-91) ∧
    (3 : Int) * (-91) + 1 = 16 * (-17) := by decide

theorem collatz_successors_odd :
    (-25 : Int) % 2 = 1 ∧ (-37 : Int) % 2 = 1 ∧
    (-55 : Int) % 2 = 1 ∧ (-41 : Int) % 2 = 1 ∧
    (-61 : Int) % 2 = 1 ∧ (-91 : Int) % 2 = 1 ∧
    (-17 : Int) % 2 = 1 := by decide

end CatalanEllipticAudit
