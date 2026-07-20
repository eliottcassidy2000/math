import Mathlib

/-!
# Pure-K fastest-wall bank quotient (THM-1277)

This module checks the arithmetic and finite combinatorial consumers of the
paper theorem.  The geometric providers are THM-1273's choice of an
`h`-active point more than `1/12` from the K/E interface, the alternating
fastest wall word, deletion-minimal interval-chain adjacency, and the fact
that every wall strictly inside K has a selected middle-owner crosser.

The module verifies the sharp wall-length conversion, the integer threshold
interface, the prefix/suffix flood layering rule, the selected/flood bank
coefficients, the distinct-owner reciprocal gain, the localized regular-run
capacity, the centered wall-residue threshold, and the exact c=140 guardrail.
It contains no covering axiom, `sorry`, `admit`, or `native_decide`.
-/

namespace LRC14.PureKFastestWallBankQuotient

/-- If an interval deeper than `1/12` contains only `P` complete fastest
teeth after its initial active tooth, the alternating `1:6` wall word gives
the exact cross-multiplied upper bound `h < (14P+16)d1`. -/
theorem wallBank_length_consumer
    (P : ℕ) (d₁ h needleLength : ℚ)
    (hh : 0 < h)
    (hdeep : 1 / 12 < needleLength)
    (hcell : needleLength ≤
      (7 * (P : ℚ) + 8) * (d₁ / (6 * h))) :
    h < (14 * (P : ℚ) + 16) * d₁ := by
  have hrearrange :
      (7 * (P : ℚ) + 8) * (d₁ / (6 * h)) =
        (((7 * (P : ℚ) + 8) * d₁) / 6) / h := by
    field_simp [ne_of_gt hh]
    ring
  rw [hrearrange] at hcell
  have hmul : needleLength * h ≤ ((7 * (P : ℚ) + 8) * d₁) / 6 :=
    (le_div_iff₀ hh).mp hcell
  have hstrict : h / 12 < needleLength * h := by
    nlinarith
  nlinarith

/-- Arithmetic interface for
`P >= floor((h-2d1)/(14d1))`: every integer threshold Q whose left endpoint
has already been reached is at most the number P of complete teeth. -/
theorem completeTooth_threshold
    (P Q : ℕ) (d₁ h : ℚ)
    (hd₁ : 0 < d₁)
    (hupper : h < (14 * (P : ℚ) + 16) * d₁)
    (hlower : (14 * (Q : ℚ) + 2) * d₁ ≤ h) :
    Q ≤ P := by
  by_contra hnot
  have hnat : P + 1 ≤ Q := by omega
  have hcast : (P : ℚ) + 1 ≤ (Q : ℚ) := by exact_mod_cast hnat
  nlinarith

/-- At `h >= 2d1`, the initial fastest-tooth remainder is no longer than
`1/12`, so a needle deeper than `1/12` crosses its first right wall. -/
theorem firstWall_forced
    (d₁ h firstDistance needleLength : ℚ)
    (hd₁ : 0 ≤ d₁) (hh : 0 < h)
    (hratio : 2 * d₁ ≤ h)
    (hfirst : firstDistance ≤ d₁ / (6 * h))
    (hdeep : 1 / 12 < needleLength) :
    firstDistance < needleLength := by
  have ha : d₁ / (6 * h) ≤ 1 / 12 := by
    rw [div_le_iff₀ (by positivity : (0 : ℚ) < 6 * h)]
    nlinarith
  linarith

/-- The pointwise extension of THM-1275: an arbitrary unselected complete
fastest tooth, including a prefix/suffix tooth, layers with a selected seam.
On their intersection the full comb has the selected seam's two low layers
plus the unselected fastest layer. -/
theorem prefixSuffixFlood_layeredMultiplicity
    (flood seam multiplicity : ℕ)
    (hLayers : 1 + flood + seam ≤ multiplicity) :
    flood + seam ≤ multiplicity - 1 := by
  omega

/-- Slow-gap normalization and the six-bin density floor turn one complete
fastest tooth into the weighted credit `c/(8h)`. -/
theorem wholeFlood_weight (c h : ℚ) (hh : h ≠ 0) :
    (7 * c / 6) * (3 / 4) * (1 / (7 * h)) = c / (8 * h) := by
  field_simp [hh]
  ring

/-- One primitive `1/(14hd5)` boundary quantum has weighted credit
`c/(16hd5)`. -/
theorem boundarySeam_weight (c h d₅ : ℚ) (hh : h ≠ 0) (hd₅ : d₅ ≠ 0) :
    (7 * c / 6) * (3 / 4) * (1 / (14 * h * d₅)) =
      c / (16 * h * d₅) := by
  field_simp [hh, hd₅]
  ring

/-- The singleton-discrepancy factor `49/6` turns the same boundary quantum
into the harmonic coefficient `7/(12hd5)`. -/
theorem boundarySeam_harmonic (h d₅ : ℚ) (hh : h ≠ 0) (hd₅ : d₅ ≠ 0) :
    (49 / 6) * (1 / (14 * h * d₅)) = 7 / (12 * h * d₅) := by
  field_simp [hh, hd₅]
  ring

/-- A whole flood dominates the phase-free two-distinct-owner seam floor.
This is the algebra which collapses `F+A=P` to one uniform P-bank. -/
theorem selectedFlood_bankCollapse
    (scale base pairFloor F A P : ℚ)
    (hscale : 0 ≤ scale) (hF : 0 ≤ F)
    (hpair : pairFloor ≤ 2)
    (hP : P = F + A) :
    base + scale * (2 * F + A * pairFloor) ≥
      base + scale * P * pairFloor := by
  have hnonneg : 0 ≤ scale * F * (2 - pairFloor) := by positivity
  rw [hP]
  nlinarith

/-- If the two lower boundary owners are distinct, the smaller is at most
`d4` and the larger at most `d5`; their reciprocal product clocks therefore
pay `(1/h)(1/d4+1/d5)`. -/
theorem distinctOwner_pairFloor
    (h jSmall jLarge d₄ d₅ : ℚ)
    (hh : 0 < h) (hjs : 0 < jSmall) (hjl : 0 < jLarge)
    (hsmall : jSmall ≤ d₄) (hlarge : jLarge ≤ d₅) :
    1 / (h * d₄) + 1 / (h * d₅) ≤
      1 / (h * jSmall) + 1 / (h * jLarge) := by
  have hpSmall : h * jSmall ≤ h * d₄ :=
    mul_le_mul_of_nonneg_left hsmall (le_of_lt hh)
  have hpLarge : h * jLarge ≤ h * d₅ :=
    mul_le_mul_of_nonneg_left hlarge (le_of_lt hh)
  have hrecSmall : 1 / (h * d₄) ≤ 1 / (h * jSmall) :=
    one_div_le_one_div_of_le (mul_pos hh hjs) hpSmall
  have hrecLarge : 1 / (h * d₅) ≤ 1 / (h * jLarge) :=
    one_div_le_one_div_of_le (mul_pos hh hjl) hpLarge
  linarith

/-- Distinct teeth of one slower comb are separated too far to meet the two
walls of a single fastest tooth.  Combined with deletion minimality, this is
why a selected complete fastest tooth has distinct boundary owners. -/
theorem sameOwner_distinctTeeth_tooFar
    (h j : ℚ) (hh : 0 < h) (hj : 0 < j) (hjh : j < h) :
    1 / (7 * h) < 6 / (7 * j) := by
  rw [div_lt_div_iff₀ (by positivity : (0 : ℚ) < 7 * h)
    (by positivity : (0 : ℚ) < 7 * j)]
  nlinarith

/-- `F` missing complete teeth and `T` multi-low turns split the selected
vertices into at most `F+T+1` regular runs.  If each has at most `e+1`
vertices, the localized block has the exact capacity used in THM-1277. -/
theorem localizedRun_capacity
    (P A F T e : ℕ)
    (hP : P = A + F)
    (hRuns : A + 1 ≤ (e + 1) * (F + T + 1)) :
    P ≤ (e + 2) * F + (e + 1) * T + e := by
  have hadd := Nat.add_le_add_right hRuns F
  calc
    P + 1 = (A + 1) + F := by omega
    _ ≤ (e + 1) * (F + T + 1) + F := hadd
    _ = ((e + 2) * F + (e + 1) * T + e) + 1 := by ring
  omega

def centeredResidue (x modulus : ℕ) : ℕ :=
  min (x % modulus) (modulus - x % modulus)

def rationalResidueDistance (x modulus : ℕ) : ℚ :=
  (centeredResidue x modulus : ℚ) / modulus

/-- Clearing the denominator at a fastest wall `z=a/(14h)` gives the exact
finite crosser test `|ja|_(14h)<h`.  The geometric equality between circle
distance at z and `rationalResidueDistance (ja) (14h)` is the paper input;
this theorem checks the strict threshold conversion. -/
theorem wallResidue_threshold
    (j a h : ℕ) (hh : 0 < h) :
    rationalResidueDistance (j * a) (14 * h) < 1 / 14 ↔
      centeredResidue (j * a) (14 * h) < h := by
  rw [rationalResidueDistance]
  have hm : (0 : ℚ) < (14 * h : ℕ) := by positivity
  rw [div_lt_iff₀ hm]
  norm_num
  exact_mod_cast (show centeredResidue (j * a) (14 * h) < h ↔
    centeredResidue (j * a) (14 * h) < h from Iff.rfl)

theorem sharp140_wall_order :
    (7476011 : ℚ) / 12938240 < 14603 / 25270 ∧
      (14603 : ℚ) / 25270 < 1133 / 1960 ∧
      (1133 : ℚ) / 1960 < 2923 / 5054 ∧
      (2923 : ℚ) / 5054 < 14617 / 25270 ∧
      (14617 : ℚ) / 25270 < 7425603 / 12837160 := by
  norm_num

theorem sharp140_wall_addresses :
    (14603 : ℚ) / 25270 = (14 * 1043 + 1) / (14 * 1805) ∧
      (2923 : ℚ) / 5054 = (14 * 1044 - 1) / (14 * 1805) ∧
      (14617 : ℚ) / 25270 = (14 * 1044 + 1) / (14 * 1805) := by
  norm_num

theorem sharp140_bridge_contains_first_two_walls :
    (2071 : ℚ) / 3584 < 14603 / 25270 ∧
      (14603 : ℚ) / 25270 < 2923 / 5054 ∧
      (2923 : ℚ) / 5054 < 2073 / 3584 := by
  norm_num

theorem sharp140_ratio_floor_cell :
    (0 : ℚ) ≤ (1805 - 2 * 254) / (14 * 254) ∧
      (1805 - 2 * 254 : ℚ) / (14 * 254) < 1 := by
  norm_num

theorem sharp140_residue_labels :
    centeredResidue (256 * 14603) (14 * 1805) = 1592 ∧
      centeredResidue (256 * 14615) (14 * 1805) = 1480 ∧
      centeredResidue (256 * 14617) (14 * 1805) = 1992 ∧
      1592 < 1805 ∧ 1480 < 1805 ∧ 1805 ≤ 1992 := by
  norm_num [centeredResidue]

#print axioms wallBank_length_consumer
#print axioms completeTooth_threshold
#print axioms firstWall_forced
#print axioms prefixSuffixFlood_layeredMultiplicity
#print axioms wholeFlood_weight
#print axioms boundarySeam_weight
#print axioms boundarySeam_harmonic
#print axioms selectedFlood_bankCollapse
#print axioms distinctOwner_pairFloor
#print axioms sameOwner_distinctTeeth_tooFar
#print axioms localizedRun_capacity
#print axioms wallResidue_threshold
#print axioms sharp140_wall_order
#print axioms sharp140_wall_addresses
#print axioms sharp140_bridge_contains_first_two_walls
#print axioms sharp140_ratio_floor_cell
#print axioms sharp140_residue_labels

end LRC14.PureKFastestWallBankQuotient
