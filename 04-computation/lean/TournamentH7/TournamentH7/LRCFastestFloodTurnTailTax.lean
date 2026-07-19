import Mathlib

/-!
# Fastest-owner flood/turn tail-tax arithmetic (THM-1275)

This module checks the arithmetic consumers of the paper theorem: the
regular-run/exception count, the private-occurrence product inequality, the
two-layer multiplicity rule, all normalization coefficients, the `7/2`
dominated-turn simplification, one literal return constant, and the six exact
protrusion values of the `c=140` local countermodel.

Deletion-minimal interval-chain extraction, the assertion that skipped
fastest teeth lie inside the carrier gap, disjointness of the physical invoice
regions, and the THM-1264 return identity are explicit paper providers.  No
claim about LRC(14) is made here.
-/

namespace LRCFastestFloodTurnTailTax

/-- If `R` regular transitions are split by `E` exceptions into at most
`E+1` runs of length at most `e`, then `K` fastest occurrences obey the
sharp block bound `K <= (e+1)(E+1)`. -/
theorem regularRun_exception_bound
    (K R E e : ℕ)
    (hK : K = R + E + 1)
    (hR : R ≤ e * (E + 1)) :
    K ≤ (e + 1) * (E + 1) := by
  calc
    K = R + E + 1 := hK
    _ ≤ e * (E + 1) + E + 1 := by omega
    _ = (e + 1) * (E + 1) := by ring

/-- Address skips `X` and consecutive-address multi-low turns `T` dominate
the number `E` of nonregular fastest-owner transitions. -/
theorem flood_turn_obligation_bound
    (K R E e X T : ℕ)
    (hK : K = R + E + 1)
    (hR : R ≤ e * (E + 1))
    (hET : E ≤ X + T) :
    K ≤ (e + 1) * (X + T + 1) := by
  have hKE : K ≤ (e + 1) * (E + 1) :=
    regularRun_exception_bound K R E e hK hR
  have hmono : (e + 1) * (E + 1) ≤ (e + 1) * (X + T + 1) := by
    exact Nat.mul_le_mul_left (e + 1) (by omega)
  exact hKE.trans hmono

/-- A skipped-fastest indicator and a selected raw-seam indicator can be
added pointwise.  On their overlap the two selected seam owners are low and
the full skipped fastest tooth supplies a third active owner. -/
theorem layeredMultiplicity
    (skipped seam multiplicity : ℕ)
    (hLayers : 1 + skipped + seam ≤ multiplicity) :
    skipped + seam ≤ multiplicity - 1 := by
  omega

/-- A normalized fastest tooth has dual mass at most `7c/(36h)`.  Hence
private mass `eta` carried by `K` such teeth gives the product form behind
`K >= ceil(36h eta/(7c))`. -/
theorem privateMass_occurrence_product
    (c h eta K : ℚ)
    (hh : 0 < h)
    (heta : eta ≤ K * (7 * c / (36 * h))) :
    36 * h * eta ≤ 7 * c * K := by
  have hscale : 0 ≤ 36 * h := by positivity
  calc
    36 * h * eta ≤ 36 * h * (K * (7 * c / (36 * h))) :=
      mul_le_mul_of_nonneg_left heta hscale
    _ = 7 * c * K := by
      field_simp [ne_of_gt hh]

/-- Slow-gap normalization `7c/6`, endpoint density floor `3/4`, and one
fastest-tooth length `1/(7h)` give exactly `c/(8h)`. -/
theorem skippedTooth_weight_coefficient (c h : ℚ) (hh : h ≠ 0) :
    (7 * c / 6) * (3 / 4) * (1 / (7 * h)) = c / (8 * h) := by
  field_simp [hh]
  ring

/-- One raw lcm seam has weighted mass at least `c/(16 lcm)`. -/
theorem rawSeam_weight_coefficient (c L : ℚ) (hL : L ≠ 0) :
    (7 * c / 6) * (3 / 4) * (1 / (14 * L)) = c / (16 * L) := by
  field_simp [hL]
  ring

/-- The same normalization turns a return seam
`Omega=(1/7)(S-6/h)` into weighted mass `(c/8)(S-6/h)`. -/
theorem turnSeam_weight_coefficient (c h S : ℚ) (hh : h ≠ 0) :
    (7 * c / 6) * (3 / 4) * ((1 / 7) * (S - 6 / h)) =
      (c / 8) * (S - 6 / h) := by
  field_simp [hh]
  ring

/-- General address return `R`: the packet bracket is
`S-(7R-1)/h`, seven times its exact raw seam mass. -/
theorem returnPacket_weight_coefficient (c h S R : ℚ) (hh : h ≠ 0) :
    (7 * c / 6) * (3 / 4) *
        ((1 / 7) * (S - (7 * R - 1) / h)) =
      (c / 8) * (S - (7 * R - 1) / h) := by
  field_simp [hh]
  ring

/-- The THM-1253 singleton-discrepancy conversion gives the unweighted
coefficient `7/6` for the same reciprocal tail expression. -/
theorem unweighted_tail_coefficient (S : ℚ) :
    (49 / 6) * ((1 / 7) * S) = (7 / 6) * S := by
  ring

/-- THM-1253's singleton-discrepancy conversion sends one lcm seam quantum
to the coefficient `7/12`. -/
theorem unweighted_lcmSeam_coefficient (L : ℚ) :
    (49 / 6) * (1 / (14 * L)) = (7 / 12) * (1 / L) := by
  ring

/-- If `h >= (7/2)d5`, any multi-low packet has reciprocal excess at least
`1/h`, since it contains at least two speeds no larger than `d5`. -/
theorem dominatedTurn_oneOverHigh
    (h d5 : ℚ)
    (hh : 0 < h)
    (hd : 0 < d5)
    (hdom : 7 * d5 ≤ 2 * h) :
    1 / h ≤ 2 / d5 - 6 / h := by
  have hfrac : (7 : ℚ) / h ≤ 2 / d5 :=
    (div_le_div_iff₀ hh hd).2 hdom
  rw [le_sub_iff_add_le]
  calc
    1 / h + 6 / h = 7 / h := by ring
    _ ≤ 2 / d5 := hfrac

/-- Abstract composition of the phase-free envelope, actual weighted excess,
and a pointwise layered tail invoice. -/
theorem weighted_invoice_consumer
    (F actual tail : ℚ)
    (hEnvelope : actual ≤ F)
    (hTail : tail ≤ actual) :
    tail ≤ F := by
  linarith

/-- The exact literal turn used by the referee has reciprocal excess
`41/780` and seam mass `41/5460`. -/
example : (1 / (5 : ℚ) + 1 / 12 - 6 / 26) = 41 / 780 := by norm_num
example : (1 / (7 : ℚ)) * (41 / 780) = 41 / 5460 := by norm_num

def rawProtrusion140 (d rho : ℚ) : ℚ :=
  1 / 2 + (7 / 6) * rho - d / 280

example : rawProtrusion140 254 (9 / 20) = 33 / 280 := by
  norm_num [rawProtrusion140]

example : rawProtrusion140 255 (1 / 8) = -89 / 336 := by
  norm_num [rawProtrusion140]

example : rawProtrusion140 256 (3 / 10) = -9 / 140 := by
  norm_num [rawProtrusion140]

example : rawProtrusion140 257 (11 / 40) = -163 / 1680 := by
  norm_num [rawProtrusion140]

example : rawProtrusion140 292 (2 / 5) = -8 / 105 := by
  norm_num [rawProtrusion140]

example : rawProtrusion140 1805 (3 / 8) = -617 / 112 := by
  norm_num [rawProtrusion140]

#print axioms regularRun_exception_bound
#print axioms flood_turn_obligation_bound
#print axioms layeredMultiplicity
#print axioms privateMass_occurrence_product
#print axioms skippedTooth_weight_coefficient
#print axioms rawSeam_weight_coefficient
#print axioms turnSeam_weight_coefficient
#print axioms returnPacket_weight_coefficient
#print axioms unweighted_tail_coefficient
#print axioms unweighted_lcmSeam_coefficient
#print axioms dominatedTurn_oneOverHigh
#print axioms weighted_invoice_consumer

end LRCFastestFloodTurnTailTax
