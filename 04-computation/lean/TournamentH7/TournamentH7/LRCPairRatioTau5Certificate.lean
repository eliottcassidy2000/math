import TournamentH7.LRCPairRatioFiniteCover

/-!
# The tau5 local-neighborhood certificate

The exact `tau5` primitive-ratio cover has 272 vertices.  Its quotient graph
has 76 isolated vertices and one 196-vertex component, 861 edges, 232
triangles, and no four-clique.  It is nevertheless 5-chromatic.  Thus a
coloring quotient cannot prove the required `K₄` bound.

Raw Zarankiewicz compression is also unsound here: the graph contains 3,440
four-cycles, and arbitrary vertex pairs have as many as thirteen common
neighbors.  Even the reciprocal two-circle cut (`ratio < 1` versus `ratio > 1`)
has 287 cross edges and 264 four-cycles.  The certificate below instead reduces
every quotient edge to the primitive
cross-product pair and checks that every local neighborhood is triangle-free.
The finite replay uses kernel reduction only; it does not use `native_decide`
or a native SAT oracle.

Tournament audit: vertices are primitive ratios and the binary relation is
membership of either directed quotient in the `tau5` cover.  Orienting all
pairs by the ambient increasing rational order gives the transitive tournament:
score multiset `0,1,...,271`, no directed cycle, 272 singleton SCCs, and one
Hamiltonian path.  That tournament forgets the quotient-edge predicate, so the
certificate retains the sparse graph instead.  It has degree
histogram
`0^76 2^22 3^2 4^42 6^20 8^28 10^30 12^20 14^10 16^4 18^8 22^4 26^2 30^2 36^2`.
This challenges the assumption that runners, residues, or fixed circle
sections form a faithful small quotient: modular residue pairs identify both
edges and nonedges, while a five-color quotient creates a `K₅` and destroys
the preserved predicate.  The reduced cross-product graph preserves exactly
the quotient edge implication needed by the LRC clique obstruction.
-/

namespace LonelyRunner
namespace LRCPairRatioTau5

open Finset SimpleGraph
open LRCB5ContinuumFloor LRCPairCovarianceKernel LRCWeightedRatioLayer
open LRCPairRatioQuotient LRCPairRatioFiniteCover

set_option maxRecDepth 100000

/-- Membership in a finite primitive-pair cover recovers the strict primitive
deficit inequality whenever the supplied fraction is already reduced. -/
theorem primitiveDeficit_above_of_mem_finiteAllowedRatios
    (threshold : ℚ) (cap first second : ℕ)
    (hsecond : 0 < second)
    (hcoprime : Nat.Coprime first second)
    (hmem : (first : ℚ) / second ∈ finiteAllowedRatios threshold cap) :
    threshold < primitiveDeficit first second := by
  rw [finiteAllowedRatios, Finset.mem_image] at hmem
  obtain ⟨pair, hpair, hratio⟩ := hmem
  rw [primitivePairCandidates, Finset.mem_biUnion] at hpair
  obtain ⟨candidateFirst, _hfirstRange, hpair⟩ := hpair
  rw [Finset.mem_image] at hpair
  obtain ⟨candidateSecond, hcandidate, hpairEq⟩ := hpair
  rw [Finset.mem_filter] at hcandidate
  rcases hcandidate with
    ⟨_hsecondRange, hcandidateFirst, hcandidateSecond,
      hcandidateCoprime, habove⟩
  subst pair
  have hratio' :
      (first : ℚ) / second = (candidateFirst : ℚ) / candidateSecond :=
    hratio.symm
  have hunique := Rat.div_int_inj
    (a := (first : ℤ)) (b := (second : ℤ))
    (c := (candidateFirst : ℤ)) (d := (candidateSecond : ℤ))
    (by exact_mod_cast hsecond) (by exact_mod_cast hcandidateSecond)
    (by simpa using hcoprime) (by simpa using hcandidateCoprime)
    (by simpa using hratio')
  have hfirstEq : first = candidateFirst := by exact_mod_cast hunique.1
  have hsecondEq : second = candidateSecond := by exact_mod_cast hunique.2
  simpa [hfirstEq, hsecondEq] using habove

set_option maxRecDepth 100000 in
def numeratorData : Array ℕ :=
  #[
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 2, 1, 1, 1, 1, 2, 1, 2, 1, 1, 1,
    1, 2, 3, 2, 1, 2, 1, 1, 1, 1, 3, 2,
    2, 2, 3, 3, 2, 3, 1, 3, 2, 1, 2, 1,
    1, 1, 3, 3, 1, 3, 2, 2, 2, 3, 3, 4,
    4, 3, 2, 3, 4, 1, 4, 3, 2, 3, 1, 2,
    1, 5, 5, 1, 4, 4, 1, 3, 3, 1, 5, 3,
    5, 5, 3, 2, 4, 4, 2, 5, 5, 5, 2, 5,
    3, 6, 3, 3, 4, 10, 11, 4, 3, 5, 5, 8,
    4, 5, 9, 10, 9, 6, 5, 11, 10, 11, 5, 11,
    17, 11, 12, 13, 15, 17, 16, 25, 17, 8, 18, 17,
    19, 9, 11, 17, 19, 19, 11, 9, 19, 12, 13, 8,
    11, 31, 31, 13, 10, 11, 23, 13, 22, 9, 23, 24,
    26, 11, 23, 25, 13, 22, 37, 38, 23, 39, 8, 25,
    26, 9, 37, 39, 10, 51, 52, 11, 23, 12, 37, 25,
    38, 51, 13, 53, 40, 27, 41, 65, 67, 52, 53, 37,
    39, 41, 65, 22, 67, 68, 23, 24, 25, 51, 26, 53,
    80, 27, 82, 55, 94, 95, 65, 67, 69, 109, 37, 38,
    39, 40, 81, 41, 83, 137, 95, 51, 52, 53, 54, 109,
    55, 123, 65, 66, 67, 68, 137, 69, 80, 81, 82, 83,
    94, 95, 96, 97, 109, 110, 111, 123, 124, 125, 137, 138,
    139, 152, 153, 166, 167, 181, 195, 209
  ]

def numerator (index : Fin 272) : ℕ := numeratorData[index.val]!

set_option maxRecDepth 100000 in
def denominatorData : Array ℕ :=
  #[
    209, 195, 181, 167, 166, 153, 152, 139, 138, 137, 125, 124,
    123, 111, 110, 109, 97, 96, 95, 94, 83, 82, 81, 80,
    69, 137, 68, 67, 66, 65, 123, 55, 109, 54, 53, 52,
    51, 95, 137, 83, 41, 81, 40, 39, 38, 37, 109, 69,
    67, 65, 95, 94, 55, 82, 27, 80, 53, 26, 51, 25,
    24, 23, 68, 67, 22, 65, 41, 39, 37, 53, 52, 67,
    65, 41, 27, 40, 53, 13, 51, 38, 25, 37, 12, 23,
    11, 52, 51, 10, 39, 37, 9, 26, 25, 8, 39, 23,
    38, 37, 22, 13, 25, 23, 11, 26, 24, 23, 9, 22,
    13, 23, 11, 10, 13, 31, 31, 11, 8, 13, 12, 19,
    9, 11, 19, 19, 17, 11, 9, 19, 17, 18, 8, 17,
    25, 16, 17, 15, 13, 12, 11, 17, 11, 5, 11, 10,
    11, 5, 6, 9, 10, 9, 5, 4, 8, 5, 5, 3,
    4, 11, 10, 4, 3, 3, 6, 3, 5, 2, 5, 5,
    5, 2, 4, 4, 2, 3, 5, 5, 3, 5, 1, 3,
    3, 1, 4, 4, 1, 5, 5, 1, 2, 1, 3, 2,
    3, 4, 1, 4, 3, 2, 3, 4, 4, 3, 3, 2,
    2, 2, 3, 1, 3, 3, 1, 1, 1, 2, 1, 2,
    3, 1, 3, 2, 3, 3, 2, 2, 2, 3, 1, 1,
    1, 1, 2, 1, 2, 3, 2, 1, 1, 1, 1, 2,
    1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1
  ]

def denominator (index : Fin 272) : ℕ := denominatorData[index.val]!

def ratio (index : Fin 272) : ℚ :=
  (numerator index : ℚ) / denominator index

def tau5Vertices : Finset ℚ := Finset.univ.image ratio

def crossNumerator (first second : Fin 272) : ℕ :=
  numerator first * denominator second

def crossDenominator (first second : Fin 272) : ℕ :=
  denominator first * numerator second

def crossGcd (first second : Fin 272) : ℕ :=
  Nat.gcd (crossNumerator first second) (crossDenominator first second)

def reducedCrossNumerator (first second : Fin 272) : ℕ :=
  crossNumerator first second / crossGcd first second

def reducedCrossDenominator (first second : Fin 272) : ℕ :=
  crossDenominator first second / crossGcd first second

set_option maxRecDepth 100000 in
theorem numerator_pos (index : Fin 272) : 0 < numerator index := by
  decide +kernel +revert

set_option maxRecDepth 100000 in
theorem denominator_pos (index : Fin 272) : 0 < denominator index := by
  decide +kernel +revert

theorem crossNumerator_pos (first second : Fin 272) :
    0 < crossNumerator first second := by
  exact Nat.mul_pos (numerator_pos first) (denominator_pos second)

theorem crossDenominator_pos (first second : Fin 272) :
    0 < crossDenominator first second := by
  exact Nat.mul_pos (denominator_pos first) (numerator_pos second)

theorem crossGcd_pos (first second : Fin 272) :
    0 < crossGcd first second := by
  exact Nat.gcd_pos_of_pos_left _ (crossNumerator_pos first second)

theorem reducedCrossNumerator_pos (first second : Fin 272) :
    0 < reducedCrossNumerator first second := by
  exact Nat.div_pos
    (Nat.le_of_dvd (crossNumerator_pos first second)
      (Nat.gcd_dvd_left _ _))
    (crossGcd_pos first second)

theorem reducedCrossDenominator_pos (first second : Fin 272) :
    0 < reducedCrossDenominator first second := by
  exact Nat.div_pos
    (Nat.le_of_dvd (crossDenominator_pos first second)
      (Nat.gcd_dvd_right _ _))
    (crossGcd_pos first second)

theorem reducedCross_coprime (first second : Fin 272) :
    Nat.Coprime (reducedCrossNumerator first second)
      (reducedCrossDenominator first second) := by
  exact Nat.coprime_div_gcd_div_gcd (crossGcd_pos first second)

theorem ratio_div_eq_reducedCross (first second : Fin 272) :
    ratio first / ratio second =
      (reducedCrossNumerator first second : ℚ) /
        reducedCrossDenominator first second := by
  have hnfirst : (numerator first : ℚ) ≠ 0 := by
    exact_mod_cast (numerator_pos first).ne'
  have hnsecond : (numerator second : ℚ) ≠ 0 := by
    exact_mod_cast (numerator_pos second).ne'
  have hdfirst : (denominator first : ℚ) ≠ 0 := by
    exact_mod_cast (denominator_pos first).ne'
  have hdsecond : (denominator second : ℚ) ≠ 0 := by
    exact_mod_cast (denominator_pos second).ne'
  have hcross : ratio first / ratio second =
      (crossNumerator first second : ℚ) /
        crossDenominator first second := by
    unfold ratio crossNumerator crossDenominator
    push_cast
    field_simp
  have hfirstFactor :
      reducedCrossNumerator first second * crossGcd first second =
        crossNumerator first second :=
    Nat.div_mul_cancel (Nat.gcd_dvd_left _ _)
  have hsecondFactor :
      reducedCrossDenominator first second * crossGcd first second =
        crossDenominator first second :=
    Nat.div_mul_cancel (Nat.gcd_dvd_right _ _)
  rw [hcross, ← hfirstFactor, ← hsecondFactor]
  push_cast
  exact mul_div_mul_right _ _ (by
    exact_mod_cast (crossGcd_pos first second).ne')

theorem reducedCross_above_of_tau5_mem (first second : Fin 272)
    (hmem : ratio first / ratio second ∈ tau5FiniteRatios) :
    tau5 < primitiveDeficit (reducedCrossNumerator first second)
      (reducedCrossDenominator first second) := by
  apply primitiveDeficit_above_of_mem_finiteAllowedRatios
    tau5 431 (reducedCrossNumerator first second)
      (reducedCrossDenominator first second)
    (reducedCrossDenominator_pos first second)
    (reducedCross_coprime first second)
  rw [← ratio_div_eq_reducedCross first second]
  simpa [tau5FiniteRatios] using hmem

/-- Arithmetic graph on the indexed cover.  Unlike a quotient-witness scan,
each adjacency test reduces only one primitive cross-product pair. -/
def indexGraph : SimpleGraph (Fin 272) where
  Adj first second := first ≠ second ∧
    (tau5 < primitiveDeficit (reducedCrossNumerator first second)
        (reducedCrossDenominator first second) ∨
      tau5 < primitiveDeficit (reducedCrossNumerator second first)
        (reducedCrossDenominator second first))
  symm := by
    intro first second hadj
    exact ⟨hadj.1.symm, hadj.2.symm⟩
  loopless := ⟨by
    intro vertex hadj
    exact hadj.1 rfl⟩

instance : DecidableRel indexGraph.Adj := by
  intro first second
  change Decidable (first ≠ second ∧ (_ ∨ _))
  infer_instance

set_option maxRecDepth 100000 in
def neighborData : Array (List ℕ) :=
  #[
    [14, 18, 37, 44, 52, 96, 102, 107, 121, 123],
    [59, 65, 67, 72, 82, 87, 90, 94],
    [],
    [],
    [],
    [22, 26, 78, 82, 120, 124],
    [18, 23, 50, 79, 82, 104, 123, 130],
    [],
    [60, 82, 83, 101],
    [],
    [80, 87],
    [42, 113],
    [53, 73, 90, 106],
    [81, 82, 89, 90],
    [0, 42, 59, 84, 87, 107, 115, 148],
    [],
    [],
    [28, 42, 82, 93, 118, 133],
    [0, 6, 31, 42, 80, 96, 119, 127, 141, 150],
    [],
    [],
    [],
    [5, 90, 147],
    [6, 31, 60, 87, 93, 111, 133, 148],
    [82, 90, 95, 101],
    [],
    [5, 42, 58, 116, 128, 151],
    [],
    [17, 52, 60, 82, 90, 98, 102, 115, 118, 138],
    [59, 67, 80, 85, 87, 92, 94, 99, 100, 103, 111, 117],
    [53, 73, 106, 120],
    [18, 23, 59, 80, 87, 102, 107, 121, 138, 144],
    [],
    [60, 82, 106, 120],
    [],
    [49, 65, 67, 82, 87, 93, 99, 103, 108, 116, 118, 130],
    [54, 62, 120, 124],
    [0, 52, 59, 96, 100, 123, 127, 164],
    [],
    [],
    [],
    [90, 106],
    [11, 14, 17, 18, 26, 52, 59, 82, 92, 104, 111, 130, 133, 143, 152, 153, 156, 158],
    [65, 70, 72, 82, 90, 91, 106, 108, 112, 118, 120, 126],
    [0, 64, 127, 169],
    [],
    [],
    [90, 95, 106, 109],
    [],
    [35, 80, 87, 88, 100, 103, 111, 112, 117, 141],
    [6, 75, 119, 167],
    [],
    [0, 28, 37, 42, 80, 87, 100, 107, 115, 121, 138, 144, 153, 175],
    [12, 30, 151, 165],
    [36, 82, 120, 147],
    [87, 93, 111, 116],
    [],
    [72, 87, 88, 93, 94, 111, 112, 116, 117, 130],
    [26, 74, 124, 155],
    [1, 14, 29, 31, 37, 42, 65, 135, 141, 148, 150, 154, 164, 177],
    [8, 23, 28, 33, 75, 87, 90, 116, 118, 126, 151, 156, 160, 170],
    [],
    [36, 75, 82, 124, 128, 151],
    [],
    [44, 82, 87, 93, 115, 121, 125, 144],
    [1, 35, 43, 59, 92, 99, 111, 117, 136, 145, 153, 181],
    [],
    [1, 29, 35, 91, 106, 108, 120, 155, 160, 184],
    [],
    [],
    [43, 99, 103, 111, 116, 151],
    [],
    [1, 43, 57, 100, 117, 141, 153, 189],
    [12, 30, 165, 181],
    [58, 147],
    [50, 60, 62, 92, 130, 143, 145, 152],
    [],
    [],
    [5, 90, 134, 189],
    [6, 98, 127, 189],
    [10, 18, 29, 31, 49, 52, 150, 154, 164, 168, 175, 184],
    [13, 181],
    [1, 5, 6, 8, 13, 17, 24, 28, 33, 35, 42, 43, 54, 62, 64, 85, 98, 106, 116, 129, 135, 146, 147, 151, 159, 160, 163, 165, 169, 170, 178, 182, 188, 192, 193, 199],
    [8, 189],
    [14, 184],
    [29, 82, 99, 108, 118, 130, 136, 171],
    [],
    [1, 10, 14, 23, 29, 31, 35, 49, 52, 55, 57, 60, 64, 91, 100, 130, 135, 150, 153, 154, 155, 156, 159, 168, 169, 172, 178, 187, 191, 204],
    [49, 57, 108, 120, 155, 160],
    [13, 189],
    [1, 12, 13, 22, 24, 28, 41, 43, 47, 60, 78, 94, 102, 129, 135, 137, 155, 162, 163, 165, 173, 176, 181, 190, 198, 206],
    [43, 67, 87, 111, 112, 116, 117, 136, 151, 165],
    [29, 42, 65, 75, 141, 154, 167, 177],
    [17, 23, 35, 55, 57, 64, 102, 111, 118, 133, 156, 159, 160, 172, 184, 189],
    [1, 29, 57, 90, 108, 112, 118, 126, 136, 160, 179, 212],
    [24, 47, 165, 181],
    [0, 18, 37, 107, 127, 171, 191, 219],
    [],
    [28, 79, 82, 111, 115, 116, 121, 142, 144, 181],
    [29, 35, 65, 70, 85, 141, 155, 160, 178, 184],
    [29, 37, 49, 52, 72, 87, 141, 154, 164, 168, 175, 186],
    [8, 24, 189, 211],
    [0, 28, 31, 90, 93, 138, 142, 184, 189, 227],
    [29, 35, 49, 70, 112, 130, 160, 171, 184, 191],
    [6, 42, 116, 126, 179, 221],
    [],
    [12, 30, 33, 41, 43, 47, 67, 82, 115, 129, 155, 163, 176, 180, 181, 189, 198, 218],
    [0, 14, 31, 52, 96, 115, 118, 125, 130, 144, 171, 191, 212, 234],
    [35, 43, 67, 85, 88, 94, 145, 151, 153, 165, 181, 189],
    [47, 181],
    [],
    [23, 29, 42, 49, 55, 57, 65, 70, 91, 93, 98, 150, 153, 154, 168, 172, 177, 178, 183, 189, 204, 211],
    [43, 49, 57, 91, 94, 103, 141, 153, 155, 178, 184, 189],
    [11, 229],
    [],
    [14, 28, 52, 64, 98, 106, 107, 141, 142, 155, 178, 184, 211, 229],
    [26, 35, 55, 57, 60, 70, 82, 91, 98, 104, 125, 133, 145, 156, 159, 165, 172, 181, 183, 184, 204, 213],
    [29, 49, 57, 65, 72, 91, 160, 171, 179, 184, 191, 212],
    [17, 28, 35, 43, 60, 85, 93, 94, 107, 146, 159, 160, 163, 184, 199, 206, 219, 229],
    [18, 50, 196, 229],
    [5, 30, 33, 36, 43, 54, 67, 88, 163, 180, 189, 201, 209, 211, 218, 245],
    [0, 31, 52, 64, 98, 160, 184, 191, 212, 253],
    [],
    [0, 6, 37, 212, 248, 257],
    [5, 36, 58, 62, 189, 197, 217, 249],
    [64, 107, 116, 138, 153, 189],
    [43, 60, 94, 104, 155, 163, 196, 206],
    [18, 37, 44, 79, 96, 164, 173, 207, 219, 240],
    [26, 62, 196, 229],
    [82, 90, 106, 156, 169, 173],
    [6, 35, 42, 57, 75, 85, 87, 103, 107, 156, 159, 171, 172, 179, 199, 212, 222, 253],
    [],
    [],
    [17, 23, 42, 93, 116, 146, 169, 219, 240, 243],
    [78, 181],
    [59, 82, 87, 90, 177, 180, 186, 206],
    [65, 85, 91, 94, 181, 184, 189, 212],
    [90, 193],
    [28, 31, 52, 102, 125, 155, 178, 229, 248, 254],
    [],
    [],
    [18, 49, 59, 72, 92, 99, 100, 112, 115, 164, 168, 184, 186, 196, 214, 229, 236, 265],
    [98, 102, 115, 165, 181, 189],
    [42, 75, 209, 245],
    [31, 52, 64, 98, 107, 175, 192, 227, 234, 253],
    [65, 75, 108, 116, 167, 177, 211, 228],
    [82, 118, 133, 155, 164, 207],
    [22, 54, 74, 82, 209, 213, 235, 266],
    [14, 23, 59, 234, 265, 271],
    [],
    [18, 59, 80, 87, 111, 173, 207, 219, 240, 271],
    [26, 53, 60, 62, 70, 82, 91, 108, 183, 204, 217, 228, 235, 238, 241, 266],
    [42, 75, 221, 253],
    [42, 52, 65, 72, 87, 108, 111, 112, 125, 164, 177, 178, 186, 211, 228, 236, 243, 254],
    [59, 80, 87, 92, 100, 111, 180, 199, 206, 214, 222, 242],
    [58, 67, 87, 88, 90, 99, 106, 112, 115, 126, 138, 146, 167, 173, 180, 189, 201, 211, 214, 216, 236, 245],
    [42, 60, 87, 93, 116, 129, 130, 164, 165, 173, 207, 219, 243, 257],
    [],
    [42, 260],
    [82, 87, 93, 116, 118, 130, 168, 177, 180, 214, 222, 228],
    [60, 67, 82, 88, 93, 94, 99, 103, 117, 118, 121, 173, 178, 180, 201, 206, 214, 216, 222, 229, 242, 248],
    [],
    [90, 224],
    [82, 90, 106, 118, 120, 126, 177, 183, 186, 204, 228, 236],
    [37, 59, 80, 100, 127, 141, 146, 153, 156, 175, 219, 240, 257, 271],
    [53, 73, 82, 90, 91, 95, 108, 116, 142, 156, 189, 204, 224, 228, 230, 238, 241, 259],
    [],
    [50, 92, 145, 155, 229, 265],
    [80, 87, 100, 111, 141, 159, 201, 222, 236, 242],
    [44, 82, 87, 129, 133, 178, 181, 240, 243, 271],
    [60, 82, 247, 263],
    [85, 96, 103, 107, 117, 130, 184, 199, 219, 222, 234, 242],
    [87, 93, 111, 116, 130, 186, 201, 206, 236, 242],
    [90, 127, 129, 150, 155, 156, 160, 189, 192, 243],
    [],
    [52, 80, 100, 144, 164, 234, 253, 271],
    [90, 106, 224, 247],
    [59, 92, 111, 135, 145, 153, 159, 163, 181, 214, 242, 270],
    [82, 87, 99, 111, 112, 115, 138, 153, 160, 169, 207, 214, 216, 236, 248, 254],
    [94, 104, 117, 130, 196, 206, 229, 242],
    [106, 120, 135, 154, 155, 159, 160, 184, 204, 228],
    [65, 73, 81, 90, 95, 98, 106, 108, 109, 116, 134, 136, 142, 169, 177, 193, 211, 224, 228, 230, 243, 247, 249, 258, 259, 270],
    [82, 258],
    [111, 116, 151, 163, 214, 222],
    [67, 80, 84, 93, 99, 102, 103, 112, 115, 116, 117, 118, 121, 136, 141, 171, 180, 207, 211, 214, 216, 219, 222, 236, 240, 242, 248, 257, 261, 270],
    [],
    [100, 135, 141, 153, 163, 172, 189, 242],
    [87, 257],
    [82, 263],
    [72, 78, 79, 83, 89, 93, 101, 102, 106, 108, 111, 112, 120, 124, 125, 136, 142, 155, 165, 173, 186, 207, 209, 217, 228, 229, 236, 238, 243, 247, 254, 258, 263, 265, 266, 270],
    [90, 258],
    [87, 96, 103, 107, 117, 121, 219, 222, 240, 242, 253, 261],
    [82, 144, 173, 265],
    [82, 137, 181, 266],
    [],
    [],
    [119, 126, 128, 141, 179, 209, 211, 221],
    [124, 213],
    [90, 106, 241, 259],
    [82, 118, 130, 154, 171, 214, 228, 270],
    [],
    [120, 155, 160, 168, 172, 228],
    [],
    [],
    [87, 111, 116, 151, 163, 165, 180, 236, 242, 270],
    [],
    [90, 118, 126, 135, 154, 160, 172, 179, 212, 228, 236, 270],
    [127, 146, 150, 156, 178, 184, 189, 227],
    [],
    [120, 143, 147, 189, 196, 235],
    [],
    [101, 111, 115, 120, 145, 153, 155, 181, 184, 196, 238, 243, 248, 263],
    [94, 107, 117, 121, 123, 130, 136, 206, 229, 234, 240, 242, 257, 270],
    [116, 147, 197, 245],
    [141, 154, 155, 159, 160, 177, 178, 183, 184, 199],
    [],
    [155, 160, 178, 184],
    [124, 151, 189, 235],
    [106, 120, 241, 259],
    [96, 118, 127, 133, 150, 156, 164, 171, 184, 191, 229, 234, 243, 271],
    [],
    [104, 152, 196, 265],
    [130, 154, 159, 160, 168, 171, 183, 184, 191, 236],
    [],
    [162, 165, 176, 181],
    [],
    [],
    [102, 144, 207, 271],
    [145, 151, 153, 159, 163, 165, 180, 181, 189, 199, 201, 206],
    [113, 115, 118, 119, 128, 138, 141, 160, 167, 179, 189, 212, 219, 245, 253, 254, 257, 260],
    [165, 181],
    [],
    [],
    [],
    [107, 144, 148, 171, 175, 212, 219, 271],
    [147, 151, 209, 217],
    [141, 153, 155, 163, 168, 172, 178, 184, 189, 204, 206, 222],
    [],
    [151, 165, 189, 211],
    [],
    [127, 133, 150, 164, 169, 184, 191, 212, 248, 253],
    [151, 165, 198, 218],
    [154, 160, 168, 171, 172, 177, 179, 184, 186, 191, 204, 212],
    [133, 153, 156, 169, 173, 181, 189, 211, 219, 254],
    [],
    [120, 143, 155, 213, 229, 266],
    [],
    [170, 176, 181, 189],
    [123, 138, 160, 178, 184, 211, 240, 265],
    [124, 181, 266],
    [],
    [],
    [],
    [121, 130, 144, 152, 175, 191, 229, 240, 265, 271],
    [138, 153, 178, 189, 229, 243],
    [],
    [],
    [123, 156, 164, 184, 187, 212, 229, 271],
    [181, 182, 189, 190],
    [165, 181, 198, 218],
    [158, 229],
    [184, 191],
    [],
    [170, 188, 189, 211],
    [],
    [141, 148, 167, 189, 192, 221, 248, 253],
    [147, 151, 189, 193, 245, 249],
    [],
    [],
    [],
    [177, 181, 184, 189, 199, 204, 206, 212],
    [148, 150, 164, 169, 175, 219, 227, 234, 253, 257]
  ]

def neighbor (index : Fin 272) : List ℕ :=
  neighborData[index.val]!

/-- Explicit sparse supergraph used by the local-neighborhood replay. -/
def tableGraph : SimpleGraph (Fin 272) where
  Adj first second := first ≠ second ∧
    (second.val ∈ neighbor first ∨ first.val ∈ neighbor second)
  symm := by
    intro first second hadj
    exact ⟨hadj.1.symm, hadj.2.symm⟩
  loopless := ⟨by
    intro vertex hadj
    exact hadj.1 rfl⟩

instance : DecidableRel tableGraph.Adj := by
  intro first second
  change Decidable (first ≠ second ∧ (_ ∨ _))
  infer_instance

def indexBlock : Fin 17 → List (Fin 272) :=
  ![[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
    [16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31],
    [32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47],
    [48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63],
    [64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79],
    [80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95],
    [96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111],
    [112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123,
      124, 125, 126, 127],
    [128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139,
      140, 141, 142, 143],
    [144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155,
      156, 157, 158, 159],
    [160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171,
      172, 173, 174, 175],
    [176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187,
      188, 189, 190, 191],
    [192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203,
      204, 205, 206, 207],
    [208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219,
      220, 221, 222, 223],
    [224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235,
      236, 237, 238, 239],
    [240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251,
      252, 253, 254, 255],
    [256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267,
      268, 269, 270, 271]]

def indexBlockIndex (index : Fin 272) : Fin 17 :=
  ⟨index.val / 16, by omega⟩

/-- Pairwise soundness check from reduced cross-products to the explicit
neighbor table. -/
def edgeTableCheckBlock (block : Fin 17) : Bool :=
  (indexBlock block).all fun first =>
    (List.finRange 272).all fun second =>
      if indexGraph.Adj first second then
        decide (tableGraph.Adj first second)
      else true

/-- Block-local `K₄` check with implications nested before later quantifiers,
so only actual table edges and common neighbors are extended. -/
def localTableCheckBlock (block : Fin 17) : Bool :=
  (indexBlock block).all fun first =>
    (List.finRange 272).all fun second =>
      if tableGraph.Adj first second then
        (List.finRange 272).all fun third =>
          if tableGraph.Adj first third then
            if tableGraph.Adj second third then
              (List.finRange 272).all fun fourth =>
                if tableGraph.Adj first fourth then
                  if tableGraph.Adj second fourth then
                    decide (¬tableGraph.Adj third fourth)
                  else true
                else true
            else true
          else true
      else true

def isolatedIndices : Finset (Fin 272) :=
  Finset.univ.filter fun index => (neighbor index).isEmpty

def degreeIndices (degree : ℕ) : Finset (Fin 272) :=
  Finset.univ.filter fun index => (neighbor index).length == degree

/-- Kernel-checked score histogram and edge-handshake count. -/
theorem quotient_graph_fingerprint :
    isolatedIndices.card = 76 ∧
      (degreeIndices 0).card = 76 ∧
      (degreeIndices 2).card = 22 ∧
      (degreeIndices 3).card = 2 ∧
      (degreeIndices 4).card = 42 ∧
      (degreeIndices 6).card = 20 ∧
      (degreeIndices 8).card = 28 ∧
      (degreeIndices 10).card = 30 ∧
      (degreeIndices 12).card = 20 ∧
      (degreeIndices 14).card = 10 ∧
      (degreeIndices 16).card = 4 ∧
      (degreeIndices 18).card = 8 ∧
      (degreeIndices 22).card = 4 ∧
      (degreeIndices 26).card = 2 ∧
      (degreeIndices 30).card = 2 ∧
      (degreeIndices 36).card = 2 ∧
      (Finset.univ.sum fun index : Fin 272 => (neighbor index).length) =
        1722 := by
  decide +kernel

def candidateListed (first second : ℕ) : Bool :=
  (List.finRange 272).any fun index =>
    numerator index == first && denominator index == second

/-- Hyperbola-shaped completeness check for the exact primitive-pair list. -/
def candidateCheck : Bool :=
  (List.range 432).all fun first =>
    (List.range (431 / first + 1)).all fun second =>
      if 0 < first ∧ 0 < second ∧ Nat.Coprime first second ∧
          tau5 < primitiveDeficit first second then
        candidateListed first second
      else true

set_option maxHeartbeats 100000000 in
theorem candidateCheck_true : candidateCheck = true := by
  decide +kernel

theorem candidateListed_of_primitive_candidate
    (first second : ℕ)
    (hfirstRange : first < 432)
    (hsecondRange : second < 431 / first + 1)
    (hfirst : 0 < first) (hsecond : 0 < second)
    (hcoprime : Nat.Coprime first second)
    (habove : tau5 < primitiveDeficit first second) :
    candidateListed first second = true := by
  have hfirstCheck := (List.all_eq_true.mp candidateCheck_true)
    first (List.mem_range.mpr hfirstRange)
  have hsecondCheck := (List.all_eq_true.mp hfirstCheck)
    second (List.mem_range.mpr hsecondRange)
  rw [if_pos ⟨hfirst, hsecond, hcoprime, habove⟩] at hsecondCheck
  exact hsecondCheck

theorem tau5FiniteRatios_subset_tau5Vertices :
    tau5FiniteRatios ⊆ tau5Vertices := by
  intro candidate hcandidate
  rw [tau5FiniteRatios, finiteAllowedRatios, Finset.mem_image] at hcandidate
  obtain ⟨pair, hpair, hratio⟩ := hcandidate
  rw [primitivePairCandidates, Finset.mem_biUnion] at hpair
  obtain ⟨candidateFirst, hfirstRange, hpair⟩ := hpair
  rw [Finset.mem_image] at hpair
  obtain ⟨candidateSecond, hcandidateData, hpairEq⟩ := hpair
  rw [Finset.mem_filter] at hcandidateData
  rcases hcandidateData with
    ⟨hsecondRange, hfirst, hsecond, hcoprime, habove⟩
  subst pair
  have hlisted := candidateListed_of_primitive_candidate
    candidateFirst candidateSecond
    (Finset.mem_range.mp hfirstRange) (Finset.mem_range.mp hsecondRange)
    hfirst hsecond hcoprime habove
  obtain ⟨index, _hindex, hmatch⟩ := List.any_eq_true.mp hlisted
  have hmatch' : numerator index = candidateFirst ∧
      denominator index = candidateSecond := by
    simpa [candidateListed] using hmatch
  rw [tau5Vertices, Finset.mem_image]
  refine ⟨index, Finset.mem_univ index, ?_⟩
  simpa [ratio, hmatch'.1, hmatch'.2] using hratio

set_option maxHeartbeats 100000000 in
theorem vertex_candidate_data (index : Fin 272) :
    numerator index < 432 ∧
      denominator index < 431 / numerator index + 1 ∧
      Nat.Coprime (numerator index) (denominator index) ∧
      tau5 < primitiveDeficit (numerator index) (denominator index) := by
  decide +kernel +revert

theorem tau5Vertices_subset_tau5FiniteRatios :
    tau5Vertices ⊆ tau5FiniteRatios := by
  intro candidate hcandidate
  rw [tau5Vertices, Finset.mem_image] at hcandidate
  obtain ⟨index, _hindex, rfl⟩ := hcandidate
  have hdata := vertex_candidate_data index
  rw [tau5FiniteRatios, finiteAllowedRatios, Finset.mem_image]
  refine ⟨(numerator index, denominator index), ?_, rfl⟩
  rw [primitivePairCandidates, Finset.mem_biUnion]
  refine ⟨numerator index, Finset.mem_range.mpr hdata.1, ?_⟩
  rw [Finset.mem_image]
  refine ⟨denominator index, ?_, rfl⟩
  rw [Finset.mem_filter]
  exact ⟨Finset.mem_range.mpr hdata.2.1,
    numerator_pos index, denominator_pos index,
    hdata.2.2.1, hdata.2.2.2⟩

theorem tau5FiniteRatios_eq_tau5Vertices :
    tau5FiniteRatios = tau5Vertices := by
  apply Finset.Subset.antisymm
  · exact tau5FiniteRatios_subset_tau5Vertices
  · exact tau5Vertices_subset_tau5FiniteRatios


#print axioms primitiveDeficit_above_of_mem_finiteAllowedRatios
#print axioms ratio_div_eq_reducedCross
#print axioms reducedCross_above_of_tau5_mem
#print axioms quotient_graph_fingerprint
#print axioms candidateCheck_true
#print axioms vertex_candidate_data
#print axioms tau5FiniteRatios_eq_tau5Vertices

end LRCPairRatioTau5
end LonelyRunner
