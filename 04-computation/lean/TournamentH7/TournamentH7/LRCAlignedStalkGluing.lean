/-
  TournamentH7.LRCAlignedStalkGluing

  The algebraic gluing layer for zero-color overlap stalks.  Zero determinant
  is transitive across a nonzero intermediate speed.  This file formalizes
  the diameter-two certificate needed by the current census carrier: a single
  zero-color anchor star forces the stalk to be complete, and the actual
  `Finset.gcd` then supplies the Bezout ledger:

      failWitness_i = (v_i / h) * r.

  The top-window estimate then forces `h*p = r*q`, and one fixed anchored
  stalk has at most `gcd(h,q)-1` active multipliers.  Selecting anchor stars
  from arbitrary connected carriers and aggregating stalks without reuse
  remain separate.

  Tournament-analysis audit: the carrier is a zero-color relation graph, not
  an orientation tournament.  The switch is `overlapDet = 0`; path closure
  preserves proportional witness data, while quotienting to a static sign
  tournament erases the zero edges that perform the gluing.  The natural tie
  path is therefore a spanning path in the zero-color graph.

  No `sorry`; no `native_decide`.
-/

import TournamentH7.LRCAlignedResonance
import TournamentH7.LRCSevenOverlapRelations
import Mathlib.RingTheory.PrincipalIdealDomain

namespace LonelyRunner
namespace LRC14Concrete

open scoped BigOperators Classical

/-- Zero overlap color is transitive across a runner of nonzero speed.  This
turns connected zero-color stalks into complete zero-color stalks by path
induction. -/
theorem overlapDet_eq_zero_trans
    (v : Fin 13 → ℤ) (q p : ℕ) (i j k : Fin 13)
    (hvj : v j ≠ 0)
    (hij : overlapDet v q p i j = 0)
    (hjk : overlapDet v q p j k = 0) :
    overlapDet v q p i k = 0 := by
  have hplucker := overlapDet_plucker v q p i j k
  rw [hjk, hij] at hplucker
  have hproduct : overlapDet v q p i k * v j = 0 := by
    nlinarith
  rcases mul_eq_zero.mp hproduct with hzero | hvzero
  · exact hzero
  · exact False.elim (hvj hvzero)

/-- A zero-color star is already a zero-color clique.  This is the finite
diameter-two form of path closure and reduces a stalk certificate from all
pair colors to one anchor spoke per vertex. -/
theorem zeroColorStalk_complete_of_anchor
    (v : Fin 13 → ℤ) (q p : ℕ) (stalk : Finset (Fin 13))
    (anchor : Fin 13) (hvanchor : v anchor ≠ 0)
    (hanchor : ∀ i ∈ stalk, overlapDet v q p i anchor = 0) :
    ∀ i ∈ stalk, ∀ j ∈ stalk, overlapDet v q p i j = 0 := by
  intro i hi j hj
  apply overlapDet_eq_zero_trans v q p i anchor j hvanchor (hanchor i hi)
  calc
    overlapDet v q p anchor j = -overlapDet v q p j anchor :=
      overlapDet_skew v q p j anchor
    _ = 0 := by rw [hanchor j hj]; simp

/-- **Finite zero-stalk factorization from a Bezout gcd ledger.**

If `h` divides every stalk speed and is expressed as an integer linear
combination of those speeds, then complete zero color forces every witness to
be `(v_i / h)` times one shared integer `r`.  In the intended application,
`h` is the positive gcd of the stalk speeds and `coeff` is a Bezout vector. -/
theorem zeroColorStalk_factorization_of_bezout
    (v : Fin 13 → ℤ) (q p : ℕ)
    (stalk : Finset (Fin 13)) (coeff : Fin 13 → ℤ) (h : ℤ)
    (hh : h ≠ 0)
    (hbezout : (∑ i ∈ stalk, coeff i * v i) = h)
    (hdiv : ∀ i ∈ stalk, h ∣ v i)
    (hzero : ∀ i ∈ stalk, ∀ j ∈ stalk,
      overlapDet v q p i j = 0) :
    ∃ r : ℤ, ∀ i ∈ stalk,
      failWitness v q p i = (v i / h) * r := by
  let r : ℤ := ∑ i ∈ stalk, coeff i * failWitness v q p i
  refine ⟨r, ?_⟩
  intro j hj
  have hcross : ∀ i ∈ stalk,
      v i * failWitness v q p j = v j * failWitness v q p i := by
    intro i hi
    have hz := hzero i hi j hj
    unfold overlapDet at hz
    omega
  have hscaled :
      h * failWitness v q p j = v j * r := by
    calc
      h * failWitness v q p j =
          (∑ i ∈ stalk, coeff i * v i) * failWitness v q p j := by
            rw [hbezout]
      _ = ∑ i ∈ stalk,
          (coeff i * v i) * failWitness v q p j := by
            rw [Finset.sum_mul]
      _ = ∑ i ∈ stalk,
          v j * (coeff i * failWitness v q p i) := by
            apply Finset.sum_congr rfl
            intro i hi
            rw [mul_assoc, hcross i hi]
            ring
      _ = v j * (∑ i ∈ stalk,
          coeff i * failWitness v q p i) := by
            rw [Finset.mul_sum]
      _ = v j * r := by rfl
  have hfactor : h * (v j / h) = v j :=
    Int.mul_ediv_cancel' (hdiv j hj)
  have hcancel :
      h * failWitness v q p j = h * ((v j / h) * r) := by
    calc
      h * failWitness v q p j = v j * r := hscaled
      _ = (h * (v j / h)) * r := by rw [hfactor]
      _ = h * ((v j / h) * r) := by ring
  exact mul_left_cancel₀ hh hcancel

/-- **Local aligned-stalk resonance law.**  The Bezout factorization, one bad
top runner, and the natural top-window inequality compose to the exact
Diophantine equation `h*p = r*q`.  This is the proof-facing bridge from a
zero-color stalk to `alignedResonantMultipliers`. -/
theorem zeroColorStalk_resonance_of_bezout
    (v : Fin 13 → ℤ) (q p : ℕ)
    (stalk : Finset (Fin 13)) (coeff : Fin 13 → ℤ) (h : ℤ)
    (hq : 0 < q) (hh : 0 < h)
    (hbezout : (∑ i ∈ stalk, coeff i * v i) = h)
    (hdiv : ∀ i ∈ stalk, h ∣ v i)
    (hzero : ∀ i ∈ stalk, ∀ j ∈ stalk,
      overlapDet v q p i j = 0)
    (top : Fin 13) (htop : top ∈ stalk) (hvtop : 0 < v top)
    (hbadtop : ¬ inBand v q p top)
    (hwindow : h * (q : ℤ) ≤ 14 * v top) :
    ∃ r : ℤ,
      (∀ i ∈ stalk, failWitness v q p i = (v i / h) * r) ∧
      h * (p : ℤ) = r * (q : ℤ) := by
  obtain ⟨r, hfactorization⟩ :=
    zeroColorStalk_factorization_of_bezout
      v q p stalk coeff h (ne_of_gt hh) hbezout hdiv hzero
  refine ⟨r, hfactorization, ?_⟩
  have htopWitness := hfactorization top htop
  have hbadBound := bad_at_witness v q p top hq hbadtop
  have htopFactor : h * (v top / h) = v top :=
    Int.mul_ediv_cancel' (hdiv top htop)
  let error : ℤ :=
    v top * (p : ℤ) - (v top / h * r) * (q : ℤ)
  have herrorBound : 14 * |error| < (q : ℤ) := by
    dsimp only [error]
    rw [← htopWitness]
    exact hbadBound
  have herrorIdentity :
      h * error = v top * (h * (p : ℤ) - r * (q : ℤ)) := by
    dsimp only [error]
    calc
      h * (v top * (p : ℤ) - (v top / h * r) * (q : ℤ)) =
          h * v top * (p : ℤ) -
            (h * (v top / h)) * r * (q : ℤ) := by ring
      _ = v top * (h * (p : ℤ) - r * (q : ℤ)) := by
        rw [htopFactor]
        ring
  have habsIdentity :
      v top * |h * (p : ℤ) - r * (q : ℤ)| = h * |error| := by
    have habs := congrArg abs herrorIdentity
    simpa [abs_mul, abs_of_pos hvtop, abs_of_pos hh] using habs.symm
  have hscaledError : h * (14 * |error|) < h * (q : ℤ) :=
    mul_lt_mul_of_pos_left herrorBound hh
  have hclose :
      14 * v top * |h * (p : ℤ) - r * (q : ℤ)| < h * (q : ℤ) := by
    calc
      14 * v top * |h * (p : ℤ) - r * (q : ℤ)| =
          14 * (v top * |h * (p : ℤ) - r * (q : ℤ)|) := by ring
      _ = 14 * (h * |error|) := by rw [habsIdentity]
      _ = h * (14 * |error|) := by ring
      _ < h * (q : ℤ) := hscaledError
  exact LRC14Grand.eq_of_aligned_resonance_closeness
    (v top) h (q : ℤ) (p : ℤ) r hvtop hh.le hwindow hclose

/-- The concrete `Finset.gcd` automatically supplies both halves of the
Bezout ledger: divisibility of every stalk speed and an integer linear
combination equal to the gcd. -/
theorem zeroColorStalk_factorization
    (v : Fin 13 → ℤ) (q p : ℕ) (stalk : Finset (Fin 13))
    (hgcd : stalk.gcd v ≠ 0)
    (hzero : ∀ i ∈ stalk, ∀ j ∈ stalk,
      overlapDet v q p i j = 0) :
    ∃ r : ℤ, ∀ i ∈ stalk,
      failWitness v q p i = (v i / stalk.gcd v) * r := by
  obtain ⟨coeff, hcoeff⟩ := Finset.gcd_eq_sum_mul stalk v
  have hbezout : (∑ i ∈ stalk, coeff i * v i) = stalk.gcd v := by
    calc
      (∑ i ∈ stalk, coeff i * v i) =
          ∑ i ∈ stalk, v i * coeff i := by
            apply Finset.sum_congr rfl
            intro i _hi
            ring
      _ = stalk.gcd v := hcoeff.symm
  have hdiv : ∀ i ∈ stalk, stalk.gcd v ∣ v i := by
    intro i hi
    exact Finset.gcd_dvd hi
  exact zeroColorStalk_factorization_of_bezout
    v q p stalk coeff (stalk.gcd v) hgcd hbezout hdiv hzero

/-- **Concrete local aligned-stalk theorem.**  A complete zero-color stalk with
positive concrete speed gcd, one bad positive top runner, and the top-window
inequality has a shared witness parameter satisfying the exact resonance law.
-/
theorem zeroColorStalk_resonance
    (v : Fin 13 → ℤ) (q p : ℕ) (stalk : Finset (Fin 13))
    (hq : 0 < q) (hgcd : 0 < stalk.gcd v)
    (hzero : ∀ i ∈ stalk, ∀ j ∈ stalk,
      overlapDet v q p i j = 0)
    (top : Fin 13) (htop : top ∈ stalk) (hvtop : 0 < v top)
    (hbadtop : ¬ inBand v q p top)
    (hwindow : stalk.gcd v * (q : ℤ) ≤ 14 * v top) :
    ∃ r : ℤ,
      (∀ i ∈ stalk,
        failWitness v q p i = (v i / stalk.gcd v) * r) ∧
      stalk.gcd v * (p : ℤ) = r * (q : ℤ) := by
  obtain ⟨coeff, hcoeff⟩ := Finset.gcd_eq_sum_mul stalk v
  have hbezout : (∑ i ∈ stalk, coeff i * v i) = stalk.gcd v := by
    calc
      (∑ i ∈ stalk, coeff i * v i) =
          ∑ i ∈ stalk, v i * coeff i := by
            apply Finset.sum_congr rfl
            intro i _hi
            ring
      _ = stalk.gcd v := hcoeff.symm
  have hdiv : ∀ i ∈ stalk, stalk.gcd v ∣ v i := by
    intro i hi
    exact Finset.gcd_dvd hi
  exact zeroColorStalk_resonance_of_bezout
    v q p stalk coeff (stalk.gcd v) hq hgcd hbezout hdiv hzero
      top htop hvtop hbadtop hwindow

/-- Multipliers supporting one fixed complete zero-color stalk.  The runner
set is fixed while multiplier activity is retained. -/
def zeroColorStalkMultipliers
    (v : Fin 13 → ℤ) (q : ℕ) (stalk : Finset (Fin 13)) : Finset ℕ :=
  (Finset.Ioo 0 q).filter fun p =>
    (∀ i ∈ stalk, ¬ inBand v q p i) ∧
    (∀ i ∈ stalk, ∀ j ∈ stalk, overlapDet v q p i j = 0)

/-- The cheaper anchor-star certificate for a fixed zero-color stalk. -/
def anchoredZeroColorStalkMultipliers
    (v : Fin 13 → ℤ) (q : ℕ) (stalk : Finset (Fin 13))
    (anchor : Fin 13) : Finset ℕ :=
  (Finset.Ioo 0 q).filter fun p =>
    (∀ i ∈ stalk, ¬ inBand v q p i) ∧
    (∀ i ∈ stalk, overlapDet v q p i anchor = 0)

/-- With a nonzero anchor speed, an anchor-star event is a complete
zero-color stalk event. -/
theorem anchoredZeroColorStalkMultipliers_subset
    (v : Fin 13 → ℤ) (q : ℕ) (stalk : Finset (Fin 13))
    (anchor : Fin 13) (hvanchor : v anchor ≠ 0) :
    anchoredZeroColorStalkMultipliers v q stalk anchor ⊆
      zeroColorStalkMultipliers v q stalk := by
  intro p hp
  unfold anchoredZeroColorStalkMultipliers at hp
  unfold zeroColorStalkMultipliers
  obtain ⟨hpWindow, hbad, hanchor⟩ := Finset.mem_filter.mp hp
  exact Finset.mem_filter.mpr
    ⟨hpWindow, hbad,
      zeroColorStalk_complete_of_anchor
        v q p stalk anchor hvanchor hanchor⟩

/-- Every event in a fixed complete zero-color stalk is an arithmetic
resonance for the concrete stalk gcd. -/
theorem zeroColorStalkMultipliers_subset_alignedResonant
    (v : Fin 13 → ℤ) (q : ℕ) (stalk : Finset (Fin 13))
    (hq : 0 < q) (hgcd : 0 < stalk.gcd v)
    (top : Fin 13) (htop : top ∈ stalk) (hvtop : 0 < v top)
    (hwindow : stalk.gcd v * (q : ℤ) ≤ 14 * v top) :
    zeroColorStalkMultipliers v q stalk ⊆
      LRC14Grand.alignedResonantMultipliers (stalk.gcd v).natAbs q := by
  intro p hp
  unfold zeroColorStalkMultipliers at hp
  obtain ⟨hpWindow, hbad, hzero⟩ := Finset.mem_filter.mp hp
  obtain ⟨r, _hfactorization, hresonance⟩ :=
    zeroColorStalk_resonance v q p stalk hq hgcd hzero
      top htop hvtop (hbad top htop) hwindow
  unfold LRC14Grand.alignedResonantMultipliers
  rw [Finset.mem_filter]
  refine ⟨hpWindow, ?_⟩
  have hdvdZ : (q : ℤ) ∣ stalk.gcd v * (p : ℤ) := by
    refine ⟨r, ?_⟩
    calc
      stalk.gcd v * (p : ℤ) = r * (q : ℤ) := hresonance
      _ = (q : ℤ) * r := by ring
  have hgcast : (((stalk.gcd v).natAbs : ℕ) : ℤ) = stalk.gcd v := by
    rw [Int.natCast_natAbs, abs_of_pos hgcd]
  have hdvdZ' : (q : ℤ) ∣
      (((stalk.gcd v).natAbs * p : ℕ) : ℤ) := by
    simpa [Int.natCast_mul, hgcast] using hdvdZ
  exact Int.natCast_dvd_natCast.mp hdvdZ'

/-- **Exact local aligned-activity bound.**  A fixed complete zero-color stalk
has at most `gcd(stalkGcd,q)-1` active multipliers in the top window. -/
theorem zeroColorStalkMultipliers_card_le
    (v : Fin 13 → ℤ) (q : ℕ) (stalk : Finset (Fin 13))
    (hq : 0 < q) (hgcd : 0 < stalk.gcd v)
    (top : Fin 13) (htop : top ∈ stalk) (hvtop : 0 < v top)
    (hwindow : stalk.gcd v * (q : ℤ) ≤ 14 * v top) :
    (zeroColorStalkMultipliers v q stalk).card ≤
      Nat.gcd (stalk.gcd v).natAbs q - 1 := by
  calc
    (zeroColorStalkMultipliers v q stalk).card ≤
        (LRC14Grand.alignedResonantMultipliers
          (stalk.gcd v).natAbs q).card :=
      Finset.card_le_card
        (zeroColorStalkMultipliers_subset_alignedResonant
          v q stalk hq hgcd top htop hvtop hwindow)
    _ = Nat.gcd (stalk.gcd v).natAbs q - 1 :=
      LRC14Grand.alignedResonantMultipliers_card
        (stalk.gcd v).natAbs q hq

/-- **Anchor-star aligned-activity bound.**  Only the colors incident to one
nonzero anchor have to vanish; transitivity supplies every omitted pair. -/
theorem anchoredZeroColorStalkMultipliers_card_le
    (v : Fin 13 → ℤ) (q : ℕ) (stalk : Finset (Fin 13))
    (hq : 0 < q) (hgcd : 0 < stalk.gcd v)
    (anchor top : Fin 13) (hvanchor : v anchor ≠ 0)
    (htop : top ∈ stalk) (hvtop : 0 < v top)
    (hwindow : stalk.gcd v * (q : ℤ) ≤ 14 * v top) :
    (anchoredZeroColorStalkMultipliers v q stalk anchor).card ≤
      Nat.gcd (stalk.gcd v).natAbs q - 1 := by
  calc
    (anchoredZeroColorStalkMultipliers v q stalk anchor).card ≤
        (zeroColorStalkMultipliers v q stalk).card :=
      Finset.card_le_card
        (anchoredZeroColorStalkMultipliers_subset
          v q stalk anchor hvanchor)
    _ ≤ Nat.gcd (stalk.gcd v).natAbs q - 1 :=
      zeroColorStalkMultipliers_card_le
        v q stalk hq hgcd top htop hvtop hwindow

/-! ## Axiom audit -/

#print axioms overlapDet_eq_zero_trans
#print axioms zeroColorStalk_complete_of_anchor
#print axioms zeroColorStalk_factorization_of_bezout
#print axioms zeroColorStalk_resonance_of_bezout
#print axioms zeroColorStalk_factorization
#print axioms zeroColorStalk_resonance
#print axioms zeroColorStalkMultipliers_subset_alignedResonant
#print axioms zeroColorStalkMultipliers_card_le
#print axioms anchoredZeroColorStalkMultipliers_card_le

end LRC14Concrete
end LonelyRunner
