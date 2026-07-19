/-
  TournamentH7.LRCSharpCombArithmetic -- arithmetic skeleton for THM-1094/1097.

  The measure-theoretic input is the sharp periodic-comb discrepancy

      |J intersect D_k| <= |J|/7 + 6/(49*k).

  This file kernel-checks the exact ratio inequalities that turn that input,
  a component-count bound, and the external finite core atlas into the
  two-comb and three-comb tails.  It also checks the arithmetic core of
  THM-1126's half-coverage/overlap reduction for the open four-comb case.  It
  checks the exact scale inequalities in THM-1128's thirteen-grid Kakeya
  rectangle gate, THM-1129's rectangle-to-needle tail arithmetic, and the
  three-range dispatch used by THM-1133.  It
  deliberately does not internalize the 66/220/495-core rational atlases or
  the interval-measure, torus-rectangle, and component-topology lemmas.

  Kernel-pure: no `sorry`, no `native_decide`, and no new axioms.
-/
import Mathlib

namespace LonelyRunner
namespace SharpCombArithmetic

/-- The uniform two-comb ratio inequality.  At the boundary `x=5/3`, its
strict margin is exactly `(122/3) * (r-1)`. -/
theorem twoComb_ratio_tail {x r : ℝ} (hx : (5 : ℝ) / 3 ≤ x) (hr : 1 < r) :
    6 * r + 29 < x * (28 * r - 7) := by
  have hcoef : 0 < 28 * r - 7 := by linarith
  have hbase : 6 * r + 29 < ((5 : ℝ) / 3) * (28 * r - 7) := by
    linarith
  have hmono : ((5 : ℝ) / 3) * (28 * r - 7) ≤ x * (28 * r - 7) :=
    mul_le_mul_of_nonneg_right hx hcoef.le
  exact hbase.trans_le hmono

/-- The uniform three-comb ratio inequality.  Here `r=k2/k1`, `s=k3/k1`,
and `q=k3/k2`, so `q*r=s`.  This is the analytic tail behind THM-1097. -/
theorem threeComb_ratio_tail {x r s q : ℝ} (hx : 7 ≤ x)
    (hr : 1 < r) (hrs : r < s) (hq : q * r = s) :
    6 * (s + q) + 37 < x * (21 * s - 7 * (1 + r)) := by
  have hr0 : 0 < r := lt_trans (by norm_num) hr
  have hs0 : 0 < s := lt_trans hr0 hrs
  have hq_eq : q = s / r := by
    exact (eq_div_iff hr0.ne').2 hq
  have hsr : s < s * r := by
    simpa only [mul_one] using mul_lt_mul_of_pos_left hr hs0
  have hq_lt_s : q < s := by
    rw [hq_eq]
    exact (div_lt_iff₀ hr0).2 hsr
  have hcoef : 0 < 21 * s - 7 * (1 + r) := by linarith
  have hbase : 6 * (s + q) + 37 < 7 * (21 * s - 7 * (1 + r)) := by
    linarith
  have hmono : 7 * (21 * s - 7 * (1 + r)) ≤
      x * (21 * s - 7 * (1 + r)) :=
    mul_le_mul_of_nonneg_right hx hcoef.le
  exact hbase.trans_le hmono

/-- The exceptional two-comb rectangle left by the coarse `5/3` tail is
strictly positive.  Integer killers satisfy `b>=a+1`; reals make the monotone
polynomial fact reusable. -/
theorem twoComb_exceptional_rectangle {a b : ℝ} (ha : 157 ≤ a)
    (hab : a + 1 ≤ b) :
    0 < 5 * a * b - 96 * b - 656 * a := by
  have hcoef : 0 ≤ 5 * a - 96 := by linarith
  have hstep :
      5 * a * (a + 1) - 96 * (a + 1) - 656 * a ≤
        5 * a * b - 96 * b - 656 * a := by
    nlinarith
  have hprod : 0 ≤ (a - 157) * (5 * a + 38) :=
    mul_nonneg (by linarith) (by linarith)
  have hdiag : 0 < 5 * a * (a + 1) - 96 * (a + 1) - 656 * a := by
    nlinarith
  exact hdiag.trans_le hstep

/-- Arithmetic form of the final-killer step.  Once an interval is longer
than `1/(7*k)`, the sharp discrepancy upper bound is strictly smaller than
the interval itself. -/
theorem sharpFinalComb_mass_lt {L k : ℝ} (hk : 0 < k)
    (hL : 1 / (7 * k) < L) :
    L / 7 + 6 / (49 * k) < L := by
  have hscale := mul_lt_mul_of_pos_left hL (show (0 : ℝ) < 6 / 7 by norm_num)
  have hfee : 6 / (49 * k) < 6 * L / 7 := by
    convert hscale using 1 <;> field_simp <;> ring
  linarith

/-- The asymptotic toothpick phase transition: after `m` removals, the density
to component-slope ratio beats the final sharp `1/7` comb threshold exactly
when `m<7/2`.  Hence the mass/component mechanism naturally reaches two and
three removals (`r=3,4`) and fails at four removals (`r=5`). -/
theorem toothpick_phase_transition {m : ℝ} (hm : 0 < m) :
    1 < (7 - m) / m ↔ m < (7 : ℝ) / 2 := by
  constructor
  · intro h
    have hraw : (1 : ℝ) * m < 7 - m :=
      (lt_div_iff₀ hm).mp (by simpa using h)
    have h' : m < 7 - m := by simpa using hraw
    linarith
  · intro h
    apply (lt_div_iff₀ hm).2
    linarith

/-- Arithmetic core of THM-1126's half-coverage gate.  The topological input
is that `b` internal danger components each cost at least `tau`, and that if
all safe components are short then their total mass is at most `(b+1)tau`.
Those two facts force the danger union to occupy at least `(ell-tau)/2`. -/
theorem halfCoverage_of_short_components {ell tau delta safeMass b : ℝ}
    (hdecomp : ell = safeMass + delta)
    (hcost : b * tau ≤ delta)
    (hshort : safeMass ≤ (b + 1) * tau) :
    (ell - tau) / 2 ≤ delta := by
  nlinarith

/-- Multi-component form of the same gate.  If the initial safe carrier has
`components` components, at most `b+components` short survivor pieces remain,
so failure of a long piece forces danger mass `(ell-components*tau)/2`. -/
theorem halfCoverage_of_short_components_multi
    {ell tau delta safeMass b components : ℝ}
    (hdecomp : ell = safeMass + delta)
    (hcost : b * tau ≤ delta)
    (hshort : safeMass ≤ (b + components) * tau) :
    (ell - components * tau) / 2 ≤ delta := by
  nlinarith

/-- The numerical last step of THM-1126's spanning-tree overlap certificate.
The four one-comb discrepancy bounds contribute `4ell/7 + endpointFee`; an
overlap tree larger than `ell/14 + endpointFee + tau/2` pushes the union below
the half-coverage obstruction. -/
theorem overlapTree_forces_below_half {ell tau endpointFee mst unionMass : ℝ}
    (hunion : unionMass ≤ 4 * ell / 7 + endpointFee - mst)
    (hoverlap : ell / 14 + endpointFee + tau / 2 < mst) :
    unionMass < (ell - tau) / 2 := by
  linarith

/-- A quadratic gap-energy certificate.  If the sum of squared component
lengths exceeds `tau` times their total length, at least one component is
strictly longer than `tau`. -/
theorem gapSquareEnergy_forces_long {ι : Type*} [Fintype ι]
    (length : ι → ℝ) {tau : ℝ} (hnonneg : ∀ i, 0 ≤ length i)
    (henergy : tau * ∑ i, length i < ∑ i, (length i) ^ 2) :
    ∃ i, tau < length i := by
  by_contra hlong
  simp only [not_exists, not_lt] at hlong
  have hterm : ∀ i, (length i) ^ 2 ≤ tau * length i := by
    intro i
    nlinarith [hnonneg i, hlong i]
  have hsum : (∑ i, (length i) ^ 2) ≤ ∑ i, tau * length i :=
    Finset.sum_le_sum fun i _ => hterm i
  rw [← Finset.mul_sum] at hsum
  exact (not_lt_of_ge hsum) henergy

/-- Discrete heart of THM-1128.  Four cyclic gaps on the thirteen-grid,
including zero gaps from repeated residues, cannot all have size at most
three. -/
theorem fourGridGaps_thirteen {g₀ g₁ g₂ g₃ : ℕ}
    (hsum : g₀ + g₁ + g₂ + g₃ = 13) :
    4 ≤ g₀ ∨ 4 ≤ g₁ ∨ 4 ≤ g₂ ∨ 4 ≤ g₃ := by
  omega

/-- Arithmetic closure for THM-1128's optimized thirteen-grid Kakeya
rectangle.  `M ≥ 24` protects the universal core window, `B ≥ 53M` pays one
full winding period, and `B ≤ k4` makes the resulting interval strictly
longer than the finest killer tooth. -/
theorem thirteenGrid_kakeya_scale {B k4 M : ℝ}
    (hB : 0 < B) (hBk : B ≤ k4) (hM24 : 24 ≤ M)
    (hscale : 53 * M ≤ B) :
    12 * (53 / (4914 * M)) ≤ (1 : ℝ) / 182 ∧
      (2809 : ℝ) / 2457 ≤ 53 * B / (2457 * M) ∧
      1 / (7 * k4) < 352 / (2457 * B) := by
  have hM : 0 < M := lt_of_lt_of_le (by norm_num) hM24
  have hk4 : 0 < k4 := lt_of_lt_of_le hB hBk
  constructor
  · calc
      12 * (53 / (4914 * M)) = (636 : ℝ) / (4914 * M) := by ring
      _ ≤ 1 / 182 := by
        apply (div_le_div_iff₀ (mul_pos (by norm_num) hM) (by norm_num)).2
        nlinarith
  constructor
  · apply (div_le_div_iff₀ (by norm_num : (0 : ℝ) < 2457)
      (mul_pos (by norm_num) hM)).2
    nlinarith
  · apply (div_lt_div_iff₀ (mul_pos (by norm_num) hk4)
      (mul_pos (by norm_num) hB)).2
    nlinarith

/-- Arithmetic part of THM-1129's rectangle-to-needle tail.  A time interval
of length at least one winding step plus one `1/7`-arc preimage can contain a
complete needle segment; positive offset span then makes that segment
strictly longer than the last killer's tooth.  The geometric containment is
the external input, while both scale comparisons are checked here. -/
theorem boundedOffset_needle_tail {d K c : ℝ}
    (hK : 0 < K) (hc : 0 < c) (hfit : 8 / (7 * K) ≤ d) :
    1 / K + 1 / (7 * K) ≤ d ∧
      1 / (7 * (K + c)) < 1 / (7 * K) := by
  constructor
  · calc
      1 / K + 1 / (7 * K) = 8 / (7 * K) := by field_simp; ring
      _ ≤ d := hfit
  · have hKc : 0 < K + c := by linarith
    exact one_div_lt_one_div_of_lt
      (mul_pos (by norm_num) hK)
      (by nlinarith : 7 * K < 7 * (K + c))

/-- Logical assembly used by THM-1133.  Once external certificates own the
bottom, exact middle, and tail ranges, there is no untracked integer row.
Keeping this tiny dispatch in the kernel prevents endpoint-language drift
between the finite ledgers and the all-scale theorem statement. -/
theorem threeRange_dispatch {Good : ℕ → Prop} {lo cut tail : ℕ}
    (hbottom : ∀ K, lo ≤ K → K < cut → Good K)
    (hmiddle : ∀ K, lo ≤ K → cut ≤ K → K < tail → Good K)
    (htail : ∀ K, lo ≤ K → tail ≤ K → Good K) :
    ∀ K, lo ≤ K → Good K := by
  intro K hlegal
  by_cases hcut : K < cut
  · exact hbottom K hlegal hcut
  by_cases htailK : K < tail
  · exact hmiddle K hlegal (Nat.le_of_not_gt hcut) htailK
  · exact htail K hlegal (Nat.le_of_not_gt htailK)

/-- Arithmetic closure for THM-1134's five-residue multiplier cone.  The
external ten-orbit lemma supplies a six-unit grid gap; these inequalities
check core stability, one complete needle winding, and the strict final-tooth
margin under `M≥80` and `B≥17M`. -/
theorem fiveGrid_kakeya_scale {B k5 M : ℝ}
    (hB : 0 < B) (hBk : B ≤ k5) (hM80 : 80 ≤ M)
    (hscale : 17 * M ≤ B) :
    12 * (10 / (273 * M)) ≤ (1 : ℝ) / 182 ∧
      1 + (67 : ℝ) / 273 ≤ B * (20 / (273 * M)) ∧
      1 / (7 * k5) < 67 / (273 * B) := by
  have hM : 0 < M := lt_of_lt_of_le (by norm_num) hM80
  have hk5 : 0 < k5 := lt_of_lt_of_le hB hBk
  constructor
  · calc
      12 * (10 / (273 * M)) = (120 : ℝ) / (273 * M) := by ring
      _ ≤ 1 / 182 := by
        apply (div_le_div_iff₀ (mul_pos (by norm_num) hM) (by norm_num)).2
        nlinarith
  constructor
  · calc
      1 + (67 : ℝ) / 273 = (340 : ℝ) / 273 := by norm_num
      _ ≤ (20 * B) / (273 * M) := by
        apply (div_le_div_iff₀ (by norm_num : (0 : ℝ) < 273)
          (mul_pos (by norm_num) hM)).2
        nlinarith
      _ = B * (20 / (273 * M)) := by ring
  · apply (div_lt_div_iff₀ (mul_pos (by norm_num) hk5)
      (mul_pos (by norm_num) hB)).2
    nlinarith

/-- Exact scale arithmetic for THM-1134's former worst step-two ray. -/
theorem r6_stepTwo_ray_scale {b : ℝ} (hb : 148 ≤ b) :
    1 + (8 : ℝ) / 49 ≤ (b + 4) * (3 / 392) ∧
      1 / (7 * (b + 8)) < 8 / (49 * (b + 4)) := by
  have hb4 : 0 < b + 4 := by linarith
  have hb8 : 0 < b + 8 := by linarith
  constructor
  · nlinarith
  · apply (div_lt_div_iff₀ (mul_pos (by norm_num) hb8)
      (mul_pos (by norm_num) hb4)).2
    nlinarith

/-- Arithmetic heart of THM-1134's separated-ratio `Q5` gate.  External
interval geometry supplies the mass lower bound and component upper bound;
strict positivity of `Q5` makes seven finest-tooth lengths times the survivor
mass exceed the component count. -/
theorem fiveComb_Q5_mass_components
    {ell K sumK reciprocal survivorMass components : ℝ}
    (hK : 0 < K)
    (hmass : 2 * ell / 7 - 6 * reciprocal / 49 ≤ survivorMass)
    (hcomponents : components ≤ ell * sumK + 47 / 7)
    (hQ5 : 0 < ell * (14 * K - 7 * sumK) -
      6 * K * reciprocal - 47) :
    components < 7 * K * survivorMass := by
  have hscaled := mul_le_mul_of_nonneg_left hmass
    (mul_nonneg (by norm_num : (0 : ℝ) ≤ 7) hK.le)
  nlinarith

#print axioms twoComb_ratio_tail
#print axioms threeComb_ratio_tail
#print axioms twoComb_exceptional_rectangle
#print axioms sharpFinalComb_mass_lt
#print axioms toothpick_phase_transition
#print axioms halfCoverage_of_short_components
#print axioms halfCoverage_of_short_components_multi
#print axioms overlapTree_forces_below_half
#print axioms gapSquareEnergy_forces_long
#print axioms fourGridGaps_thirteen
#print axioms thirteenGrid_kakeya_scale
#print axioms boundedOffset_needle_tail
#print axioms threeRange_dispatch
#print axioms fiveGrid_kakeya_scale
#print axioms r6_stepTwo_ray_scale
#print axioms fiveComb_Q5_mass_components

end SharpCombArithmetic
end LonelyRunner
