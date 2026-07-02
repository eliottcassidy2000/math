/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-02-S34)
-/
import Mathlib.MeasureTheory.Integral.IntervalIntegral.Periodic
import Mathlib.Tactic.NormNum.Prime
import Mathlib.MeasureTheory.Group.AddCircle
import Mathlib.MeasureTheory.Measure.Haar.Unique
import TournamentH7.LonelyRunnerMathlib

/-!
# The 7-commensuration lemma (project THM-602 addendum; DAG-spec row 6)

For speeds `P, Q` with `7 ∤ P` and `7 ∣ Q`, at the critical radius `r = 1/14` the overlap of
the two danger combs is **identically** `(2r)² = 1/49`, for every pair of phases — exact
independence.  This is the finite-support fact behind THM-602(C)'s deviation-table tail:
7-commensurate rows contribute zero deviation, so the partial-cycle remainder sum is a finite
table.

**Proof.**  The seven translates of `danger P ψ` by `i/7` are the preimages, under the
measure-preserving map `x ↦ P • x + ψ`, of the seven *open* balls of radius `1/14` centered at
`P·i/7` — pairwise disjoint (centers at distinct sevenths, spacing `≥ 1/7 = ` sum of radii) with
total measure `7 · (1/7) = 1`, so their union is almost all of the circle.  `danger Q φ` is
invariant under translation by `1/7` because `7 ∣ Q`.  Hence the seven translates of
`danger P ψ ∩ danger Q φ` are pairwise disjoint, all of the same measure (translation
invariance of Haar measure), and their union is `danger Q φ` up to a null set; so each has
measure `vol (danger Q φ) / 7 = (2r)/7 = (2r)²`.
-/

open MeasureTheory MeasureTheory.Measure Set Metric
open scoped ENNReal

namespace LRCCommensuration

/-- The affine runner map `x ↦ v • x + φ` on the unit circle. -/
def runnerMap (v : ℤ) (φ : UnitAddCircle) : UnitAddCircle → UnitAddCircle :=
  fun x => v • x + φ

/-- The `r`-danger set of speed `v` at phase `φ` on the unit circle: the times `x` at which
the runner of speed `v`, offset by `φ`, is within `r` of the origin. -/
def danger (v : ℤ) (φ : UnitAddCircle) (r : ℝ) : Set UnitAddCircle :=
  runnerMap v φ ⁻¹' (ball 0 r)

theorem continuous_runnerMap (v : ℤ) (φ : UnitAddCircle) : Continuous (runnerMap v φ) :=
  (continuous_zsmul v).add continuous_const

theorem measurableSet_danger (v : ℤ) (φ : UnitAddCircle) (r : ℝ) :
    MeasurableSet (danger v φ r) :=
  (isOpen_ball.preimage (continuous_runnerMap v φ)).measurableSet

theorem mem_danger {v : ℤ} {φ : UnitAddCircle} {r : ℝ} {x : UnitAddCircle} :
    x ∈ danger v φ r ↔ v • x + φ ∈ ball (0 : UnitAddCircle) r := Iff.rfl

/-- The runner map is measure-preserving for `v ≠ 0` (compact abelian group: `zsmul` is
measure-preserving by surjectivity + Haar uniqueness, and translation is Haar). -/
theorem measurePreserving_runnerMap {v : ℤ} (hv : v ≠ 0) (φ : UnitAddCircle) :
    MeasurePreserving (runnerMap v φ) volume volume := by
  have h1 : MeasurePreserving (fun x : UnitAddCircle => v • x) volume volume :=
    measurePreserving_zsmul volume hv
  have h2 : MeasurePreserving (fun x : UnitAddCircle => x + φ) volume volume :=
    measurePreserving_add_right volume φ
  exact h2.comp h1

/-- The volume of an open ball in the unit circle, subcritical case. -/
theorem volume_ball_unitAddCircle (x : UnitAddCircle) {r : ℝ} (hr : 0 < r) (hr2 : 2 * r ≤ 1) :
    volume (ball x r) = ENNReal.ofReal (2 * r) := by
  rw [← measure_congr (AddCircle.closedBall_ae_eq_ball (x := x) (ε := r)),
    AddCircle.volume_closedBall, min_eq_right hr2]

/-- The volume of a danger set: each comb has total measure exactly `2r`. -/
theorem volume_danger {v : ℤ} (hv : v ≠ 0) (φ : UnitAddCircle) {r : ℝ}
    (hr : 0 < r) (hr2 : 2 * r ≤ 1) :
    volume (danger v φ r) = ENNReal.ofReal (2 * r) := by
  rw [danger, (measurePreserving_runnerMap hv φ).measure_preimage
    measurableSet_ball.nullMeasurableSet, volume_ball_unitAddCircle _ hr hr2]

section SevenCommensuration

/-- The seventh-turn on the unit circle. -/
noncomputable def seventh : UnitAddCircle := ((1 / 7 : ℝ) : UnitAddCircle)

/-- Distinct multiples `k • (1/7)`, `0 < k < 7`, are nonzero — in fact at distance `≥ 1/7`
from `0`.  Stated with the norm for direct use in ball-disjointness. -/
theorem norm_nsmul_seventh_ge {k : ℕ} (hk0 : 0 < k) (hk7 : k < 7) :
    (1 : ℝ) / 7 ≤ ‖(k • seventh : UnitAddCircle)‖ := by
  have hval : (k • ((1 : ℝ) / 7) : ℝ) = ((k : ℝ)) / ((7 : ℕ) : ℝ) := by
    rw [nsmul_eq_mul]
    push_cast
    ring
  have hcast : (k • seventh : UnitAddCircle) = ((((k : ℝ)) / ((7 : ℕ) : ℝ) : ℝ) : UnitAddCircle) := by
    unfold seventh
    rw [← hval]
    rfl
  rw [hcast, LonelyRunner.norm_natCast_div]
  have hmod : k % 7 = k := Nat.mod_eq_of_lt hk7
  rw [hmod]
  have h1 : 1 ≤ min k (7 - k) := by omega
  have h1' : (1 : ℝ) ≤ ((min k (7 - k) : ℕ) : ℝ) := by exact_mod_cast h1
  have h77 : ((7 : ℕ) : ℝ) = 7 := by norm_num
  rw [h77]
  linarith

/-- The full turn: `7 • (1/7) = 1 = 0` on the unit circle. -/
theorem seven_nsmul_seventh : (7 : ℕ) • seventh = 0 := by
  unfold seventh
  have h1 : ((7 : ℕ) • ((1 : ℝ) / 7) : ℝ) = 1 := by rw [nsmul_eq_mul]; norm_num
  calc (7 : ℕ) • (((1 : ℝ) / 7 : ℝ) : UnitAddCircle)
      = (((7 : ℕ) • ((1 : ℝ) / 7) : ℝ) : UnitAddCircle) := rfl
    _ = ((1 : ℝ) : UnitAddCircle) := by rw [h1]
    _ = 0 := by rw [AddCircle.coe_period]

theorem seven_zsmul_seventh : (7 : ℤ) • seventh = 0 := by
  have := seven_nsmul_seventh
  rwa [← natCast_zsmul] at this

/-- Reduction of integer multiples of the seventh-turn mod 7. -/
theorem zsmul_seventh_emod (m : ℤ) : m • seventh = (m % 7) • seventh := by
  conv_lhs => rw [← Int.emod_add_mul_ediv m 7]
  rw [add_smul, mul_comm, mul_smul, seven_zsmul_seventh, smul_zero, add_zero]

/-- Nonzero residues give points at distance at least `1/7` from the origin. -/
theorem norm_zsmul_seventh_ge {m : ℤ} (hm : ¬ (7 : ℤ) ∣ m) :
    (1 : ℝ) / 7 ≤ ‖(m • seventh : UnitAddCircle)‖ := by
  rw [zsmul_seventh_emod]
  have h0 : 0 ≤ m % 7 := Int.emod_nonneg m (by norm_num)
  have h7 : m % 7 < 7 := Int.emod_lt_of_pos m (by norm_num)
  have hne : m % 7 ≠ 0 := fun h => hm (Int.dvd_of_emod_eq_zero h)
  set k : ℕ := (m % 7).toNat with hk
  have hcast : ((k : ℕ) : ℤ) = m % 7 := Int.toNat_of_nonneg h0
  have hsm : (m % 7) • seventh = (k : ℕ) • seventh := by
    rw [← hcast, natCast_zsmul]
  rw [hsm]
  exact norm_nsmul_seventh_ge (by omega) (by omega)

section MainTheorem

/-- The `i`-th seventh-turn translate. -/
noncomputable def tau (i : ℕ) : UnitAddCircle := (i : ℤ) • seventh

/-- Translating a danger set shifts its phase. -/
theorem preimage_add_danger (v : ℤ) (χ τ' : UnitAddCircle) (r : ℝ) :
    (fun x => x + τ') ⁻¹' danger v χ r = danger v (χ + v • τ') r := by
  ext x
  have harg : v • (x + τ') + χ = v • x + (χ + v • τ') := by rw [smul_add]; abel
  simp only [mem_preimage, danger, runnerMap, mem_ball, harg]

/-- A phase-shifted danger set is the preimage of a shifted ball. -/
theorem danger_eq_preimage_ball (v : ℤ) (χ τ' : UnitAddCircle) (r : ℝ) :
    danger v (χ + τ') r = runnerMap v χ ⁻¹' ball (-τ') r := by
  ext x
  have harg : v • x + (χ + τ') = (v • x + χ) - (-τ') := by abel
  simp only [danger, runnerMap, mem_preimage, mem_ball, dist_eq_norm, sub_zero, harg]

/-- Centers of distinct translates of the `P`-comb are at distance `≥ 1/7`. -/
theorem center_sep {P : ℤ} (hP7 : ¬ (7 : ℤ) ∣ P) {i j : ℕ} (hij : i < j) (hj7 : j < 7) :
    (1 : ℝ) / 7 ≤ dist (-(P • tau i)) (-(P • tau j)) := by
  rw [dist_eq_norm]
  have hval : -(P • tau i) - -(P • tau j) = (P * ((j : ℤ) - (i : ℤ))) • seventh := by
    unfold tau
    rw [smul_smul, smul_smul, neg_sub_neg, ← sub_smul]
    congr 1
    ring
  rw [hval]
  apply norm_zsmul_seventh_ge
  intro hdvd
  have h7prime : Prime (7 : ℤ) := by norm_num
  rcases h7prime.dvd_mul.mp hdvd with h | h
  · exact hP7 h
  · have h1 : (i : ℤ) < (j : ℤ) := by exact_mod_cast hij
    have h2 : (j : ℤ) < 7 := by exact_mod_cast hj7
    omega

/-- The seven translate-balls are pairwise disjoint. -/
theorem balls_disjoint {P : ℤ} (hP7 : ¬ (7 : ℤ) ∣ P) {i j : ℕ} (hi7 : i < 7) (hj7 : j < 7)
    (hij : i ≠ j) :
    Disjoint (ball (-(P • tau i)) (1 / 14 : ℝ)) (ball (-(P • tau j)) (1 / 14 : ℝ)) := by
  wlog hlt : i < j generalizing i j
  · exact (this hj7 hi7 hij.symm (by omega)).symm
  rw [Set.disjoint_left]
  intro z hzi hzj
  rw [mem_ball] at hzi hzj
  have hsep := center_sep hP7 hlt hj7
  have : dist (-(P • tau i)) (-(P • tau j)) < 1 / 7 :=
    calc dist (-(P • tau i)) (-(P • tau j))
        ≤ dist (-(P • tau i)) z + dist z (-(P • tau j)) := dist_triangle _ _ _
      _ < 1 / 14 + 1 / 14 := by rw [dist_comm (-(P • tau i)) z]; exact add_lt_add hzi hzj
      _ = 1 / 7 := by norm_num
  linarith

/-- **The 7-commensuration lemma** (project THM-602 addendum; DAG-spec row 6).  At the
critical radius `r = 1/14`, a speed `P` not divisible by `7` and a speed `Q` divisible by `7`
have danger-overlap exactly `(2r)² = 1/49`, for *every* pair of phases: 7-commensurate pairs
are exactly independent.  The seven `1/7`-translates of the `P`-comb tile the circle a.e.,
while the `Q`-comb is `1/7`-periodic; averaging over the translates gives the product. -/
theorem seven_commensuration {P Q : ℤ} (hP7 : ¬ (7 : ℤ) ∣ P) (hQ0 : Q ≠ 0) (hQ7 : (7 : ℤ) ∣ Q)
    (ψ φ : UnitAddCircle) :
    volume (danger P ψ (1 / 14) ∩ danger Q φ (1 / 14)) = ENNReal.ofReal (1 / 49) := by
  have hP0 : P ≠ 0 := fun h => hP7 (h ▸ dvd_zero 7)
  have hrpos : (0 : ℝ) < 1 / 14 := by norm_num
  have hr2 : 2 * (1 / 14 : ℝ) ≤ 1 := by norm_num
  -- (1) the Q-comb is 1/7-periodic
  have hCper : ∀ i : ℕ, (fun x => x + tau i) ⁻¹' danger Q φ (1 / 14) = danger Q φ (1 / 14) := by
    intro i
    rw [preimage_add_danger]
    have hQτ : Q • tau i = 0 := by
      obtain ⟨m, rfl⟩ := hQ7
      unfold tau
      rw [smul_smul, show (7 * m * (i : ℤ)) = (m * (i : ℤ)) * 7 by ring, ← smul_smul,
        seven_zsmul_seventh, smul_zero]
    rw [hQτ, add_zero]
  -- (2) translating carries the intersection to the shifted-phase piece
  have htrans : ∀ i : ℕ, (fun x => x + tau i) ⁻¹' (danger P ψ (1 / 14) ∩ danger Q φ (1 / 14)) =
      danger P (ψ + P • tau i) (1 / 14) ∩ danger Q φ (1 / 14) := by
    intro i
    rw [Set.preimage_inter, hCper, preimage_add_danger]
  -- (3) so every piece has the intersection's volume
  have hpiece : ∀ i : ℕ,
      volume (danger P (ψ + P • tau i) (1 / 14) ∩ danger Q φ (1 / 14)) =
      volume (danger P ψ (1 / 14) ∩ danger Q φ (1 / 14)) := by
    intro i
    rw [← htrans i]
    exact (measurePreserving_add_right volume (tau i)).measure_preimage
      (((measurableSet_danger P ψ _).inter (measurableSet_danger Q φ _)).nullMeasurableSet)
  -- (4) the union of the seven translate-pieces has full measure
  have hpre : ∀ i : ℕ, danger P (ψ + P • tau i) (1 / 14) =
      runnerMap P ψ ⁻¹' ball (-(P • tau i)) (1 / 14) := fun i =>
    danger_eq_preimage_ball P ψ (P • tau i) (1 / 14)
  have hUvol : volume (⋃ i ∈ Finset.range 7, danger P (ψ + P • tau i) (1 / 14)) = 1 := by
    have hrw : (⋃ i ∈ Finset.range 7, danger P (ψ + P • tau i) (1 / 14)) =
        runnerMap P ψ ⁻¹' (⋃ i ∈ Finset.range 7, ball (-(P • tau i)) (1 / 14)) := by
      simp only [hpre, Set.preimage_iUnion₂]
    rw [hrw, (measurePreserving_runnerMap hP0 ψ).measure_preimage]
    · rw [measure_biUnion_finset]
      · have hone : ∀ i ∈ Finset.range 7,
            volume (ball (-(P • tau i)) (1 / 14 : ℝ)) = ENNReal.ofReal (2 * (1 / 14)) :=
          fun i _ => volume_ball_unitAddCircle _ hrpos hr2
        rw [Finset.sum_congr rfl hone, Finset.sum_const, Finset.card_range, nsmul_eq_mul]
        calc ((7 : ℕ) : ℝ≥0∞) * ENNReal.ofReal (2 * (1 / 14))
            = ENNReal.ofReal 7 * ENNReal.ofReal (2 * (1 / 14)) := by norm_num
          _ = ENNReal.ofReal (7 * (2 * (1 / 14))) := (ENNReal.ofReal_mul (by norm_num)).symm
          _ = 1 := by norm_num
      · intro i hi j hj hij
        exact balls_disjoint hP7 (Finset.mem_range.mp hi) (Finset.mem_range.mp hj) hij
      · exact fun i _ => measurableSet_ball
    · exact (Finset.measurableSet_biUnion _ fun i _ => measurableSet_ball).nullMeasurableSet
  -- (5) the C-mass splits equally over the seven pieces
  have hCmeas : MeasurableSet (danger Q φ (1 / 14)) := measurableSet_danger Q φ _
  have hUmeas : MeasurableSet (⋃ i ∈ Finset.range 7, danger P (ψ + P • tau i) (1 / 14)) :=
    Finset.measurableSet_biUnion _ fun i _ => measurableSet_danger P _ _
  have hUc : volume (⋃ i ∈ Finset.range 7, danger P (ψ + P • tau i) (1 / 14))ᶜ = 0 := by
    rw [measure_compl hUmeas (measure_ne_top _ _), hUvol]
    have huniv : volume (Set.univ : Set UnitAddCircle) = 1 := by
      rw [AddCircle.measure_univ, ENNReal.ofReal_one]
    rw [huniv, tsub_self]
  have hsplit : volume (danger Q φ (1 / 14) ∩ ⋃ i ∈ Finset.range 7,
        danger P (ψ + P • tau i) (1 / 14)) = ENNReal.ofReal (2 * (1 / 14)) := by
    have h1 : volume (danger Q φ (1 / 14) ∩ ⋃ i ∈ Finset.range 7,
          danger P (ψ + P • tau i) (1 / 14)) +
        volume (danger Q φ (1 / 14) \ ⋃ i ∈ Finset.range 7,
          danger P (ψ + P • tau i) (1 / 14)) = volume (danger Q φ (1 / 14)) :=
      measure_inter_add_diff _ hUmeas
    have h2 : volume (danger Q φ (1 / 14) \ ⋃ i ∈ Finset.range 7,
        danger P (ψ + P • tau i) (1 / 14)) = 0 :=
      measure_mono_null (fun x hx => hx.2) hUc
    rw [h2, add_zero] at h1
    rw [h1, volume_danger hQ0 φ hrpos hr2]
  have hsum : volume (danger Q φ (1 / 14) ∩ ⋃ i ∈ Finset.range 7,
      danger P (ψ + P • tau i) (1 / 14)) =
      (7 : ℕ) • volume (danger P ψ (1 / 14) ∩ danger Q φ (1 / 14)) := by
    rw [Set.inter_iUnion₂, measure_biUnion_finset]
    · have heach : ∀ i ∈ Finset.range 7,
          volume (danger Q φ (1 / 14) ∩ danger P (ψ + P • tau i) (1 / 14)) =
          volume (danger P ψ (1 / 14) ∩ danger Q φ (1 / 14)) := by
        intro i _
        rw [Set.inter_comm]
        exact hpiece i
      rw [Finset.sum_congr rfl heach, Finset.sum_const, Finset.card_range]
    · intro i hi j hj hij
      have hd : Disjoint (danger P (ψ + P • tau i) (1 / 14))
          (danger P (ψ + P • tau j) (1 / 14)) := by
        rw [hpre i, hpre j]
        exact ((balls_disjoint hP7 (Finset.mem_range.mp hi) (Finset.mem_range.mp hj)
          hij).preimage _)
      exact hd.mono Set.inter_subset_right Set.inter_subset_right
    · exact fun i _ => hCmeas.inter (measurableSet_danger P _ _)
  -- (6) endgame: 7 · V = 1/7 ⟹ V = 1/49
  have hkey : ((7 : ℕ) : ℝ≥0∞) * volume (danger P ψ (1 / 14) ∩ danger Q φ (1 / 14)) =
      ENNReal.ofReal (2 * (1 / 14)) := by
    rw [← nsmul_eq_mul, ← hsum, hsplit]
  have hval : ((7 : ℕ) : ℝ≥0∞) * ENNReal.ofReal (1 / 49) = ENNReal.ofReal (2 * (1 / 14)) := by
    calc ((7 : ℕ) : ℝ≥0∞) * ENNReal.ofReal (1 / 49)
        = ENNReal.ofReal 7 * ENNReal.ofReal (1 / 49) := by norm_num
      _ = ENNReal.ofReal (7 * (1 / 49)) := (ENNReal.ofReal_mul (by norm_num)).symm
      _ = ENNReal.ofReal (2 * (1 / 14)) := by norm_num
  rw [← hval] at hkey
  exact (ENNReal.mul_right_inj (by norm_num) (by norm_num)).mp hkey

/-- Symmetric form: the commensurate speed on the left. -/
theorem seven_commensuration' {P Q : ℤ} (hP0 : P ≠ 0) (hP7 : (7 : ℤ) ∣ P) (hQ7 : ¬ (7 : ℤ) ∣ Q)
    (ψ φ : UnitAddCircle) :
    volume (danger P ψ (1 / 14) ∩ danger Q φ (1 / 14)) = ENNReal.ofReal (1 / 49) := by
  rw [Set.inter_comm]
  exact seven_commensuration hQ7 hP0 hP7 φ ψ

end MainTheorem

end SevenCommensuration

end LRCCommensuration
