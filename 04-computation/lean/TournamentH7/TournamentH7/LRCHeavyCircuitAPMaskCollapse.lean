import Mathlib
import TournamentH7.LRCContinuumTriangleCeiling

/-!
# Heavy continuum circuits collapse the beat-mask nerve (THM-1218)

THM-1203's analytic bridge says that continuum BAD mass above `60/637`
forces each of the four deletion triangles into the ratio set
`p=q`, `p=2q`, or `q=2p`.  Its exact finite core and sheared tail are
kernel-audited here at that stricter cutoff.  The remaining four-obligation
calculation forces equal consecutive gaps, so an ordered heavy quadruple is
a four-term arithmetic progression.

At the defining beat `q=b6-b5`, that progression makes all four speeds
congruent modulo `q`.  Their direct danger masks, and their pointwise
restrictions to every numerator range, are therefore identical.  In the
mixed-period application the actual master length is the divisor `L=q/d0`;
arbitrary longer ranges are not asserted to be master clocks.  The usual five
mask labels reduce to at most three distinct masks and the common-period
threshold becomes

`B3(Q) = 2 + 3*(A(Q)-1) = 6*ceil(Q/14)-4`.

The Haar-measure identification, BAD-to-deletion-event inclusions, and gcd
dilation invariance are represented by `HeavyTriangleProvider` and four
explicit domination hypotheses.  Everything after those analytic provider
hypotheses is kernel-checked without placeholders or a native evaluator.
-/

namespace LRC14.HeavyCircuitAPMaskCollapse

open scoped Classical
noncomputable section

/-! ## The strict `60/637` spectral gate -/

/-- The endpoint value of THM-1203's sheared-grid tail at total frequency
`p+q=26`.  Strictly larger mass is the heavy-circuit predicate. -/
def heavyCutoff : ℚ := 60 / 637

theorem ap_mass_above_cutoff :
    heavyCutoff < (2 : ℚ) / 21 ∧ (2 : ℚ) / 21 - heavyCutoff = 2 / 1911 := by
  norm_num [heavyCutoff]

theorem second_stratum_below_cutoff :
    (13 : ℚ) / 147 < heavyCutoff ∧ heavyCutoff - 13 / 147 = 11 / 1911 := by
  norm_num [heavyCutoff]

theorem tail_endpoint_exact :
    (3 : ℚ) / 49 + 6 / (7 * 26) = heavyCutoff := by
  norm_num [heavyCutoff]

/-- The analytic tail is at most the heavy cutoff for every `N>=26`.
At `N=26` equality is allowed, which is why THM-1218 uses strict `>` rather
than `>=`. -/
theorem tail_le_heavyCutoff (N : ℚ) (hN : 26 ≤ N) :
    (3 : ℚ) / 49 + 6 / (7 * N) ≤ heavyCutoff := by
  have hden : (0 : ℚ) < 7 * N := by positivity
  have hterm : (6 : ℚ) / (7 * N) ≤ 3 / 91 := by
    apply (div_le_iff₀ hden).2
    nlinarith
  calc
    (3 : ℚ) / 49 + 6 / (7 * N) ≤ 3 / 49 + 3 / 91 :=
      add_le_add_right hterm _
    _ = heavyCutoff := by norm_num [heavyCutoff]

/-- Denominator-cleared strict comparison for THM-1203's exact carry
numerator.  Since

`triangleMeasure = 2*S/[7*(p+q)*p*q]`,

the inequality `60/637 < triangleMeasure` is exactly `30D < 91S`. -/
def triangleHeavyAtCutoff (p q : ℕ) : Prop :=
  30 * ((p + q) * p * q) <
    91 * LonelyRunner.LRC14.ContinuumTriangleCeiling.triangleNumerator p q

instance (p q : ℕ) : Decidable (triangleHeavyAtCutoff p q) :=
  inferInstanceAs (Decidable
    (30 * ((p + q) * p * q) <
      91 * LonelyRunner.LRC14.ContinuumTriangleCeiling.triangleNumerator p q))

private theorem cleared_strict_iff (numerator denominator : ℕ)
    (hdenominator : 0 < denominator) :
    30 * denominator < 91 * numerator ↔
      (60 : ℚ) / 637 <
        (((2 * numerator : ℕ) : ℚ) / ((7 * denominator : ℕ) : ℚ)) := by
  have hdenNat : 0 < 7 * denominator := Nat.mul_pos (by norm_num) hdenominator
  have hden : (0 : ℚ) < ((7 * denominator : ℕ) : ℚ) := by
    exact_mod_cast hdenNat
  rw [lt_div_iff₀ hden]
  constructor
  · intro h
    have hcast : (((30 * denominator : ℕ) : ℚ)) <
        (((91 * numerator : ℕ) : ℚ)) := by
      exact_mod_cast h
    push_cast at hcast ⊢
    norm_num at hcast ⊢
    nlinarith
  · intro h
    have hcast : (((30 * denominator : ℕ) : ℚ)) <
        (((91 * numerator : ℕ) : ℚ)) := by
      push_cast at h ⊢
      norm_num at h ⊢
      nlinarith
    exact_mod_cast hcast

theorem triangleHeavyAtCutoff_iff_measure_gt
    (p q : ℕ) (hp : 0 < p) (hq : 0 < q) :
    triangleHeavyAtCutoff p q ↔ heavyCutoff <
      LonelyRunner.LRC14.ContinuumTriangleCeiling.triangleMeasure p q := by
  have hden : 0 < (p + q) * p * q := by positivity
  simpa only [triangleHeavyAtCutoff, heavyCutoff,
      LonelyRunner.LRC14.ContinuumTriangleCeiling.triangleMeasure,
      Nat.mul_assoc] using cleared_strict_iff
    (LonelyRunner.LRC14.ContinuumTriangleCeiling.triangleNumerator p q)
    ((p + q) * p * q) hden

/-- The exact reduced finite core pairs strictly above `60/637`. -/
def coreHeavyPairs : Finset (ℕ × ℕ) :=
  LonelyRunner.LRC14.ContinuumTriangleCeiling.corePairs.filter fun pair =>
    triangleHeavyAtCutoff pair.1 pair.2

/-- The unique reduced heavy triangle is `(1,2)`.  Ordinary kernel `decide`
checks all 99 pairs; no native evaluator is used. -/
theorem coreHeavyPairs_exact : coreHeavyPairs = {(1, 2)} := by
  set_option maxRecDepth 100000 in
    decide

theorem core_heavy_forces_ratio12 {p q : ℕ}
    (hmem : (p, q) ∈
      LonelyRunner.LRC14.ContinuumTriangleCeiling.corePairs)
    (hheavy : heavyCutoff <
      LonelyRunner.LRC14.ContinuumTriangleCeiling.triangleMeasure p q) :
    LonelyRunner.LRC14.ContinuumTriangleCeiling.Ratio12 p q := by
  have hshape := hmem
  rw [LonelyRunner.LRC14.ContinuumTriangleCeiling.corePairs] at hshape
  have hproduct := (Finset.mem_filter.mp hshape).1
  have hpBounds := Finset.mem_Icc.mp (Finset.mem_product.mp hproduct).1
  have hqBounds := Finset.mem_Icc.mp (Finset.mem_product.mp hproduct).2
  have hp : 0 < p := by omega
  have hq : 0 < q := by omega
  have hheavyNat : triangleHeavyAtCutoff p q :=
    (triangleHeavyAtCutoff_iff_measure_gt p q hp hq).2 hheavy
  have hpq : (p, q) ∈ coreHeavyPairs := Finset.mem_filter.mpr ⟨hmem, hheavyNat⟩
  rw [coreHeavyPairs_exact] at hpq
  have heq : (p, q) = (1, 2) := by simpa using hpq
  have hp1 : p = 1 := congrArg Prod.fst heq
  have hq2 : q = 2 := congrArg Prod.snd heq
  unfold LonelyRunner.LRC14.ContinuumTriangleCeiling.Ratio12
  omega

/-! ## The four-deletion provider and AP rigidity -/

/-- The precise analytic service supplied by THM-1203 after gcd reduction:
an additive-triangle mass strictly above the heavy cutoff forces ratio
`1:1`, `1:2`, or `2:1`. -/
structure HeavyTriangleProvider (J : ℕ → ℕ → ℚ) : Prop where
  heavy_forces_ratio12 : ∀ p q : ℕ,
    0 < p → 0 < q → heavyCutoff < J p q →
      LonelyRunner.LRC14.ContinuumTriangleCeiling.Ratio12 p q

/-- Once BAD is dominated by its four deletion-triangle events, strict
heavy mass forces all three consecutive gaps to agree. -/
theorem heavy_four_forces_equal_gaps
    (J : ℕ → ℕ → ℚ) (provider : HeavyTriangleProvider J)
    (μ : ℚ) (a b c : ℕ)
    (ha : 0 < a) (hb : 0 < b) (hc : 0 < c)
    (hheavy : heavyCutoff < μ)
    (hab : μ ≤ J a b) (hbc : μ ≤ J b c)
    (ha_bc : μ ≤ J a (b + c)) (hab_c : μ ≤ J (a + b) c) :
    a = b ∧ b = c := by
  apply LonelyRunner.LRC14.ContinuumTriangleCeiling.ratio12_four_triangles_rigid
    a b c ha hb hc
  · exact provider.heavy_forces_ratio12 a b ha hb (lt_of_lt_of_le hheavy hab)
  · exact provider.heavy_forces_ratio12 b c hb hc (lt_of_lt_of_le hheavy hbc)
  · exact provider.heavy_forces_ratio12 a (b + c) ha (by omega)
      (lt_of_lt_of_le hheavy ha_bc)
  · exact provider.heavy_forces_ratio12 (a + b) c (by omega) hc
      (lt_of_lt_of_le hheavy hab_c)

/-- Proof-facing form for the ordered quadruple `(bi,bj,b5,b6)`.  The four
domination hypotheses are exactly the BAD-to-three-band inclusions for its
four deletion triangles. -/
theorem heavy_ordered_quadruple_forces_AP
    (J : ℕ → ℕ → ℚ) (provider : HeavyTriangleProvider J)
    (μ : ℚ) (bi bj b5 b6 : ℕ)
    (hij : bi < bj) (hj5 : bj < b5) (h56 : b5 < b6)
    (hheavy : heavyCutoff < μ)
    (h01 : μ ≤ J (bj - bi) (b5 - bj))
    (h12 : μ ≤ J (b5 - bj) (b6 - b5))
    (h02 : μ ≤ J (bj - bi) ((b5 - bj) + (b6 - b5)))
    (h13 : μ ≤ J ((bj - bi) + (b5 - bj)) (b6 - b5)) :
    let q := b6 - b5
    bi + 2 * q = b5 ∧ bj + q = b5 ∧ b6 = b5 + q := by
  have ha : 0 < bj - bi := by omega
  have hb : 0 < b5 - bj := by omega
  have hc : 0 < b6 - b5 := by omega
  have equalGaps := heavy_four_forces_equal_gaps J provider μ
    (bj - bi) (b5 - bj) (b6 - b5) ha hb hc hheavy h01 h12 h02 h13
  dsimp
  omega

/-- The lower two values completing an AP below a fixed defining pair are
unique.  Strict ordering of the speed list then upgrades value uniqueness to
index-pair uniqueness. -/
theorem ordered_AP_completion_unique
    (q x y x' y' b5 b6 : ℕ)
    (hx : x + 2 * q = b5) (hy : y + q = b5) (h6 : b6 = b5 + q)
    (hx' : x' + 2 * q = b5) (hy' : y' + q = b5)
    (_h6' : b6 = b5 + q) :
    x = x' ∧ y = y' := by
  omega

/-! ## Direct q-masks and arbitrary clock restrictions -/

/-- Strict danger for a beat numerator `p`.  Only the residue of `speed`
modulo `q` matters. -/
def directDanger (q speed p : ℕ) : Prop :=
  14 * min ((speed * p) % q) ((q - (speed * p) % q) % q) < q

/-- The direct q-danger mask. -/
def directQMask (q speed : ℕ) : Finset (Fin q) :=
  Finset.univ.filter fun p => directDanger q speed p.val

/-- Restrict the pointwise q-danger predicate to representatives
`0,...,L-1`.  When `L=q/d0` is the mixed-period master clock and the quotient
period divides `L`, this is the corresponding lifted mask.  For arbitrary
`L` it is only a representative-range restriction, not a clock claim. -/
def numeratorRestriction (L q speed : ℕ) : Finset (Fin L) :=
  Finset.univ.filter fun p => directDanger q speed p.val

theorem directDanger_iff_of_speed_mod_eq
    {q speed speed' : ℕ} (hmod : speed % q = speed' % q) (p : ℕ) :
    directDanger q speed p ↔ directDanger q speed' p := by
  simp only [directDanger]
  have hproduct : (speed * p) % q = (speed' * p) % q := by
    simp [Nat.mul_mod, hmod]
  rw [hproduct]

theorem directQMask_eq_of_speed_mod_eq
    {q speed speed' : ℕ} (hmod : speed % q = speed' % q) :
    directQMask q speed = directQMask q speed' := by
  ext p
  simp only [directQMask, Finset.mem_filter, Finset.mem_univ, true_and]
  exact directDanger_iff_of_speed_mod_eq hmod p.val

theorem numeratorRestriction_eq_of_speed_mod_eq
    {L q speed speed' : ℕ} (hmod : speed % q = speed' % q) :
    numeratorRestriction L q speed = numeratorRestriction L q speed' := by
  ext p
  simp only [numeratorRestriction, Finset.mem_filter, Finset.mem_univ, true_and]
  exact directDanger_iff_of_speed_mod_eq hmod p.val

/-- A completed AP is one residue class modulo its defining gap. -/
theorem AP_completion_residues_coincide
    (q bi bj b5 b6 : ℕ)
    (hi : bi + 2 * q = b5) (hj : bj + q = b5) (h6 : b6 = b5 + q) :
    bi % q = bj % q ∧ bj % q = b5 % q ∧ b5 % q = b6 % q := by
  have hi5 : bi % q = b5 % q := by
    rw [← hi]
    simp
  have hj5 : bj % q = b5 % q := by
    rw [← hj]
    simp
  have h56 : b5 % q = b6 % q := by
    rw [h6]
    simp
  exact ⟨hi5.trans hj5.symm, hj5, h56⟩

/-- All four direct masks and all four representative-range restrictions
coincide.  In THM-1217, the actual master clock divides `q`; no divisibility
hypothesis is needed for this pointwise equality. -/
theorem AP_direct_and_numerator_restrictions_coincide
    (L q bi bj b5 b6 : ℕ)
    (hi : bi + 2 * q = b5) (hj : bj + q = b5) (h6 : b6 = b5 + q) :
    (directQMask q bi = directQMask q bj ∧
      directQMask q bj = directQMask q b5 ∧
      directQMask q b5 = directQMask q b6) ∧
    (numeratorRestriction L q bi = numeratorRestriction L q bj ∧
      numeratorRestriction L q bj = numeratorRestriction L q b5 ∧
      numeratorRestriction L q b5 = numeratorRestriction L q b6) := by
  obtain ⟨hij, hj5, h56⟩ := AP_completion_residues_coincide q bi bj b5 b6 hi hj h6
  exact ⟨
    ⟨directQMask_eq_of_speed_mod_eq hij,
      directQMask_eq_of_speed_mod_eq hj5,
      directQMask_eq_of_speed_mod_eq h56⟩,
    ⟨numeratorRestriction_eq_of_speed_mod_eq hij,
      numeratorRestriction_eq_of_speed_mod_eq hj5,
      numeratorRestriction_eq_of_speed_mod_eq h56⟩⟩

/-! ## Three mask classes and the improved common-Q threshold -/

/-- Natural-number ceiling of `Q/14`, repeated locally so this formalization
does not depend on the independent THM-1216 module. -/
def ceilDiv14 (Q : ℕ) : ℕ :=
  (Q + 13) / 14

/-- Cardinality of one strict danger window in a reduced period. -/
def windowCount (Q : ℕ) : ℕ :=
  2 * ceilDiv14 Q - 1

/-- Escape threshold for a union represented by `classes` common-zero masks. -/
def classThreshold (Q classes : ℕ) : ℕ :=
  2 + classes * (windowCount Q - 1)

private theorem one_le_ceilDiv14 (Q : ℕ) (hQ : 1 ≤ Q) :
    1 ≤ ceilDiv14 Q := by
  unfold ceilDiv14
  omega

/-- If three of the five representative labels are the AP mask, their union
is represented by only three masks. -/
theorem five_mask_union_collapses_to_three
    {α : Type*} [DecidableEq α]
    (s0 s1 s2 s3 s4 : Finset α) (h01 : s0 = s1) (h04 : s0 = s4) :
    s0 ∪ s1 ∪ s2 ∪ s3 ∪ s4 = s0 ∪ s2 ∪ s3 := by
  rw [← h01, ← h04]
  ext x
  simp only [Finset.mem_union]
  tauto

/-- A common point saves two units from the sum of three mask sizes. -/
theorem three_union_card_add_two_le_sum_cards
    {α : Type*} [DecidableEq α]
    (x : α) (s0 s1 s2 : Finset α)
    (hx0 : x ∈ s0) (hx1 : x ∈ s1) (hx2 : x ∈ s2) :
    (s0 ∪ s1 ∪ s2).card + 2 ≤ s0.card + s1.card + s2.card := by
  have h01 := Finset.card_union_add_card_inter s0 s1
  have h012 := Finset.card_union_add_card_inter (s0 ∪ s1) s2
  have hi01 : 1 ≤ (s0 ∩ s1).card := by
    apply Finset.card_pos.mpr
    exact ⟨x, by simp [hx0, hx1]⟩
  have hi012 : 1 ≤ ((s0 ∪ s1) ∩ s2).card := by
    apply Finset.card_pos.mpr
    exact ⟨x, by simp [hx0, hx2]⟩
  omega

/-- Three common-zero masks, each of size at most `A`, have union cardinality
at most `3A-2`. -/
theorem three_union_card_add_two_le_three_mul
    {α : Type*} [DecidableEq α]
    (x : α) (s0 s1 s2 : Finset α) (A : ℕ)
    (hx0 : x ∈ s0) (hx1 : x ∈ s1) (hx2 : x ∈ s2)
    (h0 : s0.card ≤ A) (h1 : s1.card ≤ A) (h2 : s2.card ≤ A) :
    (s0 ∪ s1 ∪ s2).card + 2 ≤ 3 * A := by
  have hsaving := three_union_card_add_two_le_sum_cards x s0 s1 s2 hx0 hx1 hx2
  omega

/-- Exact closed form of the heavy-circuit three-class threshold. -/
theorem classThreshold_three_formula (Q : ℕ) (hQ : 1 ≤ Q) :
    classThreshold Q 3 = 6 * ceilDiv14 Q - 4 := by
  unfold classThreshold windowCount
  have := one_le_ceilDiv14 Q hQ
  omega

/-- The improved three-class threshold still fits in every nontrivial
reduced period. -/
theorem classThreshold_three_le_period (Q : ℕ) (hQ : 2 ≤ Q) :
    classThreshold Q 3 ≤ Q := by
  unfold classThreshold windowCount ceilDiv14
  omega

/-- Typed finite consumer for the three-mask collapse. -/
theorem three_mask_bridge_impossible
    {α : Type*} [DecidableEq α]
    (Q : ℕ) (x : α) (block s0 s1 s2 : Finset α)
    (hx0 : x ∈ s0) (hx1 : x ∈ s1) (hx2 : x ∈ s2)
    (hblock : classThreshold Q 3 ≤ block.card)
    (h0 : s0.card ≤ windowCount Q)
    (h1 : s1.card ≤ windowCount Q)
    (h2 : s2.card ≤ windowCount Q)
    (hcover : block ⊆ s0 ∪ s1 ∪ s2) : False := by
  have hunion := three_union_card_add_two_le_three_mul x s0 s1 s2
    (windowCount Q) hx0 hx1 hx2 h0 h1 h2
  have hcard := Finset.card_le_card hcover
  unfold classThreshold at hblock
  omega

#print axioms tail_le_heavyCutoff
#print axioms coreHeavyPairs_exact
#print axioms core_heavy_forces_ratio12
#print axioms heavy_four_forces_equal_gaps
#print axioms heavy_ordered_quadruple_forces_AP
#print axioms ordered_AP_completion_unique
#print axioms AP_direct_and_numerator_restrictions_coincide
#print axioms five_mask_union_collapses_to_three
#print axioms classThreshold_three_formula
#print axioms classThreshold_three_le_period
#print axioms three_mask_bridge_impossible

end
end LRC14.HeavyCircuitAPMaskCollapse
