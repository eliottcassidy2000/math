/-
  TournamentH7.LRCMultiScaleWitness -- THE MULTI-SCALE WITNESS
  (klein-2026-07-10-S245, HYP-5950, THM-694): THM-693 iterated per scale.

  Two generic lemmas drive the recursion:

  (I) THE BAND LIFT: a speed strictly in-band at (q, a) with 0 < p <= B stays
      strictly in-band at (q*V, a*V + Delta) for ANY Delta in [0, q), once
      14*B*q < V.  Integrality supplies the margins for free (q < 14r gives
      q+1 <= 14r; 14r < 13q gives 14r <= 13q-1); the lift never inspects
      Delta -- so the next level's tuning cannot disturb lower levels.

  (II) THE CLUSTER JOIN: if the missed-class residue s = (C - e*a) % q is
      itself strictly in-band [q < 14s < 13q] and 14*e*q < V, then the new
      cluster speed V - e is strictly in-band at (q*V, a*V + (C - a*V) % q).

  (III) miss_at_next: the COARSE missed class c1 mod q* (always available --
      every cluster has <= 12 co-offsets < q*) lifts to a band-interior
      missed residue at the composite modulus:
      (c1*V2 - e*(a*V2 + D2)) % (q**V2) = sigma_e*V2 - e*D2.

  THE r = 2 ASSEMBLY (threeScale_strictlyLive): P u C2 u C1 is strictly live
  at modulus q**V2*V1 with the explicit mixed-radix multiplier
      w1 = w2*V1 + (c1*V2 - w2*V1) % (q**V2),   w2 = a*V2 + (c2 - a*V2) % q*,
  under V2 > 2184, 168f < V2 per inner co-offset, 14*e*(q**V2) < V1 and
  168e < V2 per outer co-offset, and 14*(q**V2)*V2 < V1 (quadratic scale
  separation -- the constructive price for the continuum's V1/V2 -> inf).
  General r: apply (I) to every already-placed speed and (II)+(III) to each
  new cluster, innermost out -- the same two lemmas, once per scale.

  DEMO: P = {1,2,3}, E2 = {0,1} @ V2 = 3000, E1 = {0..7} @ V1 = 2*10^9,
  (q*, a, c2, c1) = (13, 1, 2, 8): witness (q1, w1) =
  (78000000000000, 6010000020000), all thirteen residues strictly in-band.

  Kernel-pure: no native_decide, no sorry.
-/
import Mathlib
import TournamentH7.LRCTwoScaleWitness

namespace LonelyRunner
namespace MultiScaleWitness

open LRC14Grand

/-- **THE BAND LIFT**: strict in-band at `(q, a)` is preserved to
`(q·V, a·V + Δ)` for ANY `Δ ∈ [0, q)`, for speeds `0 < p ≤ B`, once
`14·B·q < V`.  The lift never inspects `Δ`. -/
theorem band_lift {q a V Δ p B : ℤ}
    (hband : q < 14 * ((p * a) % q) ∧ 14 * ((p * a) % q) < 13 * q)
    (hp : 0 < p) (hpB : p ≤ B) (hq : 0 < q)
    (hΔ0 : 0 ≤ Δ) (hΔq : Δ < q)
    (hV : 14 * B * q < V) :
    q * V < 14 * ((p * (a * V + Δ)) % (q * V)) ∧
    14 * ((p * (a * V + Δ)) % (q * V)) < 13 * (q * V) := by
  obtain ⟨hb1, hb2⟩ := hband
  set r : ℤ := (p * a) % q with hrdef
  set j : ℤ := (p * a) / q with hjdef
  have hdm : q * j + r = p * a := Int.ediv_add_emod (p * a) q
  have hr0 : 0 ≤ r := Int.emod_nonneg _ (by omega)
  have hrq : r < q := Int.emod_lt_of_pos _ hq
  have hV0 : 0 < V := by nlinarith
  have hkey : p * (a * V + Δ) = (r * V + p * Δ) + (q * V) * j := by
    linear_combination V * hdm.symm
  have hb0 : 0 ≤ r * V + p * Δ := by nlinarith
  have hbq : r * V + p * Δ < q * V := by nlinarith
  have hmod : (p * (a * V + Δ)) % (q * V) = r * V + p * Δ := by
    rw [hkey, Int.add_mul_emod_self_left, Int.emod_eq_of_lt hb0 hbq]
  rw [hmod]
  constructor
  · nlinarith
  · nlinarith

/-- **THE CLUSTER JOIN**: if the missed-class residue `(C − e·a) % q` is
strictly in-band, the new cluster speed `V − e` is strictly in-band at
`(q·V, a·V + (C − a·V) % q)`, once `14·e·q < V`. -/
theorem cluster_join {q a V C e : ℤ}
    (hmiss : q < 14 * ((C - e * a) % q) ∧ 14 * ((C - e * a) % q) < 13 * q)
    (he : 0 ≤ e) (hq : 0 < q) (hV0 : 0 < V)
    (hV : 14 * e * q < V) :
    q * V < 14 * (((V - e) * (a * V + (C - a * V) % q)) % (q * V)) ∧
    14 * (((V - e) * (a * V + (C - a * V) % q)) % (q * V)) < 13 * (q * V) := by
  obtain ⟨hm1, hm2⟩ := hmiss
  set Δ : ℤ := (C - a * V) % q with hΔdef
  have hΔ0 : 0 ≤ Δ := Int.emod_nonneg _ (by omega)
  have hΔq : Δ < q := Int.emod_lt_of_pos _ hq
  set w : ℤ := a * V + Δ with hwdef
  have hwc : w % q = C % q := by
    rw [hwdef, hΔdef, Int.add_emod, Int.emod_emod_of_dvd _ dvd_rfl,
        ← Int.add_emod]
    ring_nf
  set s : ℤ := (w - e * a) % q with hsdef
  set m : ℤ := (w - e * a) / q with hmdef
  have hdm : q * m + s = w - e * a := Int.ediv_add_emod (w - e * a) q
  have hsC : s = (C - e * a) % q := by
    rw [hsdef, Int.sub_emod, hwc, ← Int.sub_emod]
  have hs0 : 0 ≤ s := Int.emod_nonneg _ (by omega)
  have hsq : s < q := Int.emod_lt_of_pos _ hq
  have hs1 : q < 14 * s := by rw [hsC]; exact hm1
  have hs2 : 14 * s < 13 * q := by rw [hsC]; exact hm2
  have hkey : (V - e) * w = (s * V - e * Δ) + (q * V) * m := by
    linear_combination V * hdm.symm
  have hb0 : 0 ≤ s * V - e * Δ := by nlinarith
  have hbq : s * V - e * Δ < q * V := by nlinarith
  have hmod : ((V - e) * w) % (q * V) = s * V - e * Δ := by
    rw [hkey, Int.add_mul_emod_self_left, Int.emod_eq_of_lt hb0 hbq]
  rw [show ((V - e) * (a * V + (C - a * V) % q)) = (V - e) * w by rw [hwdef]]
  rw [hmod]
  constructor
  · nlinarith
  · nlinarith

/-- **The coarse missed class lifts**: `(c₁ − e·a) % q* ≠ 0` at the base
gives a band-interior missed residue at the composite modulus `q*·V₂`:
`(c₁·V₂ − e·(a·V₂ + Δ₂)) % (q*·V₂) = σ_e·V₂ − e·Δ₂`, strictly in-band. -/
theorem miss_at_next {qs a c1 V2 e Δ2 : ℤ}
    (hqs8 : 8 ≤ qs) (hqs13 : qs ≤ 13)
    (hσ : (c1 - e * a) % qs ≠ 0) (he : 0 ≤ e)
    (hΔ0 : 0 ≤ Δ2) (hΔq : Δ2 < qs)
    (hVe : 168 * e < V2) (hV2 : 2184 < V2) :
    qs * V2 < 14 * ((c1 * V2 - e * (a * V2 + Δ2)) % (qs * V2)) ∧
    14 * ((c1 * V2 - e * (a * V2 + Δ2)) % (qs * V2)) < 13 * (qs * V2) := by
  have hqs0 : (0 : ℤ) < qs := by omega
  set σ : ℤ := (c1 - e * a) % qs with hσdef
  set k : ℤ := (c1 - e * a) / qs with hkdef
  have hdm : qs * k + σ = c1 - e * a := Int.ediv_add_emod (c1 - e * a) qs
  have hσ0 : 0 ≤ σ := Int.emod_nonneg _ (by omega)
  have hσq : σ < qs := Int.emod_lt_of_pos _ hqs0
  have hσ1 : 1 ≤ σ := by
    rcases lt_or_eq_of_le hσ0 with h | h
    · omega
    · exact absurd h.symm hσ
  have hkey : c1 * V2 - e * (a * V2 + Δ2)
      = (σ * V2 - e * Δ2) + (qs * V2) * k := by
    linear_combination V2 * hdm.symm
  have hb0 : 0 ≤ σ * V2 - e * Δ2 := by nlinarith
  have hbq : σ * V2 - e * Δ2 < qs * V2 := by nlinarith
  have hmod : (c1 * V2 - e * (a * V2 + Δ2)) % (qs * V2) = σ * V2 - e * Δ2 := by
    rw [hkey, Int.add_mul_emod_self_left, Int.emod_eq_of_lt hb0 hbq]
  rw [hmod]
  constructor
  · nlinarith
  · nlinarith

/-- **THE THREE-SCALE WITNESS (THM-694, r = 2)**: `P ∪ C₂ ∪ C₁` is strictly
live at modulus `q*·V₂·V₁` with the explicit mixed-radix multiplier
`w₁ = w₂·V₁ + (c₁·V₂ − w₂·V₁) % (q*·V₂)`, `w₂ = a·V₂ + (c₂ − a·V₂) % q*`. -/
theorem threeScale_strictlyLive (v : Fin 13 → ℤ) (qs a c2 c1 V2 V1 : ℤ)
    (hqs8 : 8 ≤ qs) (hqs13 : qs ≤ 13) (ha1 : 1 ≤ a) (haq : a < qs)
    (hV2 : 2184 < V2)
    (hV1 : 14 * (qs * V2) * V2 < V1)
    (hv : ∀ i, (1 ≤ v i ∧ v i ≤ 13 ∧ (v i * a) % qs ≠ 0) ∨
          (∃ f, 0 ≤ f ∧ v i = V2 - f ∧ (c2 - f * a) % qs ≠ 0 ∧ 168 * f < V2) ∨
          (∃ e, 0 ≤ e ∧ v i = V1 - e ∧ (c1 - e * a) % qs ≠ 0 ∧ 168 * e < V2 ∧
            14 * e * (qs * V2) < V1)) :
    StrictlyLive v (qs * V2 * V1)
      ((a * V2 + (c2 - a * V2) % qs) * V1 +
        (c1 * V2 - (a * V2 + (c2 - a * V2) % qs) * V1) % (qs * V2)) := by
  have hqs0 : (0 : ℤ) < qs := by omega
  have hV20 : (0 : ℤ) < V2 := by omega
  have hq20 : (0 : ℤ) < qs * V2 := by positivity
  have hV10 : (0 : ℤ) < V1 := by nlinarith
  have hΔ20 : 0 ≤ (c2 - a * V2) % qs := Int.emod_nonneg _ (by omega)
  have hΔ2q : (c2 - a * V2) % qs < qs := Int.emod_lt_of_pos _ hqs0
  have hΔ10 : 0 ≤ (c1 * V2 - (a * V2 + (c2 - a * V2) % qs) * V1) % (qs * V2) :=
    Int.emod_nonneg _ (by positivity)
  have hΔ1q : (c1 * V2 - (a * V2 + (c2 - a * V2) % qs) * V1) % (qs * V2)
      < qs * V2 := Int.emod_lt_of_pos _ hq20
  have hw2pos : 0 < a * V2 + (c2 - a * V2) % qs := by nlinarith
  have hw2lt : a * V2 + (c2 - a * V2) % qs < qs * V2 := by nlinarith
  refine ⟨by nlinarith, ?_, ?_⟩
  · -- w₁ < q*·V₂·V₁
    have h1 : (a * V2 + (c2 - a * V2) % qs) * V1 ≤ (qs * V2 - 1) * V1 := by
      nlinarith
    nlinarith
  · intro i
    have hassoc : qs * V2 * V1 = (qs * V2) * V1 := by ring
    rw [hassoc]
    rcases hv i with ⟨h1, h13, hr⟩ | ⟨f, hf0, hvf, hc2f, hVf⟩ |
      ⟨e, he0, hve, hc1e, hVe2, hVe1⟩
    · -- small speed: level-1 band (S244) + the lift
      have hL1 := TwoScaleWitness.small_speed_band (c := c2)
        hqs8 hqs13 h1 h13 hr hV2
      exact band_lift hL1 (by omega) (by omega : v i ≤ V2) hq20 hΔ10 hΔ1q
        (by nlinarith)
    · -- inner cluster: level-1 band (S244) + the lift
      have hL1 := TwoScaleWitness.cluster_speed_band (c := c2)
        hqs8 hqs13 hf0 hc2f hVf hV2
      rw [hvf]
      exact band_lift hL1 (by omega) (by omega : V2 - f ≤ V2) hq20 hΔ10 hΔ1q
        (by nlinarith)
    · -- outer cluster: the coarse miss lifts, then the join
      have hmiss := miss_at_next (a := a) (Δ2 := (c2 - a * V2) % qs)
        hqs8 hqs13 hc1e he0 hΔ20 hΔ2q hVe2 hV2
      rw [hve]
      exact cluster_join hmiss he0 hq20 hV10 hVe1

/-- The kps composition: the three-scale family has a STRICT WITNESS at the
explicit rational time `w₁/(q*·V₂·V₁)`. -/
theorem threeScale_strictWitness (v : Fin 13 → ℤ) (qs a c2 c1 V2 V1 : ℤ)
    (hqs8 : 8 ≤ qs) (hqs13 : qs ≤ 13) (ha1 : 1 ≤ a) (haq : a < qs)
    (hV2 : 2184 < V2) (hV1 : 14 * (qs * V2) * V2 < V1)
    (hv : ∀ i, (1 ≤ v i ∧ v i ≤ 13 ∧ (v i * a) % qs ≠ 0) ∨
          (∃ f, 0 ≤ f ∧ v i = V2 - f ∧ (c2 - f * a) % qs ≠ 0 ∧ 168 * f < V2) ∨
          (∃ e, 0 ≤ e ∧ v i = V1 - e ∧ (c1 - e * a) % qs ≠ 0 ∧ 168 * e < V2 ∧
            14 * e * (qs * V2) < V1)) :
    StrictWitness v :=
  strictWitness_of_strictlyLive
    (threeScale_strictlyLive v qs a c2 c1 V2 V1 hqs8 hqs13 ha1 haq hV2 hV1 hv)

/-! ## Demo: three scales, witnessed through the theorem -/

def demoThreeScale : Fin 13 → ℤ :=
  ![1, 2, 3, 2999, 3000, 1999999993, 1999999994, 1999999995, 1999999996,
    1999999997, 1999999998, 1999999999, 2000000000]

theorem demoThreeScale_strictWitness : StrictWitness demoThreeScale := by
  refine threeScale_strictWitness demoThreeScale 13 1 2 8 3000 2000000000
    (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by norm_num) ?_
  intro i
  fin_cases i
  · exact Or.inl ⟨by decide, by decide, by decide⟩
  · exact Or.inl ⟨by decide, by decide, by decide⟩
  · exact Or.inl ⟨by decide, by decide, by decide⟩
  · exact Or.inr (Or.inl ⟨1, by decide, by decide, by decide, by decide⟩)
  · exact Or.inr (Or.inl ⟨0, by decide, by decide, by decide, by decide⟩)
  · exact Or.inr (Or.inr ⟨7, by decide, by decide, by decide, by decide,
      by decide⟩)
  · exact Or.inr (Or.inr ⟨6, by decide, by decide, by decide, by decide,
      by decide⟩)
  · exact Or.inr (Or.inr ⟨5, by decide, by decide, by decide, by decide,
      by decide⟩)
  · exact Or.inr (Or.inr ⟨4, by decide, by decide, by decide, by decide,
      by decide⟩)
  · exact Or.inr (Or.inr ⟨3, by decide, by decide, by decide, by decide,
      by decide⟩)
  · exact Or.inr (Or.inr ⟨2, by decide, by decide, by decide, by decide,
      by decide⟩)
  · exact Or.inr (Or.inr ⟨1, by decide, by decide, by decide, by decide,
      by decide⟩)
  · exact Or.inr (Or.inr ⟨0, by decide, by decide, by decide, by decide,
      by decide⟩)

end MultiScaleWitness
end LonelyRunner

-- kernel-purity audit (fleet convention)
#print axioms LonelyRunner.MultiScaleWitness.band_lift
#print axioms LonelyRunner.MultiScaleWitness.cluster_join
#print axioms LonelyRunner.MultiScaleWitness.miss_at_next
#print axioms LonelyRunner.MultiScaleWitness.threeScale_strictlyLive
#print axioms LonelyRunner.MultiScaleWitness.demoThreeScale_strictWitness
