/-
  TournamentH7.LRCRayWitness -- THE CONSTRUCTIVE RAY WITNESS
  (klein-2026-07-11-S247, HYP-5960, THM-695): ray-coupled families made
  explicit.  For S = P u {V_mid - f : f in E_mid} u {V - e : e in E_top}
  with V_mid = (a_r*V - c)/b (the rational-ray cluster of THM-688/S246),
  THE WITNESS IS THM-693's TIME AT THE b-SCALED MODULUS:

      q = q* * b * V,     w = b * (a*V + (c1 - a*V) % q*)  =  b * w_693.

  Small speeds and the top cluster INHERIT their bands by pure scaling
  ((b*x) % (b*Q) = b*(x % Q)); the ray-mid speeds u (with b*u = a_r*V - c
  - b*f) get residue

      (u*w) % (q*bV)  =  rho'_f * V - delta*(c + b*f),
      rho'_f = (a_r*(a*V + delta) - a*(c + b*f)) % (q* * b),
      delta  = (c1 - a*V) % q*,

  strictly in-band under the per-f DECIDABLE digit hypothesis
  q*b < 14*rho'_f < 13*q*b and the threshold 14*q**(c + b*f) < V.

  With mac-mini's cont.26/27 (lrc14_from_citations_only; the ruler-embedding
  de-citation narrowing THM-527-A to CERTIFICATE SUPPLY), this file extends
  the supply's constructive coverage to ray-coupled residual shapes:
  witnesses for [two-scale (THM-693)] u [multi-scale (THM-694)] u
  [ray-coupled (THIS)] are all explicit mixed-radix numbers.

  DEMO: P = {1,2,3}, ray (1,2,0) at V_mid = 5000 with E_mid = {0,1}, top
  E_top = {0,1,2,3,5,8,13,21} at V = 10000: witness (q, w) = (260000, 20002)
  -- THM-693's demo time with the doubled modulus.

  Kernel-pure: no native_decide, no sorry.
-/
import Mathlib
import TournamentH7.LRCTwoScaleWitness

namespace LonelyRunner
namespace RayWitness

open LRC14Grand

/-- Bands scale: `(b·x) % (b·Q) = b·(x % Q)`, so a strict band at modulus `Q`
is a strict band at modulus `b·Q` for the `b`-scaled value. -/
theorem scaled_band {b Q x : ℤ} (hb : 0 < b) (hQ : 0 < Q)
    (h : Q < 14 * (x % Q) ∧ 14 * (x % Q) < 13 * Q) :
    b * Q < 14 * ((b * x) % (b * Q)) ∧
    14 * ((b * x) % (b * Q)) < 13 * (b * Q) := by
  have hr0 : 0 ≤ x % Q := Int.emod_nonneg _ (by omega)
  have hrQ : x % Q < Q := Int.emod_lt_of_pos _ hQ
  have hdm : Q * (x / Q) + x % Q = x := Int.ediv_add_emod x Q
  have hkey : b * x = b * (x % Q) + (b * Q) * (x / Q) := by
    linear_combination b * hdm.symm
  have hb0 : 0 ≤ b * (x % Q) := by positivity
  have hbq : b * (x % Q) < b * Q := by nlinarith
  have hmod : (b * x) % (b * Q) = b * (x % Q) := by
    rw [hkey, Int.add_mul_emod_self_left, Int.emod_eq_of_lt hb0 hbq]
  rw [hmod]
  constructor
  · nlinarith [h.1]
  · nlinarith [h.2]

/-- **The ray-mid band**: a ray-cluster speed `u` (with `b·u = a_r·V − c −
b·f`) is strictly in-band at `(q*·b·V, b·w₆₉₃)` under the decidable digit
hypothesis on `ρ'_f` and the threshold `14·q*·(c+bf) < V`. -/
theorem ray_speed_band {qs a c1 V ar b c f u : ℤ}
    (hqs0 : 0 < qs) (hb : 0 < b)
    (hu : b * u = ar * V - c - b * f)
    (hrho : qs * b <
        14 * ((ar * (a * V + (c1 - a * V) % qs) - a * (c + b * f)) % (qs * b)) ∧
        14 * ((ar * (a * V + (c1 - a * V) % qs) - a * (c + b * f)) % (qs * b))
          < 13 * (qs * b))
    (hcf0 : 0 ≤ c + b * f)
    (hV : 14 * qs * (c + b * f) < V) (hV0 : 0 < V) :
    qs * b * V <
      14 * ((u * (b * (a * V + (c1 - a * V) % qs))) % (qs * b * V)) ∧
    14 * ((u * (b * (a * V + (c1 - a * V) % qs))) % (qs * b * V))
      < 13 * (qs * b * V) := by
  obtain ⟨hρ1, hρ2⟩ := hrho
  set δ : ℤ := (c1 - a * V) % qs with hδdef
  have hδ0 : 0 ≤ δ := Int.emod_nonneg _ (by omega)
  have hδq : δ < qs := Int.emod_lt_of_pos _ hqs0
  have hqb0 : (0 : ℤ) < qs * b := by positivity
  set ρ : ℤ := (ar * (a * V + δ) - a * (c + b * f)) % (qs * b) with hρdef
  set m : ℤ := (ar * (a * V + δ) - a * (c + b * f)) / (qs * b) with hmdef
  have hdm : (qs * b) * m + ρ = ar * (a * V + δ) - a * (c + b * f) :=
    Int.ediv_add_emod _ _
  have hρ0 : 0 ≤ ρ := Int.emod_nonneg _ (by omega)
  have hρq : ρ < qs * b := Int.emod_lt_of_pos _ hqb0
  have hkey : u * (b * (a * V + δ))
      = (ρ * V - δ * (c + b * f)) + (qs * b * V) * m := by
    linear_combination (a * V + δ) * hu - V * hdm
  have hρge : 1 ≤ ρ := by nlinarith
  have hb0 : 0 ≤ ρ * V - δ * (c + b * f) := by nlinarith
  have hbq : ρ * V - δ * (c + b * f) < qs * b * V := by nlinarith
  have hmod : (u * (b * (a * V + δ))) % (qs * b * V)
      = ρ * V - δ * (c + b * f) := by
    rw [hkey, Int.add_mul_emod_self_left, Int.emod_eq_of_lt hb0 hbq]
  rw [hmod]
  constructor
  · nlinarith
  · nlinarith

/-- **THE RAY WITNESS (THM-695)**: the ray-coupled family is strictly live
at `(q*·b·V, b·w₆₉₃)` — THM-693's witness time, scaled modulus.  Small and
top-cluster speeds inherit by scaling; ray-mid speeds by the digit
hypothesis. -/
theorem rayFamily_strictlyLive (v : Fin 13 → ℤ) (qs a c1 V ar b c : ℤ)
    (hqs8 : 8 ≤ qs) (hqs13 : qs ≤ 13) (ha1 : 1 ≤ a) (haq : a < qs)
    (hb : 0 < b) (hV : 2184 < V)
    (hv : ∀ i,
      (1 ≤ v i ∧ v i ≤ 13 ∧ (v i * a) % qs ≠ 0) ∨
      (∃ e, 0 ≤ e ∧ v i = V - e ∧ (c1 - e * a) % qs ≠ 0 ∧ 168 * e < V) ∨
      (∃ f, 0 ≤ c + b * f ∧ b * v i = ar * V - c - b * f ∧
        (qs * b < 14 * ((ar * (a * V + (c1 - a * V) % qs) - a * (c + b * f))
            % (qs * b)) ∧
         14 * ((ar * (a * V + (c1 - a * V) % qs) - a * (c + b * f))
            % (qs * b)) < 13 * (qs * b)) ∧
        14 * qs * (c + b * f) < V)) :
    StrictlyLive v (qs * b * V) (b * (a * V + (c1 - a * V) % qs)) := by
  have hqs0 : (0 : ℤ) < qs := by omega
  have hV0 : (0 : ℤ) < V := by omega
  have hδ0 : 0 ≤ (c1 - a * V) % qs := Int.emod_nonneg _ (by omega)
  have hδq : (c1 - a * V) % qs < qs := Int.emod_lt_of_pos _ hqs0
  refine ⟨?_, ?_, ?_⟩
  · have hw0 : 0 < a * V + (c1 - a * V) % qs := by nlinarith
    nlinarith [hw0, hb]
  · -- w < q:  b·w₆₉₃ < b·(q*·V) ⟸ w₆₉₃ < q*·V
    have hw : a * V + (c1 - a * V) % qs < qs * V := by nlinarith
    nlinarith [hw, hb]
  · intro i
    rcases hv i with ⟨h1, h13, hr⟩ | ⟨e, he0, hve, hs, hVe⟩ |
      ⟨f, hcf0, hu, hrho, hVf⟩
    · -- small: THM-693's band at (q*V, w₆₉₃), then scale by b
      have h693 := TwoScaleWitness.small_speed_band (c := c1)
        hqs8 hqs13 h1 h13 hr hV
      have hsc := scaled_band hb (by positivity) h693
      have e1 : b * (v i * (a * V + (c1 - a * V) % qs))
          = v i * (b * (a * V + (c1 - a * V) % qs)) := by ring
      have e2 : b * (qs * V) = qs * b * V := by ring
      rw [e1, e2] at hsc
      exact hsc
    · -- top cluster: THM-693's band, then scale
      have h693 := TwoScaleWitness.cluster_speed_band (c := c1)
        hqs8 hqs13 he0 hs hVe hV
      have hsc := scaled_band hb (by positivity) h693
      have e1 : b * ((V - e) * (a * V + (c1 - a * V) % qs))
          = (V - e) * (b * (a * V + (c1 - a * V) % qs)) := by ring
      have e2 : b * (qs * V) = qs * b * V := by ring
      rw [e1, e2] at hsc
      rw [hve]
      exact hsc
    · -- ray-mid speed
      exact ray_speed_band hqs0 hb hu hrho hcf0 hVf hV0

/-- The kps composition: the ray family has a strict witness at the explicit
rational time `b·w₆₉₃/(q*·b·V) = w₆₉₃/(q*·V)`. -/
theorem rayFamily_strictWitness (v : Fin 13 → ℤ) (qs a c1 V ar b c : ℤ)
    (hqs8 : 8 ≤ qs) (hqs13 : qs ≤ 13) (ha1 : 1 ≤ a) (haq : a < qs)
    (hb : 0 < b) (hV : 2184 < V)
    (hv : ∀ i,
      (1 ≤ v i ∧ v i ≤ 13 ∧ (v i * a) % qs ≠ 0) ∨
      (∃ e, 0 ≤ e ∧ v i = V - e ∧ (c1 - e * a) % qs ≠ 0 ∧ 168 * e < V) ∨
      (∃ f, 0 ≤ c + b * f ∧ b * v i = ar * V - c - b * f ∧
        (qs * b < 14 * ((ar * (a * V + (c1 - a * V) % qs) - a * (c + b * f))
            % (qs * b)) ∧
         14 * ((ar * (a * V + (c1 - a * V) % qs) - a * (c + b * f))
            % (qs * b)) < 13 * (qs * b)) ∧
        14 * qs * (c + b * f) < V)) :
    StrictWitness v :=
  strictWitness_of_strictlyLive
    (rayFamily_strictlyLive v qs a c1 V ar b c hqs8 hqs13 ha1 haq hb hV hv)

/-! ## Demo: the ray family, witnessed through the theorem -/

def demoRay : Fin 13 → ℤ :=
  ![1, 2, 3, 4999, 5000, 9979, 9987, 9992, 9995, 9997, 9998, 9999, 10000]

theorem demoRay_strictWitness : StrictWitness demoRay := by
  refine rayFamily_strictWitness demoRay 13 1 4 10000 1 2 0
    (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by norm_num) ?_
  intro i
  fin_cases i
  · exact Or.inl ⟨by decide, by decide, by decide⟩
  · exact Or.inl ⟨by decide, by decide, by decide⟩
  · exact Or.inl ⟨by decide, by decide, by decide⟩
  · exact Or.inr (Or.inr ⟨1, by decide, by decide, ⟨by decide, by decide⟩,
      by decide⟩)
  · exact Or.inr (Or.inr ⟨0, by decide, by decide, ⟨by decide, by decide⟩,
      by decide⟩)
  · exact Or.inr (Or.inl ⟨21, by decide, by decide, by decide, by decide⟩)
  · exact Or.inr (Or.inl ⟨13, by decide, by decide, by decide, by decide⟩)
  · exact Or.inr (Or.inl ⟨8, by decide, by decide, by decide, by decide⟩)
  · exact Or.inr (Or.inl ⟨5, by decide, by decide, by decide, by decide⟩)
  · exact Or.inr (Or.inl ⟨3, by decide, by decide, by decide, by decide⟩)
  · exact Or.inr (Or.inl ⟨2, by decide, by decide, by decide, by decide⟩)
  · exact Or.inr (Or.inl ⟨1, by decide, by decide, by decide, by decide⟩)
  · exact Or.inr (Or.inl ⟨0, by decide, by decide, by decide, by decide⟩)

end RayWitness
end LonelyRunner

-- kernel-purity audit (fleet convention)
#print axioms LonelyRunner.RayWitness.scaled_band
#print axioms LonelyRunner.RayWitness.ray_speed_band
#print axioms LonelyRunner.RayWitness.rayFamily_strictlyLive
#print axioms LonelyRunner.RayWitness.demoRay_strictWitness
