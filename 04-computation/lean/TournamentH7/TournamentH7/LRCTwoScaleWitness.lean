/-
  TournamentH7.LRCTwoScaleWitness -- THE FINITE-V TWO-SCALE WITNESS
  (klein-2026-07-10-S244, HYP-5945, THM-693): the test-point program made
  CONSTRUCTIVE at finite V.  For the two-scale family
  S_V = P u {V - e : e in E}, the EXPLICIT multiplier

      w  =  a*V + (c - a*V) % qs        at modulus  q = qs*V

  is STRICTLY LIVE, provided: qs in [8,13]; every small speed p in [1,13]
  has (p*a) % qs != 0 (the q*-test-point P-side); every co-offset e has
  (c - e*a) % qs != 0 (c = a missed class, supplied by the pigeonhole);
  and V > 2184, V > 168*e per co-offset.  The two residue computations:

      (p*w)     % (qs*V)  =  r_p*V + p*D    (r_p = (p*a) % qs  in [1, qs-1])
      ((V-e)*w) % (qs*V)  =  s_e*V - e*D    (s_e = (c - e*a) % qs, likewise)

  with D = (c - a*V) % qs in [0, qs).  Both land strictly inside
  (qs*V/14, 13*qs*V/14): the small speeds because 14*r_p >= qs + 1, the
  cluster because 14*s_e*V dominates 14*e*D once V > 168*e.

  NO LIMITS, NO TRANSFER, NO MEASURE THEORY: the witness time w/(qs*V) is
  direct, and kps's strictWitness_of_strictlyLive carries it to a strict
  witness.  This is the finite-V shadow of THM-690/691/692 -- the missed
  class c is exactly what pigeonhole_missed_class (S243) produces.

  DEMO: P = {1..5}, E = {0,1,2,3,5,8,13,21}, qs = 13, a = 1, c = 4,
  V = 10000: the witness (q, w) = (130000, 10001), all thirteen residues
  strictly in-band, via the theorem with decide-checked hypotheses.

  Kernel-pure: no native_decide, no sorry.
-/
import Mathlib
import TournamentH7.LRCStrictRuler

namespace LonelyRunner
namespace TwoScaleWitness

open LRC14Grand

/-- **The small-speed band (P-side)**: at the witness multiplier, a speed
`p ∈ [1,13]` with nonzero test-point residue lands strictly in-band at
modulus `qs·V`, for every `V > 2184`. -/
theorem small_speed_band {qs a c V p : ℤ}
    (hqs8 : 8 ≤ qs) (hqs13 : qs ≤ 13)
    (hp1 : 1 ≤ p) (hp13 : p ≤ 13)
    (hr : (p * a) % qs ≠ 0) (hV : 2184 < V) :
    qs * V < 14 * ((p * (a * V + (c - a * V) % qs)) % (qs * V)) ∧
    14 * ((p * (a * V + (c - a * V) % qs)) % (qs * V)) < 13 * (qs * V) := by
  have hqs0 : (0 : ℤ) < qs := by omega
  set Δ : ℤ := (c - a * V) % qs with hΔdef
  have hΔ0 : 0 ≤ Δ := Int.emod_nonneg _ (by omega)
  have hΔq : Δ < qs := Int.emod_lt_of_pos _ hqs0
  set r : ℤ := (p * a) % qs with hrdef
  set j : ℤ := (p * a) / qs with hjdef
  have hdm : qs * j + r = p * a := Int.ediv_add_emod (p * a) qs
  have hr0 : 0 ≤ r := Int.emod_nonneg _ (by omega)
  have hrq : r < qs := Int.emod_lt_of_pos _ hqs0
  have hr1 : 1 ≤ r := by
    rcases lt_or_eq_of_le hr0 with h | h
    · omega
    · exact absurd h.symm hr
  have hkey : p * (a * V + Δ) = (r * V + p * Δ) + (qs * V) * j := by
    linear_combination V * hdm.symm
  have hb0 : 0 ≤ r * V + p * Δ := by nlinarith
  have hbq : r * V + p * Δ < qs * V := by nlinarith
  have hmod : (p * (a * V + Δ)) % (qs * V) = r * V + p * Δ := by
    rw [hkey, Int.add_mul_emod_self_left, Int.emod_eq_of_lt hb0 hbq]
  rw [hmod]
  constructor
  · nlinarith
  · nlinarith

/-- **The cluster band (E-side)**: at the witness multiplier, a cluster speed
`V − e` with nonzero missed-class residue lands strictly in-band at modulus
`qs·V`, for every `V > 168·e`. -/
theorem cluster_speed_band {qs a c V e : ℤ}
    (hqs8 : 8 ≤ qs) (hqs13 : qs ≤ 13)
    (he0 : 0 ≤ e)
    (hs : (c - e * a) % qs ≠ 0)
    (hVe : 168 * e < V) (hV : 2184 < V) :
    qs * V < 14 * (((V - e) * (a * V + (c - a * V) % qs)) % (qs * V)) ∧
    14 * (((V - e) * (a * V + (c - a * V) % qs)) % (qs * V)) < 13 * (qs * V) := by
  have hqs0 : (0 : ℤ) < qs := by omega
  set Δ : ℤ := (c - a * V) % qs with hΔdef
  have hΔ0 : 0 ≤ Δ := Int.emod_nonneg _ (by omega)
  have hΔq : Δ < qs := Int.emod_lt_of_pos _ hqs0
  set w : ℤ := a * V + Δ with hwdef
  -- w ≡ c (mod qs), so the missed-class residue transfers to w − e·a
  have hwc : w % qs = c % qs := by
    rw [hwdef, hΔdef, Int.add_emod, Int.emod_emod_of_dvd _ dvd_rfl,
        ← Int.add_emod]
    ring_nf
  set s : ℤ := (w - e * a) % qs with hsdef
  set m : ℤ := (w - e * a) / qs with hmdef
  have hdm : qs * m + s = w - e * a := Int.ediv_add_emod (w - e * a) qs
  have hsc : s = (c - e * a) % qs := by
    rw [hsdef, Int.sub_emod, hwc, ← Int.sub_emod]
  have hs0 : 0 ≤ s := Int.emod_nonneg _ (by omega)
  have hsq : s < qs := Int.emod_lt_of_pos _ hqs0
  have hs1 : 1 ≤ s := by
    rcases lt_or_eq_of_le hs0 with h | h
    · omega
    · rw [hsc] at h
      exact absurd h.symm hs
  have hkey : (V - e) * w = (s * V - e * Δ) + (qs * V) * m := by
    linear_combination V * hdm.symm
  have hb0 : 0 ≤ s * V - e * Δ := by nlinarith
  have hbq : s * V - e * Δ < qs * V := by nlinarith
  have hmod : ((V - e) * w) % (qs * V) = s * V - e * Δ := by
    rw [hkey, Int.add_mul_emod_self_left, Int.emod_eq_of_lt hb0 hbq]
  rw [show ((V - e) * (a * V + (c - a * V) % qs)) = (V - e) * w by rw [hwdef]]
  rw [hmod]
  constructor
  · nlinarith
  · nlinarith

/-- **THE FINITE-V TWO-SCALE WITNESS (THM-693)**: for the two-scale family
(every speed either small with nonzero test-point residue, or a cluster
speed `V − e` with nonzero missed-class residue and `V > 168·e`), the
EXPLICIT multiplier `w = a·V + (c − a·V) % qs` is strictly live at modulus
`qs·V`.  No limits, no transfer — the witness is direct. -/
theorem twoScale_strictlyLive (v : Fin 13 → ℤ) (qs a c V : ℤ)
    (hqs8 : 8 ≤ qs) (hqs13 : qs ≤ 13) (ha1 : 1 ≤ a) (haq : a < qs)
    (hV : 2184 < V)
    (hv : ∀ i, (1 ≤ v i ∧ v i ≤ 13 ∧ (v i * a) % qs ≠ 0) ∨
          (∃ e, 0 ≤ e ∧ v i = V - e ∧ (c - e * a) % qs ≠ 0 ∧ 168 * e < V)) :
    StrictlyLive v (qs * V) (a * V + (c - a * V) % qs) := by
  have hqs0 : (0 : ℤ) < qs := by omega
  have hΔ0 : 0 ≤ (c - a * V) % qs := Int.emod_nonneg _ (by omega)
  have hΔq : (c - a * V) % qs < qs := Int.emod_lt_of_pos _ hqs0
  refine ⟨by nlinarith, by nlinarith, ?_⟩
  intro i
  rcases hv i with ⟨h1, h13, hr⟩ | ⟨e, he0, hve, hs, hVe⟩
  · exact small_speed_band hqs8 hqs13 h1 h13 hr hV
  · rw [hve]
    exact cluster_speed_band hqs8 hqs13 he0 hs hVe hV

/-- The composition into kps's chain: the two-scale family has a STRICT
WITNESS — the explicit rational time `w/(qs·V)`. -/
theorem twoScale_strictWitness (v : Fin 13 → ℤ) (qs a c V : ℤ)
    (hqs8 : 8 ≤ qs) (hqs13 : qs ≤ 13) (ha1 : 1 ≤ a) (haq : a < qs)
    (hV : 2184 < V)
    (hv : ∀ i, (1 ≤ v i ∧ v i ≤ 13 ∧ (v i * a) % qs ≠ 0) ∨
          (∃ e, 0 ≤ e ∧ v i = V - e ∧ (c - e * a) % qs ≠ 0 ∧ 168 * e < V)) :
    StrictWitness v :=
  strictWitness_of_strictlyLive
    (twoScale_strictlyLive v qs a c V hqs8 hqs13 ha1 haq hV hv)

/-! ## Demo: an explicit two-scale family, witnessed by the theorem

`P = {1,…,5}`, co-offsets `E = {0,1,2,3,5,8,13,21}` at `V = 10000`,
test data `qs = 13, a = 1, c = 4` (class 4 is missed: the residues
`(4 − e) mod 13` are all nonzero).  The witness: `(q, w) = (130000, 10001)`. -/

def demoTwoScale : Fin 13 → ℤ :=
  ![1, 2, 3, 4, 5, 9979, 9987, 9992, 9995, 9997, 9998, 9999, 10000]

theorem demoTwoScale_strictlyLive :
    StrictlyLive demoTwoScale 130000 10001 := by
  have h := twoScale_strictlyLive demoTwoScale 13 1 4 10000
    (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by
      intro i
      fin_cases i
      · exact Or.inl (by refine ⟨by decide, by decide, by decide⟩)
      · exact Or.inl (by refine ⟨by decide, by decide, by decide⟩)
      · exact Or.inl (by refine ⟨by decide, by decide, by decide⟩)
      · exact Or.inl (by refine ⟨by decide, by decide, by decide⟩)
      · exact Or.inl (by refine ⟨by decide, by decide, by decide⟩)
      · exact Or.inr ⟨21, by decide, by decide, by decide, by decide⟩
      · exact Or.inr ⟨13, by decide, by decide, by decide, by decide⟩
      · exact Or.inr ⟨8, by decide, by decide, by decide, by decide⟩
      · exact Or.inr ⟨5, by decide, by decide, by decide, by decide⟩
      · exact Or.inr ⟨3, by decide, by decide, by decide, by decide⟩
      · exact Or.inr ⟨2, by decide, by decide, by decide, by decide⟩
      · exact Or.inr ⟨1, by decide, by decide, by decide, by decide⟩
      · exact Or.inr ⟨0, by decide, by decide, by decide, by decide⟩)
  have : (13 : ℤ) * 10000 = 130000 := by norm_num
  have h2 : (1 : ℤ) * 10000 + (4 - 1 * 10000) % 13 = 10001 := by decide
  rw [this, h2] at h
  exact h

/-- The demo family has a strict witness at `10001/130000`. -/
theorem demoTwoScale_strictWitness : StrictWitness demoTwoScale :=
  strictWitness_of_strictlyLive demoTwoScale_strictlyLive

end TwoScaleWitness
end LonelyRunner

-- kernel-purity audit (fleet convention)
#print axioms LonelyRunner.TwoScaleWitness.small_speed_band
#print axioms LonelyRunner.TwoScaleWitness.cluster_speed_band
#print axioms LonelyRunner.TwoScaleWitness.twoScale_strictlyLive
#print axioms LonelyRunner.TwoScaleWitness.twoScale_strictWitness
#print axioms LonelyRunner.TwoScaleWitness.demoTwoScale_strictWitness
