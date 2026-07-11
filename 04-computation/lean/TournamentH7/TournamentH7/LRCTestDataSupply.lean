/-
  TournamentH7.LRCTestDataSupply -- THE TEST-DATA EXISTENCE + SUPPLY DISCHARGE
  (klein-2026-07-11-S248, HYP-5965, THM-696): the certificate-supply citation
  (mac-mini cont.27: THM527ACertificateSupply -- LRC(14) = THM-661 +
  pigeonhole + certificate existence) DISCHARGED on two-scale shapes,
  wholesale, with NO per-family data.

  (1) supply_of_strictlyLive -- THE FATTENING BRIDGE: any strictly-live
      ruler (Q, w) for speeds bounded by B yields the citation's EXACT data:

          D = 28·B·Q,  x = 28·B·w − 1,  y = 28·B·w + 1,  q = D,

      because the integer margins Q+1 <= 14r <= 13Q−1 fatten uniformly
      (28B·(14r − Q) >= 28B > 14B >= 14·v).  EVERY witness theorem
      (THM-693 two-scale, THM-694 multi-scale, THM-695 ray) now feeds
      THM527ACertificateSupply directly.

  (2) qstar_exists (P.card <= 5 gives some qs in [8,13] avoiding P: 6 > 5
      pigeonhole) + the a = 1 P-side (qs not in P forces p % qs != 0 for
      p in [1,13]) + missed_bridge ((c − x) % q != 0 iff x % q != c)
      => testData_exists: the witness data (qs, c) exists from the SHAPE
      alone.

  (3) twoScale_supply: for EVERY two-scale family -- small part P
      (card <= 5, listed in [1,13]), cluster co-offsets E with some
      qs in [8,13] avoiding P and |E| < qs (all non-top-block shapes),
      V > 2184 and 168e < V -- the citation's conclusion holds:
      exists D x y q with thirteen SafeIvStrict certificates.

  Kernel-pure: no native_decide, no sorry.
-/
import Mathlib
import TournamentH7.LRCTwoScaleWitness
import TournamentH7.LRCTestPointCore

namespace LonelyRunner
namespace TestDataSupply

open LRC14Grand MeasureTransfer

/-- **THE FATTENING BRIDGE**: a strictly-live ruler yields the certificate
data of the supply citation — the explicit interval
`[(28Bw−1)/(28BQ), (28Bw+1)/(28BQ)]` with `q = D`. -/
theorem supply_of_strictlyLive {v : Fin 13 → ℤ} {Q w B : ℤ}
    (h : StrictlyLive v Q w) (hB : ∀ i, 0 < v i ∧ v i ≤ B) (hB0 : 0 < B) :
    ∃ D x y q : ℤ, 0 < D ∧ 0 < q ∧ 0 ≤ x ∧ y < D ∧ D < (y - x) * q ∧
      ∀ i, SafeIvStrict (v i) D x y := by
  obtain ⟨hw0, hwQ, hband⟩ := h
  have hQ0 : (0 : ℤ) < Q := lt_trans hw0 hwQ
  have hD0 : (0 : ℤ) < 28 * B * Q := by positivity
  refine ⟨28 * B * Q, 28 * B * w - 1, 28 * B * w + 1, 28 * B * Q,
    hD0, hD0, by nlinarith, by nlinarith, by nlinarith, ?_⟩
  intro i
  obtain ⟨hb1, hb2⟩ := hband i
  obtain ⟨hv0, hvB⟩ := hB i
  set r : ℤ := (v i * w) % Q with hrdef
  set j : ℤ := (v i * w) / Q with hjdef
  have hdm : Q * j + r = v i * w := Int.ediv_add_emod _ _
  have h14lo : Q + 1 ≤ 14 * r := by omega
  have h14hi : 14 * r ≤ 13 * Q - 1 := by omega
  refine ⟨j, ?_, ?_⟩
  · have hval : v i * (28 * B * w - 1) - j * (28 * B * Q)
        = 28 * B * r - v i := by
      linear_combination (28 * B) * hdm.symm
    rw [hval]
    nlinarith [h14lo, hvB, hB0, hv0]
  · have hval : v i * (28 * B * w + 1) - j * (28 * B * Q)
        = 28 * B * r + v i := by
      linear_combination (28 * B) * hdm.symm
    rw [hval]
    nlinarith [h14hi, hvB, hB0, hv0]

/-- **The q\* existence**: a small part with at most five speeds leaves some
modulus in `[8,13]` untouched — six values, five speeds. -/
theorem qstar_exists (P : Finset ℤ) (hcard : P.card ≤ 5) :
    ∃ qs : ℤ, 8 ≤ qs ∧ qs ≤ 13 ∧ qs ∉ P := by
  by_contra h
  push_neg at h
  have hsub : Finset.Icc (8 : ℤ) 13 ⊆ P := by
    intro x hx
    rw [Finset.mem_Icc] at hx
    exact h x hx.1 hx.2
  have hle := Finset.card_le_card hsub
  rw [Int.card_Icc] at hle
  simp at hle
  omega

/-- The missed-class bridge: `x % q ≠ c` (the pigeonhole's output) gives the
witness theorems' hypothesis `(c − x) % q ≠ 0`. -/
theorem missed_bridge {q c x : ℤ} (hq : 0 < q) (hc0 : 0 ≤ c) (hcq : c < q)
    (h : x % q ≠ c) : (c - x) % q ≠ 0 := by
  intro h0
  apply h
  obtain ⟨k, hk⟩ := Int.dvd_of_emod_eq_zero h0
  have hx : x = c - q * k := by omega
  rw [hx, Int.sub_emod, Int.mul_emod_right, sub_zero, Int.emod_emod_of_dvd _ dvd_rfl,
    Int.emod_eq_of_lt hc0 hcq]

/-- **THE TEST-DATA EXISTENCE (the directive's lemma)**: from the shape alone
— small part of ≤ 5 speeds in `[1,13]`, some modulus in `[8,13]` avoiding it
with room for the cluster — the witness data `(qs, a = 1, c)` EXISTS. -/
theorem testData_exists (P E : Finset ℤ)
    (hPsub : ∀ p ∈ P, 1 ≤ p ∧ p ≤ 13)
    (hqsE : ∃ qs : ℤ, 8 ≤ qs ∧ qs ≤ 13 ∧ qs ∉ P ∧ (E.card : ℤ) < qs) :
    ∃ qs c : ℤ, 8 ≤ qs ∧ qs ≤ 13 ∧ 0 ≤ c ∧ c < qs ∧
      (∀ p ∈ P, (p * 1) % qs ≠ 0) ∧ (∀ e ∈ E, (c - e * 1) % qs ≠ 0) := by
  obtain ⟨qs, h8, h13, hnP, hEcard⟩ := hqsE
  obtain ⟨c, hc0, hcq, hmiss⟩ :=
    TestPointCore.pigeonhole_missed_class E qs 1 (by omega) hEcard
  refine ⟨qs, c, h8, h13, hc0, hcq, ?_, ?_⟩
  · intro p hp
    obtain ⟨hp1, hp13⟩ := hPsub p hp
    rw [mul_one]
    intro h0
    obtain ⟨k, hk⟩ := Int.dvd_of_emod_eq_zero h0
    have hpq : p = qs := by
      rcases le_or_gt k 0 with hle | hgt
      · exfalso; nlinarith [hk]
      · have hk1 : (1 : ℤ) ≤ k := hgt
        rcases eq_or_lt_of_le hk1 with heq | h2
        · rw [← heq, mul_one] at hk; exact hk
        · exfalso; nlinarith [hk]
    exact hnP (hpq ▸ hp)
  · intro e he
    exact missed_bridge (by omega) hc0 hcq (hmiss e he)

/-- **THE SUPPLY DISCHARGE ON TWO-SCALE SHAPES (THM-696)**: every two-scale
family — small part `P` (≤ 5 speeds in `[1,13]`), cluster `{V − e : e ∈ E}`
with some `qs ∈ [8,13]` avoiding `P` and `|E| < qs` (every non-top-block
shape), `V > 2184`, `168e < V` — satisfies the certificate-supply
conclusion: an interval with thirteen strict band certificates exists.
No per-family data; the citation's arithmetic content, proved wholesale. -/
theorem twoScale_supply (v : Fin 13 → ℤ) (V : ℤ) (P E : Finset ℤ)
    (hV : 2184 < V)
    (hPsub : ∀ p ∈ P, 1 ≤ p ∧ p ≤ 13)
    (hqsE : ∃ qs : ℤ, 8 ≤ qs ∧ qs ≤ 13 ∧ qs ∉ P ∧ (E.card : ℤ) < qs)
    (hE : ∀ e ∈ E, 0 ≤ e ∧ 168 * e < V)
    (hv : ∀ i, v i ∈ P ∨ ∃ e ∈ E, v i = V - e) :
    ∃ D x y q : ℤ, 0 < D ∧ 0 < q ∧ 0 ≤ x ∧ y < D ∧ D < (y - x) * q ∧
      ∀ i, SafeIvStrict (v i) D x y := by
  obtain ⟨qs, c, h8, h13, hc0, hcq, hPside, hEside⟩ :=
    testData_exists P E hPsub hqsE
  have hlive := TwoScaleWitness.twoScale_strictlyLive v qs 1 c V
    h8 h13 (by norm_num) (by omega) hV (by
      intro i
      rcases hv i with hp | ⟨e, he, hve⟩
      · exact Or.inl ⟨(hPsub _ hp).1, (hPsub _ hp).2, hPside _ hp⟩
      · exact Or.inr ⟨e, (hE e he).1, hve, hEside e he, (hE e he).2⟩)
  refine supply_of_strictlyLive hlive (fun i => ?_) (by omega : (0 : ℤ) < V)
  rcases hv i with hp | ⟨e, he, hve⟩
  · obtain ⟨h1, h13'⟩ := hPsub _ hp
    exact ⟨by omega, by omega⟩
  · obtain ⟨he0, heV⟩ := hE e he
    rw [hve]
    exact ⟨by omega, by omega⟩

end TestDataSupply
end LonelyRunner

-- kernel-purity audit (fleet convention)
#print axioms LonelyRunner.TestDataSupply.supply_of_strictlyLive
#print axioms LonelyRunner.TestDataSupply.qstar_exists
#print axioms LonelyRunner.TestDataSupply.missed_bridge
#print axioms LonelyRunner.TestDataSupply.testData_exists
#print axioms LonelyRunner.TestDataSupply.twoScale_supply
