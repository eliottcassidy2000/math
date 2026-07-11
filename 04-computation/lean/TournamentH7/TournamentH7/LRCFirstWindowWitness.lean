/-
  TournamentH7.LRCFirstWindowWitness -- THE FIRST-WINDOW WITNESS
  (klein-2026-07-11-S249, HYP-5970, THM-697): the packed TOP-BLOCK supply gap
  (S248's remainder (a)) closed by a route SIMPLER than the test points.

  THE WITNESS: t = (7j+6)/(7V) -- the multiplier w = 7j+6 puts the top
  cluster's phase at EXACTLY 6/7:

      V·(7j+6) mod 7V = 6V,   so   ((V−e)·w) mod (7V) = 6V − e·w,

  and ONE inequality per cluster speed (2·e·w < 11·V) keeps it strictly
  in-band.  Small speeds ride the OPEN FIRST WINDOW W_P = (1/(14·pmin),
  13/(14·pmax)) -- nondegenerate for EVERY P with pmax + 1 ≤ 13·pmin (all
  top-blocks [b,13]: 13(b−1) ≥ 104 of room) -- with the no-wrap band being
  the hypothesis itself.  NO test point, NO missed class: the window does
  the P-side and the 6/7-phase does the cluster.

  Coverage: e ≤ 10·pmin with the leftmost j (nearly THM-691(B)'s packed
  range).  Complementarity: P's containing {1, 13} (degenerate window) are
  exactly the NON-top-blocks that THM-693/696's test-point route covers --
  the two witnesses together cover every small part with packed clusters.

  DEMO: the k = 8 top block {9,…,13} with the CONSECUTIVE APEX cluster
  E = {0,…,7} at V = 10000 -- the extremal shape of the whole program --
  strictly live at (q, w) = (70000, 559).

  Kernel-pure: no native_decide, no sorry.
-/
import Mathlib
import TournamentH7.LRCTestDataSupply

namespace LonelyRunner
namespace FirstWindowWitness

open LRC14Grand MeasureTransfer

/-- **The 6/7-phase cluster band**: with `w = 7j+6`, the cluster residue is
`6V − e·w`, strictly in-band once `2·e·w < 11·V`. -/
theorem cluster_band {V e j : ℤ} (hV : 0 < V) (he : 0 ≤ e) (hj : 0 ≤ j)
    (hew : 2 * e * (7 * j + 6) < 11 * V) :
    7 * V < 14 * (((V - e) * (7 * j + 6)) % (7 * V)) ∧
    14 * (((V - e) * (7 * j + 6)) % (7 * V)) < 13 * (7 * V) := by
  have hw0 : (0 : ℤ) < 7 * j + 6 := by omega
  have hkey : (V - e) * (7 * j + 6)
      = (6 * V - e * (7 * j + 6)) + (7 * V) * j := by ring
  have hb0 : 0 ≤ 6 * V - e * (7 * j + 6) := by nlinarith
  have hbq : 6 * V - e * (7 * j + 6) < 7 * V := by nlinarith
  have hmod : ((V - e) * (7 * j + 6)) % (7 * V)
      = 6 * V - e * (7 * j + 6) := by
    rw [hkey, Int.add_mul_emod_self_left, Int.emod_eq_of_lt hb0 hbq]
  rw [hmod]
  constructor
  · nlinarith
  · nlinarith

/-- **THE FIRST-WINDOW WITNESS (THM-697)**: per speed, either the no-wrap
band (the small side — the window inclusion, decidable) or the packed
cluster inequality; strictly live at `(7V, 7j+6)`. -/
theorem firstWindow_strictlyLive (v : Fin 13 → ℤ) (V j : ℤ)
    (hV : 0 < V) (hj : 0 ≤ j) (hwq : 7 * j + 6 < 7 * V)
    (hv : ∀ i, (7 * V < 14 * (v i * (7 * j + 6)) ∧
                14 * (v i * (7 * j + 6)) < 13 * (7 * V)) ∨
          (∃ e, 0 ≤ e ∧ v i = V - e ∧ 2 * e * (7 * j + 6) < 11 * V)) :
    StrictlyLive v (7 * V) (7 * j + 6) := by
  refine ⟨by omega, hwq, ?_⟩
  intro i
  rcases hv i with ⟨h1, h2⟩ | ⟨e, he, hve, hew⟩
  · have h0 : 0 ≤ v i * (7 * j + 6) := by nlinarith
    have hq : v i * (7 * j + 6) < 7 * V := by nlinarith
    rw [Int.emod_eq_of_lt h0 hq]
    exact ⟨h1, h2⟩
  · rw [hve]
    exact cluster_band hV he hj hew

/-- **The leftmost-j existence**: the first `w ≡ 6 (mod 7)` past the window's
left edge exists and is bracketed — `7V < 14·pmin·w` and `w ≤ V/(2·pmin)+7`
in the scaled form `14·pmin·w ≤ 7V + 112·pmin`. -/
theorem firstWindow_j_exists (V pmin : ℤ) (hpm : 1 ≤ pmin) (hV : 0 < V) :
    ∃ j : ℤ, 0 ≤ j ∧
      7 * V < 14 * pmin * (7 * j + 6) ∧
      14 * pmin * (7 * j + 6) ≤ 7 * V + 112 * pmin := by
  set w0 : ℤ := (7 * V) / (14 * pmin) + 1 with hw0def
  have hpq0 : (0 : ℤ) < 14 * pmin := by omega
  have hdm := Int.ediv_add_emod (7 * V) (14 * pmin)
  have hm0 := Int.emod_nonneg (7 * V) hpq0.ne'
  have hmD := Int.emod_lt_of_pos (7 * V) hpq0
  have hw0lo : 7 * V < 14 * pmin * w0 := by
    rw [hw0def]
    nlinarith [hdm, hmD]
  have hw0hi : 14 * pmin * w0 ≤ 7 * V + 14 * pmin := by
    rw [hw0def]
    nlinarith [hdm, hm0]
  -- round w0 up to ≡ 6 (mod 7): w := w0 + ((6 - w0) % 7) ∈ [w0, w0 + 6]
  set d : ℤ := (6 - w0) % 7 with hddef
  have hd0 : 0 ≤ d := Int.emod_nonneg _ (by norm_num)
  have hd7 : d < 7 := Int.emod_lt_of_pos _ (by norm_num)
  have hd6 : (w0 + d) % 7 = 6 % 7 := by
    have : w0 + d = 6 - 7 * ((6 - w0) / 7) := by
      have := Int.ediv_add_emod (6 - w0) 7
      omega
    rw [this, Int.sub_emod, Int.mul_emod_right, sub_zero,
      Int.emod_emod_of_dvd _ dvd_rfl]
  obtain ⟨j, hjw⟩ : ∃ j : ℤ, w0 + d = 7 * j + 6 := by
    have h6 : (w0 + d) % 7 = 6 := by rw [hd6]; norm_num
    refine ⟨(w0 + d) / 7, ?_⟩
    have := Int.ediv_add_emod (w0 + d) 7
    omega
  have hw0pos : 0 < w0 := by nlinarith [hw0lo]
  refine ⟨j, by omega, ?_, ?_⟩
  · rw [← hjw]
    nlinarith [hw0lo, hd0]
  · rw [← hjw]
    nlinarith [hw0hi, hd7]

/-! ## Demo: the top-block apex, strictly live through the theorem -/

def demoTopBlock : Fin 13 → ℤ :=
  ![9, 10, 11, 12, 13, 9993, 9994, 9995, 9996, 9997, 9998, 9999, 10000]

theorem demoTopBlock_strictlyLive :
    StrictlyLive demoTopBlock 70000 559 := by
  have h := firstWindow_strictlyLive demoTopBlock 10000 79
    (by norm_num) (by norm_num) (by norm_num) ?_
  · have e1 : (7 : ℤ) * 10000 = 70000 := by norm_num
    have e2 : (7 : ℤ) * 79 + 6 = 559 := by norm_num
    rw [e1, e2] at h
    exact h
  intro i
  fin_cases i
  · exact Or.inl ⟨by decide, by decide⟩
  · exact Or.inl ⟨by decide, by decide⟩
  · exact Or.inl ⟨by decide, by decide⟩
  · exact Or.inl ⟨by decide, by decide⟩
  · exact Or.inl ⟨by decide, by decide⟩
  · exact Or.inr ⟨7, by decide, by decide, by decide⟩
  · exact Or.inr ⟨6, by decide, by decide, by decide⟩
  · exact Or.inr ⟨5, by decide, by decide, by decide⟩
  · exact Or.inr ⟨4, by decide, by decide, by decide⟩
  · exact Or.inr ⟨3, by decide, by decide, by decide⟩
  · exact Or.inr ⟨2, by decide, by decide, by decide⟩
  · exact Or.inr ⟨1, by decide, by decide, by decide⟩
  · exact Or.inr ⟨0, by decide, by decide, by decide⟩

/-- The top-block apex family satisfies the certificate-supply conclusion —
via the fattening bridge, from one first-window witness. -/
theorem demoTopBlock_supply :
    ∃ D x y q : ℤ, 0 < D ∧ 0 < q ∧ 0 ≤ x ∧ y < D ∧ D < (y - x) * q ∧
      ∀ i, SafeIvStrict (demoTopBlock i) D x y :=
  TestDataSupply.supply_of_strictlyLive demoTopBlock_strictlyLive
    (by intro i; fin_cases i <;> exact ⟨by decide, by decide⟩)
    (by norm_num : (0 : ℤ) < 10000)

end FirstWindowWitness
end LonelyRunner

-- kernel-purity audit (fleet convention)
#print axioms LonelyRunner.FirstWindowWitness.cluster_band
#print axioms LonelyRunner.FirstWindowWitness.firstWindow_strictlyLive
#print axioms LonelyRunner.FirstWindowWitness.firstWindow_j_exists
#print axioms LonelyRunner.FirstWindowWitness.demoTopBlock_strictlyLive
#print axioms LonelyRunner.FirstWindowWitness.demoTopBlock_supply
