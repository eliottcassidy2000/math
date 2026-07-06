/-
  TournamentH7.LRCHarmonicAP — THE AP IS THE HARMONIC (FLAT) FAMILY: vanishing
  second differences CHARACTERISE the arithmetic progression
  (kind-pasteur-2026-07-06-S30).

  The density floor is the rigidity of the AP as the unique extremal (unique M-
  minimiser / unique `safe = 0` tiler).  Across every lens the AP is the "flat"
  configuration: flat eigenvalue spectrum (Paley, THM-126, kps S29), equi-
  oscillation (mac-mini HYP-4462), min discrepancy (opus HYP-4074), maximal
  relation lattice / additive energy (opus HYP-4446).  This file gives the
  ELEMENTARY algebraic heart of all of them:

  > a family has vanishing SECOND DIFFERENCES  `v(i+2) − 2v(i+1) + v(i) = 0`
  >   ⟺  it is an arithmetic progression  `v(i) = v(0) + i·d`.

  The second-difference operator is the discrete Laplacian; its kernel is the
  "harmonic" (flat) sequences, and this theorem says the harmonic sequences are
  EXACTLY the APs.  The relations `e_i − 2 e_{i+1} + e_{i+2}` are therefore in the
  AP's relation lattice `L(AP)` (opus S112) and CHARACTERISE it — the length-3,
  coefficient-`(1,−2,1)` relations, the shortest nontrivial ones, hence (since
  `|ĥ_m| ~ 1/m`) the dominant terms in the `safe` theta-sum.  The AP maximises
  these harmonic relations; that is the additive-energy heart of `safe(AP) = 0`.

  Kernel-pure; no `sorry`, no `native_decide`.
-/
import Mathlib

namespace LonelyRunner
namespace HarmonicAP

/-- **The AP is harmonic**: an arithmetic progression `v(i) = a + i·d` has
vanishing second differences. -/
theorem ap_second_diff_zero (a d : ℤ) (i : ℕ) :
    (a + ((i : ℤ) + 2) * d) - 2 * (a + ((i : ℤ) + 1) * d) + (a + (i : ℤ) * d) = 0 := by
  ring

/-- **Harmonic ⟹ AP** (the characterisation): a sequence with vanishing second
differences everywhere is an arithmetic progression with common difference
`v(1) − v(0)`.  Proved by two-step (paired) induction — the discrete
Laplacian's kernel is exactly the linear sequences. -/
theorem ap_of_second_diff_zero (v : ℕ → ℤ)
    (h : ∀ i, v (i + 2) - 2 * v (i + 1) + v i = 0) :
    ∀ n, v n = v 0 + (n : ℤ) * (v 1 - v 0) := by
  have key : ∀ n, v n = v 0 + (n : ℤ) * (v 1 - v 0)
      ∧ v (n + 1) = v 0 + ((n : ℤ) + 1) * (v 1 - v 0) := by
    intro n
    induction n with
    | zero => exact ⟨by ring, by push_cast; ring⟩
    | succ k ih =>
      obtain ⟨ih1, ih2⟩ := ih
      refine ⟨?_, ?_⟩
      · -- v (k+1) = v0 + (k+1) d  — this is ih2
        push_cast; linarith [ih2]
      · -- v (k+2) = 2 v(k+1) − v k = v0 + (k+2) d
        have hk : v (k + 2) = 2 * v (k + 1) - v k := by linarith [h k]
        rw [hk, ih1, ih2]; push_cast; ring
  exact fun n => (key n).1

/-- **The harmonic characterisation of the AP** (both directions): a sequence has
vanishing second differences iff it is an arithmetic progression.  The kernel of
the discrete Laplacian is exactly the APs — the elementary heart of "the AP is
the flat/extremal family" behind the density floor. -/
theorem second_diff_zero_iff_ap (v : ℕ → ℤ) :
    (∀ i, v (i + 2) - 2 * v (i + 1) + v i = 0)
      ↔ ∃ a d : ℤ, ∀ n, v n = a + (n : ℤ) * d := by
  constructor
  · intro h
    exact ⟨v 0, v 1 - v 0, ap_of_second_diff_zero v h⟩
  · rintro ⟨a, d, hv⟩ i
    rw [hv (i + 2), hv (i + 1), hv i]
    push_cast; ring

#print axioms ap_second_diff_zero
#print axioms ap_of_second_diff_zero
#print axioms second_diff_zero_iff_ap

end HarmonicAP
end LonelyRunner
