/-
  TournamentH7.LRCKStratification — THE k-STRATIFICATION LEMMA CORE
  (mac-mini-2026-07-06-S1, HYP-4232; composes with kps-S17's cluster-gcd
  ladder HYP-4217 for the 3/38-cell k-reduction).

  A hypothetical 3/38-attainer's maximizer sits on a 38k-grid (grid_div_38:
  the binder's integer grid-distance 3q/38 forces 38 | q); its binders are
  all divisible by k (binder_dvd: v·m ≡ ±3k mod 38k with gcd(m,k) = 1); and
  the k-multiples' quotients inherit the mod-38 level-3 structure through
  the exact scaling identity dInt (38k) (k·x) = k · dInt 38 x (dInt_scale).
  With the cluster-gcd ladder bounding gcd-clusters, the k ≥ 2 strata are
  height-bounded relative to the non-multiples (prose: the S1 draft); the
  k = 1 base has its binding pair summing to EXACTLY 38.

  Pure integer arithmetic; no analysis.  Draft with the full composition:
  03-artifacts/drafts/k-stratification-lemma-macmini-S1.md.
-/
import Mathlib.RingTheory.Int.Basic
import Mathlib.Tactic

namespace LonelyRunner
namespace KStratification

/-- Integer circle distance at modulus q (canonical residue form). -/
def dInt (q x : ℤ) : ℤ := min (x % q) (q - x % q)

/-- THE SCALING IDENTITY: distances at modulus 38k of k-multiples are k times
    the mod-38 distances of the quotients.  (The quotient-descent engine.) -/
theorem dInt_scale (k x : ℤ) (hk : 0 < k) :
    dInt (38 * k) (k * x) = k * dInt 38 x := by
  unfold dInt
  have h1 : k * x % (38 * k) = k * (x % 38) := by
    rw [show (38 : ℤ) * k = k * 38 by ring]
    exact Int.mul_emod_mul_of_pos x 38 hk
  rw [h1, show (38 : ℤ) * k - k * (x % 38) = k * (38 - x % 38) by ring]
  rcases le_total (x % 38) (38 - x % 38) with h | h
  · rw [min_eq_left (by nlinarith), min_eq_left h]
  · rw [min_eq_right (by nlinarith), min_eq_right h]

/-- BINDER DIVISIBILITY: a runner binding at ±3k on the 38k-grid under a
    dilation coprime to k is itself divisible by k. -/
theorem binder_dvd (v m k : ℤ) (hco : IsCoprime k m)
    (h : (38 * k : ℤ) ∣ (v * m - 3 * k) ∨ (38 * k : ℤ) ∣ (v * m + 3 * k)) :
    k ∣ v := by
  have hkdvd : (k : ℤ) ∣ 38 * k := ⟨38, by ring⟩
  have hkvm : k ∣ v * m := by
    rcases h with h | h
    · have hk1 : k ∣ v * m - 3 * k := hkdvd.trans h
      have h2 : v * m = (v * m - 3 * k) + 3 * k := by ring
      rw [h2]
      exact dvd_add hk1 ⟨3, by ring⟩
    · have hk1 : k ∣ v * m + 3 * k := hkdvd.trans h
      have h2 : v * m = (v * m + 3 * k) - 3 * k := by ring
      rw [h2]
      exact dvd_sub hk1 ⟨3, by ring⟩
  exact hco.dvd_of_dvd_mul_right hkvm

/-- GRID DIVISIBILITY: if the binder's grid distance 3q/38 is an integer
    (38 ∣ 3q), then 38 ∣ q — the attaining grid is a 38k-grid. -/
theorem grid_div_38 (q : ℤ) (h : (38 : ℤ) ∣ 3 * q) : (38 : ℤ) ∣ q := by
  omega

/-- The k = 1 binding-pair anchor: a pair binding at +3 and −3 under the same
    dilation on the 38-grid sums to a multiple of 38 (with both positive and
    below 38·2, the sum is exactly 38 — the finite base of the cell attack). -/
theorem pair_sum_dvd (vi vj m : ℤ)
    (hi : (38 : ℤ) ∣ (vi * m - 3)) (hj : (38 : ℤ) ∣ (vj * m + 3))
    (hm : IsCoprime (38 : ℤ) m) :
    (38 : ℤ) ∣ (vi + vj) := by
  have hsum : (38 : ℤ) ∣ ((vi + vj) * m) := by
    have hadd := dvd_add hi hj
    have heq : (vi * m - 3) + (vj * m + 3) = (vi + vj) * m := by ring
    rwa [heq] at hadd
  exact hm.dvd_of_dvd_mul_right hsum

end KStratification
end LonelyRunner
