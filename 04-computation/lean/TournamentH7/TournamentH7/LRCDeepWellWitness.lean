/-
  TournamentH7.LRCDeepWellWitness — THE COVERING-MIN EXTREMIZER'S WITNESS, GENERAL n
  (klein-2026-07-03-S119, HYP-4065; renumbered from HYP-4062 -- collided with kind-pasteur-S37).

  The LRC open core is the covering-min lower bound `M(covering) >= n/Phi6(n)` (THM-523),
  whose unique extremizer at `n = 14` is the DEEP WELL `D_n = {1,...,n-2} U {n(n-1)}`
  (`= {1..12, 182}`), with `M(D_n) = n/Phi6(n) = 14/183` (opus: the witness is the zeta_6
  hexagonal rotation `mult-by-n` on `Z[omega]/(n-omega) = Z/Phi6(n)`, `Phi6(n)=n^2-n+1`).

  This file FORMALIZES the WITNESS side for GENERAL `n >= 3`, sorry-free: at the cyclotomic
  time `t* = n/Phi6(n)`, every runner of the deep well is at distance `>= n/Phi6(n)` from the
  integers, hence `>= n/Phi6(n) > 1/n` — so the deep well is lonely with margin
  `(n-1)/(n*Phi6(n)) > 0`.  The whole proof is the Eisenstein "tower" identity
  `n^2 = n - 1  (mod Phi6(n))`, i.e. `Phi6(n) | n^3 - n^2 + n = n * Phi6(n)`, plus elementary
  band arithmetic — pinning the covering-min TARGET value rigorously.

  Scope (honest): this is the WITNESS (upper) direction for the extremal family.  The
  covering-min LOWER bound over ALL covering families (the rigidity `M=1/n => tight locus`)
  is the remaining open crux and is NOT proved here.
-/
import Mathlib

namespace LonelyRunner.DeepWell

/-- `Phi6(n) = n^2 - n + 1`: the 6th-cyclotomic value = Eisenstein norm `N(n - omega)`
  = `|PG(2, n-1)|`.  The witness denominator `q* = Phi6(n)`. -/
def phi6 (n : ℤ) : ℤ := n ^ 2 - n + 1

theorem phi6_pos {n : ℤ} (hn : 2 ≤ n) : 0 < phi6 n := by
  unfold phi6; nlinarith

/-- The Eisenstein tower divisibility `Phi6(n) | n^3 - n^2 + n` (the factored form
  `n^3 - n^2 + n = n * Phi6(n)` IS the identity `n^2 = n - 1 (mod Phi6)`). -/
theorem phi6_dvd_tower (n : ℤ) : phi6 n ∣ (n ^ 3 - n ^ 2 + n) :=
  ⟨n, by unfold phi6; ring⟩

/-- **Band lemma.**  If `x mod q = r` with `q > 0`, `0 <= n`, `n <= r <= q - n`, and `n <= q`,
  then `x` is at distance `>= n` from EVERY integer multiple of `q`: `n <= |x - m*q|`. -/
theorem dist_multiples_ge {x q r n : ℤ} (hq : 0 < q) (hr : x % q = r) (hn0 : 0 ≤ n)
    (hlo : n ≤ r) (hhi : r ≤ q - n) (hnq : n ≤ q) :
    ∀ m : ℤ, n ≤ |x - m * q| := by
  intro m
  have hx : x = q * (x / q) + r := by
    have := Int.ediv_add_emod x q; rw [hr] at this; linarith
  set k : ℤ := x / q with hk
  have hrw : x - m * q = q * (k - m) + r := by rw [hx]; ring
  rw [hrw]
  rcases lt_trichotomy (k - m) 0 with hd | hd | hd
  · -- k - m ≤ -1 : q*(k-m)+r ≤ -q + (q-n) = -n
    have h1 : k - m ≤ -1 := by omega
    have h2 : q * (k - m) ≤ q * (-1) := mul_le_mul_of_nonneg_left h1 (le_of_lt hq)
    exact le_abs.mpr (Or.inr (by nlinarith))
  · -- k - m = 0 : |r| = r ≥ n
    rw [hd]; simp only [mul_zero, zero_add]
    exact le_abs.mpr (Or.inl hlo)
  · -- k - m ≥ 1 : q*(k-m)+r ≥ q ≥ n
    have h1 : (1 : ℤ) ≤ k - m := by omega
    have h2 : q * 1 ≤ q * (k - m) := mul_le_mul_of_nonneg_left h1 (le_of_lt hq)
    exact le_abs.mpr (Or.inl (by nlinarith))

/-- **AP-runner band.**  For `1 <= j <= n-2` and `n >= 3`, the AP runner `j` lands, under the
  witness `mult-by-n`, at residue `j*n in [n, Phi6-n]` (distance `>= n` from 0). -/
theorem ap_runner_band {n j : ℤ} (hn : 3 ≤ n) (hj1 : 1 ≤ j) (hj2 : j ≤ n - 2) :
    n ≤ (j * n) % phi6 n ∧ (j * n) % phi6 n ≤ phi6 n - n := by
  have hq : 0 < phi6 n := phi6_pos (by linarith)
  have hnn : 0 ≤ j * n := mul_nonneg (by linarith) (by linarith)
  have hlt : j * n < phi6 n := by unfold phi6; nlinarith
  have hmod : (j * n) % phi6 n = j * n := Int.emod_eq_of_lt hnn hlt
  rw [hmod]
  refine ⟨by nlinarith, ?_⟩
  unfold phi6; nlinarith

/-- **Defect-runner band.**  The pronic defect `n(n-1)` lands, under `mult-by-n`, at residue
  `Phi6-n` (distance exactly `n` from 0) — the AP tail `-n`, via `n^2 = n-1 (mod Phi6)`. -/
theorem defect_runner_band {n : ℤ} (hn : 3 ≤ n) :
    (n * (n - 1) * n) % phi6 n = phi6 n - n := by
  have hq : 0 < phi6 n := phi6_pos (by linarith)
  have hcong : (n * (n - 1) * n) % phi6 n = (phi6 n - n) % phi6 n := by
    have hdvd : phi6 n ∣ (phi6 n - n) - n * (n - 1) * n := ⟨1 - n, by unfold phi6; ring⟩
    exact Int.ModEq.symm (Int.modEq_iff_dvd.mpr hdvd) |>.symm
  rw [hcong]
  apply Int.emod_eq_of_lt
  · unfold phi6; nlinarith
  · unfold phi6; nlinarith

/-- Deep-well AP body is lonely under `mult-by-n`: distance `>= n` to every multiple of Phi6. -/
theorem ap_runner_lonely {n j : ℤ} (hn : 3 ≤ n) (hj1 : 1 ≤ j) (hj2 : j ≤ n - 2) :
    ∀ m : ℤ, n ≤ |j * n - m * phi6 n| := by
  obtain ⟨hlo, hhi⟩ := ap_runner_band hn hj1 hj2
  exact dist_multiples_ge (phi6_pos (by linarith)) rfl (by linarith) hlo hhi
    (by unfold phi6; nlinarith)

/-- Deep-well pronic defect is lonely under `mult-by-n`: distance `>= n` to every multiple. -/
theorem defect_runner_lonely {n : ℤ} (hn : 3 ≤ n) :
    ∀ m : ℤ, n ≤ |n * (n - 1) * n - m * phi6 n| := by
  refine dist_multiples_ge (phi6_pos (by linarith)) (defect_runner_band hn) (by linarith)
    ?_ ?_ ?_
  · unfold phi6; nlinarith
  · linarith
  · unfold phi6; nlinarith

/-- **The deep well is lonely at `t* = n/Phi6(n)` (real form, AP body).**  Every AP runner `j`
  is at real distance `>= n/Phi6(n)` from every integer, at the cyclotomic witness time. -/
theorem ap_runner_dist_real {n j : ℤ} (hn : 3 ≤ n) (hj1 : 1 ≤ j) (hj2 : j ≤ n - 2) (m : ℤ) :
    (n : ℝ) / phi6 n ≤ |(j : ℝ) * ((n : ℝ) / phi6 n) - m| := by
  have hqZ : 0 < phi6 n := phi6_pos (by linarith)
  have hqR : (0 : ℝ) < (phi6 n : ℝ) := by exact_mod_cast hqZ
  have hint : n ≤ |j * n - m * phi6 n| := ap_runner_lonely hn hj1 hj2 m
  have key : (n : ℝ) ≤ |(j : ℝ) * n - m * (phi6 n : ℝ)| := by
    have h : (n : ℝ) ≤ ((|j * n - m * phi6 n| : ℤ) : ℝ) := by exact_mod_cast hint
    rw [Int.cast_abs] at h; push_cast at h; exact h
  have hfac : (j : ℝ) * ((n : ℝ) / phi6 n) - m
      = ((j : ℝ) * n - m * (phi6 n : ℝ)) / (phi6 n : ℝ) := by field_simp
  rw [hfac, abs_div, abs_of_pos hqR, le_div_iff₀ hqR]
  have hmul : (n : ℝ) / phi6 n * phi6 n = (n : ℝ) := by field_simp
  rw [hmul]; exact key

/-- The witness threshold beats `1/n`: `n/Phi6(n) > 1/n` for `n >= 2` (margin `(n-1)/(n Phi6)`),
  since `Phi6(n) = n^2 - n + 1 < n^2`. -/
theorem witness_gt_threshold {n : ℤ} (hn : 2 ≤ n) :
    (1 : ℝ) / n < (n : ℝ) / phi6 n := by
  have hqZ : 0 < phi6 n := phi6_pos hn
  have hqR : (0 : ℝ) < (phi6 n : ℝ) := by exact_mod_cast hqZ
  have hnR : (0 : ℝ) < (n : ℝ) := by exact_mod_cast (show (0 : ℤ) < n by linarith)
  have hkey : (phi6 n : ℝ) < (n : ℝ) * n := by
    have : phi6 n < n * n := by unfold phi6; nlinarith
    exact_mod_cast this
  rw [div_lt_div_iff₀ hnR hqR]; nlinarith [hkey]

#print axioms dist_multiples_ge
#print axioms defect_runner_band
#print axioms ap_runner_dist_real
#print axioms witness_gt_threshold

end LonelyRunner.DeepWell
