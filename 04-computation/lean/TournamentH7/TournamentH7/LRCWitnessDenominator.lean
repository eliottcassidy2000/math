/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-06-S109)
-/
import TournamentH7.LRCGridAttainment

/-!
# The witness-denominator lever: (G) is a finite check modulo one height bound (HYP-4416)

mac-mini-S16 (HYP-4432) proved on paper: `M(S) = c/q` (lowest terms) ⟹ `q ∣ (vᵢ ± vⱼ)` for
some pair, so `q ≤ 2·max|vᵢ|` — bounding height bounds the witness denominator, and (G)
collapses to a FINITE check.  This file formalizes the keystone: at a positive maximizer, the
maximal margin, times some pair-sum `|vᵢ| + |vⱼ|` (the merge-grid denominator, THM-592), is an
INTEGER.  Hence if `M = c/q` in lowest terms then `q ∣ (|vᵢ| + |vⱼ|)`.

The proof is one line of grid attainment: at the maximizer `t*`, `(|vᵢ|+|vⱼ|)·t* = m ∈ ℤ`
(`maximizer_on_grid`), and the margin is achieved by some runner `l₀`, so
`(|vᵢ|+|vⱼ|)·distZ(v_{l₀} t*) = |v_{l₀}·m − (|vᵢ|+|vⱼ|)·round(v_{l₀} t*)| ∈ ℤ`.
-/

namespace LonelyRunner
namespace WitnessDenominator

open GridAttainment TournamentH7.LRCWitness

variable {k : ℕ} [Nonempty (Fin k)]

/-- **The witness-denominator lever (integer form).**  At a positive global maximizer, the
maximal margin times some merge-grid pair-sum `|vᵢ| + |vⱼ|` is an integer. -/
theorem pairsum_mul_margin_int (v : Fin k → ℤ) (hv : ∀ i, v i ≠ 0) (tstar : ℝ)
    (hmax : ∀ t, margin v t ≤ margin v tstar) (hM0 : 0 < margin v tstar) :
    ∃ (i j : Fin k) (a : ℤ), ((|v i| + |v j| : ℤ) : ℝ) * margin v tstar = (a : ℝ) := by
  obtain ⟨i, j, m, hgrid⟩ := maximizer_on_grid v hv tstar hmax hM0
  obtain ⟨l0, -, hl0⟩ :=
    Finset.exists_mem_eq_inf' (Finset.univ_nonempty (α := Fin k))
      (fun l => distZ ((v l : ℝ) * tstar))
  have hmarg : margin v tstar = distZ ((v l0 : ℝ) * tstar) := hl0
  refine ⟨i, j, |v l0 * m - (|v i| + |v j|) * round ((v l0 : ℝ) * tstar)|, ?_⟩
  rw [hmarg, distZ_eq_round]
  have hDnn : (0 : ℝ) ≤ ((|v i| + |v j| : ℤ) : ℝ) := by positivity
  rw [show ((|v i| + |v j| : ℤ) : ℝ)
        * |(v l0 : ℝ) * tstar - round ((v l0 : ℝ) * tstar)|
      = |((|v i| + |v j| : ℤ) : ℝ)
        * ((v l0 : ℝ) * tstar - round ((v l0 : ℝ) * tstar))|
      from by rw [abs_mul, abs_of_nonneg hDnn]]
  have hkey : ((|v i| + |v j| : ℤ) : ℝ)
        * ((v l0 : ℝ) * tstar - round ((v l0 : ℝ) * tstar))
      = ((v l0 * m - (|v i| + |v j|) * round ((v l0 : ℝ) * tstar) : ℤ) : ℝ) := by
    rw [Int.cast_sub, Int.cast_mul, Int.cast_mul]
    linear_combination (v l0 : ℝ) * hgrid
  rw [hkey, ← Int.cast_abs]

/-- **The witness-denominator lever (divisibility form).**  If the maximal margin equals
`c/q` in lowest terms, then `q` divides some merge-grid pair-sum `|vᵢ| + |vⱼ|` — hence
`q ≤ 2·max|vᵢ|`.  Bounding height bounds the witness denominator, so (G) is a FINITE check
modulo the single remaining height bound. -/
theorem denom_dvd_pairsum (v : Fin k → ℤ) (hv : ∀ i, v i ≠ 0) (tstar : ℝ)
    (hmax : ∀ t, margin v t ≤ margin v tstar) (hM0 : 0 < margin v tstar)
    (c : ℤ) (q : ℕ) (hq : 0 < q) (hcop : Nat.Coprime c.natAbs q)
    (hM : margin v tstar = (c : ℝ) / q) :
    ∃ i j : Fin k, (q : ℤ) ∣ (|v i| + |v j|) := by
  obtain ⟨i, j, a, ha⟩ := pairsum_mul_margin_int v hv tstar hmax hM0
  refine ⟨i, j, ?_⟩
  rw [hM] at ha
  have hqR : (q : ℝ) ≠ 0 := by exact_mod_cast hq.ne'
  -- (|vᵢ|+|vⱼ|)·c = a·q  in ℤ
  have hDcZ : (|v i| + |v j|) * c = a * q := by
    have : ((|v i| + |v j| : ℤ) : ℝ) * c = (a : ℝ) * q := by
      field_simp at ha; linarith [ha]
    exact_mod_cast this
  -- q ∣ (|vᵢ|+|vⱼ|)·c, and gcd(q, c) = 1, so q ∣ (|vᵢ|+|vⱼ|)
  have hdvd : (q : ℤ) ∣ (|v i| + |v j|) * c := ⟨a, by linarith [hDcZ]⟩
  have hcop' : IsCoprime (q : ℤ) c := by
    rw [Int.isCoprime_iff_gcd_eq_one]
    simpa [Int.gcd, Nat.Coprime, Nat.gcd_comm] using hcop
  exact hcop'.dvd_of_dvd_mul_right hdvd

#print axioms pairsum_mul_margin_int
#print axioms denom_dvd_pairsum

end WitnessDenominator
end LonelyRunner
