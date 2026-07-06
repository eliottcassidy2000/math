/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-06-S102)
-/
import Mathlib

/-!
# The rigidity heart of the (A) ⟸ (C) reduction (HYP-4346)

The (A) leg reduces to (C) via the rigidity "all full-support projections tight ⟹ rank ≤ 1"
(`the-A-leg-reduces-to-the-C-leg-opus-S101.md`).  Its algebraic core, formalized here:

  if two INDEPENDENT linear combinations of `u, v` are BOTH proportional to the SAME vector
  `w`, then `u` and `v` are linearly dependent.

In the rigidity, two full-support directions `(a,b) ≠ (a',b')` whose tight projections share
the same dilated-AP ORDERING both write `a·u + b·v = λ·w` and `a'·u + b'·v = λ'·w` for the
same ordering vector `w` — so by this lemma `u, v` are dependent, i.e. the 2-torus is rank 1.
The pigeonhole that produces the two same-ordering directions (finitely many orderings over
`Sym 12`, infinitely many primitive directions) is the combinatorial wrapper; this is the
2×2 solve underneath it.
-/

namespace LonelyRunner
namespace RankRigidity

/-- **The 2×2 rank-drop.**  If `a·u + b·v = λ·w` and `a'·u + b'·v = λ'·w` pointwise, with the
directions independent (`a·b' − a'·b ≠ 0`), then `u` and `v` are linearly dependent: some
`(p,q) ≠ (0,0)` has `p·u + q·v = 0` pointwise.  (Both `u` and `v` are forced proportional to
`w` by inverting the 2×2 system.) -/
theorem dep_of_two_proportional {n : ℕ} (u v w : Fin n → ℝ)
    (a b a' b' lam lam' : ℝ) (hdet : a * b' - a' * b ≠ 0)
    (h1 : ∀ i, a * u i + b * v i = lam * w i)
    (h2 : ∀ i, a' * u i + b' * v i = lam' * w i) :
    ∃ p q : ℝ, (p, q) ≠ (0, 0) ∧ ∀ i, p * u i + q * v i = 0 := by
  -- cross-multiply the 2×2 system (no division): det·u = (b'λ−bλ')·w, det·v = (aλ'−a'λ)·w
  have keyU : ∀ i, (a * b' - a' * b) * u i = (b' * lam - b * lam') * w i := by
    intro i; linear_combination b' * h1 i - b * h2 i
  have keyV : ∀ i, (a * b' - a' * b) * v i = (a * lam' - a' * lam) * w i := by
    intro i; linear_combination (-a') * h1 i + a * h2 i
  -- the dependency coefficients
  set p : ℝ := a * lam' - a' * lam with hp
  set q : ℝ := b * lam' - b' * lam with hq
  have hzero : ∀ i, p * u i + q * v i = 0 := by
    intro i
    have hdd : (a * b' - a' * b) * (p * u i + q * v i) = 0 := by
      rw [hp, hq]
      linear_combination (a * lam' - a' * lam) * keyU i + (b * lam' - b' * lam) * keyV i
    exact (mul_eq_zero.mp hdd).resolve_left hdet
  by_cases hpq0 : p = 0 ∧ q = 0
  · -- p = q = 0 forces u = 0 (det ≠ 0, since b'λ − bλ' = −q = 0); use (1, 0)
    refine ⟨1, 0, by simp, fun i => ?_⟩
    have hu0 : u i = 0 := by
      have hk := keyU i
      have hbb : b' * lam - b * lam' = 0 := by
        have hq0 : q = 0 := hpq0.2
        rw [hq] at hq0; linarith
      rw [hbb, zero_mul] at hk
      exact (mul_eq_zero.mp hk).resolve_left hdet
    rw [hu0]; ring
  · exact ⟨p, q, fun hc => hpq0 ⟨by simpa using congrArg Prod.fst hc,
      by simpa using congrArg Prod.snd hc⟩, hzero⟩

/-- **Rigidity, contrapositive-ready form.**  Two independent directions whose projections
are proportional to a common vector force `u, v` dependent — so a genuinely rank-2 pair
(no nonzero `(p,q)` with `p·u + q·v = 0`) CANNOT have two such common-proportional
directions.  This is what the pigeonhole feeds: if every full-support projection were a
dilated AP (all proportional, per shared ordering, to that ordering's vector), rank-2 fails. -/
theorem not_two_proportional_of_indep {n : ℕ} (u v : Fin n → ℝ)
    (hindep : ∀ p q : ℝ, (∀ i, p * u i + q * v i = 0) → p = 0 ∧ q = 0)
    (w : Fin n → ℝ) (a b a' b' lam lam' : ℝ) (hdet : a * b' - a' * b ≠ 0)
    (h1 : ∀ i, a * u i + b * v i = lam * w i)
    (h2 : ∀ i, a' * u i + b' * v i = lam' * w i) : False := by
  obtain ⟨p, q, hpq, hzero⟩ := dep_of_two_proportional u v w a b a' b' lam lam' hdet h1 h2
  exact hpq (Prod.ext (hindep p q hzero).1 (hindep p q hzero).2)

/-- **The pigeonhole wrapper.**  Consider the 1-parameter direction family `(1, N)`,
`N : ℤ`.  If every projection `u + N·v` is proportional to `W(cls N)` for a classifier
`cls : ℤ → S` into a FINITE type `S` (in the rigidity, `S` = signed permutations of
`{1..12}` and `W(cls N)` is the dilated-AP vector of the tight projection), then `u` and
`v` are linearly dependent.  Infinitely many `N`, finitely many classes ⟹ two share a
class ⟹ two independent directions proportional to the SAME vector ⟹ `dep_of_two_proportional`. -/
theorem dep_of_infinite_common_proportional {n : ℕ} (u v : Fin n → ℝ)
    {S : Type*} [Finite S] (cls : ℤ → S) (W : S → Fin n → ℝ) (lamf : ℤ → ℝ)
    (hprop : ∀ N : ℤ, ∀ i, u i + (N : ℝ) * v i = lamf N * W (cls N) i) :
    ∃ p q : ℝ, (p, q) ≠ (0, 0) ∧ ∀ i, p * u i + q * v i = 0 := by
  obtain ⟨N, N', hNN, hcls⟩ := Finite.exists_ne_map_eq_of_infinite cls
  refine dep_of_two_proportional u v (W (cls N)) 1 (N : ℝ) 1 (N' : ℝ)
    (lamf N) (lamf N') ?_ (fun i => ?_) (fun i => ?_)
  · -- det = 1·N' − 1·N = N' − N ≠ 0
    have : (N' : ℝ) ≠ (N : ℝ) := by exact_mod_cast (Ne.symm hNN)
    simp only [one_mul]
    intro h; apply this; linarith
  · have := hprop N i; simpa using this
  · have := hprop N' i
    rw [← hcls] at this
    simpa using this

#print axioms dep_of_two_proportional
#print axioms not_two_proportional_of_indep
#print axioms dep_of_infinite_common_proportional

end RankRigidity
end LonelyRunner
