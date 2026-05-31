/-
  TournamentH7.LonelyRunner — Lean formalization of the Lonely Runner sieve theory.

  The Lonely Runner Conjecture (reduced form): for `k` nonzero integer speeds
  `v : ι → ℤ` (`|ι| = k`) and threshold `1/n` with `n = k+1`, there is a time
  `t ∈ ℝ` with `‖v i · t‖ ≥ 1/n` for every `i`, where `‖x‖` is the distance from
  `x` to the nearest integer.

  This module formalizes the *denominator-sieve* layer of the conjecture, which
  this repo had previously only verified computationally (Python, S356–S362,
  S388) and catalogued as THM-358 / THM-360 / THM-366.  This Lean module is
  catalogued as THM-369 (the machine-checked formalization of THM-366).
  The master result `sieve_frac` proves that
  for any reduced fraction `a/q` with `q ≤ n`, if no speed is divisible by `q`
  then `a/q` is a lonely time.  It uniformly subsumes:
    · THM-360  (unit endpoint divisibility filter) — the `q = n` case;
    · the positive direction of THM-358 (initial-segment unit witnesses);
    · the "all speeds odd ⇒ t = 1/2 is lonely" antipodal tool for even `n`;
    · the CRT folding tools `t = 1/p` for `p ∣ n` used at composite `n`.
  Consequently a counterexample must contain a speed divisible by every
  `q ∈ {2,…,n}` (`counterexample_needs_all_divisors`).

  oracle-2026-05-31-S18.
-/
import Mathlib.Tactic

namespace LonelyRunner

/-- `t` is `n`-lonely for the speed family `v` iff every `v i · t` is at distance
`≥ 1/n` from every integer.  This is the standard `‖v i · t‖ ≥ 1/n` condition on
the circle `ℝ/ℤ`, written in the proof-friendly "far from every integer" form. -/
def Lonely {ι : Type*} (n : ℕ) (v : ι → ℤ) (t : ℝ) : Prop :=
  ∀ i, ∀ m : ℤ, (1 : ℝ) / n ≤ |(v i : ℝ) * t - m|

section Sieve
variable {ι : Type*}

/-- **Master sieve lemma.**  Let `q ≤ n` with `0 < q`, let `a` be coprime to `q`,
and suppose no speed `v i` is divisible by `q`.  Then the reduced fraction `a/q`
is an `n`-lonely time. -/
theorem sieve_frac (n q : ℕ) (a : ℤ) (v : ι → ℤ)
    (hqn : q ≤ n) (hq0 : 0 < q) (hcop : IsCoprime (q : ℤ) a)
    (hdiv : ∀ i, ¬ ((q : ℤ) ∣ v i)) :
    Lonely n v ((a : ℝ) / q) := by
  intro i m
  have hq0' : (0 : ℝ) < (q : ℝ) := by exact_mod_cast hq0
  have hqne : (q : ℝ) ≠ 0 := ne_of_gt hq0'
  -- Put the value over the common denominator `q`.
  have key : (v i : ℝ) * ((a : ℝ) / q) - m
      = (((v i * a - m * q : ℤ)) : ℝ) / q := by
    field_simp
    push_cast
    ring
  rw [key, abs_div, abs_of_pos hq0']
  -- The numerator is a nonzero integer: `q ∤ v i · a`, since `q ∤ v i` and `q ⟂ a`.
  have hnotdvd : ¬ ((q : ℤ) ∣ v i * a) := fun h => hdiv i (hcop.dvd_of_dvd_mul_right h)
  have hne : (v i * a - m * q) ≠ 0 := by
    intro h
    apply hnotdvd
    have hva : v i * a = m * q := by linarith [sub_eq_zero.mp h]
    exact ⟨m, by rw [hva]; ring⟩
  -- A nonzero integer has absolute value `≥ 1`.
  have h1 : (1 : ℝ) ≤ |((v i * a - m * q : ℤ) : ℝ)| := by
    rw [← Int.cast_abs]
    have : (1 : ℤ) ≤ |v i * a - m * q| := Int.one_le_abs hne
    exact_mod_cast this
  -- Chain `1/n ≤ 1/q ≤ |num|/q`.
  have hnq : (1 : ℝ) / n ≤ (1 : ℝ) / q :=
    one_div_le_one_div_of_le hq0' (by exact_mod_cast hqn)
  calc (1 : ℝ) / n ≤ (1 : ℝ) / q := hnq
    _ ≤ |((v i * a - m * q : ℤ) : ℝ)| / q := by gcongr

/-- The `a = 1` specialization: if no speed is divisible by `q` (`0 < q ≤ n`) then
`t = 1/q` is lonely.  This is the "denominator sieve". -/
theorem sieve_one_div (n q : ℕ) (v : ι → ℤ)
    (hqn : q ≤ n) (hq0 : 0 < q) (hdiv : ∀ i, ¬ ((q : ℤ) ∣ v i)) :
    Lonely n v ((1 : ℝ) / q) := by
  have h := sieve_frac n q 1 v hqn hq0 isCoprime_one_right hdiv
  simpa using h

/-- **Necessary condition for a counterexample.**  If `v` has *no* `n`-lonely
time, then for every `q` with `2 ≤ q ≤ n` some speed is divisible by `q`.  In
particular a counterexample must contain a multiple of `n` and of every prime
power `≤ n`. -/
theorem counterexample_needs_all_divisors (n : ℕ) (v : ι → ℤ)
    (hcex : ∀ t : ℝ, ¬ Lonely n v t)
    (q : ℕ) (hq2 : 2 ≤ q) (hqn : q ≤ n) :
    ∃ i, (q : ℤ) ∣ v i := by
  by_contra h
  push Not at h
  exact hcex ((1 : ℝ) / q) (sieve_one_div n q v hqn (by omega) h)

/-- **Antipodal tool (even-`n`).**  If every speed is odd and `2 ≤ n`, then
`t = 1/2` is lonely.  (Special case `q = 2`.) -/
theorem all_odd_half_lonely (n : ℕ) (v : ι → ℤ)
    (hn : 2 ≤ n) (hodd : ∀ i, ¬ (2 : ℤ) ∣ v i) :
    Lonely n v ((1 : ℝ) / 2) := by
  have h := sieve_one_div n 2 v hn (by norm_num) hodd
  simpa using h

end Sieve

/-- **Positive direction of THM-358.**  For the initial segment `v i = i+1`
(`i : Fin (n-1)`, i.e. speeds `1,…,n-1`) and any `a` coprime to `n`, the unit
fraction `a/n` is an `n`-lonely time.  Hence `{1,…,n-1}` is at best tight: it
always has the unit witnesses `a/n`. -/
theorem initial_segment_unit_lonely (n : ℕ) (hn : 1 ≤ n) (a : ℤ)
    (hcop : IsCoprime (n : ℤ) a) :
    Lonely n (fun i : Fin (n - 1) => (i : ℤ) + 1) ((a : ℝ) / n) := by
  apply sieve_frac n n a _ (le_refl n) (by omega) hcop
  intro i hdvd
  -- `n ∤ (i+1)` because `1 ≤ i+1 ≤ n-1 < n`.
  have hpos : (0 : ℤ) < (i : ℤ) + 1 := by positivity
  have hle : (n : ℤ) ≤ (i : ℤ) + 1 := Int.le_of_dvd hpos hdvd
  have hb : (i : ℕ) < n - 1 := i.isLt
  omega

/-! ### Axiom audit

Each `#print axioms` should report only the foundational
`[propext, Classical.choice, Quot.sound]` — i.e. no `sorry` and no
project-specific axioms underlie the LRC sieve theory. -/

#print axioms sieve_frac
#print axioms sieve_one_div
#print axioms counterexample_needs_all_divisors
#print axioms all_odd_half_lonely
#print axioms initial_segment_unit_lonely

end LonelyRunner
