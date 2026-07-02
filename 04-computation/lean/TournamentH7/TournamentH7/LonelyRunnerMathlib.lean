/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: kind-pasteur (LRC multi-agent project, 2026-07-01-S29)
-/
import Mathlib.Analysis.Normed.Group.AddCircle
import Mathlib.Data.ZMod.Basic
import Mathlib.NumberTheory.DiophantineApproximation.Basic

/-!
# The Lonely Runner Conjecture: definitions and elementary bounds

This file is the MATHLIB-TRACK core of the project's LRC(14) formalization: self-contained
definitions and fully proved elementary lemmas, written to mathlib conventions so the file can be
extracted as `Mathlib/NumberTheory/LonelyRunner/Basic.lean` (mathlib has no lonely-runner content
as of v4.30). Project-specific certificates stay in the other `LRC*.lean` files.

We model the track as the unit circle `UnitAddCircle = ℝ ⧸ ℤ`, with the runner at speed `v`
located at time `t` at `(v * t : UnitAddCircle)`; its distance to the stationary runner at the
origin is the quotient norm `‖(v * t : ℝ)‖ = |v t - round (v t)|`.

## Main definitions

* `LonelyRunner.IsLonelyAt S r t` : at time `t` every speed in `S` is at distance `≥ r` from 0.
* `LonelyRunner.Conjecture k` : the Lonely Runner Conjecture for `k` moving runners
  (loneliness `1/(k+1)` is achieved for every set of `k` distinct positive integer speeds).

## Main results

* `LonelyRunner.isLonelyAt_of_forall_not_dvd` : **the q-witness lemma** — if no speed in `S` is
  divisible by `q` then time `1/q` is `1/q`-lonely.  (The easy half of the project's THM-523;
  it disposes of every "non-covering" speed set at once.)
* `LonelyRunner.isLonelyAt_image_mul` : dilation invariance of loneliness (project THM-522).
* `LonelyRunner.conjecture_one`, `LonelyRunner.conjecture_two` : the cases `k = 1, 2`.
* `LonelyRunner.norm_le_of_abs_le` : the origin window `|t| ≤ r/v` is `r`-dangerous for speed `v`
  (the containment half of the project's exact pair-overlap law THM-594).

## References

* J. Wills, *Zwei Sätze über inhomogene diophantische Approximation von Irrationalzahlen* (1967).
* T. Cusick, *View-obstruction problems* (1973).
* Sungkawichai–Trakulthongchai, arXiv:2604.23906 (LRC for `k ≤ 12`).
-/

namespace LonelyRunner

open Finset

/-- `IsLonelyAt S r t` : at time `t`, every runner whose speed lies in `S` is at circle-distance
at least `r` from the stationary runner at the origin. -/
def IsLonelyAt (S : Finset ℕ) (r t : ℝ) : Prop :=
  ∀ v ∈ S, r ≤ ‖(((v : ℝ) * t : ℝ) : UnitAddCircle)‖

/-- The **Lonely Runner Conjecture** for `k` moving runners: for every set of `k` distinct
positive integer speeds there is a time at which all runners are at distance at least
`1/(k+1)` from the origin. -/
def Conjecture (k : ℕ) : Prop :=
  ∀ S : Finset ℕ, S.card = k → (∀ v ∈ S, 0 < v) → ∃ t : ℝ, IsLonelyAt S (1 / (k + 1)) t

/-- The circle-distance of `m/n` to the integers, for natural `m, n`. -/
theorem norm_natCast_div (m n : ℕ) :
    ‖(((m : ℝ) / n : ℝ) : UnitAddCircle)‖ = (min (m % n) (n - m % n) : ℕ) / (n : ℝ) := by
  rw [UnitAddCircle.norm_eq, abs_sub_round_div_natCast_eq]

/-- If `q` does not divide `v` then `v/q` is at circle-distance at least `1/q` from the
integers. -/
theorem one_div_le_norm_div {v q : ℕ} (hq : 0 < q) (h : ¬ q ∣ v) :
    (1 : ℝ) / q ≤ ‖(((v : ℝ) / q : ℝ) : UnitAddCircle)‖ := by
  rw [norm_natCast_div]
  have hmod : v % q ≠ 0 := fun h0 => h (Nat.dvd_of_mod_eq_zero h0)
  have hlt : v % q < q := Nat.mod_lt _ hq
  have h1 : 1 ≤ min (v % q) (q - v % q) := by omega
  have hq' : (0 : ℝ) < q := by exact_mod_cast hq
  gcongr
  exact_mod_cast h1

/-- **The q-witness lemma**: if no speed in `S` is divisible by `q`, then the time `1/q` is
`1/q`-lonely.  In the project's language this closes every non-covering configuration; a
counterexample to the Lonely Runner Conjecture must contain a multiple of every
`q ≤ k + 1`. -/
theorem isLonelyAt_of_forall_not_dvd (S : Finset ℕ) {q : ℕ} (hq : 0 < q)
    (h : ∀ v ∈ S, ¬ q ∣ v) : IsLonelyAt S (1 / q) (1 / q) := by
  intro v hv
  rw [mul_one_div]
  exact one_div_le_norm_div hq (h v hv)

/-- **The covering reduction**: a set of `k` speeds none of which is divisible by `k + 1`
satisfies the Lonely Runner bound.  Hence any counterexample to `Conjecture k` must contain a
multiple of `k + 1` (and, iterating the q-witness over `q ≤ k + 1`, a multiple of every such
`q` — a "covering" configuration). -/
theorem exists_isLonelyAt_of_forall_not_dvd (S : Finset ℕ) (k : ℕ)
    (h : ∀ v ∈ S, ¬ (k + 1) ∣ v) :
    ∃ t : ℝ, IsLonelyAt S (1 / ((k : ℝ) + 1)) t := by
  have hL := isLonelyAt_of_forall_not_dvd S (Nat.succ_pos k) h
  have hcast : ((k + 1 : ℕ) : ℝ) = (k : ℝ) + 1 := by push_cast; ring
  rw [hcast] at hL
  exact ⟨_, hL⟩

/-- Dilation invariance (project THM-522): multiplying every speed by `c > 0` and dividing time
by `c` preserves loneliness.  Stated as an equivalence between the dilated set at time `t` and
the original set at time `c * t`. -/
theorem isLonelyAt_image_mul (c : ℕ) (S : Finset ℕ) (r t : ℝ) :
    IsLonelyAt (S.image (c * ·)) r t ↔ IsLonelyAt S r ((c : ℝ) * t) := by
  constructor
  · intro h v hv
    have := h (c * v) (mem_image_of_mem _ hv)
    have hcast : ((c * v : ℕ) : ℝ) * t = (v : ℝ) * ((c : ℝ) * t) := by push_cast; ring
    rwa [hcast] at this
  · intro h v hv
    obtain ⟨w, hw, rfl⟩ := mem_image.mp hv
    have := h w hw
    have hcast : ((c * w : ℕ) : ℝ) * t = (w : ℝ) * ((c : ℝ) * t) := by push_cast; ring
    rwa [← hcast] at this

/-- The origin window is dangerous (the containment half of the exact pair-overlap law,
project THM-594): if `|t| ≤ r / v` then the runner at speed `v` is within `r` of the origin. -/
theorem norm_le_of_abs_le {v : ℕ} {r t : ℝ} (hv : 0 < v) (h : |t| ≤ r / v) :
    ‖(((v : ℝ) * t : ℝ) : UnitAddCircle)‖ ≤ r := by
  have hv' : (0 : ℝ) < v := by exact_mod_cast hv
  have h0 : ‖(((v : ℝ) * t : ℝ) : UnitAddCircle)‖ ≤ |(v : ℝ) * t| := by
    rw [UnitAddCircle.norm_eq]
    simpa using round_le ((v : ℝ) * t) 0
  calc ‖(((v : ℝ) * t : ℝ) : UnitAddCircle)‖ ≤ |(v : ℝ) * t| := h0
    _ = (v : ℝ) * |t| := by rw [abs_mul, abs_of_pos hv']
    _ ≤ (v : ℝ) * (r / v) := by gcongr
    _ = r := by field_simp

/-- The Lonely Runner Conjecture holds for one runner: at time `1/(2v)` the runner is exactly
opposite the origin. -/
theorem conjecture_one : Conjecture 1 := by
  intro S hcard hpos
  obtain ⟨v, rfl⟩ := card_eq_one.mp hcard
  have hv : 0 < v := hpos v (mem_singleton_self v)
  have hv' : (0 : ℝ) < v := by exact_mod_cast hv
  refine ⟨1 / (2 * v), ?_⟩
  intro w hw
  rw [mem_singleton] at hw
  subst hw
  have hval : (w : ℝ) * (1 / (2 * w)) = 1 / 2 := by field_simp
  rw [hval, UnitAddCircle.norm_eq]
  have : round ((1 : ℝ) / 2) = 1 := by
    rw [round_eq]; norm_num
  rw [this]
  norm_num

section Tightness

/-- **Tightness of the Lonely Runner bound** (Dirichlet's approximation theorem): at every time
`t`, some speed in `{1, …, k}` is within `1/(k+1)` of the origin.  So the constant `1/(k+1)` in
the conjecture cannot be improved. -/
theorem exists_norm_le_of_mem_Icc (k : ℕ) (hk : 0 < k) (t : ℝ) :
    ∃ v ∈ Finset.Icc 1 k, ‖(((v : ℝ) * t : ℝ) : UnitAddCircle)‖ ≤ 1 / ((k : ℝ) + 1) := by
  obtain ⟨j, s, hs0, hsk, h⟩ := Real.exists_int_int_abs_mul_sub_le t hk
  refine ⟨s.toNat, ?_, ?_⟩
  · rw [Finset.mem_Icc]
    omega
  · have hvs : ((s.toNat : ℕ) : ℝ) = (s : ℝ) := by
      exact_mod_cast congrArg (Int.cast : ℤ → ℝ) (Int.toNat_of_nonneg hs0.le)
    rw [hvs]
    calc ‖(((s : ℝ) * t : ℝ) : UnitAddCircle)‖
        = |(s : ℝ) * t - round ((s : ℝ) * t)| := UnitAddCircle.norm_eq
      _ ≤ |(s : ℝ) * t - j| := round_le _ j
      _ ≤ 1 / ((k : ℝ) + 1) := h

/-- The constant `1/(k+1)` in the Lonely Runner Conjecture is optimal: for the consecutive
speeds `{1, …, k}` no time is `r`-lonely for any `r > 1/(k+1)`. -/
theorem not_isLonelyAt_Icc_of_lt {k : ℕ} (hk : 0 < k) {r t : ℝ}
    (hr : 1 / ((k : ℝ) + 1) < r) : ¬ IsLonelyAt (Finset.Icc 1 k) r t := by
  intro hL
  obtain ⟨v, hv, hnorm⟩ := exists_norm_le_of_mem_Icc k hk t
  have := hL v hv
  linarith

end Tightness

section TwoRunners

/-- The coprime two-runner core: for coprime positive `a ≠ b` there is a rational time
`j/(a+b)` at which both runners are at distance at least `1/3` from the origin.  The two
distances agree because `a ≡ -b (mod a+b)`; choosing `j ≡ a⁻¹ ⌊(a+b)/2⌋ (mod a+b)` places both
runners in the middle third of the circle. -/
theorem exists_isLonelyAt_pair_coprime {a b : ℕ} (ha : 0 < a) (hb : 0 < b) (hab : a ≠ b)
    (hco : Nat.Coprime a b) :
    ∃ t : ℝ, IsLonelyAt {a, b} (1 / 3) t := by
  set n : ℕ := a + b with hn
  have hn3 : 3 ≤ n := by omega
  have hnpos : 0 < n := by omega
  haveI : NeZero n := ⟨by omega⟩
  set m : ℕ := n / 2 with hm
  have hm1 : 1 ≤ m := by omega
  have hmn : m < n := by omega
  have h3m : n ≤ 3 * m := by omega
  have h2m : 2 * m ≤ n := by omega
  -- `a` is a unit mod `n = a + b`, since `gcd (a, a + b) = gcd (a, b) = 1`
  have hcoan : Nat.Coprime a n := by
    have h1 : Nat.gcd a (b + a) = Nat.gcd a b := Nat.gcd_add_self_right a b
    unfold Nat.Coprime at hco ⊢
    rw [hn, Nat.add_comm a b, h1]
    exact hco
  -- the witness numerator
  set x : ZMod n := (a : ZMod n)⁻¹ * (m : ZMod n) with hx
  set j : ℕ := x.val with hj
  have hxj : ((j : ℕ) : ZMod n) = x := by rw [hj, ZMod.natCast_val, ZMod.cast_id]
  -- a * j ≡ m (mod n)
  have hax : (a : ZMod n) * x = (m : ZMod n) := by
    rw [hx, ← mul_assoc, ZMod.coe_mul_inv_eq_one a hcoan, one_mul]
  have hajm : (a * j) % n = m := by
    have hcast : ((a * j : ℕ) : ZMod n) = ((m : ℕ) : ZMod n) := by
      push_cast [hxj]
      exact hax
    have := (ZMod.natCast_eq_natCast_iff _ _ _).mp hcast
    unfold Nat.ModEq at this
    rwa [Nat.mod_eq_of_lt hmn] at this
  -- b * j ≡ n - m (mod n), because b ≡ -a
  have hbz : (b : ZMod n) = -(a : ZMod n) := by
    have h0 : ((a + b : ℕ) : ZMod n) = 0 := by
      rw [← hn]; exact_mod_cast ZMod.natCast_self n
    push_cast at h0
    exact eq_neg_of_add_eq_zero_right h0
  have hbjm : (b * j) % n = n - m := by
    have hnm : (((n - m : ℕ)) : ZMod n) = -(m : ZMod n) := by
      rw [Nat.cast_sub hmn.le, ZMod.natCast_self, zero_sub]
    have hcast : ((b * j : ℕ) : ZMod n) = (((n - m : ℕ)) : ZMod n) := by
      push_cast [hxj]
      rw [hbz, hnm, neg_mul, hax]
    have := (ZMod.natCast_eq_natCast_iff _ _ _).mp hcast
    unfold Nat.ModEq at this
    rwa [Nat.mod_eq_of_lt (by omega : n - m < n)] at this
  -- the witness time
  refine ⟨(j : ℝ) / n, ?_⟩
  have hn' : (0 : ℝ) < n := by exact_mod_cast hnpos
  intro v hv
  rcases Finset.mem_insert.mp hv with rfl | hv'
  · -- v = a
    have hcast : (v : ℝ) * ((j : ℝ) / n) = ((v * j : ℕ) : ℝ) / n := by push_cast; ring
    rw [hcast, norm_natCast_div, hajm]
    have hminm : min m (n - m) = m := by omega
    rw [hminm, div_le_div_iff₀ (by norm_num) hn']
    exact_mod_cast (by omega : 1 * n ≤ m * 3)
  · -- v = b
    rw [Finset.mem_singleton] at hv'
    subst hv'
    have hcast : (v : ℝ) * ((j : ℝ) / n) = ((v * j : ℕ) : ℝ) / n := by push_cast; ring
    rw [hcast, norm_natCast_div, hbjm]
    have hminm : min (n - m) (n - (n - m)) = m := by omega
    rw [hminm, div_le_div_iff₀ (by norm_num) hn']
    exact_mod_cast (by omega : 1 * n ≤ m * 3)

/-- The Lonely Runner Conjecture holds for two runners: reduce to the coprime case by dividing
out the gcd (dilation invariance), then apply the middle-third construction. -/
theorem conjecture_two : Conjecture 2 := by
  intro S hcard hpos
  obtain ⟨a, b, hab, rfl⟩ := card_eq_two.mp hcard
  have ha : 0 < a := hpos a (by simp)
  have hb : 0 < b := hpos b (by simp)
  set g : ℕ := Nat.gcd a b with hg
  have hgpos : 0 < g := Nat.gcd_pos_of_pos_left b ha
  have hg' : (g : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hgpos.ne'
  set a' : ℕ := a / g with ha'def
  set b' : ℕ := b / g with hb'def
  have hag : a = g * a' := (Nat.mul_div_cancel' (Nat.gcd_dvd_left a b)).symm
  have hbg : b = g * b' := (Nat.mul_div_cancel' (Nat.gcd_dvd_right a b)).symm
  have ha'pos : 0 < a' := Nat.div_pos (Nat.le_of_dvd ha (Nat.gcd_dvd_left a b)) hgpos
  have hb'pos : 0 < b' := Nat.div_pos (Nat.le_of_dvd hb (Nat.gcd_dvd_right a b)) hgpos
  have hab' : a' ≠ b' := by
    intro h
    exact hab (by rw [hag, hbg, h])
  have hco : Nat.Coprime a' b' := Nat.coprime_div_gcd_div_gcd hgpos
  obtain ⟨t, ht⟩ := exists_isLonelyAt_pair_coprime ha'pos hb'pos hab' hco
  refine ⟨t / g, ?_⟩
  have hr : (1 : ℝ) / ((2 : ℕ) + 1) = 1 / 3 := by norm_num
  intro v hv
  rw [hr]
  rcases Finset.mem_insert.mp hv with rfl | hv'
  · have hcast : (v : ℝ) * (t / g) = (a' : ℝ) * t := by
      rw [hag]; push_cast; field_simp
    rw [hcast]
    exact ht a' (Finset.mem_insert_self _ _)
  · rw [Finset.mem_singleton] at hv'
    subst hv'
    have hcast : (v : ℝ) * (t / g) = (b' : ℝ) * t := by
      rw [hbg]; push_cast; field_simp
    rw [hcast]
    exact ht b' (by simp)

end TwoRunners

end LonelyRunner
