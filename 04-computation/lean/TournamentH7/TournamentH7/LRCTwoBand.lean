/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-06-S97)
-/
import TournamentH7.LRCClearCert

/-!
# The two-band transport: the 1d-crux mechanism, exact (HYP-4256)

A two-band family CORE ∪ S·P is handled by transporting the PATTERN's citation witness
into the core's clear interval: with `t = (t_P + k)/S` for the right integer `k`
(existing as soon as `S·L > 1`), every top runner satisfies
`(S·p)·t = p·t_P + p·k ≡ p·t_P (mod 1)` — the pattern's margins transport EXACTLY,
with no Lipschitz loss, no measure, no free-fraction analysis. The core's margins come
from the clear interval (certified by `toothMiss` tables / `clear_of_cert` at its band).

This is height-uniformity on CRT-frozen rays in four moves: cite the core (→ J), cite
the pattern (→ t_P, margin ≥ 1/13 since |P| ≤ 12), pick k = ⌊S·a − t_P⌋ + 1, transport.
Verified end-to-end at scales to 10^12+7 (two_band_exact_opus_S97.out): top margins are
EXACTLY the pattern's, at every scale.
-/

namespace LonelyRunner
namespace TwoBand

/-- **The two-band transport.** A clear core interval `(a, a+L)` plus a pattern witness
`t_P` yields, for every scale `S` with `S·L > 1`, one time serving both bands: the core
strictly above its band, the scaled pattern at EXACTLY the pattern's margins. -/
theorem two_band_transport
    (core : List ℝ) (P : List ℤ) (S : ℤ) (hS : 0 < S)
    (a L : ℝ) (hcore : ∀ t : ℝ, a < t → t < a + L → ∀ w ∈ core, ∀ m : ℤ,
      (1 : ℝ) / 14 < |w * t - m|)
    (tP : ℝ) (htP : ∀ p ∈ P, ∀ m : ℤ, (1 : ℝ) / 13 ≤ |(p : ℝ) * tP - m|)
    (hSL : 1 < (S : ℝ) * L) :
    ∃ t : ℝ, (∀ w ∈ core, ∀ m : ℤ, (1 : ℝ) / 14 < |w * t - m|) ∧
             (∀ p ∈ P, ∀ m : ℤ, (1 : ℝ) / 13 ≤ |((S * p : ℤ) : ℝ) * t - m|) := by
  have hSR : (0 : ℝ) < (S : ℝ) := by exact_mod_cast hS
  set k : ℤ := ⌊(S : ℝ) * a - tP⌋ + 1 with hk
  set t : ℝ := (tP + k) / S with ht
  have hk1 : (S : ℝ) * a - tP < (k : ℝ) := by
    rw [hk]
    push_cast
    exact Int.lt_floor_add_one _
  have hk2 : (k : ℝ) ≤ (S : ℝ) * a - tP + 1 := by
    rw [hk]
    push_cast
    have := Int.floor_le ((S : ℝ) * a - tP)
    linarith
  have hta : a < t := by
    rw [ht, lt_div_iff₀ hSR]
    have hcomm : a * S = (S : ℝ) * a := mul_comm a S
    linarith
  have htb : t < a + L := by
    rw [ht, div_lt_iff₀ hSR]
    have hexp : (a + L) * S = (S : ℝ) * a + (S : ℝ) * L := by ring
    linarith
  refine ⟨t, fun w hw m => hcore t hta htb w hw m, fun p hp m => ?_⟩
  have key : ((S * p : ℤ) : ℝ) * t = (p : ℝ) * tP + ((p * k : ℤ) : ℝ) := by
    rw [ht]
    push_cast
    field_simp
  rw [key]
  have h2 := htP p hp (m - p * k)
  convert h2 using 2
  push_cast
  ring

/-- The `Lonely 14`-grade corollary: everything sits strictly above `1/14`. -/
theorem two_band_lonely14
    (core : List ℝ) (P : List ℤ) (S : ℤ) (hS : 0 < S)
    (a L : ℝ) (hcore : ∀ t : ℝ, a < t → t < a + L → ∀ w ∈ core, ∀ m : ℤ,
      (1 : ℝ) / 14 < |w * t - m|)
    (tP : ℝ) (htP : ∀ p ∈ P, ∀ m : ℤ, (1 : ℝ) / 13 ≤ |(p : ℝ) * tP - m|)
    (hSL : 1 < (S : ℝ) * L) :
    ∃ t : ℝ, (∀ w ∈ core, ∀ m : ℤ, (1 : ℝ) / 14 < |w * t - m|) ∧
             (∀ p ∈ P, ∀ m : ℤ, (1 : ℝ) / 14 < |((S * p : ℤ) : ℝ) * t - m|) := by
  obtain ⟨t, h1, h2⟩ := two_band_transport core P S hS a L hcore tP htP hSL
  exact ⟨t, h1, fun p hp m => lt_of_lt_of_le (by norm_num) (h2 p hp m)⟩

#print axioms two_band_transport
#print axioms two_band_lonely14

end TwoBand
end LonelyRunner
