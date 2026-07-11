/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: kind-pasteur (LRC multi-agent project, 2026-07-11-S127 cont.34)
-/
import TournamentH7.LRCDiscreteBonferroni

/-!
# `B5` is residue-periodic — the bounded-window ruler search is DIAMETER-FREE

The single LRC(14) obligation `hB5` asks each residual family for a modulus `q` with `B5 v q > 0`. A natural
strengthening is: does a **fixed, diameter-independent window** `q ∈ [8, Q]` always contain such a `q`?
(Empirically yes, with `Q ≈ 43`, even at diameter `~10⁵`.) The reason it can be diameter-free is this file's
lemma:

> **`B5_congr_mod`**: `B5 v q` depends only on the residues `v i mod q`. So for `q ≤ Q`, whether any
> `q ∈ [8, Q]` gives `B5 v q > 0` depends only on `(v i mod lcm(8..Q))` — a finite property, independent of
> the actual sizes of the speeds.

`bandCount v q p = #{i : (vᵢ·p mod q) ∉ [q/14, 13q/14]}` reads only `vᵢ·p mod q`, which is fixed by
`vᵢ mod q`; the whole `B5` alternating sum inherits the invariance.
-/

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- `bandCount` reads only the residues of `vᵢ·p mod q`. -/
lemma bandCount_congr_mod (v w : Fin 13 → ℤ) (q p : ℕ)
    (h : ∀ i, (v i * (p : ℤ)) % (q : ℤ) = (w i * (p : ℤ)) % (q : ℤ)) :
    bandCount v q p = bandCount w q p := by
  unfold bandCount
  apply congrArg Finset.card
  apply Finset.filter_congr
  intro i _
  have hi : inBand v q p i ↔ inBand w q p i := by
    unfold inBand; rw [h i]
  exact not_congr hi

/-- **`B5` is residue-periodic.**  If `vᵢ ≡ wᵢ (mod q)` for every `i`, then `B5 v q = B5 w q`.  Hence the
bounded-window ruler search `∃ q ∈ [8, Q], B5 v q > 0` depends only on `v mod lcm(8..Q)` — diameter-free. -/
theorem B5_congr_mod (v w : Fin 13 → ℤ) (q : ℕ)
    (h : ∀ i, v i ≡ w i [ZMOD (q : ℤ)]) :
    B5 v q = B5 w q := by
  have hband : ∀ p, bandCount v q p = bandCount w q p := by
    intro p
    apply bandCount_congr_mod
    intro i
    exact (h i).mul_right (p : ℤ)
  unfold B5 momentS
  simp only [hband]

/-- The bounded-window predicate: some modulus in `[8, Q]` is a live `B5`-ruler. -/
def HasWindowRuler (v : Fin 13 → ℤ) (Q : ℕ) : Prop :=
  ∃ q ∈ Finset.Icc 8 Q, 0 < B5 v q

/-- **The window predicate is residue-periodic (diameter-free).**  If `vᵢ ≡ wᵢ (mod lcm)` for a modulus
`lcm` divisible by every `q ∈ [8, Q]`, then `v` and `w` have the same window-ruler status.  So
"`HasWindowRuler _ Q`" is a property of the residue class `v mod lcm`, not of the diameter. -/
theorem hasWindowRuler_congr (v w : Fin 13 → ℤ) (Q lcm : ℕ)
    (hdvd : ∀ q ∈ Finset.Icc 8 Q, q ∣ lcm)
    (h : ∀ i, v i ≡ w i [ZMOD (lcm : ℤ)]) :
    HasWindowRuler v Q ↔ HasWindowRuler w Q := by
  have key : ∀ q ∈ Finset.Icc 8 Q, B5 v q = B5 w q := by
    intro q hq
    apply B5_congr_mod
    intro i
    have hqd : (q : ℤ) ∣ (lcm : ℤ) := Int.natCast_dvd_natCast.mpr (hdvd q hq)
    exact (h i).of_dvd hqd
  constructor
  · rintro ⟨q, hq, hpos⟩; exact ⟨q, hq, (key q hq) ▸ hpos⟩
  · rintro ⟨q, hq, hpos⟩; exact ⟨q, hq, (key q hq).symm ▸ hpos⟩

end LRC14Concrete
end LonelyRunner

#print axioms LonelyRunner.LRC14Concrete.B5_congr_mod
#print axioms LonelyRunner.LRC14Concrete.hasWindowRuler_congr
