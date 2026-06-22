/-
  TournamentH7.LRCBindingPair -- the binding-pair arithmetic core of HYP-2909
  (mac-mini S50, the forward direction of the LRC(14) tightness crux ★).

  At a primitive denominator-14 optimum `t = a/14` (a coprime to 14), tightness
  `M(S) = 1/14` forces an active runner `si` on the INCREASING side (at +1/14, i.e.
  `si·a ≡ 1 mod 14`) and an active runner `sj` on the DECREASING side (at -1/14, i.e.
  `sj·a ≡ -1 mod 14`).  Then `(si+sj)·a ≡ 0 mod 14`, and cancelling the unit `a`
  gives the **binding pair**

      14 ∣ (si + sj)        (the apex-7 antipodal residue pair, 14 = 2·7).

  This is the LRC analogue of THM-079's Move B (cycle-forcing): tightness sits at
  the apex-7 point.  Together with `LRCApex7Floor.D14_never_certifies` (a multiple
  of 14 sits on the observer at every `a/14`), a covering set cannot bind at a
  denom-14 optimum.  Sorry-free, elementary `Int.ModEq` / `IsCoprime`.
-/

import Mathlib

namespace LonelyRunner.BindingPair

/-- **Binding-pair theorem (arithmetic core of HYP-2909).** At a primitive denom-14
optimum `t = a/14` (`a` coprime to 14): if `si` is on the increasing side
(`si·a ≡ 1 [ZMOD 14]`) and `sj` on the decreasing side (`sj·a ≡ -1 [ZMOD 14]`), then
`14 ∣ (si + sj)` -- the two active runners are antipodal residues mod 14. -/
theorem binding_pair_dvd (a si sj : ℤ) (ha : IsCoprime (14 : ℤ) a)
    (hi : si * a ≡ 1 [ZMOD 14]) (hj : sj * a ≡ -1 [ZMOD 14]) :
    (14 : ℤ) ∣ (si + sj) := by
  have hsum : (si + sj) * a ≡ 0 [ZMOD 14] := by
    have h : si * a + sj * a ≡ (1 : ℤ) + (-1) [ZMOD 14] := hi.add hj
    have h0 : si * a + sj * a ≡ 0 [ZMOD 14] := by simpa using h
    have he : (si + sj) * a = si * a + sj * a := by ring
    rw [he]; exact h0
  have hdvd : (14 : ℤ) ∣ (si + sj) * a := (Int.modEq_zero_iff_dvd).mp hsum
  exact ha.dvd_of_dvd_mul_right hdvd

/-- The antipodal-residue form: the increasing runner's residue and the decreasing
runner's residue are negatives mod 14 (`si·a + sj·a ≡ 0`), i.e. the binding pair lies
in one of the 7 antipodal residue classes `{1,13},{2,12},…,{6,8},{7,7}` — the apex-7. -/
theorem binding_residues_antipodal (a si sj : ℤ)
    (hi : si * a ≡ 1 [ZMOD 14]) (hj : sj * a ≡ -1 [ZMOD 14]) :
    si * a + sj * a ≡ 0 [ZMOD 14] := by
  simpa using hi.add hj

/-! ## Axiom audit -/

#print axioms binding_pair_dvd
#print axioms binding_residues_antipodal

end LonelyRunner.BindingPair
