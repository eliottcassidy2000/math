/-
Copyright (c) 2026 The multi-agent math project contributors.
Released under Apache 2.0.
Authors: klein-S428

# Comparator challenge: GMC(2)

A **self-contained, independently auditable statement** of the Gaussian-Moments /
Mathieu-Zhao conjecture in dimension two. This file imports only Mathlib. It is
the statement half of a comparator pair, in the style of
`github.com/openai/ten-proofs/ComparatorChallenges`: the solution module must
supply a theorem of the *same name and type*, proved, using only the permitted
axioms `propext`, `Classical.choice`, `Quot.sound`.

## Why this file exists

`TournamentH7.GMC2Main.GMC2.gmc2` is already `sorry`-free and kernel-pure. But a
kernel-pure proof only certifies that the stated proposition follows; it says
nothing about whether the *statement* is the theorem anyone cares about. That
residual risk is real and has bitten this repo before (see MISTAKE-190, and the
THM-3018 scope error where `FC(n)` was stated with a homogeneity hypothesis the
actual conjecture does not have).

The risk here concentrates in one place: **is `GMC2.E` the Gaussian expectation?**
`GMC2.E` is defined as a sum over a polynomial's support against a weight `wt`.
Auditing that definition means auditing repo code.

So this challenge does not use `GMC2.E` at all. It quantifies over **every**
functional satisfying the Wick moment law. Checking the statement therefore
reduces to checking one textbook identity for a standard complex Gaussian `Z`,

```text
    E[ Z^a * conj(Z)^b ] = a! * [a = b],
```

written here with `Z = X 0` and `conj Z = X 1`. Anyone who accepts that formula
accepts this statement, without reading a line of the proof or of the repo's
definitions.

Because a functional satisfying `WickFunctional` is *unique* (additivity plus
`C`-homogeneity pin it down on monomials, and monomials span), this is equivalent
to the fixed-`E` statement -- it is not a stronger or weaker theorem, only an
independently auditable one.
-/
import Mathlib

noncomputable section

open MvPolynomial

namespace GMC2Challenge

/-- `E` is a Wick / Gaussian expectation on `ℂ[Z, W]`, where `Z = X 0` plays the
role of a standard complex Gaussian and `W = X 1` the role of its conjugate.

The three clauses are: additivity, homogeneity over scalars, and the Gaussian
moment law `E[Z^a W^b] = a! [a = b]`. Nothing here refers to the repo. -/
structure WickFunctional (E : MvPolynomial (Fin 2) ℂ → ℂ) : Prop where
  add : ∀ P Q, E (P + Q) = E P + E Q
  smul : ∀ (c : ℂ) (P), E (C c * P) = c * E P
  moment : ∀ a b : ℕ, E (X 0 ^ a * X 1 ^ b) = if a = b then (Nat.factorial a : ℂ) else 0

/-- **GMC(2).** If every central power `E(P^m)`, `m ≥ 1`, vanishes, then
`E(Q · P^m)` vanishes for all sufficiently large `m`.

Stated for an arbitrary Wick functional, so that reading this file is enough to
know what is being certified. -/
theorem gmc2 (E : MvPolynomial (Fin 2) ℂ → ℂ) (hE : WickFunctional E)
    (P Q : MvPolynomial (Fin 2) ℂ)
    (hnull : ∀ m : ℕ, 1 ≤ m → E (P ^ m) = 0) :
    ∃ N : ℕ, ∀ m ≥ N, E (Q * P ^ m) = 0 := by
  sorry

end GMC2Challenge

end
