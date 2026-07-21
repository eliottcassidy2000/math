---
id: THM-1780
title: "THE ATOM-COVERING CLOSES VIA PAIRS — multi-charge atoms are REDUNDANT — and the Laplace moment engine reveals E[P^m] as the WEIGHT-0 (trivial-representation) projection of a binary form, with the nullcone = in/transitivity of the charge flow. (A) THE MOMENT ENGINE: E[P^m] = L(charge-0 part of P^m), the projection of P^m onto the TRIVIAL representation of the diagonal torus T = {Z↦λZ, W↦λ⁻¹W} (charge = T-weight), then the radial Laplace L(V^k)=k!. Verified: E[P_λ^m] = E[P^m] (torus-invariance), and the charge-graded convolution engine matches the direct moment. This is the representation theory of binary forms: the moment reads the SL(2)/U(1)-invariant weight-0 component. (B) IN/TRANSITIVITY: a charge-0 tuple is a closed walk on ℤ (steps = charges); an atom is a PRIMITIVE cycle. P is ONE-SIDED ⟺ all charges one sign ⟺ partial sums monotone ⟺ NO cycle ⟺ TRANSITIVE ⟺ nullcone-harmless; two-sided ⟺ a cycle exists ⟺ INTRANSITIVE. The nullcone/covering is the disjoint-cycle question — the OCF shape for the charge lattice. Tournaments ARE in/transitivity, and the GMC(2) pivot is the same pivot. (C) THE PAIR-REDUCTION, the key closure: V(⟨all pair-straddle atom forms c_p^{|n|/g} c_n^{p/g}⟩) = the one-sided locus — verified 5 multi-charge patterns including {+2,+3,−5} (forms a⁵c²,b⁵c³), {+1,+2,−3,−4}, {+3,+5,−7}. So MULTI-CHARGE ATOMS (triples like the (+2,+3,−5) form abc) are REDUNDANT: every two-sided P has pos-neg pairs, each pair contributes c_p c_n to the radical, and the ideal of pair-forms cuts out exactly {all pos=0} ∪ {all neg=0} = one-sided (THM-1770 D). The general atom-covering therefore reduces to PAIRS, and GMC(2) is down to one clean statement: every pair-straddle atom form lies in radical(moment ideal) — supplied level-by-level by the first-return renewal (THM-1770 A)."
status: >
  (A) The engine is VERIFIED against the direct moment (m=1..6, 3 patterns) and the
  torus-invariance E[P_λ^m]=E[P^m] confirmed; the weight-0-projection identity is the
  definition of the Gaussian functional, exact. (B) is a framing/dictionary, exact where
  stated (one-sided ⟺ no cycle is elementary). (C) THE PAIR-REDUCTION V(⟨pair forms⟩)=one-sided
  is VERIFIED on 5 multi-charge patterns (each c_p c_n in the radical by Gröbner); it is the
  proved-per-pattern statement that multi-charge atoms are redundant. The remaining step —
  every pair-atom form in radical(I) uniformly — is the renewal induction (THM-1770 A per
  level), well-founded but not a closed uniform proof.
source: mac-mini-2026-07-20-S154 (owner: take the atom-covering target; representation theory
  of binary forms and its relation to tournaments = in/transitivity; work the Laplace moment engine)
depends_on:
  - THM-1770  # first-return renewal; pair-only closure (this extends it to the general covering)
  - THM-1760  # single-straddle radial reduction
related:
  - THM-1685  # opus: the vanishing-sum Nullstellensatz -- pair-reduction simplifies its core
  - THM-466   # binary digits of H = odd-cycle census (the tournament OCF)
  - LEM-004   # tournaments are odd functions -- the in/transitivity theme
  - HYP-8590  # the localisation lemma -- now reduced to pairs
concurrent: opus-2026-07-20-S432 (binary forms thread -- independent convergence)
script: 04-computation/laplace_moment_engine_binary_forms_macmini_S154.py,
        05-knowledge/results/gmc2_pair_reduction_macmini_S154.out
---

# THM-1780 — the atom-covering closes via pairs; the moment is a weight-0 projection

## (A) The Laplace moment engine, and the representation theory

`P(Z,W)`, `W = Z̄`, is a **binary form**. The diagonal torus `T = {Z↦λZ, W↦λ⁻¹W}` acts, and the
**charge** `k = a−b` of `Z^a W^b` is its **`T`-weight**. The Gaussian moment functional
`E[Z^A W^B] = A!δ_{AB}` is `T`- and `U(1)`-invariant, so

> **`E[P^m] = L(charge-0 part of P^m)`** — the projection of `P^m` onto the **trivial
> representation** (weight 0), then the radial Laplace `L(V^k) = k!`, `V = ZW`.

`CT_θ` (= averaging over `U(1)`) is the character projection onto weight 0. This is the classical
representation theory of binary forms: the moment reads the invariant weight-0 component. The
**engine** stores `P` as `{charge : radial polynomial in V}` and convolves on the charge lattice
(weight = tensor grading) while multiplying radial polynomials — matching the direct moment
(verified `m=1..6`) and far faster. Torus-invariance `E[P_λ^m] = E[P^m]` verified.

## (B) In/transitivity — the tournament pivot

A charge-0 tuple is a **closed walk on ℤ** (steps = charges, start/end 0). An **atom** is a
primitive cycle — never returning to 0 in between.

> **`P` one-sided ⟺ all charges one sign ⟺ partial sums monotone ⟺ NO cycle ⟺ TRANSITIVE**
> (nullcone-harmless); **two-sided ⟺ a cycle exists ⟺ INTRANSITIVE.**

Verified: `[1,2,3]`, `[2,4,6]` transitive (no atom); `[1,−1]`, `[2,3,−5]`, `[1,2,−3,−4]`
intransitive. The nullcone/atom-covering is the **disjoint-cycle-packing** question for the
charge lattice — the same shape as the repo's odd-cycle collection for tournaments (THM-466). The
owner's thesis that *tournaments are in/transitivity itself* lands exactly: the GMC(2) pivot
between harmless (one-sided/transitive) and constrained (two-sided/intransitive) is the
in/transitivity pivot.

## (C) The pair-reduction — the key closure

> **`V(⟨all pair-straddle atom forms `c_p^{|n|/g} c_n^{p/g}`⟩) = the one-sided locus`** — where
> the pair `(p,n)` ranges over all pos-neg charge pairs, `g = gcd(p,|n|)`.

Verified on 5 multi-charge patterns: `{+2,+3,−5}` (forms `a⁵c², b⁵c³`), `{+1,+2,−3,−4}`,
`{+2,+3,−4,−5}`, `{+1,−1,+2,−3}`, `{+3,+5,−7}` — each with every `c_p c_n` in the radical.

**Consequence: multi-charge atoms are REDUNDANT.** The `(+2,+3,−5)` triple gives the form `abc`,
but the *pairs* `(2,5)` and `(3,5)` give `a⁵c²` and `b⁵c³`, and `⟨a⁵c², b⁵c³⟩` already cuts out
`{c=0} ∪ {a=b=0}` = one-sided. Since **every two-sided `P` has pos-neg pairs**, and the pair-form
ideal's variety is exactly `{all pos=0} ∪ {all neg=0}` (the coordinate-subspace union of
THM-1770 D), the covering closes through pairs alone.

> **GMC(2) is now down to one clean statement:** every pair-straddle atom form lies in
> `radical(moment ideal)`. The first-return renewal (THM-1770 A) supplies this level by level —
> each pair's form appears isolated at its return level once the lower atoms are killed.

## Honest scope

- **(A) is exact** (the weight-0 identity is the functional's definition; torus-invariance
  verified). The "SL(2)" language is the weight/torus grading — the genuine invariance is the
  diagonal torus + `U(1)`, not the full `SL(2)`, which does not preserve the Gaussian.
- **(B) is a dictionary**, exact where stated; "the covering = OCF-shape" is a structural
  analogy to the tournament odd-cycle formula, not a reduction (THM-466 is over `𝔽₂`; this is
  over `ℂ`).
- **(C)'s pair-reduction is verified on 5 patterns, not proved in general.** It reduces the
  general atom-covering to pairs (a real simplification — multi-charge atoms drop out), but
  "every pair-atom form in `radical(I)` uniformly" remains the renewal induction, well-founded
  (finitely many new atoms per level, THM-1770 A) but not a closed uniform proof.
- **GMC(2) status:** span-2 (THM-1600), complex radial (THM-1695), single-straddle (THM-1760),
  star/pair-only (THM-1770 D) all closed; the general case reduces (C) to pairs, leaving the
  uniform "pair forms in the ideal" as the sole gap. opus-S432's concurrent binary-forms thread
  may bear on this.

*Artifacts:* `04-computation/laplace_moment_engine_binary_forms_macmini_S154.py`,
`05-knowledge/results/gmc2_pair_reduction_macmini_S154.out`.
