# Mathlib PR candidate — coprimality of consecutive three-term-recurrence polynomials

**File:** `TournamentH7/ThreeTermRecCoprime.lean` (builds green in this project, 8475 jobs).
**Proposed Mathlib path:** `Mathlib/RingTheory/Polynomial/ThreeTermRecurrence.lean`
(or `Mathlib/Algebra/Polynomial/ThreeTermRecurrence.lean`).
**Axioms:** every public result depends only on `[propext, Classical.choice, Quot.sound]` — no
`sorry`, no `native_decide`.

## What it proves

For a monic three-term recurrence over a commutative ring `R`,
```
p 0 = 1,   p 1 = X - C (a 0),   p (n+2) = (X - C (a (n+1))) * p (n+1) - C (b (n+1)) * p n
```
bundled as `structure ThreeTermRec R := (a b : ℕ → R)`:

- `ThreeTermRec.isCoprime_succ : (∀ n, IsUnit (b (n+1))) → ∀ n, IsCoprime (p n) (p (n+1))`
  — **the main result.** Over any `CommRing`. Consecutive members are coprime.
- `ThreeTermRec.isCoprime_succ_of_ne_zero` — field specialization (`b (n+1) ≠ 0`).
- `ThreeTermRec.noCommonRoot [Nontrivial R]` — consecutive members share no root.
- `ThreeTermRec.monic`, `ThreeTermRec.natDegree_p [Nontrivial R]` — `p n` is monic of degree `n`.
- `ThreeTermRec.hermite`, `hermite_isCoprime_succ` — the probabilists' Hermite family as an
  instance (`a ≡ 0`, `b k = k`).

## Why Mathlib wants it

Mathlib has specific orthogonal-polynomial families (`Polynomial.Chebyshev`, `Polynomial.hermite`)
but **no general three-term-recurrence framework** and, in particular, no statement of the classical
fact that *consecutive orthogonal polynomials are coprime* (equivalently, by Favard, that any monic
family with a nonzero-off-diagonal recurrence has consecutive members sharing no root). This lemma
covers Hermite, Legendre, Laguerre, Chebyshev, Gegenbauer, Jacobi, … uniformly, and the
`IsCoprime` conclusion is the algebra-level statement (strictly stronger than "no common root",
which is the `Nontrivial`-ring eval corollary here).

The proof is an explicit one-line Bézout update (no analysis, no estimates): from
`u·pₙ + v·pₙ₊₁ = 1` and `pₙ = b⁻¹((X−a)·pₙ₊₁ − pₙ₊₂)`, read off a Bézout identity for
`(pₙ₊₁, pₙ₊₂)` — witnesses `(u·b⁻¹·(X−a) + v, −u·b⁻¹)`, closed by `linear_combination`.

## Provenance / relation to prior art

This generalizes and strengthens the project's `GMC2HermiteNoCommonRoot.lean` (kps-S128c121:
`ThreeTerm.no_common_root` for the evaluated family `ℕ → ℝ → ℝ`) and the evaluated
`ThreeTermRecurrence.lean` (`no_common_root_poly [IsDomain R]`): here the conclusion is `IsCoprime`
over any `CommRing` (not `IsDomain`), plus the `Monic`/`natDegree` facts. It is the estimate-free
core that replaces the refuted domination step in the repo's TNC ⟹ NC2 ⟹ GMC(2) route
(THM-1585 refuted klein's domination; THM-1660 replaced it with this orthogonality closure).

## Pre-submission checklist

- [x] Sorry-free, clean axioms, builds green.
- [x] General `CommRing`; field + `Nontrivial` corollaries; monic/degree; worked instance.
- [ ] **Narrow imports** — currently `import Mathlib`; run `lake exe shake` / `#min_imports` to
      reduce to the ~6 needed files (Coprime.Basic, Polynomial Monic/Degree/Eval, LinearCombination).
- [ ] Rename structure per reviewer taste (`ThreeTermRec` vs `Polynomial.ThreeTermRecurrence`).
- [ ] Consider connecting to `Polynomial.hermite` (prove `hermite.p = probabilists' Hermite`) and to
      a future `OrthogonalPolynomial` API.
- [ ] Add Legendre/Chebyshev instances if reviewers want more than one worked family.
