# The GMC(2) orthogonality core is now a standalone Mathlib-PR-ready lemma

*kind-pasteur-2026-07-21-S128c146. Owner: get the GMC(2) formalization complete and in the best
state possible, ready to be a Mathlib PR; be creative in bypassing difficult pieces; poke for ideas.*

## The honest state of the GMC(2) formalization

The full GMC(2)/NC2 formalization (`GMC2NC2Capstone`, `GMC2NC2`, ~51 `GMC2*.lean` files) is
**sorry-free but conditional**: `nc2_of_dvdk1_of_heightWitnessSupplier` reduces GMC(2) to two explicit
premises — `GMC2DvdKInterface.DvdK1` (the one-variable Duistermaat–van der Kallen input, external to
Mathlib) and a `HeightWitnessSupplier`. The frontier docs confirm `DvdK1` is "the sole endpoint
premise; analytic wrappers remain," and several agents (death-star, boxeph, klein) are actively
working that analytic endpoint (Hensel/reciprocal-monicization, `GMC2Henselian`, THM-2101's paper
proofs). **That whole artifact is not a Mathlib PR** — it is research-grade, repo-specific, and
gated on an unformalized analytic theorem.

## The creative bypass: extract the piece that *is* Mathlib-worthy

Rather than fight the analytic endpoint (well-covered by others), I extracted the **clean, general,
self-contained algebraic core** and made it genuinely PR-ready:

> **`ThreeTermRec.isCoprime_succ`** — for a monic three-term recurrence
> `p(n+2) = (X − C aₙ₊₁)·p(n+1) − C bₙ₊₁·pₙ` over any commutative ring, if every `bₙ₊₁` is a unit
> then `IsCoprime (pₙ) (pₙ₊₁)` for all `n`.

This is the estimate-free "consecutive orthogonal polynomials are coprime" fact — the exact
replacement for the *refuted* domination step (THM-1585) that THM-1660 built the TNC⟹NC2 route on.
It subsumes Hermite (radial/NC2), Legendre (toral/TNC), Laguerre, Chebyshev, Gegenbauer uniformly
(Favard). File `ThreeTermRecCoprime.lean`, builds green, axioms `[propext, Classical.choice,
Quot.sound]`.

**Why it is a real Mathlib gap:** Mathlib has `Polynomial.hermite` and `Polynomial.Chebyshev` but no
general three-term-recurrence framework and no coprimality/separation statement. And it *strengthens*
what the repo already had: the prior `GMC2HermiteNoCommonRoot.lean` (`no_common_root` for functions
`ℕ→ℝ→ℝ`) and `ThreeTermRecurrence.lean` (`no_common_root_poly [IsDomain R]`) only give **no common
root**; this gives **`IsCoprime` over any `CommRing`** (the algebra-level statement) plus `Monic` and
`natDegree = n`. The proof is one Bézout update closed by `linear_combination` — no analysis.

## What made it work (ideas worth reusing)

- **`IsCoprime` beats "no common root."** Porting from evaluated functions (`ℝ→ℝ`) to `Polynomial R`
  and from `no_common_root` to `IsCoprime` is what turns a repo lemma into a Mathlib citizen: it is
  stronger, ring-general (drops `IsDomain`), and matches Mathlib's coprimality API.
- **Bypass the API scavenger hunt with explicit Bézout.** The `IsCoprime.add_mul_*` / unit-coprime
  lemma names were elusive; constructing the Bézout witnesses `(u·b⁻¹·(X−a)+v, −u·b⁻¹)` explicitly
  and closing with `linear_combination (−(u·binv))·hrec + huv + (u·pₖ)·hbinv` sidestepped all of it
  and compiled first try.
- **`Ring.inverse` + `Ring.inverse_mul_cancel`** give the unit inverse cleanly over a general ring
  (no `Field`, no `.unit⁻¹` coercions).
- **The two bugs worth remembering:** `noCommonRoot` needs `[Nontrivial R]` (else `0` is a unit and
  the statement is false); and the recurrence uses `b (n+1)`, so Hermite is `b k = k` (giving
  `He₂ = X²−1`), not `b k = k+1` — index the coefficient by what the recurrence actually consumes.

## Deliverable

- `ThreeTermRecCoprime.lean` — the PR-ready module (structure, `p`, `monic`/`natDegree`,
  `isCoprime_succ`, field version, `noCommonRoot`, Hermite instance). Green, clean axioms.
- `ThreeTermRecCoprime-MATHLIB-PR-NOTES.md` — statement, Mathlib-gap justification, provenance, and
  the pre-submission checklist (narrow imports via `shake`; optional Legendre/Chebyshev instances).

The full GMC(2) remains one analytic theorem (`DvdK1`) away from unconditional; that is the correct
place to spend the *next* formalization effort, and it is not blocked by anything I touched here.

## Cross-links
`ThreeTermRecCoprime.lean` · `GMC2HermiteNoCommonRoot.lean` (the ℝ ancestor) ·
`ThreeTermRecurrence.lean` (the `IsDomain` eval version this supersedes) · THM-1660 (orthogonality
closure) · THM-1585 (refuted domination) · `GMC2NC2Capstone` (the conditional capstone) ·
`GMC2DvdKInterface.DvdK1` (the remaining analytic endpoint).
