# HYP-2270 — The additive–multiplicative exclusion: rationality vs irrationality = the resonance/two-faces structure

**Session:** claudebox-2026-06-03-S629. **Prompt (user):** integrate the fact that not both of `π+e`, `π·e` are
rational; see rationality/irrationality like odd/even and multiplication/addition; find the deep underlying
structure with the previous work. **Threads:** HYP-2150 (two faces), HYP-2175 (Collatz), HYP-2235 (CM/norm-1),
HYP-2245 (α₁/α₂), HYP-2265 (π/3 cyclotomic).

## The exclusion (formalized)
`(sum, product) = (a+b, a·b)` are the elementary symmetric functions of a pair — the coefficients of the monic
quadratic `(X−a)(X−b) = X² − (a+b)X + ab`. **If both the sum (additive invariant) and the product (multiplicative
invariant) are rational, the pair is algebraic** (each is a Vieta root of a rational quadratic). Contrapositive: a
transcendental pair (`π, e`) has **at least one of `π+e`, `π·e` irrational** — the additive and multiplicative
invariants cannot both collapse to `ℚ`. **Formalized (math-lean):** `Math/Transcendence/SumProduct.lean`
`isAlgebraic_of_sum_prod_rat`, `not_sum_and_prod_rat_of_transcendental` (⟹ the `π,e` fact via Lindemann).

## Rational/irrational ↔ even/odd ↔ addition/multiplication — the deep structure
The exclusion is the SAME complementarity as the whole arc: a pair/system has an **additive** invariant and a
**multiplicative** invariant, and they cannot both be "tame" (rational / trivial / resonant) unless the object is
special (algebraic / collapsed / forbidden).

| thread | additive invariant (sum / e₁) | multiplicative invariant (product / e₂) |
|---|---|---|
| `π, e` | `π + e` | `π · e` (not both rational) |
| LRC two faces (HYP-2150) | resonance `Σ mᵢ vᵢ = 0` (mod `2n−1`) | doubling `x↦2x` (mod `n`) |
| Collatz (HYP-2175) | the `+1` in `3n+1` | `×3` vs the 2-adic; the `2^K=3^L` cycle |
| independence polynomial (HYP-2245) | `α₁` = #3-cycles (`e₁`) | `α₂` = disjoint pairs (`e₂`); norm-1 = `α₂=1` |
| CM / norm-1 (HYP-2235) | `ρ₁+ρ₂` | `ρ₁ρ₂ = 1 = |β|²` (conjugation-fixed) |
| parity | — | `H` odd (Rédei), the 2-adic / even-vs-odd |

Rationality is to irrationality as even is to odd as the multiplicative face is to the additive face: the **same
2-fold complementarity** read in transcendence, parity, and the resonance structure.

## The cyclotomic ↔ transcendental spectrum (the new bridge to S628)
The two extremes of the additive–multiplicative axis (verified):
- **Cyclotomic (both rational = algebraic = forbidden/resonant/collapse):** the roots of `Φ₃ = X²+X+1` (the cube
  roots of unity, the `7=Φ₃(2)` of S628) have sum `−1` and product `+1` — **both rational**, the maximally-aligned
  algebraic point. The forbidden H-values, the LRC collapse, the resonance live here (`Σmᵢvᵢ=0` IS "the symmetric
  functions are special").
- **Transcendental (at least one irrational = generic/lonely):** `π, e` — the maximally-misaligned point, where the
  additive and multiplicative invariants refuse to both be rational; loneliness / genericity lives here.

So **resonance/forbidden/collapse = the algebraic (both-symmetric-functions-rational) end; loneliness/genericity =
the transcendental (not-both) end** — one spectrum from cyclotomic to `{π,e}`, the additive–multiplicative axis. The
perspective key (free vs fixed, σ-orbit type) is this axis: the conjugation-fixed (norm-1, product rational) vs the
free/generic.

## Formalized (math-lean, sorry-free) — `Math/Transcendence/SumProduct.lean`
`isAlgebraic_of_sum_prod_rat` (the Vieta root), `not_sum_and_prod_rat_of_transcendental` (the exclusion; ⟹ `π+e`,
`π·e` not both rational).

## Open
- Make the cyclotomic↔transcendental spectrum a quantitative invariant (the "irrationality measure" of the
  symmetric functions as a loneliness/genericity gauge).
- Tie the `α₂=1` (product rational, CM/norm-1) family to the conjugation-fixed end and the forbidden Φ₃ values.
- Collatz: the `2^K=3^L` cycle as the "both-rational" coincidence (additive log-relation + multiplicative) — the
  resonance as a symmetric-function alignment, forbidden by a transcendence/linear-forms-in-logs argument.
