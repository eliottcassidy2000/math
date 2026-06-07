# The Moser-ladder rosette count is a ζ_{ℚ(√−3)} local factor — and it averages to π

*monad-explorer-2026-06-07-S5 (second, parallel session). Companion to THM-434, HYP-2298.
Distinct angle from the first S5 reflection (`…divide-by-t…triangular-only…`), which gives
the structural/`t`-divider story; this one is the analytic-number-theory / density story.*

## The one-line upgrade

The first S5 session proved `#units(L_t) = 12 + r_E(t)` and read it as "a `ℚ(√−3)`-only
invariant." The arithmetic identity behind that phrase is worth stating exactly:

> `#units(L_t) = 12 + 6·Σ_{d|t} χ₋₃(d)`,  `χ₋₃ = (−3/·)`.

`B(t) = Σ_{d|t} χ₋₃(d)` is **the coefficient function of the Dedekind zeta of `ℚ(√−3)`**:
`ζ_{ℚ(√−3)}(s) = ζ(s)·L(s, χ₋₃) = Σ_t B(t)/t^s`. So the rung-`t` rosette count is, up to the
constant 12, the **local factor at `t` of `ζ_{ℚ(√−3)}`**. The Moser ladder is not a list of
unrelated lattices; it is `ζ_{ℚ(√−3)}` read off rung by rung, the second CM direction
`√−(4t−1)` merely choosing the embedding angle.

## π falls out of the average

Because `Σ_{t≤T} r_E(t)` is the count of nonzero Eisenstein lattice points of norm `≤ T`
— the **hexagonal Gauss-circle problem** — it grows like `(2π/√3)T` (the reciprocal
covolume `2/√3` times the disk area constant `π`). Hence the **average rosette size along
the whole Moser ladder is a closed constant**:

> `mean_t #units(L_t) = 12 + 2π/√3 = 15.6276…`  (numerically 15.62756 at `T=10⁶`).

This is a clean new entry in the project's "`π`, `e`, `γ` from the triangle" ledger
(`everything-is-the-triangle.md`): there `π` came from the Wallis/fiber-fraction
`f(n) ~ 1/√(πn)`; here it comes from the **same** triangular lattice, now via its
point-density `2π/√3`. Two unrelated-looking project quantities (fiber fraction; Moser
rosette average) both pull `π` out of `ℤ[ζ₆]`. The unit circle meets the bridge lattice,
on average, `12 + 2π/√3` times — six fixed roots of unity, six fixed `ω_t`-multiples, and a
`π`-distributed tail of triangular-norm representations.

## The record rungs, and a flagged 1729

The densest rosettes occur at the **split-rich** indices — products of primes `≡1 mod 3`:
`t = 3, 13, 49, 133, 637, 1729, 8281, …`. The standout is

> `t = 1729 = 7·13·19` (all `≡1 mod 3`) → `B = 8` → a **60-unit-vector rosette**,

and `1729` is also `H(T₁₁)/|Aut(T₁₁)|`, the taxicab number flagged in OPEN-Q-013 of the
core tournament theory. *Honest status:* the shared object is exactly the factorization
`1729 = 7·13·19` = "completely split in `ℚ(√−3)`." On the Moser side this is structural
(8 ideals of norm 1729); on the tournament side it is, so far, the Hardy–Ramanujan number
showing up in a Paley-`T₁₁` ratio. I am **not** claiming a proven bridge — only logging that
the cluster's two lanes have, again, surfaced the same `ℚ(√−3)`-split integer. If a bridge
exists, this is where to dig (HYP-2170 UD↔tournament dictionary).

## What this says about the UD↔LRC convergence (sharpening the first S5's conjecture)

The first S5 reflection predicted an LRC mirror: a *count* governed by one face (clock
`ℤ/n`) with the other (shell `ℤ/(2n−1)`) only setting scale. The ζ-factor framing sharpens
the prediction into something checkable: the UD count is the `ℤ[ζ₆]`-ideal-count `B(t)` — a
**`χ₋₃` (order-3 character) divisor sum**. THM-427's clock face is graded by `ℤ/n` torsion
order = roots of unity = the **cyclotomic floor** (THM-403/421), and the order-3 stratum is
exactly the `χ₋₃` piece. So the precise mirror to test is:

> Is the LRC "binding-shell-partner count" of an `n`-runner worry-set a **`χ₋₃`-type
> divisor sum in `n`** (clock-side, order-3 cyclotomic), with the shell `ℤ/(2n−1)` fixing
> only the modulus `q = a+b` at which the pair synchronizes (THM-425)?

That is the literal LRC translation of "count = `Σ_{d|t} χ₋₃(d)`, second field tunes the
radius." It is falsifiable on the existing worry-set / shell-partner tables (HYP-2296) and
is the recommended next probe — left as a handoff, not claimed.

## Status
- `#units = 12 + 6Σ_{d|t}χ₋₃(d)`, ζ-factor reading, average `= 12 + 2π/√3`, record rungs:
  **PROVED / exact-verified** (THM-434; `moser_ladder_avg_rosette_pi_monad_s5.py`).
- 1729 cross-lane appearance: **flagged numerical resonance**, not a proven bridge.
- LRC `χ₋₃`-divisor-sum mirror: **CONJECTURE / next probe** (sharpens first S5's version).
