# The slope ladder of Pascal's triangle: one family, a ladder of growth constants, and HYP-614's φ as its s=1 rung

**Source:** kind-pasteur-2026-06-15. Dispatch (user): extend the family
`2^n / Fibonacci / construction-3` (sums of Pascal diagonals) and find the
recursion in `n`.

## The family is exactly one object: Pascal's diagonals at slope s

All three of the user's constructions are **diagonal sums of Pascal's triangle**,
indexed by a single slope parameter `s`:

> `a_s(n) = Σ_k C(n − s·k, k)`,  with the Pascal-identity recurrence
> **`a_s(n) = a_s(n−1) + a_s(n−s−1)`.**

- `s = 0` (full rows `Σ_k C(n,k)`): **`2^n`** — recurrence `a(n)=2a(n−1)` (A000079).
- `s = 1` (shallow diagonals `Σ_k C(n−k,k)`): **Fibonacci** — `a(n)=a(n−1)+a(n−2)` (A000045).
- `s = 2` (`Σ_k C(n−2k,k)`): the user's **construction-3** = **Narayana's cows**,
  `a(n)=a(n−1)+a(n−3)` (A000930): `1,1,1,2,3,4,6,9,13,19,28,…`.
- `s = 3` (next member, asked for): `a(n)=a(n−1)+a(n−4)` (A003269):
  `1,1,1,1,2,3,4,5,7,10,14,19,…`, coefficients `1 / 1 / 1 / 1 / 1+1 / 1+2 / 1+3 /
  1+4 / 1+5+1 / 1+6+3 / 1+7+6 / 1+8+10 / …`.
- general `s`: `a(n)=a(n−1)+a(n−s−1)` (the "lag-`(s+1)` Fibonacci", A003520 at s=4, …).

(All verified to n=40, recurrence and coefficient breakdowns matching the user's
`+`-expansions exactly: `04-computation/pascal_slope_family_growth_ladder_kps.py`.)

## The recursion in n, read three ways

1. **Pascal.** `C(n−sk,k) = C(n−1−sk,k) + C(n−1−sk,k−1)`; the first piece reindexes
   to `a_s(n−1)`, the second (shift `k→k−1`) to `a_s(n−s−1)`. The recurrence *is*
   Pascal's rule restricted to the slope-`s` line.
2. **Hard-core gas.** `a_s(n)` counts configurations where "site `n` is empty
   (`a_s(n−1)` ways) or occupied, forcing the previous `s` sites empty
   (`a_s(n−s−1)` ways)" — the **1-D hard-core lattice gas with exclusion radius
   `s` at fugacity 1** (binary strings with all 1s at least `s+1` apart; the
   independent sets of the `s`-th power of a path). `s=0` = no exclusion = `2^n`;
   `s=1` = no two adjacent = Fibonacci.
3. **Transfer matrix / characteristic.** The recurrence's characteristic equation
   is **`x^{s+1} = x^s + 1`**; the dominant root `β_s` is the per-site growth rate.

## The payoff: a ladder of growth constants

`β_s` = dominant root of `x^{s+1} = x^s + 1`, a strictly decreasing ladder of
algebraic "metallic/plastic" constants → 1:

| s | recurrence | sequence | β_s | name |
|---|---|---|---|---|
| 0 | a(n)=2a(n−1) | `2^n` | 2 | — |
| 1 | a(n)=a(n−1)+a(n−2) | Fibonacci | **1.61803** | golden φ (x²=x+1) |
| 2 | a(n)=a(n−1)+a(n−3) | Narayana | **1.46557** | supergolden (x³=x²+1) |
| 3 | a(n)=a(n−1)+a(n−4) | A003269 | **1.38028** | (x⁴=x³+1) |
| 4 | a(n)=a(n−1)+a(n−5) | A003520 | **1.32472** | **plastic** ρ |
| ∞ | — | — | → 1 | — |

The **plastic number** appears at `s=4` exactly: `x⁵−x⁴−1 = (x³−x−1)(x²−x+1)`, so
`ρ³=ρ+1 ⟹ ρ⁵=ρ⁴+1` (verified). The golden ratio (`s=1`) and the plastic number
(`s=4`) — the two most famous "self-similar growth" constants — are both rungs of
this one slope ladder.

The user's **auxiliary** "same-pace" sequences are the central/dominant *column*
of each diagonal: for `s=0` they are exactly the **central binomials
`C(n,⌊n/2⌋)`** (A001405) `1,1,2,3,6,10,20,35,70,126`. These grow as `β_s^n / √n`
— they carry the sub-exponential (Wallis/`√(πn)`) correction on top of the
exponential `β_s^n` that the main sequence and the dominant root share. (The exact
indexing of the `s≥1` auxiliaries the user listed is convention/transcription
sensitive; the role — the width/dominant column — is the stable content.)

## Why this lands in a tournament-parity repo

The whole family is **independence polynomials of path powers evaluated at fugacity
1**, and the repo's central object is an independence polynomial at fugacity 2:

- **HYP-614** identified `φ` (= the Dedekind regulator `R = log φ` of `ℚ(√5)`) as
  *the* growth rate controlling `H`, the Ising free energy, and the Lyapunov
  exponent. This reflection places `φ` as the **`s=1` rung of the ladder `β_s`** —
  there is a whole family of growth constants, one per exclusion radius, and the
  golden ratio is just the nearest-neighbour case. (HYP-2518: is there a
  Dedekind-regulator / arithmetic meaning for the supergolden and plastic rungs,
  paralleling HYP-614's `log φ`?)
- **THM-485** (claudebox, "two temperatures") is the **fugacity axis** of the
  `s=1` rung: `I(P_n,x) = I(P_{n−1},x) + x·I(P_{n−2},x)`, with `x=1 →` Fibonacci/φ
  (Zeckendorf) and `x=2 →` the repo's `H`/Jacobsthal. The user's family is the
  `x=1` slice across all slopes; THM-485 is the `x`-slice at `s=1`. Together they
  fill a **2-parameter grid `(s, x)`** of hard-core partition functions, with the
  OCF `H = I(Ω,2)` living at `x=2` (on the conflict graph rather than a path power).
- **"Everything is the triangle."** Pascal's triangle *is* the triangle; the slope
  `s` is the analogue of the repo's reduction Mode (the direction you read the
  staircase). Reading the one triangle at slope `s` produces the growth-constant
  ladder — the same way the repo's staircase produces √2, π, e, γ from its sides.

## The clean takeaways

1. The user's `2^n / Fib / construction-3` are `s = 0,1,2` of **one** family
   `a_s(n)=Σ_k C(n−sk,k)`, recurrence `a_s(n)=a_s(n−1)+a_s(n−s−1)`, characteristic
   `x^{s+1}=x^s+1`.
2. The growth constants form a **named ladder** `2 ⊃ φ ⊃ supergolden ⊃ … ⊃ plastic
   ⊃ … → 1` (golden at s=1, plastic at s=4).
3. The family is the **1-D hard-core gas / path-power independence polynomial at
   fugacity 1**; HYP-614's φ is its `s=1` rung and THM-485's two-temperatures is
   its `s=1` fugacity axis — the repo's `H = I(Ω,2)` is the same species at `x=2`.

Cross-links: HYP-614 (φ = Dedekind regulator = the s=1 growth rate), THM-485
(two-temperatures = the fugacity axis at s=1), THM-002 (OCF, the independence
polynomial), [[triangle_foundation]] ("everything is the triangle"; constants from
the sides), [[the-triangular-number-is-the-n4-metagraph-kps]] (the prior
Pascal/tournament thread). HYP-2518 (the arithmetic of the higher rungs).
