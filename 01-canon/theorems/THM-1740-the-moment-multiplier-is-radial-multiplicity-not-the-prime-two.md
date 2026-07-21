---
id: THM-1740
title: "THE MOMENT-COUNT MULTIPLIER IS THE RADIAL MULTIPLICITY, NOT THE PRIME 2 — correcting my own THM-1725. Challenging the 'factor 2' by the owner's prime-family lens (a constant in an equation is the x-th member of a family; ask the equation at every member) reveals the 2 was the r=2 case of M* = r·m₀. EXACT LAW (verified r=1,2,3): for the single-straddle family [Z^q, W, ZW², …] where one charge is carried by r monomials and the opposite by one, the minimal certifying moment count is M* = r·m₀ with m₀ = (max coprime charge-pair sum). r=1 ⟹ M*=m₀ (NO factor, for EVERY q prime or composite — verified q=1..7, refuting the guess that composite m₀ needs an extra level); r=2 ⟹ 2m₀; r=3 ⟹ 3m₀. The levels m₀, 2m₀, …, r·m₀ are exactly where the r same-charge coefficients get pinned (opus THM-1685's 'primitive + second level' generalized to r levels). CONSEQUENCE FOR THM-1725: its bound M* ≤ 2·max-pairs is valid ONLY for radial multiplicity ≤ 2 (all trinomials, since a two-sided trinomial's busiest charge has ≤ 2 terms) and is REFUTED at multiplicity ≥ 3 (k ≥ 4): [Z²,W,ZW²,Z²W³] (charge −1 ×3, m₀=3) needs M*=9 = 3·3 > 2·3. The corrected uniform bound is M* ≤ (max radial multiplicity)·(max-pairs). And since q+1 sweeps ALL integers and r sweeps ALL positive integers, the moment counts are the full product family {r}×{charge-return-levels} — there are no special constants, exactly as the owner said."
status: >
  EXACT LAW M* = r·m₀ VERIFIED for r = 1,2,3 (single-straddle family, exact Gröbner over ℚ);
  the r=1 case verified for m₀ = 2..8 (q = 1..7), showing NO dependence on primality of m₀.
  The CORRECTION to THM-1725 is PROVED by the explicit multiplicity-3 counterexample
  ([Z²,W,ZW²,Z²W³], M*=9 > 2·3). The general corrected bound M* ≤ (max mult)·(max-pairs) is
  CONJECTURAL (verified on the tested families); the exact-equality law is established only for
  the single-straddle family, not for general multi-charge patterns.
source: mac-mini-2026-07-20-S150 (owner: look at coprime concepts across the repo and use the
  prime-family lens to challenge assumptions — a constant like 2 is the x-th member of a
  family; read the equation at every member; there are no coincidences)
corrects: THM-1725   # the factor-2 bound is the multiplicity-2 slice only
depends_on:
  - THM-1725  # the moment-count bound whose factor this corrects
  - THM-1685  # opus: primitive + second level -- generalized here to r levels
related:
  - THM-1650  # the ESV first-return level m0
  - THM-415   # prime-modulus vanishing (the mechanism I tested and had to refine)
  - HYP-8540  # the uniform bound -- now with the multiplicity replacing the 2
---

# THM-1740 — the moment multiplier is radial multiplicity

## The challenge, and what it found

THM-1725 bounded the certifying moment count by `M* ≤ 2·max_pairs (p+|n|)/gcd`, with the `2` =
opus's `CT(m₀)+CT(2m₀)`. The owner's directive — *a constant like `2` is the `x`-th member of
a family; ask the equation at every member; there are no coincidences* — says: don't trust the
`2`. Testing the multiplier `M*/m₀` across an explicit family settles it.

**First, a refuted guess.** I predicted (from THM-415, *prime modulus = no collision, composite
= collision*) that composite `m₀` would need an extra level. **False.** For
`[Z^q, W, ZW]` (charges `+q, −1, 0`, so `m₀ = q+1`), `M* = q+1` **exactly**, for every
`q = 1..7` — prime `m₀` and composite `m₀` alike. The single-term-per-charge pattern always
saturates at the primitive level. The primality of `m₀` does **not** drive the count.

**What does: radial multiplicity.** Let one charge be carried by `r` monomials (at `r`
different radial degrees) and the opposite charge by one. Then:

| pattern | busiest-charge multiplicity `r` | `m₀` | `M*` | `M*/m₀` |
|---|---|---|---|---|
| `[Z²,W,ZW]` | 1 | 3 | 3 | 1 |
| `[Z²,W,ZW²]` | 2 | 3 | 6 | 2 |
| `[Z²,W,ZW²,Z²W³]` | 3 | 3 | 9 | 3 |
| `[Z⁴,W,ZW²]` | 2 | 5 | 10 | 2 |
| `[Z³,W,ZW²]` | 2 | 4 | 8 | 2 |

> **EXACT LAW (single straddle): `M* = r·m₀`**, `r` = number of monomials sharing the busiest
> charge, `m₀` = max coprime charge-pair sum. Verified `r = 1,2,3`.

The levels `m₀, 2m₀, …, r·m₀` are exactly where the `r` same-charge coefficients get pinned —
opus THM-1685's "primitive + second level" (`CT(m₀)+CT(2m₀)`, two coefficients) generalized to
`r` levels for `r` coefficients.

## The correction to THM-1725

THM-1725's `2·max-pairs` is the **multiplicity-2 slice only**. It held for all 132 trinomials
because a genuinely two-sided *trinomial* has its busiest charge carried by at most **2** terms
(three monomials, at least one of each sign). It is **refuted at multiplicity ≥ 3**, i.e.
`k ≥ 4`:

> `[Z²,W,ZW²,Z²W³]` (charge `−1` carried by 3 terms, `m₀ = 3`) needs `M* = 9 = 3·3 > 2·3`.

> **Corrected uniform bound: `M* ≤ (max radial multiplicity)·(max coprime charge-pair sum)`.**

So HYP-8540's factor is the **multiplicity**, not `2`. This does not change the *decidability*
of THM-1725(A) (each fixed `(k,D)` is still a finite Gröbner test, unconditionally); it corrects
the *bound* that a uniform proof must target.

## The prime-family reading — no coincidences

`m₀ = q+1` sweeps **every integer ≥ 2** as `q` runs (each realized by the explicit pattern
`[Z^q,W,ZW]`), and `r` sweeps **every positive integer**. So the moment counts `M* = r·m₀` are
the **full product family** `{r} × {return levels}` — every value is realized, and it factors
into two structural integers. The `2, 3, 7` that showed up in THM-1725 were incidental values
of `m₀` or `r` for particular patterns, not constants of the theory. The single equation
`M* = r·m₀`, read at every `(r, m₀)`, **is** the family. That is exactly the owner's point.

## Honest scope

- **Exact `M* = r·m₀` is established for the single-straddle family only** (one charge with
  multiplicity `r`, the opposite with multiplicity 1), `r ≤ 3`. For general multi-charge,
  multi-multiplicity patterns the claim is the **bound** `M* ≤ (max mult)·(max-pairs)`, verified
  on what was tested, not proved.
- **The refutation of THM-1725's factor-2 is rigorous** (explicit `M*=9` counterexample at
  `k=4`); the *replacement* bound is conjectural (HYP-8540, corrected).
- The refuted primality guess is recorded because it was the natural THM-415 prediction and its
  failure is informative: the moment count is driven by coefficient-counting (multiplicity),
  not by the arithmetic of `m₀`.
- No claim that `M*` equals `r·m₀` when several charges each carry multiple terms — cross terms
  between distinct straddles are untested here and are the obvious next case.

*Artifacts:* `04-computation/moment_multiplier_prime_family_macmini_S150.py` and the
`moment_multiplier_*_macmini_S150.out` runs.
*Credits:* the owner's prime-family directive (which produced the correction); opus THM-1685
(the primitive+second level, generalized); THM-1725 (the bound corrected).
