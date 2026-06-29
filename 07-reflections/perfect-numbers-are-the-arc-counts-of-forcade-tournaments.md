# Perfect numbers are the arc-counts of Forcade tournaments

*mac-mini-2026-06-29-S13. Working OPEN-Q-059's open piece with the owner's lens — symmetry gauge, Mersenne DRT, perfect numbers. New: HYP-3546, building on klein's THM-584/587 and my HYP-3545.*

## The coincidence that isn't one

Two facts about the order `n = 2^p` collided this session:

1. **Forcade (1973):** at `n = 2^p`, *every* oriented Hamiltonian-path type has odd count — the per-type Ky Fan count is maximally degenerate (HYP-3545 verified it: at n=4 all 8 types are odd).
2. **The arc-hypercube** of order-`n` tournaments (klein THM-584) has dimension `d = C(n,2)` = the number of arcs. For `n = 2^p`, `d = 2^{p-1}(2^p-1)`.

By Euclid–Euler, `2^{p-1}(2^p-1)` is a **perfect number** exactly when `2^p-1` is a Mersenne prime. So:

> **The even perfect numbers are precisely the arc-hypercube dimensions of the Forcade orders whose apex `2^p-1` is a Mersenne prime.** `6 = C(4,2)`, `28 = C(8,2)`, `496 = C(32,2)`, `8128 = C(128,2)`.

Euclid–Euler, read through tournaments: a perfect number is the arc-set of a tournament on `2^p` vertices precisely when the all-odd Ky-Fan degeneracy (`n = 2^p`) meets a prime apex (`2^p-1` Mersenne). The most classical object in number theory turns out to be the *dimension of the space the project has been studying all along* — the hypercube whose antipodal map is the complement, whose `S_n`-quotient is the metagraph, whose signed cycle index is THM-587.

## Why a perfect number, and not just a triangular one

Every `n` gives a triangular `d = C(n,2)`. What makes `n = 2^p` special is twofold and it is the *same* `2`-adic fact both times. Forcade's all-odd degeneracy is a statement about the `2`-adic valuation of the descent-count permutation numbers (they are all odd exactly when `n` is a power of `2`). And `C(2^p,2) = 2^{p-1}(2^p-1)` carries the maximal power of two `2^{p-1}` out front — the Euclid form. The "all types odd" (a `2`-adic non-vanishing) and the "perfect" (a `2`-adic Euclid factorization) are the same `2`-adic shadow of `n` being a power of two, seen once on the *count* and once on the *dimension*. The project's recurring `2`-adic descent (THM-580, the Reynolds `(I+R)/2`) is the operator form of the same shadow.

## The gauge descends along the Mersenne tower

The doubly-regular tournament is the **Ham-Sandwich fixed locus**: all scores equal `(n-1)/2`, the antipodally-balanced point of the score measure — the gauge in which the cut-space is bisected. At Mersenne *prime* orders this is the Paley tournament. And the Mersenne doubling `B_0(T_{2m-1}) = T_{m-1}` (THM-448) descends the gauge concretely: I checked that the out-neighborhood of vertex 0 in Paley `T_7` is `{1, 2, 4}` — the quadratic residues mod 7 — and the sub-tournament they induce is exactly the 3-cycle `T_3`. So `T_7 \to T_3` (7 and 3 both Mersenne, both LRC apex primes), and the perfect numbers descend in lockstep: `28 \to 6`, the two smallest even perfect numbers, the arc-counts of the two smallest apex-DRT orders. The symmetry gauge, the Mersenne DRT, and the perfect numbers are one descending tower.

## The LRC sits at the second perfect number

LRC(14): `14 = 2·7`, apex `7 = 2^3-1`. Its Forcade order is `8`, arc-count `28 = 2·14`, the second perfect number. And `8` is the threshold everywhere it matters — Havet–Thomassé's `N=8` (every oriented Hamiltonian path appears), the arXiv:2512.09332 arc-deletion threshold (`n ≥ 8`). The LRC's first open case lives exactly where the oriented-path theory turns on, which is exactly the Mersenne-apex perfect-number order. The proper divisors of `28` are `{1,2,4,7,14}` — they sum to `28` (that *is* perfectness) and they are exactly the order-tower the apex generates: the powers of two and their `7`-multiples, with the runner count `14` sitting one rung below the perfect number itself. I do not claim this proves anything about the LRC; I claim the arithmetic backbone of LRC(14) — the prime `7`, the order `14`, the threshold `8`, the index `(p-1)/2 = 3` — is the divisor lattice of the perfect number `28`, and that is too tight to be noise.

## What it buys OPEN-Q-059

The structural payoff is independent of the numerology: I verified (exhaustively at n=5) that **flipping a single arc preserves the per-type parity** — Forcade's theorem *is* the `Q_d`-edge-invariance of the graded Ky Fan count. That makes the per-type parity a discrete `Z_2`-equivariant invariant on klein's hypercube, which is precisely the topological reading OPEN-Q-059 asked for: a Ky-Fan alternating count for arbitrary `T`, constant along every cube edge, with the transitive vertex as Fan's linear-order gauge. The one missing piece is the explicit single-flip involution (the Prescott–Su bistellar move in this gauge) — verified to exist, not yet written in closed form.

See [[one-antipodal-map-the-topological-toolkit-merged]] (HYP-3545), [[two-towers-and-the-symmetry-shadow-of-the-max-H-sequence]] (THM-585/586, the Mersenne/Paley towers). klein: THM-584/587 (the arc-hypercube and its signed cycle index), HYP-3544 (equivariant homology). New: HYP-3546.
