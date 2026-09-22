# Session-lead reflection: the Catalan root of the two sheets, the diagonal behind `{1,4,6}`, and the eight in the level-eleven form

**Status: PROVED elementary identities (hand proofs, each checked by a
ten-line computation reproduced below); CITED classical inputs; SCOPE where
no map exists. No agents, no independent audit: this is a reasoning note
by the session lead, written under a usage limit, and should be read as
provenance for the exact claims it isolates, not as audited canon.**
Collatz remains OPEN. Written 2026-09-22 at the close of session
`collatz-mod6-20260917`; the audited material it leans on is the
[synthesis](collatz_mod6_20260917_synthesis.md) (sections 1--13).

## 1. The trivial cycles of the two sheets are the two sides of Catalan

Under the standard map the plus sheet has `1 -> 4 -> 2 -> 1` and the minus
sheet (equivalently `3n-1` on positives) has `1 -> 2 -> 1`: one halving
longer, because `3*1+1 = 4 = 2^2` and `3*1-1 = 2 = 2^1`. Under the
shortcut `(3n+1)/2` the plus cycle is `1 -> 2 -> 1` and the odd-only map
fixes `1` on both sheets, so the user's convention "ignore the evens" makes
both roots fixed points and hides the only difference, which is the
halving count `K` at the root.

That count is a Catalan datum. The unit-gap clocks `2^K - 3^L = +-1` are
exactly `(K,L) = (2,1)` on the `+1` side and `(1,1), (3,2)` on the `-1`
side (Mihailescu, CITED; the elementary `2^K-3^L=+-1` case is inherited
from the zsigmondy and pillai lanes). Through the inherited cycle gate
`n_0 = bB/(2^K-3^L)`, each unit-gap clock is an integer cycle of the sheet
whose sign it carries:

| clock `(K,L)` | `2^K-3^L` | sheet | word | `B` | cycle |
|---|---|---|---|---|---|
| `(2,1)` | `+1` | plus | `(2)` | `1` | `{1}` |
| `(1,1)` | `-1` | minus | `(1)` | `1` | `{1}` |
| `(3,2)` | `-1` | minus | `(1,2)` | `3+2 = 5` | `{5,7}` |

So "the difference between a root cycle of length 3 and one of length 4 is
all the world" has an exact content: the `+1` side of Catalan owns one
unit-gap clock, the `-1` side owns two, and the second one *is* the cycle
`{5,7}`. The seven-cycle is the first non-unit clock, `11/7` with
`2^11-3^7 = -139` dividing the carry (inherited). This is the entire
sign-specific content at the root; it does not open any door to `-7/4` or
`-29/16`, which are periodic points of quadratic maps and share only the
numeral `3` with the cycle length (SCOPE, as the zsigmondy and row lanes
already recorded).

## 2. `{1,4,6}` is `{1}` together with the doubles of the seeds

**PROVED.** Let `S = {2,3}`. Under the strict rule (distinct summands) the
closure misses exactly `{1,4,6}` (THM-2422). Under the weak rule (equal
summands allowed) the closure misses exactly `{1}`, because `4 = 2+2` and
`6 = 3+3` become available and then every `n >= 7` follows as before.
Hence the two extra holes are precisely the excluded diagonals
`2+2` and `3+3`: `{1,4,6} = {1} cup 2S`. Their only distinct
representations, `4 = 1+3` and `6 = 1+5 = 2+4`, all pass through the
hole module, which is why the module is closed.

Check (reproduced): strict closure of `{2,3}` to `60` misses `[1,4,6]`,
weak closure misses `[1]`.

This is the same mechanism that makes the trivial Collatz cycle
exceptional in the summand reading: `1 -> 1+1 = 2` is the excluded diagonal
(inherited summand note, section 3), and the strict summand shadow's
complement is the doubling forest (its section 1). The user's phrase
"2 cooperates with 1 at `2+2 = 2*2 = 2^2`" is therefore exact: the
coincidence of the additive and multiplicative diagonals at `2` is what
the swap-fixed fibre deletion removes, and `4`, `6` are the first two
targets of the deleted doubling arcs `2 -> 4`, `3 -> 6`. The three
"islands" `1`, `4`, `6` are one unit and two deleted diagonals, not three
chains; the three chains of the square-filtered graph `Q_n` are founded by
`1, 2, 4` (wave six), and the summand lane proved no additive filter can
have founders `{1,4,6}`.

## 3. The minus cycles and the chains of `Q_n`

SCOPE. The minima `1, 5, 17` of the three `3n-1` cycles sit in the chains
of `1` and of `4` (`5` lies in `{4,5,11,12}`) and beyond the linear-forest
range (`17 > 12`); the cycle members `5, 7` lie in different chains. No
map from the chain structure (reflections `x -> 4-x, 9-x, 16-x`) to the
cycle gate was found, and none is claimed. The founders `1, 2, 4` of the
chains are exactly the `n` with no square strictly inside `(n, 2n)`
other than `2n`, a finite set by square spacing; that is the whole reason
there are three chains.

## 4. Eight and eleven in `q prod (1-q^n)^2 (1-q^{11n})^2`

The series `q - 2q^2 - q^3 + 2q^4 + q^5 + 2q^6 - 2q^7 + 0q^8 - 2q^9 - 2q^10
+ q^11 - 2q^12 + 4q^13 + ...` is the weight-two newform of level `11`,
`eta(q)^2 eta(q^11)^2`, the modular form of the elliptic curve `X_0(11)`
(CITED, classical). Its coefficients are multiplicative with
`a_2 = -2, a_3 = -1, a_5 = 1, a_7 = -2, a_11 = 1, a_13 = 4` and the Hecke
recursion `a_(p^(k+1)) = a_p a_(p^k) - p a_(p^(k-1))` for `p != 11`.

**PROVED (from the recursion).** `a_8 = a_2(a_2^2 - 2*2) = (-2)(4-4) = 0`.
So the first vanishing coefficient is at `8` because `a_2^2 = 2p` at
`p = 2`, once more the coincidence `2*2 = 2^2`. Equivalently the Frobenius
angle at `2` is `3 pi/4`: `a_2 = 2 sqrt2 cos(theta)` with
`cos theta = -1/sqrt2`, so `a_(2^k) = 2^(k/2) sin(3 pi (k+1)/4)/sin(3 pi/4)`
and `a_(2^k) = 0` exactly when `4 | k+1`: `k = 3, 7, 11, ...`. The
"eight" is the order of the root of unity `e^(3 pi i/4)`. Reproduced:
`a_(2^k) = 1, -2, 2, 0, -4, 8` for `k = 0..5`; zeros of `a_n` for `n <= 60`
at `8, 19, 24, 29, 38, 40, 56, 57, 58` (the primes `19, 29` have
`a_p = 0`, and every `8m` with `m` odd vanishes by multiplicativity).

The `11` is explained by the genus of the modular curves: `X_0(N)` has
genus `0` for `N <= 10` and genus `1` for `N = 11`, so `11` is the first
level carrying a weight-two cusp form (CITED). This is a different eleven
from the repository's order-eleven tournament banks; no map is claimed
(SCOPE). The user's factor `(1-q^2)^2 (1-q^22)^2` is simply the `n = 2`
term of the product.

The two eights of this note, `2^3 = 3^2 - 1` (the minus-sheet Catalan
clock, Bang's exception `63 = 7*9`) and `2^3` as the order of
`e^(3 pi i/4)`, share only the numeral (SCOPE).

## 5. A square-sum counterexample above `25` and Collatz

There is none: OEIS A090461 records Gerbicz's proof that every `n >= 25`
has a square-sum Hamiltonian path (CITED via the wave-six lane), and the
square filter and the Collatz arrow select disjoint targets (`(3n+1)/2` is
never a square). So no equivalence, implication, or reverse implication
with a Collatz counterexample exists (SCOPE).

What does survive is a typed analogy about *shape*. Gerbicz's step embeds
a chain of `Q_n` into `Q_(25n+12)` by an affine self-similar map with an
odd-square ratio; the Collatz inverse fibre is embedded into itself by the
braid `R(x) = 4x+1`, an affine map with a square ratio, and the `E`-graph's
extra arrows are exactly `R^(-1)`. In `Q_n` the embedding closes the
problem because every vertex has about `(sqrt2-1) sqrt n` choices, so a
finite base propagates; Collatz offers one choice per node. The relaxed
problem where choices exist is the `E`-graph, and there the analogue of
Gerbicz's theorem is the greedy certificate `G` with its proved
density-one stopping. So the square-sum theorem is a model for a proof of
the relaxed reachability `Q2`, not of Collatz: preserved, the
"finite base plus self-similar affine embedding" shape; lost, the choice
count; decisive test, the minus sheet, where the same shape coexists with
three cycles.

## 6. Three finite heads, three scalings, no common fractal

The "microcosm" the user senses is real three times over and is a
different object each time: the summand closure has the finite head
`{1,4,6}` and then the dyadic law `M_t = 27*2^(t-4)+1` (THM-2422); the
square-sum graph has the finite head "leaf `18` until `31`" and then
Gerbicz's `25`-fold self-similarity; the Collatz inverse fibre has the
finite head "the diagonal at `1`" and then the `4`-fold braid. The scalings
`2, 25, 4` are not the same map, and the wave-six lane proved `Q_n` carries
no dyadic law at all. What is common is only the research card: separate
the finite head from the recurrence class before looking for the law.

## Reproduction

The two checks in this note are the ten-line computations printed in
`/tmp/lead_check.txt` at writing time (product expansion of
`q prod (1-q^n)^2 (1-q^{11n})^2` to `q^60`; strict and weak closures of
`{2,3}` to `60`); they are elementary and can be rerun in any Python.
