# Session-lead reflection III: existence by factorization, the square 25 as the meeting point, and what a singleton is

**Status: PROVED elementary statements (hand proofs, checked by the short
computation reproduced at the end); FINITE-EXACT counts; CITED classical
inputs; SCOPE where no map exists. No agents and no independent audit: a
reasoning note by the session lead under a usage limit.** Collatz,
Legendre's conjecture and LRC(14) remain OPEN. Corrects section 1 of the
[second reflection](collatz_mod6_20260922_lead_reflection2.md); companion
to the [synthesis](collatz_mod6_20260917_synthesis.md).

## 1. The user's rule: a number exists once it is a proper factor

In the strict multiplicand graph (`x -> z`, `y -> z` iff `x*y = z`,
`1 < x < y`) say that `m` *exists at time `N`* if some added vertex
`z <= N` has `m` as a factor with a distinct cofactor above `1`. **PROVED:**
`m` exists from time `2m` for `m >= 3` and from time `6` for `m = 2`
(`2*2 = 4` is the excluded diagonal), and a number that is itself never a
product (`1`, a prime, a prime square) has no other way to appear. Hence at
time `N` the numbers not yet in existence are exactly the sources (primes
and prime squares) in `(N/2, N]`, with `2` counted until `5`.

At the odd squares this gives `[5,7,9]` at `9` and `[13,17,19,23,25]` at
`25`: `11` is gone because `22 = 2*11` precedes `25`, and `15, 21` are gone
as products. This is the user's grouping, and the correction to the second
reflection (which grouped all odd sources between odd squares) is right:
the relevant primes are those in the top half `(N/2, N]`, which is exactly
Bertrand's territory (at least one prime in `(N/2, N]`, CITED classical).

From `49` on, the previous odd square exceeds the half-line
(`(2k+1)^2/2 < (2k-1)^2` iff `k >= 3`), so if each source is assigned to
the first odd square at which it is not yet in existence, the groups become
"sources strictly between consecutive odd squares, plus the square":
`[29,31,37,41,43,47,49]` (`25` belongs to its own group), then at `81`:
`53,59,61,67,71,73,79` (`49` belongs to `49`), sizes
`3, 5, 7, 7, 9, 10, ...`. The odd sizes `1, 3, 5, 7` are therefore an
accident of the first four groups; for large `k` the size is the number of
primes between `(2k-1)^2` and `(2k+1)^2` plus one, which grows like
`8k/log(4k^2)` heuristically and is at least one by Legendre's conjecture
(OPEN). FINITE-EXACT half-line groups at `9, 25, 49, 81, 121, 169, 225,
289` have sizes `3, 5, 8, 11, 14, 18, 21, 29` (the half-line count keeps the
previous square while it is still not a factor).

## 2. Twenty-five is where the three chains meet

The user's three sums `12+13`, `18+7`, `16+9` are three different edges of
the square-sum graph into `25`, and each is a structural event:

* `12 + 13 = 25` is the merge edge of `Q_13`: vertex `13` has exactly two
  square partners, `3` (square `16`) and `12` (square `25`), and they lie
  in different chains, `{1,3,6,8,10}` and `{4,5,11,12}`; this is the
  `3 -> 2` merge (wave six). Vertex `14` then joins `{2,7,9}` through
  `14 + 2 = 16` and `14 + 11 = 25`. So the square `25` (with `16`) is what
  connects all three chains. In general the odd square `(2j+1)^2` carries
  the consecutive edge `(2j^2+2j, 2j^2+2j+1)`: `4+5 = 9`, `12+13 = 25`,
  `24+25 = 49`, `40+41 = 81`.
* `18 + 7 = 25` is the scar edge: `18` has no other partner until `31`
  because `18 + 18 = 36` is the excluded diagonal, which is the whole
  reason square-sum arrangements fail for `18..22` and `24` (wave six).
* `16 + 9 = 25` is the first edge of any `Q_n` whose two endpoints are both
  squares above `1`: `3^2 + 4^2 = 5^2`, the primitive triple `3-4-5`, the
  root of the Berggren tree (THM-3756, THM-3341). It appears at `n = 16`.

Thus at `n = 25` the square-sum graph has, in one square, the merge of its
finite startup (three chains), its last scar (`18`), and the root of the
Pythagorean tree. That is an exact statement about one integer, and it is
the honest content of "the boundary region 15..25". It has no map to the
Collatz map (the Collatz target `(3n+1)/2` is never a square, wave six) and
none is claimed.

## 3. Thirteen speeds and thirteen vertices

SCOPE. `Q_13` is the unique size with exactly two components, and LRC(14)
concerns thirteen nonzero speeds. No map was found between the component
structure of `Q_n` (reflections `x -> 4-x, 9-x, 16-x`) and the lonely
runner covering problem (residue avoidance at a common time), and none is
claimed. The one shared research move is the repository's own: separate a
finite startup head (here `n <= 14`; there the sporadic tight sets) from
the recurrence class, and never promote a residue-level covering statement
to a global one without an integer height sidecar (META-PATTERNS, first
card). The user's `19 = -1 mod 5` and the window `20 +- 2` remain
numerology.

## 4. What a singleton is, and why the macrocosm is regular

Every singleton met in this session is the small solution set of a
coincidence between two operations, and the corresponding regularity of
the macrocosm is the theorem that the coincidence has no large solutions:

| singleton | coincidence | why it stops |
|---|---|---|
| summand holes `4, 6` | `2+2 = 2*2`, `3+3` are deleted diagonals | only `2S` lacks another distinct representation |
| trivial Collatz cycle `{1}` | `1+1 = 2`, `2^2 - 3 = 1` | Catalan (Mihailescu): `2^K - 3^L = 1` only at `(2,1)` |
| minus cycles `{1}`, `{5,7}` | `2 - 3 = -1`, `2^3 - 3^2 = -1` | Catalan: `3^L - 2^K = 1` only at `(1,1),(2,3)` |
| Bang's `63` | `2^3 + 1 = 3^2` | Zsigmondy: no other base-2 exception |
| multiplicand sources `p^2` | `p*p` is the diagonal | infinitely many, but each is one point |
| square-sum scars `2, 8, 18` | `2m = j^2` with one square in `(2j^2, 4j^2)` | for `j >= 4` that interval holds two squares |
| `a_8 = 0` at level `11` | `a_2^2 = 2p` at `p = 2` | the Frobenius angle `3 pi/4` fixes the period |
| Wieferich rows `1093, 3511` | `2^(p-1) = 1 mod p^2` | no third example below `6.7*10^15` (CITED) |
| the seven-cycle clock `11/7` | `2^11 - 3^7 = -139` divides the carry | Pillai/Baker: `|2^K - 3^L|` grows |

The "rough isomorphism" the user asks for is exactly this: a microcosm is
the finite solution set of an equation `f(x) = g(x)` between two growth
laws, and the macrocosm is regular because the growth laws separate. This
is a principle of classification, not a map between the singletons (SCOPE
for any such map).

It also states the honest limit of this route for Collatz. A nontrivial
cycle needs `(2^K - 3^L) | bB(w)`, which is a coincidence between an
exponential gap and a carry that is itself exponentially large, not a
coincidence at a small argument; Baker and Pillai forbid small gaps, and
the counterexample portrait shows every remaining condition is satisfiable.
Divergence is not a coincidence statement at all. So the singleton
principle explains why the finite heads are finite; it cannot, by itself,
supply the sign-specific order statement that the word-function theorem
(wave six) says a proof must contain.

## Reproduction

The checks are the short computation saved at writing time as
`/tmp/lead_check3.txt`: existence-time groups at the odd squares `9..289`
and the even squares `4, 16, 36`; the twelve pairs summing to `25`; the
consecutive-pair law `(2j^2+2j, 2j^2+2j+1)`; the edges of `13` and `14`;
the first square-square edge `(9,16,25)`. Elementary; rerun in any Python
with sympy.
