# Session-lead reflection II: every hole, source and leaf is a deleted diagonal

**Status: PROVED elementary statements (hand proofs, each checked by the
short computation reproduced at the end); FINITE-EXACT counts; CITED
classical inputs; SCOPE where no map exists. No agents and no independent
audit: a reasoning note by the session lead under a usage limit, provenance
for the exact claims it isolates, not audited canon.** Collatz and
Legendre's conjecture remain OPEN. Companion to the
[first reflection](collatz_mod6_20260922_lead_reflection.md) and the
[synthesis](collatz_mod6_20260917_synthesis.md).

## 1. The user's groups are the odd sources of the strict multiplicand graph

Take the multiplicand graph with the same strictness as the summand graph:
`x -> z` and `y -> z` iff `x*y = z` with `1 < x < y`. A number is *brought
into existence* when some pair of smaller distinct factors above `1`
multiplies to it. **PROVED:** the numbers never brought into existence are
exactly `1`, the primes, and the prime squares `p^2 = p*p`, because a
prime square's only factorization above `1` is its diagonal, while every
other composite `n` has a factor `x` with `1 < x < n/x`. Reproduced to
`60`: sources `1,2,3,4,5,7,9,11,13,17,19,23,25,29,...,49,53,59`.

The odd sources between consecutive odd squares are therefore
`[3,5,7,9]`, `[11,13,17,19,23,25]`, `[29,31,37,41,43,47,49]`, ... : the
user's groups with `11` restored. Each odd square `p^2` closes its group
because it is created only by its own diagonal `p*p`, which the strict
rule forbids. This is the exact multiplicative twin of the summand holes
`{1,4,6} = {1} cup {2+2, 3+3}` of the first reflection: in both graphs the
exceptional set is `{1}` together with excluded diagonals, finitely many
on the additive side (only `4` and `6` lack another distinct
representation) and infinitely many on the multiplicative side (every
`p^2`).

The "even number of primes" pattern is an artifact of omitting `11`: the
number of primes strictly between `(2k-1)^2` and `(2k+1)^2` is
`4,5,6,7,8,9,9,13,11,13` for `k = 1..10` (FINITE-EXACT), growing like
`8k/log(4k^2)` heuristically; that at least one prime lies between
consecutive squares is Legendre's conjecture (OPEN).

## 2. The multiplicand graph inside the summand graph plus its complement

Let `S` be the strict summand arc set (`x -> z` iff `x < z` and
`z != 2x`), `D = {u -> 2u : u >= 1}` the doubling forest, which is the
complement of `S` inside the ascending arcs (inherited summand note,
section 1), and `M` the strict multiplicand arc set above.

**PROVED.** `M` is contained in `S cup D`, and
`M cap D = D minus {1 -> 2, 2 -> 4}`. Hence

    M = (M cap S)  disjoint-union  (D minus {1 -> 2, 2 -> 4}).

Proof. An arc of `M` from the smaller factor is `x -> xy` with
`1 < x < y`, so `y >= 3` and `xy >= 3x != 2x`: it is a summand arc. An arc
from the larger factor is `y -> xy`; if `x >= 3` it is a summand arc for the
same reason, and if `x = 2` it is the doubling arc `y -> 2y` with `y >= 3`,
which is not a summand arc. Conversely every doubling arc `y -> 2y` with
`y >= 3` is the larger-factor arc of `2*y`. The two doubling arcs missing
from `M` are `1 -> 2` (the factor `1` is excluded) and `2 -> 4` (the
diagonal `2*2`). Reproduced on all arcs with target `<= 200`.

So the user's picture is exact: the multiplicand graph is a subgraph of the
summand graph together with the summand graph's complement, and the only
two arcs of that complement it does not use are the two at the coincidence
`1+1 = 2` and `2+2 = 2*2 = 4`. The number `2` is at once the doubling
multiplier that generates the complement and the summand diagonal that
creates the hole `4`.

## 3. The ray tableau is the odd/even by nonsquare/square table

Every `n` lies on the doubling ray of its odd part and on the squaring ray
`b, b^2, b^4, ...` of its nonsquare base `b` (`n = b^(2^k)`, inherited
summand note, section 2). A vertex contributes a *new* doubling ray iff it
is odd (otherwise its ray sits inside the ray of `n/2`) and a *new*
squaring ray iff it is not a square (otherwise inside the ray of
`sqrt n`). This is the whole content of the tableau:

| class | new doubling ray | new squaring ray | examples |
|---|---|---|---|
| odd nonsquare | yes | yes | `3, 5, 7, 11` |
| odd square | yes | no | `9 -> 18, 36, ...`; `25` |
| even nonsquare | no | yes | `6 -> 36`, `8 -> 64`, `10 -> 100` |
| even square | no | no | `4, 16, 36, 64, 100` |

Below `100` the classes have sizes `45, 5, 45, 5`. The vertex `1` is the
unique fixed point of squaring and doubles normally (`1 -> 2 -> 4 -> ...`),
which is why it is the root of the doubling forest but isolated in the
squaring forest. The remark "`4` adds nothing, `6, 8, 10` add only their
squares, `9` adds only its doubles, `25` caps the block" is the third
table row by row; the blocks of five are not a structural unit.

## 4. The square-sum boundary is the same diagonal

In the square-sum graph `Q_n` the window of impossibility `18..22` (and
`24`) is caused by vertex `18` having the single neighbour `7` until `31`
appears, because `18 + 18 = 36` is a square that the distinct-summand rule
excludes (wave six, PROVED). The earlier scars are the same phenomenon:
`2` is born isolated because `2 + 2 = 4`, and `8` is a leaf until `17`
because `8 + 8 = 16`; the family `m = 2j^2` has this effect exactly for
`j <= 3` (wave six). So the three finite heads of this session's objects
are one mechanism: the summand holes `4 = 2+2`, `6 = 3+3`; the trivial
Collatz arrow `1 -> 1+1`; the multiplicand sources `p*p`; and the
square-sum leaves `2+2`, `8+8`, `18+18`. Each is a deleted diagonal of the
swap-fixed operation fibre of THM-2422.

The remaining numerology of the message is typed as SCOPE: the window
`18..22 = 20 +- 2` and the preperiodic set `{-2,...,2}` of `x^2 - 2` share
an interval and nothing else; `13` being the only size with exactly two
components and `19 = -1 mod 5` are true facts without a map.

## 5. The one candidate research card

"Before reading a hole, a source, or a leaf in an incrementally built
operation graph, test whether it is a deleted diagonal `x o x`." Evidence
from four distinct threads above; counterindication: the doubling forest
itself is *made* of deleted diagonals, so the card explains exceptions,
not the generic structure. Not promoted here; recorded for a future
session with audit capacity.

## Reproduction

The checks are the short computation saved at writing time as
`/tmp/lead_check2.txt`: prime counts between odd squares for `k <= 10`;
strict multiplicand sources to `60` against `{1} cup P cup P^2`; the arc
identity `M = (M cap S) u (D minus {(1,2),(2,4)})` on all targets `<= 200`;
the four ray classes to `100`. They are elementary and rerun in any Python
with sympy.
