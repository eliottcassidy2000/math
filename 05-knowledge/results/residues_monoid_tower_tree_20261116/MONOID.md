# The multiplicative monoid of counts: hidden generators from 35 and 49, the first 32 odd numbers, and the address coordinate

Session klein-2026-11-16-monoid. Script `scripts/indecomposables.py`, receipt `indecomposables.out`.

## 1. Indecomposables (THM-5258 (B), PROVED)

In the monoid `M = (odd numbers ∖ {7, 21}, ×)` the indecomposables are exactly

- the primes `p ≠ 7`;
- the numbers `7p` for primes `p ≥ 5` — **the hidden generators**, beginning `35, 49, 77, 91, 119, 133, 161, …`;
- `49, 63, 343`.

Verified below `10^5` (11,266 indecomposables, 1,676 of them hidden; the prediction P-5 **confirmed** exactly). Since
`R`, the set of counts, is a multiplicative monoid (ordinal sums) containing every odd number ≤ 2^32 except 7, 21, the
structure above is the structure of `R` below `2^32`, and of all of `R` if Chan–Pak holds. HYP-+2880 noticed that
"35 = 5·7, 49 = 7² are atoms"; THM-4094 reduced Chan–Pak to the two lanes `p ≠ 7` and `7p`; here the complete list of
atoms follows from the two holes alone. The owner's `{35, 49}` are the two smallest hidden generators; the set has the
cardinality of the primes, hence of the natural numbers; and `63 = 7·9`, the Zsigmondy exception, is the one hidden
generator that is a multiple of 3.

## 2. The first 32 odd numbers, 3 + 1 + 28, and the address coordinate

Write `t = (H − 1)/2` (the address of hadamard's §5; for a one-apex tournament `t` is the order of the complete
conflict graph). The first 32 odd numbers are the addresses `t = 0, …, 31`; the holes are `t = 3` (7, the fourth odd
number) and `t = 10` (21, the eleventh); the block ends at `t = 31` (63). Taking the fourth out splits the block
`3 | 1 | 28`, as the owner says. Multiplication by 7 is `t ↦ 7t + 3` (`7(2t + 1) = 2(7t + 3) + 1`), so:

| source `t` | `H` | `7t + 3` | `7H` | status |
|---|---|---|---|---|
| 0 | 1 | 3 | 7 | hole |
| 1 | 3 | 10 | 21 | hole |
| 2 | 5 | 17 | 35 | hidden generator |
| 3 | 7 | 24 | 49 | hidden generator (image of a hole) |
| 4 | 9 | 31 | 63 | hidden generator (last of the first 32) |
| `(p − 1)/2`, `p` prime ≥ 5 | `p` | `7(p − 1)/2 + 3` | `7p` | hidden generator |
| 24 | 49 | 171 | 343 | hidden generator (second image of the hole 7) |
| 31 | 63 | 220 | 441 = 63·7 | decomposable (`441 = 49·9`) |

The hole at `t = 3` maps to `t = 24` (49) and then to `t = 171` (343), and stops: `7^4 = 49·49` decomposes. The hole
at `t = 10` maps to `t = 73` (147 = 49·3), decomposable. So the holes are the images of the two trivial counts 1, 3 under
`×7`, the hidden generators are the images of the primes ≥ 5 and of `7, 9, 49`, and nothing else is hidden. The
"fractal" of the owner: at scale `2^k` the hidden generators occupy the addresses `7t + 3` with `t` a prime address
below `2^k/7` — their number is `π(2^(k+1)/7)`, growing with the primes (receipt: 3, 6, 10, 21, 34, 61, 106, 192 hidden
among the first `2^5, …, 2^12` odd numbers).

## 3. Microcosm and the field with one element (typed)

- The owner's "prior micro/macrocosm sessions" is hadamard's §5 ("the 3 + 7 = 10 microcosm"): in the address
  coordinate the column `7, 21, 63` is `3, 10, 31` with `3 + 7 = 10`, `10 + 21 = 31`, i.e. `t ↦ 3t + 1 = t + H`
  (multiplication by 3 is `t ↦ 3t + 1`). Combined with §2: `×3` is `t ↦ 3t + 1`, `×7` is `t ↦ 7t + 3`, and in general
  `×(2k + 1)` is `t ↦ (2k + 1)t + k` — the address map of the whole monoid. The microcosm `3 + 7 = 10` is the
  statement that the two holes are related by `×3` (`t: 3 ↦ 10`), and the block boundary 31 is `×3` once more.
- "The field with one element": the repository's reading (CONCEPT-MAP row 29, HYP-1539 in the historical index) is
  that a tournament "over `F_1`" is a linear order, with `H = 1`. In THM-5257's terms that is the depth-first build of
  a tree; in THM-5258's terms it is the unit of the monoid. The two holes are not `F_1`-phenomena: they sit at the
  addresses `3` and `10`, not at the unit. ANALOGY only.

## 4. Converse and the Platonic count (CONVERSE.md)

The owner's "33 first appears at order six as one self-converse class plus two converse pairs" is confirmed and
placed in the full multiplicity table there; it is the same `1 + 2·2` as the five solids.
