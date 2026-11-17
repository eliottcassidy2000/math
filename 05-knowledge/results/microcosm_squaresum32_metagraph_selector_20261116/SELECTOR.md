# A selector that fills every count except 7 and 21: what this programme can state, and its exact reach

Owner: *"Aim for a selector that fills every count except 7 and 21."* A selector is a procedure that, given an odd
target `N`, names a tournament of count `N`. Every selector in the repo is exact up to a bound and conjectural beyond;
this one is assembled from the programme's proved pieces and measured honestly.

## 1. The composite selector

Given odd `N ∉ {7, 21}`:

1. **Small `N`.** If `N ≤ 2881`, take any realiser from the census of orders ≤ 9 (every odd number in `[1, 2881]` other
   than 7 and 21 is a count of some class of order ≤ 9, [MICROCOSM.md](MICROCOSM.md) §2).
2. **Near-join.** If `N = h + 2^k B_t(T)` with `k ≥ 1` for a rooted class `(T, t)` of order ≤ 10, return the near-join
   `T → M_k → z` with `z → t` (THM-5263 (I3)). This is exact for every odd `N ≤ 19,190,825` except the seven near-join
   exceptions 1, 7, 11, 13, 21, 47, 67 (bases ≤ 9 reach 398,129, THM-5265; bases ≤ 10 reach 19,190,825,
   [receipts/nj10.out](receipts/nj10.out)); 11, 13, 47, 67 are handled by step 1.
3. **Closure.** If `N = c·m` with `c` a count and `m` already selected, return the ordinal sum (THM-5264 (E) keeps
   near-joins near-joins; ordinal sums multiply counts).
4. **Seams.** If `N ± 1` is already selected at a class carrying a one-contact seam at centre `N ± 1`, flip that contact
   (THM-5267; every even centre in `[264, 10000]` has one from bases ≤ 9).
5. **Beyond the reach**, no step is exact; the two-level placements of wave 56 and the frozen-fibre bank of the
   concurrent two-cover packet reach specific primes near 10⁶ and 10¹⁰, and the repository's digital bank realises every
   odd number ≤ 10⁸ except 7 and 21 (HYP-9453's status line).

## 2. Exact reach and what fills the holes' place

The selector is exact for every odd `N ≤ 19,190,825` except 7 and 21: step 1 below 2,882, step 2 above it (step 2 misses
only the seven, all below 68). Odd misses of step 2 above its reach: 1,197 odd numbers below 2^25, all in [2^24, 2^25); by dyadic window: none below 2^24, then 1,197 of 8,388,608 in [2^24, 2^25)
([receipts/residues10.out](receipts/residues10.out) gives the residue structure: flat, as for bases ≤ 9).

## 3. What would make it fill everything

Any one of three density statements closes the selector: HYP-9453 (every odd `N` is `h + 2^k B` for some rooted
tournament), HYP-9864 (every even `c ≥ 264` is a seam centre), or THM-4094's lane reduction supplied with realisers for
every prime `p ≥ 23` and every `7p`, `p ≥ 5`. Each is verified far beyond any structural reason known for it; the
metagraph gives the ladder and the images their shape, the quadratic form gives every step an exact law, and the holes
7 and 21 are the two forced jumps. Chan–Pak OPEN.
