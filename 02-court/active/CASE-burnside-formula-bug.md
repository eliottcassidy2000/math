# Court Case: Burnside Sign-Flip Formula Bug

**Filed by**: kind-pasteur-2026-03-22-S20bq
**Status**: OPEN
**Priority**: Medium

## The Issue

The Burnside formula for counting tournament iso classes (A000568) requires tracking "sign flips" in pair orbits under permutations. My implementation in `burnside_fix_s20bq.py` gives WRONG values at n >= 4:

| n | Brute force (correct) | Formula (wrong) |
|---|----------------------|----------------|
| 3 | 2 | 2 |
| 4 | 4 | **5** |
| 5 | 12 | **15** |
| 6 | 56 | **66** |

## The Correct Approach

The brute-force Burnside (summing Fix(sigma) over all sigma in S_n, where Fix counts tournaments with sigma as automorphism) gives A000568 exactly.

The FORMULA approach (track sign flips in pair orbits) is conceptually correct:
- For each pair-orbit of sigma, count how many times sigma reverses the canonical ordering
- If total reversals is odd: the orbit contributes 0 to Fix
- If even: the orbit contributes 2 (free choice of direction)
- Fix(sigma) = product over orbits

But the implementation has a bug in the orbit-tracing code -- specifically in how orbits are closed and sign flips are accumulated.

## What Needs To Be Done

1. Fix the orbit-tracing code in `count_tournament_fixed()`
2. Verify the formula matches brute force at n=3,4,5
3. Use the corrected formula to compute A000568 at n=7+ (where brute force is infeasible)
4. Document the correct formula as a theorem

## The Known Correct Formula

From OEIS A000568 (Davis 1953):
a(n) = (1/n!) * sum over permutations sigma of S_n of 2^{c(sigma)}
where c(sigma) is the number of orbits of sigma on ORDERED pairs (i,j) with i != j
such that the orbit is "even" (doesn't reverse direction).

This is equivalent to what I'm trying to compute, just needs the correct implementation.
