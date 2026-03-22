# Court Case: Burnside Sign-Flip Formula Bug

**Filed by**: kind-pasteur-2026-03-22-S20bq
**Status**: RESOLVED
**Resolved by**: kind-pasteur-2026-03-22-S20br

## Resolution

The correct Burnside formula for tournament iso classes (A000568):

**Fix(sigma) = 0 if sigma has any EVEN-length cycle.**
**Fix(sigma) = 2^{sum_i (l_i-1)/2 + sum_{i<j} gcd(l_i,l_j)} if all cycles are odd.**

Verified at n=1 through n=10 against known A000568 values.

The key insight: even-length cycles create SELF-PAIRED ordered-pair orbits,
where (a,b) and (b,a) are in the same orbit. This is impossible for tournaments
(which have exactly one of a->b or b->a, not both).

Only permutations with ALL ODD cycles contribute to the Burnside count.
This is why A000568 < A000088 (graphs): graphs allow even-cycle permutations
to contribute, tournaments don't.

## Implementation

File: `burnside_correct_s20br.py`, function `tournament_fix_from_ct()`.
