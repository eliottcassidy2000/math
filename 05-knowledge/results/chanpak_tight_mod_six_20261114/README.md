# Absolute height windows and the arithmetic of 2 and 3

**Partial exact selection, not Chan–Pak completeness.** This
checkpoint separates two kinds of progress: one construction does
force absolute equality in explicit infinite windows; another shows
why joint local solvability at 2 and 3 can still be far from sufficient.

1. **[Complementary height windows](TIGHT_COMPLEMENT.md).**
   Pair two literal binary tails so the first-order high-prefix
   terms cancel as integers. At each center B(u,k), the mod-six
   selector keeps the error ≤4<6 and realizes every odd target in
   [B−2,B+2]. Examples are 61,63,65 and 573,575,577. Lower
   dyadic windows have a mod-four selector with error ≤2.
   A complete absolute inverse for the chart returns actual
   adjacency masks. Farther from the center a two-bit phase
   correction bounds error by the square root of the residual;
   for general odd targets there is a same-phase O(sqrt N)
   approximation, not an error less than six.
2. **[Exact inverse-Collatz transport](COLLATZ_RAYS.md).** The
   inverse exponent's parity determines the ternary branch,
   and taking the next admissible exponent literally appends 01
   to a tournament contact word. Two fixed physical rays
   together hit every odd residue modulo every 2^b·3^a.
   Nevertheless for the single prescribed attainable target 63,
   even the least scalar match on those rays grows exponentially
   in the dyadic modulus. We compute the exact first matching
   parameter at all b≥10, not an averaged obstruction.
3. **[Ternary runs and two narrow height decoders](RUN_HEIGHT.md).**
   H mod 3 is an exact function of the number of adjacent drops
   and the root's phase. This yields a sharp fixed-run height
   ceiling, and there are explicit forced collars at both
   ends of the permitted width range. At fixed middle width
   the physical remainder is smaller than the modulus; both the
   quotient and residue must still match on the same word.
4. **[Mod 3 across arbitrary homogeneous readers](READER_MOD3.md).**
   All such readers collapse to an ordinary count multiplied
   by one of two test values, selected by the source parity.
   The mechanism is special to 2 and 3, where every nonzero
   residue has square one. As with inverse Collatz, prime-sum
   sign filters, and primitive right triangles, the local
   modes have real structural content. They do not supply
   an absolute witness. In particular even mandatory divisibility
   by three does not give a physical operation dividing the
   ordinary count by three.

The constructive problem left over is to place the height windows
at a prescribed target, rather than sharpen congruence on an
unbounded trajectory. The exact current windows and collar charts
are sparse. Full coverage needs new movable source or reader
families that supply an exact quotient along with the small remainder.

Run `python3 replay_balanced_windows.py 577` for an actual
literal witness, or give another target. Run
`python3 replay_collatz_runs.py` for the independent checks of
the CRT, the fixed-target obstruction, and both word-height formulas.
