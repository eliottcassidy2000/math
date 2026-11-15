# Exact placement through marked inverses: graceful paths, square-sum, and Berggren descent

The session produced exact **target-first** decoders for new physical
subfamilies of tournament counts, and identified what prevents them
from covering all odd targets. No global Chan–Pak proof is claimed.

## Actual exact constructions

- **A ternary tree of literal +2 windows.** The square-sum edge
  absorbers at fixed aperture 18 have disjoint inverse intervals.
  Three letters give a rooted ternary tree; keeping all four gives a
  quaternary tree. On each fixed-sum depth band the tournament reader
  is injective, including across orientation and center choice.
  A prescribed count gives at most one path back to the seed.
  [Physical tree](SQUARE_TREE.md) · [discriminant inverse](DISCRIMINANT.md).
- **Variable scale without losing the inverse.** Two letters of
  aperture 8 and one of aperture 13 have disjoint *ratio* intervals,
  not just disjoint heights. This gives a denser ternary family and
  an exact list decoder with O(log N) candidate sums from the target
  (proof does not assume scalar injectivity across sums).
  [Mixed apertures](MIXED_APERTURES.md).
- **A prescribed gap and absolute endpoint.** A scaled swap of the
  two addresses cancels the bilinear and leading linear terms. Every
  positive multiple of 30 is the exact common-coordinate gap of two
  +2 windows. Fixing the gap and absolute endpoint leaves at most one
  candidate, selected by an integer square root and checked exactly.
  At fixed gap the attained endpoints are still sparse.
  [Scaled-swap inverse](SCALED_SWITCH.md).

## What does and does not transfer

**Graceful order matters, not just the set of differences.** Walecki's
ordered edge lengths let even mutually conflicting neighbors absorb a
whole interval by choosing matchings. Fibonacci weights make the
interval exponential and give the canonical finite Zeckendorf
decoder. The contract is exact additivity on legal matchings, not on
all subsets. Disjoint one-apex interval controls fail even this
weaker contract: their mixed difference is a nonzero integer, and
the touching exception cannot supply four consecutive odds.
At the polynomial clock's error scale, the strong-tournament vertex
bound also rules out implementing each switch with three private
vertices if one array must fill the whole error interval.
[Graceful decoder and obstruction](GRACEFUL.md) ·
[Fibonacci repair contract](FIBONACCI_REPAIR.md).

**Berggren's hypothesis is a physical inverse domain.** The primitive
triple tree has a unique smaller parent because its three positive
inverse cones are disjoint and exhaustive. This transfers to the
square-sum host labels using height or ratio. It does not transfer
to the natural Hamiltonian/two-slot norm profile: even a single
Berggren edge cannot join two such physical endpoints, and no
bounded-length macro alphabet fixes this for an infinite chain.
A fixed-count parabolic ray still has an exact dyadic membership
test, and primitive diagonal profiles fill the allowed joint
dyadic/triadic residues. Prime diagonal targets above 3 must be
11 mod 24. [Physical norm domain and decoder](BERGGREN.md).

**An isomorphism does not imply dense scalar coverage.** Even giving
this reader *all* square-sum address pairs from every seed and
aperture produces at most O(X^(3/4)) counts ≤X. Changing how the
reader uses the host labels is necessary; adding more branches
within the square-sum domain alone cannot close Chan–Pak.

## Next precise obligation

Build a marked family in which the prescribed absolute target
selects both the baseline and an exact compatible correction.
One concrete target is an exponential matching decoder with
interleaved or shared switches that preserve the marked interface
and the +2 options. The separated one-apex version is now ruled
out, not merely untested. Another target is a reader that retains
the mixed-aperture inverse but no longer requires the *readout*
sum to be a square.

## Replay

Run `python verify_square.py`, `python verify_discriminant.py`,
`python verify_mixed.py`, `python verify_scaled.py`,
`python verify_graceful.py`, `python verify_fibonacci.py`,
and `python verify_berggren.py` in this directory.
