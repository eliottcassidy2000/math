# The chessboard lens: rings are loneliness, the knight is the tight rider, and every tour must stall

*mac-mini, 2026-10-06, chessboard-weave session. A reflection on the owner's 8x8 prompt (rings 4, 12, 20, 28;
diagonal scaffolds {1,3,5,7,7,5,3,1} by {2,4,6,8,6,4,2}; one step or slide until obstruction; the knight's tour
weaving the two). Provenance, not truth: the statements live in
[`chessboard_weave_20261006.md`](../05-knowledge/results/chessboard_weave_20261006.md), THM-4550, THM-4551 and HYP-9210.*

## What the picture turned out to be

The owner's two structures are two `Z/2`-flavoured gradings of one board. The rings are Chebyshev shells, and a
shell is exactly a level set of the lonely-runner function `min(||x||, ||y||)` on the unit cell; the bishop, the
piece that slides diagonally until it is obstructed, measures that level exactly (`reach = 6 + 16 lambda`). The
rook sees nothing of it. So "horizontal/vertical rings" and "diagonal scaffolds" are welded by the bishop, and
the lonely runner is the theory of how deep a sliding piece gets.

The knight enters three times, each time as the extremal case:

* as a **rider** it is the tight two-speed lonely runner (`delta = 1/3`), the one direction whose line never
  reaches the centre of the 8x8 cell;
* as a **leap** on the exponent board it is the trivial Collatz cycle `1 -> 4 -> 2 -> 1` and the average
  Syracuse step;
* as a **tour** it must stall: the weave law says every closed tour makes as many moves inside the outer ring
  as inside the middle annulus, and at least one.

And the monotile paper joins through the same door: its boundary substitution replaces every rook step by the
centre-to-centre flight of a knight or a zebra. The flights that would hit a lattice corner are exactly the
colour-preserving leapers, which are exactly the perfectly lonely riders.

## The mechanism worth keeping

The weave law is a two-line argument with a general shape. A cyclic process whose every step flips one
grading (colour) and *usually* flips a second (the ring side) must sometimes fail to flip the second, unless
the two gradings are locked together on the whole space; and when the two sides have equal size, the failures
inside one side equal the failures inside the other. The 8x8 board is special only because its ring sizes
form an arithmetic progression, `4 + 28 = 12 + 20`, and because the only knight moves inside rings 0 and 3 are
the eight corner moves of the outer ring. The same argument, in Hall language (`|N(X)| = |X|` after deleting
the within-ring moves), is the obstruction genus of the tournament blocking conjecture HYP-9168: kill every
Hamiltonian object by starving a set.

## Where the analogy stops

Nothing here touches LRC(14), Collatz, or `H >= disc`. The two-dimensional rider facts are classical or
elementary; their content is vocabulary with exact meaning. The first-lonely-time object (HYP-9210) is new to the
repository and gives the cleanest local/global contrast so far: lonely-runner first passage is uniformly below
`1/2` (though its supremum creeps up to `1/2` through the progressions with the 2 removed), whereas Collatz first
descent is unbounded. The cat map destroys the rings (no shell survives a hyperbolic map), which is the same
reason the monotile's boundary is fractal.

## A prompt for the next session

Look for weave laws where the owner's prompts keep pointing: a cyclic or Hamiltonian object, two gradings, and a
size identity. The tournament version (HYP-9168's Hall genus) and the Collatz version (parity against the plus
and minus sheets) are the two untested cells.
