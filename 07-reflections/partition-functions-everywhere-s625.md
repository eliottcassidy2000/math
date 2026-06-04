# Partition functions everywhere (S625)

Hamkins assigns ordinals to positions in infinite games — infinite chess, infinite Go — and the assignment is a
recursion: a position you have already won has value zero, and otherwise the value is built from the values of the
positions you can move to, taking a minimum when it is your move and a supremum-plus-one when it is the opponent's.
The user told me to see the n+2 recursion in this and to come to see partition functions everywhere, and once I
wrote the game-value recursion next to the recursion I have been using all month, they were the same recursion.

The independence polynomial — my tournament `H`, my covering-depth generating function, the unit-distance count —
satisfies deletion-contraction: the value at a graph is the value with a vertex deleted, plus `x` times the value
with that vertex and its neighbors deleted. Condition on whether the vertex is in the independent set or not. That is
a sum over the two choices, weighted. Hamkins's game value is the same recursion with the sum replaced by a minimum
and the weight replaced by plus-one. A partition function is a sum-over-choices; a game value is a min-over-choices;
they are the same object read in two different semirings, the ordinary one and the tropical one. The ordinal game
value is the third reading, in the semiring of ordinals, and there the recursion has to be well-founded for the value
to exist at all. For a path graph the ordinary reading gives the Jacobsthal numbers, exponential; the tropical
reading gives a rank, linear. Same recursion, two growth rates, because the semiring decides whether you are counting
or ranking.

The reading that pays off is Collatz. Make it a game: you sit at a number, your move is the shortcut step, you win
when you reach one. The game value is the well-founded rank of your position, and the Collatz conjecture is exactly
the statement that every position has one — that the game is open, that there is no infinite play, that the recursion
terminates. And I already computed that rank, two sessions ago, without knowing that is what it was: the altitude
tower, the number of steps at floor one and the number of epochs at floor two. The iterated logarithms were the
ordinal. A trajectory that diverges or cycles is a position with no ordinal value, an infinite play in the open game.
Collatz is Hamkins's well-foundedness question for one specific game.

And the forbidden values fall out of the same picture, sharper than before. Twenty-one is a Jacobsthal number — it is
the independence polynomial of a six-vertex path at `x` equal to two — so it is a perfectly good partition-function
value. But twenty-one is forbidden as a tournament `H`. The resolution is that not every conflict graph is the
conflict graph of a tournament's three-cycles; the path that realizes twenty-one is not one a tournament can build.
The partition function ranges over all conflict graphs; the tournament sees only a sub-family; the gaps are the
shadow of that constraint. So the master object is the partition function, the master recursion is deletion-
contraction, and the four problems are the same recursion wearing the masks the user has been naming for a month:
the count, the gap, the rank, the game. I did not finish formalizing the recursion this session — the free-vertex
case is clean and the rest is bookkeeping — but I finally see the single object underneath, and it is, everywhere, a
partition function.
