# A scalar Collatz recurrence cannot be a uniform homogeneous lift

The earlier [four-visit obstruction](../session321_selective_returns_and_fano_wiring_20261113/FOUR_VISIT.md)
compares particular equal-count pairs with ordered cover coordinates.
The following theorem concerns every candidate scalar rule that is
supposed to act on arbitrary witnesses, without a marked sidecar.

## Classification, including several copies of the same witness

Fix a quotient tournament, mark k>=1 of its vertices, and fill each
marked block with a copy of the same nonempty tournament U. Fill the
other blocks with fixed nonempty tournaments. Write h=H(U). The
resulting count is a function of h alone, for all U, if and only if
each marked vertex belongs to a singleton SCC of the quotient.
In that case the function is exactly C h^k, for a fixed positive
integer C. In particular the only uniform one-copy scalar laws are
multiplications.

To prove this, compress maximal runs within blocks in a Hamilton
path. Each quotient itinerary contributes a product of ordered
nonempty cover counts of its factors, with nonnegative coefficients.
Every block is visited at least once. If a marked vertex is in a
singleton SCC of the quotient, no itinerary can re-enter it:
the intervening part would give a directed closed walk through that
vertex and another vertex. Thus if all marked vertices are such,
every contribution uses exactly one H(U) for each copy, proving
the monomial formula.

Conversely a marked vertex in a nonsingleton SCC lies in a directed
triangle: in its strong component there is an arc from an out-neighbor
to an in-neighbor (otherwise no path from its out-set to its in-set).
Start with that triangle traversed as a word whose first and last
vertices are the mark. Insert every other quotient vertex in some
legal gap or at an end. Such an insertion exists in any directed
word when the inserted vertex is absent: if it cannot prepend or
append, the first vertex beats it and it beats the last, so somewhere
an adjacent pair crosses from one status to the other. This yields
an itinerary with two runs at the mark and one elsewhere.

Now replace U by U→TT_m. Its H is still h. For every j>=2, its
F_j is unbounded as m grows: put one Hamilton path of U in the
first slot and distribute the transitive tail among all j slots,
making the other slots nonempty. For example reserve one tail
vertex for each other slot and freely assign the remaining vertices.
Thus the term of the displayed two-run itinerary grows without
bound, and no other term can cancel it. The output cannot be a
function of h.

A choice of context depending only on the residue of h does not
escape this: the same padding keeps h and that residue fixed.
The theorem does not prohibit scalar formulas on a restricted
marked witness class, choosing a context from the full witness,
or using a selective cut inside an old block. It explains why
native versions of 2h±1 or 3h±1 must carry more information than h
if they use these operations. Starting from finitely many seeds,
the monomial laws also cannot create a new prime factor except
from the fixed constants of the contexts.

## Which kind of padding could carry a machine?

For fixed target h, the entire strong witness set is finite, by the
[compressed bound](LOGIC_AND_BOUNDS.md). Thus reachability under
any decidable move rule on those strong equal-count witnesses is
decidable by finite search. The same is true for arbitrary
tournaments if moves are specified on their singleton-free cores,
or if a bound on Q(2) is imposed: Q_T(2)>=2^|T|.

Without either restriction, there is a genuine infinite coordinate:
ordinal singleton padding is invisible to H. A rule can read and
rewrite lengths of such runs as machine registers. But this uses
the padding as extra computational state, not the Hamilton
path-count predicate. The complementary uniform hardness
construction explains exactly how putting time in that hidden
coordinate changes the logical problem.
