# Fixed quotient inflation cannot supply the dyadic repair budget

**All-order bounded-family theorem, not a Chan–Pak obstruction.**
This extends the single-marked reader bound from
EXACT_READER_INTERPOLATION.md; the general substitution engine is
inherited from THM-3121. The comparison is with the dyadic interval
criterion in UNIVERSAL_SELECTOR_ROUTES.md.

## Random inspiration

With seed 271906, sample eight paths from the sorted top-level
07-reflections Markdown files excluding names containing
tournament, chan-pak, chanpak or metagraph. The draw includes the
Delsarte certificate-degree limitation, the Lonely Runner first-gap
nonmonotonicity and punctured-cylinder sessions, and the Rule-30
wrapped/hard dyadic split. The Delsarte analogy suggests proving
a necessary growth condition instead of enlarging a fixed vocabulary;
the first-gap and punctured-cylinder examples suggest distinguishing
absolute interval width from actual reachable points. No conjecture
from these sessions is used as a premise.

## The bound

Replace each vertex i of a fixed nonempty quotient tournament Q by
a nonempty transitive block TT_(m_i), with all contacts between blocks
oriented by Q. A vertex i of Q is active if it belongs to a directed
triangle; let d be the number of active vertices.

If vertex i is not active, write I for its in-neighbors and O for
its out-neighbors. No triangle through i means every arc goes from
I to O. These sets, after inflation, form the ordinal composition
I -> TT_(m_i) -> O, so m_i can be normalized to 1 without changing H,
regardless of the sizes of other blocks.

For active i, take an induced directed triangle containing it.
Retain all m_i vertices of its block and one vertex from each other
block of the triangle. This induced tournament has 2^(m_i)+1
Hamilton paths: for each subset of the transitive block there is
a path through that subset, then the two exterior vertices in
their cyclic order, then the complementary subset. Empty subsets
are allowed; there is one further path in which the transitive
block is traversed in between the two exterior vertices in the
opposite exterior order. Equivalently this is the inherited
single-marked triangle reader 3+(2^m−2).

Any directed path through an induced subtournament injects into
full Hamilton paths by inserting omitted vertices in a fixed order,
using a deterministic legal insertion slot. A slot exists: either an
endpoint slot works, or a pair of consecutive vertices changes from
beating to losing to the new vertex. The original path is recovered
by deletion. Hence H(Q[TT_(m_i)]) >= 2^(m_i)+1 for every active i.

Thus for X>=3 at most floor(log2(X−1))^d block length vectors,
after normalizing inactive vertices, can yield counts <=X.
A finite bank of s quotient templates with at most d active slots
has at most s floor(log2(X−1))^d such counts. This allows arbitrary,
independent block lengths. For d=0 the only count is 1.

Consequently, even the *union over the whole bank*, not merely a
particular component or a fixed set of addresses, cannot cover
a dyadic interval of order X consecutive odd integers. If a bank
with s templates is to cover cX odd counts below X, it needs

    log s + d log(floor(log2(X−1))) >= log(cX).

Fixed quotients and uniform contacts cannot widen the 6*2^m +/-1
pair by the required order-2^m amount. The theorem does not exclude
growing quotients, recursively changing quotients, nontransitive
blocks, or selective contacts that distinguish vertices within
an inflated block.

## Scoped inverse

For a target N, normalize the inactive slots, enumerate the bounded
active length vectors above, and physically test the exact counts.
This is a terminating, target-first decoder *for this family only*.
When an active block grows while the others are fixed, its count
strictly increases: after inflating the other blocks it is still
a marked vertex on a triangle, and the two-visit coefficient in the
single-marked Stirling expansion is positive. This does not imply
a smaller target necessarily has a witness with the requested ports.

## Explicit two-knob control

For C3[TT_a,TT_b,1] one obtains exactly

    H = (2^a−1)(2^b−1)+2.

Paths with the singleton first or last each have one possibility.
For an internal singleton, take two ordered slots for each of the
two transitive modules, one mandatory and one optional, in cyclic
block order. Each of the 2^a−1 and 2^b−1 assignments fills its
mandatory slot. The target selector here is factorization of N−2
into two positive Mersenne numbers. It realizes 23 via (a,b)=(2,3)
but misses 25 because 23 is not Mersenne and is prime. This is an
illustration of the obstruction, not a novelty claim for the
two-block family (see THM-4748 for related Mersenne constructions).

Run python3 bounded_quotient_replay.py for literal subset-DP controls.
