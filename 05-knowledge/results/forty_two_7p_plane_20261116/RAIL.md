# An induced 41/43 to 7-times twin rail

Two strong order-nine tournaments, rows in outgoing-mask order
and with vertex 8 marked, are:

    256 1  67  23 71 287 427 319 30
    256 1 259 279 39  79  31 127 242

Their (H,E,O,g) profiles are respectively (287,217,240,23)
and (301,217,222,5). The first contains the vertices
{0,1,2,3,5,6,8} with exactly 41 Hamiltonian paths; the second
contains the vertices {0,1,2,3,4,5,8} with exactly 43. Hence
these are induced, zero-old-reversal, two-vertex scale-by-seven
lifts specifically for these two primes.

Both roots are normalized (O=E+g) and have the same E=217.
A new near-sink loses to all old nonroot vertices and beats
the root. On normalized inputs this changes (H,E) to
(H+2E,2E) and preserves normalization. Appending k such
vertices to each witness therefore gives strong tournaments
with counts

    7*(62*2^k - 21),  7*(62*2^k - 19),  k >= 0.

Both are strong at k=0 by literal verification and stay strong
under the specified mixed attachment. At k=2 the factors
are the twin primes 227 and 229. No claim is made that the
factors are always prime or that primes of this form are
infinite. For example the second factor is divisible by
3 for every odd k. The shared E suffices for this operation,
not arbitrary ordinary-sink insertions: g differs.
The standalone check_rail.py confirms literal orientations,
strongness, induced counts, full weighted profiles and
two successive literal extensions by subset dynamic programming.
