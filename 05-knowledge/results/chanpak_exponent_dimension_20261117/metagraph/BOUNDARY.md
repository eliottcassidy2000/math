# A degree boundary gives an exact color constraint on blue edges

This is the bookkeeping rule recovered by the parallel research pass,
checked directly here. The path-frame complement and its closed/open
layer counts are inherited (THM-643, THM-4729, THM-4731 and the
[black-detour dictionary](../../black_detours_closed_blue_lines_klein_dictionary_20261020.md)).
Let the s-to-t path be the frame of T, d the outdegrees, and T' flip
all arcs except the path edges. Relative to full converse, the
interior vertices keep one incoming and one outgoing edge, and the
two ends carry exactly the two imbalances:

    d'(v) = n-1-d(v) + [v=s] - [v=t].

Consequently the centered second degree moment, with center (n-1)/2,
changes by 2-2(d(s)-d(t)). The number of transitive triples is
sum_v binom(d(v),2); the total degree is fixed. So if c3 counts directed
triangles, then

    c3(T') - c3(T) = d(s)-d(t)-1.

This law transports to the actual line quotient: c3 is invariant
under both relabeling and converse. For a blue frame, reversal
forces d(s)+d(t)=n-1. At odd orders every blue edge toggles c3 parity,
so that the blue graph on the merged classes is bipartite and
loopless. At even orders a blue loop can occur only when
d(s)-d(t)=1, a necessary condition, not a sufficient one.

Here is a separate constraint on degrees. For a self-converse class
the usual automorphism group acts freely on Hamilton orders: fixing
an ordered tuple fixes every vertex. Thus the number of tilings
presenting the class is H/|Aut|, odd by Redei. Converse-with-reversal
is an involution on these tilings; its number of fixed points is
odd. The fixed points are exactly its blue incidences, with a loop
counted twice. Therefore each connected component of the blue
graph on the self-converse classes has an even number of vertices.
This does not prove it is connected.

The same free action gives a useful arithmetic restriction not
involving shape: if a tournament has a prime count H=p>n then it
is rigid. Indeed |Aut| divides p, and p cannot divide n!. The
corresponding argument for H=7p>7n, with p>n, leaves group order
either 1 or 7. This is a direct example of why the same number can
carry different perspectives: c3 gives a pairing color, H gives
a group-order divisor, and neither encodes the other's full fibre.
