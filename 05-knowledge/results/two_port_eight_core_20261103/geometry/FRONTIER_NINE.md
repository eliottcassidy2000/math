# The first Moser-containing nine-point core has an exact two-hub relation

**Status: PROVED reduction; the resulting obstruction problem is OPEN.**
This does not assert a nine-private-point blocker exists or is impossible.
It applies only when the private nine-vertex minimum-degree-three graph
contains a seven-vertex Moser subgraph. Other nine-vertex core types are
outside its scope.

Keep the old Moser boundary and five named colours as in
[PROOF.md](PROOF.md), and retain all unit contacts. Write the private
vertices as a Moser core K together with x,y. By E3's rigid-spindle contact
bound, each of x,y has at most two unit neighbours in K. Since each has
private degree at least three, **xy is an edge and each of x,y has exactly
two neighbours in K**. Thus the first extra private layer is necessarily
the interacting two-hub configuration of H1 in
[two_port_beyond_nine_20261103/geometry/PROOF.md](../../two_port_beyond_nine_20261103/geometry/PROOF.md).

Fix any proper old boundary colouring b. Let L_v be its available private
colour lists, each of size at least three. Let C_b be the full set of
proper list colourings of K, with all contacts from K to the old boundary
included. The inherited seven-private-point extension theorem guarantees
C_b is nonempty for every such b. For f in C_b put

    A_x(f) = L_x minus { f(u) : u in N(x) intersect K },
    A_y(f) = L_y minus { f(u) : u in N(y) intersect K }.

Both sets are nonempty, because only two K-neighbours can remove colours.
The exact number of extensions over all nine private vertices is

    Z(b) = sum_(f in C_b) ( |A_x(f)| |A_y(f)| - |A_x(f) intersect A_y(f)| ).

The subtracted term removes precisely the assignments violating xy.
All terms are nonnegative integers, and a term vanishes precisely when
A_x(f)=A_y(f) is a singleton. Therefore

    Z(b)=0  iff  for every f in C_b,
                  A_x(f)=A_y(f) is a singleton.

In particular, any such blocking assignment must have |L_x|=|L_y|=3.
For every admissible f, the two K-neighbours of each hub must receive two
distinct colours inside that hub's list, and the two leftover colours
must agree. If either original hub list has size at least four, every
core colouring has a positive pair extension and blocking is impossible.

This identifies the missing relation rather than replacing it by pair
marginals: a candidate must force **all** admissible core colourings into
the identical-singleton condition, while realizing both pairs of contacts
and the edge xy by actual distinct planar unit-distance points. A zero for
one b would still be only a conditional boundary obstruction. Six-chromatic
forcing requires eliminating every proper old five-colour boundary state.
