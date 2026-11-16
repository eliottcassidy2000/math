# The factor-one reverse bridge has only two right-hand profiles

In the [rooted reverse bridge](BRIDGE.md), let E_V and O_V
be the sums over even- and odd-order paths ending at the
right root, weighted by the Hamilton count of the complement.
The low left O-coordinate has coefficient 2E_V. If this is
to be the primitive +2 dial with a unit O edit, the right
profile must have E_V=1. In that case there are only two
possibilities for the triple (H,E_V,O_V):

    (1,1,1) or (3,1,2).

Proof. Each in-neighbor of the right root contributes a
positive term to E_V from a path of order 2. Thus it has
exactly one in-neighbor w, and removing the root and w leaves
a transitive spine, since the one remaining complement count
must equal 1. A two-edge path in the spine-plus-w ending at
w would extend to a second even-order root-ending path, so
no such path exists. Every spine vertex beating w therefore
beats all other spine vertices; at most one can do so.
If none does, w and the root are the first and second
vertices of a transitive tournament. Otherwise that first
spine vertex, w, and the root form a directed triangle
preceding the rest of the spine. The weighted odd root-ending
paths are the singleton root in the first case, and the
singleton plus the order-three path in the second. That
gives the stated profiles. Their converse describes the
analogous factor-one root-start profiles.

So this bridge genuinely has a controllable factor three
(from E_V=3), but fixing E_V=1 does not leave a right-hand
four-state dial. It is not a failure of the observed
two-dial construction: that construction uses a different
one-head interface.

## What the retained-root read actually remembers

Keep the original left root x after the bridge to a right
root y. Put g=H(U−x) and K=P_V(2)/2. Write E_s,O_s for right
**start** counts (not its end counts E_V,O_V used in H).
Then the exact retained-root update is

    H' = h_V H + 2(E O_V+O E_V),
    E' = K(E+O) − g O_s,
    O' = K(E+O) − g E_s,
    g' = h_V g.

Any root-start path either stays in U or consists of a left
root-start part and one right part; it cannot use y→x, as x
is already visited. Its complementary path count factors
across the forward join after removing x. For a full
unreversed forward join, the two parities of the optional
right run each have total weighted count K. Indeed the
alternating sum over the two ordered slots vanishes by
assigning each odd-cycle/singleton component to one slot
with opposite signs. The exceptional forbidden junction
x→y eliminates precisely a singleton left part followed by
right paths starting at y, weighted by g. That swaps the
parities in the subtracted terms. Thus one must retain g:
even E+O equality does not alone synchronize the new read.

If h_V=3, its unique cyclic SCC is a triangle, with s
singleton components. P_V(2)=12·2^s, and its total weighted
root-end count E_V+O_V is 3·2^j, as follows by reading
preceding SCCs through their two-slot counts and subsequent
SCCs through H. The root SCC itself contributes 3 for a
triangle and 1 for a singleton. For s≥1 the possible j are
exactly 0,…,s+1, by placing the root and triangle; if s=0
only 0 occurs. This explains why that right factor has
limited prepared readouts rather than a new free radix digit.
