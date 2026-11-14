# Hamilton-cycle necklaces and the open-end fans of the merged metagraph

[THM-4729](../../../01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md)
counts closed presentations. The following refines those counts by exact
rotational stabilizer and describes the open ends attached to each unrooted
cycle. Directed cycles here are identified under rotation, not reversal.

## Rooted versus unrooted, including self-converse pairs

On vertices 0,...,n−1, n≥3, prescribe the ring i→i+1 mod n.
A closed tiling is exactly a tournament extending this ring. Quotienting
by rotations gives isomorphism classes of pairs (tournament, distinguished
directed Hamilton cycle). In contrast, holding the root fixed recovers the
specific closed tiling. Let d divide n, m=n/d, b=binom(n,2).

For d odd, the number F_d of closed tilings fixed by the subgroup of
rotations of order d is
$$
 F_d=2^{(b-n)/d}.
$$
For d even it is zero: the half-turn in the subgroup swaps the ends of
some arc. For odd d, no nontrivial rotation in the subgroup can fix an
unordered vertex pair, so there are b/d free-or-fixed orientation
coordinates after imposing the subgroup action. The prescribed n edges
use exactly n/d of these coordinates. Set
$$
 E_d=\sum_{d\mid e\mid n}\mu(e/d) F_e .
$$
Then E_d counts the rooted tilings of exact rotational stabilizer d,
and E_d/m counts unrooted directed-cycle classes of that type. This
global formula refines the closure identity n·HC(T)/|Aut(T)| of THM-4729:
within a fixed tournament class, each cycle orbit of stabilizer d
contributes m rooted tilings, not n. Conversely, summing rooted counts
recovers the original closure identity.

To include the converse merge, use the fixed path mirror
s(i)=n−1−i, applied as a relabeling AND reversal of every arc.
It preserves the prescribed directed ring. Call a tiling blue if it
is fixed by this operation. For d odd write
$$
 R_d=2^{(b/d+\lfloor n/2\rfloor-m-u)/2},\quad
 u=\begin{cases}1&n\text{ odd},\\2&n\text{ even}.\end{cases}
$$
For d even put R_d=0. This counts the blue tilings fixed by the same
rotation subgroup. To see this, the mirror acts on the b/d unordered-edge
orbits. An invariant edge orbit is a torsor of odd size d where the
reflection acts by j→c−j, so each such orbit contains exactly one
mirror-fixed edge. Conversely every mirror-fixed edge is in such an
orbit. The fixed edges are the floor(n/2) transposed pairs; a mirror
fixes at most one vertex, hence none of these constraints forces an
arc to be its own reversal. Nonfixed edge orbits are paired. For the
prescribed ring edges their d-rotation quotient has m edges indexed
by the start of the oriented edge, and the mirror acts by i→−2−i
mod m. This has u fixed orbits. Subtract their (m+u)/2 coordinates
from the total to obtain the exponent above.

Let A_d=sum_{d|e|n} μ(e/d)R_e be the number of blue rooted tilings
with exact period m. In a rotation orbit fixed by reversal-mirror,
the latter acts on the m roots by i→c−i, giving one fixed root when
m is odd. For m even, c is even, giving exactly two fixed roots:
otherwise a reflected representative with odd shift of the chosen
mirror would fix two opposite vertices of the original n-gon,
which is impossible for an anti-automorphism of a tournament.
Thus division of A_d by 1 (odd n) or 2 (even n) counts self-converse
cycle classes with this stabilizer; call this quotient B_d.
The number of classes of pairs modulo isomorphism *and* converse is
$$
 \frac{E_d/m+B_d}{2}.
$$
These self-converse marked cycles are not the same as self-converse
unmarked tournaments: the anti-automorphism must preserve the
particular cycle up to reversing its order. In particular, at a
self-converse metagraph node, only some cycle necklaces may contain
blue presentations. A genuinely chiral cycle necklace and its
converse have no blue roots.

## The open fan of one directed Hamilton cycle

Let Z be a directed Hamilton cycle of T, and T^Z be T with all edges
of Z reversed. This keeps every outdegree s_v unchanged. Rooting Z is
equivalent to deleting its closing edge e from the cyclic word to read
a Hamilton path. The merged-metagraph mate is isomorphic, up to
converse, to T^Z with just e reversed back: reversing the path arcs
instead of every nonpath arc differs only by taking the converse.
This follows the gauge of Proposition K in THM-4729.

The closed tilings associated to one automorphism orbit of the marked
cycle therefore form a fan of m=n/d ordinary (pre-converse-merge) line
incidences, one for each cut up to the cycle's stabilizer of order d.
The m open endpoints need not be m different unmarked tournament classes.
In the converse-merged graph a chiral necklace and its opposite give
one fan of m lines. A self-converse necklace gives (m+u)/2 lines: u blue
cuts are fixed, while the others pair. The sums below use the ordinary
m incidences (equivalently weighting a nonblue line in a self-converse
fan twice and a blue line once). Every edge occurs in an
orbit of length d on the ring. For any function phi on integer scores,
put F(T)=sum_v phi(s_v) and use T^Z with one ring edge flipped back
as the representative on the open side. Then
$$
 \sum_{\mathrm{fan}} (F(\mathrm{open})-F(T))
 =\frac1d \sum_v (\phi(s_v-1)+\phi(s_v+1)-2\phi(s_v)).
$$
Indeed before the one-edge flip, T^Z and T have the same scores.
As that edge traverses the reversed Z every vertex appears once as
its tail and once as its head. Each flip lowers the tail's score by
one and raises the head's score by one. Dividing by d takes one
representative per edge orbit. Use converse-symmetric phi to obtain
an invariant of merged nodes without this representative convention.

For the converse-invariant triangle count
c3(T)=binom(n,3)−sum_v binom(s_v,2), the discrete second difference
is −1 at every vertex. Thus the *exact sum over each closed
Hamilton-cycle fan*, not merely over the whole metagraph, is
$$
 \sum_{\mathrm{fan}} (c3(\mathrm{open})-c3(T))=-m.
$$
This implies every fan has an open mate with strictly fewer cyclic
triangles. It is a descent in cyclic triangles, not in path count,
so it supplies structural information for the enumeration problem
without a path-spectrum completeness proof.

For a nonregular tournament, along any directed cycle visiting all
vertices, the differences of consecutive scores include a positive
and a negative integer. So a fan from a nonregular class has both
a descent of at least two and a mate that does not descend. For a
regular tournament all its fan lines descend by exactly one. Iterating
the certified descent from any strong node must eventually reach a
nonstrong node because cyclic triangle count is a nonnegative
integer. This is not an assertion of connectivity between different
nonstrong nodes.

## Strong bridges count the nonstrong side from the cycle side

Let b(U) be the number of strong bridges of a strong tournament U:
an arc e is a strong bridge if U−e is not strongly connected.
In the fan of (T,Z), write U=T^Z. The open mate at cut e is
nonstrong exactly when e is a strong bridge of U. In fact if U−e
is strong then adding the reverse arc retains strong connectivity;
otherwise there is a directed cut for which e was the sole arc in
one direction, so adding its reverse cannot restore it. Every
strong bridge of U lies on every directed Hamilton cycle of U,
because such a cycle would by itself ensure connectivity after
deleting the arc if it avoided the arc. Since the reversed ring
is a Hamilton cycle of U, b(U)/d is exactly the number of
nonstrong partners of a single necklace, as line incidences rather
than distinct unmarked neighbors.

Cycle reversal (T,Z) ↦ (T^Z,Z reversed) is an equivariant
involution of marked-cycle classes that retains stabilizer order d.
Summing all incidences over ordinary (not converse-merged) closed
tilings, THM-4729 gives the exact identity
$$
 \mathrm{NS}(n)=\sum_{[U]\text{ strong}}
     b(U)\,\frac{\mathrm{HC}(U)}{|\mathrm{Aut}(U)|}.
$$
Here NS(n) is the number of nonstrong tilings, not the number of
nonstrong tournament classes. To check the coefficient, an orbit of
cycles with stabilizer d contains |Aut(U)|/d cycles, and contributes
b(U)/d incidences. The contribution at each strong class is an integer.

If U is regular of order n≥5, then every cut A,U−A has equal
traffic in both directions, because the net degree imbalance summed
over A is zero. Each way contains |A|(n−|A|)/2 ≥2 arcs, so U has no
strong bridges. Its entire fan therefore stays inside the strong
region while decreasing cyclic triangle count by one. The test is
about connectivity of the complement line's mate, not about the
number of Hamilton paths being monotone.
