# Three copies without confusing the cyclic lift and a product

For any tournament T and any chosen vertex v, replace v by a
directed triangle whose vertices have the same exterior arcs as v.
Distinguish skeletons of a Hamiltonian path after replacing each
maximal internal run by the same symbol v. Let A_j be the number
of valid skeletons with j occurrences of v and every exterior
vertex occurring exactly once; occurrences of v are never
adjacent, and ordinary exterior edges must be respected.
The j runs are ordered by their position in the path.

The directed triangle has 3 full Hamilton paths, 6 ordered covers
by two nonempty paths, and 6 ordered covers by three. Therefore
the exact substitution formula for any exterior orientation is
```
H(T[v←C3]) = 3 H(T) + 6 (A_2+A_3).
```
This follows by independently choosing the ordered path cover
inside the triangle after the skeleton is prescribed. It is the
three-vertex case of slot substitution, now with an exact
criterion for when the offset vanishes.

**The offset vanishes iff v is a singleton SCC of T.**
Two appearances of v in a skeleton require a directed cycle
through v by reading between consecutive appearances. Conversely,
if v belongs to a nonsingleton SCC, a directed Hamiltonian cycle
in that SCC can be opened at v and read with v at both ends.
Concatenate Hamilton paths of SCCs before and after in their
unique order, producing A_2>=1. For the elementary longest-cycle argument, see the
[order-seven cutoff](SEVEN_CUTOFF.md#move-the-cutoff-to-every-higher-order).
Hence in a strong source on more than one vertex the multiplication
has a strict positive excess, at least 6.

For instance starting at the directed triangle (H=3) and
inflating one vertex by a directed triangle gives a strong
H=15 (A_2=1, A_3=0), not H=9. Inflating any vertex of a
strong H=15 graph in this way cannot give H=45, since it
gives at least 51. Yet there are several strong H=45
witnesses: the four-run apex word, and the cyclic bouquet
with each of three slots a transitive pair. In that bouquet
each pair has ordered nonempty cover numbers L_1=1,L_2=2.
The cyclic run formula reads
3 L_1³ + 3 L_1² L_2 + 3 L_1 L_2² + 3 L_2³ =45.
This geometry is different from exact multiplication by
a separate triangle SCC.

In particular a scalar descent rule of dividing any multiple
of 3 by 3 cannot be justified by imagining a cyclic vertex
module, even when such a module exists. The exact input to
the inverse is H+2(A_2+A_3), not H. In the known order-eight
63 representatives the automorphism group is trivial and
there is no cyclic-triangle module: such a module would
supply a cycle of three automorphisms. The
[THM-344 class census](../../../01-canon/theorems/THM-344-n8-h63-complete-omega.md)
documents the two representatives; the companion replay
checks the displayed one directly. The new four-run proof
constructs 63 independently of the census.
