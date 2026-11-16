# All-order noncontainment from 27 to 63, and exact rewrites

An exact Collatz landing can exist with no induced extension between
the tournament counts. The old factored example was not the only
obstruction: **no count-63 tournament contains an induced count-27
tournament**, whatever the count-27 source.

## Reduction to two explicit targets

H is the number of directed Hamiltonian paths. It multiplies across
an ordinal sum of strong components. Since 7 and 21 are impossible,
a count-63 target has exactly one nonsingleton strong component and
that component has count 63. Busch's strong order bound
5^((n−1)/3)≤H bounds it by order 8 (63³<5⁸). An induced count-27
subtournament would still have count 27 after intersecting that
component, as outside vertices can only contribute factors of one.
The [complete order-8 archive](../actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md)
also captures smaller strong targets by neutral padding and
restriction. The two strong H63 classes both occur at order 8.
The independent subset path recurrence in check_27_63.py computes
all their induced counts. Neither contains 27. Thus the conclusion
is all-order, with the indicated finite atlas dependency.

## A sharp old-arc repair, on one continuing graph

Start from the previous H27 seed, outgoing masks
(62,124,120,112,32,64,17). Reverse old arc 3–6 and add one vertex
with outgoing mask 63 on the old labels. The new rows are
(62,124,120,48,32,64,153,63) and have H63. The one old reversal is
necessary over *all* choices of initial source by noncontainment,
and is sufficient by this witness. For this seven-vertex source one extra vertex is sufficient and necessary;
every H63 core itself has eight vertices. The [all-source cost
classification](CONTAINMENT.md) distinguishes outputs with and
without neutral padding.

The same physical graph can be continued, preserving all old
vertices and keeping every intermediate tournament strongly
connected. At each line, reverse the specified old pair and add
one vertex with the specified outgoing mask, using zero-based
unchanged labels for old vertices.

| New count | Old arc reversed | Added outgoing mask |
|--:|:--|--:|
| 63 | 3–6 | 63 |
| 189 | 0–2 | 24 |
| 567 | 0–6 | 288 |
| 1701 | 4–8 | 304 |
| 5103 | 3–6 | 16 |
| 15309 | 0–1 | 56 |
| 45927 | 0–1 | 8056 |
| 137781 | 3–6 | 2576 |

This is independent of the trivial ordinal-sum multiplication by
three. It gives connected changes through 137781,
but not a uniform rule picking future reversals.

For comparison, take the *factored* 27 source consisting of three
ordered C3 blocks. If no vertices may be added, the least edit cost
to reach 63 is exactly **five arc reversals**, checked by all subsets
of at most five of its 36 edges. With triangle orientation
0→1→2→0 in each block's local coordinates, reverse 0–1, 0–3, 2–5,
5–6 and 6–7. New outgoing rows
(496,509,473,465,480,396,32,320,64) give H63. edit.cpp reproduces
the complete fixed-order search, including 162 successful five-edge
sets. This cost is only for that factored graph and no new vertices.

## The 27→63 obstruction propagates up the ×3 tower

For every j≥0 take ANY count-27 tournament and adjoin j directed
triangles in ordinal sum. There is no induced enlargement to count
63·3^j. Indeed, any strong
target SCC that contains a factor 7 has count 7·3^a, a≥2 by the
two known exclusions; every other component has count a power of 3.
An entire old component lies in one target SCC. Every target SCC
whose count is 3^b can accommodate at most b old exponent units,
by induced monotonicity and the old factors 3 and 27. The old
exponent is j+3, so the 7-bearing SCC must house at least a+1 units.

If it has only old triangular components, it houses ≥a+1 triangles,
requiring ≥3(a+1) vertices. Busch excludes this for a≥2:
5^(a+2/3)>7·3^a at a=2 and the ratio grows by 5/3 per step.
If the original H27 is H9·3 and the target SCC contains the H9
component, it must also house ≥a−1 triangles. A strong H9 has order
at least 5, and the disparate old SCCs need an additional vertex
to join, so the new order is ≥5+3(a−1)+1=3a+3. Again Busch excludes
this, already at a=2. This also leaves no case for an original
source made entirely of triangles.

Finally, if the 7-bearing SCC contains a single old strong H27,
a=2 is the already-excluded H63 containment. For a≥3 it must
house ≥a−2 triangular blocks as well. The old vertices are
decomposable, so at least one new vertex is needed. If the old
core has order ≥7 the size is at least 7+3(a−2)+1=3a+2;
5^(a+1/3)>7·3^a for a≥3. For an order-6 old core, the size is at
least 3a+1; 5^a>7·3^a for a≥4. The only residual case is one
order-6 strong H27 and one triangle in an order-10 H189 target.
Test both old block orders, the two order-6 source classes
(from the same atlas), and each of 512 orientations of the
new vertex. None has H189. This exhausts the possible factorizations
of an original count 27 and proves the claimed all-order wall.

On the other hand, the explicit one-reversal, one-addition move for
the displayed strong H27 source works at EVERY j: preserve the old
triangular padding outside both sides, with ordered incidence. The
result is 27·3^j→63·3^j. The minimum possible number of old arc edits
for those particular infinitely many source graphs is exactly one.
The actual arcs have a fixed small description; mere count matching
and pure induction by induced inclusion fail.
