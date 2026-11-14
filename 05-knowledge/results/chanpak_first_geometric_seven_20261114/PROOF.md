# The first 7p count that admits a geometric rotation of order seven

This follows the [symmetry budget](
../chanpak_least_missing_jacobian_20261114/BLACK_BOX_AND_SYMMETRY.md).
All claims concern the order-seven symmetric slice; they do not
classify the rigid witnesses at the same count.

**Theorem.** The smallest Hamiltonian path count of the form 7p,
p prime, that admits an order-seven automorphism is 4151=7·593.
At equality the strong witnesses form one converse pair, at
order ten. At every order n≥10 there are exactly 2(n−9)
ordinary isomorphism classes and n−9 converse-merged classes
in this symmetric slice. Each has automorphism group exactly C7.

## Reduction to one homogeneous orbit

A permutation of order seven moves at least one block of seven
vertices. On every such orbit the induced orientation is
circulant with one forward step from each opposite pair.
The eight choices separate into two multiplier classes with
step sets {1,2,3} and {1,2,4}. The literal path DP gives
(h,f)=(175,1988) and (189,2016), where f is the ordered
two-path-cover count allowing empty slots. Independent odd
cycle packs on two disjoint seven-orbits give H≥175².
This is already greater than 4151. Thus for a candidate
at or below 4151 there is exactly one moving orbit R.
Every vertex outside it is fixed, so it sees all of R
uniformly. Collapse R to a quotient vertex v.

If v is a singleton SCC of the quotient, the inflated
tournament has H divisible by 175 or 189. Dividing by 7
leaves either a factor 25 or a factor 27, not a prime.
If its quotient SCC has exactly three vertices, it is
a directed triangle; the inflated SCC has H=h+f. This
is 2163 or 2205 and again has composite quotient by 7.
Its ordinal context only adds a multiplicative factor.

If the quotient SCC has size at least four, it has an induced
strong four-vertex subtournament containing v. For clarity,
here is the short vertex-pancyclicity argument for this step.
In a strong tournament a vertex v belongs to some triangle:
if all its in-neighbors beat all its out-neighbors, the two
sides cannot participate together in its SCC. Start with a
triangle through v. An outside vertex with mixed incidence
to that triangle can be inserted at a transition in the
directed circle. If no such vertex exists, the outside
vertices split into those beating the whole triangle and
those losing to it. Both sides are nonempty by strong
connectivity, and there must be an arc from the losing side
to the beating side. Replace any triangle vertex other than
v by those two vertices, in that order, in the old circle.
This gives the required four-cycle.

## An exact four-vertex quotient atlas, valid for any replacement module

Label the triangle v→b→a→v and append x. Its outmask r
is read over (v,a,b), with v the least significant bit.
If a marked strong order-four quotient is given, one may
choose a triangle through v, so this table is exhaustive
for the rooted quotient classes (with possible duplicates).

Every path in an inflated chart either traverses R once or
traverses it twice. A return requires at least two external
vertices between visits, because a single external vertex
cannot both beat and lose to a homogeneous module. There
are only three external vertices. A once-through word is
an ordinary Hamilton path of the quotient and contributes
h; a twice-through word contributes f−2h (both runs
nonempty). Thus one can count the quotient run words before
choosing any R. The table is

| Masks r | Once-through words | Twice-through words | H |
|---|---:|---:|---|
| 0,7 (not strong) | 3 | 1 | h+f |
| 4,5 | 5 | 2 | h+2f |
| 1,2,3,6 | 5 | 3 | 3f−h |

The replay enumerates those short words independently and
also checks every inflated tournament. For the two rotors
the strong cofactors obtained are 593,827 and 603,837
respectively. The minimum 593 is prime and attained
by masks 4,5 of the first rotor.

Induced containment makes the path count nondecreasing
by the positive odd-cycle packing formula. Thus the
four-vertex lower bound is a lower bound for the original
candidate. At equality its additional vertices must be
in singleton SCCs: any added vertex on a triangle would
supply a new positive packing, and a vertex on no triangle
belongs to no nontrivial SCC. To see the latter, its
in-neighbors then all beat its out-neighbors, with no return
from the out-side to the vertex. Therefore the strong
part is precisely one of the two order-ten cores at equality.

The two cores are converse: taking converse and exchanging
a,b changes mask 4 to 5, and the {1,2,3} rotor is isomorphic
to its converse by index negation. The cores are not
isomorphic: the uniquely moving orbit is invariant under
isomorphism, and its quotient vertex v has different
outdegrees in masks 4 and 5. Uniqueness of the moving orbit
as an isomorphism invariant follows because the whole group
is C7: the earlier symmetry budget for a 7p count excludes
any factor p=593 in the automorphism group. Remaining
singleton SCCs are distributed before and after the sole
strong component, so there are n−9 placements per core.
No class is self-converse, since the cores are chiral.

## Metagraph regression at the unique strong merged node

The mask-4 orientation is a cyclic wrap with two portals
b→a, then the ordinal source R→x, then back to b.
A directed Hamilton cycle is completely specified by
one full path of R between the two portals and x. Its
unrooted ordinary count is therefore 175. Every
automorphism fixes the portals, so its action on directed
cycles is free; there are 25 unrooted cycle-marked classes.
On each ordinary strong core, the 593 path-marked classes
split into 250 closed and **343 open**. In formulas,
175·10/7=250 and (4151−175·10)/7=343=7³.
The same counts apply to its unique converse-merged node.

In fact each of its 25 cycle fans has exactly three closed
cuts leading to nonstrong open mates. To see it, reverse
the full cycle b→a→(path on R)→x→b and let u and w be
the first and last vertices of its R path. In the reversed
tournament, b has only the outgoing edge b→x, x has
only the outgoing edge x→w, and a has only the incoming
edge u→a. Those are strong bridges. They are the only ones.
Within R a cut originally had |S|(7−|S|)/2≥3 arcs in
each direction, by regularity. Reversing a path changes
this count by at most one (net flow through a path cut
depends only on its ends). Thus R remains strong after
deleting any one internal edge. All its other crossings
with a,b,x can be bypassed: every rotor vertex beats b,
a reaches all of R except u, R is internally strong, and
the three listed edges still provide the portal closure.
More explicitly a→b or a→x can instead go through R;
an edge from R to x can instead go through b; losing one
of the parallel contacts from a to R or R to b leaves
others. By the [cycle-fan bridge criterion](
../chanpak_h9_cycle_metagraph_20261114/CYCLE_NECKLACES.md),
the nonstrong boundary therefore has 25·3=75 line
incidences. The other 175 closed lines lead to strong
partners. The replay checks this entire fan neighborhood.

At n>10 every class in this symmetric slice has singleton
SCC padding and has no Hamilton cycle. Its 593
Hamilton-path incidences per node are all open. The
generating functions for the ordinary and converse-merged
class counts across order are 2z^10/(1−z)^2 and
z^10/(1−z)^2. The labelled count at order n is
2(n−9)n!/7. These formulas do not count any rigid classes
with H=4151 and do not imply that earlier 7p counts are
unattainable by rigid tournaments.

Run `python replay.py`.
