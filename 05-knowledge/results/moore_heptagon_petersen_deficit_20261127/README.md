# The heptagon witnesses count missing Petersens

2026-11-27. Structural companion to the [21-to-63 detours](../moore_57_heptagon_63_detours_20261127/README.md).
These are necessary properties of a degree-57 Moore graph, if it exists.
The Petersen upper bound itself is known; the result here is its
exact heptagon-witness and rooted-transport interpretation.

Let P be the number of induced Petersen subgraphs and P_y the
number containing y. Around a heptagon C use the [existing
notation](../moore_57_heptagon_21_residual_cuts_20261127/README.md):
w_i is the common neighbor of v_i,v_{i+3}, z_i is the sole
single-contact neighbor of w_i, and b(C) counts sites y outside
the cycle with no cycle contact and exactly three w-neighbors.
Call (C,y) a witness. Then

    number of witnesses at y = 6(819280 - P_y),
    sum_C b(C) = 60(266266000 - P).

## Root transports measure the local deficit

At y, fix an unordered triple of its second-shell branches.
An ordering A,B,C induces the permutation L on the starting
branch by the three matching steps A->B->C->A. It has no fixed
slots by girth. Each 2-cycle of L traces a hexagon, with two slots
in each branch. The three branch parents each join an antipodal
pair of the hexagon and all meet y. This is the ten-vertex
Petersen pattern. It has no extra edges, by its internal
diameter two and girth at least five. Conversely a Petersen
containing y determines the branch triple and a transposition
of its transport. Reversal or cyclic change of branch order
does not affect the number of transpositions.

From a slot not in a transposition, make the five-edge walk
with labels A,B,C,A,B,C. Its six sites are distinct. The ends
are not adjacent (that would give a fixed point of L²) and
have different neighbors in branch B (otherwise L fixes
the first slot). Their unique common neighbor is in a fourth
branch. It cannot be y, a branch parent, or a slot in A, B
or C. The new site closes a heptagon, where the three parent
vertices of A,B,C are precisely the w-contacts of y. All
seven sites are second-shell slots, so y has no cycle contact.
Conversely the three w-contacts of a witness pair six sites
of its heptagon, covering all but one base site. Their labels
are A,B,C,A,B,C along the remaining path (the three w labels
are consecutive). Thus the six branch orders count each
witness twice, once per path direction.

If the triple transport has t transpositions, it contributes
3(56-2t) witnesses. Sum over choose(57,3) triples, and use
28 choose(57,3)=819280. Then sum over 3250 roots and use
sum_y P_y=10P. The identity does not show the deficit is
positive. The degree-57 case in [Junker's involutional model](https://home.mathematik.uni-freiburg.de/junker/preprints/mj-hs.pdf)
assumes extra conjugation relations to get his exclusion.
The published global count and percolation dichotomy are
in [Ibrahim–LaFayette–McCall, Theorem 4.18](https://ajc.maths.uq.edu.au/pdf/93/ajc_v93_p060.pdf).

## Visible reasons for a triple to percolate

In two-neighbor closure, add any vertex with two already present
neighbors, repeatedly. An independent triple without a common
neighbor determines a unique hexagon: take the common neighbor
of each pair. Each hexagon gives exactly two such triples.
By [Ibrahim–LaFayette–McCall, Proposition 4.16 and Theorem 4.18](https://ajc.maths.uq.edu.au/pdf/93/ajc_v93_p060.pdf),
each triple either closes to a Petersen or percolates to all
3250. A hexagon can be in at most one Petersen, since both
its three antipodal intermediates and their common center
are determined. Each Petersen has ten hexagons.

A witness with neighbors w_0,w_1,w_2 designates two independent
triples {v_0,v_2,v_4} and {v_1,v_3,v_5}. The first triple
builds the hexagon v_0,v_1,v_2,v_3,v_4,w_4. Its antipodal
intermediates are w_0,w_1,z_4. The only possible Petersen
center is y, already the unique common neighbor of w_0,w_1.
But y cannot neighbor z_4: then y and v_2 are both common
neighbors of z_4,w_2. For the second triple the antipodal
intermediates are w_1,w_2,z_5, and the forbidden duplicate
uses v_3 and w_0. So both triples percolate.

Each percolating independent triple of this type is designated
by exactly six witnesses. To see the upper bound, there are
three choices of which two-edge segment of its hexagon to
replace by a three-edge segment. Give the surviving four-edge
path one fixed orientation v_0,...,v_4. Now y is the forced
common neighbor of w_0,w_1, the antipodal intermediates for
v_0,v_3 and v_1,v_4. Its third w-contact is the unique common
neighbor of y and v_2. It can play the role w_2 or w_6:
in the first case v_5 is the common neighbor of w_2,v_4 and
v_6 that of v_5,v_0; in the second reverse the reconstruction.
The upper bound of six is attained by the incidence count:
there are 2,662,660,000 hexagons (the nonbacktracking length-six
trace divided by twelve, equal to n*d*(d-1)^2*(d-2)/12).
The number of percolating triples of this type is therefore
2(2662660000-10P). Six times this number equals twice
sum_C b(C), with each witness contributing exactly two.

## Closure distinguishes 21 from 63

From the heptagon, the first simultaneous round adds its
seven w sites. Round two adds the seven z sites and 14-2b(C)
distinct zero-contact sites. The fourteen pairs of w's with
no common base contact each give a site; a triple witness
accounts for three pairs at one site, subtracting two.
Thus the first three sizes are 7,14,35-2b(C), and each
canonical 21-cycle has the same eventual closure as its base.
By [the closure classification, Lemma 4.9 and Corollary 4.10](https://ajc.maths.uq.edu.au/pdf/93/ajc_v93_p060.pdf),
a proper closure containing the base must be Hoffman–Singleton:
the only other proper closures are stars, smaller Moore graphs,
and the trivial sets. The Hoffman–Singleton replay gives
7,14,35,50 on all 18,000 heptagons. In particular b(C)>0
is a sufficient certificate that the heptagon and its
21-cycles percolate.

Every 63-cycle percolates, including the marked ones in the
detour construction: even a full star has just 58 vertices,
and a smaller Moore graph has at most 50. If the 21-cycle
already percolates nothing further is needed. Otherwise,
its closure is one Hoffman–Singleton copy and adding any
vertex from the 63 extension outside that copy already
forces percolation (no proper closure can contain that copy
and another vertex). This distinction is about actual graph
closure and does not explain tournament-spectrum gaps.

Run node replay.js for the controls: counting all root triples,
hexagons, heptagon witnesses and closure rounds in the two
smaller Moore graphs. The degree-57 assertions are arithmetic.
