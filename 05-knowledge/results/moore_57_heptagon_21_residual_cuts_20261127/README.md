# A forced 7-to-21 completion and exact residual cuts in any degree-57 Moore graph

2026-11-27. This extends the [7 and 21 cycle comparison](../moore_fixed_order_spectra_20261124/THIRTEEN-ODD-ALPHABET.md)
and the [nested Moore shells](../nested_moore_shells_edge_grid_20261126/README.md).
Unlike the earlier [mixed tournament alignment obstruction](../moore_57_mixed_alignment_spectral_20261127/README.md),
the conditions here apply to **every** hypothetical degree-57 Moore graph.
They do not yield a contradiction. They give a forced graph-only search
seed and exclude separated or weakly connected residual completions.

## Every 7-cycle has two canonical 21-cycles

Write d for the Moore degree. Its graph has n=d²+1, no triangles or
squares, and precisely one common neighbor for each nonadjacent pair.
A 7-cycle has no chord (a chord would make a shorter cycle of length
at most 4). Write its vertices v_i with indices modulo 7.
The seven pairs at cyclic distance three give distinct outside common
neighbors w_i of v_i and v_{i+3}. An outside vertex cannot attach three
times: all pairwise cyclic distances would have to be at least three,
requiring circumference at least nine. All other outside vertices
either have one cycle contact (in group S_i at v_i) or none (in Y).
Every S_i has d-4 vertices.

The w_i are independent. For labels at offset ±3 they share a cycle
contact. At offsets ±1 or ±2 some pair of cycle contacts is adjacent.
A w-w edge would therefore create a triangle or square. A single
contact at v_j can neighbor w_i only if v_j is neither a contact of
w_i nor adjacent to any such contact. The only allowed column is
j=i+5. There is exactly one such single-contact neighbor z_i:
the pair (w_i,v_{i+5}) must have a common neighbor, and neither a
cycle vertex nor a double-contact vertex can serve. Moreover
z_i has one neighbor in each of S_{i+1},S_{i+2} and no other
neighbors in the single-contact groups. Its paths of at most two
edges to the other five cycle vertices are already supplied by
v_{i+5}, its two neighbors on the cycle, and w_i's contacts.
Thus z-z edges may occur only at offsets ±3.

The three disjoint classes (v_i),(w_i),(z_i) are 21 actual vertices.
The triples of edges

    v_i -- w_i -- z_i -- v_{i+5}
    v_{i+3} -- w_i -- z_i -- v_{i+5}

give two simple 21-cycles, since the induced steps on the first class
are +5 and +2 modulo 7. This pair is independent of how the
original 7-cycle is numbered or oriented: the base-to-double
incidences form a 14-cycle, whose two perfect matchings are exactly
the two choices. The original cycle edges are chords in each
21-cycle. The construction is injective on (base cycle, choice):
from a produced 21-cycle the three classes must be the three
phases modulo three. The double phase is independent. The z phase
cannot serve as a different base, because each site in either
other phase has just one neighbor in it, never two. This recovers
the original base phase and the matching used.

The actual 7-cycle count in a Moore graph is

    n d(d-1)^2 (d-2)(d-3) / 14.

For a trace proof, use nonbacktracking counts B_0=I, B_1=A,
B_2=A²-dI, B_j=A B_{j-1}-(d-1)B_{j-2} for j>=3.
The equations A²=(d-1)I-A+J and AJ=dJ give
tr B_7 = n d(d-1)^2 (d-2)^2 (multiply recursively in I,A,J).
A closed nonbacktracking walk of length seven is a cycle or a
pentagon with a single tail edge: after taking the shortest
repeated-vertex subcycle the girth and no immediate backstep
leave only these lengths. Each pentagon contributes 10(d-2)
such walks. Each edge lies on (d-1)^2 pentagons: pick the
other neighbor at each endpoint and join them by their unique
two-step route. Thus there are n d(d-1)^2/10 pentagons.
Subtracting the tails and dividing by 14 gives the formula.
[Recent primary graph paper, for the status of degree 57](https://www.mdpi.com/2227-7390/14/23/4318).

Petersen has zero 7-cycles. Hoffman–Singleton has 18,000
and hence 36,000 canonical 21-cycles. Every hypothetical
degree-57 Moore graph has 123,243,120,000 induced 7-cycles
and at least 246,486,240,000 ordinary 21-cycles.
This is a direct 7-to-21 relationship via Moore routes.
It does not identify any graph cycle census with the
[Chan–Pak Hamiltonian path spectrum](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf).

## Seven exact neighborhood packs

Write t_y for the number of w-neighbors of y in Y. Its contacts
with the base via its w-neighbors and its single-contact neighbors
partition the seven base vertices (unique routes).
Thus 2t_y + (# single-contact neighbors)=7, t_y<=3, and its
degree in Y is d-7+t_y. Counting the incidences with w_i and
pairs of such incidences gives

    |Y|=d²-7d+15,
    sum t_y=7(d-3),      sum choose(t_y,2)=14.

The latter equality holds because of the 21 unordered w-pairs,
seven already have a common cycle contact, and the remaining
14 have their unique common neighbor in Y: no single-contact
vertex touches two w's and the w's are independent.
If t_y=3, its three two-contact sets on the base must be a
matching of size three in the heptagram of distance-three
pairs. For each missed base vertex v_k, the unique such matching
consists of w_{k+1},w_{k+2},w_{k+3} (delete k from the heptagram,
then match the remaining six-site path). Two witnesses for the
same or cyclically adjacent k would share at least two w-neighbors,
contrary to girth. Thus b=#{y:t_y=3} is at most three.

For each i, put U_i=N(w_i)∩Y and L_i=N(z_i)∩Y. The closed
Y-neighborhoods of the d-3 members of U_i are disjoint, since
they have their one common neighbor w_i outside Y. Their
complement in Y is exactly the set L_i of size d-4. Indeed
the only possibilities for a short route from y to w_i are the
direct one, a two-step route via a member of U_i, or a two-step
route via z_i. The size of L_i follows because z_i has one
base, one double, and two single-contact neighbors. For matrices M,L of these seven indicators
and residual adjacency B, necessary exact tests are

    (I+B) M = J-L,
    M^T M = (d-4)I + J - A(step 3 on Z/7),
    (M^T L)_{ij} = 0 if i=j, j-i=±2, or z_i~z_j; 1 otherwise.

The cross term follows by asking whether the pair w_i,z_j
has a route in zero-contact vertices: diagonal pairs are already
adjacent; at offsets ±2 they have a common base contact; the
only other possible outside-Y intermediate is z_i. The full
boundary, after all single-contact vertices are placed, also
satisfies the usual Moore block equations.

For d=57 the residue Y has 2,865 vertices and degrees 50–53,
with census (2501-b, 350+3b, 14-3b, b) in that order. Each U_i
has 54 vertices and its closed neighborhoods cover 2,812 sites,
leaving exactly the 53 in L_i. For comparison, across all
18,000 Hoffman–Singleton base cycles the replay finds b=0,
no z-z edges, and Y consists of two 7-cycles and an isolated
vertex. The same extension arithmetic at degree 57 forbids
any such separation.

## No weakly connected residual completion: five exact cut values

Let kappa_r be the minimum size of a vertex deletion K making
Y disconnected with every remaining component having at least
r vertices. For every degree-57 candidate,

    (kappa_1, kappa_2, kappa_3, kappa_4, kappa_5)
      = (50, 98, 146, 193, 239).

**Lower bounds.** The Moore adjacency has nonprincipal maximum
eigenvalue 7. For any deletion disconnecting Y, at least one
component C has average internal degree at most 7: otherwise
a nonzero linear combination of the indicators of two
components, centered in the full graph, has quotient greater
than 7. Write c=|C|. The boundary incidence a_z at z in K
satisfies E=sum a_z>=43c and sum a_z(a_z-1)<=c(c-1) by
the no-square rule. Cauchy and monotonicity in E give

    |K| >= E²/(E+c(c-1)) >= 1849c/(c+42).

For r=1,2,3,4 this already exceeds the claimed cut value
when c>=r+1. At c=r use degree 50 and no shared external
neighbors for pairs adjacent or at distance two inside C.
A singleton costs 50; an edge costs 98; a three-site path
costs 146. A four-site tree costs at least 193: its raw
external incidence is at least 194 and at most one pair is
at internal distance three, allowing one duplicate. At r=5
the general bound exceeds 239 for c>=7; for c=6 the raw
incidences less possible duplicated pairs already give
50c-2e-choose(c,2)>=255, since e<=15.
A five-cycle costs at least 240 with no duplicated pairs.
A five-site tree costs at least 242 minus the number of
pairs at internal distance at least three, at most three
(attained on the path). It therefore costs at least 239.

**Constructing all cuts.** Let H be the graph induced by the
vertices of Y of degree 50 in Y. Its order h is at least 2498,
there are at most 367 other vertices of degree at most 53,
and its average degree a is hence at least 50-367*53/2498
>42.2. The [irregular Moore bound of Alon, Hoory and Linial](https://web.math.princeton.edu/~nalon/PDFS/ahl1.pdf)
forces a pentagon in H: girth at least six at average degree
greater than 42 would need more than 2(1+41+41²)=3446
vertices. If needed, prune vertices of degree at most one
before applying the bound; average degree does not decrease.
A consecutive path of each size 1 through 4 on this pentagon
has precisely the required external boundary size in Y.
The only duplicated external neighbor for the four-site
path is the fifth site of the pentagon, shared by its endpoints.

For r=5, H has an induced five-site path with all three long
pairs (at path distances 3, 3 and 4) having their common
neighbors in Y. Here is a counting proof. Its maximum degree
is at most 50. There are at least h a^4 oriented four-edge
walks by two applications of Cauchy to the adjacency matrix:
sum of the four-edge walk counts is ||A_H²1||² and
1^T A_H² 1 is the sum of squared degrees. At most
4h*50²*a of these walks either have a repeated site or
close an endpoint edge: for repetition use the three
possible adjacent backtrack positions (any other repeat
makes a triangle or square). Each can be counted by a
three-edge walk or three-edge star, with at most 50²
choices after the first directed edge. For endpoint closure,
choose the first three edges and then at most one common
neighbor of the start and the fourth vertex. There are no
other possible chords. So at least h(a^4-10000a) walks
are induced paths.

For a degree-50 site, its seven outside-Y neighbors are all
single-contact sites. Each has at most 53 neighbors in Y.
For special sites z_i this was shown above. An ordinary
single-contact site has only its one base contact, no
double-contact neighbors, and four neighbors in the
single-contact groups (one for each required route to the
four base sites not already reached); it has 52 neighbors
in Y. Therefore any start of a long pair has at most
7*52=364 bad endpoints, whose common neighbor is outside Y.
For fixed endpoints there are at most 50 paths of length
three (choose the first neighbor, then the middle site is
a unique common neighbor) and at most 50² paths of length
four (choose the first and last neighbors). For each of the
three long-pair positions there are at most h*364*50²
bad paths, extending the length-three path by one step
when needed. But already at a=42.2,

    a^4 - 10000a > 3*364*50² = 2730000.

This proves the desired path exists. The three duplicated
external neighbors are distinct, since one serving two long
pairs would also serve a short pair. They give boundary size
5*50-8-3=239.

Finally, deleting the external boundary of any of the five
chosen paths leaves the path as one component. Any site
outside the path and its boundary loses at most r neighbors
to the deletion: two deleted neighbors that both touch
the same path site would give two common neighbors. Those
other sites still have degree at least 50-r>=r, so none of
their components is too small. There are ample remaining
sites. This proves the five exact cut values.

The contrast with Hoffman–Singleton is a graph-level
obstruction to copying its disconnected residue, not a
nonexistence proof. A next test is to combine the seven
overlapping near-perfect neighborhood packs with the required
tight cut around a degree-50 path; both constraints now apply
to every candidate without imposing any symmetry or assuming
a smaller Moore subgraph. Run node replay.js for the
finite-graph census and the exact inequalities.
