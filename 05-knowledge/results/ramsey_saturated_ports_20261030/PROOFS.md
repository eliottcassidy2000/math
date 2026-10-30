# Saturated Ramsey neighborhoods and faithful symmetry restriction

**Status: PROVED implications, with the stated CITED Ramsey inputs.** This
packet neither constructs nor excludes a `(3,3,5;56)` coloring. Its exact
Ramsey value remains OPEN. No literature priority claim is made for the
derived implications. Local labels below are scoped to this file.

An `(s,t,u;n)` coloring is a red/green/blue coloring of every pair of an
n-element vertex set with no red `K_s`, green `K_t`, or blue `K_u`.
Automorphisms preserve each named color. Write `R(v), G(v), B(v)` for the
three neighborhoods. Edges are undirected; there is no orientation gauge.

The CITED inputs are `R(3,5)=14`, `R(3,3,4)=30`, and uniqueness of the
13-vertex graph with no triangle and no independent 5-set. That graph is

    H = Cay(Z/13, {1,5,8,12}).

The first two inputs are recorded in [Boza's introduction](https://arxiv.org/html/2603.10851v1).
Completeness of the one-graph bank is attributed to
[McKay's Ramsey database](https://users.cecs.anu.edu.au/~bdm/data/ramsey.html).
Our exact checks verify the supplied graph and its cyclic presentation;
they do not independently enumerate every graph on 13 vertices.

## P1 Saturation and the directed cross margins

Let C be an `(3,3,5;56)` coloring. At every vertex,

    (d_R,d_G,d_B) = (13,13,29).

Indeed the separate upper bounds are 13,13,29, and their sum equals 55.
The green graph inside R(v), and the red graph inside G(v), are both H.
Their complementary blue graphs have degree eight and contain a K4.

Fix v and abbreviate its ports by R,G,B. For each r in R there are exactly
four red neighbors of r in G: in the green H on R(r), the vertex v has
four green neighbors, which are precisely these vertices. Exchanging the
colors gives exactly four green neighbors in R for each g in G.

Thus the R-by-G matrix has **red row sums four and green column sums four**.
The red column sums and green row sums have average four. No pointwise
equality for those latter margins is asserted.

For a red edge uv, put `t=|G(u) intersect G(v)|`. The exact third-vertex
table, with row color from u and column color from v, is

| | R | G | B |
|---|---:|---:|---:|
| R | 0 | 4 | 8 |
| G | 4 | t | 9-t |
| B | 8 | 9-t | 12+t |

Here `0<=t<=9`. Summing t over all 364 red edges counts the 26 red edges
in each green neighborhood, so its mean is four. An average is not an
edgewise intersection parameter.

## P2 An edge budget transfer formula

More generally suppose a coloring has no red or green triangle, every
vertex has degrees `(a,a,b)`, and the opposite-color graphs on the two
a-vertex ports are d-regular. Put `n=1+2a+b`. At every vertex v,

    e_R(B(v)) = e_G(B(v)) = a(b-2a+1+3d)/2,
    e_B(B(v)) = q = b(b-1)/2 - a(b-2a+1+3d).

Consequently `nq/3` is the number of blue triangles. All displayed edge
counts must be nonnegative integers and `3` must divide `nq`.

Proof: the red cross count from R to G is ad. There are
`a(a-1-d)` red edges from R to B, by the degree of every r in R.
Summing red degrees on G gives `a(a-2d)` red edges from G to B.
Now sum red degrees on B and subtract these two cross contributions.
The remainder counts every internal red edge twice. Exchange R and G,
then subtract from `binomial(b,2)`. Finally each blue triangle appears
in three blue neighborhoods. This is a transfer of exact incidences;
it assumes no uniformity of the omitted cross margins.

For `(a,b,d)=(13,29,4)`, every 29-vertex blue port has edge counts

    (e_R,e_G,e_B) = (104,104,198).

The full coloring has exactly 3696 blue triangles, and every vertex lies
in 198 of them. The mean blue codegree along blue edges is `396/29`.
It is nonintegral, so the blue graph cannot have a constant codegree on
all its edges. In particular an edge-regular blue graph cannot be a
56-vertex candidate.

Within each 29-vertex blue port, the degree caps are `(8,8,16)`, from
`R(3,4)=9` and `R(3,3,3)=17`. Define deficits
`(8-d_R,8-d_G,16-d_B)`. They sum to four at each vertex; their total
over the port is exactly `(24,24,68)`.

The GF(16) example in P6 is a genuine positive control for the general
formula: `(a,b,d)=(5,5,2)` gives internal edge counts `(5,5,0)`.

## P3 The local graph has a two-point base

`Aut(H)={x -> ax+b : a in {1,5,8,12}, b in Z/13}` has order 52.
Every nonidentity automorphism fixes at most one vertex.

Here is an elementary upper bound, independent of a group-library result.
Fix 0. Its four neighbors are N={1,5,8,12}. Outside N and 0, the vertices
with a singleton neighborhood in N are

| a in N | private vertex f(a) |
|---:|---:|
| 1 | 2 |
| 5 | 10 |
| 8 | 3 |
| 12 | 11 |

The four other vertices have N-neighborhoods
`4:{5,12}, 6:{1,5}, 7:{8,12}, 9:{1,8}`.
Thus fixing N pointwise fixes every vertex. Two members of N are called
opposite when their only common neighbor is 0. For each a, its private
vertex is adjacent to exactly one of the four double-neighborhood vertices.
That double neighborhood contains the opposite of a and one other member
of N. Select the latter as the successor of a. This intrinsic rule yields

    1 -> 8 -> 12 -> 5 -> 1.

Every automorphism fixing 0 preserves this directed cycle, so its
stabilizer has size at most four. Translations and multiplication by the
four units displayed above supply 52 automorphisms, proving equality.
An affine map fixing two distinct points has a=1 and b=0.

H is also doubly saturated for `(3,5)`: deleting any edge produces an
independent5-set, and adding any nonedge produces a triangle. More
generally, a unique Ramsey-good isomorphism class at a fixed order is
necessarily doubly saturated, since a legal one-edge change would have
a different edge count and hence a different isomorphism type. The
native finite audit sharpens this for H: all26 edge deletions produce
exactly three independent5-sets; of its52 nonedges,26 additions create
one triangle and26 create two. This connects to the older doubly
saturated Ramsey thread while keeping H13 distinct from Paley13.

## P4 One port determines the entire rooted automorphism

For C as in P1, restriction gives an **injective** map

    Aut(C)_v -> Aut(H on R(v)).

The same holds with R and G exchanged.

To prove injectivity, suppose h fixes v and fixes R pointwise. For every
r in R, the set G(r) contains four fixed vertices in R, since the green
graph on R is H. By P3, h fixes all of G(r). Every g in G has four green
neighbors in R by P1, so h fixes G pointwise as well.

For every r in R, R(r) contains the fixed vertex v and four fixed vertices
of G. Again P3 forces h to fix all of R(r). Finally, any b in B must have
a red or green neighbor in R: otherwise it extends a blue K4 in R to a
blue K5. Thus b belongs to an already fixed R(r) or G(r). All vertices
are fixed, proving injectivity.

This is faithfulness of an action, not reconstruction of the entire
coloring from the unmarked isomorphism type H. Different global colorings
can have the same local H. A root and two distinct vertices in its red
port form a base for the color-preserving group.

## P5 The full group order divides 224

Every `(3,3,5;56)` coloring has color-preserving automorphism group of
order dividing

    224 = 4*56 = 2^5*7.

First exclude an element h of order 13. Its cycles have lengths one or
13, so it fixes a vertex v. P4 implies its restrictions on both R(v) and
G(v) are nonidentity, hence 13-cycles. Its action on B(v), of size 29,
has at least three fixed vertices.

A fixed b in B(v) has constant color to each 13-cycle. Neither constant
can be blue, since each port contains a blue K4. The two constants cannot
both be red, or both green, because the relevant total degree is only 13.
They therefore use exactly 13 red and 13 green edges. All edges from b
to vertices outside the two ports are blue. Three such fixed vertices,
v, and any further vertex in B(v) now form a blue K5, a contradiction.

By P4, every point stabilizer has order dividing 52. The preceding
exclusion and Cauchy's theorem show its order divides four. If the full
group has orbits represented by v_i, orbit-stabilizer gives

    224 = sum_i 4 |Aut(C)| / |Aut(C)_{v_i}|.

The coefficient of `|Aut(C)|` is an integer. This proves the divisibility.
Every nonidentity element of a 7-subgroup is fixed-point-free, so such a
subgroup has order at most seven. No monochromatic edge set can be a
single group orbit: its size is 364,364, or812, none dividing224.

If swapping red and green is allowed, the enlarged group has order
dividing448; its color-preserving subgroup has index at most two.

The same mechanism has a uniform version. Suppose an avoiding
`(3,3,s;1+2p+b)` coloring has degrees `(p,p,b)`, p prime, and every
opposite-color p-port is d-regular with d>=2, has a two-point base for
its automorphism group at every pair of distinct vertices, and contains
a blue `K_(s-1)`. Equivalently, no nonidentity local automorphism fixes
two vertices; mere existence of one base pair would not suffice. The argument in
P4 again makes each single-port restriction faithful. If

    b mod p >= s-2, and b >= s-1,

then no point stabilizer has an element of order p: there would be at
least s-2 fixed blue neighbors, each blue to all vertices outside the
two p-ports, producing a blue K_s with v and one further blue neighbor.
If additionally p does not divide `1+2p+b`, the full group has no element
of order p. The arithmetic clock is the residual cycle count `b mod p`,
and its carrier is the actual fixed-vertex set. At `(p,b,s)=(13,29,5)`
the residue is3. At the genuine GF(16) control `(5,5,3)`, the residue is0
and the obstruction correctly does not apply.

**Boundary:** a trivial automorphism group satisfies these restrictions.
They eliminate symmetry assumptions, not all candidate colorings. In
particular, failure of a highly symmetric search cannot establish an
upper Ramsey bound.

## P6 The three distinct roles of 17

Use GF(16)=F2[z]/(z^4+z+1). The element z has multiplicative order15.
Color the edge xy by the exponent of x+y modulo3 in this multiplicative
group. This partitions the 120 edges of K16 into three 5-regular graphs.
The subgroup `{z^(3j):0<=j<5}` is a binary five-frame: its five elements
sum to zero and any four are independent. No three sum to zero, and
the same holds for each multiplicative coset. Thus every color graph is
triangle-free. This gives `R(3,3,3)>16`. The standard bound
`R(3,3)=6` makes every colored degree at most5 in a triangle-free
three-coloring, excluding17 vertices. Hence `R(3,3,3)=17`.

Each of the three color graphs is the Clebsch graph (parameters
16,5,0,2). Exact enumeration gives 1920 automorphisms after retaining
only one color and 160 after retaining all three. The repo's
THM-4624 / finite unit-circle affine-rigidity theorem and THM-4701 /
Kuratowski-Petersen unit-circle theorem supply the older five-frame and
Petersen connections; see the paths in SOURCES.md. These are classical
structures, not newly discovered Ramsey values.

Separately, the complete `(3,6;17)` bank has seven graphs. In its stored
order, the automorphism orders are `1,2,2,2,2,8,2`, and the numbers of
vertex orbits are `17,9,9,10,9,4,9`. Consequently there are **67 rooted
types**. Passing from seven unrooted graphs to a joint attachment problem
requires retaining these root distinctions or proving an operation-specific
quotient. Completeness is CITED; the rooted counts are FINITE-EXACT and
independently checked with VF2.

The ordered two-root orbit counts are `289,145,145,149,145,47,145`:
1065 in total, or998 when the two roots must be distinct. Burnside's
fixed-point-square formula and direct pair-orbit enumeration agree.
The seven local groups all have power-of-two order; this is a complete
statement about this cited finite bank, not about global90-vertex groups.

The unique `(4,4;17)` graph is instead Paley17, with 68 edges and
automorphism order136. It is not a member of the `(3,6;17)` bank.

## P7 The next boundary for R(3,3,6)

Using the CITED bounds `R(3,6)=18` and `R(3,3,5)<=57`, an avoiding
coloring on90 vertices has one of the following profiles at each vertex:

    (16,17,56), (17,16,56), (17,17,55).

Each profile occurs an even number of times, by the handshake lemma in
each color. A 17-vertex red or green port must be one of the seven graphs
in P6; smaller ports require their own bank. A 56-vertex blue port must
satisfy P1-P5. P5 bounds the intrinsic automorphism group of that induced
port; it does not assert injectivity of a surrounding 90-vertex stabilizer.

Boza's Theorem2.4 already gives `R(3,3,6)<=91`: a hypothetical91-vertex
coloring would force degrees `(17,17,56)`. Its blue neighborhoods would
each have812 blue edges; counting blue triangles would require
`3 | 91*812`, which fails. This is a CITED result, not an improvement
produced here.

Even proving `R(3,3,5)<=56` would only force the third profile at order90;
it would not by itself exclude that order. A complete exact-value proof
still requires a matching lower construction and an exhaustive upper
certificate, including asymmetric cases.
