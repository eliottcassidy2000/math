# Seven nine-edge pieces: extend the Moore 7-to-21 construction to 63

2026-11-27. Every heptagon in a hypothetical degree-57 Moore graph
forces two canonical 21-cycles by the [previous construction](../moore_57_heptagon_21_residual_cuts_20261127/README.md).
Each of those 21-cycles extends, with all its old vertices still
in order, to a **simple 63-cycle in seven nine-edge pieces**.
The first step is canonical; the second has many choices.
The result applies to every degree-57 candidate, not just those
with a particular symmetry. It is conditional, as [degree 57 remains open](https://www.mdpi.com/2227-7390/14/23/4318).

## Why the recipe is not literal tripling

The earlier step used the unique short routes from each vertex to
the rest of an induced heptagon. Its output has chords, so there
is no claim that the same uniqueness argument applies to its
21 vertices. Also, replacing each of the 21 edges by a three-edge
path would make a square with that original edge, impossible.
Instead keep one edge in each three-edge piece, and replace the
other two by four-edge detours. The result has 7(4+4+1)=63 edges.

In either canonical 21-cycle, orient each three-edge piece from
base to double-contact to special single-contact to base, as in

    v -- w -- z -- v'.

Replace v--w and w--z and keep z--v'. The proof is that all 14
chosen original edges have simultaneous detours with disjoint
internal vertices outside the entire 21-vertex set.

## Every edge has a two-coordinate detour grid

For any edge uv in a diameter-two Moore graph of degree d, let
R=N(u)-{v} and S=N(v)-{u}. They are disjoint. A choice a in R,
b in S forces a unique common neighbor x of a,b: a and b are
not adjacent (that would make a square with uv). The common
neighbor lies outside the closed neighborhoods of u and v by
the triangle and square rules. Conversely every x outside
those closed neighborhoods has one contact in R and in S,
so the middle x identifies the cell (a,b). Thus the (d-1)²
cells are exactly the four-edge detours u--a--x--b--v
forming a pentagon with uv.

For a forbidden set F not containing u or v, write a_F and b_F
for the numbers of forbidden sites in R and S, respectively.
The number of usable detours is at least

    (d-1-a_F)(d-1-b_F) - (|F|-a_F-b_F).

In fact the first product counts the remaining rectangle, and
one then deletes just those cells in it whose middle vertices
are forbidden. This is the [earlier edge grid](../nested_moore_shells_edge_grid_20261126/README.md)
used as a pool of paths rather than an adjacency representation.

## Simultaneous choices and the 63-cycle count

Initially F holds all 19 other original vertices. The induced
21-vertex scaffold has degrees 5 at v, 3 at w, and at most 4
at z. After j detours, add the 3j new internal sites to F.
Each earlier detour has three internal sites on a two-edge path.
Any outside vertex is adjacent to at most one of those three:
two consecutive neighbors make a triangle, and two extreme
neighbors make a square. Consequently a_F,b_F rise by at most
j each, not 3j. On a v--w edge they are bounded by 4+j and
2+j; on a w--z edge they are bounded by 2+j and 3+j.
The rectangular lower bound is decreasing in either forbidden
degree over this range. Using the worse of the two edge types,
there are always at least

    f_j = (52-j)(54-j) - (13+j),   0<=j<=13

choices for the next detour. The last bound is still 1,573.
Choose an edge order once; every completion gives a different
63-cycle for this fixed marked 21-cycle. This proves the
extension for both canonical 21-cycles of every heptagon.

For unmarked counting, each 63-cycle has at most 18 preimages:
two orientations and nine choices of phase for the seven
equally spaced original base vertices. Once these are fixed,
the w and z marks are four and eight steps after each base
mark. The original heptagon on the seven base vertices is
unique if it exists, since no 7-cycle has a chord; the original
21-cycle and all 14 detours are then specified. This proves

    # simple 63-cycles >= ceil(246486240000 * product(f_j) / 18)
                       > 5.22 * 10^56.

The 246,486,240,000 distinct canonical seeds (twice the heptagon
count) were established in the preceding note. These are
ordinary cycles, not induced: among other chords every
stretched original edge is still an edge of the ambient
graph. This is not a statement about tournament path counts.

## The nearby lengths 61 and 65 can use the same seed

If all replacements are four edges, we can in fact replace
any subset of the 21 original edges: at stage j the lower
bound using degrees 5 at both ends is (52-j)²-(11+j),
positive even at j=20. This gives lengths 21,24,...,84.

To fill the two other residue classes, choose at most two
edges first for **five-edge** paths. In the grid of an edge
uv, the residual vertices (the middle x sites) have exactly
d-2 neighbors inside the grid, because each has exactly one
contact on each side. These neighbors fall in distinct columns
of the grid: two in the same column would already share their
column contact as well as x. At either of the first two
stages there are at most 23 forbidden sites. Choose an
unforbidden x with an unblocked row contact a. There are at
least 55-23-23=9 choices of a residual neighbor y with an
unblocked column contact b and y not forbidden. This gives
u--a--x--y--b--v. Each earlier four-site internal path blocks
at most two side-neighbor choices later, again by the
triangle and square rules.

For e<=2 five-edge paths followed by j four-edge paths,
the latter grid has at least (56-4-2e-j)² usable side
pairs before deleting at most 19+4e+3j forbidden middle
vertices. The difference stays positive for e+j<=20.
Thus every length 21+4e+3r with r+e<=21 occurs through the
same marked seed, in particular all lengths 27 through 86.
Length 61 takes one five-edge and twelve four-edge paths;
65 takes two and twelve. The length 63 construction above
retains the stronger equal-sized seven-piece structure.

The replay checks the grid and deletion count on all 175
edges of Hoffman–Singleton and builds a literal 42-cycle
from one canonical 21-cycle by stretching one edge per
piece. Its 50 vertices of course cannot host a 63-cycle.
The degree-57 part of the replay checks the exact arithmetic
rather than claiming to construct a degree-57 graph.
Run node replay.js.
