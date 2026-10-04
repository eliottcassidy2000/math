# Pair tiles, medial dualities, and the diagonal that remembers vertices

2026-10-04 (America/Denver).

**Status:** PROVED elementary reconstruction and counting statements;
FINITE-EXACT for the listed graph and tournament censuses; CITED classical
map background. The arrows named `F(G)` and `E(G)` in the owner's diagram
have no recovered definitions in the supplied text or targeted repository
search. They are not assigned invented meanings here. We explicitly name
the flag graph `Flag(G)` and the line graph `L(G)` when using them.

## Inherited mechanisms and scope

The closest proved route is
[the tournament/polyhedron note](camion_busch_gaps_polyhedra_collatz_20261002.md),
sections 6.1--6.4: the 24 transitive labelled four-vertex tournaments are
tetrahedral flags, the six tetrahedral edges are the vertices of its
octahedral medial, and a symmetry subgroup must be distinguished from the
full automorphism group. Its P7 example explicitly has different face
orbits under the tournament subgroup and under the full map group.

The [fifteen-position tiling note](tiling_modular_atoms_20261003.md),
sections 1--2, provides the exact loop/path/free-pair coordinate map and
warns that presentations are not isomorphism classes. The relevant canon is
[THM-1430, tiling/class/metagraph dictionary](../../01-canon/theorems/THM-1430-the-tiling-class-metagraph-dictionary-and-which-tricks-pay.md)
and [THM-4524, selfie loop gauge](../../01-canon/theorems/THM-4524-selfie-tournaments-loop-gauge-arc-parity-and-the-odd-mallows-sloane-partner.md).

The canonical hostile below is a tetrahedron: its abstract map group has
order 24, while its medial graph has 48 automorphisms. The corrected near
miss is identifying an operation's entire automorphism group with the
group transported from its input. The least-used sidecar is a cell's type:
diagonal double, vertex-star, face-boundary, or flag-adjacency color.

The live board consists of ordered pair addresses, diagonal doubles,
tournament orientation bits, embedded faces, typed flags, and symmetry
groups. The useful unification is a faithful carrier with explicit
projections, not equality of quantities because they happen to be 4, 6,
10, or 100.

For primary classical context, Hubard, del Rio Francos, Orbanic and Pisanski,
[Medial symmetry type graphs (2013)](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v20i3p29/pdf/),
sections 2--2.3, define maps, typed flag adjacency, duality, and medial
operations. They identify map automorphisms with color-preserving flag
automorphisms, and medial-map symmetries with original automorphisms and
dualities. The restricted elementary reconstructions used here are
proved below; no novelty is claimed for this map theory.

## 1. The diagonal and two triangles are a branched pair quotient

Let X be an n-element set. Transpose sends `(i,j)` to `(j,i)` in `X x X`.
Its fixed set is the diagonal `Delta={(i,i)}`. The quotient is the set
`Sym^2(X)` of two-element multisets, consisting of ordinary unordered
pairs and diagonal doubles. Thus

    n^2 = n + 2*C(n,2),
    |Sym^2(X)| = n + C(n,2) = n(n+1)/2.

Each off-diagonal quotient cell has two ordered preimages; a diagonal cell
has one. A chosen linear order places the first sheet above the diagonal
and the second below. The two triangular pictures therefore require an
ordering convention, although the quotient and its diagonal are invariant
under arbitrary vertex relabelings.

A tournament chooses exactly one member of each off-diagonal two-element
fiber. Its `C(n,2)` choices are orientation bits. Diagonal cells can carry
optional loop bits, but are not ordinary tournament edges. Reading a
diagonal address arithmetically produces a square `i^2` in a product table
or a double `2i` in a sum table; these are labels on the same address type,
not new graph edges.

**Loss control:** `1*4=2*2=4` identifies an ordinary pair and a diagonal
double after retaining only the product. The coordinate addresses preserve
the distinction. Similarly, folding an ordered pair loses its orientation
unless the sheet bit is retained. Neither numerical multiplication nor
the chosen upper-triangle convention is invariant under every relabeling
in `S_n`.

## 2. An intrinsic completion theorem for the dotted diagonal

Define `J(n,2)` on ordinary pairs: two vertices are adjacent exactly when
their supports intersect. This is the line graph `L(K_n)`. Define `D_n`
on all of `Sym^2(X)` by the same support-overlap rule. The diagonal doubles
are vertices of `D_n`, not graph loops of `D_n`.

**PROVED.** For n>=2, the uncolored abstract graph `D_n` intrinsically
recovers X and all its pair addresses; `Aut(D_n)=S_n`.

Proof: a double `{i,i}` has degree n-1; an ordinary pair `{i,j}` has degree
2n-2. The distinct degrees identify the diagonal without an external
color. The ordinary pair `{i,j}` has exactly the two diagonal neighbors
`{i,i}` and `{j,j}`. Hence a graph automorphism's permutation of diagonal
vertices determines its action on every vertex. Conversely every
permutation of X preserves support intersections. This proves both the
reconstruction and the exact automorphism group.

At n=4, the six off-diagonal tiles form an octahedron, `J(4,2)`. Its
nonedges are the three complementary-pair pairs, so its automorphism group
has order `2^3*3!=48`. Vertex relabeling contributes only `4!=24`.
The extra central involution is

    e -> X\e.

It sends the three edges through vertex 0 to the triangle on vertices
1,2,3. Thus it swaps a vertex-star and a face-boundary; it is not induced
by a permutation of the original four vertices. Adding the four diagonal
doubles gives `D_4`, with ten vertices and automorphism group of order 24.
The dotted diagonal is exactly sufficient memory to remove this extra
symmetry in the pair construction.

For n>=5 the off-diagonal graph already recovers X. A pairwise-intersecting
family of two-element sets either has a common element or lies in one
three-element triangle: start with `{a,b},{a,c}`, and a set not containing
a must be `{b,c}`. Therefore the maximal cliques of size n-1>=4 are
precisely the n vertex-stars. Their pairwise intersections recover the
original pairs, proving `Aut(J(n,2))=S_n`. At n=4 both stars and triangles
have size three, explaining the exceptional ambiguity. At n=3 the graph
is `K_3`, again with group `S_3`; at n=2 the single off-diagonal vertex
forgets the swap of its two endpoints.

## 3. Five original vertices give the requested 3, 4, 5, 10, and 100

For X of size 5 there are ten ordinary pairs. They lie in five vertex-stars
of size four, and each pair belongs to exactly two stars. A fixed pair
has six overlapping neighbors and three disjoint neighbors. The disjoint
graph is the standard pair-set model `KG(5,2)` of the Petersen graph;
the overlap graph `J(5,2)` is its complement.

Consequently the 10-by-10 table indexed by ordered pairs of ordinary
pair-tiles has the intrinsic three-part decomposition

    100 = 10 equal + 60 overlap + 30 disjoint.

After removing equality and folding transposition, its 45 distinct
unordered entries split into 30 overlaps and 15 disjoint pairs. These
relations are exactly the three `S_5` orbits on ordered pairs of tiles:
their union has size 2, 3, or 4, and bijections between equal-size unions
extend to relabelings of X.

This table is a table of **pair addresses and relations**. Identifying its
ten row labels with residues modulo 10 requires a selected bijection.
Nothing here identifies its relation classes with the values of the
modulo-10 multiplication table. Such a comparison must retain that labeling
and check an additional arithmetic predicate.

An ordinary five-vertex tournament places ten orientation bits on these
ten pairs. Forcing the path `5 -> 4 -> 3 -> 2 -> 1` fixes four, leaving
the six nonconsecutive pairs. Restoring the five diagonal loop positions
gives exactly

    15 positions = 6 free pairs + 4 path pairs + 5 doubles.

The explicit triangular coordinate bijection is already proved in the
tiling note. It transports addresses and these roles, but does not assert
that Euclidean tile adjacency equals support intersection. A path is a
retained gauge: arbitrary vertex relabelings need not preserve the selected
path presentation. The 64 binary presentations in this gauge represent
12 ordinary tournament isomorphism classes, not 64 classes or six classes.

## 4. What dual and medial graphs preserve

Now let G be a connected simple polyhedral map: vertices, edges, cyclic
face boundaries, and their incidences are part of the data. The dual
interchanges vertices and faces and keeps the edge carrier. The medial
has original edges as vertices and original corners `(v,f)` as edges.
Two original edges are adjacent in it when they meet consecutively at a
face corner. Its faces come in two remembered kinds: the vertex-stars
of G and the face-boundaries of G.

The same corner construction is symmetric in vertex and face roles,
so `M(G)=M(G*)` on the identified edge carrier. This equality concerns
maps and incidence; it does not identify their original vertices with
their original faces.

**PROVED reconstruction.** A medial map with its two face types retained
recovers G: the original vertices are the vertex-type medial faces, the
original edges are the medial vertices, and their endpoints are the two
incident vertex-type faces. The other face type recovers the original
face boundaries. Therefore

    Aut(G) = Aut(M(G), remembered face types)

as map automorphism groups. In a connected medial map its checkerboard
face coloring has exactly two choices. Forgetting the types permits an
automorphism either to preserve them or to exchange them. The latter
operation is a duality of G. Thus the embedded medial automorphism group
contains the original group with index one or two, and index two occurs
exactly when G is self-dual.

This statement is about the **embedded map**. For a naked abstract graph,
one must first check whether an automorphism preserves the relevant
embedding. For a geometric realization, one must separately check whether
an abstract map automorphism is an actual Euclidean symmetry. Neither
upgrade is included by changing terminology.

The smallest convex-polyhedral hostile has four vertices. Take a tetrahedron
with coordinates `(0,0,0),(1,0,0),(0,2,0),(0,0,3)`. Its six squared edge
lengths are all distinct: `1,4,5,9,10,13`, so any Euclidean symmetry fixes
every edge and hence every vertex. The three group orders are

    geometric tetrahedron: 1,
    abstract tetrahedral map: 24,
    abstract octahedral medial graph: 48.

The regular tetrahedron has geometric group 24, but the medial graph
still has 48 abstract automorphisms. One may choose a subgroup and compare
its vertex/edge/face orbits faithfully; allowing all symmetries of each
derived graph independently changes the equivalence relation.

**Line/medial boundary.** For the tetrahedron, `L(K_4)=M(K_4)`, because
degree three makes every two incident edges consecutive around a vertex.
The same equality holds for the cubic cube. It fails already for the
five-vertex square pyramid: opposite apex edges `{0,4}` and `{2,4}` meet
but do not share a face. Its line graph has 18 edges and its medial has
16. In particular `J(5,2)=L(K_5)` is degree six, whereas a medial map is
degree four. The n=5 pair construction must not be renamed `M(K_5)`;
K5 is also not a spherical polyhedral graph.

## 5. Typed flags give a faithful common carrier

A flag is an incident triple `(v,e,f)`. There are four flags per edge.
Define involutions `r_0,r_1,r_2` by changing only the vertex, edge, or face,
respectively. Their colored adjacency graph is `Flag(G)`; it is cubic,
and `r_0 r_2=r_2 r_0` records the four flags around one edge.

**PROVED.** Recover the original cells as components after retaining two
colors:

    vertices: colors {1,2},
    edges:    colors {0,2},
    faces:    colors {0,1}.

The components are exactly all flags containing a fixed cell; incidence
is nonempty intersection of components of different types. Hence this
carrier recovers the entire map, and its color-preserving automorphism
group is exactly the original map group. A duality exchanges colors 0
and 2. This supplies a precise meaning for unifying vertex, edge, and
face symmetry without declaring those cell types identical.

For a tetrahedron, a permutation `(a,b,c,d)` gives the flag

    ({a}, {a,b}, {a,b,c}).

The three flag moves swap adjacent positions in this permutation. The
same permutation is the transitive tournament with order `a,b,c,d`.
Reversing one arc keeps a transitive tournament precisely when it swaps
adjacent positions. Thus the flag graph is exactly the graph on the
24 transitive labelled tournaments joined by single-arc reversals; it
has 36 edges. This is the inherited truncated-octahedral connection.

All four ordinary tournament classes on four vertices are retained in the
audit:

| Kind | Sorted outdegrees | Labelled tournaments |
|---|---|---:|
| Transitive | 0,1,2,3 | 24 |
| Directed triangle plus source | 1,1,1,3 | 8 |
| Directed triangle plus sink | 0,2,2,2 | 8 |
| Strong | 1,1,2,2 | 24 |

The flag bijection selects the first class. It does not turn the four
isomorphism classes into the four vertices of a tetrahedron or identify
all 64 tournaments with its 24 flags. The existing oriented-dual construction
relates the other classes, but requires its face/orientation convention.

## Exact controls and the surviving connection

Run:

    python -X utf8 -B 04-computation/experiments/medial_pair_tiles_20261004.py
    python -X utf8 -B -O 04-computation/experiments/medial_pair_tiles_20261004.py

The [saved output](medial_pair_tiles_20261004.out) records exhaustive
unfiltered adjacency-preserving automorphism searches: diagonal-completed
pair graphs for n=2,...,6; strict pair graphs for n=3,...,6; and typed
flags/medials of tetrahedron, cube, and square pyramid. Their map/medial
group orders are `(24,48)`, `(48,48)`, `(8,16)`. Every tested abstract
medial automorphism is checked against the supplied face sets before it is
classified as type-preserving or type-exchanging.

Further controls recover V/E/F components, compare dual medials on the
same edge carrier, audit the asymmetric geometric tetrahedron, enumerate
all 64 four-vertex tournaments and all 64 five-vertex path presentations,
verify the explicit transitive-flag bijection, and count every cell in the
100-entry pair relation table. The line/medial and arithmetic-value
hostiles reject the tempting stronger identifications. All checks remain
active under optimized Python.

The surviving chain is exact: ordered pairs -> unordered pair tiles with
an orientation sheet -> diagonal-completed incidence -> selected graph
operators with remembered cell types -> typed flags. Each arrow states its
discarded information. The last two unnamed arrows of the original diagram
remain undefined until their operations are supplied.
