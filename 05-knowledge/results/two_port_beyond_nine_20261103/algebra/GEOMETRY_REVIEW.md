# Independent review of interacting hubs and the fixed-spindle extension theorem

**Verdict: PASS, with one verifier assertion hardened before this freeze.**
This review audits [geometry/PROOF.md](../geometry/PROOF.md), its exact
certificate mechanism, and its normal/optimized replay. It is an independent
mathematical and code review, not a second independent implementation of the
whole geometry engine. No six-chromatic graph is certified.

## Reviewed mathematical claims

**H1, two adjacent hubs.** The count |A_x||A_y|-|A_x intersect A_y| subtracts
exactly the equal-colour assignments from the Cartesian product of the two
allowed lists. For nonempty lists, zero forces both to be the same singleton:
if either list has at least two elements, it has a colour different from any
chosen colour at the other hub. Empty lists are handled separately. With at
most four boundary contacts, palette5 gives nonempty lists and palette6
gives lists of size at least two. The stated zero and positive-extension
criteria have the correct quantifiers.

**H2, full boundary support.** The rational/free-angle points use the actual
unit graph, not an assumed abstract adjacency. The six ports are independent.
The two forbidden four-class partitions describe exactly the two possible
orders of the private pairs' shared residual colours. Their labelled mass
at palette5 is 2*(5)_4=240, hence 15625-240=15385 extendable assignments.
The six-distinct partition is inadmissible at palette5, not a third proper
zero. The chromatic polynomial q(q-1)(q-2)^2(q-1)^4 follows directly from
the adjacent hubs, two common neighbours and four private leaves. A supplied
three-colouring prevents any inference of global six-chromatic forcing.

**H3-H4, finite bank and growth scope.** For every boundary pair, the listed
lattice candidates achieve the full Euclidean circle-intersection count
(two, one, or zero according to squared separation below, equal to, or
above4). Since the ports are distinct, this proves completeness in the
whole plane; the search is not merely restricted to lattice points. The
resulting private graph is a path plus isolated vertices. Once the two
hubs are coloured, its path leaves and isolates always extend at the
stated palettes. The iteration theorem uses the complete norm<=4 list
in Z[alpha] and applies to sequential additions with at least two already
placed contacts. Free-angle one-contact additions and unordered simultaneous
gluings are explicitly outside its scope. Unit edges may cross; no planar-
graph edge-count theorem is being smuggled into this geometric argument.

**H5, spindle incidence and complete external contact bound.** Rotation
cosine5/6 and sine sqrt11/6 makes BB' have squared length
3*((1/6)^2+11/36)=1. Each diamond forces equal colours at its two
nonadjacent tips in a three-colouring, so BB' proves the lower bound4;
the explicit four-colouring gives the matching upper bound. The exact
21-pair incidence check is indispensable because all boundary contacts
are retained in the later argument.

For a noncollinear triple, the perpendicular-bisector equations give its
unique possible centre. The independent multiplication-only criterion is
correct: circumradius1 is equivalent to

    |d|^2 |e|^2 |d-e|^2 = 4 det(d,e)^2,

with nonzero determinant. Three distinct collinear points cannot lie on a
unit circle. The 35-triple certificate exhausts the seven-point carrier;
all ten surviving centres are old vertices. As a consistency check, these
ten triples are exactly the neighbourhood triples: four centred at the
degree-four vertex A and one at each of the other six degree-three vertices.
Thus any new exterior point with three contacts would contradict one of
the exhausted triples. The argument covers all real points, rather than
only points in the computation's number field: a noncollinear triple has
at most one centre, and the exact solution already identifies it.

The bound two is attained: B and B' are unit-separated and have no common
old neighbour in the verified graph. Their two equilateral circle centres
are therefore exterior points with two old contacts; H5 excludes a third.

The field arithmetic uses the independent basis1,sqrt3,sqrt11,sqrt33.
Its multiplication correctly uses XOR for basis labels, and factors3 and11
for repeated radicals. Inversion multiplies the three nonidentity Galois
conjugates and checks that the norm denominator is nonzero and rational.
The two circle tests have different algebraic constructions, although they
share this exact field implementation. No floating comparison is used.

**H6, every boundary five-colouring extends over <=6 private points.**
The reduction is correct: H5 bounds each private vertex's old contacts by
two, so its allowed list has at least three colours. A 2-degenerate private
graph can be coloured from these lists in reverse removal order.

The small-core classification is valid. On at most five vertices, a
minimum-degree3 core is either K4 or, at order5, a graph whose complement
is a matching. Matching size0 or1 leaves K4; size2 gives K1 joined to C4.
Four equidistant planar points are impossible, and a unit-neighbour circle
graph has no C4, since its edges step by60 degrees in disjoint regular
hexagon orbits. Thus no such private core exists.

For six points, maximum degree<=3 gives at most9 edges. A degree5 vertex
has at most4 edges among its five circle neighbours. A degree4 vertex has
at most3 edges among its four circle neighbours, and the sixth point has
at most2 contacts to them by intersection of two distinct unit circles.
Both cases again give at most9 edges. A surviving private3-core must use
all six vertices and be cubic. Its complement is C6 or two disjoint C3s,
giving respectively the triangular prism or K3,3. The latter would give
two distinct points three common unit neighbours and is impossible.

The prism list argument is sound for arbitrary three-element lists, not
only a palette of five. Its first triangle has at least6 proper choices.
After matched colours are deleted, the three opposite lists have size at
least2. Hall can fail only when all three are the same two-element set S.
In that event each original opposite list is S plus its matched colour;
the matched colours are distinct. Thus the original opposite lists are
pairwise distinct, S is their pairwise intersection, and each matched
colour is uniquely determined. At most one of the >=6 first-triangle
choices fails; the other choices extend. Restricting larger lists to
three-element sublists only strengthens this existence proof.

The conclusion is therefore valid for any six or fewer distinct new
planar points with every induced unit edge included. It is specific to
the fixed seven spindle ports. A nonempty private3-core and at least seven
private points are necessary for blocking a boundary assignment on this
carrier, not sufficient for a new obstruction.

## Verifier adjustment and replay

The initial source asserted only `maximum <= 3` after constructing the
external-centre list, while the theorem states maximum2. Its earlier
per-triple assertion already required every unit centre to be an old
vertex, and frozen artifacts were empty externally, so this did not
invalidate the result. At this reviewer's request, the author hardened
the final assertion to `not outside_contacts and maximum == 2` and
refreshed the receipt. No mathematical repair was necessary.

I replayed `python3 verify.py` and `python3 -O verify.py` in the geometry
directory after that change. Both returned byte-identical output matching
the frozen receipt as JSON. The replay verified all203 boundary partitions,
all256 bank assemblies, all35 spindle triples with zero exterior triple
centres, and100000 normalized prism list choices. The latter cover10^6
unnormalized three-list choices over a five-colour palette by colour
permutation; the proof, rather than this finite palette check, establishes
the unrestricted three-list theorem. Every returned list witness is checked
against its lists and every prism edge.

Reviewed SHA256 hashes:

- geometry/PROOF.md: `1563cfaf2ef99c0b45b48378315aee61b4a78df9856a5142ca46f2258b00e166`
- geometry/verify.py: `090e00f093f75929c880f622e4bab16ecc8d238d2826ce58e461c25656346a4c`
- geometry/spindle_contacts.json: `a530071f6dbcf6491ea85870235a9d7930d71f55f8e929f1bc3a2b11013bfd74`
- geometry/receipt.json: `7d2515822ea7f001cbaf44e32f5483d52d97aeffc2d1a52e70a53d3392dbb4c4`

The complete boundary relations, fixed geometry and palette are retained.
No conditional zero is promoted to an empty global colouring relation.
