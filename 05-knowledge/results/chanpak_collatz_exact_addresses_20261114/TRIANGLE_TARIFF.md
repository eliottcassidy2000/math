# The triangle-carry obstruction to a literal Collatz lift

Take any tournament `B` with `r >= 1` nonsingleton strongly connected
components. If `T` contains `B` as an induced subtournament and
`H(T) > H(B)`, then

    H(T) - H(B) >= 2*3^(r-1).

For each r the bound is sharp. The result is an all-order constraint on
source-preserving physical constructions, not a missing-value theorem
for the unrestricted tournament spectrum.

## Proof

Use the odd-cycle-packing identity (inherited; see the repo’s [THM-002](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-002-ocf.md) and [Irving–Omar, Corollary 20](https://arxiv.org/html/2412.10572v3#S3.SS4.Thmcorollary20)): a collection of k pairwise
vertex-disjoint directed odd cycles contributes weight 2^k, including
the empty collection. The old collections contribute exactly H(B).
Consequently strict increase is equivalent to the presence of some
directed odd cycle using at least one new vertex v. In a tournament,
any vertex on a directed cycle lies on a directed triangle: shorten a
cycle v,w1,w2,... if v beats w2, and use v,w1,w2 otherwise. Thus select
one new directed triangle D, containing at least one new vertex.

If D meets at most one old nonsingleton SCC, choose any r-1 of those
SCCs disjoint from it. Each contains an old triangle. Pack D together
with any subset of the chosen disjoint old triangles. These new
collections contribute 2*3^(r-1).

Otherwise D meets exactly two old nonsingleton SCCs, at one vertex
each, and its third vertex is new. Label the old vertices x and y in
their SCC order. Then D is v→x→y→v. Within each old SCC choose an old
directed triangle through its respective vertex. Such a triangle
exists by strong connectivity and the same shortest-cycle argument.
Replace x by either the singleton x or the three-vertex directed
path beginning at x on its old triangle. Independently replace y by
either y or the three-vertex path ending at y on its old triangle.
The connecting arc between these two paths always exists, because
all cross-SCC arcs point from the x SCC to the y SCC. This gives four
distinct new odd cycles (lengths 3,5,5,7), each compatible with any
subset of r-2 triangles chosen in the remaining old nonsingleton
SCCs. Their total weight is at least 8*3^(r-2), larger than the
claimed bound. No double counting is possible: each counted
collection contains just one of these four cycles and it is the only
cycle using a new vertex.

For sharpness, take the ordinal sum of r directed triangles.
Replace the first triangle by a strong four-vertex extension with
Hamilton path count 5, e.g. rows [8,1,3,6] which contain the triangle
on vertices 0,2,3. Hamilton path counts multiply over ordinal
components, so the increase is (5-3)*3^(r-1).

A weighted, support-sensitive form follows by summing all packings in
the untouched nonsingleton SCCs. If their old counts are h_i and the
new triangle touches the index set J, the increase is at least
twice the product of h_i over i outside J for |J|≤1, or eight
times that product for |J|=2. This bound depends on which physical
SCCs the new triangle meets, not just on the scalar count.

**Equality.** A strict increase attains `2*3^(r-1)` exactly when
there is just one new directed odd cycle `D`, it is a triangle meeting
exactly one old nonsingleton SCC `C_j`, all the other nonsingleton
SCCs have `h_i=3`, and `C_j - V(D)` has no directed odd cycle. For
necessity, the zero- and two-touch cases are strictly larger, and the
one-touch bound gives `h_i=3` off the touched component. The
packings already counted exhaust the increase. A second new odd
cycle or an old odd cycle disjoint from `D` in `C_j` would give a
new packing outside that enumeration. Conversely, under the listed
conditions every new packing is precisely `D` together with an
arbitrary old packing in the other components, giving equality.

For replay run `python3 -I -B replay_triangle_tariff.py`. It directly
computes Hamilton path counts by subset-endpoint dynamic programming,
without invoking the packing formula. The computation is a control;
the proof above supplies the unbounded quantifier.

## Target-first consequence

The affine map F(h)=3h+2 is exact multiplication by 3 on the shifted
coordinate h+1, and F^k(h)=3^k(h+1)-1. It is not itself a Collatz
iteration: no dyadic division is implicit in this equality.

Suppose h=H(R)>1. Start with a disjoint ordered bank of k>=1 triangles
next to R, so the initial induced source B has H(B)=3^k h and at least
k+1 old nonsingleton SCCs (R contains a directed triangle since its
count is greater than one). No induced extension of B can give
F^k(h): the desired positive increment is 3^k-1 but the first
possible increase is at least 2*3^k. In particular freely adding
selective apex contacts, allowing several appends at once, or
burying the old components inside a new strong component does not
rescue the move. The old induced cross-arcs and old cycle packings
remain present. Changing the core, replacing an old induced arc,
or making a construction which never creates these k simultaneous
triangle SCCs is essential for escaping this particular obstruction.

## The triangle-only source also blocks the batched carry

For the bare source of r≥2 triangle SCCs, the first two positive
increment levels are separated: the increment either equals
`2*3^(r-1)`, or is at least `4*3^(r-1)`.

To see this, first suppose a new triangle touches two old SCCs,
with new vertex v. If v has a out-neighbours in the earlier SCC
and b in-neighbours in the later SCC, there are ab crossing
triangles. The originally selected triangle supplies the four
cycles counted above. There are always at least two *other* new
triangles. If ab≥3, take two other crossing triangles. If ab=2,
take the other crossing triangle and a triangle wholly within one
of the two old SCCs plus v: a or b equals 1 so that v has a
nonempty proper contact set in that old strong SCC. A nontrivial
cut of a strong SCC has crossing arcs both ways, producing the
required triangle. If ab=1 the same argument gives a separate
internal-contact triangle in each SCC. Each extra triangle gives
weight at least `2*3^(r-2)` by using old packings in untouched SCCs.
These new packs cannot repeat those already counted, since their
only cycle with a new vertex is a different triangle. Together the
contribution is at least `12*3^(r-2)=4*3^(r-1)`.

If no new triangle touches two old SCCs, then two or more new
triangles already give `4*3^(r-1)`, by counting packs with exactly
one of them and only old cycles besides. A new triangle touching
no old SCC gives `2*3^r` by itself. So below the claimed threshold
there is exactly one new triangle, touching just one old SCC.
It cannot contain two new vertices: neither could make a triangle
with two old vertices, so each has uniform contact to every old
strong SCC (a nontrivial contact cut gives a triangle). In that
case, replacing the lone old vertex on the new triangle with either
other vertex in its SCC would make another triangle.

Thus there is only one new vertex v on the unique new triangle,
and all other new vertices lie on no cycles, by the shortest-cycle
argument. Delete them without changing H. The SCC containing v
can meet only the one old SCC: if it met two, take an out-neighbour
in the earliest reachable SCC and an in-neighbour in the latest
returning SCC. They are in distinct old SCCs and make a crossing
triangle with v by ordinal orientation. This leaves a four-vertex
cyclic block containing the original triangle. It has exactly
two triangles (the old one and the unique new one) and cannot
have a longer odd cycle or disjoint pair; its count is 5. The
total count increment is therefore precisely `2*3^(r-1)`.

This also disposes of h=1: any source with `H(R)=1` and k separate
ordered triangles contains the bare `C3^k` as an induced source.
For k≥2 the desired increment `3^k-1` over that bare source is
strictly between `2*3^(k-1)` and `4*3^(k-1)`. Thus `F^k(1)`
cannot follow from retaining the k literal triangle SCCs. The k=1
case is the sharp four-vertex witness H=5 and is not excluded.
