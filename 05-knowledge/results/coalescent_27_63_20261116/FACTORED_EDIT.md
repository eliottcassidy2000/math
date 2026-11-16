# Two repairs of the three-triangle 27

The exact minimum number of old arcs to reverse when the source is
C3→C3→C3 and arbitrary new vertices are allowed is **two**, compared
with the five needed at fixed order 9. Both possible H63 strong cores
have eight directed triangles, and no two are vertex-disjoint.
Neither does ordinal padding add triangles. One flip in the old
graph would leave two disjoint old triangles, impossible in the
target. This proves the lower bound, regardless of new vertex budget.

For attainment, reverse old edges 0–1 and 3–4 (in the canonical
old graph of edit.cpp). The first two old triangles become an ordered
six-vertex transitive chain. Only the last triangle remains.
One embedding of that altered source is obtained by putting six
singleton sources ahead of the H63 core
(16,5,17,7,74,31,175,63) and identifying old labels 6,7,8
with core vertices 4,6,7 (in that order). This simple witness adds
five vertices. A smaller witness uses the first H63 core
(2,16,3,23,165,79,31,111), restricting it on old vertices to
indices {0,1,3,4,6} and prepending four singleton sources.
That restriction has two leading singleton SCCs followed by a
directed triangle. It is exactly the rewritten old graph up to
relabeling, so adds just three vertices.

**At the two-reversal optimum, three new vertices are necessary.**
Enumerate the 9-vertex restrictions of the two target cores with
arbitrary neutral padding. Such a restriction consists of a core
subset with k≤8 old vertices, with 9−k singleton old vertices
split arbitrarily before and after. This gives 76 distinct patterns
up to isomorphism. Search two edits in the fixed old source and
canonically match the 9-vertex results, tracking the largest matching
k, which is 5. Consequently at least 8−5=3 new core vertices are
needed. The program factored_edit.cpp performs precisely this
finite pattern search, independent of any bound on extra padding.

The [triangle-first code](../triangular_backbone_chanpak_20261116/CORE.md)
already assigns three independent digits to the factored 27. Every
odd cycle in a count-63 core shares a vertex with every other such
cycle, since two disjoint cycles would give two disjoint directed
triangles. Its odd-cycle-pack polynomial in a marker z is therefore
1+31z (eight 3-cycles, seventeen 5-cycles, six 7-cycles); the old one
is (1+z)^3. Both evaluate at z=2 to the desired counts, but changing
the packing rank is compulsory. This is an interface-level
explanation of the edit barrier rather than a modular one.
