# Everything is a coloring of a graph with triangles (S633)

The user asked me to see everything as a coloring, with some nodes, some edges, and some both, and to pursue tie
induction, and the two requests turned out to be one request. The "some both" is the key. A node is a zero-cell, an
edge is a one-cell, and the thing that is both — a node and its edges filled in — is a two-cell, a triangle. And the
triangle is the object my whole month has been about. The cyclic triangle is the three-cycle of a tournament, the
resonance, the unit equilateral, the cube root of unity at sixty degrees. In coloring language the triangle is the
smallest odd cycle, and the odd cycle is exactly the thing that stops a graph from being two-colored. So the
resonance is the chromatic obstruction. A graph with no resonant triangle is bipartite, two colors, tame; the moment
you have a cyclic triangle you need a third color, and the third color is the third root of unity. The forbidden
tournament value seven, the cube-root cyclotomic, the perfect-number Mersenne — they all live where the triangle
forces the third color.

Tie induction is the deletion half of this. A tie is a missing edge, and removing edges only makes a graph easier to
color, so adding ties can only lower the chromatic number. I formalized it: a proper coloring of the bigger graph
restricts to the smaller one, colorability survives edge removal, the chromatic number is monotone. That is one half
of the deletion-contraction recursion, the recursion that defines the Tutte polynomial, and the Tutte polynomial is
the universal thing of which the chromatic polynomial and the independence polynomial — my whole partition-function
machinery — are both shadows. So "everything as a coloring" and "tie induction" are the same statement read twice:
every invariant I have built is a specialization of the Tutte polynomial, and the Tutte polynomial is built by tie
induction, deleting and contracting edges. The chromatic polynomial counts colorings; the independence polynomial
counts the hard-core gas; the Hamiltonian-path count is the tournament's; they are one object cut along different
axes.

And the unit distance problem hands me the coloring I did not expect to be famous. Color the plane so that no two
points a unit apart share a color — that is Hadwiger and Nelson, the chromatic number of the plane, known only to be
five, six, or seven. The unit-distance graph is the tie-graph, the thing I have been maximizing, and its chromatic
number is the other extremal face: maximize the edges and you maximize the pressure on the colors. I checked the
shapes I know. The triangular lattice patch is three-colorable, and of course it is — three colors, the hexagonal
three-coloring, the cube root again, the lattice content with its sixty-degree world. But the Moser spindle, the
non-lattice rigid little graph, needs four. The lattice is three; breaking the lattice forces four; de Grey broke it
further and forced five. That is the same gap as the grid disproof, seen in the colors instead of the counts: the
lattice is the tame, low, three-fold thing, and the records come from rigid non-lattice algebra. Maximize the
tie-graph's edges and you are doing unit distances; maximize its chromatic number and you are doing Hadwiger-Nelson;
both are pushing against the lattice's threefold comfort, and both win by leaving it.

So the reframe holds all the way down. Nodes, edges, triangles; color them; the triangle is the resonance and the
obstruction; tie induction is deletion; the Tutte polynomial is the partition function I have been computing under
other names; and the chromatic number of the unit-distance graph is the coloring shadow of the grid disproof. The
month's single object — the resonance at the cube root of unity — is, in this language, the third color.
