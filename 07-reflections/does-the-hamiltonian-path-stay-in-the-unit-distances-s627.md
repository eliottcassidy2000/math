# Does the Hamiltonian path stay in the unit distances? (S627)

The user's mapping is a lovely idea: take a configuration of points, order them, build the transitive tournament on
that order, and then flip exactly the arcs whose endpoints are a unit apart. Every tournament has a Hamiltonian path,
guaranteed, so the question becomes whether that mandatory path can be made of flipped arcs — of unit distances. For
the first handful of points every pair is essentially at the unit, so of course it can. The user asks whether that
survives, and if not, where it breaks.

The first thing I learned is that the tournament itself is not quite the right object, because it remembers the order
I imposed. Reorder the points and the flips land on different positions and the path count `H` changes — I watched it
slide from three to five to eleven to thirty-seven as I added points in angular order, but those numbers belong to
the order as much as to the configuration. What does not depend on the order is the unit-distance graph: the points,
with an edge whenever two of them are exactly one apart. And the real question, stripped of the encoding, is whether
that graph has a Hamiltonian path. A unit path through all the points. That is what "the Hamiltonian path is part of
the unit distances" means.

So I asked the graph directly, for every compact triangular patch up to twenty points and for the Moser spindle, and
the answer was the same every time: yes, traceable, a unit Hamiltonian path exists. No flop. And once I saw why, the
"no" felt less like an accident. A compact patch of the triangular lattice snakes — you walk along a row on unit
edges, step up to the next row on a unit edge, walk back. The denser the unit-distance graph, the easier it is to
thread, and the optimal configurations are precisely the ones that pack in the most unit distances. Maximizing edges
is the same move as making the graph easy to traverse. So the flop is anti-correlated with optimality: a
configuration whose unit graph cannot be threaded has wasted its edges on a cut or a dangling piece, and a
configuration that wastes edges is not the maximum. The Hamiltonian path stays in the unit distances *because* the
configuration is optimal.

If the flop happens at all, it is at the place the grid conjecture just died. The record-holding configurations stop
being lattice patches and become sparse algebraic things — points projected from a high-dimensional number-field
lattice, unit distances coming from elements of norm one in a CM tower. Those graphs have degrees that grow, but
slowly, nowhere near half the vertices, and a sparse graph can refuse to be threaded. So the honest answer to "does
it always stay in the unit distances" is: yes for every lattice optimum I can compute, and the only place it could
turn over is the lattice-to-tower crossover, which lives at a size no one can see yet. The flop, if it exists, is the
same event as the grid conjecture's failure, read on the Hamiltonian path instead of the edge count.

The recursion the user asked about is hiding in the degenerate case. Put the points in a line at unit spacing and the
unit graph is just a path, and the number of independent sets weighted by two — the partition function, the
Hamiltonian-path weight — is the Jacobsthal sequence, one, three, five, eleven, twenty-one, climbing by the rule that
each term is the previous plus twice the one before. I formalized its closed form. And there in the fifth place sits
twenty-one, which is a perfectly good value for a path but a forbidden value for a tournament, because a path of
conflicts is not the conflict graph of any tournament's three-cycles. The one-dimensional chain realizes the value
the two-dimensional world forbids. The whole month has been this: the same partition function, counted along
different seams, allowed here and forbidden there, and the seams are unit distances, three-cycles, rotations, and
now Hamiltonian paths.
