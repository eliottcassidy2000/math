# The two tetrahedral Fano planes are dual, not equivariantly identical

Session23b. Proved in the stated scope, with independent review in the sibling session23b-duality-review directory. This note joins two already reviewed constructions with an explicit map; it is not a new definition of the Fano plane or a metric consequence of seven cospherical labels.

Let K4 be the tetrahedron's edge graph. The **edge Fano plane** has as points its six edges and a formal point omega. Its lines are the four face-edge triangles and the three triples consisting of an opposite-edge pair together with omega. The **cycle Fano plane** has as points the seven nonzero elements of the binary cycle space: four face triangles and three Hamiltonian quadrilaterals. Its lines are triples of cycle supports {A,B,A XOR B}.

There is a canonical S4-equivariant bijection from cycle points to EDGE-PLANE LINES:

- A triangular face cycle maps to the same three-edge face line.
- A Hamiltonian quadrilateral maps to the line consisting of the two omitted opposite edges and omega.

This is a projective duality: a triple of cycle points is a line in the cycle Fano plane if and only if the three corresponding edge-plane lines meet at one point.

Proof. The three quadrilateral cycles XOR to zero, so form one cycle-plane line; their images are the three opposite-pair lines, concurrent at omega. Two distinct triangular faces share exactly one tetrahedral edge e. Their XOR is the quadrilateral omitting e and its opposite edge. The corresponding three edge-plane lines all pass through e. The six unordered face pairs give the six remaining cycle-plane lines and the six remaining concurrence points, one for each tetrahedral edge. These seven cases exhaust the lines in both planes. Every definition uses only faces, opposite edges and symmetric difference, so the map commutes with all tetrahedron vertex permutations.

This is not an S4-equivariant bijection between the two POINT sets. In the edge point action, omega is fixed and the six edges form one orbit, giving orbit sizes1+6. In the cycle point action, the face cycles form an orbit of size4 and quadrilaterals an orbit of size3. An equivariant point bijection would preserve fixed points and orbit sizes, which is impossible. An arbitrary seven-element bijection exists, but does not retain the natural tetrahedral action.

An explicit coordinate version recovers the same map. Use edge points03=1,13=2,23=4,01=3,02=5,12=6 and omega=7 in F2^3. For a nonzero dual vector a define the cycle support

C_a={edge e: a dot (v_e+omega)=1}.

The four a with a dot omega=1 give the four face cycles; the remaining three give the quadrilaterals. The corresponding edge-plane line is L_a={v nonzero:a dot v=0}. The map a->C_a is linear after adjoining the zero cycle, and C_a XOR C_b=C_(a+b). Thus the cycle space is explicitly the dual vector space in this marked model. The complement in the quadrilateral rule is essential: simply identifying a cycle's support with a Fano line would give four-point 'lines' and fail.

## Relation to the tournament work

The edge-plane construction is the carrier for the seven independent cyclic Fano-line orientations and its Q7 line-flip graph. The cycle-plane construction records the supports of directed cycles in a tetrahedral tournament and then lifts them to integral circulation cones. They are genuinely related by duality, but their additional data differ: a line orientation on the former is not automatically a feasible nonnegative circulation or a directed-cycle support on the latter. The strong order4 tournament's three directed-cycle supports are a NONCOLLINEAR triple, whereas the omitted XOR cycle is not directed.

This resolves an exact compatibility question, not an analogy based on the number7. It provides a canonical bridge while proving the obstruction to replacing it by an equivariant point-to-point identification. Neither construction identifies the user's metric sphere features with the seven points, and no unit-distance or clock structure is supplied.
