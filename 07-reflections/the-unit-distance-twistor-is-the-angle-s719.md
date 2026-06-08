# The unit-distance momentum twistor is the angle (S719)

The last two sessions found that the lonely runner's hidden symmetry is multiplicative — the multiplier
group acting on residues — and that the logarithm that linearizes it, the discrete log, is the LRC's
momentum twistor. The natural next move was to ask the same question of the other problem the cluster
keeps circling, unit distances, and the answer is the cleanest possible: a unit edge is a difference vector
on the unit circle, the unit circle is the multiplicative group `U(1)` of rotations, the dual conformal
symmetry is global rotation, and the logarithm of `U(1)` is the angle. The unit-distance momentum twistor
is the angle map. Rotate the whole configuration and every edge angle shifts by the same constant — the
dual conformal symmetry, additive and linear in the twistor coordinate, verified exactly. It is the same
construction as the discrete log, one circle group standing in for another.

What the twistor buys, as always, is that the hard quantity decomposes. Resolve the unit distances by the
line-direction of each edge, and the count is just the sum over directions of how many edges point that
way. Each direction is a one-dimensional sub-problem — points spaced by a fixed unit vector form chains,
so a single direction holds at most `n-1` edges — and the total is their sum. Run this on the proven
champion, the twenty-one-point `K3 [] W7`, and its fifty-seven unit distances fall into exactly six
directions with multiplicities twelve, twelve, twelve, seven, seven, seven. The product structure is
sitting right there in the twistor: three directions carry seven, the triangle's edges replicated over the
wheel's seven vertices, and three carry twelve, the wheel's edges replicated over the triangle's three,
and `57 = 3*7 + 3*12`. The Cartesian product, seen through the angle, is two rosettes of directions
superimposed. A triangular blob of the same size manages only three directions and forty-nine distances.

That comparison is the whole principle. The count is large exactly when the edges crowd into few
directions, each nearly full — and few directions, each a clean copy of the others, is precisely a rotation
group, a set of roots of unity. So maximizing unit distances forces a multiplicative structure on the
directions, which is to say it forces the configuration into a CM lattice. This is the twistor's
explanation of the grid-disproof: the reason the dense configurations live in imaginary quadratic fields
with many modulus-one elements is that those elements ARE the rotation group whose orbit the edge-
directions must fill to make the angle-autocorrelation peak large. The arithmetic of the rotation group is
visible directly: in the Eisenstein layer where the unit is the norm-`D` vector, the number of unit
neighbours is the number of ways to write `D` as `a^2+ab+b^2`, and that representation count is the order
of the rotation group, whose log is the twistor `Z/r_Q(D)`. The basic triangular layer has six; the
`sqrt(13)` layer the AMP optima use has twelve. More directions, denser graph.

Pushed at the open frontier, the twistor gives a concrete gain and a clean account of the remaining gap.
Searching the Eisenstein layers for the densest twenty-two-point subset, the triangular layer gives
forty-nine, but the `sqrt(13)` layer — twelve directions instead of six — gives fifty-seven, an explicit
construction matching the product bound and far above the lattice baseline, confirming the thesis that the
higher-degree rotation group wins. The record sixty, and the open sixty-one, lie beyond any single lattice
patch, and the twistor says why in the same breath as the earlier sessions: a lattice patch of twenty-two
points cannot saturate enough directions at once, so reaching sixty needs a rigid, non-lattice graph that
holds more saturated directions than the lattice supplies — a non-product configuration, since
twenty-two factors only as two times eleven and products cap at fifty-seven, built on a rich rotation
group. The bleeding edge is exactly the place where you must leave the lattice but keep the rotation group.

And the two twistors turn out to be one. The lonely runner's twistor is the discrete log of the multiplier
group, sending the inversion to negation and the transversal core to a half-system that tiles the cyclic
group. The unit-distance twistor is the angle, the log of the rotation group, sending reflection to
negation and the extremal configuration to a set of edge-directions that fills a rotation orbit. In both,
the dual conformal symmetry is a multiplicative circle group, the twistor is its logarithm, the special
conformal generator is the negation, and the critical object is a difference set that fills a coset
structure of the group — a half-system tiling `Z/phi(2n-1)` for the runners, a rotation orbit filling the
CM module for the distances. The grid-disproof's supply of modulus-one CM elements is the unit-distance
mirror of the runner's cyclotomic shell tower. The amplitude lent us a single word twice, and both times
it named the logarithm of a circle.
