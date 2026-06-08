# Erdős problem 98 is the unit-distance dual; the obstruction is the shared "3" (S722)

Erdős's ninety-eighth problem asks a question that sounds unrelated to anything the cluster has done, and
turns out to be the mirror image of all of it. Take `n` points in the plane in honest general position —
no three on a line, no four on a circle — and count the distinct distances they determine; let `h(n)` be
the smallest such count. Does it grow faster than linearly? Erdős could not even prove it reaches `n`. The
best constructions, by Pach and by Erdős–Füredi–Pach, get the number of distinct distances down to
`n^{log_2 3}` and then to `n exp(c sqrt(log n))`, subquadratic, almost linear. Whether the general-position
constraint forces superlinearly many distinct distances is open.

The repo has spent six sessions on the opposite question. Unit distances ask how *many* times a single
distance can repeat — the peak of the point set's radial autocorrelation. Problem ninety-eight asks how
*few* distinct distances there can be — the support of the same autocorrelation. One maximizes a peak, the
other minimizes the support; they are two extremizations of one object, the Patterson function of the
configuration. That alone would make them cousins. What makes them the same problem in two costumes is the
number three.

Here is the mechanism. If no four of the points are concyclic, then no point can have four others at a
single distance from it, because four points equidistant from a point lie on a circle centered there, and
that is four concyclic points. So every distance, viewed from every point, has multiplicity at most three;
every distance-graph has maximum degree at most three; each distance is used at most `3n/2` times; and
since the `binom(n,2)` pairs must be shared among the distinct distances, there are at least `(n-1)/3` of
them. That `(n-1)/3` is the only lower bound anyone has, and the three in it is not a coincidence. It is
the same three the cluster met last week at the unit-distance frontier: THM-440 and S719 found that the
densest twenty-one-point unit-distance graph admits a new vertex joined to at most *three* of its points,
because you cannot place four of them on a unit circle, and that ceiling of three is exactly why `u(22)`
stalls at sixty. The local statement "no four of these on the unit circle" and the global statement "no
four of these on any circle" are the same ban at one radius and at all radii. Problem ninety-eight is the
all-radii form of the obstruction that caps `u(22)`.

Once you see that, the reason ninety-eight is hard for our toolkit becomes visible too, and it is a clean
tension rather than a vague difficulty. The way to make distinct distances few is to use a CM lattice,
where the squared distances are norms in an imaginary quadratic field and there are very few of them — this
is the whole content of S719 and of the grid-disproof. But a CM lattice is the most concyclic object
there is: a point has `r_Q(D)` neighbours on each norm-`D` circle, six on the nearest, twelve on the
`sqrt(13)` layer the optima love, always at least four. Every circle that makes the distance set small
also makes four points concyclic, which problem ninety-eight forbids outright. So the structure that
minimizes distinct distances is exactly the structure the problem bans. The extremal configurations for
ninety-eight have to be CM lattices with their concyclicity surgically removed — few distances kept, every
four-on-a-circle destroyed — which is what the Erdős–Füredi–Pach construction does and why it is delicate.
The unit-distance optimum and the ninety-eight optimum are the two poles of one autocorrelation: one wants
all the concyclicity it can get, the other none.

Three reframes the repo already owns sharpen the attack. As a coloring problem — the cluster's standing
lens — minimizing distinct distances is decomposing the complete graph into the fewest distance-classes
that are each at most three-regular and simultaneously drawable in the plane with no three points
collinear; the three-regular packing gives the `(n-1)/3` floor, and the planar realizability is the whole
difficulty. As a Sidon problem — the additive-combinatorics face of our LRC twistor and half-systems — "no
four concyclic" is a ban on coincidences, and the gap between the trivial `(n-1)/3` and Erdős's hoped `n`
is precisely the gap between the weak consequence the bound uses (each point sees each distance at most
three times) and the full ban it ignores (no four on *any* circle, including circles centered away from the
points). And as a temperature problem — last week's frame — distinct distances are the support of the
radial measure, concyclic CM lattices are the cold crystalline configurations with minimal support, and
the question `h(n)/n -> infinity` is asking whether forbidding the cold structure forces the support to
warm.

That is where to begin work in earnest. The lower bound throws away the strongest hypothesis. The unused
strength of "no four concyclic" is a global ban on circle-incidences, which is a ban on isosceles
four-tuples and perpendicular-bisector coincidences, which is a bound on the additive energy of the
distance structure. The route is to bound that energy using the zero-concyclic hypothesis and then convert
energy to support by Cauchy–Schwarz, the Elekes move, to push `h(n)` off `(n-1)/3` toward superlinear. The
autocorrelation operator from S714 is the bookkeeping: distinct distances are the size of the support of
its radial spectrum, and the concyclic ban is a bound on its off-diagonal circle-incidences. The cluster
built the operator for the peak; ninety-eight is the same operator read for the support, and the same three
that stalls the peak at twenty-two is the three that floors the support at `n/3`. The problem is to spend
the hypothesis the floor leaves on the table.
