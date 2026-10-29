# Reflection traces, geometric carriers, and encoding capacity

**PROVED elementary deductions / CITED geometric counting bound / FINITE-EXACT
examples. Same-agent audit.** Local labels L1-L2 and G1-G4 belong to this file.

## L1. Exactly what the natural LRC reflection detects

Let V be a nonempty finite set of positive integer relative velocities and
0<delta<=1/2. Define the good-time set on R/Z by

    G_delta(V) = {t : ||v t|| >= delta for every v in V}.

It is a finite union of closed intervals and isolated points, excludes 0,
and is invariant under r(t)=1-t. Because it excludes 0, its connected
components can be viewed as closed intervals inside (0,1).

**Theorem.** Reflection fixes exactly one component if every v is odd,
and fixes none otherwise. Therefore

    trace(r on H_0(G_delta(V); Q)) = [all v odd],
    number of connected components mod 2 = [all v odd].

For the good residues on the grid t=k/q, 0<=k<q,

    number of good residues mod 2 = [q even and all v odd].

**Proof.** A closed component [a,b] is fixed exactly when a+b=1,
and then it contains 1/2. At most one component does. Conversely the
component containing 1/2 must be reflection-invariant. The half-time is
good exactly when all v are odd, because its phase distances are 1/2 or 0.
All other components occur in pairs. The trace on H_0 is the number of
fixed components. On the finite grid the only reflection-fixed residues
are 0 and, when q is even, q/2. The first is never good and the second has
the same all-odd criterion. QED.

This is an exact obstruction to a particular parity strategy, not an
obstruction to LRC itself. For V=(1,2), the good sets are:

| delta | Good set | Component trace |
|---:|---|---:|
| 1/3 | {1/3,2/3} | 0 |
| 1/4 | [1/4,3/8] union [5/8,3/4] | 0 |
| 2/5 | empty | 0 |

The same trace cannot distinguish empty, isolated, and positive-length
witness sets. For V=(1,3,5), delta=1/4, the good set consists of the
points 1/4,3/4 and the interval [9/20,11/20], giving three components
and trace one. The all-odd case already has the explicit witness t=1/2.

At the LRC(14) boundary example V={1,...,13}, delta=1/14, exact interval
intersection gives six isolated good times {1,3,5,9,11,13}/14.
This is a verified tight example, not an exhaustive thirteen-speed test.

The formulas hold on the unit-period circle even when gcd(V)>1.
For an intrinsic minimal-period interpretation, first divide by gcd(V).
Dilation changes the number of repeated components, so period normalization
must remain part of the representation.

## L2. Moving the reflection centre while preserving phases does not help

A reflection r_a(t)=a-t preserves every individual function ||v t|| for all t
if and only if v a is an integer for every v in V.

Necessity follows by setting t=0. Sufficiency follows from
||v(a-t)||=||-v t|| when va is integral. Thus, for primitive V (gcd=1),
the only such reflection on R/Z is a=0. This statement preserves each
labelled distance function; a reflection allowed to permute functions
or preserve only a coarser union is a different problem.

Together L1-L2 explain why the self-converse tournament analogy does not
automatically produce an odd LRC witness count. Tournament pair classes
have extra independent converse and mask operations (P2-P3); the phase
predicate-preserving LRC reflection has no analogous freedom to change
its centre. A new parity proof would need a different object or additional
weights, boundaries, or incidences with a proved nonzero invariant.

## G1. Odd-prime theta carriers retain an arithmetic mark obstruction

Let p be an odd prime and k>=2. Form D from p internally disjoint directed
u-to-w paths of length k and p-1 isolated vertices. Then

    N=pk+1, E=pk=N-1.

Rotate the p paths by g of order p, fixing u,w and the isolates.
There are p+1 fixed vertices and no fixed carrier edges. Every free
pair orbit has a consistent orientation, so g fixes 2^a completions, with

    a = k(pk-1)/2 + (p^2-1)/2.

Any C_p-set X with pk vertices has f fixed vertices, where f is a
nonnegative multiple of p. On Tour(X), the fixed count is 2^b, where

    b = k(pk-1)/2 + (p-1)f(f-1)/(2p).

Equality would require f(f-1)=p(p+1). Its only nonnegative integer solution
is f=p+1, incompatible with p dividing f. Therefore even a weak equivariant
trijection is impossible, already upon restriction to this subgroup.

This family and its nonweak obstruction are **inherited** from THM-4835
[Euler criterion and tiers](../../../01-canon/theorems/THM-4835-euler-characteristic-criterion-for-marked-excess-minus-one-carriers-in-tournaments.md).
The calculation is restated to compare the same carriers with geometry;
it is not a new family of trijection counterexamples.

For p=3,k=3, the source fixed exponent is 16. The possible target exponents
are 12,14,22,36, corresponding to f=0,3,6,9. For p=3,k=2, source exponent
9 contrasts with target exponents 5,7,15.

## G2. Length three admits faithful unit-distance realizations for any number of paths

Here “faithful” means that edges are exactly the pairs at distance one;
vertices must be distinct. Crossings are allowed: this is not a planar
drawing requirement. These are the definitions used by
[Alon and Kupavskii](https://www.tau.ac.il/~nogaa/PDFS/ramsey10.pdf).

The p length-two paths form K_(2,p), which cannot be realized in R2 for
p>=3: two distinct centres have at most two intersections of their unit
circles. Adding isolated vertices cannot repair this obstruction.

In contrast, p length-three paths between two terminals have a faithful
unit-distance realization in R2 for every finite positive integer p.
Any prescribed finite number of isolated vertices can be added.

**Construction and proof.** Put u=(0,0), w=(2,0). For rational
t in the open interval (1/4,1/2), set

    a(t)=((1-t^2)/(1+t^2), 2t/(1+t^2)).

Then 3/5<a_x<15/17 and 0<|w-a|<2. Choose the upper of the two intersections
b(t) of the unit circles centred at a(t) and w, continuously on this interval.
Thus u-a-b-w is a unit path, and neither terminal pair nor a-w has length
one. Also b-u cannot have length one: the unit circles about u,w meet only
at (1,0), and |a-(1,0)|=1 would require a_x=1/2, outside the interval.

The point a varies injectively along an arc of the circle centred at u.
Every fixed b has at most two preimages a, since a must lie on the two
unit circles centred at u and b. The centres are distinct (b lies on the
unit circle about w and cannot be u). Thus b(t) has finite fibres.

Suppose finitely many paths have already been chosen. Coincidence of a or b
with an old vertex excludes only finitely many parameters. An unwanted
unit distance from a to an old vertex other than u excludes at most two a
positions; from b to an old vertex other than w it excludes at most two b
positions and hence finitely many parameters. The exceptions u-a and w-b
are the desired terminal edges. The other terminal distances were excluded
above. Rational parameters remain after deleting this finite set.
Inductively choose p paths. Their internal vertices are distinct and have
no unwanted unit distances. Add isolates away from the finitely many points
and unit circles. QED.

This construction gives algebraic coordinates when rational t is chosen.
The exact p=3 example uses t=1/4,1/3,1/2 (the endpoint values also satisfy
the required inequalities), followed by isolates (10,10),(20,20).
The [receipt](unit_distance.json) gives all ten coordinates and checks all
45 unordered pairs symbolically, including every nonedge. It has nine
edges and therefore belongs to the requested E=V-1 class.

Combining with G1 gives infinitely many faithfully unit-distance-realizable
carriers which are not even weakly trijective (all odd primes p, k=3).

## G3. The converse implication fails too, with a rigid eight-vertex witness

Take the undirected graph K_(2,3) on centres 0,1 and middle vertices 2,3,4,
add the component edge 5-6, and add isolate 7. Orient the arcs as

    2->0, 2->1, 0->3, 3->1, 0->4, 1->4, 5->6.

There are eight vertices and seven arcs. The (outdegree,indegree) pairs
for vertices 0 through 7 are respectively

    (2,1), (1,2), (2,0), (1,1), (0,2), (1,0), (0,1), (0,0).

They are all distinct, so the oriented automorphism group is trivial.
The free-coordinate cube has dimension binom(7,2); any coordinate bijection
to pairs of a seven-set is automatically equivariant. This is a monomial
trijection. Yet its K_(2,3) subgraph forbids a unit-distance realization
by the two-circle argument. No minimality among all possible examples
is asserted.

Consequently, geometric realizability and trijectivity are independent even
inside the E=V-1 universe. Euler characteristic alone controls neither
metric realization nor the fixed-subgraph marks needed for trijections.

A separate exact hostile repairs a historical degree assumption: place a
vertex at the origin and seven distinct rational points on its unit circle.
The rational parametrization above with t=0,1,...,6 supplies them.
The origin then has seven unit neighbors. The degree-six bound requires
a separation or triangular-lattice restriction; it is false for arbitrary
planar point sets. The same construction gives arbitrarily high degree.

## G4. A geometric analogue of sparse encoding capacity

Alon and Kupavskii's Theorem 1.2 gives, for fixed dimension d, the number
of labelled faithful unit-distance graphs on m vertices as
2^((1+o(1)) d m log2(m)). In dimension two it is at most
2^((2+o(1)) m log2(m)).
This concerns the abstract distance relation, not the continuum of its
coordinate realizations. [Primary source](https://www.tau.ac.il/~nogaa/PDFS/ramsey10.pdf).

Suppose an injective encoding sends every n-vertex tournament isomorphism
class to an unlabelled faithful planar unit-distance graph with exactly
m=m(n) vertices, and the decoder uses only the abstract graph. There are
at least 2^binom(n,2)/n! input classes. The labelled graph count bounds
the number of unlabelled output classes above, hence

    (2+o(1)) m log2(m) >= binom(n,2)-log2(n!)
                       = n^2/2-O(n log n).

It follows that

    m log2(m) >= (1/4-o(1)) n^2,
    m >= (1/8-o(1)) n^2/log2(n).

For the second inequality, if m>=n^2 it is immediate; otherwise
log2(m)<=2log2(n), which gives the displayed lower bound.
In particular m=O(n) is impossible. Fixed-length tuples multiply the
available leading exponent by their length, and recorded sidecars must
be charged separately.

This parallels the inherited sparse E=V-1 capacity obstruction in
[Session81](../session81_counterexample_families_20261028/GRAPH_CARRIERS.md)
and the sharper concurrent
[Session90 capacity theorem](../session90_even_graphs_odd_descent_20261029/capacity/CAPACITY_AND_EULER_REPAIR.md).
The underlying counting mechanisms differ: few-edge subsets versus
realisable distance relations. The obstruction concerns the same output
contract, namely lossless encoding of all isomorphism classes.

If the decoder may read arbitrary real coordinates, labels, or an unbounded
auxiliary string, this graph-only argument no longer applies. It also does
not constrain a representation designed only to answer one predicate.
No new bound for the maximum number of unit distances, no coloring theorem,
and no complexity-theoretic collapse follows.

## A targeted next frontier

The most concrete next LRC question is a weighted incidence object that
survives the three (1,2) controls above and the six-point AP13 boundary.
It must specify its operation, cancellation rule, nonzero invariant, and
map back to an actual good time. Merely replacing reflection by another
centre while retaining all phase functions is blocked by L2.

For tournament composition, compare signed spectra before searching for
a pair-coordinate codec, then retain joint stabilizer data when attaching
ports. For geometry, varying path length changes realization while leaving
the odd-prime mark obstruction alive. These are useful experiments because
their success and failure mechanisms are independent and exactly stated.
They do not identify LRC, unit distances, and tournament enumeration as
one extremal problem.
