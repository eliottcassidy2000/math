# Recovered sum/difference maps and the Petersen obstruction tests

Status: **PROVED** elementary maps and extension obstructions;
**FINITE-EXACT** constructions and hostile controls; **OPEN** universal
two-port reconstruction, the Graceful Tree Conjecture, and the requested
six-chromatic unit-distance construction. No novelty priority is claimed.

## 1. Inheritance and corrections

The closest proved mechanisms are:

- **THM-2761**, [graph-edge-sum-discriminant-codegree-factorization-and-graceful-sign-gauge](../../../../01-canon/theorems/THM-2761-graph-edge-sum-discriminant-codegree-factorization-and-graceful-sign-gauge.md):
  a bipartite sign gauge converts sums to oriented differences, while a
  separate mirror factor is required for absolute differences.
- [graceful_square_exchange_20261019.md](../../graceful_square_exchange_20261019.md):
  fixed-endpoint graceful two-edge switches are impossible; three-edge
  trades exist; exchanges among already labelled trees do not establish
  that every tree shape has a graceful labelling.
- [square_sum_ports_20261023/README.md](../../square_sum_ports_20261023/README.md):
  exact vertex-birth recurrence retains exposed endpoints and two-path
  covers. The endpoint of a Hamilton path is an essential sidecar when
  using degree-two forcing.
- **THM-4701**, [kuratowski-graphs-the-petersen-graph-and-the-finite-unit-circle-graphs](../../../../01-canon/theorems/THM-4701-kuratowski-graphs-the-petersen-graph-and-the-finite-unit-circle-graphs.md):
  the two-circle codegree bound, an induced golden unit-distance Petersen
  realization, and the pentagon two-sum model of Petersen.
- [defect_transitions_20261020/petersen.md](../../defect_transitions_20261020/petersen.md):
  cubic edge-colouring uses individually named semiedges and the three
  nonzero elements of F_2^2. Its boundary support is not a vertex-colouring
  relation or a graph-minor certificate by itself.

The canonical hostile is Petersen: it has unit geometry and a cubic core,
but its Hamilton-cycle and edge-colouring obstructions do not obstruct
vertex colouring from three-element lists. The corrected near miss is
forgetting the relation being represented: MISTAKE-643 reverses the old
meeting/skew cubic-line predicate, and MISTAKE-641 retracts an apparent
novelty claim for classical conic graphs. The least-used sidecar is the
sheet sign in a sum-to-difference conversion.

The conic identifications are already stated in the opening of Coxeter's
1983 [My Graph](https://londmathsoc.onlinelibrary.wiley.com/doi/abs/10.1112/plms/s3-46.1.117):
at q=5 the conjugacy graphs include Petersen and five disjoint triangles;
at q=7 they include the line graph of Heawood and Coxeter. That primary
page was checked. Unit-distance Petersen drawings are also classical;
the introduction of Žitnik, Horvat and Pisanski's 2012
[All generalized Petersen graphs are unit-distance graphs](https://doi.org/10.4134/JKMS.2012.49.3.475)
explicitly discusses the older twisted pentagon/pentagram drawing.
Neither fact is claimed as a new construction here.

Anchor: exact transformations which preserve the observed response.
Niche: graceful versus square-sum labelling. Wildcard: a signed pentagon
lift and the failed geometric line-graph map. The concept board retains
five coordinates: incidence predicate, endpoint/support tensor, sign sheet,
vertex versus edge colouring, and exact metric realization.

## 2. The exact incidence-matrix gauge and its boundary

Let N be the unsigned edge-by-vertex incidence matrix of a simple graph H:
row uv is e_u+e_v. Choose an orientation and let B have row e_u-e_v.
Over a field of characteristic different from two, signs sigma_v and
tau_e in{+1,-1} satisfy

    N diag(sigma)=diag(tau) B                         (1)

if and only if H is bipartite. This is a formal identity on the original
named vertex coordinates, not an equality after selecting special labels.

Proof: comparison of the two nonzero row entries forces
sigma_u=tau_e and sigma_v=-tau_e. Thus every edge joins opposite signs.
Conversely any bipartition supplies these signs. This is the exact matrix
form of the inherited THM-2761 gauge. In characteristic two, sign disappears
and N=B, so that reduction cannot recover the bipartition obstruction.

If positive labels x_u satisfy x_u+x_v=q_uv^2 on edges, the bipartite
signed labels ell_u=sigma_u x_u have absolute differences q_uv^2.
They do **not** automatically form a graceful labelling: positivity,
the bounded alphabet0,...,m, and the complete distinct difference set
1,...,m have not been preserved. The square-sum path1--3--6 becomes
labels1,-3,6, with edge lengths4 and9. No global affine rescaling turns
their ratio9/4 into the graceful two-edge ratio2.

For a candidate injection ell:V-> {0,...,m}, define oriented differences
d_e. Gracefulness is exactly

    product_(e<f)(d_e^2-d_f^2)^2 != 0.                (2)

The m nonzero absolute differences then lie in1,...,m and are distinct.
The factorization into(d_e-d_f)^2(d_e+d_f)^2 explains the missing sign
sidecar: path labels0,1,2 give oriented differences-1,+1, distinct as
oriented values but equal in absolute value. This criterion, already in
THM-2761, does not construct the injection required by the conjecture.

The native square-sum graph Q_n is not globally bipartite once n>=15:

    6--3--1--15--10--6

has sums9,4,16,25,16. The companion checks that every Q_n for n<=14 is
bipartite, so15 is the first failure. Selecting only a spanning path
removes this obstruction but also discards the remaining graph and its
future insertion choices. A gauge on that selected path is not a gauge
on Q_n.

## 3. A genuine shared transfer formula for tournaments and square sums

Let the old graph be any finite directed graph, with adjacency A. In the
squarefree algebra Z[z_v]/(z_v^2), set G=(I-ZA)^(-1)Z, e=1^tG and s=G1.
No completeness or tournament hypothesis is needed: G still enumerates
all simple directed paths with both endpoints and full support.

For two new vertices x,y, retain independent incidence vectors

    a_i=[i->x], b_i=[x->i], c_i=[i->y], d_i=[y->i],
    alpha=[x->y], beta=[y->x].

The exact number of Hamilton vertex words in the enlarged graph is

    [V]{(1+ea)(alpha+b^tGc)(1+d^ts)
         +(1+ec)(beta+d^tGa)(1+b^ts)}.                (3)

Proof: split each word into the three old runs before, between, and after
its two ports. For x before y, the middle run is either empty with weight
alpha, or is a path from an x-out-neighbour to a y-in-neighbour. The other
order gives the second term. Squarefree multiplication enforces disjoint
supports, and[V] enforces coverage. This bijection includes every empty
outer run and singleton middle run. It retains the full three-path tensor.

The tournament formula is the specialization b=1-a, d=1-c, alpha=1,
beta=0. The undirected square-sum specialization is a=b=N(x), c=d=N(y),
alpha=beta=[x+y is a square], with both directed arcs for each old edge.
For an undirected graph, reverse words pair freely since the enlarged
graph has at least two vertices; divide(3) by two for unoriented paths.
This explicitly explains the difference from tournament Hamilton parity.

For example, old labels1,3,6 form a square-sum path. Adding10 and15 gives
the five-cycle above, so(3) yields10 directed Hamilton words, or5 unoriented
paths. The earlier one-port square-sum recurrence is the corresponding
one-birth split by endpoints and two old runs. The common object is the
named support tensor; the allowed attachment rules remain different.

This gives an evaluation map, not an injectivity theorem for arbitrary
undirected graphs. Nor can the tournament's odd-count and complementary-mask
arguments be imported after replacing its relation by a symmetric graph.

## 4. A faithful twenty-point signed pentagon realization

For any symmetric unit-sum representation z_v, introduce two sheets with
points z_(v,+)=z_v and z_(v,-)=-z_v. Then

    |z_(u,+)-z_(v,-)|=|z_u+z_v|.

The abstract selected cross-sheet adjacency is the bipartite double cover
with matrix

    A_lift = [[0,A],[A,0]].                           (4)

This always preserves the chosen cross-sheet edge predicate. Faithful unit
geometry additionally requires distinct points and absence of unintended
same-sheet or diagonal unit pairs; these are separate obligations.

For the ten two-sums s_ab=zeta_5^a+zeta_5^b, a<b, the inherited unit-sum
predicate is disjointness of the two index pairs, namely Petersen. Here
all the additional obligations hold: the twenty points +/-s_ab are
distinct, and exactly30 of their190 pairs are unit distance. Those30 pairs
are exactly the opposite-sheet disjoint pairs. The companion gives every
coordinate in the integral basis1,zeta_5,zeta_5^2,zeta_5^3 and checks every
norm in Z[zeta_5], reducing by1+x+x^2+x^3+x^4. No floating tolerance is used.
This is a **FINITE-EXACT faithful unit-distance realization**.

The lift has an explicit Hamilton cycle, recorded in normal.json.
Petersen itself has no Hamilton cycle and has240 directed Hamilton paths.
Consequently this map does not preserve Hamiltonicity. It also replaces
an odd-cycle graph by a bipartite graph. The sheet and its involution let
one reconstruct the abstract quotient, but colouring or Hamiltonian
obstructions cannot be inferred from that quotient alone. This is an
exact example of information loss, not a count-based analogy.

## 5. K5 and K3,3: metric exclusions are not minor exclusions

Two distinct planar points have at most two common unit neighbours.
Subtracting the two squared-distance equations gives a line, which meets
a circle in at most two points. Hence a graph with distinct unit-distance
vertices cannot contain K_(2,3), and in particular cannot contain K_(3,3)
or K5 as a subgraph. K4 is also impossible: the two equilateral apexes
over a unit edge are sqrt(3) apart.

Use the standard Petersen edges

    01 04 05 12 16 23 27 34 38 49 57 58 68 69 79.

Its five spokes05,16,27,38,49 are disjoint connected branch sets. Contract
them and the outer cycle plus inner pentagram give every edge of K5.
Alternatively delete0, then suppress the degree-two vertices1,4,5.
The result is exactly K_(3,3) on sides{2,8,9} and{3,6,7}. Both certificates
are checked literally. Petersen nevertheless has the inherited induced
unit-distance realization in THM-4701. Thus unit-distance realizability
is not minor-closed. Planarity and unit realization require separate tests:
planar K4 already fails the metric test, while unit-distance Petersen
is nonplanar. A Kuratowski/minor certificate is not a six-colour obstruction.

## 6. Why the Petersen private core cannot block the spindle

**Lemma.** A connected subcubic graph with an induced even cycle is
colourable from every assignment of lists of size at least three.

Choose a spanning forest rooted on that cycle and colour exterior vertices
in decreasing distance from it. Each such vertex has an uncoloured parent,
so at most two neighbours are already coloured and a listed colour remains.
Each cycle vertex has at most one exterior neighbour, leaving at least two
colours. Reduce to two colours in each list. If all cycle lists are equal,
alternate the two colours. Otherwise there is an oriented cycle edge whose
initial list contains a colour missing from the terminal list. Start with
that colour and greedily follow the other direction around the cycle; at
the last vertex the starting colour is absent from its list, so closing
the cycle causes no conflict. This proves the lemma.

In the displayed Petersen labelling,0--1--2--3--8--5--0 is an induced
six-cycle. One exterior order is9,4,6,7. Therefore Petersen is3-choosable.
The construction in check.py is an explicit algorithm, tested on all
1,000 assignments obtained by independently choosing a three-subset of a
five-colour palette for each vertex class modulo3, and on the constant
three-colour assignment. The proof covers arbitrary lists, not just that
finite test family.

The fixed Moser spindle's external-contact theorem in
[the earlier geometry packet, H5-H6](../../two_port_beyond_nine_20261103/geometry/PROOF.md)
gives every private point at most two boundary neighbours. Every proper
five-colouring of that boundary consequently leaves at least three colours
at each private vertex. Thus **every such colouring extends across a
private graph which is Petersen, or any subgraph of Petersen**. The same
argument applies to any subcubic component satisfying the lemma. Extra
private edges are not silently ignored: the actual graph must satisfy the
stated hypothesis.

Petersen's snark property does not contradict this. Its six perfect
matchings have complements consisting of two five-cycles, all recorded
by exhaustive lowest-unused-vertex matching recursion. A Hamilton cycle
would have a perfect-matching complement, while a three-edge-colouring
would make the other two colours alternate on even cycles. Both are
excluded. These are edge/path predicates, distinct from vertex-list
colourability. In particular a cubic core alone is insufficient even with
ten private points; an actual failure of the boundary lists is required.

## 7. The line-graph transfer exists, but its midpoint metric map fails

Proper edge k-colouring of any graph H is exactly proper vertex k-colouring
of its line graph L(H), whose vertices are the edges of H and whose
adjacency means sharing an endpoint. This is the correct predicate-preserving
map for the Petersen edge-colouring obstruction. It changes the carrier
from ten vertices to fifteen and requires a fresh unit-distance realization.

**Midpoint obstruction.** No finite nonempty simple cubic unit-distance
graph has a unit-distance realization of its line graph obtained by taking
its edge midpoints followed by a single similarity of the plane.

Proof: if all adjacent edge midpoints had the same positive distance,
the three neighbours of each original vertex would be pairwise equidistant.
They lie on the unit circle about that vertex, so they form an equilateral
triangle whose circumcentre is the old vertex. Hence

    p_v = (1/3) sum_(u adjacent v) p_u.

On each finite connected component, every real coordinate of a harmonic
function is constant: at a maximum all neighbours must attain the same
maximum, and connectedness propagates it. Thus all p_v coincide, contrary
to the unit edges. This proves the obstruction. The argument concerns
midpoints plus one similarity only; it does not exclude an unrelated
realization of L(H).

As a finite hostile, the faithful cubic twenty-point realization above
has more than one exact squared distance between incident-edge midpoints.
Those values, multiplied by four to keep integral coefficients, are saved
in normal.json. Thus the natural edge-to-vertex map fails its metric
predicate even on a fully authenticated unit-distance carrier.

## 8. Reproduction, losses, and remaining work

Run from the repository root, replacing python3 by python3 -O for the
optimized receipts:

```text
python3 05-knowledge/results/two_port_seven_core_20261103/recovery/check.py
python3 05-knowledge/results/two_port_seven_core_20261103/recovery/transfer_check.py
```

The normal and optimized receipts agree byte-for-byte. The first checker
tests the gauge on every graph through five vertices, all190 geometric
pairs, both explicit minors, all six Petersen perfect matchings, and the
specified list controls. The second compares(3) with literal Hamilton
permutations for every digraph through two old vertices and every independent
port incidence/arc choice, and every undirected graph through three old
vertices and every undirected attachment choice. Source hashes are recorded;
all checks remain active under optimization.

The useful inherited principle is to retain products, shared labels and
the exact relation before taking scalar totals (META-PATTERNS). Here it
produces the shared two-port formula. Its counterindication is equally
concrete: a sign-sheet quotient loses odd-cycle and Hamilton predicates,
and the edge-to-vertex map loses unit geometry at the midpoint step.
Nothing here proves a graceful labelling for an arbitrary tree, a universal
two-port decoder, or a six-chromatic unit-distance graph. The surviving
geometric target needs an actual private list obstruction with every
unit contact retained; Petersen and its simple midpoint-line-graph route
are now explicitly excluded from the proposed shortcuts.
