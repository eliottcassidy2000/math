# Six-colour forcing: the boundary, stress, and geometry must all survive

**PROVED:** a planar framework with independent pinned rigidity rows
extends every proper four-colouring of its boundary. Four is sharp.
A minimal failure to extend a prescribed six-colouring carries at least
`h+ceil(a/2)` independent pinned stresses, where h is its number of
exterior vertices and a its number of port edges. The exact sharper
bound below retains repeated forbidden colours.

**PROVED:** the smallest blocker banks for the inherited palette-doubling
method have exactly two forms. Geometry eliminates one form when the
ports have distinct colours. **FINITE-EXACT:** the stated small controls.
Neither result supplies a six-chromatic unit-distance graph or settles
Hadwiger–Nelson. The recent primary
[Alon–Bucic–Davies introduction](https://arxiv.org/html/2609.36220v1#S1)
still states the Euclidean bounds `5<=chi(R²)<=7`; its theorem about
typical norms does not apply to the Euclidean norm.

## 1. Inherit the actual obstruction

The closest proved mechanism is
[THM-4614, three-colourable contexts](../../../01-canon/theorems/THM-4614-three-colourable-contexts-cannot-change-the-six-colour-port-relation-in-q-zeta55.md).
At its **fixed eight ports** in `Q(zeta55)`, there are only nine possible
exterior double-port vertices and no exterior triple-port vertex. Every
proper six-colouring extends across every three-colourable exterior.
The port-0 exception and the residue-prime scope were repaired in
MISTAKE-593; neither is discarded here.

The later [interface profile](../six_colour_interface_profile_20261015f.md)
and [component-palette matching](../six_colour_component_pairs_20261016a.md)
retain which exterior classes touch the ports and their independent
two-class component flips. They give stronger extension tests for some
four-colourable exteriors. Merely escaping three-colourability is not a
forcing certificate.

The least-used relevant sidecar is the **pinned** rigidity stress space.
The [unit-distance carry packet](../unit_distance_carry_20261019/README.md)
already proves that an unpinned stress-free framework is four-colourable
and that first arithmetic carries must annihilate every stress. Section 2
below strengthens its colouring conclusion to arbitrary prescribed
boundary colourings and quantifies the cost of boundary forcing.

The canonical hostile is changing the consumer: a forbidden boundary
assignment is much weaker than a large ordinary chromatic number.
[THM-4661, odd-degree carriers and equality testers](../../../01-canon/theorems/THM-4661-odd-degree-carriers-a-failed-six-colour-exterior-has-at-least-twenty-vertices-and-equality-testers-are-seven-colour-certificates.md)
proves that a nonvacuous six-colour equality tester would already imply
`chi(R²)=7`. That is a different, stronger target than finding a
six-chromatic unit graph. Its boundary and geometry assumptions matter.

Concept board: fixed boundary relations; paired palettes; pinned rigidity;
carry compatibility; literal tournament moves. Anchor: a valid geometric
forcing criterion. Niche: equality in the palette union bound. Wildcard:
whether quotient returns retain enough data to lift actual unit edges.

## 2. Independent pinned rows preserve every four-colouring

Let G be a finite simple graph faithfully placed at distinct points
`x_v` of the Euclidean plane. Unit lengths are allowed but are not needed
for this theorem. Mark a set B of ports, and put `H=V(G)-B`.
The **pinned rigidity matrix** `R_B` has two columns for each v in H
and one row for every edge with at least one endpoint in H. Its row for
uv is the differential of squared edge length divided by two:

    (x_u-x_v) dot (h_u-h_v),       h_b=0 for b in B.       (1)

Edges internal to B have no row. A pinned stress is a vector in the left
kernel of `R_B`; equilibrium is imposed at exterior vertices, while the
ports may carry reaction forces. This is not the unpinned stress space.

**Theorem S1 (PROVED).** If the rows of `R_B` are independent, every
proper four-colouring of `G[B]` extends to G. More precisely, every
nonempty `U subset H` has a vertex v with

    deg_(G[U])(v) + |N_G(v) intersect B| <= 3.            (2)

**Proof.** Let u=|U|, let e count edges internal to U, and let a count
edges from U to B. The corresponding e+a rows remain independent after
discarding the unused exterior columns, so `e+a<=2u`.
If a>0, the sum of the left sides of (2) is

    2e+a <= 4u-a < 4u.

Thus one left side is at most three. If a=0 and u>=2, the two translations
and infinitesimal rotation are independent kernel vectors on the distinct
points of U. Therefore `e<=2u-3`, giving the same strict average bound.
If u=1 and a=0, its degree is zero. Delete a vertex satisfying (2), recurse,
and colour the deleted vertices in reverse order. At most three already
coloured neighbours forbid colours at each step. QED.

The same proof works over a field of odd characteristic with a
nondegenerate planar quadratic form and distinct positions. The
infinitesimal rotation is the corresponding invertible skew-adjoint
linear map. This permits exact good-prime rank certificates for algebraic
Euclidean frameworks; a modular dependence alone would not prove a
characteristic-zero dependence.

**Sharpness.** The Moser spindle with a single vertex pinned has eleven
independent rigidity rows, but no three-colouring. Its classical exact
coordinates, inherited from [the spindle interface note](../spindle_tournament_interfaces_20261011g.md),
are independently checked here in `Q(s,t)`, `s²=-3,t²=-11`. Reduction at
17 sends `s->4theta,t->6theta,theta²=3` and certifies rank eleven after
one vertex is pinned. Thus four cannot be replaced by three, even with a
nonempty boundary. The graph has explicit four-colourings.

This is an extension theorem, not a claim that all four-colourable
frameworks have independent rows. Its contrapositive is the usable gate:
any actual failure to extend four or more prescribed colours requires a
nonzero pinned stress. For more than four colours, use the same at-most-three
forbidden-neighbours greedy step; do not recolour or restrict the prescribed
boundary to four colours.

## 3. A quantitative stress cost for a failed six-colouring

Fix a proper q-colouring beta of the ports, with q>=1. Let H be a vertex-minimal
nonempty exterior on which beta fails to extend. Write

    h=|H|,    e=|E(H)|,
    a_v=|N(v) intersect B|,
    t_v=|beta(N(v) intersect B)|,
    a=sum a_v,    t=sum t_v.

The available list at v has size `q-t_v`. Every proper induced subgraph
of H is list-colourable, so

    deg_H(v)>=q-t_v,       2e>=qh-t.                     (3)

Let `sigma_B=dim ker(R_B^T)` on precisely this exterior and its incident
port edges. Since `rank R_B<=2h`,

    sigma_B >= e+a-2h
             >= ceil(((q-4)h+2a-t)/2).                 (4)

For q=6, this gives the two useful inequalities

    sigma_B >= h+ceil((2a-t)/2) >= h+ceil(a/2).           (5)

Repeated port colours decrease t and strengthen the bound. If only one
geometric port touches H, rotation about it adds one to the right side
of (4). If no port touches H, translations and rotation add three.
These refinements use actual pinned positions, not the number of colours
seen at the boundary.

Equation (5) is relevant to an arithmetic lift with prescribed port lifts
held fixed: every one of these pinned stresses must annihilate the first
carry right-hand side, including those fixed port data. If the ports may
move, the applicable matrix and stress space change. Independent pinned
rows cannot evade the colouring test because they already guarantee
extension. Passing the stress equations only removes a first-order
obstruction; higher carries, global realization, and colour forcing
remain separate duties.

**A sharp hostile.** Put one exterior vertex at 0 and take the six ports

    p_j=((1-j²)/(1+j²), 2j/(1+j²)),       j=1,...,6.       (6)

They are distinct rational points on the unit circle. No two ports are a
unit distance apart: otherwise they and 0 would form an equilateral
triangle with rational coordinates, impossible because its area is
`sqrt(3)/4`, whereas a rational-coordinate triangle has rational area.
The full unit graph is exactly `K_(1,6)`.

Give the ports six different colours. The central list is empty, so this
prescribed six-colouring does not extend. Here `h=1,a=t=6`, the pinned
matrix has rank two, and its stress dimension is four, attaining (5).
Yet the whole graph is two-colourable. Thus nonzero or even large pinned
stress is necessary for the relevant failure, never sufficient for a
high chromatic number. Changing the ports can manufacture trivial
boundary forcing while evading the actual inherited question.

The unit diamond with its two tips pinned to different colours similarly
attains (4) at q=3: `h=2,a=t=4`, five pinned rows, rank four, one stress.
Both interiors are forced to the same remaining colour although adjacent.

## 4. Classify equality in the palette-doubling bound

Let k>=3, take a palette of n=2k colours, and fix a proper k-colouring of
the exterior. Assume each exterior vertex has at most two port neighbours.
An option partitions the n colours into k labelled pairs, one pair per
exterior class. A **blocker type** `(i,{a,b})` says that some vertex in
class i sees the two distinct boundary colours a,b. Duplicate types have
the same effect and are counted once. Let F_i be class i's pair family.

There are `n!/2^k` options, and every blocker removes exactly the fraction
`1/binom(n,2)`. The inherited union bound therefore guarantees an option
with fewer than `binom(n,2)` blockers.

**Theorem P1 (PROVED).** A bank of exactly `binom(2k,2)` distinct blocker
types covers every option if and only if it has one of these forms:

1. one class forbids every colour pair, all other classes forbid none;
2. for one colour c, every class forbids all `2k-1` pairs containing c.

At k=3 these are exactly nine labelled banks: three concentrated banks
and six common-star banks. The two shapes were already controls in the
original doubling note; their exhaustive classification is the addition.

**Proof.** At this threshold the sum of blocker sizes equals the number
of options. Covering everything forces a disjoint cover. Two blocker
types in different classes overlap exactly when their pairs are disjoint:
those two pairs can be completed to a labelled matching. Consequently
every two nonempty families F_i,F_j must be cross-intersecting.

For two nonempty cross-intersecting edge families A,B on n>=6 vertices,

    |A|+|B| <= 2(n-1).                                  (7)

Here is an elementary proof. A singleton family has at most `2n-3`
partners. If neither family is a singleton and one contains disjoint
edges, the other has at most four edges, the four pairs joining their
endpoints. If that other family also has disjoint edges, both have at
most four. Otherwise it is a star or triangle: a two-edge star has
`n` common intersecting edges, a star of three or four edges has `n-1`,
and a triangle has three. All resulting sums are at most `2(n-1)`.
If neither family contains disjoint edges, each lies in a star or a
triangle, and each has size at most n-1. This proves (7).

If exactly r>=2 families are nonempty, sum (7) over their pairs to get
`sum |F_i|<=r(n-1)`. The required total is `k(n-1)`, so r=k and equality
holds in every pair inequality. Since k>=3, every family has size n-1.
A family with disjoint edges would give its partner at most four members,
contradicting n-1>=5. Thus every family is a full star. Cross-intersection
forces their centres to agree. The case r=1 is the concentrated bank.
Conversely a concentrated bank blocks the assigned pair of its class,
while a common-star bank blocks the unique pair containing its centre.
In each case exactly one blocker hits each option. QED.

**Palette failure is not colouring failure.** A concentrated bank is
realized abstractly by independent exterior vertices, one for each pair
of prescribed port colours. Each vertex can choose one of the other
colours, so the whole boundary colouring extends although this palette
compiler fails. Conversely, the common-star bank at six colours can be
made into an actual failed extension: use `K_(5,5,5)` outside six ports;
the j-th vertex of each part sees ports of colours 0 and j, j=1,...,5.
Each part needs at least two of the five colours other than 0, and
different parts need disjoint sets. This is impossible. The whole graph
is nevertheless four-colourable. Neither construction is claimed to be
a faithful unit-distance realization.

## 5. Geometry rejects the common-star equality case

**Corollary P2 (PROVED).** In a faithful planar unit-distance graph,
suppose there are exactly 2k ports with 2k distinct prescribed colours,
k>=3, and every exterior vertex has at most two port neighbours. If the
exterior is k-colourable and has at most `binom(2k,2)` vertices with two
port neighbours, a failure of the palette method at the threshold can
only be concentrated in one class, with exactly one blocker for each
port pair.

Indeed a common-star bank would require k different exterior vertices
adjacent to the same pair of ports, one in each exterior class. Two
distinct Euclidean points have at most two common unit neighbours, so
k>=3 is impossible. Vertices with more than two port neighbours must be
excluded as in P1; the rational star of (6) shows why this guard matters.

At six colours a useful sufficient extension test follows: with six
distinct-coloured ports, no exterior triple-port vertices, and at most
15 double-port vertices, **any** supplied proper three-colouring which
places blockers in at least two classes gives a valid pair-palette lift.
Fewer than 15 blockers is the inherited unconditional case. For the fixed
eight Q(zeta55) ports the stronger nine-point census already wins, so
this corollary is intended for changed interfaces, not a relabelling of
that old result. Repeated boundary colours or more ports do not inherit
the two-common-neighbour argument for colour pairs.

## 6. What transports to tournament moves, and what does not

An orientation of every pair on n vertices forgets to `K_n` regardless of
its Hamiltonian paths, blue reversals, or K moves. Thus the ordinary
chromatic number after this forgetful map is always n. For example the
transitive tournament and the strong four-vertex tournament both forget
to `K_4`, while the former is pure-blue and the latter is mixed. Therefore
no invariant of that undirected graph alone can detect the colour of a
metagraph move or supply a unit-distance forcing certificate.

The useful common structure is a certificate interface, with different
consumers. Tournament detours retain the full directed path and exact
labelled endpoint; boundary colour compilers retain the prescribed port
partition and lift; arithmetic geometry retains the actual coordinates
and stress/carry pairing. There is no identified graph isomorphism or
predicate-preserving map between these three objects.

For a complement-edge game on an undirected graph, toggling a pair can
change a colouring constraint, but its legality as a unit edge needs
coordinates. A formal return or an Euler edge count alone supplies no
such lift: the inherited straight p-cycle in a residue unit graph has
the expected closed walk, yet its first carry has nonzero trace and
forbids a section at precision two. A `(p,q)` torus model should retain
that kind of return/lift coordinate before transporting its trajectory
length to a Euclidean forcing claim. No particular game rules are
silently assigned here.

The most useful next search is therefore a **fixed-interface** candidate
which survives the palette screens, together with its full pinned stress
space and compatible carry data. Stress-free families are now excluded
even as boundary-forcing devices at four colours. Stressful candidates
still require an exact full boundary relation and faithful geometric
realization. None of these criteria supplies a six-chromatic unit graph.

## 7. The incoming nineteen-point gate attains the stress bound

The latest [geometric palette gate](../kernels_games_detours_20261020f/forcing.md),
commit411c1fd35a, has the seven-vertex Moser spindle as its interior H and
twelve ports. Nine ports form three matching covers of the six noncentral
vertices; three more ports are private unit neighbours of the central
vertex. Its exact contact certificate gives

```
h=7,  e=11,  a=t=21,  every exterior vertex sees three distinct ports.
```

For the displayed proper six-colour boundary assignment, the three forbidden
colours at each interior vertex are3,4,5. Every list is{0,1,2}; the spindle
has no such colouring. Deleting any interior vertex permits a three-colouring,
so this is a vertex-minimal failed exterior, as required in section3.
The full nineteen-point unit graph nevertheless has chromatic number four.

The stress count is exact, not merely bounded. At any interior vertex,
three distinct ports lie on its unit circle. Their difference vectors
cannot all be collinear, because a line through the centre meets that
circle in at most two points. Select two nonparallel port rows per vertex.
Those fourteen rows form a block-diagonal nonsingular minor of the pinned
matrix, so its rank is14. Port-port edges do not contribute rows. Hence

```
sigma_B = (11+21)-14 = 18
        = h+ceil(a/2).
```

This attains the six-colour stress lower bound with a four-chromatic
interior, strengthening the one-vertex rational-star sharpness example.
The incoming proof separately shows that removing any one of the three
central ports leaves an eighteen-point graph whose entire proper
six-colour boundary relation is full. Both claims retain the actual
boundary graph, which contains two edges, rather than treating the ports
as an independent set.

[incoming_gate_stress.py](incoming_gate_stress.py) imports the pinned exact
coordinate/contact receipt, checks seven nonzero port determinants,
certifies all seven vertex-deletion colourings, and replays the numerical
stress equality. It does not duplicate the incoming global geometry census.
The remaining construction problem is still a geometric palette synchronizer
that supplies the prescribed input relation; identifying equally coloured
ports would change the geometry.

## 8. Reproduction and finite scope

[six_colour_check.py](six_colour_check.py) uses no repository modules.
It checks exact rational star distances and rank; every distance of the
seven-point spindle in the biquadratic field; its faithful residue graph;
all 128 port subsets; and all proper four-colour boundary assignments
for the 26 subsets with modularly certified independent pinned rows.
All 737 such assignments extend. There are 384 labelled four-colourings
of the spindle and no three-colouring. The single-pin spindle and the
pinned diamond give the two sharp controls.

It also examines all 32,768 families of pairs on six colours, recovering
the inherited 36 extremal first families for (7) and exactly six
five-by-five star families. All nine minimum banks partition the 90
options. The concentrated extension and abstract tripartite failure
are independently checked; the latter's explicit four-colouring and K4
certify its ordinary chromatic number four.

    python 05-knowledge/results/strong_blue_games_20261020f/six_colour_check.py
    python -O 05-knowledge/results/strong_blue_games_20261020f/six_colour_check.py

These finite controls audit the stated proofs, not Hadwiger–Nelson or
all-order black connectivity. META-PATTERNS used: retain operation and
carrier; recompute essentiality per consumer; preserve the selected
boundary; compose witness fibres before scalar projection. No new
universal research rule is promoted from this one lane.
