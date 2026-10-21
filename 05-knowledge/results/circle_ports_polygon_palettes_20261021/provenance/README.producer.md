# Circle phases, physical ports, and a six-colour forcing obstruction

2026-10-21. **PROVED scoped elementary statements; FINITE-EXACT controls.**
No new chromatic lower bound and no claim of literature priority. The
Hadwiger–Nelson and LRC(14) problems remain open. All radius statements below
refer to the Euclidean unit of distance, not to a normalized drawing.

## 1. Inheritance and the exact next question

The closest inherited mechanism is
`05-knowledge/results/zero_heavy_robustness_20261021a/palettes.md`, especially
its section 5: to add a centre and force a sixth colour, every proper
five-colouring of the source must use all five colours on the centre's
actual neighbours. Its root-only gates have arbitrarily large conditional
repair radius while their whole unit graphs remain four-colourable.

The canonical hostile is therefore a prescribed boundary obstruction read
as an ordinary chromatic obstruction. The corrected near miss is
`05-knowledge/results/strong_blue_games_20261020f/six_colour.md`: pinned
stress is necessary for its boundary-extension failure, not sufficient for
large chromatic number. Its six-coloured leaf-star is already a warning.
The least-used sidecar here is the full set of available colours at every
physical circle vertex, rather than one palette cardinality.

The board is: **radius; angular phase; unit-edge equality; off-circle
contacts; prescribed ports; runner clock**. Anchor: a necessary geometric
condition for circle saturation. Niche: exact phased regular polygons.
Wildcard: a single nonresonant radius destroys nonlinear clock transport.

## 2. Exact graph on one circle

Let C_R be the circle of radius R centred at o. For R<1/2 there are no
unit edges. At R=1/2 unit edges form antipodal pairs. For R>1/2 set

    alpha_R = arcsin(1/(2R))/pi,       0<alpha_R<1/2.

In angular turns, two vertices are a unit distance apart precisely when
their difference is +alpha_R or -alpha_R modulo 1.

**Circle classification.** If alpha_R=p/q is reduced, every connected
component of the whole-circle unit graph is a q-cycle. A finite subset
induces paths and possibly complete q-cycles. If alpha_R is irrational,
the whole components are double rays, and every finite subset is a forest
of paths. Consequently every finite circle unit graph is at most
three-colourable; odd cycles are the only obstruction to two colours.

Proof: each angular orbit of addition by alpha_R is a cyclic group of
order q, or an infinite copy of Z. The unit edges are exactly consecutive
steps in that orbit. No other edges exist. The q=2 case was separated to
avoid double-counting antipodal edges.

For the particular radius **R=1**, alpha_R=1/6. Every component of a
finite circle graph is a path or a C6. Thus it is bipartite, even if the
vertex set contains very many overlapping regular polygons. Adding the
centre uses at most a third colour. At any other radius the centre has no
unit edges to the circle, so circle plus centre is still at most
three-colourable. An odd cycle on another radius does not contradict the
special unit-circle statement.

### 2.1 Regular polygons with the actual phases retained

For P_n(theta)={theta+i/n: 0<=i<n}, a unit chord exists on C_R iff
n alpha_R is an integer. If alpha_R=p/q is reduced and q divides n,
the induced unit graph is n/q disjoint C_q's (with matching convention at
q=2). Otherwise it is empty. In particular, a radius-one regular n-gon
has any unit chords iff 6 divides n.

For two labelled polygon lists P_n(theta), P_m(phi) on the same circle,
put L=lcm(n,m), g=gcd(n,m), and Delta=phi-theta. Their labelled vertex
pairs collide iff L Delta is an integer. The number of colliding pairs
is then g. For each sign separately, there are exactly g unit-contact
pairs iff L(Delta-alpha_R) or L(Delta+alpha_R), respectively, is an
integer, and none otherwise. For alpha_R<1/2 the two signed pair sets
are disjoint. At alpha_R=1/2 count them once.

Proof: (i,j) -> j/m-i/n maps Z/n x Z/m onto the L-torsion subgroup,
with every fibre of size nm/L=g. Colliding labels must be identified
before constructing the physical graph. The cross-list pair counts are
not undirected edge counts after that identification.

## 3. The positive result: preserve the exterior ports and extend

**Lemma (constructive 2-choosability).** A finite union of paths and even
cycles admits a proper colouring from arbitrary vertex lists of size at
least two. Paths are greedy. On an even cycle, if every list is identical,
alternate two colours. Otherwise choose neighbouring u,v and a colour in
L(u) minus L(v), reversing the order if necessary. Give that colour to u,
then greedily colour the path from u around the cycle to v. The final
vertex v cannot conflict with u because u's colour is absent from L(v).

Let G be a finite faithfully placed unit-distance graph with vertex set
the disjoint union {o}, N, and X, where o is the centre, N is its actual
unit-circle neighbour set, and X=V(G) minus N[o] is the entire off-circle
exterior. Include every actual edge; in particular o
has no neighbours in X. Fix a proper q-colouring beta of G[X], q>=3.
These are fixed physical ports, not a quotient palette. For z in N let

    F_z = { beta(x): x in X and |x-z|=1 }.

**Theorem (circle-port extension).** A prescribed colour c for o extends
beta to the full graph whenever

    |F_z union {c}| <= q-2             for every z in N.       (1)

Indeed the exact available list is

    L_c(z) = {0,...,q-1} minus (F_z union {c}),

and section 2 and the lemma colour N. All X vertices, their colours,
their identities, and all unit edges remain unchanged.

In particular, |F_z|<=q-3 for every z guarantees extension for **every**
centre colour c. At q=5 this yields two useful certificates:

1. Every rim vertex has at most two distinct forbidden exterior colours.
2. Every F_z has size at most three, and the intersection of all size-three
   F_z is nonempty. Any c in that intersection works. With no size-three
   sets, every c works.

These conditions are sufficient, not necessary. Actual lists of size one
can still have a colouring. The checker retains them and solves the
finite path/C6 constraints directly rather than declaring failure unsafe.

**Port boundary.** The theorem leaves N free to recolour. If a rim vertex
is itself a prescribed port, replace its list by the corresponding
singleton (or the empty list on a conflict), and use exact list checking.
Do not apply the size-two theorem to those singleton inputs. A fixed
five-coloured star boundary can fail at the centre even though the graph
is bipartite. This is not a counterexample: its prescribed ports are on
the circle, not the free vertices N in this theorem.

### 3.1 Consequences for genuine five-colour saturation

If the off-circle unit graph G[X] is three-colourable, colour it with
three colours, the bipartite circle with two fresh colours, and o with an
exterior colour. Thus G is five-colourable. A genuine six-chromatic centre
construction needs an off-circle core of chromatic number at least four.

More pointedly, if G[X] is five-colourable and every z in N has at most
two exterior neighbours, (1) constructs a five-colouring of G. Therefore,
in every vertex-critical six-chromatic unit-distance graph, for **each**
root o there is a neighbour z having at least three neighbours outside
the closed neighbourhood N[o]. To prove this, G-o is five-colourable,
so its restriction to X is an admissible beta; a failure of the conclusion
would invoke the preceding extension. This is a necessary condition, not
a sufficient construction or a new general order bound.

The critical-graph consequence also holds for noninduced unit-distance
graphs: use their graph neighbourhood and define F_z using only listed
graph edges zx in E(G), not every geometric unit contact. The same proof
applies. The neighbourhood graph is still a subgraph of disjoint
C6's. In that convention, points outside N[o] need not be geometrically
off-circle; the full-induced convention above is what makes them so.

The chosen-centre threshold is sharp: take two adjacent rim vertices on
the unit circle and attach three independent off-circle leaves to each,
coloured 0,1,2. Prescribe centre colour 3 with five colours available.
Both rim lists are then {4}, so no extension exists. Prescribing centre
colour 0 instead repairs it. The checker gives an exact induced Euclidean
realization in Q(sqrt(3)); the underlying graph is only three-chromatic.

## 4. Runner clocks: what the unit graph really observes

For integer runner speeds v_i and a reduced rational clock t=a/q on the
radius-one circle, positions have angular phase a v_i/q. A unit edge
between runners exists iff

    6 divides q and a(v_i-v_j) = +/-q/6 mod q.             (2)

Collisions are v_i=v_j mod q. At the tight AP clock t=1/14 for speeds
0,1,...,13, the radius-one graph is empty, although the runner separation
is exactly 1/14. At t=1/12 the full AP gives two C6 components. Equations
(2), including the actual clock and phase, are the exact bridge; a
tournament score or polygon count alone cannot recover them.

For any N-runner configuration, loneliness at threshold 1/N is equivalent
to all observer-runner Euclidean chord distances being at least
2 sin(pi/N) on the unit circle. Rescaling to radius
1/(2 sin(pi/N)) makes this a distance-at-least-one statement. The unit
graph remembers only equality, not which side of the threshold a
nonedge lies on. At the tight AP clock its rescaled graph is merely C_N,
two-colourable when N=14. There is no six-colour forcing consequence.

## 5. Finite resonance versus all radii

The following observation was contributed by the parent research lane.
Let H be an orientation-preserving circle homeomorphism preserving every
unit edge on C_R, with R>1/2 and alpha=alpha_R. A degree-one increasing
lift satisfies

    H(x+alpha)-H(x) in {alpha,1-alpha}.

The left side is continuous in x and lies strictly between zero and one,
so it is constant. Writing H(x)=x+h(x) with h periodic and integrating
over one period shows the constant is alpha. Therefore H commutes with
rotation by alpha.

If alpha=p/q is rational, this is equivalent to commuting with rotation
by 1/q. Nonlinear such H exist: the lift
H(x)=x+epsilon sin(2 pi q x) is increasing when 0<|epsilon|<1/(2 pi q).
If alpha is irrational, continuity and
density force h to be constant, so H is a rotation. Radius R=2 already
has irrational alpha: cos(2 pi alpha)=7/8; were alpha rational, twice
this cosine would be both a rational number and an algebraic integer,
hence an integer, contrary to 7/4.

Thus a nonlinear clock compatible with a finite common period can retain
unit-circle contacts, but cannot preserve all unit contacts even on this
single additional radius. This is an exact obstruction to promoting a
resonant angular construction to a plane-wide unit-distance symmetry.

## 6. Connection contracts and tests

| Source | Target and map | Preserved | Lost / required sidecar | Cheapest hostile |
|---|---|---|---|---|
| Labelled regular polygon lists on C_R | phase differences modulo 1 | exact unit contacts and collisions | radius R, phases, label-to-physical-vertex map | n=14 at R=1 versus its side-unit radius |
| Integer runners at actual t | circle embedding exp(2 pi i v_i t) | all phases and actual chord lengths | inequalities and observer; equality graph alone loses them | AP14 graph empty at radius 1 |
| Fixed off-circle beta and actual N | exact lists L_c(z) | all fixed port colours and all cross-edges | rim prescriptions must become singleton lists | prescribed five-colour star |
| Finite resonant clock homeomorphism | radial projection to C_1 | contacts at +/-1/6 when clock commutes with 1/6 | physical polygon chords and other radii | radius 2 forces a rotation |

There is no intrinsic binary tournament relation in unit-edge equality.
Orienting all nonedges would add data unrelated to colouring. The faithful
carrier here is the graph of angular orbit steps with literal list states.

## 7. Primary literature and scope

- Alon–Bucic–Davies, *Hadwiger–Nelson Problem for Typical Norms*,
  [arXiv:2609.36220v1](https://arxiv.org/html/2609.36220v1), submitted
  2026-09-28; checked 2026-10-21. Its introduction retains Euclidean
  bounds 5 through 7. Theorem 1 gives a 2d upper bound for typical norms,
  not for the Euclidean norm. Theorem 3 is a genuine shifted matrix
  avoidance theorem: for 1<=k<n, full-spark n rows in R^k guarantee
  supremum clearance at least k/(2n), attained for rational rows. It neither proves LRC(14) nor turns
  this circle carrier into a high-chromatic graph.
- De Grey, *The chromatic number of the plane is at least 5*,
  [arXiv:1804.02385v3](https://arxiv.org/abs/1804.02385v3), 2018-05-30,
  supplies the established finite non-four-colourable unit-distance graph.
  The elementary statements above do not improve that bound.

## Reproduction and remaining frontier

Run `python3 check.py` and `python3 -O check.py` from this directory.
There are no third-party packages or repository imports. Both must emit
the identical JSON receipt. Exact rational phase arithmetic certifies the
geometric phase predicates by the analytic chord identity. An independent
backtracker verifies list witnesses; exact Q(sqrt(3)) arithmetic verifies
every distance in the sharpness example. No floating tolerance is used.

`python3 verify.py` checks every frozen local file hash and reruns both
modes against the saved receipts; it also works from another directory.
See `AUDIT.md` for the independent review and its evidence boundaries.

The open forcing target is an actual off-circle context which constrains
every five-colouring so that the centre sees all five colours. The present
theorem excludes broad free-rim low-contact constructions; it does not
construct that context. No counts here are an unbounded graph census.
