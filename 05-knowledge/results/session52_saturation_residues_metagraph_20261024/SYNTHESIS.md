# Session 52 synthesis: saturation, incidence, and what survives a quotient

24 October 2026. This is a map of the certified results and remaining
obligations, not a claim that the requested open problems have been solved.

## 1. What the repeated word actually records

The word `3, 2, 3, 5, 1, 1, 2` is reproduced exactly, but only after specifying the
observable on each side. The Fortunate observable is equality of the
primorial offset F_n with the next indexed prime p_(n+1), for n=1..17.
The graph observable is equality of the minimum vertex-disjoint path-cover
number pc(G_N) with the number of connected components c(G_N), for N=10..26.
Both yield `11100111000001011`. It is not a comparison of raw numerical gaps,
and the first three graph values do not assert whole-graph Hamiltonian-path existence.

Both graph-window end runs are truncated. At the very next pair, n=18 and
N=27, the arithmetic equality fails and the graph equality holds. Thus the
right common language is **lower-bound saturation**, not an established
map from prime offsets to graph paths. This is a useful interpretation of
the observation, but an interpretation selected after seeing the word is
not statistical evidence of a hidden universal correspondence.

The genuine general theorem is on graph growth: adding a vertex changes pc
by -1, 0, or 1. A decrease occurs exactly when two distinct paths in some
minimum old cover have endpoints adjacent to the new vertex. If that vertex
meets r old components, the excess pc-c changes by change(pc)-1+r.
Component merging and endpoint repair are separate mechanisms. The exact
endpoint witnesses explain transitions that the Boolean word alone hides.

See [the theorem](fortunate-square-pattern/THEOREM.md) and
[independent audit](audit-fortunate-square/AUDIT.md).

## 2. A count-generating carrier that retains the actual path

For even n=2 p, almost-regular tournaments have scores p-1 and p. Native C
reverses all arcs outside the fixed Hamiltonian path. Its score change at
the first and last vertex differs from its change at an interior vertex.
Consequently an almost-regular source stays in that score layer precisely
when its path starts at score p and ends at score p-1.

These internal ports correspond, C-equivariantly, to regular tournaments on
n+1 vertices containing one specified Hamilton cycle, or to Eulerian
orientations of K_(n+1) with that cycle deleted. This is a structural reason
for a coefficient formula, not just a fit to a numerical sequence. The
generic score partition product was already THM-554/555; this session
extracts endpoint and reflection refinements, the cycle carrier, and an
order-twelve extension with independent coefficient algorithms.

At n=12 the exact internal C-line count is 1,606,912,987,280. Its blue part
is 3,401,824 and its black part is 1,606,909,585,456. These do not count
isomorphism classes, distinguish loops, or determine which two classes
each line joins. The native lines have not been divided again by reflection.

See [balanced counts](balanced-metagraph/THEOREM.md) and
[audit](audit-balanced-metagraph/AUDIT.md).

## 3. Why A000568 and full wiring are different obligations

The earlier repository already contains a fast ordinary-tournament count
algorithm, including a reported run through order 200. Recomputing Burnside
or its signed cycle index is a useful independent baseline, not a new
prediction or a new fast enumeration method.

The signed pair-cycle polynomial P_n(t) has two exact specializations:
P_n(1) counts ordinary classes, and P_n(-1) counts self-converse classes.
Their half-sum counts converse-merged nodes. Its Fourier coefficients give
the spectrum of the **single-arc-flip weighted class quotient**, not the
native simultaneous non-path complement graph. At n=3 those operators
already have spectra {3,-1} and {1,-1}, respectively.

At n=4, two presentations of one source class have native complements in
different converse-merged target classes. Thus the bare source class cannot
determine a native target. No improved scalar node-count formula removes
this missing incidence information. A proposed class-level wiring rule must
either retain the Hamiltonian-path presentation or correctly aggregate its
entire fibre, with multiplicities and endpoint data.

The correct self-converse cycle filter was already established in THM-283:
at most one fixed vertex, with every other permutation cycle length equal
to 2 modulo 4. This session supplies independent checks and flags a stale
historical draft that omitted those restrictions. Separate narrow precision
repairs distinguish character traces from fixed-subspace dimensions and
center the cyclic-triangle eigenfunction by binom(n, 3)/4. A graded trace
does not supply a differential or a homology theory.

The known order-eleven count cannot be substituted into a native identity
and called a native prediction. In particular, 9,941,438,418 obtained from
11*A000568(11)+152690 is a backsolved consistency target. Independently
evaluating that native upper trace remains a separate obligation.

See [counting lenses](counting-lenses/README.md),
[audit](audit-counting-lenses/AUDIT.md), and
[precision notes](counting-precision/PRECISION_NOTES.md).

## 4. A complete finite blue graph, not just a line total

Completeness uses the frozen inherited almost-regular atlas; the blue
reflection cube and every recorded incidence are independently checked.
The order-ten calculation retains every blue port and its literal source
path and native-C partner. On all 229 self-converse almost-regular classes,
the induced blue multigraph has 11,452 edges, including 24 loops. It is
connected, has no bridges or articulation vertices, and has cycle rank
11,224. There are 18,795 outgoing blue lines to 786 distinct boundary classes.
The boundary target score profile is (3, 4, 4, 4, 4, 5, 5, 5, 5, 6).

Ignoring loops and parallel multiplicity for graph distance, the unordered
vertex-pair distance histogram is 9,219 at distance 1, 16,886 at distance 2,
and exactly one pair at distance 3. That exceptional pair consists of atlas
indices 63 and 86, each of simple degree 17 and automorphism order 5. Their
total blue presentation counts are 37 and 35, with 20 and 18 boundary edges.
These are exact finite distinctions, not an identification with a prime
family or a claim that equal degrees explain the obstruction.

There is an additional structural feature: the inherited 20 positive-error
balanced nodes form a connected dominating set. Their induced graph has 90
edges; 1,857 edges join this core to the other 209 balanced nodes; 9,505 lie
entirely outside the core. Every outside node is one step from the core.
The four isolated positive-error nodes in the earlier larger study are
nonbalanced and are not among these 229 nodes. Thus this result must not be
described as connecting those four components.

A spanning tree with actual C-word witnesses and an independent ear
certificate make connectivity and two-connectivity checkable without
trusting a visualization. This is a finite structural result. Its general
extension would require a proof that an analogous core dominates every
admissible class, with actual compatible ports, not merely enough total
edges. The black graph and full all-size wiring remain outside this claim.

See [blue wiring](balanced-blue-wiring/README.md) and
[independent audit](audit-balanced-blue-wiring/AUDIT.md).

### An all-size existence route, with its source boundary retained

The next probe found a relevant existing result, rather than evidence that
the finite graph must extrapolate. Farrugia's survey, section 2.30, printed
page 59, attributes to Salvi-Zagaglia (1980) a symmetric-cycle theorem for a
regular self-converse tournament with a specified anti-involution. For
order at least 7 its stated cycle-length range includes the Hamilton length
through every invariant arc. The original 1980 paper was not accessible in
this session; this is an explicitly secondary-source lead, not a new proof
or a claim to have inspected the original theorem.

The translation itself is proved: extend the anti-involution over the
canonical regular completion by fixing the new vertex z, and delete z from
a symmetric Hamilton cycle. The remaining path is symmetric, starts high,
and ends low. Conversely every such path closes through z. Small completion
orders 3 and 5 are checked directly, with order 7 as another finite control.
Under the quoted cycle theorem, every specified anti-involution therefore
admits an internal blue path, with each invariant paired arc attainable as
its middle arc. The packet conservatively retains CONDITIONAL status for
this all-size consequence pending primary-source verification.

Existence of one or several internal ports at each class does not imply
connectivity, two-connectivity, or compatible global wiring. Those remain
distinct obligations even if the cited existence theorem is accepted.
See the [existence probe](blue-existence-probe/PROBE.md) and its
[translation audit](audit-blue-existence/AUDIT.md).

## 5. The real fifteen-object bridge, and its limit

Fifteen perfect matchings on six labels are the fifteen maximal isotropic
planes of a binary symplectic four-space. This carrier was already present
in prior work. The current certificate makes the six-pentad construction
explicit: the graph of matchings sharing an edge is isomorphic to KG(6, 2).
Transporting the action requires the outer automorphism of S6. An ordinary
transposition fixes seven two-element subsets but only three matchings,
so equal cardinalities cannot establish ordinary equivariance.

For the Bolza Jacobian in characteristic outside {2,3,5}, over an
algebraic closure or a sufficiently split field, the cited fifteen Richelot
kernels split into
automorphism orbits of sizes 6, 4, 4, 1. The size-six orbit has an elliptic-square
codomain with j=8000; the singleton is the loop. Four weighted edge classes
need not be four distinct destination vertices. Characteristic exceptions
and the ground-field qualification remain essential.

This is a genuine route from a previously studied curve to fifteen kernel
choices and to the j=8000 CM geometry appearing in the m=8 MRB residue.
It is not a bijection with Ogg's fifteen primes and does not identify the
MRB constant or its small defect. At odd good primes, the elliptic curve's supersingular
reduction is governed by the character (-8/p), not by membership in Ogg's
set and not by reduction modulo 15. The primes 17 and 37 separate those
predicates explicitly.

The exact mod 15 unit-filtered alternating MRB series is
M(1)-M(3)-M(5)+M(15), not 8/15 of the original constant. Its mask has period 30
because parity remains part of the state. Removing parity changes a
convergent alternating problem into a divergent positive one.

See [fifteen kernels](fifteen-kernels/README.md),
[audit](audit-fifteen-kernels/AUDIT.md), and
[MRB/CM analysis](mrb-supersingular/RESEARCH.md).

## 6. What the unification is, and is not

Across these examples, the successful common practice is to retain the
data a quotient discards: path endpoints, the chosen observer path,
orientation signs around permutation orbits, actual stabilizers, and
kernel incidences. This is a methodological synthesis, not a newly proved
universal mathematical equivalence.

The marked regularization makes the warning concrete. A known Type-I
construction lifts each source arc to a disjoint directed four-cycle and
each source Hamiltonian path to a distinguished alternating cycle. Keeping
the center preserves the source class. Forgetting it can identify different
classes, including a converse-merged collision already at source order 5.
The source native operation also differs from the full target operation on
5(n-1) edges. A carrier embedding is not automatically an operation-preserving
or guard-preserving embedding.

See [regularization](marked-regularization/THEOREM.md) and
[audit](audit-marked-regularization/AUDIT.md).

The next useful all-size result would explain when these retained local
incidences admit compatible global assignments, rather than merely matching
their totals. Neither the finite saturation word, the Eulerian coefficient
formula, nor the fifteen-object geometry by itself supplies that theorem.
