# What actually connects—and what only looks similar

## Geometry needs three different questions

A functional map supplies at most one outgoing transition per vertex.
For a finite simple graph with n vertices, m edges and t tree components,
compressing to one such successor system loses exactly m-n+t edges.
Subdividing edges does not repair that deficit. For k layers, the exact test
is |E(W)|<=k|W| on every vertex set W, not just a global edge count.
Applied to all plane unit-distance graphs below the first u(n)>kn crossover,
this becomes an exact threshold theorem. The layers are arbitrary partial
successor systems; they are not automatically Collatz maps.

At the same time, generic relative rotations realize Cartesian products
faithfully at unit distance. The K3 powers have arbitrarily large m/n=d
while needing only three colors. Dense unit-distance configurations therefore
supply no automatic Hadwiger–Nelson lower bound. A chosen honeycomb routing
is a subdivision with branch/route/clock sidecars, not a dense extremizer;
its induced honeycomb completion may add unused edges, still within a
bipartite degree-three graph.

## Twenty needs its arithmetic map and its source

Write s(n)=sigma(n)-n, where proper divisors include1. Then
s(20)=1+2+4+5+10=22. Removing1 gives21, not the standard aliquot sum.
The complete total-sum fibre sigma^-1(42)={20,26,41} connects sources with
very different aliquot sums22,16,1. At fixed sigma, s(n)=sigma-n is an
anti-diagonal, so forgetting the source erases the excess and forward orbit.
A sum coincidence is a fibre relation, not a conjugacy of dynamics.

The exact20 tail is20->22->14->10->8->7->1->0. Here7=1+2+4=s(8)
uses the geometric-series divisor structure. Separately, the marked Collatz
cycle1->4->2->1 is multiplication by4 modulo7, inverse to doubling on the
Paley7 triangle. In the specified F8 chart, residue doubling becomes
Frobenius and the marked triple has rank2, not a normal basis. This is a
small marked correspondence; it does not identify aliquot dynamics with
Collatz or prove a global descent modulo7.

## Connection contracts

| Source -> target | Preserved predicate | Lost information / required sidecar | Decisive hostile |
|---|---|---|---|
| Graph -> k successor layers | Every edge, iff all subset capacities hold | Layer and endpoint-slot assignment; no arithmetic map identity | A dense induced subset despite a harmless global average |
| Graph -> edge subdivision | Pseudoforest deletion defect and branch degrees | Original vertices, route lengths and transition clocks | Degree>3 cannot be routed as a honeycomb subdivision |
| Cartesian point sums -> unit-distance product | Every edge/nonedge after forbidden angles are avoided | Relative rotation and collision exclusions | K3 powers: many edges, always three colors |
| Divisor list -> sigma or s value | The named sum | Source, prime exponents, multiplicity and omitted1 convention |20 and26 have equal sigma, unequal s |
| Marked residue exponent -> F8 power | Doubling/Frobenius action | Marking and linear rank | The chosen orbit is rank2; another generator has rank3 |

The common lesson is operation-specific sufficiency, not a cross-field
isomorphism. A number such as20,24,25 or26 changes meaning when the source,
map or invariant changes. Curvature and square-path connections require
their own preservation proofs and are staged separately.
