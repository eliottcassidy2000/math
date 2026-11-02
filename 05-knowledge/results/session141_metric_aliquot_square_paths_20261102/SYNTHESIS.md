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
map or invariant changes. Curvature and square-path connections below supply their own preservation
proofs rather than borrowing authority from this analogy.


## Triangle inequality is the packing mechanism—not the obstruction to curvature

If every source edge has ambient length at most L, triangle inequality gives
`d_X(f(u),f(v))<=L*d_G(u,v)` and places a radius-r source ball inside an
ambient radius-Lr ball. It gives no reverse inequality. With image separation
delta and multiplicity at most K, disjoint radius-delta/2 disks yield

    (N/K)[cosh(kappa*delta/2)-1]
        <= cosh(kappa*L*r+kappa*delta/2)-1.

Euclidean disk capacity is quadratic; constant negative-curvature capacity
is exponential. Hence a binary tree requires kappa*L>=log2. A positive
relative gap between child intervals repairs the missing lower-distance
bound: the explicit q<1/2 construction is proper, noncrossing, separated and
quasi-isometrically embedded, with every edge of any prescribed length
satisfying the strict bound. The critical equality case remains unresolved.
At q=1/2, distant cross-root vertices remain at bounded ambient distance;
the inherited dyadic placement is radially efficient but not quasi-isometric.
The repaired family need not have its minimum-distance pairs exactly on edges.

For native Collatz, clock matters. The inherited five-step binary supply
implies only kappa*L>=log2/5 under uniform separation and bounded fibres.
Odd macro edges can share even native tails:511 odd states at depth8 expand
to2001 distinct native states, not2811 separate subdivision vertices.
The source union and first-hit clock are the needed sidecars. None of these
statements constructs a full native hyperbolic Collatz embedding.

## Square paths retain joint compatibility, not just pair counts

For a full labelled tournament T, H_G(T) counts only Hamilton paths whose
edges pass the square-sum guard. Its polynomial is not unrestricted H(T).
For a compatible path family S, the union of its prescribed directed arcs
has exactly `n-b+|J|` arcs, where b is the number of maximal common directed
blocks and J is the set of distinct ordered block joins. Opposite block
joins use different original endpoint pairs unless both blocks are
singletons; only in that singleton case are they contradictory directions.
This proves both the compatibility criterion and the joint weight.

The inherited n25 response is

    2^300 +20*2^276(z-1)+2^276(z-1)^2.

There are20 directed paths (10 up to reversal), four compatible pairs and
no triples. At n26 the complete bank has24 paths,12 compatible pairs and
two triples. Eight pairs prescribe27 arcs, two prescribe29 and two prescribe30;
each triple prescribes29. This yields

    2^296(536870561+320z+29z^2+2z^3).

The n25 matching-only simplification fails because n26 has compatible triples.
The common-block theorem, rather than a guessed analogy, supplies the exact
extension. At n24 the path bank is empty and the polynomial is2^276.
Finite maxima at n25..30 are2,3,5,5,3,3: no monotonicity is claimed.
A source-preserving five-block directed carrier explains the n26 triple,
but arbitrary abstract blow-ups are not automatically square-sum labelled.

## Binary content explains the normalized powers of two

Integer translation z=w+1 is unimodular on polynomial coefficients, so the
coefficient gcd equals the gcd of the binomial moments. Reversal pairs every
nonempty compatible subset with a distinct one of the same union size.
For a nonempty reversal-closed bank at n>=2 this guarantees divisibility by
`2^(N-Umax+1)`, where Umax is the largest union-arc count. Joint multiplicities
can force an additional factor. The exact contents are2^276 at n25 and
2^296 at n26: reversal alone misses one factor at25 but is sharp at26.
Here Umax grows26->30 and N grows300->325. The exponent jump20 is therefore
`25-4-1`:25 additional tournament pairs, four extra maximally prescribed
arcs, and the loss of the one-bit surplus over the reversal bound. This is
a derived identity of the proved response data, not a connection to aliquot20.

## Four exact appearances of20, and one false proposed role

| Role | Defining object and equation | What20 means |
|---|---|---|
| Arithmetic source | sigma(20)=42; s(20)-20=2 | A labelled excess-two integer, not its divisor sum |
| Discrete curvature | Closed spherical triangulation, all degrees5; sum(6-deg)=12 | Euler/incidence force V=12,E=30,F=20 |
| Square-path bank | G25 on labels1..25, complete directed Hamilton bank |20 orientations of10 reversal classes |
| Binary-content exponent jump | Exact contents2^276 and2^296 |20=25-4-1 from changed free-orientation and reversal data |
| Radius21 count difference |704 source states versus694 honeycomb sites | The difference is10, not20 |

These are exact but different statements. The sphere role uses an actual
curvature budget; the path and binary-content roles use joint overlap and
reversal; the aliquot role uses source factorization. No map identifying their counted objects has been supplied.
Matching numerals are recorded as prompts, not promoted into a theorem.


## A graph-theoretic bridge, not merely a shared lesson

The [common-block / successor-loss corollary](CONNECTIONS.md) applies the
functional defect formula to the actual compatible directed path union.
For b common blocks and |J| distinct ordered joins, it proves
`beta=|J|-b+1` and `epsilon=max(0,|J|-b)` when reorientation is allowed.
The n25 pair has cycle rank2 and deletion loss1; the n26 triple has cycle
rank4 and loss3. The contracted quotient must retain parallel edges from
opposite joins between nontrivial blocks. This is a precise cross-lane
operation preserving a named invariant; it does not assert that square-sum
labels make an associated Euclidean unit-distance realization.


Keeping the path arrows adds a sink sidecar: the loss is m-n+z, not merely
the undirected defect. Displayed families starting18 have costs1/3;
reversing them makes18 a shared terminal and raises costs to2/4. Reversal
therefore preserves response weights and undirected defect but can change
fixed-successor cost. This is the exact lost-coordinate mechanism, not a
numerical1–2–4 pattern.
