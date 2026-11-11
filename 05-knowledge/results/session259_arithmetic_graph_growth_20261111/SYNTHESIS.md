# Two growing arithmetic graphs: four clocks and a common diagonal mechanism

**Session259, 2026-11-11. COMPLETE mathematical synthesis.**
The central advance is a usable distinction between **numeric birth**, **parent
arrival**, **recursive generation**, and **topological lifetime**. The same
integer can play very different roles on those clocks. Doubles and squares
are not merely visual analogies: they arise as fixed points of the same swap,
as the two boundaries of a sum–product graph, and as exact birth/death scalings
in multiplicative topology.

The [definition ledger](README.md) preserves the original paired-parent
conventions. Full proofs are in the four linked lane notes. New proofs received
a [focused independent review](review/REVIEW.md); existing repository claims
were inputs, not the subject of a new audit. Finite receipts support only their
explicit ranges. No literature-priority claim or numbered canon promotion is
made here.

## 1. What actually changes when integer n is added?

Use distinct positive summands and distinct nonunit factors. At birth n:

    summand pairs       floor((n-1)/2),
    multiplicand pairs  (tau(n)-2-[n is a square])/2, n>=2.

Weak variants restore equal parents. They add one pair at n=2a and n=a²,
respectively. Both statements are inherited from
[THM-2422](../../../01-canon/theorems/THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md):
the involution (a,b)<->(b,a) pairs every nonfixed witness, leaving precisely
those diagonals. The multiplicative analogue of even/odd is **square/nonsquare**,
not another period-two count.

A summand birth is a matching on the older nodes, with the half-node unmatched
at even n. A multiplicand birth is a matching on proper divisors, with the
square root unmatched at squares. The parent pairs, not just their two arcs,
are the faithful event records.

The [node atlas](BIRTH_ATLAS.csv) records every n from 1 through 128: strict and
weak pair counts, diagonal indicators, weak and strict topological births and
deaths. For a quick reading:

| n | Strict sum pairs | Strict proper factor pairs | What the extra topology sees |
|---|---:|---:|---|
| 4 | 1 | 0 | Strict square atom begins an H0 interval [4,8) |
| 6 | 2 | 1 | Weak component born at3 merges now |
| 9 | 4 | 0 | Strict square atom begins [9,18) |
| 12 | 5 | 2 | Strict square interaction births a loop [12,24) |
| 15 | 7 | 1 | Weak squarefree interaction births a loop [15,30) |
| 16 | 7 | 1 | Square diagonal event, but no reduced homology event |
| 27 | 13 | 1 | Strict cube interaction births a loop [27,54) |
| 30 | 14 | 3 | The weak loop born at15 dies |
| 60 | 29 | 5 | Strict interaction births a two-dimensional class [60,120) |

“H0 interval” means an extra component relative to the essential component.
Filled triangles matter: a loop in the ordinary undirected graph need not be
a homology class in its order complex.

### The old arc-shadow weld, now kept separate

Let S be the strict additive **arc shadow**, M the strict nonunit multiplicative
arc shadow, and D={u->2u:u>=1}. Then

    M = (M intersect S) disjoint union {u->2u:u>=3}.         (W)

Proof: S contains every increasing arc x->z except z=2x. Such a missing arc
belongs to M exactly when its other factor2 is nonunit and distinct from x,
which excludes x=1,2 and permits all x>=3. Thus the only multiplicative arcs
invisible to the additive shadow are the surviving doubling arcs. The square
removal in M and the double removal in S meet exceptionally at2->4.

This is an elementary arc-level connection, not a statement about paired
closures, the transport graph below, or equality of operation fibres. Retaining
that distinction is what lets the old dichotomy join the new growth theory.

## 2. The product is the last missing summand, not just another operation

[Summand proof, Sections 2–4](summand/RESEARCH.md).
For arbitrary seeds 0<a<b, the exact distinct-parent closure is

    Cl({a,b}) = {a,b} union (a+b + <a,b>),

where <a,b>={ia+jb:i,j>=0}. The first mixed event is a+b; every nonseed
construction contains both seed types, and all positive coefficients can be
built by repeatedly adding either smaller seed.

For coprime 2<=a<b this yields

    last missing integer = ab,
    conductor = ab+1,
    number of missing positive integers = (ab+a+b-5)/2.

For seeds3 and5 the holes are1,2,4,6,7,9,10,12,15. The first mixed child8 and
last hole15 are the sum and product of the **same pair**. Their difference7
is the final gap of the ordinary semigroup where equal-parent addition is
allowed. Thus diagonal deletion shifts the classical conductor boundary by
a+b. For noncoprime seeds, normalize by gcd and work in its lattice; off-lattice
nodes never appear. Seed1 is a separate explicitly solved boundary.

Finite startup holes leave a permanent witness-storage scar. If the hole set
has size h, sum q and maximum F, then for every N>=2F,

    total legal strict sum events through N
      = floor((N-1)^2/4) - hN + q + h(h+1)/2.

The full record is quadratic, while one selected pair per nonseed output gives
a linear-size existence DAG. Compression preserves constructibility and loses
alternative witnesses, fault tolerance and generation information.

## 3. Recursive growth can hide parity for exponentially long frontiers

[Summand proof, Sections 5–6](summand/RESEARCH.md).
Starting with {2,b}, b odd and at least3, let every generation add all legal
sums of currently available distinct values. Let M_t be the maximum after t
generations, and T=(b+5)/2. Exactly for 3<=t<=T,

    M_t=(3b+4)2^(t-3)+2,
    opposite-parity deficit = b-2(t-3).

The three largest values of the leading parity are M_t,M_t-2,M_t-4. At T the
opposite parity finally reaches M_T-1. Thereafter three consecutive top values
persist and M_(t+1)=2M_t-1.

For {2,101}, the static closure's last hole is only202, yet the moving upper
frontier does not enter its consecutive regime until generation53, at

    M_53=307*2^50+2.

There is no contradiction: final membership and time of construction are
different invariants. This is an exact, arbitrarily long parity transient,
not a finite-data conjecture. It does not assert that every value between the
static conductor and moving frontier has appeared by generation T.

The reusable coordinate is the deficit set D={M-x:x available}. With
 g=min(D\{0}), its exact update is

    D'={d+e-g:d<e in D} union {M-g+d:d in D}.

Exponential absolute growth becomes a subtractive small-gap evolution. The
absolute maximum M must remain as a sidecar for the retained old vertices.

## 4. Multiplicand topology: squarefree births, doubling deaths, square repairs

[Multiplicand proof, Sections 2–4 and 8–9](multiplicand/RESEARCH.md).
Take the unit-deleted divisibility order on2,...,N, and fill every comparable
chain as a simplex. Its radical map n->rad(n) retracts the poset naturally to its squarefree
subposet. The resulting order complex is the barycentric subdivision of the
complex of prime sets whose product is at most N. This removes exponent
multiplicity while preserving the entire filtered homology.

The complete weak reduced barcode is

    [m,2m) in degree omega(m)-1,
    for every odd squarefree m>1.

For m=15 the boundary on prime vertices2,3,5 is completed when the edge3–5
arrives at15, and filled when the triangle arrives at30. Every bar has exactly
one doubling of life, or length log2 on a logarithmic time axis. Restricting to
an allowed prime set with least member q replaces the death clock2 by q.

This separates two dichotomies that the fibre count alone does not distinguish:

- **Square/nonsquare:** whether the parent-swap involution has a fixed point.
- **Squarefree/nonsquarefree:** whether a weak divisor birth changes topology.

For instance12 is nonsquare yet topologically silent in the weak order complex.
Deleting equal-factor events changes the answer in a precise way. Transitive
closure repairs every removed a->a² relation with composite a, but cannot
repair p->p² for prime p. These missing relations make prime squares new atoms.
The **complete strict-reachability barcode** consists of the weak bars plus

    [p²,2p²) in degree0,
    [p²m,2p²m) in degree omega(m)
       for every prime p and odd squarefree m>1.

The extra bars are suspended, prime-square-rescaled copies of the weak ones.
The repeated prime is permitted in m, explaining [27,54). Extra isolated
components alone would miss the first new loop [12,24) and the higher class
[60,120).

At most one reduced barcode birth or death happens at each integer. Birth
values have 2-adic valuation0 or2; their doubled death values have valuation1
or3. Consequently every multiple of16 is homologically silent for this strict
filtration, even though it may add many ordinary factor edges. This is a new
precise odd/even versus square/squarefree switchboard.

The signed count also gives the exact identities

    reduced Euler(weak_N) = -Mertens(N),
    reduced Euler(strict_N)
      = -Mertens(N) + sum_(p²<=N) Mertens(floor(N/p²)).

They translate cancellation into topology but do not establish a new bound for
the Mertens function or a statement about zeta zeros. Raw strict-arc clique
persistence remains distinct from the solved reachability-order persistence.

## 5. A causal graph puts doubling and squaring on opposite boundaries

[Transport proof](transport/RESEARCH.md).
Keep each strict nonunit pair as an edge

    a+b -> ab,             2<=a<b.

This graph is lossless on unordered parent pairs: s=a+b and p=ab recover them
as roots of X²-sX+p. Its positive-square discriminant is (b-a)². Every edge
points forward in numerical time; including unit parents would introduce
backward edges, and including2,2 would introduce the exceptional loop4->4.

For a source n, products appear from the first time2n-4 up to a nearly square
last time n²/4. Exact strict maxima are (n²-1)/4 for odd n and n²/4-1 for even n;
the weak even variant restores the final square. Thus the whole summand fibre
has a **linear doubling boundary** and a **quadratic squaring boundary**.

Shift to u=n-4, centering at the unique binary product-sum fixed point4.
The smaller-parent-a lane becomes

    u -> a*u-(a-2)².

The first three lanes are2u,3u-1 and4u-4. Identifying doubling rays exposes a
minus-Collatz branch, but not a deterministic description of the entire
quotient and not plus-Collatz convergence. The factor-four lane proves a more
useful global fact immediately:

    u -> 4(u-1) <- 2(u-1) <- u-1,        u>=5.

Hence the underlying undirected transport graph is connected, with explicit
O(log u)-length paths to u=1. Directed reachability is a separate unresolved
problem: a prime source such as7 cannot be reached from5.

The exact number of summand pairs whose products have not arrived by N is

    floor((N-1)^2/4)-(N-2)
      -(sum_(k<=N)tau(k)-2N+2-floor(sqrt(N)))/2.

It grows as N²/4+O(N log N): dense additive promises accumulate much faster
than sparse multiplicative completions. This is a storage/scheduling statement,
not a prime-existence theorem.

## 6. Old dichotomies become higher-arity rules

[Dichotomy proof](dichotomies/RESEARCH.md).
The binary swap mechanism extends along prime cyclic rotation. For prime ell,
ordered ell-parent counts modulo ell retain only constant tuples:

    addition detects n=ell*a,
    multiplication detects n=a^ell.

Doubles/squares are ell=2. This is the additive/multiplicative form of the
same Frobenius rule. Composite arity has smaller nonconstant orbits and does
not satisfy this naive extension. Collision strata must be retained.

For the lossless binary chart j=(s,p), scaling parents by c is
D_c(s,p)=(cs,c²p), while squaring them is Q(s,p)=(s²-2p,p²). They satisfy

    Q D_c = D_(c²) Q.

This is an exact noncommuting growth relation, not a visual analogy. Higher
powers follow the Newton recurrence. The chart ceases to be complete for
three parents: (1,6,6) and (2,2,9) share sum13 and product36 but have different
parent-birth times6 and9. Reconstructing a higher-arity parent tuple needs the full elementary-symmetric
coordinates or the tuple itself; a collision partition alone does not suffice.

The earlier blue/black and self-converse involution stories share the
fixed/free-orbit mechanism. Transferring a growth theorem additionally needs
a map that commutes with the involution and retains birth labels. Equal parity
counts alone do not create a bijection of the underlying graphs. Likewise the
log map turns products into sums on log(N), not on an evenly spaced integer
lattice.

## 7. What is genuinely finished, and what deserves the next session?

Finished here:

- Exact two-seed strict summand closure, conductor and witness economics.
- An all-order family with exactly delayed parity mixing.
- Full weak and strict-reachability multiplicand persistence, with the explicit
  birth/death basis and all-order square correction.
- A causal sum–product graph, its complete local birth clocks, weak connectivity,
  and exact backlog count.
- Prime-arity fixed-point and lossless binary transformation rules.

Most useful next questions:

1. Does every coprime pair of seeds eventually acquire three consecutive top
   values in synchronous distinct-sum growth? Find its exact stabilization law,
   then separate that from filling the whole moving interval.
2. What is the birth schedule and eventual size of the large weak component
   in each finite causal transport graph? Classify directed descendants of5
   without mistaking weak connectivity for constructibility.
3. What survives if one uses the raw strict-edge clique complex, rather than
   repairing missing relations by transitive closure? The present barcode
   explicitly belongs to the latter object.
4. Can marked factor multiplicity be added to the radical complex without
   destroying its simple death clock? The present quotient deliberately forgets
   exponents except for the explicit strict square correction.
5. For thin parent sets such as twin-center ranks, which parts of the matching
   and diagonal laws survive, and which need genuinely new correlation input?
   The old all-target twin-center ancestry question remains open.

The recent Chan–Pak work is relevant as a reminder that arithmetic on output
values does not automatically lift to a lawful operation on its combinatorial
objects. This session supplies faithful arithmetic carriers; it does not
supply that missing tournament lift or solve Chan–Pak completeness.

## Reproduction and scope

From repository root, run each lane's script:

- summand/experiment.py: 120 seed pairs through target500, deletion checks and
  seven synchronous families through generation8.
- multiplicand/experiment.py: weak barcode through300, graph laws through120.
- multiplicand/strict_experiment.py: all178 strict bars born through300,
  including56 square corrections, plus Euler/event-word controls.
- dichotomies/scout.py: chart, power, arity and insertion controls.
- transport/check.py: 61,752 pair events with sum<=500 and2,000 explicit paths.
- birth_atlas.py: atlas through128 and unique-event/silence controls through20,000.

Paths are relative to this packet. Normal and optimized Python are compared
for the final scripts. Finite checks are modest on purpose: they test the new
constructions, while all-order conclusions rest on the stated proofs. See
[review](review/REVIEW.md), per-lane receipts and the final SHA256 manifest.
