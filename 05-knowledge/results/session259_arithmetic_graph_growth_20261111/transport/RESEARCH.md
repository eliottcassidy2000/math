# Causal transport from summand births to multiplicand births

**Status: PROVED elementary statements below; scoped finite controls separate.**
This is a new carrier built from the inherited parent fibres, not a replacement
for their definitions. Claims of novelty mean new to this session's synthesis,
not priority claims about the literature.

## 1. A single graph retaining both operations

Let T have nodes n>=5 and an edge

    s=a+b -> p=ab,             2 <= a < b.

The same pair is first seen by summation at s and later by multiplication at p.
Since

    p-s=(a-1)(b-1)-1 >= 1,

this is a locally finite DAG. Distinct unordered parent pairs give distinct
edges: from (s,p), the parents are the two integer roots of X²-sX+p. Equivalently,
Delta=s²-4p must be a positive square, sqrt(Delta) must have s's parity, and
(s-sqrt(Delta))/2>=2. Thus T loses only parent order, which both operations
already forget. The discriminant restores the gap b-a.

Including equal parents adds edges 2a->a², a>=3; enlarging the node set
to include4 also adds the exceptional loop4->4.
Including unit parents introduces backward edges n+1->n and destroys the
causal DAG. Unit deletion is therefore structural, not cosmetic.

## 2. Exact one-node growth and completion clocks

At birth n, every incoming transport edge is supplied by a strict proper
factor pair of n. No old source acquires an edge to an old target. Thus

    indeg_T(n) = (tau(n)-2-[n is square])/2,
    outdeg_T(n) = floor((n-1)/2)-1,                         n>=5.

The indegree-zero sources are primes and prime squares. Every source n>=5
has outgoing edges; there are no sinks. The out-neighbours of n are precisely

    a(n-a),             2 <= a < n/2.

They are strictly increasing with a. Consequently the first and last future
completion times of a summand birth are

    L(n)=2n-4,
    U(n)=floor((n²-1)/4)       if n is odd,
         n²/4-1              if n is even.

For even n, the weak variant adds exactly the square n²/4 at the final time.
The uniform linear first-completion boundary and the near-square last-completion
boundary are an exact doubles/squares dichotomy across the *whole* fibre,
not just a statement about its diagonal. For example birth10 promises products
16,21,24; weak birth10 adds25. Birth11 promises18,24,28,30.

**Proof.** These are the inherited factor count and the strictly concave
quadratic a(n-a); its forward difference is n-2a-1>0 in the indicated range.

For N>=2 let B(N) count pairs already seen by summation but not yet by
multiplication: a+b<=N<ab, 2<=a<b. Let D(N)=sum_{k<=N}tau(k). Then

    B(N)=floor((N-1)²/4)-(N-2)
         -(D(N)-2N+2-floor(sqrt(N)))/2.                    (1)

The first term minus N-2 counts all nonunit strict summand pairs; the last
fraction counts all proper strict product pairs. Every product pair counted
has its sum already born. This proves (1) without an asymptotic approximation.
Elementary D(N)<=N(1+log N) gives

    B(N)=N²/4+O(N log N).

Thus the backlog of uncompleted pair correspondences is quadratic, although
the number of product completions is O(N log N). This is a capacity contrast,
not an obstruction to any number-theoretic existence conjecture.

## 3. Shift by the unique product-sum fixed point

The inherited binary product-sum equality has the unique solution2+2=2*2=4
(THM-361). Shift all transport nodes by this distinguished value:

    u=n-4 >=1.

An edge using the smaller parent a becomes exactly

    u -> a*u-(a-2)²,             2<=a<(u+4)/2.             (2)

Three especially useful lanes are

    a=2: u -> 2u,               u>=1;
    a=3: u -> 3u-1,             u>=3;
    a=4: u -> 4u-4,             u>=5.

Equation (2) is obtained by writing ab=a(n-a) and subtracting4. The fixed
point4 simultaneously removes the additive intercept from the first lane and
turns all the others into square defects. It is therefore a natural coordinate,
not an arbitrary pattern fit.

After identifying the a=2 doubling rays, the a=3 edge at an odd u>=3 projects
to

    u -> oddpart(3u-1).

This is a genuine **minus-Collatz** edge. The quotient also has edges induced
from even representatives, so the full quotient is a multivalued relation,
not just the deterministic minus-Collatz map. The a=3 branch at u=1 is absent
(the two parents would not be distinct). No claim about plus-Collatz convergence
or convergence of the minus system follows.

## 4. Factor four dissolves the weak-connectivity obstruction

**Theorem.** The underlying undirected graph of T is connected. In shifted
coordinates every u>=1 connects to1 by O(log u) edges, with all visited values
at most max(8,4u).

**Proof.** The doubling edges connect1,2,4. The a=3 edge3->8 joins3 to that
component. For every u>=5, the a=4 edge and two reversed doubling edges give

    u -> 4(u-1) <- 2(u-1) <- u-1.                         (3)

Hence all positive u belong to one component by induction. For the logarithmic
bound, halve even u. For odd u>=5 apply (3), then halve u-1. Within at most four
undirected steps the current u drops by at least a factor2. The base u<=4 has
paths of length at most4. All excursion values are at most4u initially and
only decrease as this algorithm proceeds. This proves the theorem.

This is a useful boundary on a tempting analogy: directed arithmetic ancestry
may be subtle, but weak connectivity here is completely elementary because
parent4 supplies a local decrement after dyadic identification. Connectivity
does not imply directed reachability from5. For example prime7 has indegree0
and is never reachable from5, though the two are weakly connected.

For numerical cutoffs, all nodes

    5<=n<=floor((N-4)/4)+4,             N>=12,              (4)

lie in the component of5 in T restricted to nodes<=N. The explicit paths
above stay below max(12,4(n-4)+4), which proves this sufficient bound.
The limit graph is connected while each finite prefix has newly arriving
prime/prime-square sources. “Connected eventually” must therefore specify
fixed old nodes, not all nodes in each finite graph.

## 5. Seed closure makes the transport edge an ancestry certificate

The summand lane proves, for 0<a<b,

    Cl_strict({a,b})={a,b} union (a+b+<a,b>),

where <a,b> is the ordinary nonnegative additive semigroup. When a,b>=2
are coprime, its last missing positive integer is ab. Therefore the transport
edge

    a+b -> ab

joins the first mixed-parent birth to the last missing target in the full
strict summand closure. The latency ab-a-b is exactly the classical two-generator
semigroup Frobenius number, and the transported endpoint itself is absent
from strict closure. With equal parents allowed, the closure is simply
<a,b> minus{0}; its final hole is ab-a-b. Thus deleting the diagonal shifts
this particular conductor boundary by the full sum a+b. These are precise
conductor statements, not an identification of generation time with value.

Example seeds3,5: first mixed target8, strict final hole15, weak final hole7.
The graph edge8->15 carries exactly those three numbers: sum, product, delay.
For noncoprime seeds, pass to the gcd lattice; ab itself is no longer the
last missing multiple (the last one is ab/gcd(a,b) when both reduced seeds>=2).

## 6. Open problems with useful separation

1. Classify directed descendants of a fixed source, first for5. No claim of
   cofinite descendants is possible because infinitely many prime sources
   remain unreachable.
2. Determine the size and birth/merger schedule of the giant weak component
   in finite T_N. The explicit old-node connected core is a lower bound,
   not a density-one theorem.
3. Keep a marked parent pair while combining transport edges. A node-level
   walk changes parent witnesses at every step; it is not an expression tree
   in either original operation without an extra compatibility condition.
4. Quantify the distribution of completion lags ab-a-b under a fixed summand
   birth. Its extremal boundaries are proved here; refined modular patterns
   can be tested without assuming independence of the two operations.

## Connection contract

Source: strict nonunit unordered parent pairs. Target: edges (sum,product).
Map: {a,b}->(a+b,ab). Preserved: pair identity, gap, both targets, causal delay,
all diagonal/unit conventions. Lost: parent order only. Sidecar: orientation
if a noncommutative operation is later attached; witness compatibility for
paths. Tests:4's loop, birth10's missing25,7's directed nonarrival from5,
and the coprime seeds3,5. No old theorem has been re-audited for this packet.
