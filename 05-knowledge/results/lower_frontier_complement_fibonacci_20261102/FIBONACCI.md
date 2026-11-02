# Fibonacci transfers and marked packing boundaries

Status: **PROVED** elementary matching identities; inherited source status
for the packing correspondence; **OPEN** extension to complete tournament
metagraph classification.

## A literal Fibonacci object inside the stability certificate

For an undirected graph G define its matching partition function
Z_G(t)=sum_k m_k(G)t^k. On the path with n vertices,

    Z_n(t)=Z_(n-1)(t)+t Z_(n-2)(t),  Z_0=Z_1=1.          (1)

Partition matchings by whether the final edge is used. At t=1 this is
Z_n(1)=F_(n+1). At the variance endpoint t=-1/4, the characteristic
polynomial has a repeated root 1/2, and

    Z_n(-1/4)=(n+1)/2^n.                               (2)

The exact two-state transfer is [[1,t],[1,0]]. At t=-1/4 its difference
from (1/2)I is nonzero and square-zero. Dropping that nilpotent coordinate
loses the factor n+1 even though the eigenvalue remains unchanged.
This is an elementary instance of the same information loss recorded in
the earlier F243 transfer: a clock or eigenvalue does not determine the
observable response. It is not a claim that the two algebras are isomorphic.

The order-ten stability polynomial is built from these very matching
coefficients. At ten vertices,

    P=V-M22+M222-M2222=1-Z_G(-q)-M22222

in multivariate notation, where Z_G(-q) weights each selected edge by
-q_e. The last term corrects the full perfect-matching coefficient.
This is a concrete source-to-target map: the uncertainty graph goes to
its matching polynomial; disjoint uncertain-edge incidence is preserved;
orientations and cycle contributions are lost and must be supplied by
the joint triangle/cycle terms in [the proof](STABILITY_10.md).

For more general graphs, boundary states replace the two path states.
The order-ten vertex-addition formula retains deletion matching counts
and the outside-edge counts Qij. These are the joint boundary responses
needed to evaluate the next neighbourhood. The resulting graph polynomial
is quadratic in that neighbourhood, avoiding repeated full pattern counts.

## What the older square packing work actually contributes

The inherited
[eleven-square routing packet](../eleven_squares_routing_bridges_20261007.md)
uses the supplied *The Fibonacci Geometry of Eleven Squares* as an
accepted source for its specified contact family. It retains

    Q=[[0,1],[1,1]], Q^2=Q+I,
    period lattice (I+wQ)(4I-Q) Z^2,

together with a canonical mark and the two clearance inequalities. Equal
unmarked tori can have different marks and enclosing side lengths. A
quarter-turn transports the frame AND the marks/clearance functionals.
This note uses that previously recorded scoped premise; it does not
re-audit the supplied packing paper or extend its optimum to every packing.

The transfer Q at t=1 is conjugate by swapping coordinates to the path
matching transfer in (1). That is an actual integer matrix relation.
At t=-1/4 the repeated-root degeneration changes its algebra; it is not
the same golden-ratio dynamics. The marks in the packing correspondence
and the endpoint occupancy in matching composition play the same specified
role: they retain information queried by the next operation. This is an
operation-level comparison, not a map of individual squares to tournament
vertices or a proof of packing optimality from Fibonacci counts.

The earlier separator correction remains decisive. Two edges glued at
one vertex form P3, whose independent-set polynomial is 1+3t+t^2.
Dividing the product (1+2t)^2 by 1+t does NOT give that polynomial.
Composition needs separate boundary-occupied and boundary-free channels.
See [the corrected boundary packet](../tournament_future_boundary_20261020/README.md).
Similarly, current square-path legality needs returning versus closed
leaf states once actual complements are enabled.

## A synthesis with a precise limit

The useful common construction is a boundary response with its operation
language: multiply compatible internal contributions, retain joint boundary
conditions, then minimize only after closing under the next legal operations.
The path matching transfer has two occupancy coordinates. The square25
complement observer has 33 minimal states. The five-vertex identity couples
variance, triangle response and a square remainder instead of summarizing
them separately. The packing model carries a mark and clearance guards.

These are different carriers produced by the same mathematical test.
Incoming [native Collatz joint ports](../collatz_hexagonal_growth_20261102/JOINT_PORTS.md)
give another literal reconstruction: an odd inverse tree with edge valuations
determines its shared halving tails. If two branches enter their odd LCA q
with valuations a,b, their first native meeting is 2^min(a,b) q, with both
clocks shortened by min(a,b). Thus the sidecar is sufficient because an
explicit reconstruction map exists. Independent edge subdivision would
lose the shared-tail identification. This is a comparison of source-preserving
composition, not an equivalence of Collatz convergence and tournament stability.

Complete metagraph classification would additionally need all legal paths,
full target relations, vertex transport and symmetry weights. In particular,
the square15 sources with full H=1 and H=9 forbid using the 33-state
restricted construction as an unrestricted enumeration recurrence.

Next concrete problems: prove the repaired order-eleven graph inequality
with its necessary wedge charge; classify path banks whose compatibility
graph is a forest rather than a matching; and find boundary transfers whose
guard-stable quotient remains small under an enlarged actual-path alphabet.
