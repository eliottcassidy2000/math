# Source-weighted routing after finite feasible coverage

**PROVED conditional inequality**, using the FINITE-EXACT native graph and
trees. This does not prove scale-uniform congestion or Collatz entry.

The primary census records a rooted breadth-first tree of depth at most5 in
every six-connected fibre. The independent sparse-completion checker builds
trees of depth at most5 in every margin4/5 fibre using only the1944 small moves
and220 selected bridges. Complement gives the remaining grades. Therefore
every pair of states in a full supplied fibre can be joined through a root
in at most10 events, each changing at most8 slots and22 actual owner clocks:
total work at most220. The largest fibre has20150 states. These are upper
bounds, not claims of minimum diameter or least-work paths.

Let pi be the ACTUAL positive normalized state law on one fibre F, and let P
be any reversible kernel with pi-stationary edge conductance

    Q(x,y)=pi(x)P(x,y)>0

on each edge of the selected tree. For a tree edge e, delete it and denote one
side by A_e. Route the demand pi(x)pi(y) along the unique tree path. For either
orientation of e the unweighted demand is exactly

    pi(A_e)(1-pi(A_e)).

Every routed path has at most10 edges, so the length-weighted congestion is

    Lambda <= max_e 10 pi(A_e)(1-pi(A_e)) / Q(e).

For clarity, apply Cauchy--Schwarz to
f(x)-f(y)=sum of edge increments on the path, multiply by pi(x)pi(y), and sum.
The left side divided by2 is Var_pi(f); the right side is at most Lambda times
the Dirichlet form of P. Thus the spectral gap is at least1/Lambda. The same
cut demand formula provides an exact test for a proposed source-weighted
allocation; uniform cardinality is not substituted for actual pi mass.

This bound necessarily retains small conductances. A two-state fibre with
transition probability epsilon in both directions has gap2epsilon although
its graph is connected. No connectivity theorem can remove that obstruction.
Likewise moving most probability behind a small cut can worsen the weighted
bound without changing any charge or normalized resource.

The physical electrical multiplicity in the native packet remains2, and all
source owners/weights of those edge records remain unchanged. The statistical
law pi in this paragraph is a different object. Arbitrarily changing native
electrical source weights would change the conservation equations and is not
licensed by this routing result.

## Reusable move and boundary

Separate three tasks: certify a complete native fibre, connect it with feasible
rules, and compare the actual measure across those rules. The present finite
closure solves the first two for one inventory and supplies the exact
weighted-cut interface for the third. The earlier optimal product-flow theorem
in `balanced_native_20261023/PRODUCT_FLOW.md` already shows why refresh cost
can grow quadratically even on an easy cube. The current finite constant220
cannot be exported to larger inventories without a new supply/depth theorem.

For Collatz the retained source, refinement weights and integer event clock
have the same bookkeeping roles, but no map carrying these tournament moves
to positive ROOT certificates has been established. Similarly a geometric
colouring application must retain the relation of extendible boundary colour
assignments; connectivity of a quotient or preservation of a graph minor is
insufficient. These are typed transfer tests, not claimed equivalences.

