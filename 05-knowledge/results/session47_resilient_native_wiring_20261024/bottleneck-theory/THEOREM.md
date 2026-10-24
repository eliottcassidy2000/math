# Exact bottleneck hierarchies for guarded native moves

Status: PROVED conditional all-size framework; FINITE-EXACT independent controls.
This is a synthesis of the actual-state obstruction with classical weighted
connectivity, not a new minimum-spanning-tree theorem. The compatibility input
is pinned in `SOURCEPINS.json`; the native guard conventions are not weakened.

## 1. What is being measured

Fix a finite set F of complete feasible assignments in one invariant fibre.
A state includes the actual source/path records, occupied normalized resources,
shared owners, fixed background, margins and all required original nodewise
current guards. Inactive helper records are included. Helpers required to
return have their original records at the two specified endpoints, not
necessarily at intermediate states. No formal charge label creates a resource.

Let E be the legal **reversible atomic events** between these states. In the
packet model an event is a simultaneous replacement between feasible arrivals;
it is not a claimed sequence of microscopic legal owner changes. Give each
edge its actual nonnegative symmetric cost c(e), such as changed slots or
changed owner-clock occurrences. An absent edge has cost infinity. The cost
model and transition rules are fixed before optimization.

For x,y in F define

    beta(x,y) = min over x--y paths P of max over e in P of c(e).

The empty path has cost 0; unreachable pairs have beta=infinity. Finiteness
means every finite minimum is attained by a simple path. This is **peak event
cost**, not the number of events, total work, or direct endpoint distance.
For b>=0 let G_b contain exactly the edges of cost at most b.

The conventions matter. Reversible capacity-one borrowing fits this graph.
A temporal prohibition on ever returning to a released port does not: legality
then depends on history and need not be reversible. The undirected/MST results
below cannot be applied to that different policy by forgetting its history.
Likewise preserving all nodewise harmonic currents is not the same as
preserving every separately marked ordinary port phase.

## 2. Hierarchy, cuts and tree certificates

**Theorem 1 (exact bottleneck hierarchy).**

1. x and y are in the same component of G_b iff beta(x,y)<=b.
2. beta is an extended pseudoultrametric:

       beta(x,z) <= max(beta(x,y), beta(y,z)).

   It is symmetric, nonnegative and zero on the diagonal. Distinct states may
   have zero distance if zero-cost moves exist. On a connected graph with
   strictly positive edge costs it is an ordinary ultrametric.
3. If A contains x but not y, put

       kappa(A) = min { c(u,v): u in A, v outside A, {u,v} in E },

   with the empty minimum infinity. For distinct x,y,

       beta(x,y) = max over x--y separating A of kappa(A).       (1)

4. On a connected F, let T be any minimum-total-weight spanning tree. Then

       beta(x,y) = max { c(e): e on the x--y path in T }.         (2)

   Consequently the sharp budget connecting **all** states of F is

       B(F) = max_(x,y) beta(x,y)
            = max_(nonempty proper A) kappa(A)
            = min_(spanning trees S) max_(e in S) c(e)
            = max_(e in T) c(e).                               (3)

   For a singleton, B=0. For a disconnected graph with at least two states,
   B=infinity and the tree minimum is over an empty collection. Apply the tree
   formula within each component, using a minimum spanning forest.

Proof. The first statement is the definition of threshold connectivity.
Concatenating two paths gives the strong triangle inequality. Every x--y path
crosses every separating cut, proving beta>=kappa(A). If beta=t is finite and
positive, take A to be the vertices reachable from x using edges of cost
strictly less than t. It excludes y and has no crossing edge below t, proving
kappa(A)>=t. If t=0 every kappa(A)>=0, and the opposite bound already gives
equality. If t=infinity, take the finite-edge connected component of x; it
has no crossing edge. This proves (1).

For (2), let e be a heaviest edge on the tree path. Removing e separates x
and y. A crossing edge of weight less than c(e) would replace e in T and
strictly lower its total weight. Thus this cut has kappa>=c(e), while the tree
path gives beta<=c(e). Formula (3) follows: a threshold graph is connected
exactly when it contains a spanning tree. Equivalently every nontrivial cut
has an edge within the budget. QED.

There is a useful maximality characterization. Among all extended
pseudoultrametrics u on F satisfying u(x,y)<=c(x,y) on legal edges, beta is the
largest pointwise. The strong triangle inequality along any path gives
u(x,y)<=max_e c(e); minimize over paths. The direct edge gives beta<=c.
The hierarchy therefore records the exact transitive closure of inexpensive
native moves, not a guessed geometric distance between charge labels.

The MST/single-linkage connection is classical. The bibliographic anchor is
J. C. Gower and G. J. S. Ross, *Minimum Spanning Trees and Single Linkage
Cluster Analysis*, Applied Statistics 18(1), 54--64 (1969),
[doi:10.2307/2346439](https://doi.org/10.2307/2346439).
The publisher's [issue record](https://academic.oup.com/jrsssc/issue/18/1)
verifies authors, title and publication details. The proofs above are supplied
here; no claim depends on access to the article's full text.

### What constitutes a sharp native certificate

For a proposed integer budget C:

- UPPER: give a spanning tree in each complete feasible fibre with every edge
  cost at most C, or another exact connectivity certificate with that bound.
- LOWER: give one nontrivial cut in one fibre and verify that **every** legal
  crossing event has cost at least C.
- AUTHENTICATION: certify the full feasible state set, all costs and transitions,
  actual resources, guards and supply boundaries. A cut in a selected subgraph
  is not a lower bound for a larger state graph unless its other crossings are
  also excluded or checked.

If the upper and lower match, C is exact. For integer weights, failure at C-1
is equivalent to a cut with no crossing edge of cost below C. This concerns a
minimum crossing edge, not the sum of crossing weights in a maximum-flow cut.
For a finite collection of nonempty invariant fibres, the universal threshold
is max_F B(F), not the connectivity of their disjoint union.

## 3. Contracting neutral fibres without losing a budget

The general graph principle is simple: at a fixed budget b, partition the full
state set into blocks that are each connected by b-budget paths **inside the
block**. Join two blocks when some crossing edge costs at most b. The quotient
and full graph have exactly the same connected components after pulling the
quotient partition back to the original states.

Necessity projects a full route. For sufficiency, route inside the current
block to the chosen endpoint of each quotient edge, cross it, and continue.
The witnesses used by two successive quotient edges need not be the same
state. Internal routing is the fact that permits changing witnesses.

Here is an explicit version retaining the complete compatibility sidecar.
Use the authenticated binary-record model from the pinned neutral-fibre note.
Old-old pairwise constraints become a finite implication poset P after
strongly connected components are contracted. Selecting component x means
selecting all its records. Write w_x>=0 for its actual change cost and q_x in
an abelian group H for its **complete** retained additive change. H can have
integer factors; no scalar checksum substitutes for unproved current tests.
The selected old components form an upset U.

A complete helper configuration h has charge a_h, exact change cost d(h,g)
and compatibility interval [L_h,R_h] of old upsets. These intervals compile
actual pairwise old/helper predicates; unencoded higher-arity constraints
invalidate this reduction. Full feasible states satisfy

    L_h subset U subset R_h,       a_h + sum_(x in U) q_x = 0.   (4)

Fix b and retain a visible set V containing every charged component q_x!=0
and every heavy neutral component w_x>b. One may retain more. Set N=P\V.
A projected state alpha=(h,T), T subset V, forces T on and V\T off. Let

    L = up(L_h union T),
    R = P \ down((P\R_h) union (V\T)).

Discard alpha unless L subset R and a_h+sum_(x in T)q_x=0. For retained alpha,
L and R have the same visible part T. Set

    l_alpha = L intersect N,        r_alpha = R intersect N.

Its full inverse image consists exactly of T union W, where W is an upset of
the induced poset on N with l_alpha subset W subset r_alpha. All x in N have
zero charge and w_x<=b.

**Internal routing proof.** For two such upsets W,W', remove minimal elements
of W\W' and then add maximal elements of W'\W. No removed component is
required and no added component is forbidden. Removing a minimal element
preserves upward closure; adding a maximal missing element does too. These
operations stay inside the interval. Each differing component changes once,
and each event costs at most b. Thus every retained projected fibre is
b-connected without pretending that all its owner changes cost zero.

For two projected fibres on this same N, define the boundary gap

    delta(alpha,beta)
      = sum_(x in l_alpha\r_beta) w_x
        + sum_(x in l_beta\r_alpha) w_x.                      (5)

Its value is the exact minimum weighted Hamming cost between interval members.
Indeed, every component forced on in one and forbidden in the other must
change. The bound is attained by the two upsets

    W_alpha = l_alpha union (l_beta intersect r_alpha),
    W_beta  = l_beta  union (l_alpha intersect r_beta).        (6)

Both belong to their intervals, and their respective set differences are
exactly the two sets in (5). Their unions and intersections are upsets, so the
witness construction respects the implication order, not merely the bits.

**Theorem 2 (exact weighted compatibility quotient).** Give projected vertices
the crossing cost

    c_b(alpha,beta) = d(h,g)
       + sum_(x in T symmetric_difference T') w_x
       + delta(alpha,beta).                                 (7)

A full pair of endpoints is b-connected iff their projected vertices are
connected using quotient edges of cost at most b.

Proof. Any full crossing has the helper/visible costs in (7) and neutral cost
at least (5). Conversely, (6) realizes (7) exactly; both endpoint assignments
are feasible. The internal routing proof lets one reach the appropriate
boundary witness before each crossing and the desired record at the end.
This is the general contraction principle above. QED.

A k-edge quotient route lifts using at most (k+1)|N|+k events. This is an upper
bound, not an optimal event count. Internal adjustment can add substantial
total work. Exact bottleneck compression does not preserve shortest-route
length or minimum total work.

## 4. A hierarchy needs budget-aware compression

Set V_b={x:q_x!=0 or w_x>b}. As b increases these visible sets decrease.
Each budget gives an exact quotient by Theorem 2. Pull all its connected
components back to the same original full-state set F. The resulting
partitions are exactly the nested components of G_b, so they preserve the
entire hierarchy of Theorem 1. At b<=b', forgetting the additional neutral
coordinates defines the compatible vertex map (h,T)->(h,T intersect V_b');
the two full-state projections commute with it. Connected components merge,
never split.

This does **not** license treating one aggressively compressed quotient as an
exact weighted graph at every smaller budget. For a fixed partition whose
fibres are internally connected with bottleneck at most r, write beta_Q for
the minimax distance of its exact minimum-crossing-cost graph. Then

    beta_Q(p(x),p(y)) <= beta_F(x,y)
                    <= max(r, beta_Q(p(x),p(y))).             (8)

The left inequality projects paths; the right lifts them. Hence the quotient
is exact at every threshold b>=r and for every pair whose beta_Q>=r. Below r,
internal neutral barriers may be invisible. Retaining all heavy neutral
components is a sufficient uniform way to avoid this error. It is not a
claim that every heavy component must be retained in every special instance:
a forced singleton fibre has no internal barrier to preserve.

### Three decisive controls

1. **Dropping heavy neutral components creates false connectivity.** One
   neutral implication component has weight 2, one helper configuration allows
   it either off or on, and b=1. The full graph has two states separated by
   cost 2. Dropping that component identifies them at one quotient vertex,
   falsely declaring them connected. Keeping it visible fixes the error.
   Electrical neutrality alone is insufficient even without helper changes.

2. **Interval intersection is too restrictive.** One neutral component has
   weight 1. Helper h forces it off, helper g forces it on, d(h,g)=1, and
   b=2. The intervals are disjoint. Nevertheless their exact crossing cost is
   1+1=2, a valid simultaneous event. Requiring a common neutral record would
   lose that edge and invent a separation.

3. **Zero boundary gap need not mean intersection.** With weight zero, the
   off-only and on-only intervals are disjoint but their gap is zero. All
   costs may be nonnegative, so the theory uses a pseudoultrametric rather
   than silently assuming positive definiteness.

## 5. Boundary gap is not a metric, and that is not a defect

On a one-component neutral poset of weight 1, let

    I={off},       J={off,on},       K={on}.

The exact gaps satisfy delta(I,J)=delta(J,K)=0 but delta(I,K)=1. Therefore
the triangle inequality fails. Distinct overlapping intervals also have zero
gap, so identity of indiscernibles fails even with strictly positive component
weights. These are minimum set-to-set gaps, not Hausdorff distances and not
metrics on the interval space.

The first zero-gap crossing uses J's off witness; the second uses J's on
witness. Changing between those witnesses costs 1. Neutral routing makes that
change legal at budget at least 1; it does not make it free. At budget below 1
the neutral component must remain visible, splitting the alleged intermediary.
Thus the correct object is the **budgeted state graph plus internal routing**,
not the raw gap alone. Its minimax closure is the pseudoultrametric of Theorem
1. A nonmetric edge-cost table is perfectly valid input to that theorem.

## 6. Helper supply: what is monotone, and what is not

Let one actual supply embed into another by fixing the added helpers at a
specified compatible baseline. Require that every old feasible state and
old event embeds with the same costs and all fixed-background predicates.
Then for fixed original endpoints x,y,

    beta_larger(i(x),i(y)) <= beta_smaller(x,y).               (9)

Every old route is still available. The theorem does not assert the required
embedding merely because the larger supply has more formal charge labels.
Removing helper options, while preserving the same endpoints and unchanged
costs on the surviving induced graph, gives the opposite inequality.

The universal worst case over **all endpoints of the enlarged supply** need
not decrease. Two old vertices joined by cost 1 have B=1; add a new vertex
whose two incident costs are 9 and the enlarged B is 9. The old pair still
has bottleneck 1. This metric three-state control prevents confusing easier
routing for fixed obligations with a smaller worst case over new obligations.

For a specified family R of allowed supply failures and endpoints surviving
each one, the sharp robust budget is max_(D in R) beta_D(x_D,y_D). Each beta_D
has its own exact cut lower bounds and path upper bounds. A supply packing or
charge calculation can help prove an upper bound, but it must authenticate
its actual compatible records and all intermediate states.

The [supply-to-route corollary](SUPPLY_TO_ROUTE.md) makes this separation
quantitative: disjoint representations survive bounded forced-U exclusions,
but a route also needs the actual adjacent neutral boundary costs.

## 7. Why a number 22 is a weighted statement, not 2 times 11

There is no established structural 11-threshold assumed here. Neither the
prime 11, a suggestive {2,3,11} analogy, nor an integer used to display a group
coordinate yields a native bottleneck 22.

If two cost functions on the **same graph** satisfy a*s(e)<=c(e)<=A*s(e), then

    a*beta_s(x,y) <= beta_c(x,y) <= A*beta_s(x,y).             (10)

This follows by taking the path maximum and then the path minimum. It proves
only bounds unless the scaling is exact or an independent optimal route/cut
closes the gap. As an abstract two-state control, eleven coupled changed
slots with all costs 2 have clock bottleneck 22; replacing one slot cost by 3
leaves the slot bottleneck 11 but changes the clock bottleneck to 23. These
controls do not authenticate new native paths.

An eight-slot event with costs (2,2,2,2,3,3,4,4) also costs 22. To prove that
22 is a sharp endpoint or universal threshold still requires the actual
weighted cut and matching upper certificate from Section 2. Counting the
eight changed slots or noticing a factorization of 22 does not do that work.

In particular, a packed parity label 11 in (C2)^4 is not the prime 11 as a
structural invariant. On coordinates (b0,b1,b2,b3), the invertible relabeling

    (b0,b1,b2,b3) -> (b0,b1 XOR b0,b2,b3 XOR b0)

maps binary 1011 to 0001. Applied consistently to every parity charge, it
preserves all charge equations and therefore the feasible graph and its
clock costs. A cyclic coordinate can be left unchanged. Thus a charge label
written (123,11) can equally be written (123,1); a genuinely certified peak
cost 22 is unchanged. No preserved mathematical map from an unrelated
appearance of 11 to this weighted cut has been supplied.

## 8. All-size content and limitations

Theorems 1--2 hold for every finite authenticated instance satisfying their
hypotheses; they do not assume order eight, a bounded number of slots or a
uniform helper density. They give a rule for certifying a number, not a
universal numerical bound independent of the native supply.

The inherited silent-kernel family provides a concrete all-size interpretation.
Its m+1 genuine slots have exactly three shared-owner diagonal phases; every
nontrivial phase change affects all m+1 slots and all three final C3 clocks in
each slot. The full three-state graph therefore has every edge cost m+1 in
slot units or 3(m+1) in clock units. One singleton-versus-two cut gives the
lower bound and any two edges give a spanning-tree upper bound. Its neutral
charge cannot reduce either value. Native path/resource authentication is
in the pinned inherited note, not newly replayed or transported between
network orders here. The small checker verifies only this stated diagonal
record model, not the inherited tournament construction.

Likewise an integer generating set for the charge kernel does not provide a
spanning tree in an actual bounded-capacity fibre. It omits precisely the
applicability information that determines the cut. The session-46 two-state
four-selected-slot fibre is an inherited example of this distinction; this
note does not rerun that native census.

The cohesive obstruction hierarchy is consequently:

1. actual source/path legality and guards determine the feasible state set;
2. compatibility and capacity determine which supplied labels coexist;
3. actual event weights determine the edge filtration;
4. neutral-fibre routing allows only justified contractions;
5. cuts identify the exact remaining bottleneck, and weighted trees certify
   its attainability.

This closes a conditional certification framework. It does not prove missing
helper supply, a universal small move budget, preservation of every ordinary
port phase, strict no-return buffered routes, or a new classification of
all tournament metagraphs.

## 9. Reproduction

Run `python check.py` and `python -O check.py` in this directory. The checker
uses explicit exceptions, not Python assertions. It independently compares
minimax paths, threshold components, all cuts and spanning-tree certificates;
compiles and enumerates small compatibility quotients; and executes the
counterexamples. `receipt-normal.json` and `receipt-optimized.json` report the
exact finite scopes. `SOURCEPINS.json` authenticates the inherited text inputs;
`MANIFEST.json` pins the delivered packet. Tests do not substitute for the
all-size proofs above.
