# Relational hole obstructions and robust witnessed sections

Status: **PROVED conditional all-size finite-record theorems**, with standalone
finite exact controls. No new all-order native supply, native census, global
observer, or physical electrical-code equivalence is asserted here.

## 1. Inheritance and the genuinely new boundary

The unrestricted binary repair identity and its hole-difference test are
inherited from `inputs/RECEIVER_BANKS.txt`, Section 3. The two different
branch-conditioned supports required by a receiver-reaching relation are
already explicit in `inputs/SESSION48_RECEIVER_REPAIRS.txt`, Section 5.
Full-tuple boundary descent, paid neutral fibres, and the conditional native
ternary-relation theorem are inherited from
`inputs/SESSION48_RELATIONAL_SECTIONS.txt`. The full-record separator and
paired-cost recurrences are inherited from `inputs/SEPARATOR_SECTIONS.txt`.
All seven archived source files are byte-identical Git blobs at the exact
commit and paths in `SOURCEPINS.json`; Markdown is archived as plain text.

The extensions proved here are:

1. A multi-branch **joint difference obstruction**, retaining the actual
   relation-conditioned domains, grades, and full witness choices.
2. Exact robustness against specified deletion atoms, with weighted vertex
   cover only for its justified whole-assignment special case.
3. An exact Hall-cut characterization of coverage under arbitrary relational
   tables with specified row/column counts. This gives a nontrivial family
   of arbitrary correlations for which coverage can be universally proved.
4. A simultaneous slot/clock Pareto interface and a component-aware exact
   lift, including explicit history expansion when temporal nonreturn is
   requested. This refines the inherited scalar-budget interface rather than
   claiming charge coverage alone proves routing.

The graph, cut, and flow arguments are standard finite combinatorics; no
claim of literature priority is made. Their role here is to retain the
sidecars that the native wiring problem cannot discard.

## 2. Full boundary, branches, and complete guards

Partition finitely many labelled slots into a frozen boundary B, receiver X,
and a repair block Z. A boundary b is its **complete actual record tuple**,
not merely its charge, target classes, coordinate envelopes, or unary
projections. Each record retains source, target, actual path, kernel lift,
normalized capacity unit, actual inverse/converse partners, all original
pins, and owner clocks. All other fixed data are part of b or an explicitly
retained parameter. Constants contributed by b are subtracted below.

Let Q be a finite (hence nonempty) abelian charge group and let g be an integer grade vector.
The group must encode the entire requested physical guard via a separately
authenticated equivalence. Alternatively retain the physical vector itself;
then the elementary union identities apply without finite-Q complements.
Exact populations stay in g, never in an accidentally reduced congruence.
Charges and clock occurrences belong to uniquely owned slots. Shared slots
are not counted a second time at a separator.

For each actual repair full assignment z, let c(z) be its charge and d(z)
its grade. Let E_b consist of every pair (x,z) that truly glues at b and
passes **all** factors: arbitrary relational pins of their full scopes,
capacity one for normalized units, background reservations, physical row
conditions not encoded in Q, and all other requested guards. For total grade
k write

    E_(b,k) = {(x,z) in E_b : g(x)+d(z)=k}.

No product assumption about E_b is made. Empty branch domains are permitted.
The same charge may have many full records, all of which remain distinguishable.
A repair block is reassignment of declared actual slots, not creation of
resources. One full kernel lift is a record, never an additional unit.

## 3. Multi-branch hole-difference theorem

At fixed b,k define the **branch-conditioned** support and holes

    A_z = {q(x) : (x,z) in E_(b,k)},       H_z = Q \ A_z.

Then the total support S and its holes H are exactly

    S = union_z (c(z)+A_z),
    H = intersection_z (c(z)+H_z).                         (1)

The intersection over an empty repair menu is Q. A section for q in S
consists of a branch z and a full x with (x,z) in E_(b,k) and
q(x)=q-c(z). Store that pair, not a charge-only representative.

**Proof.** Each allowed total assignment has exactly one repair restriction
z, and the residual receiver charge belongs to A_z. Conversely every member
of A_z has an actual compatible full witness by its definition. Taking the
complement in Q gives the second identity. Grades have already been imposed
inside E_(b,k), so no population is manufactured. QED.

For a nonempty menu z_0,...,z_m form the joint difference set

    D(H_0,...,H_m)
      = {(h_0-h_1,...,h_0-h_m) : h_j in H_j for every j}.

**Joint obstruction theorem.** Full coverage holds if and only if

    (c(z_1)-c(z_0),...,c(z_m)-c(z_0)) not in D(H_0,...,H_m). (2)

A forbidden vector retains a literal hole tuple h_0,...,h_m, and the missing
total charge is q=c(z_0)+h_0. Indeed the difference equations give
c(z_j)+h_j=q for every j. The converse follows by setting h_j=q-c(z_j).
For a single branch, D is the singleton empty tuple exactly when H_0 is
nonempty; the same criterion still applies. QED.

This is a joint obstruction, not a collection of pairwise differences.
In Q=C3 let all c(z)=0 and take holes {0,1}, {1,2}, {0,2}. Every pair of
translated hole sets intersects, but their triple intersection is empty.
All three branches repair the bank although no two do. The three branches
may be three full lifts of one resource: they are alternative arrivals,
not three simultaneous capacity units. No resource-resilience follows.

For two free compatible binary records with differences 0 and lambda and
grades 0 and 1, (1) becomes the inherited

    H'_k = H^U_k intersect (lambda+H^V_(k-1)).

Only when both branch restrictions leave the same old domain does this
reduce to the old H_k intersect (lambda+H_(k-1)) formula. Even equally
charged repair candidates can impose different factors, so their H_z must
be recomputed. An abstract permitted label is not evidence of a native
record having that label and those compatibilities.

### 3.1 Exact whole-branch deletion margin

Fix a nonempty requested charge set K subset Q. Let
C_q={z : q-c(z) in A_z}. If a new restriction only deletes whole repair
branches and leaves all other branch domains unchanged, it destroys q
exactly when it deletes every member of C_q. With unit deletion cost the
least number of branches that can destroy any requested charge is

    rho = min_(requested q) |C_q|.                          (3)

Thus all requested charges survive every deletion of at most t branches
if and only if t<rho. For nonnegative branch-deletion weights, replace
|C_q| by their sum. A zero value correctly detects an already missing q.
The statement concerns the specified deletion operation; an arbitrary new
relation inside surviving branches is not a whole-branch deletion.

## 4. Witness bans, resources, and exact covering obstructions

Fix b,k,q and retain the full witness set

    W_q = {(x,z) in E_(b,k) : q(x)+c(z)=q}.

Let T be a declared finite set of deletion atoms. For each actual witness w
record D(w) subset T: the atoms whose activation invalidates w. The model
requires that after choosing atoms A subset T, w survives exactly when
A intersect D(w) is empty; mutations of records or factors outside this
monotone deletion model need fresh witness sets.

**Witness-cover theorem.** A ban A destroys q if and only if A hits every
D(w), w in W_q. For nonnegative atom weights, the exact least destructive
cost is the minimum weight of such a hitting set, with infinity if no such
set exists. Hence all requested charges survive every ban of total cost at
most T_0 if and only if each minimum destructive cost exceeds T_0.

**Proof.** A charge survives exactly when one of its complete witnesses
survives. Negating this existential statement gives precisely the hitting
condition. The finite minimum is attained whenever a destructive set exists.
If W_q is empty, the empty set destroys it at cost zero. QED.

This covers deleting a normalized resource, a supplied record, a specific
branch, or a stated local option, provided its actual witness-incidence set
is retained. A normalized resource can invalidate many kernel-lift records
and many multislotted assignments. It is therefore normally a **grouped
atom**, not one bipartite vertex.

### 4.1 The justified bipartite vertex-cover case

When atoms are exactly whole receiver assignments x and whole repair
assignments z, each witness is an edge of the bipartite graph G_q=(X,Z,W_q).
Destroying q is exactly a vertex cover. With weights a_x,b_z its least cost
is the s-t minimum cut of

    s -> x, capacity a_x;
    x -> z for (x,z) in W_q, capacity M;
    z -> t, capacity b_z,

where M exceeds the sum of all finite vertex weights. A minimum cut cannot
cut an M-edge; its left vertices on the t side and right vertices on the s
side form the cover, and conversely every cover gives such a cut. This also
proves the equality without mistaking one charge representative for all
compatible full records.

For separate unweighted limits a,b, robustness has the exact rectangular
form: q can be destroyed if and only if there is U subset X with

    |U| <= a and |N_q(X\U)| <= b.                           (4)

After deleting U, precisely all remaining right neighbors must be deleted.
This is not equivalent to testing only the minimum *total* cover size.

Example of the resource quotient hazard: three disjoint witness edges have
whole-assignment vertex-cover number three. If the three repair endpoints
are three lifts of the same normalized resource rho, banning rho destroys
all three edges at cost one. Each individual arrival used rho once, so the
witness graph was capacity valid. The invalid step would be counting its
three alternative records as three independently bannable resource units.

These are certificates for **one requested charge at a time**. They do not
allocate several requested charge witnesses simultaneously, pack ports,
or provide disjoint witness inventories for many concurrent assemblies.
Such a problem requires its own joint capacity factors.

## 5. Arbitrary relational tables with prescribed marginals

Unrestricted arbitrary predicates cannot guarantee a nontrivial charge
range: the singleton predicate selecting one feasible full assignment
leaves just its charge and grade. This remains true even when every unary
projection or every unpinned receiver was previously onto.

There are, however, useful precisely delimited families of arbitrary
correlations. Fix b,k and the full compatibility graph E=E_(b,k). Specify
nonnegative integer row counts d_x and column counts e_z with equal total
M. Consider **all relation tables** R subset E satisfying

    degree_R(x)=d_x,       degree_R(z)=e_z.                 (5)

The row/column counts refer to numbers of allowed full pairs in the
predicate. They are not physical resource multiplicities or simultaneous
allocations. First verify this family is nonempty; universal assertions
over an empty family are not reported as coverage.

For a requested charge q let E^notq=E\W_q. For U subset X, V subset Z put

    deficit_q(U,V)
      = sum_(x in U) d_x - sum_(z in V) e_z
        - |E^notq intersect (U x (Z\V))|.                 (6)

**Relational Hall-cut theorem.** Assuming (5) is feasible on E, every
relation R satisfying (5) retains q if and only if there exist U,V with

    deficit_q(U,V)>0.                                    (7)

Thus all requested charges are retained by every such arbitrary relation
if and only if a cut certificate (7) exists separately for each q.

**Proof.** A relation avoiding q is an integral b-matching on E^notq with
these exact degrees. Use capacities d_x on s->x, capacity one on each
allowed x->z edge, and e_z on z->t. A flow of value M saturates every outer
edge and determines the required subset of allowed pairs; every such subset
gives that flow. Integer capacities admit an integral maximum flow: augment
along residual s-t paths by the integer bottleneck until none remains.
When no augmenting path remains, the residual-reachable cut has capacity
equal to the flow. Any cut can be written {s} union U union V, with capacity

    M - sum_U d_x + |E^notq intersect (U x (Z\V))| + sum_V e_z.

It is less than M exactly when (7) holds. Hence a positive deficit forbids
every q-avoiding relation, and absence of one produces an actual q-avoiding
relation via an integral flow. QED.

This is an exact alternative: either a witnessed unavoidable-charge cut or
a fully specified allowed relation that removes q. The same test on E first
certifies that the declared quota family itself is nonempty. Upper/lower
quota variants can be handled by the corresponding circulation construction,
but only exact quotas (5) are asserted and implemented here.

### 5.1 Positive and hostile full-table controls

Let X and Z each have four actual records, all pairs compatible at one grade.
In C2 let their charges be respectively (0,1,1,1) and (0,0,1,1). Set every
row and column count equal to one, so the relation is any of the 24 perfect
matching tables. Every such relation retains both charges:

- Avoiding charge 1 permits only equal-charge pairs. The three left charge-1
  records cannot fit into the two right charge-1 records. In (6), take those
  three left records as U and those two right records as V: deficit is one.
- Avoiding charge 0 permits only opposite-charge pairs. The same three left
  charge-1 records cannot fit into the two right charge-0 records, again a
  deficit-one cut.

This is genuine robustness to arbitrary non-product correlations within the
declared quota class. In each resulting table every fixed repair branch has
exactly one receiver record, so that conditioned branch is not charge-complete.
With two records on each side and charges (0,1),(0,1), both perfect-matching
relations have full unary projections, but one retains only charge 0 and
the other only charge 1. Full projections alone are not the quota theorem.

### 5.2 Sharp all-size binary marginal certificate

Suppose every full pair is compatible at one fixed total grade, with n>=1
receiver records and n repair records. In Q=C2 let a and b be the respective
numbers of charge-zero records. Let R range over all degree-one relation
tables, equivalently all bijective pairings of the **actual named records**.
Then the sharp minimum numbers of retained witnesses of the two charges are

    min_R |W_0 intersect R| = |n-a-b|,
    min_R |W_1 intersect R| = |a-b|.                        (8)

In particular every such arbitrary relation retains both charges if and
only if a!=b and a+b!=n. This is a scoped all-size sufficient-and-necessary
family, not universal coverage under unrestricted arbitrary predicates.

Proof: let t count zero-zero pairs. The four type counts are

    00: t,  01: a-t,  10: b-t,  11: n-a-b+t.

The exact possible t values are all integers from max(0,a+b-n) to min(a,b).
Every such choice is attained by matching the named records within these
four nonnegative type counts. Thus the charge-zero count is n-a-b+2t, whose
minimum is |n-a-b| and whose maximum is n-|a-b|. The complementary charge-one
count gives the second equality. Both extrema are attained. QED.

The exact marked-table multiplicities are retained as well. For each allowed
t define

    N_t = a! (n-a)! b! (n-b)! /
          [t! (a-t)! (b-t)! (n-a-b+t)!].

Choose which left records are paired with each right charge class, then
bijections to the named right records within each class. This gives N_t
actual tables. Consequently their two-charge counting polynomial is exactly

    sum_t N_t u^(n-a-b+2t) v^(a+b-2t),

and sum_t N_t=n!. These count relation tables on actual named records, not
unmarked native hosts, normalized resources, or simultaneous allocations.

Consequently, after selecting any such relation, deleting at most
min(|n-a-b|,|a-b|)-1 of its allowed **pair witnesses** still leaves both
charges whenever this bound is nonnegative. This is a whole-pair deletion
bound, not a resource-ban bound: one grouped resource deletion may erase
many of those witnesses. The n=4,a=1,b=2 example above has sharp minima one
for each charge, and no promised tolerance of even one arbitrary pair ban.

Requiring merely a nonempty set of neighbors at each record, with no upper
counts, is weaker. A q-avoiding relation with specified degree **lower**
bounds exists exactly when the whole graph E^notq meets those lower bounds:
selecting all its edges already suffices. A Hall obstruction cannot be
invented by silently imposing a matching requirement on an arbitrary
full-projection predicate.

## 6. Sections, exact simultaneous costs, and boundary lifting

All preceding results concern existence. For a current full assignment u
and a desired boundary/grade/charge section F, retain the finite witnessed
Pareto set

    P(u,F) = minimal pairs (s(u,v), c(u,v)), v in F,         (9)

where s is the number of changed labelled slots and c is the sum of actual
per-slot owner-clock disagreements. If the boundary changes, its slots are
included. Helper, outside and separately owned control changes must also be
included, exactly once. Paths, kernel lifts and literal endpoint records
are not replaced by convenient charge representatives.

**Exact source-fixed event criterion.** Under the atomic feasible-arrival
contract and an additional transition predicate if one is declared, an event
from u into F with budgets S,T exists exactly when (9), restricted to legal
pairs, contains a point (s,c) with s<=S,c<=T. Retain its actual endpoint v.
There is no claim that a legal atomic packet decomposes into legal individual
metagraph microedges.

Minimizing the coordinates separately is invalid. From binary source 000,
accept target records 100 and 011 with actual switch costs (10,1,1).
Their event costs are (1,10) and (2,2). The separate minima are 1 slot and
2 clocks, but there is no event satisfying both those budgets.

For arbitrary finite relational factors, the inherited paired separator
compiler can retain these Pareto pairs rather than one scalar minimum.
Combine disjoint owned costs by vector addition and discard dominated pairs
only at equal full separator tuples and equal retained charges/grades.
A dominated record may be discarded for a **single source-fixed query**;
forgetting it globally is not sound when later compatibility or history
can distinguish it. This is a finite exact compiler, not a polynomial bound
on its Pareto table size or on arbitrary separator width.

### 6.1 Component-aware boundary criterion

Fix S,T and assume the declared transition predicate is symmetric. Form the
undirected graph on all accepted full assignments with the legal atomic
budget edges just defined. Partition vertices by the desired full
boundary projection. **Split every projected fibre into its actual internal
budget-connected components.** Make a quotient vertex for each such
component, and put an edge when a pair of its full witnesses makes a legal
budget edge. Then two literal endpoints connect in the original graph if
and only if their component vertices connect in the quotient.

Proof: project a route for necessity. Conversely, before each quotient edge,
route internally to its retained source witness, cross to its retained target
witness, and continue. The final internal route reaches the exact prescribed
endpoint. Internal routes are paid events, not free record changes. QED.

For a directed transition predicate, split each projected fibre into its
strongly connected components and retain directed quotient edges instead.
The same necessity and sufficiency proof gives directed reachability. Weakly
connected fibre components do not suffice: they might not reach the selected
source witness of the next crossing.

If every fibre is internally connected, splitting is unnecessary. Otherwise
it is essential: allowed binary states 000,100,111,011 at slot budget one
have components {000,100} and {111,011}. Group 100 and111 into one projected
fibre B, and the other states separately into A,C. The unsplit quotient has
A-B-C, although no full route crosses. The minimum crossings chose different
records of B. Existence of all charge sections never repairs this defect.

### 6.2 Temporal no-reuse is an additional state constraint

Capacity one at every arrival is already inside E. It does not imply that a
released normalized unit is never reacquired later or that a departed
per-slot full record is never revisited.

For a strict no-reacquisition contract, expand the state to (u,L), where L
is the set of previously released normalized units. An event u->v is legal
only if no newly occupied unit is in L; update

    L' = L union (occupied(u) \ occupied(v)).

This counts release at the declared simultaneous-arrival boundary; a different
microscopic release convention needs its own transition rule. For per-slot
full-record nonreturn, instead retain each slot's departed records, forbid
arrival at one, and add the old record whenever that slot changes. Both
finite directed history-expanded graphs give exact reachability by ordinary
directed search.
The state can grow exponentially; this is not an efficient all-size routing
claim. Quotient splitting must retain the relevant history, not just u.

If every slot keeps one fixed distinct normalized resource while changing
kernel lifts, no resource is released, so strict resource no-reacquisition
holds automatically. Full-record nonreturn can nevertheless fail. The
inherited ternary code generated by (1,0,0,1,1),(0,1,1,0,1) is connected at
slot budget three, but 00000 cannot reach (1,2,2,1,0) without returning its
fifth record; direct budget four does suffice. The checker verifies this
nine-state distinction with an explicit history search.

## 7. Physical and native-use contract

A native user of these theorems must supply and authenticate E, q, g and
costs. The present packet does not infer those data from the abstract graph.
In the inherited order-eight application, Q=C128 x C2^4 is authorized only
by the original complete paired dyadic theorem, together with the literal
factor-two source and target rows, target diagonals, and fixed integer
populations. `inputs/ELECTRICAL_CODES.txt` documents that inherited gate.
It is not preservation of odd-primary characters, exact integer X currents,
individual marked-port phases, or a new carrier's electrical bank.

All requested native records and kernel lifts remain distinct until an
explicit quotient is proven harmless for the question. Normalize resources
before capacity checks; retain the actual inverse/converse closures without
multiplying their capacity. Converse transport applies to the **converse
input**: it transports factors, boundary tuples, guard labels, grades and
costs as authenticated. A relation need not be converse-invariant. Thus a
route has a converse route in the transported instance when the actual
converse bijection preserves its transition contract; no unverified
symmetry of the unchanged pinned instance is assumed.

A charge-complete receiver gives conditional all-size migration only when
its sections are compatible with every complete external frame used by the
route and the true frame transitions and boundary costs are supplied. The
inherited receiver migration proof can then be invoked, with the actual
sections from (1), the quota theorem, or a separator compiler. None of these
static witnesses proves that the hypotheses hold at every tournament order,
that arbitrary receiver-wide predicates preserve coverage, that one global
observer realizes all local records, or that ports may be reused.

The standalone checker exhausts the finite universes recorded in its receipt,
uses explicit exceptions rather than removable assertions, and runs normally
and under Python -O. These controls check the stated identities and hostile
models; the all-size statements rest on the proofs above, not extrapolation.
