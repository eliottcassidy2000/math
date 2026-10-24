# Compatibility intervals and exact neutral-fibre lifting

Status: PROVED, scoped all-size full-record compiler and move-lifting theorem.
FINITE-EXACT controls are in `check.py` and `receipt.json`. Integer trade span
alone is not used as a route certificate.

## 1. Inheritance and the new target

The closest mechanism is
`native_packet_decomposition_20261023/interface/ENDPOINT_FILTERS.md`: old binary
endpoint records with pairwise compatibility become an implication poset,
while exact per-owner currents and margins remain a joint additive factor.
The same session's `COMPENSATOR_LAYERS.md` retains genuine helper slots and
their translated charge levels. The present note supplies the missing
compatibility sidecar when helpers touch changing old pins or resources.

The canonical hostile is its algebra lane's vacancy path: a finite charge
group does not bound the necessary move when resource precedence is retained.
The corrected near miss is `balanced_native_20261023/SIX_TRADE_LATTICE.md`:
the 1,944 four/six-slot trades generate the entire integer charge kernel, but
this does not establish applicability inside a binary fibre. The least-used
sidecar here is a charge-neutral component's actual owner-clock weight. A
zero electrical charge is not permission to alter that component for free.

Concept board: actual shared pins; implication components; helper charge;
normalized resource capacity; neutral kernel records; atomic move cost.
Anchor: exact compatibility lifting. Niche: eliminating neutral components
without losing bounded connectivity. Wildcard: native C3 kernels whose
electrical silence conceals an arbitrarily expensive coupled owner change.

## 2. Actual records and the old implication order

Keep a fixed background and finitely many old slots, each restricted to two
supplied full records A_i,B_i. These records retain the literal source,
Hamiltonian path, normalized F/R resource reservoir, inverse/converse partners,
kernel labels, every fixed local guard, shared owners and integer clocks.
Both alternatives must pass the requested local guards, including any chosen
original per-slot phase. No hidden kernel lift is an extra resource.

Assume the all-A and all-B old configurations pass their pairwise old-old
constraints. Each pair forbids at most the mixed choices10 and01. Compile
these to implications and contract their strongly connected components. Let P
be the resulting finite poset, with x<=y meaning that selecting x requires y.
Thus selected components form an upset. Component x has fixed nonnegative
change cost w_x: for changed-slot cost it is the number of slots in x; for
clock-occurrence cost it is the sum of the two-record clock distances.

Let H be an abelian group, allowing Z factors. Component x has charge q_x in H,
the sum of its actual record differences. For electrical guards H contains the
ENTIRE retained per-owner row-current vector, including degree-diagonal terms
and all required original characters. Coincident physical source/target rows
are combined before testing. Exact margins are coordinates. A scalar source
checksum is allowed only after a separate complete equivalence theorem.

Supply a finite set of complete helper configurations h. Each includes all
helper slots, even inactive ones, passes its internal and background guards,
and has charge a_h relative to the declared baseline. Let d(h,g) be its exact
changed-slot or clock-occurrence cost. Old/helper interactions must be pairwise
record predicates. Higher-arity constraints require another retained factor;
the theorem does not erase them. Helpers are genuine supplied records and
resources, not formal labels carrying a convenient charge.

An atomic packet move is replacement of any set of records between two
feasible arrivals; its cost is the sum of their actual changed-slot or
clock-occurrence costs. The theorem concerns this declared packet graph. It
does not turn these updates into single edges of the original tournament
metagraph or assert a microscopic path through illegal intermediate records.

## 3. Each helper has an exact compatibility interval

For a fixed helper h, test both old records against the complete helper
configuration. At each old slot the allowed choices form a subset of{0,1}.
If a slot allows neither value, h has no old completion. Otherwise combine
the allowed sets over each implication component. Let F_h be the components
forced1 and Z_h those forced0. Any component forced both makes h infeasible.

Write up(S) for all successors of S, including S, and down(S) for all
predecessors. Define

    L_h = up(F_h),       R_h = P \ down(Z_h).

**Compatibility-interval theorem.** The old configurations compatible with h
are exactly the upsets U with L_h subset U subset R_h. This interval is
nonempty iff L_h subset R_h.

Proof. An upset containing F_h contains up(F_h). An upset avoiding Z_h cannot
contain any predecessor of Z_h. Conversely these two inclusions satisfy every
unary restriction and every old implication. Both L_h and R_h are upsets;
when L_h subset R_h either endpoint is a witness. QED.

This is a full iff for the pairwise compatibility carrier, BEFORE imposing
the separate joint charge equation

    a_h + sum_(x in U) q_x = 0 in H.                         (1)

It is not a 2-SAT encoding of that electrical equation. Intersecting several
sets of compatible fixed observations simply unions their forced1 and forced0
sets before taking closure. For actual helper alternatives one does not
intersect their observations: they occupy different arrival states.

## 4. Fix the visible coordinates and eliminate neutral ones

Fix an atomic budget b>=0. Retain as visible every component with q_x!=0 and
every component with w_x>b. One may retain additional components. Write V
for this visible set and N=P\V. Thus every x in N has q_x=0 and w_x<=b.

For a helper h and selected visible set T subset V, force T on and V\T off:

    L = up(F_h union T),
    R = P \ down(Z_h union (V\T)).

Discard the state alpha=(h,T) if L is not a subset of R, or if

    a_h + sum_(x in T) q_x != 0.                            (2)

Otherwise L intersect V = R intersect V = T. Put

    L_alpha = L intersect N,       R_alpha = R intersect N.

**Neutral-fibre theorem.** The full feasible records projecting to alpha are
exactly U=T union W, where W is an upset of the induced poset on N and

    L_alpha subset W subset R_alpha.                       (3)

In particular this fibre is connected by changing one neutral implication
component at a time, with every such move costing at most b and preserving
the complete charge equation.

Proof. Restriction of an upset gives (3). Conversely an upset W of N in these
bounds yields an upset T union W of P. A successor of x in T is already
forced into L. If x in W has a visible successor outside T, x would lie in
down(V\T), contrary to x in R. Successors in N are handled by upward closure
of W. All neutrality and charge assertions now follow from q_x=0 on N.

To route between two fibre upsets W,W', first remove minimal elements of
W\W', then add maximal elements of W'\W. Each step stays an upset between
the same bounds. No required element is removed and no forbidden one added.
Each differing component changes once, giving exactly |W symmetric_difference
W'| events and the claimed per-event cost. QED.

Retaining heavy neutral components is essential: electrical neutrality does
not let an indivisible shared-pin component of cost>b move inside a b-route.

## 5. Exact boundary distance between two neutral fibres

For two nonempty intervals alpha,beta of upsets on the SAME N, define

    delta_N(alpha,beta)
      = sum_(x in L_alpha\R_beta) w_x
        + sum_(x in L_beta\R_alpha) w_x.                    (4)

**Boundary-distance lemma.** This is the exact minimum weighted Hamming
distance between their members. It is attained by

    W_alpha = L_alpha union (L_beta intersect R_alpha),
    W_beta  = L_beta  union (L_alpha intersect R_beta).       (5)

Proof. A component forced on in one interval and forbidden in the other
must change, giving the lower bound (4). Both sets in (5) are upsets in N
and lie in the required intervals. Their two set differences are exactly
L_alpha\R_beta and L_beta\R_alpha, respectively. Thus they attain (4).
Nonnegative weights suffice; they need not be equal or strictly positive.
QED.

In particular the two fibres have a common member iff their forced bounds
satisfy L_alpha union L_beta subset R_alpha intersect R_beta. The weighted
distance can be zero without intersection if some mandatory changes have
zero cost; do not conflate these statements.

## 6. The quotient preserves every bounded-connectivity question

Make one quotient vertex for every retained feasible alpha=(h,T). Join
alpha=(h,T) and beta=(g,T') when

    c(alpha,beta) = d(h,g)
       + sum_(x in T symmetric_difference T') w_x
       + delta_N(alpha,beta) <= b.                          (6)

**Exact lifting theorem.** Two supplied feasible full assignments lie in the
same b-budget packet component iff their quotient vertices lie in the same
component of this quotient graph.

Proof of necessity. Project any actual move. The helper and visible parts
cost the first two terms of (6), while the neutral part costs at least (4).
Thus every projected change is a quotient edge or stays at one vertex.

Proof of sufficiency. For a quotient edge, first use the connected neutral
fibre to reach W_alpha in (5). Replace the helper configuration, visible
components and mandatory neutral differences atomically, arriving at W_beta.
Both arrivals are feasible by construction, and the cost is exactly (6).
Continue along the quotient route, using neutral fibre moves before each
edge, then route inside the last fibre to the prescribed final record.
All moves cost at most b. Full sidecars recover every actual owner, resource,
kernel, chosen phase and clock at each arrival. QED.

A quotient path with k edges lifts with at most (k+1)|N|+k events; each fibre
adjustment changes each relevant neutral component at most once. No minimum
event-count assertion is made. The quotient has at most |helpers|*2^|V|
vertices, even when its full fibres contain exponentially many neutral
assignments. Poset closures are precomputed once; each fibre and edge uses
set operations and its exact charge test. Discovery/authentication of helper
records and the complexity of the remaining visible charge states are not
removed.

This is a quantitative lifting criterion, rather than integer generation:
the helper/visible graph must use the precise boundary cost (4). Replacing it
by zero, or checking only that each translated charge level is nonempty,
can create a nonexistent route. The theorem also allows the neutral values
to change in the SAME event as helper/visible values; requiring literal
interval intersection would wrongly reject these admissible coupled moves.

## 7. Native all-size control: silent kernels can force a large event

Use the inherited actual family S_r=T_r -> C3 from
`native_packet_decomposition_20261023/workspace/ONE_BIT_WORKSPACE.md` §3.
For each prefix bit string a in{0,1}^r, its complete prefix-pinned menu is

    R(a,k)=P_a followed by (8r+k,8r+(k+1 mod3),8r+(k+2 mod3)), k in Z/3.

All three records consume the SAME normalized resource, indexed by a. Distinct
a consume distinct resources, and every record has the same target class.
The latter is authenticated by the inherited endpoint-fixed partner maps and
literal final C3 automorphisms; it is not inferred from a formal charge code.

Choose m+1 distinct prefix keys, possible when 2^r>=m+1. Pin every prefix
owner in each slot. Impose equality, across all slots, of the actual owner
at clock8r, the first final-C3 owner. Each local menu still has exactly three
records, but the full accepted assembly has exactly the three diagonal
phases k=0,1,2. Every target-dependent per-slot phase and every selected
class-incidence/current difference is constant. The actual shared owner
nevertheless forces all m+1 slots to change together, costing m+1 changed
slots or 3(m+1) clock occurrences. All three normalized resource lifts remain
one unit in each slot; no capacity is borrowed from their multiplicity.

Treat one slot as the helper and restrict attention to endpoint phases0 and1.
The old m slots form one neutral implication component of cost m in slot
units. The helper0 fibre forces it off, and helper1 forces it on. Formula(4)
gives neutral distance m, and the helper costs1, so (6) gives the exact bridge
cost m+1. For budget b=m the neutral component is light and each fibre is
individually connected (indeed a singleton), yet their bridge is forbidden.
The third available kernel phase cannot evade this lower bound: the retained
literal owner equality forces every nontrivial phase change to affect all
slots. This is an all-size native control of the compatibility mechanism,
not an assertion that a charge-changing helper family has been constructed.

Shared owner equality here uses the same actual source S_r and the same
clock; no unrelated local owner names are silently identified. Converse
transport uses reflected clock n-1-8r=2 and the corresponding literal owner.
The source order changes with r, so no original order-eight atlas harmonic
bank is transported to a different atlas by numerical coincidence.

## 8. Map, scope and reproduction

Source: actual binary old records, complete helper configurations, their
pairwise shared-resource predicates and full additive guards. Target: a finite
visible/helper graph with interval boundary weights. Preserved: EXACT
b-budget packet connectivity. Lost in the quotient: neutral owner choices.
Required sidecar: their implication order, forced bounds and full records,
which reconstruct an explicit route. Cheapest decisive test: two singleton
neutral fibres forced to opposite values by helpers; ignoring boundary cost
would invent a helper switch.

The result does not supply missing helpers, prove every charged fibre
connected, or infer binary feasibility from a lattice basis. It does provide
an all-size certificate once the actual compatible helper supply is given,
and it retains nontrivial shared-pin components instead of splitting them
into fictitious independent moves.

Run from the repository root:

    python 05-knowledge/results/compatible_charged_helpers_20261023/compatibility/check.py
    python -O 05-knowledge/results/compatible_charged_helpers_20261023/compatibility/check.py

The exact checker compares literal full-state budget components with this
quotient on its stated finite universes; checks interval membership and
distance by independent enumeration; and reconstructs the native C3 controls.
Explicit exceptions keep every check active under optimization.
