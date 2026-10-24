# Equal-rank implication components and paid boundary routing

Status: PROVED, scoped conditional all-size routing theorem, independently
reviewed. No unrestricted native wiring claim or external novelty claim.

## Interface and inheritance

Use the full-record interface of
`../../compatible_charged_helpers_20261023/compatibility/NEUTRAL_FIBRE_LIFTING.md`.
All old binary records, normalized resource choices, source kernels, literal
converse records, shared owners and actual clock indices are retained.
Old-old compatibility is an implication poset P after contracting strongly
connected components. This requires both endpoint assemblies to be compatible;
ties or missing record predicates are not fabricated as tournament edges.

Fix visible components V containing every nonneutral charge component.
Let N=P\V be genuinely neutral for EVERY retained additive guard, including
exact margins. Write P_V for the induced order on V, including comparisons
whose paths pass through N. Neutral components have nonnegative costs w_x.
Do not remove a heavy neutral component when claiming a smaller event budget.

Rank |T|=k is a declared condition on the selected visible upset. It follows
from a fixed exact integer margin when every visible component has the SAME
NONZERO integer margin increment and the helper contribution to that margin
is fixed. Without these hypotheses, equal margin does not imply equal rank.
Additional nonconstant margin coordinates must
remain in the charge equation and the section test. This equal-rank
hypothesis is substantive: unequal component sizes need not admit exchanges.

For each rank-k upset T of P_V supply a nonempty set A(T) of complete helper
configurations h. Each (T,h) passes the full charge equation and has nonempty
neutral completion interval [L_(T,h),R_(T,h)]. Compatibility of the helpers
with both one another and the unchanged background is authenticated. In
particular a formal onto checksum is insufficient. A(T) may depend on T.

The closest mechanism is the earlier exact neutral boundary formula. The
hostile is a compulsory neutral shared-owner component of arbitrary weight.
The corrected near miss is integer charge generation without binary routes.
The new ingredient is a core route forced by the implication order itself,
so quotient connectivity need not be assumed independently.

## 1. Equal-rank upset exchange

**Lemma.** For any finite poset and any two upsets A,B with |A|=|B|, there is
a route of exactly |A\B| exchanges, each deleting one element of A\B and
adding one of B\A, preserving the upset property throughout. Each differing
component changes exactly once; this is a shortest route in two-component
exchanges, not a weighted minimax claim.

Proof. If A differs from B, choose an element i minimal in A\B. No element
of A lies strictly below i: such an element in A\B contradicts minimality;
such an element in A intersect B would force i into the upset B. Therefore
A\{i} is an upset. Choose j maximal in B\A. Every strict successor of j
lies in B and, by maximality, in A. It is not i, because i is outside B.
Thus (A\{i}) union {j} is an upset of the same size. Its difference from B
is smaller by one on each side. Iterate. QED.

Here x<=y means selection of x forces selection of y. One must remove a
minimal selected element and add a maximal missing one, not reverse them.

Every upset T of P_V has an old completion before helper constraints: up_P(T)
restricts to T on V. Helpers can remove that completion, which is precisely
why the nonemptiness assumption on A(T) is retained separately.

## 2. Exact cost of a lifted exchange

For accepted alpha=(T,h) and beta=(T',g), let

    delta(alpha,beta)
      = w(L_alpha\R_beta) + w(L_beta\R_alpha),

where L,R here are restricted to the neutral set N. Define

    c(alpha,beta)=d(h,g)+w(T symmetric_difference T')+delta(alpha,beta).

The inherited boundary lemma proves that c is the EXACT cheapest one-event
bridge between these two neutral fibres. It is attained by the explicit
two bounds in that lemma; helper resource multiplicities are not capacities.

**Constructive routing theorem.** Choose any exchange route T_0,...,T_m
from the lemma and any h_j in A(T_j), using the prescribed full helper
records at its two ends; these prescribed records must belong to A(T_0)
and A(T_m). Put

    B_route=max(max_(x in N) w_x, max_(0<=j<m) c((T_j,h_j),(T_(j+1),h_(j+1)))).

Then the prescribed full endpoint records are joined by legal packet events
of cost at most B_route. There are at most (m+1)|N|+m events, and all retained
additive guards hold at every arrival. If T_0=T_m but the helper records
differ, use a single helper-only quotient edge, counted as m=1 in the bound,
with the same visible set at both ends. If the helper records also agree,
only the neutral-fibre route is needed, with m=0. Empty maxima are zero.

Proof. First move inside the current neutral fibre to the explicit minimizer
for the next boundary crossing. Its components are neutral and can be
changed one at a time in upset order. Cross with the complete helper/visible
change and the exact forced neutral difference in one atomic event. Repeat,
then restore the prescribed final neutral record. The retained full records
authenticate every arrival; the equal-rank skeleton terminates after m
exchanges. QED.

This is stronger than using a worst-case charge alphabet to guess a block:
the route has a monotone visible distance and an exact recorded boundary bill.
It remains conditional on actual compatible sections at every required T.
It concerns simultaneous packets between feasible arrivals, not an ordering
of single-slot changes through potentially infeasible intermediate states.

## 3. A uniform sufficient condition that can be composed

If helper distances are at most D, each two-visible-component exchange costs
at most W, and all section pairs across an exchange (also same-T endpoint
resets) have neutral boundary cost at most C, then

    B=max(max_(x in N) w_x, D+W+C)

suffices for the whole rank-k guarded family. A same-T reset has visible cost
zero. D and C may instead be bounded on chosen sections plus the endpoint
resets; the theorem does not require arbitrary unused helper choices.

A simple universally valid bound is C<=sum_(x in N) w_x; useful improvement
requires the actual forced bounds. Arbitrarily large native neutral
components rule out a universal constant C without a structural hypothesis.
For a truly shielded fixed-occupancy r helper bank, D is bounded by the sum
of its largest 2min(r,h-r) record distances. Correlated control records belong
in the visible or helper cost; they are never free neutral lifts.

The separator compiler and protected-reservoir construction in this packet
are independent ways of certifying the section hypothesis. This theorem
provides their terminating core route whenever the equal-margin hypothesis
holds. It neither infers native supply at every tournament order nor
transports the order-eight harmonic bank to a new atlas.

## 4. Necessary hypotheses and small hostiles

1. With unequal margin weights (1,1,2) on a three-element antichain, margin2
   has states {first,second} and {third}. Their only bridge changes THREE
   components. Equal weighted margin does not imply a two-component route.
2. In a three-element chain, deleting a maximal selected element instead of
   a minimal one can leave a non-upset; the orientation matters.
3. Two helper configurations can force the same neutral component to opposite
   values. Even zero charge requires paying its entire weight. The previous
   native C3 example realizes this cost with actual kernel records.
4. A helper may return to its own reserved resource later, with every event
   counted, but each arrival has capacity one. Strict temporal nonreturn
   forbids moving and later restoring that helper and is not asserted here.

## 5. Check and map contract

`check.py` enumerates all naturally labelled posets through its stated bound,
their complete upsets and every equal-rank endpoint pair. It checks the route
against independent literal upset membership, distance and changed-coordinate
counts, and tests neutral boundary minimizers by independent enumeration.
No assertion statement is needed for correctness under Python optimization.

Map: full compatible packet -> visible upset plus complete helper record.
Preserved: charge, margin, shared-pin feasibility and bounded packet routing.
Lost: neutral owner choices. Sidecar: induced implication order, forced
interval bounds, weights and literal endpoint records. Decisive hostile:
one mandatory neutral flip of weight exceeding the purported budget.
