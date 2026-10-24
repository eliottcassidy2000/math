# Repairing pinned charge holes from the other population

Status: **PROVED conditional all-size constructions.** Finite controls below
are abstract full-record models. Native applications require the original
record, electrical-code and resource authentication separately.

## 1. Inheritance and change of direction

The incoming `three_target570_20261024/THEOREM.md` uses a graded reservoir
while allowing core/helper populations to vary. Its `SHARED-PIN.md` shows
why a pin changes the available grade interval; `TWELVE-PIN.md` keeps all
twelve owner equations. Our preceding
`guarded_component_completion_20261024/theory/SEPARATOR_SECTIONS.md` retains
actual boundary records and computes exact sections and bridge costs.

The present operation keeps the TWO populations separately fixed. A helper
whose pinned charge image has holes need not itself be repaired: a receiver
inside the core can pay its changing charge. The new sufficient conditions
below construct the receiver from incomplete pieces, and then give a
terminating route. Neither step assumes connectivity of the final charged
fibre.

Concept board: charge holes; exact grade; shared owner boundary; protected
receiver; full-record refresh; endpoint cost. Anchor: repair externally pinned
families. Niche: complementary incomplete supplies. Wildcard: a charge-onto
union whose every relevant fixed grade fails.

## 2. Full records and the retained boundary

Let Q be a finite abelian group. It must encode the required physical rows,
not only a convenient statistic. Exact integer margins are kept as separate
grades. For three targets, separate target-owned charges are retained unless
an authenticated equivalence eliminates them.

Every option means a fixed FULL record, with actual source, target, path,
kernel lift, normalized resource, inverse/converse records, indexed pins and
clocks. Two different kernel lifts are not two capacity units. A boundary b
specifies every shared record/value needed to glue the interiors. All omitted
cross-factor relations must be justified; disjoint complete resource stores
are one sufficient justification for absence of resource collisions.

Charges and occurrence costs belong to uniquely owned slots. A shared owner
is not automatically a uniquely counted clock: the cost model here counts
the declared per-slot owner-clock occurrences. Distinct-global-owner costs
require a different ownership allocation.

## 3. Completing two incomplete graded pieces

At a FIXED actual boundary b and grades r,s, suppose two disjoint slot blocks
have charge sets A and B. Retain one authenticated witness a_hat per a in A,
and b_hat per b in B. The two sets need not be subgroups or cosets. Let D_q
bound the number of incompatible witness pairs of the form

    (a_hat, b_hat),  a+b=q.

This bound checks ALL cross-pins, resource collisions and higher-arity factors
against the retained boundary. In particular, D_q=0 is justified only when
every such witness pair really glues.

**Deficiency-repair lemma.** If

    |A|+|B| > |Q| + D_q                                      (1)

then these two incomplete pieces supply q at exact combined grade r+s.
There are at least |A|+|B|-|Q|-D_q compatible charge decompositions.
If (1) holds for every q, they give a complete section of Q.

**Proof.** Candidate first charges form A intersect (q-B), of size at least
|A|+|B|-|Q|. Each first charge determines exactly one second charge and the
chosen pair of full witnesses. At most D_q pairs fail. A surviving pair
glues and has the required grade and charge. Enumerating A and looking up
q-a in B constructs it. QED.

This is an elementary pigeonhole operation, not a claim of a new additive
combinatorics theorem. Its useful sidecar is the actual boundary and the
incompatibility count on witnessed pairs. It extends the earlier coset-only
sufficient construction to arbitrary supports.

For boundary-dependent grade intervals I_A=[l_A,u_A], I_B=[l_B,u_B], suppose
the same inequality holds for every allowed b, grade pair and requested q.
Every m in [l_A+l_B,u_A+u_B] has r+s=m in these intervals. Choosing such r,s
then applying the lemma builds a **conditional complete graded receiver**.
It is enough to certify one suitable split per m; all splits are not needed.

An extra restriction invalidating at most e_A chosen witnesses in the first
block and e_B in the second is tolerated whenever

    |A|+|B|-e_A-e_B > |Q|+D,

where D bounds remaining cross-incompatibilities. These are counts of lost
CHARGES in the retained witness maps, not numbers of pinned slots. One pin
can invalidate many witnesses, so a pin-to-loss estimate must be supplied.

The strict boundary is sharp. In C4, A=B={0,2} has |A|+|B|=|Q| and misses
both odd charges. With Q=C5, A=Q and B={0,1}, deleting the two candidate
pairs for charge0 gives |A|+|B|=|Q|+D_0 but still loses charge0. Conversely
A={0,1,2}, B={0,2,4} are both incomplete in C5 but supply every charge.
Giving each block three binary slots with these labels makes grades1 AND2
have three distinct charges. Their independently glued union is complete
at every total grade2..4, although neither block alone is complete there.

Boundary records cannot be averaged away. Binary slots x,y in C2 are each
charge-onto before imposing x=y; afterward x+y=0. Conditioning on the
actual shared value gives singleton supports, and (1) correctly fails.

### An exact difference-set test for one existing repair slot

At boundary b let A_(b,r) be the current grade-r support and H_(b,r)=Q minus
A_(b,r) its holes. Reassign a compatible actual binary slot from the outside
core into this receiver, with U-to-V difference label lambda. Its new
support is A_(b,r) union (lambda+A_(b,r-1)), so its new holes are exactly

    H_(b,r) intersect (lambda+H_(b,r-1)).

**Single-slot repair criterion.** The slot repairs every requested boundary
and grade if and only if its label avoids

    F = union_(b,r) [H_(b,r)-H_(b,r-1)].                    (4)

Indeed a surviving hole q is precisely a pair q in H_(b,r), p in
H_(b,r-1) with q-p=lambda. Each forbidden label retains an explicit
boundary/grade/hole witness. Each permitted label constructs sections by
trying the actual U case and then the actual V case. Empty hole sets
contribute nothing. This replaces a candidate-slot search by a difference-set
computation when the hole lists are short.

The criterion is exact for the TWO declared full records and their compatible
gluing, not for hypothetical records sharing their charge label. It requires
exact grades: omitting the r-1 shift manufactures population. Reassignment
changes which slots are called receiver/outside and occupies no extra unit.

## 4. Receiver-bank migration with both populations fixed

Partition the CORE into a receiver R of h binary slots and an outside set E
of e binary slots. Each slot has two fixed full records; a selected bit means
its V record. Subtract the all-U charge as a constant. Let lambda(I) be the
sum of selected slot labels. The core V population k is fixed throughout.

Let H be a connected finite graph of actual HELPER configurations. Every
vertex satisfies the declared helper population, all internal factors and
all frozen external pins. Its edges are explicit full-record transitions.
The graph ignores total charge, but not any of these predicates. Its charge
is eta(y). For example, an unrestricted fixed-grade free helper set has its
ordinary Johnson graph; a pinned helper graph may be much smaller.

A frame is (J,y), where J is an outside subset and y is a helper vertex.
Assume every outside subset is permitted as a frame and that, for every
frame, every grade r in [l,u] and every z in Q, there is an ACTUAL compatible
receiver record I of grade r with lambda(I)=z. All cross-factors between
I,J,y are checked. This conditional section may be constructed by Section3
or the existing separator compiler; it need not choose the same I at two
frames. No section is required at grades outside [l,u]. Assume 0<=l<=u<=h.

**Receiver migration theorem.** For

    l <= k <= u+e,

every nonempty fibre

    |I|+|J|=k,  lambda(I)+lambda(J)+eta(y)=q                 (2)

is connected by the following explicit full-record events. The helper
population remains fixed, and the core population remains k at EVERY event.

Set r0=max(l,k-e). It lies in [l,u], and 0<=k-r0<=e.

1. If |I|<l, remove l-|I| selected outside slots and refresh R to a grade-l
   section paying the new debt. There are enough selected outside slots
   because k>=l. If |I|>u, add |I|-u outside slots and refresh R to grade u;
   there is enough outside space because k<=u+e.
2. While receiver grade r>r0, add one outside slot and refresh R at r-1.
   The outside set cannot fill before r0 because k-r0<=e.
3. Now connect the two normalized outside sets by ordinary Johnson exchanges
   at grade k-r0. Hold the helper fixed and refresh R to pay every exchange.
4. Hold the outside set fixed and follow a path in H to the required helper
   endpoint, refreshing R after each helper transition.
5. Refresh to the exact receiver record at the end of the second endpoint's
   normalization, if needed, and reverse that normalization.

**Proof.** Each refresh grade lies in [l,u], and the section hypothesis gives
the exact residual q-lambda(J)-eta(y) at the actual new frame. Hence every
arrival satisfies (2), all full-record factors and both fixed populations.
Steps1--2 terminate. Johnson exchange decreases the outside symmetric
difference by two at every step. The helper path terminates by its given
connected skeleton. Step5 uses the retained literal second endpoint, not a
replacement representative. Each move is reversible on the given records.
This proves connectivity. QED.

Normalizing either endpoint takes at most u-l+1 events. If d_E is the number
of Johnson exchanges and d_H the chosen helper-path length, a safe bound is
2(u-l+1)+d_E+d_H+1 events. The last +1 allows a receiver-only reset even
when both skeleton distances are zero. Normalization events of zero length
are omitted.

This differs from assuming charged-fibre connectivity: H may have severe
charge holes; only its uncharged compatible skeleton is used. In Section3
even the receiver may be assembled from pieces with charge holes.

## 5. Exact event costs and source-specific choices

Let c_i be the actual occurrence cost of switching the two fixed records of
slot i. For receiver subsets I,I' of grades r,s, the exact event contribution
is sum_(i in I symmetric-difference I') c_i, and

    |I symmetric-difference I'| <= min(r+s, 2h-r-s).         (3)

The right side is attained by minimum possible overlap, so it is the sharp
grade-only support bound. Let W_R(t) be the sum of the min(t,h) largest
receiver costs (W_R(0)=0), with W_E similarly clamped at e. Then:

* entering from r<l costs at most W_R(r+l)+W_E(l-r);
* entering from r>u costs at most W_R(2h-r-u)+W_E(r-u);
* a consecutive-grade change r to r-1 costs at most
  W_R(min(2r-1,2h-2r+1))+max_(i in E)c_i;
* a normalized outside exchange costs at most
  W_R(2min(r0,h-r0))+the actual two outside switching costs;
* a helper edge y to y' costs at most
  W_R(2min(r0,h-r0))+the ACTUAL full-record cost of that edge.

Receiver-only refresh uses the same receiver term. The support bounds are
the corresponding bounds with all c_i=1, plus changed helper slots. These
are upper bounds for every pair of witnessing sections; hence they compose.
An exact source-fixed minimum section can improve them. Independent minima
between two *fibres* cannot be concatenated without keeping the middle
record or proving within-fibre connectivity. The earlier nonmetric hostile
remains applicable.

The outside E event must pay any attached forced controls as part of its
actual frame-edge cost; the displayed W_E bound applies only when switching
those E records is the entire outside work. Likewise helper shared-owner
components may be large. Width or charge rank does not erase their cost.

## 6. Complements, obstructions and physical scope

Complementing only the core is a coordinate change, NOT a charge-preserving
event. If lambda_all is the charge of all core labels and U is the set of
unselected core slots, (2) becomes

    -lambda(U)+eta(y)=q-lambda_all.

Thus the core labels are NEGATED, the helper charge is unchanged, and the
new core grade is h+e-k. Graded coverage survives negation/complement with
the appropriate grade change. Dropping this minus sign is invalid.

Charge coverage with margins forgotten is insufficient. A single binary
receiver labelled1 in C2 covers both charges in its union of grades, but at
fixed grade0 supplies only0. It cannot pay a helper debt1 while preserving
the core grade. At l or u, removing a pinned option may destroy completeness;
the new conditional band must be recomputed, not padded with imaginary slots.

The native incoming bank's code uses Y0..3, whereas this session's original
paired-current adapter uses Y1..4. The inherited relation
X mod2=Y0+Y4 mod2 in
`injective_native_blocks_20261023/electrical/ELECTRICAL_CODES.md` equation(3)
gives the explicit additive adapter

    (x,y0,y1,y2,y3) -> (x,y1,y2,y3,(x mod2) xor y0).

It is invertible and needs no margin correction. The packet independently
checks this adapter against all original physical potentials. One must apply
the map or recompute; silently identifying the displayed four-bit words
would still be incorrect.

The source is a compatible finite full-record supply; the target is a route
in one fixed-margin, fixed-charge fibre. The map keeps actual records at all
arrivals, charge and margins; it forgets no ports. The charge-table projection
does forget record costs and resources, so its section witnesses and factor
tables are indispensable. This proves no all-order native supply, no carrier
transfer, and no external conjecture. It permits later reuse of released
capacity and does not claim temporal nonreturn.

## 7. Native repair labels and up to ten fixed core records

The authenticated tables in `../erasure/receipt.json`, independently rebuilt
by `../erasure/audit.py` from original electrical potentials, cover every one
of the1140 triples frozen in the twenty cheap core sources. Their remaining
seventeen-source banks have only SIX cases with a hole at grade8 or9. Four
have no complete grade at all; two have complete grades7 and10 but incomplete
8 and9. This is a genuine nonmonotone grade obstruction.

`repair_labels.py` applies (4) to their exact holes at grades7..10. All other
triple erasures have complete grades8 and9, so contribute no forbidden label
for requested grades8,9,10. The union F has exactly32 labels. Every one of
the EIGHT old outside-core sources avoids F:

| Existing repair source | Paired label | Actual switching cost |
|---:|---:|---:|
|1532|725|4|
|1841|1592|3|
|1923|176|4|
|2087|1720|3|
|2553|1672|4|
|2835|1347|3|
|3015|487|3|
|3360|1173|3|

Thus ANY surviving one of these sources, reassigned from the outside into
the receiver, repairs every triple erasure at all grades8..10. This conclusion
uses actual original charge vectors and path costs. It does not assume that
every abstract allowed repair label has an available native source.

**Flexible finite-native corollary.** Keep external shared owner4, helper
population four and the original core population fixed. Freeze any THREE
cheap core records and any additional a of the eight outside-core records,
where 0<=a<=7. Each frozen option may be U or V, independently. At least one
outside source p remains free. Set e=7-a and c=its actual switching cost.
Every nonempty fibre connects at20 changed slots/85 primary slot clocks per
event. For free-core V count t in [8,10+e], every residual Q charge exists
and the stronger20/(49+c) bound holds:52 clocks if c=3 and53 if c=4.

**Proof.** Subtract the frozen count and charge as constants; frozen full
records never move. The free core has18+e slots, with repaired receiver R
of size18 and outside E of size e. Its receiver has seventeen costs2 and
one cost c. Apply the receiver theorem with l=8,u=10. The actual free helper
pairs form their Johnson graph, with each edge costing at most15 and the
externally pinned tree unchanged. Old core frames and resources are private,
as authenticated in `../native/CORE_PAYS_BOUNDARY.md`, so these unary pins
introduce no new cross-factor.

At canonical grade r0=max(8,t-e), receiver refresh costs at most34+c. The
helper step therefore costs at most49+c. Outside exchanges cost at most42+c;
consecutive-grade changes at most36+c; entry at most46+c (a safe bound).
Each changes at most20 slots. The exact grade formula computes tighter
source-dependent costs; for c=3,e=7 the entry bound is47. These estimates
hold for every selected section, so endpoints compose and restore exactly.

Outside [8,10+e], either t<=7 or (18+e)-t<=7. Two feasible endpoints then
differ on at most14 free core slots (cost<=56) and four free helper slots
(cost<=29). Their direct simultaneous event costs at most85, support18.
Complementing the free core is used only to bound distance, not as an
unguarded event. Boundary charge coverage is not asserted. QED.

The a=0,p=1841 case recovers the20/52 interior and20/85 global corollary
proved in `../erasure/ERASURE.md`; the new difference certificate removes
dependence on that one repair source. The most extreme allowed pattern has
ten frozen core cells (three cheap and seven outside), leaving an eighteen-
slot receiver. Its complete band is8..10.

“Freeze” means a fixed actual U/V record, enforceable at a clock distinguishing
its two paths, not an arbitrary new shared-owner relation or deletion of an
occupied slot. If all eight outside sources are frozen, or further relational
pins restrict the repaired receiver, this corollary supplies no conclusion.

## 8. Exact controls

`check.py` independently checks the deficiency threshold on all subsets of
small cyclic groups, enumerates graded incomplete-block examples, and runs
the terminating receiver compiler on every ordered pair of feasible states
in declared small conditional full-record models. It retains the actual
source at every refresh, verifies every boundary pin, both populations,
charge, exact switching cost and the grade-only bound. Positive and sharp
hostile controls are explicit. Run normal and optimized Python; the receipt
is a frozen exact experiment, not evidence beyond the theorem's quantifiers.
