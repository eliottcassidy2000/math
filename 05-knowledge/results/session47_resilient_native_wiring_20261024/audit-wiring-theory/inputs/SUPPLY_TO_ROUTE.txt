# From helper coverage to a costed route: the missing interface condition

Status: PROVED conditional all-size corollary. No new native helper supply is
asserted. The source/path/resource authentication and the actual compensation
charge of each helper set remain separate inputs.

## 1. Arrival coverage alone

At a proposed visible/core arrival i, suppose an actual helper bank has
p pairwise disjoint sets S_i,1,...,S_i,p, each of cardinality r, each producing
the required compensating charge. Selecting S means choosing the specified V
record at exactly those helper slots; **every other helper has its specified U
record**, not an omitted record. All choices preserve the required helper
margin. Each complete configuration has been authenticated internally,
against the fixed background, and with its inverse/converse and resources.

Suppose the only additional restriction relevant to choosing among these
candidates is a set Z_i of at most t helper slots forced U, and every candidate
avoiding Z_i is feasible at that arrival. In particular there is no unstated
forced-V requirement, higher-arity factor, resource collision or owner
predicate outside this stated criterion. If p>=t+1, some S_i,j avoids Z_i:
each forbidden slot can meet at most one member of a disjoint packing.

This is an arrival-existence certificate. It is not yet a low-cost edge or
route certificate. Different arrivals may require different members of their
packings, and their old neutral compatibility intervals need not overlap.
For r=4, nine disjoint representations suffice against any eight forced-U
slots; twelve suffice against any eleven. These are pigeonhole consequences
**conditional on such authenticated packings**, not assertions that a
particular supply has either packing count.

## 2. A sufficient route criterion

Fix a budget b and the visibility rule from Theorem 2 in `THEOREM.md`: every
charged component and every neutral component heavier than b is visible.
Suppose T_0,...,T_k is the intended core/visible sequence, and one can choose
at each step a complete r-selected helper configuration h_i so that:

1. The actual projected state alpha_i=(h_i,T_i) has a nonempty compatibility
   interval and satisfies every retained additive equation.
2. Consecutive visible changes cost at most K in the chosen cost units.
3. Consecutive chosen neutral fibres have exact boundary gap at most D.
4. Consecutive complete helper configurations cost at most H.
5. K+H+D<=b. Endpoint helper configurations are the required initial and final
   records, not merely configurations of equal charge.

Then every pair of prescribed full endpoint records in alpha_0 and alpha_k
has a legal b-budget route with exactly those helper endpoint records.

Proof. Formula (7) of the theorem gives each quotient edge cost at most
K+H+D<=b. All omitted neutral components are b-light and charge zero, so the
neutral-fibre routing theorem lifts the projected path to the prescribed
full endpoints. QED.

Conditions 1--5 can be checked for a selected sequence; a stronger uniform
bound over every allowed adjacent pair makes arbitrary compatible choices
safe. Merely finding a nonempty set of helper states at each arrival does
not provide either bound.

When each helper has two actual records, every h_i selects exactly r helper
V records, and all remaining helper records are U,

    changed helper slots = |S_i symmetric_difference S_(i+1)| <= 2r.

Thus H=2r is valid in changed-slot units. If every helper's actual U/V change
costs at most W owner-clock occurrences, H=2rW is valid in clock units. Use
actual visible costs K and neutral gaps D in the **same units**. As examples
of this conditional arithmetic, r=4 gives at most eight changed helper slots;
W=2 gives at most sixteen helper clocks. A two-slot visible change of four
clocks would therefore cost at most ten slots or twenty clocks **before** any
mandatory neutral boundary changes. This is not a universal cost-20 claim.

## 3. Why the boundary clause cannot be removed

Let one neutral implication component have weight m>=1, and take b=m. One
helper arrival forces it off and the other forces it on. Suppose the helper
change costs 1 and all other visible costs are zero. Both translated
compensation levels are populated, each compatible neutral fibre is a
singleton, and the neutral component is light enough to omit at budget b.
Nevertheless the exact crossing cost is 1+m>b. Ignoring the boundary term
would wrongly certify that route. Repeating each arrival option in a large
supply does not alter this forced old-record change.

This identifies the division of labor: disjoint helper packings can certify
survival of an actual supply after bounded exclusions; interval boundaries
certify compatibility of transitions between surviving arrivals; weighted
cuts and trees certify the optimal peak cost of the resulting complete state
graph. None of these certificates substitutes for another.
