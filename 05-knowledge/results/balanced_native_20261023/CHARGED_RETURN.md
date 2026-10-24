# One genuine helper repairs the new eighteen-slot circuit

**FINITE-EXACT / PROVED from complete finite domains.** Add source1638 to the new eighteen-source circuit, with its U path at both endpoints. The complete nineteen-slot pinned and margin-restricted fibre has42 states. All42 are connected by events changing at most six slots. Budget four has23 components and separates the marked endpoints. Thus six is the sharp changed-slot budget for this complete inventory.

The helper is real source/path capacity. It is not another lift of an already occupied tile. All38 U/V alternatives in this inventory have distinct normalized resources; their closed reservoir has152 ports. Every arrival occupies19 resources, or76 closed ports. The single raw extra path at source1923 still leads to1559 and remains excluded by the exact zero target count. Shared pin roles and all private fixed owner clocks remain unchanged.

## Route, current and clock

Use the A,B endpoints in NATIVE18.md, with source1638 initially and finally at U=1965. A shortest route in the complete six-slot move graph has four atomic events, each changing six slots. Their owner-clock costs are14,15,15,13, totaling57. All full arrival records are retained in route18.normal.json.

At each of the three intermediate arrivals the helper is at V=2234. Relative to the initial packet, its U-owned X current change is-104, while the old eighteen slots contribute-152. The total is-256, hence zero modulo256. The V-owned changes are the negatives, and every source-row change vanishes. The four binary character debts also vanish. The old contribution includes its diagonal degree term: when one old V slot is replaced by the helper, ignoring the changed old target count would report the wrong physical current.

At the terminal packet the helper returns to U, and the old slots' U-current difference is-256. Thus this route retains every original dyadic harmonic bank, via the exponent256 theorem. It does NOT retain the stronger, different requirement that the same raw integer X current remain exactly zero. Indeed its endpoints already differ by-256, so that stronger requirement would forbid the requested return before any route search.

Four events are shortest at budget six: the eighteen old choices must all change, and any active helper must leave and return, requiring at least twenty changed-slot occurrences. Three events could supply at most eighteen. The saved route attains four events. This argument concerns event count;57 is the cost of this recorded shortest-event route, not a claimed globally minimal owner-clock cost.

## Complete helper comparison

Among the ten remaining sources of the supplied28-source pool, three individual helpers connect these endpoints at budget six; eight do so at budget eight. Source1638 is the first successful budget-six helper in source order. With all ten helpers supplied but at most one at V at a time, there are exactly217 states. Their budget-six graph is connected; the budget-four graph still separates the endpoints. Thus no choice among these ten translated single-active layers, even with handoffs, repairs the four-slot obstruction.

For the chosen single helper the42-state domain is the complete nineteen-slot binary fibre, since at most one helper can be active automatically. The layer enumeration is independently checked by all binom(19,9)=92378 target subsets. `route18_audit.py` also enumerates all8! owner permutations per source and directly computes each selected packet's physical Laplacian action on X and64Y0,...,64Y3. It checks every target and source row at every arrival, including diagonal contributions; it does not trust only a scalar total or the producer's charge table.

## A stronger four-slot obstruction, even with all helpers active

There are exactly27 possible charge-safe four-slot trades in the full28-source binary supply, found from all378 distinct source pairs. Their signed incidence matrix has rational rank23, hence a five-dimensional space of linear invariants. The following integer coefficients define one such invariant f(I)=sum_(s in I) c_s on the set I of sources assigned to V:

|c_s|Sources|
|---:|---|
|-3|1691,3000,3369|
|-2|1841,2176,2581,2945|
|-1|1480,2553,3360|
|0|1394,1457,1549,1921,1923,2071,2087,2185,2835,2905|
|1|176,1532,2837|
|2|1477|
|3|1425,1638,1974,3015|

Every one of the27 four-slot trade vectors has exact dot product zero with c. The initial and terminal values are-3 and-1. Distinct source labels rule out two-slot trades, and retained target margins require even support. Therefore budgets four or five cannot connect the marked endpoints anywhere in the full28-source binary supply, even with arbitrarily many of its helpers active. The six-slot route already constructed proves the sharp budget six in that larger context as well.

The concurrent packet's earlier small invariant takes value-6 at BOTH new endpoints; it does not detect this obstruction. `four_move.py` supplies a different exact invariant and retains the full rank computation. Its small coefficients were found with an optional integer optimizer, but the saved verifier uses only exact integer products; no numerical infeasibility or coefficient-minimality claim is made. This is a concrete reason to retain an invariant space rather than one successful scalar summary.

## Integration with the concurrent charged-layer result

Concurrent commit `a793424e20`, `native_packet_decomposition_20261023`, repairs a DIFFERENT sixteen-source circuit with two helpers1532/1923. Its eighteen-source enlarged fibre has16 exact-integer-current states and a five-event route of work76. Our eighteen-source *primitive* packet and nineteen-source *repaired* packet have different source sets and different current scopes. The shared numbers do not identify the objects.

The concurrent translated-layer theorem is the mechanism used here: an active helper h supplies the complete old subset layer with cardinality k-1 and charge c-ell(h), together with its actual resource and margin. It moves outside the original endpoint-mixture interval while keeping the enlarged physical guard. This explains why a genuine charged helper succeeds where an electrically neutral resource bit cannot.

The next all-size question is not whether abstract charge translations exist. It is whether enough compatible native layers exist with retained pins and distinct capacities to bridge every surviving small-move component. The three-target fibre in MULTITARGET.md gives a second precise supply target.
