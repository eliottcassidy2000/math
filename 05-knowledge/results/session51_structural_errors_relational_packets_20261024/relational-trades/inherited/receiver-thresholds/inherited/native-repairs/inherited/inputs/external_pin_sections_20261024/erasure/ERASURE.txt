# Receiver sections after additional external pins

Status: **FINITE-EXACT universal census on the stated twenty native slots**,
with the conditional routing consequences proved below. Arbitrary relational
pins and arbitrary tournament orders are not included.

## The receiver bank and its actual frozen records

The twenty two-clock core sources are

    176,1394,1425,1457,1477,1480,1549,1638,1691,1921,
    1974,2071,2176,2185,2581,2837,2905,2945,3000,3369.

Every source has the two fixed full U/V records of the inherited core menu.
Keep the actual XOR-tree external owner pinned to4, so helpers368,776,3139
remain U,V,V. Freezing a receiver slot here means fixing its entire U or V
record. It does not delete its occupied resource or its current contribution.
Its fixed charge and V occupancy are subtracted as constants before computing
the remaining variable section. The old core owner frames are private, so
these unary restrictions do not create unrecorded shared-owner conflicts.

A literal indexed-owner pin can implement either frozen bit: choose a clock
at which its two fixed paths differ and require the owner in the chosen path.
The full inherited target menu contains exactly those two records. No kernel
mark is exchanged while the target bit stays fixed.

The charge group is Q=Z/128 x F2^4, with X modulo128 and radical coordinates
Y1,...,Y4 modulo2. This census rebuilds the tables in those coordinates rather
than copying the incoming Y0,...,Y3 display.

## Exact conditional bands

For every set of e distinct frozen receiver slots, the following grades of
the remaining free slots cover ALL2048 charges:

| Frozen slots e | Free slots | Grades complete for every frozen set |
|---:|---:|---|
|0|20|5 through15|
|1|19|6 through13|
|2|18|7 through11|
|3|17|No common complete grade|

These are statements about every one of the1,351 erasure sets of size at
most three. Each of the2^e assignments to the frozen records merely translates
the required charge and grade by its fixed contribution, so the section
statements hold for all those assignments without identifying them.

Exactly four three-slot sets leave NO complete grade:

    {176,1425,2071}, {1425,2071,2837},
    {1477,1691,2905}, {1974,2071,2945}.

At grade8 their missing charge sets are respectively

    {1130}, {1008,1424,1602}, {7}, {1332}.

These integers encode charges, not source-class names. Two other three-slot
deletions have complete grades7 and10 but incomplete grades8 and9. Thus
coverage is genuinely nonmonotone in occupancy. The frozen receipt records
both sets and their holes; padding an existing witness is not a valid proof
of an intervening grade.

## Reassign one existing source to repair all three-slot cases

The old core source1841 is outside the twenty-slot bank. It has two genuine
full records, switch cost3, and paired charge1592. Reassign it to the receiver
bank; it remains the SAME occupied source with the SAME two reserved options.

**Repair theorem.** For every choice of three frozen cheap core slots, the
remaining seventeen together with1841 form an eighteen-slot receiver whose
grades8,9,10 each cover all of Q. All1,140 triple choices are checked, including
the six sets with deficient central sections. No port is added, duplicated,
or released by this reassignment.

Mechanism: if C_r(q) counts grade-r subsets of the seventeen-slot bank, the
new count is exactly

    C_r(q) + C_(r-1)(q-lambda_1841).

The two summands mean that the actual1841 slot uses U or V. They are disjoint
full-record cases. At grades8 through10 their sum is positive for every q,
even where either older section has a hole. This is completion by two
translated adjacent grades, not a fictitious copy of a missing resource.

## Routing consequence and exact boundaries

With three core records frozen, twenty-five core slots are free. The repaired
receiver uses eighteen, leaving seven outside core slots. Let t be the free
core V population, keeping both the frozen contribution and helper V
population four fixed.

For8<=t<=17 the receiver migration theorem applies with (h,e,l,u)=(18,7,8,10).
Its conditions hold in this native host: all receiver and outside core owners
are private, the tree stays pinned, and every record uses its own reserved
resource. The helper skeleton is the Johnson graph on pairs of the51 free
helpers. An edge switches two helpers costing at most15 clocks.

The receiver's total possible switching cost is17*2+3=37. Hence its refresh
plus a helper edge costs at most52 clocks and changes at most20 slots.
Entry into the complete band costs at most49 clocks; consecutive-grade and
outside-exchange events cost at most39 and45. These are valid upper bounds
for every choice of section, so exact endpoint restoration composes. The
result is **20-slot/52-clock connectivity throughout the interior band**.

For t<8, two feasible endpoints differ on at most fourteen free core slots
and four helpers. Their direct full-record replacement costs at most
14*4+29=85 clocks and changes at most18 slots. For t>17 use the same distance
bound on the at most seven unselected free core records. This is a distance
argument, not an unguarded complementation event. Thus **every nonempty fibre
under these three additional unary pins connects at20 slots/85 clocks**.
Full charge coverage is guaranteed in the interior band; it is not asserted
at all boundary grades. These budgets are upper bounds, not optima.

The receiver repair requires1841 to remain free. Further pins involving that
source, shared relations between surviving receiver slots, or additional
target menus require new conditional tables. This is three-slot unary-pin
resilience, not arbitrary shared-pin resilience.

## Reproduction and independent calculation

`check.py` builds integer graded multiplicities and removes a whole slot by
the exact recurrence D_r(q)=C_r(q)-D_(r-1)(q-lambda). It checks nonnegativity,
the total2^(20-e), complementation and fresh-count controls. The receipt
retains every coverage profile, all exceptional central holes and repairs.

`audit.py` imports no producer code and reads the original electrical
potentials. A separate Boolean group-algebra dynamic program computes each
erased bank afresh with bit permutations and cyclic shifts. It agrees on
ALL1,351 deletion profiles and checks ALL1,140 eighteen-slot repairs.
Run both scripts normally and with Python `-O`.
