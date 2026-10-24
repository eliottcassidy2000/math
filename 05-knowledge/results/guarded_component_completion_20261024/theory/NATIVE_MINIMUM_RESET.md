# Exact paid resets inside two native shared-pin components

Status: **PROVED finite dynamic-programming mechanism; FINITE-EXACT native
certificate.** The optimum below is for one fixed core route and its full
supplied menus. It is not a uniform103-clock theorem or an optimum over other
core routes, extra records, or different source orders.

## 1. The complete actual state space

Use the independently authenticated family `twelve_two_components` in
`../native/NATIVE_COMPONENTS.md` and `../native/receipt.json`. Every one of the
original28 core slots and53 helper slots remains occupied. Twelve helper bits
are forced by actual shared-clock equations on core sources176 and1532;
helper2948 is fixed U. The protected40-helper bank has exactly three V slots.
Each control pattern contains six V slots, so helper occupancy is exactly9.

The ten controls attached to176 cost61 primary owner-clock occurrences when
that root changes. The two attached to1532 cost11. Changing both requires72;
this mandatory boundary work cannot be removed by choosing a cheaper charge
representative. All paths, embeddings, pins, normalized resource identities
and original electrical coordinates remain those of the native certificate.

There are exactly C(40,3)=9,880 possible reservoir triples. Every triple is a
genuinely compatible native choice: its complete option reservoir is disjoint
from the controls, core and fixed background. Its charge is the authenticated
Q=C128 x (C2)^4 code, with the declared margins fixed. No unlisted cross-cell
predicate is being dropped. Target U is1965 and target V is2234.

## 2. Source-specific cheapest sections

For an ACTUAL current reservoir triple R, define

    d_R(q)=min_{J subset bank, |J|=3, Sigma(J)=q}
                 sum_(s in R symmetric_difference J) w_s,       (1)

where w_s is the literal Hamming distance of the two supplied path records.
The records at U count as occupied choices, not absent slots.

The min-plus recurrence processes one slot at a time, retaining selected
count k, full charge q, minimum paid cost, number of minimizers, and an actual
witness mask. Choosing bit b adds its contribution to k,q and costs
w_s*[b differs from current bit]. States with k>3 can be removed because all
future count increments are nonnegative. At the end the k=3 table is (1).
This is the empty-separator source-fixed specialization of
`SEPARATOR_SECTIONS.md`, retaining an exact integer margin as well as Q.

**Correctness.** Every choice of the processed slots has one last bit and a
unique previous state. Taking the minimum therefore loses no feasible
candidate. Each slot is charged once, including unchanged U records. Ties
are counted as distinct full assignments, never as extra physical capacity.
Induction gives the exact cost and complete minimizing multiplicity at every
state; backpointers give a literal full record assignment. QED.

The checker independently compares all2048 outputs with direct evaluation of
ALL9,880 triples, including the count and selected representative of every
minimum. It does this for the five distinct source triples encountered in
the optimized route. Their worst required-charge costs are respectively
40,40,39,41,38. These are source-specific figures, not a worst-case theorem
over every possible current triple; the inherited unrestricted reset bound
remains42 for this reservoir.

## 3. A fixed route can be optimized without losing interior records

Retain the native core route

    {1394,2185} -> {176,2185} -> {1532,2185}
      -> {176,1532} -> {176,1841} -> {1477,1841}.

The required total encoded charge is134. The initial and final reservoir
triple is prescribed to be {650,676,760}. At each arrival the control bits
are fixed by the actual core bits. The complete residual-charge reservoir
layers have sizes

    1, 3, 8, 3, 3, 1.

There are consequently216 full candidate trajectories, with actual triple
records retained as vertices. Between consecutive layers the exact event cost
is the core-record cost plus forced control cost plus reservoir symmetric-
difference cost. All arrivals satisfy their shared pins and full charge.

A first dynamic-programming pass minimizes the maximum event cost. A second
pass minimizes total work subject to that optimal peak cap. Each pass retains
the actual triple at the intermediate layer. This avoids the invalid operation
of joining minima that choose different records of the same intermediate
fibre. It also avoids a different trap: keeping only a locally lexicographically
best `(peak,total)` label can discard a lower-total history whose larger peak
will be hidden by a later mandatory event. Fixing the optimal cap first makes
the second objective additive and sound.

Independent enumeration of all216 trajectories confirms the exact optimum:

    minimum peak = 103 primary owner-clock occurrences;
    minimum total among peak-103 routes = 415 occurrences.

One attaining reservoir sequence is

    {650,676,760}
    {862,1254,2393}
    {827,1254,2393}
    {1254,1400,1948}
    {650,1254,1391}
    {650,676,760}.

Its five events have costs103,89,91,43,89 and changed-slot supports
18,16,16,8,16. The inherited demonstrated route had peak114 and total433;
the uniform theorem's bound120 is unchanged. Greedily choosing the cheapest
next full section also attains peak103 on this example, but uses total427.
The global finite optimization saves another12 occurrences while returning
every literal helper endpoint record.

## 4. Guard, boundary and scope verification

The checker reads each shared-clock equation from the frozen native receipt
and evaluates it on the actual embedded paths at every optimized arrival.
It does not substitute equality of charges for equality of owners. Each
mandatory control flip is charged using its actual source path distance.

For every arrival it reconstructs every selected weight-two source-target
edge and all six original potential-current rows. Both source and target
diagonal contributions are included. The full original dyadic conclusion
uses the authenticated paired-code equivalence; six finite potential tests
alone are not declared an independent all-level basis. The physical guard
remains on the original order-eight carrier, not the ambient owner host.

The shared source/partner host construction and complete domain/resource
authentication are inherited from the pinned native and supply certificates.
Only reservoir records already in that complete reservation are selected.
There is no simultaneous resource reuse, missing slot, hidden kernel choice,
or released old pin. Later revisiting an own reserved alternative is allowed.
The final prescribed helper configuration is compatible with the final root
bits; the construction would not promise an inconsistent return.

Source-specific one-event sections can be used causally once a destination
is revealed. The fixed-route optimum additionally knows the entire six-arrival
core route in advance. These are distinct guarantees. No result here turns
the atomic packet events into single native metagraph edges or proves that
all compatible larger supplies admit the same clock bound.

Run:

    python -B 05-knowledge/results/guarded_component_completion_20261024/theory/check_reset.py
    python -O -B 05-knowledge/results/guarded_component_completion_20261024/theory/check_reset.py

The receipt retains every source-specific cost and minimizing witness, the
optimized full states and event ledgers, complete layer sizes and input hashes.
