# One-shot absorption and small-support rewiring

Status: **PROVED elementary reductions; FINITE-EXACT census through order 7;
OPEN unbounded coverage.** This is a research packet, not a Chan--Pak proof
or a canon promotion. The relevant statistic is the complete Hamiltonian-path
count H(T), not walks, rooted starts, or a selected path subfamily.

## The change of question

Ask for a small residual repair at a selected approximate witness. Do not
require that the repaired output be ready to execute the same operation again.
Two separate obligations remain: place an approximate witness near the target,
and certify a repair budget for that actual witness. Large subset-sum sets,
common scalar endpoints, and an arbitrary long interval do not settle either.

The most encouraging new finite fact is that reversing at most **three arcs
at once**, while changing H by 0 or 2 in absolute value, connects exactly the
consecutive blocks of the complete spectrum through seven vertices. This
includes every odd count from 23 through 61 in one actual component. The
restriction to one-arc moves stops the 23-component at 25. This changes the
operation, rather than searching longer for the same ready birth.

## Inheritance and concept board

The closest proved mechanism is
[THM-5122, shared-core exchange intervals and negative-inertia capacity](../../../../01-canon/theorems/THM-5122-shared-core-exchange-intervals-and-negative-inertia-capacity.md):
the affine Brown criterion and the unit-Lipschitz interacting cube criterion
are already proved. Neither is new here. The canonical hostile example is
[THM-5117, physical cut switches and primitive sparse cubes](../../../../01-canon/theorems/THM-5117-physical-cut-switches-reversal-cancellation-and-primitive-sparse-cubes.md):
arbitrarily large actual affine cubes have coefficient gcd 2 but minimum gap
50*37^(k-2). The corrected near miss is treating an independent collection of
present modules as a sum of independent corrections:
[THM-5124, additive rigidity on presence cubes](../../../../01-canon/theorems/THM-5124-additive-rigidity-of-the-path-count-on-presence-cubes.md)
rules out three active additive presence directions. The least-used relevant
sidecar is the identity of the actual tournament at an interval endpoint.

| Concept | Preserved predicate | Lost datum / obstruction | Test here |
|---|---|---|---|
| Absorption | Final exact count | Output need not be reusable | Bounded-jump lemma below |
| Additive completeness | Residual integer coverage | Actual simultaneous tournament realization | Inherited gcd-2 hostile cube |
| Reconfiguration | Actual tournament, not count quotient | Equal H does not identify vertices | 25-fibre separating witness |
| Locality budget | At most r arc changes per step | No reason r stays bounded as n grows | r=1,2,3 census |
| Scale | Maximum attainable count in an anchored component | Arbitrary interval location is uncontrolled | First missing count after 23 |

The earlier [THM-4960 twelve-bit moving-core intervals](../../../../01-canon/theorems/THM-4960-twelve-bit-moving-core-consecutive-odd-intervals.md)
and [THM-5000 coprime interface windows](../../../../01-canon/theorems/THM-5000-coprime-interface-composition-consecutive-odd-windows.md)
are positive interval constructions, but their exhibited baselines grow
exponentially while the window width grows linearly in their repetition
parameters. That parameter mismatch is why the present questions retain
target placement and maximum gap, not merely cardinality or gcd.

## A1. Bounded-jump terminal absorption

**PROVED.** Let K be a connected graph whose vertices are actual tournaments
and whose edges are certified operations. Suppose |H(U)-H(V)|<=2D on every
edge, for an integer D>=1. Suppose each vertex U admits terminal witnesses
A(U,j) with H(A(U,j))=H(U)+2j for 0<=j<=D-1. Then every odd integer between
any two attained counts of K is attained by a tournament.

For a target m, walk from a vertex below m to one at or above m. If m occurs,
stop. Otherwise the first crossing edge U,V has H(U)<m<H(V), and parity gives
1<=(m-H(U))/2<=D-1. Select A(U,(m-H(U))/2). No edge after this repair is needed.
Neither the repair output's readiness nor its rooted signature occurs in the
conclusion. All approximate vertices must carry the stated repair guarantee;
finding one absorber somewhere is insufficient. The same proof permits a
different D at each crossing edge, with the corresponding lower-end budget.
This is an elementary crossing argument proved here; no external absorbing
theorem is being invoked as a count-preserving theorem.

An unbounded family K_t with a common attained lower value M and maxima
tending to infinity would therefore cover every odd m>=M, **if** its residual
repair hypotheses hold. This is a sufficient certificate, not an equivalence
or a claim that the hypotheses currently hold for tournaments.

The operation graph may use rewiring or deletion. A birth-only nested chain
is a strictly stronger demand and is not required by this lemma. In particular,
the root session's persistent cyclic-vertex obstruction to a nested +2 chain
does not apply to the one-shot conclusion.

## A2. Correct gluing condition for interval components

**PROVED.** If every operation changes an odd-valued count by 0 or 2, the
image of each connected component is an odd interval. This follows by taking
a path between minimum and maximum and crossing each intermediate value.
The assertion applies to infinite components too, with intervals possibly
unbounded, by considering a finite path between any two chosen vertices.

Two components with a common scalar count do not thereby become one
component. A new count-preserving edge between actual vertices in the two
components merges them and realizes the union of their overlapping intervals.
The operation can alter every secondary signature unless a later composition
explicitly needs one. Thus terminal fibre transport is a weaker requirement
than rooted-context equivalence.

This formulation is elementary; the contribution of the census is to measure
the omitted transport and exhibit a small repair that supplies it.

## A3. Exact finite universe and findings

The universe contains every labelled tournament of orders 1 through 7:
1,2,8,64,1024,32768,2097152 tournaments, respectively; total 2131019.
We identify all isomorphisms using adjacent vertex transpositions and allow
source/sink births and their inverses while staying within the order cap.
Relabelling before deletion permits deleting any source or sink. At each
order, an allowed arc operation simultaneously toggles the indicated arcs;
its endpoints, not intermediate individual reversals, must satisfy the count
condition. There are no rooted, strong, or ready filters.

The H-count is computed by exact subset dynamic programming. Independent
permutation enumeration checks every labelled tournament through order 5 and
sampled codes at orders 6 and 7 (1361 checks). Every candidate relabelling and
neutral birth is also checked to preserve H. All errors raise exceptions;
the checker does not rely on assertions. Triple neighbourhoods are enumerated
from every isomorphism class, not one representative of a larger component.
The isomorphism-class totals are 1,1,2,4,12,56,456, agreeing with the known
small-order tournament census (the program verifies the numbers internally
only via the exact orbit enumeration).

| Moves beyond relabelling and neutral births/deletions | Components through n=7 | Component containing H=23 |
|---|---:|---|
| Single arc, absolute count change <=2 | 146 | [23,25] |
| Above, plus count-neutral two-arc changes | 40 | [23,35] |
| At most two arcs, absolute count change <=2 | 17 | [23,59] |
| At most three arcs, absolute count change <=2 | 10 | [23,61] |

In the last row there is exactly one component for each consecutive block of
the spectrum:

    [1,5], [9,19], [23,61], [65,105], [109,117],
    [121,147], [151,159], {171}, {175}, {189}.

These are odd intervals. In particular there is no duplicate component with
the same scalar interval in that row. At n<=6, two-arc steps already join
[23,33]; 35 is absent at that order. At n<=7 the three-arc 23-component stops
at 61 because **63 is absent in the entire finite universe**, not because of
a remaining reconfiguration obstruction at that boundary.

Relaxing the allowed *single*-arc count jump alone behaves differently:
absolute jumps <=4, <=6, or <=8 give a 23-component with maximum 45; jumps
<=16 reach 159. Such components can skip values and cannot directly establish
coverage. This is a hostile control for conflating large connected hulls with
unit-step interval coverage.

### Explicit 25-fibre repair

Rows encode the out-neighbour set as a bit mask, with vertices numbered 0..6.

    A = [44,17,2,6,13,30,63], H(A)=25
    B = [110,16,2,6,13,30,62], H(B)=25.

B reverses precisely the arcs {0,1} and {0,6} of A. In the single-arc graph
through seven vertices, A lies in a component with interval [23,25], while B
lies in one with interval [25,31]. The intermediate counts for reversing just
one of these arcs are 15 and 49. Thus neither individual step is neutral or
unit-sized, even though the simultaneous repair is exactly neutral.
Before the repair, A consists of a strong six-vertex component and singleton
source 6. After it, B is strong on all seven vertices: the move absorbs a
neutral singleton into the strong core while preserving H=25. This is beyond
an ordinal factor rearrangement. Two components of H=5 would need at least
eight vertices in total. The independent checker verifies the SCCs; the
parallel fibre audit also finds that the ordered two-bin path-cover count
changes from 608 to 548, so future response is not preserved.

The mixed two-bit response is 25-15-49+25=-14. The individual responses -10
and +24 are cancelled by this interaction. The repair works **because** its
interaction is retained, and is not evidence for independent additive ports.
Ordinal sums L->A->R and L->B->R give an infinite family of exact neutral
two-arc repairs, since both counts are 25 H(L)H(R). This generalization is
proved by concatenating the uniquely ordered component paths. Its interval
steps scale by H(L)H(R), so it supplies no new claim of unit-scale coverage.

### Why neutral transport alone stops at 35

Even if all equal-count tournaments were freely identified, a single-arc
unit move cannot take 35 to 37 in this finite universe: there is no such edge
anywhere. All single-arc component hulls are in the receipt; none contains
both values. This is an edge-existence obstruction, not merely a split fibre.
Allowing a two-arc unit change removes that obstruction and reaches 59.
There is likewise no two-arc 59-to-61 edge anywhere through seven vertices.
The three-arc witness below crosses that final local boundary:

    codes 32897 -> 37065, toggled bits 3,6,12,

with bit ordering (0,1),(0,2),(1,2),(0,3),... . The independent checker prints
the rows and counts. All three-arc neighbourhoods of all 456 types are checked,
so the resulting [23,61] interval is a finite exhaustive result.

## Outside mechanisms and their exact transfer boundary

**CITED, method connection only.** The absorbing framework replaces a prepared
small configuration after the bulk object is chosen. For example, Treglown
and Zhao define a matching absorber by requiring perfect matchings both with
and without the prospective leftover set; this is an actual realizability
predicate, not just divisibility. The corresponding tournament proposal is to
reserve a residual-count capability of the *actual approximate witness*.
The matching existence theorem does not imply a theorem about the number of
Hamiltonian paths. [Primary paper, Section 2.2](https://web.mat.bham.ac.uk/A.Treglown/exactmatchingFINAL.pdf).

**CITED, method connection only.** Conlon--Fox--Pham build interval coverage
from modular subset-sum expansion, density in intervals, and a final interval
combination. This suggests measuring residue diversity, largest gaps and
placement separately. The map from a set of repair weights to tournament
outputs exists only after all mixed responses are certified; THM-5124 and
the explicit -14 interaction show why it cannot be assumed. [Primary paper,
Introduction and Section 2](https://arxiv.org/html/2104.14766v1).

## Ranked questions with stopping tests

1. **Boundary repair, highest value.** For the complete order-n spectrum,
   let R(n) be the smallest arc radius for which unit-count rewiring plus
   neutral padding connects each maximal odd interval of the spectrum.
   The present finite statement is R(7)=3 (and R(6)=2). Is R(n) uniformly
   bounded, or can a bounded rule find a repair without already having a
   witness for the target count? The next cheap test is
   graph-isomorphism representatives at order 8, with the newly appearing
   63-fibre kept explicitly. No bounded-radius claim is made here.
2. **Growth, independent obligation.** Does the first absent odd count above
   23 tend to infinity with n? This is itself the remaining Chan--Pak
   coverage assertion, not a consequence of the connectivity census.
   Seek a construction that raises this boundary, with terminal rewiring
   only where needed. Proving connectivity without growth would not suffice.
3. **Budgeted absorption.** Produce approximate witnesses whose largest
   gap is at most twice their certified terminal correction radius plus 2.
   Measure both at the target scale. A single family of exponential
   baselines and linear windows fails this test before any new enumeration.
4. **Robust mixed-response reservoirs.** Seek signed exchange directions
   with controlled interaction and residue-diverse effective weights.
   First test simultaneous realizability and the normalized gap; do not
   restart a search for independent active presence modules, already ruled
   out. A failure of affine cancellation can still be useful, as the
   neutral 25 repair demonstrates.

**Calibration of locality.** An existential O((log M)^2) arc budget between
chosen witnesses for *already attained* adjacent counts at most M is automatic
from the inherited core bound. Delete singleton SCCs, pad the two cores to
K(M) vertices, and reverse all differing arcs in a single batch; there are at
most binom(K(M),2) of them. This already assumes both witnesses. Consequently,
a useful repair theorem needs a stronger feature: a uniform small support,
a preserved external boundary, or constructive access without knowing the
target witness. Neither this automatic bound nor R(n)-connectivity creates
an absent count.

The census classifies the specified finite operation graphs completely. It
does not classify arbitrary orders, prove a uniform local rewrite theorem,
prove Chan--Pak, or establish an asymptotic theorem for any prime sequence.

## Reproduction

From this directory, with a C++17 compiler and Python:

    g++ -O3 -std=c++17 corridor.cpp -o corridor.exe
    ./corridor.exe 7 1 1 > corridor.json
    ./corridor.exe 7 1 2 > corridor-neutral2.json
    ./corridor.exe 7 1 2 1 > corridor-pair-unit.json
    ./corridor.exe 7 1 3 1 > corridor-triple-unit.json
    g++ -O0 -std=c++17 corridor.cpp -o corridor-unoptimized.exe
    ./corridor-unoptimized.exe 7 1 3 1 > corridor-triple-unit.unoptimized.json
    python check.py

The additional jump2, jump3, jump4 and jump8 receipts use arguments `7 D 1`.
The optimized/unoptimized comparison is recorded in check.json. Executables
are disposable outputs and are not part of the research packet.
