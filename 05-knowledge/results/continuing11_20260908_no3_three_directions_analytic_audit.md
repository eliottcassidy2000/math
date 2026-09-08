# Independent analytic referee: third-direction deletion bonus

**Status: PASS / PROVED analytical implications within the stated sampling
and original-cell-loss scope.** No mathematical repair is requested. This
referee read the complete frozen proof and the inherited two-direction
avoidance/concentration proof and audit. It did not import or execute the
producer and does not certify the new finite census; the root's independent
finite referee handles that separate obligation.

The exact audited candidate is
`continuing11_20260908_no3_three_directions.md`, raw SHA256

    8354bc56f6c8dc663d76a6b728c7bb6fbe3df1ae18b6373210d9d4adc7936f4c.

The inherited proof and audit are
`05-knowledge/results/continuing10_20260907_no3_repair.md` and
`05-knowledge/results/continuing10_20260907_no3_repair_audit.md`.
Their already audited uniform conditional two-family mean is used with its
original coefficient gamma. The new positive coefficient is

    delta_3 = exp(-416/45)/4.

## 1. Deterministic deletion and Boolean dual

Let U contain the original cells whose original difference and sum lines
both have occupancy at most two. If a deletion set clears the two old
families, restoring its deleted cells from U cannot create an old violation:
originally safe lines remain safe, and originally overfull lines contain no
restored cell. Thus any three-family deletion has at least tau_2 deletions
outside U. Since the slope-two lines partition U, it also deletes at least
sum_q (|U intersect q|-2)_+ cells inside U. These are disjoint charges, so

    tau_3 >= tau_2 + J_3 >= tau_2 + I_3.

This argument requires no stochastic independence and applies to every
finite distinct-cell board. It avoids precisely the double-counting error
in the unrestricted third-direction excess.

The two-family flow has unit cell edges and vertex capacities two. A cut
with old source-side vertices A,H has the displayed capacity; the cells
outside N((D minus A) union H) are exactly E(A,S minus H). Subtracting the
cut from |B| therefore gives c(Q)=|N(Q)|-2|Q|, with every Q represented.
This validates the exact Boolean dual for tau_2.

For a fixed old line selection Q, adding a slope-two line exposes exactly
its previously uncovered cells and costs two. Distinct new lines expose
disjoint cell sets, so optimizing each choice separately proves the exact
Boolean compiler in equation (4). Removing a selected old line of degree
at most two never decreases c(Q). Hence an optimal old selection can avoid
all such lines; its covered cells are disjoint from U. This gives

    tau_3 >= tau_LP >= beta_3 >= tau_2+J_3.

The fractional inequality follows because any feasible fractional retained
set has total retained mass at most 2|Q| in N(Q). It does not assume
integrality of the three-family hypergraph.

For the displayed three-triple triangle, every pair intersects but no cell
lies in all three constraints. Thus two integer deletions are necessary
and sufficient. Half a deletion at each shared corner is feasible; half
weight on each triple gives the matching fractional lower bound 3/2. The
Boolean values of one, two, and three overfull lines are 1,1,0; low-degree
line vertices can again be removed without decreasing the objective. Thus
beta_3=1<3/2<tau_3=2. The determinant-minus-two minor records the precise
failed integrality mechanism. The n<=4 geometric dimension boundary is
correct because a slope-two grid line then contains at most two cells.

## 2. Eligibility, forced companions and uniform avoidance

Fix any simple bipartite two-regular skeleton and any physical row order;
only the column labels are random. For targets (r,q+2r), their difference
and sum labels are q+r and q+3r. Distinct target rows therefore have
distinct labels within each old family and distinct target columns.

Each row has two source-column neighbours and each such column has a
different companion row v. A selected row r excludes at most the six
other-row addresses v,2r-v,(2r+v)/3 over its two neighbours. For a line
with L>=3 available rows, choose the offending directed pair and then the
third row: at most 6L(L-2) unordered triples can fail eligibility, with
overcounting harmless. There are O(n) lines and L=O(n), so the discarded
count is O(n^3) with an absolute constant independent of the skeleton or
row order.

Eligibility forbids source-column sharing between selected rows. Thus
every one of the eight intended source-edge choices uses three distinct
columns and has assignment probability exactly 1/(n)_3. A forced companion
(v,q+2r) lies on the target slope-two line only if v=r, on another target
d-line only if its row is 2r-v, and on another target s-line only if its
row is (2r+v)/3. The same-row cases also force v=r. All are excluded, so
every forced companion lies outside the complete seven-line union.

After these assignments, the residual permutation is uniform on n-3
columns and labels. Each free source column has two row neighbours, each
meeting at most seven target lines; its forbidden degree is at most 14.
A physical column likewise meets at most seven target rows, each having
two source-column neighbours, giving column degree at most 14. Each grid
cell on the line union accounts for at most two forbidden entries. Hence
m<=2S<=13n+1, including all overlaps and discarded assigned labels only as
harmless overcounts.

The inherited bounded-degree permutation-avoidance proof applies directly
with fixed Delta=14 and N=n-3. Its factorial-moment collision estimate and
odd/even Bonferroni truncations give exp(-m/N)+o(1), uniformly in all these
matrices. Since S=O(n) uniformly,

    exp(-m/(n-3))+o(1) >= exp(-2S/n)-o(1).

The change of denominator contributes only O(1/n). No independence of
matrix entries or occupied lines has entered this deduction.

Avoidance gives exactly the three target cells on q and singleton old
diagonals through each. The eight assignments for one triple are mutually
exclusive because a target label cannot be assigned to two different
source columns. Different triples on the same q-line are mutually
exclusive because the event prescribes exactly its three occupied cells.
Events on different lines need not be independent: their expectations are
simply added. The eligibility condition is sufficient, not necessary for
a realized isolated triple; the report retains this distinction.

## 3. Triple normalization and the Jensen coefficient

For a triple with extreme row distance h, there are n-h choices for the
first row, h-1 for the interior row and n-2h for the first column. This
proves the exact potential-triple formula and its n^4/32 leading term.
Deleting O(n^3) ineligible triples does not change this coefficient.

The scaled available row interval at q is
[max(0,-q/2), min(1,(1-q)/2)]. It has the stated L(q). The opposite
diagonal lengths through a point are 1-|r+q| and 1-|3r+q-1|, giving the
stated d(r,q). Splitting the row integral at these two absolute-value
zeros validates the three polynomial pieces of H(q). Ordinary polynomial
integration gives the reported values

    integral L^3 = 3/16,
    integral L^4 = 7/80,
    integral L^2 H = 187/720.

The factor three in the mean is essential. The q-line length occurs once
per triple, while each target row's two diagonal lengths occur in exactly
binom(L_q-1,2) triples. Dividing by the total triple count therefore gives

    mean(2S/n) ->
    2 [integral L^4 + 3 integral L^2 H] / integral L^3
    = 416/45.

All costs are at most 13+O(1/n), so removing the uniform O(n^3) exceptional
set changes this average by O(1/n). Convexity of exp(-a) applies to this
deterministic list of geometric costs. It therefore gives exp(-416/45)
as a uniform lower limit for the average avoidance envelope. Finally,
eight assignments times n^4/32 triples divided by (n)_3 contributes n/4.
The uniform avoidance error contributes o(n), proving

    liminf_n inf_(G,rho) E_sigma I_3/n >= exp(-416/45)/4.

This is an analytic asymptotic conclusion, not an extrapolation from the
finite event bank. Combining the deterministic inequality with the
inherited uniform conditional mean yields gamma+delta_3. No finite
threshold or optimality assertion follows or is claimed.

## 4. Concentration and adaptive output scope

For equal-size boards differing by s original points, intersecting an
optimal retained subset with the other board loses at most s points.
Reversing this comparison proves |tau_3(B)-tau_3(B')|<=s. Feasibility in
three families is hereditary, which is all this proof needs. A row or
column transposition changes at most four original points of a two-regular
board. Exposing n-1 column labels and coupling completions by one
transposition gives conditional increment ranges at most four. The
bounded-range exponential argument therefore yields the exact displayed
denominator 8(n-1).

The expectation used in this tail bound is conditional on the fixed
skeleton and row order. Its uniform asymptotic lower bound gives the
rate -(gamma+delta_3-kappa)^2/8 for any fixed nonnegative loss fraction
kappa below gamma+delta_3. Uniformity also permits the inherited uncolored
saturated-board mixture, where the columns remain conditionally uniform.
An arbitrary changed starting distribution is not covered automatically.

For any final no-three-in-line set T, the intersection B with T meets all
three selected line constraints. Hence |B minus T|>=tau_3(B) pointwise,
regardless of adaptive deletions, insertions, moves, or auxiliary
randomness. This is an obstruction to retaining too many original cells.
It is neither a running-time bound for unrestricted repair nor an
extremal nonexistence theorem. The report preserves both limitations.

The analytic audit is complete. The separate finite census, named boards,
residual-permutation event bank and execution pins remain the responsibility
of their frozen producer and the root's independent finite controls.
