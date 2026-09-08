# Independent audit: a third direction forces extra original-cell repair

**Verdict: ACCEPTED analytic theorem and FINITE-EXACT controls.** The
deterministic increment, uniform asymptotic increment and concentration
apply to every fixed simple two-regular row/column skeleton with every
fixed row order, under a uniformly random column permutation. This
does not prove an extremal no-three-in-line nonexistence theorem.

The producer is
[the three-direction report](continuing11_20260908_no3_three_directions.md).
The previous two-direction expectation constant and uniform permutation
avoidance lemma are proved in
[the repair report](continuing10_20260907_no3_repair.md). This audit
rederives the new argument and uses no producer code.

## 1. The actual additive cost and the lost integrality

For a set B of original cells let tau2 be the minimum deletion count
leaving at most two cells on each diagonal of slopes +1 and -1. Define
tau3 similarly after also imposing slope2, whose lines are c-2r=q.
Let U consist of cells whose two original opposite diagonals each
have occupancy at most two. Put

    J3=sum_q (|U intersect {c-2r=q}|-2)_+.

If D is any three-direction deletion set, D minus U still clears both
old directions: no originally overfull old line contains a cell of U.
Thus |D minus U|>=tau2. Every slope2 line needs at least its displayed
number of U deletions, independently of deletions outside U. These lines
partition U, so |D intersect U|>=J3. Therefore tau3>=tau2+J3.
The argument needs the ORIGINAL old occupancies; recomputing a safe set
after an optimal old repair would allow overlapping charges.

For any selected collection Q of lines, let N(Q) be its union of original
cells and c(Q)=|N(Q)|-2|Q|. Any feasible retained set has at most 2|Q|
points in N(Q), giving a Boolean dual lower bound. For two directions its
maximum is tau2 by the inherited unit-cell-capacity bipartite flow theorem.
For three directions write beta3 for the maximum over all three families.
Fixing the selected old lines, the slope2 lines have disjoint contributions
outside their union. Optimizing each independently gives the exact identity

    beta3=max_Qold [c(Qold)+sum_q (|q minus N(Qold)|-2)_+].

Old lines of size at most two never improve the objective, so an old
maximizer may use only overfull lines. Its union avoids U. Consequently

    tau3 >= beta3 >= tau2+J3.

The first inequality need not be equality. On the n=5 board with zero-based
row pairs (03,12,04,34,12), the three overfull triples form an actual
triangle: beta3=1, the deletion LP has optimum3/2, and tau3=2. Giving
weight1/2 to cell indices0,5,6 (row-major order) is primal feasible;
giving weight1/2 to each of the three line constraints is dual feasible.
Both costs are3/2. The corresponding cell-incidence minor has determinant
of absolute value2. This identifies the failed implication: the third
direction destroys the integral bipartite structure, even for a genuine
two-regular board. The compiled Boolean dual remains a valid lower bound.

## 2. Uniform eligibility of isolated triples

Let I3 count slope2 lines with exactly three original cells, each of whose
two old diagonals is a singleton. Then I3<=J3. Fix three distinct target
rows r_i on a line c=2r+q. Each row has two source-column choices.
Exclude target pairs sharing a source column. If the other row incident to
a chosen source column at r is v, its forced cell will be (v,2r+q).
It cannot lie on the target slope2 line, since v differs from r. It hits
another chosen cell's old diagonal only if that row is

    2r-v, or (2r+v)/3.

Exclude these pairs for BOTH source-column choices at every target row.
Each row has at most two companion rows and at most four further forbidden
positions. Thus there are O(n) forbidden ordered pairs, O(n^2) triples
per line and O(n^3) excluded triples over all slope2 lines, uniformly in
the fixed skeleton and row order. All eight assignments on an eligible
triple are valid, use distinct source columns and have probability1/(n)_3.
They are disjoint because the singleton requirement determines the unique
original source column mapping to each chosen target column.

Condition on one assignment. The remaining permutation has size N=n-3.
Forbid every remaining original cell on the target line and the six old
diagonals of the triple. The fixed companion cells already avoid them.
Each remaining source column has two rows, so the forbidden matrix has
row degree at most14. For each target column each of the seven lines has
at most one row, giving column degree at most14. If M is its number of
forbidden entries, then M is at most twice the sum of the seven physical
line lengths. Removing fixed source/target columns only decreases M.

The inherited bounded-degree rook-count/Bonferroni lemma gives avoidance
probability exp(-M/N)+o(1), uniformly for all these matrices. Its hypotheses
are paid independently of the skeleton's cycle lengths or row order.
Multiplying its uniform error by O(n^4)/(n)_3 gives o(n), so the ensuing
expectation bound is uniform, not merely true for a typical row order.

## 3. Exact seven-line average and the Jensen gain

Scale the square to [0,1]^2 and q to [-2,1]. Let I_q be the allowed row
interval for c=2r+q and L(q) its length. The sum of the two old diagonal
lengths through (r,2r+q) is

    d(r,q)=2-|r+q|-|3r+q-1|.

Writing J(q)=integral over I_q of d(r,q), direct integration gives

| q interval | L(q) | J(q) |
|---|---|---|
| [-2,-1] | (q+2)/2 | (q+2)(q+5)/6 |
| [-1,0] | 1/2 | -(q-1)(q+2)/3 |
| [0,1] | (1-q)/2 | (q-4)(q-1)/6 |

The independent code splits the crossing of absolute-value breakpoints
at q=-1/2 before integration. It obtains exactly

    integral L^3=3/16, integral L^4=7/80,
    integral L^2 J=187/720.

The count of all target triples is asymptotic to n^4 times
integral L^3/6=n^4/32. Under this triple-weighted distribution the mean
upper bound for M/N tends to

    2 (integral L^4+3 integral L^2 J)/(integral L^3)=416/45.

All upper costs are bounded by13+o(1); removing the O(n^3) ineligible
triples does not change either limiting average or leading count.
Jensen's inequality for the convex exponential now yields

    liminf E I3/n >= delta3=exp(-416/45)/4.

The factor1/4 is8/32: all eight source choices are included once.
This improves the crude maximum-cost bound exp(-13)/4. With the inherited
two-direction constant gamma2, linearity of expectation and the
deterministic inequality prove

    liminf E tau3/n >= gamma2+delta3
                        =0.3286222170740135... .

The proof is a lower bound; it does not claim that Jensen's inequality
is sharp or that the two random deletion costs are independent.

## 4. Adaptive scope, tail and independent universe

For any no-three-in-line set T, B intersect T is a three-direction feasible
subset of the original board B. Hence |B minus T|>=tau3(B), even when T
depends on B and arbitrary new cells are added. A column transposition
removes at most four original cells and adds at most four. The maximum
retained feasible subset size changes by at most four, and therefore so
does tau3. The permutation exposure martingale has conditional range
at most four at each of n-1 steps, giving

    P(tau3<=k) <= exp(-(E tau3-k)_+^2/[8(n-1)]).

This gives an exponentially small probability for sublinear original-cell
repair, uniformly for every fixed row order. An exp(-c n) bound does not
pay a union bound over n! row orders or all skeletons, and does not settle
the extremal no-three-in-line problem.

The independent executable enumerates every row-pair product and filters
column degrees, obtaining all2137 boards for n=2..5 without the producer's
recursive deficit pruning. It solves each deletion problem by increasing
cardinality subset search, separately computes both Boolean duals, and
checks the complete inequality chain. There are142 boards with strict
tau3>tau2 at n=5 and no isolated event below n=6. Five named controls include
the n=6 board (01,13,05,23,24,45), with (tau2,tau3,J3,I3)=(3,4,1,1),
and a board refuting naive addition of full slope2 excess. No complete
n=6 census is claimed. Normal and optimized runs pass4307 always-active
gates with identical raw LF output and regenerated certificate bytes.
